import openmm.app as app
import openmm as mm
import openmm.unit as unit
import numpy as np

from oqp.openqp import OPENQP
import oqp


def unpack_lower_tri_single(packed_atom, nbf):
    """
    Unpack a single lower-triangular packed array (length nbf*(nbf+1)/2)
    into a full symmetric (nbf x nbf) matrix.
    """
    packed_atom = np.asarray(packed_atom)
    nbf_tri = nbf * (nbf + 1) // 2
    if packed_atom.size != nbf_tri:
        raise ValueError(f"Size mismatch: got {packed_atom.size}, expected {nbf_tri}")
    full = np.zeros((nbf, nbf), dtype=packed_atom.dtype)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            val = packed_atom[idx]
            full[i, j] = val
            full[j, i] = val
            idx += 1
    return full

def unpack_lower_tri_multi(packed, nbf, natm):
    """
    packed: Matrix in packed lower form for all atoms.
            Layout: natm blocks, each of length nbf*(nbf+1)/2.
    Returns: full[natm, nbf, nbf] – one full matrix per atom.
    """
    packed = np.asarray(packed)
    nbf_tri = nbf * (nbf + 1) // 2
    flat = packed.ravel(order="C")
    if flat.size != natm * nbf_tri:
        raise ValueError(f"Size mismatch: got {flat.size}, expected {natm * nbf_tri}")
    # Shape: (natm, nbf_tri), each row = packed lower-tri for one atom
    packed_by_atom = flat.reshape((natm, nbf_tri))
    full_all = np.zeros((natm, nbf, nbf), dtype=packed.dtype)
    for a in range(natm):
        full_all[a] = unpack_lower_tri_single(packed_by_atom[a], nbf)

    return full_all

def pack_lower_tri_single(full):
    """
    Pack a single symmetric (nbf x nbf) matrix into its lower-triangular
    packed form (length nbf*(nbf+1)/2).
    """
    full = np.asarray(full)
    if full.shape[0] != full.shape[1]:
        raise ValueError("Matrix must be square")

    nbf = full.shape[0]
    nbf_tri = nbf * (nbf + 1) // 2
    packed = np.zeros(nbf_tri, dtype=full.dtype)

    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            packed[idx] = full[i, j]
            idx += 1

    return packed


class OpenQpQMMM:
    """
    High-level QM/MM driver using OpenMM (MM) + OpenQP (QM).

    Usage
    -----
    qmmm = OpenQpQMMM(
        positions=init_positions,
        topology=topology,
        forcefield=forcefield,
        qm_atoms=[0, 1, 2],     # QM atom indices
        oqp_cfg=oqp_cfg_dict,   # OpenQP configuration dict
        Cutoff=app.NoCutoff,    # or app.PME / app.CutoffNonPeriodic
        Embedding='electrostatic'  # or 'mechanical'
    )

    total_E, total_F = qmmm.compute(current_positions)
    # or simply:
    total_E, total_F = qmmm(current_positions)
    """

    def __init__(
        self,
        positions,
        topology,
        forcefield,
        qm_atoms,
        oqp_cfg,
        Cutoff=app.NoCutoff,
        Embedding='mechanical',
    ):
        # Store basic stuff
        self.positions = positions
        self.topology = topology
        self.forcefield = forcefield
        self.qm_atoms = np.array(qm_atoms, dtype=int)
        self.Cutoff = Cutoff
        self.Embedding = Embedding

        # Store *base* OpenQP configuration (user-defined)
        # We will inject "input.system" at each call
        self.oqp_cfg_base = oqp_cfg  # shallow copy is enough

        # Build MM systems (periodic, Ewald, original, etc.)
        self.prepare_mm()

    # --- Internal helpers -------------------------------------------------

    def _build_xyz_string(self):
        """
        Build the 'atom string' for OpenQP:
        'H 0.000000 0.000000 0.000000; O 0.958000 0.000000 0.000000; ...'
        """
        xyz_atoms = []
        for atom in self.topology.atoms():
            at_index = atom.index
            if at_index in self.qm_atoms:
                sym = atom.element.symbol
                x = self.positions[at_index][0].value_in_unit(unit.angstrom)
                y = self.positions[at_index][1].value_in_unit(unit.angstrom)
                z = self.positions[at_index][2].value_in_unit(unit.angstrom)
                xyz_atoms.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
        return '; '.join(xyz_atoms)

    def forces_qm_openqp(self, potmm=None, potqm=None):
        """
        QM energy + forces with OpenQP, using self.oqp_cfg_base and
        the current positions. This is basically your forces_qm()
        but driven by a user configuration dict.
        """

        # Build system string
        xyz_atoms = self._build_xyz_string()

        # --- OpenQP setup --------------------------------------------------
        self.oqp_cfg_base["input.system"] = xyz_atoms
        op = OPENQP(self.oqp_cfg_base, True)
        op.sp._prep_guess()

        op.mol.data["OQP::POTMM"] = potmm
        op.mol.data["OQP::POTQM"] = potqm
        oqp.espf_op_corr(op.mol)
        espf_op_corr = op.mol.data["OQP::ESPF_CORR"]

        # --- Modify 1e Hamiltonian with ESPF term -------------------------
        basis = op.mol.data.get_basis()
        nat = op.mol.data["natom"]
        nbf = basis["nbf"]

        espf_op_corr_f = unpack_lower_tri_multi(espf_op_corr, nbf, nat)

        if potmm is not None:
            hcore = op.mol.get_hcore()
            hcore_full = unpack_lower_tri_single(hcore, nbf)
            # add ESPF contribution
            hcore_full += np.einsum("ijk,i->jk", espf_op_corr_f, potmm)
            hcore = pack_lower_tri_single(hcore_full)
            op.mol.set_hcore(hcore)

        # --- QM energy & charges ------------------------------------------
        op.sp.scf()
        self.eqm = op.mol.get_scf_energy()

        oqp.form_esp_charges(op.mol)
        self.pchg_qm = op.mol.data["OQP::partial_charges"]

        # Correct potmm to add QM–QM interactions
        if potqm is not None and potmm is not None:
            # Make a local copy to avoid in-place surprises
            potmm -= np.einsum(
                "ij,j->i", potqm, self.pchg_qm - op.mol.get_atoms("charge")
            )
            op.mol.data["OQP::POTMM"] = potmm

        # Subtract QM/MM energy from eqm (OpenMM will add it later)
        if potmm is not None:
            self.eqm -= np.dot(
                self.pchg_qm - op.mol.get_atoms("charge"),
                potmm
            )

        # --- Gradients: pure QM + ESPF contribution -----------------------
        oqp.hf_gradient(op.mol)
        grad = op.mol.get_grad()
        gqm = np.array(grad.copy()).reshape(
            (1, op.mol.get_atoms("natom"), 3)
        )

        oqp.grad_esp_qmmm(op.mol)
        gqm += op.mol.data["OQP::ESPF_GRAD"]

        # --- Unit conversion to OpenMM conventions ------------------------
        self.eqm *= 2625.5 * unit.kilojoule_per_mole
        self.gqm = gqm[0]*49578.9  # already dimensionless, factor to kJ/mol/nm

        return self.eqm, self.gqm, self.pchg_qm

    #----------------------------------------------------------------------
    #-----------         MM FORCE (INCLUDING CLASSIC QM/MM)       ---------
    #----------------------------------------------------------------------
    def forces_mm(self, pchg_qm):

        system = self.mm_systems["sys0"]
        simulation = self.mm_systems["sim0"]

        # Replace QM partial charges
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        nonbonded = forces['NonbondedForce']
        for i in range(nonbonded.getNumParticles()):
            if i in self.qm_atoms:
               charge, sigma, epsilon = nonbonded.getParticleParameters(i)
               charge = pchg_qm[np.where(self.qm_atoms == i)[0][0]]*unit.elementary_charge
               nonbonded.setParticleParameters(i, charge, sigma, epsilon)
        for i in range(nonbonded.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = nonbonded.getExceptionParameters(i)
            if p1 in self.qm_atoms and p2 not in self.qm_atoms:
               charge1 = pchg_qm[np.where(self.qm_atoms == p1)[0][0]]*unit.elementary_charge
               charge2, sigma2, epsilon2 = nonbonded.getParticleParameters(p2)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
            elif p2 in self.qm_atoms and p1 not in self.qm_atoms:
               charge2 = pchg_qm[np.where(self.qm_atoms == p2)[0][0]]*unit.elementary_charge
               charge1, sigma1, epsilon1 = nonbonded.getParticleParameters(p1)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
        nonbonded.updateParametersInContext(simulation.context)

        state = simulation.context.getState(getEnergy=True,getForces=True)

        emm = state.getPotentialEnergy()
        gmm = state.getForces(asNumpy=True)

        return emm, gmm


    # --- Public API -------------------------------------------------------
    def compute_force(self, positions, topology,mm_systems,qm_atoms):

        """
        Compute total QM/MM energy and forces for a given set of positions.

        Parameters
        ----------
        positions : openmm.Vec3 array with units
            Current atomic coordinates.

        Returns
        -------
        total_energy : openmm.unit.Quantity
            Total QM + MM energy (same units as OpenMM, kJ/mol).
        total_forces : np.ndarray (n_atoms, 3) with units of kJ/mol/nm
            Total forces on all atoms (QM + MM combined).
        """
        self.positions = positions
        self.topology = topology
        self.mm_systems = mm_systems
        self.qm_atoms = qm_atoms

        # 1. Electrostatic potential (if needed)
        potmm = potqm = None
        if self.Embedding == "electrostatic":
            potmm, potqm = self.electrostatic_potential()

        # 2. QM part with OpenQP
        eqm, gqm, pchg_qm = self.forces_qm_openqp(
            potmm=potmm,
            potqm=potqm
        )

        # 3. MM+classic QM/MM with updated QM charges
        emm, gmm = self.forces_mm(pchg_qm)

        # 4. Combine energies and forces
        total_energy = eqm + emm
        total_forces = gmm.copy()

        # subtract QM gradient from the corresponding atoms
        for k, i in enumerate(self.qm_atoms):
            total_forces[i] = total_forces[i] - gqm[k]

        return total_energy, total_forces

    def prepare_mm(self):

       system=self.forcefield.createSystem(self.topology,nonbondedMethod=self.Cutoff,constraints=None,rigidWater=False)
       int0=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
       simulation=app.Simulation(self.topology, system, int0)
       simulation.context.setPositions(self.positions)

       #Deactivate QM-QM interactions
       for f in system.getForces():
       #Deactivating non-bonded terms between QM atoms (but keep QM-MM)
          if isinstance(f, mm.NonbondedForce):
             for p1 in self.qm_atoms:
                for p2 in self.qm_atoms:
                   if p1 != p2: f.addException(p1,p2,0,0,0,replace=True)
             f.updateParametersInContext(simulation.context)
          elif isinstance(f, mm.CustomNonbondedForce):
             for p1 in self.qm_atoms:
               for p2 in self.qm_atoms:
                  if p1 != p2: f.addExclusion(p1,p2)
             f.updateParametersInContext(simulation.context)
       #Deactivating QM bonded terms
          elif isinstance(f, mm.HarmonicBondForce):
             for i in range(f.getNumBonds()):
                p1, p2, length, k = f.getBondParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms)
                if exclude: f.setBondParameters(i, p1, p2, length, 0)
             f.updateParametersInContext(simulation.context)
          elif isinstance(f, mm.CustomBondForce):
             for i in range(f.getNumBonds()):
                p1, p2, parameters = f.getBondParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms)
                if exclude: f.setBondParameters(i, p1, p2, (0,0))
       #Deactivating QM angle terms
          elif isinstance(f, mm.HarmonicAngleForce):
             for i in range(f.getNumAngles()):
                p1, p2, p3, angle, k = f.getAngleParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms and p3 in self.qm_atoms)
                if exclude: f.setAngleParameters(i, p1, p2, p3, angle, 0)
             f.updateParametersInContext(simulation.context)
          elif isinstance(f, mm.CustomAngleForce):
             for i in range(f.getNumAngles()):
                p1, p2, p3, parameters = f.getAngleParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms and p3 in self.qm_atoms)
                if exclude: f.setAngleParameters(i, p1, p2, p3, (0, 0))
             f.updateParametersInContext(simulation.context)
       #Deactivating QM torsion terms
          elif isinstance(f, mm.PeriodicTorsionForce):
             for i in range(f.getNumTorsions()):
                p1, p2, p3, p4, periodicity, phase, k = f.getTorsionParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms and p3 in self.qm_atoms and p4 in self.qm_atoms)
                if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, 0)
             f.updateParametersInContext(simulation.context)
          elif isinstance(f, mm.CustomTorsionForce):
             for i in range(f.getNumTorsions()):
                p1, p2, p3, p4, parameters = f.getTorsionParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms and p3 in self.qm_atoms and p4 in self.qm_atoms)
                if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, (0,0,0))
             f.updateParametersInContext(simulation.context)
          elif isinstance(f, (mm.CMAPTorsionForce)):
             for i in range(f.getNumTorsions()):
                cmap, p1, p2, p3, p4, q1, q2, q3, q4 = f.getTorsionParameters(i)
                exclude = (p1 in self.qm_atoms and p2 in self.qm_atoms and p3 in self.qm_atoms and p4 in self.qm_atoms)
                exclude = exclude and (q1 in self.qm_atoms and q2 in self.qm_atoms and q3 in self.qm_atoms and q4 in self.qm_atoms)
                if exclude: f.setMapParameters(i,cmap.size,0)
             if f.getNumTorsions() != 0: f.updateParametersInContext(simulation.context)
       #Exception, unless CMMotionRemover
          else:
             if not isinstance(f, mm.CMMotionRemover): exit(f"Force not found")

       if self.Cutoff is not app.NoCutoff:
          sysew=self.forcefield.createSystem(self.topology,nonbondedMethod=app.Ewald,constraints=None,rigidWater=False)
          intew=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
          simew=app.Simulation(self.topology, sysew, intew)
          simew.context.setPositions(self.positions)

          sysor=self.forcefield.createSystem(self.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
          intor=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
          simor=app.Simulation(self.topology, sysor, intor)
          simor.context.setPositions(self.positions)
       else:
          sysew = simew = sysor = simor = None

       self.mm_systems = {
       "sys0": system,
       "sim0": simulation,
       "sysew": sysew,
       "simew": simew,
       "sysor": sysor,
       "simor": simor,
       }

#Parameters to choose
    def electrostatic_potential(self):

       syspbc=self.mm_systems["sys0"]
       simpbc=self.mm_systems["sim0"]

    #Focus on non-bonded interactions
       forces = { force.__class__.__name__ : force for force in syspbc.getForces() }
       nonbonded = forces['NonbondedForce']
       if nonbonded is None: ValueError(f"Non-bonded interactions are not present, what shall I do?")

    #######################################################################
    #  1. Compute MM potential (needs QM-QM contributions to be removed)  #
    #######################################################################
       potmm=np.zeros((len(self.qm_atoms)))
       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)
       nonbonded.updateParametersInContext(simpbc.context)
       state = simpbc.context.getState(getEnergy=True)
       e_pbc_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, sigma, epsilon)
           nonbonded.updateParametersInContext(simpbc.context)
           state = simpbc.context.getState(getEnergy=True)
           e_pbc_qm_charge=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potmm[i]=e_pbc_qm_charge-e_pbc_no_qm_charge
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)

       if self.Cutoff == app.NoCutoff: return potmm,None

    #######################################################################
    #                 2. Compute QM pair potential                        #
    #######################################################################
       potqm = np.zeros((len(self.qm_atoms), len(self.qm_atoms)))

    # Create an Ewald system
       sysew=self.mm_systems["sysew"]
       simew=self.mm_systems["simew"]
       forces = { force.__class__.__name__ : force for force in sysew.getForces() }
       nonbondedew = forces['NonbondedForce']

    # Create a non-periodic system
       sysor=self.mm_systems["sysor"]
       simor=self.mm_systems["simor"]
       forcesor = { force.__class__.__name__ : force for force in sysor.getForces() }
       nonbondedor = forcesor['NonbondedForce']

    # Create a system with no charges, both QM and MM (Ewald & Original)
       for i in range(sysew.getNumParticles()):
           nonbondedew.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
       nonbondedew.updateParametersInContext(simew.context)
       nonbondedor.updateParametersInContext(simor.context)

       state = simew.context.getState(getEnergy=True)
       e_ew_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       state = simor.context.getState(getEnergy=True)
       e_or_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

    #######################################################################
    # 2.1. Compute the QM-QM diagonal potential and remove QM-QM from MM  #
    #         Note 1: the 1/2 factor needs to be corrected later          #
    #         Note 2: here the PME/Ew method			      #
    #######################################################################
       for i in range(len(self.qm_atoms)):
           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedew.updateParametersInContext(simew.context)
           state = simew.context.getState(getEnergy=True)
           e_ew_qm_chargei = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potqm[i,i] = e_ew_qm_chargei - e_ew_no_qm_charge
           potmm[i] -= potqm[i,i] #Remove QM-QM interactions from MM potential
           potqm[i,i] *= 2.0
           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

    #######################################################################
    #        2.2. Compute the QM-QM off-diagonal potential                #
    #######################################################################
       for i in range(len(self.qm_atoms)):

           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)

           for j in range(i+1,len(self.qm_atoms)):

               nonbondedew.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedew.updateParametersInContext(simew.context)
               state = simew.context.getState(getEnergy=True)
               e_ew_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

               nonbondedor.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.updateParametersInContext(simor.context)
               state = simor.context.getState(getEnergy=True)
               e_or_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
               ecorr = e_or_qm_chargeij - e_or_no_qm_charge

               potij = e_ew_qm_chargeij - e_ew_no_qm_charge - ecorr - 0.5*(potqm[j,j] + potqm[i,i])
               potqm[i,j] = potqm[j,i] = potij

               nonbondedew.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)

           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

       return potmm,potqm


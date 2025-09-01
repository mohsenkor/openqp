module trah
  use precision, only: dp
  character(len=*), parameter :: module_name = "trah"
  public :: calc_g_k
  public :: calc_h_op
  public :: compute_energy
  public :: rotate_orbs_trah

contains

  subroutine calc_g_h(grad, h_diag, fock_ao, mo_a, nbf, nocc_a, scftype, mo_b, nocc_b)
!    use types, only: information
    use mathlib, only: pack_matrix, unpack_matrix
    implicit none
!    type(information), intent(inout) :: infos
    integer, intent(in) :: scftype, nbf, nocc_a
    integer, intent(in), optional :: nocc_b
    real(dp), intent(out) :: grad(:)
    real(dp), intent(out) :: h_diag(:)
    real(dp), intent(in)  :: fock_ao(:,:)      ! packed AO fock
    real(dp), intent(inout)  :: mo_a(nbf,nbf)     ! MO coefficient matrix
    real(dp), intent(inout), optional  :: mo_b(nbf,nbf)     ! MO coefficient matrix

    integer :: nvir, i, a, k
    real(dp), allocatable :: work1(:,:), work2(:,:), work3(:,:)

    ! dimensions
    nvir = nbf - nocc_a

    allocate(work1(nbf,nbf))
    allocate(work2(nbf,nbf))
    allocate(work3(nbf,nbf))

!    if (size(grad)   /= nocc_a*nvir) stop 'grad wrong size'
!    if (size(h_diag)/= nocc_a*nvir) stop 'h_diag wrong size'
    select case(scftype)
    case (1)
      ! unpack the AO fock into a full (nbf×nbf) matrix
      call unpack_matrix(fock_ao(:,1), work1)
      work2 = 0.0_dp
      work3 = 0.0_dp
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo_a, nbf, &
                 0.0_dp, work2, nbf)

      ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo_a, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)

      k = 0
      do i = nocc_a +1, nbf
        do a = 1, nocc_a
          k= k+1
          grad(k) =  2 * work3(i,a)
        end do
      end do

      k = 0
      do i = nocc_a +1, nbf
        do a = 1, nocc_a
          k= k+1
          h_diag(k) = 2.0_dp*( work3(i,i) - work3(a,a) )
        end do
      end do
! UHF
    case(2)
      ! unpack the AO Beta fock matrix
      call unpack_matrix(fock_ao(:,1), work1)
      work2 = 0.0_dp
      work3 = 0.0_dp
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo_a, nbf, &
                 0.0_dp, work2, nbf)

      ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo_a, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)

      k = 0
      do i = nocc_a +1, nbf
        do a = 1, nocc_a
          k= k+1
          grad(k) =  2 * work3(i,a)
        end do
      end do

      k = 0
      do i = nocc_a +1, nbf
        do a = 1, nocc_a
          k= k+1
          h_diag(k) = 2.0_dp*( work3(i,i) - work3(a,a) )
        end do
      end do
      ! unpack the Beta fock matrix
      call unpack_matrix(fock_ao(:,2), work1)
      work2 = 0.0_dp
      work3 = 0.0_dp
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo_b, nbf, &
                 0.0_dp, work2, nbf)

      ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo_b, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)

      k = (nbf-nocc_a) * nocc_a
      do i = nocc_b +1, nbf
        do a = 1, nocc_b
          k= k+1
          grad(k) =  2 * work3(i,a)
        end do
      end do

      k = (nbf-nocc_a) * nocc_a
      do i = nocc_b +1, nbf
        do a = 1, nocc_b
          k= k+1
          h_diag(k) = 2.0_dp*( work3(i,i) - work3(a,a) )
        end do
      end do
    end select
    deallocate(work1, work2, work3)

  end subroutine calc_g_h

   !> @brief Computes the two-electron part (Coulomb and exchange) of the Fock matrix.
  !> @detail Forms the Coulomb (J) and exchange (K) contributions to the Fock matrix
  !>         using two-electron integrals,
  !>         with optional scaling of the exchange term for hybrid DFT methods.
  !> @param[in] basis Basis set information.
  !> @param[in] d Density matrices (triangular format).
  !> @param[inout] f Fock matrices to be updated (triangular format).
  !> @param[in] scalefactor Optional scaling factor for exchange (default = 1.0).
  !> @param[inout] infos System information.
  subroutine fock_jk(basis, d, f, scalefactor, infos)
    use precision, only: dp
    use io_constants, only: iw
    use util, only: measure_time
    use basis_tools, only: basis_set
    use types, only: information
    use int2_compute, only: int2_compute_t, int2_fock_data_t, &
                            int2_rhf_data_t, int2_urohf_data_t

    implicit none

    type(basis_set), intent(in) :: basis
    type(information), intent(inout) :: infos
    real(kind=dp), optional, intent(in) :: scalefactor

    integer :: i, ii, nf, nschwz
    real(kind=dp) :: scalef
    real(kind=dp), target, intent(in) :: d(:,:)
    real(kind=dp), intent(inout) :: f(:,:)

    type(int2_compute_t) :: int2_driver
    class(int2_fock_data_t), allocatable :: int2_data

    ! Initial Settings
    scalef = 1.0d0
    if (present(scalefactor)) scalef = scalefactor



    ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_data = int2_rhf_data_t(nfocks=1, d=d, scale_exchange=scalefactor)


    ! Constructing two electron Fock matrix
    call int2_driver%run(int2_data)
    nschwz = int2_driver%skipped

    ! Scaling (everything except diagonal is halved)
    f =  0.5 * int2_data%f(:,:,1)
    do nf = 1, ubound(f,2)
      ii = 0
      do i = 1, basis%nbf
         ii = ii + i
         f(ii,nf) = 2*f(ii,nf)
      end do
    end do

    call int2_driver%clean()

  end subroutine fock_jk

  subroutine calc_h_op(infos,fock_ao, x, x2, mo, mo_b)
    use types, only: information
    use basis_tools,  only: basis_set
    use mathlib, only: pack_matrix, unpack_matrix
    implicit none

    class(information), intent(inout), target :: infos
    real(dp), intent(in)                 :: x(:)    ! length nocc*nvir
    real(dp), intent(out)                :: x2(:)   ! same length
    real(dp), pointer, intent(in)        :: fock_ao(:,:)
    real(dp), intent(inout)              :: mo(:, :)
    real(dp), optional, intent(inout)    :: mo_b(:,:)

    integer :: nbf, nocc_a,nocc_b, nvir_a, nvir_b, i, a, k
    real(dp), allocatable :: foo(:,:), fvv(:,:), xmat(:,:), x2mat(:,:),&
             foo_b(:,:), fvv_b(:,:), xmat_b(:,:),x2mat_b(:,:),&
             v(:,:), work2(:,:), work1(:,:), &
            work3(:,:)
    real(dp), allocatable :: pfock(:,:)
    real(dp), allocatable :: dm(:,:), dm_tri(:,:)
    real(kind=dp) :: scalefactor
    type(basis_set), pointer :: basis
    logical :: is_dft

    integer :: scf_type                 ! Type of SCF calculation
    integer, parameter :: scf_rhf  = 1  ! Restricted Hartree-Fock
    integer, parameter :: scf_uhf  = 2  ! Unrestricted Hartree-Fock
    integer, parameter :: scf_rohf = 3  ! Restricted Open-Shell Hartree-Fock
    integer:: nfocks

    ! dims
    basis => infos%basis
    nbf = basis%nbf
    nocc_a = infos%mol_prop%nelec_a
    nocc_b = infos%mol_prop%nelec_b
    nvir_a = nbf - nocc_a
    nvir_b = nbf - nocc_b

    is_dft = infos%control%hamilton >= 20
    if (is_dft) then
      scalefactor = infos%dft%HFscale
    else
      scalefactor = 1.0_dp
    end if
    basis => infos%basis

    allocate(foo(nocc_a,nocc_a), fvv(nvir_a,nvir_a))
    allocate(work3(nbf,nbf), work2(nbf,nbf), work1(nbf,nbf))
    allocate(xmat(nvir_a,nocc_a), x2mat(nvir_a,nocc_a))
    allocate(dm(nbf,nbf), v(nbf,nbf))
    allocate(pfock(nbf*(nbf+1)/2,2),dm_tri(nbf*(nbf+1)/2,2))

    select case (infos%control%scftype)
    case (1)
      scf_type = scf_rhf
      nfocks = 1
    case (2)
      scf_type = scf_uhf
      nfocks = 2
      allocate(foo_b(nocc_b,nocc_b), fvv_b(nvir_b,nvir_b))
      allocate(xmat_b(nvir_b,nocc_b), x2mat_b(nvir_a,nocc_b))

    case (3)
      scf_type = scf_rohf
      nfocks = 2
    end select
    select case (scf_type)
    case(scf_rhf)
      k = 0
      do i = 1, nvir_a
        do a = 1, nocc_a
          k= k+1
          xmat(i,a) = x(k)
        end do
      end do
      call unpack_matrix(fock_ao(:,1), work1)

      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo, nbf, &
                 0.0_dp, work2, nbf)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)
       foo = work3(1:nocc_a,1:nocc_a)
       fvv = work3(nocc_a+1:nbf,nocc_a+1:nbf)


      call dgemm('N','N', nvir_a, nocc_a, nvir_a, &
                 1.0_dp, fvv,  nvir_a, &
                         xmat, nvir_a, &
                 0.0_dp, x2mat, nvir_a)
      call dgemm('N','N', nvir_a, nocc_a, nocc_a, &
                -1.0_dp, xmat, nvir_a, &
                          foo, nocc_a, &
                 1.0_dp, x2mat, nvir_a)

      dm = 0.0_dp
      work2 = 0
      call dgemm('N','N', nbf, nocc_a, nvir_a, &
                 2.0_dp, mo(:, nocc_a+1:nbf), nbf, &
                         xmat,             nvir_a, &
                0.0_dp, work2,          nbf)
      call dgemm('N','T', nbf, nbf, nocc_a, &
                 1.0_dp, work2,          nbf, &
                         mo(:, 1:nocc_a),   nbf, &
                 0.0_dp, work3,          nbf)
      do i = 1, nbf
        do a = 1, nbf
          dm(i,a) = work3(i,a) + work3(a,i)
        end do
      end do
      call pack_matrix(dm,dm_tri(:,1))
      dm_tri(:,2) = dm_tri(:,1)
      call fock_jk(infos%basis, dm_tri, pfock, scalefactor, infos)
      call unpack_matrix(pfock(:,1), v)
      work2 = 0

      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo, nbf, &
                         v ,             nbf, &
                 0.0_dp, work2,          nbf)
      work3 = 0
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work2, nbf, &
                         mo, nbf, &
                 0.0_dp, work3,  nbf)
      x2mat = x2mat + work3(nocc_a+1:,1:nocc_a)

      k = 0
      do i = 1, nvir_a
        do a = 1, nocc_a
          k= k+1
          x2(k) = 2*x2mat(i,a)
        end do
      end do
      !#################### UHF
    case (scf_uhf)
! alpha
      k = 0
      do i = 1, nvir_a
        do a = 1, nocc_a
          k= k+1
          xmat(i,a) = x(k)
        end do
      end do
      call unpack_matrix(fock_ao(:,1), work1)

      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo, nbf, &
                 0.0_dp, work2, nbf)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)
       foo = work3(1:nocc_a,1:nocc_a)
       fvv = work3(nocc_a+1:nbf,nocc_a+1:nbf)


      call dgemm('N','N', nvir_a, nocc_a, nvir_a, &
                 1.0_dp, fvv,  nvir_a, &
                         xmat, nvir_a, &
                 0.0_dp, x2mat, nvir_a)
      call dgemm('N','N', nvir_a, nocc_a, nocc_a, &
                -1.0_dp, xmat, nvir_a, &
                          foo, nocc_a, &
                 1.0_dp, x2mat, nvir_a)

      dm = 0.0_dp
      work2 = 0
      call dgemm('N','N', nbf, nocc_a, nvir_a, &
                 2.0_dp, mo(:, nocc_a+1:nbf), nbf, &
                         xmat,             nvir_a, &
                0.0_dp, work2,          nbf)
      call dgemm('N','T', nbf, nbf, nocc_a, &
                 1.0_dp, work2,          nbf, &
                         mo(:, 1:nocc_a),   nbf, &
                 0.0_dp, work3,          nbf)
      do i = 1, nbf
        do a = 1, nbf
          dm(i,a) = work3(i,a) + work3(a,i)
        end do
      end do
      call pack_matrix(dm,dm_tri(:,1))

! beta
      k = nvir_a*nocc_a
      do i = 1, nvir_b
        do a = 1, nocc_b
          k= k+1
          xmat_b(i,a) = x(k)
        end do
      end do
      call unpack_matrix(fock_ao(:,2), work1)

      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work1, nbf,      &
                         mo_b, nbf, &
                 0.0_dp, work2, nbf)
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo_b, nbf, &
                         work2,         nbf, &
                 0.0_dp, work3,nbf)
      foo_b = work3(1:nocc_b,1:nocc_b)
      fvv_b = work3(nocc_b+1:nbf,nocc_b+1:nbf)

      call dgemm('N','N', nvir_b, nocc_b, nvir_b, &
                 1.0_dp, fvv_b,  nvir_b, &
                         xmat_b, nvir_b, &
                 0.0_dp, x2mat_b, nvir_b)
      call dgemm('N','N', nvir_b, nocc_b, nocc_b, &
                -1.0_dp, xmat_b, nvir_b, &
                          foo_b, nocc_b, &
                 1.0_dp, x2mat_b, nvir_b)

      dm = 0.0_dp
      work2 = 0
      call dgemm('N','N', nbf, nocc_b, nvir_b, &
                 2.0_dp, mo_b(:, nocc_b+1:nbf), nbf, &
                         xmat_b,             nvir_b, &
                0.0_dp, work2,          nbf)
      call dgemm('N','T', nbf, nbf, nocc_b, &
                 1.0_dp, work2,          nbf, &
                         mo_b(:, 1:nocc_b),   nbf, &
                 0.0_dp, work3,          nbf)
      do i = 1, nbf
        do a = 1, nbf
          dm(i,a) = work3(i,a) + work3(a,i)
        end do
      end do

      call pack_matrix(dm,dm_tri(:,2))
! end of dm calculation
      call fock_jk(infos%basis, dm_tri, pfock, scalefactor, infos)
! alpha x2mat
      call unpack_matrix(pfock(:,1), v)
      work2 = 0
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo, nbf, &
                         v ,             nbf, &
                 0.0_dp, work2,          nbf)
      work3 = 0
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work2, nbf, &
                         mo, nbf, &
                 0.0_dp, work3,  nbf)
      x2mat = x2mat + work3(nocc_a+1:,1:nocc_a)

! beta x2mat
      call unpack_matrix(pfock(:,2), v)
      work2 = 0
      call dgemm('T','N', nbf, nbf, nbf, &
                 1.0_dp, mo, nbf, &
                         v ,             nbf, &
                 0.0_dp, work2,          nbf)
      work3 = 0
      call dgemm('N','N', nbf, nbf, nbf, &
                 1.0_dp, work2, nbf, &
                         mo, nbf, &
                 0.0_dp, work3,  nbf)
      x2mat_b = x2mat_b + work3(nocc_b+1:,1:nocc_b)

      k = 0
      do i = 1, nvir_b
        do a = 1, nocc_b
          k= k+1
          x2(k) = 2*x2mat(i,a)
        end do
      end do
    end select

    deallocate(foo,fvv,work1,work2,work3,dm,v,xmat,x2mat)

  end subroutine calc_h_op

  subroutine rotate_orbs_trah(infos, step, nbf,nocc_a, mo)
    use types, only : information
    use oqp_tagarray_driver
    implicit none

    class(information), intent(inout), target :: infos
    real(kind=dp), intent(in)        :: step(:)
    integer, intent(in)              :: nbf,nocc_a
    real(kind=dp), intent(inout)     :: mo(:,:)
    real(kind=dp), allocatable       :: work_1(:,:), work_2(:,:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    integer            :: i, idx
    logical :: second_term

    allocate(work_1(nbf,nbf),work_2(nbf,nbf))
    work_1 = 0
    work_2 = 0
    idx = 0
    second_term = .true.
    if (infos%control%scftype == 3) then! ROHF
      second_term = .false.
    end if
    call exp_scaling(work_1, step, idx, nocc_a, nbf, second_term)

    call orthonormalize(work_1, nbf)
    call dgemm('N','N', nbf, nbf, nbf, 1.0_dp, mo, nbf, work_1, nbf, 0.0_dp, work_2, nbf)
    mo = work_2
 !   mo_a = mo


  contains

    subroutine exp_scaling(G, step, idx, nocc_a, nbf, second_term)
      use matrix_expm_mod, only: expm
      real(kind=dp), intent(out)     :: G(:,:)
      real(kind=dp), intent(in)      :: step(:)
      integer, intent(inout)         :: idx
      integer, intent(in)            :: nocc_a, nbf

      real(kind=dp), allocatable     :: K(:,:), K2(:,:)
      integer                        :: occ, virt, istart, i
      logical, intent(inout)         :: second_term
      integer :: info

      allocate(K(nbf,nbf), source=0.0_dp)
      allocate(K2(nbf,nbf), source=0.0_dp)

      istart = nocc_a +1 !merge(nocc_b+1, nocc_a+1, occ <= nocc_b)
      do virt = istart, nbf
        do occ = 1, nocc_a
          idx = idx + 1
          K(virt, occ) =  step(idx)
          K(occ, virt) = -step(idx)
        end do
      end do

      G = 0.0_dp
      do i = 1, nbf
        G(i,i) = 1.0_dp
      end do
      G = G + K
      if (second_term) then
        call dgemm('N','T', nbf, nbf, nbf, 1.0_dp, K, nbf, K, nbf, 0.0_dp, K2, nbf)
        G = G + 0.5_dp * K2
      end if
!      call expm(nbf, K, nbf, G, nbf, info)

      deallocate(K, K2)
    end subroutine exp_scaling

    subroutine orthonormalize(G, nbf)
      real(kind=dp), intent(inout) :: G(:,:)
      integer, intent(in)          :: nbf
      integer                      :: i, j
      real(kind=dp)                :: norm, dot

      do i = 1, nbf
        norm = sqrt(dot_product(G(:,i), G(:,i)))
        call dscal(nbf, 1.0_dp / norm, G(1,i), 1)
        if (i == nbf) cycle
        do j = i+1, nbf
          dot = dot_product(G(:,i), G(:,j))
          call daxpy(nbf, -dot, G(1,i), 1, G(1,j), 1)
        end do
      end do
    end subroutine orthonormalize
  end subroutine rotate_orbs_trah

  function compute_energy(infos) result(etot)
    use types, only: information
    implicit none
    type(information), intent(in):: infos
    real(dp)  :: etot
    etot = infos%mol_energy%energy
  end function

  subroutine get_fock(basis, infos, molgrid, fock_ao, mo_a_in, mo_b_in)
    use precision, only: dp
    use oqp_tagarray_driver
    use types, only: information
    use int2_compute, only: int2_compute_t, int2_fock_data_t, int2_rhf_data_t, int2_urohf_data_t
    use dft, only: dftexcor
    use mod_dft_molgrid, only: dft_grid_t
    use basis_tools, only: basis_set
    use util, only: e_charge_repulsion
    use messages, only: show_message, WITH_ABORT
    use mathlib, only: traceprod_sym_packed, matrix_invsqrt, pack_matrix, unpack_matrix
    use guess, only: get_ab_initio_density
    implicit none

    type(basis_set), intent(in) :: basis
    type(information), target, intent(inout) :: infos
    type(dft_grid_t), intent(in) :: molgrid
    real(dp), pointer, intent(inout)        :: fock_ao(:,:)
    real(dp), intent(inout), optional :: mo_a_in(:,:)
    real(dp), intent(inout), optional :: mo_b_in(:,:)

    integer :: nbf, nbf_tri, nfocks, scf_type, nelec, nelec_a, nelec_b
    integer :: i, ii, ok
    real(kind=dp) :: vshift, scalefactor
    logical :: is_dft
    type(int2_compute_t) :: int2_driver
    class(int2_fock_data_t), allocatable :: int2_data

    !==============================================================================
    ! Matrices and Vectors for SCF Calculation
    !==============================================================================
    real(kind=dp), allocatable, target :: smat_full(:,:)  ! Full overlap matrix
    real(kind=dp), allocatable, target :: pdmat(:,:)  ! Density matrices in triangular format
    real(kind=dp), allocatable, target :: pfock(:,:)  ! Fock matrices in triangular format
    real(kind=dp), allocatable, target :: rohf_bak(:,:)  ! Backup for ROHF Fock
    real(kind=dp), allocatable, target :: dold(:,:)  ! Old density for incremental builds
    real(kind=dp), allocatable, target :: fold(:,:)  ! Old Fock for incremental builds
    real(kind=dp), allocatable :: pfxc(:,:)  ! DFT exchange-correlation matrix
    real(kind=dp), allocatable :: qmat(:,:)  ! Orthogonalization matrix
    real(kind=dp), allocatable, target :: work1(:,:)  ! Work matrix 1
    real(kind=dp), allocatable, target :: work2(:,:)  ! Work matrix 2

    !==============================================================================
    ! Energy Components
    !==============================================================================
    real(kind=dp) :: ehf      ! Electronic energy (HF part)
    real(kind=dp) :: ehf1     ! One-electron energy
    real(kind=dp) :: nenergy  ! Nuclear repulsion energy
    real(kind=dp) :: etot     ! Total SCF energy
    real(kind=dp) :: e_old    ! Energy from previous iteration
    real(kind=dp) :: psinrm   ! Wavefunction normalization
    real(kind=dp) :: vne      ! Nucleus-electron potential energy
    real(kind=dp) :: vnn      ! Nucleus-nucleus potential energy
    real(kind=dp) :: vee      ! Electron-electron potential energy
    real(kind=dp) :: vtot     ! Total potential energy
    real(kind=dp) :: virial   ! Virial ratio (V/T)
    real(kind=dp) :: tkin     ! Kinetic energy
    real(kind=dp) :: eexc     ! Exchange-correlation energy for DFT
    real(kind=dp) :: totele   ! Total electron density for DFT
    real(kind=dp) :: totkin   ! Total kinetic energy for DFT

    !==============================================================================
    ! Tag Arrays for Accessing Data
    !==============================================================================
    real(kind=dp), contiguous, pointer :: smat(:), hcore(:), tmat(:), &
                                          fock_a(:), fock_b(:), &
                                          dmat_a(:), dmat_b(:), &
                                          mo_energy_b(:), mo_energy_a(:), &
                                          mo_a(:,:), mo_b(:,:)
    character(len=*), parameter :: tags_general(3) = &
      (/ character(len=80) :: OQP_SM, OQP_TM, OQP_Hcore /)
    character(len=*), parameter :: tags_alpha(4) = &
      (/ character(len=80) :: OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(4) = &
      (/ character(len=80) :: OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    ! SCF type codes
    integer, parameter :: scf_rhf = 1, scf_uhf = 2, scf_rohf = 3

    ! 1. Get SCF type and electron counts
    select case (infos%control%scftype)
    case (1)
      scf_type = scf_rhf
      nfocks = 1
    case (2)
      scf_type = scf_uhf
      nfocks = 2
    case (3)
      scf_type = scf_rohf
      nfocks = 2
    end select
    nelec   = infos%mol_prop%nelec
    nelec_a = infos%mol_prop%nelec_a
    nelec_b = infos%mol_prop%nelec_b
    nbf     = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    ! 2. DFT and scale
    is_dft = infos%control%hamilton >= 20
    if (is_dft) then
      scalefactor = infos%dft%HFscale
    else
      scalefactor = 1.0_dp
    end if

    !==============================================================================
    ! Retrieve Tag Arrays and Allocate Memory
    !==============================================================================
    ! Get general tag arrays
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_TM, tmat)

    ! Get alpha-spin tag arrays
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    ! Get beta-spin tag arrays if needed
    if (nfocks > 1) then
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    end if

    ! Unpack overlap matrix to full for ROHF/level shift
    allocate(smat_full(nbf, nbf))
    call unpack_matrix(smat, smat_full, nbf, 'U')
    if (present(mo_a_in)) then
       mo_a = mo_a_in
    end if
    if (present(mo_b_in)) then
       mo_b = mo_b_in
    end if
    ! 4. Compute Nuclear-Nuclear Repulsion Energy
    nenergy = e_charge_repulsion(infos%atoms%xyz, infos%atoms%zn - infos%basis%ecp_zn_num)

    ! 5. Allocate Fock and density arrays
    allocate(pfock(nbf_tri, nfocks), pdmat(nbf_tri, nfocks),&
            rohf_bak(nbf_tri, nfocks), work1(nbf, nbf), work2(nbf, nbf))
    allocate(pfxc(nbf_tri, nfocks), &
         stat=ok, &
         source=0.0_dp)
    pfock = 0.0_dp
    pdmat = 0.0_dp

    ! 6. Set density matrices
    pdmat(:,1) = dmat_a
    if (nfocks > 1) pdmat(:,2) = dmat_b

    ! 8. --- FOCK BUILD ---
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    select case (scf_type)
    case (scf_rhf)
      allocate(int2_rhf_data_t :: int2_data)
      int2_data = int2_rhf_data_t(nfocks=1, d=pdmat, scale_exchange=scalefactor)
    case (scf_uhf, scf_rohf)
      allocate(int2_urohf_data_t :: int2_data)
      int2_data = int2_urohf_data_t(nfocks=2, d=pdmat, scale_exchange=scalefactor)
    end select
    call int2_driver%run(int2_data)

    ! 9. --- SCALE TWO-ELECTRON TERMS ---
    pfock(:,:) = 0.5_dp * int2_data%f(:,:,1)
    ii = 0
    do i = 1, nbf
      ii = ii + i
      pfock(ii,1:nfocks) = 2.0_dp * pfock(ii,1:nfocks)
    end do

    ehf = 0.0_dp
    ehf1 = 0.0_dp

    ! 10. --- ADD CORE HAMILTONIAN ---
    do i = 1, nfocks
      pfock(:,i) = pfock(:,i) + hcore
    end do

    ! Compute one and two-electron energies
    do i = 1, nfocks
      ehf1 = ehf1 + traceprod_sym_packed(pdmat(:,i), hcore, nbf)
      ehf = ehf + traceprod_sym_packed(pdmat(:,i), pfock(:,i), nbf)
    end do

    ! Total HF energy = 0.5*(E1 + E2) (to avoid double-counting)
    ehf = 0.5_dp * (ehf + ehf1)
    etot = ehf + nenergy


    !----------------------------------------------------------------------------
    ! Compute DFT Exchange-Correlation Contribution (if DFT)
    !----------------------------------------------------------------------------
    if (is_dft) then
      if (scf_type == scf_rhf) then
        call dftexcor(basis, molgrid, 1, pfxc, pfxc, mo_a, mo_a, &
                      nbf, nbf_tri, eexc, totele, totkin, infos)
      else if (scf_type == scf_uhf) then
        call dftexcor(basis, molgrid, 2, pfxc(:,1), pfxc(:,2), mo_a, mo_b, &
                      nbf, nbf_tri, eexc, totele, totkin, infos)
      else if (scf_type == scf_rohf) then
        ! ROHF does not have MO_B. So we copy MO_A to MO_B.
        mo_b = mo_a
        call dftexcor(basis, molgrid, 2, pfxc(:,1), pfxc(:,2), mo_a, mo_b, &
                      nbf, nbf_tri, eexc, totele, totkin, infos)
      end if

      ! Add XC contribution to Fock and total energy
      pfock = pfock + pfxc
      etot = etot + eexc
    end if


    if (scf_type == scf_rohf) then
      call form_rohf_fock_b(pfock(:,1), pfock(:,2), mo_a, smat_full, nelec_a, nelec_b, nbf, vshift, work1, work2)
      pdmat(:,1) = pdmat(:,1) + pdmat(:,2)
    end if
    call get_ab_initio_density(pdmat(:,1),mo_a,pdmat(:,2),mo_b,infos,basis)

    select case (scf_type)
    case (scf_rhf)
      fock_a = pfock(:,1)
      dmat_a = pdmat(:,1)
      fock_ao(:,1) = pfock(:,1)
    case (scf_uhf)
      fock_a = pfock(:,1)
      fock_b = pfock(:,2)
      dmat_a = pdmat(:,1)
      dmat_b = pdmat(:,2)
      fock_ao(:,1) = pfock(:,1)
      fock_ao(:,2) = pfock(:,2)
    case (scf_rohf)
      fock_a = rohf_bak(:,1)
      dmat_a = pdmat(:,1) - pdmat(:,2)
      dmat_b = pdmat(:,2)
      mo_b = mo_a
      mo_energy_b = mo_energy_a
      fock_ao(:,1) = rohf_bak(:,1)
      fock_ao(:,2) = rohf_bak(:,2)
    end select

    fock_ao(:,1) = pfock(:,1)
    infos%mol_energy%energy = etot

  end subroutine get_fock
  !> @brief Forms the ROHF Fock matrix in the MO basis using the Guest-Saunders method.
  !> @detail Transforms alpha and beta Fock matrices from the AO basis to the MO basis,
  !>         constructs the ROHF Fock matrix following the Guest-Saunders approach,
  !>         and optionally applies a level shift to virtual orbitals.
  !>         Reference: M. F. Guest, V. Saunders. Mol. Phys. 28, 819 (1974).
  !> @author Konstantin Komarov, 2023
  !> @param[inout] fock_a_ao Alpha Fock matrix in AO basis (triangular format).
  !> @param[inout] fock_b_ao Beta Fock matrix in AO basis (triangular format).
  !> @param[in] mo_a Alpha MO coefficients.
  !> @param[in] smat_full Full overlap matrix.
  !> @param[in] nocca Number of occupied alpha orbitals.
  !> @param[in] noccb Number of occupied beta orbitals.
  !> @param[in] nbf Number of basis functions.
  !> @param[in] vshift Level shift parameter for virtual orbitals.
  !> @param[inout] work1 Work array 1 (nbf x nbf).
  !> @param[inout] work2 Work array 2 (nbf x nbf).
  subroutine form_rohf_fock_b(fock_a_ao, fock_b_ao, &
                            mo_a, smat_full, &
                            nocca, noccb, nbf, vshift, &
                            work1, work2)
    use precision, only: dp
    use mathlib, only: orthogonal_transform_sym, &
                       orthogonal_transform2, &
                       unpack_matrix, &
                       pack_matrix

    implicit none

    real(kind=dp), intent(inout), dimension(:) :: fock_a_ao
    real(kind=dp), intent(inout), dimension(:) :: fock_b_ao
    real(kind=dp), intent(in), dimension(:,:) :: mo_a
    real(kind=dp), intent(in), dimension(:,:) :: smat_full
    real(kind=dp), intent(inout), dimension(:,:) :: work1
    real(kind=dp), intent(inout), dimension(:,:) :: work2
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: vshift

    real(kind=dp), allocatable, dimension(:) :: fock_mo
    real(kind=dp), allocatable, dimension(:,:) :: &
          work_matrix, fock, fock_a, fock_b
    real(kind=dp) :: acc, aoo, avv, bcc, boo, bvv
    integer :: i, nbf_tri

    acc = 0.5_dp; aoo = 0.5_dp; avv = 0.5_dp
    bcc = 0.5_dp; boo = 0.5_dp; bvv = 0.5_dp
    nbf_tri = nbf * (nbf + 1) / 2

    ! Allocate full matrices
    allocate(work_matrix(nbf, nbf), &
             fock(nbf, nbf), &
             fock_mo(nbf_tri), &
             fock_a(nbf, nbf), &
             fock_b(nbf, nbf), &
             source=0.0_dp)

    ! Transform alpha and beta Fock matrices to MO basis
    call orthogonal_transform_sym(nbf, nbf, fock_a_ao, mo_a, nbf, fock_mo)
    fock_a_ao(:nbf_tri) = fock_mo(:nbf_tri)

    call orthogonal_transform_sym(nbf, nbf, fock_b_ao, mo_a, nbf, fock_mo)
    fock_b_ao(:nbf_tri) = fock_mo(:nbf_tri)

    ! Unpack triangular matrices to full matrices
    call unpack_matrix(fock_a_ao, fock_a)
    call unpack_matrix(fock_b_ao, fock_b)

    ! Construct ROHF Fock matrix in MO basis using Guest-Saunders method
    associate ( na => nocca &
              , nb => noccb &
      )
      fock(1:nb, 1:nb) = acc * fock_a(1:nb, 1:nb) &
                       + bcc * fock_b(1:nb, 1:nb)
      fock(nb+1:na, nb+1:na) = aoo * fock_a(nb+1:na, nb+1:na) &
                             + boo * fock_b(nb+1:na, nb+1:na)
      fock(na+1:nbf, na+1:nbf) = avv * fock_a(na+1:nbf, na+1:nbf) &
                               + bvv * fock_b(na+1:nbf, na+1:nbf)
      fock(1:nb, nb+1:na) = fock_b(1:nb, nb+1:na)
      fock(nb+1:na, 1:nb) = fock_b(nb+1:na, 1:nb)
      fock(1:nb, na+1:nbf) = 0.5_dp * (fock_a(1:nb, na+1:nbf) &
                                     + fock_b(1:nb, na+1:nbf))
      fock(na+1:nbf, 1:nb) = 0.5_dp * (fock_a(na+1:nbf, 1:nb) &
                                     + fock_b(na+1:nbf, 1:nb))
      fock(nb+1:na, na+1:nbf) = fock_a(nb+1:na, na+1:nbf)
      fock(na+1:nbf, nb+1:na) = fock_a(na+1:nbf, nb+1:na)

      ! Apply Vshift to the diagonal
      do i = nb+1, na
        fock(i,i) = fock(i,i) + vshift * 0.5_dp
      end do
      do i = na+1, nbf
        fock(i,i) = fock(i,i) + vshift
      end do
    end associate

    ! Back-transform ROHF Fock matrix to AO basis
    call dsymm('l', 'u', nbf, nbf, &
               1.0_dp, smat_full, nbf, &
                       mo_a, nbf, &
               0.0_dp, work1, nbf)
    call orthogonal_transform2('t', nbf, nbf, work1, nbf, fock, nbf, &
                               work_matrix, nbf, work2)
    call pack_matrix(work_matrix, fock_a_ao)

    deallocate(work_matrix, fock, fock_mo, fock_a, fock_b)

  end subroutine form_rohf_fock_b


end module trah

!==============================================================================
! Module: trah_opentrustregion_interface
! A simple interface to call OpenTrustRegion solver from OpenQP without needing
! opentrustregion_solver_class
!==============================================================================
module otr_interface

  use, intrinsic :: iso_fortran_env, only: int32
  use opentrustregion, only: solver, update_orbs_type, obj_func_type, hess_x_type, logger_type
  use mathlib,           only: unpack_matrix
  use trah,        only: calc_g_h, calc_h_op, compute_energy,  rotate_orbs_trah, get_fock
  use precision, only: dp
  use types, only:information
  use mod_dft_molgrid, only: dft_grid_t
  use basis_tools,      only: basis_set
  use guess, only:get_ab_initio_density
  implicit none

  ! Module-level state for callbacks
  class(information), pointer :: infos       ! OpenQP information object
  type(dft_grid_t), pointer :: molgrid
  real(dp),         pointer :: fock_ao(:,:)  ! Packaged AO Fock matrix
  real(dp),         pointer :: mo_a(:,:)     ! MO alpha coefficient matrix
  real(dp),         pointer :: mo_b(:,:)     ! MO beta coefficient matrix
  integer                   :: nbf, nocc_a, nocc_b, nvir_a, nvir_b, nfock
  integer(int32)            :: n_param      ! number of parameters = nocc*nvir
  integer(int32)            :: max_iter     ! macro iteration limit
  real(dp)                  :: conv_tol     ! convergence tolerance
  integer(int32)            :: verbose      ! verbosity level

contains

  subroutine init_trah_solver(infos_in, molgrid_in, fock_in, mo_a_in,&
                  mo_b_in)
    class(information), intent(inout), target :: infos_in
    type(dft_grid_t), intent(in), target :: molgrid_in
    real(dp), intent(inout),  target :: fock_in(:,:)
    real(dp), intent(inout),  target :: mo_a_in(:,:)
    real(dp), intent(inout),  target :: mo_b_in(:,:)
    type(basis_set), pointer :: basis
    integer :: nvir
!    character(len=*), intent(in)          :: print_level
    ! Initialize module state
    infos   => infos_in
    fock_ao => fock_in
    molgrid => molgrid_in
    mo_a    => mo_a_in
    mo_b    => mo_b_in

    basis => infos_in%basis
    nbf = basis%nbf
    nocc_a  = infos%mol_prop%nelec_a
    nocc_b  = infos%mol_prop%nelec_b
    nvir_a    = nbf - nocc_a
    nvir_b    = nbf - nocc_b

    select case (infos%control%scftype)
    case (1) !RHF
      n_param = int(nocc_a * nvir_a, kind=int32)
      nfock = 1
    case (2) !UHF
      n_param = int((nocc_a * nvir_a) +(nocc_b * nvir_b), kind=int32)
      nfock = 2
    case (3) !ROHF
      nvir = min(nvir_a, nvir_b)
      n_param = int((nocc_a * nvir_a) + (nocc_b * nvir_a), kind=int32)
      nfock = 2
    end select

    max_iter = int(infos%control%maxit, kind=int32)
    conv_tol = infos%control%conv
    verbose  = int(3, kind=int32)
    select case (infos%control%scftype)
    case (1)
      call get_fock(basis, infos_in, molgrid, fock_ao,mo_a)
    case (2)
      call get_fock(basis, infos_in, molgrid, fock_ao,mo_a,mo_b)
    end select 

  end subroutine init_trah_solver

  subroutine run_trah_solver()
    logical(kind=4) :: error
    procedure(update_orbs_type), pointer :: p_update
    procedure(obj_func_type),   pointer :: p_obj
    procedure(logger_type),     pointer :: p_log

    ! Bind callbacks
    p_update => update_orbs
    p_obj    => obj_func
    p_log    => logger

    ! Call OpenTrustRegion solver
    call solver(p_update, p_obj, n_param, error, &
!                conv_tol=conv_tol, n_macro=max_iter, &
                verbose=verbose)! , logger=p_log)
    if (error) then
      write(*,*) 'OpenTrustRegion solver failed.'
    end if
  end subroutine run_trah_solver

  subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr)
    real(dp), intent(in)                     :: kappa(:)
    real(dp), intent(out)                    :: func
    real(dp), intent(out)                    :: grad(:), h_diag(:)
    procedure(hess_x_type), pointer, intent(out) :: hess_x_funptr
    type(basis_set), pointer :: basis
    basis => infos%basis
    ! Rotate orbitals
    select case (infos%control%scftype)
    case (1)
      call rotate_orbs_trah(infos, kappa, nbf, nocc_a, mo_a)
      call get_fock(basis, infos, molgrid, fock_ao, mo_a)
   ! Gradient & Hessian diag
      call calc_g_h(grad, h_diag, fock_ao, mo_a, nbf, nocc_a, infos%control%scftype)
    case (2)
      call rotate_orbs_trah(infos, kappa(1:nocc_a*nvir_a), nbf, nocc_a, mo_a)
      call rotate_orbs_trah(infos, kappa(nocc_a*nvir_a+1:), nbf, nocc_b, mo_b)
      call get_fock(basis, infos, molgrid, fock_ao, mo_a, mo_b)
      call calc_g_h(grad, h_diag, fock_ao, mo_a, nbf, nocc_a, infos%control%scftype, mo_b, nocc_b)
    end select
    func = compute_energy(infos)
    hess_x_funptr => hess_x_cb
    h_diag = 2.0_dp * h_diag
    grad = 2.0_dp * grad
  end subroutine update_orbs


  function hess_x_cb(x) result(hx)
    real(dp), intent(in) :: x(:)
    real(dp)             :: hx(size(x))
    select case (infos%control%scftype)
    case (1)
      call calc_h_op(infos, fock_ao, x, hx, mo_a)
    case (2)
      call calc_h_op(infos, fock_ao, x, hx, mo_a, mo_b)
    end select
    hx = 2.0_dp * hx
  end function hess_x_cb


  function obj_func(kappa) result(val)
    real(dp), intent(in) :: kappa(:)
    real(dp)             :: val
    real(dp), allocatable :: mo_tmp_a(:,:), mo_tmp_b(:,:),dens(:,:)
    type(basis_set), pointer :: basis
    basis => infos%basis
    allocate(mo_tmp_a(nbf,nbf))
    allocate(dens(nbf*(nbf+1)/2,nfock))
    select case(infos%control%scftype)
    case (1)
      mo_tmp_a = mo_a
      call rotate_orbs_trah(infos, kappa, nbf, nocc_a, mo_a)
      call get_fock(basis, infos, molgrid, fock_ao, mo_a)
      val = compute_energy(infos)
      mo_a = mo_tmp_a
    case (2)
      allocate(mo_tmp_b(nbf,nbf))
      mo_tmp_a = mo_a
      mo_tmp_b = mo_b
      call rotate_orbs_trah(infos, kappa(1:nvir_a*nocc_a), nbf, nocc_a, mo_a)
      call rotate_orbs_trah(infos, kappa(nvir_a*nocc_a+1:), nbf, nocc_b, mo_b)
      call get_ab_initio_density(dens(:,1),mo_a,dens(:,2),mo_b,infos,basis)
      call get_fock(basis, infos, molgrid, fock_ao, mo_a, mo_b,dens)
      val = compute_energy(infos)
      mo_a = mo_tmp_a
      mo_b = mo_tmp_b
      deallocate(mo_tmp_b)
    end select
    deallocate(mo_tmp_a,dens)
  end function obj_func


  subroutine logger(message)
    use io_constants, only: IW
    implicit none
    character(*), intent(in) :: message
    write(IW, "(A)") trim(message)
  end subroutine 

end module otr_interface

module trah
  use precision, only: dp
  character(len=*), parameter :: module_name = "trah"
  public :: calc_g_k
  public :: calc_h_op
  public :: compute_energy
  public :: rotate_orbs_trah

contains

  subroutine calc_g_h(grad, h_diag, fock_ao, mo, nbf, nocc)
    use mathlib, only: pack_matrix, unpack_matrix
    implicit none
    integer, intent(in) :: nbf, nocc
    real(dp), intent(out) :: grad(:)
    real(dp), intent(out) :: h_diag(:)
    real(dp), intent(in)  :: fock_ao(:)      ! packed AO fock
    real(dp), intent(inout)  :: mo(nbf,nbf)     ! MO coefficient matrix
!    real(dp), allocatable :: grad(:)
!    real(dp), allocatable :: h_diag(:)

    integer :: nvir, i, a, k
    real(dp), allocatable :: foo(:,:), fvv(:,:)
    real(dp), allocatable :: work1(:,:), work2(:,:), work3(:,:)

    ! dimensions
    nvir = nbf - nocc

    allocate(foo(nocc,nocc))
    allocate(fvv(nvir,nvir))
    allocate(work1(nbf,nbf))
    allocate(work2(nbf,nbf))
    allocate(work3(nbf,nbf))
!    allocate(grad(nocc*nvir))
!    allocate(h_diag(nocc*nvir))
    if (size(grad)   /= nocc*nvir) stop 'grad wrong size'
    if (size(h_diag)/= nocc*nvir) stop 'h_diag wrong size'

    ! unpack the AO fock into a full (nbf×nbf) matrix
    call unpack_matrix(fock_ao, work1)
    !----------------------------------------------------------------
    ! Initialize fock_ao (13×13)
!    !----------------------------------------------------------------
!    work1 = reshape((/ &
!      -2.05297329d+01, -5.12554089d+00, -3.63750957d+00, -7.69087435d-18,  2.08166817d-17, -3.10181652d-02, -1.01872724d-18,  2.33103467d-18, -5.08017286d-03, -7.76993049d-01, -1.51167851d+00, -7.76993049d-01, -1.51167851d+00, &
!      -5.12554089d+00, -2.12539314d+00, -2.11565100d+00,  9.89881715d-18, -1.73472348d-17, -9.28013564d-02,  1.07187122d-17, -7.80625564d-18, -7.99286942d-02, -8.09999644d-01, -1.15639865d+00, -8.09999644d-01, -1.15639865d+00, &
!      -3.63750957d+00, -2.11565100d+00, -1.67093694d+00, -5.46704089d-18, -2.94902991d-17, -1.23552132d-01, -6.01333577d-18, -1.73472348d-17, -2.04463839d-01, -7.17852599d-01, -9.56671584d-01, -7.17852599d-01, -9.56671584d-01, &
!      -7.69087435d-18,  9.89881715d-18, -5.46704089d-18, -1.76471020d-01,  6.35997945d-17, -7.09239599d-19, -6.57913299d-01,  7.21131869d-17, -2.77852637d-17, -2.58558417d-17, -3.96161958d-17, -1.58725169d-17, -3.76764130d-18, &
!       2.08166817d-17, -1.73472348d-17, -2.94902991d-17,  6.35997945d-17, -7.67989714d-02, -3.77953317d-17,  9.08747298d-17, -6.88096400d-01, -1.21430643d-17,  4.77757288d-01,  2.49871627d-01, -4.77757288d-01, -2.49871627d-01, &
!      -3.10181652d-02, -9.28013564d-02, -1.23552132d-01, -7.09239599d-19, -3.77953317d-17, -1.58833933d-01, -3.53491507d-18,  2.60208521d-17, -6.94009791d-01, -3.90159710d-01, -2.61815588d-01, -3.90159710d-01, -2.61815588d-01, &
!      -1.01872724d-18,  1.07187122d-17, -6.01333577d-18, -6.57913299d-01,  9.08747298d-17, -3.53491507d-18,  1.62158567d-02,  8.74313656d-17, -2.22193561d-17, -4.39479830d-17, -4.48946333d-17,  9.19224786d-18,  5.76911160d-18, &
!       2.33103467d-18, -7.80625564d-18, -1.73472348d-17,  7.21131869d-17, -6.88096400d-01,  2.60208521d-17,  8.74313656d-17, -2.08737204d-01,  1.56125113d-17,  3.31259385d-01,  1.29688274d-01, -3.31259385d-01, -1.29688274d-01, &
!      -5.08017286d-03, -7.99286942d-02, -2.04463839d-01, -2.77852637d-17, -1.21430643d-17, -6.94009791d-01, -2.22193561d-17,  1.56125113d-17, -1.17854236d-01, -2.67761024d-01, -2.12305595d-01, -2.67761024d-01, -2.12305595d-01, &
!      -7.76993049d-01, -8.09999644d-01, -7.17852599d-01, -2.58558417d-17,  4.77757288d-01, -3.90159710d-01, -4.39479830d-17,  3.31259385d-01, -2.67761024d-01, -2.52870265d-01, -6.51857992d-01, -1.77347607d-01, -3.74762076d-01, &
!      -1.51167851d+00, -1.15639865d+00, -9.56671584d-01, -3.96161958d-17,  2.49871627d-01, -2.61815588d-01, -4.48946333d-17,  1.29688274d-01, -2.12305595d-01, -6.51857992d-01, -5.91224421d-01, -3.74762076d-01, -5.30630294d-01, &
!      -7.76993049d-01, -8.09999644d-01, -7.17852599d-01, -1.58725169d-17, -4.77757288d-01, -3.90159710d-01,  9.19224786d-18, -3.31259385d-01, -2.67761024d-01, -1.77347607d-01, -3.74762076d-01, -2.52870265d-01, -6.51857992d-01, &
!      -1.51167851d+00, -1.15639865d+00, -9.56671584d-01, -3.76764130d-18, -2.49871627d-01, -2.61815588d-01,  5.76911160d-18, -1.29688274d-01, -2.12305595d-01, -3.74762076d-01, -5.30630294d-01, -6.51857992d-01, -5.91224421d-01  /), (/13,13/))
!
!    !----------------------------------------------------------------
!    ! Initialize mo_coeff0 (13×13)
!    !----------------------------------------------------------------
!    mo = reshape((/ &
!       9.95706989d-01, -2.11589279d-01, -1.66948851d-16,  7.63844969d-02,  1.40972527d-16, -9.04076480d-02,  3.10017110d-15,  1.23048202d-15,  5.32034281d-02,  2.44660535d-16, -6.38119630d-02, -5.91239204d-16,  5.50652084d-02, &
!       2.23514367d-02,  4.67713904d-01,  8.65748312d-16, -1.90973008d-01, -3.62977450d-16,  1.22323369d-01, -4.46385980d-15, -1.04813563d-14, -2.83020887d-01, -1.33851438d-15,  3.88708047d-01,  7.64037207d-15, -1.68077550d+00, &
!      -8.07967478d-03,  4.68129768d-01, -8.75448116d-16, -2.70477112d-01, -5.35585665d-16,  1.18147474d+00, -3.69843335d-14,  1.25529073d-14,  2.06937018d-01,  8.94582031d-16,  1.24635420d-02, -1.06053210d-14,  2.76358718d+00, &
!       6.21172186d-19,  8.33058317d-18, -3.91583901d-16, -1.09237882d-15,  6.53918573d-01,  1.18250517d-16,  2.62978170d-16, -2.06012770d-16,  3.77260075d-15, -9.53119518d-01, -5.37979039d-16,  9.11916972d-16, -7.94957768d-17, &
!      -1.80785178d-18, -1.31187374d-16,  4.81056586d-01,  1.07160450d-15,  3.04172129d-16,  9.14032059d-15,  3.59097413d-01,  1.21429542d-01, -1.03793883d-14, -1.21748643d-15,  5.52565695d-15, -1.03721085d+00, -2.64241966d-15, &
!       2.59421528d-03,  1.23913580d-01, -9.62917673d-16,  5.76343394d-01,  9.34199707d-16,  2.36402846d-01, -7.81472370d-15, -1.22077155d-14, -6.06424994d-01, -2.24729236d-15, -7.54379725d-01,  4.88490946d-16, -1.86054311d-01, &
!      -2.49008082d-19,  1.96918010d-18, -5.20304807d-16, -8.22148066d-16,  4.96633571d-01,  1.51215398d-17, -3.78705978d-17,  3.22623595d-16, -4.21993025d-15,  1.04374394d+00,  4.86023553d-16, -1.32110152d-15, -9.12704677d-18, &
!       2.26274282d-18,  7.27321017d-16,  2.72638707d-01, -5.10168713d-17, -9.53202400d-17,  2.79515896d-14,  8.37849052d-01,  6.97385450d-01, -6.16145308d-15,  1.73856622d-15, -6.14605649d-15,  1.55732280d+00,  2.26726584d-15, &
!      -2.00766990d-03,  7.06586478d-02, -1.28327870d-15,  4.06296573d-01,  7.00128288d-16,  4.52194838d-01, -1.49161098d-14,  5.00478694d-15,  1.76192831d-01,  2.39295927d-16,  1.21742351d+00,  1.91288732d-15,  8.40301895d-01, &
!      -7.99155806d-05,  1.27275153d-01, -2.52567084d-01,  1.17004427d-01, -8.13341731d-17, -8.15029772d-02,  4.64004123d-02,  9.80298930d-01,  8.76317193d-01,  3.49115755d-15, -4.91868874d-01, -6.47755149d-02, -4.60513476d-01, &
!       1.97486378d-03,  7.03797424d-03, -1.65134060d-01,  6.33046350d-02,  1.97076684d-17, -9.77741575d-01,  1.37096333d+00, -5.27952150d-01, -5.74377108d-01, -1.09679257d-15,  1.96434745d-02,  9.01823806d-01, -6.15951017d-01, &
!      -7.99155806d-05,  1.27275153d-01,  2.52567084d-01,  1.17004427d-01,  5.30782588d-16, -8.15029772d-02, -4.64004123d-02, -9.80298930d-01,  8.76317193d-01,  4.19721860d-15, -4.91868874d-01,  6.47755149d-02, -4.60513476d-01, &
!       1.97486378d-03,  7.03797424d-03,  1.65134060d-01,  6.33046350d-02,  3.26603876d-16, -9.77741575d-01, -1.37096333d+00,  5.27952150d-01, -5.74377108d-01, -3.73786870d-15,  1.96434745d-02, -9.01823806d-01, -6.15951017d-01  /), (/13,13/))
!
!    work1=transpose(work1)
!    mo=transpose(mo)
!    print *, "fock_ao", work1(:,1)
!    print *, "mo", mo(:,1)
    !—— compute gradient block: g_{i,a} = 2 * (C_i^T F C_a) ——
    ! work2(μ,a) = ∑_ν F_ao(μ,ν) · C(ν,nocc+a)
!    call dgemm('N','N', nbf, nvir, nbf, &
!               1.0_dp, work1, nbf,      &
!                       mo(:,nocc+1:nbf), nbf, &
!               0.0_dp, work2, nbf)
!
!    ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
!    call dgemm('T','N', nocc, nvir, nbf, &
!               1.0_dp, mo(:,1:nocc), nbf, &
!                       work2,         nbf, &
!               0.0_dp, work3(1:nocc, 1:nvir),nocc)
    print *, "fock_ao",work1
    print *, "in grad mo", mo
    work2 = 0.0_dp
    work3 = 0.0_dp
    call dgemm('N','N', nbf, nbf, nbf, &
               1.0_dp, work1, nbf,      &
                       mo, nbf, &
               0.0_dp, work2, nbf)

    ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
    call dgemm('T','N', nbf, nbf, nbf, &
               1.0_dp, mo, nbf, &
                       work2,         nbf, &
               0.0_dp, work3,nbf)
    print *, "work3",work3
!print *, "work3", work3
!    k = 0
!    do a = 1, nvir
!      do i = 1, nocc
!        k = k + 1
!        grad(k) = 2.0_dp * work3(i,a)
!      end do
!    end do

    k = 0
    do i = nocc +1, nbf
      do a = 1, nocc
        k= k+1
        grad(k) =  2 * work3(i,a)
      end do
    end do
!print *, "grad",grad
    !—— compute Hessian diagonal: h_{i,a;i,a} = 2*(F_{aa} - F_{ii}) ——
    ! build occupied block foo_occ and virtual block fvv_virt
    ! foo_occ(i,j) = ∑_μν C(μ,i)^T F_ao(μ,ν) C(ν,j)
!    call dgemm('N','N', nbf, nocc, nbf, &
!               1.0_dp, work1, nbf,      &
!                       mo(:,1:nocc), nbf, &
!               0.0_dp, work2,        nbf)
!    call dgemm('T','N', nocc, nocc, nbf, &
!               1.0_dp, mo(:,1:nocc), nbf, &
!                       work2,          nbf, &
!               0.0_dp, foo,          nocc)
!
!    ! fvv(a,b) = ∑_μν C(μ,nocc+a)^T F_ao(μ,ν) C(ν,nocc+b)
!    ! (re‑reuse work2)
!    call dgemm('N','N', nbf, nvir, nbf, &
!               1.0_dp, work1, nbf,      &
!                       mo(:,nocc+1:nbf), nbf, &
!               0.0_dp, work2,        nbf)
!    call dgemm('T','N', nvir, nvir, nbf, &
!               1.0_dp, mo(:,nocc+1:nbf), nbf, &
!                       work2,              nbf, &
!               0.0_dp, fvv,              nvir)

!    k = 0
!    do a = 1, nvir
!      do i = 1, nocc
!        k = k + 1
!        h_diag(k) = 2.0_dp*( fvv(a,a) - foo(i,i) )
!      end do
!    end do

    k = 0
    do i = nocc +1, nbf
      do a = 1, nocc
        k= k+1
        h_diag(k) = 2.0_dp*( work3(i,i) - work3(a,a) )
      end do
    end do
!    print *, "h_diag", h_diag
    ! clean up
    deallocate(foo, fvv, work1, work2, work3)

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

    call measure_time(print_total=1, log_unit=IW)

!    write(IW,"(/3x,'Form Two-Electron J and K Fock')")

    ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_data = int2_rhf_data_t(nfocks=1, d=d, scale_exchange=scalefactor)

!    call flush(IW)

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

  subroutine calc_h_op(infos,fock_ao, x, x2, mo, nbf ,nocc)
    use types, only: information
    use basis_tools,  only: basis_set
    use mathlib, only: pack_matrix, unpack_matrix
    implicit none

    class(information), intent(inout), target :: infos
    real(dp), intent(in)                 :: x(:)    ! length nocc*nvir
!    real(dp), allocatable                 :: x(:), x2(:)
    real(dp), intent(out)                :: x2(:)   ! same length
    real(dp), pointer, intent(in)        :: fock_ao(:)
    real(dp), intent(inout)                 :: mo(:, :)

    integer :: nbf, nocc, nvir, i, a, k
    real(dp), allocatable :: foo(:,:), fvv(:,:), xmat(:,:),&
            x2mat(:,:), v(:,:), work2(:,:), work1(:,:), &
            work3(:,:),temp(:,:)
    real(dp), allocatable :: pfock(:,:)
    real(dp), allocatable :: dm(:,:), dm_tri(:,:)
    real(kind=dp) :: scalefactor
    type(basis_set), pointer :: basis
    logical :: is_dft

    ! dims
    nvir = nbf - nocc

    is_dft = infos%control%hamilton >= 20
    if (is_dft) then
      scalefactor = infos%dft%HFscale
    else
      scalefactor = 1.0_dp
    end if
    basis => infos%basis

    allocate(foo(nocc,nocc), fvv(nvir,nvir))
    allocate(work3(nbf,nbf), work2(nbf,nbf), work1(nbf,nbf))
    allocate(xmat(nvir,nocc), x2mat(nvir,nocc), temp(nvir,nocc))
    allocate(dm(nbf,nbf), v(nbf,nbf))
!    allocate(x(nocc*nvir),x2(nocc*nvir))
    allocate(pfock(nbf*(nbf+1)/2,2),dm_tri(nbf*(nbf+1)/2,2))
    k = 0
    do i = 1, nvir
      do a = 1, nocc
        k= k+1
        xmat(i,a) = x(k)
      end do
    end do
!print *, "xmat",transpose(xmat)
!    xmat(1:nvir, 1:nocc) = reshape(x, [nvir, nocc])
    call unpack_matrix(fock_ao, work1)
!    work1 = reshape((/ &
!      -2.05297329d+01, -5.12554089d+00, -3.63750957d+00, -7.69087435d-18,  2.08166817d-17, -3.10181652d-02, -1.01872724d-18,  2.33103467d-18, -5.08017286d-03, -7.76993049d-01, -1.51167851d+00, -7.76993049d-01, -1.51167851d+00, &
!      -5.12554089d+00, -2.12539314d+00, -2.11565100d+00,  9.89881715d-18, -1.73472348d-17, -9.28013564d-02,  1.07187122d-17, -7.80625564d-18, -7.99286942d-02, -8.09999644d-01, -1.15639865d+00, -8.09999644d-01, -1.15639865d+00, &
!      -3.63750957d+00, -2.11565100d+00, -1.67093694d+00, -5.46704089d-18, -2.94902991d-17, -1.23552132d-01, -6.01333577d-18, -1.73472348d-17, -2.04463839d-01, -7.17852599d-01, -9.56671584d-01, -7.17852599d-01, -9.56671584d-01, &
!      -7.69087435d-18,  9.89881715d-18, -5.46704089d-18, -1.76471020d-01,  6.35997945d-17, -7.09239599d-19, -6.57913299d-01,  7.21131869d-17, -2.77852637d-17, -2.58558417d-17, -3.96161958d-17, -1.58725169d-17, -3.76764130d-18, &
!       2.08166817d-17, -1.73472348d-17, -2.94902991d-17,  6.35997945d-17, -7.67989714d-02, -3.77953317d-17,  9.08747298d-17, -6.88096400d-01, -1.21430643d-17,  4.77757288d-01,  2.49871627d-01, -4.77757288d-01, -2.49871627d-01, &
!      -3.10181652d-02, -9.28013564d-02, -1.23552132d-01, -7.09239599d-19, -3.77953317d-17, -1.58833933d-01, -3.53491507d-18,  2.60208521d-17, -6.94009791d-01, -3.90159710d-01, -2.61815588d-01, -3.90159710d-01, -2.61815588d-01, &
!      -1.01872724d-18,  1.07187122d-17, -6.01333577d-18, -6.57913299d-01,  9.08747298d-17, -3.53491507d-18,  1.62158567d-02,  8.74313656d-17, -2.22193561d-17, -4.39479830d-17, -4.48946333d-17,  9.19224786d-18,  5.76911160d-18, &
!       2.33103467d-18, -7.80625564d-18, -1.73472348d-17,  7.21131869d-17, -6.88096400d-01,  2.60208521d-17,  8.74313656d-17, -2.08737204d-01,  1.56125113d-17,  3.31259385d-01,  1.29688274d-01, -3.31259385d-01, -1.29688274d-01, &
!      -5.08017286d-03, -7.99286942d-02, -2.04463839d-01, -2.77852637d-17, -1.21430643d-17, -6.94009791d-01, -2.22193561d-17,  1.56125113d-17, -1.17854236d-01, -2.67761024d-01, -2.12305595d-01, -2.67761024d-01, -2.12305595d-01, &
!      -7.76993049d-01, -8.09999644d-01, -7.17852599d-01, -2.58558417d-17,  4.77757288d-01, -3.90159710d-01, -4.39479830d-17,  3.31259385d-01, -2.67761024d-01, -2.52870265d-01, -6.51857992d-01, -1.77347607d-01, -3.74762076d-01, &
!      -1.51167851d+00, -1.15639865d+00, -9.56671584d-01, -3.96161958d-17,  2.49871627d-01, -2.61815588d-01, -4.48946333d-17,  1.29688274d-01, -2.12305595d-01, -6.51857992d-01, -5.91224421d-01, -3.74762076d-01, -5.30630294d-01, &
!      -7.76993049d-01, -8.09999644d-01, -7.17852599d-01, -1.58725169d-17, -4.77757288d-01, -3.90159710d-01,  9.19224786d-18, -3.31259385d-01, -2.67761024d-01, -1.77347607d-01, -3.74762076d-01, -2.52870265d-01, -6.51857992d-01, &
!      -1.51167851d+00, -1.15639865d+00, -9.56671584d-01, -3.76764130d-18, -2.49871627d-01, -2.61815588d-01,  5.76911160d-18, -1.29688274d-01, -2.12305595d-01, -3.74762076d-01, -5.30630294d-01, -6.51857992d-01, -5.91224421d-01  /), (/13,13/))

    call dgemm('N','N', nbf, nbf, nbf, &
               1.0_dp, work1, nbf,      &
                       mo, nbf, &
               0.0_dp, work2, nbf)

    ! foo(i,a) = ∑_μ C(μ,i)^T · work2(μ,a)
    call dgemm('T','N', nbf, nbf, nbf, &
               1.0_dp, mo, nbf, &
                       work2,         nbf, &
               0.0_dp, work3,nbf)
     foo = work3(1:nocc,1:nocc)
     fvv = work3(nocc+1:nbf,nocc+1:nbf)

!     print *,"fvv", fvv
!    call dgemm('N','N', nbf, nocc, nbf, &
!               1.0_dp, work1, nbf, &
!                       mo(:,1:nocc),   nbf, &
!               0.0_dp, work2,        nbf)
!    call dgemm('T','N', nocc, nocc, nbf, &
!               1.0_dp, mo(:,1:nocc), nbf, &
!                       work2,          nbf, &
!               0.0_dp, foo,          nocc)
!    work2 = 0_dp
!    call dgemm('N','N', nbf, nvir, nbf, &
!               1.0_dp, work1, nbf, &
!                       mo(:,nocc+1:nbf),nbf, &
!               0.0_dp, work2,        nbf)
!    call dgemm('T','N', nvir, nvir, nbf, &
!               1.0_dp, mo(:,nocc+1:nbf), nbf, &
!                       work2,            nbf, &
!               0.0_dp, fvv,            nvir)

    call dgemm('N','N', nvir, nocc, nvir, &
               1.0_dp, fvv,  nvir, &
                       xmat, nvir, &
               0.0_dp, x2mat, nvir)
    call dgemm('N','N', nvir, nocc, nocc, &
              -1.0_dp, xmat, nvir, &
                        foo, nocc, &
               1.0_dp, x2mat, nvir) 
    ! Transpose temp(nocc,nvir) to x2mat(nvir,nocc), then do:
!    x2mat = x2mat + transpose(temp)

!    print *,"x2mat",transpose(x2mat)

!    call dgemm('N','N', nvir, nocc, nocc, &
!              -1.0_dp, xmat,  nvir, &
!                        foo,   nocc, &
!               1.0_dp, x2mat, nvir)

    dm = 0.0_dp
    work2 = 0
    call dgemm('N','N', nbf, nocc, nvir, &
               2.0_dp, mo(:, nocc+1:nbf), nbf, &
                       xmat,             nvir, &
               0.0_dp, work2,          nbf)
!    print *, "mo(:, nocc+1:nbf)"
    ! 2) work3 = work2 * orbo^T
    call dgemm('N','T', nbf, nbf, nocc, &
               1.0_dp, work2,          nbf, &
                       mo(:, 1:nocc),   nbf, &
               0.0_dp, work3,          nbf)
    ! 3) dm1 = work3 + work3^T
!    dm = work3
    do i = 1, nbf
      do a = 1, nbf
        dm(i,a) = work3(i,a) + work3(a,i)
      end do
    end do
!dm = 1.0_dp
print *, "dmat", dm
    call pack_matrix(dm,dm_tri(:,1))
    dm_tri(:,1) = dm_tri(:,1)
    call fock_jk(infos%basis, dm_tri, pfock, scalefactor, infos)
    call unpack_matrix(pfock(:,1), v)
!v = reshape((/ -0.29858673, -0.23395133,  0.0423584 , -0.16794323,&
! -0.23395133, -0.44040524, -0.16794323, -0.41066484,&
!  0.0423584 , -0.16794323, -0.29858673, -0.23395133,&
! -0.16794323, -0.41066484, -0.23395133, -0.44040524/), (/4,4/))
!print *, "x2 before 2e", x2mat
print *, "v", v   
    work2 = 0
print *, "orbv",mo(1, nocc+1:nbf)
!    call dgemm('T','N', nvir, nbf, nbf, &
!               1.0_dp, mo(:, nocc+1:nbf), nbf, &
!                       v ,             nbf, &
!               0.0_dp, work2,          nbf)
!    ! 2) work3 = work2 * orbo^T
!    work3 = 0
!    call dgemm('N','N', nvir, nocc, nbf, &
!               1.0_dp, work2(1:nvir,nbf), nbf, &
!                       mo(:, 1:nocc), nbf, &
!               0.0_dp, work3,  nbf)

    call dgemm('T','N', nbf, nbf, nbf, &
               1.0_dp, mo, nbf, &
                       v ,             nbf, &
               0.0_dp, work2,          nbf)
    ! 2) work3 = work2 * orbo^T
    work3 = 0
    call dgemm('N','N', nbf, nbf, nbf, &
               1.0_dp, work2, nbf, &
                       mo, nbf, &
               0.0_dp, work3,  nbf)
    x2mat = x2mat + work3(nocc+1:,1:nocc)
print *, "x2 afet 2e", work3(nocc+1:,1:nocc)
    
    k = 0
    do i = 1, nvir
      do a = 1, nocc
        k= k+1
        x2(k) = 2*x2mat(i,a)
      end do
    end do
print *, "x2k", x2
print *, "==========================="

!    call tdhf_sigma(infos, mo,x, v)



    ! clean up
    deallocate(foo,fvv,work1,work2,work3,dm,v,xmat,x2mat)

  end subroutine calc_h_op

  subroutine rotate_orbs_trah(infos, molgrid ,step, nocc_a, nocc_b, mo, fock_ao)

    use types, only: information
    use basis_tools,      only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use oqp_tagarray_driver
    implicit none


    class(information), intent(inout), target :: infos
    type(dft_grid_t), intent(in) :: molgrid
    real(kind=dp), intent(in)        :: step(:)
    integer, intent(in)              :: nocc_a, nocc_b
    real(kind=dp), intent(inout)     :: mo(:,:)
    real(kind=dp), allocatable       :: work_1(:,:), work_2(:,:)
    real(dp), pointer, intent(inout)        :: fock_ao(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:)
    character(len=*), parameter :: tags_alpha(1) = &
            (/ character(len=80) :: OQP_VEC_MO_A /)
    integer            :: nbf, i, idx
    logical :: second_term
    type(basis_set), pointer :: basis

!    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    basis => infos%basis
    nbf   = basis%nbf
    allocate(work_1(nbf,nbf),work_2(nbf,nbf))
    work_1 = 0
    work_2 = 0
    idx = 0
    print *,"mo in orbital befor rot",mo
    print *, "step", step
    second_term = .true.
    if (infos%control%scftype == 3) then! ROHF
      second_term = .false.
    end if
    call exp_scaling(work_1, step, idx, nocc_a, nocc_b, nbf, second_term)

    print *, "work_1", work_1 
    call orthonormalize(work_1, nbf)
    call dgemm('N','N', nbf, nbf, nbf, 1.0_dp, mo, nbf, work_1, nbf, 0.0_dp, work_2, nbf)
    mo = work_2
!    mo_a = mo
    print *,"mo in orbital after rot",mo
    call get_fock(basis, infos, molgrid, fock_ao, mo)
!    print *, "fock_ao",fock_ao
!   here we need to calculate fock matrix again


  contains

    subroutine exp_scaling(G, step, idx, nocc_a, nocc_b, nbf, second_term)
      use matrix_expm_mod, only: expm
      real(kind=dp), intent(out)     :: G(:,:)
      real(kind=dp), intent(in)      :: step(:)
      integer, intent(inout)         :: idx
      integer, intent(in)            :: nocc_a, nocc_b, nbf

      real(kind=dp), allocatable     :: K(:,:), K2(:,:)
      integer                        :: occ, virt, istart, i
      logical, intent(inout)         :: second_term
      integer :: info

      allocate(K(nbf,nbf), source=0.0_dp)
      allocate(K2(nbf,nbf), source=0.0_dp)

!      do occ = 1, nocc_a
!        istart = merge(nocc_b+1, nocc_a+1, occ <= nocc_b)
!        do virt = istart, nbf
!          idx = idx + 1
!          K(virt, occ) =  step(idx)
!          K(occ, virt) = -step(idx)
!        end do
!      end do

      istart = nocc_a +1 !merge(nocc_b+1, nocc_a+1, occ <= nocc_b)
      do virt = istart, nbf
        do occ = 1, nocc_a
          idx = idx + 1
          K(virt, occ) =  step(idx)
          K(occ, virt) = -step(idx)
        end do
      end do
!    do i = 1, nvir
!      do a = 1, nocc
!        k= k+1
!        x2(k) = 2*x2mat(i,a)
!      end do
!    end do
!      do occ = 1, nocc_a
!        istart = merge(nocc_b+1, nocc_a+1, occ <= nocc_b)
!        do virt = istart, nbf
!          idx = idx + 1
!          K(virt, occ) =  step(idx)
!          K(occ, virt) = -step(idx)
!        end do
!      end do
      
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

  subroutine get_fock(basis, infos, molgrid, fock_ao, mo_in)
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
    real(dp), pointer, intent(in)        :: fock_ao(:)
    real(dp), target, intent(in), optional :: mo_in(:,:)

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
    if (present(mo_in)) then
       mo_a = mo_in
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
    case (scf_uhf)
      fock_a = pfock(:,1)
      fock_b = pfock(:,2)
      dmat_a = pdmat(:,1)
      dmat_b = pdmat(:,2)
    case (scf_rohf)
      fock_a = rohf_bak(:,1)
!      call mo_to_ao(fock_b, pfock(:,2), smat_full, mo_a, nbf, nbf, work1, work2)
      dmat_a = pdmat(:,1) - pdmat(:,2)
      dmat_b = pdmat(:,2)
      mo_b = mo_a
      mo_energy_b = mo_energy_a
    end select

    fock_ao = pfock(:,1)
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

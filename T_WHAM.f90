!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the calculation of the free energy from paralell tempering simulations !
! Written by                                                                        !
!                                        Ye Mei                                     !
!                                      09/06/2014                                   !
!                           East China Normal University                            !
! Reference:                                                                        !
!      J. Chem. Theory Comput., 3, 26-41 (2007)                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module T_WHAM
  use precision_m
  implicit none
  
  integer(kind=4), public :: NumR   ! number of replicas
  integer(kind=4), public :: NumS   ! number of snapshots
  integer(kind=4), public :: NumT   ! number of temperatures
  
  private
  integer(kind=4) :: IndexR ! index of replicas
  integer(kind=4) :: IndexS ! index of snapshot
  integer(kind=4) :: IndexB ! index of energy bin
  integer(kind=4) :: IndexT ! index of temperature

  integer(kind=4), public :: NumB = 100  ! number of bins
  real(kind=fp_kind), parameter :: TOLERANCE = 1.0e-5 ! convergence criterion for free energy 
  integer(kind=4), parameter :: MaxITS = 1000  ! max number of iterations

  real(kind=fp_kind), public :: target_temperature = 300.d0

  real(kind=fp_kind), allocatable, public :: U_nk(:,:)      ! potential energy of snapshot n in replica k
  real(kind=fp_kind), allocatable, public :: T_nk(:,:)      ! temperature of snapshot n in replica 
  integer(kind=4), allocatable, public :: state_nk(:,:)     ! temperature state of snapshot n in replica k
  integer(kind=4), allocatable :: bin_nk(:,:)               ! bin affiliation of snapshot n in replica k
  real(kind=fp_kind), allocatable :: U_m(:)                 ! potential energy at the center of the bin
  real(kind=fp_kind), allocatable :: beta_l(:)              ! inverse temperature
  integer(kind=4),    allocatable :: N_kl(:,:)              ! number of visits to temperature l by replica k
  integer(kind=4),    allocatable :: H_mk(:,:)              ! histogram of energy bin m by replica k
  integer(kind=4),    allocatable :: H_m(:)                 ! histogram of energy bin m by all the replicas
  real(kind=fp_kind), allocatable :: g_mk(:,:)              ! statistical inefficiency of energy bin m from replica k
  real(kind=fp_kind), allocatable :: delta2H_m(:)           ! statistical uncertainty of histogram of energy bin m in all replicas
  real(kind=fp_kind), allocatable :: delta2H_mk(:,:)        ! statistical uncertainty of histogram of energy bin m in replica k
  real(kind=fp_kind), allocatable :: f_l(:)                 ! dimensionless free energy of temperature index l
  real(kind=fp_kind), allocatable :: omega_mk(:,:)          ! density of states for energy bin m in replica k
  real(kind=fp_kind), allocatable :: omega_m(:)             ! density of states for energy bin m in all replicas
  real(kind=fp_kind), allocatable :: delta2omega_m(:)       ! statistical uncertainty of density of states for energy bin m in all replicas
  real(kind=fp_kind), allocatable :: delta2omega_mk(:,:)    ! statistical uncertainty of density of states for energy bin m in replica k
  real(kind=fp_kind), allocatable, public :: weight_nk(:,:) ! weighting factor for snapshot n of replica k

  real(kind=fp_kind) :: deltaU ! width of the energy bin

  integer(kind=4), public :: NumRCbin  ! number of bins for the reaction coordinate

  integer(kind=4), allocatable, public :: iRCbin_nk(:,:) ! RC bin affiliation of snapshot n in replica k

  real(kind=fp_kind), allocatable, public :: pmf(:) ! potential of mean force along the reaction coordinates

  public :: WHAM_initialize, WHAM_finalize, WHAM_computePMF

contains

  subroutine WHAM_initialize()
    use precision_m
    implicit none

    NumT=NumR

    allocate(state_nk(NumS,NumR))
    allocate(bin_nk(NumS,NumR))
    allocate(U_m(NumB))
    allocate(beta_l(NumT))
    allocate(N_kl(NumR,NumT))
    allocate(H_mk(NumB,NumR))
    allocate(H_m(NumB))
    allocate(g_mk(NumB,NumR))
    allocate(delta2H_m(NumB))
    allocate(delta2H_mk(NumB,NumR))
    allocate(f_l(NumT))
    allocate(omega_mk(NumB,NumR))
    allocate(omega_m(NumB))
    allocate(delta2omega_m(NumB))
    allocate(delta2omega_mk(NumB,NumR))
    allocate(weight_nk(NumS,NumR))
    allocate(iRCbin_nk(NumS,NumR))

    call construct_beta_l

    call construct_Um

    call assign_snapshots_to_state

    call construct_Histogram

    call construct_Nkl

    call construct_gmk

    call calculateStatUncertaintyHistogram

    call calculateDensityOfState

    call calculateStatUncertaintyOfDensityOfState

    call calculateWeight

  end subroutine WHAM_initialize

  subroutine WHAM_finalize()
    use precision_m
    implicit none
    if(allocated(U_nk)) deallocate(U_nk)
    if(allocated(T_nk)) deallocate(T_nk)
    if(allocated(state_nk)) deallocate(state_nk)
    if(allocated(bin_nk)) deallocate(bin_nk)
    if(allocated(U_m)) deallocate(U_m)
    if(allocated(beta_l)) deallocate(beta_l)
    if(allocated(N_kl)) deallocate(N_kl)
    if(allocated(H_mk)) deallocate(H_mk)
    if(allocated(H_m)) deallocate(H_m)
    if(allocated(g_mk)) deallocate(g_mk)
    if(allocated(delta2H_m)) deallocate(delta2H_m)
    if(allocated(delta2H_mk)) deallocate(delta2H_mk)
    if(allocated(f_l)) deallocate(f_l)
    if(allocated(omega_mk)) deallocate(omega_mk)
    if(allocated(omega_m)) deallocate(omega_m)
    if(allocated(delta2omega_m)) deallocate(delta2omega_m)
    if(allocated(delta2omega_mk)) deallocate(delta2omega_mk)
    if(allocated(weight_nk)) deallocate(weight_nk)
    if(allocated(iRCbin_nk)) deallocate(iRCbin_nk)
  end subroutine WHAM_finalize
  
  subroutine WHAM_computePMF
    use precision_m
    use constant
    implicit none
    pmf=0.0d0
    do IndexR = 1, NumR
      do IndexS = 1, NumS
        pmf(iRCbin_nk(IndexS,IndexR)) = pmf(iRCbin_nk(IndexS,IndexR)) + weight_nk(IndexS,IndexR)
      end do
    end do 
    pmf = -kB*target_temperature*log(pmf + 1.D-7)
    pmf = pmf - minval(pmf)
  end subroutine WHAM_computePMF

  subroutine construct_beta_l
    use precision_m
    use constant
    implicit none
    beta_l = T_nk(1,:)
    call sort(NumT, beta_l)
    beta_l = 1 / (beta_l*kB)
  end subroutine construct_beta_l

  subroutine construct_Um
    use precision_m
    implicit none
    real(kind=fp_kind) :: Umax, Umin
    real(kind=fp_kind) :: U_shift
    integer(kind=4) :: i, j, k
    U_shift = minval(U_nk)
    U_nk = U_nk - U_shift
    Umax = maxval(U_nk)
    Umin = minval(U_nk)
    deltaU = ( Umax - Umin ) / ( NumB - 1 )   
    forall(i=1:NumB)
      U_m(i) = Umin + deltaU * (i-1) 
    end forall
  end subroutine construct_Um

  subroutine assign_snapshots_to_state
    use precision_m
    implicit none
    integer(kind=4) :: ilocation
    real(kind=fp_kind) :: Tempr(NumT)
    Tempr = T_nk(1,:)
    call sort(NumT, Tempr)
    state_nk = 0
    do IndexS = 1, NumS
      do IndexR = 1, NumR
        do ilocation = 1, NumT
          if(abs(T_nk(IndexS,IndexR)-Tempr(ilocation))<1.D-6)then
            state_nk(IndexS,IndexR) = ilocation
          end if
        end do
      end do
!      write(99,'(<NumR>I4)')state_nk(IndexS,:)
    end do
  end subroutine assign_snapshots_to_state

  subroutine construct_Histogram
    use precision_m
    implicit none
    integer(kind=4) :: binID
    H_mk = 0
    do IndexR = 1, NumR
      do IndexS = 1, NumS
        binID = floor((U_nk(IndexS,IndexR)-(U_m(1)-deltaU/2.d0))/deltaU) + 1
        bin_nk(IndexS,IndexR) = binID
        H_mk(binID, IndexR) = H_mk(binID, IndexR) + 1
      end do
    end do
!   write(100,'(<NumR>I8)')((H_mk(IndexB,IndexR),IndexR=1,NumR), IndexB=1,NumB)
    forall (IndexB=1:NumB)
      H_m(IndexB)=sum(H_mk(IndexB,:))
    end forall
!    write(33,'(<NumB>I8)')H_m
  end subroutine construct_Histogram

  subroutine construct_Nkl
    use precision_m
    implicit none
    N_kl = 0
    do IndexR = 1, NumR
      do IndexS = 1, NumS
        N_kl(IndexR,state_nk(IndexS,IndexR)) = N_kl(IndexR,state_nk(IndexS,IndexR)) + 1
      end do
    end do
!    write(22,'(<NumT>I8)')N_kl
  end subroutine construct_Nkl

  subroutine construct_gmk
    use precision_m
    implicit none
    real(kind=fp_kind) :: psi_avg(NumB,NumR) ! average of the occupacy in bin m by snapshots in replica K.
                                             ! psi_avg = psi^2_avg 
    integer(kind=4) :: dt
    integer(kind=4) :: ibin1, ibin2
    integer(kind=4) :: C_mk(NumB,NumR)
    integer(kind=4) :: ibinarray(NumS)
    real(kind=fp_kind) :: corr(NumS)
    integer(kind=4) :: increment
    print*,' Calculating statistical inefficiency'
    psi_avg = 0.d0
    do IndexR = 1, NumR
      do IndexB = 1, NumB
        do IndexS = 1, NumS
          if(bin_nk(IndexS,IndexR) == IndexB) psi_avg(IndexB,IndexR) = psi_avg(IndexB,IndexR) + 1
        end do
      end do
    end do 
    psi_avg = psi_avg / NumS

    g_mk=1.d0
    dt = 1
    increment = 1
    do while (dt < NumS-1)
      print*, 'dT = ', dt
      C_mk = 0.d0
      do IndexR = 1, NumR
        do IndexS = 1, NumS - dt
          ibin1 = bin_nk(IndexS     , IndexR)
          ibin2 = bin_nk(IndexS + dt, IndexR)
          if(ibin1 == ibin2) C_mk(ibin1,IndexR) = C_mk(ibin1,IndexR) + 1
        end do
      end do
      C_mk = C_mk / (NumS - dt)
      where(psi_avg > 1.D-8) C_mk = (C_mk - psi_avg**2) / (psi_avg - psi_avg**2)
      g_mk = g_mk + 2.d0 * C_mk * (1.d0 - real(dt)/real(NumS)) * increment
      if( all(C_mk < 0.d0) ) exit

      dt = dt + increment
      increment = increment + 1
    end do
!    write(*,'(<NumR>F8.3)')g_mk
  end subroutine construct_gmk

  subroutine calculateStatUncertaintyHistogram
    use precision_m
    implicit none
    real(kind=fp_kind) :: H_mk_mean
    do IndexB = 1, NumB
      do IndexR = 1, NumR
        ! TO BE UPDATED
      end do
    end do
  end subroutine calculateStatUncertaintyHistogram

  subroutine calculateDensityOfState
    use precision_m
    implicit none
    real(kind=fp_kind) :: numerator(NumR), denominator(NumR)
    integer(kind=4) :: converged
    real(kind=fp_kind) :: f_l_old(NumT)
    real(kind=fp_kind) :: rmsd
    integer(kind=4) :: iITS
    converged = 0
    f_l=0.d0
    iITS = 0 
    do while (converged == 0)
      iITS = iITS + 1
      if(iITS >MaxITS) then
        write(*,*)'Error: Maximal number of iterations has been reached.'
        write(*,*)'Quit.'
        stop
      end if

      f_l_old=f_l

      do IndexB = 1, NumB
        do IndexR = 1, NumR
          numerator(IndexR) = H_mk(IndexB, IndexR)
          denominator(IndexR) = 0.d0
          do IndexT = 1, NumT
            denominator(IndexR) = denominator(IndexR) + N_kl(IndexR,IndexT)*deltaU*exp(f_l(IndexT)-U_m(IndexB)*beta_l(IndexT))
          end do
          omega_mk(IndexB,IndexR) = numerator(IndexR)/denominator(IndexR)
        end do
        omega_m(IndexB) = sum(numerator(:)/g_mk(IndexB,:)) / sum(denominator(:)/g_mk(IndexB,:))
      end do
   
      f_l = 0.d0
      do IndexT = 1, NumT
        do IndexB = 1, NumB
          f_l(IndexT) = f_l(IndexT) + omega_m(IndexB)*deltaU*exp(-U_m(IndexB)*beta_l(IndexT))
        end do
        f_l(IndexT) = - log(f_l(IndexT))
      end do

      f_l = f_l - minval(f_l)

      converged=1
      do IndexT = 1, NumT
        if(abs(f_l(IndexT)-f_l_old(IndexT))>TOLERANCE)converged=0
      end do
      write(*,'(A,I5,A,ES10.3)')'Iteration ',iITS, ': RMSD of dimensionless free energy:', rmsd(NumT, f_l, f_l_old)
    end do
!    write(44,'(<NumB>F10.7)')omega_m
!    write(44,'(F10.7)')sum(omega_m)
  end subroutine calculateDensityOfState

  subroutine calculateStatUncertaintyOfDensityOfState
    use precision_m
    real(kind=fp_kind) :: denominator(NumR)
    do IndexB = 1, NumB
      denominator = 0.d0
      do IndexR = 1, NumR
        do IndexT = 1, NumT
          denominator(IndexR) = denominator(IndexR) + N_kl(IndexR,IndexT)*deltaU*exp(f_l(IndexT)-beta_l(IndexT)*U_m(IndexB))
        end do
        denominator(IndexR) = denominator(IndexR)/g_mk(IndexB,IndexR)
        delta2omega_mk(IndexB,IndexR)=omega_m(IndexB)/denominator(IndexR)
      end do
      delta2omega_m(IndexB)=omega_m(IndexB)/sum(denominator(:))
    end do   
  end subroutine calculateStatUncertaintyOfDensityOfState

  subroutine calculateWeight
    use precision_m
    use constant
    implicit none
    real(kind=fp_kind) :: beta_target
    real(kind=fp_kind) :: weight_sum
    beta_target = 1.d0 / (kB*target_temperature)
    weight_nk = 0.d0
    do IndexS = 1, NumS
      do IndexR = 1, NumR
        do IndexB = 1, NumB
          if(bin_nk(IndexS,IndexR) == IndexB)then
            weight_nk(IndexS,IndexR) = weight_nk(IndexS,IndexR) + omega_m(IndexB) * exp(-beta_target*U_m(IndexB)) / (H_m(IndexB)+1.D-9)
          end if
!          write(54,'(2E12.5,I8,E12.5)')omega_m(IndexB), U_m(IndexB), H_m(IndexB), exp(-beta_target*U_m(IndexB))
        end do
      end do
    end do
    weight_sum = sum(weight_nk)
    weight_nk = weight_nk / weight_sum
!    write(66,'(<NumR>E12.3)')((weight_nk(IndexS,IndexR), IndexR=1,NumR), IndexS=1,NumS)
!    print*,'sum of weight_nk: ', sum(weight_nk)
  end subroutine calculateWeight
end module T_WHAM


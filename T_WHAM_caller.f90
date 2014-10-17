program WHAMcaller
  use precision_m
  use T_WHAM
  implicit none
interface
  subroutine GetSimulData(nreplicas,nsamples,pote,temperature)
    use precision_m
    implicit none
    integer(kind=4), intent(out) :: nreplicas
    integer(kind=4), intent(out) :: nsamples
    real(kind=fp_kind), allocatable, intent(out) :: pote(:,:)
    real(kind=fp_kind), allocatable, intent(out) :: temperature(:,:)
  end subroutine GetSimulData

end interface
  real(kind=fp_kind), allocatable :: rc1(:,:)
  real(kind=fp_kind), allocatable :: rc2(:,:)
  integer(kind=4) :: iostatus
  integer(kind=4) :: IndexS, IndexR
  print*, 'get simulation data'
  call GetSimulData(NumR,NumS,U_nk,T_nk)

! you can customize NumB before WHAM_initialize
!  NumB = 200

  call WHAM_initialize
  allocate(rc1(NumS,NumR))
  allocate(rc2(NumS,NumR))
  print*,'Reading reaction coordinate 1'
  open(11,file='rmsd.dat',status='OLD',iostat=iostatus)
  read(11,'(<NumR>F10.6)')((rc1(IndexS,IndexR),IndexR=1,NumR),IndexS=1,NumS)
  close(11)
  print*,'Reading reaction coordinate 2'
  open(11,file='rg.dat',status='OLD',iostat=iostatus)
  read(11,'(<NumR>F10.6)')((rc2(IndexS,IndexR),IndexR=1,NumR),IndexS=1,NumS)
  close(11)
  call getpmf_2D(rc1,rc2)
  call WHAM_finalize
  deallocate(rc1,rc2)

end program WHAMcaller

subroutine getpmf_2D(rc1,rc2)
  use precision_m
  use constant
  use T_WHAM
  implicit none
  real(kind=fp_kind) :: rc1(NumS,NumR), rc2(NumS,NumR)
  real(kind=fp_kind) :: rc1min, rc1max
  real(kind=fp_kind) :: rc2min, rc2max
  integer(kind=4) :: NumBinRC1=201, NumBinRC2=201
  real(kind=fp_kind) :: deltaRC1, deltaRC2
  integer(kind=4) :: RC1binID, RC2binID, RCID
  integer(kind=4) :: IndexS, IndexR
  real(kind=fp_kind), allocatable :: pmf_2D(:,:)
  rc1min = floor(minval(rc1))
  rc1max = ceiling(maxval(rc1))
  rc2min = floor(minval(rc2))
  rc2max = ceiling(maxval(rc2))
  deltaRC1 = ( rc1max - rc1min) / ( NumBinRC1 - 1 )
  deltaRC2 = ( rc2max - rc2min) / ( NumBinRC2 - 1 )
  print*,'RC1 :', rc1min, rc1max, deltaRC1
  print*,'RC2 :', rc2min, rc2max, deltaRC2
  print*,'Calculating histogram of (RC1, RC2)'
  NumRCbin = NumBinRC1*NumBinRC2
  do IndexR = 1, NumR
    do IndexS = 1, NumS
      RC1binID = floor( ( rc1( IndexS, IndexR ) - ( rc1min - deltaRC1/2.d0 ) ) / deltaRC1 ) + 1
      RC2binID = floor( ( rc2( IndexS, IndexR ) - ( rc2min - deltaRC2/2.d0 ) ) / deltaRC2 ) + 1
      RCID = ( RC2binID - 1 ) * NumBinRC1 + RC1binID
      iRCbin_nk( IndexS, IndexR ) = RCID
    end do
  end do

  allocate(pmf_2D(NumBinRC1,NumBinRC2))
  allocate(pmf(NumRCbin))
  call WHAM_computePMF
  pmf_2D = reshape(pmf, (/NumBinRC1, NumBinRC2/)) 
  deallocate(pmf)
 
  do RC1binID = 1, NumBinRC1
    do RC2binID = 1, NumBinRC2
      write(88,'(3F10.4)')rc1min+RC1binID*deltaRC1, rc2min+RC2binID*deltaRC2, pmf_2D(RC1binID,RC2binID)
    end do
    write(88,*)
  end do
  deallocate(pmf_2D)
end subroutine getpmf_2D

subroutine GetSimulData(nreplicas,nsamples,pote,temperature)
  use precision_m
  implicit none
  integer(kind=4), intent(out) :: nreplicas
  integer(kind=4), intent(out) :: nsamples
  real(kind=fp_kind), allocatable, intent(out) :: pote(:,:)
  real(kind=fp_kind), allocatable, intent(out) :: temperature(:,:)
  character(len=120) :: datafile
  integer(kind=4) :: iostatus
  integer(kind=4) :: i, j, k, m, n
  read*,datafile
  open(11,file=datafile,status='OLD',iostat=iostatus)
  read(11,*)nreplicas,nsamples
  print*,'Number of replicas:', nreplicas
  print*,'Number of samples:', nsamples
  if(allocated(pote))deallocate(pote)
  allocate(pote(nsamples,nreplicas))
  if(allocated(temperature))deallocate(temperature)
  allocate(temperature(nsamples,nreplicas))
  print*,'Reading potential energy from data file'
  read(11,'(<nreplicas>F8.2)')((pote(i,j),j=1,nreplicas),i=1,nsamples)
  read(11,*)
  print*,'Reading temperature information from data file'
  read(11,'(<nreplicas>F8.2)')((temperature(i,j),j=1,nreplicas),i=1,nsamples)
  close(11)
end subroutine GetSimulData

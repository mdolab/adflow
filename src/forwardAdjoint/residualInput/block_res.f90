! This is a super-combined function that combines the original
! functionality of: 
! Pressure Computation
! timeStep
! applyAllBCs
! initRes
! residual 

! The real difference between this and the original modules is that it
! it only operates on a single block at a time and as such the nominal
! block/sps loop is outside the calculation. This routine is suitable
! for forward mode AD with Tapenade

subroutine block_res(nn, sps, useSpatial, useExtra) 

  !alpha, beta, liftIndex,

  use blockPointers       
  use flowVarRefState     
  use inputPhysics 
  use inputTimeSpectral
  use section
  use monitor
  use iteration
  use NKSolverVars, only : resSUm
  implicit none

  ! Input Arguments:
  integer(kind=intType), intent(in) :: nn, sps
  !real(kind=realType), intent(in) :: alpha, beta
  !integer(kind=intType), intent(in) :: liftIndex
  logical, intent(in) :: useSpatial, useExtra

  ! Working Variables
  real(kind=realType) :: gm1,v2
  integer(kind=intType) :: i,j,k,sps2,mm,l
  real(kind=realType), dimension(nSections) :: t
  logical :: useOldCoor

  useOldCoor = .False.

  ! Set pointers to input/output variables
  w  => flowDoms(nn,1,sps)%w
  dw => flowDoms(nn,1,sps)%dw
  x  => flowDoms(nn,1,sps)%x

  ! ------------------------------------------------
  !        Additional 'Extra' Components
  ! ------------------------------------------------

  if (useExtra) then
!      call adjustInflowAngle(alpha,beta,liftIndex)
!      call referenceState_mod()
!      call setFlowInfinityState()
  end if

  ! ------------------------------------------------
  !        Additional Spatial Components
  ! ------------------------------------------------

  if (useSpatial .or. useExtra) then

     call xhalo_block(nn,1,sps)
     call metric_block(nn,1,sps)

     t = timeUnsteadyRestart
  
     if(equationMode == timeSpectral) then
        do mm=1,nSections
           t(mm) = t(mm) + (sps-1)*sections(mm)%timePeriod &
                /         real(nTimeIntervalsSpectral,realType)
        enddo
     endif

     !call gridVelocitiesFineLevel_block(useOldCoor, t, sps) ! Required for TS
     !call normalVelocities_block(sps) ! Required for TS
     !call slipVelocitiesFineLevel(.false., t, mm) !required for wall Functions
  end if

  ! ------------------------------------------------
  !        Normal Residual Computation
  ! ------------------------------------------------

  ! Compute the pressures
  gm1 = gammaConstant - one

  ! Compute P 
  do k=0,kb
     do j=0,jb
        do i=0,ib
           v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 + w(i,j,k,ivz)**2
           p(i,j,k) = gm1*(w(i,j,k,irhoE) - half*w(i,j,k,irho)*v2)
           p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
        enddo
     enddo
  enddo

  ! Compute Laminar/eddy viscosity if required
  call computeLamViscosity
  !call computeEddyViscosity # Required for turblence models

  !  Apply all BC's
  call applyAllBC_block(.True.)
  
  ! Compute skin_friction Velocity
  !call computeUtau_block ! Required for turblence
  
  ! Compute time step and spectral radius
  call timeStep_block(.false.)
  
  if( equations == RANSEquations ) then
     ! Initialize only the Turblent Variables
     call initres_block(nt1MG, nMGVar,nn,sps) 
     call turbResidual_block
  endif
  
  ! Next initialize residual for flow variables. The is the only place
  ! where there is an n^2 dependance
  do sps2 = 1,nTimeIntervalsSpectral
     dw => flowDoms(nn,1,sps2)%dw
     call initRes_block(1,nwf, nn, sps2)
  end do

  ! Reset dw pointer to sps instance
  dw => flowDoms(nn,1,sps)%dw

  !  Actual residual calc
  call residual_block

  ! Divide through by the volume
  do sps2 = 1,nTimeIntervalsSpectral
     ! Set dw and vol to looping sps2 instance
     dw => flowDoms(nn,1,sps2)%dw
     vol => flowDoms(nn,currentLevel,sps2)%vol

     do l=1,nw
        do k=2,kl
           do j=2,jl
              do i=2,il
                 dw(i,j,k,l) = dw(i,j,k,l) / vol(i,j,k) * resSum(l)
              end do
           end do
        end do
     end do
  end do
  
  ! Reset dw and vol to sps instance
  dw => flowDoms(nn,1,sps)%dw
  vol => flowDoms(nn,currentLevel,sps)%vol

end subroutine block_res

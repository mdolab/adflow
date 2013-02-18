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

subroutine block_res(nn, sps, useSpatial, useForces, &
     alpha, beta, liftIndex, &
     Force, Moment, Lift, Drag, cForce, cMoment, CL, CD)

  use blockPointers       
  use flowVarRefState     
  use inputPhysics 
  use inputTimeSpectral
  use section
  use monitor
  use iteration
  use costFunctions
  implicit none

  ! Input Arguments:
  integer(kind=intType), intent(in) :: nn, sps
  logical, intent(in) :: useSpatial, useForces
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) :: liftIndex

  ! Output Arguments: Note: Cannot put intent(out) since reverse mode
  ! may NOT compute these values and then compilation will fail
  real(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
  real(kind=realType) :: Lift, Drag, CL, CD

  ! Working Variables
  real(kind=realType) :: gm1, v2, fact
  integer(kind=intType) :: i, j, k, sps2, mm, l
  real(kind=realType), dimension(nSections) :: t
  real(kind=realType), dimension(3) :: cFp, cFv, cMp, cMv
  real(kind=realType) :: yplusMax, scaleDim
  logical :: useOldCoor
  
  useOldCoor = .False.

  ! Set pointers to input/output variables
  w  => flowDoms(nn, currentLevel, sps)%w
  dw => flowDoms(nn, 1           , sps)%dw
  x  => flowDoms(nn, currentLevel, sps)%x

  ! ------------------------------------------------
  !        Additional 'Extra' Components
  ! ------------------------------------------------ 

  call adjustInflowAngle(alpha, beta, liftIndex)
  call referenceState
  call setFlowInfinityState

  ! ------------------------------------------------
  !        Additional Spatial Components
  ! ------------------------------------------------
  if (useSpatial) then

     call xhalo_block
     call metric_block
     ! -------------------------------------
     ! These functions are required for TS
     ! --------------------------------------

     t = timeUnsteadyRestart
     if(equationMode == timeSpectral) then
        do mm=1,nSections
           t(mm) = t(mm) + (sps-1)*sections(mm)%timePeriod &
                /         real(nTimeIntervalsSpectral,realType)
        enddo
     endif
     
     call gridVelocitiesFineLevel_block(useOldCoor, t, sps) ! Required for TS
     call normalVelocities_block(sps) ! Required for TS
  end if
  
  ! ------------------------------------------------
  !        Normal Residual Computation
  ! ------------------------------------------------

  ! Compute the pressures
  gm1 = gammaConstant - one

  ! Compute P 
  do k=0, kb
     do j=0, jb
        do i=0, ib
           v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
           p(i, j, k) = gm1*(w(i, j, k, irhoE) - half*w( i, j, k, irho)*v2)
           p(i, j, k) = max(p(i, j, k), 1.e-4_realType*pInfCorr)
        enddo
     enddo
  enddo

  ! Compute Laminar/eddy viscosity if required
  call computeLamViscosity
  call computeEddyViscosity 

  !  Apply all BC's
#ifndef TAPENADE_REVERSE
  call applyAllBC_block(.True.)
#endif  
  ! Compute skin_friction Velocity (only for wall Functions)
#ifndef TAPENADE_REVERSE
  call computeUtau_block
#endif
  ! Compute time step and spectral radius
  call timeStep_block(.false.)

  ! -------------------------------
  ! The forward ADjoint is NOT currently setup for RANS equations
  !   if( equations == RANSEquations ) then
  !      ! Initialize only the Turblent Variables
  !      call initres_block(nt1MG, nMGVar,nn,sps) 
  !      call turbResidual_block
  !   endif
  ! -------------------------------  

  ! -------------------------------

  ! Next initialize residual for flow variables. The is the only place
  ! where there is an n^2 dependance

  ! sps here is the on-spectral instance
  call initRes_block(1, nwf, nn, sps)
  
  !  Actual residual calc
  call residual_block

  ! Divide through by the volume
  do sps2 = 1,nTimeIntervalsSpectral
     ! Set dw and vol to looping sps2 instance
     dw => flowDoms(nn, 1, sps2)%dw
     vol => flowDoms(nn, currentLevel, sps2)%vol

     do l=1, nw
        do k=2, kl
           do j=2, jl
              do i=2, il
                 dw(i, j, k, l) = dw(i, j, k, l) / vol(i, j, k) 
              end do
           end do
        end do
     end do
  end do
  
  ! Reset dw and vol to sps instance
  dw => flowDoms(nn,1,sps)%dw
  vol => flowDoms(nn,currentLevel,sps)%vol

  ! We are now done with the residuals, we move on to the forces and moments

  ! This routine compute Force, Moment, Lift, Drag, and the
  ! coefficients of the values

  if (useForces) then
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax)
     scaleDim = pRef/pInf

     ! Sum pressure and viscous contributions
     cForce = cFp + cFv
     cMoment = cMp + cMv
     
     ! Get Lift coef and Drag coef
     CD =  cForce(1)*dragDirection(1) &
          + cForce(2)*dragDirection(2) &
          + cForce(3)*dragDirection(3)
     
     CL =  cForce(1)*liftDirection(1) &
          + cForce(2)*liftDirection(2) &
          + cForce(3)*liftDirection(3)
     
     ! Divide by fact to get the forces, Lift and Drag back
     
     fact = two/(gammaInf*pInf*MachCoef*MachCoef &
          *surfaceRef*LRef*LRef*scaleDim)
     Force = cForce / fact
     Lift  = CL / fact
     Drag  = CD / fact
     
     ! Moment factor has an extra lengthRef
     fact = fact/(lengthRef*LRef)
     
     Moment = cMoment / fact
  else
     Force = zero
     Moment = zero
     cForce = zero
     cMoment = zero
     Lift = zero
     Drag = zero
     CD = zero
     CD = zero
  end if

  call getCostFuncMat(alpha, beta, liftIndex)

end subroutine block_res

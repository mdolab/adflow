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
  use inputADjoint
  use diffSizes
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
  integer(kind=intType) :: i, j, k, sps2, mm, l, ii, ll, jj, lEnd
  integer(kind=intType) :: nState
  real(kind=realType), dimension(nSections) :: t
  real(kind=realType), dimension(3) :: cFp, cFv, cMp, cMv
  real(kind=realType) :: yplusMax, scaleDim, tmp
  logical :: useOldCoor
  real(kind=realType), pointer, dimension(:,:,:,:) :: wsp
  real(kind=realType), pointer, dimension(:,:,:) :: volsp
  useOldCoor = .False.

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  ! Set pointers to input/output variables
  w  => flowDoms(nn, currentLevel, sps)%w
  dw => flowDoms(nn, 1           , sps)%dw
  x  => flowDoms(nn, currentLevel, sps)%x
  vol=> flowDoms(nn, currentLevel, sps)%vol 

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
  call applyAllBC_block(.True.)

  ! Compute skin_friction Velocity (only for wall Functions)
! #ifndef TAPENADE_REVERSE
!   call computeUtau_block
! #endif
  ! Compute time step and spectral radius
  call timeStep_block(.false.)

  ! -------------------------------
  ! Compute turbulence residual for RANS equations
  if( equations == RANSEquations) then
     ! Initialize only the Turblent Variables
     call unsteadyTurbSpectral_block(itu1, itu2, nn, sps)
     
     select case (turbModel)
        
     case (spalartAllmaras)
        !call determineDistance2(1, sps)
        call sa_block(.true.)
        
     case default
        call terminate("turbResidual", & 
             "Only SA turbulence adjoint implemented")
        
     end select
     
  endif

  ! -------------------------------  

  ! Next initialize residual for flow variables. The is the only place
  ! where there is an n^2 dependance. There are issues with
  ! initRes. So only the necesary timespectral code has been copied
  ! here. See initres for more information and comments.

  ! sps here is the on-spectral instance
  if (nTimeIntervalsSpectral == 1) then
     dw(:,:,:,1:nwf) = zero
  else
     ! Zero dw on all spectral instances
     spectralLoop1: do sps2=1,nTimeIntervalsSpectral
        flowDoms(nn, 1, sps2)%dw(:,:,:,1:nwf) = zero
     end do spectralLoop1

     spectralLoop2: do sps2=1,nTimeIntervalsSpectral
        jj = sectionID    

        timeLoopFine: do mm=1,nTimeIntervalsSpectral
           ii    =  3*(mm-1)

           varLoopFine: do l=1,nwf
              if(l == ivx .or. l == ivy .or. l == ivz) then
                 if(l == ivx) ll = 3*sps2 - 2
                 if(l == ivy) ll = 3*sps2 - 1
                 if(l == ivz) ll = 3*sps2
                 do k=2,kl
                    do j=2,jl
                       do i=2,il
                          tmp = dvector(jj,ll,ii+1)*flowDoms(nn,1,mm)%w(i,j,k,ivx) &
                               + dvector(jj,ll,ii+2)*flowDoms(nn,1,mm)%w(i,j,k,ivy) &
                               + dvector(jj,ll,ii+3)*flowDoms(nn,1,mm)%w(i,j,k,ivz)
                          flowDoms(nn, 1, sps2)%dw(i,j,k,l) = &
                               flowDoms(nn, 1, sps2)%dw(i,j,k,l) + &
                               tmp*flowDoms(nn,1,mm)%vol(i,j,k)*&
                               flowDoms(nn,1,mm)%w(i,j,k,irho)
                       enddo
                    enddo
                 enddo
              else
                 do k=2,kl
                    do j=2,jl
                       do i=2,il
                          ! This is: dw = dw + dscalar*vol*w
                          flowDoms(nn, 1, sps2)%dw(i,j,k,l) = &
                               flowDoms(nn, 1, sps2)%dw(i,j,k,l) + &
                               dscalar(jj,sps2,mm)*flowDoms(nn,1,mm)%vol(i,j,k)*&
                               flowDoms(nn,1,mm)%w(i,j,k,l) 
                       enddo
                    enddo
                 enddo
              endif
           end do varLoopFine
        end do timeLoopFine
     end do spectralLoop2
  end if

  !  Actual residual calc
  call residual_block

  ! Divide through by the volume
  do sps2 = 1,nTimeIntervalsSpectral
     do l=1, nState
        do k=2, kl
           do j=2, jl
              do i=2, il
                 flowDoms(nn, 1, sps2)%dw(i, j, k, l)  = & 
                      flowDoms(nn, 1, sps2)%dw(i, j, k, l)  / &
                      flowDoms(nn, currentLevel, sps2)%vol(i, j, k)
              end do
           end do
        end do
     end do
  end do

  ! We are now done with the residuals, we move on to the forces and
  ! moments

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

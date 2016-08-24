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

subroutine block_res(nn, sps, useSpatial, alpha, beta, liftIndex, &
     frozenTurb)
  use BCRoutines
  use blockPointers       
  use flowVarRefState     
  use inputPhysics 
  use inputIteration
  use inputTimeSpectral
  use section
  use monitor
  use iteration
  use diffSizes
  use costFunctions
  use wallDistanceData
  use inputDiscretization 
  use saModule
  use inputUnsteady
  use utils, only : terminate
  implicit none

  ! Input Arguments:
  integer(kind=intType), intent(in) :: nn, sps
  logical, intent(in) :: useSpatial, frozenTurb
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) :: liftIndex

  ! Output Variables
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment
  real(kind=realType) :: sepSensor, Cavitation, sepSensorAvg(3)
  
  ! Working Variables
  real(kind=realType) :: gm1, v2, fact, tmp
  integer(kind=intType) :: i, j, k, sps2, mm, l, ii, ll, jj, m
  integer(kind=intType) :: nState
  real(kind=realType), dimension(nSections) :: t
  logical :: useOldCoor
  real(kind=realType), dimension(3) :: cFp, cFv, cMp, cMv
  real(kind=realType) :: yplusMax, scaleDim, oneOverDt
  useOldCoor = .False.

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurb ) then
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

  ! ------------------------------------------------
  !        Additional Spatial Components
  ! ------------------------------------------------
  if (useSpatial) then

     call volume_block
     call metric_block
     call boundaryNormals

#ifdef TAPENADE_REVERSE
     if (equations == RANSEquations .and. useApproxWallDistance) then 
        call updateWallDistancesQuickly(nn, 1, sps)
     end if
#endif

#ifndef TAPENADE_REVERSE
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
     call slipVelocitiesFineLevel_block(useOldCoor, t, sps)
#endif
  end if

  ! ------------------------------------------------
  !        Normal Residual Computation
  ! ------------------------------------------------

  ! Compute the pressures
  call computePressureSimple

  ! Compute Laminar/eddy viscosity if required
  call computeLamViscosity
  call computeEddyViscosity 
  
  call applyAllBC_block(.True.)

  if (equations == RANSequations) then 
     call bcTurbTreatment
     call applyAllTurbBCThisBLock(.True.)
  end if

  ! Compute skin_friction Velocity (only for wall Functions)
! #ifndef TAPENADE_REVERSE
!   call computeUtau_block
! #endif
  ! Compute time step and spectral radius
  call timeStep_block(.false.)

  spectralLoop0: do sps2=1,nTimeIntervalsSpectral
     flowDoms(nn, 1, sps2)%dw(:,:,:,:) = zero
  end do spectralLoop0
  
  ! -------------------------------
  ! Compute turbulence residual for RANS equations
  if( equations == RANSEquations) then

     ! ! Initialize only the Turblent Variables
     ! call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

     select case (turbModel)
        
     case (spalartAllmaras)
        call sa_block(.true.)
     !case (menterSST)
        ! Not implemented yet
        !call SST_block(.true.)
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

  !call initres_block(1, nwf, nn, sps)

  if (equationMode == Steady) then 
     dw(:,:,:,1:nwf) = zero
  else if (equationMode == timeSpectral) then 
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
  else if (equationMode == unsteady) then 
     
     ! Assume only MD or BDF types

     ! Store the inverse of the physical nonDimensional
     ! time step a bit easier.
     
     oneOverDt = timeRef/deltaT

     ! Ground level of the multigrid cycle. Initialize the
     ! owned cells to the unsteady source term. First the
     ! term for the current time level. Note that in w the
     ! velocities are stored and not the momentum variables.
     ! Therefore the if-statement is present to correct this.

     do l=1,nw
        if(l == ivx .or. l == ivy .or. l == ivz) then
           ! Momentum variables.
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    flowDoms(nn, 1, sps)%dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                         * w(i,j,k,l)*w(i,j,k,irho)
                 enddo
              enddo
           enddo
        else
           ! Non-momentum variables, for which the variable
           ! to be solved is stored; for the flow equations this
           ! is the conservative variable, for the turbulent
           ! equations the primitive variable.
           
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    flowDoms(nn, 1, sps)%dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                         * w(i,j,k,l)
                 enddo
              enddo
           enddo
        end if
     end do

     ! The terms from the older time levels. Here the
     ! conservative variables are stored. In case of a
     ! deforming mesh, also the old volumes must be taken.
     
     deformingTest: if( deforming_Grid ) then
        
        ! Mesh is deforming and thus the volumes can change.
        ! Use the old volumes as well.
        
        do m=1,nOldLevels
           do l=1,nw
              do k=2,kl
                 do j=2,jl
                    do i=2,il
                       flowDoms(nn, 1, sps)%dw(i,j,k,l) = flowDoms(nn, 1, sps)%dw(i,j,k,l)                 &
                            + coefTime(m)*volOld(m,i,j,k) &
                            * wOld(m,i,j,k,l)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     else deformingTest
        ! Rigid mesh. The volumes remain constant.
           
        do m=1,nOldLevels
           do l=1,nw
              do k=2,kl
                 do j=2,jl
                    do i=2,il
                       flowDoms(nn, 1, sps)%dw(i,j,k,l) = flowDoms(nn, 1, sps)%dw(i,j,k,l)            &
                            + coefTime(m)*vol(i,j,k) &
                            * wOld(m,i,j,k,l)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     end if deformingTest

     ! Multiply the time derivative by the inverse of the
     ! time step to obtain the true time derivative.
     ! This is done after the summation has been done, because
     ! otherwise you run into finite accuracy problems for
     ! very small time steps.
        
     do l=1,nw
        do k=2,kl
           do j=2,jl
              do i=2,il
                 flowDoms(nn, 1, sps)%dw(i,j,k,l) = oneOverDt*flowDoms(nn, 1, sps)%dw(i,j,k,l)
              enddo
           enddo
        enddo
     enddo
  end if

  !  Actual residual calc
  call residual_block

  ! Divide through by the reference volume
  do sps2 = 1,nTimeIntervalsSpectral
     do l=1, nwf
        do k=2, kl
           do j=2, jl
              do i=2, il
                 flowDoms(nn, 1, sps2)%dw(i, j, k, l)  = & 
                      flowDoms(nn, 1, sps2)%dw(i, j, k, l)  / &
                      flowDoms(nn, currentLevel, sps2)%volref(i, j, k)
              end do
           end do
        end do
     end do

     ! Treat the turblent residual with the scaling factor on the
     ! residual
     do l=nt1,nState
        do k=2, kl
           do j=2, jl
              do i=2, il
                 flowDoms(nn, 1, sps2)%dw(i, j, k, l)  = & 
                      flowDoms(nn, 1, sps2)%dw(i, j, k, l)  / &
                      flowDoms(nn, currentLevel, sps2)%volref(i, j, k)*turbResScale(l-nt1+1)
              end do
           end do
        end do
     end do
  end do

  call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, &
       sepSensorAvg, Cavitation)

  ! Convert back to actual forces. Note that even though we use
  ! MachCoef, Lref, and surfaceRef here, they are NOT differented,
  ! since F doesn't actually depend on them. Ideally we would just get
  ! the raw forces and moment form forcesAndMoments. 
  force = zero
  moment = zero
  scaleDim = pRef/pInf
  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
       *surfaceRef*LRef*LRef*scaleDim)
  do sps2 = 1,nTimeIntervalsSpectral
     force(:, sps2) = (cFp + cFV)/fact
  end do

  fact = fact/(lengthRef*LRef)
  
  do sps2 = 1,nTimeIntervalsSpectral
     moment(:, sps2) = (cMp + cMV)/fact
  end do

#ifndef USE_COMPLEX  
  call getCostFunction2(force, moment, sepSensor, sepSensorAvg, &
       Cavitation, alpha, beta, liftIndex)
#endif
end subroutine block_res

subroutine resScale
  
  use blockPointers
  use flowVarRefState     
  use inputIteration

  implicit none

  ! Local Variables
  integer(kind=intType) :: i, j, k, l
  real(kind=realType) :: ovol
  
 ! Divide through by the reference volume

  do k=2, kl
     do j=2, jl
        do i=2, il
           oVol = one/volref(i,j,k)
           do l=1, nwf
              dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l))* ovol
           end do
           do l=nt1, nt2
              dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l))* ovol * turbResScale(l-nt1+1)
           end do
        end do
     end do
  end do
end subroutine resScale

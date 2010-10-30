!
!      ******************************************************************
!      *                                                                *
!      * File:          solverUnsteadyExplicitRK.F90                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-13-2006                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine solverUnsteadyExplicitRK
!
!      ******************************************************************
!      *                                                                *
!      * solverUnsteadyExplicitRK solves the unsteady equations using   *
!      * the explicit Runge-Kutta schemes for the multigrid level       *
!      * groundLevel.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use flowVarRefState
       use inputPhysics
       use inputUnsteady
       use iteration
       use killSignals
       use monitor
       implicit none
!
!      Local parameter.
!
       integer(kind=intType), parameter :: nWriteConvHeader = 50
!
!      Local variables.
!
       integer(kind=intType) :: iter, nTimeSteps, stage
       integer(kind=intType) :: i, j, k, l, m, nn

       real(kind=realType) :: tmp, timeUnsteadyOld

       character(len=7) :: numberString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set rkStage to 0. This variable is only relevant if Runge-Kutta
       ! smoothers are used as an iterative algorithm, but rkStage needs
       ! to be initialized.

       rkStage = 0

       ! Initializations of the write and signal parameters.

       writeVolume  = .false.
       writeSurface = .false.
       writeGrid    = .false.
       globalSignal = noSignal
       localSignal  = noSignal

       ! Initialize timeStepUnsteady to 0 and set the number of
       ! time steps depending on the grid level.

       timeStepUnsteady = 0

       nTimeSteps = nTimeStepsCoarse
       if(groundLevel == 1) nTimeSteps = nTimeStepsFine

       ! Write a message to stdout about how many time steps will
       ! be taken on the current grid. Only processor 0 does this.

       if(myID == 0) then
         write(numberString, "(i6)") nTimeSteps
         numberString = adjustl(numberString)

         print "(a)", "#"
         print 100, groundLevel, trim(numberString)
         print "(a)", "#"
 100     format("# Grid",1x,i1,": Performing",1x,a,1x, &
                "explicit Runge Kutta time steps.")

         ! Also write the convergence header. Technically this is
         ! not really a convergence history for explicit RK, but
         ! the routine is called this way.

         call convergenceHeader
       endif

       ! Determine and write the initial convergence info. Again not
       ! really convergence for explicit RK.

       call convergenceInfo
!
!      ******************************************************************
!      *                                                                *
!      * Loop over the number of time steps to be computed.             *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the unsteady time.

       timeUnsteadyOld = timeUnsteady
       timeUnsteady    = timeUnsteadyOld + gammaRKUnsteady(1)*deltaT

       timeLoop: do iter=1,nTimeSteps

         ! Rewrite the convergence header after a certain number of
         ! iterations. Only processor 0 does this.

         if(myID == 0 .and. mod(iter,nWriteConvHeader) == 0) &
           call convergenceHeader

         ! Loop over the number of RK stages.

         stageLoop: do stage=1,nRKStagesUnsteady

           ! Compute the residual. Note that the turbulent residual must
           ! be computed first, because it uses the array of the flow
           ! field residuals as temporary buffers.

           if(equations == RANSEquations) then
             call initres(nt1, nt2)
             call turbResidual
           endif

           call initres(1_intType, nwf)
           call residual

           ! Loop over the number of domains.

           domainLoop: do nn=1,nDom

             ! Set the pointers to this domain. Note that there is only
             ! one time instance for a time accurate computation.

             call setPointers(nn,currentLevel,1)

             ! Step 1: Store in dw the conservative update for the
             !         flow field variables. This means a multiplication
             !         with deltaT/vol. Also store in w the conservative
             !         flow field variables.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   tmp = deltaT/vol(i,j,k)
                   do l=1,nw
                     dw(i,j,k,l) = tmp*dw(i,j,k,l)
                   enddo

                   w(i,j,k,ivx) = w(i,j,k,ivx)*w(i,j,k,irho)
                   w(i,j,k,ivy) = w(i,j,k,ivy)*w(i,j,k,irho)
                   w(i,j,k,ivz) = w(i,j,k,ivz)*w(i,j,k,irho)

                 enddo
               enddo
             enddo

             ! Step 2: Compute the new conservative flow field variables
             !         and primitive turbulence variables.

             do m=1,(stage-1)
               tmp = betaRKUnsteady(stage,m)
               if(tmp /= zero) then
                 do l=1,nw
                   do k=2,kl
                     do j=2,jl
                       do i=2,il
                         w(i,j,k,l) = w(i,j,k,l) - tmp*dwOldRK(m,i,j,k,l)
                       enddo
                     enddo
                   enddo
                 enddo
               endif
             enddo

             tmp = betaRKUnsteady(stage,stage)
             do l=1,nw
               do k=2,kl
                 do j=2,jl
                   do i=2,il
                     w(i,j,k,l) = w(i,j,k,l) - tmp*dw(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo

             ! Step 3. Convert the conservative variables back to
             !         primitive variables for this block. Compute the
             !         laminar and eddy viscosities as well.

             ! Convert the momentum into velocities.
             ! Also store the total energy in the pressure, such that
             ! the routine computePressure can be used.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   tmp          = one/w(i,j,k,irho)
                   w(i,j,k,ivx) = w(i,j,k,ivx)*tmp
                   w(i,j,k,ivy) = w(i,j,k,ivy)*tmp
                   w(i,j,k,ivz) = w(i,j,k,ivz)*tmp
                   p(i,j,k)     = w(i,j,k,irhoE)
                 enddo
               enddo
             enddo

             ! Compute the pressure.

             call computePressure(2_intType, il, 2_intType, jl, &
                                  2_intType, kl, 0_intType)

             ! Swap the pressure and total energy, because
             ! computePressure stores the pressure in the position of
             ! the total energy.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   tmp            = p(i,j,k)
                   p(i,j,k)       = w(i,j,k,irhoE)
                   w(i,j,k,irhoE) = tmp
                 enddo
               enddo
             enddo

             ! Compute the laminar and eddy viscosities.

             call computeLamViscosity
             call computeEddyViscosity

             ! Step 4. Store dw in dwOldRK if this is not the last
             !         RK stage.

             if(stage < nRKStagesUnsteady) then
               do l=1,nw
                 do k=2,kl
                   do j=2,jl
                     do i=2,il
                       dwOldRK(stage,i,j,k,l) = dw(i,j,k,l)
                     enddo
                   enddo
                 enddo
               enddo
             endif

           enddo domainLoop

           ! Exchange the pressure if the pressure must be exchanged
           ! early. Only the first halo's are needed, thus whalo1 is
           ! called.

           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
                       .false., .false.)

           ! Apply the boundary conditions to all blocks.

           call applyAllBC(.true.)
           if(equations == RANSEquations) call applyAllTurbBC(.true.)

           ! Exchange the halo data. As we are on the fine grid
           ! the second halo is needed.

           call whalo2(currentLevel, 1_intType, nw, .true., &
                       .true., .true.)

           ! Determine the time and initialize the geometrical and
           ! boundary info for next stage, if needed.

           if(stage < nRKStagesUnsteady) then
             timeUnsteady = timeUnsteadyOld &
                          + gammaRKUnsteady(stage+1)*deltaT

             call initStageRK(stage+1)
           endif

         enddo stageLoop

         ! Increment timeStepUnsteady and update
         ! timeUnsteady with the current time step.

         timeStepUnsteady = timeStepUnsteady + 1
         timeUnsteady     = timeUnsteadyOld  + deltaT

         ! Write the convergence info.

         call convergenceInfo

         ! Determine the time and initialize the geometrical and
         ! boundary info for the next time step, such that the writing
         ! of grid files with moving geometries is done correctly.

         timeUnsteadyOld = timeUnsteady
         timeUnsteady    = timeUnsteadyOld + gammaRKUnsteady(1)*deltaT

         call initStageRK(1_intType)

         ! Update PV3 solution if PV3 support is required.

#ifdef USE_PV3
         call pv_update(real(timeStepUnsteady,realPV3Type))
#endif

         ! Determine whether or not solution files must be written.

         call checkWriteUnsteadyInLoop

         ! Exit the loop if the corresponding kill signal
         ! has been received.

         if(globalSignal == signalWriteQuit) exit

       enddo timeLoop

       ! Determine whether or not the final solution must be written.

       call checkWriteUnsteadyEndLoop

       end subroutine solverUnsteadyExplicitRK

       !=================================================================

       subroutine initStageRK(stage)
!
!      ******************************************************************
!      *                                                                *
!      * initStageRK performs the initialization tasks for the          *
!      * Runge-Kutta schemes in unsteady mode.                          *
!      *                                                                *
!      ******************************************************************
!
       use inputMotion
       use inputUnsteady
       use iteration
       use monitor
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: stage
!
!      Local variables.
!
       integer(kind=intType) :: nn
       real(kind=realType)   :: fact

       real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If the grid is changing a whole lot of geometric
       ! info must be adapted.

       testChanging: if(changing_Grid .or. gridMotionSpecified) then

         ! Determine the time step relative to the previous stage.

         if(stage == 1) then
           fact = deltaT*(one - gammaRKUnsteady(nRKStagesUnsteady))
         else
           fact = deltaT*(gammaRKUnsteady(stage) &
                -         gammaRKUnsteady(stage-1))
         endif

         ! Set the time for all sections; also store their
         ! time step. They are the same for all sections, but all
         ! of them should be updated because of consistency.

         do nn=1,nSections
           tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
           deltaTSec(nn) = fact
         enddo

         ! Advance the coordinates 1 time step.

         call updateCoorFineMesh(deltaTSec, 1_intType)

         ! Adapt the geometric info on all grid levels needed for the
         ! current ground level and multigrid cycle.
         ! The wall distance only needs to be recomputed when the
         ! grid is changing; not when a rigid body motion is
         ! specified. Furthermore, the user can choose not to update
         ! the wall distance, because he may know a priori that the
         ! changes in geometry happen quite far away from the boundary
         ! layer. This is accomplished via updateWallDistanceUnsteady.
         ! Also note changingOverset can be true only when the grid
         ! is changing.

         call updateCoordinatesAllLevels
         if(changing_Grid .and. updateWallDistanceUnsteady) &
           call updateWallDistanceAllLevels

         if(changingOverset) call updateOversetAllLevels
         call updateSlidingAllLevels
         call updateMetricsAllLevels

         ! Update the rotation matrices of the faces. Only needed
         ! on the finest grid level.

         call faceRotationMatrices(currentLevel, .false.)

         ! Determine the velocities of the cell centers and faces
         ! for the current ground level. Note that the spectral mode
         ! is always 1 for unsteady mode.

         call gridVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

         ! Determine the new slip velocities on the viscous walls.

         call slipVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

       endif testChanging

       ! Determine the prescribed boundary condition data for the
       ! boundary subfaces for the (currently) finest grid.

       call setBCDataFineGrid(.false.)

       ! Non-dimensionalize the boundary data.

       call nonDimBoundData

       ! Set the primitive variables to the free stream values for
       ! the supersonic inflow faces for which this data has not
       ! been prescribed.

       call setSupersonicInletFreeStream

       ! Set the turbulent quantities to the free stream values for
       ! the inflow faces for which this data has not been prescribed.
       ! These can be both subsonic and supersonic inflow faces.

       call setInletFreestreamTurb

       end subroutine initStageRK

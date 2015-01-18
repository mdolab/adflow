!
!      ******************************************************************
!      *                                                                *
!      * File:          convergenceInfo.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-15-2003                                      *
!      * Last modified: 11-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine convergenceInfo
!
!      ******************************************************************
!      *                                                                *
!      * convergenceInfo computes and writes the convergence info to    *
!      * standard output. In spectral mode a convergence history for    *
!      * every mode is written.                                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use cgnsNames
       use inputIO
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use monitor
       use iteration
       use killSignals
       use NKsolverVars
       use flowVarRefState  ! eran-massf 
       use bleedFlows       ! eran-massf 
       use couplerParam     ! eran_idendifyname
       implicit none
!
!      Local variables.
!
       integer :: ierr, iConvStdout

       integer(kind=intType) :: sps, nn, mm, iConv
       integer(kind=intType) :: fTempMon=87 ! !eran-tempmon

       real(kind=realType) :: hdiffMax, MachMax
       real(kind=realType) :: eddyvisMax, yplusMax, sepSensor, Cavitation

       real(kind=realType) :: L2ConvThisLevel
       real(kind=realType) :: L2ConvThisLevelRel

       real(kind=realType), dimension(3) :: cfp, cfv, cmp, cmv

       integer(kind=intType) :: tempCurrentLevel,tempMGStartLevel

       logical :: nanOccurred, writeIterations
       logical :: relNotConv,absNotConv
!
!      Function definition
!
       logical :: myIsNAN
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!
!----------eran-cbd  for CBD run- just go the CBD printout -----
!
       chCBD : if(componentsBreakDown)then
          call componentsBreakDownPrintout(0)
          return
       end if chCBD
!
! -------eran-cbd ---------------------------------------------
!

       ! Determine whether or not the iterations must be written.

       writeIterations = .true.
       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) writeIterations = .false.

       ! Initialize converged to .true. This will be corrected to
       ! .false. later on if at least one of the spectral modes is not
       ! converged yet. Initialize nanOccurred to .false.

       converged   = .true.
       nanOccurred = .false.

       ! Set the L2 norm for convergence for this level.

       L2ConvThisLevel = L2ConvCoarse
       if(groundLevel == 1) L2ConvThisLevel = L2Conv

       L2ConvThisLevelRel = L2ConvRel


       ! Set the value of nIterCur, depending on groundLevel.

       nIterCur = iterTot
       if(groundLevel == 1) nIterCur = nIterCur + nIterOld

       ! Set the value of iConv, the place in the convergence array.
       ! On the finest level this is nIterCur. On the coarser grids this
       ! would be a logical choice as well. However it is theoretically
       ! possible (although highly unlikely) that more iterations are
       ! prescribed on the coarser grids than on the fine. As only the
       ! fine grid convergence history is stored, the dimensions of the
       ! convergence arrays are based on the number of fine grid
       ! iterations. To avoid any possible overflow this counter is set
       ! to either 0 or 1 on the coarse grids. The convergence history
       ! is overwritten anyway by the fine mesh.

       iConv = nIterCur
       if(groundLevel > 1) iConv = min(iConv,1_intType)

       ! Set the value of iConvStdout. For a steady and a spectral
       ! compution this is nIterCur; for an unsteady computation it is
       ! iterTot, the iteration number for the current time step.

       iConvStdout = nIterCur
       if(equationMode == unsteady) iConvStdout = iterTot

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Initialize the local monitoring variables to zero.

          monLoc = zero

         ! Loop over the blocks.

          domains: do nn=1,nDom

           ! Set the pointers for this block.

             call setPointers(nn, groundLevel, sps)

           ! Compute the forces and moments for this block.

             call forcesAndMoments(cfp, cfv, cmp, cmv, yplusMax, sepSensor, Cavitation)


           ! Determine the maximum values of the monitoring variables
           ! of this block.

             call maxHdiffMach(hdiffMax, MachMax)
             call maxEddyv(eddyvisMax)

           ! Loop over the number of monitoring variables.
           nMonitoringVar: do mm=1,nMon

             ! Determine the monitoring variable and act accordingly.

              select case (monNames(mm))
                
              case ('totalR')
                 call sumAllResiduals(mm)

                case (cgnsL2resRho)
                   call sumResiduals(irho, mm)

                case (cgnsL2resMomx)
                   call sumResiduals(imx, mm)

               case (cgnsL2resMomy)
                 call sumResiduals(imy, mm)

               case (cgnsL2resMomz)
                 call sumResiduals(imz, mm)

               case (cgnsL2resRhoe)
                 call sumResiduals(irhoE, mm)

               case (cgnsL2resNu, cgnsL2resK)
                 call sumResiduals(itu1, mm)

               case (cgnsL2resOmega, cgnsL2resTau, cgnsL2resEpsilon)
                 call sumResiduals(itu2, mm)

               case (cgnsL2resV2)
                 call sumResiduals(itu3, mm)

               case (cgnsL2resF)
                 call sumResiduals(itu4, mm)

               case (cgnsCl)
                 monLoc(mm) = monLoc(mm)                         &
                            + (cfp(1) + cfv(1))*liftDirection(1) &
                            + (cfp(2) + cfv(2))*liftDirection(2) &
                            + (cfp(3) + cfv(3))*liftDirection(3)

                case (cgnsClp)
                   monLoc(mm) = monLoc(mm) + cfp(1)*liftDirection(1) &
                        +              cfp(2)*liftDirection(2) &
                        +              cfp(3)*liftDirection(3)

                case (cgnsClv)
                   monLoc(mm) = monLoc(mm) + cfv(1)*liftDirection(1) &
                        +              cfv(2)*liftDirection(2) &
                        +              cfv(3)*liftDirection(3)

               case (cgnsCd)
                 monLoc(mm) = monLoc(mm)                         &
                            + (cfp(1) + cfv(1))*dragDirection(1) &
                            + (cfp(2) + cfv(2))*dragDirection(2) &
                            + (cfp(3) + cfv(3))*dragDirection(3)



                case (cgnsCdp)
                   monLoc(mm) = monLoc(mm) + cfp(1)*dragDirection(1) &
                        +              cfp(2)*dragDirection(2) &
                        +              cfp(3)*dragDirection(3)

                case (cgnsCdv)
                   monLoc(mm) = monLoc(mm) + cfv(1)*dragDirection(1) &
                        +              cfv(2)*dragDirection(2) &
                        +              cfv(3)*dragDirection(3)

                case (cgnsCfx)
                   monLoc(mm) = monLoc(mm) + cfp(1) + cfv(1)

                case (cgnsCfy)
                   monLoc(mm) = monLoc(mm) + cfp(2) + cfv(2)

                case (cgnsCfz)
                   monLoc(mm) = monLoc(mm) + cfp(3) + cfv(3)
                   
                case (cgnsCmx)
                   monLoc(mm) = monLoc(mm) + cmp(1) + cmv(1)

                case (cgnsCmy)
                   monLoc(mm) = monLoc(mm) + cmp(2) + cmv(2)

                case (cgnsCmz)
                   monLoc(mm) = monLoc(mm) + cmp(3) + cmv(3)

                case (cgnsHdiffMax)
                   monLoc(mm) = max(monLoc(mm), hdiffMax)

                case (cgnsMachMax)
                   monLoc(mm) = max(monLoc(mm), MachMax)

                case (cgnsYplusMax)
                   monLoc(mm) = max(monLoc(mm), yplusMax)

                case (cgnsEddyMax)
                   monLoc(mm) = max(monLoc(mm), eddyvisMax)
                   
                case (cgnsSepSensor)
                   monLoc(mm) = monLoc(mm) + sepSensor
                case (cgnsCavitation)
                   monLoc(mm) = monLoc(mm) + Cavitation

                end select ! monNames(mm)

             end do nMonitoringVar
          end do domains

         ! Determine the global sum of the summation monitoring
         ! variables. The sum only needs to be known on processor 0.

         if(nMonSum > 0) &
           call mpi_reduce(monLoc, monGlob, nMonSum, sumb_real, &
                           mpi_sum, 0, SUmb_comm_world, ierr)

         ! Idem for the maximum monitoring variables.
#ifndef USE_COMPLEX
          if(nMonMax > 0) &
               call mpi_reduce(monLoc(nMonSum+1), monGlob(nMonSum+1), &
               nMonMax, sumb_real, mpi_max, 0,        &
               SUmb_comm_world, ierr)
#else
          if (nMonMax < 0) & 
               monGlob(nMonSum+1) = zero
#endif
!
! ---- eran-massf ---
!
          GroundLevelCh : if(groundLevel == 1 )then

             if(nOutflowSubsonic + nOutflowBleeds + nInflowSubsonic > 0 )then
!
! calculate mass flux through bcOutflowSubsonic boundaries
!
                call calMassFlux(massFluxG)
 
             end if ! n BCs

          end if GroundLevelCh
         
! ----eran-massf end --
!
         ! Write the convergence info; only processor 0 does this.

          testRootProc: if(myID == 0) then

             if(.not.standAloneMode .and. printIterations )then
!
! ---------- eran_idendifyname ------  For CHIMPS run print identifier infront each line (code-name)
!
                write(*,'(a15)',advance="no")trim(codeName)
             end if
!----------end eran_idendifyname ---

           ! The variables which must always be written.

            if(printIterations) then
             write(*,"(1x,i6,2x)",advance="no") groundLevel
             
             if(equationMode == unsteady) then

                  write(*,"(i6,1x)",advance="no") timeStepUnsteady + &
                       nTimeStepsRestart
                  write(*,"(e12.5,1x)",advance="no") timeUnsteady + &
                       timeUnsteadyRestart
                  
               else if(equationMode == timeSpectral) then
                  
                  write(*,"(i8,3x)",advance="no") sps
                  
               endif


               if( writeIterations ) &
                    write(*,"(i6,1x)",advance="no") iConvStdout
               if( showCPU ) &
                    write(*,"(e12.5,1x)",advance="no") mpi_wtime() - t0Solver
            end if
           ! Loop over the number of monitoring values.

             do mm=1,nMon

             ! The residual variables must be corrected.

                select case (monNames(mm))

               case (cgnsL2resRho,  cgnsL2resMomx,    &
                     cgnsL2resMomy, cgnsL2resMomz,    &
                     cgnsL2resRhoe, cgnsL2resNu,      &
                     cgnsL2resK,    cgnsL2resOmega,   &
                     cgnsL2resTau,  cgnsL2resEpsilon, &
                     cgnsL2resV2,   cgnsL2resF        )
                  monGlob(mm) = sqrt(monGlob(mm)/nCellGlobal(groundLevel))
                  if( myIsNAN(monGlob(mm)) ) nanOccurred = .true.                  
               case ('totalR')
                  monGlob(mm) = sqrt(monGlob(mm))
                  if( myIsNAN(monGlob(mm)) ) nanOccurred = .true.                  
               end select

             ! Write the convergence info to stdout.
             if ( printIterations) then
#ifndef USE_COMPLEX
                write(*,"(e24.16,1x)",advance="no") monGlob(mm)
#else
                write(*,"(2e24.16,1x)",advance="no") monGlob(mm)
#endif
             end if
          enddo

!
             !--- eran-massf ---

             if(nOutflowSubsonic + nOutflowBleeds + nInflowSubsonic > 0 )& 
                  write(*,"(e12.5,1x)",advance="no")massFluxG
! 
             ! ----eran-massf end
!

           ! Write the carriage return.
          if (printIterations) then
             print "(1x)"
          end if
           ! Store the convergence info in convArray, if desired.

           if( storeConvInnerIter ) then
             do mm=1,nMon
               convArray(iConv,sps,mm) = monGlob(mm)
             enddo
           endif

           ! Determine whether or not the solution is converged.
           ! A distinction must be made between unsteady mode and the
           ! other modes, because in unsteady mode the convergence of
           ! the inner iterations is not necessarily stored.

             select case (equationMode)
             case (steady, timeSpectral)
                if (.not. coeffConvCheck) then

                   ! Steady or time spectral mode. The convergence histories
                   ! are stored and this info can be used. The logical
                   ! converged is set to .false. if the density residual
                   ! has not converged yet.
                   if (iterTot > minIterNum) then 
                      if(convArray(iConv,sps,1) > L2ConvThisLevel*convArray(0,sps,1)) then
                         absNotConv = .True.
                      else
                         absNotConv = .False.
                      end if
                      
                      if(fromPython) then
                         if (convArray(iConv,sps,1) > L2ConvThisLevelRel*convArray(1,sps,1)) then
                            relNotConv = .True.
                         else
                            relNotConv = .False.
                         end if
                      else
                         relNotConv = .True.
                      end if
                      
                      if (absNotConv .and. relNotConv) then ! Not converged if the absCheck is True and the rel Check is true.
                         converged = .False.
                      end if
                   else
                      converged = .False.
                   end if
                else
!----eran-coeffConv starts
                   
                   convergenceQuality = 0 ! that is no convergence
                   converged          = .false. 
                   
                   if (iterTot >= minIterNum) then
                      
                      if (iterTot == minIterNum)then
                         write(*,*)'#***************************************************************'
                         write(*,*)'# Note: at step ',iterTot,&
                              ' Starting to test for convergence'
                         write(*,*)'#***************************************************************'
                      end if
                      
                      if(convArray(iConv,sps,1) <= L2ConvThisLevel*convArray(0,sps,1)) then
                         converged = .true.
                         convergenceQuality = 10
                      end if
                      
                      if (epsCoefConv > zero   .and.&
                           (groundLevel == 1 .and. converged .eqv. .false.)  )then
                         ! !
                         ! ! ---- Check if coefficients reached a cconstant value
                         ! !
                         call coeffConvergenceCheck(iConv,iterTot,sps)
                         if(convergenceQuality > 0) then
                            converged = .true.
                            write(*,*)&
                                 'convegenceInfo: Coefficients convergence criterion reached'
                         end if
                      end if ! epsCoefConv > zero
                      
                      if(converged)then
                         select case (convergenceQuality)
                         case(10)
                            write(*,*)'Convergence: Residual < Convergence criterion'
                         case(6)
                            write(*,*)'Coefficient uniform (up to criterion) in ',ConvCheckWindowSize,&
                                 ' iterations'
                         case(4)
                            write(*,*)'Coefficient uniform (up to criterion in ',10*ConvCheckWindowSize,&
                                 ' iterations'
                         case(2)
                            write(*,*)'Coefficient uniform (up to criterion in ',100*ConvCheckWindowSize,&
                                 ' iterations'
                         end select
                      end if ! converged
                   end if ! iterTot >= minIterNum
                   ! ! ------- end eran-coeffConv
                end if

             !===========================================================

             case (unsteady)
             
               ! Unsteady mode. The array convArray may not be present
               ! and therefore something else must be done.
               ! First determine the position in the array timeDataArray
               ! that can be used. For the coarser grids this is 1,
               ! because the time evolution is overwritten. For the fine
               ! mesh the actual position is determined.
             
                nn = 1
                if(groundLevel == 1) &
                     nn = timeStepUnsteady + nTimeStepsRestart
                nn = max(nn,1_intType)

               ! Make a distinction between the time integration
               ! schemes.


                select case(timeIntegrationScheme)

                case (explicitRK)

                   ! Explicit scheme. Simply store the data in the
                   ! convergence arrays.

                   timeArray(nn) = timeUnsteady + timeUnsteadyRestart

                   ! For explicit schemes the residuals are not
                   ! monitored and therefore the monitoring variables
                   ! can simply be copied.

                   do mm=1,nMon
                      timeDataArray(nn,mm) = monGlob(mm)
                   enddo

                 !=======================================================

               case (BDF,implicitRK)

                   ! An implicit scheme is used and therefore an
                   ! iterative algorithm within every time step.
                   ! The array convArray may not be present and
                   ! therefore
                   ! something else must be done. First determine the
                   ! position in the array timeDataArray that can be
                   ! used.
                   ! For the coarser grids this is 1, because the time
                   ! evolution is overwritten. For the fine mesh the
                   ! actual position is determined.

                   ! Determine the situation we have here.

                   testInitUnsteady: if(iterTot == 0) then

                     ! This is the initialization phase for this time step.
                     ! Simply copy monGlob into monRef, store the value
                     ! of the physical time and set converged to .false..


                      do mm=1,nMon
                         monRef(mm) = monGlob(mm)
                      enddo

                      timeArray(nn) = timeUnsteady + timeUnsteadyRestart
                     converged     = .false.

                   else testInitUnsteady

                     ! Iteration for this time step. Store the relative
                     ! convergence for the residual compared to the start
                     ! of this time step; an absolute norm does not give
                     ! any information here. For all other monitoring 
                     ! variables the current value is stored.

                      do mm=1,nMon
                         select case (monNames(mm))
                            
                         case (cgnsL2resRho,  cgnsL2resMomx,    &
                               cgnsL2resMomy, cgnsL2resMomz,    &
                               cgnsL2resRhoE, cgnsL2resNu,      &
                               cgnsL2resK,    cgnsL2resOmega,   &
                               cgnsL2resTau,  cgnsL2resEpsilon, &
                               cgnsL2resV2,   cgnsL2resF)

                          timeDataArray(nn,mm) = monGlob(mm) &
                                               / max(monRef(mm),eps)

                         !===============================================

                         case default


                            timeDataArray(nn,mm) = monGlob(mm)

                       end select
                     enddo

                     ! Set the logical converged to .false. if the density
                     ! residual has not converged yet.

                     if(timeDataArray(nn,1) > L2ConvThisLevel) &
                       converged = .false.

                   endif testInitUnsteady
                end select ! temporal integration scheme
!
! add temporal monitoring option for unsteady cases  ! eran-tempmon starts
!
                if(nIterCur == nCycles .or. converged .eqv. .true.)then

                   write(fTempMon,"(i6,1x)",advance="no") timeStepUnsteady + &
                        nTimeStepsRestart
                   write(fTempMon,"(e12.5,1x)",advance="no") timeUnsteady + &
                        timeUnsteadyRestart
                   do mm=1,nMon
                      select case (monNames(mm))
                      case(cgnsCl,cgnsClp,cgnsClv, cgnsCd, cgnsCdp, cgnsCdv ,&
                           cgnsCfx, cgnsCfy, cgnsCfz, cgnsCmx, cgnsCmy, cgnsCmz)

                         write(fTempMon,"(e12.5,1x)",advance="no") monGlob(mm)

                      end select ! monitored vars
                   end do ! mm
                   
                   if(nOutflowSubsonic + nOutflowBleeds > 0 )& 
                        write(fTempMon,"(e12.5,1x)",advance="no")massFluxG

                   
                   write(fTempMon,"(1x)") !  eran-tempmon
                   
                end if ! temp-mon

! -----------------eran-tempmon ends

             end select ! unsteady

          endif testRootProc
       enddo spectralLoop

       ! Check if a NaN occured. If so the computation is terminated,
       ! such that possible solution files written earlier are not
       ! corrupted.
       
       call mpi_bcast(nanOccurred, 1, MPI_LOGICAL, 0, SUmb_comm_world, ierr)
       if( nanOccurred )then
          !reset flow and exit!
          print *,'Nan occured in Convergence Info on proc:',myid
          !call initflow
          ! Initialize the flow field to uniform flow.
          tempMGStartLevel = mgStartLevel
          tempCurrentLevel = currentLevel
          
          mgStartLevel = 1
          currentlevel = 1

          call setUniformFlow

          mgStartLevel = tempMGStartLevel
          currentLevel = tempCurrentLevel
  
          call terminate("convergenceInfo", &
               "A NaN occurred during the computation.")

          ! in a normal computation, code will simply exit.
          ! in a python based computation, code will set 
          ! routinedFailed to .True. and return to the 
          ! python level...
          return
       endif
       if((fromPython).and. (nIterCur==nCycles))then
          
          !Check to see if residuals are diverging or stalled for python
          select case (equationMode)
             
          case (steady, timeSpectral)
             
             ! Steady or time spectral mode. The convergence histories
             ! are stored and this info can be used. If the residuals 
             ! are diverging the, logical routineFailed in killSignals
             ! is set to true and the progress is halted.
             !only check on root porcessor
             if (myID==0)then

                ! If we made it to ncycles, check to see if we're
                ! "close" to being converged. 
                routineFailed = .False.
                do sps = 1,nTimeIntervalsSpectral
                   if(convArray(iConv,sps,1) > &
                        maxL2DeviationFactor * L2ConvThisLevel) then 
                      routineFailed = .True.
                      exit
                   end if
                enddo

             endif
             call mpi_bcast(routineFailed, 1, MPI_LOGICAL, 0, SUmb_comm_world, ierr)

          case(unsteady)
             !print *,'divergence check for unsteady not implemented...'
             !return
          end select
       else
          routineFailed = .False.
       end if

       ! Determine whether or not the solution is considered
       ! converged.  This info is only known at processor 0 and must
       ! therefore be broadcast to the other processors. MPI supports
       ! the communication of logicals, but this is done with
       ! integers.  Therefore send an integer and retrieve the
       ! information from it.  Remember that the logical converged was
       ! initialized to .true.

       mm = 0
       if( converged ) mm = 1

       call mpi_bcast(mm, 1, sumb_integer, 0, SUmb_comm_world, ierr)
       if(mm == 0) converged = .false.


       end subroutine convergenceInfo

!      ==================================================================

       subroutine sumResiduals(nn, mm)
!
!      ******************************************************************
!      *                                                                *
!      * sumResiduals adds the sum of the residuals squared at          *
!      * position nn to the array monLoc at position mm. It is assumed  *
!      * that the arrays of blockPointers already point to the correct  *
!      * block.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use monitor
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, mm
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of owned cells of this block and
       ! accumulate the residual.

       do k=2,kl
         do j=2,jl
           do i=2,il
             monLoc(mm) = monLoc(mm) + (dw(i,j,k,nn)/vol(i,j,k))**2
           enddo
         enddo
       enddo

       end subroutine sumResiduals

      subroutine sumAllResiduals(mm)
!
!      ******************************************************************
!      *                                                                *
!      * sumAllResiduals adds the sum of the ALL residuals squared at   *
!      * to monLoc at position mm.                                      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use monitor
       use flowvarrefstate
       use inputIteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: mm
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l
       real(kind=realType) :: state_sum,ovv

!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of owned cells of this block and
       ! accumulate the residual.

       do k=2,kl
         do j=2,jl
           do i=2,il
              state_sum = 0.0
              ovv = one/vol(i,j,k)
              do l=1,nwf
                 state_sum = state_sum + (dw(i,j,k,l)*ovv)**2
              end do
              do l=nt1,nt2
                 state_sum = state_sum + (dw(i,j,k,l)*ovv*turbResScale)**2
              end do
              monLoc(mm) = monLoc(mm) + state_sum

           enddo
         enddo
       enddo

     end subroutine sumAllResiduals

     subroutine printCurrentR()
       use communication 
       use block
       use blockPointers
       use flowvarrefstate
       use inputtimespectral
       implicit none

       integer(kind=intType) :: nn,sps,i,j,k,l,ierr
       real(kind=realType) :: state_sum,rnorm,ovv

       state_sum = 0.0
       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom
            call setPointers(nn,1,sps)
            do k=2,kl
               do j=2,jl
                  do i=2,il
                     ovv = 1/vol(i,j,k)
                     do l=1,nw
                        state_sum = state_sum + (dw(i,j,k,l)*ovv)**2
                     end do
                  end do
               end do
            end do
         end do domains
      end do spectralLoop
      call mpi_reduce(state_sum,rnorm,1,sumb_real,mpi_sum,0,&
           SUmb_comm_world, ierr)

      if (myid==0) then
         print *,'Current R:',sqrt(rnorm)
      end if

     end subroutine printCurrentR

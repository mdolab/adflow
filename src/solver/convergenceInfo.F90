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
  use nksolvervars
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr, iConvStdout

  integer(kind=intType) :: sps, nn, mm, iConv

  real(kind=realType) :: hdiffMax, MachMax
  real(kind=realType) :: eddyvisMax, yplusMax, sepSensor, Cavitation
  real(kind=realType) :: sepSensorAvg(3)

  real(kind=realType) :: L2ConvThisLevel
  real(kind=realType), dimension(3) :: cfp, cfv, cmp, cmv
  logical :: nanOccurred, writeIterations
  logical :: absConv, relConv
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
  ! Determine whether or not the iterations must be written.

  writeIterations = .true.
  if(equationMode          == unsteady .and. &
       timeIntegrationScheme == explicitRK) writeIterations = .false.

  ! Initializations
  converged   = .False.
  nanOccurred = .false.

  ! Set the L2 norm for convergence for this level.

  L2ConvThisLevel = L2ConvCoarse
  if(groundLevel == 1) L2ConvThisLevel = L2Conv

  ! Loop over the number of spectral solutions.

  spectralLoop: do sps=1,nTimeIntervalsSpectral

     ! Initialize the local monitoring variables to zero.

     monLoc = zero

     ! Loop over the blocks.

     domains: do nn=1,nDom

        ! Set the pointers for this block.

        call setPointers(nn, groundLevel, sps)

        ! Compute the forces and moments for this block.

        call forcesAndMoments(cfp, cfv, cmp, cmv, yplusMax, sepSensor, &
             sepSensorAvg, Cavitation)

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

     ! Add the corrections from zipper meshes from proc 0
     ! ----------------------------------------------------------
     call forcesAndMomentsZipper(cfp, cfv, cmp, cmv)

     ! Save lift and drag contributions from zipper meshes
     if (myid == 0) then
        write(4000,"(1x,i6,2x)",advance="no") groundLevel
        write(4000,"(i6,1x)",advance="no") iterTot
        write(4000,"(i6,1x)",advance="no") approxTotalIts

        ! Lift
        write(4000,"(e24.16,1x)",advance="no") &
          (cfp(1) + cfv(1))*liftDirection(1) &
        + (cfp(2) + cfv(2))*liftDirection(2) &
        + (cfp(3) + cfv(3))*liftDirection(3)

        ! Drag
        write(4000,"(e24.16,1x)",advance="no") &
          (cfp(1) + cfv(1))*dragDirection(1) &
        + (cfp(2) + cfv(2))*dragDirection(2) &
        + (cfp(3) + cfv(3))*dragDirection(3)
        write(4000,"(1x)")
     end if
     
     ! Loop over the number of monitoring variables.
     nMonitoringVarZip: do mm=1,nMon

        ! Determine the monitoring variable and act accordingly.

        select case (monNames(mm))

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

        end select ! monNames(mm)

     end do nMonitoringVarZip
     ! --- End corrections from zipper meshes from proc 0 ---

     ! Determine the global sum of the summation monitoring
     ! variables. This is an all reduce since every processor needs to
     ! know the residual to make the same descisions. 

     if(nMonSum > 0) &
          call mpi_allreduce(monLoc, monGlob, nMonSum, sumb_real, &
          mpi_sum, SUmb_comm_world, ierr)

     ! Idem for the maximum monitoring variables.
#ifndef USE_COMPLEX
     if(nMonMax > 0) &
          call mpi_allreduce(monLoc(nMonSum+1), monGlob(nMonSum+1), &
          nMonMax, sumb_real, mpi_max, SUmb_comm_world, ierr)
#else
     if (nMonMax < 0) & 
          monGlob(nMonSum+1) = zero
#endif

     ! Write the convergence info; only processor 0 does this.

     testRootProc: if(myID == 0) then

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

           if( writeIterations ) then 
              write(*,"(i6,1x)",advance="no") iterTot
              write(*,"(i6,1x)",advance="no") approxTotalIts
              write(*,"(a,1x)", advance="no") iterType
              write(*,"(e10.2,1x)",advance="no") CFLMonitor
           end if

           if( showCPU ) &
                write(*,"(e12.5,1x)",advance="no") mpi_wtime() - t0Solver
        end if
     end if testRootProc

     ! Loop over the number of monitoring values.
     variableLoop: do mm=1, nMon

        ! The residual variables must be corrected.

        select case (monNames(mm))

        case (cgnsL2resRho,  cgnsL2resMomx,    &
             cgnsL2resMomy, cgnsL2resMomz,    &
             cgnsL2resRhoe, cgnsL2resNu,      &
             cgnsL2resK,    cgnsL2resOmega,   &
             cgnsL2resTau,  cgnsL2resEpsilon, &
             cgnsL2resV2,   cgnsL2resF        )
           monGlob(mm) = sqrt(monGlob(mm)/nCellGlobal(groundLevel))
           if (monNames(mm) == cgnsL2resRho) then 
              rhoRes = monGlob(mm)
           end if
        case ('totalR')
           monGlob(mm) = sqrt(monGlob(mm))
           if( myIsNAN(monGlob(mm)) ) nanOccurred = .true.
           totalR = monGlob(mm)
        end select

        if( myIsNAN(monGlob(mm)) ) nanOccurred = .true.                  

        if (myid == 0 .and. printIterations) then 
           ! Write the convergence info to stdout.
#ifndef USE_COMPLEX
           write(*,"(e24.16,1x)",advance="no") monGlob(mm)
#else
           write(*,"(2e24.16,1x)",advance="no") monGlob(mm)
#endif
        end if

        ! Store the convergence info in convArray, if desired.
        if( storeConvInnerIter ) then
           convArray(iterTot, sps, mm) = monGlob(mm)
        endif
     end do variableLoop

     if (myid == 0 .and. printIterations) then 
        ! Write the carriage return.
        print "(1x)"
     end if

     ! Determine whether or not the solution is converged.
     ! A distinction must be made between unsteady mode and the
     ! other modes, because in unsteady mode the convergence of
     ! the inner iterations is not necessarily stored.

     select case (equationMode)
     case (steady, timeSpectral)

        ! Steady or time spectral mode. The convergence histories
        ! are stored and this info can be used. The logical
        ! converged is set to .false. if the density residual
        ! has not converged yet.

        absConv = .False.
        relConv = .False. 

        ! We make a split here based on if we are operating on a
        ! coarse grid or the finest. On the coarse grid levels, we
        ! l2convThisLevel refers to the relative convergence of
        ! rhoRes. 

        if (currentLevel /= 1) then 
           if (rhoRes < L2ConvThisLevel * rhoResStart) then 
              relConv = .True. 
           end if
        else
           ! We are on the fine level. All checking is done using
           ! the total residual. 

           ! Absolute convergence check
           if (totalR < L2ConvThisLevel * totalR0) then 
              absConv = .True. 
           end if

           ! Relative check only done on finest level
           if (totalR < L2ConvRel*totalRStart) then 
              relConv = .True. 
           end if

           ! If the totla number of iterations is less than the
           ! RKReset, don't check the residual
           if (iterTot < minIterNum .and. rkreset) then 
              relConv = .False.
              absConv = .False.
           end if
        end if

        ! Combine the two flags.
        if (absConv .or. relConv)  then 
           converged = .True.
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

        case (MD,BDF,implicitRK)

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

              if(timeDataArray(nn,1) < L2ConvThisLevel) &
                   converged = .True.

           endif testInitUnsteady
        end select ! temporal integration scheme

     end select ! unsteady

  end do spectralLoop


  if( nanOccurred )then

     if (myid == 0) then 
        print *,'Nan occured in Convergence Info on proc:',myid
     end if

     routineFailed = .True.

     call returnFail("convergenceInfo", &
          "A NaN occurred during the computation.")

     ! in a normal computation, code will simply exit.
     ! in a python based computation, code will set 
     ! routinedFailed to .True. and return to the 
     ! python level...
     return
  endif

  ! ! If we are at the max iteration limit but the residual is
  ! ! *close*, ie within maxL2DeviationFactor we say that's fine

  ! if(fromPython .and. groundLevel ==1 .and. nIterCur >=nCycles) then

  !    !Check to see if residuals are diverging or stalled for python
  !    select case (equationMode)

  !    case (steady, timeSpectral)

  !       ! Steady or time spectral mode. The convergence histories
  !       ! are stored and this info can be used. If the residuals 
  !       ! are diverging the, logical routineFailed in killSignals
  !       ! is set to true and the progress is halted.
  !       !only check on root porcessor
  !       if (myID == 0) then

  !          ! If we made it to ncycles, check to see if we're
  !          ! "close" to being converged. 
  !          do sps = 1,nTimeIntervalsSpectral         
  !             if (totalR < maxL2DeviationFactor * totalR0 * L2ConvThisLevel) then 
  !                routineFailed = .False.
  !             else
  !                routineFailed = .True.
  !             end if
  !          enddo
  !       endif

  !       call mpi_bcast(routineFailed, 1, MPI_LOGICAL, 0, SUmb_comm_world, ierr)

  !    end select
  ! end if

  ! Determine whether or not the solution is considered
  ! converged.  This info is only known at processor 0 and must
  ! therefore be broadcast to the other processors.

  call mpi_bcast(converged, 1, MPI_LOGICAL, 0, SUmb_comm_world, ierr)

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
              ! l-nt1+1 will index the turbResScale properly
              state_sum = state_sum + (dw(i,j,k,l)*ovv*turbResScale(l-nt1+1))**2
           end do
           monLoc(mm) = monLoc(mm) + state_sum
        enddo
     enddo
  enddo

end subroutine sumAllResiduals

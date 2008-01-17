!
!      ******************************************************************
!      *                                                                *
!      * File:          sumb_coupler_m.F90                              *
!      * Author:        Seonghyeon Hahn                                 *
!      * Starting date: 02-16-2006                                      *
!      * Last modified: 05-11-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       module sumb_coupler_m
!
!      ******************************************************************
!      *                                                                *
!      * This module contains all the API subroutines for coupled       *
!      * simulations using CHIMPS.                                      *
!      *                                                                *
!      ******************************************************************
!
       implicit none

       interface sumb_getParam
         module procedure sumb_getParamInt
         module procedure sumb_getParamIntRk1
         module procedure sumb_getParamReal
         module procedure sumb_getParamRealRk1
         module procedure sumb_getParamStr
         module procedure sumb_getParamLog
       end interface

       interface sumb_setParam
         module procedure sumb_setParamInt
         module procedure sumb_setParamIntRk1
         module procedure sumb_setParamReal
         module procedure sumb_setParamRealRk1
         module procedure sumb_setParamStr
         module procedure sumb_setParamLog
       end interface

       contains

       subroutine sumb_initialize(paramName, comm)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_initialize initializes this particular instance of SUmb   *
!      * in case of coupled simulation.                                 *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use couplerParam
       use inputIO
       use iteration
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments
!
       character(len=*), intent(in) :: paramName
       integer, intent(in) :: comm
!
!      Local variables
!
       integer :: ierr
       integer :: lenParam
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       standAloneMode  = .false.
       deforming_Grid  = .false.
       SUmb_comm_world = comm

       ! Determine the rank and number of processors inside the group
       ! defined by SUmb_comm_world.

       call mpi_comm_rank(SUmb_comm_world, myID,  ierr)
       call mpi_comm_size(SUmb_comm_world, nProc, ierr)

       ! Set the name of the parameter file.

       lenParam = len(paramName)
       paramFile = paramName(:lenParam)
       paramFile = adjustl(paramFile)
       paramFile = trim(paramFile)

       ! Write a message to stdout with information how the
       ! executable was built.

       call writeIntroMessageCpl

       ! Allocate the memory for sendRequests and recvRequests.

       allocate(sendRequests(nProc), &
                recvRequests(nProc), stat=ierr)
       if(ierr /= 0)                &
         call terminate("sumb_initialize", &
                        "Memory allocation failure for sendRequests &
                        &and recvRequests")

       ! Read the parameter file.

       call readParamFile

       ! Partition the blocks and read the grid.

       call partitionAndReadGrid

       ! Perform the preprocessing task.

       call preprocessing

       ! Initialize of the flow variables.

       call initFlow

       end subroutine sumb_initialize

!========================================================================

       subroutine writeIntroMessageCpl
!
!      ******************************************************************
!      *                                                                *
!      * writeIntroMessageCpl writes a message to stdout with           *
!      * information how the executable was built, e.g. whether single  *
!      * or double precision is used for the integers and reals, etc.   *
!      * To avoid a messy output only processor 0 prints this info.     *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       implicit none
!
!      Local variables
!
       character(len=7) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return if this is not processor 0.

       if(myID > 0) return

       ! I'm processor 0. Write the info to stdout.

       print "(a)", "#"
       print "(a)", "# SUmb, multiblock structured flow solver"
       print "(a)", "#"
       print "(a)", "# This code solves the 3D RANS, laminar NS or &
                    &Euler equations"
       print "(a)", "# on multiblock structured hexahedral grids."

       if( SU_MPI_isSequential ) then
         print "(a)", "# This is a sequential executable compiled &
                      &with the following options:"
       else
         write(integerString,"(i7)") nProc
         integerString = adjustl(integerString)
         print "(3a)", "# This is a parallel executable running on ", &
                       trim(integerString), " processors."
         print "(a)", "# It has been compiled with the &
                      &following options:"
       endif

       if( debug ) then
         print "(a)", "# - Debug mode."
       else
         print "(a)", "# - Optimized mode."
       endif

#ifdef USE_LONG_INT
       print "(a)", "# - Size of standard integers: 8 bytes."
#else
       print "(a)", "# - Size of standard integers: 4 bytes."
#endif

#ifdef USE_SINGLE_PRECISION
       print "(a)", "# - Size of standard floating point types: &
                    &4 bytes."

#elif  USE_QUADRUPLE_PRECISION
       print "(a)", "# - Size of standard floating point types: &
                    &16 bytes."
#else
       print "(a)", "# - Size of standard floating point types: &
                    &8 bytes."
#endif

#ifdef USE_PV3
       print "(a)", "# - With pV3 support"
#else
       print "(a)", "# - Without pV3 support"
#endif

#ifdef USE_NO_CGNS
       print "(a)", "# - Without cgns support"
#else
       print "(a)", "# - With cgns support"
#endif

       if(.not. SU_MPI_isSequential) then
         if( SU_MPI_noMPIO ) then
           print "(a)", "# - Without parallel IO support"
         else
           print "(a)", "# - With parallel IO support"
         endif
       endif

#ifdef USE_NO_SIGNALS
       print "(a)", "# - Without support for signals."
#else
       print "(a)", "# - With support for signals."
#endif

       print "(a)", "#"

       end subroutine writeIntroMessageCpl

!      ==================================================================

       subroutine sumb_runIteration(nIterCpl)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_runIteration runs the actual multigrid cycles or time     *
!      * marching in case of coupled simulation. Here, nIterCpl is the  *
!      * number of unsteady physical time steps or of multigrid cycles  *
!      * on the finest mesh in unsteady and steady simulations,         *
!      * respectively.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use precision
       use communication
       use constants
       use inputDiscretization
       use inputIteration
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use killSignals
       use iteration
       use monitor
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nIterCpl
!
!      Local variables.
!
       real(kind=realType), dimension(nSections) :: dtAdvance
       integer(kind=intType) :: i, j, k, kk, ll, mm, nn
       integer :: ierr
       real(kind=cgnsRealType), dimension(:),     &
                                allocatable :: tmpTimeArray
       real(kind=cgnsRealType), dimension(:,:),   &
                                allocatable :: tmpTimeDataArray
       real(kind=cgnsRealType), dimension(:,:,:), &
                                allocatable :: tmpConvArray
       integer(kind=intType) :: iter, nTimeSteps
!
!      Function definitions.
!
       logical          :: EulerWallsPresent
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if((.not. restart) .and. &
          nIterOld == 0 .and. nIterCur == 0 .and. iterTot == 0) then

         ! No iteration has been done.

         ! If unsteady, reset the number of unsteady time steps.
         ! If steady, reset the number of multigrid cycles on the
         ! finest mesh.

         if(equationMode == unsteady .and. nIterCpl /= nTimeStepsFine) &
           nTimeStepsFine = nIterCpl
         if(equationMode == steady .and. nIterCpl /= nCycles) &
           nCycles = nIterCpl

         ! Reallocate the arrays for convergence history and time data
         ! with a new size.

         if(myID == 0) then
           if(equationMode == unsteady) then
             if(allocated(timeArray)) then
               deallocate(timeArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeArray")
             endif

             if(allocated(timeDataArray)) then
               deallocate(timeDataArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeDataArray")
             endif

             call allocTimeArrays(nTimeStepsFine)
           endif

           if(storeConvInnerIter) then
             if(allocated(convArray)) then
               deallocate(convArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for convArray")
             endif

             kk = nTimeStepsFine
             if(equationMode == steady) kk = 1
             call allocConvArrays(kk*nCycles)

           endif
         endif

       else if(nIterCur == 0 .and. iterTot == 0) then

         ! Started from a restart file but no new iteration.

         ! If unsteady, reset the number of unsteady time steps.
         ! If steady, reset the number of multigrid cycles on the
         ! finest mesh.

         if(equationMode == unsteady .and. nIterCpl /= nTimeStepsFine) &
           nTimeStepsFine = nIterCpl
         if(equationMode == steady .and. nIterCpl /= nCycles) &
           nCycles = nIterCpl

         ! Reallocate the arrays for convergence history and time data
         ! with a new size, after storing values from the restart file.

         if(myID == 0) then
           if(equationMode == unsteady) then
             if(allocated(tmpTimeArray)) then
               deallocate(tmpTimeArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpTimeArray")
             endif

             allocate(tmpTimeArray(nTimeStepsRestart), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpTimeArray")

             do i=1,nTimeStepsRestart
               tmpTimeArray(i) = timeArray(i)
             enddo

             if(allocated(timeArray)) then
               deallocate(timeArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeArray")
             endif

             if(allocated(tmpTimeDataArray)) then
               deallocate(tmpTimeDataArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpTimeDataArray")
             endif

             nn = ubound(timeDataArray,2)
             allocate(tmpTimeDataArray(nTimeStepsRestart,nn), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpTimeDataArray")

             do j=1,nn
               do i=1,nTimeStepsRestart
                 tmpTimeDataArray(i,j) = timeDataArray(i,j)
               enddo
             enddo

             if(allocated(timeDataArray)) then
               deallocate(timeDataArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeDataArray")
             endif

             call allocTimeArrays(nTimeStepsRestart+nTimeStepsFine)

             do i=1,nTimeStepsRestart
               timeArray(i) = tmpTimeArray(i)
             enddo

             do j=1,nn
               do i=1,nTimeStepsRestart
                 timeDataArray(i,j) = tmpTimeDataArray(i,j)
               enddo
             enddo

             deallocate(tmpTimeArray, tmpTimeDataArray, stat=ierr)
             if(ierr /= 0)                                 &
               call terminate("sumb_runIteration",         &
                              "Memory deallocation failure &
                              &for tmpTimeArray and tmpTime&
                              &DataArray")
           endif

           if(storeConvInnerIter) then
             if(allocated(tmpConvArray)) then
               deallocate(tmpConvArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpConvArray")
             endif

             mm = ubound(convArray,2)
             nn = ubound(convArray,3)
             allocate(tmpConvArray(0:nIterOld,mm,nn), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpConvArray")

             do k=1,nn
               do j=1,mm
                 do i=0,nIterOld
                   tmpConvArray(i,j,k) = convArray(i,j,k)
                 enddo
               enddo
             enddo

             if(allocated(convArray)) then
               deallocate(convArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for convArray")
             endif

             kk = nTimeStepsFine
             if(equationMode == steady) kk = 1
             call allocConvArrays(nIterOld+kk*nCycles)

             do k=1,nn
               do j=1,mm
                 do i=0,nIterOld
                   ConvArray(i,j,k) = tmpConvArray(i,j,k)
                 enddo
               enddo
             enddo

             deallocate(tmpConvArray, stat=ierr)
             if(ierr /= 0)                                 &
               call terminate("sumb_runIteration",         &
                              "Memory deallocation failure &
                              &for tmpConvArray")

           endif
         endif

       else

         ! More unsteady time steps or multigrid iterations
         ! in the same session.

         ! If unsteady, reset the number of unsteady time steps.
         ! If steady, reset the number of multigrid cycles on the
         ! finest mesh.

         if(equationMode == unsteady .and. nIterCpl /= nTimeStepsFine) &
           nTimeStepsFine = nIterCpl
         if(equationMode == steady .and. nIterCpl /= nCycles) &
           nCycles = nIterCpl

         ! Reallocate the arrays for convergence history and time data
         ! with a new size, after storing values from previous sessions.

         if(myID == 0) then
           if(equationMode == unsteady) then
             if(allocated(tmpTimeArray)) then
               deallocate(tmpTimeArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpTimeArray")
             endif

             ll = ubound(timeArray, 1)
             allocate(tmpTimeArray(ll), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpTimeArray")

             do i=1,ll
               tmpTimeArray(i) = timeArray(i)
             enddo

             if(allocated(timeArray)) then
               deallocate(timeArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeArray")
             endif

             if(allocated(tmpTimeDataArray)) then
               deallocate(tmpTimeDataArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpTimeDataArray")
             endif

             mm = ubound(timeDataArray,1)
             nn = ubound(timeDataArray,2)
             allocate(tmpTimeDataArray(mm,nn), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpTimeDataArray")

             do j=1,nn
               do i=1,mm
                 tmpTimeDataArray(i,j) = timeDataArray(i,j)
               enddo
             enddo

             if(allocated(timeDataArray)) then
               deallocate(timeDataArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for timeDataArray")
             endif

             call allocTimeArrays(mm+nTimeStepsFine)

             do i=1,ll
               timeArray(i) = tmpTimeArray(i)
             enddo

             do j=1,nn
               do i=1,mm
                 timeDataArray(i,j) = tmpTimeDataArray(i,j)
               enddo
             enddo

             deallocate(tmpTimeArray, tmpTimeDataArray, stat=ierr)
             if(ierr /= 0)                                 &
               call terminate("sumb_runIteration",         &
                              "Memory deallocation failure &
                              &for tmpTimeArray and tmpTime&
                              &DataArray")
           endif

           if(storeConvInnerIter) then
             if(allocated(tmpConvArray)) then
               deallocate(tmpConvArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for tmpConvArray")
             endif

             ll = ubound(convArray,1)
             mm = ubound(convArray,2)
             nn = ubound(convArray,3)
             allocate(tmpConvArray(0:ll,mm,nn), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("sumb_runIteration",       &
                              "Memory allocation failure &
                              &for tmpConvArray")

             do k=1,nn
               do j=1,mm
                 do i=0,ll
                   tmpConvArray(i,j,k) = convArray(i,j,k)
                 enddo
               enddo
             enddo

             if(allocated(convArray)) then
               deallocate(convArray, stat=ierr)
               if(ierr /= 0)                                 &
                 call terminate("sumb_runIteration",         &
                                "Memory deallocation failure &
                                &for convArray")
             endif

             kk = nTimeStepsFine
             if(equationMode == steady) kk = 1
             call allocConvArrays(ll+kk*nCycles)

             do k=1,nn
               do j=1,mm
                 do i=0,ll
                   ConvArray(i,j,k) = tmpConvArray(i,j,k)
                 enddo
               enddo
             enddo

             deallocate(tmpConvArray, stat=ierr)
             if(ierr /= 0)                                 &
               call terminate("sumb_runIteration",         &
                              "Memory deallocation failure &
                              &for tmpConvArray")

           endif
         endif

         ! Reinitialize iteration variables.

         mgStartlevel = 1
         nIterOld = nIterCur
         nIterCur = 0
         iterTot = 0

         ! Update the number of unsteady time steps at the restart.

         nTimeStepsRestart = nTimeStepsRestart + timeStepUnsteady

         ! Reinitialize the time step at the current run.

         timeStepUnsteady = 0

         ! Update the time at the restart.

         timeUnsteadyRestart = timeUnsteadyRestart + timeUnsteady

         ! Reinitialize the time at the current run.

         timeUnsteady = 0.

       endif

       ! If the normal momentum equation should be used to determine
       ! the pressure in the halo for inviscid walls, find out if there
       ! actually are inviscid walls. If so, set the logical
       ! exchangePressureEarly to .true.; otherwise set it to .false.

       if(wallBcTreatment == normalMomentum) then
         exchangePressureEarly = EulerWallsPresent()
       else
         exchangePressureEarly = .false.
       endif

       ! Connect the kill signals with the appropriate functions.
       ! Initialize localSignal for safety.
       ! Only if signalling is supported.

#ifndef USE_NO_SIGNALS
       localSignal = noSignal
       call connect_signals
#endif

       ! Determine the reference time for the solver.

       t0Solver = mpi_wtime()

       ! Set timeUnsteady to zero; this is amount of time simulated
       ! in unsteady mode.

       timeUnsteady = zero

       ! Initialize PV3 routines if PV3 support is required and has not
       ! been called before. In order to make sure that the overset 
       ! iblanking works properly with pV3, set groundLevel to 1 so that
       ! the iblank arrays are allocated with the maximum size.

#ifdef USE_PV3
       if (.not. PV3Initialized) then
         groundLevel = 1
         call initializePV3
         PV3Initialized = .true.
       end if
#endif

       ! Solve either the steady or the unsteady equations for this
       ! grid level. The time spectral method can be considered as
       ! a special kind of steady mode.

       select case (equationMode)
         case (steady, timeSpectral)
           call solverSteady

         case (unsteady)
           select case (timeIntegrationScheme)
               case (BDF)
                 call solverUnsteadyBDF

               case (explicitRK)
                 call solverUnsteadyExplicitRK

               case (implicitRK)
                 call solverUnsteadyImplicitRK
             end select

       end select

       end subroutine sumb_runIteration

!      ==================================================================

       subroutine sumb_writeVolumeSolutionFile
!
!      ******************************************************************
!      *                                                                *
!      * sumb_writeVolumeSolutionFile writes the current state of the   *
!      * volume flow solution to a file.                                *
!      *                                                                *
!      ******************************************************************
!
       use inputMotion
       use iteration
       use monitor
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if(changing_Grid .or. gridMotionSpecified) then
         writeGrid = .true.
       else
         writeGrid = .false.
       endif
       writeVolume = .true.
       writeSurface = .false.

       call writeSol

       end subroutine sumb_writeVolumeSolutionFile

!      ==================================================================

       subroutine sumb_writeSurfaceSolutionFile
!
!      ******************************************************************
!      *                                                                *
!      * sumb_writeSurfaceSolutionFile writes the current state of the  *
!      * surface flow solution to a file.                               *
!      *                                                                *
!      ******************************************************************
!
       use monitor
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       writeGrid = .false.
       writeVolume = .false.
       writeSurface = .true.

       call writeSol

       end subroutine sumb_writeSurfaceSolutionFile

!      ==================================================================

       subroutine sumb_finalize
!
!      ******************************************************************
!      *                                                                *
!      * sumb_finalize makes this particular instance of SUmb exit      *
!      * gracefully in case of coupled simulation.                      *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! First part to release the memory.

       call releaseMemoryPart1

       ! Check if for the time spectral method additional solution
       ! files must be written.

       if(equationMode == timeSpectral) then

         if( writeUnsteadyRestartSpectral ) &
           call writeUnsteadyFromSpectral

         if(writeUnsteadyVolSpectral .or. &
            writeUnsteadySurfSpectral)    &
           call writeInterpolFromSpectral

       endif

       ! Second part to release the memory.

       call releaseMemoryPart2

       ! Write the parameters used for this run to stdout.

       call writeInputParam

       end subroutine sumb_finalize

!      ==================================================================

       subroutine sumb_getPointSize(nPoints, geomName, includeHalosIn)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getPointSize provides the coupler with the number of      *
!      * points (x,y,z) belonging to the family with the given name,    *
!      * where information from other solvers is required.              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use blockPointers
       use couplerParam
       use cgnsGrid
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(in) :: geomName
       integer(kind=intType), intent(out) :: nPoints
       logical, intent(in), optional      :: includeHalosIn
!
!      Local variables.
!
       integer(intType) :: nLevels, level, nn, j, mm
       integer(intType) :: iNum, jNum, kNum
       character(len=maxCGNSNameLen) :: famName
       character(len=maxCplNameLen) :: trimGeomName

       logical :: includeHalos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not to include the cell halo's of the
       ! subface.

       if( present(includeHalosIn) ) then
         includeHalos = includeHalosIn
       else
         includeHalos = .false.
       endif


       nPoints = 0

       nLevels = 1
       if(cplGetCoarseSol) nLevels = ubound(flowDoms,2)

       levelLoop: do level=1,nLevels
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! About the time-spectral method, only the zero mode is counted.

           call setPointers(nn, level, 1)

           bocosLoop: do j=1,nBocos
             if(BCType(j) == DomainInterfaceAll    .or. &
                BCType(j) == DomainInterfaceRhoUVW .or. &
                BCType(j) == DomainInterfaceP      .or. &
                BCType(j) == DomainInterfaceRho    .or. &
                BCType(j) == DomainInterfaceTotal) then

               mm = groupNum(j)

               famName = cgnsFamilies(mm)%familyName
               call convertToLowerCase(famName)
               famName = adjustl(famName)
               famName = trim(famName)

               trimGeomName = geomName
               call convertToLowerCase(trimGeomName)
               trimGeomName = adjustl(trimGeomName)
               trimGeomName = trim(trimGeomName)

               if(famName == trimGeomName) then

                 if( includeHalos ) then
                   iNum = iabs(icEnd(j)-icBeg(j)) + 1
                   jNum = iabs(jcEnd(j)-jcBeg(j)) + 1
                   kNum = iabs(kcEnd(j)-kcBeg(j)) + 1
                 else
                   iNum = iabs(inEnd(j)-inBeg(j))
                   jNum = iabs(jnEnd(j)-jnBeg(j))
                   kNum = iabs(knEnd(j)-knBeg(j))
 		 endif

                 select case (BCFaceId(j))
                   case (iMin, iMax)
                     iNum = 1

                   case (jMin, jMax)
                     jNum = 1

                   case (kMin, kMax)
                     kNum = 1

                 end select

                 nPoints = nPoints + iNum*jNum*kNum
               endif
             endif
           enddo bocosLoop
         enddo domainLoop
       enddo levelLoop

       end subroutine sumb_getPointSize

!      ==================================================================

       subroutine sumb_getPointGeom(xyz, geomName, includeHalosIn)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getPointGeom provides the coupler with the local list of  *
!      * (x,y,z) at which information from other solvers is required.   *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use blockPointers
       use constants
       use couplerParam
       use cgnsGrid
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: xyz
       character(len=*), intent(in) :: geomName
       logical, intent(in), optional :: includeHalosIn
!
!      Local variables.
!
       integer(kind=intType) :: nCnt, nLevels, level, nn, nBFs, mm
       integer(kind=intType) :: iiI, iiF, jjI, jjF, kkI, kkF       
       integer(kind=intType) :: i, j, k, l, m, n
       character(len=maxCGNSNameLen) :: famName
       character(len=maxCplNameLen) :: trimGeomName

       logical :: includeHalos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not to include the cell halo's of the
       ! subface.

       if( present(includeHalosIn) ) then
         includeHalos = includeHalosIn
       else
         includeHalos = .false.
       endif

       nCnt = 0

       nLevels = 1
       if(cplGetCoarseSol) nLevels = ubound(flowDoms,2)

       levelLoop: do level=1,nLevels
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! About the time-spectral method, only the zero mode is counted.

           call setPointers(nn, level, 1)

           bocosLoop: do nBFs=1,nBocos
             if(BCType(nBFs) == DomainInterfaceAll    .or. &
                BCType(nBFs) == DomainInterfaceRhoUVW .or. &
                BCType(nBFs) == DomainInterfaceP      .or. &
                BCType(nBFs) == DomainInterfaceRho    .or. &
                BCType(nBFs) == DomainInterfaceTotal) then

               mm = groupNum(nBFs)

               famName = cgnsFamilies(mm)%familyName
               call convertToLowerCase(famName)
               famName = adjustl(famName)
               famName = trim(famName)

               trimGeomName = geomName
               call convertToLowerCase(trimGeomName)
               trimGeomName = adjustl(trimGeomName)
               trimGeomName = trim(trimGeomName)

               if(famName == trimGeomName) then

                 ! Find the cell centers.

		 if( includeHalos ) then
                   iiI = icBeg(nBFs)
                   iiF = icEnd(nBFs)
                   jjI = jcBeg(nBFs)
                   jjF = jcEnd(nBFs)
                   kkI = kcBeg(nBFs)
                   kkF = kcEnd(nBFs)
		 else
                   iiI = inBeg(nBFs) + 1
                   iiF = inEnd(nBFs)
                   jjI = jnBeg(nBFs) + 1
                   jjF = jnEnd(nBFs)
                   kkI = knBeg(nBFs) + 1
                   kkF = knEnd(nBFs)
                 endif

                 select case (BCFaceId(nBFs))
                   case (iMin)
                     iiI = 1
                     iiF = 1

                   case(iMax)
                     iiI = ie
                     iiF = ie

                   case(jMin)
                     jjI = 1
                     jjF = 1

                   case(jMax)
                     jjI = je
                     jjF = je

                   case(kMin)
                     kkI = 1
                     kkF = 1

                   case(kMax)
                     kkI = ke
                     kkF = ke

                 end select

                 do k=kkI,kkF
                   n = k - 1
                   do j=jjI,jjF
                     m = j - 1
                     do i=iiI,iiF
                       l = i - 1

                       nCnt = nCnt + 1

                       xyz(1,nCnt) = &
                           eighth*(x(i,j,k,1) + x(i,m,k,1) &
                         +         x(i,m,n,1) + x(i,j,n,1) &
                         +         x(l,j,k,1) + x(l,m,k,1) &
                         +         x(l,m,n,1) + x(l,j,n,1))
                       xyz(2,nCnt) = &
                           eighth*(x(i,j,k,2) + x(i,m,k,2) &
                         +         x(i,m,n,2) + x(i,j,n,2) &
                         +         x(l,j,k,2) + x(l,m,k,2) &
                         +         x(l,m,n,2) + x(l,j,n,2))
                       xyz(3,nCnt) = &
                           eighth*(x(i,j,k,3) + x(i,m,k,3) &
                         +         x(i,m,n,3) + x(i,j,n,3) &
                         +         x(l,j,k,3) + x(l,m,k,3) &
                         +         x(l,m,n,3) + x(l,j,n,3))
                     enddo
                   enddo
                 enddo
               endif
             endif
           enddo bocosLoop
         enddo domainLoop
       enddo levelLoop

       end subroutine sumb_getPointGeom

!      ==================================================================

       subroutine sumb_getPointData(data, dataNames, nDataNames, &
                                    geomName, includeHalosIn)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getPointData provides the coupler with the requested flow *
!      * variables at interface points.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use blockPointers
       use cgnsGrid
       use cgnsNames
       use communication
       use constants
       use couplerParam
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use iteration

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
       logical, intent(in), optional :: includeHalosIn
!
!      Local variables.
!
       integer(kind=intType) :: nCnt, nLevels, level, nn, nBFs, mm
       integer(kind=intType) :: iiI, iiF, jjI, jjF, kkI, kkF
       integer(kind=intType) :: nfI, nfO
       integer(kind=intType) :: i, j, k, ll
       real(kind=realType) :: vecMag

       character(len=maxCGNSNameLen) :: famName
       character(len=maxCplNameLen) :: trimGeomName, tmpName

       logical :: includeHalos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not to include the cell halo's of the
       ! subface.

       if( present(includeHalosIn) ) then
         includeHalos = includeHalosIn
       else
         includeHalos = .false.
       endif

       nCnt = 0

       nLevels = 1
       if(cplGetCoarseSol) nLevels = ubound(flowDoms,2)

       levelLoop: do level=1,nLevels
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! About the time-spectral method, only the zero mode is counted.

           call setPointers(nn, level, 1)

           bocosLoop: do nBFs=1,nBocos
             if(BCType(nBFs) == DomainInterfaceAll    .or. &
                BCType(nBFs) == DomainInterfaceRhoUVW .or. &
                BCType(nBFs) == DomainInterfaceP      .or. &
                BCType(nBFs) == DomainInterfaceRho    .or. &
                BCType(nBFs) == DomainInterfaceTotal) then

               mm = groupNum(nBFs)

               famName = cgnsFamilies(mm)%familyName
               call convertToLowerCase(famName)
               famName = adjustl(famName)
               famName = trim(famName)

               trimGeomName = geomName
               call convertToLowerCase(trimGeomName)
               trimGeomName = adjustl(trimGeomName)
               trimGeomName = trim(trimGeomName)

               if(famName == trimGeomName) then

                 ! Find the cell centers.

		 if( includeHalos ) then
                   iiI = icBeg(nBFs)
                   iiF = icEnd(nBFs)
                   jjI = jcBeg(nBFs)
                   jjF = jcEnd(nBFs)
                   kkI = kcBeg(nBFs)
                   kkF = kcEnd(nBFs)
		 else
                   iiI = inBeg(nBFs) + 1
                   iiF = inEnd(nBFs)
                   jjI = jnBeg(nBFs) + 1
                   jjF = jnEnd(nBFs)
                   kkI = knBeg(nBFs) + 1
                   kkF = knEnd(nBFs)
                 endif

                 select case (BCFaceId(nBFs))
                   case (iMin)
                     iiI = 1
                     iiF = 1

                   case(iMax)
                     iiI = ie
                     iiF = ie

                   case(jMin)
                     jjI = 1
                     jjF = 1

                   case(jMax)
                     jjI = je
                     jjF = je

                   case(kMin)
                     kkI = 1
                     kkF = 1

                   case(kMax)
                     kkI = ke
                     kkF = ke

                 end select

                 do k=kkI,kkF
                   do j=jjI,jjF
                     do i=iiI,iiF

                       nCnt = nCnt + 1

                       select case (BCFaceID(nBFs))
                         case (iMin, iMax)
                           nfO = k
                           nfI = j

                         case (jMin, jMax)
                           nfO = k
                           nfI = i

                         case (kMin, kMax)
                           nfO = j
                           nfI = i

                       end select

                       do ll=1,nDataNames
                         tmpName = dataNames(ll)
                       ! call convertToLowerCase(tmpName)
                         tmpName = adjustl(tmpName)
                         tmpName = trim(tmpName)
                         if(tmpName(1:3) == 'AVG' .or. &
                            tmpName(1:3) == 'Avg' .or. &
                            tmpName(1:3) == 'avg')     &
                           tmpName = tmpName(4:)

                         if     (tmpName == cgnsDensity) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             data(ll,nCnt) = BCData(nBFs)%rho(nfI,nfO)
                           else
                             data(ll,nCnt) = w(i,j,k,irho)
                           endif

                         else if(tmpName == cgnsVelX) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             data(ll,nCnt) = BCData(nBFs)%velX(nfI,nfO)
                           else
                             data(ll,nCnt) = w(i,j,k,ivx)
                           endif

                         else if(tmpName == cgnsVelY) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             data(ll,nCnt) = BCData(nBFs)%velY(nfI,nfO)
                           else
                             data(ll,nCnt) = w(i,j,k,ivy)
                           endif

                         else if(tmpName == cgnsVelZ) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             data(ll,nCnt) =  BCData(nBFs)%velZ(nfI,nfO)
                           else
                             data(ll,nCnt) = w(i,j,k,ivz)
                           endif

                         else if(tmpName == cgnsTurbSANu .or. &
                                 tmpName == cgnsTurbK) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%turbInlet(nfI,nfO,itu1)
                           else
                             data(ll,nCnt) = w(i,j,k,itu1)
                           endif

                         else if(tmpName == cgnsTurbOmega .or. &
                                 tmpName == cgnsTurbTau   .or. &
                                 tmpName == cgnsTurbEpsilon) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%turbInlet(nfI,nfO,itu2)
                           else
                             data(ll,nCnt) = w(i,j,k,itu2)
                           endif

                         else if(tmpName == cgnsTurbV2) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%turbInlet(nfI,nfO,itu3)
                           else
                             data(ll,nCnt) = w(i,j,k,itu3)
                           endif

                         else if(tmpName == cgnsTurbF) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%turbInlet(nfI,nfO,itu4)
                           else
                             data(ll,nCnt) = w(i,j,k,itu4)
                           endif

                         else if(tmpName == cgnsPressure) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceP) then
                             data(ll,nCnt) = BCData(nBFs)%ps(nfI,nfO)
                           else
                             data(ll,nCnt) = p(i,j,k)
                           endif

                         else if(tmpName == cgnsPTot) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%ptInlet(nfI,nfO)
                           else
                             call computePtot(w(i,j,k,irho), w(i,j,k,ivx), &
                                              w(i,j,k,ivy) , w(i,j,k,ivz), &
                                              p(i,j,k), data(ll,nCnt), 1)
                           endif

                         else if(tmpName == cgnsTTot) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%ttInlet(nfI,nfO)
                           else
                             call computeTtot(w(i,j,k,irho), w(i,j,k,ivx), &
                                              w(i,j,k,ivy) , w(i,j,k,ivz), &
                                              p(i,j,k), data(ll,nCnt), 1)
                           endif

                         else if(tmpName == cgnsVelVecX) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%flowXDirInlet(nfI,nfO)
                           else
                             vecMag = sqrt(w(i,j,k,ivx)**2 &
                                         + w(i,j,k,ivy)**2 &
                                         + w(i,j,k,ivz)**2)
                             if(vecMag > eps) then
                               data(ll,nCnt) = w(i,j,k,ivx)/vecMag
                             else
                               data(ll,nCnt) = 0.
                             endif
                           endif

                         else if(tmpName == cgnsVelVecY) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%flowYDirInlet(nfI,nfO)
                           else
                             vecMag = sqrt(w(i,j,k,ivx)**2 &
                                         + w(i,j,k,ivy)**2 &
                                         + w(i,j,k,ivz)**2)
                             if(vecMag > eps) then
                               data(ll,nCnt) = w(i,j,k,ivy)/vecMag
                             else
                               data(ll,nCnt) = 0.
                             endif
                           endif

                         else if(tmpName == cgnsVelVecZ) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             data(ll,nCnt) = &
                             BCData(nBFs)%flowZDirInlet(nfI,nfO)
                           else
                             vecMag = sqrt(w(i,j,k,ivx)**2 &
                                         + w(i,j,k,ivy)**2 &
                                         + w(i,j,k,ivz)**2)
                             if(vecMag > eps) then
                               data(ll,nCnt) = w(i,j,k,ivz)/vecMag
                             else
                               data(ll,nCnt) = 0.
                             endif
                           endif

                         else
                           ! Write the message only once.

                           if(i == iiI .and. j == jjI .and. k == kkI) then
                             write(*,*) "********* SUmb point: geomName =", &
                                        geomName, ", myID =", myID," *********"
                             write(*,*) tmpName, ": This variable cannot be &
                                       &provided from the present SUmb interface."
                             write(*,*) 'Zero is provided instead.'
                           endif
                           data(ll,nCnt) = 0.
                         endif

                       enddo
                     enddo
                   enddo
                 enddo
               endif
             endif
           enddo bocosLoop
         enddo domainLoop
       enddo levelLoop

       end subroutine sumb_getPointData

!      ==================================================================

       subroutine sumb_setPointData(data, dataNames, nDataNames, &
                                    geomName, includeHalosIn)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setPointData gets the interpolated solution at the        *
!      * interfaces from the coupler. It will be stored in the arrays   *
!      * under flowDoms%BCData.                                         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use blockPointers
       use cgnsGrid
       use cgnsNames
       use constants
       use couplerParam
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use iteration

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(in) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
       logical, intent(in), optional :: includeHalosIn
!
!      Local variables.
!
       integer(kind=intType) :: nCnt, nLevels, level, nn, nBFs, mm
       integer(kind=intType) :: iiI, iiF, jjI, jjF, kkI, kkF
       integer(kind=intType) :: nfI, nfO
       integer(kind=intType) :: i, j, k, ll
       integer(kind=intType) :: nfIB, nfIBH, nfIE, nfIEH,           &
                                nfOB, nfOBH, nfOE, nfOEH,           &
                                nfIH, nfI1, nfI2, nfOH, nfO1, nfO2, &
                                iiH, jjH, kkH
       integer(kind=intType) :: ndir, ibe
       integer(kind=intType) :: kco, lco
       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd, iiMax, jjMax
       integer(kind=intType) :: levm1
       integer(kind=intType), dimension(:,:), pointer :: iFine, jFine
       integer :: ierr

       real(kind=realType), allocatable, dimension(:,:,:) :: &
                            xxmIH, xxmI1, xxmI2, xxmOH, xxmO1, xxmO2
       real(kind=realType), allocatable, dimension(:,:) :: dsI1, dsI2, &
                                                           dsO1, dsO2
       real(kind=realType) :: factor1, factor2
       real(kind=realType) :: var
       real(kind=realType), dimension(3) :: dir

       character(len=maxCGNSNameLen) :: famName
       character(len=maxCplNameLen) :: trimGeomName, tmpName

       logical :: includeHalos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not to include the cell halo's of the
       ! subface.

       if( present(includeHalosIn) ) then
         includeHalos = includeHalosIn
       else
         includeHalos = .false.
       endif

       nCnt = 0

       nLevels = 1
       if(cplGetCoarseSol) nLevels = ubound(flowDoms,2)

       levelLoop: do level=1,nLevels
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! About the time-spectral method, only the zero mode is counted.

           call setPointers(nn, level, 1)

           bocosLoop: do nBFs=1,nBocos
             if(BCType(nBFs) == DomainInterfaceAll    .or. &
                BCType(nBFs) == DomainInterfaceRhoUVW .or. &
                BCType(nBFs) == DomainInterfaceP      .or. &
                BCType(nBFs) == DomainInterfaceRho    .or. &
                BCType(nBFs) == DomainInterfaceTotal) then

               mm = groupNum(nBFs)

               famName = cgnsFamilies(mm)%familyName
               call convertToLowerCase(famName)
               famName = adjustl(famName)
               famName = trim(famName)

               trimGeomName = geomName
               call convertToLowerCase(trimGeomName)
               trimGeomName = adjustl(trimGeomName)
               trimGeomName = trim(trimGeomName)

               if(famName == trimGeomName) then

                 ! Find the cell centers.

		 if( includeHalos ) then
                   iiI = icBeg(nBFs)
                   iiF = icEnd(nBFs)
                   jjI = jcBeg(nBFs)
                   jjF = jcEnd(nBFs)
                   kkI = kcBeg(nBFs)
                   kkF = kcEnd(nBFs)
		 else
                   iiI = inBeg(nBFs) + 1
                   iiF = inEnd(nBFs)
                   jjI = jnBeg(nBFs) + 1
                   jjF = jnEnd(nBFs)
                   kkI = knBeg(nBFs) + 1
                   kkF = knEnd(nBFs)
                 endif

                 select case (BCFaceId(nBFs))
                   case (iMin)
                     iiI = 1
                     iiF = 1

                   case(iMax)
                     iiI = ie
                     iiF = ie

                   case(jMin)
                     jjI = 1
                     jjF = 1

                   case(jMax)
                     jjI = je
                     jjF = je

                   case(kMin)
                     kkI = 1
                     kkF = 1

                   case(kMax)
                     kkI = ke
                     kkF = ke

                 end select

                 do k=kkI,kkF
                   do j=jjI,jjF
                     do i=iiI,iiF

                       nCnt = nCnt + 1

                       select case (BCFaceID(nBFs))
                         case (iMin, iMax)
                           nfO = k
                           nfI = j

                         case (jMin, jMax)
                           nfO = k
                           nfI = i

                         case (kMin, kMax)
                           nfO = j
                           nfI = i

                       end select

                       do ll=1,nDataNames
                         tmpName = dataNames(ll)
                       ! call convertToLowerCase(tmpName)
                         tmpName = adjustl(tmpName)
                         tmpName = trim(tmpName)
                         if(tmpName(1:3) == 'AVG' .or. &
                            tmpName(1:3) == 'Avg' .or. &
                            tmpName(1:3) == 'avg')     &
                           tmpName = tmpName(4:)

                         if     (tmpName == cgnsDensity) then
                           if(BCType(nBFs) == DomainInterfaceAll    .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             BCData(nBFs)%rho(nfI,nfO)  = data(ll,nCnt)
                           else
                             w(i,j,k,irho) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsVelX) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             BCData(nBFs)%velX(nfI,nfO) = data(ll,nCnt)
                           else if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowXDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,ivx) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsVelY) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             BCData(nBFs)%velY(nfI,nfO) = data(ll,nCnt)
                           else if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowYDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,ivy) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsVelZ) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceRhoUVW) then
                             BCData(nBFs)%velZ(nfI,nfO) = data(ll,nCnt)
                           else if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowZDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,ivz) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsTurbSANu .or. &
                                 tmpName == cgnsTurbK) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             BCData(nBFs)%turbInlet(nfI,nfO,itu1) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,itu1) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsTurbOmega .or. &
                                 tmpName == cgnsTurbTau   .or. &
                                 tmpName == cgnsTurbEpsilon) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             BCData(nBFs)%turbInlet(nfI,nfO,itu2) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,itu2) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsTurbV2) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             BCData(nBFs)%turbInlet(nfI,nfO,itu3) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,itu3) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsTurbF) then
                           if(BCType(nBFs) /= DomainInterfaceP) then
                             BCData(nBFs)%turbInlet(nfI,nfO,itu4) = &
                             data(ll,nCnt)
                           else
                             w(i,j,k,itu4) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsPressure) then
                           if(BCType(nBFs) == DomainInterfaceAll .or. &
                              BCType(nBFs) == DomainInterfaceP) then
                             BCData(nBFs)%ps(nfI,nfO) = data(ll,nCnt)
                           else
                             p(i,j,k) = data(ll,nCnt)
                           endif

                         else if(tmpName == cgnsPTot) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%ptInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             ! Write the message only once.

                             if(i == iiI .and. j == jjI .and. k == kkI) then
                               write(*,*) 'Warning in sumb_setPointData : pTot is &
                                          &obtained although it is not necessary.'
                               write(*,*) 'Information is ignored.'
                             endif
                           endif

                         else if(tmpName == cgnsTTot) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%ttInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             ! Write the message only once.

                             if(i == iiI .and. j == jjI .and. k == kkI) then
                               write(*,*) 'Warning in sumb_setPointData : TTot is &
                                          &obtained although it is not necessary.'
                               write(*,*) 'Information is ignored.'
                             endif
                           endif

                         else if(tmpName == cgnsVelVecX) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowXDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             ! Write the message only once.

                             if(i == iiI .and. j == jjI .and. k == kkI) then
                               write(*,*) 'Warning in sumb_setPointData : VelVecX is &
                                          &obtained although it is not necessary.'
                               write(*,*) 'Information is ignored.'
                             endif
                           endif

                         else if(tmpName == cgnsVelVecY) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowYDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             ! Write the message only once.

                             if(i == iiI .and. j == jjI .and. k == kkI) then
                               write(*,*) 'Warning in sumb_setPointData : VelVecY is &
                                          &obtained although it is not necessary.'
                               write(*,*) 'Information is ignored.'
                             endif
                           endif

                         else if(tmpName == cgnsVelVecZ) then
                           if(BCType(nBFs) == DomainInterfaceTotal) then
                             BCData(nBFs)%flowZDirInlet(nfI,nfO) = &
                             data(ll,nCnt)
                           else
                             ! Write the message only once.

                             if(i == iiI .and. j == jjI .and. k == kkI) then
                               write(*,*) 'Warning in sumb_setPointData : VelVecZ is &
                                          &obtained although it is not necessary.'
                               write(*,*) 'Information is ignored.'
                             endif
                           endif

                         endif
                       enddo

                       ! Compute some additional info in case of the total condition.

                       if(BCType(nBFs) == DomainInterfaceTotal) then

                         ! Compute the total enthalpy from the given
                         ! total temperature.

                         call computeHtot(BCData(nBFs)%ttInlet(nfI,nfO), &
                                          BCData(nBFs)%htInlet(nfI,nfO))

                         ! Determine the unit vector of the flow direction.

                         dir(1) = BCData(nBFs)%flowXdirInlet(nfI,nfO)
                         dir(2) = BCData(nBFs)%flowYdirInlet(nfI,nfO)
                         dir(3) = BCData(nBFs)%flowZdirInlet(nfI,nfO)

                         var = one/max(eps,sqrt(dir(1)**2 + dir(2)**2 &
                                              + dir(3)**2))

                         BCData(nBFs)%flowXdirInlet(nfI,nfO) = var*dir(1)
                         BCData(nBFs)%flowYdirInlet(nfI,nfO) = var*dir(2)
                         BCData(nBFs)%flowZdirInlet(nfI,nfO) = var*dir(3)

                       endif
                     enddo
                   enddo
                 enddo

                 ! Extrapolate values for halos at the interface if these
                 ! values were not requested.

		 if(.not. includeHalos ) then

                   select case (BCFaceID(nBFs))
                     case (iMin, iMax)
                       nfOB  = kkI
                       nfOBH = kkI - 1
                       nfOE  = kkF
                       nfOEH = kkF + 1

                       nfIB  = jjI
                       nfIBH = jjI - 1
                       nfIE  = jjF
                       nfIEH = jjF + 1

                       iiH = iiI

                       allocate(xxmIH(nfIBH:nfIEH,2,3), xxmI1(nfIBH:nfIEH,2,3), &
                                xxmI2(nfIBH:nfIEH,2,3), dsI1(nfIBH:nfIEH,2),    &
                                dsI2(nfIBH:nfIEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmI and dsI")

                       allocate(xxmOH(nfOBH:nfOEH,2,3), xxmO1(nfOBH:nfOEH,2,3), &
                                xxmO2(nfOBH:nfOEH,2,3), dsO1(nfOBH:nfOEH,2),    &
                                dsO2(nfOBH:nfOEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmO and dsO")

                       do ndir=1,3
                         do nfI=nfIBH,nfIEH
                           xxmIH(nfI,1,ndir) = eighth*( x(iiH  ,nfI  ,nfOBH-1,ndir) &
                                                      + x(iiH  ,nfI  ,nfOBH  ,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH-1,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH  ,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH-1,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH  ,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH-1,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH  ,ndir) )
                           xxmIH(nfI,2,ndir) = eighth*( x(iiH  ,nfI  ,nfOEH-1,ndir) &
                                                      + x(iiH  ,nfI  ,nfOEH  ,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH-1,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH  ,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH-1,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH  ,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH-1,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH  ,ndir) )

                           xxmI1(nfI,1,ndir) = eighth*( x(iiH  ,nfI  ,nfOBH  ,ndir) &
                                                      + x(iiH  ,nfI  ,nfOBH+1,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH  ,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH+1,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH  ,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH+1,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH  ,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH+1,ndir) )
                           xxmI1(nfI,2,ndir) = eighth*( x(iiH  ,nfI  ,nfOEH-2,ndir) &
                                                      + x(iiH  ,nfI  ,nfOEH-1,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH-2,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH-1,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH-2,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH-1,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH-2,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH-1,ndir) )

                           xxmI2(nfI,1,ndir) = eighth*( x(iiH  ,nfI  ,nfOBH+1,ndir) &
                                                      + x(iiH  ,nfI  ,nfOBH+2,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH+1,ndir) &
                                                      + x(iiH  ,nfI-1,nfOBH+2,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH+1,ndir) &
                                                      + x(iiH-1,nfI  ,nfOBH+2,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH+1,ndir) &
                                                      + x(iiH-1,nfI-1,nfOBH+2,ndir) )
                           xxmI2(nfI,2,ndir) = eighth*( x(iiH  ,nfI  ,nfOEH-3,ndir) &
                                                      + x(iiH  ,nfI  ,nfOEH-2,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH-3,ndir) &
                                                      + x(iiH  ,nfI-1,nfOEH-2,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH-3,ndir) &
                                                      + x(iiH-1,nfI  ,nfOEH-2,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH-3,ndir) &
                                                      + x(iiH-1,nfI-1,nfOEH-2,ndir) )
                         enddo

                         do nfO=nfOBH,nfOEH
                           xxmOH(nfO,1,ndir) = eighth*( x(iiH  ,nfIBH-1,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH  ,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH-1,nfO-1,ndir) &
                                                      + x(iiH  ,nfIBH  ,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH-1,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH  ,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH-1,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH  ,nfO-1,ndir) )
                           xxmOH(nfO,2,ndir) = eighth*( x(iiH  ,nfIEH-1,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH  ,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH-1,nfO-1,ndir) &
                                                      + x(iiH  ,nfIEH  ,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH-1,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH  ,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH-1,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH  ,nfO-1,ndir) )

                           xxmO1(nfO,1,ndir) = eighth*( x(iiH  ,nfIBH  ,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH+1,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH  ,nfO-1,ndir) &
                                                      + x(iiH  ,nfIBH+1,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH  ,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH+1,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH  ,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH+1,nfO-1,ndir) )
                           xxmO1(nfO,2,ndir) = eighth*( x(iiH  ,nfIEH-2,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH-1,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH-2,nfO-1,ndir) &
                                                      + x(iiH  ,nfIEH-1,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH-2,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH-1,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH-2,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH-1,nfO-1,ndir) )

                           xxmO2(nfO,1,ndir) = eighth*( x(iiH  ,nfIBH+1,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH+2,nfO  ,ndir) &
                                                      + x(iiH  ,nfIBH+1,nfO-1,ndir) &
                                                      + x(iiH  ,nfIBH+2,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH+1,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH+2,nfO  ,ndir) &
                                                      + x(iiH-1,nfIBH+1,nfO-1,ndir) &
                                                      + x(iiH-1,nfIBH+2,nfO-1,ndir) )
                           xxmO2(nfO,2,ndir) = eighth*( x(iiH  ,nfIEH-3,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH-2,nfO  ,ndir) &
                                                      + x(iiH  ,nfIEH-3,nfO-1,ndir) &
                                                      + x(iiH  ,nfIEH-2,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH-3,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH-2,nfO  ,ndir) &
                                                      + x(iiH-1,nfIEH-3,nfO-1,ndir) &
                                                      + x(iiH-1,nfIEH-2,nfO-1,ndir) )
                         enddo
                       enddo

                       do ibe=1,2
                         do nfI=nfIBH,nfIEH
                           dsI1(nfI,ibe) = sqrt( (xxmIH(nfI,ibe,1) - xxmI1(nfI,ibe,1))**2 &
                                               + (xxmIH(nfI,ibe,2) - xxmI1(nfI,ibe,2))**2 &
                                               + (xxmIH(nfI,ibe,3) - xxmI1(nfI,ibe,3))**2 )
                           dsI2(nfI,ibe) = sqrt( (xxmI1(nfI,ibe,1) - xxmI2(nfI,ibe,1))**2 &
                                               + (xxmI1(nfI,ibe,2) - xxmI2(nfI,ibe,2))**2 &
                                               + (xxmI1(nfI,ibe,3) - xxmI2(nfI,ibe,3))**2 )
                         enddo

                         do nfO=nfOBH,nfOEH
                           dsO1(nfO,ibe) = sqrt( (xxmOH(nfO,ibe,1) - xxmO1(nfO,ibe,1))**2 &
                                               + (xxmOH(nfO,ibe,2) - xxmO1(nfO,ibe,2))**2 &
                                               + (xxmOH(nfO,ibe,3) - xxmO1(nfO,ibe,3))**2 )
                           dsO2(nfO,ibe) = sqrt( (xxmO1(nfO,ibe,1) - xxmO2(nfO,ibe,1))**2 &
                                               + (xxmO1(nfO,ibe,2) - xxmO2(nfO,ibe,2))**2 &
                                               + (xxmO1(nfO,ibe,3) - xxmO2(nfO,ibe,3))**2 )
                         enddo
                       enddo

                       deallocate(xxmIH, xxmI1, xxmI2, xxmOH, xxmO1, xxmO2, stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory deallocation failure for xxmI and xxmO")

                     case (jMin, jMax)
                       nfOB  = kkI
                       nfOBH = kkI - 1
                       nfOE  = kkF
                       nfOEH = kkF + 1

                       nfIB  = iiI
                       nfIBH = iiI - 1
                       nfIE  = iiF
                       nfIEH = iiF + 1

                       jjH = jjI

                       allocate(xxmIH(nfIBH:nfIEH,2,3), xxmI1(nfIBH:nfIEH,2,3), &
                                xxmI2(nfIBH:nfIEH,2,3), dsI1(nfIBH:nfIEH,2),    &
                                dsI2(nfIBH:nfIEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmI and dsI")

                       allocate(xxmOH(nfOBH:nfOEH,2,3), xxmO1(nfOBH:nfOEH,2,3), &
                                xxmO2(nfOBH:nfOEH,2,3), dsO1(nfOBH:nfOEH,2),    &
                                dsO2(nfOBH:nfOEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmO and dsO")

                       do ndir=1,3
                         do nfI=nfIBH,nfIEH
                           xxmIH(nfI,1,ndir) = eighth*( x(nfI  ,jjH  ,nfOBH-1,ndir) &
                                                      + x(nfI  ,jjH  ,nfOBH  ,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH-1,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH  ,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH-1,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH  ,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH-1,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH  ,ndir) )
                           xxmIH(nfI,2,ndir) = eighth*( x(nfI  ,jjH  ,nfOEH-1,ndir) &
                                                      + x(nfI  ,jjH  ,nfOEH  ,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH-1,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH  ,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH-1,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH  ,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH-1,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH  ,ndir) )

                           xxmI1(nfI,1,ndir) = eighth*( x(nfI  ,jjH  ,nfOBH  ,ndir) &
                                                      + x(nfI  ,jjH  ,nfOBH+1,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH  ,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH+1,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH  ,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH+1,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH  ,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH+1,ndir) )
                           xxmI1(nfI,2,ndir) = eighth*( x(nfI  ,jjH  ,nfOEH-2,ndir) &
                                                      + x(nfI  ,jjH  ,nfOEH-1,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH-2,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH-1,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH-2,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH-1,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH-2,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH-1,ndir) )

                           xxmI2(nfI,1,ndir) = eighth*( x(nfI  ,jjH  ,nfOBH+1,ndir) &
                                                      + x(nfI  ,jjH  ,nfOBH+2,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH+1,ndir) &
                                                      + x(nfI-1,jjH  ,nfOBH+2,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH+1,ndir) &
                                                      + x(nfI  ,jjH-1,nfOBH+2,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH+1,ndir) &
                                                      + x(nfI-1,jjH-1,nfOBH+2,ndir) )
                           xxmI2(nfI,2,ndir) = eighth*( x(nfI  ,jjH  ,nfOEH-3,ndir) &
                                                      + x(nfI  ,jjH  ,nfOEH-2,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH-3,ndir) &
                                                      + x(nfI-1,jjH  ,nfOEH-2,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH-3,ndir) &
                                                      + x(nfI  ,jjH-1,nfOEH-2,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH-3,ndir) &
                                                      + x(nfI-1,jjH-1,nfOEH-2,ndir) )
                         enddo

                         do nfO=nfOBH,nfOEH
                           xxmOH(nfO,1,ndir) = eighth*( x(nfIBH-1,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH  ,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH-1,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH  ,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH-1,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH  ,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH-1,jjH-1,nfO-1,ndir) &
                                                      + x(nfIBH  ,jjH-1,nfO-1,ndir) )
                           xxmOH(nfO,2,ndir) = eighth*( x(nfIEH-1,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH  ,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH-1,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH  ,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH-1,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH  ,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH-1,jjH-1,nfO-1,ndir) &
                                                      + x(nfIEH  ,jjH-1,nfO-1,ndir) )

                           xxmO1(nfO,1,ndir) = eighth*( x(nfIBH  ,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH+1,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH  ,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH+1,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH  ,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH+1,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH  ,jjH-1,nfO-1,ndir) &
                                                      + x(nfIBH+1,jjH-1,nfO-1,ndir) )
                           xxmO1(nfO,2,ndir) = eighth*( x(nfIEH-2,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH-1,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH-2,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH-1,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH-2,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH-1,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH-2,jjH-1,nfO-1,ndir) &
                                                      + x(nfIEH-1,jjH-1,nfO-1,ndir) )

                           xxmO2(nfO,1,ndir) = eighth*( x(nfIBH+1,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH+2,jjH  ,nfO  ,ndir) &
                                                      + x(nfIBH+1,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH+2,jjH  ,nfO-1,ndir) &
                                                      + x(nfIBH+1,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH+2,jjH-1,nfO  ,ndir) &
                                                      + x(nfIBH+1,jjH-1,nfO-1,ndir) &
                                                      + x(nfIBH+2,jjH-1,nfO-1,ndir) )
                           xxmO2(nfO,2,ndir) = eighth*( x(nfIEH-3,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH-2,jjH  ,nfO  ,ndir) &
                                                      + x(nfIEH-3,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH-2,jjH  ,nfO-1,ndir) &
                                                      + x(nfIEH-3,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH-2,jjH-1,nfO  ,ndir) &
                                                      + x(nfIEH-3,jjH-1,nfO-1,ndir) &
                                                      + x(nfIEH-2,jjH-1,nfO-1,ndir) )
                         enddo
                       enddo

                       do ibe=1,2
                         do nfI=nfIBH,nfIEH
                           dsI1(nfI,ibe) = sqrt( (xxmIH(nfI,ibe,1) - xxmI1(nfI,ibe,1))**2 &
                                               + (xxmIH(nfI,ibe,2) - xxmI1(nfI,ibe,2))**2 &
                                               + (xxmIH(nfI,ibe,3) - xxmI1(nfI,ibe,3))**2 )
                           dsI2(nfI,ibe) = sqrt( (xxmI1(nfI,ibe,1) - xxmI2(nfI,ibe,1))**2 &
                                               + (xxmI1(nfI,ibe,2) - xxmI2(nfI,ibe,2))**2 &
                                               + (xxmI1(nfI,ibe,3) - xxmI2(nfI,ibe,3))**2 )
                         enddo

                         do nfO=nfOBH,nfOEH
                           dsO1(nfO,ibe) = sqrt( (xxmOH(nfO,ibe,1) - xxmO1(nfO,ibe,1))**2 &
                                               + (xxmOH(nfO,ibe,2) - xxmO1(nfO,ibe,2))**2 &
                                               + (xxmOH(nfO,ibe,3) - xxmO1(nfO,ibe,3))**2 )
                           dsO2(nfO,ibe) = sqrt( (xxmO1(nfO,ibe,1) - xxmO2(nfO,ibe,1))**2 &
                                               + (xxmO1(nfO,ibe,2) - xxmO2(nfO,ibe,2))**2 &
                                               + (xxmO1(nfO,ibe,3) - xxmO2(nfO,ibe,3))**2 )
                         enddo
                       enddo

                       deallocate(xxmIH, xxmI1, xxmI2, xxmOH, xxmO1, xxmO2, stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory deallocation failure for xxmI and xxmO")

                     case (kMin, kMax)
                       nfOB  = jjI
                       nfOBH = jjI - 1
                       nfOE  = jjF
                       nfOEH = jjF + 1
                       nfIB  = iiI
                       nfIBH = iiI - 1
                       nfIE  = iiF
                       nfIEH = iiF + 1

                       kkH = kkI

                       allocate(xxmIH(nfIBH:nfIEH,2,3), xxmI1(nfIBH:nfIEH,2,3), &
                                xxmI2(nfIBH:nfIEH,2,3), dsI1(nfIBH:nfIEH,2),    &
                                dsI2(nfIBH:nfIEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmI and dsI")

                       allocate(xxmOH(nfOBH:nfOEH,2,3), xxmO1(nfOBH:nfOEH,2,3), &
                                xxmO2(nfOBH:nfOEH,2,3), dsO1(nfOBH:nfOEH,2),    &
                                dsO2(nfOBH:nfOEH,2), stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory allocation failure for xxmO and dsO")

                       do ndir=1,3
                         do nfI=nfIBH,nfIEH
                           xxmIH(nfI,1,ndir) = eighth*( x(nfI  ,nfOBH-1,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH  ,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH-1,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH  ,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH-1,kkH-1,ndir) &
                                                      + x(nfI  ,nfOBH  ,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH-1,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH  ,kkH-1,ndir) )
                           xxmIH(nfI,2,ndir) = eighth*( x(nfI  ,nfOEH-1,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH  ,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH-1,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH  ,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH-1,kkH-1,ndir) &
                                                      + x(nfI  ,nfOEH  ,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH-1,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH  ,kkH-1,ndir) )

                           xxmI1(nfI,1,ndir) = eighth*( x(nfI  ,nfOBH  ,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH+1,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH  ,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH+1,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH  ,kkH-1,ndir) &
                                                      + x(nfI  ,nfOBH+1,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH  ,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH+1,kkH-1,ndir) )
                           xxmI1(nfI,2,ndir) = eighth*( x(nfI  ,nfOEH-2,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH-1,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH-2,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH-1,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH-2,kkH-1,ndir) &
                                                      + x(nfI  ,nfOEH-1,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH-2,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH-1,kkH-1,ndir) )

                           xxmI2(nfI,1,ndir) = eighth*( x(nfI  ,nfOBH+1,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH+2,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH+1,kkH  ,ndir) &
                                                      + x(nfI-1,nfOBH+2,kkH  ,ndir) &
                                                      + x(nfI  ,nfOBH+1,kkH-1,ndir) &
                                                      + x(nfI  ,nfOBH+2,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH+1,kkH-1,ndir) &
                                                      + x(nfI-1,nfOBH+2,kkH-1,ndir) )
                           xxmI2(nfI,2,ndir) = eighth*( x(nfI  ,nfOEH-3,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH-2,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH-3,kkH  ,ndir) &
                                                      + x(nfI-1,nfOEH-2,kkH  ,ndir) &
                                                      + x(nfI  ,nfOEH-3,kkH-1,ndir) &
                                                      + x(nfI  ,nfOEH-2,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH-3,kkH-1,ndir) &
                                                      + x(nfI-1,nfOEH-2,kkH-1,ndir) )
                         enddo

                         do nfO=nfOBH,nfOEH
                           xxmOH(nfO,1,ndir) = eighth*( x(nfIBH-1,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH  ,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH-1,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH  ,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH-1,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH  ,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH-1,nfO-1,kkH-1,ndir) &
                                                      + x(nfIBH  ,nfO-1,kkH-1,ndir) )
                           xxmOH(nfO,2,ndir) = eighth*( x(nfIEH-1,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH  ,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH-1,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH  ,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH-1,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH  ,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH-1,nfO-1,kkH-1,ndir) &
                                                      + x(nfIEH  ,nfO-1,kkH-1,ndir) )

                           xxmO1(nfO,1,ndir) = eighth*( x(nfIBH  ,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH+1,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH  ,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH+1,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH  ,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH+1,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH  ,nfO-1,kkH-1,ndir) &
                                                      + x(nfIBH+1,nfO-1,kkH-1,ndir) )
                           xxmO1(nfO,2,ndir) = eighth*( x(nfIEH-2,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH-1,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH-2,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH-1,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH-2,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH-1,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH-2,nfO-1,kkH-1,ndir) &
                                                      + x(nfIEH-1,nfO-1,kkH-1,ndir) )

                           xxmO2(nfO,1,ndir) = eighth*( x(nfIBH+1,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH+2,nfO  ,kkH  ,ndir) &
                                                      + x(nfIBH+1,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH+2,nfO-1,kkH  ,ndir) &
                                                      + x(nfIBH+1,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH+2,nfO  ,kkH-1,ndir) &
                                                      + x(nfIBH+1,nfO-1,kkH-1,ndir) &
                                                      + x(nfIBH+2,nfO-1,kkH-1,ndir) )
                           xxmO2(nfO,2,ndir) = eighth*( x(nfIEH-3,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH-2,nfO  ,kkH  ,ndir) &
                                                      + x(nfIEH-3,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH-2,nfO-1,kkH  ,ndir) &
                                                      + x(nfIEH-3,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH-2,nfO  ,kkH-1,ndir) &
                                                      + x(nfIEH-3,nfO-1,kkH-1,ndir) &
                                                      + x(nfIEH-2,nfO-1,kkH-1,ndir) )
                         enddo
                       enddo

                       do ibe=1,2
                         do nfI=nfIBH,nfIEH
                           dsI1(nfI,ibe) = sqrt( (xxmIH(nfI,ibe,1) - xxmI1(nfI,ibe,1))**2 &
                                               + (xxmIH(nfI,ibe,2) - xxmI1(nfI,ibe,2))**2 &
                                               + (xxmIH(nfI,ibe,3) - xxmI1(nfI,ibe,3))**2 )
                           dsI2(nfI,ibe) = sqrt( (xxmI1(nfI,ibe,1) - xxmI2(nfI,ibe,1))**2 &
                                               + (xxmI1(nfI,ibe,2) - xxmI2(nfI,ibe,2))**2 &
                                               + (xxmI1(nfI,ibe,3) - xxmI2(nfI,ibe,3))**2 )
                         enddo

                         do nfO=nfOBH,nfOEH
                           dsO1(nfO,ibe) = sqrt( (xxmOH(nfO,ibe,1) - xxmO1(nfO,ibe,1))**2 &
                                               + (xxmOH(nfO,ibe,2) - xxmO1(nfO,ibe,2))**2 &
                                               + (xxmOH(nfO,ibe,3) - xxmO1(nfO,ibe,3))**2 )
                           dsO2(nfO,ibe) = sqrt( (xxmO1(nfO,ibe,1) - xxmO2(nfO,ibe,1))**2 &
                                               + (xxmO1(nfO,ibe,2) - xxmO2(nfO,ibe,2))**2 &
                                               + (xxmO1(nfO,ibe,3) - xxmO2(nfO,ibe,3))**2 )
                         enddo
                       enddo

                       deallocate(xxmIH, xxmI1, xxmI2, xxmOH, xxmO1, xxmO2, stat=ierr)
                       if(ierr /= 0)                         &
                         call terminate("sumb_setPointData", &
                                        "Memory deallocation failure for xxmI and xxmO")

                   end select

                   do ibe=1,2
                     if(ibe == 1) then
                       nfIH = nfIBH
                       nfI1 = nfIB
                       nfI2 = nfIB + 1

                       nfOH = nfOBH
                       nfO1 = nfOB
                       nfO2 = nfOB + 1
                     else if(ibe == 2) then
                       nfIH = nfIEH
                       nfI1 = nfIE
                       nfI2 = nfIE - 1

                       nfOH = nfOEH
                       nfO1 = nfOE
                       nfO2 = nfOE - 1
                     endif

                     do nfI=nfIB,nfIE
                       factor1 =  (dsI1(nfI,ibe) + dsI2(nfI,ibe))/dsI2(nfI,ibe)
                       factor2 = - dsI1(nfI,ibe)/dsI2(nfI,ibe)

                       if(associated(BCData(nBFs)%rho)) then
                         BCData(nBFs)%rho(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%rho(nfI,nfO1) &
                         + factor2*BCData(nBFs)%rho(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velX)) then
                         BCData(nBFs)%velX(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velX(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velX(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velY)) then
                         BCData(nBFs)%velY(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velY(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velY(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velZ)) then
                         BCData(nBFs)%velZ(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velZ(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velZ(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowXDirInlet)) then
                         BCData(nBFs)%flowXDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowXDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowXDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowYDirInlet)) then
                         BCData(nBFs)%flowYDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowYDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowYDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowZDirInlet)) then
                         BCData(nBFs)%flowZDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowZDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowZDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ps)) then
                         BCData(nBFs)%ps(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ps(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ps(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ptInlet)) then
                         BCData(nBFs)%ptInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ptInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ptInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ttInlet)) then
                         BCData(nBFs)%ttInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ttInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ttInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%htInlet)) then
                         BCData(nBFs)%htInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%htInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%htInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%turbInlet)) then
                         BCData(nBFs)%turbInlet(nfI,nfOH,:) =         &
                           factor1*BCData(nBFs)%turbInlet(nfI,nfO1,:) &
                         + factor2*BCData(nBFs)%turbInlet(nfI,nfO2,:)
                       endif
                     enddo

                     do nfO=nfOB,nfOE
                       factor1 =  (dsO1(nfO,ibe) + dsO2(nfO,ibe))/dsO2(nfO,ibe)
                       factor2 = - dsO1(nfO,ibe)/dsO2(nfO,ibe)

                       if(associated(BCData(nBFs)%rho)) then
                         BCData(nBFs)%rho(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%rho(nfI1,nfO) &
                         + factor2*BCData(nBFs)%rho(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%velX)) then
                         BCData(nBFs)%velX(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%velX(nfI1,nfO) &
                         + factor2*BCData(nBFs)%velX(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%velY)) then
                         BCData(nBFs)%velY(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%velY(nfI1,nfO) &
                         + factor2*BCData(nBFs)%velY(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%velZ)) then
                         BCData(nBFs)%velZ(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%velZ(nfI1,nfO) &
                         + factor2*BCData(nBFs)%velZ(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%flowXDirInlet)) then
                         BCData(nBFs)%flowXDirInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%flowXDirInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%flowXDirInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%flowYDirInlet)) then
                         BCData(nBFs)%flowYDirInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%flowYDirInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%flowYDirInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%flowZDirInlet)) then
                         BCData(nBFs)%flowZDirInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%flowZDirInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%flowZDirInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%ps)) then
                         BCData(nBFs)%ps(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%ps(nfI1,nfO) &
                         + factor2*BCData(nBFs)%ps(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%ptInlet)) then
                         BCData(nBFs)%ptInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%ptInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%ptInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%ttInlet)) then
                         BCData(nBFs)%ttInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%ttInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%ttInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%htInlet)) then
                         BCData(nBFs)%htInlet(nfIH,nfO) =         &
                           factor1*BCData(nBFs)%htInlet(nfI1,nfO) &
                         + factor2*BCData(nBFs)%htInlet(nfI2,nfO)
                       endif

                       if(associated(BCData(nBFs)%turbInlet)) then
                         BCData(nBFs)%turbInlet(nfIH,nfO,:) =         &
                           factor1*BCData(nBFs)%turbInlet(nfI1,nfO,:) &
                         + factor2*BCData(nBFs)%turbInlet(nfI2,nfO,:)
                       endif
                     enddo
                   enddo

                   ! Corner points.

                   do ibe=1,2
                     if(ibe == 1) then
                       nfOH = nfOBH
                       nfO1 = nfOB
                       nfO2 = nfOB + 1
                     else if(ibe == 2) then
                       nfOH = nfOEH
                       nfO1 = nfOE
                       nfO2 = nfOE - 1
                     endif

                     do nfI=nfIBH,nfIEH,nfIEH-nfIBH
                       factor1 =  (dsI1(nfI,ibe) + dsI2(nfI,ibe))/dsI2(nfI,ibe)
                       factor2 = - dsI1(nfI,ibe)/dsI2(nfI,ibe)

                       if(associated(BCData(nBFs)%rho)) then
                         BCData(nBFs)%rho(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%rho(nfI,nfO1) &
                         + factor2*BCData(nBFs)%rho(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velX)) then
                         BCData(nBFs)%velX(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velX(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velX(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velY)) then
                         BCData(nBFs)%velY(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velY(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velY(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%velZ)) then
                         BCData(nBFs)%velZ(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%velZ(nfI,nfO1) &
                         + factor2*BCData(nBFs)%velZ(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowXDirInlet)) then
                         BCData(nBFs)%flowXDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowXDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowXDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowYDirInlet)) then
                         BCData(nBFs)%flowYDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowYDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowYDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%flowZDirInlet)) then
                         BCData(nBFs)%flowZDirInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%flowZDirInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%flowZDirInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ps)) then
                         BCData(nBFs)%ps(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ps(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ps(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ptInlet)) then
                         BCData(nBFs)%ptInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ptInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ptInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%ttInlet)) then
                         BCData(nBFs)%ttInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%ttInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%ttInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%htInlet)) then
                         BCData(nBFs)%htInlet(nfI,nfOH) =         &
                           factor1*BCData(nBFs)%htInlet(nfI,nfO1) &
                         + factor2*BCData(nBFs)%htInlet(nfI,nfO2)
                       endif

                       if(associated(BCData(nBFs)%turbInlet)) then
                         BCData(nBFs)%turbInlet(nfI,nfOH,:) =         &
                           factor1*BCData(nBFs)%turbInlet(nfI,nfO1,:) &
                         + factor2*BCData(nBFs)%turbInlet(nfI,nfO2,:)
                       endif
                     enddo
                   enddo

                   deallocate(dsI1, dsI2, dsO1, dsO2, stat=ierr)
                   if(ierr /= 0)                         &
                     call terminate("sumb_setPointData", &
                                    "Memory deallocation failure for ds")

                 endif
               endif
             endif
           enddo bocosLoop
         enddo domainLoop
       enddo levelLoop

       currentLevel = mgStartlevel
       groundLevel  = mgStartlevel

       call applyAllBC(.true.)

       ! Obtain the coarse-level solutions at the interfaces
       ! by interpolation from the fine-level one, if they
       ! are not provided from the coupler. This is simply
       ! a copy of the procedure performed in the subroutine
       ! setBCDataCoarseGrid.

       if(.not. cplGetCoarseSol) then

         nLevels = ubound(flowDoms,2)

         ! Loop over the coarser grid levels. It is assumed that the
         ! bc data of the finest level is set correctly.

         intpolLevelLoop: do level=2,nLevels

           ! Store the fine grid level a bit easier.

           levm1 = level - 1

           ! Loop over the number of local blocks.

           intpolDomainLoop: do i=1,nDom

             ! Set the pointers to the coarse block.
             ! About the time-spectral method, only the zero mode is considered.

             call setPointers(i, level, 1)

             ! Loop over the boundary subfaces and interpolate the
             ! prescribed boundary data for this grid level.

               intpolBocosLoop: do j=1,nBocos

               ! Determine the block face on which the subface is
               ! located and set some multigrid variables accordingly.

               select case (BCFaceID(j))
                 case (iMin,iMax)
                   iiMax = jl; jjMax = kl
                   iFine => mgJFine; jFine => mgKFine

                 case (jMin,jMax)
                   iiMax = il; jjMax = kl
                   iFine => mgIFine; jFine => mgKFine

                 case (kMin,kMax)
                   iiMax = il; jjMax = jl
                   iFine => mgIFine; jFine => mgJFine

               end select

               ! Abbreviate the size of the subface a bit easier.

               iBeg = BCData(j)%icBeg; iEnd = BCData(j)%icEnd
               jBeg = BCData(j)%jcBeg; jEnd = BCData(j)%jcEnd

               if(equations == RANSEquations)                    &
                 call interpolateBCVecData(BCData(j)%turbInlet,  &
                       flowDoms(i,levm1,1)%BCData(j)%turbInlet,  &
                       nt1, nt2)

               if(BCType(j) == DomainInterfaceAll) then
                 call interpolateBCData(BCData(j)%rho,           &
                    flowDoms(i,levm1,1)%BCData(j)%rho)
                 call interpolateBCData(BCData(j)%velX,          &
                    flowDoms(i,levm1,1)%BCData(j)%velX)
                 call interpolateBCData(BCData(j)%velY,          &
                    flowDoms(i,levm1,1)%BCData(j)%velY)
                 call interpolateBCData(BCData(j)%velZ,          &
                    flowDoms(i,levm1,1)%BCData(j)%velZ)
                 call interpolateBCData(BCData(j)%ps,            &
                    flowDoms(i,levm1,1)%BCData(j)%ps)

               else if(BCType(j) == DomainInterfaceRhoUVW) then
                 call interpolateBCData(BCData(j)%rho,           &
                    flowDoms(i,levm1,1)%BCData(j)%rho)
                 call interpolateBCData(BCData(j)%velX,          &
                    flowDoms(i,levm1,1)%BCData(j)%velX)
                 call interpolateBCData(BCData(j)%velY,          &
                    flowDoms(i,levm1,1)%BCData(j)%velY)
                 call interpolateBCData(BCData(j)%velZ,          &
                    flowDoms(i,levm1,1)%BCData(j)%velZ)

               else if(BCType(j) == DomainInterfaceP) then
                 call interpolateBCData(BCData(j)%ps,            &
                    flowDoms(i,levm1,1)%BCData(j)%ps)

               else if(BCType(j) == DomainInterfaceTotal) then
                 call interpolateBCData(BCData(j)%ptInlet,       &
                    flowDoms(i,levm1,1)%BCData(j)%ptInlet)
                 call interpolateBCData(BCData(j)%ttInlet,       &
                    flowDoms(i,levm1,1)%BCData(j)%ttInlet)
                 call interpolateBCData(BCData(j)%flowXDirInlet, &
                    flowDoms(i,levm1,1)%BCData(j)%flowXDirInlet)
                 call interpolateBCData(BCData(j)%flowYDirInlet, &
                    flowDoms(i,levm1,1)%BCData(j)%flowYDirInlet)
                 call interpolateBCData(BCData(j)%flowZDirInlet, &
                    flowDoms(i,levm1,1)%BCData(j)%flowZDirInlet)

                 ! Some additional variables should be
                 ! computed/corrected for total conditions.

                 do lco=jBeg,jEnd
                   do kco=iBeg,iEnd

                     ! Compute the total enthalpy.

                     call computeHtot(BCData(j)%ttInlet(kco,lco), &
                                      BCData(j)%htInlet(kco,lco))

                     ! Flow direction.

                     dir(1) = BCData(j)%flowXdirInlet(kco,lco)
                     dir(2) = BCData(j)%flowYdirInlet(kco,lco)
                     dir(3) = BCData(j)%flowZdirInlet(kco,lco)

                     var = one/max(eps,sqrt(dir(1)**2 + dir(2)**2 &
                         +                  dir(3)**2))

                     BCData(j)%flowXdirInlet(kco,lco) = var*dir(1)
                     BCData(j)%flowYdirInlet(kco,lco) = var*dir(2)
                     BCData(j)%flowZdirInlet(kco,lco) = var*dir(3)

                   enddo
                 enddo
               endif

             enddo intpolBocosLoop
           enddo intpolDomainLoop
         enddo intpolLevelLoop
       endif

       !=================================================================

       contains

         !===============================================================

         subroutine interpolateBCData(varCoarse, varFine)
!
!        ****************************************************************
!        *                                                              *
!        * interpolateBCData interpolates the given data array from     *
!        * the fine to the coarse grid. Of course only if the fine      *
!        * array is associated with some data.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(:,:), pointer :: varCoarse
         real(kind=realType), dimension(:,:), pointer :: varFine
!
!        Local variables.
!
         integer(kind=intType) :: i, j, if1, if2, jf1, jf2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check if varFine is associated to data. If not return.

         if(.not. associated(varFine)) return

         ! Loop over the faces of the given subface.
         ! First the j-direction.

         do j=jBeg,jEnd

           ! Determine the two children in this direction. Take care of
           ! the halo's, as this info is only available for owned cells.

           if(j < 2) then
             jf1 = 1; jf2 = 1
           else if(j > jjMax) then
             jf1 = jFine(jjMax,2) +1; jf2 = jf1
           else
             jf1 = jFine(j,1); jf2 = jFine(j,2)
           endif

           ! Loop in the i-direction.

           do i=iBeg,iEnd

             ! Determine the two children in this direction.
             ! Same story as in j-direction.

             if(i < 2) then
               if1 = 1; if2 = 1
             else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
             else
               if1 = iFine(i,1); if2 = iFine(i,2)
             endif

             ! Compute the coarse grid data as the average of the
             ! 4 fine grid values.

             varCoarse(i,j) = fourth*(varFine(if1,jf1) &
                             +        varFine(if2,jf1) &
                             +        varFine(if1,jf2) &
                             +        varFine(if2,jf2))
           enddo
         enddo

         end subroutine interpolateBCData

         !===============================================================

         subroutine interpolateBCVecData(varCoarse, varFine, &
                                         nStart, nEnd)
!
!        ****************************************************************
!        *                                                              *
!        * interpolateBCVecData interpolates the given data array       *
!        * from the fine to the coarse grid. Of course only if the fine *
!        * array is associated with some data.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nStart, nEnd

         real(kind=realType), dimension(:,:,:), pointer :: varCoarse
         real(kind=realType), dimension(:,:,:), pointer :: varFine
!
!        Local variables.
!
         integer(kind=intType) :: nn, i, j, if1, if2, jf1, jf2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check if varFine is associated to data. if not return.

         if(.not. associated(varFine)) return

         ! Loop over the faces of the given subface.
         ! First the j-direction.

         do j=jBeg,jEnd

           ! Determine the two children in this direction. Take care of
           ! the halo's, as this info is only available for owned cells.

           if(j < 2) then
             jf1 = 1; jf2 = 1
           else if(j > jjMax) then
             jf1 = jFine(jjMax,2) +1; jf2 = jf1
           else
             jf1 = jFine(j,1); jf2 = jFine(j,2)
           endif

           ! Loop in the i-direction.

           do i=iBeg,iEnd

             ! Determine the two children in this direction.
             ! Same story as in j-direction.

             if(i < 2) then
               if1 = 1; if2 = 1
             else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
             else
               if1 = iFine(i,1); if2 = iFine(i,2)
             endif

             ! Compute the coarse grid data as the average of the
             ! 4 fine grid values.

             do nn=nStart,nEnd
               varCoarse(i,j,nn) = fourth*(varFine(if1,jf1,nn) &
                                  +        varFine(if2,jf1,nn) &
                                  +        varFine(if1,jf2,nn) &
                                  +        varFine(if2,jf2,nn))
             enddo
           enddo
         enddo

         end subroutine interpolateBCVecData

       end subroutine sumb_setPointData

!      ==================================================================

       subroutine sumb_getMeshSize(nNodes, nTetra, nPyra, nPrism, &
                                   nHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshSize provides the coupler with the number of       *
!      * vertices in the local domain. Note that duplicate              *
!      * interpolations are done when local domains are passed based on *
!      * the vertex positions. However, our agreement is to pass        *
!      * the coordinates and flow variables at vertices to the coupler. *
!      * Also note that only the value of nHexa is nonzero for SUmb.    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: nNodes, nTetra, nPyra, &
                                             nPrism, nHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: nn
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNodes = 0
       nTetra = 0
       nPyra  = 0
       nPrism = 0
       nHexa  = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           ! Passing the number of cell vertices to the coupler.

           nHexa  = nHexa  + nx*ny*nz
           nNodes = nNodes + (nx + 1)*(ny + 1)*(nz + 1)
         enddo domainLoop
       endif

       nTetraTrue = nTetra
       nPyraTrue  = nPyra
       nPrismTrue = nPrism
       nHexaTrue  = nHexa
       nNodesTrue = nNodes

       nTetraAlloc = max(1,nTetra)
       nPyraAlloc  = max(1,nPyra)
       nPrismAlloc = max(1,nPrism)
       nHexaAlloc  = max(1,nHexa)
       nNodesAlloc = max(1,nNodes)

       end subroutine sumb_getMeshSize

!      ==================================================================

       subroutine sumb_getMeshSizeC(nNodes, nTetra, nPyra, nPrism, &
                                    nHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshSizeC provides the coupler with the number of      *
!      * cell centers in the local domain. This is the counterpart of   *
!      * sumb_getMeshSize for cell centers.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: nNodes, nTetra, nPyra, &
                                             nPrism, nHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: nn
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNodes = 0
       nTetra = 0
       nPyra  = 0
       nPrism = 0
       nHexa  = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           ! Passing the number of cell centers to the coupler.

           nHexa  = nHexa  + (nx - 1)*(ny - 1)*(nz - 1)
           nNodes = nNodes + nx*ny*nz
         enddo domainLoop
       endif

       nTetraTrue = nTetra
       nPyraTrue  = nPyra
       nPrismTrue = nPrism
       nHexaTrue  = nHexa
       nNodesTrue = nNodes

       nTetraAlloc = max(1,nTetra)
       nPyraAlloc  = max(1,nPyra)
       nPrismAlloc = max(1,nPrism)
       nHexaAlloc  = max(1,nHexa)
       nNodesAlloc = max(1,nNodes)

       end subroutine sumb_getMeshSizeC

!      ==================================================================

       subroutine sumb_getMeshSizeCH(nNodes, nTetra, nPyra, nPrism, &
                                     nHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshSizeCH provides the coupler with the number of     *
!      * cell centers, including the 1st-level halos, in the local      *
!      * domain. This is the counterpart of sumb_getMeshSize for        *
!      * cell centers + 1st-level halos.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: nNodes, nTetra, nPyra, &
                                             nPrism, nHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: nn
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNodes = 0
       nTetra = 0
       nPyra  = 0
       nPrism = 0
       nHexa  = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           ! Passing the number of cell centers to the coupler.

           nHexa  = nHexa  + (nx + 1)*(ny + 1)*(nz + 1)
           nNodes = nNodes + (nx + 2)*(ny + 2)*(nz + 2)
         enddo domainLoop
       endif

       nTetraTrue = nTetra
       nPyraTrue  = nPyra
       nPrismTrue = nPrism
       nHexaTrue  = nHexa
       nNodesTrue = nNodes

       nTetraAlloc = max(1,nTetra)
       nPyraAlloc  = max(1,nPyra)
       nPrismAlloc = max(1,nPrism)
       nHexaAlloc  = max(1,nHexa)
       nNodesAlloc = max(1,nNodes)

       end subroutine sumb_getMeshSizeCH

!      ==================================================================

       subroutine sumb_getMeshGeom(xyz, connTetra, connPyra, connPrism, &
                                   connHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshGeom provides the coupler with coordinates of all  *
!      * node points and connectivity of cells.                         *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: xyz
       integer(kind=intType), dimension(:,:), intent(out) :: connTetra
       integer(kind=intType), dimension(:,:), intent(out) :: connPyra
       integer(kind=intType), dimension(:,:), intent(out) :: connPrism
       integer(kind=intType), dimension(:,:), intent(out) :: connHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nNod, nCel, nn, i, j, k, l, m, n
       integer(kind=intType), allocatable, dimension(:,:,:) :: indMap
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       xyz = 0.
       connTetra = 0
       connPyra  = 0
       connPrism = 0
       connHexa  = 0

       nNod = 0
       nCel = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           allocate(indMap(il,jl,kl), stat=ierr)
           if(ierr /= 0)                        &
             call terminate("sumb_getMeshGeom", &
                            "Memory allocation failure for indMap")

           ! Passing the cell vertices to the coupler.

           do k=1,kl
             do j=1,jl
               do i=1,il
                 nNod = nNod + 1

                 xyz(1,nNod) = x(i,j,k,1)
                 xyz(2,nNod) = x(i,j,k,2)
                 xyz(3,nNod) = x(i,j,k,3)
                 indMap(i,j,k) = nNod
               enddo
             enddo
           enddo

           ! Creating the hexahedral connectivity

           do k=1,nz
             n = k + 1
             do j=1,ny
               m = j + 1
               do i=1,nx
                 l = i + 1

                 nCel = nCel + 1
                 connHexa(1,nCel) = indMap(i,j,k)
                 connHexa(2,nCel) = indMap(l,j,k)
                 connHexa(3,nCel) = indMap(l,m,k)
                 connHexa(4,nCel) = indMap(i,m,k)
                 connHexa(5,nCel) = indMap(i,j,n)
                 connHexa(6,nCel) = indMap(l,j,n)
                 connHexa(7,nCel) = indMap(l,m,n)
                 connHexa(8,nCel) = indMap(i,m,n)
               enddo
             enddo
           enddo

           deallocate(indMap, stat=ierr)
           if(ierr /= 0)                        &
             call terminate("sumb_getMeshGeom", &
                            "Memory deallocation failure for indMap")
         enddo domainLoop
       endif

       end subroutine sumb_getMeshGeom

!      ==================================================================

       subroutine sumb_getMeshGeomC(xyz, connTetra, connPyra, connPrism, &
                                    connHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshGeomC provides the coupler with coordinates of all *
!      * cell-center points and connectivity of cells. This is the      *
!      * counterpart of sumb_getMeshGeom for cell centers.              *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use constants
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: xyz
       integer(kind=intType), dimension(:,:), intent(out) :: connTetra
       integer(kind=intType), dimension(:,:), intent(out) :: connPyra
       integer(kind=intType), dimension(:,:), intent(out) :: connPrism
       integer(kind=intType), dimension(:,:), intent(out) :: connHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nCel1, nCel2, nn, i, j, k, l, m, n
       integer(kind=intType), allocatable, dimension(:,:,:) :: indMap
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       xyz = 0.
       connTetra = 0
       connPyra  = 0
       connPrism = 0
       connHexa  = 0

       nCel1 = 0
       nCel2 = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           allocate(indMap(nx,ny,nz), stat=ierr)
           if(ierr /= 0)                         &
             call terminate("sumb_getMeshGeomC", &
                            "Memory allocation failure for indMap")

           ! Passing the cell centers to the coupler.

           do k=1,nz
             do j=1,ny
               do i=1,nx
                 nCel1 = nCel1 + 1

                 xyz(1,nCel1) = eighth*( x(i  ,j  ,k  ,1) &
                                       + x(i  ,j  ,k+1,1) &
                                       + x(i  ,j+1,k  ,1) &
                                       + x(i  ,j+1,k+1,1) &
                                       + x(i+1,j  ,k  ,1) &
                                       + x(i+1,j  ,k+1,1) &
                                       + x(i+1,j+1,k  ,1) &
                                       + x(i+1,j+1,k+1,1) )
                 xyz(2,nCel1) = eighth*( x(i  ,j  ,k  ,2) &
                                       + x(i  ,j  ,k+1,2) &
                                       + x(i  ,j+1,k  ,2) &
                                       + x(i  ,j+1,k+1,2) &
                                       + x(i+1,j  ,k  ,2) &
                                       + x(i+1,j  ,k+1,2) &
                                       + x(i+1,j+1,k  ,2) &
                                       + x(i+1,j+1,k+1,2) )
                 xyz(3,nCel1) = eighth*( x(i  ,j  ,k  ,3) &
                                       + x(i  ,j  ,k+1,3) &
                                       + x(i  ,j+1,k  ,3) &
                                       + x(i  ,j+1,k+1,3) &
                                       + x(i+1,j  ,k  ,3) &
                                       + x(i+1,j  ,k+1,3) &
                                       + x(i+1,j+1,k  ,3) &
                                       + x(i+1,j+1,k+1,3) )
                 indMap(i,j,k) = nCel1
               enddo
             enddo
           enddo

           ! Creating the hexahedral connectivity

           do k=1,(nz - 1)
             n = k + 1
             do j=1,(ny - 1)
               m = j + 1
               do i=1,(nx - 1)
                 l = i + 1

                 nCel2 = nCel2 + 1
                 connHexa(1,nCel2) = indMap(i,j,k)
                 connHexa(2,nCel2) = indMap(l,j,k)
                 connHexa(3,nCel2) = indMap(l,m,k)
                 connHexa(4,nCel2) = indMap(i,m,k)
                 connHexa(5,nCel2) = indMap(i,j,n)
                 connHexa(6,nCel2) = indMap(l,j,n)
                 connHexa(7,nCel2) = indMap(l,m,n)
                 connHexa(8,nCel2) = indMap(i,m,n)
               enddo
             enddo
           enddo

           deallocate(indMap, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("sumb_getMeshGeomC", &
                            "Memory deallocation failure for indMap")
         enddo domainLoop
       endif

       end subroutine sumb_getMeshGeomC

!      ==================================================================

       subroutine sumb_getMeshGeomCH(xyz, connTetra, connPyra, connPrism, &
                                     connHexa, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshGeomCH provides the coupler with coordinates of    *
!      * all cell-center and 1st-level halo points and connectivity of  *
!      * cells. This is the counterpart of sumb_getMeshGeom for cell    *
!      * centers + 1st-level halos.                                     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use constants
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: xyz
       integer(kind=intType), dimension(:,:), intent(out) :: connTetra
       integer(kind=intType), dimension(:,:), intent(out) :: connPyra
       integer(kind=intType), dimension(:,:), intent(out) :: connPrism
       integer(kind=intType), dimension(:,:), intent(out) :: connHexa
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nCel1, nCel2, nn, i, j, k, l, m, n
       integer(kind=intType), allocatable, dimension(:,:,:) :: indMap
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       xyz = 0.
       connTetra = 0
       connPyra  = 0
       connPrism = 0
       connHexa  = 0

       nCel1 = 0
       nCel2 = 0

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName == trimGeomName) then
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block.
           ! Only the finest level is considered.

           call setPointers(nn, 1, 1)

           allocate(indMap(nx+2,ny+2,nz+2), stat=ierr)
           if(ierr /= 0)                          &
             call terminate("sumb_getMeshGeomCH", &
                            "Memory allocation failure for indMap")

           ! Passing the cell centers, including the 1st-level halos,
           ! to the coupler.

           do k=1,nz+2
             do j=1,ny+2
               do i=1,nx+2
                 nCel1 = nCel1 + 1

                 xyz(1,nCel1) = eighth*( x(i-1,j-1,k-1,1) &
                                       + x(i-1,j-1,k  ,1) &
                                       + x(i-1,j  ,k-1,1) &
                                       + x(i-1,j  ,k  ,1) &
                                       + x(i  ,j-1,k-1,1) &
                                       + x(i  ,j-1,k  ,1) &
                                       + x(i  ,j  ,k-1,1) &
                                       + x(i  ,j  ,k  ,1) )
                 xyz(2,nCel1) = eighth*( x(i-1,j-1,k-1,2) &
                                       + x(i-1,j-1,k  ,2) &
                                       + x(i-1,j  ,k-1,2) &
                                       + x(i-1,j  ,k  ,2) &
                                       + x(i  ,j-1,k-1,2) &
                                       + x(i  ,j-1,k  ,2) &
                                       + x(i  ,j  ,k-1,2) &
                                       + x(i  ,j  ,k  ,2) )
                 xyz(3,nCel1) = eighth*( x(i-1,j-1,k-1,3) &
                                       + x(i-1,j-1,k  ,3) &
                                       + x(i-1,j  ,k-1,3) &
                                       + x(i-1,j  ,k  ,3) &
                                       + x(i  ,j-1,k-1,3) &
                                       + x(i  ,j-1,k  ,3) &
                                       + x(i  ,j  ,k-1,3) &
                                       + x(i  ,j  ,k  ,3) )
                 indMap(i,j,k) = nCel1
               enddo
             enddo
           enddo

           ! Creating the hexahedral connectivity

           do k=1,(nz + 1)
             n = k + 1
             do j=1,(ny + 1)
               m = j + 1
               do i=1,(nx + 1)
                 l = i + 1

                 nCel2 = nCel2 + 1
                 connHexa(1,nCel2) = indMap(i,j,k)
                 connHexa(2,nCel2) = indMap(l,j,k)
                 connHexa(3,nCel2) = indMap(l,m,k)
                 connHexa(4,nCel2) = indMap(i,m,k)
                 connHexa(5,nCel2) = indMap(i,j,n)
                 connHexa(6,nCel2) = indMap(l,j,n)
                 connHexa(7,nCel2) = indMap(l,m,n)
                 connHexa(8,nCel2) = indMap(i,m,n)
               enddo
             enddo
           enddo

           deallocate(indMap, stat=ierr)
           if(ierr /= 0)                          &
             call terminate("sumb_getMeshGeomCH", &
                            "Memory deallocation failure for indMap")
         enddo domainLoop
       endif

       end subroutine sumb_getMeshGeomCH

!      ==================================================================

       subroutine sumb_getMeshData(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshData provides the coupler with the requested flow  *
!      * fields.                                                        *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: kk, ll, nNod, nn, ic, jc, kc, &
                                ip, jp, kp, mm, ietcVar, iwReq
       integer(kind=intType), dimension(nDataNames) :: iwMapReq
     ! character(len=3) :: prefix, prefixTurb
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
       real(kind=realType) :: fac1, fac2, fac3, fac4, &
                              fac5, fac6, fac7, fac8
       real(kind=realType) :: volC
       real(kind=realType) :: sumbScal1, sumbScal2, sumbScal3, &
                              sumbScal4, sumbScal5, sumbScal6, &
                              sumbScal7, sumbScal8, sumbVecMag
       real(kind=realType) :: sumbVecX1, sumbVecX2, sumbVecX3, &
                              sumbVecX4, sumbVecX5, sumbVecX6, &
                              sumbVecX7, sumbVecX8, sumbVecXC
       real(kind=realType) :: sumbVecY1, sumbVecY2, sumbVecY3, &
                              sumbVecY4, sumbVecY5, sumbVecY6, &
                              sumbVecY7, sumbVecY8, sumbVecYC
       real(kind=realType) :: sumbVecZ1, sumbVecZ2, sumbVecZ3, &
                              sumbVecZ4, sumbVecZ5, sumbVecZ6, &
                              sumbVecZ7, sumbVecZ8, sumbVecZC
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Prefix to tell the averaged variables from instantaneous ones.

     ! prefix = 'Avg'
     ! prefixTurb = 'Avg'

       ! First, define a list of flow variables which this particular
       ! instance of SUmb can provide to the coupler.
       ! More variables can be added here.

       if(.not. allocated(dataNamesSUmb)) then
         nDataSUmb = nw + 6
         allocate(dataNamesSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("sumb_getMeshData", &
                          "Memory allocation failure for dataNamesSUmb")

         if(allocated(iwSUmb)) then
           deallocate(iwSUmb, stat=ierr)
           if(ierr /= 0)                        &
             call terminate("sumb_getMeshData", &
                            "Memory deallocation failure for iwSUmb")
         endif

         allocate(iwSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("sumb_getMeshData", &
                          "Memory allocation failure for iwSUmb")

         dataNamesSUmb(1) = cgnsDensity
         iwSUmb(1) = irho

         dataNamesSUmb(2) = cgnsVelX
         iwSUmb(2) = ivx

         dataNamesSUmb(3) = cgnsVelY
         iwSUmb(3) = ivy

         dataNamesSUmb(4) = cgnsVelZ
         iwSUmb(4) = ivz

         dataNamesSUmb(5) = cgnsEnergy
         iwSUmb(5) = irhoE

         ietcVar = 6

         if(equations == RANSEquations) then
           select case (turbModel)

             case (spalartAllmaras, spalartAllmarasEdwards)
               dataNamesSUmb(6) = cgnsTurbSANu
               iwSUmb(6) = itu1

               ietcVar = 7

             case (komegaWilcox, komegaModified, menterSST)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbOmega
               iwSUmb(7) = itu2

               ietcVar = 8

             case (ktau)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbTau
               iwSUmb(7) = itu2

               ietcVar = 8

             case (v2f)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbEpsilon
               iwSUmb(7) = itu2

               dataNamesSUmb(8) = cgnsTurbV2
               iwSUmb(8) = itu3

               dataNamesSUmb(9) = cgnsTurbF
               iwSUmb(9) = itu4

               ietcVar = 10

           end select
         endif

         dataNamesSUmb(ietcVar) = cgnsPressure
         iwSUmb(ietcVar) = 0

         dataNamesSUmb(ietcVar+1) = cgnsPTot
         iwSUmb(ietcVar+1) = 0

         dataNamesSUmb(ietcVar+2) = cgnsTTot
         iwSUmb(ietcVar+2) = 0

         dataNamesSUmb(ietcVar+3) = cgnsVelVecX
         iwSUmb(ietcVar+3) = 0

         dataNamesSUmb(ietcVar+4) = cgnsVelVecY
         iwSUmb(ietcVar+4) = 0

         dataNamesSUmb(ietcVar+5) = cgnsVelVecZ
         iwSUmb(ietcVar+5) = 0
       endif

       ! Initialize the array data.
       data = 0.

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       dataLoop: do ll=1,nDataNames
         tmpName = dataNames(ll)
       ! call convertToLowerCase(tmpName)
         tmpName = adjustl(tmpName)
         tmpName = trim(tmpName)
         if(tmpName(1:3) == 'AVG' .or. &
            tmpName(1:3) == 'Avg' .or. &
            tmpName(1:3) == 'avg')     &
           tmpName = tmpName(4:)

         do kk=1,nDataSUmb
         ! call convertToLowerCase(tmpName1(kk))
           dataNamesSUmb(kk) = adjustl(dataNamesSUmb(kk))
           dataNamesSUmb(kk) = trim(dataNamesSUmb(kk))
           if(tmpName == dataNamesSUmb(kk)) then
             iwMapReq(ll) = iwSUmb(kk)
             cycle dataLoop
           endif
         enddo

         write(*,*) "************ SUmb mesh: geomName =", geomName, &
                    ", myID =", myID," ************"
         write(*,*) tmpName, ": This variable cannot be provided from &
                   &the present SUmb domain."
         call terminate("sumb_getMeshData","Program terminated!")
       enddo dataLoop

       nNod = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         ! Passing the variables at vertices to the coupler.

         do kc=1,kl
           kp = kc + 1
           do jc=1,jl
             jp = jc + 1
             do ic=1,il
               ip = ic + 1

               volC = vol(ic,jc,kc) + vol(ic,jc,kp) &
                    + vol(ic,jp,kc) + vol(ic,jp,kp) &
                    + vol(ip,jc,kc) + vol(ip,jc,kp) &
                    + vol(ip,jp,kc) + vol(ip,jp,kp)

               fac1 = vol(ic,jc,kc) / volC
               fac2 = vol(ic,jc,kp) / volC
               fac3 = vol(ic,jp,kc) / volC
               fac4 = vol(ic,jp,kp) / volC
               fac5 = vol(ip,jc,kc) / volC
               fac6 = vol(ip,jc,kp) / volC
               fac7 = vol(ip,jp,kc) / volC
               fac8 = vol(ip,jp,kp) / volC

               nNod = nNod + 1

               do mm=1,nDataNames
                 iwReq = iwMapReq(mm)

                 if(iwReq == 0) then
                   select case(dataNames(mm))
                     case(cgnsPressure)
                       sumbScal1 = p(ip,jp,kp)
                       sumbScal2 = p(ip,jp,kc)
                       sumbScal3 = p(ip,jc,kp)
                       sumbScal4 = p(ip,jc,kc)
                       sumbScal5 = p(ic,jp,kp)
                       sumbScal6 = p(ic,jp,kc)
                       sumbScal7 = p(ic,jc,kp)
                       sumbScal8 = p(ic,jc,kc)

                     case(cgnsPTot)
                       call computePtot(w(ip,jp,kp,1), w(ip,jp,kp,2), &
                                        w(ip,jp,kp,3), w(ip,jp,kp,4), &
                                        p(ip,jp,kp), sumbScal1, 1)
                       call computePtot(w(ip,jp,kc,1), w(ip,jp,kc,2), &
                                        w(ip,jp,kc,3), w(ip,jp,kc,4), &
                                        p(ip,jp,kc), sumbScal2, 1)
                       call computePtot(w(ip,jc,kp,1), w(ip,jc,kp,2), &
                                        w(ip,jc,kp,3), w(ip,jc,kp,4), &
                                        p(ip,jc,kp), sumbScal3, 1)
                       call computePtot(w(ip,jc,kc,1), w(ip,jc,kc,2), &
                                        w(ip,jc,kc,3), w(ip,jc,kc,4), &
                                        p(ip,jc,kc), sumbScal4, 1)
                       call computePtot(w(ic,jp,kp,1), w(ic,jp,kp,2), &
                                        w(ic,jp,kp,3), w(ic,jp,kp,4), &
                                        p(ic,jp,kp), sumbScal5, 1)
                       call computePtot(w(ic,jp,kc,1), w(ic,jp,kc,2), &
                                        w(ic,jp,kc,3), w(ic,jp,kc,4), &
                                        p(ic,jp,kc), sumbScal6, 1)
                       call computePtot(w(ic,jc,kp,1), w(ic,jc,kp,2), &
                                        w(ic,jc,kp,3), w(ic,jc,kp,4), &
                                        p(ic,jc,kp), sumbScal7, 1)
                       call computePtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal8, 1)

                     case(cgnsTTot)
                       call computeTtot(w(ip,jp,kp,1), w(ip,jp,kp,2), &
                                        w(ip,jp,kp,3), w(ip,jp,kp,4), &
                                        p(ip,jp,kp), sumbScal1, 1)
                       call computeTtot(w(ip,jp,kc,1), w(ip,jp,kc,2), &
                                        w(ip,jp,kc,3), w(ip,jp,kc,4), &
                                        p(ip,jp,kc), sumbScal2, 1)
                       call computeTtot(w(ip,jc,kp,1), w(ip,jc,kp,2), &
                                        w(ip,jc,kp,3), w(ip,jc,kp,4), &
                                        p(ip,jc,kp), sumbScal3, 1)
                       call computeTtot(w(ip,jc,kc,1), w(ip,jc,kc,2), &
                                        w(ip,jc,kc,3), w(ip,jc,kc,4), &
                                        p(ip,jc,kc), sumbScal4, 1)
                       call computeTtot(w(ic,jp,kp,1), w(ic,jp,kp,2), &
                                        w(ic,jp,kp,3), w(ic,jp,kp,4), &
                                        p(ic,jp,kp), sumbScal5, 1)
                       call computeTtot(w(ic,jp,kc,1), w(ic,jp,kc,2), &
                                        w(ic,jp,kc,3), w(ic,jp,kc,4), &
                                        p(ic,jp,kc), sumbScal6, 1)
                       call computeTtot(w(ic,jc,kp,1), w(ic,jc,kp,2), &
                                        w(ic,jc,kp,3), w(ic,jc,kp,4), &
                                        p(ic,jc,kp), sumbScal7, 1)
                       call computeTtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal8, 1)

                     case (cgnsVelVecX, cgnsVelVecY, cgnsVelVecZ)
                       sumbVecX1 = w(ip,jp,kp,2)
                       sumbVecX2 = w(ip,jp,kc,2)
                       sumbVecX3 = w(ip,jc,kp,2)
                       sumbVecX4 = w(ip,jc,kc,2)
                       sumbVecX5 = w(ic,jp,kp,2)
                       sumbVecX6 = w(ic,jp,kc,2)
                       sumbVecX7 = w(ic,jc,kp,2)
                       sumbVecX8 = w(ic,jc,kc,2)

                       sumbVecY1 = w(ip,jp,kp,3)
                       sumbVecY2 = w(ip,jp,kc,3)
                       sumbVecY3 = w(ip,jc,kp,3)
                       sumbVecY4 = w(ip,jc,kc,3)
                       sumbVecY5 = w(ic,jp,kp,3)
                       sumbVecY6 = w(ic,jp,kc,3)
                       sumbVecY7 = w(ic,jc,kp,3)
                       sumbVecY8 = w(ic,jc,kc,3)

                       sumbVecZ1 = w(ip,jp,kp,4)
                       sumbVecZ2 = w(ip,jp,kc,4)
                       sumbVecZ3 = w(ip,jc,kp,4)
                       sumbVecZ4 = w(ip,jc,kc,4)
                       sumbVecZ5 = w(ic,jp,kp,4)
                       sumbVecZ6 = w(ic,jp,kc,4)
                       sumbVecZ7 = w(ic,jc,kp,4)
                       sumbVecZ8 = w(ic,jc,kc,4)

                       sumbVecXC = fac1*sumbVecX1 &
                                 + fac2*sumbVecX2 &
                                 + fac3*sumbVecX3 &
                                 + fac4*sumbVecX4 &
                                 + fac5*sumbVecX5 &
                                 + fac6*sumbVecX6 &
                                 + fac7*sumbVecX7 &
                                 + fac8*sumbVecX8

                       sumbVecYC = fac1*sumbVecY1 &
                                 + fac2*sumbVecY2 &
                                 + fac3*sumbVecY3 &
                                 + fac4*sumbVecY4 &
                                 + fac5*sumbVecY5 &
                                 + fac6*sumbVecY6 &
                                 + fac7*sumbVecY7 &
                                 + fac8*sumbVecY8

                       sumbVecZC = fac1*sumbVecZ1 &
                                 + fac2*sumbVecZ2 &
                                 + fac3*sumbVecZ3 &
                                 + fac4*sumbVecZ4 &
                                 + fac5*sumbVecZ5 &
                                 + fac6*sumbVecZ6 &
                                 + fac7*sumbVecZ7 &
                                 + fac8*sumbVecZ8

                   end select

                   select case (dataNames(mm))
                     case (cgnsPressure, cgnsPTot, cgnsTTot)
                       data(mm,nNod) = fac1*sumbScal1 &
                                     + fac2*sumbScal2 &
                                     + fac3*sumbScal3 &
                                     + fac4*sumbScal4 &
                                     + fac5*sumbScal5 &
                                     + fac6*sumbScal6 &
                                     + fac7*sumbScal7 &
                                     + fac8*sumbScal8

                     case(cgnsVelVecX)
                       sumbVecMag = sqrt(sumbVecXC**2 &
                                       + sumbVecYC**2 &
                                       + sumbVecZC**2)

                       if(sumbVecMag > eps) then
                         data(mm,nNod) = sumbVecXC/sumbVecMag
                       else
                         data(mm,nNod) = 0.
                       endif

                     case(cgnsVelVecY)
                       sumbVecMag = sqrt(sumbVecXC**2 &
                                       + sumbVecYC**2 &
                                       + sumbVecZC**2)

                       if(sumbVecMag > eps) then
                         data(mm,nNod) = sumbVecYC/sumbVecMag
                       else
                         data(mm,nNod) = 0.
                       endif

                     case(cgnsVelVecZ)
                       sumbVecMag = sqrt(sumbVecXC**2 &
                                       + sumbVecYC**2 &
                                       + sumbVecZC**2)

                       if(sumbVecMag > eps) then
                         data(mm,nNod) = sumbVecZC/sumbVecMag
                       else
                         data(mm,nNod) = 0.
                       endif

                   end select

                 else
                   data(mm,nNod) = fac1*w(ip,jp,kp,iwReq) &
                                 + fac2*w(ip,jp,kc,iwReq) &
                                 + fac3*w(ip,jc,kp,iwReq) &
                                 + fac4*w(ip,jc,kc,iwReq) &
                                 + fac5*w(ic,jp,kp,iwReq) &
                                 + fac6*w(ic,jp,kc,iwReq) &
                                 + fac7*w(ic,jc,kp,iwReq) &
                                 + fac8*w(ic,jc,kc,iwReq)
                 endif
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop

       end subroutine sumb_getMeshData

!      ==================================================================

       subroutine sumb_getMeshDataC(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshDataC provides the coupler with the requested flow *
!      * fields. This is the counterpart of sumb_getMeshData for cell   *
!      * centers.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: kk, ll, nCel, nn, ic, jc, kc, &
                                mm, ietcVar, iwReq
       integer(kind=intType), dimension(nDataNames) :: iwMapReq
     ! character(len=3) :: prefix, prefixTurb
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
       real(kind=realType) :: sumbScal, sumbVecMag
       real(kind=realType) :: sumbVecX, sumbVecY, sumbVecZ
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Prefix to tell the averaged variables from instantaneous ones.

     ! prefix = 'Avg'
     ! prefixTurb = 'Avg'

       ! First, define a list of flow variables which this particular
       ! instance of SUmb can provide to the coupler.
       ! More variables can be added here.

       if(.not. allocated(dataNamesSUmb)) then
         nDataSUmb = nw + 6
         allocate(dataNamesSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("sumb_getMeshDataC", &
                          "Memory allocation failure for dataNamesSUmb")

         if(allocated(iwSUmb)) then
           deallocate(iwSUmb, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("sumb_getMeshDataC", &
                            "Memory deallocation failure for iwSUmb")
         endif

         allocate(iwSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("sumb_getMeshDataC", &
                          "Memory allocation failure for iwSUmb")

         dataNamesSUmb(1) = cgnsDensity
         iwSUmb(1) = irho

         dataNamesSUmb(2) = cgnsVelX
         iwSUmb(2) = ivx

         dataNamesSUmb(3) = cgnsVelY
         iwSUmb(3) = ivy

         dataNamesSUmb(4) = cgnsVelZ
         iwSUmb(4) = ivz

         dataNamesSUmb(5) = cgnsEnergy
         iwSUmb(5) = irhoE

         ietcVar = 6

         if(equations == RANSEquations) then
           select case (turbModel)

             case (spalartAllmaras, spalartAllmarasEdwards)
               dataNamesSUmb(6) = cgnsTurbSANu
               iwSUmb(6) = itu1

               ietcVar = 7

             case (komegaWilcox, komegaModified, menterSST)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbOmega
               iwSUmb(7) = itu2

               ietcVar = 8

             case (ktau)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbTau
               iwSUmb(7) = itu2

               ietcVar = 8

             case (v2f)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbEpsilon
               iwSUmb(7) = itu2

               dataNamesSUmb(8) = cgnsTurbV2
               iwSUmb(8) = itu3

               dataNamesSUmb(9) = cgnsTurbF
               iwSUmb(9) = itu4

               ietcVar = 10

           end select
         endif

         dataNamesSUmb(ietcVar) = cgnsPressure
         iwSUmb(ietcVar) = 0

         dataNamesSUmb(ietcVar+1) = cgnsPTot
         iwSUmb(ietcVar+1) = 0

         dataNamesSUmb(ietcVar+2) = cgnsTTot
         iwSUmb(ietcVar+2) = 0

         dataNamesSUmb(ietcVar+3) = cgnsVelVecX
         iwSUmb(ietcVar+3) = 0

         dataNamesSUmb(ietcVar+4) = cgnsVelVecY
         iwSUmb(ietcVar+4) = 0

         dataNamesSUmb(ietcVar+5) = cgnsVelVecZ
         iwSUmb(ietcVar+5) = 0
       endif

       ! Initialize the array data.
       data = 0.

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       dataLoop: do ll=1,nDataNames
         tmpName = dataNames(ll)
       ! call convertToLowerCase(tmpName)
         tmpName = adjustl(tmpName)
         tmpName = trim(tmpName)
         if(tmpName(1:3) == 'AVG' .or. &
            tmpName(1:3) == 'Avg' .or. &
            tmpName(1:3) == 'avg')     &
           tmpName = tmpName(4:)

         do kk=1,nDataSUmb
         ! call convertToLowerCase(tmpName1(kk))
           dataNamesSUmb(kk) = adjustl(dataNamesSUmb(kk))
           dataNamesSUmb(kk) = trim(dataNamesSUmb(kk))
           if(tmpName == dataNamesSUmb(kk)) then
             iwMapReq(ll) = iwSUmb(kk)
             cycle dataLoop
           endif
         enddo

         write(*,*) "************ SUmb mesh: geomName =", geomName, &
                    ", myID =", myID," ************"
         write(*,*) tmpName, ": This variable cannot be provided from &
                   &the present SUmb domain."
         call terminate("sumb_getMeshDataC","Program terminated!")
       enddo dataLoop

       nCel = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         ! Passing the variables at cell centers to the coupler.

         do kc=2,kl
           do jc=2,jl
             do ic=2,il

               nCel = nCel + 1

               do mm=1,nDataNames
                 iwReq = iwMapReq(mm)

                 if(iwReq == 0) then
                   select case(dataNames(mm))
                     case(cgnsPressure)
                       sumbScal = p(ic,jc,kc)

                     case(cgnsPTot)
                       call computePtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal, 1)

                     case(cgnsTTot)
                       call computeTtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal, 1)

                     case (cgnsVelVecX, cgnsVelVecY, cgnsVelVecZ)
                       sumbVecX = w(ic,jc,kc,2)
                       sumbVecY = w(ic,jc,kc,3)
                       sumbVecZ = w(ic,jc,kc,4)

                   end select

                   select case (dataNames(mm))
                     case (cgnsPressure, cgnsPTot, cgnsTTot)
                       data(mm,nCel) = sumbScal

                     case(cgnsVelVecX)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecX/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                     case(cgnsVelVecY)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecY/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                     case(cgnsVelVecZ)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecZ/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                   end select

                 else
                   data(mm,nCel) = w(ic,jc,kc,iwReq)
                 endif
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop

       end subroutine sumb_getMeshDataC

!      ==================================================================

       subroutine sumb_getMeshDataCH(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getMeshDataCH provides the coupler with the requested     *
!      * flow fields. This is the counterpart of sumb_getMeshData for   *
!      * cell centers + 1st-level halos.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(out) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: kk, ll, nCel, nn, ic, jc, kc, &
                                mm, ietcVar, iwReq
       integer(kind=intType), dimension(nDataNames) :: iwMapReq
     ! character(len=3) :: prefix, prefixTurb
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
       real(kind=realType) :: sumbScal, sumbVecMag
       real(kind=realType) :: sumbVecX, sumbVecY, sumbVecZ
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Prefix to tell the averaged variables from instantaneous ones.

     ! prefix = 'Avg'
     ! prefixTurb = 'Avg'

       ! First, define a list of flow variables which this particular
       ! instance of SUmb can provide to the coupler.
       ! More variables can be added here.

       if(.not. allocated(dataNamesSUmb)) then
         nDataSUmb = nw + 6
         allocate(dataNamesSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("sumb_getMeshDataC", &
                          "Memory allocation failure for dataNamesSUmb")

         if(allocated(iwSUmb)) then
           deallocate(iwSUmb, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("sumb_getMeshDataC", &
                            "Memory deallocation failure for iwSUmb")
         endif

         allocate(iwSUmb(nDataSUmb), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("sumb_getMeshDataC", &
                          "Memory allocation failure for iwSUmb")

         dataNamesSUmb(1) = cgnsDensity
         iwSUmb(1) = irho

         dataNamesSUmb(2) = cgnsVelX
         iwSUmb(2) = ivx

         dataNamesSUmb(3) = cgnsVelY
         iwSUmb(3) = ivy

         dataNamesSUmb(4) = cgnsVelZ
         iwSUmb(4) = ivz

         dataNamesSUmb(5) = cgnsEnergy
         iwSUmb(5) = irhoE

         ietcVar = 6

         if(equations == RANSEquations) then
           select case (turbModel)

             case (spalartAllmaras, spalartAllmarasEdwards)
               dataNamesSUmb(6) = cgnsTurbSANu
               iwSUmb(6) = itu1

               ietcVar = 7

             case (komegaWilcox, komegaModified, menterSST)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbOmega
               iwSUmb(7) = itu2

               ietcVar = 8

             case (ktau)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbTau
               iwSUmb(7) = itu2

               ietcVar = 8

             case (v2f)
               dataNamesSUmb(6) = cgnsTurbK
               iwSUmb(6) = itu1

               dataNamesSUmb(7) = cgnsTurbEpsilon
               iwSUmb(7) = itu2

               dataNamesSUmb(8) = cgnsTurbV2
               iwSUmb(8) = itu3

               dataNamesSUmb(9) = cgnsTurbF
               iwSUmb(9) = itu4

               ietcVar = 10

           end select
         endif

         dataNamesSUmb(ietcVar) = cgnsPressure
         iwSUmb(ietcVar) = 0

         dataNamesSUmb(ietcVar+1) = cgnsPTot
         iwSUmb(ietcVar+1) = 0

         dataNamesSUmb(ietcVar+2) = cgnsTTot
         iwSUmb(ietcVar+2) = 0

         dataNamesSUmb(ietcVar+3) = cgnsVelVecX
         iwSUmb(ietcVar+3) = 0

         dataNamesSUmb(ietcVar+4) = cgnsVelVecY
         iwSUmb(ietcVar+4) = 0

         dataNamesSUmb(ietcVar+5) = cgnsVelVecZ
         iwSUmb(ietcVar+5) = 0
       endif

       ! Initialize the array data.
       data = 0.

       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       dataLoop: do ll=1,nDataNames
         tmpName = dataNames(ll)
       ! call convertToLowerCase(tmpName)
         tmpName = adjustl(tmpName)
         tmpName = trim(tmpName)
         if(tmpName(1:3) == 'AVG' .or. &
            tmpName(1:3) == 'Avg' .or. &
            tmpName(1:3) == 'avg')     &
           tmpName = tmpName(4:)

         do kk=1,nDataSUmb
         ! call convertToLowerCase(tmpName1(kk))
           dataNamesSUmb(kk) = adjustl(dataNamesSUmb(kk))
           dataNamesSUmb(kk) = trim(dataNamesSUmb(kk))
           if(tmpName == dataNamesSUmb(kk)) then
             iwMapReq(ll) = iwSUmb(kk)
             cycle dataLoop
           endif
         enddo

         write(*,*) "************ SUmb mesh: geomName =", geomName, &
                    ", myID =", myID," ************"
         write(*,*) tmpName, ": This variable cannot be provided from &
                   &the present SUmb domain."
         call terminate("sumb_getMeshDataCH","Program terminated!")
       enddo dataLoop

       nCel = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         ! Passing the variables at cell centers to the coupler.

         do kc=1,kl+1
           do jc=1,jl+1
             do ic=1,il+1

               nCel = nCel + 1

               do mm=1,nDataNames
                 iwReq = iwMapReq(mm)

                 if(iwReq == 0) then
                   select case(dataNames(mm))
                     case(cgnsPressure)
                       sumbScal = p(ic,jc,kc)

                     case(cgnsPTot)
                       call computePtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal, 1)

                     case(cgnsTTot)
                       call computeTtot(w(ic,jc,kc,1), w(ic,jc,kc,2), &
                                        w(ic,jc,kc,3), w(ic,jc,kc,4), &
                                        p(ic,jc,kc), sumbScal, 1)

                     case (cgnsVelVecX, cgnsVelVecY, cgnsVelVecZ)
                       sumbVecX = w(ic,jc,kc,2)
                       sumbVecY = w(ic,jc,kc,3)
                       sumbVecZ = w(ic,jc,kc,4)

                   end select

                   select case (dataNames(mm))
                     case (cgnsPressure, cgnsPTot, cgnsTTot)
                       data(mm,nCel) = sumbScal

                     case(cgnsVelVecX)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecX/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                     case(cgnsVelVecY)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecY/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                     case(cgnsVelVecZ)
                       sumbVecMag = sqrt(sumbVecX**2 &
                                       + sumbVecY**2 &
                                       + sumbVecZ**2)

                       if(sumbVecMag > eps) then
                         data(mm,nCel) = sumbVecZ/sumbVecMag
                       else
                         data(mm,nCel) = 0.
                       endif

                   end select

                 else
                   data(mm,nCel) = w(ic,jc,kc,iwReq)
                 endif
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop

       end subroutine sumb_getMeshDataCH

!      ==================================================================

       subroutine sumb_setMeshData(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setMeshData gets the interpolated solution at the domain  *
!      * mesh points from the coupler.                                  *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(in) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l, m, n, ll, mm, nn
       integer(kind=intType) :: nNod, nCel
       integer(kind=intType), dimension(8) :: iNod
       integer(kind=intType), allocatable, dimension(:,:,:) :: indMap
       integer(kind=intType), allocatable, dimension(:,:) :: tmpConnHexa
       real(kind=realType) :: tmpData
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       ! First, restore the connectivity information.

       nNod = 0
       nCel = 0

       if(allocated(tmpConnHexa)) then
         deallocate(tmpConnHexa, stat=ierr)
         if(ierr /= 0)                        &
           call terminate("sumb_setMeshData", &
                          "Memory deallocation failure for tmpConnHexa")
       endif

       allocate(tmpConnHexa(8,nHexaAlloc), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("sumb_setMeshData", &
                        "Memory allocation failure for tmpConnHexa")

       domainLoop1: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         allocate(indMap(il,jl,kl), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("sumb_setMeshData", &
                          "Memory allocation failure for indMap")

         ! Define a map between the one- and three-dimensional
         ! indices.

         do k=1,kl
           do j=1,jl
             do i=1,il
               nNod = nNod + 1
               indMap(i,j,k) = nNod
             enddo
           enddo
         enddo

         ! Creating the hexahedral connectivity

         do k=1,nz
           n = k + 1
           do j=1,ny
             m = j + 1
             do i=1,nx
               l = i + 1

               nCel = nCel + 1
               tmpConnHexa(1,nCel) = indMap(i,j,k)
               tmpConnHexa(2,nCel) = indMap(l,j,k)
               tmpConnHexa(3,nCel) = indMap(l,m,k)
               tmpConnHexa(4,nCel) = indMap(i,m,k)
               tmpConnHexa(5,nCel) = indMap(i,j,n)
               tmpConnHexa(6,nCel) = indMap(l,j,n)
               tmpConnHexa(7,nCel) = indMap(l,m,n)
               tmpConnHexa(8,nCel) = indMap(i,m,n)
             enddo
           enddo
         enddo

         deallocate(indMap, stat=ierr)
         if(ierr /= 0)                        &
           call terminate("sumb_setMeshData", &
                          "Memory deallocation failure for indMap")
       enddo domainLoop1

       ! Reinitialize the cell counter.

       nCel = 0

       domainLoop2: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         do k=1,nz
           n = k + 1
           do j=1,ny
             m = j + 1
             do i=1,nx
               l = i + 1

               nCel = nCel + 1

               do ll=1,8
                 iNod(ll) = tmpConnHexa(ll,nCel)
               enddo

               do mm=1,nDataNames
                 tmpName = dataNames(mm)
               ! call convertToLowerCase(tmpName)
                 tmpName = adjustl(tmpName)
                 tmpName = trim(tmpName)
                 if(tmpName(1:3) == 'AVG' .or. &
                    tmpName(1:3) == 'Avg' .or. &
                    tmpName(1:3) == 'avg')     &
                   tmpName = tmpName(4:)
                 tmpData = eighth*( data(mm,iNod(1)) &
                                  + data(mm,iNod(2)) &
                                  + data(mm,iNod(3)) &
                                  + data(mm,iNod(4)) &
                                  + data(mm,iNod(5)) &
                                  + data(mm,iNod(6)) &
                                  + data(mm,iNod(7)) &
                                  + data(mm,iNod(8)) )

                 select case(tmpName)
                   case(cgnsDensity)
                     w(l,m,n,irho)  = tmpData
                   case(cgnsVelX)
                     w(l,m,n,ivx)   = tmpData
                   case(cgnsVelY)
                     w(l,m,n,ivy)   = tmpData
                   case(cgnsVelZ)
                     w(l,m,n,ivz)   = tmpData
                   case(cgnsEnergy)
                     w(l,m,n,irhoE) = tmpData
                   case (cgnsTurbSANu, cgnsTurbK)
                     w(l,m,n,itu1)  = tmpData
                   case (cgnsTurbOmega, cgnsTurbTau, &
                         cgnsTurbEpsilon)
                     w(l,m,n,itu2)  = tmpData
                   case (cgnsTurbV2)
                     w(l,m,n,itu3)  = tmpData
                   case (cgnsTurbF)
                     w(l,m,n,itu4)  = tmpData
                   case(cgnsPressure)
                     p(l,m,n)       = tmpData
                   case default
                     ! Write the message only once.

                     if(i == 1 .and. j == 1 .and. k == 1) then
                       write(*,*) 'Setting ', tmpName, ' as a volume &
                                  &variable is not ready yet.'
                       write(*,*) 'Information is ignored.'
                     endif
                 end select
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop2

       ! Deallocate the connectivity array.

       deallocate(tmpConnHexa, stat=ierr)
       if(ierr /= 0)                        &
         call terminate("sumb_setMeshData", &
                        "Memory deallocation failure for tmpConnHexa")

       end subroutine sumb_setMeshData

!      ==================================================================

       subroutine sumb_setMeshDataC(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setMeshDataC gets the interpolated solution at the domain *
!      * mesh points from the coupler. This is the counterpart of       *
!      * sumb_setMeshData for cell centers.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(in) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, mm, nn
       integer(kind=intType) :: nCel
       real(kind=realType) :: tmpData
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       ! Initialize the cell counter.

       nCel = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         do k=2,kl
           do j=2,jl
             do i=2,il

               nCel = nCel + 1

               do mm=1,nDataNames
                 tmpName = dataNames(mm)
               ! call convertToLowerCase(tmpName)
                 tmpName = adjustl(tmpName)
                 tmpName = trim(tmpName)
                 if(tmpName(1:3) == 'AVG' .or. &
                    tmpName(1:3) == 'Avg' .or. &
                    tmpName(1:3) == 'avg')     &
                   tmpName = tmpName(4:)
                 tmpData = data(mm,nCel)

                 select case(tmpName)
                   case(cgnsDensity)
                     w(i,j,k,irho)  = tmpData
                   case(cgnsVelX)
                     w(i,j,k,ivx)   = tmpData
                   case(cgnsVelY)
                     w(i,j,k,ivy)   = tmpData
                   case(cgnsVelZ)
                     w(i,j,k,ivz)   = tmpData
                   case(cgnsEnergy)
                     w(i,j,k,irhoE) = tmpData
                   case (cgnsTurbSANu, cgnsTurbK)
                     w(i,j,k,itu1)  = tmpData
                   case (cgnsTurbOmega, cgnsTurbTau, &
                         cgnsTurbEpsilon)
                     w(i,j,k,itu2)  = tmpData
                   case (cgnsTurbV2)
                     w(i,j,k,itu3)  = tmpData
                   case (cgnsTurbF)
                     w(i,j,k,itu4)  = tmpData
                   case(cgnsPressure)
                     p(i,j,k)       = tmpData
                   case default
                     ! Write the message only once.

                     if(i == 2 .and. j == 2 .and. k == 2) then
                       write(*,*) 'Setting ', tmpName, ' as a volume &
                                  &variable is not ready yet.'
                       write(*,*) 'Information is ignored.'
                     endif
                 end select
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop

       end subroutine sumb_setMeshDataC

!      ==================================================================

       subroutine sumb_setMeshDataCH(data, dataNames, nDataNames, geomName)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setMeshDataCH gets the interpolated solution at the       *
!      * domain mesh points from the coupler. This is the counterpart   *
!      * of sumb_setMeshData for cell centers + 1st-level halos.        *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use communication
       use constants
       use couplerParam
       use cgnsNames
       use flowVarRefState
       use inputPhysics

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:,:), intent(in) :: data
       integer(kind=intType), intent(in) :: nDataNames
       character(len=maxCplNameLen), dimension(nDataNames), &
                                     intent(in) :: dataNames
       character(len=*), intent(in) :: geomName
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, mm, nn
       integer(kind=intType) :: nCel
       real(kind=realType) :: tmpData
       character(len=maxCplNameLen) :: trimCodeName, trimGeomName
       character(len=maxCplNameLen) :: tmpName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       trimCodeName = codeName
       call convertToLowerCase(trimCodeName)
       trimCodeName = adjustl(trimCodeName)
       trimCodeName = trim(trimCodeName)

       trimGeomName = geomName
       call convertToLowerCase(trimGeomName)
       trimGeomName = adjustl(trimGeomName)
       trimGeomName = trim(trimGeomName)

       if(trimCodeName /= trimGeomName) return

       ! Initialize the cell counter.

       nCel = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers to this block.
         ! Only the finest level is considered.

         call setPointers(nn, 1, 1)

         do k=1,kl+1
           do j=1,jl+1
             do i=1,il+1

               nCel = nCel + 1

               do mm=1,nDataNames
                 tmpName = dataNames(mm)
               ! call convertToLowerCase(tmpName)
                 tmpName = adjustl(tmpName)
                 tmpName = trim(tmpName)
                 if(tmpName(1:3) == 'AVG' .or. &
                    tmpName(1:3) == 'Avg' .or. &
                    tmpName(1:3) == 'avg')     &
                   tmpName = tmpName(4:)
                 tmpData = data(mm,nCel)

                 select case(tmpName)
                   case(cgnsDensity)
                     w(i,j,k,irho)  = tmpData
                   case(cgnsVelX)
                     w(i,j,k,ivx)   = tmpData
                   case(cgnsVelY)
                     w(i,j,k,ivy)   = tmpData
                   case(cgnsVelZ)
                     w(i,j,k,ivz)   = tmpData
                   case(cgnsEnergy)
                     w(i,j,k,irhoE) = tmpData
                   case (cgnsTurbSANu, cgnsTurbK)
                     w(i,j,k,itu1)  = tmpData
                   case (cgnsTurbOmega, cgnsTurbTau, &
                         cgnsTurbEpsilon)
                     w(i,j,k,itu2)  = tmpData
                   case (cgnsTurbV2)
                     w(i,j,k,itu3)  = tmpData
                   case (cgnsTurbF)
                     w(i,j,k,itu4)  = tmpData
                   case(cgnsPressure)
                     p(i,j,k)       = tmpData
                   case default
                     ! Write the message only once.

                     if(i == 1 .and. j == 1 .and. k == 1) then
                       write(*,*) 'Setting ', tmpName, ' as a volume &
                                  &variable is not ready yet.'
                       write(*,*) 'Information is ignored.'
                     endif
                 end select
               enddo
             enddo
           enddo
         enddo
       enddo domainLoop

       end subroutine sumb_setMeshDataCH

!      ==================================================================

       subroutine sumb_getParamInt(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamInt provides the current value of the integer     *
!      * parameter with the given name to the coupler.                  *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("discretization scheme")
           value = spaceDiscr

         case ("discretization scheme coarse grid")
           value = spaceDiscrCoarse

         case ("order turbulent equations")
           value = orderTurb

         case ("riemann solver")
           value = riemann

         case ("riemann solver coarse grid")
           value = riemannCoarse

         case ("limiter")
           value = limiter

         case ("preconditioner")
           value = precond

         case ("wall boundary treatment")
           value = wallBCTreatment

         case ("outflow boundary treatment")
           value = outflowTreatment

         case ("non-matching block to block treatment")
           value = nonMatchTreatment
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("file format read")
           value = fileFormatRead

         case ("file format write")
           value = fileFormatWrite

         case ("write precision grid")
           value = precisionGrid

         case ("write precision solution")
           value = precisionSol
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("smoother")
           value = smoother

         case ("treatment turbulent equations")
           value = turbTreatment

         case ("turbulent smoother")
           value = turbSmoother

         case ("turbulent relaxation")
           value = turbRelax

         case ("residual averaging")
           value = resAveraging

         case ("treatment boundary multigrid corrections")
           value = mgBoundCorr

         case ("number of multigrid cycles")
           value = nCycles

         case ("number of single grid startup iterations")
           value = nsgStartup

         case ("save every")
           value = nSaveVolume

         case ("save surface every")
           value = nSaveSurface

         case ("number of runge kutta stages")
           value = nRKStages

         case ("number of multigrid cycles coarse grid")
           value = nCyclesCoarse

         case ("multigrid start level")
           value = mgStartlevel
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Polynomial rotation parameters.

         case ("degree polynomial x-rotation")
           value = degreePolXRot

         case ("degree polynomial y-rotation")
           value = degreePolYRot

         case ("degree polynomial z-rotation")
           value = degreePolZRot

         ! Fourier rotation parameters.

         case ("degree fourier x-rotation")
           value = degreeFourXRot

         case ("degree fourier y-rotation")
           value = degreeFourYRot

         case ("degree fourier z-rotation")
           value = degreeFourZRot
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("equations")
           value = equations

         case ("mode")
           value = equationMode

         case ("flow type")
           value = flowType

         case ("cp model")
           value = cpModel

         case ("turbulence model")
           value = turbModel

         case ("v2f version (n1 or n6)")
           value = rvfN

         case ("turbulence production term")
           value = turbProd
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("number time intervals spectral")
           value = nTimeIntervalsSpectral

         case ("number of unsteady solution files")
           value = nUnsteadySolSpectral
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("time accuracy unsteady")
           value = timeAccuracy

         case ("number of unsteady time steps coarse grid")
           value = nTimeStepsCoarse

         case ("number of unsteady time steps fine grid")
           value = nTimeStepsFine
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("overset interpolation type")
           value = oversetInterpType

         case ("overset interpolation type coarse grid")
           value = oversetInterpTypeCoarse

         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_getParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamInt

!      ==================================================================

       subroutine sumb_getParamIntRk1(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamIntRk1 provides the current value of the          *
!      * parameter defined as the rank-1 integer array having the given *
!      * name to the coupler.                                           *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), dimension(:), intent(inout) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamIntRk1

!      ==================================================================

       subroutine sumb_getParamReal(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamReal provides the current value of the real       *
!      * parameter with the given name to the coupler.                  *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(out) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos, p1, p2
       integer(kind=intType) :: iRotDeg

       character (len=maxStringLen)   :: keyword, tmpName
       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)

       if(keyword(1:34) == "polynomial coefficients x-rotation" .or. &
          keyword(1:34) == "polynomial coefficients y-rotation" .or. &
          keyword(1:34) == "polynomial coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for coefficients of &
                                 &polynomial rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_getParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:34)
         keyword = tmpName

       else if(keyword(1:38) == "fourier cosine coefficients x-rotation" .or. &
               keyword(1:38) == "fourier cosine coefficients y-rotation" .or. &
               keyword(1:38) == "fourier cosine coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for Fourier cosine &
                                 &coefficients of rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_getParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:38)
         keyword = tmpName

       else if(keyword(1:36) == "fourier sine coefficients x-rotation" .or. &
               keyword(1:36) == "fourier sine coefficients y-rotation" .or. &
               keyword(1:36) == "fourier sine coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for Fourier sine &
                                 &coefficients of rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_getParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:36)
         keyword = tmpName
       endif
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("vis2")
           value = vis2

         case ("vis4")
           value = vis4

         case ("vis2 coarse grid")
           value = vis2Coarse

         case ("exponent dissipation scaling")
           value = adis

         case ("kappa interpolation value")
           value = kappaCoef
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("residual averaging smoothing parameter")
           value = smoop

         case ("cfl number")
           value = cfl

         case ("alpha turbulent dd-adi")
           value = alfaTurb

         case ("beta turbulent dd-adi")
           value = betaTurb

         case ("relative l2 norm for convergence")
           value = L2Conv

         case ("cfl number coarse grid")
           value = cflCoarse

         case ("relative l2 norm for convergence coarse grid")
           value = L2ConvCoarse

         case ("restriction relaxation factor")
           value = fcoll
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Rotation point.

         case ("rotation point body x")
           value = rotPoint(1)

         case ("rotation point body y")
           value = rotPoint(2)

         case ("rotation point body z")
           value = rotPoint(3)

         ! Polynomial rotation parameters.

         case ("polynomial coefficients x-rotation")
           value = coefPolXRot(iRotDeg)

         case ("polynomial coefficients y-rotation")
           value = coefPolYRot(iRotDeg)

         case ("polynomial coefficients z-rotation")
           value = coefPolZRot(iRotDeg)

         ! Fourier rotation parameters.

         case ("omega fourier x-rotation")
           value = omegaFourXRot

         case ("omega fourier y-rotation")
           value = omegaFourYRot

         case ("omega fourier z-rotation")
           value = omegaFourZRot

         case ("fourier cosine coefficients x-rotation")
           value = cosCoefFourXRot(iRotDeg)

         case ("fourier sine coefficients x-rotation")
           value = sinCoefFourXRot(iRotDeg)

         case ("fourier cosine coefficients y-rotation")
           value = cosCoefFourYRot(iRotDeg)

         case ("fourier sine coefficients y-rotation")
           value = sinCoefFourYRot(iRotDeg)

         case ("fourier cosine coefficients z-rotation")
           value = cosCoefFourZRot(iRotDeg)

         case ("fourier sine coefficients z-rotation")
           value = sinCoefFourZRot(iRotDeg)
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable load imbalance")
           value = loadImbalance
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("offset from wall in wall functions")
           value = wallOffset

         case ("max ratio k-prod/dest")
           value = pklim

         case ("mach")
           value = Mach

         case ("mach for coefficients")
           value = MachCoef

         case ("reynolds")
           value = Reynolds

         case ("free stream velocity direction x")
           value = velDirFreestream(1)

         case ("free stream velocity direction y")
           value = velDirFreestream(2)

         case ("free stream velocity direction z")
           value = velDirFreestream(3)

         case ("lift direction x")
           value = liftDirection(1)

         case ("lift direction y")
           value = liftDirection(2)

         case ("lift direction z")
           value = liftDirection(3)

         case ("reynolds length (in meter)")
           value = ReynoldsLength

         case ("free stream temperature (in k)")
           value = tempFreestream

         case ("constant specific heat ratio")
           value = gammaConstant

         case ("gas constant (j/(kg k))")
           value = RGasDim

         case ("prandtl number")
           value = prandtl

         case ("turbulent prandtl number")
           value = prandtlTurb

         case ("free stream eddy viscosity ratio")
           value = eddyVisInfRatio

         case ("free stream turbulent intensity")
           value = turbIntensityInf

         case ("reference surface")
           value = surfaceRef

         case ("reference length")
           value = lengthRef

         case ("moment reference point x")
           value = pointRef(1)

         case ("moment reference point y")
           value = pointRef(2)

         case ("moment reference point z")
           value = pointRef(3)
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("time step (in sec) for unsteady restart")
           value = dtUnsteadyRestartSpectral
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("unsteady time step (in sec)")
           value = deltaT
!
!        ****************************************************************
!        *                                                              *
!        * Reference state values.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("reference pressure (in pa)")
           value = pRef

         case ("reference density (in kg/m^3)")
           value = rhoRef

         case ("reference temperature (in k)")
           value = TRef

         case ("conversion factor grid units to meter")
           value = LRef
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("mach for initialization")
           value = MachIni

         case ("pressure for initialization")
           value = pIni

         case ("density for initialization")
           value = rhoIni

         case ("velocity direction x for initialization")
           value = velDirIni(1)

         case ("velocity direction y for initialization")
           value = velDirIni(2)

         case ("velocity direction z for initialization")
           value = velDirIni(3)
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable donor quality")
           value = allowableDonorQuality
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_getParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamReal

!      ==================================================================

       subroutine sumb_getParamRealRk1(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamRealRk1 provides the current value of the         *
!      * parameter defined as the rank-1 real array having the given    *
!      * name to the coupler.                                           *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:), intent(out) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("rotation point body (x,y,z)")
           value(1) = rotPoint(1)
           value(2) = rotPoint(2)
           value(3) = rotPoint(3)
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("free stream velocity direction")
           value(1) = velDirFreestream(1)
           value(2) = velDirFreestream(2)
           value(3) = velDirFreestream(3)

         case ("lift direction")
           value(1) = liftDirection(1)
           value(2) = liftDirection(2)
           value(3) = liftDirection(3)

         case ("moment reference point (x,y,z)")
           value(1) = pointRef(1)
           value(2) = pointRef(2)
           value(3) = pointRef(3)
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("velocity direction for initialization")
           value(1) = velDirIni(1)
           value(2) = velDirIni(2)
           value(3) = velDirIni(3)
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamRealRk1

!      ==================================================================

       subroutine sumb_getParamStr(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamStr provides the current value of the string      *
!      * parameter with the given name to the coupler. Note that        *
!      * sumb_getParamStr can handle any parameter once it is in the    *
!      * string format.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       character (len=maxCplNameLen), intent(out) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in string.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * The parameters to monitor the convergence, the surface and   *
!        * extra volume variables to written to the solution files.     *
!        *                                                              *
!        ****************************************************************
!
         case ("monitoring variables")
           call setMonitorVariablesString(value)

         case ("surface output variables")
           call setSurfaceVariablesString(value)

         case ("volume output variables")
           call setVolumeVariablesString(value)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("discretization scheme")
           select case (spaceDiscr)
             case (dissScalar)
               value = "central plus scalar dissipation"
             case (dissMatrix)
               value = "central plus matrix dissipation"
             case (dissCusp)
               value = "central plus cusp dissipation"
             case (upwind)
               value = "upwind"
           end select

         case ("discretization scheme coarse grid")
           select case (spaceDiscrCoarse)
             case (dissScalar)
               value = "central plus scalar dissipation"
             case (dissMatrix)
               value = "central plus matrix dissipation"
             case (dissCusp)
               value = "central plus cusp dissipation"
             case (upwind)
               value = "upwind"
           end select

         case ("order turbulent equations")
           select case (orderTurb)
             case (firstOrder)
               value = "first order"
             case (secondOrder)
               value = "second order"
           end select

         case ("riemann solver")
           select case (riemann)
             case (Roe)
               value = "roe"
             case (vanLeer)
               value = "van leer"
             case (ausmdv)
               value = "ausmdv"
           end select

         case ("riemann solver coarse grid")
           select case (riemannCoarse)
             case (Roe)
               value = "roe"
             case (vanLeer)
               value = "van leer"
             case (ausmdv)
               value = "ausmdv"
           end select

         case ("limiter")
           select case (limiter)
             case (firstOrder)
               value = "first order"
             case (noLimiter)
               value = "no limiter"
             case (vanAlbeda)
               value = "van albeda"
             case (minmod)
               value = "minmod"
           end select

         case ("preconditioner")
           select case (precond)
             case (noPrecond)
               value = "no preconditioner"
             case (Turkel)
               value = "turkel"
             case (ChoiMerkle)
               value = "choi merkle"
           end select

         case ("wall boundary treatment")
           select case (wallBcTreatment)
             case (constantPressure)
               value = "constant pressure"
             case (linExtrapolPressure)
               value = "linear extrapolation pressure"
             case (quadExtrapolPressure)
               value = "quadratic extrapolation pressure"
             case (normalMomentum)
               value = "normal momentum"
           end select

         case ("outflow boundary treatment")
           select case (outflowTreatment)
             case (constantExtrapol)
               value = "constant extrapolation"
             case (linExtrapol)
               value = "linear extrapolation"
           end select

         case ("non-matching block to block treatment")
           select case (nonMatchTreatment)
             case (NonConservative)
               value = "nonconservative"
             case (Conservative)
               value = "conservative"
           end select

         case ("vortex correction")
           call setYesNo(value, keyword, vortexCorr)

         case ("vis2")
           write(value,*) vis2

         case ("vis4")
           write(value,*) vis4

         case ("vis2 coarse grid")
           write(value,*) vis2Coarse

         case ("directional dissipation scaling")
           call setYesNo(value, keyword, dirScaling)

         case ("exponent dissipation scaling")
           write(value,*) adis

         case ("total enthalpy scaling inlet")
           call setYesNo(value, keyword, hScalingInlet)

         case ("kappa interpolation value")
           write(value,*) kappaCoef
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("file format read")
           select case (fileFormatRead)
             case (cgnsFormat)
               value = "cgns"
             case (plot3DFormat)
               value = "plot3d"
           end select

         case ("file format write")
           select case (fileFormatWrite)
             case (cgnsFormat)
               value = "cgns"
             case (plot3DFormat)
               value = "plot3d"
           end select

         case ("grid file")
           value = gridFile

         case ("plot3d connectivity file")
           value = plot3DConnFile

         case ("restart file")
           value = restartFile

         case ("restart")
           call setYesNo(value, keyword, restart)

         case ("check nondimensionalization")
           call setYesNo(value, keyword, checkRestartSol)

         case ("new grid file")
            value = newGridFile

         case ("solution file")
           value = solFile

         case ("surface solution file")
           value = surfaceSolFile

         case ("write precision grid")
           select case (precisionGrid)
             case (precisionSingle)
               value = "single"
             case (precisionDouble)
               value = "double"
           end select

         case ("write precision solution")
           select case (precisionSol)
             case (precisionSingle)
               value = "single"
             case (precisionDouble)
               value = "double"
           end select

         case ("rind layer in solution files")
           call setYesNo(value, keyword, storeRindLayer)

         case ("automatic parameter update")
           call setYesNo(value, keyword, autoParameterUpdate)

         case ("write coordinates in meter")
           call setYesNo(value, keyword, writeCoorMeter)

         case ("cp curve fit file")
           value = cpFile

         case ("store convergence inner iterations")
           call setYesNo(value, keyword, storeConvInnerIter)
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("smoother")
           select case (smoother)
             case (RungeKutta)
               value = "runge kutta"
             case (nlLusgs)
               value = "nonlinear lusgs"
             case (nlLusgsLine)
               value = "nonlinear lusgs line"
           end select

         case ("treatment turbulent equations")
           select case (turbTreatment)
             case (segregated)
               value = "segregated"
             case (coupled)
               value = "coupled"
           end select

         case ("turbulent smoother")
           select case (turbSmoother)
             case (gmres)
               value = "gmres"
             case (adi)
               value = "adi"
           end select

         case ("freeze turbulent source terms in mg")
           call setYesNo(value, keyword, freezeTurbSource)

         case ("turbulent relaxation")
           select case (turbRelax)
             case (turbRelaxExplicit)
               value = "explicit"
             case (turbRelaxImplicit)
               value = "implicit"
           end select

         case ("residual averaging")
           select case (resAveraging)
             case (noResAveraging)
               value = "no"
             case (alwaysResAveraging)
               value = "all stages"
             case (alternateResAveraging)
               value = "alternate stages"
           end select

         case ("residual averaging smoothing parameter")
           write(value,*) smoop

         case ("number of multigrid cycles")
           write(value,*) nCycles

         case ("number of single grid startup iterations")
           write(value,*) nsgStartup

         case ("save every")
           write(value,*) nSaveVolume

         case ("save surface every")
           write(value,*) nSaveSurface

         case ("number of runge kutta stages")
           write(value,*) nRKStages

         case ("cfl number")
           write(value,*) cfl

         case ("alpha turbulent dd-adi")
           write(value,*) alfaTurb

         case ("beta turbulent dd-adi")
           write(value,*) betaTurb

         case ("relative l2 norm for convergence")
           write(value,*) L2Conv

         case ("number of multigrid cycles coarse grid")
           write(value,*) nCyclesCoarse

         case ("cfl number coarse grid")
           write(value,*) cflCoarse

         case ("relative l2 norm for convergence coarse grid")
           write(value,*) L2ConvCoarse

         case ("treatment boundary multigrid corrections")
           select case (mgBoundCorr)
             case (bcDirichlet0)
               value = "zero dirichlet"
             case (bcNeumann)
               value = "neumann"
           end select

         case ("restriction relaxation factor")
           write(value,*) fcoll

         case ("multigrid start level")
           write(value,*) mgStartlevel

!        case ("multigrid cycle strategy")
!          value = mgDescription
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Rotation point.

         case ("rotation point body (x,y,z)")
           write(value,*) rotPoint(1), rotPoint(2), rotPoint(3)

         ! Polynomial rotation parameters.

         case ("degree polynomial x-rotation")
           write(value,*) degreePolXRot

         case ("degree polynomial y-rotation")
           write(value,*) degreePolYRot

         case ("degree polynomial z-rotation")
           write(value,*) degreePolZRot

         case ("polynomial coefficients x-rotation")
           call writeMotionCoef(value, 0_intType, degreePolXRot, &
                                coefPolXRot)

         case ("polynomial coefficients y-rotation")
           call writeMotionCoef(value, 0_intType, degreePolYRot, &
                                coefPolYRot)

         case ("polynomial coefficients z-rotation")
           call writeMotionCoef(value, 0_intType, degreePolZRot, &
                                coefPolZRot)

         ! Fourier rotation parameters.

         case ("degree fourier x-rotation")
           write(value,*) degreeFourXRot

         case ("degree fourier y-rotation")
           write(value,*) degreeFourYRot

         case ("degree fourier z-rotation")
           write(value,*) degreeFourZRot

         case ("omega fourier x-rotation")
           write(value,*) omegaFourXRot

         case ("omega fourier y-rotation")
           write(value,*) omegaFourYRot

         case ("omega fourier z-rotation")
           write(value,*) omegaFourZRot

         case ("fourier cosine coefficients x-rotation")
           call writeMotionCoef(value, 0_intType, degreeFourXRot, &
                                cosCoefFourXRot)

         case ("fourier sine coefficients x-rotation")
           call writeMotionCoef(value, 1_intType, degreeFourXRot, &
                                sinCoefFourXRot)

         case ("fourier cosine coefficients y-rotation")
           call writeMotionCoef(value, 0_intType, degreeFourYRot, &
                                cosCoefFourYRot)

         case ("fourier sine coefficients y-rotation")
           call writeMotionCoef(value, 1_intType, degreeFourYRot, &
                                sinCoefFourYRot)

         case ("fourier cosine coefficients z-rotation")
           call writeMotionCoef(value, 0_intType, degreeFourZRot, &
                                cosCoefFourZRot)

         case ("fourier sine coefficients z-rotation")
           call writeMotionCoef(value, 1_intType, degreeFourZRot, &
                                sinCoefFourZRot)
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable load imbalance")
           write(value,*) loadImbalance

         case ("split blocks for load balance")
           call setYesNo(value, keyword, splitBlocks)
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("equations")
           select case (equations)
             case (EulerEquations)
               value = "euler"
             case (NSEquations)
               value = "laminar ns"
             case (RANSEquations)
               value = "rans"
           end select

         case ("mode")
           select case (equationMode)
             case (steady)
               value = "steady"
             case (unsteady)
               value = "unsteady"
             case (timeSpectral)
               value = "time spectral"
           end select

         case ("flow type")
           select case (flowType)
             case (internalFlow)
               value = "internal flow"
             case (externalFlow)
               value = "external flow"
           end select

         case ("cp model")
           select case (cpModel)
             case (cpConstant)
               value = "constant"

             case (cpTempCurveFits)
               value = "temperature curve fits"
           end select

         case ("turbulence model")
           select case (turbModel)
             case (baldwinLomax)
               value = "baldwin lomax"
             case (spalartAllmaras)
               value = "spalart allmaras"
             case (spalartAllmarasEdwards)
               value = "spalart allmaras edwards"
             case (komegaWilcox)
               value = "komega wilcox"
             case (komegaModified)
               value = "komega modified"
             case (ktau)
               value = "ktau"
             case (menterSST)
               value = "menter sst"
             case (v2f)
               value = "v2f"
           end select

         case ("v2f version (n1 or n6)")
           write(value,*) rvfN

         case ("v2f with upper bound")
           call setYesNo(value, keyword, rvfB)

         case ("turbulence production term")
           select case (turbProd)
             case (strain)
               value = "strain"
             case (vorticity)
               value = "vorticity"
             case (katoLaunder)
               value = "kato-launder"
           end select

         case ("use wall functions")
           call setYesNo(value, keyword, wallFunctions)

         case ("offset from wall in wall functions")
           write(value,*) wallOffset

         case ("max ratio k-prod/dest")
           write(value,*) pklim

         case ("mach")
           write(value,*) Mach

         case ("mach for coefficients")
           write(value,*) MachCoef

         case ("reynolds")
           write(value,*) Reynolds

         case ("free stream velocity direction")
           write(value,*) velDirFreestream(1), velDirFreestream(2), &
                          velDirFreestream(3)

         case ("lift direction")
           write(value,*) liftDirection(1), liftDirection(2), &
                          liftDirection(3)

         case ("reynolds length (in meter)")
           write(value,*) ReynoldsLength

         case ("free stream temperature (in k)")
           write(value,*) tempFreestream

         case ("constant specific heat ratio")
           write(value,*) gammaConstant

         case ("gas constant (j/(kg k))")
           write(value,*) RGasDim

         case ("prandtl number")
           write(value,*) prandtl

         case ("turbulent prandtl number")
           write(value,*) prandtlTurb

         case ("free stream eddy viscosity ratio")
           write(value,*) eddyVisInfRatio

         case ("free stream turbulent intensity")
           write(value,*) turbIntensityInf

         case ("reference surface")
           write(value,*) surfaceRef

         case ("reference length")
           write(value,*) lengthRef

         case ("moment reference point x")
           write(value,*) pointRef(1)

         case ("moment reference point y")
           write(value,*) pointRef(2)

         case ("moment reference point z")
           write(value,*) pointRef(3)

         case ("moment reference point (x,y,z)")
           write(value,*) pointRef(1), pointRef(2), pointRef(3)
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("number time intervals spectral")
           write(value,*) nTimeIntervalsSpectral

         case ("write file for unsteady restart")
           call setYesNo(value, keyword, writeUnsteadyRestartSpectral)

         case ("time step (in sec) for unsteady restart")
           write(value,*) dtUnsteadyRestartSpectral

         case ("write unsteady volume solution files")
           call setYesNo(value, keyword, writeUnsteadyVolSpectral)

         case ("write unsteady surface solution files")
           call setYesNo(value, keyword, writeUnsteadySurfSpectral)

         case ("number of unsteady solution files")
           write(value,*) nUnsteadySolSpectral
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("time accuracy unsteady")
           select case (timeAccuracy)
             case (firstOrder)
               value = "first"
             case (secondOrder)
               value = "second"
             case (thirdOrder)
               value = "third"
           end select

         case ("number of unsteady time steps coarse grid")
           write(value,*) nTimeStepsCoarse

         case ("number of unsteady time steps fine grid")
           write(value,*) nTimeStepsFine

         case ("unsteady time step (in sec)")
           write(value,*) deltaT

         case ("update wall distance unsteady mode")
           call setYesNo(value, keyword, updateWallDistanceUnsteady)
!
!        ****************************************************************
!        *                                                              *
!        * Visualization parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("pv3 visualization only")
           call setYesNo(value, keyword, PV3VisOnly)
!
!        ****************************************************************
!        *                                                              *
!        * Reference state values.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("reference pressure (in pa)")
           write(value,*) pRef

         case ("reference density (in kg/m^3)")
           write(value,*) rhoRef

         case ("reference temperature (in k)")
           write(value,*) TRef

         case ("conversion factor grid units to meter")
           write(value,*) LRef
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("code name")
           value = codeName

         case ("get coarse-level sol")
           call setYesNo(value, keyword, cplGetCoarseSol)

         case ("mach for initialization")
           write(value,*) MachIni

         case ("pressure for initialization")
           write(value,*) pIni

         case ("density for initialization")
           write(value,*) rhoIni

         case ("velocity direction for initialization")
           write(value,*) velDirIni(1), velDirIni(2), velDirIni(3)
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("input overset donors are guesses")
           call setYesNo(value, keyword, oversetDonorsAreGuesses)

         case ("average restricted residual for blanks")
           call setYesNo(value, keyword, avgRestrictResForBlanks)

         case ("overset interpolation type")
           select case (oversetInterpType)
             case (trilinear)
               value = "trilinear"
           end select

         case ("overset interpolation type coarse grid")
           select case (oversetInterpTypeCoarse)
             case (trilinear)
               value = "trilinear"
           end select

         case ("allowable donor quality")
           write(value,*) allowableDonorQuality
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamStr

!      ==================================================================

       subroutine setYesNo(value, keyword, logVal)
!
!      ******************************************************************
!      *                                                                *
!      * setYesNo sets the output string to either "yes" or "no"        *
!      * according to the given logical value (.True. or .False.).      *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
!
!      Subroutine arguments.
!
       character (len=*), intent(out) :: value
       character (len=*), intent(in)  :: keyword
       logical, intent(in)  :: logVal
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if(      logVal) value = "yes"
       if(.not. logVal) value = "no"

       end subroutine setYesNo

!      ==================================================================

       subroutine writeMotionCoef(string, start, end, coef)
!
!      ******************************************************************
!      *                                                                *
!      * writeMotionCoef writes the coefficients from start to end onto *
!      * the given string. These coefficients correspond to the         *
!      * description of the rigid body motion and are either polynomial *
!      * or a fourier series. In both cases it is assumed that the      *
!      * number of coefficients is specified before the actual          *
!      * coefficients are written.                                      *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: start, end
       character(len=*),      intent(out) :: string

       real(kind=realType), dimension(start:*), intent(in) :: coef
!
!      Local variables.
!
       integer :: pos

       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if end >= 0, i.e. if the order of the polynomial/fourier
       ! series is specified before this subroutine is called.

       if(end < 0) then

         ! Processor 0 prints the error message, while the others
         ! wait to get killed.

         if(myID == 0)                       &
           call terminate("writeMotionCoef", &
                          "Order of motion coefficients not known yet &
                          &when the coefficients are requested")
         call mpi_barrier(SUmb_comm_world, pos)
       endif

       ! Write all the coefficients onto the string.

       write(string,*) (coef(i), i=start,end)

       ! The length of the string should be nonzero. If no data
       ! are specified, print a warning to indicate this.

       if(myID == 0 .and. len_trim(string) <= 0) then
         print "(a)", "#"
         print "(a)", "#*==================== !!! Warning !!! &
                      &======================"
         print "(a)", "#* No coefficient written for a &
                      &certain rigid body motion"
         print "(a)", "#* Information is ignored."
         print "(a)", "#*=======================================&
                      &===================="
         print "(a)", "#"
       endif

       end subroutine writeMotionCoef

!      ==================================================================

       subroutine setMonitorVariablesString(variables)
!
!      ******************************************************************
!      *                                                                *
!      * setMonitorVariablesString writes the names of the variables    *
!      * to be monitored during the convergence onto the output string. *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use monitor
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(out) :: variables
!
!      Local variables.
!
       integer :: pos
       integer(kind=intType) :: iMon
       logical :: resMomWritten
       character(len=maxCplNameLen) :: tmpString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       resMomWritten = .False.
       pos = 0

       if(showCPU) then
         write(variables,*) "cpu_"
         pos = pos + 4
       endif

       ! Loop to write the info onto the string variables.

       do iMon=1,nMon
         select case (monNames(iMon))
           case (cgnsL2resRho)
             tmpString = variables(1:pos)//'resrho_'
             pos = pos + 7
             variables = tmpString

           case (cgnsL2resMomx, cgnsL2resMomy, cgnsL2resMomz)
             if(.not. resMomWritten) then
               tmpString = variables(1:pos)//'resmom_'
               pos = pos + 7
               resMomWritten = .True.
             endif

           case (cgnsL2resRhoe)
             tmpString = variables(1:pos)//'resrhoe_'
             pos = pos + 8
             variables = tmpString

           case (cgnsCl)
             tmpString = variables(1:pos)//'cl_'
             pos = pos + 3
             variables = tmpString

           case (cgnsClp)
             tmpString = variables(1:pos)//'clp_'
             pos = pos + 4
             variables = tmpString

           case (cgnsClv)
             tmpString = variables(1:pos)//'clv_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCd)
             tmpString = variables(1:pos)//'cd_'
             pos = pos + 3
             variables = tmpString

           case (cgnsCdp)
             tmpString = variables(1:pos)//'cdp_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCdv)
             tmpString = variables(1:pos)//'cdv_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCfx)
             tmpString = variables(1:pos)//'cfx_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCfy)
             tmpString = variables(1:pos)//'cfy_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCfz)
             tmpString = variables(1:pos)//'cfz_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCmx)
             tmpString = variables(1:pos)//'cmx_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCmy)
             tmpString = variables(1:pos)//'cmy_'
             pos = pos + 4
             variables = tmpString

           case (cgnsCmz)
             tmpString = variables(1:pos)//'cmz_'
             pos = pos + 4
             variables = tmpString

           case (cgnsHdiffMax)
             tmpString = variables(1:pos)//'hdiff_'
             pos = pos + 6
             variables = tmpString

           case (cgnsMachMax)
             tmpString = variables(1:pos)//'mach_'
             pos = pos + 5
             variables = tmpString

           case (cgnsYplusMax)
             tmpString = variables(1:pos)//'yplus_'
             pos = pos + 6
             variables = tmpString

           case (cgnsEddyMax)
             tmpString = variables(1:pos)//'eddyv_'
             pos = pos + 6
             variables = tmpString

         end select
       enddo

       ! Remove the trailing '_'.

       tmpString = variables(1:pos-1)
       variables = tmpString

       end subroutine setMonitorVariablesString

!      ==================================================================

       subroutine setSurfaceVariablesString(variables)
!
!      ******************************************************************
!      *                                                                *
!      * setSurfaceVariablesString writes the names of the surface      *
!      * variables to be written onto the output string.                *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use extraOutput
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(out) :: variables
!
!      Local variables.
!
       integer :: pos
       character(len=maxCplNameLen) :: tmpString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       pos = 0

       ! A chain of if statements to write the info onto the string
       ! variables.

       if(surfWriteRho) then
         tmpString = variables(1:pos)//'rho_'
         pos = pos + 4
         variables = tmpString
       endif

       if(surfWriteP) then
         tmpString = variables(1:pos)//'p_'
         pos = pos + 2
         variables = tmpString
       endif

       if(surfWriteTemp) then
         tmpString = variables(1:pos)//'temp_'
         pos = pos + 5
         variables = tmpString
       endif

       if(surfWriteVx) then
         tmpString = variables(1:pos)//'vx_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWriteVy) then
         tmpString = variables(1:pos)//'vy_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWriteVz) then
         tmpString = variables(1:pos)//'vz_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWriteCp) then
         tmpString = variables(1:pos)//'cp_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWritePtotloss) then
         tmpString = variables(1:pos)//'ptloss_'
         pos = pos + 7
         variables = tmpString
       endif

       if(surfWriteMach) then
         tmpString = variables(1:pos)//'mach_'
         pos = pos + 5
         variables = tmpString
       endif

       if(surfWriteCf) then
         tmpString = variables(1:pos)//'cf_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWriteCh) then
         tmpString = variables(1:pos)//'ch_'
         pos = pos + 3
         variables = tmpString
       endif

       if(surfWriteYplus) then
         tmpString = variables(1:pos)//'yplus_'
         pos = pos + 6
         variables = tmpString
       endif

       if(surfWriteCfx) then
         tmpString = variables(1:pos)//'cfx_'
         pos = pos + 4
         variables = tmpString
       endif

       if(surfWriteCfy) then
         tmpString = variables(1:pos)//'cfy_'
         pos = pos + 4
         variables = tmpString
       endif

       if(surfWriteCfz) then
         tmpString = variables(1:pos)//'cfz_'
         pos = pos + 4
         variables = tmpString
       endif

       if(surfWriteBlank) then
         tmpString = variables(1:pos)//'blank_'
         pos = pos + 5
         variables = tmpString
       endif

       ! Remove the trailing '_'.

       tmpString = variables(1:pos-1)
       variables = tmpString

       end subroutine setSurfaceVariablesString

!      ==================================================================

       subroutine setVolumeVariablesString(variables)
!
!      ******************************************************************
!      *                                                                *
!      * setVolumeVariablesString writes the names of the volume        *
!      * variables to be written onto the output string.                *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use extraOutput
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(out) :: variables
!
!      Local variables.
!
       integer :: pos
       character(len=maxCplNameLen) :: tmpString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       pos = 0

       ! A chain of if statements to write the info onto the string
       ! variables.

       if(volWriteMx) then
         tmpString = variables(1:pos)//'mx_'
         pos = pos + 3
         variables = tmpString
       endif

       if(volWriteMy) then
         tmpString = variables(1:pos)//'my_'
         pos = pos + 3
         variables = tmpString
       endif

       if(volWriteMz) then
         tmpString = variables(1:pos)//'mz_'
         pos = pos + 3
         variables = tmpString
       endif

       if(volWriteRhoe) then
         tmpString = variables(1:pos)//'rhoe_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteTemp) then
         tmpString = variables(1:pos)//'temp_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteVort) then
         tmpString = variables(1:pos)//'vort_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteVortx) then
         tmpString = variables(1:pos)//'vortx_'
         pos = pos + 6
         variables = tmpString
       endif

       if(volWriteVorty) then
         tmpString = variables(1:pos)//'vorty_'
         pos = pos + 6
         variables = tmpString
       endif

       if(volWriteVortz) then
         tmpString = variables(1:pos)//'vortz_'
         pos = pos + 6
         variables = tmpString
       endif

       if(volWriteCp) then
         tmpString = variables(1:pos)//'cp_'
         pos = pos + 3
         variables = tmpString
       endif

       if(volWriteMach) then
         tmpString = variables(1:pos)//'mach_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteMachTurb) then
         tmpString = variables(1:pos)//'macht_'
         pos = pos + 6
         variables = tmpString
       endif

       if(volWritePtotloss) then
         tmpString = variables(1:pos)//'ptloss_'
         pos = pos + 7
         variables = tmpString
       endif

       if(volWriteEddyVis) then
         tmpString = variables(1:pos)//'eddy_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteRatioEddyVis) then
         tmpString = variables(1:pos)//'eddyratio_'
         pos = pos + 10
         variables = tmpString
       endif

       if(volWriteDist) then
         tmpString = variables(1:pos)//'dist_'
         pos = pos + 5
         variables = tmpString
       endif

       if(volWriteResRho) then
         tmpString = variables(1:pos)//'resrho_'
         pos = pos + 7
         variables = tmpString
       endif

       if(volWriteResMom) then
         tmpString = variables(1:pos)//'resmom_'
         pos = pos + 7
         variables = tmpString
       endif

       if(volWriteResRhoe) then
         tmpString = variables(1:pos)//'resrhoe_'
         pos = pos + 8
         variables = tmpString
       endif

       if(volWriteResTurb) then
         tmpString = variables(1:pos)//'resturb_'
         pos = pos + 8
         variables = tmpString
       endif

       if(volWriteBlank) then
         tmpString = variables(1:pos)//'blank_'
         pos = pos + 6
         variables = tmpString
       endif

       ! Remove the trailing '_'.

       tmpString = variables(1:pos-1)
       variables = tmpString
 
       end subroutine setVolumeVariablesString

!      ==================================================================

       subroutine sumb_getParamLog(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_getParamLog sets the logical parameter with the given     *
!      * name to the given value.                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(out) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("vortex correction")
           value = vortexCorr

         case ("directional dissipation scaling")
           value = dirScaling

         case ("total enthalpy scaling inlet")
           value = hScalingInlet
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("restart")
           value = restart

         case ("check nondimensionalization")
           value = checkRestartSol

         case ("rind layer in solution files")
           value = storeRindLayer

         case ("automatic parameter update")
           value = autoParameterUpdate

         case ("write coordinates in meter")
           value = writeCoorMeter

         case ("store convergence inner iterations")
           value = storeConvInnerIter
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("freeze turbulent source terms in mg")
           value = freezeTurbSource
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("split blocks for load balance")
           value = splitBlocks
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("v2f with upper bound")
           value = rvfB

         case ("use wall functions")
           value = wallFunctions
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("write file for unsteady restart")
           value = writeUnsteadyRestartSpectral

         case ("write unsteady volume solution files")
           value = writeUnsteadyVolSpectral

         case ("write unsteady surface solution files")
           value = writeUnsteadySurfSpectral
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("update wall distance unsteady mode")
           value = updateWallDistanceUnsteady
!
!        ****************************************************************
!        *                                                              *
!        * Visualization parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("pv3 visualization only")
           value = PV3VisOnly
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("get coarse-level sol")
           value = cplGetCoarseSol
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("input overset donors are guesses")
           value = oversetDonorsAreGuesses

         case ("average restricted residual for blanks")
           value = avgRestrictResForBlanks
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_getParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_getParamLog

!      ==================================================================

       subroutine sumb_setParamInt(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamInt sets the integer parameter with the given     *
!      * name to the given value.                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos, ierr
       integer(kind=intType) :: nn

       character (len=maxStringLen)   :: keyword
       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("discretization scheme")
           if(value /= dissScalar .and. value /= dissMatrix .and. &
              value /= dissCusp   .and. value /= upwind) then
             write(errorMessage,*) "Unknown ", trim(keyword), &
                                   ", ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           spaceDiscr = value

         case ("discretization scheme coarse grid")
           if(value /= dissScalar .and. value /= dissMatrix .and. &
              value /= dissCusp   .and. value /= upwind) then
             write(errorMessage,*) "Unknown ", trim(keyword), &
                                   ", ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           spaceDiscrCoarse = value

         case ("order turbulent equations")
           if(value /= firstOrder .and. value /= secondOrder) then
             write(errorMessage,*) "Order turbulent equations must &
                                    &be first order or second order, &
                                    &not ", value
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           orderTurb = value

         case ("riemann solver")
           if(value /= Roe .and. value /= vanLeer .and. &
              value /= ausmdv) then
             write(errorMessage,*) "Unknown ", trim(keyword), &
                                   ", ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           riemann = value

         case ("riemann solver coarse grid")
           if(value /= Roe .and. value /= vanLeer .and. &
              value /= ausmdv) then
             write(errorMessage,*) "Unknown ", trim(keyword), &
                                   ", ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           riemannCoarse = value

         case ("limiter")
           if(value /= firstOrder .and. value /= noLimiter .and. &
              value /= vanAlbeda  .and. value /= minmod) then
             write(errorMessage,*) "Unknown limiter, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           limiter = value

         case ("preconditioner")
           if(value /= noPrecond .and. value /= Turkel .and. &
              value /= ChoiMerkle) then
             write(errorMessage,*) "Unknown preconditioner, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           precond = value

         case ("wall boundary treatment")
           if(value /= constantPressure     .and. &
              value /= linExtrapolPressure  .and. &
              value /= quadExtrapolPressure .and. &
              value /= normalMomentum) then
             write(errorMessage,*) "Unknown wall boundary &
                                   &treatment, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           wallBCTreatment = value

         case ("outflow boundary treatment")
           if(value /= constantExtrapol .and. &
              value /= linExtrapol) then
             write(errorMessage,*) "Unknown outflow boundary &
                                   &treatment, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           outflowTreatment = value

         case ("non-matching block to block treatment")
           if(value /= NonConservative .and. &
              value /= Conservative) then
             write(errorMessage,*) "Unknown non-matching block to &
                                   &block treatment, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           nonMatchTreatment = value
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("file format read")
           if(value /= cgnsFormat .and. value /= plot3DFormat) then
             write(errorMessage,*) "Unknown file format, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           fileFormatRead = value

         case ("file format write")
           if(value /= cgnsFormat .and. value /= plot3DFormat) then
             write(errorMessage,*) "Unknown file format, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           fileFormatWrite = value

         case ("write precision grid")
           if(value /= precisionSingle .and. &
              value /= precisionDouble) then
             write(errorMessage,*) "Unknown write precision &
                                   &grid, ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           precisionGrid = value

         case ("write precision solution")
           if(value /= precisionSingle .and. &
              value /= precisionDouble) then
             write(errorMessage,*) "Unknown write precision &
                                   &solution, ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           precisionSol = value
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("smoother")
           if(value /= RungeKutta .and. value /= nlLusgs .and. &
              value /= nlLusgsLine) then
             write(errorMessage,*) "Unknown smoother, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           smoother = value

         case ("treatment turbulent equations")
           if(value /= segregated .and. value /= coupled) then
             write(errorMessage,*) "Unknown treatment turbulent &
                                   &equations, ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           turbTreatment = value

         case ("turbulent smoother")
           if(value /= gmres .and. value /= adi) then
             write(errorMessage,*) "Unknown turbulent smoother, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           turbSmoother = value

         case ("turbulent relaxation")
           if(value /= turbRelaxExplicit .and. &
              value /= turbRelaxImplicit) then
             write(errorMessage,*) "Unknown turbulent relaxation, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           turbRelax = value

         case ("residual averaging")
           if(value /= noResAveraging     .and. &
              value /= alwaysResAveraging .and. &
              value /= alternateResAveraging) then
             write(errorMessage,*) "Unknown residual averaging, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           resAveraging = value

         case ("treatment boundary multigrid corrections")
           if(value /= bcDirichlet0 .and. &
              value /= bcNeumann) then
             write(errorMessage,*) "Unknown treatment boundary &
                                   &multigrid corrections, ",  &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           mgBoundCorr = value

         case ("number of multigrid cycles")
           nCycles = value

         case ("number of single grid startup iterations")
           nsgStartup = value

         case ("save every")
           nSaveVolume = value

         case ("save surface every")
           nSaveSurface = value

         case ("number of runge kutta stages")
           nRKStages = value

         case ("number of multigrid cycles coarse grid")
           nCyclesCoarse = value

         case ("multigrid start level")
           mgStartlevel = value
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Polynomial rotation parameters.

         case ("degree polynomial x-rotation")
           degreePolXRot = value

           nn = max(degreePolXRot,0_intType)
           if(allocated(coefPolXRot)) then
             deallocate(coefPolXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolXRot")
           endif
           allocate(coefPolXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolXRot")

         case ("degree polynomial y-rotation")
           degreePolYRot = value

           nn = max(degreePolYRot,0_intType)
           if(allocated(coefPolYRot)) then
             deallocate(coefPolYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolYRot")
           endif
           allocate(coefPolYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolYRot")

         case ("degree polynomial z-rotation")
           degreePolZRot = value

           nn = max(degreePolZRot,0_intType)
           if(allocated(coefPolZRot)) then
             deallocate(coefPolZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolZRot")
           endif
           allocate(coefPolZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolZRot")

         ! Fourier rotation parameters.

         case ("degree fourier x-rotation")
           degreeFourXRot = value

           nn = max(degreeFourXRot,0_intType)
           if(allocated(cosCoefFourXRot)) then
             deallocate(cosCoefFourXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourXRot")
           endif
           allocate(cosCoefFourXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourXRot")

           nn = max(degreeFourXRot,1_intType)
           if(allocated(sinCoefFourXRot)) then
             deallocate(sinCoefFourXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourXRot")
           endif
           allocate(sinCoefFourXRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourXRot")

         case ("degree fourier y-rotation")
           degreeFourYRot = value

           nn = max(degreeFourYRot,0_intType)
           if(allocated(cosCoefFourYRot)) then
             deallocate(cosCoefFourYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourYRot")
           endif
           allocate(cosCoefFourYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourYRot")

           nn = max(degreeFourYRot,1_intType)
           if(allocated(sinCoefFourYRot)) then
             deallocate(sinCoefFourYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourYRot")
           endif
           allocate(sinCoefFourYRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourYRot")

         case ("degree fourier z-rotation")
           degreeFourZRot = value

           nn = max(degreeFourZRot,0_intType)
           if(allocated(cosCoefFourZRot)) then
             deallocate(cosCoefFourZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourZRot")
           endif
           allocate(cosCoefFourZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourZRot")

           nn = max(degreeFourZRot,1_intType)
           if(allocated(sinCoefFourZRot)) then
             deallocate(sinCoefFourZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourZRot")
           endif
           allocate(sinCoefFourZRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourZRot")
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("equations")
           if(value /= EulerEquations .and. &
              value /= NSEquations    .and. &
              value /= RANSEquations) then
             write(errorMessage,*) "Unknown equations, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           equations = value

         case ("mode")
           if(value /= steady .and. value /= unsteady .and. &
              value /= timeSpectral) then
             write(errorMessage,*) "Unknown mode, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           equationMode = value

         case ("flow type")
           if(value /= internalFlow .and. &
              value /= externalFlow .and. &
              value /= timeSpectral) then
             write(errorMessage,*) "Unknown flow type, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           flowType = value

         case ("cp model")
           if(value /= cpConstant .and. value /= cpTempCurveFits) then
             write(errorMessage,*) "Unknown Cp model, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           cpModel = value

         case ("turbulence model")
           if(value /= baldwinLomax           .and. &
              value /= spalartAllmaras        .and. &
              value /= spalartAllmarasEdwards .and. &
              value /= komegaWilcox           .and. &
              value /= komegaModified         .and. &
              value /= ktau                   .and. &
              value /= menterSST              .and. &
              value /= v2f) then
             write(errorMessage,*) "Unknown turbulence model, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           turbModel = value

         case ("v2f version (n1 or n6)")
           if(value /= 1 .and. value /= 6) then
             write(errorMessage,*) "v2f version must be either &
                                   &1 or 6, not ", value
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           rvfN = value

         case ("turbulence production term")
           if(value /= strain           .and. &
              value /= vorticity        .and. &
              value /= katoLaunder .and. &
              value /= komegaWilcox           .and. &
              value /= komegaModified         .and. &
              value /= ktau                   .and. &
              value /= menterSST              .and. &
              value /= v2f) then
             write(errorMessage,*) "Unknown turbulence production &
                                   &term, ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           turbProd = value
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("number time intervals spectral")
           nTimeIntervalsSpectral = value

         case ("number of unsteady solution files")
           nUnsteadySolSpectral = value
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("time accuracy unsteady")
           if(value /= firstOrder .and. value /= secondOrder .and. &
              value /= thirdOrder) then
             write(errorMessage,*) "Unknown time accuracy unsteady, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           timeAccuracy = value

         case ("number of unsteady time steps coarse grid")
           nTimeStepsCoarse = value

         case ("number of unsteady time steps fine grid")
           nTimeStepsFine = value
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("overset interpolation type")
           if(value /= trilinear) then
             write(errorMessage,*) "Unknown overset interpolation &
                                   &type, ", value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           oversetInterpType = value

         case ("overset interpolation type coarse grid")
           if(value /= trilinear) then
             write(errorMessage,*) "Unknown overset interpolation &
                                   &type coarse grid, ", &
                                    value, ", specified"
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

           oversetInterpTypeCoarse = value

         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_setParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamInt

!      ==================================================================

       subroutine sumb_setParamIntRk1(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamIntRk1 sets the parameter defined as the rank-1   *
!      * integer array with the given name to the given value.          *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), dimension(:), intent(in) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamIntRk1

!      ==================================================================

       subroutine sumb_setParamReal(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamReal sets the real parameter with the given       *
!      * name to the given value.                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(in) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos, p1, p2
       integer(kind=intType) :: iRotDeg

       character (len=maxStringLen)   :: keyword, tmpName
       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)

       if(keyword(1:34) == "polynomial coefficients x-rotation" .or. &
          keyword(1:34) == "polynomial coefficients y-rotation" .or. &
          keyword(1:34) == "polynomial coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for coefficients of &
                                 &polynomial rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_setParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:34)
         keyword = tmpName

       else if(keyword(1:38) == "fourier cosine coefficients x-rotation" .or. &
               keyword(1:38) == "fourier cosine coefficients y-rotation" .or. &
               keyword(1:38) == "fourier cosine coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for Fourier cosine &
                                 &coefficients of rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_setParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:38)
         keyword = tmpName

       else if(keyword(1:36) == "fourier sine coefficients x-rotation" .or. &
               keyword(1:36) == "fourier sine coefficients y-rotation" .or. &
               keyword(1:36) == "fourier sine coefficients z-rotation") then
         p1 = index(keyword, "(")
         p2 = index(keyword, ")")
         if(p2 <= p1+1) then
           write(errorMessage,*) "Degree indices for Fourier sine &
                                 &coefficients of rotation wrongly specified"
           if(myID == 0) &
             call terminate("sumb_setParam", errorMessage)
           call mpi_barrier(SUmb_comm_world, pos)
         endif
         tmpName = keyword(p1+1:p2-1)
         read(tmpName,*) iRotDeg
         tmpName = keyword(1:36)
         keyword = tmpName
       endif
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("vis2")
           vis2 = value

         case ("vis4")
           vis4 = value

         case ("vis2 coarse grid")
           vis2Coarse = value

         case ("exponent dissipation scaling")
           adis = abs(value)

         case ("kappa interpolation value")
           kappaCoef = value
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("residual averaging smoothing parameter")
           smoop = value

         case ("cfl number")
           cfl = value

         case ("alpha turbulent dd-adi")
           alfaTurb = max(1.e-10_realType,min(value,0.99_realType))

         case ("beta turbulent dd-adi")
           betaTurb = max(1.e-10_realType,min(value,0.99_realType))

         case ("relative l2 norm for convergence")
           L2Conv = value

         case ("cfl number coarse grid")
           cflCoarse = value

         case ("relative l2 norm for convergence coarse grid")
           L2ConvCoarse = value

         case ("restriction relaxation factor")
           fcoll = min(value,one)
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Rotation point.

         case ("rotation point body x")
           rotPoint(1) = value

         case ("rotation point body y")
           rotPoint(2) = value

         case ("rotation point body z")
           rotPoint(3) = value

         ! Polynomial rotation parameters.
         ! It is assumed that polynomial degrees are already reset
         ! and the arrays for coefficients are already reallocated.

         case ("polynomial coefficients x-rotation")
           coefPolXRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("polynomial coefficients y-rotation")
           coefPolYRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("polynomial coefficients z-rotation")
           coefPolZRot(iRotDeg) = value
           gridMotionSpecified = .true.

         ! Fourier rotation parameters.
         ! It is assumed that Fourier degrees are already reset
         ! and the arrays for coefficients are already reallocated.

         case ("omega fourier x-rotation")
           omegaFourXRot = value

         case ("omega fourier y-rotation")
           omegaFourYRot = value

         case ("omega fourier z-rotation")
           omegaFourZRot = value

         case ("fourier cosine coefficients x-rotation")
           cosCoefFourXRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("fourier sine coefficients x-rotation")
           sinCoefFourXRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("fourier cosine coefficients y-rotation")
           cosCoefFourYRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("fourier sine coefficients y-rotation")
           sinCoefFourYRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("fourier cosine coefficients z-rotation")
           cosCoefFourZRot(iRotDeg) = value
           gridMotionSpecified = .true.

         case ("fourier sine coefficients z-rotation")
           sinCoefFourZRot(iRotDeg) = value
           gridMotionSpecified = .true.
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable load imbalance")
           loadImbalance = value
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("offset from wall in wall functions")
           wallOffset = max(zero, value)

         case ("max ratio k-prod/dest")
           pklim = value
           if(pklim <= zero) pklim = 20.0_realType

         case ("mach")
           Mach = value

         case ("mach for coefficients")
           MachCoef = value

         case ("reynolds")
           Reynolds = value

         case ("free stream velocity direction x")
           velDirFreestream(1) = value

         case ("free stream velocity direction y")
           velDirFreestream(2) = value

         case ("free stream velocity direction z")
           velDirFreestream(3) = value

         case ("lift direction x")
           liftDirection(1) = value
         ! liftDirSpecified = .true.

         case ("lift direction y")
           liftDirection(2) = value
         ! liftDirSpecified = .true.

         case ("lift direction z")
           liftDirection(3) = value
         ! liftDirSpecified = .true.

         case ("reynolds length (in meter)")
           ReynoldsLength = value

         case ("free stream temperature (in k)")
           tempFreestream = value

         case ("constant specific heat ratio")
           gammaConstant = value

         case ("gas constant (j/(kg k))")
           RGasDim = value

         case ("prandtl number")
           prandtl = value

         case ("turbulent prandtl number")
           prandtlTurb = value

         case ("free stream eddy viscosity ratio")
           eddyVisInfRatio = value

         case ("free stream turbulent intensity")
           turbIntensityInf = value
           if(turbIntensityInf < 0.0) turbIntensityInf = 0.001_realType

         case ("reference surface")
           surfaceRef = value

         case ("reference length")
           lengthRef = value

         case ("moment reference point x")
           pointRef(1) = value

         case ("moment reference point y")
           pointRef(2) = value

         case ("moment reference point z")
           pointRef(3) = value
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("time step (in sec) for unsteady restart")
           dtUnsteadyRestartSpectral = value
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("unsteady time step (in sec)")
           deltaT = value
!
!        ****************************************************************
!        *                                                              *
!        * Reference state values.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("reference pressure (in pa)")
           pRef = value

         case ("reference density (in kg/m^3)")
           rhoRef = value

         case ("reference temperature (in k)")
           TRef = value

         case ("conversion factor grid units to meter")
           LRef = value
           LRefSpecified = .true.
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("mach for initialization")
           MachIni = value

         case ("pressure for initialization")
           pIni = value

         case ("density for initialization")
           rhoIni = value

         case ("velocity direction x for initialization")
           velDirIni(1) = value

         case ("velocity direction y for initialization")
           velDirIni(2) = value

         case ("velocity direction z for initialization")
           velDirIni(3) = value
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable donor quality")
           allowableDonorQuality = value
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_setParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamReal

!      ==================================================================

       subroutine sumb_setParamRealRk1(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamRealRk1 sets the parameter defined as the rank-1  *
!      * real array with the given name to the given value.             *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:), intent(in) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("rotation point body (x,y,z)")
           rotPoint(1) = value(1)
           rotPoint(2) = value(2)
           rotPoint(3) = value(3)
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("free stream velocity direction")
           velDirFreestream(1) = value(1)
           velDirFreestream(2) = value(2)
           velDirFreestream(3) = value(3)

         case ("lift direction")
           liftDirection(1) = value(1)
           liftDirection(2) = value(2)
           liftDirection(3) = value(3)
         ! liftDirSpecified = .true.

         case ("moment reference point (x,y,z)")
           pointRef(1) = value(1)
           pointRef(2) = value(2)
           pointRef(3) = value(3)
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("velocity direction for initialization")
           velDirIni(1) = value(1)
           velDirIni(2) = value(2)
           velDirIni(3) = value(3)
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamRealRk1

!      ==================================================================

       subroutine sumb_setParamStr(valueIn, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamStr sets the string parameter with the given      *
!      * name to the given value. Note that sumb_setParamStr can        *
!      * handle any parameter once it is transferred to SUmb in the     *
!      * string format.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       character (len=*), intent(in) :: valueIn
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos, ierr
       integer(kind=intType) :: nn

       character (len=maxStringLen)   :: keyword, value
       character (len=2*maxStringLen) :: errorMessage
!
!      Function definition.
!
       logical               :: checkYesNo
       integer(kind=intType) :: determineDiscretization
       integer(kind=intType) :: determineFileFormat
       integer(kind=intType) :: determineRiemann
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       value = valueIn
       call replaceTabsAndReturns(value)
       value = adjustl(value)
       value = trim(value)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0 .or. len_trim(value) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * The parameters to monitor the convergence, the surface and   *
!        * extra volume variables to written to the solution files.     *
!        *                                                              *
!        ****************************************************************
!
         case ("monitoring variables")
           call monitorVariables(value)

         case ("surface output variables")
           call surfaceVariables(value)

         case ("volume output variables")
           call volumeVariables(value)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("discretization scheme")
           spaceDiscr = determineDiscretization(value, keyword)

         case ("discretization scheme coarse grid")
           spaceDiscrCoarse = determineDiscretization(value, keyword)

         case ("order turbulent equations")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("first order")
               orderTurb = firstOrder
             case ("second order")
               orderTurb = secondOrder
             case default
               write(errorMessage,*) "Order turbulent equations must &
                                      &be first order or second order, &
                                      &not ", trim(value)
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("riemann solver")
           riemann = determineRiemann(value, keyword)

         case ("riemann solver coarse grid")
           riemannCoarse = determineRiemann(value, keyword)

         case ("limiter")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("first order")
               limiter = firstOrder
             case ("no limiter")
               limiter = noLimiter
             case ("van albeda")
               limiter = vanAlbeda
             case ("minmod")
               limiter = minmod
             case default
               write(errorMessage,*) "Unknown limiter, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("preconditioner")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("no preconditioner")
               precond = noPrecond
             case ("turkel")
               precond = Turkel
             case ("choi merkle")
               precond = ChoiMerkle
             case default
               write(errorMessage,*) "Unknown preconditioner, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("wall boundary treatment")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("constant pressure")
               wallBcTreatment = constantPressure
             case ("linear extrapolation pressure")
               wallBcTreatment = linExtrapolPressure
             case ("quadratic extrapolation pressure")
               wallBcTreatment = quadExtrapolPressure
             case ("normal momentum")
               wallBcTreatment = normalMomentum
             case default
               write(errorMessage,*) "Unknown wall boundary &
                                      &treatment, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("outflow boundary treatment")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("constant extrapolation")
               outflowTreatment = constantExtrapol
             case ("linear extrapolation")
               outflowTreatment = linExtrapol
             case default
               write(errorMessage,*) "Unknown outflow boundary &
                                      &treatment, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("non-matching block to block treatment")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("nonconservative")
               nonMatchTreatment = NonConservative
             case ("conservative")
               nonMatchTreatment = Conservative
             case default
               write(errorMessage,*) "Unknown non-matching block to &
                                      &block treatment, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("vortex correction")
           vortexCorr = checkYesNo(value, keyword)

         case ("vis2")
           read(value,*) vis2

         case ("vis4")
           read(value,*) vis4

         case ("vis2 coarse grid")
           read(value,*) vis2Coarse

         case ("directional dissipation scaling")
           dirScaling = checkYesNo(value, keyword)

         case ("exponent dissipation scaling")
           read(value,*) adis
           adis = abs(adis)

         case ("total enthalpy scaling inlet")
           hScalingInlet = checkYesNo(value, keyword)

         case ("kappa interpolation value")
           read(value,*) kappaCoef
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("file format read")

           fileFormatRead = determineFileFormat(value, keyword)

         case ("file format write")
           fileFormatWrite = determineFileFormat(value, keyword)

         case ("grid file")
           gridFile = value

         case ("plot3d connectivity file")
           plot3DConnFile = value

         case ("restart file")
           restartFile = value

         case ("restart")
           restart = checkYesNo(value, keyword)

         case ("check nondimensionalization")
           checkRestartSol = checkYesNo(value, keyword)

         case ("new grid file")
            newGridFile = value

         case ("solution file")
           solFile = value

         case ("surface solution file")
           surfaceSolFile = value

         case ("write precision grid")
           select case (value)
             case ("single")
               precisionGrid = precisionSingle
             case ("double")
               precisionGrid = precisionDouble
             case default
               write(errorMessage,*) "Unknown write precision &
                                     &grid, ", trim(value), &
                                     ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("write precision solution")
           select case (value)
             case ("single")
               precisionSol = precisionSingle
             case ("double")
               precisionSol = precisionDouble
             case default
               write(errorMessage,*) "Unknown write precision &
                                     &solution, ", trim(value), &
                                     ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("rind layer in solution files")
           storeRindLayer = checkYesNo(value, keyword)

         case ("automatic parameter update")
           autoParameterUpdate = checkYesNo(value, keyword)

         case ("write coordinates in meter")
           writeCoorMeter = checkYesNo(value, keyword)

         case ("cp curve fit file")
           cpFile = value

         case ("store convergence inner iterations")
           storeConvInnerIter = checkYesNo(value, keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("smoother")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("runge kutta")
               smoother = RungeKutta
             case ("nonlinear lusgs")
               smoother = nlLusgs
             case ("nonlinear lusgs line")
               smoother = nlLusgsLine
             case default
               write(errorMessage,*) "Unknown smoother, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("treatment turbulent equations")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("segregated")
               turbTreatment = segregated
             case ("coupled")
               turbTreatment = coupled
             case default
               write(errorMessage,*) "Unknown treatment turbulent &
                                      &equations, ", trim(value),  &
                                      ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("turbulent smoother")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("gmres")
               turbSmoother = gmres
             case ("adi")
               turbSmoother = adi
             case default
               write(errorMessage,*) "Unknown turbulent smoother, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("freeze turbulent source terms in mg")
           freezeTurbSource = checkYesNo(value, keyword)

         case ("turbulent relaxation")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("explicit")
               turbRelax = turbRelaxExplicit
             case ("implicit")
               turbRelax = turbRelaxImplicit
             case default
               write(errorMessage,*) "Unknown turbulent relaxation, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("residual averaging")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("no")
               resAveraging = noResAveraging
             case ("all stages")
               resAveraging = alwaysResAveraging
             case ("alternate stages")
               resAveraging = alternateResAveraging
             case default
               write(errorMessage,*) "Unknown residual averaging, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("residual averaging smoothing parameter")
           read(value,*) smoop

         case ("number of multigrid cycles")
           read(value,*) nCycles

         case ("number of single grid startup iterations")
           read(value,*) nsgStartup

         case ("save every")
           read(value,*) nSaveVolume

         case ("save surface every")
           read(value,*) nSaveSurface

         case ("number of runge kutta stages")
           read(value,*) nRKStages

         case ("cfl number")
           read(value,*) cfl

         case ("alpha turbulent dd-adi")
           read(value,*) alfaTurb
           alfaTurb = max(1.e-10_realType,min(alfaTurb,0.99_realType))

         case ("beta turbulent dd-adi")
           read(value,*) betaTurb
           betaTurb = max(1.e-10_realType,min(betaTurb,0.99_realType))

         case ("relative l2 norm for convergence")
           read(value,*) L2Conv

         case ("number of multigrid cycles coarse grid")
           read(value,*) nCyclesCoarse

         case ("cfl number coarse grid")
           read(value,*) cflCoarse

         case ("relative l2 norm for convergence coarse grid")
           read(value,*) L2ConvCoarse

         case ("treatment boundary multigrid corrections")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("zero dirichlet")
               mgBoundCorr = bcDirichlet0
             case ("neumann")
               mgBoundCorr = bcNeumann
             case default
               write(errorMessage,*) "Unknown treatment boundary &
                                      &multigrid corrections, ",  &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("restriction relaxation factor")
           read(value,*) fcoll
           fcoll = min(fcoll,one)

         case ("multigrid start level")
           read(value,*) mgStartlevel

       ! case ("multigrid cycle strategy")
       !   mgDescription = value
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Rotation point.

         case ("rotation point body (x,y,z)")
           read(value,*) rotPoint(1), rotPoint(2), rotPoint(3)

         ! Polynomial rotation parameters.

         case ("degree polynomial x-rotation")
           read(value,*) degreePolXRot

         case ("degree polynomial y-rotation")
           read(value,*) degreePolYRot

         case ("degree polynomial z-rotation")
           read(value,*) degreePolZRot

         case ("polynomial coefficients x-rotation")
           nn = max(degreePolXRot,0_intType)
           if(allocated(coefPolXRot)) then
             deallocate(coefPolXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolXRot")
           endif
           allocate(coefPolXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolXRot")

           call readMotionCoef(value, 0_intType, degreePolXRot, &
                               coefPolXRot)

           gridMotionSpecified = .true.

         case ("polynomial coefficients y-rotation")
           nn = max(degreePolYRot,0_intType)
           if(allocated(coefPolYRot)) then
             deallocate(coefPolYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolYRot")
           endif
           allocate(coefPolYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolYRot")

           call readMotionCoef(value, 0_intType, degreePolYRot, &
                               coefPolYRot)

           gridMotionSpecified = .true.

         case ("polynomial coefficients z-rotation")
           nn = max(degreePolZRot,0_intType)
           if(allocated(coefPolZRot)) then
             deallocate(coefPolZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for coefPolZRot")
           endif
           allocate(coefPolZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &coefPolZRot")

           call readMotionCoef(value, 0_intType, degreePolZRot, &
                               coefPolZRot)

           gridMotionSpecified = .true.

         ! Fourier rotation parameters.

         case ("degree fourier x-rotation")
           read(value,*) degreeFourXRot

         case ("degree fourier y-rotation")
           read(value,*) degreeFourYRot

         case ("degree fourier z-rotation")
           read(value,*) degreeFourZRot

         case ("omega fourier x-rotation")
           read(value,*) omegaFourXRot

         case ("omega fourier y-rotation")
           read(value,*) omegaFourYRot

         case ("omega fourier z-rotation")
           read(value,*) omegaFourZRot

         case ("fourier cosine coefficients x-rotation")
           nn = max(degreeFourXRot,0_intType)
           if(allocated(cosCoefFourXRot)) then
             deallocate(cosCoefFourXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourXRot")
           endif
           allocate(cosCoefFourXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourXRot")

           call readMotionCoef(value, 0_intType, degreeFourXRot, &
                               cosCoefFourXRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients x-rotation")
           nn = max(degreeFourXRot,1_intType)
           if(allocated(sinCoefFourXRot)) then
             deallocate(sinCoefFourXRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourXRot")
           endif
           allocate(sinCoefFourXRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourXRot")

           call readMotionCoef(value, 1_intType, degreeFourXRot, &
                               sinCoefFourXRot)

           gridMotionSpecified = .true.

         case ("fourier cosine coefficients y-rotation")
           nn = max(degreeFourYRot,0_intType)
           if(allocated(cosCoefFourYRot)) then
             deallocate(cosCoefFourYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourYRot")
           endif
           allocate(cosCoefFourYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourYRot")

           call readMotionCoef(value, 0_intType, degreeFourYRot, &
                               cosCoefFourYRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients y-rotation")
           nn = max(degreeFourYRot,1_intType)
           if(allocated(sinCoefFourYRot)) then
             deallocate(sinCoefFourYRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourYRot")
           endif
           allocate(sinCoefFourYRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourYRot")

           call readMotionCoef(value, 1_intType, degreeFourYRot, &
                               sinCoefFourYRot)

           gridMotionSpecified = .true.

         case ("fourier cosine coefficients z-rotation")
           nn = max(degreeFourZRot,0_intType)
           if(allocated(cosCoefFourZRot)) then
             deallocate(cosCoefFourZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for cosCoefFourZRot")
           endif
           allocate(cosCoefFourZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &cosCoefFourZRot")

           call readMotionCoef(value, 0_intType, degreeFourZRot, &
                               cosCoefFourZRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients z-rotation")
           nn = max(degreeFourZRot,1_intType)
           if(allocated(sinCoefFourZRot)) then
             deallocate(sinCoefFourZRot, stat=ierr)
             if(ierr /=0) call terminate("sumb_setParam", &
                                         "Memory deallocation failure &
                                         &for sinCoefFourZRot")
           endif
           allocate(sinCoefFourZRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("sumb_setParam", &
                                        "Memory allocation failure for &
                                        &sinCoefFourZRot")

           call readMotionCoef(value, 1_intType, degreeFourZRot, &
                               sinCoefFourZRot)

           gridMotionSpecified = .true.
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("allowable load imbalance")
           read(value,*) loadImbalance

         case ("split blocks for load balance")
           splitBlocks = checkYesNo(value, keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("equations")

           ! Convert value to lower case and check the three options.

           call convertToLowerCase(value)

           select case (value)
             case ("euler")
               equations = EulerEquations
             case ("laminar ns")
               equations = NSEquations
             case ("rans")
               equations = RANSEquations
             case default
               write(errorMessage,*) "Unknown equations, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("mode")

           ! Convert value to lower case and check the two options.

           call convertToLowerCase(value)

           select case (value)
             case ("steady")
               equationMode = steady
             case ("unsteady")
               equationMode = unsteady
             case ("time spectral")
               equationMode = timeSpectral
             case default
               write(errorMessage,*) "Unknown mode, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("flow type")

           ! Convert value to lower case and check the two options.

           call convertToLowerCase(value)

           select case (value)

             case ("internal flow")
               flowType = internalFlow
             case ("external flow")
               flowType = externalFlow
             case default
               write(errorMessage,*) "Unknown flow type, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("cp model")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)

             case ("constant")
               cpModel = cpConstant

             case ("temperature curve fits")
               cpModel = cpTempCurveFits

             case default
               write(errorMessage,*) "Unknown Cp model, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("turbulence model")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("baldwin lomax")
               turbModel = baldwinLomax
             case ("spalart allmaras")
               turbModel = spalartAllmaras
             case ("spalart allmaras edwards")
               turbModel = spalartAllmarasEdwards
             case ("komega wilcox")
               turbModel = komegaWilcox
             case ("komega modified")
               turbModel = komegaModified
             case ("ktau")
               turbModel = ktau
             case ("menter sst")
               turbModel = menterSST
             case ("v2f")
               turbModel = v2f
             case default
               write(errorMessage,*) "Unknown turbulence model, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("v2f version (n1 or n6)")
           read(value,*) rvfN
           if(rvfN /= 1 .and. rvfN /= 6) then
             write(errorMessage,*) "v2f version must be either &
                                    &1 or 6, not ", trim(value)
             if(myID == 0) &
               call terminate("sumb_setParam", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)
           endif

         case ("v2f with upper bound")
           rvfB = checkYesNo(value, keyword)

         case ("turbulence production term")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("strain")
               turbProd = strain
             case ("vorticity")
               turbProd = vorticity
             case ("kato-launder")
               turbProd = katoLaunder
             case default
               write(errorMessage,*) "Unknown turbulence production &
                                      &term, ", trim(value), &
                                      ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("use wall functions")
           wallFunctions = checkYesNo(value, keyword)

         case ("offset from wall in wall functions")
           read(value,*) wallOffset
           wallOffset = max(zero, wallOffset)

         case ("max ratio k-prod/dest")
           read(value,*) pklim
           if(pklim <= zero) pklim = 20.0_realType

         case ("mach")
           read(value,*) Mach

         case ("mach for coefficients")
           read(value,*) MachCoef

         case ("reynolds")
           read(value,*) Reynolds

         case ("free stream velocity direction")
           read(value,*) velDirFreestream(1), velDirFreestream(2), &
                         velDirFreestream(3)

         case ("lift direction")
           read(value,*) liftDirection(1), liftDirection(2), &
                         liftDirection(3)
         ! liftDirSpecified = .true.

         case ("reynolds length (in meter)")
           read(value,*) ReynoldsLength

         case ("free stream temperature (in k)")
           read(value,*) tempFreestream

         case ("constant specific heat ratio")
           read(value,*) gammaConstant

         case ("gas constant (j/(kg k))")
           read(value,*) RGasDim

         case ("prandtl number")
           read(value,*) prandtl

         case ("turbulent prandtl number")
           read(value,*) prandtlTurb

         case ("free stream eddy viscosity ratio")
           read(value,*) eddyVisInfRatio

         case ("free stream turbulent intensity")
           read(value,*) turbIntensityInf
           if(turbIntensityInf < 0.0) turbIntensityInf = 0.001_realType

         case ("reference surface")
           read(value,*) surfaceRef

         case ("reference length")
           read(value,*) lengthRef

         case ("moment reference point x")
           read(value,*) pointRef(1)

         case ("moment reference point y")
           read(value,*) pointRef(2)

         case ("moment reference point z")
           read(value,*) pointRef(3)

         case ("moment reference point (x,y,z)")
           read(value,*) pointRef(1), pointRef(2), pointRef(3)
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("number time intervals spectral")
           read(value,*) nTimeIntervalsSpectral

         case ("write file for unsteady restart")
           writeUnsteadyRestartSpectral = checkYesNo(value, keyword)

         case ("time step (in sec) for unsteady restart")
           read(value,*) dtUnsteadyRestartSpectral

         case ("write unsteady volume solution files")
           writeUnsteadyVolSpectral = checkYesNo(value, keyword)

         case ("write unsteady surface solution files")
           writeUnsteadySurfSpectral = checkYesNo(value, keyword)

         case ("number of unsteady solution files")
           read(value,*) nUnsteadySolSpectral
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("time accuracy unsteady")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("first")
               timeAccuracy = firstOrder
             case ("second")
               timeAccuracy = secondOrder
             case ("third")
               timeAccuracy = thirdOrder
             case default
               write(errorMessage,*) "Unknown time accuracy unsteady, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("number of unsteady time steps coarse grid")
           read(value,*) nTimeStepsCoarse

         case ("number of unsteady time steps fine grid")
           read(value,*) nTimeStepsFine

         case ("unsteady time step (in sec)")
           read(value,*) deltaT

         case ("update wall distance unsteady mode")
           updateWallDistanceUnsteady = checkYesNo(value, keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Visualization parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("pv3 visualization only")
           PV3VisOnly = checkYesNo(value, keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Reference state values.                                      *
!        *                                                              *
!        ****************************************************************
!
         case ("reference pressure (in pa)")
           read(value,*) pRef

         case ("reference density (in kg/m^3)")
           read(value,*) rhoRef

         case ("reference temperature (in k)")
           read(value,*) TRef

         case ("conversion factor grid units to meter")
           read(value,*) LRef
           LRefSpecified = .true.
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("code name")
           codeName = value

         case ("get coarse-level sol")
           cplGetCoarseSol = checkYesNo(value, keyword)

         case ("mach for initialization")
           read(value,*) MachIni

         case ("pressure for initialization")
           read(value,*) pIni

         case ("density for initialization")
           read(value,*) rhoIni

         case ("velocity direction for initialization")
           read(value,*) velDirIni(1), velDirIni(2), velDirIni(3)
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("input overset donors are guesses")
           oversetDonorsAreGuesses = checkYesNo(value, keyword)

         case ("average restricted residual for blanks")
           avgRestrictResForBlanks = checkYesNo(value, keyword)

         case ("overset interpolation type")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("trilinear")
               oversetInterpType = trilinear
             case default
               write(errorMessage,*) "Unknown overset interpolation &
                                     &type, ", &
                                     trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("overset interpolation type coarse grid")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("trilinear")
               oversetInterpTypeCoarse = trilinear
             case default
               write(errorMessage,*) "Unknown overset interpolation &
                                     &type coarse grid, ", &
                                     trim(value), ", specified"
               if(myID == 0) &
                 call terminate("sumb_setParam", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("allowable donor quality")
           read(value,*) allowableDonorQuality
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in the input file"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamStr

!      ==================================================================

       logical function checkYesNo(value, keyword)
!
!      ******************************************************************
!      *                                                                *
!      * CheckYesNo checks whether the given first string is either     *
!      * yes (.True.) or no (.False.). If the string is neither of      *
!      * these two options an error message is printed.                 *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       implicit none
!
!      Function arguments.
!
       character (len=*), intent(inout) :: value
       character (len=*), intent(in)    :: keyword
!
!      Local variables.
!
       integer :: error

       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Create a lower case version of value

       call convertToLowerCase(value)

       ! Determine the case we are having here.

       select case (value)
         case ("yes")
           checkYesNo = .true.
         case ("no")
           checkYesNo = .false.
         case default
           write(errorMessage,*) trim(keyword), " must be yes or no, &
                                 &not ", trim(value)
           if(myID == 0) &
             call terminate("checkYesNo", errorMessage)
           call mpi_barrier(SUmb_comm_world, error)
       end select

       end function checkYesNo

!      ==================================================================

       function determineDiscretization(value, keyword)
!
!      ******************************************************************
!      *                                                                *
!      * DetermineDiscretization determines the discretization stored   *
!      * in the string value. If it does not match with one of the      *
!      * keywords stored an error message is printed and the program is *
!      * stopped.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use inputDiscretization
       implicit none
!
!      Function type
!
       integer(kind=intType) :: determineDiscretization
!
!      Function arguments.
!
       character (len=*), intent(inout) :: value
       character (len=*), intent(in)    :: keyword
!
!      Local variables.
!
       integer :: error

       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Create a lower case version of value

       call convertToLowerCase(value)

       ! Determine the case we are having here.

       select case (value)
         case ("central plus scalar dissipation")
           determineDiscretization = dissScalar
         case ("central plus matrix dissipation")
           determineDiscretization = dissMatrix
         case ("central plus cusp dissipation")
           determineDiscretization = dissCusp
         case ("upwind")
           determineDiscretization = upwind
         case default
           write(errorMessage,*) "Unknown ", trim(keyword), &
                                 ", ", trim(value), ", specified"
           if(myID == 0) &
             call terminate("determineDicretization", errorMessage)
           call mpi_barrier(SUmb_comm_world, error)
       end select

       end function determineDiscretization

!      ==================================================================

       function determineFileFormat(value, keyword)
!
!      ******************************************************************
!      *                                                                *
!      * determineFileFormat determines the file format stored in       *
!      * value. If it does not match with one of the expected values an *
!      * error message is printed and the program is stopped.           *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       implicit none
!
!      Function type
!
       integer(kind=intType) :: determineFileFormat
!
!      Function arguments
!
       character (len=*), intent(inout) :: value
       character (len=*), intent(in)    :: keyword
!
!      Local variables
!
       integer :: error

       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Create a lower case version of value and check the options.

       call convertToLowerCase(value)

       select case (value)
         case ("cgns")
           determineFileFormat = cgnsFormat
         case ("plot3d")
           determineFileFormat = plot3DFormat
         case default
           write(errorMessage,*) "Unknown ", trim(keyword), &
                                 ", ", trim(value), ", specified"
           if(myID == 0) &
             call terminate("determineFileFormat", errorMessage)
           call mpi_barrier(SUmb_comm_world, error)
       end select

       end function determineFileFormat

!      ==================================================================

       function determineRiemann(value, keyword)
!
!      ******************************************************************
!      *                                                                *
!      * DetermineRiemann determines the riemann solver stored in the   *
!      * string value. If it does not match with one of the keywords    *
!      * stored an error message is printed and the program is stopped. *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use inputDiscretization
       implicit none
!
!      Function type
!
       integer(kind=intType) :: determineRiemann
!
!      Function arguments.
!
       character (len=*), intent(inout) :: value
       character (len=*), intent(in)    :: keyword
!
!      Local variables.
!
       integer :: error

       character (len=2*maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Create a lower case version of value

       call convertToLowerCase(value)

       ! Determine the case we are having here.

       select case (value)
         case ("roe")
           determineRiemann = Roe
         case ("van leer")
           determineRiemann = vanLeer
         case ("ausmdv")
           determineRiemann = ausmdv
         case default
           write(errorMessage,*) "Unknown ", trim(keyword), &
                                 ", ", trim(value), ", specified"
           if(myID == 0) &
             call terminate("determineRiemann", errorMessage)
           call mpi_barrier(SUmb_comm_world, error)
       end select

       end function determineRiemann

!      ==================================================================

       subroutine readMotionCoef(string, start, end, coef)
!
!      ******************************************************************
!      *                                                                *
!      * ReadMotionCoef reads the coefficients start to end from the    *
!      * given string. These coefficients correspond to the description *
!      * of the rigid body motion and are either polynomial or a        *
!      * fourier series. In both cases it is assumed that the number of *
!      * coefficients is specified before the actual coefficients are   *
!      * read.                                                          *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)   :: start, end
       character(len=*),      intent(inout) :: string

       real(kind=realType), dimension(start:*), intent(out) :: coef
!
!      Local variables.
!
       integer :: pos

       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if end >= 0, i.e. if the order of the polynomial/fourier
       ! series is specified before this subroutine is called.

       if(end < 0) then

         ! Processor 0 prints the error message, while the others
         ! wait to get killed.

         if(myID == 0)                       &
           call terminate("readMotionCoef", &
                          "Order of motion coefficients not known yet &
                          &when the coefficients are specified")
         call mpi_barrier(SUmb_comm_world, pos)
       endif

       ! Loop over the number of coefficients to be read.

       do i=start,end

         ! Check if the string still contains data. If not processor
         ! zero prints and error message, while the others wait to get
         ! killed.

         if(len_trim(string) == 0) then
           if(myID == 0)                       &
             call terminate("readMotionCoef", &
                            "Not enough coefficients specified")
           call mpi_barrier(SUmb_comm_world, pos)
         endif

         ! Read the i-th coefficient from the string.

         read(string,*) coef(i)

         ! Remove this coefficient from the string.

         pos = index(string, " ")
         if(pos > 0) then
           string = string(pos:)
           string = adjustl(string)
           string = trim(string)
         else
           string = ""
         endif

       enddo

       ! The length of the string should be zero. If not too much data
       ! is specified. Print a warning to indicate this.

       if(myID == 0 .and. len_trim(string) > 0) then
         print "(a)", "#"
         print "(a)", "#*==================== !!! Warning !!! &
                      &======================"
         print "(a)", "#* Too many coefficients specified for a &
                      &certain rigid body motion"
         print "(a)", "#* Information is ignored."
         print "(a)", "#*=======================================&
                      &===================="
         print "(a)", "#"
       endif

       end subroutine readMotionCoef

!      ==================================================================

       subroutine sumb_setParamLog(value, name)
!
!      ******************************************************************
!      *                                                                *
!      * sumb_setParamLog sets the logical parameter with the given     *
!      * name to the given value.                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       use couplerParam
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: value
       character (len=*), intent(in) :: name
!
!      Local variables.
!
       integer :: pos
       integer(kind=intType) :: nn

       character (len=maxStringLen)   :: keyword
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace all the tab and return characters by spaces and
       ! get rid of the leading and trailing spaces in name.

       keyword = name
       call replaceTabsAndReturns(keyword)
       keyword = adjustl(keyword)
       keyword = trim(keyword)

       ! In case this is an empty string, return.

       if(len_trim(keyword) == 0) return

       ! Convert keyword to lower case, such that a comparison can be
       ! made with the predefined keywords.

       call convertToLowerCase(keyword)
!
!      ******************************************************************
!      *                                                                *
!      * All the initialization stuff has been done.                    *
!      * Now search for keyword in the set of keywords for this code.   *
!      *                                                                *
!      ******************************************************************
!
       select case(keyword)
!
!        ****************************************************************
!        *                                                              *
!        * Discretization parameters.                                   *
!        *                                                              *
!        ****************************************************************
!
         case ("vortex correction")
           vortexCorr = value

         case ("directional dissipation scaling")
           dirScaling = value

         case ("total enthalpy scaling inlet")
           hScalingInlet = value
!
!        ****************************************************************
!        *                                                              *
!        * IO parameters.                                               *
!        *                                                              *
!        ****************************************************************
!
         case ("restart")
           restart = value

         case ("check nondimensionalization")
           checkRestartSol = value

         case ("rind layer in solution files")
           storeRindLayer = value

         case ("automatic parameter update")
           autoParameterUpdate = value

         case ("write coordinates in meter")
           writeCoorMeter = value

         case ("store convergence inner iterations")
           storeConvInnerIter = value
!
!        ****************************************************************
!        *                                                              *
!        * Iteration parameters.                                        *
!        *                                                              *
!        ****************************************************************
!
         case ("freeze turbulent source terms in mg")
           freezeTurbSource = value
!
!        ****************************************************************
!        *                                                              *
!        * Parallel or load balance parameters.                         *
!        *                                                              *
!        ****************************************************************
!
         case ("split blocks for load balance")
           splitBlocks = value
!
!        ****************************************************************
!        *                                                              *
!        * Physics parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("v2f with upper bound")
           rvfB = value

         case ("use wall functions")
           wallFunctions = value
!
!        ****************************************************************
!        *                                                              *
!        * Time spectral parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("write file for unsteady restart")
           writeUnsteadyRestartSpectral = value

         case ("write unsteady volume solution files")
           writeUnsteadyVolSpectral = value

         case ("write unsteady surface solution files")
           writeUnsteadySurfSpectral = value
!
!        ****************************************************************
!        *                                                              *
!        * Unsteady parameters.                                         *
!        *                                                              *
!        ****************************************************************
!
         case ("update wall distance unsteady mode")
           updateWallDistanceUnsteady = value
!
!        ****************************************************************
!        *                                                              *
!        * Visualization parameters.                                    *
!        *                                                              *
!        ****************************************************************
!
         case ("pv3 visualization only")
           PV3VisOnly = value
!
!        ****************************************************************
!        *                                                              *
!        * Coupler parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("get coarse-level sol")
           cplGetCoarseSol = value
!
!        ****************************************************************
!        *                                                              *
!        * Overset parameters.                                          *
!        *                                                              *
!        ****************************************************************
!
         case ("input overset donors are guesses")
           oversetDonorsAreGuesses = value

         case ("average restricted residual for blanks")
           avgRestrictResForBlanks = value
!
!        ****************************************************************
!        *                                                              *
!        * The keyword does not correspond to one of the keywords for   *
!        * the input parameters. It is possible that this is either a   *
!        * family property, which is overwritten or parameters for the  *
!        * level 0 turbine cooling  model. If the keyword does not      *
!        * belong to the above mentioned category processor 0 prints    *
!        * a warning message.                                           *
!        *                                                              *
!        ****************************************************************
!
         case default

           pos = index(keyword, "family")
           if(pos == 0) pos = index(keyword, "cooling plane")

           if(pos == 0 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(3a)", "#* Unknown keyword, ", trim(keyword), &
                           ", encountered in sumb_setParam"
             print "(a)", "#* Information is ignored."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

       end select

       end subroutine sumb_setParamLog

       end module sumb_coupler_m

!
!      ******************************************************************
!      *                                                                *
!      * File:          checkInputParam.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-13-2002                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkInputParam
!
!      ******************************************************************
!      *                                                                *
!      * checkInputParam checks if all necessary data has been          *
!      * specified. If some key data is missing an error message will   *
!      * be printed and the program will exit. Key data depends on the  *
!      * case to be solved. E.g. for the Navier Stokes equations it is  *
!      * necessary to specify the Reynolds number, but for Euler this   *
!      * can be omitted.                                                *
!      *                                                                *
!      * Furthermore warnings are printed in case parameters have been  *
!      * specified that are ignored, e.g. Mach number for internal flow *
!      * computations.                                                  *
!      *                                                                *
!      * Note that only processor 0 prints warning and error messages,  *
!      * such that the output does not become messy.                    *
!      *                                                                *
!      ******************************************************************
!
       use allInputParam
       use communication
       use constants
       use couplerParam
       use flowVarRefState
       use iteration
       use monitor
       use localMG
       implicit none
!
!      Local variables
!
       integer :: ierr

       integer(kind=intType) :: nn, oldSolWrittenSize

       real(kind=realType) :: vecLength, dot

       logical :: gridPrecisionWarning, solPrecisionWarning
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Discretization parameters. Check if the key parameters have    *
!      * been specified and set some coarse grid parameters in case     *
!      * these have not been specified.                                 *
!      *                                                                *
!      ******************************************************************
!
       if(spaceDiscr == none) then
         if(myID == 0)                       &
           call terminate("checkInputParam", &
                          "Discretization scheme not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(spaceDiscrCoarse == none) spaceDiscrCoarse = spaceDiscr

       if(riemannCoarse == none) riemannCoarse = riemann

       ! Set dirScaling to .false. if a scheme other than scalar
       ! dissipation is used.

       if(spaceDiscr /= dissScalar) dirScaling = .false.

       ! Determine whether or not the spectral radIi are needed for
       ! the flux computations.

       radiiNeededFine = .false.
       if(spaceDiscr == dissScalar) radiiNeededFine = .true.

       radiiNeededCoarse = .false.
       if(spaceDiscrCoarse == dissScalar) radiiNeededCoarse = .true.
!
!      ******************************************************************
!      *                                                                *
!      * IO parameters. Check if the grid file has been specified and,  *
!      * if needed, the plot3D connectivity file. Possibly correct the  *
!      * value of restart. Note that restart got the default value of   *
!      * .true. in case no restart file has been specified it is now    *
!      * set to false. Set the names of the solution files if not       *
!      * specified and check if a cp curve fit file has been specified  *
!      * if curve fits must be used.                                    *
!      *                                                                *
!      * If the code has been compiled without cgns check that the file *
!      * format is not cgns.                                            *
!      *                                                                *
!      * Overwrite storeConvInnerIter to .true. if this is not an       *
!      * unsteady computation.                                          *
!      *                                                                *
!      ******************************************************************
!
       if(gridFile == "") then
         if(myID == 0) &
           call terminate("checkInputParam", "Grid file not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(fileFormatRead  == noFormat) fileFormatRead  = fileFormatWrite
       if(fileFormatWrite == noFormat) fileFormatWrite = fileFormatRead

#ifdef USE_NO_CGNS
       if(fileFormatRead  == noFormat) fileFormatRead  = plot3DFormat
       if(fileFormatWrite == noFormat) fileFormatWrite = plot3DFormat
#else
       if(fileFormatRead  == noFormat) fileFormatRead  = cgnsFormat
       if(fileFormatWrite == noFormat) fileFormatWrite = cgnsFormat
#endif

       if(fileFormatRead  == plot3DFormat .or. &
          fileFormatWrite == plot3DFormat) then
         if(plot3DConnFile == "") then
           if(myID == 0)                       &
             call terminate("checkInputParam", &
                            "plot3D connectivity file not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       endif

       if(restartFile == "") restart = .false.

       if(newGridFile == "") then
         select case(fileFormatWrite)
           case (cgnsFormat)
             newGridFile = "NewGrid.cgns"
           case (plot3DFormat)
             newGridFile = "NewGrid.xyz"
         end select
       endif

       if(solFile == "") then
         select case(fileFormatWrite)
           case (cgnsFormat)
             solFile = "SolSUmb.cgns"
           case (plot3DFormat)
             solFile = "SolSUmb.q"
         end select
       endif

       if(surfaceSolFile == "") &
         surfaceSolFile = trim(solfile)//"Surface"

       if(cpModel == cpTempCurveFits .and. cpFile == "") then
         if(myID == 0)                        &
           call terminate("checkInputParam", &
                          "Cp curve fit file not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

#ifdef USE_NO_CGNS

       if(fileFormatRead  == cgnsFormat .or. &
          fileFormatWrite == cgnsFormat) then
         if(myID == 0)                       &
           call terminate("checkInputParam", &
                          "cgns support disabled during compile time")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

#endif

       if(equationMode == unsteady) then
         if(timeIntegrationScheme == explicitRK) &
           storeConvInnerIter = .false.
       else
         storeConvInnerIter = .true.
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Iteration parameters. Check if the key parameters have been    *
!      * specified and set some coarse grid parameters in case these    *
!      * have not been specified.                                       *
!      *                                                                *
!      ******************************************************************
!
       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) then
         smoother = none
       else
         if(smoother == none) then
           if(myID == 0) &
             call terminate("checkInputParam", "Smoother not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(ncycles < 0) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Number of multigrid cycles not or wrongly &
                            &specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(cfl < zero) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "cfl number not or wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(l2Conv <= zero .or. L2Conv >= one) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Relative L2 norm for convergence must be a &
                            & number between 0 and 1.")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(l2ConvCoarse <= zero .or. L2ConvCoarse >= one) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Relative L2 norm for convergence coarse grid &
                            &must be a number between 0 and 1.")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Grid motion parameters. These can only be specified for an     *
!      * external flow problem.                                         *
!      *                                                                *
!      ******************************************************************
!
       if(flowType == internalFlow .and. gridMotionSpecified) then
         if(myID == 0) &
           call terminate("checkInputParam", &
                          "Grid motion specified for an internal flow; &
                          &this is not possible")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Physics parameters. Check if the key parameters have been      *
!      * specified and set the unit vector for the free-stream velocity.*
!      *                                                                *
!      ******************************************************************
!
       if(equations == none) then
         if(myID == 0) &
           call terminate("checkInputParam", "Equations not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(equationMode == none) then
         if(myID == 0) &
           call terminate("checkInputParam", "Mode not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(flowType == none) then
         if(myID == 0) &
           call terminate("checkInputParam", "Flow type not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(Mach < zero .and. flowType == externalFlow) then
         if(myID == 0)                        &
           call terminate("checkInputParam", &
                          "Mach not or wrongly specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       if(equations == NSEquations .or. equations == RANSEquations) then
         if(Reynolds < zero .and. flowType == externalFlow) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Reynolds not or wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       endif

       if(equations == RANSEquations .and. turbModel == none) then
         if(myID == 0)                        &
           call terminate("checkInputParam", &
                          "Turbulence model not specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Create a unit vector for the free stream velocity. It is checked
       ! if the vector specified is a valid one. If not processor 0 prints
       ! an error message. Only for external flows.

       if(flowType == externalFlow) then
         vecLength = sqrt(velDirFreestream(1)*velDirFreestream(1) &
                   +      velDirFreestream(2)*velDirFreestream(2) &
                   +      velDirFreestream(3)*velDirFreestream(3))
         if(vecLength < eps) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Free stream velocity direction wrongly &
                            &specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         vecLength = one/vecLength
         velDirFreestream(1) = velDirFreestream(1)*vecLength
         velDirFreestream(2) = velDirFreestream(2)*vecLength
         velDirFreestream(3) = velDirFreestream(3)*vecLength
       else
         ! Internal flow; simply reset the velocity direction. The value
         ! will be determined later from the inflow boundary conditions.

         velDirFreestream(1) = one
         velDirFreestream(2) = zero
         velDirFreestream(3) = zero
       endif

       ! Set the drag direction to the velocity direction.

       dragDirection = velDirFreestream

       ! Check the lift direction if it was specified for an external
       ! flow. Otherwise set the default direction.

       if(liftDirSpecified .and. flowType == externalFlow) then

         ! Create a unit vector. Perform the same check as for
         ! for the free stream velocity direction.

         vecLength = sqrt(liftDirection(1)*liftDirection(1) &
                   +      liftDirection(2)*liftDirection(2) &
                   +      liftDirection(3)*liftDirection(3))
         if(vecLength < eps) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Lift direction wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         vecLength = one/vecLength
         liftDirection(1) = liftDirection(1)*vecLength
         liftDirection(2) = liftDirection(2)*vecLength
         liftDirection(3) = liftDirection(3)*vecLength

         ! Check the orthogonality with the drag direction.

         dot = liftDirection(1)*dragDirection(1) &
             + liftDirection(2)*dragDirection(2) &
             + liftDirection(3)*dragDirection(3)

         if(abs(dot) > 1.e-3_realType) then
           if(myID == 0)                       &
             call terminate("checkInputParam", &
                            "Lift direction not orthogonal to &
                            &free-stream")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

       else

         ! Lift direction not specified. Set the default direction.
         ! It will have a zero component in the y-direction and a positive
         ! one in the z-direction.

         liftDirection(1) = -dragDirection(3)
         liftDirection(2) =  zero
         liftDirection(3) =  dragDirection(1)

         if(liftDirection(3) < zero) then
           liftDirection(1) = -liftDirection(1)
           liftDirection(3) = -liftDirection(3)
         endif
       endif

       ! Normalize the velocity direction used for the initialization
       ! of a multi-disciplinary computation.

       vecLength = sqrt(velDirIni(1)*velDirIni(1) &
                 +      velDirIni(2)*velDirIni(2) &
                 +      velDirIni(3)*velDirIni(3))
       if(vecLength < eps) then
         if(myID == 0)                       &
           call terminate("checkInputParam", &
                          "Velocity direction for initialization &
                          &wrongly specified")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       vecLength = one/vecLength
       velDirIni(1) = velDirIni(1)*vecLength
       velDirIni(2) = velDirIni(2)*vecLength
       velDirIni(3) = velDirIni(3)*vecLength

       ! Set the Mach number for the coefficients equal to the Mach
       ! number if it was not specified. For internal flow field this
       ! will again be changed in initFlo.

       if(MachCoef < zero) MachCoef = Mach
!
!      ******************************************************************
!      *                                                                *
!      * Time spectral parameters. They only need to be specified for a *
!      * time spectral computation.                                     *
!      *                                                                *
!      ******************************************************************
!
       testSpectral: if(equationMode == timeSpectral) then

         ! Check if the number of time intervals was specified.

         if(nTimeIntervalsSpectral < 0) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Number time intervals spectral not or &
                            &wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! If an unsteady restart solution file must be written, check
         ! if the corresponding time step has been specified.

         if( writeUnsteadyRestartSpectral ) then
           if(dtUnsteadyRestartSpectral <= zero) then
             if(myID == 0)                        &
               call terminate("checkInputParam", &
                              "Time step (in sec) for unsteady restart &
                              &not or wrongly specified.")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif
         endif

         ! If solution files (for postprocessing) must be written,
         ! check if the number has been specified.

         if( writeUnsteadyVolSpectral .or. &
             writeUnsteadySurfSpectral) then
           if(nunsteadySolSpectral <= 0) then
             if(myID == 0)                        &
               call terminate("checkInputParam", &
                              "Number of unsteady solution files &
                              &not or wrongly specified.")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif
         endif

       else testSpectral

         ! No spectral method. Set nTimeIntervalsSpectral to 1.

         nTimeIntervalsSpectral = 1

       endif testSpectral
!
!      ******************************************************************
!      *                                                                *
!      * Unsteady parameters. They only need to be specified for an     *
!      * unsteady computation.                                          *
!      *                                                                *
!      ******************************************************************
!
       testUnsteady: if(equationMode == unsteady) then

         ! Physical time step parameters.

         if(nTimeStepsFine < 0) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Number of unsteady time steps fine grid &
                            &not or wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(nTimeStepsCoarse < 0) nTimeStepsCoarse = nTimeStepsFine

         if(deltaT < 0) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Unsteady time step (in sec) &
                            &not or wrongly specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Check if the rigid body rotation parameters are consistent.
         ! The polynomial rotation coefficients.

         if(degreePolXRot >= 0 .and. &
            .not. allocated(coefPolXRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Polynomial coefficients x-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreePolYRot >= 0 .and. &
            .not. allocated(coefPolYRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Polynomial coefficients y-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreePolZRot >= 0 .and. &
            .not. allocated(coefPolZRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Polynomial coefficients z-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! The fourier rotation coefficients.

         if(degreeFourXRot >= 0 .and. &
            .not. allocated(cosCoefFourXRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier cosine coefficients x-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreeFourXRot >= 1 .and. &
            .not. allocated(sinCoefFourXRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier sine coefficients x-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreeFourYRot >= 0 .and. &
            .not. allocated(cosCoefFourYRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier cosine coefficients y-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreeFourYRot >= 1 .and. &
            .not. allocated(sinCoefFourYRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier sine coefficients y-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreeFourZRot >= 0 .and. &
            .not. allocated(cosCoefFourZRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier cosine coefficients z-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         if(degreeFourZRot >= 1 .and. &
            .not. allocated(sinCoefFourZRot)) then
           if(myID == 0)                        &
             call terminate("checkInputParam", &
                            "Fourier sine coefficients z-rotation &
                            &not specified")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

       endif testUnsteady
!
!      ******************************************************************
!      *                                                                *
!      *                       Warning messages.                        *
!      *                                                                *
!      ******************************************************************
!
       ! Check for an invisid problem if the Reynolds number is specified.
       ! If so, print a Warning that this info is ignored.

       if(myID == 0 .and. equations == EulerEquations .and. &
          Reynolds > zero) then

         print "(a)", "#"
         print "(a)", "#                      Warning"
         print "(a)", "# Reynolds number specified for the Euler &
                      &equations."
         print "(a)", "# This information is ignored."
         print "(a)", "#"

       endif

       ! Check if the Mach and Reynolds number are specified for an
       ! internal flow problem. If so, print a Warning message that this
       ! info is ignored.

       if(flowType == internalFlow) then

         ! Check whether a viscous or an inviscid problem is to be solved.
         ! For an inviscid problem you do not want to mention that the
         ! Reynolds number is ignored, because this has already been
         ! taken care of.

         if((equations == NSEquations .or.     &
             equations == RANSEquations) .and. &
            Mach > zero .and. Reynolds > zero) then

           ! Viscous problem, where both the Mach and Reynolds were
           ! specified. Processor 0 prints the Warning.

           if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print "(a)", "# Mach and Reynolds number specified &
                          &for an internal flow problem."
             print "(a)", "# This information is ignored."
             print "(a)", "#"
           endif

         else if(Mach > zero) then

           ! The Mach number has been specified. Processor 0 prints
           ! a Warning.

           if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print "(a)", "# Mach number specified for an internal &
                          &flow problem."
             print "(a)", "# This information is ignored."
             print "(a)", "#"
           endif

         endif

       endif

       ! For a steady computation possible specified rigid body
       ! rotation info is ignored. Processor 0 will print the Warning.

       if(degreePolXRot  >= 0 .or. degreePolYRot  >= 0 .or. &
          degreePolZRot  >= 0 .or. degreeFourXRot >= 0 .or. &
          degreeFourYRot >= 0 .or. degreeFourZRot >= 0) then

         if(equationMode == steady .and. myID == 0) then
           print "(a)", "#"
           print "(a)", "#                      Warning"
           print "(a)", "# Rigid body rotation info specified for &
                        &a steady computation."
           print "(a)", "# This information is ignored."
           print "(a)", "#"
         endif
       endif

       ! Print warning messages if the precision to be written
       ! is larger than the precision used in the computation.

       gridPrecisionWarning = .false.
       solPrecisionWarning  = .false.

#ifdef USE_SINGLE_PRECISION
       if(precisionGrid == precisionDouble) gridPrecisionWarning = .true.
       if(precisionSol  == precisionDouble) solPrecisionWarning = .true.
#endif

       if(gridPrecisionWarning .and. myID == 0) then
         print "(a)", "#"
         print "(a)", "#                      Warning"
         print "(a)", "# Precision of the grid file to write is &
                      &bigger than used in the computation."
         print "(a)", "# This does not make sense and is a waste &
                      &of disk space"
         print "(a)", "#"
       endif

       if(solPrecisionWarning .and. myID == 0) then
         print "(a)", "#"
         print "(a)", "#                      Warning"
         print "(a)", "# Precision of the solution file to write is &
                      &bigger than used in the computation."
         print "(a)", "# This does not make sense and is a waste &
                      &of disk space"
         print "(a)", "#"
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Wall functions can only be used if the RANS equations are to   *
!      * be solved. If no wall functions are used the wall offset is    *
!      * set to zero.                                                   *
!      *                                                                *
!      ******************************************************************
!
       if(equations /= RANSEquations) wallFunctions = .false.
       if(.not. wallFunctions) wallOffset = zero
!
!      ******************************************************************
!      *                                                                *
!      * Check whether or not the wall distance is needed for the       *
!      * turbulence model.                                              *
!      *                                                                *
!      ******************************************************************
!
       if(equations == RANSEquations) then

         ! RANS simulation. Determine if the turbulence model is
         ! wall distance free. Note that updateWallDistanceUnsteady is
         ! NOT overruled, because this is just the case for which this
         ! parameter was intended.

         select case (turbModel)
           case (komegaWilcox, komegaModified, ktau)

             ! Wall distance free turbulence models.

             wallDistanceNeeded = .false.

           !=============================================================

           case default

             ! The turbulence model needs the wall distance

             wallDistanceNeeded = .true.

         end select

       else

         ! Laminar or inviscid computation. Simply initialize the
         ! logicals for the wall distance to .false.

         wallDistanceNeeded         = .false.
         updateWallDistanceUnsteady = .false.

       endif
!
!      ******************************************************************
!      *                                                                *
!      * Parallelization parameters. Set the minimum load imbalance to  *
!      * 3 percent to avoid any problems.                               *
!      *                                                                *
!      ******************************************************************
!
       loadImbalance = max(loadImbalance, 0.03_realType)
!
!      ******************************************************************
!      *                                                                *
!      * Some default parameters, which depend on other parameters.     *
!      * Only if these have not been specified of course.               *
!      *                                                                *
!      ******************************************************************
!
       if(nsgStartup < 0)    nsgStartup    = 0
       if(ncyclesCoarse < 0) nCyclesCoarse = nCycles
       if(cflCoarse < zero)  cflCoarse     = cfl
       if(betaTurb  < zero)  betaTurb      = alfaTurb

       if(turbRelax == turbRelaxNotDefined) then
         turbRelax = turbRelaxImplicit
         if(turbModel == v2f) turbRelax = turbRelaxExplicit
       endif

       ! V2f should only be solved with explicit underrelaxation.

       if(equations == RANSEquations .and. turbModel == v2f .and. &
          turbRelax == turbRelaxImplicit) then

         turbRelax = turbRelaxExplicit

         if(myID == 0) then
           print "(a)", "#"
           print "(a)", "#                      Warning"
           print "(a)", "# Implicit underrelaxation specified for &
                        &the v2f model."
           print "(a)", "# This is overwritten to explicit &
                        &underrelaxation."
           print "(a)", "#"
         endif

       endif

       if(nsaveVolume <= 0) then
         select case (equationMode)
           case (steady, timeSpectral)
             nSaveVolume = nCycles + nCyclesCoarse + nsgStartup + 1

           case (unsteady)
             nSaveVolume = nTimeStepsFine + nTimeStepsCoarse &
                          + nTimeStepsRestart + 1
         end select
       endif

       if(nsaveSurface <= 0) nSaveSurface  = nSaveVolume

       if(eddyVisInfRatio < zero) then

         ! Default value depends on the turbulence model.

         select case (turbModel)

           case (spalartAllmaras, spalartAllmarasEdwards)
             eddyVisInfRatio = 0.009_realType

           case default
             eddyVisInfRatio = 0.1_realType

         end select
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of old grid levels needed for the BDF     *
!      * time integration of unsteady problems and allocate the memory  *
!      * for the coefficients. The actual values are not yet set,       *
!      * because in the first (and possibly second) time step a reduced *
!      * order must be used, because the older states are not available *
!      * yet. Also allocate the memory for the logicals to indicate     *
!      * whether or not old solutions have been written.                *
!      * If a Runge Kutta scheme must be used for the time integration, *
!      * either explicit or implicit, a separate routine is called to   *
!      * set all the necessary variables.                               *
!      *                                                                *
!      ******************************************************************
!
       select case (timeIntegrationScheme)
         case (BDF)

           ! First check if the accuracy is okay.

           if(timeAccuracy > thirdOrder) then
             if(myID == 0) then
               print "(a)", "#"
               print "(a)", "#                      Warning"
               print "(a)", "# Maximum third order possible for BDF."
               print "(a)", "# Order has been reduced to third."
               print "(a)", "#"
             endif

             timeAccuracy = thirdOrder
           endif

           ! Determine the accuracy and set nOldLevels accordingly.

           select case (timeAccuracy)
             case (firstOrder)
               nOldLevels = 1

             case (secondOrder)
               nOldLevels = 2

             case (thirdOrder)
               nOldLevels = 3
           end select

           ! Allocate the memory for coefTime.

           allocate(coefTime(0:nOldLevels), stat=ierr)
           if(ierr /= 0)                       &
             call terminate("checkInputParam", &
                            "Memory allocation error for coefTime")

         !===============================================================

         case (explicitRK)
           nOldLevels = 1
           call setStageCoeffExplicitRK

         case (implicitRK)
           nOldLevels = 1
           call setStageCoeffImplicitRK

       end select

       ! Set the logicals whether or not the old solutions have been
       ! written. Note that this is only used for the second and
       ! higher order BDF schemes. However it is allocated with a
       ! minimum size of 1 to avoid problems.

       oldSolWrittenSize = max(nOldLevels-1_intType, 1_intType)
       allocate(oldSolWritten(oldSolWrittenSize), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("checkInputParam", &
                        "Memory allocation error for oldSolWritten")

       do nn=1,oldSolWrittenSize
         oldSolWritten(nn) = .false.
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the values of the runge kutta parameters, depending  *
!      * on the number of stages specified.                             *
!      *                                                                *
!      ******************************************************************
!
       ! Limit the number of stages between 1 and 6 and allocate the
       ! memory.

       nRKStages = min(6_intType,max(1_intType,nRKStages))

       allocate(etaRk(nRKStages), cdisRK(nRKStages), stat=ierr)
       if(ierr /= 0) &
         call terminate("checkInputParam", &
                        "Memory allocation error for etaRK and cdisRK")

       ! Determine the case we are having here.

       select case (nRKStages)
         case (1_intType)
           etaRK(1) = one

           cdisRK(1) = one

         case (2_intType)
           etaRK(1) = 0.2222_realType
           etaRK(2) = one

           cdisRK(1) = one
           cdisRK(2) = one

         case (3_intType)
           etaRK(1) = 0.2846_realType
           etaRK(2) = 0.6067_realType
           etaRK(3) = one

           cdisRK(1) = one
           cdisRK(2) = one
           cdisRK(3) = one

         case (4_intType)
           etaRK(1) = 0.33333333_realType
           etaRK(2) = 0.26666667_realType
           etaRK(3) = 0.55555555_realType
           etaRK(4) = one

           cdisRK(1) = one
           cdisRK(2) = half
           cdisRK(3) = zero
           cdisRK(4) = zero

         case (5_intType)
           etaRK(1) = fourth
           etaRK(2) = 0.16666667_realType
           etaRK(3) = 0.37500000_realType
           etaRK(4) = half
           etaRK(5) = one

           cdisRK(1) = one
           cdisRK(2) = zero
           cdisRK(3) = 0.56_realType
           cdisRK(4) = zero
           cdisRK(5) = 0.44_realType

         case (6_intType)
           etaRK(1) = 0.0722_realType
           etaRK(2) = 0.1421_realType
           etaRK(3) = 0.2268_realType
           etaRK(4) = 0.3425_realType
           etaRK(5) = 0.5349_realType
           etaRK(6) = one

           cdisRK(1) = one
           cdisRK(2) = one
           cdisRK(3) = one
           cdisRK(4) = one
           cdisRK(5) = one
           cdisRK(6) = one
       end select
!
!      ******************************************************************
!      *                                                                *
!      * To avoid any problems later on, allocate the memory for the    *
!      * rigid body motion parameters if these values were not present  *
!      * in the parameter file.                                         *
!      *                                                                *
!      ******************************************************************
!
       if(.not. allocated(coefPolXRot) ) then
         allocate(coefPolXRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for coefPolXRot")
         coefPolXRot = zero
       endif

       if(.not. allocated(coefPolYRot) ) then
         allocate(coefPolYRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for coefPolYRot")
         coefPolYRot = zero
       endif

       if(.not. allocated(coefPolZRot) ) then
         allocate(coefPolZRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for coefPolZRot")
         coefPolZRot = zero
       endif

       if(.not. allocated(cosCoefFourXRot) ) then
         allocate(cosCoefFourXRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &cosCoefFourXRot")
         cosCoefFourXRot = zero
       endif

       if(.not. allocated(sinCoefFourXRot) ) then
         allocate(sinCoefFourXRot(1), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &sinCoefFourXRot")
         sinCoefFourXRot = zero
       endif

       if(.not. allocated(cosCoefFourYRot) ) then
         allocate(cosCoefFourYRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &cosCoefFourYRot")
         cosCoefFourYRot = zero
       endif

       if(.not. allocated(sinCoefFourYRot) ) then
         allocate(sinCoefFourYRot(1), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &sinCoefFourYRot")
         sinCoefFourYRot = zero
       endif

       if(.not. allocated(cosCoefFourZRot) ) then
         allocate(cosCoefFourZRot(0:0), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &cosCoefFourZRot")
         cosCoefFourZRot = zero
       endif

       if(.not. allocated(sinCoefFourZRot) ) then
         allocate(sinCoefFourZRot(1), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("checkInputParam", &
                          "Memory allocation failure for &
                          &sinCoefFourZRot")
         sinCoefFourZRot = zero
       endif

       end subroutine checkInputParam

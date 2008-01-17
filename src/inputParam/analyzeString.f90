!
!      ******************************************************************
!      *                                                                *
!      * File:          analyzeString.f90                               *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-12-2002                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine analyzeString(string)
!
!      ******************************************************************
!      *                                                                *
!      * analyzeString extracts a possible input parameter from the     *
!      * given string. If the string is a candidate to carry an input   *
!      * parameter, the keyword part is converted to lower case, such   *
!      * that the specified keyword is case insensitive. If an unknown  *
!      * keyword is encountered an warning message is printed.          *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use flowVarRefState
       use allInputParam
       use localMG
       use couplerParam
       use monitor
       implicit none
!
!      Subroutine argument.
!
       character (len=*), intent(inout) :: string
!
!      Local variables
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
       ! get rid of the leading and trailing spaces in string.

       call replaceTabsAndReturns(string)
       string = adjustl(string)
       string = trim(string)

       ! In case this is an empty string, return.

       if(len_trim(string) == 0) return

       ! In case the first character is a comment sign, return.

       if(string(:1) == "#") return

       ! Find a possible comment sign somewhere in the string. If present
       ! the info following the comment sign is ignored.

       pos = index(string, "#")
       if(pos > 0) then
         string = string(:pos-1)
         string = trim(string)
       endif

       ! Search for the : in the string. If not present, return.

       pos = index(string, ":")
       if(pos == 0) return

       ! Create the strings keyword and value and get rid of the leading
       ! and trailing spaces. As this operation has already been applied
       ! for string, only a trim needs to be done for keyword.

       keyword = string(:pos-1)
       keyword = trim(keyword)

       value = string(pos+1:)
       value = adjustl(value)
       value = trim(value)

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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("number additional turbulence iterations")
           read(value,*) nSubIterTurb

         case ("update bleeds every")
           read(value,*) nUpdateBleeds

         case ("relaxation factor bleed boundary conditions")
           read(value,*) relaxBleeds

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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("restriction relaxation factor")
           read(value,*) fcoll
           fcoll = min(fcoll,one)

         case ("multigrid start level")
           read(value,*) mgStartlevel

         case ("multigrid cycle strategy")
           mgDescription = value
!
!        ****************************************************************
!        *                                                              *
!        * Grid motion parameters.                                      *
!        *                                                              *
!        ****************************************************************
!
         ! Translation parameters.


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
           allocate(coefPolXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &coefPolXRot")

           call readMotionCoef(value, 0_intType, degreePolXRot, &
                               coefPolXRot)

           gridMotionSpecified = .true.

         case ("polynomial coefficients y-rotation")
           nn = max(degreePolYRot,0_intType)
           allocate(coefPolYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &coefPolYRot")


           call readMotionCoef(value, 0_intType, degreePolYRot, &
                               coefPolYRot)

           gridMotionSpecified = .true.

         case ("polynomial coefficients z-rotation")
           nn = max(degreePolZRot,0_intType)
           allocate(coefPolZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
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
           allocate(cosCoefFourXRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &cosCoefFourXRot")

           call readMotionCoef(value, 0_intType, degreeFourXRot, &
                               cosCoefFourXRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients x-rotation")
           nn = max(degreeFourXRot,1_intType)
           allocate(sinCoefFourXRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &sinCoefFourXRot")

           call readMotionCoef(value, 1_intType, degreeFourXRot, &
                               sinCoefFourXRot)

           gridMotionSpecified = .true.

         case ("fourier cosine coefficients y-rotation")
           nn = max(degreeFourYRot,0_intType)
           allocate(cosCoefFourYRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &cosCoefFourYRot")

           call readMotionCoef(value, 0_intType, degreeFourYRot, &
                               cosCoefFourYRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients y-rotation")
           nn = max(degreeFourYRot,1_intType)
           allocate(sinCoefFourYRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &sinCoefFourYRot")

           call readMotionCoef(value, 1_intType, degreeFourYRot, &
                               sinCoefFourYRot)

           gridMotionSpecified = .true.

         case ("fourier cosine coefficients z-rotation")
           nn = max(degreeFourZRot,0_intType)
           allocate(cosCoefFourZRot(0:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
                                        "Memory allocation failure for &
                                        &cosCoefFourZRot")

           call readMotionCoef(value, 0_intType, degreeFourZRot, &
                               cosCoefFourZRot)

           gridMotionSpecified = .true.

         case ("fourier sine coefficients z-rotation")
           nn = max(degreeFourZRot,1_intType)
           allocate(sinCoefFourZRot(1:nn), stat=ierr)
           if(ierr /= 0) call terminate("analyzeString", &
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("v2f version (n1 or n6)")
           read(value,*) rvfN
           if(rvfN /= 1 .and. rvfN /= 6) then
             write(errorMessage,*) "v2f version must be either &
                                    &1 or 6, not ", trim(value)
             if(myID == 0) &
               call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
           liftDirSpecified = .true.

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
         case ("time integration scheme")

           ! Convert value to lower case and check the options.

           call convertToLowerCase(value)

           select case (value)
             case ("bdf")
               timeIntegrationScheme = BDF
             case ("explicit runge-kutta")
               timeIntegrationScheme = explicitRK
             case ("implicit runge-kutta")
               timeIntegrationScheme = implicitRK
             case default
               write(errorMessage,*) "Unknown time integration &
                                     &scheme, ", trim(value),  &
                                     ", specified"
               if(myID == 0) &
                 call terminate("analyzeString", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

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
             case ("fourth")
               timeAccuracy = fourthOrder
             case ("fifth")
               timeAccuracy = fifthOrder
             case default
               write(errorMessage,*) "Unknown time accuracy unsteady, ", &
                                      trim(value), ", specified"
               if(myID == 0) &
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
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
                 call terminate("analyzeString", errorMessage)
               call mpi_barrier(SUmb_comm_world, pos)
           end select

         case ("allowable donor quality")
           read(value,*) allowableDonorQuality

!
!        ****************************************************************
!        *                                                              *
!        * Monitoring of the mass flow of the sliding interfaces.       *
!        *                                                              *
!        ****************************************************************
!
         case ("monitor massflow sliding interfaces")
           monMassSliding = checkYesNo(value, keyword)
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

       end subroutine analyzeString

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
!      Subroutine arguments                                             *
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

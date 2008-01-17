!
!      ******************************************************************
!      *                                                                *
!      * File:          setDefaultValues.F90                            *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 09-19-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setDefaultValues
!
!      ******************************************************************
!      *                                                                *
!      * setDefaultValues sets the default values for the input         *
!      * parameters where-ever possible. The parameters that must be    *
!      * set by the user are initialized such a check can be performed  *
!      * later.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use iteration
       use monitor
       use allInputParam
       use localMG
       use couplerParam
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize monitoring the turbulent residuals as well as the
       ! monitoring of mass flow of the sliding interfaces to .false.

       monDturb       = .false.
       monMassSliding = .false.

       ! Initialize the logicals to check whether or not monitoring,
       ! surface output and volume output variables were specified to
       ! .false.

       monitorSpecified    = .false.
       surfaceOutSpecified = .false.
       volumeOutSpecified  = .false.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the discretization parameters.      *
!      *                                                                *
!      ******************************************************************
!
       spaceDiscr = none                  ! Serves as a check later on.
       orderTurb  = firstOrder            ! First order discretization.
                                          ! Of turbulent advective terms.
       riemann     = Roe
       limiter     = noLimiter            ! No limiter in upwind schemes.
       precond     = noPrecond            ! No preconditioning.

       wallBCTreatment = normalMomentum    ! Normal momentum equation is
                                           ! Used to determine ghost
                                           ! cell pressure.
       outflowTreatment = constantExtrapol ! Constant extrapolation at
                                           ! outflow boundaries.

       spaceDiscrCoarse = none             ! Serves as a check. If nothing
       riemannCoarse    = none             ! is specified the fine grid
                                           ! parameter is taken.

       nonMatchTreatment = NonConservative ! Non conservative treatment
                                           ! of non-matching block to
                                           ! block boundaries.

       vortexCorr = .false.                ! No vortex correction is
                                           ! applied.

       vis2       = half
       vis4       = one/64.0_realType
       vis2Coarse = half

       dirScaling = .true.                 ! Apply isotropic directional
       adis       = two*third              ! scaling in the artificial
                                           ! dissipation schemes.

       hScalingInlet = .false.             ! No total enthalpy scaling.

       kappaCoef = third
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the IO-parameters.                  *
!      *                                                                *
!      ******************************************************************
!
       fileFormatRead  = NoFormat   ! Serves as a check later on.
       fileFormatWrite = NoFormat   ! Idem.

       gridFile       = ""          ! Serves as a check later on.
       plot3DConnFile = ""          ! Idem.
       restartFile    = ""          ! Idem.

       restart         = .true.     ! This will be corrected later if no
                                    ! restart file is specified.
       checkRestartSol = .true.     ! Restart solution is checked for
                                    ! correct nonDimensionalization.

       newGridFile = ""             ! This will be corrected later on
       solFile     = ""             ! if nothing is specified. The
                                    ! default names depend on the
                                    ! format used

       surfaceSolFile = ""          ! This will be corrected later if no
                                    ! surface solution file is specified.

       storeRindLayer = .false.     ! No halo cells in solution files.

       autoParameterUpdate = .true. ! Update the input parameter file
                                    ! when a restart file is written.
       writeCoorMeter = .false.     ! Use original coordinate units
                                    ! when writing solution files.

       cpFile = ""                  ! Serves as a check later on.

       storeConvInnerIter = .false. ! Do not store the convergence of
                                    ! the inner iterations in unsteady
                                    ! mode.

#ifdef USE_SINGLE_PRECISION
       precisionGrid = precisionSingle   ! Default IO precision depends
       precisionSol  = precisionSingle   ! on the default floating
#else                                    ! point type used. Note that
       precisionGrid = precisionDouble   ! for quadrupole precision the
       precisionSol  = precisionDouble   ! IO takes place in double
#endif                                   ! precision.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the iteration parameters.           *
!      *                                                                *
!      ******************************************************************
!
       nCycles       = -1    ! Serves as a check later on.
       nsgStartup    =  0    ! No single grid startup iterations.
       nSubIterTurb  =  0    ! No additional turbulent subiterations.
       nUpdateBleeds = 50    ! Update the bleeds every 50 iterations.

       nSaveVolume  = 0      ! Only save at the end of the computation.
       nSaveSurface = 0

       smoother  = none
       nRKStages = 5

       resAveraging = noResAveraging ! No residual averaging.
       smoop        = 1.5_realType

       turbTreatment     = segregated     ! Segregated solver for the
                                          ! turbulent equations
       turbSmoother      = adi            ! solved using an adi scheme.
       freezeTurbSource = .true.          ! Freeze the coarse grid source
                                          ! terms for a coupled solver.
       turbRelax = turbRelaxNotDefined    ! Will be set later, depending
                                          ! on the turbulence model.

       cfl = -one                         ! Serves as a check later on.

       relaxBleeds = 0.1_realType         ! Relaxation factor for the
                                          ! bleed boundary conditions.

       alfaTurb =  0.8_realType
       betaTurb = -one                    ! Serves as a check later on.

       L2Conv        = 1.e-6_realType     ! Six orders of magnitude for
                                          ! convergence.
       L2ConvCoarse = 1.e-2_realType      ! Only two on coarse grids in
                                          ! full mg.

       nCyclesCoarse = -1             ! If these parameters are not
       cflCoarse     = -one           ! specified the corresponding fine
                                      ! grid values are taken.

       fcoll = one       ! No relaxation when restricting the residuals.

       mgBoundCorr = bcDirichlet0 ! Zero out the boundary halo's for
                                  ! the multigrid corrections.

       mgStartlevel = -1    ! Start at the coarsest grid of the mg cycle
                            ! when no restart is performed.
       mgDescription = "sg" ! Single grid computation.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the motion parameters,              *
!      * i.e. no motion.                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Translation data.


       ! Rotation data.

       rotPoint = zero

       degreePolXRot = -1      ! -1, because the start index is 0.
       degreePolYRot = -1
       degreePolZRot = -1

       degreeFourXRot = -1     ! -1, because the start index is 0,
                               ! at least of the cosine part.
       degreeFourYRot = -1
       degreeFourZRot = -1

       omegaFourXRot = zero
       omegaFourYRot = zero
       omegaFourZRot = zero

       ! The logical to determine whether or not a motion is specified.
       ! Initialize it to .false.

       gridMotionSpecified = .false.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the overset parameters.             *
!      *                                                                *
!      ******************************************************************
!
       oversetDonorsAreGuesses = .false.
       avgRestrictResforBlanks = .false.
       oversetInterpType       = TriLinear
       oversetInterpTypeCoarse = TriLinear
       allowableDonorQuality   = one
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the parallel parameters.            *
!      *                                                                *
!      ******************************************************************
!
       loadImbalance = 0.1_realType  ! Allow 10 percent load imbalance.
       splitBlocks   = .true.        ! Allow the splitting of blocks to
                                     ! obtain a better load balancing.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the physics parameters.             *
!      *                                                                *
!      ******************************************************************
!
       equations     = none       ! These are parameters that must be
       equationMode = none        ! specified. If not, the program
       flowType     = none        ! exits.
       turbModel    = none

       cpModel = cpConstant       ! Constant cp.

       turbProd = strain          ! Strain is used in the production
                                  ! term of transport turbulence models.

       wallFunctions = .false.    ! No wall functions used.

       Mach     = -one            ! Both parameters must be specified
       Reynolds = -one            ! for external flows. The -1. serves
                                  ! as a check later on.

       MachCoef = -one            ! If not specified MachCoef will
                                  ! be set to Mach.

       velDirFreestream(1) = one  ! Free stream velocity
       velDirFreestream(2) = zero ! is specified in the
       velDirFreestream(3) = zero ! x-axis direction.

       liftDirSpecified = .false. ! Lift direction not specified.

       ReynoldsLength = one
       tempFreestream = 288.15_realType
       gammaConstant  = 1.4_realType
       RGasDim        = 287.87_realType

       prandtl     = 0.72_realType
       prandtlTurb = 0.90_realType
       pklim       = 20.0_realType
       wallOffset  = zero

       rvfN = 1                             ! Version 1 of the v2f
                                            ! model is used.
       rvfB = .true.                        ! An upper bound is used
                                            ! in the v2f scales.
       eddyVisInfRatio = -one               ! Default value depends on
                                            ! the turbulence model.
       turbIntensityInf = 0.001_realType

       surfaceRef = one
       lengthRef  = one

       pointRef(1) = zero
       pointRef(2) = zero
       pointRef(3) = zero
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the time spectral parameters.       *
!      *                                                                *
!      ******************************************************************
!
       nTimeIntervalsSpectral   = -1   ! Serves as a check later on.

       nUnsteadySolSpectral = -1       ! Serves as a check later on.

       writeUnsteadyVolSpectral  = .false. ! No writing of the files
       writeUnsteadySurfSpectral = .false. ! for postprocessing.

       writeUnsteadyRestartSpectral = .false. ! No writing of an unsteady
                                              ! mode restart file.

       dtUnsteadyRestartSpectral    = -one    ! Is checked later on.
!
!      ******************************************************************
!      *                                                                *
!      * Set the default values for the unsteady parameters.            *
!      *                                                                *
!      ******************************************************************
!
       timeAccuracy = secondOrder  ! Second order time accuracy.

       nTimeStepsCoarse = -1       ! Serves as a check later on.
       nTimeStepsFine   = -1       ! Serves as a check later on.

       deltaT = -one               ! Serves as a check later on.

       updateWallDistanceUnsteady = .true.  ! This default value is
                                            ! overruled for models that
                                            ! are wall distance free.
!
!      ******************************************************************
!      *                                                                *
!      * Visualization parameters.                                      *
!      *                                                                *
!      ******************************************************************
!
       PV3VisOnly = .false.    ! Perform an actual computation.
!
!      ******************************************************************
!      *                                                                *
!      * The reference state variables. Set them to -1, such that they  *
!      * can be checked later on.                                       *
!      *                                                                *
!      ******************************************************************
!
       pRef   = -one
       rhoRef = -one
       TRef   = -one
!
!      ******************************************************************
!      *                                                                *
!      * The conversion factor of the grid units to meters. Default 1.  *
!      *                                                                *
!      ******************************************************************
!
       LRef           = one
       LRefSpecified = .false.
!
!      ******************************************************************
!      *                                                                *
!      * Initialization of some unsteady restart parameters. These will *
!      * be overwritten when an actual unsteady restart is performed.   *
!      *                                                                *
!      ******************************************************************
!
       nOldSolAvail        = 1
       nTimeStepsRestart   = 0
       timeUnsteadyRestart = zero
!
!      ******************************************************************
!      *                                                                *
!      * Variables needed for the writing of grid and solution files.   *
!      *                                                                *
!      ******************************************************************
!
       timeSpectralGridsNotWritten = .true.
!
!      ******************************************************************
!      *                                                                *
!      * Coupler parameters.                                            *
!      *                                                                *
!      ******************************************************************
!
       codeName        = "SUmb"
       cplGetCoarseSol = .false.

       ! Parameters used for the flow field initialization if not
       ! enough boundary data is present in multidisciplinary mode.

       MachIni =      0.5_realType
       pIni    = 101325.0_realType
       rhoIni  =      1.2_realType

       ! Velocity direction in the free stream.

       velDirIni(1) = one
       velDirIni(2) = zero
       velDirIni(3) = zero

       end subroutine setDefaultValues

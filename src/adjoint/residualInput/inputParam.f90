!
!      ******************************************************************
!      *                                                                *
!      * File:          inputParam.f90                                  *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module accuracy
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
       integer(kind=intType), parameter :: firstOrder  = 1_intType, &
                                           secondOrder = 2_intType, &
                                           thirdOrder  = 3_intType, &
                                           fourthOrder = 4_intType, &
                                           fifthOrder  = 5_intType
       end module accuracy

!      ==================================================================

       module inputDiscretization
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to the discretization of    *
!      * the governing equations, i.e. scheme parameters, time accuracy *
!      * (in case of an unsteady computation) and preconditioning info. *
!      *                                                                *
!      ******************************************************************
!
       use accuracy
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: dissScalar = 1_intType,  &
                                           dissMatrix = 2_intType,  &
                                           dissCusp   = 3_intType,  &
                                           upwind     = 9_intType
       integer(kind=intType), parameter :: Roe     = 1_intType,     &
                                           vanLeer = 2_intType,     &
                                           ausmdv  = 3_intType
       integer(kind=intType), parameter :: noLimiter  = 2_intType,  &
                                           vanAlbeda  = 3_intType,  &
                                           minmod     = 4_intType
       integer(kind=intType), parameter :: noPrecond  = 1_intType,  &
                                           Turkel     = 2_intType,  &
                                           ChoiMerkle = 3_intType
       integer(kind=intType), parameter ::                          &
                                  constantPressure     = 1_intType, &
                                  linExtrapolPressure  = 2_intType, &
                                  quadExtrapolPressure = 3_intType, &
                                  normalMomentum       = 4_intType

       integer(kind=intType), parameter ::                      &
                                  constantExtrapol = 1_intType, &
                                  linExtrapol      = 2_intType

       integer(kind=intType), parameter :: NonConservative = 1_intType, &
                                           Conservative    = 2_intType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the discretization input parameters.             *
!      *                                                                *
!      ******************************************************************
!
       ! spaceDiscr:        Fine grid discretization.
       ! spaceDiscrCoarse:  Coarse grid discretization.
       ! orderTurb:         Order of the discretization of the advective
       !                    terms of the turbulent transport equations.
       !                    Possibilities are 1st and 2nd order.
       ! riemann:           Fine grid riemann solver, upwind schemes only.
       ! riemannCoarse:     Idem, but on the coarse grids.
       ! limiter:           Limiter, upwind schemes only.
       ! precond:           Preconditioner.
       ! wallBCTreatment:   Wall boundary condition treatment.
       ! outflowTreatment:  Treatment of the outflow boundaries. Either
       !                    constantExtrapol or linExtrapol.
       ! nonMatchTreatment: Treatment of the non-matching block
       !                    boundaries. Either NonConservative or
       !                    Conservative.
       ! vis2:              Coefficient of the second order dissipation.
       ! vis4:              Coefficient of the fourth order dissipation.
       ! vis2Coarse:        Coefficient of the second order dissipation
       !                    on the coarser grids in the mg cycle. On the
       !                    coarser grids a first order scheme is used.
       ! adis:              Exponent for directional scaling of the
       !                    dissipation. adis == 0: no directional scaling,
       !                                 adis == 1: isotropic dissipation.
       ! kappaCoef:         Coefficient in the upwind reconstruction
       !                    schemes, both linear and nonlinear.
       ! vortexCorr:        Whether or not a vortex correction must be
       !                    applied. Steady flow only.
       ! dirScaling:        Whether or not directional scaling must be
       !                    applied.
       ! hScalingInlet:     Whether or not the outgoing Riemann invariant
       !                    must be scaled for a subsonic inlet. May be
       !                    needed for stability when strong total
       !                    temperature gradients are present.
       ! radiiNeededFine:   Whether or not the spectral radii are needed
       !                    to compute the fluxes of the fine grid.
       ! radiiNeededCoarse: Idem for the coarse grid.

       integer(kind=intType) :: spaceDiscr, spaceDiscrCoarse
       integer(kind=intType) :: orderTurb, limiter
       integer(kind=intType) :: riemann, riemannCoarse, precond
       integer(kind=intType) :: wallBCTreatment, outflowTreatment
       integer(kind=intType) :: nonMatchTreatment

       real(kind=realType) :: vis2, vis4, vis2Coarse, adis
       real(kind=realType) :: kappaCoef

       logical :: vortexCorr, dirScaling, hScalingInlet
       logical :: radiiNeededFine, radiiNeededCoarse

       end module inputDiscretization

!      ==================================================================

       module inputIO
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to io issues, like file     *
!      * names and corresponding info.                                  *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: noFormat     = 0_intType, &
                                           cgnsFormat   = 1_intType, &
                                           plot3DFormat = 2_intType

       integer(kind=intType), parameter :: precisionSingle = 1_intType, &
                                           precisionDouble = 2_intType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the IO input parameters.                         *
!      *                                                                *
!      ******************************************************************
!
       ! paramFile:           Parameter file, command line argument.
       ! fileFormatRead:      What file format for the grid and
       !                      possibly solution is used when reading.
       !                      Options are cgnsFormat and plot3DFormat.
       ! fileFormatWrite:     Idem, but then for writing.
       ! firstWrite:          Whether or not this is the first time a
       !                      solution is written. Needed when different
       !                      file formats are used for reading and
       !                      writing.
       ! gridFile:            Grid file.
       ! plot3DConnFile:      Connectivity file for the grid if a
       !                      plot3D grid file is used.
       ! newGridFile:         File to which the changed grid is
       !                      written. Needed for moving and/or
       !                      deforming geometries.
       ! restartFile:         Restart solution file; for cgns this
       !                      could be the same as the grid file, but
       !                      not necesarrily.
       ! solFile:             Solution file; for cgns this could be the
       !                      same as the grid or restart file, but not
       !                      necesarrily.
       ! surfaceSolFile:      Surface solution file.
       ! cpFile:              File which contains the curve fits for cp.
       ! precisionGrid:       Precision of the grid file to be written.
       !                      Possibilities are precisionSingle and
       !                      precisionDouble.
       ! precisionSol:        Idem for the solution file(s).
       ! storeRindLayer:      Whether or not to store 1 layer of rind
       !                      (halo) cells in the solution file.
       ! restart:             Whether or not continue from a previous
       !                      computation.
       ! checkRestartSol:     Whether or not the solution in the restart
       !                      file must be checked for correct
       !                      nondimensionalization.
       ! autoParameterUpdate: Whether or not the parameter file must be
       !                      updated automatically. After a restart file
       !                      is written, such that a restart can be made
       !                      without editing the parameter file.
       ! writeCoorMeter:      Whether or not the coordinates in the
       !                      solution files must be written in meters.
       !                      If not, the original units are used.
       ! storeConvInnerIter:  Whether or not to store the convergence of
       !                      the inner iterations for unsteady mode.
       !                      On systems with a limited amount of memory
       !                      the storage of this info could be a
       !                      bottleneck for memory.

       integer(kind=intType) :: fileFormatRead, fileFormatWrite
       integer(kind=intType) :: precisionGrid, precisionSol

       character(len=maxStringLen) :: paramFile, gridFile
       character(len=maxStringLen) :: plot3DConnFile, newGridFile
       character(len=maxStringLen) :: restartFile, solFile
       character(len=maxStringLen) :: surfaceSolFile, cpFile

       logical :: storeRindLayer, restart, checkRestartSol
       logical :: autoParameterUpdate, writeCoorMeter
       logical :: storeConvInnerIter

       logical :: firstWrite = .true.

       end module inputIO

!      ==================================================================

       module inputIteration
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to the iteration process,   *
!      * i.e. multigrid parameters, cfl numbers, smoothers and          *
!      * convergence.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: RungeKutta  = 1_intType,  &
                                           nlLusgs     = 2_intType,  &
                                           nlLusgsLine = 3_intType
       integer(kind=intType), parameter :: segregated = 1_intType,   &
                                           coupled    = 2_intType
       integer(kind=intType), parameter :: gmres = 1_intType,        &
                                           adi   = 2_intType
       integer(kind=intType), parameter :: bcDirichlet0 = 0_intType, &
                                           bcNeumann    = 1_intType
       integer(kind=intType), parameter ::                           &
                                  noResAveraging        = 0_intType, &
                                  alwaysResAveraging    = 1_intType, &
                                  alternateResAveraging = 2_intType
       integer(kind=intType), parameter :: &
                                   turbRelaxNotDefined = 0_intType,  &
                                   turbRelaxExplicit   = 1_intType,  &
                                   turbRelaxImplicit   = 2_intType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the iteration input parameters.                  *
!      *                                                                *
!      ******************************************************************
!
       ! nCycles:          Maximum number of multigrid cycles.
       ! nCyclesCoarse:    Idem, but on the coarse grids in full multigrid.
       ! nSaveVolume:      Number of fine grid cycles after which a volume
       !                   solution file is written.
       ! nSaveSurface:     Number of fine grid cycles after which a
       !                   surface solution file is written.
       ! nsgStartup:       Number of single grid iterations, before
       !                   switching to multigrid. Could be useful for
       !                   supersonic problems with strong shocks.
       ! nSubIterTurb:     Number of turbulent subiterations when using
       !                   a segregated approach for the turbulence.
       ! nUpdateBleeds:    Number of iterations after which the bleed
       !                   boundary conditions must be updated.
       ! smoother:         Smoother to be used.
       ! nRKStages:        Number of stages in the runge kutta scheme.
       ! turbTreatment:    Treatment of the turbulent transport equations;
       !                   either segregated or coupled.
       ! turbSmoother:     Smoother to use in case a segregated solver
       !                   is to be used.
       ! turbRelax:        What kind of turbulent relaxation to use.
       !                   Either turbRelaxExplicit or
       !                   turbRelaxImplicit.
       ! resAveraging:     What kind of residual averaging to use.
       ! freezeTurbSource: Whether or not the turbulent source terms must
       !                   be frozen on the coarser grid levels; only if
       !                   a coupled solver is to be used.
       ! mgBoundCorr:      Treatment of the boundary halo's for the
       !                   multigrid corrections. Either dirichlet0,
       !                   set the corrections to zero, or neumann.
       ! mgStartlevel:     Grid level on which the multigrid must be
       !                   started in the full mg cycle. In case a restart
       !                   is specified this info is overruled and the
       !                   start level is the finest grid.
       ! nMGSteps:         Number of steps in the array cycleStrategy.
       ! nMGLevels:        Number of levels in the multigrid. This info
       !                   is derived from the cycle strategy.
       ! cycleStrategy:    Array which describes the mg cycle.
       ! cfl:              Cfl number on the fine grid.
       ! cflCoarse:        Idem, but on the coarse grids.
       ! alfaTurb:         Relaxation factor in turbulent dd-adi smoother.
       ! betaTurb:         Relaxation factor in vf dd-adi smoother.
       ! relaxBleeds:      Relaxation coefficient for the update
       !                   of the bleed boundary condition.
       ! smoop:            Coefficient in the implicit smoothing.
       ! fcoll:            Relaxation factor for the restricted residuals.
       ! L2Conv:           Relative L2 norm of the density residuals for
       !                   which the computation is assumed converged.
       ! L2ConvCoarse:     Idem, but on the coarse grids during full mg.
       ! etaRk:            Coefficients in the runge kutta scheme. The
       !                   values depend on the number of stages specified.
       ! cdisRk:           Dissipative coefficients in the runge kutta
       !                   scheme. The values depend on the number of
       !                   stages specified.

       integer(kind=intType) :: nCycles, nCyclesCoarse
       integer(kind=intType) :: nSaveVolume, nSaveSurface
       integer(kind=intType) :: nsgStartup, smoother, nRKStages
       integer(kind=intType) :: nSubIterTurb, nUpdateBleeds
       integer(kind=intType) :: resAveraging
       integer(kind=intType) :: turbTreatment, turbSmoother, turbRelax
       integer(kind=intType) :: mgBoundCorr, mgStartlevel
       integer(kind=intType) :: nMGSteps, nMGLevels

       integer(kind=intType), allocatable, dimension(:) :: cycleStrategy

       real(kind=realType) :: cfl, cflCoarse, fcoll, smoop
       real(kind=realType) :: alfaTurb, betaTurb
       real(kind=realType) :: L2Conv, L2ConvCoarse
       real(kind=realType) :: relaxBleeds

       real(kind=realType), allocatable, dimension(:) :: etaRK, cdisRK

       logical :: freezeTurbSource

       end module inputIteration

!      ==================================================================

       module inputMotion
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to the rigid body motion of *
!      * the entire mesh, i.e. translation and rotation.                *
!      * These parameters can only be specified for an external flow    *
!      * computation.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! rotPoint(3): Rotation point of the rigid body rotation.

       real(kind=realType), dimension(3) :: rotPoint

       ! degreePolXRot: Degree of the x-rotation polynomial.
       ! degreePolYRot: Degree of the y-rotation polynomial.
       ! degreePolZRot: Degree of the z-rotation polynomial.

       integer(kind=intType) :: degreePolXRot
       integer(kind=intType) :: degreePolYRot
       integer(kind=intType) :: degreePolZRot

       ! coefPolXRot(0:): coefficients of the x-rotation polynomial.
       ! coefPolYRot(0:): coefficients of the y-rotation polynomial.
       ! coefPolZRot(0:): coefficients of the z-rotation polynomial.

       real(kind=realType), dimension(:), allocatable :: coefPolXRot
       real(kind=realType), dimension(:), allocatable :: coefPolYRot
       real(kind=realType), dimension(:), allocatable :: coefPolZRot

       ! degreeFourXRot: Degree of the x-rotation fourier series.
       ! degreeFourYRot: Degree of the y-rotation fourier series.
       ! degreeFourZRot: Degree of the z-rotation fourier series.

       integer(kind=intType) :: degreeFourXRot
       integer(kind=intType) :: degreeFourYRot
       integer(kind=intType) :: degreeFourZRot

       ! omegaFourXRot: Fourier frequency of the x-rotation; the
       !                   period of the motion is 2*pi/omega.
       ! omegaFourYRot: Fourier frequency of the y-rotation.
       ! omegaFourZRot: Fourier frequency of the z-rotation.

       real(kind=realType) :: omegaFourXRot
       real(kind=realType) :: omegaFourYRot
       real(kind=realType) :: omegaFourZRot

       ! cosCoefFourXRot(0:): cosine coefficients of the
       !                      x-rotation fourier series.
       ! cosCoefFourYRot(0:): cosine coefficients of the
       !                      y-rotation fourier series.
       ! cosCoefFourZRot(0:): cosine coefficients of the
       !                      z-rotation fourier series.

       real(kind=realType), dimension(:), allocatable :: cosCoefFourXRot
       real(kind=realType), dimension(:), allocatable :: cosCoefFourYRot
       real(kind=realType), dimension(:), allocatable :: cosCoefFourZRot

       ! sinCoefFourXRot(1:): sine coefficients of the
       !                      x-rotation fourier series.
       ! sinCoefFourYRot(1:): sine coefficients of the
       !                      y-rotation fourier series.
       ! sinCoefFourZRot(1:): sine coefficients of the
       !                      z-rotation fourier series.

       real(kind=realType), dimension(:), allocatable :: sinCoefFourXRot
       real(kind=realType), dimension(:), allocatable :: sinCoefFourYRot
       real(kind=realType), dimension(:), allocatable :: sinCoefFourZRot

       ! gridMotionSpecified: Whether or not a rigid body motion of
       !                      the grid has been specified.

       logical :: gridMotionSpecified

       end module inputMotion

!      ==================================================================

       module inputParallel
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to the parallelization.     *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! loadImbalance: Allowable load imbalance
       ! splitBlocks:   Whether or not blocks can be split to improve
       !                the load balance.

       real(realType) :: loadImbalance
       logical        :: splitBlocks

       end module inputParallel

!      ==================================================================

       module inputPhysics
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to the physics of the flow, *
!      * like governing equations, mode of the equations, turbulence    *
!      * model and free stream conditions.                              *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: EulerEquations = 1_intType,  &
                                           NSEquations    = 2_intType,  &
                                           RANSEquations  = 3_intType
       integer(kind=intType), parameter :: steady        = 1_intType,   &
                                           unsteady      = 2_intType,   &
                                           timeSpectral  = 3_intType
       integer(kind=intType), parameter :: internalFlow = 1_intType,    &
                                           externalFlow = 2_intType
       integer(kind=intType), parameter :: cpConstant      = 1_intType, &
                                           cpTempCurveFits = 2_intType
       integer(kind=intType), parameter ::                              &
                                  baldwinLomax           =  1_intType,  &
                                  spalartAllmaras        =  2_intType,  &
                                  spalartAllmarasEdwards =  3_intType,  &
                                  komegaWilcox           =  4_intType,  &
                                  komegaModified         =  5_intType,  &
                                  ktau                   =  6_intType,  &
                                  menterSST              =  7_intType,  &
                                  v2f                    = 10_intType
       integer(kind=intType), parameter :: strain       = 1_intType,    &
                                           vorticity    = 2_intType,    &
                                           katoLaunder  = 3_intType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the physics input parameters.                    *
!      *                                                                *
!      ******************************************************************
!
       ! equations:           Governing equations to be solved.
       ! equationMode:        Mode of the equations, steady, unsteady
       !                      or timeSpectral.
       ! flowType:            Type of flow, internal or external.
       ! cpModel:             Which cp model, constant or function of
       !                      temperature via curve fits.
       ! turbModel:           Turbulence model.
       ! turbProd:            Which production term to use in the transport
       !                      turbulence equations, strain, vorticity or
       !                      kato-launder.
       ! rvfN:                Determines the version of v2f turbulence model.
       ! rvfB:                Whether or not to solve v2f with an
       !                      upper bound.
       ! wallFunctions:       Whether or not to use wall functions.
       ! wallDistanceNeeded:  Whether or not the wall distance is needed
       !                      for the turbulence model in a RANS problem.
       ! Mach:                Free stream Mach number.
       ! MachCoef:            Mach number used to compute coefficients;
       !                      only relevant for translating geometries.
       ! velDirFreestream(3): Direction of the free-stream velocity.
       !                      Internally this vector is scaled to a unit
       !                      vector, so there is no need to specify a
       !                      unit vector. Specifying this vector solves
       !                      the problem of angle of attack and yaw angle
       !                      definition as well as the direction of the
       !                      axis (e.g. y- or z-axis in spanwise direction).
       ! liftDirection(3):    Direction vector for the lift.
       ! dragDirection(3):    Direction vector for the drag.
       ! Reynolds:            Reynolds number.
       ! ReynoldsLength:      Length used to compute the Reynolds number.
       ! tempFreestream:      Free stream temperature in Kelvin.
       ! gammaConstant:       Constant specific heat ratio.
       ! RGasDim:             Gas constant in S.I. units.
       ! Prandtl:             Prandtl number.
       ! PrandtlTurb:         Turbulent prandtl number.
       ! pklim:               Limiter for the production of k, the production
       !                      is limited to pklim times the destruction.
       ! wallOffset:          Offset from the wall when wall functions
       !                      are used.
       ! eddyVisInfRatio:     Free stream value of the eddy viscosity.
       ! turbIntensityInf:    Free stream value of the turbulent intensity.
       ! surfaceRef:          Reference area for the force and moments
       !                      computation.
       ! lengthRef:           Reference length for the moments computation.
       ! pointRef(3):         Moment reference point.

       integer(kind=intType) :: equations, equationMode, flowType
       integer(kind=intType) :: turbModel, cpModel, turbProd
       integer(kind=intType) :: rvfN
       logical               :: rvfB

       logical :: wallFunctions, wallDistanceNeeded

       real(kind=realType) :: Mach, MachCoef
       real(kind=realType) :: Reynolds, ReynoldsLength
       real(kind=realType) :: tempFreestream, gammaConstant, RGasDim
       real(kind=realType) :: Prandtl, PrandtlTurb, pklim, wallOffset
       real(kind=realType) :: eddyVisInfRatio, turbIntensityInf
       real(kind=realType) :: surfaceRef, lengthRef

       real(kind=realType), dimension(3) :: velDirFreestream
       real(kind=realType), dimension(3) :: liftDirection
       real(kind=realType), dimension(3) :: dragDirection
       real(kind=realType), dimension(3) :: pointRef

       end module inputPhysics

!      ==================================================================

       module inputTimeSpectral
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters for time spectral problems.                   *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! nTimeIntervalsSpectral: Number of time instances used.

       integer(kind=intType) :: nTimeIntervalsSpectral

       ! dscalar(:,:,:): Matrix for the time derivatices of scalar
       !                 quantities; different for every section to
       !                 allow for different periodic angles.
       !                 The second and third dimension equal the
       !                 number of time intervals.
       ! dvector(:,:,:): Matrices for the time derivatives of vector
       !                 quantities; different for every section to
       !                 allow for different periodic angles and for
       !                 sector periodicity.
       !                 The second and third dimension equal 3 times
       !                 the number of time intervals.

       real(kind=realType), dimension(:,:,:), allocatable :: dscalar
       real(kind=realType), dimension(:,:,:), allocatable :: dvector

       ! writeUnsteadyRestartSpectral: Whether or not a restart file
       !                               must be written, which is
       !                               capable to do a restart in
       !                               unsteady mode.
       ! dtUnsteadyRestartSpectral:    The corresponding time step.

       real(kind=realType) :: dtUnsteadyRestartSpectral
       logical ::             writeUnsteadyRestartSpectral

       ! writeUnsteadyVolSpectral:  Whether or not the corresponding
       !                            unsteady volume solution files
       !                            must be written after the
       !                            computation.
       ! writeUnsteadySurfSpectral: Idem for the surface solution
       !                            files.
       ! nUnsteadySolSpectral:      The corresponding number of
       !                            unsteady solutions to be created.

       integer(kind=intType) :: nUnsteadySolSpectral
       logical ::               writeUnsteadyVolSpectral
       logical ::               writeUnsteadySurfSpectral

       ! rotMatrixSpectral(:,3,3):  The corresponding rotation matrices
       !                            for the velocity. No rotation
       !                            point is needed, because only the
       !                            velocities need to be transformed.
       !                            The matrix stored is the one used
       !                            when the upper bound of the mode
       !                            number is exceeded; for the lower
       !                            bound the inverse (== transpose)
       !                            must be used. The 1st dimension
       !                            is the number of sections.

       real(kind=realType), dimension(:,:,:), allocatable :: &
                                                   rotMatrixSpectral

       end module inputTimeSpectral

!      ==================================================================

       module inputUnsteady
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters for unsteady problems.                        *
!      *                                                                *
!      ******************************************************************
!
       use accuracy
       implicit none
       save

       ! Definition of the parameters for the time integration scheme.

       integer(kind=intType), parameter :: BDF        = 1_intType, &
                                           explicitRK = 2_intType, &
                                           implicitRK = 3_intType

       ! timeIntegrationScheme: Time integration scheme to be used for
       !                        unsteady problems. Possibilities are
       !                        Backward difference schemes, explicit
       !                        RungeKutta schemes and implicit
       !                        RungeKutta schemes.

       integer(kind=intType) :: timeIntegrationScheme

       ! timeAccuracy:     Accuracy of the time integrator for unsteady
       !                   problems. Possibilities are 1st, 2nd and 3rd
       !                   order accurate schemes.
       ! nTimeStepsCoarse: Number of time steps on the coarse mesh;
       !                   only relevant for periodic problems for
       !                   which a full mg can be used.
       ! nTimeStepsFine:   Number of time steps on the fine mesh.
       ! deltaT:           Physical time step in seconds.

       integer(kind=intType) :: timeAccuracy
       integer(kind=intType) :: nTimeStepsCoarse, nTimeStepsFine

       real(kind=realType) :: deltaT

       ! nRKStagesUnsteady:   Number of stages used in the Runge-Kutta
       !                      schemes for a time accurate computation.
       ! betaRKUnsteady(:,:): Matrix with the Runge-Kutta coefficients
       !                      for the residuals.
       ! gammaRKUnsteady(:):  Vector with the time portion of the
       !                      Runge-Kutta stages.

       integer(kind=intType) :: nRKStagesUnsteady

       real(kind=realType), dimension(:,:), allocatable :: betaRKUnsteady
       real(kind=realType), dimension(:),   allocatable :: gammaRKUnsteady

       ! nOldGridRead: Number of old grid levels read from the grid
       !               files. Needed only for a consistent restart
       !               on the deforming meshes.

       integer(kind=intType) :: nOldGridRead

       ! updateWallDistanceUnsteady: Whether or not to update the wall
       !                             distance in unsteady mode. For a
       !                             RANS simulation on a changing grid
       !                             this should be done if the
       !                             turbulence model requires the wall
       !                             distance. However, the user may
       !                             overrule this if he thinks it is
       !                             not necessary.

       logical :: updateWallDistanceUnsteady

       end module inputUnsteady

!      ==================================================================

       module inputVisualization
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters for visualization.                            *
!      *                                                                *
!      ******************************************************************
!
       implicit none
       save

       ! PV3VisOnly: Whether or not to run in visualization mode only.

       logical :: PV3VisOnly

       end module inputVisualization

!      ==================================================================

       module inputOverset
!
!      ******************************************************************
!      *                                                                *
!      * Input parameters which are related to verset grid assembly and *
!      * interpolation procedures.                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which make the code more         *
!      * readable. The actual values of this parameters are arbitrary;  *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: TriLinear    = 1_intType

       integer(kind=intType), dimension(1), parameter :: &
                                                nDonorWeights = (/ 8 /)
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the overset input parameters.                    *
!      *                                                                *
!      ******************************************************************
!
       ! oversetDonorsAreGuesses: Whether or not the input overset donors
       !                          should be treated as guesses, which
       !                          causes the interpolants to be ignored
       !                          and are determined automatically.
       ! avgRestrictResforBlanks: Whether or not to amplify or average
       !                          the restricted residual in multigrid
       !                          to account for the fact that a coarse
       !                          unblanked cell may contain blanked
       !                          cells on the next finer level.
       ! oversetInterpType:       Type of interpolation to use on the
       !                          fine grid level.
       ! oversetInterpTypeCoarse: Idem for the coarse levels.
       ! allowableDonorQuality:   The cut-off value for the quality of
       !                          a donor stencil when searches are
       !                          performed.

       logical :: oversetDonorsAreGuesses, avgRestrictResforBlanks

       integer(kind=intType) :: oversetInterpType
       integer(kind=intType) :: oversetInterpTypeCoarse

       real(kind=realType) :: allowableDonorQuality

       end module inputOverset

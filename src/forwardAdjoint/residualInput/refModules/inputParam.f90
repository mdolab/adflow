!
!      ******************************************************************
!      *                                                                *
!      * File:          inputParam.f90                                  *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                C.A.(Sandy) Mader                               *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 09-17-2009                                      *
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
       integer(kind=intType), parameter :: firstOrder  = 1, &
                                           secondOrder = 2, &
                                           thirdOrder  = 3, &
                                           fourthOrder = 4, &
                                           fifthOrder  = 5
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
       integer(kind=intType), parameter :: dissScalar = 1,  &
                                           dissMatrix = 2,  &
                                           dissCusp   = 3,  &
                                           upwind     = 9
       integer(kind=intType), parameter :: Roe     = 1,     &
                                           vanLeer = 2,     &
                                           ausmdv  = 3
       integer(kind=intType), parameter :: noLimiter  = 2,  &
                                           vanAlbeda  = 3,  &
                                           minmod     = 4
       integer(kind=intType), parameter :: noPrecond  = 1,  &
                                           Turkel     = 2,  &
                                           ChoiMerkle = 3
       integer(kind=intType), parameter ::                          &
                                  constantPressure     = 1, &
                                  linExtrapolPressure  = 2, &
                                  quadExtrapolPressure = 3, &
                                  normalMomentum       = 4

       integer(kind=intType), parameter ::                      &
                                  constantExtrapol = 1, &
                                  linExtrapol      = 2

       integer(kind=intType), parameter :: NonConservative = 1, &
                                           Conservative    = 2
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
       ! lumpedDiss :       logical factor for determining whether or not
       !                    lumped dissipation is used for preconditioner
       ! sigma    : Scaling parameter for dissipation lumping in approximate
       !            precondtioner

       integer(kind=intType) :: spaceDiscr, spaceDiscrCoarse
       integer(kind=intType) :: orderTurb, limiter
       integer(kind=intType) :: riemann, riemannCoarse, precond
       integer(kind=intType) :: wallBCTreatment, outflowTreatment
       integer(kind=intType) :: nonMatchTreatment

       real(kind=realType) :: vis2, vis4, vis2Coarse, adis
       real(kind=realType) :: kappaCoef

       real(kind=realType) :: vis2b, vis4b, kappaCoefb,adisb

       logical :: vortexCorr, dirScaling, hScalingInlet
       logical :: radiiNeededFine, radiiNeededCoarse

       logical :: lumpedDiss
       real(kind=realType) :: sigma

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
       integer(kind=intType), parameter :: noFormat     = 0, &
                                           cgnsFormat   = 1, &
                                           plot3DFormat = 2

       integer(kind=intType), parameter :: precisionSingle = 1, &
                                           precisionDouble = 2
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
       integer(kind=intType), parameter :: RungeKutta  = 1,  &
                                           nlLusgs     = 2,  &
                                           nlLusgsLine = 3
       integer(kind=intType), parameter :: segregated = 1,   &
                                           coupled    = 2
       integer(kind=intType), parameter :: gmres = 1,        &
                                           adi   = 2
       integer(kind=intType), parameter :: bcDirichlet0 = 0, &
                                           bcNeumann    = 1
       integer(kind=intType), parameter ::                           &
                                  noResAveraging        = 0, &
                                  alwaysResAveraging    = 1, &
                                  alternateResAveraging = 2
       integer(kind=intType), parameter :: &
                                   turbRelaxNotDefined = 0,  &
                                   turbRelaxExplicit   = 1,  &
                                   turbRelaxImplicit   = 2
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
       ! printIterations: If True, iterations are printed to stdout
       integer(kind=intType) :: nCycles, nCyclesCoarse
       integer(kind=intType) :: nSaveVolume, nSaveSurface
       integer(kind=intType) :: nsgStartup, smoother, nRKStages
       integer(kind=intType) :: nSubIterTurb, nUpdateBleeds
       integer(kind=intType) :: resAveraging
       integer(kind=intType) :: turbTreatment, turbSmoother, turbRelax
       integer(kind=intType) :: mgBoundCorr, mgStartlevel
       integer(kind=intType) :: nMGSteps, nMGLevels

       integer(kind=intType), allocatable, dimension(:) :: cycleStrategy
       integer(kind=intType) :: miniterNum
       real(kind=realType) :: cfl, cflCoarse, fcoll, smoop
       real(kind=realType) :: alfaTurb, betaTurb
       real(kind=realType) :: L2Conv, L2ConvCoarse
       real(kind=realType) :: L2ConvRel
       real(kind=realType) :: maxL2DeviationFactor
       real(kind=realType) :: relaxBleeds
       real(kind=realtype) :: epscoefconv
       integer(kind=inttype) :: convcheckwindowsize

       real(kind=realType), allocatable, dimension(:) :: etaRK, cdisRK,cdisrkb

       logical :: freezeTurbSource
       logical :: printIterations

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
       real(kind=realType), dimension(3) :: rotPointd

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

       real(kind=realType) :: omegaFourXRot,omegaFourXRotb
       real(kind=realType) :: omegaFourYRot,omegaFourYRotb
       real(kind=realType) :: omegaFourZRot,omegaFourZRotb

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

       ! degreePolAlpha: Degree of the Alpha polynomial.
   
       integer(kind=intType) :: degreePolAlpha
      
       ! coefPolAlpha(0:): coefficients of the Alpha polynomial.
     
       real(kind=realType), dimension(:), allocatable :: coefPolAlpha
       real(kind=realType), dimension(:), allocatable :: coefPolAlphab
 
       ! degreeFourAlpha: Degree of the Alpha fourier series.
      
       integer(kind=intType) :: degreeFourAlpha

       ! omegaFourAlpha: Fourier frequency of the Alpha; the
       !                   period of the motion is 2*pi/omega.
   
       real(kind=realType) :: omegaFourAlpha,omegafouralphab
  
       ! cosCoefFourAlpha(0:): cosine coefficients of the
       !                      x-rotation fourier series.

       real(kind=realType), dimension(:), allocatable :: cosCoefFourAlpha
       real(kind=realType), dimension(:), allocatable :: cosCoefFourAlphab
       
       ! sinCoefFourAlpha(1:): sine coefficients of the
       !                      Alpha fourier series.

       real(kind=realType), dimension(:), allocatable :: sinCoefFourAlpha
       real(kind=realType), dimension(:), allocatable :: sinCoefFourAlphab

       ! degreePolXRot: Degree of the Beta polynomial.
  
       integer(kind=intType) :: degreePolBeta  

       ! coefPolXRot(0:): coefficients of the Beta polynomial.
 
       real(kind=realType), dimension(:), allocatable :: coefPolBeta
       real(kind=realType), dimension(:), allocatable :: coefPolBetab

       ! degreeFourBeta: Degree of the Beta fourier series.
   
       integer(kind=intType) :: degreeFourBeta

       ! omegaFourBeta: Fourier frequency of the Beta; the
       !                   period of the motion is 2*pi/omega.
 
       real(kind=realType) :: omegaFourBeta,omegafourbetab

       ! cosCoefFourBeta(0:): cosine coefficients of the
       !                      Beta fourier series.
    
       real(kind=realType), dimension(:), allocatable :: cosCoefFourBeta
       real(kind=realType), dimension(:), allocatable :: cosCoefFourBetab

       ! sinCoefFourBeta(1:): sine coefficients of the
       !                      Beta fourier series.
 
       real(kind=realType), dimension(:), allocatable :: sinCoefFourBeta
       real(kind=realType), dimension(:), allocatable :: sinCoefFourBetab

       ! degreePolMach: Degree of the Mach polynomial.
 
       integer(kind=intType) :: degreePolMach

       ! coefPolMach(0:): coefficients of the Mach polynomial.

       real(kind=realType), dimension(:), allocatable :: coefPolMach
       real(kind=realType), dimension(:), allocatable :: coefPolMachb

       ! degreeFourMach: Degree of the Mach fourier series.
 
       integer(kind=intType) :: degreeFourMach

       ! omegaFourMach: Fourier frequency of the Mach Number; the
       !                   period of the motion is 2*pi/omega.
 
       real(kind=realType) :: omegaFourMach,omegafourmachb

       ! cosCoefFourMach(0:): cosine coefficients of the
       !                      Mach Number fourier series.
 
       real(kind=realType), dimension(:), allocatable :: cosCoefFourMach
       real(kind=realType), dimension(:), allocatable :: cosCoefFourMachb

       ! sinCoefFourMach(1:): sine coefficients of the
       !                      Mach Number fourier series.
  
       real(kind=realType), dimension(:), allocatable :: sinCoefFourMach
       real(kind=realType), dimension(:), allocatable :: sinCoefFourMachb

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
       integer(kind=intType), parameter :: EulerEquations = 1,  &
                                           NSEquations    = 2,  &
                                           RANSEquations  = 3
       integer(kind=intType), parameter :: steady        = 1,   &
                                           unsteady      = 2,   &
                                           timeSpectral  = 3
       integer(kind=intType), parameter :: internalFlow = 1,    &
                                           externalFlow = 2
       integer(kind=intType), parameter :: cpConstant      = 1, &
                                           cpTempCurveFits = 2
       integer(kind=intType), parameter ::                              &
                                  baldwinLomax           =  1,  &
                                  spalartAllmaras        =  2,  &
                                  spalartAllmarasEdwards =  3,  &
                                  komegaWilcox           =  4,  &
                                  komegaModified         =  5,  &
                                  ktau                   =  6,  &
                                  menterSST              =  7,  &
                                  v2f                    = 10
       integer(kind=intType), parameter :: strain       = 1,    &
                                           vorticity    = 2,    &
                                           katoLaunder  = 3
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
       ! MachGrid:            Mach number of the Mesh. Used in stability 
       !                      derivative calculations. Specified as the
       !                      negative of the desired freestream Mach number.
       !                      When this option is set, set Mach = 0.0...
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
       ! pointRefEC(3):       Elastic center. Bending moment refernce point

       integer(kind=intType) :: equations, equationMode, flowType
       integer(kind=intType) :: turbModel, cpModel, turbProd
       integer(kind=intType) :: rvfN
       logical               :: rvfB

       logical :: wallFunctions, wallDistanceNeeded

       real(kind=realType) :: Mach, MachCoef,MachGrid
       !AD derivative values
       real(kind=realType) :: Machd, MachCoefd,MachGridd
       real(kind=realType) :: Reynolds, ReynoldsLength
       real(kind=realType) :: tempFreestream, gammaConstant, RGasDim
       real(kind=realType) :: Prandtl, PrandtlTurb, pklim, wallOffset
       real(kind=realType) :: eddyVisInfRatio, turbIntensityInf
       real(kind=realType) :: surfaceRef, lengthRef
       !bending moment derivative
       real(kind=realType) :: lengthrefb

       real(kind=realType), dimension(3) :: velDirFreestream,velDirFreestreamd
       real(kind=realType), dimension(3) :: liftDirection,liftDirectiond
       real(kind=realType), dimension(3) :: dragDirection,dragDirectiond
       real(kind=realType), dimension(3) :: pointRef
       !bending moment derivative
       real(kind=realType), dimension(3) :: pointRefb
       real(kind=realType), dimension(3) :: pointRefEC

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

       integer(kind=intType), parameter :: BDF        = 1, &
                                           explicitRK = 2, &
                                           implicitRK = 3, &
                                           MD         = 4
       
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
       integer(kind=intType), parameter :: TriLinear    = 1

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

       module inputADjoint
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters ADjoint.                         *
!      * The actual values of this parameters are arbitrary;            *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       !definition of parameters for the ADjoint Linear Solver

       integer(kind=intType), parameter :: PETSCBICGStab = 1, &
                                           PETSCGMRES    = 2, &
                                           PETSCCG       = 3, &
                                           PETSCFGMRES   = 4
       
       !Definitions for Global Preconditioners
       integer(kind=intType), parameter :: BlockJacobi      = 1, &
                                           Jacobi           = 2, &
                                           AdditiveSchwartz = 3

       !Definitions for local Preconditioners
       integer(kind=intType), parameter :: ILU       = 1, &
                                           ICC       = 2, &
                                           LU        = 3, &
                                           Cholesky  = 4

       !Definitions for matrix ordering method
       integer(kind=intType), parameter :: Natural               = 1, &
                                           ReverseCuthillMckee   = 2, &
                                           NestedDissection      = 3, &
                                           OnewayDissection      = 4, &
                                           QuotientMinimumDegree = 5

       !Definitions for Preconditioner Side
       integer(kind=intType), parameter :: Left  = 1, &
                                           Right = 2
       
       !Definitions for matrix ordering method
       integer(kind=intType), parameter :: Normal = 1, &
                                           RowMax = 2, &
                                           RowSum = 3, &
                                           RowAbs = 4
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the adjoint input parameters.                    *
!      *                                                                *
!      ******************************************************************
!
       ! solveADjoint : Whether or not the adjoint equations should be solved
       ! Monitor      : Whether or not to enable the monitor for the KSP 
       !                contexts.
       ! ApproxPC     : Whether or not to use the approximate jacobian 
       !                preconditioner
       ! restartADjoint: Whether or not we want to restart the adjoint 
       !                 from the previous solution
       ! useDiagTSPC   : Whether or not the off time instance terms are
       !                 included in the TS preconditioner.
       logical :: solveADjoint, setMonitor, ApproxPC,restartADjoint,useDiagTSPC

       ! ADjointSolverType: Type of linear solver for the ADjoint
       ! PreCondType      : Type of Preconditioner to use
       ! ScaleType        : Type of Scaling to use in Jacobi Preconditioner
       ! Matrix Ordering  : Type of matrix ordering to use
       integer(kind=intType)::ADjointSolverType,PreCondType,ScaleType,&
            MatrixOrdering
       
       ! PCSide  : Right or Left Preconditioning
       ! LocalPC : Local preconditioning type 
       integer(kind=intType):: PCSide,LocalPCType

       ! FillLevel     : Number of levels of fill for the ILU local PC
       ! Overlap       : Amount of overlap in the ASM PC
       integer(kind=intType):: FillLevel, Overlap

       ! adjRelTol     : Relative tolerance
       ! adjAbsTol     : Absolute tolerance
       ! adjDivTol     : Relative tolerance increase to divergence
       ! adjMaxIter    : Maximum number of iterations
       ! adjRestart    : Maximum number of steps before restart
       !                 It has a high impact on the required memory!
       ! adjMonStep    : Convergence monitor step
       
       real(kind=realType)    :: adjRelTol  
       real(kind=realType)    :: adjAbsTol
       real(kind=realType)    :: adjRelTolRel
       real(kind=realType)    :: adjDivTol  
       integer(kind=intType)  :: adjMaxIter 
       integer(kind=intType)  :: adjRestart 
       integer(kind=intType)  :: adjMonStep 

       
       real(kind=realType)    :: sigmab
       logical :: printTiming
       logical :: finitedifferencepc
       integer(kind=intType) :: subKSPSubspaceSize
     end module inputADjoint

     module inputTSStabDeriv
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters for Time Spectral stability      *
!      * derivatives.                                                   *
!      * The actual values of this parameters are arbitrary;            *
!      * in the code always the symbolic names are (should be) used.    *
!      *                                                                *
!      *******************************************m***********************
!

       ! TSStability : Whether or not the TS stability derivatives should
       !               be computed
       logical:: TSStability,TSAlphaMode,TSBetaMode,TSpMode,&
            TSqMode,TSrMode,TSAltitudeMode,TSMachMode

       ! useWindAxis : whether to rotate around the wind axis or the body
       !               axis...
       logical:: useWindAxis

     end module inputTSStabDeriv

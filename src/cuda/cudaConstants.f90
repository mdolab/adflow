module cudaFlowVarRefState
    use cudafor
    use precision, only: realType, intType
    implicit none
    save 

    ! nw:       Total number of independent variables including the
    !           turbulent variables.
    ! nwf:      Number of flow variables. For perfect gas computations
    !           this is 5.
    ! nwt:      Number of turbulent variables, nwt = nw - nwf
    ! nt1, nt2: Initial and final indices for turbulence variables

    integer(kind=intType),constant :: nw, nwf, nwt, nt1, nt2

    ! pRef:     Reference pressure (in Pa) used to nondimensionalize
    !           the flow equations.
    ! rhoRef:   Reference density (in kg/m^3) used to
    !           nondimensionalize the flow equations.
    ! TRef:     Reference temperature (in K) used to nondimensionalize
    !           the flow equations.
    ! muRef:    Scale for the viscosity,
    !           muRef = rhoRef*sqrt(pRef/rhoRef); there is also a
    !           reference length in the nondimensionalization of mu,
    !           but this is 1.0, because all the coordinates are
    !           converted to meters.
    ! timeRef:  time scale; needed for a correct
    !           nondimensionalization of unsteady problems.
    !           timeRef = sqrt(rhoRef/pRef); for the reference
    !           length, see the comments for muRef.
    ! uRef:     velocity scale;
    !           uRef = sqrt(pRef/rhoRef);
    ! hRef:     enthalpy scale;
    !           hRef = pRef/rhoRef;

    real(kind=realType),constant :: pRef, rhoRef, TRef
    real(kind=realType),constant :: muRef, timeRef, uRef, href

    ! LRef:          Conversion factor of the length unit of the
    !                grid to meter. e.g. if the grid is in mm.,
    !                LRef = 1.e-3.
    ! LRefSpecified: Whether or not a conversion factor is specified
    !                in the input file.

    real(kind=realType),constant :: LRef
    logical,constant :: LRefSpecified

    ! pInfDim:   Free stream pressure in Pa.
    ! rhoInfDim: Free stream density in kg/m^3.
    ! muDim:     Free stream molecular viscosity in kg/(m s)

    real(kind=realType),constant :: pInfDim, rhoInfDim, muDim, TinfDim, muInfDim
    !AD derivative values

    ! wInf(nw): Nondimensional free stream state vector.
    !           Variables stored are rho, u, v, w and rhoE.
    ! pInf:     Nondimensional free stream pressure.
    ! pInfCorr: Nondimensional free stream pressure, corrected for
    !           a possible presence of 2/3 rhok.
    ! rhoInf:   Nondimensional free stream density.
    ! uInf:     Nondimensional free stream velocity
    ! muInf:    Nondimensional free stream viscosity.
    ! RGas:     Nondimensional gas constant.
    ! gammaInf: Free stream specific heat ratio.

    real(kind=realType),constant :: rhoInf, uInf, pInf, pInfCorr
    real(kind=realType),constant :: RGas, muInf, gammaInf
    
    real(kind=realType), dimension(10),constant :: wInf

    ! viscous:   whether or not this is a viscous computation.
    ! kPresent:  whether or not a turbulent kinetic energy is present
    !            in the turbulence model.
    ! eddyModel: whether or not the turbulence model is an eddy
    !            viscosity model.
    logical,constant :: kPresent, eddyModel, viscous

end module cudaFlowVarRefState

module cudaInputDiscretization
    !
    !       Input parameters which are related to the discretization of
    !       the governing equations, i.e. scheme parameters, time accuracy
    !       (in case of an unsteady computation) and preconditioning info.
    !
    use constants, only: intType, realType
    implicit none
    save
    !
    !       Definition of the discretization input parameters.
    !
    ! spaceDiscr:             Fine grid discretization.
    ! spaceDiscrCoarse:       Coarse grid discretization.
    ! orderTurb:              Order of the discretization of the advective
    !                         terms of the turbulent transport equations.
    !                         Possibilities are 1st and 2nd order.
    ! riemann:                Fine grid riemann solver, upwind schemes only.
    ! riemannCoarse:          Idem, but on the coarse grids.
    ! limiter:                Limiter, upwind schemes only.
    ! precond:                Preconditioner.
    ! eulerWallBCTreatment:   Wall boundary condition treatment for inviscid
    !                         simulations.
    ! viscWallBCTreatment:    Wall boundary condition treatment for viscous
    !                         simulations.
    ! outflowTreatment:       Treatment of the outflow boundaries. Either
    !                         constantExtrapol or linExtrapol.
    ! nonMatchTreatment:      Treatment of the non-matching block
    !                         boundaries. Either NonConservative or
    !                         Conservative.
    ! vis2:                   Coefficient of the second order dissipation.
    ! vis4:                   Coefficient of the fourth order dissipation.
    ! vis2Coarse:             Coefficient of the second order dissipation
    !                         on the coarser grids in the mg cycle. On the
    !                         coarser grids a first order scheme is used.
    ! adis:                   Exponent for directional scaling of the
    !                         dissipation. adis == 0: no directional scaling,
    !                                      adis == 1: isotropic dissipation.
    ! kappaCoef:              Coefficient in the upwind reconstruction
    !                         schemes, both linear and nonlinear.
    ! vortexCorr:             Whether or not a vortex correction must be
    !                         applied. Steady flow only.
    ! dirScaling:             Whether or not directional scaling must be
    !                         applied.
    ! hScalingInlet:          Whether or not the outgoing Riemann invariant
    !                         must be scaled for a subsonic inlet. May be
    !                         needed for stability when strong total
    !                         temperature gradients are present.
    ! radiiNeededFine:        Whether or not the spectral radii are needed
    !                         to compute the fluxes of the fine grid.
    ! radiiNeededCoarse:      Idem for the coarse grid.
    ! lumpedDiss :            logical factor for determining whether or not
    !                         lumped dissipation is used for preconditioner
    ! approxSA:               Determines if the approximate source terms form
    !                         the SA model is used.
    ! sigma      :            Scaling parameter for dissipation lumping in
    !                         approximateprecondtioner
    ! useApproxWallDistance : logical to determine if the user wants to
    !                         use the fast approximate wall distance
    !                         computations. Typically only used for
    !                         repeated calls when the wall distance would
    !                         not have changed significantly
    ! updateWallAssociation : Logical to determine if the full wall distance
    !                         assocation is to be performed on the next
    !                         wall distance calculation. This is only
    !                         significant when useApproxWallDistance is
    !                         set to True. This allows the user to
    !                         reassociate the face a cell is associated
    !                         with.
    ! lowspeedpreconditoner:  Whether or not to use low-speed precondioner

    integer(kind=intType),constant :: spaceDiscr, spaceDiscrCoarse
    integer(kind=intType),constant :: orderTurb, limiter
    integer(kind=intType),constant :: riemann, riemannCoarse, precond
    integer(kind=intType),constant :: eulerWallBCTreatment, viscWallBCTreatment, outflowTreatment
    integer(kind=intType),constant :: nonMatchTreatment

    real(kind=realType),constant :: vis2, vis4, vis2Coarse, adis
    real(kind=realType),constant :: acousticScaleFactor
    real(kind=realType),constant :: kappaCoef
    logical,constant :: lumpedDiss
    logical,constant :: approxSA
    real(kind=realType),constant :: sigma
    logical,constant :: useBlockettes


    logical,constant :: vortexCorr, dirScaling, hScalingInlet
    logical,constant :: radiiNeededFine, radiiNeededCoarse

    logical,constant :: useApproxWallDistance
    logical,constant :: lowSpeedPreconditioner
end module cudaInputDiscretization

module cudaInputIteration
    !
    !       Input parameters which are related to the iteration process,
    !       i.e. multigrid parameters, cfl numbers, smoothers and
    !       convergence.
    !
    use constants
    implicit none
    save
    !
    !       Definition of the iteration input parameters.
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
    !                   a decoupled approach for the turbulence.
    ! nUpdateBleeds:    Number of iterations after which the bleed
    !                   boundary conditions must be updated.
    ! smoother:         Smoother to be used.
    ! nRKStages:        Number of stages in the runge kutta scheme.
    ! nSubiterations:   Maximum number of subiterations used in
    !                   DADI.
    ! turbTreatment:    Treatment of the turbulent transport equations;
    !                   either decoupled or coupled.
    ! turbSmoother:     Smoother to use in case a decoupled solver
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
    ! cfllimit          Limit used to determine how much residuals are smoothed
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
    ! turbresscale: Scaling factor for turbulent residual. Necessary for
    !            NKsolver with RANS. Only tested on SA.
    ! iterType : String used for specifying which type of iteration was taken
    !
    ! Definition of the string, which stores the multigrid cycling
    ! strategy.
    !

    integer(kind=intType),constant :: nCycles, nCyclesCoarse
    integer(kind=intType),constant :: nSaveVolume, nSaveSurface
    integer(kind=intType),constant :: nsgStartup, smoother, nRKStages
    integer(kind=intType),constant :: nSubiterations
    integer(kind=intType),constant :: nSubIterTurb, nUpdateBleeds
    integer(kind=intType),constant :: resAveraging
    real(kind=realType),constant :: CFLLimit
    integer(kind=intType),constant :: turbTreatment, turbSmoother, turbRelax
    integer(kind=intType),constant :: mgBoundCorr, mgStartlevel
    integer(kind=intType),constant :: nMGSteps, nMGLevels
    real(kind=realType),constant :: timeLimit
    integer(kind=intType), allocatable, dimension(:) :: cycleStrategy
    integer(kind=intType),constant :: miniterNum
    real(kind=realType),constant :: cfl, cflCoarse, fcoll, smoop
    real(kind=realType),constant :: alfaTurb, betaTurb
    real(kind=realType),constant :: L2Conv, L2ConvCoarse
    real(kind=realType),constant :: L2ConvRel
    real(kind=realType),constant :: maxL2DeviationFactor
    real(kind=realType),constant :: relaxBleeds
    real(kind=realtype),constant :: epscoefconv
    integer(kind=inttype),constant :: convcheckwindowsize
    real(kind=realType), allocatable, dimension(:) :: etaRK, cdisRK
    character(len=maxStringLen),constant :: mgDescription
    logical,constant :: rkReset
    logical,constant :: useLinResMonitor
    logical,constant :: freezeTurbSource
    logical,constant :: printIterations
    logical,constant :: printWarnings
    logical,constant :: printNegativeVolumes
    real(kind=realType),constant, dimension(4) :: turbResScale
    logical,constant :: useDissContinuation
    real(kind=realType),constant :: dissContMagnitude, dissContMidpoint, dissContSharpness

end module cudaInputIteration

module cudaInputPhysics
    !
    !       Input parameters which are related to the physics of the flow,
    !       like governing equations, mode of the equations, turbulence
    !       model and free stream conditions.
    !
    use precision
    implicit none
    save

    !       Definition of the physics input parameters.
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
    ! useQCR:              Determines if the QCR term is applied to the shear tensor computation
    !                      when considering turbulence model effects
    ! useRotationSA:       Determines if we will use rotation correction (SA model only)
    ! useft2SA:            Determines if we will use the ft2 term (SA model only)
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
    ! SSuthDim:            Sutherlands law temperature (SI Units)
    ! muSuthDim:           Reference viscosity at reference temperature for Sutherlands law (SI Units)
    ! TSuthDim:            Reference temperature for Sutherlands law (SI Units)
    ! momentAxis(3,2)      Axis about which to calculate a moment, provided as 2 points in 3-D
    ! cavitationnumber     Negative Cp value that triggers the traditional
    !                      step-function based cavitation sensor.
    ! cpmin_rho            The rho parameter used with the KS-based cavitation sensor.
    ! cpmin_family         The cpmin for a given surface family that does not use
    !                      KS-aggregation, but rather an exact min computation.

    integer(kind=intType),constant :: equations, equationMode, flowType
    integer(kind=intType),constant :: turbModel, cpModel, turbProd
    integer(kind=intType),constant :: rvfN
    logical,constant :: rvfB
    logical,constant :: useQCR, useRotationSA, useft2SA

    logical,constant :: wallFunctions, wallDistanceNeeded

    real(kind=realType),constant :: alpha, beta
    integer(kind=intType),constant :: liftIndex
    real(kind=realType),constant :: Mach, MachCoef, MachGrid
    real(kind=realType),constant :: Reynolds, ReynoldsLength
    real(kind=realType),constant :: gammaConstant, RGasDim
    real(kind=realType),constant :: Prandtl, PrandtlTurb, pklim, wallOffset, wallDistCutoff
    real(kind=realType),constant :: eddyVisInfRatio, turbIntensityInf
    real(kind=realType),constant :: surfaceRef, lengthRef
    real(kind=realType),constant, dimension(3) :: velDirFreestream
    real(kind=realType),constant, dimension(3) :: liftDirection
    real(kind=realType),constant, dimension(3) :: dragDirection
    real(kind=realType),constant, dimension(3) :: pointRef
    real(kind=realType),constant, dimension(3, 2) :: momentAxis
    real(kind=realType),constant :: SSuthDim, muSuthDim, TSuthDim
    real(kind=realType),constant :: cavitationnumber
    real(kind=realType),constant :: cpmin_rho
    real(kind=realType), dimension(:), allocatable :: cpmin_family

    real(kind=realType),constant, dimension(3) :: pointRefEC

    ! Return forces as tractions instead of forces:
    logical,constant:: forcesAsTractions

end module cudaInputPhysics

module cudaParamTurb
    use precision, only: realType, intType
    implicit none
    save
    real(kind=realType), constant :: rvfLimitK, rvfLimitE, rvfCl
    real(kind=realType), constant :: rvfCmu

    !
!       Variables to store the parameters for the wall functions fits.
!       As these variables depend on the turbulence model they are set
!       during runtime. Allocatables are used, because the number of
!       fits could be different for the different models.
!       The curve is divided in a number of intervals and is
!       constructed such that both the function and the derivatives
!       are continuous. Consequently cubic polynomials are used.
!
    ! nFit:               Number of intervals of the curve.
    ! ypT(0:nFit):        y+ values at the interval boundaries.
    ! reT(0:nFit):        Reynolds number at the interval
    !                     boundaries, where the Reynolds number is
    !                     defined with the local velocity and the
    !                     wall distance.
    ! up0(nFit):          Coefficient 0 in the fit for the
    !                     nondimensional tangential velocity as a
    !                     function of the Reynolds number.
    ! up1(nFit):          Idem for coefficient 1.
    ! up2(nFit):          Idem for coefficient 2.
    ! up3(nFit):          Idem for coefficient 3.
    ! tup0(nFit,nt1:nt2): Coefficient 0 in the fit for the
    !                     nondimensional turbulence variables as a
    !                     function of y+.
    ! tup1(nFit,nt1:nt2): Idem for coefficient 1.
    ! tup2(nFit,nt1:nt2): Idem for coefficient 2.
    ! tup3(nFit,nt1:nt2): Idem for coefficient 3.
    ! tuLogFit(nt1:nt2):  Whether or not the logarithm of the variable
    !                     has been fitted.
    integer(kind=intType), constant :: nFit
    real(kind=realType), dimension(:), allocatable, device :: ypT, reT
    real(kind=realType), dimension(:), allocatable, device :: up0, up1
    real(kind=realType), dimension(:), allocatable, device :: up2
    
    real(kind=realType), dimension(:, :), allocatable,device :: tup0, tup1
    real(kind=realType), dimension(:, :), allocatable,device :: tup2, tup3

end module cudaParamTurb


module cudaIteration
    use precision, only: realType, intType
    implicit none
    save
    ! groundLevel:  Current ground level of the computation. Needed
    !               to determine what kind of action must be
    !               undertaken. E.G. On the coarse grids no solution
    !               will be written.
    ! currentLevel: MG level at which the compution currently resides.

    integer(kind=intType),device :: groundLevel, currentLevel
    !rFil : coefficient to control the fraction of the dissipation
    !        residual of the previous runge-kutta stage.
    real(kind=realType), device :: rFil
    ! Variables for monitoring the residuals
    real(kind=realType),device :: totalR0, totalR
end module cudaIteration

module cudaTurbMod
    implicit none
    save
    ! secondOrd:  whether or not a second order discretization for
    !             the advective terms must be used.
    logical,device :: secondOrd
end module cudaTurbMod
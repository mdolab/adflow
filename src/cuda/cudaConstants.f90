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
    contains
    subroutine cudaCopyFlowVarRefState
        use flowVarRefState, only: h_nw=> nw, h_nwf=> nwf, h_nwt=> nwt, h_nt1=> nt1, h_nt2=> nt2, &
                                 h_pRef=> pRef, h_rhoRef=> rhoRef, h_TRef=> TRef, h_muRef=> muRef, &
                                    h_timeRef=> timeRef, h_uRef=> uRef, h_href=> href, &
                                    h_LRef=> LRef, h_LRefSpecified=> LRefSpecified, &
                                    h_pInfDim=> pInfDim, h_rhoInfDim=> rhoInfDim, h_muDim=> muDim, &
                                    h_TinfDim=> TinfDim, h_muInfDim=> muInfDim, &
                                    h_wInf=> wInf, h_pInf=> pInf, h_pInfCorr=> pInfCorr, &
                                    h_rhoInf=> rhoInf, h_uInf=> uInf, h_muInf=> muInf, &
                                    h_RGas=> RGas, h_gammaInf=> gammaInf, &
                                    h_viscous=> viscous, h_kPresent=> kPresent, h_eddyModel=> eddyModel
        
        
        implicit none
        !copy data to gpu
        nw = h_nw
        nwf = h_nwf
        nwt = h_nwt
        nt1 = h_nt1
        nt2 = h_nt2
        pRef = h_pRef
        rhoRef = h_rhoRef
        TRef = h_TRef
        muRef = h_muRef
        timeRef = h_timeRef
        uRef = h_uRef
        href = h_href
        LRef = h_LRef
        LRefSpecified = h_LRefSpecified
        pInfDim = h_pInfDim
        rhoInfDim = h_rhoInfDim
        muDim = h_muDim
        TinfDim = h_TinfDim
        muInfDim = h_muInfDim
        wInf(1:h_nw) = h_wInf(1:h_nw)
        pInf = h_pInf
        pInfCorr = h_pInfCorr
        rhoInf = h_rhoInf
        uInf = h_uInf
        muInf = h_muInf
        RGas = h_RGas
        gammaInf = h_gammaInf
        viscous = h_viscous
        kPresent = h_kPresent
        eddyModel = h_eddyModel


    end subroutine cudaCopyFlowVarRefState

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
    contains
    subroutine cudaCopyInputDiscretization
        use inputDiscretization, only: h_spaceDiscr=>spaceDiscr, h_spaceDiscrCoarse=>spaceDiscrCoarse, &
                                        h_orderTurb=>orderTurb, h_riemann=>riemann, &
                                        h_riemannCoarse=>riemannCoarse, h_limiter=>limiter, &
                                        h_precond=>precond, h_eulerWallBCTreatment=>eulerWallBCTreatment, &
                                        h_viscWallBCTreatment=>viscWallBCTreatment, h_outflowTreatment=>outflowTreatment, &
                                        h_nonMatchTreatment=>nonMatchTreatment, h_vis2=>vis2, &
                                        h_vis4=>vis4, h_vis2Coarse=>vis2Coarse, h_adis=>adis, &
                                        h_kappaCoef=>kappaCoef, h_vortexCorr=>vortexCorr, &
                                        h_dirScaling=>dirScaling, h_hScalingInlet=>hScalingInlet, &
                                        h_radiiNeededFine=>radiiNeededFine, h_radiiNeededCoarse=>radiiNeededCoarse, &
                                        h_lumpedDiss=>lumpedDiss, h_approxSA=>approxSA, h_sigma=>sigma, &
                                        h_useApproxWallDistance=>useApproxWallDistance, &
                                        h_lowSpeedPreconditioner=>lowSpeedPreconditioner, &
                                        h_acousticScaleFactor=>acousticScaleFactor
        implicit none
        spaceDiscr = h_spaceDiscr
        spaceDiscrCoarse = h_spaceDiscrCoarse
        orderTurb = h_orderTurb
        riemann = h_riemann
        riemannCoarse = h_riemannCoarse
        limiter = h_limiter
        precond = h_precond
        eulerWallBCTreatment = h_eulerWallBCTreatment
        viscWallBCTreatment = h_viscWallBCTreatment
        outflowTreatment = h_outflowTreatment
        nonMatchTreatment = h_nonMatchTreatment
        vis2 = h_vis2
        vis4 = h_vis4
        vis2Coarse = h_vis2Coarse
        adis = h_adis
        kappaCoef = h_kappaCoef
        vortexCorr = h_vortexCorr
        dirScaling = h_dirScaling
        hScalingInlet = h_hScalingInlet
        radiiNeededFine = h_radiiNeededFine
        radiiNeededCoarse = h_radiiNeededCoarse
        lumpedDiss = h_lumpedDiss
        approxSA = h_approxSA
        sigma = h_sigma
        useApproxWallDistance = h_useApproxWallDistance
        lowSpeedPreconditioner = h_lowSpeedPreconditioner
        acousticScaleFactor = h_acousticScaleFactor
    end subroutine cudaCopyInputDiscretization

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

    contains
    subroutine cudaCopyInputIteration
        use inputIteration, only: h_nCycles=>nCycles, h_nCyclesCoarse=>nCyclesCoarse, h_nSaveVolume=>nSaveVolume, h_nSaveSurface=>nSaveSurface,&
                                 h_nsgStartup=>nsgStartup, h_smoother=>smoother, h_nRKStages=>nRKStages, h_nSubiterations=>nSubiterations, &
                                    h_nSubIterTurb=>nSubIterTurb, h_nUpdateBleeds=>nUpdateBleeds, h_resAveraging=>resAveraging,  h_CFLLimit=>CFLLimit, &
                                    h_turbTreatment=>turbTreatment, h_turbSmoother=>turbSmoother, h_turbRelax=>turbRelax, &
                                    h_mgBoundCorr=>mgBoundCorr, h_mgStartlevel=>mgStartlevel, h_nMGSteps=>nMGSteps, &
                                    h_nMGLevels=>nMGLevels, h_timelimit=>timeLimit,h_cycleStrategy=>cycleStrategy, h_miniterNum=>miniterNum, &
                                    h_cfl=>cfl, h_cflCoarse=>cflCoarse, h_fcoll=>fcoll, &
                                    h_smoop=>smoop, h_alfaTurb=>alfaTurb, h_betaTurb=>betaTurb, h_L2Conv=>L2Conv, &
                                    h_L2ConvCoarse=>L2ConvCoarse, h_L2ConvRel=>L2ConvRel, h_maxL2DeviationFactor=>maxL2DeviationFactor, &

                                    h_relaxBleeds=>relaxBleeds, h_epscoefconv=>epscoefconv, h_convcheckwindowsize=>convcheckwindowsize, &
                                    h_etaRK=>etaRK, h_cdisRK=>cdisRK, h_mgDescription=>mgDescription, h_rkReset=>rkReset, &
                                    h_useLinResMonitor=>useLinResMonitor, h_freezeTurbSource=>freezeTurbSource, &

                                    h_printIterations=>printIterations, h_printWarnings=>printWarnings, &
                                    h_printNegativeVolumes=>printNegativeVolumes, h_turbResScale=>turbResScale, &
                                    h_useDissContinuation=>useDissContinuation, h_dissContMagnitude=>dissContMagnitude, &
                                    h_dissContMidpoint=>dissContMidpoint, h_dissContSharpness=>dissContSharpness

        implicit none
        nCycles = h_nCycles
        nCyclesCoarse = h_nCyclesCoarse
        nSaveVolume = h_nSaveVolume
        nSaveSurface = h_nSaveSurface
        nsgStartup = h_nsgStartup
        smoother = h_smoother
        nRKStages = h_nRKStages
        nSubiterations = h_nSubiterations
        nSubIterTurb = h_nSubIterTurb
        nUpdateBleeds = h_nUpdateBleeds
        resAveraging = h_resAveraging
        CFLLimit = h_CFLLimit
        turbTreatment = h_turbTreatment
        turbSmoother = h_turbSmoother
        turbRelax = h_turbRelax
        mgBoundCorr = h_mgBoundCorr
        mgStartlevel = h_mgStartlevel
        nMGSteps = h_nMGSteps
        nMGLevels = h_nMGLevels
        timeLimit = h_timelimit
        !TODO allocate this
        if (.not. allocated(cycleStrategy)) then
            allocate(cycleStrategy(nMGSteps))
        end if
        cycleStrategy(1:nMGSteps) = h_cycleStrategy(1:nMGSteps)
        miniterNum = h_miniterNum
        cfl = h_cfl
        cflCoarse = h_cflCoarse
        fcoll = h_fcoll
        smoop = h_smoop
        alfaTurb = h_alfaTurb
        betaTurb = h_betaTurb
        L2Conv = h_L2Conv
        L2ConvCoarse = h_L2ConvCoarse
        L2ConvRel = h_L2ConvRel
        maxL2DeviationFactor = h_maxL2DeviationFactor
        relaxBleeds = h_relaxBleeds
        epscoefconv = h_epscoefconv
        convcheckwindowsize = h_convcheckwindowsize
        !TO DO ALLOCATE THIS DATA
        if (.not. allocated(etaRk)) then
            allocate(etaRk(nRKStages))
        end if
        etaRK(1:nRKStages) = h_etaRK(1:nRKStages)
        if (.not. allocated(cdisRK)) then
            allocate(cdisRK(nRKStages))
        end if
        cdisRK(1:nRKStages) = h_cdisRK(1:nRKStages)
        mgDescription = h_mgDescription
        rkReset = h_rkReset
        useLinResMonitor = h_useLinResMonitor
        freezeTurbSource = h_freezeTurbSource
        printIterations = h_printIterations
        printWarnings = h_printWarnings
        printNegativeVolumes = h_printNegativeVolumes
        turbResScale(1:4) = h_turbResScale(1:4)
        useDissContinuation = h_useDissContinuation
        dissContMagnitude = h_dissContMagnitude
        dissContMidpoint = h_dissContMidpoint
        dissContSharpness = h_dissContSharpness
    end subroutine cudaCopyInputIteration
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
    contains
    subroutine cudaCopyInputPhysics
        use inputPhysics, only: h_equations=>equations, h_equationMode=>equationMode, h_flowType=>flowType, &
                                h_turbModel=>turbModel, h_cpModel=>cpModel, h_turbProd=>turbProd, &
                                h_rvfN=>rvfN, h_rvfB=>rvfB, h_useQCR=>useQCR, h_useRotationSA=>useRotationSA, &
                                h_useft2SA=>useft2SA, h_wallFunctions=>wallFunctions, h_wallDistanceNeeded=>wallDistanceNeeded, &
                                h_alpha=>alpha, h_beta=>beta, h_liftIndex=>liftIndex, h_Mach=>Mach, &
                                h_MachCoef=>MachCoef, h_MachGrid=>MachGrid, h_Reynolds=>Reynolds, &
                                h_ReynoldsLength=>ReynoldsLength, h_gammaConstant=>gammaConstant, &
                                h_RGasDim=>RGasDim, h_Prandtl=>Prandtl, h_PrandtlTurb=>PrandtlTurb, &
                                h_pklim=>pklim, h_wallOffset=>wallOffset, h_wallDistCutoff=>wallDistCutoff, &
                                h_eddyVisInfRatio=>eddyVisInfRatio, h_turbIntensityInf=>turbIntensityInf, &
                                h_surfaceRef=>surfaceRef, h_lengthRef=>lengthRef, h_velDirFreestream=>velDirFreestream, &
                                h_liftDirection=>liftDirection, h_dragDirection=>dragDirection, h_pointRef=>pointRef, &
                                h_momentAxis=>momentAxis, h_SSuthDim=>SSuthDim, h_muSuthDim=>muSuthDim, &
                                h_TSuthDim=>TSuthDim, h_cavitationnumber=>cavitationnumber, h_cpmin_rho=>cpmin_rho, &
                                h_cpmin_family=>cpmin_family, h_pointRefEC=>pointRefEC, h_forcesAsTractions=>forcesAsTractions

        implicit none
        equations = h_equations
        equationMode = h_equationMode
        flowType = h_flowType
        turbModel = h_turbModel
        cpModel = h_cpModel
        turbProd = h_turbProd
        rvfN = h_rvfN
        rvfB = h_rvfB
        useQCR = h_useQCR
        useRotationSA = h_useRotationSA
        useft2SA = h_useft2SA
        wallFunctions = h_wallFunctions
        wallDistanceNeeded = h_wallDistanceNeeded
        alpha = h_alpha
        beta = h_beta
        liftIndex = h_liftIndex
        Mach = h_Mach
        MachCoef = h_MachCoef
        MachGrid = h_MachGrid
        Reynolds = h_Reynolds
        ReynoldsLength = h_ReynoldsLength
        gammaConstant = h_gammaConstant
        RGasDim = h_RGasDim
        Prandtl = h_Prandtl
        PrandtlTurb = h_PrandtlTurb
        pklim = h_pklim
        wallOffset = h_wallOffset
        wallDistCutoff = h_wallDistCutoff
        eddyVisInfRatio = h_eddyVisInfRatio
        turbIntensityInf = h_turbIntensityInf
        surfaceRef = h_surfaceRef
        lengthRef = h_lengthRef
        velDirFreestream(1:3) = h_velDirFreestream(1:3)
        liftDirection(1:3) = h_liftDirection(1:3)
        dragDirection(1:3) = h_dragDirection(1:3)
        pointRef(1:3) = h_pointRef(1:3)
        momentAxis(1:3, 1:2) = h_momentAxis(1:3, 1:2)
        SSuthDim = h_SSuthDim
        muSuthDim = h_muSuthDim
        TSuthDim = h_TSuthDim
        cavitationnumber = h_cavitationnumber
        cpmin_rho = h_cpmin_rho
        !TODO allocate this data
        if (.not. allocated(cpmin_family)) then
            allocate(cpmin_family(1:size(h_cpmin_family)))
        end if
        cpmin_family(1:size(h_cpmin_family)) = h_cpmin_family(1:size(h_cpmin_family))
        pointRefEC(1:3) = h_pointRefEC(1:3)
        forcesAsTractions = h_forcesAsTractions
    end subroutine cudaCopyInputPhysics

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
    logical, dimension(:), allocatable,device :: tuLogFit
    contains
    subroutine cudaCopyParamTurb
        use paramTurb, only: h_rvfLimitK=>rvflimitK, h_rvfLimitE=>rvfLimitE, h_rvfCl=>rvfCl, &
                            h_rvfCmu=>rvfCmu, h_nFit=>nFit, h_ypT=>ypT, h_reT=>reT, &
                            h_up0=>up0, h_up1=>up1, h_up2=>up2, h_tup0=>tup0, &
                            h_tup1=>tup1, h_tup2=>tup2, h_tup3=>tup3, h_tuLogFit=>tuLogFit
        
        use flowVarRefState, only: nt1, nt2
        implicit none
        h_rvfLimitK = rvfLimitK
        h_rvfLimitE = rvfLimitE
        h_rvfCl = rvfCl
        h_rvfCmu = rvfCmu
        h_nFit = nFit
        !TOOD allocate these and fix copy for proper dimensions
        if (.not. allocated(ypT)) then
            allocate(ypT(0:h_nFit))
        end if
        ypT(0:h_nFit) = h_ypT(0:h_nFit)
        
        if (.not. allocated(reT)) then
            allocate(reT(0:h_nFit))
        end if
        reT(0:h_nFit) = h_reT(0:h_nFit)

        if(.not. allocated(up0)) then
            allocate(up0(h_nFit))
        end if
        up0(1:h_nFit) = h_up0(1:h_nFit)

        if(.not. allocated(up1)) then
            allocate(up1(h_nFit))
        end if
        up1(1:h_nFit) = h_up1(1:h_nFit)
        
        if(.not. allocated(up2)) then
            allocate(up2(h_nFit))
        end if
        up2(1:h_nFit) = h_up2(1:h_nFit)

        if (.not. allocated(tup0)) then
            allocate(tup0(1:h_nFit, nt1:nt2))
        end if
        tup0(1:h_nFit, nt1:nt2) = h_tup0(1:h_nFit, nt1:nt2)

        if (.not. allocated(tup1)) then
            allocate(tup1(1:h_nFit, nt1:nt2))
        end if
        tup1(1:h_nFit, nt1:nt2) = h_tup1(1:h_nFit, nt1:nt2)

        if (.not. allocated(tup2)) then
            allocate(tup2(1:h_nFit, nt1:nt2))
        end if
        tup2(1:h_nFit, nt1:nt2) = h_tup2(1:h_nFit, nt1:nt2)

        if (.not. allocated(tup3)) then
            allocate(tup3(1:h_nFit, nt1:nt2))
        end if
        tup3(1:h_nFit, nt1:nt2) = h_tup3(1:h_nFit, nt1:nt2)
        
        if (.not. allocated(tuLogFit)) then
            allocate(tuLogFit(nt1:nt2))
        end if
        tuLogFit(nt1:nt2) = h_tuLogFit(nt1:nt2)
    end subroutine cudaCopyParamTurb
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
    contains
    subroutine cudaCopyIteration
        use iteration, only: h_groundLevel=>groundLevel, h_currentLevel=>currentLevel, &
                             h_rFil=>rFil, h_totalR0=>totalR0, h_totalR=>totalR
        implicit none
        groundLevel = h_groundLevel
        currentLevel = h_currentLevel
        rFil = h_rFil
        totalR0 = h_totalR0
        totalR = h_totalR
    end subroutine cudaCopyIteration

end module cudaIteration

module cudaTurbMod
    implicit none
    save
    ! secondOrd:  whether or not a second order discretization for
    !             the advective terms must be used.
    logical,device :: secondOrd
    contains
    subroutine cudaCopyTurbMod
        use turbMod, only: h_secondOrd=>secondOrd
        implicit none
        secondOrd = h_secondOrd
    end subroutine cudaCopyTurbMod
end module cudaTurbMod
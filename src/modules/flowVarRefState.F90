module flowVarRefState
!
!       Module that contains information about the reference state as
!       well as the nondimensional free stream state.
!
    use constants, only: intType, realType
    implicit none
    save

    ! nw:       Total number of independent variables including the
    !           turbulent variables.
    ! nwf:      Number of flow variables. For perfect gas computations
    !           this is 5.
    ! nwt:      Number of turbulent variables, nwt = nw - nwf
    ! nt1, nt2: Initial and final indices for turbulence variables

    integer(kind=intType) :: nw, nwf, nwt, nt1, nt2

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

    real(kind=realType) :: pRef, rhoRef, TRef
    real(kind=realType) :: muRef, timeRef, uRef, href

    ! LRef:          Conversion factor of the length unit of the
    !                grid to meter. e.g. if the grid is in mm.,
    !                LRef = 1.e-3.
    ! LRefSpecified: Whether or not a conversion factor is specified
    !                in the input file.

    real(kind=realType) :: LRef
    logical :: LRefSpecified

    ! pInfDim:   Free stream pressure in Pa.
    ! rhoInfDim: Free stream density in kg/m^3.
    ! muDim:     Free stream molecular viscosity in kg/(m s)

    real(kind=realType) :: pInfDim, rhoInfDim, muDim, TinfDim, muInfDim
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

    real(kind=realType) :: rhoInf, uInf, pInf, pInfCorr
    real(kind=realType) :: RGas, muInf, gammaInf
#ifdef USE_TAPENADE
    real(kind=realType), dimension(10) :: wInf
#else
    real(kind=realType), dimension(:), allocatable :: wInf
#endif

#ifndef USE_TAPENADE
    REAL(kind=realtype) :: prefd, rhorefd, trefd, uRefd, Hrefd
    REAL(kind=realtype) :: murefd, timerefd
    REAL(kind=realtype) :: pinfdimd, rhoinfdimd, tinfdimd
    real(kind=realtype) :: mudimd, muinfdimd
    REAL(kind=realtype) :: pinfdimb, rhoinfdimb
    REAL(kind=realtype) :: rhoinfd, uinfd, pinfd, pinfcorrd
    REAL(kind=realtype) :: rgasd, muinfd, gammainfd
    real(kind=realType), dimension(:), allocatable :: wInfd, wInfb
#endif
    ! viscous:   whether or not this is a viscous computation.
    ! kPresent:  whether or not a turbulent kinetic energy is present
    !            in the turbulence model.
    ! eddyModel: whether or not the turbulence model is an eddy
    !            viscosity model.

    logical :: kPresent, eddyModel, viscous

end module flowVarRefState

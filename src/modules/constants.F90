module constants
  
  ! Define all constants used in the code. This is the *only* module 
  ! that is allowed to be imported without an 'only' qualifier. 

  use precision
  implicit none
  save

  ! Maximum numbers of characters in a string and in a cgns name

  integer, parameter :: maxStringLen   = 256
  integer, parameter :: maxCGNSNameLen =  32

  ! Numerical constants

  real(kind=realType), parameter :: pi    = 3.1415926535897931_realType
  real(kind=realType), parameter :: eps   = 1.e-25_realType
  real(kind=realType), parameter :: large = 1.e+37_realType

  ! Constants to define the porosity values

  integer(kind=porType), parameter :: noFlux     = -1_porType
  integer(kind=porType), parameter :: boundFlux  =  0_porType
  integer(kind=porType), parameter :: normalFlux =  1_porType

  ! Indices in the array of independent variables

  integer, parameter :: irho    =  1  ! Density
  integer, parameter :: ivx     =  2  ! x-Velocity
  integer, parameter :: ivy     =  3  ! y-velocity
  integer, parameter :: ivz     =  4  ! z-Velocity
  integer, parameter :: irhoE   =  5  ! Energy

  integer, parameter :: itu1    =  6  ! Turbulent kinetic energy,
  ! SA viscosity
  integer, parameter :: itu2    =  7  ! Dissipation rate, time scale 
  integer, parameter :: itu3    =  8  ! Scalar V2 
  integer, parameter :: itu4    =  9  ! Scalar F2 
  integer, parameter :: itu5    = 10  ! Eddy-viscosity used for
  ! wall functions.

  ! Parameters to indicate the position in the work array dw for
  ! turbulence models.

  integer, parameter :: idvt    = 1  ! Tmp RHS storage; at max a
  ! 2x2 subsystem is solved.
  integer, parameter :: ivort   = 3  ! Tmp vort storage
  integer, parameter :: istrain = 3  ! Tmp strain storage
  integer, parameter :: iprod   = 3  ! Tmp prod storage
  integer, parameter :: icd     = 4  ! Tmp cross term storage
  integer, parameter :: if1SST  = 5  ! Tmp F1 (for SST) storage
  integer, parameter :: isct    = 4  ! Tmp time scale (for v2f) storage
  integer, parameter :: iscl2   = 5  ! Tmp length scale (for v2f) storage
  integer, parameter :: iqq     = 6  ! Central jacobian storage

  ! Indices in the array of conservative flow residuals for the
  ! momentum variables.

  integer, parameter :: imx   = ivx ! x-Momentum
  integer, parameter :: imy   = ivy ! y-Momentum
  integer, parameter :: imz   = ivz ! z-Momentum

  ! Floating point parameters.

  real(kind=realType), parameter :: zero  = 0.0_realType
  real(kind=realType), parameter :: one   = 1.0_realType
  real(kind=realType), parameter :: two   = 2.0_realType
  real(kind=realType), parameter :: three = 3.0_realType
  real(kind=realType), parameter :: four  = 4.0_realType
  real(kind=realType), parameter :: five  = 5.0_realType
  real(kind=realType), parameter :: six   = 6.0_realType
  real(kind=realType), parameter :: eight = 8.0_realType

  real(kind=realType), parameter :: half   = 0.5_realType
  real(kind=realType), parameter :: third  = one/three
  real(kind=realType), parameter :: fourth = 0.25_realType
  real(kind=realType), parameter :: sixth  = one/six
  real(kind=realType), parameter :: eighth = 0.125_realType
  real(kind=realType), parameter :: threefourth = 0.75_realType
  real(kind=realType), parameter :: sqrtthree = 1.7320508075688772_realType

  ! String constants
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

  ! Threshold parameter for real types; the value depends
  ! whether single or double precision is used.

#ifdef USE_SINGLE_PRECISION
  real(kind=realType), parameter :: thresholdReal = 1.e-5_realType
#else
  real(kind=realType), parameter :: thresholdReal = 1.e-10_realType
#endif

  ! Definition of the tab and carriage return character.
#ifndef USE_TAPENADE
  character(len=1), parameter :: tabChar = achar(9)
  character(len=1), parameter :: retChar = achar(13)
#endif

  integer(kind=intType), parameter :: &
       EulerEquations = 1,  &
       NSEquations    = 2,  &
       RANSEquations  = 3

  integer(kind=intType), parameter :: &
       steady        = 1,   &
       unsteady      = 2,   &
       timeSpectral  = 3

  integer(kind=intType), parameter :: &
       internalFlow = 1,    &
       externalFlow = 2

  integer(kind=intType), parameter :: &
       cpConstant      = 1, &
       cpTempCurveFits = 2

  integer(kind=intType), parameter ::                              &
       spalartAllmaras        =  2,  &
       spalartAllmarasEdwards =  3,  &
       komegaWilcox           =  4,  &
       komegaModified         =  5,  &
       ktau                   =  6,  &
       menterSST              =  7,  &
       v2f                    = 10

  integer(kind=intType), parameter :: &
       strain       = 1,    &
       vorticity    = 2,    &
       katoLaunder  = 3

  integer(kind=intType), parameter :: &
       firstOrder  = 1, &
       secondOrder = 2, &
       thirdOrder  = 3, &
       fourthOrder = 4, &
       fifthOrder  = 5

  integer(kind=intType), parameter :: &
       dissScalar = 1,  &
       dissMatrix = 2,  &
       dissCusp   = 3,  &
       upwind     = 9

  integer(kind=intType), parameter :: &
       Roe     = 1,     &
       vanLeer = 2,     &
       ausmdv  = 3

  integer(kind=intType), parameter :: &
       noLimiter  = 2,  &
       vanAlbeda  = 3,  &
       minmod     = 4

  integer(kind=intType), parameter :: &
       noPrecond  = 1,  &
       Turkel     = 2,  &
       ChoiMerkle = 3

  integer(kind=intType), parameter ::                          &
       constantPressure     = 1, &
       linExtrapolPressure  = 2, &
       quadExtrapolPressure = 3, &
       normalMomentum       = 4

  integer(kind=intType), parameter ::  &
       constantExtrapol = 1, &
       linExtrapol      = 2

  integer(kind=intType), parameter :: &
       NonConservative = 1, &
       Conservative    = 2

  integer(kind=intType), parameter :: &
       precisionSingle = 1, &
       precisionDouble = 2

  ! Definition of the parameters for the time integration scheme.
  integer(kind=intType), parameter :: &
       BDF        = 1, &
       explicitRK = 2, &
       implicitRK = 3, &
       MD         = 4

  ! Line search parameters
  integer(kind=intType), parameter :: &
       noLineSearch = 0_intType, &
       cubicLineSearch = 1_intType, &
       nonMonotoneLineSearch = 2_intType

  integer(kind=intType), parameter :: &
       RungeKutta  = 1,  &
       DADI        = 2,  &
       nlLusgs     = 3,  &
       nlLusgsLine = 4

  integer(kind=intType), parameter :: &
       segregated = 1,   &
       coupled    = 2
  integer(kind=intType), parameter :: &
       gmres = 1,        &
       adi   = 2

  integer(kind=intType), parameter :: &
       bcDirichlet0 = 0, &
       bcNeumann0    = 1

  integer(kind=intType), parameter :: &
       noResAveraging        = 0, &
       alwaysResAveraging    = 1, &
       alternateResAveraging = 2

  integer(kind=intType), parameter :: &
       turbRelaxNotDefined = 0,  &
       turbRelaxExplicit   = 1,  &
       turbRelaxImplicit   = 2

  ! Parameters used for coarsening definition.
  integer(kind=porType), parameter :: &
       leftStarted  = -1_porType, &
       regular      =  0_porType, &
       rightStarted =  1_porType

  ! Parameters used for subsonic inlet bc treatment.
  integer(kind=intType), parameter :: &
       noSubInlet      = 0, &
       totalConditions = 1, &
       massFlow        = 2
  !
  integer, parameter :: adtSurfaceADT = 1
  integer, parameter :: adtVolumeADT  = 2
  integer(kind=intType), parameter :: nCoorMaxLowerLimit = 100000

  integer(kind=adtElementType), parameter :: adtTriangle      = 1
  integer(kind=adtElementType), parameter :: adtQuadrilateral = 2
  integer(kind=adtElementType), parameter :: adtTetrahedron   = 3
  integer(kind=adtElementType), parameter :: adtPyramid       = 4
  integer(kind=adtElementType), parameter :: adtPrism         = 5
  integer(kind=adtElementType), parameter :: adtHexahedron    = 6
  
  ! BCDefinitions
  integer(kind=intType), parameter :: BCNull                =   0
  integer(kind=intType), parameter :: Symm                  =  -1
  integer(kind=intType), parameter :: SymmPolar             =  -2
  integer(kind=intType), parameter :: NSWallAdiabatic       =  -3
  integer(kind=intType), parameter :: NSWallIsothermal      =  -4
  integer(kind=intType), parameter :: EulerWall             =  -5
  integer(kind=intType), parameter :: FarField              =  -6
  integer(kind=intType), parameter :: SupersonicInflow      =  -7
  integer(kind=intType), parameter :: SubsonicInflow        =  -8
  integer(kind=intType), parameter :: SupersonicOutflow     =  -9
  integer(kind=intType), parameter :: SubsonicOutflow       = -10
  integer(kind=intType), parameter :: MassBleedInflow       = -11
  integer(kind=intType), parameter :: MassBleedOutflow      = -12
  integer(kind=intType), parameter :: mDot                  = -13
  integer(kind=intType), parameter :: Thrust                = -14
  integer(kind=intType), parameter :: Extrap                = -15
  integer(kind=intType), parameter :: B2BMatch              = -16
  integer(kind=intType), parameter :: B2BMismatch           = -17
  integer(kind=intType), parameter :: SlidingInterface      = -18
  integer(kind=intType), parameter :: OversetOuterBound     = -19
  integer(kind=intType), parameter :: DomainInterfaceAll    = -20
  integer(kind=intType), parameter :: DomainInterfaceRhoUVW = -21
  integer(kind=intType), parameter :: DomainInterfaceP      = -22
  integer(kind=intType), parameter :: DomainInterfaceRho    = -23
  integer(kind=intType), parameter :: DomainInterfaceTotal  = -24
  integer(kind=intType), parameter :: BCNotValid            = -25
  !
  !      Number of actual boundary conditions supported by the code
  !      This number refers to bocos, not flow-through BCs
  !      Edit this number when additional boundary conditions are
  !      supported
  !
  integer(kind=intType), parameter :: nBCs = 24

    !Block faces on which boundary conditions may be imposed
  integer(kind=intType), parameter :: iMin = 1
  integer(kind=intType), parameter :: iMax = 2
  integer(kind=intType), parameter :: jMin = 3
  integer(kind=intType), parameter :: jMax = 4
  integer(kind=intType), parameter :: kMin = 5
  integer(kind=intType), parameter :: kMax = 6

  integer(kind=intType) :: myIntStack(32)
  integer(kind=intType) :: myIntPtr = 0

  ! BC specific input variable counts
  real(kind=realType), parameter :: nbcVarSubsonicInflow = 17
  real(kind=realType), parameter :: nbcVarSubsonicOutflow = 1
  real(kind=realType), parameter :: nbcVarSupersonicInflow = 7
  real(kind=realType), parameter :: nbcVarIsothermalWall = 1


end module constants

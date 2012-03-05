!
!      ******************************************************************
!      *                                                                *
!      * File:          su_cgns.F90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-16-2005                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module su_cgns
!
!      ******************************************************************
!      *                                                                *
!      * Module that contains the definition of the cgns parameters.    *
!      * Depending on the compiler flags either the file cgnslib_f.h is *
!      * included or the functionality is faked by just defining the    *
!      * parameters.                                                    *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       ! CGNS is not used. define the parameters.

       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Modes for cgns file                                            *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: MODE_READ   = 0
       integer, parameter :: MODE_WRITE  = 1
       integer, parameter :: MODE_CLOSED = 2
       integer, parameter :: MODE_MODIFY = 3
!
!      ******************************************************************
!      *                                                                *
!      * Some error codes                                               *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: ALL_OK         = 0
       integer, parameter :: ERROR          = 1
       integer, parameter :: NODE_NOT_FOUND = 2
       integer, parameter :: INCORRECT_PATH = 3

       integer, parameter :: CG_OK             = 0
       integer, parameter :: CG_ERROR          = 1
       integer, parameter :: CG_NODE_NOT_FOUND = 2
       integer, parameter :: CG_INCORRECT_PATH = 3
!
!      ******************************************************************
!      *                                                                *
!      * Dimensional units                                              *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Null = 0
       integer, parameter :: UserDefined = 1

       integer, parameter :: Kilogram  = 2
       integer, parameter :: Gram      = 3
       integer, parameter :: Slug      = 4
       integer, parameter :: PoundMass = 5

       character(len=32), dimension(0:5), PARAMETER :: &
            MassUnitsName = (/"Null       ",           &
                              "UserDefined",           &
                              "Kilogram   ",           &
                              "Gram       ",           &
                              "Slug       ",           &
                              "PoundMass  "/)

       integer, parameter :: Meter      = 2
       integer, parameter :: Centimeter = 3
       integer, parameter :: Millimeter = 4
       integer, parameter :: Foot       = 5
       integer, parameter :: Inch       = 6

       character(len=32), dimension(0:6), PARAMETER :: &
            LengthUnitsName = (/"Null       ",         &
                                "UserDefined",         &
                                "Meter      ",         &
                                "Centimeter ",         &
                                "Millimeter ",         &
                                "Foot       ",         &
                                "Inch       "/)

       integer, parameter :: Second = 2

       character(len=32), dimension(0:2), PARAMETER :: &
            TimeUnitsName = (/"Null       ",           &
                              "UserDefined",           &
                              "Second     "/)

       integer, parameter :: Kelvin     = 2
       integer, parameter :: Celcius    = 3
       integer, parameter :: Rankine    = 4
       integer, parameter :: Fahrenheit = 5

       character(len=32), dimension(0:5), PARAMETER :: &
            TemperatureUnitsName = (/"Null       ",    &
                                     "UserDefined",    &
                                     "Kelvin     ",    &
                                     "Celcius    ",    &
                                     "Rankine    ",    &
                                     "Fahrenheit "/)

       integer, parameter :: Degree = 2
       integer, parameter :: Radian = 3

       character(len=32), dimension(0:3), PARAMETER :: &
            AngleUnitsName = (/"Null       ",          &
                               "UserDefined",          &
                               "Degree     ",          &
                               "Radian     "/)

       integer, parameter :: Ampere     = 2
       integer, parameter :: Abampere   = 3
       integer, parameter :: Statampere = 4
       integer, parameter :: Edison     = 5
       integer, parameter :: auCurrent  = 6

       character(len=32), dimension(0:6), PARAMETER ::  &
            ElectricCurrentUnitsName = (/"Null       ", &
                                         "UserDefined", &
                                         "Ampere     ", &
                                         "Abampere   ", &
                                         "Statampere ", &
                                         "Edison     ", &
                                         "a.u.       "/)

       integer, parameter :: Mole               = 2
       integer, parameter :: Entities           = 3
       integer, parameter :: StandardCubicFoot  = 4
       integer, parameter :: StandardCubicMeter = 5

       character(len=32), dimension(0:5), PARAMETER ::         &
            SubstanceAmountUnitsName = (/"Null              ", &
                                         "UserDefined       ", &
                                         "Mole              ", &
                                         "Entities          ", &
                                         "StandardCubicFoot ", &
                                         "StandardCubicMeter"/)

       integer, parameter :: Candela = 2
       integer, parameter :: Candle  = 3
       integer, parameter :: Carcel  = 4
       integer, parameter :: Hefner  = 5
       integer, parameter :: Violle  = 6

       character(len=32), dimension(0:6), PARAMETER ::    &
            LuminousIntensityUnitsName = (/"Null       ", &
                                           "UserDefined", &
                                           "Candela    ", &
                                           "Candle     ", &
                                           "Carcel     ", &
                                           "Hefner     ", &
                                           "Violle     "/)
!
!      ******************************************************************
!      *                                                                *
!      * Data class                                                     *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Dimensional                    = 2
       integer, parameter :: NormalizedByDimensional        = 3
       integer, parameter :: NormalizedByUnknownDimensional = 4
       integer, parameter :: NondimensionalParameter        = 5
       integer, parameter :: DimensionlessConstant          = 6

       character(len=32), dimension(0:6), PARAMETER ::          &
            DataClassName = (/"Null                          ", &
                              "UserDefined                   ", &
                              "Dimensional                   ", &
                              "NormalizedByDimensional       ", &
                              "NormalizedByUnknownDimensional", &
                              "NondimensionalParameter       ", &
                              "DimensionlessConstant         "/)
!
!      ******************************************************************
!      *                                                                *
!      * Grid location                                                  *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Vertex      = 2
       integer, parameter :: CellCenter  = 3
       integer, parameter :: FaceCenter  = 4
       integer, parameter :: IFaceCenter = 5
       integer, parameter :: JFaceCenter = 6
       integer, parameter :: KFaceCenter = 7
       integer, parameter :: EdgeCenter  = 8

       character(len=32), dimension(0:8), PARAMETER :: &
            GridLocationName = (/"Null       ",        &
                                 "UserDefined",        &
                                 "Vertex     ",        &
                                 "CellCenter ",        &
                                 "FaceCenter ",        &
                                 "IFaceCenter",        &
                                 "JFaceCenter",        &
                                 "KFaceCenter",        &
                                 "EdgeCenter "/)
!
!      ******************************************************************
!      *                                                                *
!      * Grid connectivity types                                        *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Overset      = 2
       integer, parameter :: Abutting     = 3
       integer, parameter :: Abutting1to1 = 4

       character(len=32), dimension(0:4), PARAMETER ::   &
            GridConnectivityTypeName = (/"Null        ", &
                                         "UserDefined ", &
                                         "Overset     ", &
                                         "Abutting    ", &
                                         "Abutting1to1"/)
!
!      ******************************************************************
!      *                                                                *
!      * Point set types                                                *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: PointList       = 2
       integer, parameter :: PointListDonor  = 3
       integer, parameter :: PointRange      = 4
       integer, parameter :: PointRangeDonor = 5
       integer, parameter :: ElementRange    = 6
       integer, parameter :: ElementList     = 7
       integer, parameter :: CellListDonor   = 8

       character(len=32), dimension(0:8), PARAMETER :: &
            PointSetTypeName = (/"Null           ",    &
                                 "UserDefined    ",    &
                                 "PointList      ",    &
                                 "PointListDonor ",    &
                                 "PointRange     ",    &
                                 "PointRangeDonor",    &
                                 "ElementRange   ",    &
                                 "ElementList    ",    &
                                 "CellListDonor  "/)
!
!      ******************************************************************
!      *                                                                *
!      * Governing equations and physical models types                  *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: FullPotential             = 2
       integer, parameter :: Euler                     = 3
       integer, parameter :: NSLaminar                 = 4
       integer, parameter :: NSTurbulent               = 5
       integer, parameter :: NSLaminarIncompressible   = 6
       integer, parameter :: NSTurbulentIncompressible = 7

       character(len=32), dimension(0:7), PARAMETER ::                  &
            GoverningEquationsTypeName = (/"Null                     ", &
                                           "UserDefined              ", &
                                           "FullPotential            ", &
                                           "Euler                    ", &
                                           "NSLaminar                ", &
                                           "NSTurbulent              ", &
                                           "NSLaminarIncompressible  ", &
                                           "NSTurbulentIncompressible"/)

       integer, parameter :: Ideal                       =  2
       integer, parameter :: VanderWaals                 =  3
       integer, parameter :: Constant                    =  4
       integer, parameter :: PowerLaw                    =  5
       integer, parameter :: SutherlandLaw               =  6
       integer, parameter :: ConstantPrandtl             =  7
       integer, parameter :: EddyViscosity               =  8
       integer, parameter :: ReynoldsStress              =  9
       integer, parameter :: ReynoldsStressAlgebraic     = 10
       integer, parameter :: Algebraic_BaldwinLomax      = 11
       integer, parameter :: Algebraic_CebeciSmith       = 12
       integer, parameter :: HalfEquation_JohnsonKing    = 13
       integer, parameter :: OneEquation_BaldwinBarth    = 14
       integer, parameter :: OneEquation_SpalartAllmaras = 15
       integer, parameter :: TwoEquation_JonesLaunder    = 16
       integer, parameter :: TwoEquation_MenterSST       = 17
       integer, parameter :: TwoEquation_Wilcox          = 18
       integer, parameter :: CaloricallyPerfect          = 19
       integer, parameter :: ThermallyPerfect            = 20
       integer, parameter :: ConstantDensity             = 21
       integer, parameter :: RedlichKwong                = 22
       integer, parameter :: Frozen                      = 23
       integer, parameter :: ThermalEquilib              = 24
       integer, parameter :: ThermalNonequilib           = 25
       integer, parameter :: ChemicalEquilibCurveFit     = 26
       integer, parameter :: ChemicalEquilibMinimization = 27
       integer, parameter :: ChemicalNonequilib          = 28
       integer, parameter :: EMElectricField             = 29
       integer, parameter :: EMMagneticField             = 30
       integer, parameter :: EMConductivity              = 31
       integer, parameter :: Voltage                     = 32
       integer, parameter :: Interpolated                = 33
       integer, parameter :: Equilibrium_LinRessler      = 34
       integer, parameter :: Chemistry_LinRessler        = 35

       character(len=32), dimension(0:35), PARAMETER ::      &
            ModelTypeName = (/"Null                       ", &
                              "UserDefined                ", &
                              "Ideal                      ", &
                              "VanderWaals                ", &
                              "Constant                   ", &
                              "PowerLaw                   ", &
                              "SutherlandLaw              ", &
                              "ConstantPrandtl            ", &
                              "EddyViscosity              ", &
                              "ReynoldsStress             ", &
                              "ReynoldsStressAlgebraic    ", &
                              "Algebraic_BaldwinLomax     ", &
                              "Algebraic_CebeciSmith      ", &
                              "HalfEquation_JohnsonKing   ", &
                              "OneEquation_BaldwinBarth   ", &
                              "OneEquation_SpalartAllmaras", &
                              "TwoEquation_JonesLaunder   ", &
                              "TwoEquation_MenterSST      ", &
                              "TwoEquation_Wilcox         ", &
                              "CaloricallyPerfect         ", &
                              "ThermallyPerfect           ", &
                              "ConstantDensity            ", &
                              "RedlichKwong               ", &
                              "Frozen                     ", &
                              "ThermalEquilib             ", &
                              "ThermalNonequilib          ", &
                              "ChemicalEquilibCurveFit    ", &
                              "ChemicalEquilibMinimization", &
                              "ChemicalNonequilib         ", &
                              "EMElectricField            ", &
                              "EMMagneticField            ", &
                              "EMConductivity             ", &
                              "Voltage                    ", &
                              "Interpolated               ", &
                              "Equilibrium_LinRessler     ", &
                              "Chemistry_LinRessler       "/)
!
!      ******************************************************************
!      *                                                                *
!      * Boundary condition types                                       *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: BCAxisymmetricWedge     =  2
       integer, parameter :: BCDegenerateLine        =  3
       integer, parameter :: BCDegeneratePoint       =  4
       integer, parameter :: BCDirichlet             =  5
       integer, parameter :: BCExtrapolate           =  6
       integer, parameter :: BCFarfield              =  7
       integer, parameter :: BCGeneral               =  8
       integer, parameter :: BCInflow                =  9
       integer, parameter :: BCInflowSubsonic        = 10
       integer, parameter :: BCInflowSupersonic      = 11
       integer, parameter :: BCNeumann               = 12
       integer, parameter :: BCOutflow               = 13
       integer, parameter :: BCOutflowSubsonic       = 14
       integer, parameter :: BCOutflowSupersonic     = 15
       integer, parameter :: BCSymmetryPlane         = 16
       integer, parameter :: BCSymmetryPolar         = 17
       integer, parameter :: BCTunnelInflow          = 18
       integer, parameter :: BCTunnelOutflow         = 19
       integer, parameter :: BCWall                  = 20
       integer, parameter :: BCWallInviscid          = 21
       integer, parameter :: BCWallViscous           = 22
       integer, parameter :: BCWallViscousHeatFlux   = 23
       integer, parameter :: BCWallViscousIsothermal = 24
       integer, parameter :: FamilySpecified         = 25

       character(len=32), dimension(0:25), PARAMETER :: &
            BCTypeName = (/"Null                   ",   &
                           "UserDefined            ",   &
                           "BCAxisymmetricWedge    ",   &
                           "BCDegenerateLine       ",   &
                           "BCDegeneratePoint      ",   &
                           "BCDirichlet            ",   &
                           "BCExtrapolate          ",   &
                           "BCFarfield             ",   &
                           "BCGeneral              ",   &
                           "BCInflow               ",   &
                           "BCInflowSubsonic       ",   &
                           "BCInflowSupersonic     ",   &
                           "BCNeumann              ",   &
                           "BCOutflow              ",   &
                           "BCOutflowSubsonic      ",   &
                           "BCOutflowSupersonic    ",   &
                           "BCSymmetryPlane        ",   &
                           "BCSymmetryPolar        ",   &
                           "BCTunnelInflow         ",   &
                           "BCTunnelOutflow        ",   &
                           "BCWall                 ",   &
                           "BCWallInviscid         ",   &
                           "BCWallViscous          ",   &
                           "BCWallViscousHeatFlux  ",   &
                           "BCWallViscousIsothermal",   &
                           "FamilySpecified        "/)
!
!      ******************************************************************
!      *                                                                *
!      * Data types                                                     *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Integer    = 2
       integer, parameter :: RealSingle = 3
       integer, parameter :: RealDouble = 4
       integer, parameter :: Character  = 5

       character(len=32), dimension(0:5), PARAMETER :: &
            DataTypeName = (/"Null       ",            &
                             "UserDefined",            &
                             "Integer    ",            &
                             "RealSingle ",            &
                             "RealDouble ",            &
                             "Character  "/)
!
!      ******************************************************************
!      *                                                                *
!      * BCData_t types                                                 *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Dirichlet = 2
       integer, parameter :: Neumann   = 3

       character(len=32), dimension(0:3), PARAMETER :: &
            BCDataTypeName = (/"Null       ",          &
                               "UserDefined",          &
                               "Dirichlet  ",          &
                               "Neumann    "/)
!
!      ******************************************************************
!      *                                                                *
!      * Element types                                                  *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: NODE     =  2
       integer, parameter :: BAR_2    =  3
       integer, parameter :: BAR_3    =  4
       integer, parameter :: TRI_3    =  5
       integer, parameter :: TRI_6    =  6
       integer, parameter :: QUAD_4   =  7
       integer, parameter :: QUAD_8   =  8
       integer, parameter :: QUAD_9   =  9
       integer, parameter :: TETRA_4  = 10
       integer, parameter :: TETRA_10 = 11
       integer, parameter :: PYRA_5   = 12
       integer, parameter :: PYRA_14  = 13
       integer, parameter :: PENTA_6  = 14
       integer, parameter :: PENTA_15 = 15
       integer, parameter :: PENTA_18 = 16
       integer, parameter :: HEXA_8   = 17
       integer, parameter :: HEXA_20  = 18
       integer, parameter :: HEXA_27  = 19
       integer, parameter :: MIXED    = 20
       integer, parameter :: NGON_n   = 21

       character(len=32), dimension(0:21), PARAMETER :: &
            ElementTypeName = (/"Null       ",          &
                                "UserDefined",          &
                                "NODE       ",          &
                                "BAR_2      ",          &
                                "BAR_3      ",          &
                                "TRI_3      ",          &
                                "TRI_6      ",          &
                                "QUAD_4     ",          &
                                "QUAD_8     ",          &
                                "QUAD_9     ",          &
                                "TETRA_4    ",          &
                                "TETRA_10   ",          &
                                "PYRA_5     ",          &
                                "PYRA_14    ",          &
                                "PENTA_6    ",          &
                                "PENTA_15   ",          &
                                "PENTA_18   ",          &
                                "HEXA_8     ",          &
                                "HEXA_20    ",          &
                                "HEXA_27    ",          &
                                "MIXED      ",          &
                                "NGON_n     "/)
!
!      ******************************************************************
!      *                                                                *
!      * Zone types                                                     *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Structured   =  2
       integer, parameter :: Unstructured =  3

       character(len=32), dimension(0:3), PARAMETER :: &
            ZoneTypeName = (/"Null        ",           &
                             "UserDefined ",           &
                             "Structured  ",           &
                             "Unstructured"/)
!
!      ******************************************************************
!      *                                                                *
!      * Rigid grid motion types                                        *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: ConstantRate = 2
       integer, parameter :: VariableRate = 3

       character(len=32), dimension(0:3), PARAMETER ::  &
            RigidGridMotionTypeName = (/"Null        ", &
                                        "UserDefined ", &
                                        "ConstantRate", &
                                        "VariableRate"/)
!
!      ******************************************************************
!      *                                                                *
!      * Arbitrary grid motion types                                    *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: NonDeformingGrid = 2
       integer, parameter :: DeformingGrid    = 3

       character(len=32), dimension(0:3), PARAMETER ::          &
            ArbitraryGridMotionTypeName = (/"Null            ", &
                                            "UserDefined     ", &
                                            "NonDeformingGrid", &
                                            "DeformingGrid   "/)
!
!      ******************************************************************
!      *                                                                *
!      * Simulation type                                                *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: TimeAccurate    = 2
       integer, parameter :: NonTimeAccurate = 3

       character(len=32), dimension(0:3), PARAMETER :: &
            SimulationTypeName = (/"Null           ",  &
                                   "UserDefined    ",  &
                                   "TimeAccurate   ",  &
                                   "NonTimeAccurate"/)
!
!      ******************************************************************
!      *                                                                *
!      * BC Property types                                              *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: Generic = 2

       character(len=32), dimension(0:2), PARAMETER :: &
            WallFunctionTypeName = (/"Null       ",    &
                                     "UserDefined",    &
                                     "Generic    "/)

       integer, parameter :: BleedArea   = 2
       integer, parameter :: CaptureArea = 3

       character(len=32), dimension(0:3), PARAMETER :: &
            AreaTypeName = (/"Null       ",            &
                             "UserDefined",            &
                             "BleedArea  ",            &
                             "CaptureArea"/)
!
!      ******************************************************************
!      *                                                                *
!      * Grid connectivity property types                               *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: AverageAll             = 2
       integer, parameter :: AverageCircumferential = 3
       integer, parameter :: AverageRadial          = 4
       integer, parameter :: AverageI               = 5
       integer, parameter :: AverageJ               = 6
       integer, parameter :: AverageK               = 7

       character(len=32), dimension(0:7), PARAMETER ::            &
           AverageInterfaceTypeName = (/"Null                  ", &
                                        "UserDefined           ", &
                                        "AverageAll            ", &
                                        "AverageCircumferential", &
                                        "AverageRadial         ", &
                                        "AverageI              ", &
                                        "AverageJ              ", &
                                        "AverageK              "/)

#else
       ! CGNS is used. simply include the file cgnslib_f.h.

       implicit none
       save
       include "cgnslib_f.h"

#endif

       end module su_cgns

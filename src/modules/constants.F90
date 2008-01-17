!
!      ******************************************************************
!      *                                                                *
!      * File:          constants.F90                                   *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      * Starting date: 12-09-2002                                      *
!      * Last modified: 07-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module constants
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the constants used in the code.                  *
!      *                                                                *
!      ******************************************************************
!
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

       ! Constant in Sutherland's law; SI-units.

       real(kind=realType), parameter :: SSuthDim  = 110.55_realType
       real(kind=realType), parameter :: muSuthDim = 1.716e-5_realType
       real(kind=realType), parameter :: TSuthDim  = 273.15_realType

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

       ! Threshold parameter for real types; the value depends
       ! whether single or double precision is used.

#ifdef USE_SINGLE_PRECISION
       real(kind=realType), parameter :: thresholdReal = 1.e-5_realType
#else
       real(kind=realType), parameter :: thresholdReal = 1.e-10_realType
#endif

       ! Definition of the tab and carriage return character.

       character(len=1), parameter :: tabChar = achar(9)
       character(len=1), parameter :: retChar = achar(13)

       end module constants

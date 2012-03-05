!
!      ******************************************************************
!      *                                                                *
!      * File:          BCTypes.f90                                     *
!      * Author:        Edwin van der Weide, Juan J. Alonso,            *
!      *                Steve Repsher, Seonghyeon Hahn                  *
!      * Starting date: 01-21-2003                                      *
!      * Last modified: 11-04-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module BCTypes
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the definition of the all boundary        *
!      * condition types available including internal block boundaries. *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Supported boundary condition types. These symbolic names should
!      be used instead of their actual numbers.
!
       integer(kind=intType), parameter :: BCNull                =   0_intType
       integer(kind=intType), parameter :: Symm                  =  -1_intType
       integer(kind=intType), parameter :: SymmPolar             =  -2_intType
       integer(kind=intType), parameter :: NSWallAdiabatic       =  -3_intType
       integer(kind=intType), parameter :: NSWallIsothermal      =  -4_intType
       integer(kind=intType), parameter :: EulerWall             =  -5_intType
       integer(kind=intType), parameter :: FarField              =  -6_intType
       integer(kind=intType), parameter :: SupersonicInflow      =  -7_intType
       integer(kind=intType), parameter :: SubsonicInflow        =  -8_intType
       integer(kind=intType), parameter :: SupersonicOutflow     =  -9_intType
       integer(kind=intType), parameter :: SubsonicOutflow       = -10_intType
       integer(kind=intType), parameter :: MassBleedInflow       = -11_intType
       integer(kind=intType), parameter :: MassBleedOutflow      = -12_intType
       integer(kind=intType), parameter :: mDot                  = -13_intType
       integer(kind=intType), parameter :: Thrust                = -14_intType
       integer(kind=intType), parameter :: Extrap                = -15_intType
       integer(kind=intType), parameter :: B2BMatch              = -16_intType
       integer(kind=intType), parameter :: B2BMismatch           = -17_intType
       integer(kind=intType), parameter :: SlidingInterface      = -18_intType
       integer(kind=intType), parameter :: OversetOuterBound     = -19_intType
       integer(kind=intType), parameter :: DomainInterfaceAll    = -20_intType
       integer(kind=intType), parameter :: DomainInterfaceRhoUVW = -21_intType
       integer(kind=intType), parameter :: DomainInterfaceP      = -22_intType
       integer(kind=intType), parameter :: DomainInterfaceRho    = -23_intType
       integer(kind=intType), parameter :: DomainInterfaceTotal  = -24_intType
       integer(kind=intType), parameter :: BCNotValid            = -25_intType
!
!      Number of actual boundary conditions supported by the code
!      This number refers to bocos, not flow-through BCs
!      Edit this number when additional boundary conditions are
!      supported
!
       integer(kind=intType), parameter :: nBCs = 24_intType
!
!      block faces on which boundary conditions may be imposed
!
       integer(kind=intType), parameter :: iMin = 1_intType
       integer(kind=intType), parameter :: iMax = 2_intType
       integer(kind=intType), parameter :: jMin = 3_intType
       integer(kind=intType), parameter :: jMax = 4_intType
       integer(kind=intType), parameter :: kMin = 5_intType
       integer(kind=intType), parameter :: kMax = 6_intType
!
       end module BCTypes

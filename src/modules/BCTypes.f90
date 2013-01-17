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
!
!      block faces on which boundary conditions may be imposed
!
       integer(kind=intType), parameter :: iMin = 1
       integer(kind=intType), parameter :: iMax = 2
       integer(kind=intType), parameter :: jMin = 3
       integer(kind=intType), parameter :: jMax = 4
       integer(kind=intType), parameter :: kMin = 5
       integer(kind=intType), parameter :: kMax = 6
!
       end module BCTypes

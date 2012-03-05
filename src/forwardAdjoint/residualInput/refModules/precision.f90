!
!      ******************************************************************
!      *                                                                *
!      * File:          precision.F90                                   *
!      * Author:        Gaetan Kenway                                   *
!      *                                                                *
!      ******************************************************************
!
       module precision
!
!      ******************************************************************
!      *                                                                *
!      * A modified precision module for AD rouitines                   *
!      *                                                                *
!      ******************************************************************

       implicit none
       save

       ! ADjoint routines are only known to work with 4 byte integers:

       integer(kind=4), private :: dummyInt
       integer, parameter       :: sizeOfInteger = 4

       ! ADjoint routines are only known to work with 8 byte reals:
       
       real(kind=8), private :: dummyReal
       integer, parameter    :: sizeOfReal = 8

       real(kind=8), private :: dummyCGNSReal

!
!      ******************************************************************
!      *                                                                *
!      * Definition of the porosity type. As this is only a flag to     *
!      * indicate whether or not fluxes must be computed, an integer1   *
!      * is perfectly okay.                                             *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=1), private :: dummyPor

!      ******************************************************************
!      *                                                                *
!      * Definition of the kind parameters for the integer and real     *
!      * types.                                                         *
!      *                                                                *
!      ******************************************************************
!
       integer, parameter :: intType      = kind(dummyInt)
       integer, parameter :: porType      = kind(dummyPor)
       integer, parameter :: realType     = kind(dummyReal)

       end module precision

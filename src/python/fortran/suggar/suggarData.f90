!
!      ******************************************************************
!      *                                                                *
!      * File:          suggarData.f90                                  *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 09-02-2005                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module suggarData
!
!      ******************************************************************
!      *                                                                *
!      * Variables used in the routines for interfacing with SUGGAR.    *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Type definition to store the way the original CGNS zones are   *
!      * split among the processes.                                     *
!      *                                                                *
!      ******************************************************************
!
       type splitZoneType

         ! nSubBlocks:             # of subblocks into which the original
         !                         zone is split.
         ! proc(nSubBlocks):       Process number for each subblock.
         ! localBlock(nSubBlocks): Local block number for each subblock.
         ! iBegOr(nSubBlocks):     The nodal range in the original
         ! etc.                    block of these subblocks.

         integer(kind=intType) :: nSubBlocks

         integer(kind=intType), dimension(:), pointer :: proc
         integer(kind=intType), dimension(:), pointer :: localBlock

         integer(kind=intType), dimension(:), pointer :: iBegOr, iEndOr
         integer(kind=intType), dimension(:), pointer :: jBegOr, jEndOr
         integer(kind=intType), dimension(:), pointer :: kBegOr, kEndOr

       end type splitZoneType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the variables.                                   *
!      *                                                                *
!      ******************************************************************
!
       ! nZones:                Same as cgnsNDom.
       ! zoneNames(nZones):     Sorted list of all zone names.
       ! unsortedZone(nZones):  Corresponding zone number for each of
       !                        the zones in the sorted list.
       ! splitData(nZones):     Split type defining the way the zones
       !                        were split among the processes.

       integer(kind=intType) :: nZones
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                              zoneNames
       integer(kind=intType), dimension(:), allocatable :: unsortedZone
       type(splitZoneType),  dimension(:), allocatable :: splitData

       end module suggarData

!
!      ******************************************************************
!      *                                                                *
!      * File:          coarse1to1Subface.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-05-2003                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module coarse1to1Subface
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived datatype                *
!      * coarse1to1SubfaceType, which is used to determine the 1 to 1   *
!      * block boundaries for the coarser grids.                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype.                        *
!      *                                                                *
!      ******************************************************************
!
       type coarse1to1SubfaceType

         ! Nodal range in the three coordinates directions for the
         ! coarse grid subface.

         integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd

         ! Processor and block id of the neighboring block.

         integer(kind=intType) :: neighProc, neighBlock

         ! Number of points in the three coordinate directions for the
         ! coarse grid donor subface.

         integer(kind=intType) :: ndi, ndj, ndk

         ! Corresponding i, j and k indices of the fine grid donor block
         ! for each of the coarse grid subface lines.

         integer(kind=intType), dimension(:), pointer :: idfine
         integer(kind=intType), dimension(:), pointer :: jdfine
         integer(kind=intType), dimension(:), pointer :: kdfine

       end type coarse1to1SubfaceType

       ! Number of 1 to 1 fine grid subfaces on this processor.

       integer(kind=intType) :: nSubface1to1

       ! Array of 1 to 1 subfaces.

       type(coarse1to1SubfaceType), dimension(:), allocatable :: subface1to1

       end module coarse1to1Subface

!      ==================================================================

       module coarseningInfo
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived datatype                *
!      * coarseningInfoType, which stores for a given block the grid    *
!      * lines to keep for the coarse grid.                             *
!      *                                                                *
!      ******************************************************************
!
       type coarseningInfoType

         ! Logical, which indicate whether or not a fine grid 1 to 1
         ! block boundary is still a 1 to 1 block boundary on the
         ! coarse grid.

         logical, dimension(:), pointer :: coarseIs1to1

       end type coarseningInfoType

       ! Array to store the info for all the blocks.

       type(coarseningInfoType), dimension(:), allocatable :: coarseInfo

       end module coarseningInfo

!
!      ******************************************************************
!      *                                                                *
!      * File:          boundaryList.f90                                *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-06-2005                                      *
!      * Last modified: 08-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module boundaryList
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains temporary variables to create the   *
!      * list of overset boundary cells and their donor info. It also   *
!      * contains some variables that assist in the coarse boundary     *
!      * creation.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

!      ******************************************************************
!      *                                                                *
!      * The definition of the variables for the lists. Note this is a  *
!      * repeat of the derived type declared in the preprocessing local *
!      * halo list module, and is here so that some routines from that  *
!      * directory may be used with no confusion.                       *
!      *                                                                *
!      ******************************************************************
!
       type haloListType

         ! myBlock        : local block id of the halo.
         ! myI, myJ, myK  : i,j,k indices of the halo.
         ! donorProc      : processor where donor is stored. In case
         !                  the halo is a boundary halo, donorProc
         !                  is set to -1.
         ! donorBlock     : block id of the donor. In case the halo
         !                  is a boundary halo donorBlock is set
         !                  to the corresponding boundary condition.
         ! dI, dJ, dK     : i,j,k indices of the donor.
         ! levOfInd       : level of indirectness.
         ! interp(..)     : interpolants for the donor stencil; only
         !                  allocated for lists requiring this info.

         integer(kind=intType) :: myBlock
         integer(kind=intType) :: myI, myJ, myK
         integer(kind=intType) :: donorProc, donorBlock
         integer(kind=intType) :: dI, dJ, dK
         integer(kind=intType) :: levOfInd

         real(kind=realType), dimension(:), pointer :: interp

       end type haloListType

       integer(kind=intType) :: nHaloOver
       type(haloListType), dimension(:), allocatable :: oversetHalo
!
!      ******************************************************************
!      *                                                                *
!      * Definition and variable to help in estimating the donor info   *
!      * for a coarse grid boundary from the the next finier grid.      *
!      *                                                                *
!      ******************************************************************
!
       type blockBndryType
 
         ! istartOverHole       - starting index into boundary data
         !                        which the boundary cell contains
         !                        no finer grid boundary. This is an
         !                        index into ibndry in flowDoms.
         ! nearestBndry(        - index into bndry info for those
         !   istartOverHole:)     cells added over top of fine holes.
         !                        Also an index into ibndry but this
         !                        index is < istartOverHole.
 
         integer(kind=intType) :: istartOverHole
         integer(kind=intType), pointer :: nearestBndry(:)

       end type blockBndryType

       type(blockBndryType), dimension(:), allocatable :: blockBndry
!
!      ******************************************************************
!      *                                                                *
!      * Other variables in this module.                                *
!      *                                                                *
!      ******************************************************************
!
       ! fringeSize:      The required fringe size in each indical
       !                  direction in order to compute the residual for
       !                  a cell (e.g. pure central differencing had a 
       !                  fringeSize = 1).
       ! donorInfo(5,..): Array to store donor estimates (processor,
       !                  block, and 3 indices) for 1-to-1 halo cells
       !                  that are flagged for telling the connecting
       !                  block to add to its fringe.
       ! nDonorInfo:      Current number of donor info sets stored in
       !                  the donorInfo array.

       integer(kind=intType) :: fringeSize, nDonorInfo

       integer(kind=intType), dimension(:,:), pointer :: donorInfo

       end module boundaryList

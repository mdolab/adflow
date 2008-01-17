!
!      ******************************************************************
!      *                                                                *
!      * File:          indirectHalo.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-31-2003                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module indirectHalo
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived data type used to       *
!      * determine the indirect halo's as well as an array of this type.*
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualIndirectHaloType
       private :: lessIndirectHaloType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type indirectHaloType.      *
!      *                                                                *
!      ******************************************************************
!
       type indirectHaloType

         ! myBlock      : Local block ID of the halo.
         ! myI, myJ, myK: i,j,k indices of the halo.
         ! myDirectHalo : Index in the haloListType where the
         !                corresponding direct halo is stored.
         ! levOfInd     : Level of indirectness.
         ! donorProc    : Processor where donor of the direct halo is
         !                stored. In case this halo is a boundary
         !                halo, donorProc is set to -1.

         integer(kind=intType) :: myBlock
         integer(kind=intType) :: myI, myJ, myK
         integer(kind=intType) :: myDirectHalo
         integer(kind=intType) :: levOfInd
         integer(kind=intType) :: donorProc

       end type indirectHaloType

       ! Interface for the extension of the operators <= and <.
       ! These are needed for the sorting of indirectHaloType.

       interface operator(<=)
         module procedure lessEqualIndirectHaloType
       end interface

       interface operator(<)
         module procedure lessIndirectHaloType
       end interface

       ! nIndHalo         : Number of indirect halo's to be treated.
       ! indHalo(nIndHalo): The indirect halo's.

       integer(kind=intType) :: nIndHalo
       type(indirectHaloType), dimension(:), allocatable :: indHalo

       ! nLevOfInd               : Number of levels of indirectness
       ! nHaloPerLev(0:nLevOfInd): Number of indirect halo's per level
       !                           of indirectness; stored in
       !                           cumulative storage format.
       ! nHaloPerProc(0:nProc)   : Number of indirect halo's per
       !                           processor for a given level of
       !                           indirectness; cumulative storage.
       !                           nHaloPerProc(0) is not 0,
       !                           because of the presence of boundary
       !                           halo's, which get proc ID -1.

       integer(kind=intType) :: nLevOfInd
       integer(kind=intType), dimension(:), allocatable :: nHaloPerLev
       integer(kind=intType), dimension(:), allocatable :: nHaloPerProc

       contains
!
!        ****************************************************************
!        *                                                              *
!        * Functions to simulate the operators <= and <.                *
!        *                                                              *
!        ****************************************************************
!
         logical function lessEqualIndirectHaloType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 <= g2 and .false.         *
!        * otherwise. The comparison is firstly based on the level of   *
!        * indirectness followed by the donor processor, the            *
!        * corresponding direct halo, my block ID and finally the i, j  *
!        * and k indices.                                               *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(indirectHaloType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the level of indirectness. If not equal, set
         ! lessEqual appropriately and return.

         if(g1%levOfInd < g2%levOfInd) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%levOfInd > g2%levOfInd) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! Compare the donor processors.

         if(g1%donorProc < g2%donorProc) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%donorProc > g2%donorProc) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! Compare the direct halo.

         if(g1%myDirectHalo < g2%myDirectHalo) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%myDirectHalo > g2%myDirectHalo) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! Compare my block ID.

         if(g1%myBlock < g2%myBlock) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%myBlock > g2%myBlock) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! Finally compare the halo indices. Start with k.

         if(g1%myK < g2%myK) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%myK > g2%myK) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! The j index.

         if(g1%myJ < g2%myJ) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%myJ > g2%myJ) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! The i index.

         if(g1%myI < g2%myI) then
           lessEqualIndirectHaloType = .true.
           return
         else if(g1%myI > g2%myI) then
           lessEqualIndirectHaloType = .false.
           return
         endif

         ! Both entities are identical. So set lessEqual to .true.

         lessEqualIndirectHaloType = .true.

         end function lessEqualIndirectHaloType

!        ================================================================

         logical function lessIndirectHaloType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. If g1 < g2 and .false.          *
!        * otherwise. It is basically the same as the lessEqual         *
!        * function, except that the equality is now considered as      *
!        * .false.                                                      *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(indirectHaloType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the level of indirectness. If not equal, set
         ! lessIndirectHaloType appropriately and return.

         if(g1%levOfInd < g2%levOfInd) then
           lessIndirectHaloType = .true.
           return
         else if(g1%levOfInd > g2%levOfInd) then
           lessIndirectHaloType = .false.
           return
         endif

         ! Compare the donor processors.

         if(g1%donorProc < g2%donorProc) then
           lessIndirectHaloType = .true.
           return
         else if(g1%donorProc > g2%donorProc) then
           lessIndirectHaloType = .false.
           return
         endif

         ! Compare the direct halo.

         if(g1%myDirectHalo < g2%myDirectHalo) then
           lessIndirectHaloType = .true.
           return
         else if(g1%myDirectHalo > g2%myDirectHalo) then
           lessIndirectHaloType = .false.
           return
         endif

         ! Compare my block id.

         if(g1%myBlock < g2%myBlock) then
           lessIndirectHaloType = .true.
           return
         else if(g1%myBlock > g2%myBlock) then
           lessIndirectHaloType = .false.
           return
         endif

         ! Finally compare the halo indices. Start with k.

         if(g1%myK < g2%myK) then
           lessIndirectHaloType = .true.
           return
         else if(g1%myK > g2%myK) then
           lessIndirectHaloType = .false.
           return
         endif

         ! The j index.

         if(g1%myJ < g2%myJ) then
           lessIndirectHaloType = .true.
           return
         else if(g1%myJ > g2%myJ) then
           lessIndirectHaloType = .false.
           return
         endif

         ! The i index.

         if(g1%myI < g2%myI) then
           lessIndirectHaloType = .true.
           return
         else if(g1%myI > g2%myI) then
           lessIndirectHaloType = .false.
           return
         endif

         ! Both entities are identical.
         ! So set lessIndirectHaloType to .false.

         lessIndirectHaloType = .false.

         end function lessIndirectHaloType

       end module indirectHalo

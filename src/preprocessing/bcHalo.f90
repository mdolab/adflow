!
!      ******************************************************************
!      *                                                                *
!      * File:          bcHalo.f90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-18-2003                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module bcHalo
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived datatype bcHaloType,    *
!      * which is used to determine the boundary condition for an       *
!      * indirect halo when the nearest direct halo's are all boundary  *
!      * halo's.                                                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualBCHaloType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype.                        *
!      *                                                                *
!      ******************************************************************
!
       type bcHaloType

         ! directHalo: Index in the haloListType where the
         !             corresponding direct halo is stored.
         ! BC:         Corresponding boundary condition.

         integer(kind=intType) :: directHalo, BC

       end type bcHaloType

       ! Interface for the extension of the operator <= needed for the
       ! sorting of bcHaloType. Note that the = operator does not
       ! need to be defined, because bcHaloType only contains
       ! primitive types.

       interface operator(<=)
         module procedure lessEqualBCHaloType
       end interface

       contains
!
         logical function lessEqualBCHaloType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * Function to simulate the operator <= for bcHaloType.         *
!        * It first compares the boundary condition. If equal the index *
!        * of the direct halo is compared, although this is not really  *
!        * important.                                                   *
!        * LessEqual returns .true. if g1 <= g2 and .false. otherwise.  *
!        *                                                              *
!        ****************************************************************
!
         use BCTypes
         implicit none
!
!        Function arguments.
!
         type(bcHaloType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! First compare the boundary conditions. Note that the sequence
         ! in BCTypes is such that the most important BC has the
         ! highest number.

         if(g1%BC < g2%BC) then
           lessEqualBCHaloType = .true.
           return
         else if(g1%BC > g2%BC) then
           lessEqualBCHaloType = .false.
           return
         endif

         ! Boundary conditions are equal. Just compare the index.

         if(g1%directHalo < g2%directHalo) then
           lessEqualBCHaloType = .true.
           return
         else if(g1%directHalo > g2%directHalo) then
           lessEqualBCHaloType = .false.
           return
         endif

         ! g1 and g2 are equal. Return .true.

         lessEqualBCHaloType = .true.

         end function lessEqualBCHaloType

!        ================================================================

         subroutine sortBCHaloType(bcHaloArray, nn)
!
!        ****************************************************************
!        *                                                              *
!        * SortBCHaloType sorts the given number of BCHalo's in         *
!        * increasing order. Note that this routine is called sort and  *
!        * not qsort, because only an insertion sort is done here. The  *
!        * reason is that nn <= 3 and thus an insertion sort is okay.   *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments
!
         integer(kind=intType), intent(in) :: nn
         type(bcHaloType), dimension(*), intent(inout) :: bcHaloArray
!
!        Local variables.
!
         integer(kind=intType) :: i, j

         type(bcHaloType) :: a
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         do j=1,nn
           a = bcHaloArray(j)
           do i=(j-1),1,-1
             if(bcHaloArray(i) <= a) exit
             bcHaloArray(i+1) = bcHaloArray(i)
           enddo
           bcHaloArray(i+1) = a
         enddo

         end subroutine sortBCHaloType

       end module bcHalo

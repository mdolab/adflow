!
!      ******************************************************************
!      *                                                                *
!      * File:          updateComm.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-14-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module updateComm
!
!      ******************************************************************
!      *                                                                *
!      * Local module to store the derived data types for the temporary *
!      * storage of the interpolation weights and for sorting them.     *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualSortedDonorComm
       private :: lessSortedDonorComm
       private :: equalSortedDonorComm

       ! Definition of the derived data type updateCommType

       type updateCommType

         ! ncopy1st:            Number of entities used in the
         !                      interpolation for 1st level halo cells.
         ! ncopy2nd:            Idem for the 2nd level halo cells.
         ! indBuf(ncopy2nd):    Corresponding index in the input buffer.
         !                      First the 1st halo cells will be stored,
         !                      followed by the second halo's.
         ! block(ncopy2nd):     Local block id of the donor cell used in
         !                      the interpolation.
         ! indices(ncopy2nd,3): Indices in the block of that donor cell.
         ! weight(ncopy2nd):    Interpolation weights.

         integer(kind=intType) :: ncopy1st, ncopy2nd

         integer(kind=intType), pointer, dimension(:)   :: indBuf
         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices
         real(kind=realType),   pointer, dimension(:)   :: weight

       end type updateCommType

       ! Definition of the derived data type sortedDonorCommType

       type sortedDonorCommType

         ! haloLevel:  Level of the corresponding halo cell,
         !             either 1 or 2.
         ! block:      Block ID of the donor cell.
         ! indices(3): Indices in the block of the donor cell.

         integer(kind=intType) ::               haloLevel, block
         integer(kind=intType), dimension(3) :: indices

       end type sortedDonorCommType

       ! Interface for the extension of the operators <=, < and ==.
       ! These are needed for the sorting and searching of
       ! sortedDonorCommType. Note that the = operator does not need
       ! to be defined, because sortedDonorCommType only contains
       ! primitive types.

       interface operator(<=)
         module procedure lessEqualSortedDonorComm
       end interface

       interface operator(<)
         module procedure lessSortedDonorComm
       end interface

       interface operator(==)
         module procedure equalSortedDonorComm
       end interface

       !=================================================================

       contains

         !===============================================================

         logical function lessEqualSortedDonorComm(g1,g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessEqualSortedDonorComm returns .true. If g1 <= g2 and      *
!        * .false. otherwise. The comparison is firstly based on the    *
!        * halo level, followed by the block ID and finally the indices.*
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(sortedDonorCommType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the halo level.

         if(g1%haloLevel < g2%haloLevel) then
           lessEqualSortedDonorComm = .true.
           return
         else if(g1%haloLevel > g2%haloLevel) then
           lessEqualSortedDonorComm = .false.
           return
         endif

         ! Halo levels are identical. Compare the block ID's.

         if(g1%block < g2%block) then
           lessEqualSortedDonorComm = .true.
           return
         else if(g1%block > g2%block) then
           lessEqualSortedDonorComm = .false.
           return
         endif

         ! Block ID's are also identical. Compare the indices.

         if(g1%indices(1) < g2%indices(1)) then
           lessEqualSortedDonorComm = .true.
           return
         else if(g1%indices(1) > g2%indices(1)) then
           lessEqualSortedDonorComm = .false.
           return
         endif

         if(g1%indices(2) < g2%indices(2)) then
           lessEqualSortedDonorComm = .true.
           return
         else if(g1%indices(2) > g2%indices(2)) then
           lessEqualSortedDonorComm = .false.
           return
         endif

         if(g1%indices(3) < g2%indices(3)) then
           lessEqualSortedDonorComm = .true.
           return
         else if(g1%indices(3) > g2%indices(3)) then
           lessEqualSortedDonorComm = .false.
           return
         endif

         ! g1 is identical to g2. Return .true.

         lessEqualSortedDonorComm = .true.

         end function lessEqualSortedDonorComm

         !===============================================================

         logical function lessSortedDonorComm(g1,g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessSortedDonorComm returns .true. If g1 < g2 and .false.    *
!        * otherwise. The comparison is firstly based on the halo       *
!        * level, followed by the block ID and finally the indices.     *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(sortedDonorCommType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the halo level.

         if(g1%haloLevel < g2%haloLevel) then
           lessSortedDonorComm = .true.
           return
         else if(g1%haloLevel > g2%haloLevel) then
           lessSortedDonorComm = .false.
           return
         endif

         ! Halo levels are identical. Compare the block ID's.

         if(g1%block < g2%block) then
           lessSortedDonorComm = .true.
           return
         else if(g1%block > g2%block) then
           lessSortedDonorComm = .false.
           return
         endif

         ! Block ID's are also identical. Compare the indices.

         if(g1%indices(1) < g2%indices(1)) then
           lessSortedDonorComm = .true.
           return
         else if(g1%indices(1) > g2%indices(1)) then
           lessSortedDonorComm = .false.
           return
         endif

         if(g1%indices(2) < g2%indices(2)) then
           lessSortedDonorComm = .true.
           return
         else if(g1%indices(2) > g2%indices(2)) then
           lessSortedDonorComm = .false.
           return
         endif

         if(g1%indices(3) < g2%indices(3)) then
           lessSortedDonorComm = .true.
           return
         else if(g1%indices(3) > g2%indices(3)) then
           lessSortedDonorComm = .false.
           return
         endif

         ! g1 is identical to g2. Return .false.

         lessSortedDonorComm = .false.

         end function lessSortedDonorComm

         !===============================================================

         logical function equalSortedDonorComm(g1,g2)
!
!        ****************************************************************
!        *                                                              *
!        * EqualSortedDonorComm returns .true. if g1 == g2 and          *
!        * .false. otherwise. All member variables must be identical    *
!        * for the equality to hold.                                    *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(sortedDonorCommType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize equalSortedDonorComm to .False.

         equalSortedDonorComm = .false.

         ! Check the member variables, which must all be identical.

         if(g1%haloLevel  /= g2%haloLevel)  return
         if(g1%block      /= g2%block)      return
         if(g1%indices(1) /= g2%indices(1)) return
         if(g1%indices(2) /= g2%indices(2)) return
         if(g1%indices(3) /= g2%indices(3)) return

         ! g1 is identical to g2. Return .true.

         equalSortedDonorComm = .true.

         end function equalSortedDonorComm

       end module updateComm

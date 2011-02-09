!
!      ******************************************************************
!      *                                                                *
!      * File:          periodicInfo.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-10-2003                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module periodicInfo
!
!      ******************************************************************
!      *                                                                *
!      * Local module that contains derived datatypes as well as arrays *
!      * of these derived datatypes to store information related to     *
!      * periodicity.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessCGNSPeriodicType
       private :: equalCGNSPeriodicType
       private :: lessPeriodicSubfacesHaloT
       private :: lessEqualPeriodicSubfacesHaloT
       private :: equalPeriodicSubfacesHaloT

       ! Definition of the derived data type for storing the periodic
       ! faces of the cgns grid a bit easier.

       type cgnsPeriodicType

         ! cgnsBlock   :: the block ID in the cgns grid.
         ! cgnsSubface :: the suface ID in this block.

         integer(kind=intType) :: cgnsBlock, cgnsSubface

       end type cgnsPeriodicType

       ! Interface for the extension of the operators < and ==.
       ! These are needed for the sorting and searching of
       ! cgnsPeriodicType.

       interface operator(<)
         module procedure lessCGNSPeriodicType
       end interface

       interface operator(==)
         module procedure equalCGNSPeriodicType
       end interface

       ! nPeriodicGlobal :: Total number of periodic faces in cgns grid.
       ! periodicGlobal  :: The corresponding faces.

       integer(kind=intType) :: nPeriodicGlobal
       type(cgnsPeriodicType), dimension(:), allocatable :: periodicGlobal

       ! Definition of the derived data type to store the periodic
       ! subfaces that are crossed when the halo and the donor are
       ! connected. The direction is from the halo to the donor.

       type periodicSubfacesHaloType

         ! internalHalo:      Whether or not the halo is an internal
         !                    halo, i.e. it is stored on the
         !                    same processor as the donor.
         ! indexInHaloList:   The corresponding index in either the
         !                    node or cell halo list.
         ! nPeriodicSubfaces: Number of periodic subfaces that are
         !                    crossed. This is at most the level of
         !                    indirectness of the halo.
         ! periodicSubfaces:  The corresponding subfaces ID's according
         !                    to the sequence defined in periodicGlobal.

         logical               :: internalHalo
         integer(kind=intType) :: indexInHaloList
         integer(kind=intType) :: nPeriodicSubfaces
         integer(kind=intType), dimension(:), pointer :: periodicSubfaces

       end type periodicSubfacesHaloType

       ! Interface for the extension of the operators <, <= and ==.
       ! This is needed for the sorting and comparing of variables of
       ! the type periodicSubfacesHaloType

       interface operator(<)
         module procedure lessPeriodicSubfacesHaloT
       end interface

       interface operator(<=)
         module procedure lessEqualPeriodicSubfacesHaloT
       end interface

       interface operator(==)
         module procedure equalPeriodicSubfacesHaloT
       end interface

       !=================================================================

       contains

       !=================================================================
!
!        ****************************************************************
!        *                                                              *
!        * Functions to simulate the operators < and ==.                *
!        *                                                              *
!        ****************************************************************
!
         logical function lessCGNSPeriodicType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessCGNSPeriodicType returns .true. if g1 is considered      *
!        * smaller than g2. This comparison is first based on the block *
!        * ID followed by the subface id.                               *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(cgnsPeriodicType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the block ID. If not equal set lessCGNSPeriodicType
         ! accordingly.

         if(g1%cgnsBlock < g2%cgnsBlock) then
           lessCGNSPeriodicType = .true.
           return
         else if(g1%cgnsBlock > g2%cgnsBlock) then
           lessCGNSPeriodicType = .false.
           return
         endif

         ! Block ID's are identical. Compare the subfaces.

         if(g1%cgnsSubface < g2%cgnsSubface) then
           lessCGNSPeriodicType = .true.
           return
         else if(g1%cgnsSubface > g2%cgnsSubface) then
           lessCGNSPeriodicType = .false.
           return
         endif

         ! Both objects are identical.
         ! Set lessCGNSPeriodicType to .false.

         lessCGNSPeriodicType = .false.

         end function lessCGNSPeriodicType

!        ================================================================

         logical function equalCGNSPeriodicType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * equalCGNSPeriodicType returns .true. if g1 is considered     *
!        * equal to g2, i.e. both the block and subface ID must match,  *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(cgnsPeriodicType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         equalCGNSPeriodicType = .false.
         if(g1%cgnsBlock   == g2%cgnsBlock .and. &
            g1%cgnsSubface == g2%cgnsSubface)    &
           equalCGNSPeriodicType = .true.

         end function equalCGNSPeriodicType

!        ================================================================

         logical function lessPeriodicSubfacesHaloT(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessPeriodicSubfacesHaloT returns .true. if g1 is            *
!        * considered smaller than g2.                                  *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(periodicSubfacesHaloType), intent(in) :: g1, g2
!
!        Local variables.
!
         integer(kind=intType) :: nn, i1, i2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! First compare whether or not both g1 and g2 are internal
         ! halo's. Fortran does not allow a direct comparison of
         ! logicals and therefore the integers i1 and i2 are used.

         i1 = 1; if( g1%internalHalo ) i1 = 0
         i2 = 1; if( g2%internalHalo ) i2 = 0

         if(i1 < i2) then
           lessPeriodicSubfacesHaloT = .true.
           return
         else if(i1 > i2) then
           lessPeriodicSubfacesHaloT = .false.
           return
         endif

         ! Compare the number of periodic subfaces.

         if(g1%nPeriodicSubfaces < g2%nPeriodicSubfaces) then
           lessPeriodicSubfacesHaloT = .true.
           return
         else if(g1%nPeriodicSubfaces > g2%nPeriodicSubfaces) then
           lessPeriodicSubfacesHaloT = .false.
           return
         endif

         ! The number of periodic subfaces is the same. Compare the
         ! subfaces themselves. It is assumed that the subfaces are
         ! sorted in increading order. This can be done, because the
         ! periodic transformations are commuting matrices.

         do nn=1,g1%nPeriodicSubfaces
           if(g1%periodicSubfaces(nn) < g2%periodicSubfaces(nn)) then
             lessPeriodicSubfacesHaloT = .true.
             return
           else if(g1%periodicSubfaces(nn) > g2%periodicSubfaces(nn)) then
             lessPeriodicSubfacesHaloT = .false.
             return
           endif
         enddo

         ! The periodic subfaces are identical as well. Compare the
         ! indices in the list.

         if(g1%indexInHaloList < g2%indexInHaloList) then
           lessPeriodicSubfacesHaloT = .true.
           return
         else if(g1%indexInHaloList > g2%indexInHaloList) then
           lessPeriodicSubfacesHaloT = .false.
           return
         endif

         ! Both objects are the same. Return .false.

         lessPeriodicSubfacesHaloT = .false.

         end function lessPeriodicSubfacesHaloT

!        ================================================================

         logical function lessEqualPeriodicSubfacesHaloT(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessEqualPeriodicSubfacesHaloT returns .true. if g1 is       *
!        * considered smaller than or equal to g2.                      *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(periodicSubfacesHaloType), intent(in) :: g1, g2
!
!        Local variables.
!
         integer(kind=intType) :: nn, i1, i2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! First compare whether or not both g1 and g2 are internal
         ! halo's. Fortran does not allow a direct comparison of
         ! logicals and therefore the integers i1 and i2 are used.

         i1 = 1; if( g1%internalHalo ) i1 = 0
         i2 = 1; if( g2%internalHalo ) i2 = 0

         if(i1 < i2) then
           lessEqualPeriodicSubfacesHaloT = .true.
           return
         else if(i1 > i2) then
           lessEqualPeriodicSubfacesHaloT = .false.
           return
         endif

         ! Compare the number of periodic subfaces.

         if(g1%nPeriodicSubfaces < g2%nPeriodicSubfaces) then
           lessEqualPeriodicSubfacesHaloT = .true.
           return
         else if(g1%nPeriodicSubfaces > g2%nPeriodicSubfaces) then
           lessEqualPeriodicSubfacesHaloT = .false.
           return
         endif

         ! The number of periodic subfaces is the same. Compare the
         ! subfaces themselves. It is assumed that the subfaces are
         ! sorted in increading order. This can be done, because the
         ! periodic transformations are commuting matrices.

         do nn=1,g1%nPeriodicSubfaces
           if(g1%periodicSubfaces(nn) < g2%periodicSubfaces(nn)) then
             lessEqualPeriodicSubfacesHaloT = .true.
             return
           else if(g1%periodicSubfaces(nn) > g2%periodicSubfaces(nn)) then
             lessEqualPeriodicSubfacesHaloT = .false.
             return
           endif
         enddo

         ! The periodic subfaces are identical as well. Compare the
         ! indices in the list.

         if(g1%indexInHaloList < g2%indexInHaloList) then
           lessEqualPeriodicSubfacesHaloT = .true.
           return
         else if(g1%indexInHaloList > g2%indexInHaloList) then
           lessEqualPeriodicSubfacesHaloT = .false.
           return
         endif

         ! Both objects are the same. Return .true.

         lessEqualPeriodicSubfacesHaloT = .true.

         end function lessEqualPeriodicSubfacesHaloT

!        ================================================================

         logical function equalPeriodicSubfacesHaloT(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * equalPeriodicSubfacesHaloT returns .true. if g1 is           *
!        * considered equal to g2. The equal operator is only used to   *
!        * find the different number of periodic transformations in     *
!        * determinePeriodicData. Hence only the periodic subfaces of   *
!        * the halo's are compared and g1 and g2 are considered equal   *
!        * if the subfaces are equal, even if other member variables    *
!        * differ.                                                      *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(periodicSubfacesHaloType), intent(in) :: g1, g2
!
!        Local variables.
!
         integer(kind=intType) :: nn
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         if(g1%nPeriodicSubfaces /= g2%nPeriodicSubfaces) then
           equalPeriodicSubfacesHaloT = .false.
           return
         endif

         do nn=1,g1%nPeriodicSubfaces
           if(g1%periodicSubfaces(nn) /= g2%periodicSubfaces(nn)) then
             equalPeriodicSubfacesHaloT = .false.
             return
           endif
         enddo

         equalPeriodicSubfacesHaloT = .true.

         end function equalPeriodicSubfacesHaloT

       end module periodicInfo

!
!      ******************************************************************
!      *                                                                *
!      * File:          sortSubfaces.f90                                *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 07-23-2003                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine sortSubfaces(oldSubfaceID, blockID)
!
!      ******************************************************************
!      *                                                                *
!      * sortSubfaces sorts the boundary subfaces of the given block    *
!      * such that viscous subfaces are numbered first, followed by     *
!      * inviscid, etc.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use partitionMod
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), dimension(*), intent(out) :: oldSubfaceID
       type(distributionBlockType), intent(in)          :: blockID
!
!      Local variables.
!
       integer(kind=intType) :: i, ii, nDiff

       integer(kind=intType), dimension(blockID%nBocos) :: bcPrior
       integer(kind=intType), dimension(blockID%nBocos) :: bcPriorSort

       integer(kind=intType), dimension(0:blockID%nBocos) :: mult
!
!      Function definition.
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary subfaces and determine the priorities.
       ! Store the priorities in both bcPrior and bcPriorSort.

       do i=1,blockID%nBocos

         select case (blockID%BCType(i))

           case (NSWallAdiabatic)
             bcPrior(i) = 1

           case (NSWallIsothermal)
             bcPrior(i) = 2

           case (EulerWall)
             bcPrior(i) = 3

           case (Symm)
             bcPrior(i) = 4

           case (SymmPolar)
             bcPrior(i) = 5

           case (FarField)
             bcPrior(i) = 6

           case (SupersonicInflow)
             bcPrior(i) = 7

           case (SubsonicInflow)
             bcPrior(i) = 8

           case (SupersonicOutflow)
             bcPrior(i) = 9

           case (SubsonicOutflow)
             bcPrior(i) = 10

           case (MassBleedInflow)
             bcPrior(i) = 11

           case (MassBleedOutflow)
             bcPrior(i) = 12

           case (mDot)
             bcPrior(i) = 13

           case (Thrust)
             bcPrior(i) = 14

           case (Extrap)
             bcPrior(i) = 15

           case (SlidingInterface)
             bcPrior(i) = 19

           case (OversetOuterBound)
             bcPrior(i) = 20

           case (DomainInterfaceAll)
             bcPrior(i) = 21

           case (DomainInterfaceRhoUVW)
             bcPrior(i) = 22

           case (DomainInterfaceP)
             bcPrior(i) = 23

           case (DomainInterfaceRho)
             bcPrior(i) = 24

           case (DomainInterfaceTotal)
             bcPrior(i) = 25

         end select

         bcPriorSort(i) = bcPrior(i)

       enddo

       ! Sort bcPriorSort in increasing order.

       call qsortIntegers(bcPriorSort, blockID%nBocos)

       ! Get rid of the multiple entries and store the multiplicity in
       ! cumulative storage format. nDiff contains the number of
       ! different boundary conditions for this block. The initialization
       ! with the min function is necessary to be able to treat blocks
       ! without a boundary condition correctly.

       nDiff       = min(1_intType, blockID%nBocos)
       mult(0)     = 0
       mult(nDiff) = 1

       do i=2,blockID%nBocos
         if(bcPriorSort(i) == bcPriorSort(nDiff)) then
           mult(nDiff) = mult(nDiff) + 1
         else
           nDiff              = nDiff + 1
           mult(nDiff)        = mult(nDiff-1) + 1
           bcPriorSort(nDiff) = bcPriorSort(i)
         endif
       enddo

       ! Determine the old subface ID by searching in the sorted
       ! priorities.

       do i=1,blockID%nBocos

         ii = bsearchIntegers(bcPrior(i), bcPriorSort, nDiff)

         ! Update the lower boundary in the multiplicity for this
         ! boundary condition and store the entry a bit easier.

         mult(ii-1) = mult(ii-1) + 1
         ii         = mult(ii-1)

         ! Store the mapping of the new to old subface id.

         oldSubfaceID(ii) = i

       enddo

       end subroutine sortSubfaces

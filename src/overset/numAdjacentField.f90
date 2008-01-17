!
!      ******************************************************************
!      *                                                                *
!      * File:          numAdjacentField.f90                            *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 10-04-2005                                      *
!      * Last modified: 10-04-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function numAdjacentField(level, sps, blk, i, j, k)
!
!      ******************************************************************
!      *                                                                *
!      * numAdjacentField returns the number of field or non-orphan     *
!      * fringe cells that share a face with the input cell by using    *
!      * the iblank values of the neighbors.                            *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Function type.
!
       integer(kind=intType) :: numAdjacentField
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, blk, i, j, k
!
!      Local variables.
!
       integer(kind=intType) :: l, m, ni, nj, nk

       integer(kind=intType), dimension(3) :: del
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the result to 0.

       numAdjacentField = 0

       ! Loop over the 3 coordinate directions, and for both the
       ! positive and negative direction set the delta vector to be a
       ! unit vector in that direction.

       directionLoop: do m = 1,3
         plusMinusLoop: do l = -1,1,2

           del    = 0
           del(m) = l

           ! Compute the neighbor indices and skip if it is outside the
           ! boundaries of the block.

           ni = i + del(1)
           nj = j + del(2)
           nk = k + del(3)

           if (ni < 0 .or. ni > flowDoms(blk,level,sps)%ib .or. &
               nj < 0 .or. nj > flowDoms(blk,level,sps)%jb .or. &
               nk < 0 .or. nk > flowDoms(blk,level,sps)%kb) cycle

           ! If the neighbor is a field (iblank = 1) or non-orphan fringe
           ! cell (iblank >= 10) then increment the result.

           if (flowDoms(blk,level,sps)%iblank(ni,nj,nk) == 1 .or. &
               flowDoms(blk,level,sps)%iblank(ni,nj,nk) >= 10) then
             numAdjacentField = numAdjacentField + 1
           end if

         end do plusMinusLoop
       end do directionLoop

       end function numAdjacentField

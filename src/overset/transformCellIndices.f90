!
!      ******************************************************************
!      *                                                                *
!      * File:          transformCellIndices.F90                        *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-20-2005                                      *
!      * Last modified: 04-20-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine transformCellIndices(mm, ind)
!
!      ******************************************************************
!      *                                                                *
!      * TransformCellIndices transforms the given direct halo cell     *
!      * indices into the corresponding donor cell indices using the    *
!      * subface info given by index mm and the block currently set in  *
!      * blockPointers. This is essentially a replica of what is done   *
!      * in the preprocessing routine determineFaceHalos, except for    *
!      * a single cell.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType),               intent(in)    :: mm
       integer(kind=intType), dimension(3), intent(inout) :: ind
!
!      Local variables.
!
       integer(kind=intType) :: i

       integer(kind=intType), dimension(3) :: myOffset, donorOffset
       integer(kind=intType), dimension(3) :: step

       integer(kind=intType), dimension(3,3) :: tmat
       integer(kind=intType), dimension(3,2) :: myCellRange
       integer(kind=intType), dimension(3,2) :: myNodeRange
       integer(kind=intType), dimension(3,2) :: donorCellRange
!
!      Function definition.
!
       integer(kind=intType) :: delta
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the cell and nodal range for the halo's of this
       ! subface as well as the direction normal to the subface.

       call haloRanges(mm, myOffset, myCellRange, &
                        myNodeRange, step)

       ! Determine the complete transformation matrix from the
       ! given shorthand.

       i = sign(1_intType,l1(mm))
       tmat(1,1) = i * delta(l1(mm),1_intType)
       tmat(2,1) = i * delta(l1(mm),2_intType)
       tmat(3,1) = i * delta(l1(mm),3_intType)

       i = sign(1_intType,l2(mm))
       tmat(1,2) = i * delta(l2(mm),1_intType)
       tmat(2,2) = i * delta(l2(mm),2_intType)
       tmat(3,2) = i * delta(l2(mm),3_intType)

       i = sign(1_intType,l3(mm))
       tmat(1,3) = i * delta(l3(mm),1_intType)
       tmat(2,3) = i * delta(l3(mm),2_intType)
       tmat(3,3) = i * delta(l3(mm),3_intType)

       ! Determine the offset of the donor block.

       donorOffset = matmul(tmat, myOffset)

       ! First determine the cell range of the donor block on
       ! the subface. This equals the nodal range, except that 1 is
       ! added to the smallest index. As it is possible that the
       ! index is running negatively, this should be taken into account.

       donorCellRange(1,1) = dinBeg(mm)
       donorCellRange(2,1) = djnBeg(mm)
       donorCellRange(3,1) = dknBeg(mm)

       donorCellRange(1,2) = dinEnd(mm)
       donorCellRange(2,2) = djnEnd(mm)
       donorCellRange(3,2) = dknEnd(mm)

       ! The loop to add 1 to the lowest index and to correct the
       ! index corresponding to the face we are on.

       do i=1,3
         if(donorCellRange(i,1) == donorCellRange(i,2)) then

           ! If the face corresponds to a min face, indicated by
           ! donorCellRange(i,1) == 1 then 1 must be added;
           ! otherwise nothing needs to be done.

           if(donorCellRange(i,1) == 1) then
             donorCellRange(i,1) = 2
             donorCellRange(i,2) = 2
           endif

         else if(donorCellRange(i,1) > donorCellRange(i,2)) then
           donorCellRange(i,2) = donorCellRange(i,2) + 1
         else
           donorCellRange(i,1) = donorCellRange(i,1) + 1
         endif
       enddo

       ! Determine the indices of the donor point by applying
       ! the transformation matrix to i,j,k.

       step = ind - myCellRange(:,1)

       ind = donorCellRange(:,1) + matmul(tmat, step)

       end subroutine transformCellIndices

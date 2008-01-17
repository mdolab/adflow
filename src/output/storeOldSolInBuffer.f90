!
!      ******************************************************************
!      *                                                                *
!      * File:          storeOldSolInBuffer.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-08-2004                                      *
!      * Last modified: 06-29-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine storeOldSolInBuffer(buffer, ind, wID, &
                                      iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
!
!      ******************************************************************
!      *                                                                *
!      * storeOldSolInBuffer stores the given range of the wID'th       *
!      * conservative variable of an old solution in buffer. Needed for *
!      * a time accurate restart. It is assumed that the variables in   *
!      * blockPointers already point to the correct block.              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: ind, wID
       integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType), intent(in) :: kBeg, kEnd

       real(kind=realType), dimension(*), intent(out) :: buffer
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, nOld, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the index in wOld a bit easier.

       nOld = ind - 1

       ! Loop over the cell range of the block and copy the wID'th
       ! variable in buffer.

       nn = 0

       do k=kBeg,kEnd
         do j=jBeg,jEnd
           do i=iBeg,iEnd
             nn = nn + 1; buffer(nn) = wOld(nOld,i,j,k,wID)
           enddo
         enddo
       enddo

       end subroutine storeOldSolInBuffer

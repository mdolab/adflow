!
!      ******************************************************************
!      *                                                                *
!      * File:          copyDataBuf.90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-12-2005                                      *
!      * Last modified: 10-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine copyDataBufSinglePrecision(val, buffer,      &
                                             iBeg, jBeg, kBeg, &
                                             iEnd, jEnd, kEnd, subRange)
!
!      ******************************************************************
!      *                                                                *
!      * copyDataBufSinglePrecision stores the given 1D buffer into the *
!      * subrange of the 3D single precision val array.                 *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBeg, jBeg, kBeg
       integer(kind=intType), intent(in) :: iEnd, jEnd, kEnd
       integer(kind=intType), dimension(3,2), intent(in) :: subRange

       real(kind=realType), dimension(*), intent(in) :: buffer
       real(kind=4), dimension(iBeg:iEnd,jBeg:jEnd,kBeg:kEnd), &
                                                    intent(inout) :: val
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the subrange into val.

       ll = 0
       do k=subRange(3,1), subRange(3,2)
         do j=subRange(2,1), subRange(2,2)
           do i=subRange(1,1), subRange(1,2)
             ll = ll + 1
             val(i,j,k) = buffer(ll)
           enddo
         enddo
       enddo

       end subroutine copyDataBufSinglePrecision

!      ==================================================================

       subroutine copyDataBufDoublePrecision(val, buffer,      &
                                             iBeg, jBeg, kBeg, &
                                             iEnd, jEnd, kEnd, subRange)
!
!      ******************************************************************
!      *                                                                *
!      * copyDataBufDoublePrecision stores the given 1D buffer into the *
!      * subrange of the 3D double precision val array.                 *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBeg, jBeg, kBeg
       integer(kind=intType), intent(in) :: iEnd, jEnd, kEnd
       integer(kind=intType), dimension(3,2), intent(in) :: subRange

       real(kind=realType), dimension(*), intent(in) :: buffer
       real(kind=8), dimension(iBeg:iEnd,jBeg:jEnd,kBeg:kEnd), &
                                                    intent(inout) :: val
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the subrange into val.

       ll = 0
       do k=subRange(3,1), subRange(3,2)
         do j=subRange(2,1), subRange(2,2)
           do i=subRange(1,1), subRange(1,2)
             ll = ll + 1
             val(i,j,k) = buffer(ll)
           enddo
         enddo
       enddo

       end subroutine copyDataBufDoublePrecision

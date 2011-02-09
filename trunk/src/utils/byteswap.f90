!
!      ******************************************************************
!      *                                                                *
!      * File:          byteswap.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-16-2005                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine byteswap(buffer, baseSize, nn)
!
!      ******************************************************************
!      *                                                                *
!      * byteswap swaps the byte sequence of the given buffer. The size *
!      * of each entity is given by baseSize and the number of          *
!      * elements by nn. This means that buffer contains nn*baseSize    *
!      * bytes.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       character, dimension(*), intent(inout) :: buffer
       integer(kind=intType),  intent(in)     :: baseSize, nn
!
!      Local variables.
!
       integer(kind=intType) :: ii, jj, kk, i, j
       character             :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the half of baseSize in kk.

       kk = baseSize/2

       ! Loop over the number of elements.

       do j=1,nn

         ! Initialize ii and jj, which are used to store the indices of
         ! the bytes.

         ii = (j-1)*baseSize
         jj = ii + baseSize +1

         ! Loop over half the number of bytes in the entity.

         do i=1,kk

           ! Store the indices of the bytes to be swapped
           ! and swap them.

           ii = ii +1
           jj = jj -1

           tmp        = buffer(jj)
           buffer(jj) = buffer(ii)
           buffer(ii) = tmp
         enddo
       enddo

       end subroutine byteswap

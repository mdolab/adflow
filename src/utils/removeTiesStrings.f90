!
!      ******************************************************************
!      *                                                                *
!      * File:          removeTiesStrings.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-28-2005                                      *
!      * Last modified: 04-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine removeTiesStrings(array, mSizes, nn)
!
!      ******************************************************************
!      *                                                                *
!      * removeTiesStrings removes the ties in the given sorted         *
!      * array and stores the number of ties per entry in cumulative    *
!      * storage format in mSizes. On entry nn contains the number of   *
!      * elements in array; on exit the number of different elements in *
!      * array.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), dimension(*), intent(out) :: mSizes
       integer(kind=intType), intent(inout)             :: nn

       character(len=*), dimension(*), intent(inout) :: array
!
!      Local variables
!
       integer(kind=intType) :: i, ii
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the first element of mSizes to 0.

       mSizes(1) = 0

       ! If the size of the array is 0, return

       if(nn == 0) return

       ! nn is at least 1, so both mSizes(2) and ii (the number of
       ! different entities in array) can be set to 1.

       mSizes(2) = 1
       ii        = 1

       ! Loop over the rest of the array.

       do i=2,nn

         if(array(i) /= array(ii)) then

           ! The array elements differ. Update the number of different
           ! entities, ii, store the value of the array element at the
           ! appropriate position and initialize the new value of mSizes
           ! to the previous value plus 1. Rememember that mSizes is
           ! stored in cumulative storage format.

           ii = ii +1
           array(ii) = array(i)
           mSizes(ii+1) = mSizes(ii) +1

         else

           ! Both array elements are the same. Simply update mSizes.

           mSizes(ii+1) = mSizes(ii+1) + 1
         endif

       enddo

       ! Set nn to ii, the number of different entities in array.

       nn = ii

       end subroutine removeTiesStrings

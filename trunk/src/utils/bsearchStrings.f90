!
!      ******************************************************************
!      *                                                                *
!      * File:          bsearchStrings.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-19-2002                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function bsearchStrings(key, base, nn)
!
!      ******************************************************************
!      *                                                                *
!      * bsearchStrings returns the index in base where key is stored.  *
!      * A binary search algorithm is used here, so it is assumed that  *
!      * base is sorted in increasing order. In case key appears more   *
!      * than once in base, the result is arbitrary. If key is not      *
!      * found, a zero is returned.                                     *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Function type
!
       integer(kind=intType) :: bsearchStrings
!
!      Function arguments.
!
       character(len=*), intent(in)               :: key
       character(len=*), dimension(*), intent(in) :: base
       integer(kind=intType)                      :: nn
!
!      Local variables.
!
       integer(kind=intType) :: ii, pos, start
       logical               :: entryFound
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize some values.

       start       = 1
       ii          = nn
       entryFound = .false.

       ! Binary search to find key.

       do
         ! Condition for breaking the loop

         if(ii == 0) exit

         ! Determine the position in the array to compare.

         pos = start + ii/2

         ! In case this is the entry, break the search loop.

         if(base(pos) == key) then
           entryFound = .true.
           exit
         endif

         ! In case the search key is larger than the current position,
         ! only parts to the right must be searched. Remember that base
         ! is sorted in increasing order. Nothing needs to be done if the
         ! key is smaller than the current element.

         if(key > base(pos)) then
           start = pos +1
           ii    = ii -1
         endif

         ! Modify ii for the next branch to search.

         ii = ii/2
       enddo

       ! Set bsearchStrings. This depends whether the key was found.

       if( entryFound ) then
         bsearchStrings = pos
       else
         bsearchStrings = 0
       endif

       end function bsearchStrings

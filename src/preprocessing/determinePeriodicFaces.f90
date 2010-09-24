!
!      ******************************************************************
!      *                                                                *
!      * File:          determinePeriodicFaces.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-10-2003                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determinePeriodicFaces
!
!      ******************************************************************
!      *                                                                *
!      * determinePeriodicFaces determines and stores the number of     *
!      * periodic faces present in the complete mesh. The sequence of   *
!      * storing the data is such that the array periodicGlobal is      *
!      * sorted with the definition of the < operator for this datatype.*
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use periodicInfo
       implicit none
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn, ii, i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of periodic faces present in the cgns grid.

       nPeriodicGlobal = 0
       do nn=1,cgnsNDom
         do i=1,cgnsDoms(nn)%n1to1
           if( cgnsDoms(nn)%conn1to1(i)%periodic ) &
             nPeriodicGlobal = nPeriodicGlobal + 1
         enddo
       enddo

       ! Allocate the memory for periodicGlobal.

       allocate(periodicGlobal(nPeriodicGlobal), stat=ierr)
       if(ierr /= 0)                              &
         call terminate("determinePeriodicFaces", &
                        "Memory allocation failure for periodicGlobal")

       ! Repeat the loop over the faces of the cgns grid and store the
       ! periodic faces.

       ii = 0
       do nn=1,cgnsNDom
         do i=1,cgnsDoms(nn)%n1to1
           if( cgnsDoms(nn)%conn1to1(i)%periodic ) then

             ii = ii + 1
             periodicGlobal(ii)%cgnsBlock   = nn
             periodicGlobal(ii)%cgnsSubface = i

           endif
         enddo
       enddo

       end subroutine determinePeriodicFaces

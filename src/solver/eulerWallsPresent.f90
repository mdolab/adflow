!
!      ******************************************************************
!      *                                                                *
!      * File:          EulerWallsPresent.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-08-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       logical function EulerWallsPresent()
!
!      ******************************************************************
!      *                                                                *
!      * eulerWallsPresent determines whether or not inviscid walls     *
!      * are present in the whole grid. It first determines if these    *
!      * walls are present locally and performs an allReduce            *
!      * afterwards.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use BCTypes
       use communication
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn, i
       integer               :: ierr
       logical               :: localEulerWalls
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize localEulerWalls to .false. and loop over the
       ! boundary subfaces of the blocks to see if Euler walls are
       ! present on this processor. As the info is the same for all
       ! spectral solutions, only the 1st needs to be considered.

       localEulerWalls = .false.
       do nn=1,nDom
         do i=1,flowDoms(nn,1,1)%nBocos
           if(flowDoms(nn,1,1)%BCType(i) == EulerWall) &
             localEulerWalls = .true.
         enddo
       enddo

       ! Set i to 1 if Euler walls are present locally and to 0
       ! otherwise. Determine the maximum over all processors
       ! and set EulerWallsPresent accordingly.

       i = 0
       if( localEulerWalls ) i = 1
       call mpi_allreduce(i, nn, 1, sumb_integer, mpi_max, &
                          SUmb_comm_world, ierr)

       if(nn == 0) then
         EulerWallsPresent = .false.
       else
         EulerWallsPresent = .true.
       endif

       end function EulerWallsPresent

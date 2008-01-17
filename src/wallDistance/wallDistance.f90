!
!      ******************************************************************
!      *                                                                *
!      * File:          wallDistance.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-27-2003                                      *
!      * Last modified: 11-03-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine wallDistance(level, allocMem)
!
!      ******************************************************************
!      *                                                                *
!      * wallDistance computes the distances of the cell centers to     *
!      * the nearest viscous wall. An adt type of method is used, which *
!      * guarantees to find the minimum distance to the wall. Possible  *
!      * periodic transformations are taken into account, such that     *
!      * also in case of a periodic problem the correct distance is     *
!      * computed; the nearest wall point may lie in a periodic domain. *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use inputPhysics
       use inputTimeSpectral
       use viscSurface
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
       logical, intent(in)               :: allocMem
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: sps

       double precision :: t0

       character(len=3) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if the RANS equations are solved. If not, the wall
       ! distance is not needed and a return can be made.

       if(equations /= RANSEquations) return

       ! If the turbulence model is wall distance free just compute the
       ! normal spacing of the first cell and store this in the wall
       ! distance. It may be needed for the boundary conditions and
       ! the monitoring of the y+. Return afterwards.

       if(.not. wallDistanceNeeded) then

         ! Loop over the number of spectral solutions, initialize the
         ! distance and compute the initial normal spacing.

         do sps=1,nTimeIntervalsSpectral
           call initWallDistance(level, sps, allocMem)
           call computeNormalSpacing(level, sps)
         enddo

         ! And return.

         return
       endif

       ! Write a message to stdout to indicate that the wall distance
       ! computation starts for the given level.

       write(integerString,"(i2)") level
       integerString = adjustl(integerString)
       if(myID == 0) then
         print 101, trim(integerString)
       endif
 101   format("# Start wall distances level",1X,A)

       ! Store the start time.

       t0 = mpi_wtime()

       ! Release temporarily some memory such that the overall memory
       ! requirement is not dictated by this routine. What memory is
       ! released depends on allocMem. If this is .True. It means that
       ! this is the first time the wall distances are computed, i.e. in
       ! the preprocessing phase. Then only send and receive buffers are
       ! released. If allocMem is .False., this means that the wall
       ! distances are computed in a moving mesh computation and
       ! consequently some more memory can be released.

       if( allocMem ) then
         deallocate(sendBuffer, recvBuffer, stat=ierr)
         if(ierr /= 0)                    &
           call terminate("wallDistance", &
                          "Deallocation error for communication buffers")
       else
         call deallocateTempMemory(.false.)
       endif

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Initialize the wall distances.

         call initWallDistance(level, sps, allocMem)

         ! Build the viscous surface mesh.

         call viscousSurfaceMesh(level, sps)

         ! If there are no viscous faces, processor 0 prints a warning
         ! and the wall distances are not computed.

         if(nquadViscGlob == 0) then

           if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#               Warning!!!!"
             print "(a)", "# No viscous boundary found. Wall &
                          &distances are set to infinity"
             print "(a)", "#"
           endif

         else

           ! Determine the wall distances for the owned cells.

           call determineDistance(level, sps)

         endif

       enddo spectralLoop

       ! Allocate the temporarily released memory again. For more info
       ! see the comments at the beginning of this routine.

       if( allocMem ) then
         allocate(sendBuffer(sendBufferSize), &
                  recvBuffer(recvBufferSize), stat=ierr)
         if(ierr /= 0)                    &
           call terminate("wallDistance", &
                          "Memory allocation failure for comm buffers")
       else
         call allocateTempMemory(.false.)
       endif

       ! Synchronize the processors.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! Write a message to stdout with the amount of time it
       ! took to compute the distances.

       if(myID == 0) then
         print 102, trim(integerString)
         print 103, mpi_wtime() - t0
         print "(a)", "#"
       endif
 102   format("# End wall distances level",1X,A)
 103   format("# Wall clock time:",E12.5," sec.")

       end subroutine wallDistance

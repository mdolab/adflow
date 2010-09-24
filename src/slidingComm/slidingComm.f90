!
!      ******************************************************************
!      *                                                                *
!      * File:          slidingComm.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine slidingComm(level, firstTime)
!
!      ******************************************************************
!      *                                                                *
!      * slidingComm determines the communication pattern for the       *
!      * sliding mesh interfaces on the given grid level.               *
!      * In case of a steady computation the mixing plane approach is   *
!      * used.                                                          *
!      *                                                                *
!      ******************************************************************
!
       use block
       use commSliding
       use communication
       use inputPhysics
       use inputTimeSpectral
       use interfaceGroups
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
       logical, intent(in)               :: firstTime
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Release the memory of the entire communication pattern if this
       ! is not the first time this routine is called and initialize
       ! it again afterwards.

       if(.not. firstTime) call releaseMemSliding(level)
       call initMemSliding(level)

       ! Perform the initialization for the mixing planes.

       call initMemMixingPlane(level)

       ! Release temporarily some memory such that the overall memory
       ! requirement is not dictated by this routine. What memory is
       ! released depends on firstTime. If firstTime is .True. This
       ! means that the routine is called in the preprocessing phase
       ! and thus only the send and receive buffer can be released.
       ! Otherwise this routine is called in a moving mesh computation
       ! and some more memory can be released.

       if( firstTime ) then
         deallocate(sendBuffer, recvBuffer, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingComm", &
                          "Deallocation error for communication buffers")
       else
         call deallocateTempMemory(.false.)
       endif

       ! Determine for every cell and the 1st level halo cells what type
       ! of halo cells they are. This info is identical for all spectral
       ! solutions and therefore it is determined only once.

       call haloType(level)

       ! Determine the situation we have here.

       select case(equationMode)

         case (steady)

           ! Steady state computation. The mixing plane assumption
           ! is used.

           ! Loop over the number of sliding interface groups, such
           ! that a processor works at most on 1 interface at the
           ! same time.

           do ii=1,nInterfaceGroups

             ! Determine the communication/interpolation pattern
             ! for the interface in this group.

             call mixingPlaneComm(level,ii)

             ! Block until all processors have reached this point.

             call mpi_barrier(SUmb_comm_world, ierr)

           enddo

           ! Although no sliding mesh interpolation needs to be
           ! performed, the arrays need to be initialized.

           sps = 1
           call cumulativeNSendReceives(commSlidingCell_1st(level,sps))
           call cumulativeNSendReceives(commSlidingCell_2nd(level,sps))

           ! Determine the new size of the communication buffers.

           call setBufferSizes(level, sps, .false., .true., .false.)

         !===============================================================

         case (unsteady, timeSpectral)

           ! An unsteady computation. The sliding mesh communication
           ! pattern must be determined.

           ! Loop over the number of spectral solutions.

           spectralLoop: do sps=1,nTimeIntervalsSpectral

             ! Loop over the number of sliding interface groups, such
             ! that a processor works at most on 1 interface at the
             ! same time.

             do ii=1,nInterfaceGroups

               ! Update the communication pattern for the sliding mesh
               ! interface in this group.

               call slidingMesh(level,sps,ii)

               ! Block until all processors have reached this point.

               call mpi_barrier(SUmb_comm_world, ierr)

             enddo

             ! Determine the cumulative variants of the send and receive
             ! arrays in the external communication pattern for both the
             ! 1st and 2nd level cell halo's.

             call cumulativeNSendReceives(commSlidingCell_1st(level,sps))
             call cumulativeNSendReceives(commSlidingCell_2nd(level,sps))

             ! Determine the new size of the communication buffers.

             call setBufferSizes(level, sps, .false., .true., .false.)

           enddo spectralLoop

       end select

       ! Release the memory of the halo info.

       do ii=1,nDom
         deallocate(donorDoms(ii)%haloInfo, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingComm", &
                          "Deallocation error for haloInfo")
       enddo

       deallocate(donorDoms, stat=ierr)
       if(ierr /= 0)                   &
         call terminate("slidingComm", &
                        "Deallocation error for donorDoms")

       ! Allocate the temporarily released memory again. For more info
       ! see the comments at the beginning of this routine.

       if( firstTime ) then
         allocate(sendBuffer(sendBufferSize), &
                  recvBuffer(recvBufferSize), stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingComm", &
                          "Memory allocation failure for comm buffers")
       else
         call allocateTempMemory(.false.)
       endif

       end subroutine slidingComm

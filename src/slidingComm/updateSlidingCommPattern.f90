!
!      ******************************************************************
!      *                                                                *
!      * File:          updateSlidingCommPattern.f90                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateSlidingCommPattern(level, sps, color)
!
!      ******************************************************************
!      *                                                                *
!      * updateSlidingCommPattern updates the sliding mesh              *
!      * communication patterns for the 1st and 2nd level cell halo's   *
!      * for the active sliding interface.                              *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use interfaceGroups
       use localSubfacesMod
       use thisSlide
       use updateComm
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, color
!
!      Local variables.
!
       integer :: ierr
       integer :: nProcSlide, myIDSlide, commSlide
       integer :: size, procID

       integer, dimension(mpi_status_size) :: status
       integer, dimension(0:myInterfaces(color)%nProcSlide-1) :: trueProc

       integer(kind=intType) :: i
       integer(kind=intType) :: nn, mm, ii, jj, kk, ll
       integer(kind=intType) :: nMessagesSend, nMessagesRecv

       integer(kind=intType), &
        dimension(0:myInterfaces(color)%nProcSlide) :: &
           nDirToProc, nDirFromProc

       integer(kind=intType), &
        dimension(0:myInterfaces(color)%nProcSlide-1) :: &
           countInt, countReal, countIndices

       integer(kind=intType), dimension(:), allocatable :: intSend
       integer(kind=intType), dimension(:), allocatable :: tmpInt
       integer(kind=intType), dimension(:), allocatable :: indicesSend

       real(kind=realType), dimension(:), allocatable :: realSend
       real(kind=realType), dimension(:), allocatable :: tmpReal

       type(updateCommType), &
        dimension(0:myInterfaces(color)%nProcSlide-1) :: localInterpol
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Abbreviate the number of processors, my processor id and the
       ! communicator for the subgroup of this sliding interface.

       nProcSlide = myInterfaces(color)%nProcSlide
       myIDSlide  = myInterfaces(color)%myIDSlide
       commSlide  = myInterfaces(color)%commSlide

       ! Initialize the number of entities in localInterpol to 0.

       do nn=0,(nProcSlide-1)
         localInterpol(nn)%nCopy1st = 0
         localInterpol(nn)%nCopy2nd = 0
       enddo

       ! Make sure that the processor ID's of the SUmb_comm_world are
       ! known for all processors in the group commSlide.

       call mpi_allgather(myID, 1, mpi_integer, trueProc, 1, &
                         mpi_integer, commSlide, ierr)

       ! Synchronize the processors of this group. This is done, because
       ! wild cards are used in the receives. It is probably not
       ! necessary, but it is just as safe to do it.

       call mpi_barrier(commSlide, ierr)

       ! Determine the number of cells to be sent to other processors
       ! for which a corresponding donor was searched in the surface
       ! grid.

       nDirToProc = 0
       call getNCell2Proc(nMySubfaces1, mySubfaces1)
       call getNCell2Proc(nMySubfaces2, mySubfaces2)

       ! Determine the number of cells I will receive from every
       ! processor.

       nDirFromProc(0) = 0

       call mpi_alltoall(nDirToProc(1),   1, sumb_integer, &
                         nDirFromProc(1), 1, sumb_integer, &
                         commSlide, ierr)

       ! Put all the arrays nDirToProc and nDirFromProc in cumulative
       ! storage format.

       do i=1,nProcSlide
         nDirToProc(i)   = nDirToProc(i)   + nDirToProc(i-1)
         nDirFromProc(i) = nDirFromProc(i) + nDirFromProc(i-1)
       enddo

       ! Allocate the memory for both the send buffers and for the
       ! integer array to store some additional info.

       nn = 2*nDirToProc(nProcSlide)
       mm = 2*nDirToProc(nProcSlide)
       ii = 4*nDirToProc(nProcSlide)

       allocate(intSend(nn), realSend(mm), indicesSend(ii), &
                stat=ierr)
       if(ierr /= 0)                                &
         call terminate("updateSlidingCommPattern", &
                        "Memory allocation failure for intSend, etc.")

       ! Initialize the counters, which store the addresses for the
       ! buffers just allocated.

       do i=0,(nProcSlide-1)
         countInt(i)     = 2*nDirToProc(i)
         countReal(i)    = 2*nDirToProc(i)
         countIndices(i) = 4*nDirToProc(i)
       enddo

       ! Fill the send buffers intSend and realSend and the array
       ! indicesSend. The latter stores the indices of the halo cell
       ! as well as the info whether it is 1st or 2nd halo.
       ! The last argument, 1 or -1, is there to indicate whether or not
       ! to negate the donorDualQ ID's. In this way a distinction can be
       ! made between subfaces from part 1 and from part 2. If negative
       ! the donor quad is part of part 1.

       call dirIntAndRealBuffer(nMySubfaces1, mySubfaces1,  1_intType)
       call dirIntAndRealBuffer(nMySubfaces2, mySubfaces2, -1_intType)

       ! Release the memory of the arrays of mySubfaces, which are
       ! not needed anymore.

       do nn=1,nMySubfaces1
         deallocate(mySubfaces1(nn)%searchQuad, &
                    mySubfaces1(nn)%donorDualQ, &
                    mySubfaces1(nn)%donorProc,  &
                    mySubfaces1(nn)%u,          &
                    mySubfaces1(nn)%v, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation error for mySubfaces1.")
       enddo

       do nn=1,nMySubfaces2
         deallocate(mySubfaces2(nn)%searchQuad, &
                    mySubfaces2(nn)%donorDualQ, &
                    mySubfaces2(nn)%donorProc,  &
                    mySubfaces2(nn)%u,          &
                    mySubfaces2(nn)%v, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation error for mySubfaces2.")
       enddo

       ! Send the data I have to send to other processors.
       ! Use nonblocking sends to avoid deadlock.

       nMessagesSend = 0
       sendLoop1: do nn=0,(nProcSlide-1)

         ! Do not send a message to myself.

         testMyself1: if(nn /= myIDSlide) then

           ! Test if something must be sent to processor nn.

           if(nDirToProc(nn+1) > nDirToProc(nn)) then

             ! Update nMessagesSend. Note that nMessagesSend actually
             ! contains half the number of messages send; the combination
             ! of the integer and real messages is interpreted as
             ! 1 message.

             nMessagesSend = nMessagesSend + 1

             ! Determine the starting index in intSend and the size of
             ! the message.

             i  = 2*nDirToProc(nn)
             mm = 2*nDirToProc(nn+1) - i
             i  = i + 1

             ! Send the integer information.

             size   = mm
             procID = nn

             call mpi_isend(intSend(i), size, sumb_integer, &
                            procID, procID, commSlide,     &
                            sendRequests(nMessagesSend), ierr)

             ! Determine the starting index in realSend and the size of
             ! the message.

             i  = 2*nDirToProc(nn)
             mm = 2*nDirToProc(nn+1) - i
             i  = i + 1

             ! Send the real information.

             size    = mm
             procID = nn

             call mpi_isend(realSend(i), size, sumb_real,   &
                            procID, procID+1, commSlide,   &
                            recvRequests(nMessagesSend), ierr)
           endif

         endif testMyself1

       enddo sendLoop1

       ! Check if local data must be interpolated.

       i = myIDSlide
       if(nDirToProc(i+1) > nDirToProc(i)) then

         ! There is local data to be interpolated. Determine the starting
         ! indices in intSend and realSend for the local data, determine
         ! the number of local halo cells and call the volume
         ! interpolation routine.

         nn = 2*nDirToProc(i) + 1
         mm = 2*nDirToProc(i) + 1

         ii = nDirToProc(i+1) - nDirToProc(i)

         call volumeInterpol(intSend(nn), realSend(mm), &
                             localInterpol(i), ii, level, sps)
       endif

       ! Determine the number of messages I will receive. Remember that
       ! no message from myself is received.

       nMessagesRecv = 0
       do nn=0,(nProcSlide-1)
         if(nn /= myIDSlide) then
           if(ndirFromProc(nn+1) > ndirFromProc(nn)) &
             nMessagesRecv = nMessagesRecv + 1
         endif
       enddo

       ! Loop over the number of messages I must receive.

       recvLoop1: do nn=1,nMessagesRecv

         ! Block until a message arrives and determine the processor ID.

         call mpi_probe(mpi_any_source, myIDSlide, commSlide, &
                       status, ierr)
         procID = status(mpi_source)
         i      = procID

         ! Determine the number of halo cells coming from that processor.

         ii = ndirFromProc(i+1) - ndirFromProc(i)

         ! Allocate the memory for the receive buffers.

         allocate(tmpInt(2*ii), tmpReal(2*ii), stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Memory allocation failure for tmpInt &
                          &and tmpReal.")

         ! Receive both messages. Use blocking sends, because they will
         ! be needed right away and the messages have already arrived.

         size = 2*ii
         call mpi_recv(tmpInt, size, sumb_integer, procID, &
                       myIDSlide, commSlide, status, ierr)

         size = 2*ii
         call mpi_recv(tmpReal, size, sumb_real, procID, &
                       myIDSlide+1, commSlide, status, ierr)

         ! Determine the interpolation for these halo cells.

         call volumeInterpol(tmpInt, tmpReal, localInterpol(i), &
                             ii, level, sps)

         ! Release the memory of tmpInt and tmpReal.

         deallocate(tmpInt, tmpReal, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation failure for tmpInt &
                          &and tmpReal.")
       enddo recvLoop1

       ! Release the memory of the surface grids of the sliding mesh
       ! interfaces which is still allocated. Note that this memory
       ! was always allocated, albeit with size zero if not needed.

       deallocate(subface1, quadID1, &
                  subface2, quadID2, stat=ierr)
       if(ierr /= 0) call terminate("updateSlidingCommPattern",  &
                                    "Deallocation error for subface1, &
                                    &quadID1, subface2 and quadID2.")

       ! Deallocate the memory of the local subfaces which has not been
       ! deallocated yet. This data structure is no longer needed.

       do nn=1,nMySubfaces1
         deallocate(mySubfaces1(nn)%indHalo1, &
                    mySubfaces1(nn)%indHalo2, &
                    mySubfaces1(nn)%connDual, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation error for the remaining arrays &
                          &mySubfaces1.")
       enddo

       do nn=1,nMySubfaces2
         deallocate(mySubfaces2(nn)%indHalo1, &
                    mySubfaces2(nn)%indHalo2, &
                    mySubfaces2(nn)%connDual, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation error for he remaining arrays &
                          &mySubfaces2.")
       enddo

       deallocate(mySubfaces1, mySubfaces2, stat=ierr)
       if(ierr /= 0)                                &
         call terminate("updateSlidingCommPattern", &
                        "Deallocation error for mySubfaces1 &
                        &and mySubfaces2.")

       ! Complete the nonblocking sends.

       size = nMessagesSend
       do nn=1,nMessagesSend
         call mpi_waitany(size, sendRequests, procID, status, ierr)
         call mpi_waitany(size, recvRequests, procID, status, ierr)
       enddo

       ! Release the memory of intSend and realSend.

       deallocate(intSend, realSend, stat=ierr)
       if(ierr /= 0)                                &
         call terminate("updateSlidingCommPattern", &
                        "Deallocation failure for intSend &
                        &and realSend.")

       ! Determine the new size of the integer and real send buffers.

       jj = 0
       do nn=0,(nProcSlide-1)
         if(nn /= myIDSlide) jj = jj + localInterpol(nn)%nCopy2nd
       enddo

       ii = 2*jj + 4*nMessagesRecv

       ! Allocate the memory for the buffers intSend and realSend.

       allocate(intSend(ii), realSend(jj), stat=ierr)
       if(ierr /= 0)                                &
         call terminate("updateSlidingCommPattern", &
                        "Memory allocation failure for intSend &
                        &and realSend.")

       ! Send the interpolation data back to the correct processors.
       ! Again use nonblocking sends to avoid deadlock.

       ii = 1
       jj = 1
       kk = 0
       sendLoop2: do i=0,(nProcSlide-1)

         ! Do not send a message to myself.

         testMyself2: if(i /= myIDSlide) then

           ! Test if something must be sent to processor nn.

           if(localInterpol(i)%nCopy2nd > 0) then

             ! Store the update of the interpolation information
             ! in the buffers.

             ll = trueProc(i)
             call updateInterpolSendBuf(intSend(ii), realSend(jj), &
                                        localInterpol(i), level,   &
                                        sps, ll)

             ! Send the integer data.

             mm     = localInterpol(i)%nCopy2nd
             nn     = 2*mm + 4
             size   = nn
             procID = i
             kk      = kk + 1
             call mpi_isend(intSend(ii), size, sumb_integer, &
                            procID, procID+2, commSlide,    &
                            sendRequests(kk), ierr)

             ! Send the real data.

             size = mm
             call mpi_isend(realSend(jj), size, sumb_real, &
                            procID, procID+3, commSlide,  &
                            recvRequests(kk), ierr)

             ! Update the indices ii and jj.

             ii = ii + nn
             jj = jj + mm

           endif
         endif testMyself2
       enddo sendLoop2

       ! Check if local interpolation data must be stored.

       i = myIDSlide
       if(nDirToProc(i+1) > nDirToProc(i)) then

         ! Store the start index in indicesSend and call the
         ! routine updateLocalCommSlide.

         nn = 4*nDirToProc(i) + 1

         call updateLocalCommSlide(indicesSend(nn), &
                                   localInterpol(i), level, sps)

       endif

       ! Loop over the number of double messages I must receive.

       recvLoop2: do nn=1,nMessagesSend

         ! Block until a message arrives and determine the processor id
         ! and the size of the integer message.

         call mpi_probe(mpi_any_source, myIDSlide+2, commSlide, &
                       status, ierr)
         procID = status(mpi_source)

         call mpi_get_count(status, sumb_integer, size, ierr)
         ii = size

         ! Allocate the memory for tmpInt and receive the message.

         allocate(tmpInt(ii), stat=ierr)
         if(ierr /= 0) &
           call terminate("updateSlidingCommPattern", &
                          "Memory allocation failure for tmpInt.")

         call mpi_recv(tmpInt, size, sumb_integer, procID, &
                       myIDSlide+2, commSlide, status, ierr)

         ! Determine the size of the corresponding real message.

         call mpi_probe(procID, myIDSlide+3, commSlide, &
                        status, ierr)
         call mpi_get_count(status, sumb_real, size, ierr)
         ii = size

         ! Allocate the memory for tmpReal and receive the message.

         allocate(tmpReal(ii), stat=ierr)
         if(ierr /= 0) &
           call terminate("updateSlidingCommPattern", &
                          "Memory allocation failure for tmpReal.")

         call mpi_recv(tmpReal, size, sumb_real, procID, &
                       myIDSlide+3, commSlide, status, ierr)

         ! Determine the starting index in indicesSend and call the
         ! routine updateExternalCommSlide.

         i = procID
         mm = 4*nDirToProc(i) + 1

         ll = trueProc(i)
         call updateExternalCommSlide(indicesSend(mm), tmpInt, &
                                      tmpReal, level, sps, ll)

         ! Release the memory of the receive buffers.

         deallocate(tmpInt, tmpReal, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("updateSlidingCommPattern", &
                          "Deallocation error for tmpInt and tmpReal.")

       enddo recvLoop2

       ! Complete the nonblocking sends.

       size = nMessagesRecv
       do nn=1,nMessagesRecv
         call mpi_waitany(size, sendRequests, procID, status, ierr)
         call mpi_waitany(size, recvRequests, procID, status, ierr)
       enddo

       ! Release the memory for intSend, realSend and indicesSend.

       deallocate(intSend, realSend, indicesSend, stat=ierr)
       if(ierr /= 0)                                &
         call terminate("updateSlidingCommPattern", &
                        "Deallocation error for intSend, &
                        &realSend and indicesSend.")

       ! Synchronize the processors. Just to be sure, because wild
       ! cards have been used in the receives.

       call mpi_barrier(commSlide, ierr)

       !=================================================================

       contains

         !===============================================================

         subroutine getNCell2Proc(nMySubfaces, mySubfaces)
!
!        ****************************************************************
!        *                                                              *
!        * getNCell2Proc increments the number of cells sent to the     *
!        * processors for interpolation.                                *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nMySubfaces
         type(localSubfaceType), dimension(nMySubfaces), &
                                                intent(in) :: mySubfaces
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the subfaces and the quadrilateral faces of these
         ! subfaces.

         do nn=1,nMySubfaces
           do mm=1,mySubfaces(nn)%nQuad

             ! Check if a halo should be constructed at all.

             if( mySubfaces(nn)%searchQuad(mm) ) then

               ! Either a direct halo or a corner halo. In the former
               ! case 2 cells will be sent for interpolation, in the
               ! latter only 1. A corner halo is indicated by a -1 for
               ! the first index of the second level halo cell. The
               ! offset + 1 in the processor ID is present, because
               ! nDirToProc will be put in cumulative storage format
               ! later on.

               ii = mySubfaces(nn)%donorProc(mm) + 1
               if(mySubfaces(nn)%indHalo2(mm,1) == -1) then
                 nDirToProc(ii) = nDirToProc(ii) + 1
               else
                 nDirToProc(ii) = nDirToProc(ii) + 2
               endif

             endif
           enddo
         enddo

         end subroutine getNcell2Proc

         !===============================================================

         subroutine dirIntAndRealBuffer(nMySubfaces, mySubfaces, &
                                        signSubface)
!
!        ****************************************************************
!        *                                                              *
!        * dirIntAndRealBuffer stores the information in                *
!        * mySubfaces in the correct place in the buffers intSend,      *
!        * realSend and indicesSend. In this routine only info of the   *
!        * direct halo cells is stored.                                 *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nMySubfaces
         integer(kind=intType), intent(in) :: signSubface

         type(localSubfaceType), dimension(nMySubfaces), &
                                               intent(in) :: mySubfaces
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the subfaces and the quadrilateral faces of these
         ! subfaces. First store the info of the halo's that must be
         ! interpolated.

         do nn=1,nMySubfaces
           do mm=1,mySubfaces(nn)%nQuad

             ! Check if the halo('s) should be interpolated.

             testInterpol: if( mySubfaces(nn)%searchQuad(mm) ) then

               ! The 1st halo must always be constructed. Determine the
               ! processor to which the info must be sent.

               jj = mySubfaces(nn)%donorProc(mm)

               ! Store the halo location. Update the corresponding entry
               ! in the counter countIndices.

               i = countIndices(jj)
               countIndices(jj) = countIndices(jj) + 4

               indicesSend(i+1) = mySubfaces(nn)%blockID
               indicesSend(i+2) = mySubfaces(nn)%indHalo1(mm,1)
               indicesSend(i+3) = mySubfaces(nn)%indHalo1(mm,2)
               indicesSend(i+4) = mySubfaces(nn)%indHalo1(mm,3)

               ! Store the information for intSend and update the
               ! corresponding entry in the counter countInt.
               ! Note that the donorDualQ is multiplied by signSubface,
               ! such that a distinction can be made between part 1 and
               ! part 2 of the interface.

               i  = countInt(jj)
               countInt(jj) = countInt(jj) + 2

               intSend(i+1) = signSubface &
                            * mySubfaces(nn)%donorDualQ(mm)
               intSend(i+2) = 1

               ! Store the information for realSend, i.e. the
               ! interpolation weights u and v. Update the corresponding
               ! entry in the counter countReal.

               i = countReal(jj)
               countReal(jj) = countReal(jj) + 2

               realSend(i+1) = mySubfaces(nn)%u(mm)
               realSend(i+2) = mySubfaces(nn)%v(mm)

               ! Test if the second halo cell must be interpolated for
               ! this face. This is indicated by a nonnegative value of
               ! the first index of indHalo2.

               test2nd: if(mySubfaces(nn)%indHalo2(mm,1) > -1) then

                 ! The 2nd level halo must also be interpolated. Store
                 ! the same info as for the 1st level halo.

                 ! Store the information for the local halo location.

                 i = countIndices(jj)
                 countIndices(jj) = countIndices(jj) + 4

                 indicesSend(i+1) = mySubfaces(nn)%blockID
                 indicesSend(i+2) = mySubfaces(nn)%indHalo2(mm,1)
                 indicesSend(i+3) = mySubfaces(nn)%indHalo2(mm,2)
                 indicesSend(i+4) = mySubfaces(nn)%indHalo2(mm,3)

                 ! The info stored in the buffer intSend.

                 i  = countInt(jj)
                 countInt(jj) = countInt(jj) + 2

                 intSend(i+1) = signSubface &
                              * mySubfaces(nn)%donorDualQ(mm)
                 intSend(i+2) = 2

                 ! The info stored in realSend.

                 i = countReal(jj)
                 countReal(jj) = countReal(jj) + 2

                 realSend(i+1) = mySubfaces(nn)%u(mm)
                 realSend(i+2) = mySubfaces(nn)%v(mm)

               endif test2nd

             endif testInterpol
           enddo
         enddo

         end subroutine dirIntAndRealBuffer

       end subroutine updateSlidingCommPattern

!
!      ******************************************************************
!      *                                                                *
!      * File:          check1to1Subfaces.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-25-2004                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine check1to1Subfaces
!
!      ******************************************************************
!      *                                                                *
!      * check1to1Subfaces checks if the 1 to 1 internal subfaces,      *
!      * including the periodic ones, match up to a certain tolerance.  *
!      * If not, a warning will be printed. The computation is not      *
!      * terminated, because sometimes gaps are introduced on purpose,  *
!      * e.g. near a wing tip in an H-topology in spanwise direction.   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use communication
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer :: ierr, procId, size

       integer, dimension(mpi_status_size) :: status
       integer, dimension(nProc)           :: sizeMessage
       integer, dimension(nProc)           :: recvcounts, displs

       integer(kind=intType) :: i, j, k, ll1, ll2, ll3, ic, jc, kc
       integer(kind=intType) :: sps, nn, mm, ll, ii, proc, nFCheck
       integer(kind=intType) :: stepI, stepJ, stepK
       integer(kind=intType) :: nMessagesSend, nMessagesReceive
       integer(kind=intType) :: nBad, nBadGlobal

       integer(kind=intType), dimension(3,3) :: tMat

       integer(kind=intType), dimension(0:nProc) :: nFSend, nCoor
       integer(kind=intType), dimension(nProc)   :: nFCount, nCCount

       integer(kind=intType), dimension(:,:), allocatable :: intBuf
       integer(kind=intType), dimension(:,:), allocatable :: intRecv
       integer(kind=intType), dimension(:,:), allocatable :: badSubfaces
       integer(kind=intType), dimension(:,:), allocatable :: badGlobal

       real(kind=realType), dimension(:),   allocatable :: badDist
       real(kind=realType), dimension(:),   allocatable :: badDistGlobal
       real(kind=realType), dimension(:,:), allocatable :: realBuf
       real(kind=realType), dimension(:,:), allocatable :: realRecv

       character(len=7) :: intString
!
!      Function definitions.
!
       integer(kind=intType) :: delta
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local number of faces that must be sent to other *
!      * processors, including to myself. Also determine the number of  *
!      * faces I have to check.                                         *
!      *                                                                *
!      ******************************************************************
!
       nFCheck = 0
       nFSend  = 0
       nCoor   = 0

       ! Loop over the local blocks and its subfaces of the number of
       ! spectral solutions to be checked.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           call setPointers(nn, 1_intType, sps)

           ! Loop over the 1 to 1 subfaces.

           do mm=1,n1to1

             ! Add the offset of nBocos to mm, such that the entries
             ! in the arrays corresponds to this 1 to 1 subface.

             ll = mm + nBocos

             ! Update nFSend and nCoor of the correct processor and
             ! update nFCheck.

             ii = (abs(inEnd(ll) - inBeg(ll)) + 1) &
                * (abs(jnEnd(ll) - jnBeg(ll)) + 1) &
                * (abs(knEnd(ll) - knBeg(ll)) + 1)
             ll = neighProc(ll) + 1

             nFSend(ll) = nFSend(ll) + 1
             nCoor(ll)  = nCoor(ll)  + ii

             nFCheck = nFCheck + 1

           enddo
         enddo
       enddo

       ! Put nFSend and nCoor in cumulative storage format. Note
       ! that nFSend and nCoor start at index 0. Store the starting
       ! value in nFCount and nCCount respectively.

       do nn=1,nProc
         nFCount(nn) = nFSend(nn-1)
         nCCount(nn) = nCoor(nn-1)

         nFSend(nn) = nFSend(nn) + nFSend(nn-1)
         nCoor(nn)  = nCoor(nn)  + nCoor(nn-1)
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the integer and real buffers to store the subface    *
!      * information to be communicated.                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the integer and real buffers to store
       ! the information of the subfaces to be communicated.

       nn = nFSend(nProc)
       mm = nCoor(nProc)

       allocate(intBuf(10,nn), realBuf(3,mm), stat=ierr)
       if(ierr /= 0) &
         call terminate("check1to1Subfaces", &
                        "Memory allocation failure for intBuf and &
                        &realBuf")

       ! Repeat the loop over the number of local blocks and its subfaces
       ! of the number of spectral solutions to be checked.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           call setPointers(nn, 1_intType, sps)

           ! Loop over the 1 to 1 subfaces.

           do mm=1,n1to1

             ! Add the offset of nBocos to mm, such that the entries
             ! in the arrays corresponds to this 1 to 1 subface.

             ll = mm + nBocos

             ! Store the donor range, the local block number, the
             ! spectral solution and global block ID in intBuf.

             proc          = neighProc(ll) + 1
             nFCount(proc) = nFCount(proc) + 1
             ii            = nFCount(proc)

             intBuf( 1,ii) = dinBeg(ll)
             intBuf( 2,ii) = dinEnd(ll)
             intBuf( 3,ii) = djnBeg(ll)
             intBuf( 4,ii) = djnEnd(ll)
             intBuf( 5,ii) = dknBeg(ll)
             intBuf( 6,ii) = dknEnd(ll)

             intBuf( 7,ii) = neighBlock(ll)
             intBuf( 8,ii) = sps

             intBuf( 9,ii) = nbkGlobal
             intBuf(10,ii) = cgnsSubface(ll)

             ! Determine whether the subface has positive or negative
             ! running indices. To be sure that the correct sequence is
             ! stored in realBuf, the loop must be performed over the
             ! donor range (i, j, and k indices could be swapped and
             ! therefore stored wrongly in the 1D buffer).

             stepI = 1; if(dinEnd(ll) < dinBeg(ll)) stepI = -1
             stepJ = 1; if(djnEnd(ll) < djnBeg(ll)) stepJ = -1
             stepK = 1; if(dknEnd(ll) < dknBeg(ll)) stepK = -1

             ! Determine the transformation matrix between the donor
             ! and the current face. As the information stored in l1,
             ! l2 and l3 is for the transformation matrix to the donor
             ! face, the transpose must be taken.

             ll1 = l1(ll); ll2 = L2(ll); ll3 = l3(ll)

             tMat(1,1) = sign(1_intType,ll1) * delta(ll1,1_intType)
             tMat(1,2) = sign(1_intType,ll1) * delta(ll1,2_intType)
             tMat(1,3) = sign(1_intType,ll1) * delta(ll1,3_intType)

             tMat(2,1) = sign(1_intType,ll2) * delta(ll2,1_intType)
             tMat(2,2) = sign(1_intType,ll2) * delta(ll2,2_intType)
             tMat(2,3) = sign(1_intType,ll2) * delta(ll2,3_intType)

             tMat(3,1) = sign(1_intType,ll3) * delta(ll3,1_intType)
             tMat(3,2) = sign(1_intType,ll3) * delta(ll3,2_intType)
             tMat(3,3) = sign(1_intType,ll3) * delta(ll3,3_intType)

             ! Store the coordinates in realBuf by looping over the
             ! points of the subface.

             ii = nCCount(proc)

             do k=dknBeg(ll), dknEnd(ll), stepK
               do j=djnBeg(ll), djnEnd(ll), stepJ
                 do i=dinBeg(ll), dinEnd(ll), stepI

                   ! Determine the nodal indices in the current block.

                   ll1 = i - dinBeg(ll)
                   ll2 = j - djnBeg(ll)
                   ll3 = k - dknBeg(ll)

                   ic = inBeg(ll) &
                      + tMat(1,1)*ll1 + tMat(1,2)*ll2 + tMat(1,3)*ll3
                   jc = jnBeg(ll) &
                      + tMat(2,1)*ll1 + tMat(2,2)*ll2 + tMat(2,3)*ll3
                   kc = knBeg(ll) &
                      + tMat(3,1)*ll1 + tMat(3,2)*ll2 + tMat(3,3)*ll3

                   ! Store the coordinates in the buffer.

                   ii = ii + 1
                   realBuf(1,ii) = x(ic,jc,kc,1)
                   realBuf(2,ii) = x(ic,jc,kc,2)
                   realBuf(3,ii) = x(ic,jc,kc,3)
                 enddo
               enddo
             enddo

             ! Set ii to the number of points stored for this subface.
             ! This number is needed when for the periodic correction.

             ii = ii - nCCount(proc)

             ! Check if a periodic correction must be applied and if so
             ! call the routine to do so for the coordinates of this
             ! subface. Note that internally created subfaces due to
             ! block splitting must be excluded from this test.

             k = cgnsSubface(ll)
             if(k > 0) then
               if( cgnsDoms(nbkGlobal)%conn1to1(k)%periodic ) then

                 j = nCCount(proc) + 1
                 call periodicTransformSubface(realBuf(1,j), ii,     &
                     cgnsDoms(nbkGlobal)%conn1to1(k)%rotationCenter, &
                     cgnsDoms(nbkGlobal)%conn1to1(k)%rotationAngles, &
                     cgnsDoms(nbkGlobal)%conn1to1(k)%translation)
               endif
             endif

             ! Update the counter nCCount(proc).

             nCCount(proc) = nCCount(proc) + ii

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of messages I will send and receive.      *
!      * The term message in this routine means a the complete subface  *
!      * info, i.e. a combination of an integer and real message.       *
!      *                                                                *
!      ******************************************************************
!
       ! Fill the array nCCount with 0's and 1's. A 0 indicates that no
       ! message is sent to the corresponding processor.

       do nn=1,nProc
         nCCount(nn) = 0
         if(nFSend(nn) > nFSend(nn-1)) nCCount(nn) = 1
       enddo

       ! Make sure that no message is sent to myself. An offset of + 1
       ! must be added, because the processor numbers start at 0 and
       ! nCCount at 1.

       nCCount(myID+1) = 0

       ! Determine the number of message I have to receive.

       sizeMessage = 1
       call mpi_reduce_scatter(nCCount, nMessagesReceive,          &
                               sizeMessage, sumb_integer, mpi_sum, &
                               SUmb_comm_world, ierr)

       ! Send the data I have to send. Do not send a message to myself.
       ! That is handled separately. As nonblocking sends must be used
       ! to avoid deadlock and two messages are sent to a processor,
       ! the sendRequests (stored in the module communication) are used
       ! for the integer messages and the recvRequests (same module)
       ! are used for the real messages.

       ii = 0
       do nn=1,nProc

         ! Check if something must be sent to this processor.

         if(nCCount(nn) == 1) then

           ! Send the integer buffer.

           ii     = ii + 1
           procID = nn - 1
           size   = 10*(nFSend(nn) - nFSend(nn-1))
           mm     = nFSend(nn-1) + 1

           call mpi_isend(intBuf(1,mm), size, sumb_integer, procID,  &
                          procID, SUmb_comm_world, sendRequests(ii), &
                          ierr)

           ! Send the real buffer.

           size = 3*(nCoor(nn) - nCoor(nn-1))
           mm   = nCoor(nn-1) + 1

           call mpi_isend(realBuf(1,mm), size, sumb_real, procID, &
                          procID+1, SUmb_comm_world,              &
                          recvRequests(ii), ierr)

         endif
       enddo

       ! Store the number of messages sent.

       nMessagesSend = ii
!
!      ******************************************************************
!      *                                                                *
!      * Check the coordinates of the subfaces which should have been   *
!      * sent to myself.                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize nBad to 0 and allocate the memory to store possible
       ! bad subfaces.

       nBad = 0
       allocate(badSubfaces(4,nFCheck), badDist(nFCheck), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("check1to1Subfaces", &
                        "Memory allocation failure for badSubfaces &
                        &and badDist")

       ! Determine the number of local subfaces that must be checked.
       ! If any present, check them.

       nn = nFSend(myID+1) - nFSend(myID)
       if(nn > 0) then
         ii = nFSend(myID) + 1
         mm = nCoor(myID)  + 1
         call checkSubfaceCoor(intBuf(1,ii), realBuf(1,mm), nn, &
                               nBad, badSubfaces, badDist,      &
                               nTimeIntervalsSpectral)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Check the coordinates of the subfaces which are received from  *
!      * other processors.                                              *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of messages to be received.

       do ii=1,nMessagesReceive

         ! Wait until an integer message arrives and determine the
         ! source and size of the message.

         call mpi_probe(mpi_any_source, myID, SUmb_comm_world, &
                       status, ierr)

         procID = status(mpi_source)
         call mpi_get_count(status, sumb_integer, size, ierr)

         ! Check in debug mode that the incoming message is of
         ! correct size.

         if( debug ) then
           if(size == mpi_undefined .or. mod(size,10) /= 0) &
             call terminate("check1to1Subfaces",            &
                            "Unexpected size of integer message")
         endif

         ! Determine the number of subfaces this message contains and
         ! allocate the memory for the integer receive buffer.

         nn = size/10
         allocate(intRecv(10,nn), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("check1to1Subfaces", &
                          "Memory allocation failure for intRecv")

         ! Receive the integer buffer. Blocking receives can be used,
         ! because the message has already arrived.

         call mpi_recv(intRecv, size, sumb_integer, procID, &
                       myID, SUmb_comm_world, status, ierr)

         ! Probe for the corresponding real buffer and determine its
         ! size.

         call mpi_probe(procID, myID+1, SUmb_comm_world, &
                       status, ierr)
         call mpi_get_count(status, sumb_real, size, ierr)

         ! Check in debug mode that the incoming message is of
         ! correct size.

         if( debug ) then
           if(size == mpi_undefined .or. mod(size,3) /= 0) &
             call terminate("check1to1Subfaces",           &
                            "Unexpected size of real message")
         endif

         ! Determine the total number of coordinates in the message and
         ! allocate the memory for the real receive buffer.

         mm = size/3
         allocate(realRecv(3,mm), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("check1to1Subfaces", &
                          "Memory allocation failure for realRecv")

         ! Receive the real buffer. Blocking receives can be used,
         ! because the message has already arrived.

         call mpi_recv(realRecv, size, sumb_real, procID, &
                       myID+1, SUmb_comm_world, status, ierr)

         ! Check the subfaces stored in these messages.

         call checkSubfaceCoor(intRecv, realRecv, nn, nBad, & 
                               badSubfaces, badDist,        &
                               nTimeIntervalsSpectral)

         ! Release the memory of the integer and real receive buffer.

         deallocate(intRecv, realRecv, stat=ierr)
         if(ierr /= 0)                         &
           call terminate("check1to1Subfaces", &
                          "Deallocation failure for intRecv &
                          &and realRecv")
       enddo

       ! Complete the nonblocking sends.

       size = nMessagesSend
       do nn=1,nMessagesSend
         call mpi_waitany(size, sendRequests, procID, status, ierr)
         call mpi_waitany(size, recvRequests, procID, status, ierr)
       enddo

       ! Deallocate the memory for the integer and real buffers.

       deallocate(intBuf, realBuf, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("check1to1Subfaces", &
                        "Deallocation failure for intBuf and realBuf")
!
!      ******************************************************************
!      *                                                                *
!      * Determine the global number of bad subfaces and gather this    *
!      * information.                                                   *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the global number of bad subfaces.

       call mpi_allgather(nBad, 1, sumb_integer, nCCount, 1, &
                          sumb_integer, SUmb_comm_world, ierr)

       ! Determine the global number of bad subfaces and the arrays
       ! recvcounts and displs needed for the call to allgatherv

       nBadGlobal    = nCCount(1)
       recvcounts(1) = nCCount(1)
       displs(1)     = 0

       do nn=2,nProc
         nBadGlobal     = nBadGlobal + nCCount(nn)
         recvcounts(nn) = nCCount(nn)
         displs(nn)     = displs(nn-1) + recvcounts(nn-1)
       enddo

       ! Allocate the memory to store the global bad surfaces.

       allocate(badGlobal(4,nBadGlobal), badDistGlobal(nBadGlobal), &
                stat=ierr)
       if(ierr /= 0)                         &
         call terminate("check1to1Subfaces", &
                        "Memory allocation failure for badGlobal &
                        &and badDistGlobal")

       ! Gather the data. First the distance info.

       size = nBad
       call mpi_allgatherv(badDist, size, sumb_real, badDistGlobal, &
                           recvcounts, displs, sumb_real,           &
                           SUmb_comm_world, ierr)

       ! And the integer info. Multiply recvcounts and displs
       ! by 4, because 4 integers are received.

       do nn=1,nProc
         recvcounts(nn) = 4*recvcounts(nn)
         displs(nn)     = 4*displs(nn)
       enddo

       size = 4*nBad
       call mpi_allgatherv(badSubfaces, size, sumb_integer, badGlobal, &
                           recvcounts, displs, sumb_integer,           &
                           SUmb_comm_world, ierr)

       ! Sort the bad subfaces and get rid of the multiple entries.

       call sortBadEntities(nBadGlobal, badGlobal, badDistGlobal, .true.)

       ! Check for the presence of any internally created subfaces.
       ! This only occurs when something goes wrong in the block
       ! splitting and therefore the program is terminated.

       do nn=1,nBadGlobal
         if(badGlobal(2,nn) == 0) then
           if(myID == 0)                         &
             call terminate("check1to1Subfaces", &
                            "Non-matching internally created &
                            &face found.")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       enddo

       ! Print the bad subfaces, if present. Only processor 0 performs
       ! this task.

       if(myID == 0 .and. nBadGlobal > 0) then

         write(intString,"(i6)") nBadGlobal
         intString = adjustl(intString)

         print "(a)",  "#"
         print 101
         print 102, trim(intString)
         print 103
         print 104
         print "(a)", "#"

         do nn=1,nBadGlobal
           i = badGlobal(1,nn)
           j = badGlobal(2,nn)

           write(intString,"(i6)") badGlobal(3,nn)
           intString = adjustl(intString)

           ! Write a different error message if more than one grid has
           ! been read.

           if(nTimeIntervalsSpectral > 1) then

             if(badGlobal(4,nn) == 1) then
               print 105, trim(intString),                           &
                          trim(cgnsDoms(i)%zoneName),                &
                          trim(cgnsDoms(i)%conn1to1(j)%connectName), &
                          trim(cgnsDoms(i)%conn1to1(j)%donorName),   &
                          badDistGlobal(nn)
             else
               print 106, trim(intString),                           &
                          trim(cgnsDoms(i)%zoneName),                &
                          trim(cgnsDoms(i)%conn1to1(j)%connectName), &
                          trim(cgnsDoms(i)%conn1to1(j)%donorName),   &
                          badDistGlobal(nn)
             endif

           else

             if(badGlobal(4,nn) == 1) then
               print 107, trim(cgnsDoms(i)%zoneName),                &
                          trim(cgnsDoms(i)%conn1to1(j)%connectName), &
                          trim(cgnsDoms(i)%conn1to1(j)%donorName),   &
                          badDistGlobal(nn)
             else
               print 108, trim(cgnsDoms(i)%zoneName),                &
                          trim(cgnsDoms(i)%conn1to1(j)%connectName), &
                          trim(cgnsDoms(i)%conn1to1(j)%donorName),   &
                          badDistGlobal(nn)
             endif

           endif
         enddo

         print "(a)", "#"

 101     format("#                      Warning")
 102     format("# Found ",a," one to one subfaces which do not &
                &coincide.")
 103     format("# Computation continues, but be aware of it.")
 104     format("# List of nonmatching one to one subfaces.")
 105     format("# Spectral grid ",a, ", zone ",a, ", periodic &
                &subface ", a, " does not match donor ", a,". &
                &Maximum deviation: ", e12.5)
 106     format("# Spectral grid ",a, ", zone ",a, ", subface ", &
                a, " does not match donor ", a,". &
                &Maximum deviation: ", e12.5)
 107     format("# Zone ",a, ", periodic subface ", a, " does not &
                &match donor ", a,". Maximum deviation: ", e12.5)

 108     format("# Zone ",a, ", subface ", a, " does not match &
                &donor ", a,". Maximum deviation: ", e12.5)

       endif

       ! Deallocate the memory of to store the bad subfaces.

       deallocate(badSubfaces, badGlobal, badDist, badDistGlobal, &
                  stat=ierr)
       if(ierr /= 0)                         &
         call terminate("check1to1Subfaces", &
                        "Deallocation failure for badSubfaces, &
                        &badGlobal, badDist and badDistGlobal")

       ! Synchronize the processors, because wild cards have been used.

       call mpi_barrier(SUmb_comm_world, ierr)

       end subroutine check1to1Subfaces

!      ==================================================================

       subroutine periodicTransformSubface(coor, nn, rotCenter, &
                                           rotAngles, translation)
!
!      ******************************************************************
!      *                                                                *
!      * periodicTransformSubface transforms the given set of           *
!      * coordinates using the periodic transformation defined by       *
!      * rotCenter, rotAngles and translation.                          *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn

       real(kind=realType), dimension(3), intent(in) :: rotCenter
       real(kind=realType), dimension(3), intent(in) :: rotAngles
       real(kind=realType), dimension(3), intent(in) :: translation

       real(kind=realType), dimension(3,nn), intent(inout) :: coor
!
!      Local variables.
!
       integer(kind=intType) :: i

       real(kind=realType) :: cosTheta, cosPhi, cosPsi
       real(kind=realType) :: sinTheta, sinPhi, sinPsi
       real(kind=realType) :: dx, dy, dz

       real(kind=realType), dimension(3)   :: trans
       real(kind=realType), dimension(3,3) :: rotMatrix
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Construct from the given rotation angles the rotation matrix
       ! from the current coordinates to the donor coordinates.
       ! Note that the sequence of rotation is first rotation around the
       ! x-axis, followed by rotation around the y-axis and finally
       ! rotation around the z-axis.

       cosTheta = cos(rotAngles(1)); sinTheta = sin(rotAngles(1))
       cosPhi   = cos(rotAngles(2)); sinPhi   = sin(rotAngles(2))
       cosPsi   = cos(rotAngles(3)); sinPsi   = sin(rotAngles(3))

       rotMatrix(1,1) =  cosPhi*cosPsi
       rotMatrix(2,1) =  cosPhi*sinPsi
       rotMatrix(3,1) = -sinPhi

       rotMatrix(1,2) = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi
       rotMatrix(2,2) = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi
       rotMatrix(3,2) = sinTheta*cosPhi

       rotMatrix(1,3) = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi
       rotMatrix(2,3) = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi
       rotMatrix(3,3) = cosTheta*cosPhi

       ! Store the translation plus the rotation center in trans.

       trans(1) = rotCenter(1) + translation(1)
       trans(2) = rotCenter(2) + translation(2)
       trans(3) = rotCenter(3) + translation(3)

       ! Loop over the number of coordinates to be corrected.

       do i=1,nn

         ! Determine the relative position w.R.T. The rotation center.

         dx = coor(1,i) - rotCenter(1)
         dy = coor(2,i) - rotCenter(2)
         dz = coor(3,i) - rotCenter(3)

         ! Determine the coordinates after the transformation.

         coor(1,i) = rotMatrix(1,1)*dx + rotMatrix(1,2)*dy &
                   + rotMatrix(1,3)*dz + trans(1)
         coor(2,i) = rotMatrix(2,1)*dx + rotMatrix(2,2)*dy &
                   + rotMatrix(2,3)*dz + trans(2)
         coor(3,i) = rotMatrix(3,1)*dx + rotMatrix(3,2)*dy &
                   + rotMatrix(3,3)*dz + trans(3)

       enddo

       end subroutine periodicTransformSubface

!      ==================================================================

       subroutine checkSubfaceCoor(subfaceInfo, coor, nFace,   &
                                   nBad, badSubfaces, badDist, &
                                   nSpectral)
!
!      ******************************************************************
!      *                                                                *
!      * checkSubfaceCoor checks if the coordinates of the subfaces     *
!      * defined in subfaceInfo and coor match the coordinates stored   *
!      * in flowDoms.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: nFace, nSpectral
       integer(kind=intType), intent(inout) :: nBad

       integer(kind=intType), dimension(10,*), intent(in) :: subfaceInfo
       integer(kind=intType), dimension( 4,*), intent(inout) :: badSubfaces

       real(kind=realType), dimension(*),   intent(inout) :: badDist
       real(kind=realType), dimension(3,*), intent(in)    :: coor
!
!      Local parameter.
!
       real(kind=realType), parameter :: relTol = 0.05_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, mm, nn
       integer(kind=intType) :: stepI, stepJ, stepK
       integer(kind=intType) :: im1, ip1, jm1, jp1, km1, kp1
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
       integer(kind=intType) :: blockID, sps, globalDonID, subfaceID

       real(kind=realType) :: dist2, dist2I, dist2J, dist2K, tol
       real(kind=realType) :: factI, factJ, factK, diffMax

       character(len=maxStringLen) :: errorMessage
       character(len=7)            :: intString

       logical :: badFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the counter ii for the coordinates and loop over the
       ! number of subfaces to be checked.

       ii = 0
       subfaceLoop: do nn=1,nFace

         ! Store the integer info a bit easier.

         iBeg = subfaceInfo(1,nn)
         iEnd = subfaceInfo(2,nn)
         jBeg = subfaceInfo(3,nn)
         jEnd = subfaceInfo(4,nn)
         kBeg = subfaceInfo(5,nn)
         kEnd = subfaceInfo(6,nn)

         blockID     = subfaceInfo( 7,nn)
         sps         = subfaceInfo( 8,nn)
         globalDonID = subfaceInfo( 9,nn)
         subfaceID   = subfaceInfo(10,nn)

         ! Set the pointers to the block to be searched.

         call setPointers(blockID, 1_intType, sps)

         ! Find the matching 1 to 1 subface.

         do mm=1,n1to1

           ! Add the offset of nBocos to nn such that it contains the
           ! correct index of the arrays.

           jj = mm + nBocos

           ! If this is the subface exit the loop.

           if(min(iBeg,iEnd) == min(inBeg(jj),inEnd(jj)) .and. &
              max(iBeg,iEnd) == max(inBeg(jj),inEnd(jj)) .and. &
              min(jBeg,jEnd) == min(jnBeg(jj),jnEnd(jj)) .and. &
              max(jBeg,jEnd) == max(jnBeg(jj),jnEnd(jj)) .and. &
              min(kBeg,kEnd) == min(knBeg(jj),knEnd(jj)) .and. &
              max(kBeg,kEnd) == max(knBeg(jj),knEnd(jj))) exit

         enddo

         ! Check if a subface was found. If not terminate.

         if(mm > n1to1) then

           if(nSpectral > 1) then
             write(intString,"(i6)") sps
             intString = adjustl(intString)
             write(errorMessage,301) trim(intString), &
               trim(cgnsDoms(globalDonID)%zoneName),  &
               trim(cgnsDoms(globalDonID)%conn1to1(subfaceID)%connectName)
 301         format("Spectral grid ",a,", Zone ",a, "Connectivity ", a, &
                    ": 1 to 1 subface not found. Something is &
                    &seriously wrong with the zone connectivity.")
           else
             write(errorMessage,302) &
               trim(cgnsDoms(globalDonID)%zoneName), &
               trim(cgnsDoms(globalDonID)%conn1to1(subfaceID)%connectName)
 302         format("Zone ",a, "Connectivity", a, &
                    ": 1 to 1 subface not found. Something is &
                    &seriously wrong with the zone connectivity.")
           endif

           call terminate("checkSubfaceCoor", errorMessage)

         endif

         ! Determine whether positive or negative running indices must
         ! be used for the subface. This depends on the sequence stored
         ! in the coordinate buffer and not of the sequence of the
         ! 1 to 1 subface mm.

         stepI = 1; if(iEnd < iBeg) stepI = -1
         stepJ = 1; if(jEnd < jBeg) stepJ = -1
         stepK = 1; if(kEnd < kBeg) stepK = -1

         ! Initialize badFace to .false. to indicate that it is
         ! a correct subface. Set the maximum difference to zero.

         badFace = .false.
         diffMax = zero

         ! Loop over the coordinates of the subface to see if they
         ! match of to a certain tolerance.

         do k=kBeg,kEnd,stepK

           ! Determine the indices to the left and to the right.
           ! Do not exceed the block boundary. Determine the
           ! scaling factor of the tolerance accordingly.

           km1 = max(1_intType,k-1_intType)
           kp1 = min(kl,       k+1_intType)

           factK = one
           if(km1 == k) factK = factK*four
           if(kp1 == k) factK = factK*four

           ! Loop in j-direction.

           do j=jBeg,jEnd,stepJ

             ! Determine the neighbors to the left and right and factJ.

             jm1 = max(1_intType,j-1_intType)
             jp1 = min(jl,       j+1_intType)

             factJ = one
             if(jm1 == j) factJ = factJ*four
             if(jp1 == j) factJ = factJ*four

             ! Loop in i-direction.

             do i=iBeg,iEnd,stepI

               ! Determine the neighbors to the left and right and factI.

               im1 = max(1_intType,i-1_intType)
               ip1 = min(il,       i+1_intType)

               factI = one
               if(im1 == i) factI = factI*four
               if(ip1 == i) factI = factI*four

               ! Determine the distances squared in i, j and k direction.
               ! The reason for squared is that some square roots are
               ! avoided this way.

               dist2I = (x(ip1,j,k,1) - x(im1,j,k,1))**2 &
                      + (x(ip1,j,k,2) - x(im1,j,k,2))**2 &
                      + (x(ip1,j,k,3) - x(im1,j,k,3))**2

               dist2J = (x(i,jp1,k,1) - x(i,jm1,k,1))**2 &
                      + (x(i,jp1,k,2) - x(i,jm1,k,2))**2 &
                      + (x(i,jp1,k,3) - x(i,jm1,k,3))**2

               dist2K = (x(i,j,kp1,1) - x(i,j,km1,1))**2 &
                      + (x(i,j,kp1,2) - x(i,j,km1,2))**2 &
                      + (x(i,j,kp1,3) - x(i,j,km1,3))**2

               ! Make sure that singular lines are excluded.

               if(dist2I < eps) dist2I = large
               if(dist2J < eps) dist2J = large
               if(dist2K < eps) dist2K = large

               ! Determine the tolerance for identical points.
               ! Note the multiplication with the square of the relative
               ! tolerance, because the distances are squared as well.

               tol = factK*dist2K
               tol = min(tol, factJ*dist2J)
               tol = min(tol, factI*dist2I)

               tol = tol*relTol*relTol

               ! Update the counter for the coordinate in the buffer
               ! and determine the distance squared between the points
               ! that should be identical.

               ii = ii + 1
               dist2 = (x(i,j,k,1) - coor(1,ii))**2 &
                     + (x(i,j,k,2) - coor(2,ii))**2 &
                     + (x(i,j,k,3) - coor(3,ii))**2

               ! Flag the subface to bad if the nodes do not coincide
               ! within the given tolerance. Store the difference if
               ! it is larger than the currently stored value.

               if(dist2 > tol) then
                 badFace = .true.
                 diffMax = max(diffMax, dist2)
               endif

             enddo
           enddo
         enddo

         ! Store this subface if it is not matching. The following data
         ! is stored: global block number of the donor, the subface ID
         ! on this block, the spectral solution, whether or not this
         ! is a periodic face and the maximum difference.

         if( badFace ) then
           nBad = nBad + 1

           badSubfaces(1,nBad) = globalDonID
           badSubfaces(2,nBad) = subfaceID
           badSubfaces(3,nBad) = sps

           badSubfaces(4,nBad) = 0
           if(subfaceID > 0) then
             if( cgnsDoms(globalDonID)%conn1to1(subfaceID)%periodic ) &
               badSubfaces(4,nBad) = 1
           endif

           badDist(nBad) = sqrt(diffMax)

        !  if(subfaceID == 0 .and. myID == 1) then
        !    write(*,*) "myID: ", myID, badDist(nBad)
        !    write(*,*) blockID, mm, globalDonID, nbkGlobal
        !    write(*,"(6I4)") iBegor, iEndor, jBegor, jEndor, kBegor, kEndor

        !    jj = mm + nBocos
        !    write(*,"(6I4)") inBeg(jj),inEnd(jj), jnBeg(jj),jnEnd(jj), knBeg(jj),knEnd(jj)
        !  endif
         endif

       end do subfaceLoop

       end subroutine checkSubfaceCoor

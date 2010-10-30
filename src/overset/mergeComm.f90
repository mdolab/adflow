!
!      ******************************************************************
!      *                                                                *
!      * File:          mergeComm.f90                                   *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-25-2005                                      *
!      * Last modified: 08-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mergeComm(commPattern1, internal1, &
                            commPattern2, internal2)

!      ******************************************************************
!      *                                                                *
!      * MergeComm combines the 2 communication patterns into a single  *
!      * pattern. In particular, the 2nd internal and external patterns *
!      * are merged into the 1st one.                                   *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       type(commType), intent(inout) :: commPattern1
       type(commType), intent(in)    :: commPattern2

       type(internalCommType), intent(inout) :: internal1
       type(internalCommType), intent(in)    :: internal2
!
!      Local variables.

       integer               :: ierr
       integer(kind=intType) :: nn, i, k, j1, j2, nInterp

       logical :: freeMem

       type(sendCommListType), dimension(:), pointer :: tmpSendList
       type(recvCommListType), dimension(:), pointer :: tmpRecvList
!
!      Interfaces
!
       interface
         subroutine reallocateInteger(intArray, newSize, oldSize, &
                                       alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize, oldSize
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger

         !===============================================================

         subroutine reallocateInteger2(intArray,            &
                                        newSize1, newSize2, &
                                        oldSize1, oldSize2, &
                                        alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:,:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger2

         !===============================================================

         subroutine reallocateReal2(realArray,           &
                                     newSize1, newSize2, &
                                     oldSize1, oldSize2, &
                                     alwaysFreeMem)
           use precision
           implicit none

           real(kind=realType), dimension(:,:), pointer :: realArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateReal2
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the number of interpolants in this communication pattern.
       ! Note this is the same for all so just use the internal part.

       nInterp = ubound(internal1%donorInterp,2)

       ! Reallocate the internal communication lists.

       nn              = internal1%ncopy
       internal1%ncopy = internal1%ncopy + internal2%ncopy

       call reallocateInteger(internal1%haloBlock, &
                              internal1%ncopy, nn, .true.)
       call reallocateInteger(internal1%donorBlock, &
                              internal1%ncopy, nn, .true.)
       call reallocateInteger2(internal1%haloIndices, &
                               internal1%ncopy, 3_intType, &
                               nn,              3_intType, .true.)
       call reallocateInteger2(internal1%donorIndices, &
                               internal1%ncopy, 3_intType, &
                               nn,              3_intType, .true.)
       call reallocateReal2(internal1%donorInterp, &
                            internal1%ncopy, nInterp, &
                            nn,              nInterp, .true.)

       ! Tack 2nd onto end of 1st internal.

       do j2 = 1,internal2%ncopy
         j1 = nn + j2
         internal1%haloBlock(j1)      = internal2%haloBlock(j2)
         internal1%haloIndices(j1,:)  = internal2%haloIndices(j2,:)
         internal1%donorBlock(j1)     = internal2%donorBlock(j2)
         internal1%donorIndices(j1,:) = internal2%donorIndices(j2,:)
         internal1%donorInterp(j1,:)  = internal2%donorInterp(j2,:)
       end do

!      ==================================================================

       ! Loop over the sending lists for the 2nd.

       sends: do i = 1,commPattern2%nProcSend

         ! Find the index for this processor in the 1st.

         k = commPattern1%indexSendProc(commPattern2%sendProc(i))

         createSendList: if (k == 0) then

           ! Add another processor and list.

           k = commPattern1%nProcSend + 1

           call reallocateInteger(commPattern1%sendProc, k, &
                                  commPattern1%nProcSend, .false.)
           call reallocateInteger(commPattern1%nsend, k, &
                                  commPattern1%nProcSend, .false.)
           call reallocateInteger(commPattern1%nsendCum, k+1, &
                                  commPattern1%nProcSend+1, .false.)

           ! The reallocation of sendList. This is slightly more
           ! complicated, because it is a derived datatype.

           tmpSendList => commPattern1%sendList

           allocate(commPattern1%sendList(k), stat=ierr)
           if(ierr /= 0)                  &
             call terminate("mergeComm", &
                            "Memory allocation failure for sendList.")

           ! Set the pointers back to the original list.

           do nn = 1,commPattern1%nProcSend
             commPattern1%sendList(nn)%block   => tmpSendList(nn)%block
             commPattern1%sendList(nn)%indices => tmpSendList(nn)%indices
             commPattern1%sendList(nn)%interp  => tmpSendList(nn)%interp
           end do

           ! Set the send size to 0 and bullify.

           commPattern1%nsend(k)     = 0
           commPattern1%nsendCum(k) = commPattern1%nsendCum(k-1)

           nullify(commPattern1%sendList(k)%block,    &
                   commPattern1%sendList(k)%indices,  &
                   commPattern1%sendList(k)%interp)

           ! Release the memory of temporary list.

           deallocate(tmpSendList, stat=ierr)
           if(ierr /= 0)                  &
             call terminate("mergeComm", &
                            "Deallocation error for tmpSendList.")

           ! Modify the processor parameters for the pattern.

           commPattern1%indexSendProc(commPattern2%sendProc(i)) = k

           commPattern1%sendProc(k) = commPattern2%sendProc(i)
           commPattern1%nProcSend   = k

         end if createSendList

         ! Whether or not to free the old list memory depends on if it
         ! is currently allocated. It will be null if the list was just
         ! added above for example.

         freeMem = associated(commPattern1%sendList(k)%block)

         ! Reallocate the list size.

         nn                      = commPattern1%nsend(k)
         commPattern1%nsend(k) = nn + commPattern2%nsend(i)

         call reallocateInteger(commPattern1%sendList(k)%block, &
                                commPattern1%nsend(k), nn, freeMem)
         call reallocateInteger2(commPattern1%sendList(k)%indices, &
                                 commPattern1%nsend(k), 3_intType, &
                                 nn, 3_intType, freeMem)
         call reallocateReal2(commPattern1%sendList(k)%interp,  &
                              commPattern1%nsend(k), nInterp, &
                              nn, nInterp, freeMem)

         ! Tack 2nd onto end of 1st list for this proc.

         do j2 = 1,commPattern2%nsend(i)
           j1 = nn + j2

           commPattern1%sendList(k)%block(j1) = &
           commPattern2%sendList(i)%block(j2)

           commPattern1%sendList(k)%indices(j1,:) = &
           commPattern2%sendList(i)%indices(j2,:)

           commPattern1%sendList(k)%interp(j1,:) = &
           commPattern2%sendList(i)%interp(j2,:)
         end do

         ! Modify the cumulative sizes.

         commPattern1%nsendCum(k:) = commPattern1%nsendCum(k:) + &
                                        commPattern2%nsend(i)

       end do sends

!      ==================================================================

       ! Loop over the receiving lists for the 2nd.

       receives: do i = 1,commPattern2%nProcRecv

         ! Find the index for this processor in the 1st.

         k = commPattern1%indexRecvProc(commPattern2%recvProc(i))

         createRecvList: if (k == 0) then

           ! Add another processor and list.

           k = commPattern1%nProcRecv + 1

           call reallocateInteger(commPattern1%recvProc, k, &
                                  commPattern1%nProcRecv, .false.)
           call reallocateInteger(commPattern1%nrecv, k, &
                                  commPattern1%nProcRecv, .false.)
           call reallocateInteger(commPattern1%nrecvCum, k+1, &
                                  commPattern1%nProcRecv+1, .false.)

           ! The reallocation of recvList. This is slightly more
           ! complicated, because it is a derived datatype.

           tmpRecvList => commPattern1%recvList

           allocate(commPattern1%recvList(k), stat=ierr)
           if(ierr /= 0)                  &
             call terminate("mergeComm", &
                            "Memory allocation failure for recvList.")

           ! Set the pointers back to the original list.

           do nn = 1,commPattern1%nProcRecv

             commPattern1%recvList(nn)%block   => tmpRecvList(nn)%block
             commPattern1%recvList(nn)%indices => tmpRecvList(nn)%indices
           end do

           ! Set the recv size to 0 and bullify.

           commPattern1%nrecv(k)     = 0
           commPattern1%nrecvCum(k) = commPattern1%nrecvCum(k-1)

           nullify(commPattern1%recvList(k)%block,    &
                   commPattern1%recvList(k)%indices)

           ! Release the memory of temporary list.

           deallocate(tmpRecvList, stat=ierr)
           if(ierr /= 0)                  &
             call terminate("mergeComm", &
                            "Deallocation error for tmpRecvList.")

           ! Modify the processor parameters for the pattern.

           commPattern1%indexRecvProc(commPattern2%recvProc(i)) = k

           commPattern1%recvProc(k) = commPattern2%recvProc(i)
           commPattern1%nProcRecv   = k

         end if createRecvList

         ! Whether or not to free the old list memory depends on if it
         ! is currently allocated. It will be null if the list was just
         ! added above for example.

         freeMem = associated(commPattern1%recvList(k)%block)

         ! Reallocate the list size.

         nn                      = commPattern1%nrecv(k)
         commPattern1%nrecv(k) = nn + commPattern2%nrecv(i)

         call reallocateInteger(commPattern1%recvList(k)%block, &
                                commPattern1%nrecv(k), nn, freeMem)
         call reallocateInteger2(commPattern1%recvList(k)%indices, &
                                 commPattern1%nrecv(k), 3_intType, &
                                 nn, 3_intType, freeMem)

         ! Tack 2nd onto end of 1st list for this proc.

         do j2 = 1,commPattern2%nrecv(i)
           j1 = nn + j2

           commPattern1%recvList(k)%block(j1) = &
           commPattern2%recvList(i)%block(j2)

           commPattern1%recvList(k)%indices(j1,:) = &
           commPattern2%recvList(i)%indices(j2,:)
         end do

         ! Modify the cumulative sizes.

         commPattern1%nrecvCum(k:) = commPattern1%nrecvCum(k:) + &
                                        commPattern2%nrecv(i)

       end do receives

       end subroutine mergeComm

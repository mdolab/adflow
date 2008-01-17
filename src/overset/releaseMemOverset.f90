!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseMemOverset.f90                           *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-07-2005                                      *
!      * Last modified: 10-19-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseMemOverset(level, sps, relBlockData)
!
!      ******************************************************************
!      *                                                                *
!      * releaseMemOverset releases the memory of the overset mesh      *
!      * communication pattern for the given grid level. This must be   *
!      * done, because every new time step the communication pattern    *
!      * changes and must be recomputed.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps

       logical, intent(in) :: relBlockData
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Release the memory for the internal and external communication 
       ! patterns.

       call releaseCommPattern(commPatternOverset(level,sps))
       call releaseInternalComm(internalOverset(level,sps))

       ! Loop over the blocks if their data is to be released.

       if (relBlockData) then
         do nn=1,nDom

           deallocate(flowDoms(nn,level,sps)%ibndry,         &
                      flowDoms(nn,level,sps)%neighProcOver,  &
                      flowDoms(nn,level,sps)%neighBlockOver, &
                      flowDoms(nn,level,sps)%idonor,         &
                      flowDoms(nn,level,sps)%overint, stat=ierr)
           if (ierr /= 0) &
             call terminate("releaseMemOverset", &
                            "Deallocation failure for block data")

         end do
       end if

       end subroutine releaseMemOverset

!      ==================================================================

       subroutine releaseCommPattern(commPattern)

!      ******************************************************************
!      *                                                                *
!      * ReleaseCommPattern releases the memory of the given argument   *
!      * which should have type commPattern defined in the module       *
!      * communication.
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       type(commType), intent(inout) :: commPattern
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       do nn=1,commPattern%nProcSend

         deallocate(commPattern%sendList(nn)%block,   &
                    commPattern%sendList(nn)%indices, &
                    commPattern%sendList(nn)%interp,  &
                    stat=ierr)
         if(ierr /= 0) call terminate("releaseCommPattern", &
                                      "Deallocation error for sendList")
       enddo

       do nn=1,commPattern%nProcRecv

         deallocate(commPattern%recvList(nn)%block,   &
                    commPattern%recvList(nn)%indices, &
                    stat=ierr)
         if(ierr /= 0) call terminate("releaseCommPattern", &
                                      "Deallocation error for recvList")

       enddo

       deallocate(commPattern%sendProc,       &
                  commPattern%recvProc,       &
                  commPattern%indexSendProc, &
                  commPattern%indexRecvProc, &
                  commPattern%nsend,           &
                  commPattern%nrecv,           &
                  commPattern%nsendCum,       &
                  commPattern%nrecvCum,       &
                  commPattern%sendList,       &
                  commPattern%recvList,       &
                  stat=ierr)
       if(ierr /= 0)                                 &
         call terminate("releaseCommPattern", &
                        "Deallocation error for commPattern")

       end subroutine releaseCommPattern

!      ==================================================================

       subroutine releaseInternalComm(internal)
!
!      ******************************************************************
!      *                                                                *
!      * ReleaseInternalComm releases memory of the given argument      *
!      * which should have type internalComm defined in the module      *
!      * communication.
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       type(internalCommType), intent(inout) :: internal
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       deallocate(internal%donorBlock,   &
                  internal%donorIndices, &
                  internal%donorInterp,  &
                  internal%haloBlock,    &
                  internal%haloIndices,  &
                  stat=ierr)
       if(ierr /= 0)                                &
         call terminate("releaseInternalComm", &
                        "Deallocation error for internal")

       end subroutine releaseInternalComm

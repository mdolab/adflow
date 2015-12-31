!
!      ******************************************************************
!      *                                                                *
!      * File:          setBufferSizes.f90                              *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 02-21-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBufferSizes(level, sps, determine1to1Buf,    &
                                             determineSlidingBuf, & 
                                             determineOversetBuf)
!
!      ******************************************************************
!      *                                                                *
!      * setBufferSizes determines the size of the send and receive     *
!      * buffers for this grid level. After that the maximum value of   *
!      * these sizes and the currently stored value is taken, such that *
!      * for all mg levels the same buffer can be used. Normally the    *
!      * size on the finest grid should be enough, but it is just as    *
!      * safe to check on all mg levels. A distinction is made between  *
!      * 1 to 1 and sliding mesh communication, because these happen    *
!      * consecutively and not simultaneously. Consequently the actual  *
!      * buffer size is the maximum of the two and not the sum.         *
!      * For steady state computations a mixing plane boundary          *
!      * condition is used instead of a sliding mesh. However, both     *
!      * communication patterns are allocated and initialized.          *
!      * Therefore the maximum of the two can be taken without checking *
!      * the situation we are dealing with.                             *
!      *                                                                *
!      ******************************************************************
!
       use commMixing
       use commSliding
       use communication
       use flowVarRefState
       use inputPhysics
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
       logical, intent(in) :: determine1to1Buf, determineSlidingBuf
       logical, intent(in) :: determineOversetBuf
!
!      Local variables.
!
       integer(kind=intType) :: i
       integer(kind=intType) :: sendSize, recvSize, nVarComm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the maximum number of variables to be communicated.

       nVarComm = nw + 1
       if(cpModel == cpTempCurveFits) nVarComm = nVarComm + 1
       if( viscous )   nVarComm = nVarComm + 1
       if( eddyModel ) nVarComm = nVarComm + 1

       ! Check if the 1 to 1 communication must be considered.

       if( determine1to1Buf ) then

         ! Store the send and receive buffer sizes needed for the nodal
         ! exchange. Determine the maximum for the number of send and
         ! receive processors.

         i = commPatternNode_1st(level)%nProcSend
         sendSize = commPatternNode_1st(level)%nsendCum(i)

         i = commPatternNode_1st(level)%nProcRecv
         recvSize = commPatternNode_1st(level)%nrecvCum(i)

         ! Determine the buffer sizes for the 2nd level cell exchange and
         ! set the size for this processor to the maximum needed. Note
         ! that it is not needed to test the 1st level cell halo, because
         ! it is entirely incorporated in the 2nd level.
         ! Determine the maximum for the number of send and receive
         ! processors as well.

         i = commPatternCell_2nd(level)%nProcSend
         sendSize = max(sendSize, &
                        commPatternCell_2nd(level)%nsendCum(i))

         i = commPatternCell_2nd(level)%nProcRecv
         recvSize = max(recvSize, &
                        commPatternCell_2nd(level)%nrecvCum(i))

         ! Multiply sendSize and recvSize with the number of variables to
         ! be communicated.

         sendSize = sendSize*nVarComm
         recvSize = recvSize*nVarComm

         ! Store the maximum of the current values and the old values
         ! in sendBufferSize1to1 and recvBufferSize1to1.

         sendBufferSize_1to1 = max(sendBufferSize_1to1, sendSize)
         recvBufferSize_1to1 = max(recvBufferSize_1to1, recvSize)

       endif

       ! Check if the sliding mesh communication must be considered.

       if( determineSlidingBuf ) then

         ! Only the second level cell halo communication pattern needs
         ! to be considered. Note that there is no nodal communication
         ! pattern for sliding meshes.

         i = commSlidingCell_2nd(level,sps)%nProcSend
         sendSize = commSlidingCell_2nd(level,sps)%nsendCum(i)

         i = commSlidingCell_2nd(level,sps)%nProcRecv
         recvSize = commSlidingCell_2nd(level,sps)%nrecvCum(i)

         ! Multiply sendSize and recvSize with the number of variables to
         ! be communicated.

         sendSize = sendSize*nVarComm
         recvSize = recvSize*nVarComm

         ! Store the maximum of the current values and the old values
         ! in sendBufferSizeSlide and recvBufferSizeSlide.

         sendBufferSizeSlide = max(sendBufferSizeSlide, sendSize)
         recvBufferSizeSlide = max(recvBufferSizeSlide, recvSize)

         ! Take possible mixing plane boundaries into account.

         sendSize = 0
         do i=1,nInterfaceGroups
           sendSize = max(sendSize,                            &
                          commPatternMixing(level,i,1)%nInter, &
                          commPatternMixing(level,i,2)%nInter)
         enddo

         sendSize = sendSize*nVarComm
         recvSize = sendSize

         ! Store the maximum value in sendBufferSizeSlide and
         ! recvBufferSizeSlide.

         sendBufferSizeSlide = max(sendBufferSizeSlide, sendSize)
         recvBufferSizeSlide = max(recvBufferSizeSlide, recvSize)

       endif

       ! Check if the overset communication must be considered.

       if( determineOversetBuf ) then

         ! Same deal for the overset communication.

         i = commPatternOverset(level,sps)%nProcSend
         sendSize = commPatternOverset(level,sps)%nsendCum(i)

         i = commPatternOverset(level,sps)%nProcRecv
         recvSize = commPatternOverset(level,sps)%nrecvCum(i)

         ! Multiply sendSize and recvSize with the number of variables to
         ! be communicated.

         sendSize = sendSize*nVarComm
         recvSize = recvSize*nVarComm

         ! Store the maximum of the current values and the old values.

         sendBufferSizeOver = max(sendBufferSizeOver, sendSize)
         recvBufferSizeOver = max(recvBufferSizeOver, recvSize)

       endif

       ! Take the maximum for of all the buffers to
       ! obtain the actual size to be allocated.

       sendBufferSize = max(sendBufferSize_1to1, &
                            sendBufferSizeOver,  &
                            sendBufferSizeSlide)
       recvBufferSize = max(recvBufferSize_1to1, &
                            recvBufferSizeOver,  &
                            recvBufferSizeSlide)

       end subroutine setBufferSizes

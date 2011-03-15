!
!      ******************************************************************
!      *                                                                *
!      * File:          commSliding.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-25-2003                                      *
!      * Last modified: 03-21-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module commSliding
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the communication pattern, the            *
!      * interpolation coefficients and the transformation matrices for *
!      * the sliding mesh interfaces for all the grid levels.           *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the data type slidingCommListType, which is  *
!      * is identical to comm_list_type. However it is duplicated,      *
!      * because in the future it may be different.                     *
!      *                                                                *
!      ******************************************************************
!
       type slidingCommListType

         ! block(..):     Local block ID to which the cell/node belongs.
         !                The dimension is equal to the number of entities
         !                to be communicated with the particular processor.
         ! indices(..,3): i, j and k indices of the data to be communicated.
         !                for the first dimension, see block.

         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices

       end type slidingCommListType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the data type interpolHaloListType, which    *
!      * stores the indices and block IDs as well as the interpolation  *
!      * weights of the halo's to be constructed by interpolation and   *
!      * the indices in the receive buffer. This datatype is intended   *
!      * for external communication via buffers and should not be used  *
!      * for internal communication.                                    *
!      *                                                                *
!      ******************************************************************
!
       type interpolHaloListType

         ! nCopy:             Number of entities to copy from the receive
         !                    buffer. This number can be larger than the
         !                    actual size of the buffer, because entities
         !                    may be needed more than once.
         ! indRecv(nCopy):    Indices in the receive buffer.
         ! block(nCopy):      Local block id to which the entity belongs.
         ! indices(nCopy,..): Indices of the data to be constructed. The
         !                    second dimension is 3 in case a halo of a
         !                    normal block is constructed and 2 if a
         !                    halo of conservative interface_type is
         !                    constructed.
         ! weight(nCopy):     Interpolation weights.

         integer(kind=intType) :: nCopy
         integer(kind=intType), pointer, dimension(:)   :: indRecv
         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices
         real(kind=realType),   pointer, dimension(:)   :: weight

       end type interpolHaloListType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type slidingCommType,       *
!      * which stores the communication pattern for sliding mesh        *
!      * interfaces for a certain halo type.                            *
!      *                                                                *
!      ******************************************************************
!
       type slidingCommType

         ! nProcSend:           # of procs, to whom this proc will send.
         ! nProcRecv:           # of procs, from whom this proc will 
         !                      receive.
         ! sendProc(nProcSend): Send processor numbers.
         ! recvProc(nProcRecv): Receive processor numbers.

         integer(kind=intType)  :: nProcSend, nProcRecv
         integer(kind=intType), pointer, dimension(:) :: sendProc
         integer(kind=intType), pointer, dimension(:) :: recvProc

         ! nSend(nProcSend): # of entities to send to other processors.
         ! nRecv(nProcRecv): # of entities to receive from other
         !                   processors.

         integer(kind=intType), pointer, dimension(:) :: nSend, nRecv

         ! nSendCum(0:nProcSend): Cumulative version of nSend.
         ! nRecvCum(0:nProcRecv): Cumulative version of nRecv.

         integer(kind=intType), pointer, dimension(:) :: nSendCum
         integer(kind=intType), pointer, dimension(:) :: nRecvCum

         ! sendList(nProcSend): Indices and block IDs to send to these
         !                      processors.

         type(slidingCommListType), pointer, dimension(:) :: sendList

         ! recvList(nProcRecv): Halo indices and block ids as well as
         !                      the interpolation weights to construct
         !                      the halo's.

         type(interpolHaloListType), pointer, dimension(:) :: recvList

       end type slidingCommType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the data type internalSlidingCommType,       *
!      * which stores the memory to memory copy on this processors for  *
!      * sliding mesh interfaces for a certain halo type.               *
!      *                                                                *
!      ******************************************************************
!
       type internalSlidingCommType

         ! nSlidingHalos:           Number of sliding halos stored on
         !                          this processor.
         ! slidingHaloList:         Indices and block id's of the
         !                          sliding mesh halo's.
         ! rotIndex(nSlidingHalos): Which rotation matrix to be used
         !                          for the velocities; 0 means that
         !                          velocities do not need to be rotated.

         integer(kind=intType)     :: nSlidingHalos
         type(slidingCommListType) :: slidingHaloList

         integer(kind=intType), dimension(:), pointer :: rotIndex

         ! nCopy:         Number of entities to copy locally.
         ! donorList:     Indices and block id's of the donors.
         ! haloList:      Indices and block id's of the halo's.
         ! weight(nCopy): Interpolation weights.

         integer(kind=intType)     :: nCopy
         type(slidingCommListType) :: donorList, haloList

         real(kind=realType), dimension(:), pointer :: weight

       end type internalSlidingCommType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! nRotSliding:        The number of rotation matrices needed to
       !                     compute the velocities in the sliding mesh
       !                     halo's on this processor. This value is 0
       !                     if there are no periodic boundaries. the
       !                     unity transformation is not included.
       ! rotSliding(..,3,3): The actual rotation matrices. The first
       !                     dimension equals nRotSliding.

       integer(kind=intType) :: nRotSliding
       real(kind=realType), dimension(:,:,:), allocatable :: rotSliding

       ! commSlidingCell_1st(:,:): The communication pattern for the
       !                           1st level nonconservative sliding
       !                           mesh cell halo's on the multiple
       !                           grids. Dimensions are (nLevel,
       !                           nTimeIntervalsSpectral)
       ! commSlidingCell_2nd(:,:): The communication pattern for the
       !                           2nd level nonconservative sliding
       !                           mesh cell halo's on the multiple
       !                           grids. Same dimensions as 1st level.

       type(slidingCommType), allocatable, dimension(:,:) :: &
                                            commSlidingCell_1st
       type(slidingCommType), allocatable, dimension(:,:) :: &
                                            commSlidingCell_2nd

       ! intSlidingCell_1st(:,:): Memory to memory copies for 1st
       !                          level nonconservative sliding mesh
       !                          cell halo's on the multiple grids.
       !                          dimensions are
       !                          (nLevel,nTimeIntervalsSpectral)
       ! intSlidingCell_2nd(:,:): Memory to memory copies for 2nd
       !                          level nonconservative sliding mesh
       !                          cell halo's on the multiple grids.
       !                          Same dimensions as 1st level.

       type(internalSlidingCommType), allocatable, dimension(:,:) :: &
                                                     intSlidingCell_1st
       type(internalSlidingCommType), allocatable, dimension(:,:) :: &
                                                     intSlidingCell_2nd

       ! sendBufferSizeSlide: Size of the send buffer needed to
       !                      perform all sliding mesh communication.
       ! recvBufferSizeSlide: Idem for the receive buffer.

       integer(kind=intType) :: sendBufferSizeSlide
       integer(kind=intType) :: recvBuffersizeSlide

       end module commSliding

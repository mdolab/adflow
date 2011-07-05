!
!      ******************************************************************
!      *                                                                *
!      * File:          communication.F90                               *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module communication
!
!      ******************************************************************
!      *                                                                *
!      * Contains the variable definition of the processor number,      *
!      * myID and the number of processors, nProc, which belong to the  *
!      * group defined by the communicator SUmb_comm_world. The range   *
!      * of processor numbers is <0..Nproc-1>, i.e. the numbering       *
!      * starts at 0. This is done for compatibility with MPI.          *
!      * Furthermore this module contains the communication pattern for *
!      * all the multigrid levels.                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type commListType, which    *
!      * stores the i,j and k indices as well as the block id of the    *
!      * data to be communicated. Send lists may contain interpolants   *
!      * since the indices may refer to a stencil, while the receive    *
!      * list does not. All interpolations should be done on the send   *
!      * side to keep message sizes to a minimum.                       *
!      *                                                                *
!      ******************************************************************
!
       type sendCommListType

         ! block(..):     Local block id to which the cell/node belongs.
         !                The dimension is equal to the number of entities
         !                to be communicated with the particular processor.
         ! indices(..,3): I, j and k indices of the data to be communicated.
         !                For the first dimension, see block.
         ! interp(..,..): Interpolants for indices that represent a cell
         !                stencil (allocated only when needed, e.g. send
         !                list for overset communication).
         !                For the first dimension, see block.

         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices
         real(kind=realType),   pointer, dimension(:,:) :: interp

       end type sendCommListType
 
       type recvCommListType

         ! block(..):     Local block id to which the cell/node belongs.
         !                The dimension is equal to the number of entities
         !                to be communicated with the particular processor.
         ! indices(..,3): I, j and k indices of the data to be communicated.
         !                For the first dimension, see block.

         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices

       end type recvCommListType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type periodicDataType,      *
!      * which stores the rotation matrix, the rotation center and the  *
!      * translation vector of the periodic transformation, as well as  *
!      * the halos to which this transformation must be applied.        *
!      *                                                                *
!      ******************************************************************
!
       type periodicDataType

         ! rotMatrix(3,3): Rotation matrix.
         ! rotCenter(3):   Coordinates of center of rotation.
         ! translation(3):  Translation vector.

         real(kind=realType), dimension(3,3) :: rotMatrix
         real(kind=realType), dimension(3)   :: rotCenter, translation

         ! nHalos:            # of halos to which this periodic
         !                    transformation must be applied.
         ! block(nHalos):     Local block id to which the halos belong.
         ! indices(nHalos,3): I, j and k indices of the halos in
         !                    this block.

         integer(kind=intType) :: nHalos

         integer(kind=intType), pointer, dimension(:)   :: block
         integer(kind=intType), pointer, dimension(:,:) :: indices

       end type periodicDataType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type commType, which        *
!      * stores the communication pattern for a certain halo type for a *
!      * certain grid level.                                            *
!      *                                                                *
!      ******************************************************************
!
       type commType

         ! nProcSend:           # of procs, to whom this proc will send
         ! nProcRecv:           # of procs, from whom this proc will
         !                      receive.
         ! sendProc(nProcSend): Send processor numbers.
         ! recvProc(nProcRecv): Receive processor numbers.

         integer(kind=intType)  :: nProcSend, nProcRecv
         integer(kind=intType), pointer, dimension(:) :: sendProc, recvProc

         ! nsend(nProcSend): # of entities to send to other processors.
         ! nrecv(nProcRecv): # of entities to receive from other processors.

         integer(kind=intType), pointer, dimension(:) :: nsend, nrecv

         ! nsendCum(0:NprocSend): cumulative version of nsend.
         ! nrecvCum(0:NprocRecv): cumulative version of nrecv.

         integer(kind=intType), pointer, dimension(:) :: nsendCum
         integer(kind=intType), pointer, dimension(:) :: nrecvCum

         ! indexSendProc(0:Nproc-1): index of the processors in sendProc.
         !                           If nothing is to be sent to a
         !                           processor this value is 0.
         ! indexRecvProc(0:Nproc-1): index of the processors in recvProc.
         !                           If nothing is to be received from a
         !                           processor this value is 0.

         integer(kind=intType), pointer, dimension(:) :: indexSendProc
         integer(kind=intType), pointer, dimension(:) :: indexRecvProc

         ! sendList(nProcSend): Indices and block ids to send to these
         !                      processors.
         ! recvList(nProcRecv): Indices and block ids to receive from
         !                      these processors.

         type(sendCommListType), pointer, dimension(:) :: sendList
         type(recvCommListType), pointer, dimension(:) :: recvList

         ! nPeriodic:               # of periodic data arrays.
         ! periodicData(nPeriodic): Periodic data and entities to which
         !                          the transformation must be applied.

         integer(kind=intType) :: nPeriodic
         type(periodicDataType), pointer, dimension(:) :: periodicData

       end type commType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type internalCommType,      *
!      * which stores the memory to memory copy on this processor for a *
!      * certain halo type.                                             *
!      *                                                                *
!      ******************************************************************
!
       type internalCommType

         ! ncopy: # of entities to copy internally.

         ! donorBlock(ncopy):     Local block id of the donor cell.
         ! donorIndices(ncopy,3): The indices of the donor cell.
         ! donorInterp(ncopy,3):  Interpolants of the donor stencil.
         !                        Only allocated when needed, e.g.
         !                        overset communication).

         ! haloBlock(ncopy):     Local block id of the halo cell.
         ! haloIndices(ncopy,3): The indices of the halo cell.

         integer(kind=intType) :: ncopy

         integer(kind=intType), pointer, dimension(:)   :: donorBlock
         integer(kind=intType), pointer, dimension(:,:) :: donorIndices
         real(kind=realType),   pointer, dimension(:,:) :: donorInterp

         integer(kind=intType), pointer, dimension(:)   :: haloBlock
         integer(kind=intType), pointer, dimension(:,:) :: haloIndices

         ! nPeriodic:               # of periodic data arrays.
         ! periodicData(nPeriodic): Periodic data and entities to which
         !                          the transformation must be applied.

         integer(kind=intType) :: nPeriodic
         type(periodicDataType), pointer, dimension(:) :: periodicData

       end type internalCommType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! SUmb_comm_world: The communicator of this processor group.
       ! myID:            My processor number in SUmb_comm_world.
       ! nProc:           The number of processors in SUmb_comm_world.

       integer :: SUmb_comm_world, myID, nProc

       ! commPatternCell_1st(nLevel): The communication pattern for 1st
       !                              level cell halo's on the multiple
       !                              grids.
       ! commPatternCell_2nd(nLevel): The communication pattern for 2nd
       !                              level cell halo's (including the
       !                              1st level) on the multiple grids.
       ! commPatternNode_1st(nLevel): The communication pattern for 1st
       !                              level node halo's on the multiple
       !                              grids.
       ! CommPatternOverset(nLevel,   the communication pattern for
       !                      nsps):  Overset halos on the multiple
       !                              grids.

       type(commType), allocatable, dimension(:)   :: commPatternCell_1st
       type(commType), allocatable, dimension(:)   :: commPatternCell_2nd
       type(commType), allocatable, dimension(:)   :: commPatternNode_1st
       type(commType), allocatable, dimension(:,:) :: commPatternOverset

       ! internalCell_1st(nLevel): Memory to memory copies for 1st level
       !                           cell halo's on the multiple grids.
       ! internalCell_2nd(nLevel): Memory to memory copies for 2nd level
       !                           cell halo's on the multiple grids.
       ! internalNode_1st(nLevel): Memory to memory copies for 1st level
       !                           node halo's on the multiple grids.
       ! InternalOverset(nLevel,   internal communication for overset
       !                  nsps):   Halos on the multiple grids.

       type(internalCommType), allocatable, dimension(:)   :: internalCell_1st
       type(internalCommType), allocatable, dimension(:)   :: internalCell_2nd
       type(internalCommType), allocatable, dimension(:)   :: internalNode_1st
       type(internalCommType), allocatable, dimension(:,:) :: internalOverset

       ! sendBufferSize_1to1: Size of the send buffer needed to perform
       !                      all 1 to 1 communication.
       ! recvBufferSize_1to1: Idem for the receive buffer.
       ! sendBufferSizeOver:  Size of the send buffer needed to perform
       !                      all overset communication.
       ! recvBufferSizeOver:  Idem for the receive buffer.
       ! sendBufferSize:      Size of the send buffer to perform all
       !                      possible communication.
       ! recvBufferSize:      Idem for the receive buffer.

       integer(kind=intType) :: sendBufferSize_1to1, sendBufferSize
       integer(kind=intType) :: recvBufferSize_1to1, recvBufferSize
       integer(kind=intType) :: sendBufferSizeOver, recvBufferSizeOver

       ! sendBuffer:   Buffer used to store the info to be send during
       !               a nonblocking communication.
       ! recvBuffer:   Buffer used to store the info to be received
       !               during a nonblocking communication.
       ! sendRequests: Array of requests for the nonblocking sends.
       ! recvRequests: Array of requests for the nonblocking receives.

       real(kind=realType), allocatable, dimension(:) :: sendBuffer
       real(kind=realType), allocatable, dimension(:) :: recvBuffer

       integer, allocatable, dimension(:) :: sendRequests, recvRequests

       end module communication

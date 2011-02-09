!
!      ******************************************************************
!      *                                                                *
!      * File:          IOModule.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-02-2005                                      *
!      * Last modified: 10-31-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module IOModule
!
!      ******************************************************************
!      *                                                                *
!      * Constants and variables used in the IO routines.               *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the parameters.                                  *
!      *                                                                *
!      ******************************************************************
!
       ! The maximum amount of data a processor can read in one call
       ! to su_file_read. This value is 2Gbyte.
 
       integer(kind=mpi_offset_kind), parameter :: &
                                  maxSizeIO = 2147483648_mpi_offset_kind

       ! Definition of the possibilities for storing data.

       integer(kind=intType), parameter :: nodeData          = 1
       integer(kind=intType), parameter :: cellDataNoHalo    = 2
       integer(kind=intType), parameter :: cellDataPlusHalo  = 3

       ! Definition of the possible Plot3D block formats.

       integer(kind=intType), parameter :: P3D_SingleBlock = 1
       integer(kind=intType), parameter :: P3D_MultiBlock  = 2
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the derived datatype IOType, which is used to    *
!      * to make the IO as general as needed.                           *
!      *                                                                *
!      ******************************************************************
!
       type IOType

         ! pointerOffset: offset due to the usage of a pointer to a
         !                subarray. The Fortran standard is such that
         !                the starting indices of the pointer array
         !                is 1, no matter what the original starting
         !                index is. This can lead to a shift in the
         !                indices.
         ! w:             The variable(s) to be read/written.

         integer(kind=intType) :: pointerOffset
         real(kind=realType), dimension(:,:,:,:), pointer :: w

       end type IOType
!
!      ******************************************************************
!      *                                                                *
!      * Definition P3D_IOPartType, which is used to distinguish in the *
!      * IO buffer the physical blocks.                                 *
!      *                                                                *
!      ******************************************************************
!
       type P3D_IOPartType

         ! offsetIO:       Offset in bytes relative to the start of the
         !                 local IO buffer.
         ! offsetBuffer:   Offset to be used for the buffer when actual
         !                 data is extracted from the IO buffer. On some
         !                 platforms it is needed to apply a shift for
         !                 offsetIO if this value is not an integer
         !                 multiple of the size of the real type.
         ! nItemsTotal:    Number of items read/written for this part.
         ! nItemsLocal:    Number of items that should be copied to/from
         !                 the local buffers.
         ! nItemsNonLocal: Number of items that should be communicated
         !                 to/from other processors.

         integer(kind=intType) :: offsetIO, offsetBuffer
         integer(kind=intType) :: nItemsTotal
         integer(kind=intType) :: nItemsLocal, nItemsNonLocal

         ! blockID(nItemsLocal):   The block ID's of the local items.
         ! indices(nItemsLocal,3): The corresponding i, j, k indices.
         ! indexW(nItemsLocal):    The index in the w variables of IOVar
         !                         assuming that the starting index is
         !                         0. The user must define the offset in
         !                         the argument list of the actual
         !                         routines that perform the IO, i.e.
         !                         readPlot3DVar and writePlot3DVar.
         ! posLocal(nItemsLocal):  The corresponding positions in
         !                         the buffer.

         integer(kind=intType), dimension(:),   pointer :: blockID
         integer(kind=intType), dimension(:),   pointer :: indexW
         integer(kind=intType), dimension(:),   pointer :: posLocal
         integer(kind=intType), dimension(:,:), pointer :: indices

         ! posComm(nItemsNonLocal):     The index in the communication
         !                              buffer.
         ! posNonLocal(nItemsNonLocal): The corresponding positions
         !                              in the IO buffer.

         integer(kind=intType), dimension(:), pointer :: posComm
         integer(kind=intType), dimension(:), pointer :: posNonLocal

       end type P3D_IOPartType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the variables used for both CGNS and Plot3D.     *
!      *                                                                *
!      ******************************************************************
!
       ! IOVar(nDom,nIOFiles): Array of the derived datatype IOType to
       !                       facilitate a general IO implementation.

       type(IOType), dimension(:,:), allocatable :: IOVar
!
!      ******************************************************************
!      *                                                                *
!      * Plot3D variables (only relevant if Plot3D is used).            *
!      *                                                                *
!      ******************************************************************
!
       ! fh:              File handler.
       ! P3D_BlockFormat: What type of grid is stored, single block
       !                  or multi-block.
       ! P3D_Precision:   Precision of the floating point type, 
       !                  possibilities are precisionSingle and
       !                  precisionDouble.
       ! P3D_DataStorage: Type of storage for the data to be
       !                  read/written.  Possibilities are nodeData,
       !                  cellDataNoHalo and cellDataPlusHalo.
       ! P3D_Offset:      Offset in bytes from the beginning of the
       !                  file where the IO must be performed.
       ! P3D_ByteSwap:    Whether or not byte swapping must be applied.
       ! P3D_iblank:      Whether or not an iblanking array is present.
       ! P3D_nVar:        Number of variables to be read/written during
       !                  one call to the IO routines. E.g. for the
       !                  coordinates this value is 3, for solution
       !                  variables it is 1.

       integer :: fh

       integer(kind=intType) :: P3D_BlockFormat, P3D_Precision
       integer(kind=intType) :: P3D_DataStorage, P3D_nVar

       integer(kind=mpi_offset_kind) :: P3D_Offset

       logical :: P3D_ByteSwap, P3D_iblank
!
!      ******************************************************************
!      *                                                                *
!      * Plot3D variables only used in the routines that perform the    *
!      * actual IO.                                                     *
!      *                                                                *
!      ******************************************************************
!
       ! P3D_mySizeIO:   Number of bytes I have to read/write. This is
       !                 a normal integer, because MPI requires this.
       !                 This means that 1 processor cannot read/write
       !                 a buffer that is larger than 2 GByte. I don't
       !                 think that is a real restriction.
       ! P3D_myOffset:   Relative offset to P3D_Offset where I have to
       !                 read/write data.
       ! P3D_nIOParts:   Number of IO parts in which the local IO buffer
       !                 is split. This corresponds to the number of
       !                 global block (parts) I'm responsible for.
       ! P3D_IOParts(:): The corresponding array of the derived datatype
       !                 to store the necessary data.

       integer :: P3D_mySizeIO

       integer(kind=mpi_offset_kind) :: P3D_myOffset

       integer(kind=intType) :: P3D_nIOParts

       type(P3D_IOPartType), dimension(:), allocatable :: P3D_IOParts

       ! P3D_nProcSend:                 Number of processors to
       !                                which I have to send data.
       ! P3D_nProcRecv:                 Number of processors from
       !                                which I receive data.
       ! P3D_sendSize(0:P3D_nProcSend): Size of the send messages in
       !                                bytes, cumulative storage format.
       ! P3D_recvSize(0:P3D_nProcRecv): Size of the receive messages in
       !                                bytes, cumulative storage format.
       ! P3D_procSend(P3D_nProcSend):   Corresponding processors.
       ! P3D_procRecv(P3D_nProcRecv):   Processor ID's from which I
       !                                receive data.
       ! P3D_commPart:                  Derived data type to copy
       !                                data to/from IOVar from/to the
       !                                communication buffer.

       integer(kind=intType) :: P3D_nProcSend, P3D_nProcRecv

       integer(kind=intType), dimension(:), allocatable :: P3D_procSend
       integer(kind=intType), dimension(:), allocatable :: P3D_procRecv
       integer(kind=intType), dimension(:), allocatable :: P3D_sendSize
       integer(kind=intType), dimension(:), allocatable :: P3D_recvSize

       type(P3D_IOPartType) :: P3D_commPart

       ! P3D_nRecordIntegersWrite:   The number of integer records I
       !                             am responsible for to write.
       ! P3D_recordIntegersWrite(:): The corresponding sizes in bytes.
       ! P3D_recordPosition(:):      The position in the local write
       !                             buffer where these integers should
       !                             be stored.

       integer(kind=intType) :: P3D_nRecordIntegersWrite

       integer(kind=intType), dimension(:), allocatable :: &
                                                 P3D_recordIntegersWrite
       integer(kind=intType), dimension(:), allocatable :: &
                                                      P3D_recordPosition

       end module IOModule

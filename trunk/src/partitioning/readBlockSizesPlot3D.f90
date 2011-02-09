!
!      ******************************************************************
!      *                                                                *
!      * File:          readBlockSizesPlot3D.F90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-16-2005                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readBlockSizesPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readBlockSizesPlot3D determines the actual format of the       *
!      * plot3D file and reads the number of blocks and their sizes     *
!      * from this file. The data is stored in the module cgnsGrid.     *
!      * This name comes from historic reasons, because the plot3D      *
!      * option was added later on; in the early versions only CGNS     *
!      * format was supported.                                          *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use IOModule
       use partitionMod
       implicit none
!
!      Local variables
!
       integer :: sizeBuf, ierr

       integer, dimension(mpi_status_size) :: status

       integer(kind=intPLOT3DType) :: intPlot3D
       integer(kind=intRecordPLOT3DType), dimension(2) :: recordSize

       integer(kind=intPLOT3DType), dimension(:,:), allocatable :: &
                                                            intPlot3DBuf
       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record

       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Allocate the memory for the variables to store the information
       ! on the byte swapping and the block format.

       allocate(blockFormatGrids(nGridsRead), &
                byteSwapGrids(nGridsRead), stat=ierr)
       if(ierr /= 0)                             &
          call terminate("readBlockSizesPlot3D", &
                         "Memory allocation failure for &
                         &blockFormatGrids and byteSwapGrids")
!
!      ******************************************************************
!      *                                                                *
!      * Loop over the grids to be read to perform some initialization  *
!      * tasks.                                                         *
!      *                                                                *
!      ******************************************************************
!
       gridInitLoop: do mm=1,nGridsRead

         ! Open the grid file for reading and check if it went okay.

         call mpi_file_open(SUmb_comm_world, gridFiles(mm), &
                           mpi_mode_rdonly, mpi_info_null, &
                           fh, ierr)

         if(ierr /= mpi_success) then
           write(errorMessage,"(3a)") "plot3D grid file ", &
                                      trim(gridFiles(mm)), " not found"
           if(myID == 0) &
             call terminate("readBlockSizesPlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Store the file handler for later purposes.

         fileIDs(mm) = fh
!
!        ****************************************************************
!        *                                                              *
!        * Determine the format of the plot3D file. It is assumed that  *
!        * this file is an unformatted fortran file. At the moment the  *
!        * following options are recognized:                            *
!        * - whether or not byte swapping must be applied.              *
!        * - single or multi block format.                              *
!        *                                                              *
!        * The precision of the grid file and whether or not iblanks    *
!        * are stored is determined when the coordinates are read in    *
!        * the routine readGridPlot3D.                                  *
!        *                                                              *
!        ****************************************************************
!
         ! Read the size of the first fortran record. With this
         ! information it can be determined whether or not byte swapping
         ! must be applied and the block format of the file.
         ! All the processors read it, because that is the most efficient
         ! way of doing it. A collective read is nonsense here.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)

         ! Determine whether or not byte swapping must be applied.
         ! It is stored for later purposes.

         if(recordSize(1) ==    sizeP3D_Int .or. &
            recordSize(1) ==  3*sizeP3D_Int) then
           P3D_ByteSwap = .false.
         else
           call byteswap(recordSize, sizeP3D_Record, 1_intType)
           P3D_ByteSwap = .true.
         endif

         byteSwapGrids(mm) = P3D_ByteSwap

         ! Determine the block format. If no valid format is found an
         ! error message is printed and the program terminates.
         ! The block format is stored for later purposes.

         if(recordSize(1) == sizeP3D_Int) then
           P3D_BlockFormat = P3D_MultiBlock
         else if(recordSize(1) ==  3*sizeP3D_Int) then
           P3D_BlockFormat = P3D_SingleBlock
         else
           write(errorMessage,"(3a)") "Unknown plot3D format for grid &
                                       &file ", trim(gridFiles(mm)), "."
           if(myID == 0) &
             call terminate("readBlockSizesPlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         blockFormatGrids(mm) = P3D_BlockFormat
!
!        ****************************************************************
!        *                                                              *
!        * Determine the number of blocks and the their sizes.          *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the number of blocks. For a single block format this
         ! of course is 1; for a multi-block format the number is read.

         testMultiBlock: if(P3D_BlockFormat == P3D_MultiBlock) then

           call mpi_file_read(fh, intPlot3D, 1, sumb_integerPLOT3D, &
                             status, ierr)
           if( P3D_ByteSwap ) &
             call byteswap(intPlot3D, sizeP3D_Int, 1_intType)

           ! Read the closing integer of the 1st record and the opening
           ! integer of the 2nd record.

           call mpi_file_read(fh, recordSize, 2, &
                              sumb_integerRecordPLOT3D, status, ierr)

         else testMultiBlock

           intPlot3D = 1

         endif testMultiBlock

         ! Store the number of blocks if this is the first grid to be
         ! read and check if the number of blocks is identical.

         if(mm == 1) cgnsNDom = intPlot3D

         if(intPlot3D /= cgnsNDom) then
           write(errorMessage,"(4a)") "plot3D grid file ", &
                                      trim(gridFiles(mm)), &
                                      ": Different number of &
                                      &blocks than in file ", &
                                      trim(gridFiles(1))
           if(myID == 0) &
             call terminate("readBlockSizesPlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

       enddo gridInitLoop
!
!      ******************************************************************
!      *                                                                *
!      * Determine and check the block sizes.                           *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for cgnsDoms, the derived data type to
       ! store the info of the grid in the file, and for the read buffer
       ! to do the actual reading. Nullify the pointers of cgnsDoms.

       allocate(cgnsDoms(cgnsNDom), intPlot3DBuf(3,cgnsNDom), &
                stat=ierr)
       if(ierr /= 0)                            &
         call terminate("readBlockSizesPlot3D", &
                        "Memory allocation failure for cgnsDoms and &
                        &intPlot3DBuf")

       do nn=1,cgnsNDom
         call nullifyCGNSDomPointers(nn)
       enddo

       ! Loop over the grid files to read the block sizes.

       gridSizeLoop: do mm=1,nGridsRead

         ! Easier storage of the file pointer and the byte swapping.

         fh = fileIDs(mm)
         P3D_ByteSwap = byteSwapGrids(mm)

         ! Read the block sizes. Allocate the memory for the
         ! read buffer, read it and apply the byte swapping if needed.

         nn      = 3*cgnsNDom
         sizeBuf = nn
         call mpi_file_read(fh, intPlot3DBuf, sizeBuf, &
                           sumb_integerPLOT3D, status, ierr)

         if( P3D_ByteSwap ) &
           call byteswap(intPlot3DBuf, sizeP3D_Int, nn)

         ! Copy the data into cgnsDoms if this is the first grid.

         if(mm == 1) then
           do nn=1,cgnsNDom
             cgnsDoms(nn)%il = intPlot3DBuf(1,nn)
             cgnsDoms(nn)%jl = intPlot3DBuf(2,nn)
             cgnsDoms(nn)%kl = intPlot3DBuf(3,nn)
           enddo
         endif

         ! Check if the block sizes are the same.

         do nn=1,cgnsNDom
           if(cgnsDoms(nn)%il /= intPlot3DBuf(1,nn) .or. &
              cgnsDoms(nn)%jl /= intPlot3DBuf(2,nn) .or. &
              cgnsDoms(nn)%kl /= intPlot3DBuf(3,nn)) then
             write(errorMessage,"(4a)") "plot3D grid file ", &
                                        trim(gridFiles(mm)), &
                                        ": Different block sizes &
                                        &than in file ", &
                                        trim(gridFiles(1))
             if(myID == 0) &
               call terminate("readBlockSizesPlot3D", errorMessage)
             call mpi_barrier(SUmb_comm_world, ierr)
           endif
         enddo

       enddo gridSizeLoop

       ! Release the memory of intPlot3DBuf again.

       deallocate(intPlot3DBuf, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("readBlockSizesPlot3D", &
                        "Deallocation failure for intPlot3DBuf.")

       ! Determine the number of cells in i, j and k-direction for
       ! every block.

       do nn=1,cgnsNDom
         cgnsDoms(nn)%nx = cgnsDoms(nn)%il - 1
         cgnsDoms(nn)%ny = cgnsDoms(nn)%jl - 1
         cgnsDoms(nn)%nz = cgnsDoms(nn)%kl - 1
       enddo

       end subroutine readBlockSizesPlot3D

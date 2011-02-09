!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DGridFile.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2005                                      *
!      * Last modified: 10-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DGridFile
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DGridFile writes the grid file(s) in Plot3D format   *
!      * using parallel IO.                                             *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use flowVarRefState
       use inputIO
       use IOModule
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr, amode, sizeHeader

       integer, dimension(mpi_status_size) :: status

       integer(kind=intRecordPLOT3DType) :: sizeRecord

       integer(kind=intPLOT3DType), dimension(3) :: intPlot3DBuf

       integer(kind=intType) :: ii, nn, mm
       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record

       integer(kind=intType), dimension(cgnsNDom) :: cgnsBlockIDs

       real(kind=realType) :: LRefInv

       character(len=maxStringLen) :: errorMessage

       character, dimension(:), allocatable :: bufHeader
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number and names of the grid files.
       ! Also set the pointers for IOVar needed for the general
       ! treatment of the IO.

       call gridFileNamesWrite

       ! Return immediately if no grids have to be written.

       if(nGridsToWrite == 0) return

       ! Set the variables needed for the writing of the grid files.
       ! These are the same for all the grid files to be written.

       P3D_Precision   = precisionGrid
       P3D_iblank      = .false.
       P3D_DataStorage = nodeData
       P3D_nVar        = 3

       ! All the grid files to be written will be of identical format.
       ! Therefore the preparations can be done outside the loop
       ! over the number of files. All blocks are written.

       do ii=1,cgnsNDom
         cgnsBlockIDs(ii) = ii
       enddo

       call prepareWritePlot3D(cgnsNDom, cgnsBlockIDs)

       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Write a message that the grid file(s) is (are) written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "# Writing grid file(s) ..."
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the header of the Plot3D file and store it in        *
!      * bufHeader.                                                     *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the size, in bytes, of the header.

       sizeHeader = 2*sizeP3D_Record + 3*cgnsNDom*sizeP3D_Int
       if(P3D_BlockFormat == P3D_MultiBlock) &
         sizeHeader = sizeHeader + 2*sizeP3D_Record + sizeP3D_Int

       ! Allocate the memory for the file header.

       allocate(bufHeader(sizeHeader), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writePlot3DGridFile", &
                        "Memory allocation failure for bufHeader")

       ! Write the header to the character buffer. First the number of
       ! blocks in case of multi block format.

       ii = 1

       testMB: if(P3D_BlockFormat == P3D_MultiBlock) then

         ! The opening integer of the record.

         sizeRecord = sizeP3D_Int
         if( P3D_ByteSwap ) &
           call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
         call writeVarsToBuffer(bufHeader(ii), sizeRecord, &
                                sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! The number of blocks.

         intPlot3DBuf(1) = cgnsNDom
         if( P3D_ByteSwap ) &
           call byteswap(intPlot3DBuf, sizeP3D_Int, 1_intType)
         call writeVarsToBuffer(bufHeader(ii), intPlot3DBuf, sizeP3D_Int)

         ii = ii + sizeP3D_Int

         ! The closing integer of the record.

         call writeVarsToBuffer(bufHeader(ii), sizeRecord, &
                                sizeP3D_Record)
         ii = ii + sizeP3D_Record

       endif testMB

       ! The dimensions of the blocks.
       ! First the opening integer of the record.

       sizeRecord = 3*cgnsNDom*sizeP3D_Int
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(bufHeader(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the block sizes to the character buffer.

       do nn=1,cgnsNDom
         intPlot3DBuf(1) = cgnsDoms(nn)%il
         intPlot3DBuf(2) = cgnsDoms(nn)%jl
         intPlot3DBuf(3) = cgnsDoms(nn)%kl

         if( P3D_ByteSwap ) &
           call byteswap(intPlot3DBuf, sizeP3D_Int, 3_intType)

         mm = 3*sizeP3D_Int
         call writeVarsToBuffer(bufHeader(ii), intPlot3DBuf, mm)
         ii = ii + mm
       enddo

       ! The closing integer of the record.

       call writeVarsToBuffer(bufHeader(ii), sizeRecord, sizeP3D_Record)
!
!      ******************************************************************
!      *                                                                *
!      * Loop over the number of grid files to be written and write     *
!      * the coordinates.                                               *
!      *                                                                *
!      ******************************************************************
!
       gridLoop: do nn=1,nGridsToWrite

         ! Open the grid file for writing. There is no guarantee that an
         ! old file will be clobbered, so wipe out any previous output
         ! file. If the file cannot be opened, print an error message
         ! and exit.

         if( myID == 0) &
           call mpi_file_delete(gridFileNames(nn), mpi_info_null, ierr)

         amode = mpi_mode_create + mpi_mode_wronly
         call mpi_file_open(SUmb_comm_world, gridFileNames(nn), &
                           amode, mpi_info_null, fh, ierr)

         if(ierr /= mpi_success) then
           write(errorMessage,*) "File ", trim(gridFileNames(nn)), &
                                 " could not be opened for writing"
           if(myID == 0) &
             call terminate("writePlot3DGridFile", errorMessage)

           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Only processor 0 does the writing of the header, because it
         ! is a limited amount of data.

         if(myID == 0) &
           call mpi_file_write(fh, bufHeader, sizeHeader, mpi_character, &
                              status, ierr)

         ! Compute the multiplication factor to obtain the original
         ! coordinates. Note that LRef is corrected to 1.0 when the
         ! coordinates should be written in meters. This happens when
         ! the grid is read. It is assumed that all blocks have the
         ! same scaling factor, i.e. the original blocks in the grid
         ! file were specified in the same dimensional unit.

         LRefInv = one/LRef

         ! All the pointers to IOVar have already been set.
         ! So just call the general routine to perform the writing.

         P3D_Offset = sizeHeader

         call writePlot3DVar(nn, LRefInv, 1_intType)

         ! Close the grid file.

         call mpi_file_close(fh, ierr)

       enddo gridLoop

       ! Release the memory of the variables needed for the writing.

       call releaseMemIOPlot3D

       ! Deallocate the memory of bufHeader and IOVar.

       deallocate(bufHeader, IOVar, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writePlot3DGridFile", &
                        "Deallocation failure for bufHeader and IOVar")

       ! Write a message that the grid file(s) have been written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "# Grid file(s) written"
         print "(a)", "#"
       endif

       end subroutine writePlot3DGridFile

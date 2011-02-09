!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DZoneFile.f90                         *
!      * Author:        Steve Repsher, Edwin van der Weide              *
!      * Starting date: 09-01-2005                                      *
!      * Last modified: 10-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DZoneFile(izone, zoneFileName, sps, &
                                      P3D_ByteSwap_SUGGAR)
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DZoneFile writes a grid file in Plot3D format for    *
!      * the CGNS zone given by izone and spectral mode sps to the file *
!      * zoneFileName using parallel IO.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use inputIO
       use IOModule
       use suggarData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: izone, sps

       logical, intent(in) :: P3D_ByteSwap_SUGGAR

       character(len=*), intent(in) :: zoneFileName
!
!      Local variables.
!
       integer :: ierr, amode, sizeHeader
       integer :: P3D_BlockFormat_SUGGAR, P3D_RealPrecision_SUGGAR

       integer, dimension(mpi_status_size) :: status

       integer(kind=intRecordPLOT3DType) :: sizeRecord

       integer(kind=intPLOT3DType), dimension(3) :: intPlot3DBuf

       integer(kind=intType) :: ii, nn, mm
       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record

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
       ! Hard set the Plot3D file format and precision for this routine.

       P3D_BlockFormat_SUGGAR   = P3D_MultiBlock
       P3D_RealPrecision_SUGGAR = precisionDouble

       ! Set some variables needed for the writing of the grid file.

       P3D_BlockFormat = P3D_BlockFormat_SUGGAR
       P3D_Precision   = P3D_RealPrecision_SUGGAR
       P3D_DataStorage = nodeData
       P3D_iblank      = .false.
       P3D_ByteSwap    = P3D_ByteSwap_SUGGAR
       P3D_nVar        = 3

       ! Do the preparations for the writing. Only 1 block, izone,
       ! is to be written.

       call prepareWritePlot3D(1_intType, izone)

       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Write a message that the file is being written.

       if (myID == 0) then
         print "( a)", "#"
         print "(3a)", "# Writing Plot3D grid file for zone: ", &
                      trim(adjustl(cgnsDoms(izone)%zoneName)), "..."
       end if
!
!      ******************************************************************
!      *                                                                *
!      * Determine the header of the Plot3D file and store it in        *
!      * bufHeader.                                                     *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the size, in bytes, of the header.

       sizeHeader = 2*sizeP3D_Record + 3*1*sizeP3D_Int
       if(P3D_BlockFormat == P3D_MultiBlock) &
         sizeHeader = sizeHeader + 2*sizeP3D_Record + sizeP3D_Int

       ! Allocate the memory for the file header.

       allocate(bufHeader(sizeHeader), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writePlot3DZoneFile", &
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

         intPlot3DBuf(1) = 1
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

       sizeRecord = 3*1*sizeP3D_Int
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(bufHeader(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the block size to the character buffer.

       intPlot3DBuf(1) = cgnsDoms(izone)%il
       intPlot3DBuf(2) = cgnsDoms(izone)%jl
       intPlot3DBuf(3) = cgnsDoms(izone)%kl

       if( P3D_ByteSwap ) &
         call byteswap(intPlot3DBuf, sizeP3D_Int, 3_intType)

       mm = 3*sizeP3D_Int
       call writeVarsToBuffer(bufHeader(ii), intPlot3DBuf, mm)
       ii = ii + mm

       ! The closing integer of the record.

       call writeVarsToBuffer(bufHeader(ii), sizeRecord, sizeP3D_Record)
!
!      ******************************************************************
!      *                                                                *
!      * The writing of the grid file.                                  *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for IOVar, the data structure which enables
       ! the general IO treatment.

       allocate(IOVar(nDom,1), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writePlot3DZoneFile", &
                        "Memory allocation failure for IOVar")

       ! Set the pointers of IOVar to the coordinates such that the
       ! general writing routine can be used.

       do nn=1,nDom
         IOVar(nn,1)%pointerOffset = 0
         IOVar(nn,1)%w => flowDoms(nn,1,sps)%x(1:,1:,1:,:)
       enddo

       ! Open the grid file for writing. There is no guarantee that an
       ! old file will be clobbered, so wipe out any previous output
       ! file. If the file cannot be opened, print an error message
       ! and exit.

       if( myID == 0) &
         call mpi_file_delete(zoneFileName, mpi_info_null, ierr)

       amode = mpi_mode_create + mpi_mode_wronly
       call mpi_file_open(SUmb_comm_world, zoneFileName, &
                          amode, mpi_info_null, fh, ierr)

       if(ierr /= mpi_success) then
         write(errorMessage,*) "File ", trim(zoneFileName), &
                               " could not be opened for writing"
         if(myID == 0) &
           call terminate("writePlot3DZoneFile", errorMessage)

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
       ! the grid is read.

       LRefInv = one/cgnsDoms(izone)%LRef

       ! Call the general routine to perform the writing.

       P3D_Offset = sizeHeader

       call writePlot3DVar(1_intType, LRefInv, 1_intType)

       ! Close the grid file.

       call mpi_file_close(fh, ierr)

       ! Release the memory of the variables needed for the writing.

       call releaseMemIOPlot3D

       ! Deallocate the memory of bufHeader and IOVar.

       deallocate(bufHeader, IOVar, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writePlot3DZoneFile", &
                        "Deallocation failure for bufHeader and IOVar")

       ! Write a message that the file has been written.

       if (myID == 0) then
         print "(a)", "# Plot3D grid file written"
         print "( a)", "#"
       endif

       end subroutine writePlot3DZoneFile

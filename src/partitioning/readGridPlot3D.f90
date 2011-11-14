!
!      ******************************************************************
!      *                                                                *
!      * File:          readGridPlot3D.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-22-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readGridPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readGridPlot3D reads the coordinates for the blocks or block   *
!      * parts to be stored on this processor from the plot3D file(s).  *
!      * These file has already been opened in the subroutine           *
!      * readBlockSizesPlot3D.                                          *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use flowVarRefState
       use inputIO
       use IOModule
       use partitionMod
       use su_cgns
       implicit none
!
!      Local variables.
!
       integer :: ierr
       integer, dimension(mpi_status_size) :: status

       integer(kind=intRecordPLOT3DType) :: recordSize

       integer(kind=intType) :: nn, mm, sizeP3D_Record

       integer(kind=mpi_offset_kind) :: sizeHeader

       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the size of the integer record type to sizeP3D_Record.
       ! The reason is that the variable nBytesPerRecordIntPLOT3D is
       ! of type integer and not integer(kind=intType).

       sizeP3D_Record = nBytesPerRecordIntPLOT3D
!
!      ******************************************************************
!      *                                                                *
!      * The reading of the coordinates. A loop is performed over the   *
!      * desired number of grids to be read.                            *
!      *                                                                *
!      ******************************************************************
!
       nGridLoop: do nn=1,nGridsRead

         ! Easier storage of the file pointer and the variables needed
         ! to do the actual reading.

         fh = fileIDs(nn)

         P3D_ByteSwap    = byteSwapGrids(nn)
         P3D_BlockFormat = blockFormatGrids(nn)
         P3D_DataStorage = nodeData
         P3D_nVar        = 3

         ! Determine the size of the file header in bytes. This is used
         ! to put the file pointer at the correct location, i.e. where
         ! the record of the coordinates of the first block starts. 
         ! The size of the header depends on the type of PLOT3D format.

         sizeHeader = 3*cgnsNDom*nBytesPerIntPLOT3D &
                    + 2*nBytesPerRecordIntPLOT3D
         if(P3D_BlockFormat == P3D_MultiBlock) &
           sizeHeader = sizeHeader + nBytesPerIntPLOT3D &
                      + 2*nBytesPerRecordIntPLOT3D

         ! Jump to the position in the file just after the file header.

         call mpi_file_seek(fh, sizeHeader, mpi_seek_set, ierr)
!
!        ****************************************************************
!        *                                                              *
!        * Determine the precision of the coordinates and whether or    *
!        * not an iblanking array is stored.                            *
!        *                                                              *
!        ****************************************************************
!
         ! Read the leading record integer.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)
         if( P3D_ByteSwap ) &
           call byteswap(recordSize, sizeP3D_Record, 1_intType)

         ! Determine the situation. If the situation does not correspond
         ! to an expected one, print an error message and exit.

         mm = cgnsDoms(1)%il * cgnsDoms(1)%jl * cgnsDoms(1)%kl

         if(recordSize == 3*4*mm) then
           P3D_Precision = precisionSingle
           P3D_iblank    = .false.
         else if(recordSize == 3*8*mm) then
           P3D_Precision = precisionDouble
           P3D_iblank    = .false.
         else if(recordSize == (3*4+nBytesPerIntPLOT3D)*mm) then
           P3D_Precision = precisionSingle
           P3D_iblank    = .true.
         else if(recordSize == (3*8+nBytesPerIntPLOT3D)*mm) then
           P3D_Precision = precisionDouble
           P3D_iblank    = .true.
         else
           write(errorMessage,"(3a)") "Unknown plot3D format for grid &
                                      &file ", trim(gridFiles(nn)), "."
           if(myID == 0) &
             call terminate("readGridPlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! All the pointers to IOVar have already been set.
         ! So just call the general routines to perform the reading.

         P3D_Offset = sizeHeader

         call prepareReadPlot3D(.false.)
         call readPlot3DVar(nn, LRef, 1_intType)
         call releaseMemIOPlot3D

         ! Close the grid file.

         call mpi_file_close(fh, ierr)

       enddo nGridLoop

       ! If the coordinates in the solution files must be written in
       ! meters, correct this info for all cgns blocks.

       if( writeCoorMeter ) then
         do nn=1,cgnsNDom
           cgnsDoms(nn)%len = Meter
           cgnsDoms(nn)%gridUnitsSpecified = .true.
           cgnsDoms(nn)%LRef = one
         enddo

         LRef = one
       endif

       end subroutine readGridPlot3D

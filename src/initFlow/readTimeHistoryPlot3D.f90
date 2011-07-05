!
!      ******************************************************************
!      *                                                                *
!      * File:          readTimeHistoryPlot3D.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-26-2005                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTimeHistoryPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readTimeHistoryPlot3D attempts to read the time history of an  *
!      * unsteady computation from the given Plot3D restart file.       *
!      * If present it will be stored in the arrays timeArray and       *
!      * timeDataArray, for which memory is allocated.                  *
!      * It is assumed that the file pointer is already positioned at   *
!      * the correct location.                                          *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use inputIO
       use inputUnsteady
       use IOModule
       use monitor
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr, sizeRead
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record, sizeP3D_Real
       integer(kind=intType) :: i, j, ii, nConv

       integer(kind=intRecordPLOT3DType) :: recordSize
       integer(kind=mpi_offset_kind) :: disp

       integer(kind=intPLOT3DType), dimension(2) :: intPlot3DBuf

       integer(kind=intType), dimension(:), allocatable :: sortedNames2Or

       character(len=maxCGNSNameLen) :: message
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                               convNames, convNamesSorted

       real(kind=4), dimension(:), allocatable :: timeBuf4
       real(kind=8), dimension(:), allocatable :: timeBuf8

       logical :: historyPresent, allHistInfo
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
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

       ! Set the size of the floating point type depending on the
       ! precision used.

       select case (P3D_Precision)
         case (precisionSingle)
           sizeP3D_Real = 4
         case (precisionDouble)
           sizeP3D_Real = 8
       end select

       ! Read the size of the record and determine whether or not a
       ! convergence history is present.

       historyPresent = .true.

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)

       if(ierr == MPI_SUCCESS) then

         ! The reading was successfull. Check the record size.

         if( P3D_ByteSwap ) &
           call byteswap(recordSize, sizeP3D_Record, 1_intType)

         if(recordSize == (maxCGNSNameLen + 2*sizeP3D_Int)) then

           ! A history is present. Check if it is a convergence history.

           sizeRead = maxCGNSNameLen
           call mpi_file_read(fh, message, sizeRead, mpi_character, &
                              status, ierr)

           if(message /= TimeHistory) historyPresent = .false.

         else

           ! Don't know what this is, but it is not a time history.

           historyPresent = .false.
         endif

       else

         ! Reading was not okay.

         historyPresent = .false.

       endif

       ! Check if a time history is not present.

       checkHistoryPresent: if(.not. historyPresent) then

         ! No time history present. Set nTimeStepsRestart and
         ! timeUnsteadyRestart to zero, allocate the memory for the
         ! time history of the monitoring variables, print a warning
         ! and return.

         nTimeStepsRestart   = 0
         timeUnsteadyRestart = zero

         call allocTimeArrays(nTimeStepsFine)

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# No time history found in restart file."
         print "(a)", "# Starting at timestep 1 on the finest level."
         print "(a)", "#"

         return

       endif checkHistoryPresent

       ! Time history is present. Determine the number of old time
       ! levels and the number of monitoring variables stored.
       ! Read the terminating integer of the record as well.

       call mpi_file_read(fh, intPlot3DBuf, 2, sumb_integerPLOT3D, &
                          status, ierr)
       if( P3D_ByteSwap ) &
         call byteswap(intPlot3DBuf, sizeP3D_Int, 2_intType)

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)

       nConv             = intPlot3DBuf(1)
       nTimeStepsRestart = intPlot3DBuf(2)

       ! Allocate the memory for the names of the variables whose
       ! convergence histories are stored in the file as well as the
       ! mapping from sorted to original numbering.

       allocate(convNames(nConv), convNamesSorted(nConv), &
                sortedNames2Or(nConv), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("readTimeHistoryPlot3D", &
                        "Memory allocation failure for convNames, &
                        &convNamesSorted and sortedNames2Or")

       ! Allocate the memory for the buffer to store the convergence
       ! data. This depends on the precision used in the file.

       select case (P3D_Precision)
         case (precisionSingle)
           allocate(timeBuf4(nTimeStepsRestart), stat=ierr)

         case (precisionDouble)
           allocate(timeBuf8(nTimeStepsRestart), stat=ierr)
       end select

       if(ierr /= 0)                             &
         call terminate("readTimeHistoryPlot3D", &
                        "Memory allocation failure for either timeBuf4 &
                        &or timeBuf8")

       ! Read the names of the convergence variables.

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)
       if( P3D_ByteSwap ) &
         call byteswap(recordSize, sizeP3D_Record, 1_intType)

       sizeRead = nConv*maxCGNSNameLen
       if(recordSize /= sizeRead)                &
         call terminate("readTimeHistoryPlot3D", &
                        "Unexpected size for the convergence names")

       call mpi_file_read(fh, convNames, sizeRead, mpi_character, &
                          status, ierr)

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)

       ! Create a sorted version of the convergence names and determine
       ! the mapping from the sorted to the original names.

       do j=1,nConv
         convNamesSorted(j) = convNames(j)
       enddo

       call qsortStrings(convNamesSorted, nConv)

       do j=1,nConv
         ii = bsearchStrings(convNames(j), convNamesSorted, nConv)
         sortedNames2Or(ii) = j
       enddo

       ! Determine the total number of time levels and allocate the
       ! memory for the time history arrays.

       j = nTimeStepsRestart + nTimeStepsFine
       call allocTimeArrays(j)

       ! Read the old time values.

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)
       if( P3D_ByteSwap ) &
         call byteswap(recordSize, sizeP3D_Record, 1_intType)

       sizeRead = nTimeStepsRestart*sizeP3D_Real
       if(recordSize /= sizeRead)                &
         call terminate("readTimeHistoryPlot3D", &
                        "Unexpected size for the convergence names")

       select case (P3D_Precision)
         case (precisionSingle)

           sizeRead = nTimeStepsRestart
           call mpi_file_read(fh, timeBuf4, sizeRead, mpi_real4, &
                              status, ierr)
           if( P3D_ByteSwap ) &
             call byteswap(timeBuf4, sizeP3D_Real, nTimeStepsRestart)

           do i=1,nTimeStepsRestart
             timeArray = timeBuf4(i)
           enddo

         !===============================================================

         case (precisionDouble)

           sizeRead = nTimeStepsRestart
           call mpi_file_read(fh, timeBuf8, sizeRead, mpi_real8, &
                              status, ierr)
           if( P3D_ByteSwap ) &
             call byteswap(timeBuf8, sizeP3D_Real, nTimeStepsRestart)

           do i=1,nTimeStepsRestart
             timeArray = timeBuf8(i)
           enddo

       end select

       ! Set the value of timeUnsteadyRestart to the last value in
       ! timeArray.

       timeUnsteadyRestart = timeArray(nTimeStepsRestart)

       ! Initialize allHistInfo to .true. and perform the loop over
       ! the number of monitoring variables.

       allHistInfo = .true.

       monReadLoop: do j=1,nMon

         ! Search for the monitoring name in the sorted convergence
         ! names present in the restart file and check if the name
         ! is present.

         ii = bsearchStrings(monNames(j), convNamesSorted, nConv)

         testMonPresent: if(ii == 0) then

           ! Name not present in the restart file. Set allHistInfo
           ! to .false. and the corresponding entries in timeDataArray
           ! to zero

           allHistInfo = .false.
           do i=1,nTimeStepsRestart
             timeDataArray(i,j) = zero
           enddo

         else testMonPresent

           ! Name is present in the restart file. Determine the index
           ! in the original array and jump to the place in the file
           ! where the time history is stored.

           ii = sortedNames2Or(ii)

           disp = sizeHeader       + nVar*sizeVolumeSol             &
                + sizeConvHistory                                   &
                + 2*sizeP3D_Record + maxCGNSNameLen + 2*sizeP3D_Int &
                + 2*sizeP3D_Record + nConv*maxCGNSNameLen           &
                + sizeP3D_Record   + (ii-1)*nTimeStepsRestart*sizeP3D_Real

           call mpi_file_seek(fh, disp, mpi_seek_set, ierr)

           ! Read the convergence history, depending on the precision,
           ! and store it in the correct place of convArray.

           select case (P3D_Precision)
             case (precisionSingle)

               sizeRead = nTimeStepsRestart
               call mpi_file_read(fh, timeBuf4, sizeRead, mpi_real4, &
                                  status, ierr)
               if( P3D_ByteSwap ) &
                 call byteswap(timeBuf4, sizeP3D_Real, nTimeStepsRestart)

               do i=1,nTimeStepsRestart
                 timeDataArray(i,j) = timeBuf4(i)
               enddo

             !===========================================================

             case (precisionDouble)

               sizeRead = nTimeStepsRestart
               call mpi_file_read(fh, timeBuf8, sizeRead, mpi_real8, &
                                  status, ierr)
               if( P3D_ByteSwap ) &
                 call byteswap(timeBuf8, sizeP3D_Real, nTimeStepsRestart)

               do i=1,nTimeStepsRestart
                 timeDataArray(i,j) = timeBuf8(i)
               enddo

           end select

         endif testMonPresent

       enddo monReadLoop

       ! Print a warning in case not all the time history could
       ! be retrieved from the restart file.

       if(.not. allHistInfo) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# Not all the time history could be &
                      &retrieved from the restart file."
         print "(a)", "# Missing information is initialized to zero."
         print "(a)", "#"

       endif

       ! Release the memory of the locally allocated variables.

       deallocate(convNames, convNamesSorted, sortedNames2Or, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("readTimeHistoryPlot3D", &
                        "Deallocation failure for convNames, &
                        &convNamesSorted and sortedNames2Or")

       select case (P3D_Precision)
         case (precisionSingle)
           deallocate(timeBuf4, stat=ierr)

         case (precisionDouble)
           deallocate(timeBuf8, stat=ierr)
       end select

       if(ierr /= 0)                             &
         call terminate("readTimeHistoryPlot3D", &
                        "Deallocation failure for either timeBuf4 &
                        &or timeBuf8")

       end subroutine readTimeHistoryPlot3D

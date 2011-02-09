!
!      ******************************************************************
!      *                                                                *
!      * File:          readConvHistoryPlot3D.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-26-2005                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readConvHistoryPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readConvHistoryPlot3D attempts to read the convergence history *
!      * from the Plot3D restart file. If present it will be stored     *
!      * in the array convArray, for which memory is allocated.         *
!      * It is assumed that the file pointer is already positioned at   *
!      * the correct location.                                          *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use inputIO
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
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
       integer(kind=intType) :: i, j, ii, nn, sol, nConv

       integer(kind=intType), dimension(:), allocatable :: sortedNames2Or

       integer(kind=intRecordPLOT3DType) :: recordSize
       integer(kind=mpi_offset_kind) :: disp

       integer(kind=intPLOT3DType), dimension(2) :: intPlot3DBuf

       real(kind=4), dimension(:), allocatable :: convBuf4
       real(kind=8), dimension(:), allocatable :: convBuf8

       character(len=8)              :: intString
       character(len=maxCGNSNameLen) :: message
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                               convNames, convNamesSorted

       logical :: convPresent, allConvInfo
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
       ! Initialize sizeConvHistory to 0.

       sizeConvHistory = 0

       ! Check if we should return immediately for time spectral.

       if(equationMode == timeSpectral) then
         if(solID > 1 .and. nSolsRead /= nTimeIntervalsSpectral) return
       endif

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

       convPresent = .true.

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

           if(message /= ConvHistory) convPresent = .false.

         else

           ! Don't know what this is, but it is not a convergence
           ! history.

           convPresent = .false.
         endif

       else

         ! Reading was not okay.

         convPresent = .false.

       endif

       ! Store the total number of iterations to be performed in nn.
       ! For unsteady computations the entire convergence history
       ! will be stored, if desired.

       nn = nsgStartup + nCycles
       if(equationMode == unsteady) nn = nTimeStepsFine*nn

       ! Check if a convergence history is not present.

       checkConvPresent: if(.not. convPresent) then

         ! Convergence history is not present. Set nIterOld to 0 and
         ! allocate the memory for the convergence histories if this
         ! is the 1st solution to be stored. Otherwise do a test for
         ! consistency.

         if(solID == 1) then

           nIterOld = 0
           call allocConvArrays(nn)

           ! Print a warning if this is a steady or a time spectral
           ! computation. For an unsteady computation it could have
           ! have intention that the convergence history was not stored.

           if(equationMode == steady .or. &
              equationMode == timeSpectral) then

             print "(a)", "#"
             print "(a)", "#                 Warning"
             print "(a)", "# No convergence info found in restart file."
             print "(a)", "# Starting at iteration 0 on the finest &
                          &level."
             print "(a)", "#"

           endif

         else if(nIterOld /= 0) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Inconsistent convergence histories for the &
                        &time spectral mode."
           print "(a)", "# Starting at iteration 0 on the finest level."
           print "(a)", "#"

           nIterOld = 0

         endif

         ! Reset the file pointer, such that a possible time history
         ! will be read correctly.

         disp = sizeHeader + nVar*sizeVolumeSol
         call mpi_file_seek(fh, disp, mpi_seek_set, ierr)

         ! Return, because the convergence history is not present.

         return

       endif checkConvPresent

       ! Convergence history is present. Determine the number of old
       ! iterations and the number of monitoring variables stored.
       ! Read the terminating integer of the record as well.

       call mpi_file_read(fh, intPlot3DBuf, 2, sumb_integerPLOT3D, &
                          status, ierr)
       if( P3D_ByteSwap ) &
         call byteswap(intPlot3DBuf, sizeP3D_Int, 2_intType)

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)

       ! Determine the size of the entire convergence history.

       sizeConvHistory = 6*sizeP3D_Record + maxCGNSNameLen              &
                       + 2*sizeP3D_Int + intPlot3DBuf(1)*maxCGNSNameLen &
                       + intPlot3DBuf(1)*intPlot3DBuf(2)*sizeP3D_Real

       ! If the convergence history should not be stored, unsteady mode,
       ! determine the size of the remaining part of the convergence
       ! history and jump over it and return.

       if(.not. storeConvInnerIter) then

         disp = 4*sizeP3D_Record + intPlot3DBuf(1)*maxCGNSNameLen &
              + intPlot3DBuf(1)*intPlot3DBuf(2)*sizeP3D_Real

         call mpi_file_seek(fh, disp, mpi_seek_cur, ierr)
         return

       endif

       ! Check for time spectral mode, solID > 1, if the number of
       ! iterations is consistent. If not, set nIterOld to zero, and
       ! return. No need to allocate the memory, because that has
       ! already been done.

       if(solID > 1 .and. nIterOld /= (intPlot3DBuf(2) - 1)) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# Inconsistent convergence histories for the &
                      &time spectral mode."
         print "(a)", "# Starting at iteration 0 on the finest level."
         print "(a)", "#"

         nIterOld = 0
         return

       endif

       ! Store the number of convergence histories and iterations. Note
       ! that 1 is substracted from nIterOld, because also the reference
       ! L2 norm is stored.

       nConv    = intPlot3DBuf(1)
       nIterOld = intPlot3DBuf(2) - 1

       ! Allocate the memory for the convergence arrays if this is
       ! the first solution. The memory for all other spectral modes are
       ! allocated here as well, which explains the test.

       nn = nn + nIterOld
       call allocConvArrays(nn)

       ! Allocate the memory for the names of the variables whose
       ! convergence histories are stored in the file as well as the
       ! mapping from sorted to original numbering.

       allocate(convNames(nConv), convNamesSorted(nConv), &
                sortedNames2Or(nConv), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("readConvHistoryPlot3D", &
                        "Memory allocation failure for convNames, &
                        &convNamesSorted and sortedNames2Or")

       ! Allocate the memory for the buffer to store the convergence
       ! data. This depends on the precision used in the file.

       select case (P3D_Precision)
         case (precisionSingle)
           allocate(convBuf4(0:nIterOld), stat=ierr)

         case (precisionDouble)
           allocate(convBuf8(0:nIterOld), stat=ierr)
       end select

       if(ierr /= 0)                             &
         call terminate("readConvHistoryPlot3D", &
                        "Memory allocation failure for either convBuf4 &
                        &or convBuf8")

       ! Read the names of the convergence variables.

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)
       if( P3D_ByteSwap ) &
         call byteswap(recordSize, sizeP3D_Record, 1_intType)

       sizeRead = nConv*maxCGNSNameLen
       if(recordSize /= sizeRead)                &
         call terminate("readConvHistoryPlot3D", &
                        "Unexpected size for the convergence names")

       call mpi_file_read(fh, convNames, sizeRead, mpi_character, &
                          status, ierr)

       call mpi_file_read(fh, recordSize, 1, sumb_integerRecordPLOT3D, &
                          status, ierr)

       ! Create a sorted version of the convergence names and determine
       ! the mapping from the sorted to the original names.

       do nn=1,nConv
         convNamesSorted(nn) = convNames(nn)
       enddo

       call qsortStrings(convNamesSorted, nConv)

       do nn=1,nConv
         ii = bsearchStrings(convNames(nn), convNamesSorted, nConv)
         sortedNames2Or(ii) = nn
       enddo

       ! Check for the L2 norm of the density residual. This must be
       ! present, otherwise the reading of the convergence file does
       ! not make sense and the rest of the convergence info will be
       ! ignored.

       ii = bsearchStrings(cgnsL2resRho, convNamesSorted, nConv)

       testDensityRes: if(ii == 0 .and. solID == 1) then

         ! Density residual not present for the 1st solution. This means
         ! that reading the convergence info is completely pointless.
         ! Print a warning for the steady mode and the time spectral mode.

         if(equationMode == steady .or. &
            equationMode == timeSpectral) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# No convergence info for the density &
                        &residual found in the restart file."
           print "(a)", "# The rest of the convergence history is &
                        &ignored."
           print "(a)", "# Starting at iteration 0 on the &
                        &finest level."
           print "(a)", "#"

         endif

         ! Set nIterOld to 0 and jump to the place where the memory of
         ! the local variables is released.

         nIterOld = 0
         goto 99

       endif testDensityRes

       ! Density residual is present. Print the restarting info
       ! if this is solution 1 and this is either a steady or a time
       ! spectral computation.

       if(solID == 1 .and. (equationMode == steady .or. &
                            equationMode == timeSpectral)) then

         write(intString,"(i7)") nIterOld
         intString = adjustl(intString)

         print "(a)", "#"
         print 120, trim(intString)
         print "(a)", "#"
 120     format("# Restarting at iteration",1X,A,".")
       endif

       ! Initialize allConvInfo to .true.

       allConvInfo = .true.

       ! Loop over all the variables that must be monitored.

       monReadLoop: do j=1,nMon

         ! Search for the monitoring name in the sorted convergence
         ! names present in the restart file and check if the name
         ! is present.

         ii = bsearchStrings(monNames(j), convNamesSorted, nConv)

         testMonPresent: if(ii == 0) then

           ! Name not present in the restart file. Set allConvInfo
           ! to .false. and the corresponding entries in convArray
           ! either to zero or to the values of solution 1.

           allConvInfo = .false.

           if(solID == 1) then
             do i=0,nIterOld
               convArray(i,solID,j) = zero
             enddo
           else
             do i=0,nIterOld
               convArray(i,solID,j) = convArray(i,1,j)
             enddo
           endif

         else testMonPresent

           ! Name is present in the restart file. Determine the index
           ! in the original array and jump to the place in the file
           ! where the convergence is stored.

           ii = sortedNames2Or(ii)

           disp = sizeHeader       + nVar*sizeVolumeSol             &
                + 2*sizeP3D_Record + maxCGNSNameLen + 2*sizeP3D_Int &
                + 2*sizeP3D_Record + nConv*maxCGNSNameLen           &
                + sizeP3D_Record   + (ii-1)*(nIterOld+1)*sizeP3D_Real

           call mpi_file_seek(fh, disp, mpi_seek_set, ierr)

           ! Read the convergence history, depending on the precision,
           ! and store it in the correct place of convArray.

           select case (P3D_Precision)
             case (precisionSingle)

               nn       = nIterOld + 1
               sizeRead = nn
               call mpi_file_read(fh, convBuf4, sizeRead, mpi_real4, &
                                  status, ierr)
               if( P3D_ByteSwap ) &
                 call byteswap(convBuf4, sizeP3D_Real, nn)

               do i=0,nIterOld
                 convArray(i,solID,j) = convBuf4(i)
               enddo

             !===========================================================

             case (precisionDouble)

               nn       = nIterOld + 1
               sizeRead = nn
               call mpi_file_read(fh, convBuf8, sizeRead, mpi_real8, &
                                  status, ierr)
               if( P3D_ByteSwap ) &
                 call byteswap(convBuf8, sizeP3D_Real, nn)

               do i=0,nIterOld
                 convArray(i,solID,j) = convBuf8(i)
               enddo
           end select

         endif testMonPresent

       enddo monReadLoop

       ! Print a warning in case not all the convergence info could
       ! be retrieved from the restart file. The error message depends
       ! on the situation.

       checkAllConvInfo: if(.not. allConvInfo) then

         if(equationMode == timeSpectral .and. solID > 1) then

           ! Not the first spectral solution in time spectral mode.
           ! Inform that the data which is not present is set to the
           ! data of the 1st spectral solution.

           write(intString,"(i7)") solID
           intString = adjustl(intString)

           print "(a)", "#"
           print 130
           print 140, trim(intString)
           print 160
           print "(a)", "#"

         else

           ! 1st solution. Missing information was set to zero.

           print "(a)", "#"
           print 130
           print 150
           print 170
           print "(a)", "#"

         endif

       endif checkAllConvInfo

       ! Format statements for the missing convergence info.

 130   format("#                 Warning")
 140   format("# Not all the convergence info could be retrieved &
              &from the restart file for solution",1X,A,".")
 150   format("# Not all the convergence info could be retrieved &
              &from the restart file.")
 160   format("# Missing information is copied from the first &
              &solution.")
 170   format("# Missing information is initialized to zero.")


       ! Continue statement needed when no convergence information
       ! of the density residual was present.

 99    continue

       ! Release the memory of the locally allocated variables.

       deallocate(convNames, convNamesSorted, sortedNames2Or, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("readConvHistoryPlot3D", &
                        "Deallocation failure for convNames, &
                        &convNamesSorted and sortedNames2Or")

       select case (P3D_Precision)
         case (precisionSingle)
           deallocate(convBuf4, stat=ierr)

         case (precisionDouble)
           deallocate(convBuf8, stat=ierr)
       end select

       if(ierr /= 0)                             &
         call terminate("readConvHistoryPlot3D", &
                        "Deallocation failure for either convBuf4 &
                        &or convBuf8")

       ! Check for he spectral mode if the number of convergence
       ! histories equals the number of spectral solutions. If not
       ! print a warning and initialize them.

       if(equationMode == timeSpectral) then
         if(nSolsRead /= nTimeIntervalsSpectral) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Inconsistent number of convergence histories &
                        &found for the time spectral method."
           print "(a)", "# Copying the data of the first solution to &
                &all spectral solutions."
           print "(a)", "#"

           do sol=2,nTimeIntervalsSpectral
             do j=1,nMon
               do i=0,nIterOld
                 convArray(i,sol,j) = convArray(i,1,j)
               enddo
             enddo
           enddo

         endif
       endif

       end subroutine readConvHistoryPlot3D

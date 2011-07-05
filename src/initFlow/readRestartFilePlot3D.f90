!
!      ******************************************************************
!      *                                                                *
!      * File:          readRestartFilePlot3D.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-30-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readRestartFilePlot3D(halosRead)
!
!      ******************************************************************
!      *                                                                *
!      * readRestartFilePlot3D reads a fine grid solution from a        *
!      * restart file previously created by this solver when using the  *
!      * PLOT3D grid format. The solution file is not in PLOT3D format, *
!      * because it is too restrictive. Instead an internally used      *
!      * unformatted format is used.                                    *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use flowVarRefState
       use inputIO
       use inputPhysics
       use IOModule
       use monitor
       use restartMod
       implicit none
!
!      Subroutine arguments
!
       logical, intent(out) :: halosRead
!
!      Local variables.
!
       integer :: ierr, sizeRead

       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: ii, nn, mm
       integer(kind=intType) :: sizeP3D_Real, sizeP3D_Int, sizeP3D_Record

       integer(kind=intRecordPLOT3DType) :: recordSize
       integer(kind=intPLOT3DType), dimension(2) :: intPlot3DBuf

       integer(kind=intPLOT3DType), dimension(:,:), allocatable :: &
                                                                blockDim

       integer(kind=mpi_offset_kind) :: disp

       real(kind=4), dimension(3) :: refState4
       real(kind=8), dimension(3) :: refState8

       character(len=7)            :: integerString
       character(len=maxStringLen) :: errorMessage
       character(len=maxCGNSNameLen), dimension(:), allocatable :: namesVar

       logical :: includeHalos
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
       ! Initialize halosRead to .true. This will be overwritten if
       ! halo data cannot be read.

       halosRead = .true.

       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Set the iblanking to .false., because this can only be present
       ! in a grid file. Also indicate that 1 variable at a time is
       ! read, because this is a solution file.

       P3D_iblank = .false.
       P3D_nVar   = 1
!
!      ******************************************************************
!      *                                                                *
!      * Loop over the number of files to be read and perform the       *
!      * initialization tasks for the open file.                        *
!      *                                                                *
!      ******************************************************************
!
       fileLoop: do solID=1,nSolsRead

         ! Open the file for reading. Check if it went okay.

         call mpi_file_open(SUmb_comm_world, solFiles(solID), &
                           mpi_mode_rdonly, mpi_info_null,   &
                           fh, ierr)

         if(ierr /= mpi_success) then
           write(errorMessage,"(3a)") "Solution file ",      &
                                      trim(solFiles(solID)), &
                                      " not found"
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Read the size of the first Fortran record and determine
         ! whether or not byte swapping must be used for this file.
         ! All the processors read it, because that is the most efficient
         ! way of doing it. A collective read is nonsense here.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)

         if(recordSize ==   sizeP3D_Int .or. &
            recordSize == 2*sizeP3D_Int) then
           P3D_ByteSwap = .false.
         else
           call byteswap(recordSize, sizeP3D_Record, 1_intType)
           P3D_ByteSwap = .true.
         endif

         ! Determine what type of block format is used and read the
         ! corresponding info. Again all processors do the reading.

         if(recordSize == sizeP3D_Int) then
           P3D_BlockFormat = P3D_SingleBlock

           intPlot3DBuf(1) = 1
           call mpi_file_read(fh, intPlot3DBuf(2), 1, &
                              sumb_integerPLOT3D, status, ierr)
           if( P3D_ByteSwap ) &
             call byteswap(intPlot3DBuf(2), sizeP3D_Int, 1_intType)

         else if(recordSize == 2*sizeP3D_Int) then
           P3D_BlockFormat = P3D_MultiBlock

           call mpi_file_read(fh, intPlot3DBuf, 2, &
                              sumb_integerPLOT3D, status, ierr)
           if( P3D_ByteSwap ) &
             call byteswap(intPlot3DBuf, sizeP3D_Int, 2_intType)

         else
           write(errorMessage,"(3a)") "Unknown solution format for &
                                      &file ", &
                                      trim(solFiles(solID)), "."
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Check the number of blocks and read the closing integer
         ! of the record.

         if(intPlot3DBuf(1) /= cgnsNDom) then
           write(errorMessage,"(3a)") "Number of blocks in solution &
                                      &file ", trim(solFiles(solID)), &
                                      " differs from the number in &
                                      &the grid file."
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         nVar = intPlot3DBuf(2)

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)
!
!        ****************************************************************
!        *                                                              *
!        * Read the reference state used in the previous computation    *
!        * and determine the scale factors accordingly. These scale     *
!        * factors will only differ from 1.0 if the reference state is  *
!        * different for the previous and the current computation.      *
!        *                                                              *
!        ****************************************************************
!
         ! Read the opening integer of the record and determine the
         ! precision of the floating point type.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)
         if( P3D_ByteSwap ) &
           call byteswap(recordSize, sizeP3D_Record, 1_intType)

         if(recordSize == 4*3) then
           P3D_Precision = precisionSingle
           sizeP3D_Real = 4
         else if(recordSize == 8*3) then
           P3D_Precision = precisionDouble
           sizeP3D_Real = 8
         else
           write(errorMessage,"(3a)") "Unexpected record size for the &
                                      &reference state in file ", &
                                      trim(solFiles(solID)), "."
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Read the reference state, depending on the precision of the
         ! solution file, and store the values for the moment in rhoScale
         ! and pScale. The reference temperature is not important for a
         ! consistent restart.

         select case (P3D_Precision)
           case (precisionSingle)

             call mpi_file_read(fh, refState4, 3, mpi_real4, status, ierr)
             if( P3D_ByteSwap ) &
               call byteswap(refState4, sizeP3D_Real, 3_intType)

             rhoScale = refState4(1)
             pScale   = refState4(2)

           !=============================================================

           case (precisionDouble)

             call mpi_file_read(fh, refState8, 3, mpi_real8, status, ierr)
             if( P3D_ByteSwap ) &
               call byteswap(refState8, sizeP3D_Real, 3_intType)

             rhoScale = refState8(1)
             pScale   = refState8(2)

         end select

         ! Read the closing integer of this record.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)

         ! Determine the scale factors for this restart.

         rhoScale = rhoScale/rhoRef
         pScale   = pScale/pRef
         velScale = sqrt(pScale/rhoScale)
         muScale  = sqrt(pScale*rhoScale)
!
!        ****************************************************************
!        *                                                              *
!        * Determine the variable names and the block dimensions.       *
!        *                                                              *
!        ****************************************************************
!
         ! Allocate the memory the variable names, its sorted version
         ! and the mapping from the sorted to the original names.

         allocate(namesVar(nVar), varNames(nVar), sorted2Or(nVar), &
                  stat=ierr)
         if(ierr /= 0)                             &
           call terminate("readRestartFilePlot3D", &
                          "Memory allocation failure for namesVar, etc.")

         ! Read the record size, check it and read the variable names.
         ! For the latter it is not needed to do the byte swapping.
         ! Again all the processors do the reading.
         ! Also read the closing integer of the record.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)
         if( P3D_ByteSwap ) &
           call byteswap(recordSize, sizeP3D_Record, 1_intType)

         sizeRead = nVar*maxCGNSNameLen
         if(recordSize /= sizeRead) then
           write(errorMessage,"(3a)") "Unexpected record size for the &
                                      &variable names in solution &
                                      &file ", &
                                      trim(solFiles(solID)), "."
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         call mpi_file_read(fh, namesVar, sizeRead, mpi_character, &
                           status, ierr)

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)

         ! Copy namesVar in varNames and sort the latter.

         do nn=1,nVar
           varNames(nn)  = namesVar(nn)
           sorted2Or(nn) = 0
         enddo

         mm = nVar
         call qsortStrings(varNames, mm)

         ! Create the mapping from the sorted to the original names.

         do nn=1,nVar
           ii = bsearchStrings(namesVar(nn), varNames, mm)
           if(sorted2Or(ii) /= 0) then
             write(errorMessage,"(5a)") "Variable ",        &
                                        trim(namesVar(nn)), &
                                        "is stored more than once in &
                                        &solution file ",   &
                                        trim(solFiles(solID)), "."
             if(myID == 0) &
               call terminate("readRestartFilePlot3D", errorMessage)
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           sorted2Or(ii) = nn
         enddo

         ! Allocate the memory for dimensions of the blocks.

         allocate(blockDim(3,cgnsNDom), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("readRestartFilePlot3D", &
                          "Memory allocation failure for blockDim")

         ! Check the size of the record and read the dimensions of
         ! the blocks. Again done by all processors.

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)
         if( P3D_ByteSwap ) &
           call byteswap(recordSize, sizeP3D_Record, 1_intType)

         sizeRead = 3*cgnsNDom
         if(recordSize /= sizeRead*sizeP3D_Int)  then
           write(errorMessage,"(3a)") "Unexpected record size for the &
                                      &block dimensions in solution &
                                      &file ", &
                                      trim(solFiles(solID)), "."
           if(myID == 0) &
             call terminate("readRestartFilePlot3D", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         call mpi_file_read(fh, blockDim, sizeRead, &
                            sumb_integerPLOT3D, status, ierr)
         if( P3D_ByteSwap ) then
            nn = 3*cgnsNDom
           call byteswap(blockDim, sizeP3D_Int, nn)
         endif

         call mpi_file_read(fh, recordSize, 1, &
                            sumb_integerRecordPLOT3D, status, ierr)

         ! Check if halo values are stored and if the block dimensions
         ! are okay. Remember that this solution file stores the 
         ! cell-centered values. Set the variable P3D_DataStorage
         ! accordingly. If this is not an older state for an unsteady
         ! computation correct the value of halosRead if no halo cells
         ! are present. Also set the value of includeHalos, which is
         ! an additional argument needed when halos are present in the
         ! file, but should not be read.

         if(blockDim(1,1) == cgnsDoms(1)%nx) then
           ii = -1;
           P3D_DataStorage = cellDataNoHalo
           includeHalos    = .false.
           if(solID == 1 .or. equationMode == timeSpectral) &
             halosRead = .false.
         else
           ii = 1;
           P3D_DataStorage = cellDataPlusHalo
           includeHalos    = .false.
           if(solID == 1 .or. equationMode == timeSpectral) &
             includeHalos = .true.
         endif

         do nn=1,cgnsNDom
           if(cgnsDoms(nn)%il /= blockDim(1,nn) - ii .or. &
              cgnsDoms(nn)%jl /= blockDim(2,nn) - ii .or. &
              cgnsDoms(nn)%kl /= blockDim(3,nn) - ii) then
             write(errorMessage,"(3a)") "Wrong block dimensions in &
                                        &solution file ", &
                                        trim(solFiles(solID)), "."
             if(myID == 0) &
               call terminate("readRestartFilePlot3D", errorMessage)
             call mpi_barrier(SUmb_comm_world, ierr)
           endif
         enddo

         ! Determine the size of the header part of the file.
         ! This depends on the format of the solution file.

         sizeHeader = 2*sizeP3D_Record + sizeP3D_Int            &
                    + 2*sizeP3D_Record + 3*sizeP3D_Real         &
                    + 2*sizeP3D_Record + nVar*maxCGNSNameLen    &
                    + 2*sizeP3D_Record + 3*cgnsNDom*sizeP3D_Int

         if(P3D_BlockFormat == P3D_MultiBlock) &
           sizeHeader = sizeHeader + sizeP3D_Int

         ! Determine the size of one volume solution, i.e. the
         ! number of bytes to store one variable.

         sizeVolumeSol = 0

         do nn=1,cgnsNDom
           sizeVolumeSol = sizeVolumeSol + 2*sizeP3D_Record &
                         + blockDim(1,nn) * blockDim(2,nn)  &
                         * blockDim(3,nn) * sizeP3D_Real
         enddo

         ! Release the memory of blockDim.

         deallocate(blockDim, stat=ierr)
         if(ierr /= 0)                             &
           call terminate("readRestartFilePlot3D", &
                          "Deallocation failure for blockDim")
!
!        ****************************************************************
!        *                                                              *
!        * Reading of the solution variables.                           *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the necessary variables for the reading.

         call prepareReadPlot3D(includeHalos)

         ! Read the density and the turbulence variables.

         call readDensityPlot3D
         call readTurbvarPlot3D

         ! Read the other variables, depending on the situation.

         testPrim: if(solID == 1 .or. equationMode == timeSpectral) then

           ! Either the first solution or time spectral mode. Read
           ! the primitive variables from the restart file.

           call readXvelocityPlot3D
           call readYvelocityPlot3D
           call readZvelocityPlot3D
           call readPressurePlot3D(includeHalos)

         else testPrim

           ! Old solution in unsteady mode. Read the conservative
           ! variables.

           call readXmomentumPlot3D
           call readYmomentumPlot3D
           call readZmomentumPlot3D
           call readEnergyPlot3D

         endif testPrim

         ! Release the memory of the help variables needed to
         ! perform the reading.

         call releaseMemIOPlot3D

         ! Release the memory of namesVar, its sored version and the
         ! mapping between the two.

         deallocate(namesVar, varNames, sorted2Or, stat=ierr)
         if(ierr /= 0)                             &
           call terminate("readRestartFilePlot3D", &
                          "Deallocation failure for namesVar, etc.")
!
!        ****************************************************************
!        *                                                              *
!        * Reading of the convergence histories and possibly time       *
!        * histories. Only processor 0 performs the reading.            *
!        *                                                              *
!        ****************************************************************
!
         ! Set the file pointer at the place where the convergence
         ! histories should be stored.

         disp = sizeHeader + nVar*sizeVolumeSol
         call mpi_file_seek(fh, disp, mpi_seek_set, ierr)

         ! Read the convergence histories and possibly the time
         ! histories. Only done by processor 0. For unsteady mode these
         ! histories are only read for the current solution, i.e.
         ! solID == 1.

         if(myID == 0) then

           select case (equationMode)
             case (steady, timeSpectral)
               call readConvHistoryPlot3D

             case (unsteady)
               if(solID == 1) then
                 call readConvHistoryPlot3D
                 call readTimeHistoryPlot3D
               endif
           end select

         endif

         ! Close the file.

         call mpi_file_close(fh, ierr)

       enddo fileLoop

       ! Broadcast nTimeStepsRestart and timeUnsteadyRestart to all
       ! processors. These values are needed to perform a consistent
       ! unsteady restart.

       call mpi_bcast(nTimeStepsRestart, 1, sumb_integer, 0, &
                      SUmb_comm_world, ierr)
       call mpi_bcast(timeUnsteadyRestart, 1, sumb_real, 0, &
                      SUmb_comm_world, ierr)

       ! Write a message about the time step number for which is
       ! restarted.

       if(equationMode == unsteady .and. myID == 0) then
         write(integerString,"(i7)") nTimeStepsRestart+1
         integerString = adjustl(integerString)

         print "(a)", "#"
         print 110, trim(integerString)
         print "(a)", "#"
 110     format("# Restarting at time step",1X,A,".")
       endif

       end subroutine readRestartFilePlot3D

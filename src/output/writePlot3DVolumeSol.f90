!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DVolumeSol.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-21-2005                                      *
!      * Last modified: 10-31-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DVolumeSol
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DVolumeSol writes the solution file(s) when the      *
!      * plot3D format is used for the grids. The solution file written *
!      * in this routine is not a plot3D solution file. It is just an   *
!      * internally used format.                                        *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use cgnsNames
       use communication
       use flowVarRefState
       use inputIO
       use inputPhysics
       use IOModule
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr
       integer :: sizeFirstRecord, sizeHeader
       integer :: amode

       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: ii, jj, mm, nn, ll, var
       integer(kind=intType) :: nVolSolvar, nVolDiscrVar, nVarWritten
       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record, sizeP3D_Real

       integer(kind=intType), dimension(nDom) :: iBeg, jBeg, kBeg
       integer(kind=intType), dimension(nDom) :: iEnd, jEnd, kEnd

       integer(kind=intType), dimension(cgnsNDom) :: cgnsBlockIDs

       integer(kind=intRecordPLOT3DType) :: sizeRecord
       integer(kind=intPLOT3DType), dimension(3) :: intPlot3DBuf

       integer(kind=mpi_offset_kind) :: sizeVolumeSol, disp

       real(kind=realType) :: dummyBuf

       real(kind=4), dimension(3) :: buf4
       real(kind=8), dimension(3) :: buf8

       character(len=maxStringLen) :: errorMessage

       character, dimension(:), allocatable :: writeBuf
       character(len=maxCGNSNameLen), &
                                   dimension(:), allocatable :: solNames

       logical :: unsteadyHigherSol, writeRindLayer
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number and names of the solution files.
       ! Also set the pointers for IOVar needed for the general
       ! treatment of the IO.

       call volSolFileNamesWrite

       ! Return immediately if no solution files have to be written.

       if(nVolSolToWrite == 0) return

       ! All blocks will be written. Set cgnsBlockIDs accordingly.

       do ii=1,cgnsNDom
         cgnsBlockIDs(ii) = ii
       enddo

       ! Set the variables needed for the writing of the grid files.
       ! These are the same for all the grid files to be written.

       P3D_Precision = precisionSol
       P3D_iblank    = .false.
       P3D_nVar      = 1

       ! Set the size of sizeP3D_Real, depending on P3D_Precision.

       select case (P3D_Precision)
         case (precisionSingle)
           sizeP3D_Real = 4
         case (precisionDouble)
           sizeP3D_Real = 8
       end select

       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Write a message that solution file(s) is (are) written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "# Writing volume solution file(s) ..."
       endif
!
!      ******************************************************************
!      *                                                                *
!      * A lot of initializations.                                      *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of variables to be written to the volume
       ! solution file as well as the cgns names. In unsteady mode these
       ! names are overwritten for the older solutions.

       call numberOfVolSolVariables(nVolSolvar, nVolDiscrVar)
       allocate(solNames(nVolSolvar+nVolDiscrVar), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writePlot3DVolumeSol", &
                        "Memory allocation failure for solNames")
       call volSolNames(solNames)

       ! Determine the size of the section with the reference
       ! state, block dimensions and variable names.

       nn = 3*sizeP3D_Real         + 2*sizeP3D_Record  &
          + 3*cgnsNDom*sizeP3D_Int + 2*sizeP3D_Record  &
          + (nVolSolvar + nVolDiscrVar)*maxCGNSNameLen &
          + 2*sizeP3D_Record

       ! Allocate the memory for writeBuf.

       allocate(writeBuf(nn), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writePlot3DVolumeSol", &
                        "Memory allocation failure for writeBuf.")
!
!      ******************************************************************
!      *                                                                *
!      * Loop over the number of solution files to be written and write *
!      * the solutions.                                                 *
!      *                                                                *
!      ******************************************************************
!
       nVarWritten       = nVolSolvar + nVolDiscrVar
       unsteadyHigherSol = .false.
       writeRindLayer    = storeRindLayer

       solLoop: do nn=1,nVolSolToWrite

         ! Open the solution file for writing. There is no guarantee that
         ! an old file will be clobbered, so wipe out any previous output
         ! file. If the file cannot be opened, print an error message
         ! and exit.

         if( myID == 0) &
           call mpi_file_delete(volSolFileNames(nn), mpi_info_null, ierr)

         amode = mpi_mode_create + mpi_mode_wronly
         call mpi_file_open(SUmb_comm_world, volSolFileNames(nn), &
                            amode, mpi_info_null, fh, ierr)

         if(ierr /= mpi_success) then
           write(errorMessage,*) "File ", trim(volSolFileNames(nn)), &
                                 " could not be opened for writing"
           if(myID == 0) &
             call terminate("writePlot3DVolumeSol", errorMessage)

           call mpi_barrier(SUmb_comm_world, ierr)
         endif
!
!        ****************************************************************
!        *                                                              *
!        * File header.                                                 *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the case we are having here.

         select case (P3D_BlockFormat)
           case (P3D_SingleBlock)

             ! Single block format, only the number of variables need
             ! to be written. First the opening integer of the record.

             ii = 1
             sizeRecord = sizeP3D_Int
             if( P3D_ByteSwap ) &
               call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
             call writeVarsToBuffer(writeBuf(ii), sizeRecord, &
                                    sizeP3D_Record)
             ii = ii + sizeP3D_Record

             ! The actual number of variables.

             intPlot3DBuf(1) = nVarWritten
             if( P3D_ByteSwap ) &
               call byteswap(intPlot3DBuf, sizeP3D_Int, 1_intType)

             call writeVarsToBuffer(writeBuf(ii), intPlot3DBuf, &
                                    sizeP3D_Int)
             ii = ii + sizeP3D_Int

           !=============================================================

           case (P3D_MultiBlock)

             ! Multiblock format. Both the number of blocks and the
             ! number of variables must be written. First the opening
             ! integer of the record.

             ii = 1
             sizeRecord = 2*sizeP3D_Int
             if( P3D_ByteSwap ) &
               call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
             call writeVarsToBuffer(writeBuf(ii), sizeRecord, &
                                    sizeP3D_Record)
             ii = ii + sizeP3D_Record

             ! The actual number of blocks and the number of variables. 

             intPlot3DBuf(1) = cgnsNDom
             intPlot3DBuf(2) = nVarWritten
             if( P3D_ByteSwap ) &
               call byteswap(intPlot3DBuf, sizeP3D_Int, 2_intType)

             jj = 2*sizeP3D_Int
             call writeVarsToBuffer(writeBuf(ii), intPlot3DBuf, jj)
             ii = ii + jj

         end select

         ! The closing integer of the record.

         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! Only processor 0 does the writing of this record.

         sizeFirstRecord = ii - 1
         if(myID == 0) &
           call mpi_file_write(fh, writeBuf, sizeFirstRecord, &
                               mpi_character, status, ierr)
!
!        ****************************************************************
!        *                                                              *
!        * Write the reference state, the names of the variables and    *
!        * the cell dimensions of every block.                          *
!        *                                                              *
!        ****************************************************************
!
         ! The size of the record of the reference state.

         ii = 1
         sizeRecord = 3*sizeP3D_Real
         if( P3D_ByteSwap ) &
           call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! Write the reference state to writeBuf. Make a distinction
         ! between single and double precision.

         jj = 3*sizeP3D_Real

         select case (P3D_Precision)
           case (precisionSingle)

             ! Single precision. Write the reference density, pressure
             ! and temperature to writeBuf. Do not update the counter ii
             ! yet. This is done after the possible byte swapping.

             buf4(1) = rhoRef
             buf4(2) = pRef
             buf4(3) = TRef

             call writeVarsToBuffer(writeBuf(ii), buf4, jj)

           !=============================================================

           case (precisionDouble)

             ! Double precision. Write the reference density, pressure
             ! and temperature to writeBuf. Do not update the counter ii
             ! yet. This is done after the possible byte swapping.

             buf8(1) = rhoRef
             buf8(2) = pRef
             buf8(3) = TRef

             call writeVarsToBuffer(writeBuf(ii), buf8, jj)
           
         end select

         ! Apply the byte swapping, if needed, and update the counter
         ! ii afterwards.

         if( P3D_ByteSwap ) &
           call byteswap(writeBuf(ii), sizeP3D_Real, 3_intType)
         ii = ii + jj

         ! The closing integer of this record.

         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! The size of the record of the variable names.

         sizeRecord = maxCGNSNameLen*nVarWritten
         if( P3D_ByteSwap ) &
           call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! The names of the variables. Note that for characters
         ! no byte swapping needs to be applied.

         jj = maxCGNSNameLen
         do mm=1,nVarWritten
           call writeVarsToBuffer(writeBuf(ii), solNames(mm), jj)
           ii = ii + jj
         enddo

         ! The closing integer of the record.

         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! Write the cell dimensions of every block; include the halo
         ! values if desired. First the opening integer of the record.

         sizeRecord = 3*cgnsNDom*sizeP3D_Int
         if( P3D_ByteSwap ) &
           call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! Loop over the number of blocks and store the sizes in the
         ! write buffer. The possible byte swapping is handled
         ! afterwards.

         jj = ii

         if( writeRindLayer ) then

           ! Halo layers must be written. Useful only for postprocessing.

           P3D_DataStorage = cellDataPlusHalo

           do mm=1,cgnsNDom
             intPlot3DBuf(1) = cgnsDoms(mm)%il + 1
             intPlot3DBuf(2) = cgnsDoms(mm)%jl + 1
             intPlot3DBuf(3) = cgnsDoms(mm)%kl + 1

             call writeVarsToBuffer(writeBuf(ii), intPlot3DBuf, &
                                    3*sizeP3D_Int)
             ii = ii + 3*sizeP3D_Int
           enddo

         else

           ! Only the internal cell centered values are written.

           P3D_DataStorage = cellDataNoHalo

           do mm=1,cgnsNDom
             intPlot3DBuf(1) = cgnsDoms(mm)%il - 1
             intPlot3DBuf(2) = cgnsDoms(mm)%jl - 1
             intPlot3DBuf(3) = cgnsDoms(mm)%kl - 1

             call writeVarsToBuffer(writeBuf(ii), intPlot3DBuf, &
                                    3*sizeP3D_Int)
             ii = ii + 3*sizeP3D_Int
           enddo

         endif

         ! Apply the byte swapping, if needed.

         mm = 3*cgnsNDom
         if( P3D_ByteSwap ) &
           call byteswap(writeBuf(jj), sizeP3D_Int, mm)

         ! Store the closing integer of the record in writeBuf.

         call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
         ii = ii + sizeP3D_Record

         ! Write the records for the reference state, variable names and
         ! block dimensions; only processor 0 performs the writing.

         sizeHeader = ii - 1
         if(myID == 0)                                  &
           call mpi_file_write(fh, writeBuf, sizeHeader, &
                               mpi_character, status, ierr)
!
!        ****************************************************************
!        *                                                              *
!        * Some additional initializations, which depend on the         *
!        * settings of the current solution file.                       *
!        *                                                              *
!        ****************************************************************
! 
         ! Determine the cell index range to be stored for every block.
         ! This depends whether or not a halo layer must be written
         ! at all and if the current computational block contains
         ! boundaries of the original block.

         do mm=1,nDom
           iBeg(mm) = 2; iEnd(mm) = flowDoms(mm,1,1)%il
           jBeg(mm) = 2; jEnd(mm) = flowDoms(mm,1,1)%jl
           kBeg(mm) = 2; kEnd(mm) = flowDoms(mm,1,1)%kl

           ll = flowDoms(mm,1,1)%cgnsBlockID

           if( writeRindLayer ) then
             if(flowDoms(mm,1,1)%iBegor == 1) iBeg(mm) = 1
             if(flowDoms(mm,1,1)%jBegor == 1) jBeg(mm) = 1
             if(flowDoms(mm,1,1)%kBegor == 1) kBeg(mm) = 1

             if(flowDoms(mm,1,1)%iEndOr == cgnsDoms(ll)%il) &
               iEnd(mm) = flowDoms(mm,1,1)%ie
             if(flowDoms(mm,1,1)%jEndOr == cgnsDoms(ll)%jl) &
               jEnd(mm) = flowDoms(mm,1,1)%je
             if(flowDoms(mm,1,1)%kEndOr == cgnsDoms(ll)%kl) &
               kEnd(mm) = flowDoms(mm,1,1)%ke
           endif
         enddo

         ! Determine the size of one volume solution, i.e. the
         ! number of bytes to store one variable.

         sizeVolumeSol = 0

         if( writeRindLayer ) then
           do mm=1,cgnsNDom
             sizeVolumeSol = sizeVolumeSol + 2*sizeP3D_Record          &
                           + (cgnsDoms(mm)%il+1) * (cgnsDoms(mm)%jl+1) &
                           * (cgnsDoms(mm)%kl+1) * sizeP3D_Real
           enddo
         else
           do mm=1,cgnsNDom
             sizeVolumeSol = sizeVolumeSol + 2*sizeP3D_Record          &
                           + (cgnsDoms(mm)%il-1) * (cgnsDoms(mm)%jl-1) &
                           * (cgnsDoms(mm)%kl-1) * sizeP3D_Real
           enddo
         endif
!
!        ****************************************************************
!        *                                                              *
!        * Writing of the solution variables. Loop over the variables   *
!        * to be written. This means that for a certain variable all    *
!        * blocks are written before the next variable is written.      *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the necessary variables for the writing.
         ! All blocks are written.

         call prepareWritePlot3D(cgnsNDom, cgnsBlockIDs)

         ! Starting of the loop.

         varLoop: do var=1,nVarWritten

           ! If this is not an unsteady higher solution the variable
           ! to be written must be copied into IOVar for all domains.
           ! For the unsteady higher solution it already points to
           ! the correct entries.

           testCopyData: if(.not. unsteadyHigherSol) then

             ! Data must be copied. loop over the number of blocks.

             do mm=1,nDom

               ! Set the pointers to this block. As this is either the
               ! first solution or a time spectral computation when this
               ! part is reached, nn should be used for time spectral
               ! index.

               call setPointers(mm, 1_intType, nn)

               ! Store the variable to be written in IOVar.

               call storeSolInBuffer(dummyBuf, .false., solNames(var), &
                                     iBeg(mm), iEnd(mm), jBeg(mm),     &
                                     jEnd(mm), kBeg(mm), kEnd(mm))
             enddo
           endif testCopyData

           ! Compute the offset and write the variable to file.
           ! A distinction must be made between an unsteady higher
           ! solution and other cases, because in the former case IOVar
           ! points to the actual entry of wOld, while for the latter
           ! it is used as a buffer.

           P3D_Offset = sizeFirstRecord + sizeHeader &
                      + (var-1)*sizeVolumeSol

           if( unsteadyHigherSol ) then
             call writePlot3DVar(nn, one, var)
           else
             call writePlot3DVar(nn, one, 1_intType)
           endif

         enddo varLoop

         ! Release the memory of the variables needed for the writing.

         call releaseMemIOPlot3D
!
!        ****************************************************************
!        *                                                              *
!        * Writing of the convergence histories.                        *
!        *                                                              *
!        ****************************************************************
!
         ! Set the file pointer to the position after the solution.

         disp = sizeFirstRecord + sizeHeader + nVarWritten*sizeVolumeSol
         call mpi_file_seek(fh, disp, mpi_seek_set, ierr)

         ! Write the convergence history and possibly the time
         ! history. Only done by processor 0.

         if(myID == 0) then

           select case (equationMode)
             case (steady, timeSpectral)
               call writePlot3DConvInfo(nn)

             case (unsteady)
               if(nn == 1) then
                 call writePlot3DConvInfo(nn)
                 call writePlot3DTimeHistory
               endif

           end select

         endif

         ! Close the grid file.

         call mpi_file_close(fh, ierr)

         ! Synchronize the processors, just to be sure.

         call mpi_barrier(SUmb_comm_world, ierr)

         ! If this is an unsteady computation reset a couple of things
         ! for the next solution to write.

         if(equationMode == unsteady) then

           ! Set the variables, such that they corresponds to a
           ! higher unsteady solution.

           nVarWritten       = nw
           unsteadyHigherSol = .true.
           writeRindLayer    = .false.

           ! Change the name of the flow variables to the conservative
           ! ones. Note that the turbulence variables do not need to
           ! be changed; the values given in the routine volSolNames
           ! are the correct ones.

           solNames(irho)  = cgnsDensity
           solNames(imx)   = cgnsMomx
           solNames(imy)   = cgnsMomy
           solNames(imz)   = cgnsMomz
           solNames(irhoE) = cgnsEnergy

         endif

       enddo solLoop

       ! Deallocate the memory of IOVar. Note that the first entry
       ! is used as a temporary buffer.

       do nn=1,nDom
         deallocate(IOVar(nn,1)%w, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("writeCGNSVolumeSol", &
                          "Deallocation error for IOVar%w")
       enddo

       ! Deallocate IOVar itself as well as writeBuf and solNames.

       deallocate(IOVar, writeBuf, solNames, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
                        "Deallocation error for IOVar, writeBuf &
                        &and solNames")

       ! Write a message that the solution filei(s) have been written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "# Volume solution file(s) written"
         print "(a)", "#"
       endif

       end subroutine writePlot3DVolumeSol


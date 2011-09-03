!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DTimeHistory.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-19-2005                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DTimeHistory
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DTimeHistory writes the time history for an unsteady *
!      * computation to the Plot3D file indicated by fh.                *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use inputIO
       use IOModule
       use monitor
       implicit none
!
!      Local variables.
!
       integer :: ierr, sizeWrite
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nn, mm, ii, jj
       integer(kind=intType) :: sizeP3D_Int, sizeP3D_Record, sizeP3D_Real
       integer(kind=intType) :: sizeBuf, nStepsTotal

       integer(kind=intRecordPLOT3DType) :: sizeRecord
       integer(kind=intPLOT3DType), dimension(2) :: intPlot3DBuf

       real(kind=4), dimension(:), allocatable :: buf4
       real(kind=8), dimension(:), allocatable :: buf8

       character(len=maxCGNSNameLen) :: message
       character, dimension(:), allocatable :: writeBuf
!
!-- eran-tconv
       logical :: tconv_file_exist
       integer :: nnr
       character :: formatString*30
!-- eran-tconv
!
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the total number of time steps in nStepsTotal.

       nStepsTotal = timeStepUnsteady + nTimeStepsRestart

       select case (precisionSol)
         case (precisionSingle)
           sizeP3D_Real = 4
           allocate(buf4(nStepsTotal), stat=ierr)

         case (precisionDouble)
           sizeP3D_Real = 8
           allocate(buf8(nStepsTotal), stat=ierr)
       end select

       if(ierr /= 0)                              &
         call terminate("writePlot3DTimeHistory", &
                        "Memory allocation failure for either buf4 &
                        &or buf8")

       ! Copy the size of the integer types in the Plot3D file in
       ! sizeP3D_Int and sizeP3D_Record. The reason is that the
       ! variables nBytesPerIntPLOT3D and nBytesPerRecordIntPLOT3D are
       ! of type integer and not integer(kind=intType).

       sizeP3D_Int    = nBytesPerIntPLOT3D
       sizeP3D_Record = nBytesPerRecordIntPLOT3D

       ! Determine the size of the write buffer and allocate the memory.

       sizeBuf = nStepsTotal*sizeP3D_Real + 2*sizeP3D_Record

       mm      = nMon*maxCGNSNameLen + 4*sizeP3D_Record &
               + maxCGNSNameLen      + 2*sizeP3D_Int
       sizeBuf = max(sizeBuf, mm)

       allocate(writeBuf(sizeBuf), stat=ierr)
       if(ierr /= 0)                              &
         call terminate("writePlot3DTimeHistory", &
                        "Memory allocation failure for writeBuf")

       ! Write the convergence header to the write buffer.
       ! First a message that this is the convergence info and
       ! the number of iterations and monitoring variables.

       ii = 1
       sizeRecord = maxCGNSNameLen + 2*sizeP3D_Int
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       message = TimeHistory
       jj = maxCGNSNameLen
       call writeVarsToBuffer(writeBuf(ii), message, jj)
       ii = ii + jj

       intPlot3DBuf(1) = nMon
       intPlot3DBuf(2) = nStepsTotal
       if( P3D_ByteSwap ) &
         call byteswap(intPlot3DBuf, sizeP3D_Int, 2_intType)

       jj = 2*sizeP3D_Int
       call writeVarsToBuffer(writeBuf(ii), intPlot3DBuf, jj)
       ii = ii + jj

       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the names of the monitoring variables.

       sizeRecord = nMon*maxCGNSNameLen
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       jj = maxCGNSNameLen
       do mm=1,nMon
         call writeVarsToBuffer(writeBuf(ii), monNames(mm), jj)
         ii = ii + jj
       enddo

       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the buffer to file.

       sizeWrite = ii - 1
       call mpi_file_write(fh, writeBuf, sizeWrite, &
                           mpi_character, status, ierr)

       ! Determine the size of the record for the time values
       ! and write that value to writeBuf.

       ii = 1
       sizeRecord = nStepsTotal*sizeP3D_Real
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the time values to writeBuf, depending on the
       ! precision to be used.

       select case (precisionSol)
         case (precisionSingle)

           do nn=1,nStepsTotal
             buf4(nn) = timeArray(nn)
           enddo

           if( P3D_ByteSwap ) &
             call byteswap(buf4, sizeP3D_Real, nStepsTotal)

           jj = nStepsTotal*sizeP3D_Real
           call writeVarsToBuffer(writeBuf(ii), buf4, jj)
           ii = ii + jj

         !===============================================================

         case (precisionDouble)

           do nn=1,nStepsTotal
             buf8(nn) = timeArray(nn)
           enddo

           if( P3D_ByteSwap ) &
             call byteswap(buf8, sizeP3D_Real, nStepsTotal)

           jj = nStepsTotal*sizeP3D_Real
           call writeVarsToBuffer(writeBuf(ii), buf8, jj)
           ii = ii + jj

       end select

       ! Write the closing record for the time values.

       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Write the buffer to file.

       sizeWrite = ii - 1
       call mpi_file_write(fh, writeBuf, sizeWrite, &
                           mpi_character, status, ierr)

       ! Determine the size of the record for the time history
       ! and write that value to writeBuf.

       ii = 1
       sizeRecord = nStepsTotal*nMon*sizeP3D_Real
       if( P3D_ByteSwap ) &
         call byteswap(sizeRecord, sizeP3D_Record, 1_intType)
       call writeVarsToBuffer(writeBuf(ii), sizeRecord, sizeP3D_Record)
       ii = ii + sizeP3D_Record

       ! Loop over the convergence histories.

       monitorLoop: do mm=1,nMon

         ! Copy the convergence history into either buf4 or buf8.
         ! Apply the byte swapping and store it in writeBuf.

         select case (precisionSol)
           case (precisionSingle)

             do nn=1,nStepsTotal
               buf4(nn) = timeDataArray(nn,mm)
             enddo

             if( P3D_ByteSwap ) &
               call byteswap(buf4, sizeP3D_Real, nStepsTotal)

             jj = nStepsTotal*sizeP3D_Real
             call writeVarsToBuffer(writeBuf(ii), buf4, jj)
             ii = ii + jj

           !=============================================================

           case (precisionDouble)

             do nn=1,nStepsTotal
               buf8(nn) = timeDataArray(nn,mm)
             enddo

             if( P3D_ByteSwap ) &
               call byteswap(buf8, sizeP3D_Real, nStepsTotal)

             jj = nStepsTotal*sizeP3D_Real
             call writeVarsToBuffer(writeBuf(ii), buf8, jj)
             ii = ii + jj

         end select

         ! Write the closing integer of the record to the buffer
         ! if this is the last monitoring variable to be written.

         if(mm == nMon) then
           call writeVarsToBuffer(writeBuf(ii), sizeRecord, &
                                  sizeP3D_Record)
           ii = ii + sizeP3D_Record
         endif

         ! Write the buffer to file.

         sizeWrite = ii - 1
         call mpi_file_write(fh, writeBuf, sizeWrite, &
                             mpi_character, status, ierr)

         ! Set ii to 1 for the next variable.

         ii = 1

       enddo monitorLoop
!
!--eran-tconv starts
! 
       write(*,*)'# writing time monitoring on su_tconv.dat'
       if(nMon <= 8)then
          write(formatString,'("(i8,1x,",i1,"(1pg14.7,1x))")')nMon+1
       else
          write(formatString,'("(i8,1x,",i2,"(1pg14.7,1x))")')nMon+1
       end if
!
       inquire(file='su_tconv.dat',exist=tconv_file_exist)

       if(tconv_file_exist)then
          open(unit=89,file='su_tconv.dat',form='formatted',status='old')
! skip the 2 header lines
          do ii=1,3
             read(89,*,iostat=ierr)
             if(ierr > 0)then
                close(89)
                go to 111
             end if
          end do ! ii
!
! read convergence data of previous run
!
          do nn=1,nTimeStepsRestart
             read(89,*,iostat=ierr)nnr,timeArray(nn),&
                  (timeDataArray(nn,mm),mm=1,nMon)
             if(ierr > 0)then
                close(89)
                go to 111
             end if
          end do! nn
             rewind(89)
       else ! tconv_file_exist
111       continue
          if(nTimeStepsRestart > 0)then
             write(*,*)'# WARNING: This is a restart run ',&
               'but something wrong with old  su_tconv.dat.'
             write(*,*)'Monitoring will be corrupted due to a bug in PL3D data'
          end if
          
             open(unit=89,file='su_tconv.dat',form='formatted',status='new')

       end if ! tconv_file_exist
!
       write(89,*)'Title: "Time Convergence file of SUMB"'
       write(89,&
            '("VARIABLES = ",1h","Nstep",1h",1x,1h","Time",1h",1x,50(1h",a,1h",1x))')& 
            (monNames(mm)(1:len_trim(monNames(mm))),mm=1,nMon)
       write(89,*)'-------------'
       
       do nn=1,nStepsTotal
          write(89,formatString)nn,timeArray(nn),(timeDataArray(nn,mm),mm=1,nMon)
       end do ! nn
       close(89)
!
!-- end eran-tconv
!
       ! Release the memory of the temporary buffers.

       select case (precisionSol)
         case (precisionSingle)
           deallocate(buf4, stat=ierr)

         case (precisionDouble)
           deallocate(buf8, stat=ierr)
       end select

       if(ierr /= 0)                              &
         call terminate("writePlot3DTimeHistory", &
                        "Deallocation failure for either buf4 or buf8")

       deallocate(writeBuf, stat=ierr)
       if(ierr /= 0)                              &
         call terminate("writePlot3DTimeHistory", &
                        "Deallocation failure for writeBuf")

       end subroutine writePlot3DTimeHistory

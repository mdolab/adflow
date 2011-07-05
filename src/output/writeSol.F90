!
!      ******************************************************************
!      *                                                                *
!      * File:          writeSol.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-14-2003                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeSol
!
!      ******************************************************************
!      *                                                                *
!      * writeSol controls the writing of a new grid file, a volume     *
!      * solution file and a surface solution file. And if needed the   *
!      * parameter file is updated to the new situation.                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use extraOutput
       use flowVarRefState
       use inputIO
       use inputTimeSpectral
       use killSignals
       use monitor
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If residuals must be written to the volume solution and if halo
       ! values must be stored, exchange the data here, because in the
       ! next call the communication buffers are deleted.

       if(writeVolume     .and. storeRindLayer .and.  &
          (volWriteResRho  .or. volWriteResMom .or.   &
           volWriteResRhoE .or. volWriteResTurb))     &
         call resHalo1(1_intType, 1_intType, nw)

       ! Temporary deallocate some memory, such that writeSolution
       ! is not a memory killer.

       call deallocateTempMemory(.true.)

       ! If blanks are to be written then change the iblank values on
       ! the overset boundary such that output data there isn't blanked.

       if(writeGrid .or. volWriteBlank .or. surfWriteBlank) then
         do nn=1,nTimeIntervalsSpectral
           do i=1,nDom
             call setPointers(i, 1_intType, nn)
             call changeIblanks(.false., 1_intType)
           enddo
         enddo
       endif

       ! If the reading format is different from the writing format
       ! and a volume solution file must be written automatically
       ! a grid file is written as well if this is the first time
       ! a solution file is written.

       writePlot3DConn    = .false.
       writeFormatInParam = .false.

       if((fileFormatWrite /= fileFormatRead) .and. &
          writeVolume .and. firstWrite) then
         writeGrid          = .true.
         writePlot3DConn    = .true.
         writeFormatInParam = .true.
         firstWrite         = .false.
       endif

       ! Write the files. The routines called depend on the IO
       ! format used.

       select case(fileFormatWrite)
         case (cgnsFormat)
           call setHelpVariablesWriting
           call writeCGNSGridFile
           call writeCGNSVolumeSol
           call writeCGNSSurfaceSol
           call releaseHelpVariablesWriting

         case (plot3DFormat)
           call writePlot3DGridFile
           call writePlot3DConnFile
           call writePlot3DVolumeSol
           call writePlot3DSurfaceSol
       end select

       ! Update the parameter file, if needed.

       call updateParamfile

       ! Release the memory of the file names.

       deallocate(gridFileNames, volSolFileNames, &
                  surfSolFileNames, stat=ierr)
       if(ierr /= 0)                &
         call terminate("writeSol", &
                        "Deallocation failure for the file names.")

       ! If blanks were written then change the iblank values on the
       ! overset boundary back to 0.

       if(writeGrid .or. volWriteBlank .or. surfWriteBlank) then
         do nn=1,nTimeIntervalsSpectral
           do i=1,nDom
             call setPointers(i, 1_intType, nn)
             call changeIblanks(.false., 0_intType)
           enddo
         enddo
       endif

       ! Allocate the memory again that was deallocated in the beginning
       ! of this routine.

       call allocateTempMemory(.true.)

#ifndef USE_NO_SIGNALS

       ! It is possible that a kill signal was sent during the writing.
       ! Therefore determine the global signal as the maximum of the
       ! local ones.

       call mpi_allreduce(localSignal, globalSignal, 1, sumb_integer, &
                          mpi_max, SUmb_comm_world, ierr)
#endif

       end subroutine writeSol

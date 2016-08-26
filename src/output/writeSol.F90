!
!       File:          writeSol.f90                                    
!       Author:        Edwin van der Weide, Steve Repsher              
!       Starting date: 03-14-2003                                      
!       Last modified: 03-29-2006                                      
!
       subroutine writeSol
!
!       writeSol controls the writing of a new grid file, a volume     
!       solution file and a surface solution file. And if needed the   
!       parameter file is updated to the new situation.                
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
       use utils, only : terminate, deallocateTempMemory, allocateTempMemory
       use haloExchange, only : resHalo1
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, nn

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

       ! If the reading format is different from the writing format
       ! and a volume solution file must be written automatically
       ! a grid file is written as well if this is the first time
       ! a solution file is written.

       writeFormatInParam = .false.

       ! Write the files. The routines called depend on the IO
       ! format used.

       call setHelpVariablesWriting
       call writeCGNSGridFile
       call writeCGNSVolumeSol
       call writeCGNSSurfaceSol
       call releaseHelpVariablesWriting

       ! Update the parameter file, if needed.

       call updateParamfile

       ! Release the memory of the file names.

       deallocate(gridFileNames, volSolFileNames, &
                  surfSolFileNames, stat=ierr)
       if(ierr /= 0)                &
         call terminate("writeSol", &
                        "Deallocation failure for the file names.")

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

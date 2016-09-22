
           subroutine writeSol
!
!       writeSol controls the writing of a new grid file, a volume     
!       solution file and a surface solution file.            
!
       use constants
       use extraOutput
       use communication, only : adflow_comm_world
       use monitor, only : writeVolume
       use inputIO, only : storeRindLayer
       use killSignals, only : localSignal, globalSignal
       use flowVarRefState, only : nw
       use outputMod, only : setHelpVariablesWriting, &
            releaseHelpVariablesWriting, gridFileNames, volSolFileNames, &
            surfSolFileNames
       use writeCGNSGrid, only : writeCGNSGridFile
       use writeCGNSVolume, only : writeCGNSVolumeSol
       use writeCGNSSurface, only : writeCGNSSurfaceSol
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

       ! Write the files. The routines called depend on the IO
       ! format used.

       call setHelpVariablesWriting
       call writeCGNSGridFile
       call writeCGNSVolumeSol
       call writeCGNSSurfaceSol
       call releaseHelpVariablesWriting

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

       call mpi_allreduce(localSignal, globalSignal, 1, adflow_integer, &
                          mpi_max, ADflow_comm_world, ierr)
#endif

       end subroutine writeSol



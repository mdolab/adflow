!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSVolumeSol.F90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-11-2003                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSVolumeSol
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSVolumeSol and its subroutines write the cell          *
!      * centered CGNS solution file(s).                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("writeCGNSVolumeSol", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use block
       use cgnsGrid
       use communication
       use inputPhysics
       use IOModule
       use su_cgns
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
       integer(kind=intType) :: nVolSolvar, nVolDiscrVar

       character(len=maxCGNSNameLen), &
                                   dimension(:), allocatable :: solNames
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

       ! Return immediately if no volume solution files must be written.

       if(nVolSolToWrite == 0) return

       ! Write a message that the solution file(s) are being written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "# Writing volume solution file(s) ..."
       endif

       ! Open the CGNS file(s), the convergence info and if needed the
       ! time accurate data for an unsteady computation. This is only
       ! done by processor 0.

       if(myID == 0) then
         call openCGNSVolumeSol
         call writeCGNSConvInfo

         if(equationMode == unsteady) call writeCGNSTimeHistory
       endif

       ! Determine the number of variables to be written to the volume
       ! solution file(s) as well as the CGNS names.

       call numberOfVolSolVariables(nVolSolvar, nVolDiscrVar)
       allocate(solNames(nVolSolvar+nVolDiscrVar), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
                        "Memory allocation failure for solNames")
       call volSolNames(solNames)

       ! Loop over the number of CGNS blocks and write the cell centered
       ! solution(s) of this block.

       do nn=1,cgnsNDom
         call writeSolCGNSZone(nn, nVolSolvar, nVolDiscrVar, solNames)
       enddo

       ! Close the cgns file(s). Only processor 0 does this.

       if(myID == 0) then
         do nn=1,nVolSolToWrite
           call cg_close_f(fileIDs(nn), ierr)
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSVolumeSol", &
                            "Something wrong when calling cg_close_f")
         enddo

         ! Deallocate the memory of fileIDs and cgnsBases.
         ! These are only allocated on processor 0.

         deallocate(fileIDs, cgnsBases, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("writeCGNSVolumeSol", &
                          "Deallocation error for fileIDs &
                          &and cgnsBases.")
       endif

       ! Deallocate the memory of solNames.

       deallocate(solNames, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
                        "Deallocation error for solNames.")

       ! Deallocate the memory of IOVar. Note that the first entry
       ! is used as a temporary buffer.

       do nn=1,nDom
         deallocate(IOVar(nn,1)%w, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("writeCGNSVolumeSol", &
                          "Deallocation error for IOVar%w")
       enddo

       deallocate(IOVar, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
                        "Deallocation error for IOVar")

       ! Wait until all processors (especially processor 0) reach
       ! this point.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! Write a message that the solution file has been written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "# Volume solution file(s) written"
         print "(a)", "#"
       endif

#endif

       end subroutine writeCGNSVolumeSol

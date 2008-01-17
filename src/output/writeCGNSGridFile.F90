!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSGridFile.F90                           *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 01-20-2004                                      *
!      * Last modified: 03-30-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSGridFile
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSGridFile and its subroutines write the CGNS grid      *
!      * file(s). Typically this is needed when the coordinates have    *
!      * changed due to moving parts, deformation or both.              *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("writeCGNSGridFile", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use cgnsGrid
       use communication
       use IOModule
       use monitor
       use su_cgns
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer, dimension(cgnsNDom) :: cgnsZone

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number and names of the grid files.
       ! Also set the pointers for IOVar needed for the general
       ! treatment of the IO.

       call gridFileNamesWrite

       ! Return immediately if no grids have to be written.

       if(nGridsToWrite == 0) return

       ! Allocate the memory for the fileIDs and the bases.

       allocate(fileIDs(nGridsToWrite), cgnsBases(nGridsToWrite), &
                stat=ierr)
       if(ierr /= 0)                         &
         call terminate("writeCGNSGridFile", &
                        "Memory allocation failure for fileIDs and &
                        &cgnsBases")

       ! Write a message that the grid file(s) are being written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "# Writing grid file(s) ..."
       endif

       ! All grid information is stored on all processors, with the
       ! exception of data which can vary in time or use a large amount
       ! of memory, including the coordinates and the overset holes and
       ! connectivities. Processor 0 writes this information first as a
       ! frame for each file.

       if(myID == 0) then
         do nn=1,nGridsToWrite
           call writeCGNSGridFrame(cgnsZone, nn)
         enddo
       endif

       ! Loop over the number of cgns blocks and write the coordinates
       ! and the overset data one zone at a time to conserve memory.

       do nn=1,cgnsNDom
         call writeCoorCGNSZone(nn, cgnsZone(nn))
         call writeOversetCGNSZone(nn, cgnsZone(nn))
       enddo

       ! Check if the solution must be written in a different file.
       ! If so the files must be closed the memory for the fileIDs
       ! and bases must be released.

       testGridOnly: if(useLinksInCGNS .or. (.not. writeVolume)) then

         ! Processor 0 closes the files.

         if(myID == 0) then
           do nn=1,nGridsToWrite
             call cg_close_f(fileIDs(nn), ierr)
             if(ierr /= all_ok)                    &
               call terminate("writeCGNSGridFile", &
                              "Something wrong when calling cg_close_f")
           enddo
         endif

         ! Release the memory of fileIDs and cgnsBases.

         deallocate(fileIDs, cgnsBases, stat=ierr)
         if(ierr /= 0)                         &
           call terminate("writeCGNSGridFile", &
                          "Deallocation failure for fileIDs and &
                          &cgnsBases")

       endif testGridOnly

       ! Releases the memory of IOVar.

       deallocate(IOVar, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("writeCGNSGridFile", &
                        "Deallocation failure for IOVar")

       ! Wait until all processors (especially processor 0) reach
       ! this point.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! Write a message that the grid file has been written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "# Grid file(s) written"
         print "(a)", "#"
       endif

#endif

       end subroutine writeCGNSGridFile

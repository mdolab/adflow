!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCgnsSurfaceSol.F90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-15-2003                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSSurfaceSol
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSSurfaceSol and its subroutines write the surface      *
!      * solution file(s). The unknowns are stored in the center of the *
!      * surface quadrilaterals.                                        *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("writeCGNSSurfaceSol", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use cgnsGrid
       use communication
       use su_cgns
       use outputMod
       implicit none
!
!      Local parameter, the cell dimension.
!
       integer, parameter :: celldim = 2
!
!      Local variables.
!
       integer :: cgnsInd, ierr

       integer(kind=intType) :: nn, mm, ll
       integer(kind=intType) :: nSolVar, nZonesWritten

       character(len=maxStringLen) :: errorMessage

       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                                solNames
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number and names of the solution files.

       call surfSolFileNamesWrite

       ! Return immediately if no surface solution files must
       ! be written.

       if(nSurfSolToWrite == 0) return

       ! Write a message that the solution file(s) are being written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "# Writing surface solution file(s) ..."
       endif

       ! Allocate the memory for the fileIDs and the bases.

       allocate(fileIDs(nSurfSolToWrite), cgnsBases(nSurfSolToWrite), &
                stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
                        "Memory allocation failure for fileIDs &
                        &and cgnsBases")

       ! Open the cgns file(s) and write the header. This is only done
       ! by processor 0.

       testRootProc: if(myID == 0) then

         ! Loop over the number of surface solution files to write.

         solLoop: do nn=1,nSurfSolToWrite

           ! Open the cgns file for writing and check if it went okay.
           ! Store the file index for later purposes.

           call cg_open_f(surfSolFileNames(nn), mode_write, cgnsInd, &
                          ierr)
           if(ierr /= all_ok) then
             write(errorMessage,101) trim(surfSolFileNames(nn))
 101         format("File",1X,A,1X,"could not be opened by cgns for &
                    &writing")

             call terminate("writeCGNSSurfaceSol", errorMessage)
           endif

           fileIDs(nn) = cgnsInd

           ! Create the base.

           call cg_base_write_f(cgnsInd, "BaseSurfaceSol", celldim, &
                                cgnsPhysdim, cgnsBases(nn), ierr)
           if(ierr /= all_ok)                      &
             call terminate("writeCGNSSurfaceSol", &
                            "Something wrong when calling &
                            &cg_base_write_f")

           ! Write the header in the cgns file.

           call writeCGNSHeader(cgnsInd, cgnsBases(nn))

         enddo solLoop

       endif testRootProc

       ! Determine the number of variables to be written to the surface
       ! solution file as well as the cgns names.

       call numberOfSurfSolVariables(nSolVar)
       allocate(solNames(nSolVar), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
                        "Memory allocation failure for solNames")
       call surfSolNames(solNames)

       ! Loop over the number of cgns blocks and its boundary subfaces
       ! and write the cell centered surface solution of the subface.

       nZonesWritten = 0
       zoneLoop: do nn=1,cgnsNDom

         ! Determine the number of blocks on this processor that belong
         ! to this cgns block.

         mm = nblocksCGNSblock(nn) - nblocksCGNSblock(nn-1)

         ! Loop over the number of boundary subfaces of the original
         ! cgns block and write the cell centered surface solution to
         ! the cgns surface file.

         do ll=1,cgnsDoms(nn)%nBocos

           ! Only write the solution to file if this is a true subface.

           if( cgnsDoms(nn)%bocoInfo(ll)%actualFace )                 &
             call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
                                       nZonesWritten, .false.)
         enddo

         ! Loop over the number of internal block boundaries of the
         ! original grid and write the periodic boundaries.

         do ll=1,cgnsDoms(nn)%n1to1

           ! Only periodic boundaries are written; check for this.

           if( cgnsDoms(nn)%conn1to1(ll)%periodic )                   &
             call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
                                       nZonesWritten, .true.)
         enddo

       enddo zoneLoop

       ! Close the cgns file(s). Only processor 0 does this.

       if(myID == 0) then
         do nn=1,nSurfSolToWrite
           call cg_close_f(fileIDs(nn), ierr)
           if(ierr /= all_ok)                      &
             call terminate("writeCGNSSurfaceSol", &
                            "Something wrong when calling cg_close_f")
         enddo
       endif

       ! Deallocate the memory of solNames, fileIDs and cgnsBases.

       deallocate(solNames, fileIDs, cgnsBases, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
                        "Deallocation error for solNames, fileIDs &
                        &and cgnsBases")

       ! Wait until all processors (especially processor 0) reach
       ! this point.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! Write a message that the solution file(s) have been written.
       ! Of course only processor 0 does this.

       if(myID == 0) then
         print "(a)", "# Surface solution file(s) written"
         print "(a)", "#"
       endif

#endif

       end subroutine writeCGNSSurfaceSol

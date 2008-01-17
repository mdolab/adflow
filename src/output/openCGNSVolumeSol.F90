!
!      ******************************************************************
!      *                                                                *
!      * File:          openCGNSVolumeSol.F90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-11-2003                                      *
!      * Last modified: 07-03-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine openCGNSVolumeSol
!
!      ******************************************************************
!      *                                                                *
!      * openCGNSVolumeSol opens the cgns solution file(s) if needed.   *
!      * If opened the files are opened either for writing or for       *
!      * modification. When the grid file(s) have been written, these   *
!      * files are still open and nothing needs to be done.             *
!      * Only processor 0 performs this task.                           *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use monitor
       use su_cgns
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn

       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("openCGNSVolumeSol", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else

       ! Check if grid files have been written and if the solution must
       ! be stored in the same file. If so only the CGNS header must
       ! be written. Do this and return immediately.

       if(writeGrid .and. (.not. useLinksInCGNS)) then
         do nn=1,nVolSolToWrite
           call writeCGNSHeader(fileIDs(nn), cgnsBases(nn))
         enddo

         return
       endif

       ! Solution files must be created. Allocate the memory for
       ! the fileIDs and the bases.

       allocate(fileIDs(nVolSolToWrite), cgnsBases(nVolSolToWrite), &
                stat=ierr)
       if(ierr /= 0)                         &
         call terminate("openCGNSVolumeSol", &
                        "Memory allocation failure for fileIDs and &
                        &cgnsBases")

       ! Determine the situation we are having here.

       testForLinks: if( useLinksInCGNS ) then

         ! Links are used to the coordinates of the zones. This means
         ! that the files must be opened in write mode.

         do nn=1,nVolSolToWrite
           call cg_open_f(volSolFileNames(nn), mode_write, &
                          fileIDs(nn), ierr)
           if(ierr /= all_ok) then
             write(errorMessage,*) "File ", trim(volSolFileNames(nn)), &
                                   " could not be opened by CGNS&
                                   & for writing"
             call terminate("openCGNSVolumeSol", errorMessage)
           endif

           ! Create the base.

           call cg_base_write_f(fileIDs(nn), cgnsBaseName, cgnsCelldim, &
                                cgnsPhysdim, cgnsBases(nn), ierr)
           if(ierr /= all_ok)                    &
             call terminate("openCGNSVolumeSol", &
                            "Something wrong when calling &
                            &cg_base_write_f")
         enddo

       else testForLinks

         ! Solutions must be written in the same file(s) as the grid.
         ! As the grid file(s) are not written during the current call
         ! to writeSol, these are old files and must therefore be opened
         ! in modify mode.

         do nn=1,nVolSolToWrite
           call cg_open_f(volSolFileNames(nn), mode_modify, &
                          fileIDs(nn), ierr)
           if(ierr /= all_ok) then
             write(errorMessage,*) "File ", trim(volSolFileNames(nn)), &
                                   " could not be opened by CGNS&
                                   & for writing"
             call terminate("openCGNSVolumeSol", errorMessage)
           endif

           ! Simply set the base IDs to 1.

           cgnsBases(nn) = 1
         enddo

       endif testForLinks

       ! Write the CGNS header.

       do nn=1,nVolSolToWrite
         call writeCGNSHeader(fileIDs(nn), cgnsBases(nn))
       enddo

#endif

       end subroutine openCGNSVolumeSol

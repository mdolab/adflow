!
!      ******************************************************************
!      *                                                                *
!      * File:          gridFileNamesWrite.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-11-2005                                      *
!      * Last modified: 11-02-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine gridFileNamesWrite
!
!      ******************************************************************
!      *                                                                *
!      * gridFileNamesWrite determines the names and number of grid     *
!      * files to be written. Furthermore, it sets the pointers for     *
!      * IOVar to make a general treatment of the writing possible.     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use IOModule
       use iteration
       use monitor
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, kk, nAvail

       character(len=7) :: intString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the names and number of grid files to be written.    *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       select case (equationMode)

         case (steady)

           ! Steady state computation. Allocate the memory for the
           ! grid file names. Even if no file needs to be written the
           ! memory is still allocated because the name is always set.

           allocate(gridFileNames(1), stat=ierr)
           if(ierr /= 0)                          &
             call terminate("gridFileNamesWrite", &
                            "Memory allocation failure for grid &
                            &file names")

           ! Set the number of grid files to be written to either 0 or 1
           ! and set the name accordingly. The name is always set,
           ! because it may be needed for a link in the solution file.

           if( writeGrid ) then
             nGridsToWrite    = 1
             gridFileNames(1) = newGridFile
           else
             nGridsToWrite    = 0
             gridFileNames(1) = gridFile
           endif

         !===============================================================

         case (unsteady)

           ! Unsteady computation. For a consistent restart for a
           ! deforming mesh computation nOldLevels grids must be written.
           ! First determine the number of available solutions.

           nAvail = timeStepUnsteady + nTimeStepsRestart + 1
           nAvail = min(nAvail,nOldLevels)

           ! Allocate the memory for the file names. Note that this is
           ! an upper boundary. It is possible that less files need
           ! to be written.

           allocate(gridFileNames(nAvail), stat=ierr)
           if(ierr /= 0)                          &
             call terminate("gridFileNamesWrite", &
                            "Memory allocation failure for &
                            &gridFileNames")

           ! Set the names of the files.

           do nn=1,nAvail
             write(intString,"(i7)") timeStepUnsteady + &
                                     nTimeStepsRestart + 1 - nn
             intString = adjustl(intString)

             gridFileNames(nn) = trim(newGridFile)//"&
                                 &Timestep"//trim(intString)
           enddo

           ! Determine the number of grid files to be written.
           ! This depends on quite a few things.

           if( writeGrid ) then

             ! Initialize nGridsToWrite to 1. This may change
             ! when the mesh is deforming.

             nGridsToWrite = 1

             if( deforming_Grid ) then

               ! Grids deform during the computation. Check if the
               ! older grids must be written.

               do nn=1,(nAvail-1)
                 if(.not. oldSolWritten(nn) ) then
                   nGridsToWrite = nGridsToWrite + 1
                   gridFileNames(nGridsToWrite) = gridFileNames(nn+1)
                 endif
               enddo

             endif

           else

             ! No grids need to be written. Correct the grid file name
             ! to the original grid file.

             nGridsToWrite = 0
             gridFileNames(1) = gridFile

           endif

         !===============================================================

         case (timeSpectral)

           ! Time spectral computation. Allocate the file names.
           ! Again this is an upper bound.

           allocate(gridFileNames(nTimeIntervalsSpectral), stat=ierr)
           if(ierr /= 0)                          &
             call terminate("gridFileNamesWrite", &
                            "Memory allocation failure for &
                            &gridFileNames")

           ! Set the names of the files.

           do nn=1,nTimeIntervalsSpectral
             write(intString,"(i7)") nn
             intString = adjustl(intString)

             gridFileNames(nn) = trim(newGridFile)//"&
                                 &Spectral"//trim(intString)
           enddo

           ! Set the number of grid files to be written.
           ! This depends on quite a few things.

           if( writeGrid ) then

             ! Need some additional checks.

             if( deforming_Grid ) then

               ! Grids deform during the computation and thus
               ! they must be written.

               nGridsToWrite = nTimeIntervalsSpectral

             else if( timeSpectralGridsNotWritten ) then

               ! Grids do not deform, but the time spectral grids have
               ! not been written earlier. So write the grids and set
               ! timeSpectralGridsNotWritten to .false.

               nGridsToWrite = nTimeIntervalsSpectral
               timeSpectralGridsNotWritten = .false.

             else

               ! Although indicated that the grids must be written,
               ! this is not necessary, because they have already been
               ! written earlier and they have not changed.

               nGridsToWrite = 0

             endif

           else

             ! It is not needed to write the grid files.

             nGridsToWrite = 0

           endif

       end select
!
!      ******************************************************************
!      *                                                                *
!      * Determine whether or not to use links in CGNS.                 *
!      *                                                                *
!      ******************************************************************
!
       if( writeGrid ) then

         ! Grid file(s) will be written. Compare the (base) names of the
         ! grid and solution files and set useLinksInCGNS accordingly.

         if(newGridFile == solFile) then
           useLinksInCGNS = .false.
         else
           useLinksInCGNS = .true.
         endif

       else

         ! Grid file(s) will not be written. Compare the (base) names of
         ! the original grid and solution files and set useLinksInCGNS
         ! accordingly.

         if(gridFile == solFile) then
           useLinksInCGNS = .false.
         else
           useLinksInCGNS = .true.
         endif

       endif
!
!      ******************************************************************
!      *                                                                *
!      * Set the pointers for IOVar if grid files need to be written.   *
!      *                                                                *
!      ******************************************************************
!
       testGridsToWrite: if(nGridsToWrite > 0) then

         ! Allocate the memory for IOVar.

         allocate(IOVar(nDom,nGridsToWrite), stat=ierr)
         if(ierr /= 0)                          &
           call terminate("gridFileNamesWrite", &
                          "Memory allocation failure for IOVar")

         ! Set the pointer w of IOVar to the correct coordinates.

         select case(equationMode)

           case (steady, timeSpectral)

             ! Steady state or time spectral computation. Simply set the
             ! pointers to the current coordinates.

             do nn=1,nDom
               do mm=1,nGridsToWrite
                 IOVar(nn,mm)%pointerOffset = 0
                 IOVar(nn,mm)%w => flowDoms(nn,1,mm)%x(1:,1:,1:,:)
               enddo
             enddo

           !=============================================================

           case (unsteady)

             ! Unsteady computation. First coordinates to be written
             ! are the current coordinates.

             do nn=1,nDom
               IOVar(nn,1)%pointerOffset = 0
               IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)
             enddo

             ! It is possible (for a case with deforming meshes) that
             ! older coordinates need to be written as well.

             if( deforming_Grid ) then
               kk = 1
               do mm=1,(nAvail-1)
                 if(.not. oldSolWritten(mm) ) then
                   kk = kk + 1
                   do nn=1,nDom
                     IOVar(nn,kk)%pointerOffset = 0
                     IOVar(nn,kk)%w => flowDoms(nn,1,1)%xOld(mm,1:,1:,1:,:)
                   enddo
                 endif
               enddo
             endif

         end select

       endif testGridsToWrite

       end subroutine gridFileNamesWrite

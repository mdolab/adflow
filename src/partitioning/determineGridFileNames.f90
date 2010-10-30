!
!      ******************************************************************
!      *                                                                *
!      * File:          determineGridFileNames.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-21-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineGridFileNames
!
!      ******************************************************************
!      *                                                                *
!      * determineGridFileNames determines the number and names of the  *
!      * files that contain the grids. For steady computations only one *
!      * file must be present no matter if a restart is performed or    *
!      * not. For unsteady the situation is a little more complicated.  *
!      * If no restart is performed only one file must be present. If a *
!      * restart is performed in unsteady or time spectral mode and a   *
!      * rigid body motion is prescribed again only one grid file is    *
!      * required; however for a consistent restart with deforming      *
!      * meshes the grids in the past must be read as well. If this is  *
!      * not possible only a first order restart can be made in         *
!      * unsteady mode and some kind of interpolation is used for the   *
!      * time spectral method.                                          *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use iteration
       use partitionMod
       implicit none
!
!      Local variables
!
       integer :: ierr

       integer(kind=intType) :: ii, nn, restartID

       character(len=7)            :: integerString
       character(len=maxStringLen) :: tmpName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialization of nOldGridRead and interpolSpectral.

       nOldGridRead     = 1
       interpolSpectral = .true.

       ! Determine the desired number of files from which grids at
       ! certain time levels should be read.  This depends on the
       ! equation mode we have to solve for. Also set the corresponding
       ! file names.

       select case(equationMode)

         case (steady)

           ! Steady computation. Only one grid needs to be read.

           nGridsRead = 1
           allocate(fileIDs(nGridsRead), &
                    gridFiles(nGridsRead), stat=ierr)
           if(ierr /= 0)                              &
             call terminate("determineGridFileNames", &
                            "Memory allocation failure for fileIDs &
                            &and gridFiles")

           gridFiles(1) = gridFile

         !===============================================================

         case (unsteady)

           ! Unsteady computation. A further check is required.

           testMultipleUnsteady: if(deforming_Grid .and. restart) then

             ! A restart is made with deforming meshes. For a consistent
             ! restart nOldLevels grids must be read. First determine
             ! the prefix of the grid file and the time step number
             ! from which a restart should be made.

             ii = len_trim(gridFile)
             do
               if(gridFile(ii:ii) < "0" .or. gridFile(ii:ii) > "9") exit
               ii = ii - 1
             enddo

             ! If the last characters of the file name do not contain a
             ! number, the grid file does not come from a previous 
             ! unsteady deforming mesh computation and therefore only
             ! one grid will be read.

             if(ii == len_trim(gridFile)) then

               nGridsRead = 1
               allocate(fileIDs(nGridsRead), &
                        gridFiles(nGridsRead), stat=ierr)
               if(ierr /= 0)                              &
                 call terminate("determineGridFileNames", &
                                "Memory allocation failure for fileIDs &
                                &and gridFiles")

               gridFiles(1) = gridFile

             else

               ! Read the integer number from the last characters
               ! of the grid file.

               read(gridFile(ii+1:),*) restartID

               ! Allocate the memory for the file names and set them.

               nGridsRead = nOldLevels
               allocate(fileIDs(nGridsRead), &
                        gridFiles(nGridsRead), stat=ierr)
               if(ierr /= 0)                              &
                 call terminate("determineGridFileNames", &
                                "Memory allocation failure for fileIDs &
                                &and gridFiles")

               do nn=1,nGridsRead
                 write(integerString,"(i6)") restartID - nn + 1
                 integerString = adjustl(integerString)
                 gridFiles(nn) = gridFile(:ii)//trim(integerString)
               enddo

             endif

           else testMultipleUnsteady

             ! The computation either starts from scratch or an unsteady
             ! restart for a rigid grid (possibly moving) is made. In
             ! all cases only one grid file is needed.

             nGridsRead = 1
             allocate(fileIDs(nGridsRead), &
                      gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
               call terminate("determineGridFileNames", &
                              "Memory allocation failure for fileIDs &
                              &and gridFiles")

             gridFiles(1) = gridFile

           endif testMultipleUnsteady

           ! Check if the files can be opened.

           do nn=1,nGridsRead
             open(unit=21,file=gridFiles(nn),status="old",iostat=ierr)
             if(ierr /= 0) exit
             close(unit=21)
           enddo

           ! Possibly correct nGridsRead and set nOldGridRead.
           ! If nOldGridRead == 0, i.e. not a valid grid is found,
           ! print an error message and terminate.

           nGridsRead   = nn - 1
           nOldGridRead = nGridsRead

           if(nOldGridRead == 0) then
             if(myID == 0)                              &
               call terminate("determineGridFileNames", &
                              "Grid file(s) could not be opened")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

         !===============================================================

         case (timeSpectral)

           ! Time spectral computation. A further check is required.

           testMultipleTS: if(deforming_Grid .and. restart) then

             ! A restart is made with deforming meshes. For a consistent
             ! restart multiple grids must be read. First determine the
             ! the prefix of the grid file from which a restart should
             ! be made.

             ii = len_trim(gridFile)
             do
               if(gridFile(ii:ii) < "0" .or. gridFile(ii:ii) > "9") exit
               ii = ii - 1
             enddo

             ! If the last characters of the file name do not contain a
             ! number, the grid file does not come from a previous 
             ! time spectral deforming mesh computation and therefore
             ! only one grid will be read.

             if(ii == len_trim(gridFile)) then

               nGridsRead = 1
               allocate(fileIDs(nGridsRead), &
                        gridFiles(nGridsRead), stat=ierr)
               if(ierr /= 0)                              &
                 call terminate("determineGridFileNames", &
                                "Memory allocation failure for fileIDs &
                                &and gridFiles")

               gridFiles(1) = gridFile

             else

               ! Loop to find out how many time instances were used in
               ! the previous computation from which a restart is made.

               nn = 0
               do
                 nn = nn + 1
                 write(integerString,"(i6)") nn
                 integerString = adjustl(integerString)
                 tmpName = gridFile(:ii)//trim(integerString)

                 open(unit=21,file=tmpName,status="old",iostat=ierr)
                 if(ierr /= 0) exit
                 close(unit=21)
               enddo

               nn = nn - 1

               ! Take care of the exceptional situation that nn == 0.
               ! This happens when the restart file ends at with an
               ! integer, but does not correspond to a time spectral
               ! solution. Allocate the memory.

               nGridsRead = max(nn, 1_intType)
               allocate(fileIDs(nGridsRead), &
                        gridFiles(nGridsRead), stat=ierr)
               if(ierr /= 0)                              &
                 call terminate("determineGridFileNames", &
                                "Memory allocation failure for fileIDs &
                                &and gridFiles")

               if(nn == 0) then
                 gridFiles(1) = gridFile
               else
                 do nn=1,nGridsRead
                   write(integerString,"(i6)") nn
                   integerString = adjustl(integerString)
                   gridFiles(nn) = gridFile(:ii)//trim(integerString)
                 enddo
               endif

               ! Check whether or not the coordinates must be interpolated,
               ! i.e. check if nGridsRead == nTimeIntervalsSpectral.

               if(nGridsRead == nTimeIntervalsSpectral) &
                 interpolSpectral = .false.

             endif

           else testMultipleTS

             ! The computation either starts from scratch or a
             ! restart for a rigid grid (possibly moving) is made. In
             ! all cases only one grid file is needed.

             nGridsRead = 1
             allocate(fileIDs(nGridsRead), &
                      gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
               call terminate("determineGridFileNames", &
                              "Memory allocation failure for fileIDs &
                              &and gridFiles")

             gridFiles(1) = gridFile

           endif testMultipleTS

       end select

       end subroutine determineGridFileNames

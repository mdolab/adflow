!
!      ******************************************************************
!      *                                                                *
!      * File:          determineSolFileNames.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-05-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineSolFileNames
!
!      ******************************************************************
!      *                                                                *
!      * determineSolFileNames determines the number and names of the   *
!      * files that contain the solutions. For steady computations only *
!      * one file must be present. For unsteady the situation is a      *
!      * little more complicated. It is attempted to read as many       *
!      * solutions as needed for a consistent restart. If not possible  *
!      * as many as possible solutions are read. For an unsteady        *
!      * computation the order will be reduced; for time spectral mode  *
!      * the solution will be interpolated.                             *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use iteration
       use restartMod
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
       ! Initialize copySpectral and interpolSpectral to .false.

       copySpectral     = .false.
       interpolSpectral = .false.

       ! Determine the desired number of files to be read. This depends
       ! on the equation mode we have to solve for. Also set the
       ! corresponding file names.

       select case(equationMode)

         case (steady)

           ! Steady computation. Only one solution needs to be read.

           nSolsRead = 1
           allocate(solFiles(nSolsRead), stat=ierr)
           if(ierr /= 0)                              &
             call terminate("determineGridFileNames", &
                            "Memory allocation failure for solFiles")

           solFiles(1) = restartFile

         !===============================================================

         case (unsteady)

           ! Unsteady computation. For a consistent restart nOldLevels
           ! solutions must be read. First determine the prefix of the
           ! restart file and the time step number from which a restart
           ! should be made.

           ii = len_trim(restartFile)
           do
             if(restartFile(ii:ii) < "0" .or. &
                restartFile(ii:ii) > "9") exit
             ii = ii - 1
           enddo

           ! If the last characters of the file name do not contain a
           ! number, the restart file does not come from a previous
           ! unsteady computation and therefore only one solution
           ! will be read.

           if(ii == len_trim(restartFile)) then

             nSolsRead = 1
             allocate(solFiles(nSolsRead), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("determineSolFileNames", &
                              "Memory allocation failure for solFiles")

             solFiles(1) = restartFile

           else

             ! Read the integer number from the last characters
             ! of the restart file.

             read(restartFile(ii+1:),*) restartID

             ! Allocate the memory for the file names and set them.

             nSolsRead = nOldLevels
             allocate(solFiles(nSolsRead), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("determineSolFileNames", &
                              "Memory allocation failure for solFiles")

             do nn=1,nSolsRead
               write(integerString,"(i6)") restartID - nn + 1
               integerString = adjustl(integerString)
               solFiles(nn) = restartFile(:ii)//trim(integerString)
             enddo

           endif

           ! Check if the files can be opened.

           do nn=1,nSolsRead
             open(unit=21,file=solFiles(nn),status="old",iostat=ierr)
             if(ierr /= 0) exit
             close(unit=21)
           enddo

           ! Possibly correct nSolsRead. Take nOldGridRead if
           ! a deforming mesh computation is performed.
           ! If nSolsRead == 0, i.e. not a valid solution is found,
           ! print an error message and terminate.

           nSolsRead = nn - 1
           if( deforming_Grid ) &
             nSolsRead = min(nSolsRead, nOldGridRead)

           if(nSolsRead == 0) then
             if(myID == 0)                             &
               call terminate("determineSolFileNames", &
                              "Solution file(s) could not be opened")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           ! Set nOldSolAvail to nSolsRead and check if a consistent
           ! restart can be made. If not, processor 0 prints a warning.

           nOldSolAvail = nSolsRead
           if(myID == 0 .and. nOldSolAvail < nOldLevels) then

             print "(a)", "#"
             print "(a)", "#                 Warning"
             print "(a)", "# Not enough data found for a consistent &
                          &time accurate restart."
             print "(a)", "# Order is reduced in the first time steps &
                          &until enough data is available again."
             print "(a)", "#"

           endif

           ! Set the logicals oldSolWritten, such that nothing is
           ! overwritten when solution files are dumped for this
           ! computation.

           ii = min(nSolsRead,nOldLevels-1)
           do nn=1,ii
             oldSolWritten(nn) = .true.
           enddo

         !===============================================================

         case (timeSpectral)

           ! Time spectral computation. For a consistent restart
           ! nTimeIntervalsSpectral solutions must be read. First 
           ! determine the prefix of the restart file.

           ii = len_trim(restartFile)
           do
             if(restartFile(ii:ii) < "0" .or. &
                restartFile(ii:ii) > "9") exit
             ii = ii - 1
           enddo

           ! If the last characters of the file name do not contain a
           ! number, the restart file does not come from a previous
           ! time spectral computation and therefore only one solution
           ! will be read.

           if(ii == len_trim(restartFile)) then

             nSolsRead = 1
             allocate(solFiles(nSolsRead), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("determineSolFileNames", &
                              "Memory allocation failure for solFiles")

             solfiles(1) = restartFile

           else

             ! Loop to find out how many time instances were used in
             ! the previous computation from which a restart is made.

             nn = 0
             do
               nn = nn + 1
               write(integerString,"(i6)") nn
               integerString = adjustl(integerString)
               tmpName = restartFile(:ii)//trim(integerString)

               open(unit=21,file=tmpName,status="old",iostat=ierr)
               if(ierr /= 0) exit
               close(unit=21)
             enddo

             ! Subtract 1 to obtain the correct number of files.

             nn = nn - 1

             ! Take care of the expeptional situation that nn == 0.
             ! This happens when the restart file ends at with an
             ! integer, but does not correspond to a time spectral
             ! solution. Allocate the memory.

             nSolsRead = max(nn, 1_intType)
             allocate(solFiles(nSolsRead), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("determineSolFileNames", &
                              "Memory allocation failure for solFiles")

             if(nn == 0) then
               solFiles(1) = restartFile
             else
               do nn=1,nSolsRead
                 write(integerString,"(i6)") nn
                 integerString = adjustl(integerString)
                 solFiles(nn) = restartFile(:ii)//trim(integerString)
               enddo
             endif

           endif

           ! Check whether or not the spectral solution must be copied
           ! or interpolated.

           if(nSolsRead == 1) then
             copySpectral = .true.
           else if(nSolsRead /= nTimeIntervalsSpectral) then
             interpolSpectral = .true.
           endif

       end select

       end subroutine determineSolFileNames

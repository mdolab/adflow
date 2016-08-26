       subroutine determineSolFileNames
!
!       determineSolFileNames determines the number and names of the   
!       files that contain the solutions. For steady computations only 
!       one file must be present. For unsteady the situation is a      
!       little more complicated. It is attempted to read as many       
!       solutions as needed for a consistent restart. If not possible  
!       as many as possible solutions are read. For an unsteady        
!       computation the order will be reduced; for time spectral mode  
!       the solution will be interpolated.                             
!
       use constants
       use communication, only : myID
       use inputIO, only : restartFiles
       use inputPhysics, only : equationMode
       use inputTimeSpectral, only : nTimeIntervalsSpectral
       use iteration, only : nOldSolAvail, oldSolWritten, nOldLevels
       use restartMod, only : solFiles, copySpectral, interpolSpectral,&
             nSolsRead
       use utils, only : terminate
       implicit none
!
!      Local variables
!
       integer :: ierr

       integer(kind=intType) :: ii, nn

       character(len=7)            :: integerString
       character(len=maxStringLen) :: tmpName

       ! Initialize copySpectral and interpolSpectral to .false.

       copySpectral     = .false.
       interpolSpectral = .false.

       ! Determine the desired number of files to be read. This depends
       ! on the equation mode we have to solve for. Also set the
       ! corresponding file names.

       select case(equationMode)

         case (steady)

           ! Steady computation. Only one solution needs to be read.
           ! In case a list of restart files were provided in the python
           ! script, we force a read from only the first solution file.

           nSolsRead = 1
           allocate(solFiles(nSolsRead), stat=ierr)
           if(ierr /= 0)                              &
             call terminate("determineSolFileNames", &
                            "Memory allocation failure for solFiles")

           solFiles(1) = restartFiles(1)

           ! Check if the files can be opened, exit if that fails.
           call checkSolFileNames()

         !===============================================================

         case (unsteady)

           ! Unsteady computation. For a consistent restart nOldLevels
           ! solutions must be read. All restart files are provided explicitly
           ! from python script.
           call setSolFileNames()

           ! Check if the files can be opened, exit if that fails.
           call checkSolFileNames()

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
           ! determine the the restart files.
           call setSolFileNames()

           ! Check if the files can be opened, exit if that fails.
           call checkSolFileNames()

           ! Check whether or not the spectral solution must be copied
           ! or interpolated.

           if(nSolsRead == 1) then
             copySpectral = .true.
           else if(nSolsRead /= nTimeIntervalsSpectral) then
             interpolSpectral = .true.
           endif

       end select

       end subroutine determineSolFileNames



       subroutine setSolFileNames
!
!       setSolFileNames allocates and set the solution files that      
!       will be read and loaded in the restart                         
!
       use communication
       use restartMod
       use inputIO 
       use utils, only : terminate
       implicit none
!
!      Local variables
!
       integer :: ierr
       integer(kind=intType) :: nn


       ! The length of the array provided gives the number of nSolsRead.
       nSolsRead = SIZE(restartFiles,1)

       ! Allocate the memory for the file names and set them.
       allocate(solFiles(nSolsRead), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("determineSolFileNames", &
                      "Memory allocation failure for solFiles")

       do nn=1,nSolsRead
         solFiles(nn) = restartFiles(nn)
       enddo

       end subroutine setSolFileNames



       subroutine checkSolFileNames
!
!       checkSolFileNames will check if the provided restart files     
!       are readable on disk. If not readable return fail, if readable 
!       message will be printed to let the user know that restart file 
!       will be tried to read.                                         
!
       use communication
       use restartMod
       use utils, only : terminate
       implicit none
!
!      Local variables
!
       integer :: ierr
       character(len=maxStringLen) :: errorMessage
       integer(kind=intType) :: nn

       do nn=1,nSolsRead
         open(unit=21,file=solFiles(nn),status="old",iostat=ierr)
         if(ierr /= 0) then
           write(errorMessage,*) "Restart file ", trim(solFiles(nn)), &
                                  " could not be opened for reading"
           call terminate("checkSolFileNames", errorMessage)
           exit
         end if
         close(unit=21)

         if(myID == 0) then
           print "(a)", "#"
           write (*,100) trim(solFiles(nn))
           print "(a)", "#"
           100 format("# Found restart file: ", A, 1X)
         end if
       enddo

       end subroutine checkSolFileNames

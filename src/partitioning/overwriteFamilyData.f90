!
!      ******************************************************************
!      *                                                                *
!      * File:          overwriteFamilyData.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-11-2004                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine overwriteFamilyData(sortedFamName, famID)
!
!      ******************************************************************
!      *                                                                *
!      * OverwriteFamilyData overwrites family data specified in the    *
!      * cgns file with the data specified in the parameter file. This  *
!      * can be the rotation information or boundary condition          *
!      * information of the families. The corresponding values          *
!      * specified in the grid file, if any, are ignored later on.      *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputIO
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), dimension(*), intent(in) :: sortedFamName
       integer(kind=intType), dimension(*), intent(in) :: famID
!
!      Local parameter.
!
       integer, parameter :: readUnit = 32
!
!      Local variables.
!
       integer :: ios, pos, ierr

       integer(kind=intType) :: nn

       character(len=2*maxStringLen) :: errorMessage
       character(len=maxStringLen)   :: keyword, value
       character(len=512)            :: string
       character(len=15)             :: familyTest
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize rotatingFrameSpecified to .false. for all families
       ! and initialize the center and rates to zero to avoid
       ! uninitialized data. Also set the number of boundary data sets
       ! to 0.

       do nn=1,cgnsNFamilies
         cgnsFamilies(nn)%rotatingFrameSpecified = .false.
         cgnsFamilies(nn)%rotCenter = zero
         cgnsFamilies(nn)%rotRate   = zero
         cgnsFamilies(nn)%nDataSet  = 0
       enddo

       ! Open the parameter file for reading. This should normally be
       ! no problem, because the input parameters have already been read
       ! from it. But a check does not really hurt.

       open(unit=readUnit, file=paramFile, status="old", &
            action="read", iostat=ios)

       ! Print an error message if the file could not be opened.
       ! The message is only printed by processor 0, while the others
       ! wait to get killed.

       if(ios /= 0) then

         if(myID == 0) then
           write(errorMessage,*) "Parameter file ", trim(paramFile), &
                                 " not found anymore."
           call terminate("overwriteFamilyData", errorMessage)
         endif

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Loop to read the data

       dataLoop: do

         ! Read a string from the file. In case the end of the file
         ! has been reached, exit the loop.

         read(unit=readUnit, fmt="(a512)", iostat=ios) string
         if(ios /= 0) exit dataLoop

         ! Replace all the tab and return characters by spaces.

         call replaceTabsAndReturns(string)

         ! Get rid of the leading and trailing spaces in string.

         string = adjustl(string)
         string = trim(string)

         ! In case this is an empty string or if the first character
         ! is a comment sign, continue with the next line.

         if(len_trim(string) == 0) cycle
         if(string(:1) == "#") cycle

         ! Find a possible comment sign somewhere in the string.
         ! If present the info following the comment sign is ignored.

         pos = index(string, "#")
         if(pos > 0) then
           string = string(:pos-1)
           string = trim(string)
         endif

         ! Search for the : in the string. If not present, continue
         ! with the next line.

         pos = index(string, ":")
         if(pos == 0) cycle

         ! Create the strings keyword and value and get rid of the
         ! leading and trailing spaces of keyword. As this operation has
         ! already been applied for string, only a trim needs to be done.

         keyword = string(:pos-1)
         keyword = trim(keyword)

         value = string(pos+1:)

         ! Extract the substring needed for the family test and
         ! remove that part from the string keyword.

         familyTest = keyword(:15)
         call convertToLowerCase(familyTest)

         keyword = keyword(16:)
         keyword = adjustl(keyword)
         keyword = trim(keyword)

         ! Check whether the first 15 characters of keyword are
         ! "rotating family" or "boundary family". If so call the
         ! corresponding routine to extract that info.

         if(familyTest == "rotating family") &
           call overwriteRotatingRate(sortedFamName, famID, &
                                      keyword, value)

         if(familyTest == "boundary family") &
           call overwriteBoundaryInfo(sortedFamName, famID, &
                                      keyword, value, readUnit)

       enddo dataLoop

       ! Close the file

       close(unit=readUnit)

       end subroutine overwriteFamilyData

!========================================================================

       subroutine overwriteRotatingRate(sortedFamName, famID, &
                                        famName, value)
!
!      ******************************************************************
!      *                                                                *
!      * overwriteRotatingRate overwrites the (possible) rotating       *
!      * rate specified in the cgns file by the family values specified *
!      * in the input parameter file. This data is stored in the string *
!      * value.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(in) :: famName, value

       character(len=*), dimension(*), intent(in) :: sortedFamName
       integer(kind=intType), dimension(*), intent(in) :: famID
!
!      Local variables
!
       integer(kind=intType) :: nn
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Search the name in sortedFamName and check if it was found.

       nn = bsearchStrings(famName, sortedFamName, cgnsNFamilies)

       testFound: if(nn > 0) then

         ! Family name is found. Extract the info. It is assumed
         ! that the rotation center is specified in the same units
         ! as the grid and the rotation rate is radians per second.

         nn = famID(nn)
         cgnsFamilies(nn)%rotatingFrameSpecified = .true.

         read(value,*) cgnsFamilies(nn)%rotCenter(1), &
                       cgnsFamilies(nn)%rotCenter(2), &
                       cgnsFamilies(nn)%rotCenter(3), &
                       cgnsFamilies(nn)%rotRate(1),   &
                       cgnsFamilies(nn)%rotRate(2),   &
                       cgnsFamilies(nn)%rotRate(3)

         ! Set changing_Grid to .true.

         changing_Grid = .true.

       else testFound

         ! Family name not present. Processor 0 prints a warning.

         if(myID == 0) then

           print "(a)",  "#"
           print "(a)",  "#*==================== !!! Warning !!! &
                         &======================"
           print "(3a)", "#* Unknown rotating family, ", &
                         trim(famName),                  &
                         ", encountered in the input file"
           print "(a)",  "#* Information is ignored."
           print "(a)",  "#*=====================================&
                         &======================"
           print "(a)",  "#"

         endif

       endif testFound

       end subroutine overwriteRotatingRate

!========================================================================

       subroutine overwriteBoundaryInfo(sortedFamName, famID, &
                                        famName, value, readUnit)
!
!      ******************************************************************
!      *                                                                *
!      * overwriteBoundaryInfo overwrites the (possible) boundary       *
!      * conditions specified in the cgns file by the family values     *
!      * specified in the input parameter file.                         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use su_cgns
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(in)    :: famName
       character(len=*), intent(inout) :: value
       integer, intent(in)             :: readUnit

       character(len=*), dimension(*), intent(in) :: sortedFamName
       integer(kind=intType), dimension(*), intent(in) :: famID
!
!      Local variables
!
       integer :: ierr, pos, ios

       integer(kind=intType) :: i, j, nn, nVar, nDim, nTot
       integer(kind=intType), dimension(3) :: nPoints

       character(len=512)            :: string
       character(len=maxCGNSNameLen) :: varName

       type(cgnsBCDatasetType), pointer, dimension(:) :: dataSet
       type(cgnsBCDataArray),   pointer, dimension(:) :: BCData
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
       logical               :: dataIsDirichlet
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Search the name in sortedFamName and check if it was found.

       nn = bsearchStrings(famName, sortedFamName, cgnsNFamilies)

       testFound: if(nn > 0) then

         ! Family name is found. Determine the index in the array.

         nn = famID(nn)

         ! Get rid of the leading and trailing spaces in value and
         ! store it in string

         value  = adjustl(value)
         value  = trim(value)
         string = value

         ! Convert value to lower case and check if a colon is present.

         call convertToLowerCase(value)
         pos = index(value, ":")

         testColon: if(pos > 0) then

           ! Test for the possibility to monitor the mass flow on
           ! a per family basis.

           if(value(:pos-1) == "monitor mass flow") then

             ! Determine whether or not the mass flow must be monitored
             ! for this family.

             value = value(pos+1:)
             value  = adjustl(value)

             select case (value)
               case ("yes")
                 cgnsFamilies(nn)%monitorMassFlow = .true.
               case ("no")
                 cgnsFamilies(nn)%monitorMassFlow = .false.
               case default
                 write(string,"(3a)") "Boundary family, ", &
                                      trim(famName),       &
                                      ": monitor mass flow should &
                                      &be either yes or no in the &
                                      &parameter file."
                 if(myID == 0) &
                   call terminate("overwriteBoundaryInfo", string)
                 call mpi_barrier(SUmb_comm_world, ierr)
             end select

             ! Information has been extracted, so make a return.

             return
           endif

         endif testColon

         ! No logical information for this family. Determine the number
         ! of variables specified. Their names are stored in the
         ! string value.

         nVar = 0
         do
           ! Condition to exit the loop.

           if(len_trim(value) == 0) exit

           ! Update nVar and determine the next space in the string.
           ! If not found exit the loop.

           nVar = nVar + 1
           pos = index(value, " ")
           if(pos == 0) exit

           ! Skip the current variable name and set the string value
           ! for the next one. Get rid of the spaces.

           value = value(pos:)
           value = adjustl(value)
           value = trim(value)

         enddo

         ! Return if no variables have been specified.

         if(nVar == 0) return

         ! Store the original string again in value for later purposes.

         value = string

         ! Find the line, which contains the number of points.

         do
           read(unit=readUnit, fmt="(a512)", iostat=ios) string
           if(ios /= 0) then
             if(myID == 0) then
               write(string,100) trim(cgnsFamilies(nn)%familyName)
               call terminate("overwriteBoundaryInfo", string)
             endif
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           ! Replace all the tab and return characters by spaces
           ! and get rid of the spaces and exit the loop if it is
           ! not a comment or an empty line.

           call replaceTabsAndReturns(string)
           string = adjustl(string)
           string = trim(string)

           if(len_trim(string) == 0) cycle
           if(string(:1) /= "#") exit
         enddo

         ! Find a possible comment sign somewhere in the string.
         ! If present the info following the comment sign is ignored.

         pos = index(string, "#")
         if(pos > 0) then
           string = string(:pos-1)
           string = trim(string)
         endif

         ! Find the "=" in string and adapt string such that it only
         ! contains the part after the "=" sign.

         pos = index(string, "=")
         if(pos == 0) then
           if(myID == 0) then
             write(string,100) trim(cgnsFamilies(nn)%familyName)
             call terminate("overwriteBoundaryInfo", string)
           endif
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         string = string(pos+1:)

         ! Determine the number of dimensions for which data is
         ! specified. This should be between 1 and 3. Determine in same
         ! loop the total number of points specified.

         nDim = 0
         nTot = 1
         do i = 1,3

           ! Get rid of the leading spaces and check if this is
           ! still a nonzero string. If not, exit the loop.

           string = adjustl(string)
           string = trim(string)
           if(len_trim(string) == 0) exit

           ! Update nDim and read its number of data points and update
           ! the total number of points specified.

           nDim = nDim + 1
           read(string,*) nPoints(nDim)
           nTot = nTot*nPoints(nDim)

           ! Skip the current variable name and set the string value for
           ! the next one. If no next variable is found, exit the loop.

           pos = index(string, " ")
           if(pos == 0) exit
           string = string(pos:)

         enddo

         ! Check if nDim > 0. If not write an error message.

         if(nDim == 0) then
           if(myID == 0) then
             write(string,100) trim(cgnsFamilies(nn)%familyName)
             call terminate("overwriteBoundaryInfo", string)
           endif
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Allocate the memory for the boundary data sets and set
         ! the pointer afterwards to improve readability.

         cgnsFamilies(nn)%nDataSet = nVar

         allocate(cgnsFamilies(nn)%dataSet(nVar), stat=ierr)
         if(ierr /= 0) &
           call terminate("overwriteBoundaryInfo", &
                          "Memory allocation failure for dataSet")

         dataSet => cgnsFamilies(nn)%dataSet

         ! Loop over the number of variables prescribed to store their
         ! names, allocate the memory, etc.

         do i=1,nVar

           ! Initialize some data for this data set.

           dataSet(i)%datasetName      = "ReadFromParamFile"
           dataSet(i)%BCType           = cgnsFamilies(nn)%BCTypeCGNS
           dataSet(i)%nDirichletArrays = 0
           dataSet(i)%nNeumannArrays   = 0

           nullify(dataSet(i)%dirichletArrays)
           nullify(dataSet(i)%neumannArrays)

           ! Determine the next space in value and set the name of
           ! the variable. If a space is found remove the current name
           ! from the string value.

           pos = index(value, " ")
           if(pos == 0) then
             varName = value
           else
             varName = value(:pos-1)

             value = value(pos:)
             value = adjustl(value)
             value = trim(value)
           endif

           ! Determine if this is a Dirichlet or a Neumann type of
           ! boundary condition, allocate the corresponding array and
           ! set the pointer BCData to this array.

           if( dataIsDirichlet(varName) ) then

             allocate(dataSet(i)%dirichletArrays(1), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("overwriteBoundaryInfo", &
                              "Memory allocation failure for &
                              &dirichletArrays")

             dataSet(i)%nDirichletArrays = 1
             BCData => dataSet(i)%dirichletArrays

           else

             allocate(dataSet(i)%neumannArrays(1), stat=ierr)
             if(ierr /= 0)                               &
               call terminate("overwriteBoundaryInfo", &
                              "Memory allocation failure for &
                              &neumannArrays")

             dataSet(i)%nNeumannArrays = 1
             BCData => dataSet(i)%neumannArrays

           endif

           ! Set the units for the data specified to SI-units,
           ! because data specified in the param file must be
           ! given in si-units.

           BCData(1)%mass  = Kilogram
           BCData(1)%len   = Meter
           BCData(1)%time  = Second
           BCData(1)%temp  = Kelvin
           BCData(1)%angle = Radian

           ! Store the name of the variable in BCData

           BCData(1)%arrayName = varName

           ! Set the number of dimensions and the number of points in
           ! each dimension.

           BCData(1)%nDimensions = nDim

           do j=1,nDim
             BCData(1)%dataDim(j) = nPoints(j)
           enddo

           ! Allocate the memory for the data array.

           allocate(BCData(1)%dataArr(nTot), stat=ierr)
           if(ierr /= 0) &
           call terminate("overwriteBoundaryInfo", &
                          "Memory allocation failure for dataArr")
         enddo

         ! Loop over the total number of data points which must be read.

         loopDataPoints: do j=1,nTot

           ! Find the next line with data points.

           do
             ! Read the next line.

             read(unit=readUnit, fmt="(a512)", iostat=ios) string
             if(ios /= 0) then
               if(myID == 0) then
                 write(string,100) trim(cgnsFamilies(nn)%familyName)
                 call terminate("overwriteBoundaryInfo", string)
               endif
               call mpi_barrier(SUmb_comm_world, ierr)
             endif

             ! Replace all the tab and return characters by spaces and
             ! get rid of the leading and trailing spaces in string.
             ! The line is found when this is not a comment and not
             ! an empty line.

             call replaceTabsAndReturns(string)
             string = adjustl(string)
             string = trim(string)

             if(len_trim(string) == 0) cycle
             if(string(:1) /= "#") exit
           enddo

           ! Find a possible comment sign somewhere in the string.
           ! If present the info following the comment sign is ignored.

           pos = index(string, "#")
           if(pos > 0) then
             string = string(:pos-1)
             string = trim(string)
           endif

           ! Loop to read the variables for this data point.

           do i=1,nVar

             ! Find out whether this is a dirichlet or a neumann
             ! boundary and store the data accordingly.

             if(dataSet(i)%nDirichletArrays > 0) then
               read(string,*,iostat=ios) &
                 dataSet(i)%dirichletArrays(1)%dataArr(j)
             else
               read(string,*,iostat=ios) &
                 dataSet(i)%neumannArrays(1)%dataArr(j)
             endif

             ! Print an error message if the data could not be read.

             if(ios /= 0) then
               if(myID == 0) then
                 write(string,101) trim(cgnsFamilies(nn)%familyName)
                 call terminate("overwriteBoundaryInfo", string)
               endif
               call mpi_barrier(SUmb_comm_world, ierr)
             endif

             ! If this is the last variable to be read exit the loop,
             ! because there is no info left on this line.

             if(i == nVar) exit

             ! Find the next space in the string and skip the current
             ! entry as well as the spaces in between. If the next string
             ! is not found an error message is printed.

             pos = index(string, " ")
             if(pos == 0) then
               if(myID == 0) then
                 write(string,100) trim(cgnsFamilies(nn)%familyName)
                 call terminate("overwriteBoundaryInfo", string)
               endif
               call mpi_barrier(SUmb_comm_world, ierr)
             endif

             string = string(pos:)
             string = adjustl(string)
             string = trim(string)
           enddo

         enddo loopDataPoints

       else testFound

         ! Family name not present. Processor 0 prints a warning.

         if(myID == 0) then

           print "(a)",  "#"
           print "(a)",  "#*==================== !!! Warning !!! &
                         &======================"
           print "(3a)", "#* Unknown boundary family, ", &
                         trim(famName),                  &
                         ", encountered in the input file"
           print "(a)",  "#* Information is ignored."
           print "(a)",  "#*=====================================&
                         &======================"
           print "(a)", "#"

         endif

       endif testFound

!      Format statements.

 100   format("Wrong format for overwriting boundary condition &
              &info in parameter file for family",1X,A)
 101   format("Number of points specified is larger than the actual &
              &data present for overwriting boundary condition info &
              &in parameter file for family",1X,A)

       end subroutine overwriteBoundaryInfo

!========================================================================

       logical function dataIsDirichlet(varName)
!
!      ******************************************************************
!      *                                                                *
!      * DataIsDirichlet determines whether or not the data given by    *
!      * varName is dirichlet data. If not then it is assumed to be     *
!      * neumann data.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(in) :: varName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize dataIsDirichlet to .true., because there are many
       ! more dirichlet names than neumann names.

       dataIsDirichlet = .true.

       ! Set dataIsDirichlet to .false. for the neumann type of data.
       ! The only one I can think of is heat flux, but i don't know
       ! how that is specified at the moment.

       end function dataIsDirichlet

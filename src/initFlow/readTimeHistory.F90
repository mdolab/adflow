!
!      ******************************************************************
!      *                                                                *
!      * File:          readTimeHistory.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-20-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTimeHistory(fileIDs)
!
!      ******************************************************************
!      *                                                                *
!      * readTimeHistory attempts to read the time history of an        *
!      * unsteady computation from the given cgns restart file.         *
!      * If present it will be stored in the arrays timeArray and       *
!      * timeDataArray, for which memory is allocated.                  *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use su_cgns
       use inputUnsteady
       use monitor
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer, dimension(nSolsRead), intent(in) :: fileIDs

#ifdef USE_NO_CGNS

       call terminate("readTimeHistory", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS, dummyInt
       integer :: i, nConv, nDim, nSize

       integer(kind=intType) :: j, ii, nn

       integer(kind=intType), dimension(:), allocatable :: ind

       character(len=maxCGNSNameLen) :: cgnsName
       character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
                                                 convNames, tmpNames

       logical :: allConvInfo
!
!      Function definitions.
!
       integer               :: setCGNSRealType
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the file ID and the base a bit easier. Note that the time
       ! history only needs to be present in the first solution file.

       cgnsInd  = fileIDs(1)
       cgnsBase = 1

       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Check if the time history is present by trying to read it.

       call cg_biter_read_f(cgnsInd, cgnsBase, cgnsName, &
                            dummyInt, ierr)
       if(ierr /= all_ok) then

         ! No time history present. Set nTimeStepsRestart and
         ! timeUnsteadyRestart to zero, allocate the memory for the
         ! time history of the monitoring variables, print a warning
         ! and return.

         nTimeStepsRestart   = 0
         timeUnsteadyRestart = zero

         call allocTimeArrays(nTimeStepsFine)

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# No time history found in restart file."
         print "(a)", "# Starting at timestep 1 on the finest level."
         print "(a)", "#"

         return

       endif

       ! Store the number of old time levels.

       nTimeStepsRestart = dummyInt

       ! Go to the place in the cgns file where the time history
       ! should be stored.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                      "BaseIterativeData_t", 1, "end")
       if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
                        "Something wrong when calling cg_goto_f")

       ! Find out how many convergence variables are stored.

       call cg_narrays_f(nConv, ierr)
       if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
                        "Something wrong when calling cg_narrays_f")

       ! Allocate the memory for convNames, tmpNames and ind.

       allocate(convNames(nConv), tmpNames(nConv), ind(nConv), &
                stat=ierr)
       if(ierr /= 0)                       &
         call terminate("readTimeHistory", &
                        "Memory allocation failure for convNames, etc.")

       ! Read the names of the convergence variables. Store them in
       ! tmpNames as well. Furthermore check the dimension of the
       ! data stored.

       do i=1,nConv
         call cg_array_info_f(i, convNames(i), dummyInt, nDim, &
                              nSize, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readConvHistory", &
                          "Something wrong when calling cg_array_info_f")

         if(nDim /= 1) then
           print "(a)", "#"
           print "(a)", "#                 Warning"
           print 100, trim(convNames(i))
           print "(a)", "# Information is ignored."
           print "(a)", "#"
 100       format("# Dimension of time history for",1X,A,1X, &
                  "is not 1.")

           ! Screw up the string such that it does not correspond to
           ! a legal name. It is appended, because it is important that
           ! all strings differ.

           convNames(i) = convNames(i)//"#$@&^!#$%!"
         endif

         if(nSize /= nTimeStepsRestart) then
           print "(a)", "#"
           print "(a)", "#                 Warning"
           print 110, trim(convNames(i))
           print "(a)", "# Displayed information might be incorrect."
           print "(a)", "#"
 110       format("# Inconsistent time history for",1X,A,".")
         endif

         ! Copy the name in tmpNames for the sorting.

         tmpNames(i) = convNames(i)
       enddo

       ! Sort convNames in increasing order.

       nn = nConv
       call qsortStrings(convNames, nn)

       ! Find the numbers for the just sorted convergence names.

       do i=1,nConv
         ii      = bsearchStrings(tmpNames(i), convNames, nn)
         ind(ii) = i
       enddo

       ! Find out whether the old time values are present.
       ! If not the time history stored will be ignored.

       ii = bsearchStrings(cgnsTimeValue, convNames, nn)
       if(ii == 0) then
         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# No time values found in the time history &
                      &in the restart file."
         print "(a)", "# The rest of the time history is ignored."
         print "(a)", "# Starting at timestep 1 on the finest level."
         print "(a)", "#"

         ! Set nTimeStepsRestart and timeUnsteadyRestart to 0.

         nTimeStepsRestart   = 0
         timeUnsteadyRestart = zero
       endif

       ! Determine the total number of time levels and allocate the
       ! memory for the time history arrays.

       j = nTimeStepsRestart + nTimeStepsFine
       call allocTimeArrays(j)

       ! If the time values are not present, jump to the place where
       ! the memory of the variables used in this routine is released.

       if(ii == 0) goto 99

       ! Read the time values.

       i = ind(ii)
       call cg_array_read_as_f(i, realTypeCGNS, timeArray, ierr)
       if(ierr /= all_ok)                  &
         call terminate("readTimeHistory", &
                        "Something wrong when calling &
                        &cg_array_read_as_f")

       ! Set the value of timeUnsteadyRestart to the last value in
       ! timeArray.

       timeUnsteadyRestart = timeArray(nTimeStepsRestart)

       ! Initialize allConvInfo to .true. and perform the loop over
       ! the number of monitoring variables.

       allConvInfo = .true.

       do j=1,nMon

         ! Search for the monitoring name in the sorted
         ! convergence names present in the restart file.

         ii = bsearchStrings(monNames(j), convNames, nn)

         ! Check if the name was found.

         if(ii == 0) then

           ! Name not present in the restart file. Set allConvInfo
           ! to .false. and the corresponding entries in timeDataArray
           ! to zero

           allConvInfo = .false.
           do i=1,nTimeStepsRestart
             timeDataArray(i,j) = zero
           enddo

         else

           ! Name is present in the restart file. Read the corresponding
           ! time history.

           i = ind(ii)
           call cg_array_read_as_f(i, realTypeCGNS, &
                                   timeDataArray(1,j), ierr)
           if(ierr /= all_ok)                  &
             call terminate("readTimeHistory", &
                            "Something wrong when calling &
                            &cg_array_read_as_f")
         endif

       enddo

       ! Print a warning in case not all the time history could
       ! be retrieved from the restart file.

       if(.not. allConvInfo) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# Not all the time history could be &
                      &retrieved from the restart file."
         print "(a)", "# Missing information is initialized to zero."
         print "(a)", "#"

       endif

 99    continue

       ! Release the memory of the variables allocated in this routine.

       deallocate(convNames, tmpNames, ind, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("readTimeHistory", &
                        "Deallocation error for convNames, etc.")

#endif

       end subroutine readTimeHistory

!
!      ******************************************************************
!      *                                                                *
!      * File:          readConvHistory.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-24-2003                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readConvHistory(fileIDs)
!
!      ******************************************************************
!      *                                                                *
!      * readConvHistory attempts to read the convergence history       *
!      * from the CGNS restart file(s). If present it will be stored    *
!      * in the array convArray, for which memory is allocated.         *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use su_cgns
       use inputIO
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use monitor
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer, dimension(nSolsRead), intent(in) :: fileIDs

#ifdef USE_NO_CGNS
       call terminate("readConvHistory", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS, dummyInt
       integer :: i, nConv, nConvHistories, nDim, nSize, sol, base

       integer(kind=intType) :: j, ii, nn

       integer(kind=intType), dimension(:), allocatable :: ind

       character(len=1024) :: string
       character(len=8)    :: intString

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
       ! Return immediately if the convergence history (of the inner
       ! iterations) does not need to be stored. This logical can
       ! only be .false. for an unsteady computation.

       if(.not. storeConvInnerIter) then
         nIterOld = 0
         return
       endif

       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Store the total number of iterations to be performed in nn.
       ! For unsteady computations the entire convergence history
       ! will be stored.

       nn = nsgStartup + nCycles
       if(equationMode == unsteady) nn = nTimeStepsFine*nn

       ! Determine the number of convergence histories to be read.
       ! For the spectral mode it is attempted to read the convergence
       ! histories of all spectral solutions. If this is not possible,
       ! only 1 convergence history will be read, which will then be
       ! copied to all solutions.

       nConvHistories = 1
       if(equationMode == timeSpectral) then
         nConvHistories = nTimeIntervalsSpectral
         if(nSolsRead /= nConvHistories) nConvHistories = 1
       endif

       ! Loop over the number of convergence histories to determine
       ! nIterOld.

       nIterOld = 0
       do sol=1,nConvHistories

         ! Store the file index and base a bit easier and go to the
         ! position in the CGNS file where the convergence history
         ! should be stored.

         cgnsInd = fileIDs(sol)
         base    = 1

         call cg_goto_f(cgnsInd, base, ierr, "end")
         if(ierr /= all_ok)                  &
           call terminate("readConvHistory", &
                          "Something wrong when calling cg_goto_f")

         ! Check if a convergence history is present by trying to read
         ! a convergence history node.

         call cg_convergence_read_f(i, string, ierr)

         if(ierr /= all_ok) then

           ! No convergence info present. Determine the situation
           ! we have here. If this is not the 1st convergence history,
           ! i.e. spectral mode, simply set nConvHistories to 1.
           ! In all cases, exit the do-loop.

           if(sol > 1) nConvHistories = 1
           exit

         endif

         ! Convergence history is present. If this is the first solution
         ! set nIterOld to the value just read. Otherwise, spectral mode,
         ! compare the given number and if it is not equal to the 1st
         ! base set nConvHistories to 1 and exit the do loop.

         if(sol == 1) then
           nIterOld = i
         else if(nIterOld /= i) then
           nConvHistories = 1
           exit
         endif

       enddo

       ! If no convergence history is present, allocate the memory for
       ! the convergence arrays, print a warning (if this a steady or a
       ! spectral computation) and return.

       if(nIterOld == 0) then
         call allocConvArrays(nn)

         if(equationMode == steady .or. &
            equationMode == timeSpectral) then
           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# No convergence info found in restart file."
           print "(a)", "# Starting at iteration 0 on the finest level."
           print "(a)", "#"
         endif

         return
       endif

       ! Substract 1 from nIterOld, because the numbering starts at 0
       ! in this code. Allocate the memory for the convergence arrays.

       nIterOld = nIterOld - 1
       nn       = nn + nIterOld
       call allocConvArrays(nn)

       ! Loop over the number of solutions to read the convergence.

       solLoop: do sol=1,nConvHistories

         ! Go to the node where the convergence history is stored.
 
         cgnsInd = fileIDs(sol)
         base    = 1
         call cg_goto_f(cgnsInd, base, ierr, &
                        "ConvergenceHistory_t", 1, "end")
         if(ierr /= all_ok)                  &
           call terminate("readConvHistory", &
                          "Something wrong when calling cg_goto_f")

         ! Find out how many convergence variables are stored.

         call cg_narrays_f(nConv, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readConvHistory", &
                          "Something wrong when calling cg_narrays_f")

         ! Allocate the memory for convNames, tmpNames and ind.

         allocate(convNames(nConv), tmpNames(nConv), ind(nConv), &
                  stat=ierr)
         if(ierr /= 0)                       &
           call terminate("readConvHistory", &
                          "Memory allocation failure for convNames, &
                          &etc.")

         ! Read the names of the convergence variables. Store them in
         ! tmpNames as well. Furthermore check the dimension of the
         ! data stored.

         do i=1,nConv
           call cg_array_info_f(i, convNames(i), dummyInt, nDim, &
                                nSize, ierr)
           if(ierr /= all_ok)                  &
             call terminate("readConvHistory", &
                            "Something wrong when calling &
                            &cg_array_info_f")

           if(nDim /= 1) then
             print "(a)", "#"
             print "(a)", "#                 Warning"
             print 100, trim(convNames(i))
             print "(a)", "# Information is ignored."
             print "(a)", "#"
 100         format("# Dimension of convergence info for",1X,A,1X, &
                    "is not 1.")

             ! Screw up the string such that it does not correspond to
             ! a legal name. It is appended, because it is important that
             ! all strings differ.

             convNames(i) = convNames(i)//"#$@&^!#$%!"
           endif

           if((nSize-1) /= nIterOld) then
             print "(a)", "#"
             print "(a)", "#                 Warning"
             print 110, trim(convNames(i))
             print "(a)", "# Displayed information might be incorrect."
             print "(a)", "#"
 110         format("# Inconsistent convergence info for",1X,A,".")
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

         ! Check for the L2 norm of the density residual. This must be
         ! present, otherwise the reading of the convergence file does
         ! not make sense and the rest of the convergence info will be
         ! ignored.

         ii = bsearchStrings(cgnsL2resRho, convNames, nn)

         densityResFound: if(ii == 0) then

           ! Determine the situation we have here.

           if(sol == 1) then

             ! 1st solution. This means that reading the convergence info
             ! is completely pointless. Print a warning for the steady
             ! mode and the time spectral mode.

             if(equationMode == steady .or. &
                equationMode == timeSpectral) then

               print "(a)", "#"
               print "(a)", "#                 Warning"
               print "(a)", "# No convergence info for the density &
                            &residual found in the restart file."
               print "(a)", "# The rest of the convergence history is &
                            &ignored."
               print "(a)", "# Starting at iteration 0 on the &
                            &finest level."
               print "(a)", "#"

             endif

             ! Set nIterOld to 0 and nConvHistories to 1 and jump to
             ! the place where the memory of the variables used in the
             ! loop over the solutions is released. The setting of
             ! nConvHistories is only there to exit the current loop
             ! over the solutions.

             nIterOld = 0
             nConvHistories = 1
             goto 99

           else

             ! No density residuals found for higher solutions.
             ! Set nConvHistories to 1 such that the loop over the
             ! solutions is not executed any further and jump to the
             ! place where the memory used in this loop is released.

             nConvHistories = 1
             goto 99

           endif

         else densityResFound

           ! Density residual is present. Print the restarting info
           ! if this is solution 1 and this is either a steady or a time
           ! spectral computation.

           if(sol == 1 .and. (equationMode == steady .or. &
                              equationMode == timeSpectral)) then

             write(intString,"(i7)") nIterOld
             intString = adjustl(intString)

             print "(a)", "#"
             print 120, trim(intString)
             print "(a)", "#"
 120         format("# Restarting at iteration",1X,A,".")
           endif

         endif densityResFound

         ! Initialize allConvInfo to .true.

         allConvInfo = .true.

         ! Loop over all the variables that must be monitored.

         do j=1,nMon

           ! Search for the monitoring name in the sorted
           ! convergence names present in the restart file.

           ii = bsearchStrings(monNames(j), convNames, nn)

           ! Check if the name was found.

           if(ii == 0) then

             ! Name not present in the restart file. Set allConvInfo
             ! to .false. and the corresponding entries in convArray
             ! either to zero or to the values of solution 1.

             allConvInfo = .false.

             if(sol == 1) then
               do i=0,nIterOld
                 convArray(i,sol,j) = zero
               enddo
             else
               do i=0,nIterOld
                 convArray(i,sol,j) = convArray(i,1,j)
               enddo
             endif

           else

             ! Name is present in the restart file. Read the corresponding
             ! convergence history.

             i = ind(ii)
             call cg_array_read_as_f(i, realTypeCGNS, &
                                     convArray(0,sol,j), ierr)
             if(ierr /= all_ok)                  &
               call terminate("readConvHistory", &
                              "Something wrong when calling &
                              &cg_array_read_as_f")
           endif

         enddo

         ! Print a warning in case not all the convergence info could
         ! be retrieved from the restart file. The error message depends
         ! on the situation.

         if(.not. allConvInfo) then

           if(equationMode == timeSpectral .and. sol > 1) then

             ! Not the first spectral solution in time spectral mode.
             ! Inform that the data which is not present is set to the
             ! data of the 1st spectral solution.

             write(intString,"(i7)") sol
             intString = adjustl(intString)

             print "(a)", "#"
             print 130
             print 140, trim(intString)
             print 160
             print "(a)", "#"

           else

             ! 1st solution. Missing information was set to zero.

             print "(a)", "#"
             print 130
             print 150
             print 170
             print "(a)", "#"

           endif

         endif

         ! Format statements for the missing convergence info.

 130     format("#                 Warning")
 140     format("# Not all the convergence info could be retrieved &
                &from the restart file for solution",1X,A,".")
 150     format("# Not all the convergence info could be retrieved &
                &from the restart file.")
 160     format("# Missing information is copied from the first &
                &solution.")
 170     format("# Missing information is initialized to zero.")


         ! Continue statement needed when no convergence information
         ! of the density residual was present.

 99      continue

         ! Release the memory of the variables allocated in this routine.

         deallocate(convNames, tmpNames, ind, stat=ierr)
         if(ierr /= 0)                       &
           call terminate("readConvHistory", &
                          "Deallocation error for convNames, etc.")

       enddo solLoop

       ! Check for the spectral mode if the number of convergence
       ! histories equals the number of spectral solutions. If not
       ! print a warning and initialize them

       if(equationMode == timeSpectral .and. &
          nConvHistories == 1 .and. nIterOld > 0) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# Inconsistent number of convergence histories &
                      &found for the time spectral method."
         print "(a)", "# Copying the data of the first solution to &
                      &all spectral solutions."
         print "(a)", "#"

         ! Loop over the number of spectral solutions, starting at 2,
         ! and copy the convergence data from the 1st solution.

         do sol=2,nTimeIntervalsSpectral
           do j=1,nMon
             do i=0,nIterOld
               convArray(i,sol,j) = convArray(i,1,j)
             enddo
           enddo
         enddo

       endif

#endif

       end subroutine readConvHistory

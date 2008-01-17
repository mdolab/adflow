!
!      ******************************************************************
!      *                                                                *
!      * File:          extractMGInfo.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-13-2002                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extractMGInfo
!
!      ******************************************************************
!      *                                                                *
!      * extractMgInfo creates the integer array cycleStrategy from     *
!      * the string describing the multigrid strategy. This string      *
!      * either contains a predefined strategy, like sg, 2v, 4w, etc.,  *
!      * or a combination of -1's, 0's and 1's, which defines a user    *
!      * defined strategy. The integers -1, 0 and 1 have the following  *
!      * meaning:  0 -> perform an iteration step on the current grid.  *
!      *           1 -> go to next coarser grid.                        *
!      *          -1 -> go to next finer grid.                          *
!      * For a valid cycling strategy the sum of the elements of the    *
!      * array should be 0.                                             *
!      *                                                                *
!      ******************************************************************
!
       use inputIO
       use inputIteration
       use inputPhysics
       use inputUnsteady
       use communication
       use constants
       use localMG
       implicit none
!
!      Local variables
!
       integer                      :: stringLen, error
       integer(kind=intType)        :: i, ii, nMinus, nn
       character (len=maxStringLen) :: errorMessage
!
!      Function definitions.
!
       logical          :: digitsOnlyInString
       integer(intType) :: computeNstepsWcycle
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! For an unsteady computation using explicit Runge-Kutta schemes
       ! overrule mgDescription to sg.

       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) mgDescription = "sg"

       ! Create a lower case version of mgDescription and determine the
       ! length of the string. Note that if the string contains a user
       ! defined cycling strategy the contents (0's, 1's, -1's) is not
       ! changed by the call to convertToLowerCase

       mgDescription = trim(mgDescription)
       call convertToLowerCase(mgDescription)
       stringLen = len_trim(mgDescription)

       ! Check for predefined cycling strategies.

       if(mgDescription == "sg") then

         ! Single grid computation. Set the values for the parameters
         ! nMGSteps (number of steps in the cycle strategy) and
         ! cycleStrategy (the cycle strategy itself).

         nMGSteps = 1
         allocate(cycleStrategy(1), stat=error)
         if(error /= 0) call terminate("extractMgInfo", &
                                       "allocation error for 1 integer")
         cycleStrategy(1) = 0

         ! Set the number of grid levels needed by the multigrid to 1.

         nMGLevels = 1

       else if(mgDescription(stringLen:stringLen) == "v") then

         ! Must be a v-cycle. The rest of the string should only contain
         ! digits. Check this.

         if(.not. digitsOnlyInString(mgDescription(:stringLen-1))) then
           write(errorMessage,*) "Invalid cycle strategy, ", &
                                  mgDescription(:stringLen), ", specified"
           if(myID == 0) call terminate("extractMgInfo", errorMessage)
         endif

         ! Read the number of levels in the cycle.

         read(mgDescription(:stringLen-1),*) nMGLevels

         ! Determine the number of steps in cycleStrategy and allocate
         ! the memory for it.

         nMGSteps = 4*nMGLevels - 4
         allocate(cycleStrategy(nMGSteps), stat=error)
         if(error /= 0) then
           write(errorMessage,*) "Allocation error for", nMGSteps, &
                                  "integers for the v-cycle ",       &
                                  mgDescription(:stringLen)
           call terminate("extractMgInfo", errorMessage)
         endif

         ! Set the values of cycleStrategy.

         ii = 1
         do i=1,(nMGLevels-1)
           cycleStrategy(ii)   = 0
           cycleStrategy(ii+1) = 1
           ii = ii+2
         enddo

         do i=1,(nMGLevels-1)
           cycleStrategy(ii)   =  0
           cycleStrategy(ii+1) = -1
           ii = ii+2
         enddo

       else if(mgDescription(stringLen:stringLen) == "w") then

         ! Must be a w-cycle. The rest of the string should only contain
         ! digits. Check this.

         if(.not. digitsOnlyInString(mgDescription(:stringLen-1))) then
           write(errorMessage,*) "Invalid cycle strategy, ", &
                                  mgDescription(:stringLen), ", specified"
           if(myID == 0) call terminate("extractMgInfo", errorMessage)
         endif

         ! Read the number of levels in the cycle.

         read(mgDescription(:stringLen-1),*) nMGLevels

         ! Determine the number of steps in cycleStrategy and allocate
         ! the memory for it.

         nMGSteps = computeNstepsWcycle(nMGLevels)
         allocate(cycleStrategy(nMGSteps), stat=error)
         if(error /= 0) then
           write(errorMessage,*) "Allocation error for", nMGSteps, &
                                  "integers for the w-cycle ",       &
                                  mgDescription(:stringLen)
           call terminate("extractMgInfo", errorMessage)
         endif

         ! Set the values of cycleStrategy.

         ii = 1
         call setEntriesWcycle(ii, nMGLevels)

       else

         ! The string must be a collection of 0's, -1's and 1's to
         ! describe the cycle strategy. Get rid of the internal spaces
         ! first and determine the amount of -'s.

         ii = 0
         nMinus = 0
         do i=1,stringLen
           if(mgDescription(i:i) /= " ") then
             ii = ii+1
             mgDescription(ii:ii) = mgDescription(i:i)
             if(mgDescription(ii:ii) == "-") nMinus = nMinus+1
           endif
         enddo
         stringLen = ii

         ! Determine the number of steps in the cycle strategy and
         ! allocate the memory for it.

         nMGSteps = ii - nMinus
         allocate(cycleStrategy(nMGSteps), stat=error)
         if(error /= 0) then
           write(errorMessage,*) "Allocation error for", nMGSteps, &
                                  "integers for the cycle strategy"
           call terminate("extractMgInfo", errorMessage)
         endif

         ! Determine the entries for cycleStrategy.

         i = 1
         nn = 1
         do
           ii = i
           if(mgDescription(i:i) == "-") i = i+1

           ! Determine the case we are having here.

           select case (mgDescription(ii:i))
             case ("0")
               cycleStrategy(nn) = 0
             case ("1")
               cycleStrategy(nn) = 1
             case ("-1")
               cycleStrategy(nn) = -1
             case default
               write(errorMessage, *) "Invalid character, ", &
                                      mgDescription(ii:i), &
                                      ", in the string describing &
                                      &cycling strategy"
               if(myID == 0) call terminate("extractMgInfo", errorMessage)
           end select

           ! Update i and nn

           i  = i  + 1
           nn = nn + 1

           ! Exit the do loop in case i is larger than stringLen.

           if(i > stringLen) exit
         enddo

         ! Check if the string specified is valid and determine the
         ! maximum grid level needed in the cycle.

         nn = 0
         nMGLevels = 0
         do i=1,nMGSteps
           nn = nn + cycleStrategy(i)
           nMGLevels = max(nn, nMGLevels)
         enddo
         nMGLevels = nMGLevels + 1

         if(nn /= 0 .and. myID == 0) &
           call terminate("extractMgInfo", &
                     "sum of coefficients in cycle strategy is not 0")
       endif

       ! Correct the value of mgStartlevel in case a nonpositive number
       ! has been specified. In that case it is set to -1, the default
       ! value.

       if(mgStartlevel <= 0) mgStartlevel = -1

       ! Determine the value of mgStartlevel. This parameter might be
       ! specified in the parameter file, but it is checked here for
       ! consistency. If mgStartlevel has not been specified it is
       ! either set to the coarsest level in the mg cycle (starting
       ! from free stream) or to the finest level (restart).

       if(mgStartlevel == -1) then

         ! Value has not been specified. Default value is set, see
         ! the comments above.

         if( restart ) then
           mgStartlevel = 1
         else
           mgStartlevel = nMGLevels
         endif

       else

         ! Value has been specified. Correct this value for a restart
         ! in case not the finest grid has been specified. Processor 0
         ! prints a warning.

         if( restart ) then

           if(mgStartlevel /= 1 .and. myID == 0) then
             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print 100, mgStartlevel
 100         format("#* Multigrid start level specified is ",i1,".")
             print "(a)", "#* This is in conflict with the restart."
             print "(a)", "#* Therefore Multigrid start level is &
                          &set to 1."
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"
           endif

           ! Correct the start level.

           mgStartlevel = 1

         endif

       endif

       end subroutine extractMGInfo

!      ==================================================================

       logical function digitsOnlyInString(string)
!
!      ******************************************************************
!      *                                                                *
!      * digitsOnlyInString checks whether the given string contains    *
!      * digits only or if other character types are present. In the    *
!      * former case the function returns .True., otherwise .False.     *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine argument                                              *
!
       character (len=*), intent(in) :: string
!
!      Local variables
!
       integer :: i, stringLen
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize digitsOnlyInString to .True.

       digitsOnlyInString = .true.

       ! Determine the length of the string.

       stringLen = len_trim(string)

       ! Loop over the elements of the string and check if they are digits.

       do i=1,stringLen
         if(string(i:i) < "0" .or. string(i:i) > "9") &
           digitsOnlyInString = .false.
       enddo

       end function digitsOnlyInString

!      ==================================================================

       recursive function computeNstepsWcycle(nLevels) result(nSteps)
!
!      ******************************************************************
!      *                                                                *
!      * computeNstepsWcycle is recursive function, which determines    *
!      * the number of entries of a w-cycle of a given level.           *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use communication
       implicit none
!
!      Result variable
!
       integer(kind=intType) :: nSteps
!
!      Function argument
!
       integer(kind=intType), intent(in) :: nLevels
!
!      Local variables
!
       character (len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the case we are having here. For nLevels is less
       ! than 2 an error message is printed, in case nLevels is 2
       ! the recursion is broken and otherwise a recursive call is made.
 
       if(nLevels < 2) then
         write(errorMessage,*) "Wrong value of nLevels", nLevels
         if(myID == 0) call terminate("computeNstepsWcycle", errorMessage)
       else if(nLevels == 2) then
         nSteps = 4
       else
         nSteps = 4 + 2*computeNstepsWcycle(nLevels-1)
       endif

       end function computeNstepsWcycle

!      ==================================================================

       recursive subroutine setEntriesWcycle(counter, nLevels)
!
!      ******************************************************************
!      *                                                                *
!      * setEntriesWcycle is a recursive subroutine, which actually     *
!      * fills the entries of cycleStrategy for a w-cycle.              *
!      *                                                                *
!      ******************************************************************
!
       use inputIteration
       use communication
       use constants
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType), intent(inout) :: counter
       integer(kind=intType), intent(in)    :: nLevels
!
!      Local variables
!
       character (len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the case we are having here. For nLevels is less
       ! than 2 an error message is printed, in case nLevels is 2
       ! the recursion is broken and otherwise a recursive call is made.
 
       if(nLevels < 2) then

         write(errorMessage,*) "Wrong value of nLevels", nLevels
         if(myID == 0) call terminate("setEntriesWcycle", errorMessage)

       else if(nLevels == 2) then

         cycleStrategy(counter)   =  0
         cycleStrategy(counter+1) =  1
         cycleStrategy(counter+2) =  0
         cycleStrategy(counter+3) = -1

         counter = counter + 4

       else

         cycleStrategy(counter)   =  0
         cycleStrategy(counter+1) =  1
         counter = counter + 2

         call setEntriesWcycle(counter, nLevels-1)
         call setEntriesWcycle(counter, nLevels-1)

         cycleStrategy(counter)   =  0
         cycleStrategy(counter+1) = -1
         counter = counter + 2

       endif

       end subroutine setEntriesWcycle

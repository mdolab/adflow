!
!      ******************************************************************
!      *                                                                *
!      * File:          getSortedVarNumbers.F90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-18-2003                                      *
!      * Last modified: 10-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine getSortedVarNumbers
!
!      ******************************************************************
!      *                                                                *
!      * getSortedVarNumbers reads the names of variables stored in     *
!      * the given solution node of the cgns file, indicated by         *
!      * cgnsInd, cgnsBase and cgnsZone. Afterwards the variable        *
!      * names are sorted in increasing order, such that they can be    *
!      * used in a binary search. Their original variable number and    *
!      * type is stored.                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("getSortedVarNumbers", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use constants
       use su_cgns
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: i, ierr
       integer, dimension(:), allocatable :: tmpTypes

       integer(kind=intType) :: nn, ii

       integer(kind=intType), dimension(:), allocatable :: varNumbers

       character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
                                                                tmpNames
!
!      Function definition.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of solution variables stored.

       call cg_nfields_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                         nVar, ierr)
       if(ierr /= all_ok)                      &
         call terminate("getSortedVarNumbers", &
                        "Something wrong when calling cg_nfield_f")

       ! Allocate the memory for varnames, vartypes and varnumber

       allocate(varNames(nVar), varTypes(nVar), varNumbers(nVar), &
                stat=ierr)
       if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
                        "Memory allocation failure for varNames, etc.")

       ! Loop over the number of variables and store their names and
       ! types.

       do i=1,nVar
         call cg_field_info_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                              i, varTypes(i), varNames(i), ierr)
         if(ierr /= 0)                           &
           call terminate("getSortedVarNumbers", &
                          "Something wrong when calling cg_field_info_f")
       enddo

       ! Allocate the memory for tmpTypes and tmpNames and initialize
       ! their values.

       allocate(tmpTypes(nVar), tmpNames(nVar), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
                        "Memory allocation failure for tmp variables")

       do i=1,nVar
         tmpTypes(i) = varTypes(i)
         tmpNames(i) = varNames(i)
       enddo

       ! Sort varNames in increasing order.

       nn = nVar
       call qsortStrings(varNames, nn)

       ! Initialize varNumbers to -1. This serves as a check during
       ! the search.

       do i=1,nVar
         varNumbers(i) = -1
       enddo

       ! Find the original types and numbers for the just sorted
       ! variable names.

       do i=1,nVar
         ii = bsearchStrings(tmpNames(i), varNames, nn)

         ! Check if the variable number is not already taken. If this is
         ! the case, this means that two identical var names are present.

         if(varNumbers(ii) /= -1)                &
           call terminate("getSortedVarNumbers", &
                          "Error occurs only when two identical &
                          &variable names are present")

         ! And set the variable number and type.

         varNumbers(ii) = i
         varTypes(ii)   = tmpTypes(i)
       enddo

       ! Release the memory of varNumbers, tmpNames and tmpTypes.

       deallocate(varNumbers, tmpTypes, tmpNames, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
                        "Deallocation error for tmp variables")

#endif

       end subroutine getSortedVarNumbers

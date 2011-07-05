!
!      ******************************************************************
!      *                                                                *
!      * File:          scaleFactors.F90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-28-2003                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine scaleFactors(fileIDs)
!
!      ******************************************************************
!      *                                                                *
!      * scaleFactors determines the scale factors for the density,     *
!      * pressure and velocity from either the reference state in the   *
!      * given cgns file or they are simply set to 1.0; the latter      *
!      * occurs if the input parameter checkRestartSol is .false.       *
!      * If no reference state is present checkRestartSol is .true.     *
!      * An error message will be printed and the program terminates.   *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use communication
       use flowVarRefState
       use inputIO
       use su_cgns
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer, dimension(nSolsRead), intent(in) :: fileIDs

#ifdef USE_NO_CGNS
       call terminate("scaleFactors", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS, typeCGNS
       integer :: i, nsize, nDim, nRef

       integer(kind=intType) :: ii, nn

       integer(kind=intType), dimension(:), allocatable :: ind

       real(kind=cgnsRealType) :: tmpScale

       character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
                                                     refNames, tmpNames

       logical :: muScalePresent
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
       ! Store the file ID and the base a bit easier. Note that the
       ! reference state only needs to be present in the first file.

       cgnsInd  = fileIDs(1)
       cgnsBase = 1

       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Initialize the scale factors to 1.0, i.e. assume that the
       ! correct non-dimensional solution is stored in the restart file.

       rhoScale = one
       velScale = one
       pScale   = one
       muScale  = one

       ! Go to the place in the cgns file where the reference state
       ! should be stored.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "end")
       if(ierr /= all_ok)                &
         call terminate("scaleFactors", &
                        "Something wrong when calling cg_goto_f")

       ! Try to determine the size of the string describing the reference
       ! state. If the error code is not all_ok, this means that no
       ! reference node is present.

       call cg_state_size_f(nsize, ierr)
       if(ierr /= all_ok) then

         ! Reference state does not exist. Check if the restart solution
         ! must be checked. If not, return; otherwise print an error
         ! message and terminate the execution. This error message is
         ! only printed by processor 0 to avoid a messy output.

         if(.not. checkRestartSol) return

         if(myId == 0)                    &
           call terminate("scaleFactors", &
                          "Reference state not presented in restart &
                          &file. Scaling factors cannot be determined.")

         ! The other processors will wait until they are killed.

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Go to the reference state node.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                      "ReferenceState_t", 1, "end")
       if(ierr /= all_ok)                &
         call terminate("scaleFactors", &
                        "Something wrong when calling cg_goto_f")

       ! Found out how many reference variables are stored.

       call cg_narrays_f(nRef, ierr)
       if(ierr /= all_ok)               &
         call terminate("scaleFactors", &
                        "Something wrong when calling cg_narrays_f")

       ! Allocate the memory for refNames, tmpNames and ind.

       allocate(refNames(nRef), tmpNames(nRef), ind(nRef), &
                stat=ierr)
       if(ierr /= 0)                     &
         call terminate("scaleFactors", &
                        "Memory allocation failure for refNames, etc.")

       ! Read the names of the reference variables. Store them in
       ! tmpNames as well.

       do i=1,nRef
         call cg_array_info_f(i, refNames(i), typeCGNS, nDim, &
                              nsize, ierr)
         if(ierr /= all_ok)                &
           call terminate("scaleFactors", &
                          "Something wrong when calling cg_array_info_f")

         ! Check the dimension and the size of the array.
         ! Both should be 1. If not, screw up the name such that it
         ! will never be found in the search later on.

         if(nDim /= 1 .or. nsize /= 1) &
           refNames(i) = refNames(i)//"#$@&^!#$%!"

         ! And copy it in tmpNames.

         tmpNames(i) = refNames(i)
       enddo

       ! Sort refNames in increasing order.

       nn = nRef
       call qsortStrings(refNames, nn)

       ! Find the numbers for the just sorted reference names.

       do i=1,nRef
         ii      = bsearchStrings(tmpNames(i), refNames, nn)
         ind(ii) = i
       enddo

       ! Determine the scale factors if these must be determined.

       if( checkRestartSol ) then

         ! Read the reference density from the restart file.

         ii = bsearchStrings(cgnsDensity, refNames, nn)
         if(ii == 0) then
           if(myId == 0)                    &
             call terminate("scaleFactors", &
                            "No reference density found in restart file")

           ! The other processors will wait until they are killed.

           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         i = ind(ii)
         call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
         if(ierr /= all_ok)               &
           call terminate("scaleFactors", &
                          "Something wrong when calling &
                          &cg_array_read_as_f")
         rhoScale = tmpScale

         ! Read the reference pressure from the restart file.

         ii = bsearchStrings(cgnsPressure, refNames, nn)
         if(ii == 0) then
           if(myId == 0)                    &
             call terminate("scaleFactors", &
                            "No reference pressure found in &
                            &restart file")

           ! The other processors will wait until they are killed.

           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         i = ind(ii)
         call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
         if(ierr /= all_ok)                &
           call terminate("scaleFactors", &
                          "Something wrong when calling &
                          &cg_array_read_as_f")
         pScale = tmpScale

         ! Read the reference velocity from the restart file.

         ii = bsearchStrings(cgnsVelocity, refNames, nn)
         if(ii == 0) then
           if(myId == 0)                    &
             call terminate("scaleFactors", &
                            "No reference velocity found in &
                            &restart file")

           ! The other processors will wait until they are killed.

           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         i = ind(ii)
         call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
         if(ierr /= all_ok)                &
           call terminate("scaleFactors", &
                          "Something wrong when calling &
                          &cg_array_read_as_f")
         velScale = tmpScale

         ! Set muScalePresent to .true. to indicate that it is present
         ! and read or construct the reference molecular viscosity.

         muScalePresent = .true.

         ii = bsearchStrings(cgnsViscMol, refNames, nn)
         if(ii > 0) then

           ! Scale is present; read the value.

           i = ind(ii)
           call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
           if(ierr /= all_ok)                &
             call terminate("scaleFactors", &
                            "Something wrong when calling &
                            &cg_array_read_as_f")
           muScale = tmpScale

         else

           ! Try to read the kinematic viscosity.

           ii = bsearchStrings(cgnsViscKin, refNames, nn)
           if(ii > 0) then

             ! Scale is present; read the value and multiply it by the
             ! density.

             i = ind(ii)
             call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
             if(ierr /= all_ok)                &
               call terminate("scaleFactors", &
                              "Something wrong when calling &
                              &cg_array_read_as_f")

             muScale = tmpScale*rhoScale

           else

             ! Final possibility. Try to read the length scale.

             ii = bsearchStrings(cgnsLength, refNames, nn)
             if(ii > 0) then

               ! Scale is present; read the value and create the
               ! reference viscosity.

               i = ind(ii)
               call cg_array_read_as_f(i, realTypeCGNS, tmpScale, ierr)
               if(ierr /= all_ok)               &
                 call terminate("scaleFactors", &
                                "Something wrong when calling &
                                &cg_array_read_as_f")

               muScale = tmpScale*sqrt(pScale*rhoScale)

             else

               ! Set the logical muScalePresent to .false.

               muScalePresent = .false.

             endif
           endif
         endif

         ! Create the correct scaling factors for density, pressure,
         ! velocity and possibly viscosity.

         rhoScale = rhoScale/rhoRef
         pScale   = pScale/pRef
         velScale = velScale/sqrt(pRef/rhoRef)

         if( muScalePresent ) muScale = muScale/muRef

       endif

       ! Release the memory of refNames, tmpNames and ind.

       deallocate(refNames, tmpNames, ind, stat=ierr)
       if(ierr /= 0)                    &
         call terminate("scaleFactors", &
                        "Deallocation error for convNames, etc.")

#endif

       end subroutine scaleFactors

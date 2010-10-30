!
!      ******************************************************************
!      *                                                                *
!      * File:          setChangingOverset.f90                          *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 06-03-2005                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setChangingOverset
!
!      ******************************************************************
!      *                                                                *
!      * setChangingOverset determines the logical value of variable    *
!      * changingOverset which is stored in the module iteration. It is *
!      * only relevant for unsteady or timeSpectral mode. In either,    *
!      * case the value is set to true when the two blocks of any input *
!      * overset connectivities are moving with respect to one another. *
!      * Also, if the grids are deforming, the assumption is made that  *
!      * the overset connectivities are changing in time.               *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputPhysics
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i

       real(kind=realType), dimension(3) :: diff
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the value to false. If the equation mode is steady,
       ! there is no need to continue.

       changingOverset = .false.
       if (equationMode == steady) return

       ! Loop over the domains.

       domains: do nn = 1,cgnsNDom

         ! If the grid is deforming, then it is assumed that overset
         ! changings will be required, but only if there are overset
         ! connectivities. Since this is currently an all or nothing
         ! process, exit the domain loop as soon as one is found.

         if (deforming_Grid .and. cgnsDoms(nn)%nOverset > 0) then
           changingOverset = .true.
           exit
         end if

         ! Loop over the overset connectivities for this block.

         conns: do i = 1,cgnsDoms(nn)%nOverset

           ! Store the index of the donor block easier.

           mm = cgnsDoms(nn)%connOver(i)%donorBlock

           ! If motion has been specified for both this block and the
           ! donor, then compare the rotation centers and rate. As there
           ! could be small rounding differences in the conversion math,
           ! compute the difference to compare. Also, if rotation was
           ! specified for one and not the other, then relative motion
           ! between them is assumed.

           if (cgnsDoms(nn)%rotatingFrameSpecified .and. &
               cgnsDoms(mm)%rotatingFrameSpecified) then

             diff = abs(cgnsDoms(nn)%rotCenter - cgnsDoms(mm)%rotCenter)
             if (any(diff > eps)) changingOverset = .true.

             diff = abs(cgnsDoms(nn)%rotRate - cgnsDoms(mm)%rotRate)
             if (any(diff > eps)) changingOverset = .true.

           else if ((      cgnsDoms(nn)%rotatingFrameSpecified .and. &
                     .not. cgnsDoms(mm)%rotatingFrameSpecified) .or. &
                    (.not. cgnsDoms(nn)%rotatingFrameSpecified .and. &
                           cgnsDoms(mm)%rotatingFrameSpecified)) then

             changingOverset = .true.

           end if
         end do conns
       end do domains

       ! To guarantee that all processors have come up with the same 
       ! logical value, perform a summation reduction. This way the sum
       ! should be 0 for all false, or nProc for all true.

       nn = 0
       if (changingOverset) nn = 1

       call mpi_allreduce(nn, mm, 1, sumb_integer, mpi_sum, &
                          SUmb_comm_world, ierr)

       if (mm /= 0 .and. mm /= nProc) then

         if (myID == 0) &
           call terminate("setChangingOverset", &
                          "There was a discrepancy in the value of &
                          &changingOverset among the processes.")

         call mpi_barrier(SUmb_comm_world, ierr)
       end if

       ! If the code was compiled in stand-alone mode, there is no 
       ! support yet for updating the overset assembly so terminate.

       if (standAloneMode .and. changingOverset) &
         call terminate("setChangingOverset", &
                        "Cannot currently update overset connectivity &
                        &in stand-alone mode")

       end subroutine setChangingOverset

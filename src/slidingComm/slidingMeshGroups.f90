!
!      ******************************************************************
!      *                                                                *
!      * File:          slidingMeshGroups.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine slidingMeshGroups
!
!      ******************************************************************
!      *                                                                *
!      * slidingMeshGroups determines the coloring for the sliding      *
!      * mesh interfaces, such that as many interfaces are treated      *
!      * simultaneously, but every processor works at only one at a     *
!      * time. Furthermore when intersecting sliding mesh interfaces    *
!      * are present, the sequence of treating them is important. The   *
!      * one with the lower id is treated first. This is just a choice  *
!      * used in this code; nothing special.                            *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use interfaceGroups
       use tmpSliding
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, nn, mm, jj, kk, nSliding

       logical :: canBeStored

       integer(kind=intType), dimension(0:nProc-1)    :: flag
       integer(kind=intType), dimension(cgnsNSliding) :: nFaceMax
       integer(kind=intType), dimension(cgnsNSliding) :: interfaceIDs
       integer(kind=intType), dimension(cgnsNSliding) :: sortedColor

       integer(kind=intType), dimension(cgnsNSliding+1) :: mult

       logical, dimension(cgnsNSliding) :: interfaceNotDistributed

       type(tmpSlidingType), dimension(cgnsNSliding) :: distrSliding
!
!      Function definition.
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of sections, i.e. grid parts between
       ! sliding mesh interfaces, in the grid.

       call determineSections

       ! Determine the number of intersections between the sliding mesh
       ! interfaces and store for the highest number the ID of the
       ! lowest sliding mesh interface.

       call slidingIntersections

       ! Determine the processors that contribute in each sliding
       ! mesh interface.

       call procSliding(distrSliding)

       ! Determine a sorted version of the maximum number of faces,
       ! which occurs on a processor.

       do mm=1,cgnsNSliding
         nFaceMax(mm) = distrSliding(mm)%nFaceMax
       enddo

       ! Sort nFaceMax in increasing order.

       call qsortIntegers(nFaceMax,cgnsNSliding)

       ! Get rid of the multiple entries and store the multiplicity.

       nn   = min(1_intType,cgnsNSliding)
       mult = 1
       do mm=2,cgnsNSliding
         if(nFaceMax(mm) == nFaceMax(nn)) then
           mult(nn+1) = mult(nn+1) + 1
         else
           nn = nn + 1
           nFaceMax(nn) = nFaceMax(mm)
         endif
       enddo

       ! Put mult in cumulative storage format.

       mult(1) = 0
       do mm=1,nn
         mult(mm+1) = mult(mm+1) + mult(mm)
       enddo

       ! Determine the interface ID's belonging to the sorted
       ! version of nFaceMax.

       do mm=1,cgnsNSliding
         ii = bsearchIntegers(distrSliding(mm)%nFaceMax, &
                              nFaceMax, nn)
         mult(ii) = mult(ii) + 1
         interfaceIDs(mult(ii)) = mm
       enddo

       ! Initialize flag to 0. Flag is used as an indication whether or
       ! not a processor is already in use for the given color.
       ! Furthermore use mult as an indication whether or not a sliding
       ! mesh interface has already been distributed; 0 means that it
       ! has not been treated yet. The same info is stored in
       ! interfaceNotDistributed, but both are needed, because
       ! mult uses the sorted numbering and interfaceNotDistributed
       ! the original numbering.

       flag = 0
       mult = 0
       interfaceNotDistributed = .true.

       ! Determine the coloring of the sliding mesh interfaces.

       kk = 0
       nInterfaceGroups = 0
       nSliding = cgnsNSliding

       distribution: do

         ! Condition to exit the loop.

         if(nSliding == 0) exit

         ! Update the number of groups, i.e. colors.

         nInterfaceGroups = nInterfaceGroups + 1

         ! Loop to put as many sliding mesh interfaces for this color.

         color: do nn=1,cgnsNSliding

           ! Check if this interface has not been distributed yet.

           distributionCheck: if(mult(nn) == 0) then

             ! Store the interface id in mm and check if this interface
             ! can be distributed inside this color.

             ! First check if all the processors are still available.

             mm = interfaceIDs(nn)
             canBeStored = .true.

             do ii=1,distrSliding(mm)%nProcs1
               jj = distrSliding(mm)%procs1(ii)
               if(flag(jj) == nInterfaceGroups) canBeStored = .false.
             enddo

             do ii=1,distrSliding(mm)%nProcs2
               jj = distrSliding(mm)%procs2(ii)
               if(flag(jj) == nInterfaceGroups) canBeStored = .false.
             enddo

             ! Check if this interface intersects with lower numbered
             ! sliding interfaces. If so, these interfaces must be
             ! treated first. That is the convention used in this code.

             do ii=(nIntersecting(mm-1)+1),nIntersecting(mm)
               jj = intersecting(ii)
               if( interfaceNotDistributed(jj) ) canBeStored = .false.
             enddo

             ! If this interface can be stored for this color, do so.

             storeCheck: if( canBeStored ) then

               ! Update the counter kk, store the color in sortedColor
               ! and the index kk in mult(nn). Here a nonzero value of
               ! mult(nn) indicates that this interface has been treated.
               ! Later this value, or better its inverse, is used to
               ! store the coloring in myInterfaces

               kk = kk + 1
               sortedColor(kk) = nInterfaceGroups
               mult(nn) = kk

               ! Set interfaceNotDistributed(mm) to .false. to
               ! indicate that the interface has been distributed.

               interfaceNotDistributed(mm) = .false.

               ! Flag the processors of this interface to the current
               ! color to indicate that they are in use.

               do ii=1,distrSliding(mm)%nProcs1
                 jj = distrSliding(mm)%procs1(ii)
                 flag(jj) = nInterfaceGroups
               enddo

               do ii=1,distrSliding(mm)%nProcs2
                 jj = distrSliding(mm)%procs2(ii)
                 flag(jj) = nInterfaceGroups
               enddo

               ! Decrease nSliding.

               nSliding = nSliding - 1

             endif storeCheck
           endif distributionCheck

         enddo color
       enddo distribution

       ! Release the memory of nIntersecting and intersecting.

       deallocate(nIntersecting, intersecting, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("slidingMeshGroups", &
                        "Deallocation error for nIntersecting &
                        &and intersecting")

       ! In the determination of myInterfaces the inverse is needed
       ! of what is currently stored in mult. Use nFaceMax as a
       ! temporary buffer.

       do nn=1,cgnsNSliding
         nFaceMax(mult(nn)) = nn
       enddo

       do nn=1,cgnsNSliding
         mult(nn) = nFaceMax(nn)
       enddo

       ! Allocate the memory for myInterfaces.

       allocate(myInterfaces(nInterfaceGroups), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("slidingMeshGroups", &
                        "Memory allocation failure for myInterfaces")

       ! Determine for each color whether or not this processor
       ! participates in a search.

       nn = 1
       nGroups: do mm=1,nInterfaceGroups

         ! Find the starting index for this color. Note that sortedColor
         ! has been constructed for this purpose in the distribution loop.

         do
           if(sortedColor(nn) == mm) exit
           nn = nn + 1
         enddo

         ! Initialize procContributes to .false.

         myInterfaces(mm)%procContributes = .false.

         ! Loop to find out whether or not this processor must do a
         ! search for this color.

         do
           ! Conditions to exit this loop.

           if(myInterfaces(mm)%procContributes) exit
           if(nn > cgnsNSliding)                exit
           if(sortedColor(nn) > mm)             exit

           ! Store the ID of the sliding interface in ii, for which
           ! the mapping stored in mult is used.

           ii = interfaceIDs(mult(nn))

           ! Check if this processor participates in the search
           ! for this interface.

           do jj=1,distrSliding(ii)%nProcs1
             if(distrSliding(ii)%procs1(jj) == myID) &
               myInterfaces(mm)%procContributes = .true.
           enddo

           do jj=1,distrSliding(ii)%nProcs2
             if(distrSliding(ii)%procs2(jj) == myID) &
               myInterfaces(mm)%procContributes = .true.
           enddo

           ! Store the info in myInterfaces if this processor
           ! participates for this sliding mesh interface.

           if( myInterfaces(mm)%procContributes ) then

             ! Copy the scalar data.

             myInterfaces(mm)%globalSlideID = ii
             myInterfaces(mm)%nProcs1 = distrSliding(ii)%nProcs1
             myInterfaces(mm)%nProcs2 = distrSliding(ii)%nProcs2

             ! Allocate the memory for the processor id's.

             allocate(myInterfaces(mm)%procs1(myInterfaces(mm)%nProcs1), &
                      myInterfaces(mm)%procs2(myInterfaces(mm)%nProcs2), &
                      stat=ierr)
             if(ierr /= 0)                         &
               call terminate("slidingMeshGroups", &
                              "Memory allocation failure for procIDs")

             ! Copy the processor numbers.

             do jj=1,myInterfaces(mm)%nProcs1
               myInterfaces(mm)%procs1(jj) = distrSliding(ii)%procs1(jj)
             enddo

             do jj=1,myInterfaces(mm)%nProcs2
               myInterfaces(mm)%procs2(jj) = distrSliding(ii)%procs2(jj)
             enddo

           endif

           ! Update nn.

           nn = nn + 1

         enddo

       enddo nGroups

       ! Release the memory of distrSliding.

       do mm=1,cgnsNSliding
         deallocate(distrSliding(mm)%procs1, &
                    distrSliding(mm)%procs2, stat=ierr)
         if(ierr /= 0)                         &
           call terminate("slidingMeshGroups", &
                          "Deallocation failure for procs1 and procs2")
       enddo

       ! Loop over the number of interface groups and perform some
       ! initialization.

       do mm=1,nInterfaceGroups

         ! Do the initializations for this color.

         call initThisSlide(mm)

         ! Determine the rotational periodicity, i.e. the number of
         ! grid slices in a complete rotation, for both sides of
         ! this interface.

         myInterfaces(mm)%nSlices1 = 1
         myInterfaces(mm)%nSlices2 = 1

         if( myInterfaces(mm)%procContributes ) then

           call determineNSlices( myInterfaces(mm)%nSlices1,      &
                                 -myInterfaces(mm)%globalSlideID, &
                                  myInterfaces(mm)%commSlide)
 
           call determineNSlices(myInterfaces(mm)%nSlices2,      &
                                 myInterfaces(mm)%globalSlideID, &
                                 myInterfaces(mm)%commSlide)
         endif

         ! Block until all processors have reached this point.

         call mpi_barrier(SUmb_comm_world, ierr)

       enddo

       ! Determine the corresponding rotation matrices.

       call determineRotationMatrices

       end subroutine slidingMeshGroups

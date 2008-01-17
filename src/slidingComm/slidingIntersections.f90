!
!      ******************************************************************
!      *                                                                *
!      * File:          slidingIntersections.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-22-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module intersect
!
!      ******************************************************************
!      *                                                                *
!      * Local module to store the for every sliding mesh interface its *
!      * intersections.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! Definition of the derived data type to store the info.

       type intersectType

         ! nIntersections:  Number of intersections with lower numbered
         !                  sliding mesh interfaces.
         ! intersections(): The corresponding id's.

         integer(kind=intType) :: nIntersections

         integer(kind=intType), dimension(:), pointer :: intersections

       end type intersectType

       end module intersect

!      ==================================================================

       subroutine slidingIntersections
!
!      ******************************************************************
!      *                                                                *
!      * slidingIntersections determines for every sliding mesh         *
!      * interface the number of intersections as well as the sliding   *
!      * mesh ID's with lower numbered sliding mesh interfaces.         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use intersect
       use tmpSliding
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, nn, mm, il, jl, kl

       type(intersectType), dimension(cgnsNSliding) :: slideInt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the number of intersections for the sliding mesh
       ! interfaces to 0 and do an initial allocation of the
       ! intersections.

       do nn=1,cgnsNSliding
         slideInt(nn)%nIntersections = 0
         allocate(slideInt(nn)%intersections(0), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("slidingIntersections", &
                          "Memory allocation failure for intersections")
       enddo

       ! Loop over the number of blocks of the entire original grid.

       do nn=1,cgnsNDom

         ! Store the number of nodes a bit easier.

         il = cgnsDoms(nn)%il
         jl = cgnsDoms(nn)%jl
         kl = cgnsDoms(nn)%kl

         ! Update the intersections using the info of this block.

         call intersectionThisBlock

       enddo

       ! Determine the total number of intersections; store it in ii.

       ii = 0
       do nn=1,cgnsNSliding
         ii = ii + slideInt(nn)%nIntersections
       enddo

       ! Allocate the memory to store the final information.

       allocate(nIntersecting(0:cgnsNSliding), intersecting(ii), &
                stat=ierr)
       if(ierr /= 0)                            &
         call terminate("slidingIntersections", &
                        "Memory allocation failure for nIntersecting &
                        &and intersecting.")

       ! Copy the info from slideInt.

       ii = 0
       nIntersecting(0) = 0

       do nn=1,cgnsNSliding
         do mm=1,slideInt(nn)%nIntersections
           ii = ii + 1
           intersecting(ii) = slideInt(nn)%intersections(mm)
         enddo
         nIntersecting(nn) = ii
       enddo

       ! Release the memory of intersections inside the derived datatype.

       do nn=1,cgnsNSliding
         deallocate(slideInt(nn)%intersections, stat=ierr)
         if(ierr /= 0) &
           call terminate("slidingIntersections", &
                          "Deallocation error for intersections.")
       enddo

       !=================================================================

       contains

         !===============================================================

         subroutine intersectionThisBlock
!
!        ****************************************************************
!        *                                                              *
!        * IntersectionThisBlock determines the intersections for       *
!        * block nn in the original cgns grid and stores this info in   *
!        * the array slideInt.                                          *
!        *                                                              *
!        ****************************************************************
!
         use su_cgns
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: i, j, k
         integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
         integer(kind=intType) :: slideID

         integer(kind=intType), dimension(il) :: j1k1, jlk1, j1kl, jlkl
         integer(kind=intType), dimension(jl) :: i1k1, ilk1, i1kl, ilkl
         integer(kind=intType), dimension(kl) :: i1j1, ilj1, i1jl, iljl
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize the edge indicators of the block to 0.

         do i=1,il
           j1k1(i) = 0; jlk1(i) = 0; j1kl(i) = 0; jlkl(i) = 0
         enddo

         do j=1,jl
           i1k1(j) = 0; ilk1(j) = 0; i1kl(j) = 0; ilkl(j) = 0
         enddo

         do k=1,kl
           i1j1(k) = 0; ilj1(k) = 0; i1jl(k) = 0; iljl(k) = 0
         enddo

         ! Loop over the number of boundary subfaces of this block.

         do mm=1,cgnsDoms(nn)%nBocos

           ! Test if we are having a sliding mesh subface here.

           if(cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
              cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             ! Subface on a sliding mesh interface. Determine the global
             ! id of the interface via the boco info. Note that the
             ! absolute value must be taken, because the slidingID can
             ! be both positive and negative, depending on the side.

             slideID = abs(cgnsDoms(nn)%bocoInfo(mm)%slidingID)

             ! Copy the nodal range in iBeg, iEnd, etc. Make sure that
             ! iBeg contains the lowest value and iEnd the highest.

             iBeg = min(cgnsDoms(nn)%bocoInfo(mm)%iBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%iEnd)
             jBeg = min(cgnsDoms(nn)%bocoInfo(mm)%jBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%jEnd)
             kBeg = min(cgnsDoms(nn)%bocoInfo(mm)%kBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%kEnd)

             iEnd = max(cgnsDoms(nn)%bocoInfo(mm)%iBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%iEnd)
             jEnd = max(cgnsDoms(nn)%bocoInfo(mm)%jBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%jEnd)
             kEnd = max(cgnsDoms(nn)%bocoInfo(mm)%kBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%kEnd)

             ! Determine the block face on which the subface is located.

             if(iBeg == iEnd) then

               ! Either iMin or iMax. Determine this and set the value
               ! of the 4 edge indicators to the maximum of the current
               ! value and slideID. In this case always the highest id
               ! is taken in case of intersecting interfaces.

               if(iBeg == 1) then
                 do j=jBeg,jEnd
                   i1k1(j) = max(slideID, i1k1(j))
                   i1kl(j) = max(slideID, i1kl(j))
                 enddo

                 do k=kBeg,kEnd
                   i1j1(k) = max(slideID, i1j1(k))
                   i1jl(k) = max(slideID, i1jl(k))
                 enddo
               else
                 do j=jBeg,jEnd
                   ilk1(j) = max(slideID, ilk1(j))
                   ilkl(j) = max(slideID, ilkl(j))
                 enddo

                 do k=kBeg,kEnd
                   ilj1(k) = max(slideID, ilj1(k))
                   iljl(k) = max(slideID, iljl(k))
                 enddo
               endif

             else if(jBeg == jEnd) then

               ! Either jMin or jMax. Same story as for iMin/iMax.

               if(jBeg == 1) then
                 do i=iBeg,iEnd
                   j1k1(i) = max(slideID, j1k1(i))
                   j1kl(i) = max(slideID, j1kl(i))
                 enddo

                 do k=kBeg,kEnd
                   i1j1(k) = max(slideID, i1j1(k))
                   ilj1(k) = max(slideID, ilj1(k))
                 enddo
               else
                 do i=iBeg,iEnd
                   jlk1(i) = max(slideID, jlk1(i))
                   jlkl(i) = max(slideID, jlkl(i))
                 enddo

                 do k=kBeg,kEnd
                   i1jl(k) = max(slideID, i1jl(k))
                   iljl(k) = max(slideID, iljl(k))
                 enddo
               endif

             else   ! no need to test; this has been done before.

               ! Either kMin or kMax. Same story as for iMin/iMax.

               if(kBeg == 1) then
                 do i=iBeg,iEnd
                   j1k1(i) = max(slideID, j1k1(i))
                   jlk1(i) = max(slideID, jlk1(i))
                 enddo

                 do j=jBeg,jEnd
                   i1k1(j) = max(slideID, i1k1(j))
                   ilk1(j) = max(slideID, ilk1(j))
                 enddo
               else
                 do i=iBeg,iEnd
                   j1kl(i) = max(slideID, j1kl(i))
                   jlkl(i) = max(slideID, jlkl(i))
                 enddo

                 do j=jBeg,jEnd
                   i1kl(j) = max(slideID, i1kl(j))
                   ilkl(j) = max(slideID, ilkl(j))
                 enddo
               endif

             endif

           endif
         enddo

         ! Repeat the loop over the sliding mesh subfaces, but now
         ! store the intersection info.

         do mm=1,cgnsDoms(nn)%nBocos

           ! Test if we are having a sliding mesh subface here.

           if(cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
              cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             ! Subface on a sliding mesh interface. Determine the global
             ! ID of the interface.

             slideID = abs(cgnsDoms(nn)%bocoInfo(mm)%slidingID)

             ! Copy the nodal range in iBeg, iEnd, etc. Make sure that
             ! iBeg contains the lowest value and iEnd the highest.

             iBeg = min(cgnsDoms(nn)%bocoInfo(mm)%iBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%iEnd)
             jBeg = min(cgnsDoms(nn)%bocoInfo(mm)%jBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%jEnd)
             kBeg = min(cgnsDoms(nn)%bocoInfo(mm)%kBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%kEnd)

             iEnd = max(cgnsDoms(nn)%bocoInfo(mm)%iBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%iEnd)
             jEnd = max(cgnsDoms(nn)%bocoInfo(mm)%jBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%jEnd)
             kEnd = max(cgnsDoms(nn)%bocoInfo(mm)%kBeg, &
                        cgnsDoms(nn)%bocoInfo(mm)%kEnd)

             ! Determine the block face on which the subface is located.

             if(iBeg == iEnd) then

               ! Either iMin or iMax. Check if the edges are flagged
               ! with a different number. If so, there is an intersection
               ! where slideID is the smallest.

               if(iBeg == 1) then
                 do j=jBeg,jEnd
                   if(i1k1(j) /= slideID) &
                         call updateInt(slideInt, i1k1(j), slideID)
                   if(i1kl(j) /= slideID) &
                         call updateInt(slideInt, i1kl(j), slideID)
                 enddo

                 do k=kBeg,kEnd
                   if(i1j1(k) /= slideID) &
                         call updateInt(slideInt, i1j1(k), slideID)
                   if(i1jl(k) /= slideID) &
                         call updateInt(slideInt, i1jl(k), slideID)
                 enddo
               else
                 do j=jBeg,jEnd
                   if(ilk1(j) /= slideID) &
                         call updateInt(slideInt, ilk1(j), slideID)
                   if(ilkl(j) /= slideID) &
                         call updateInt(slideInt, ilkl(j), slideID)
                 enddo

                 do k=kBeg,kEnd
                   if(ilj1(k) /= slideID) &
                         call updateInt(slideInt, ilj1(k), slideID)
                   if(iljl(k) /= slideID) &
                         call updateInt(slideInt, iljl(k), slideID)
                 enddo
               endif

             else if(jBeg == jEnd) then

               ! Either jMin or jMax. Same story as for iMin/iMax.

               if(jBeg == 1) then
                 do i=iBeg,iEnd
                   if(j1k1(i) /= slideID) &
                         call updateInt(slideInt, j1k1(i), slideID)
                   if(j1kl(i) /= slideID) &
                         call updateInt(slideInt, j1kl(i), slideID)
                 enddo

                 do k=kBeg,kEnd
                   if(i1j1(k) /= slideID) &
                         call updateInt(slideInt, i1j1(k), slideID)
                   if(ilj1(k) /= slideID) &
                         call updateInt(slideInt, ilj1(k), slideID)
                 enddo
               else
                 do i=iBeg,iEnd
                   if(jlk1(i) /= slideID) &
                         call updateInt(slideInt, jlk1(i), slideID)
                   if(jlkl(i) /= slideID) &
                         call updateInt(slideInt, jlkl(i), slideID)
                 enddo

                 do k=kBeg,kEnd
                   if(i1jl(k) /= slideID) &
                         call updateInt(slideInt, i1jl(k), slideID)
                   if(iljl(k) /= slideID) &
                         call updateInt(slideInt, iljl(k), slideID)
                 enddo
               endif

             else   ! No need to test; this has been done before.

               ! Either kMin or kMax. Same story as for iMin/iMax.

               if(kBeg == 1) then
                 do i=iBeg,iEnd
                   if(j1k1(i) /= slideID) &
                         call updateInt(slideInt, j1k1(i), slideID)
                   if(jlk1(i) /= slideID) &
                         call updateInt(slideInt, jlk1(i), slideID)
                 enddo

                 do j=jBeg,jEnd
                   if(i1k1(j) /= slideID) &
                         call updateInt(slideInt, i1k1(j), slideID)
                   if(ilk1(j) /= slideID) &
                         call updateInt(slideInt, ilk1(j), slideID)
                 enddo
               else
                 do i=iBeg,iEnd
                   if(j1kl(i) /= slideID) &
                         call updateInt(slideInt, j1kl(i), slideID)
                   if(jlkl(i) /= slideID) &
                         call updateInt(slideInt, jlkl(i), slideID)
                 enddo

                 do j=kBeg,kEnd
                   if(i1kl(j) /= slideID) &
                         call updateInt(slideInt, i1kl(j), slideID)
                   if(ilkl(j) /= slideID) &
                         call updateInt(slideInt, ilkl(j), slideID)
                 enddo
               endif

             endif

           endif
         enddo

         end subroutine intersectionThisBlock

       end subroutine slidingIntersections

!      ==================================================================

       subroutine updateInt(slideInt, high, low)
!
!      ******************************************************************
!      *                                                                *
!      * updateInt checks whether or not the intersection between the   *
!      * sliding mesh interfaces high and low is already stored. If not *
!      * it is stored in slideInt.                                      *
!      *                                                                *
!      ******************************************************************
!
       use intersect
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: high, low
       type(intersectType), dimension(*), intent(inout) :: slideInt
!
!      Local variables.
!
       integer(kind=intType) :: l, lOld
!
!      Interfaces
!
       interface
         subroutine reallocateInteger(intArray, newSize, oldSize, &
                                      alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize, oldSize
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if the intersection is already stored.

       do l=1,slideInt(high)%nIntersections
         if(slideInt(high)%intersections(l) == low) exit
       enddo

       ! If the intersection has not been stored yet, store it now.

       if(l > slideInt(high)%nIntersections) then

         lOld = slideInt(high)%nIntersections
         l     = lOld + 1
         slideInt(high)%nIntersections = l

         call reallocateInteger(slideInt(high)%intersections, l, &
                                lOld, .true.)

         slideInt(high)%intersections(l) = low

       endif

       end subroutine updateInt

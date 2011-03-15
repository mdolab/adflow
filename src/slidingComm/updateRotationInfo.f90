!
!      ******************************************************************
!      *                                                                *
!      * File:          updateRotationInfo.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-10-2003                                      *
!      * Last modified: 04-08-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateRotationInfo(level, sps, color)
!
!      ******************************************************************
!      *                                                                *
!      * updateRotationInfo determines for every cell which is part     *
!      * of the sliding mesh interface (halo's included) the rotation   *
!      * matrix needed to obtain the correct cartesian velocities.      *
!      * The connectivity needed for this is also used for              *
!      * initialization in the communication routines whaloSliding.     *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       use interfaceGroups
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, color
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ii
       integer(kind=intType) :: n1st, n2nd, n1stOld, n2ndOld
       integer(kind=intType) :: offsetColor, nSlices

       integer(kind=intType), dimension(:), pointer :: block1, block2
       integer(kind=intType), dimension(:), pointer :: rotInd1, rotInd2
       integer(kind=intType), dimension(:,:), pointer :: ind1, ind2
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

         !===============================================================

         subroutine reallocateInteger2(intArray,           &
                                       newSize1, newSize2, &
                                       oldSize1, oldSize2, &
                                       alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:,:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger2
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the local number of first and second level halo cells
       ! for this sliding interface.

       n1st = 0
       n2nd = 0

       ! First the subfaces of part 1 of the interface.

       do nn=1,nMySubfaces1
         do mm=1,mySubfaces1(nn)%nQuad

           ! Update the number of 1st and 2nd halo's if a donor must be
           ! sought for this face. If this is a direct halo n2nd must be
           ! updated by 2; a direct halo is indicated by a value > -1
           ! for the first index of the second halo.

           if( mySubfaces1(nn)%searchQuad(mm) ) then

             n1st = n1st + 1
             n2nd = n2nd + 1
             if(mySubfaces1(nn)%indHalo2(mm,1) > -1) n2nd = n2nd + 1

           endif

         enddo
       enddo

       ! And for the subfaces of part 2.

       do nn=1,nMySubfaces2
         do mm=1,mySubfaces2(nn)%nQuad

           if( mySubfaces2(nn)%searchQuad(mm) ) then

             n1st = n1st + 1
             n2nd = n2nd + 1
             if(mySubfaces2(nn)%indHalo2(mm,1) > -1) n2nd = n2nd + 1

           endif

         enddo
       enddo

       ! Store the old value of nSlidingHalos for both the 1st and
       ! 2nd level halo's and update the values.

       n1stOld = intSlidingCell_1st(level,sps)%nSlidingHalos
       n2ndOld = intSlidingCell_2nd(level,sps)%nSlidingHalos

       n1st = n1st + n1stOld
       n2nd = n2nd + n2ndOld

       intSlidingCell_1st(level,sps)%nSlidingHalos = n1st
       intSlidingCell_2nd(level,sps)%nSlidingHalos = n2nd

       ! Reallocate the memory to store both the block indices and
       ! the rotation indices for both the 1st and 2nd level halo's.

       call reallocateInteger(                                       &
                intSlidingCell_1st(level,sps)%slidingHaloList%block, &
                n1st, n1stOld, .true.)
       call reallocateInteger2(                                        &
                intSlidingCell_1st(level,sps)%slidingHaloList%indices, &
                n1st, 3_intType, n1stOld, 3_intType, .true.)
       call reallocateInteger(intSlidingCell_1st(level,sps)%rotIndex, &
                              n1st, n1stOld, .true.)

       call reallocateInteger(                                       &
                intSlidingCell_2nd(level,sps)%slidingHaloList%block, &
                n2nd, n2ndOld, .true.)
       call reallocateInteger2(                                        &
                intSlidingCell_2nd(level,sps)%slidingHaloList%indices, &
                n2nd, 3_intType, n2ndOld, 3_intType, .true.)
       call reallocateInteger(intSlidingCell_2nd(level,sps)%rotIndex, &
                              n2nd, n2ndOld, .true.)

       ! Set a couple of pointers to make the code more readable.

       block1 => intSlidingCell_1st(level,sps)%slidingHaloList%block
       block2 => intSlidingCell_2nd(level,sps)%slidingHaloList%block

       ind1 => intSlidingCell_1st(level,sps)%slidingHaloList%indices
       ind2 => intSlidingCell_2nd(level,sps)%slidingHaloList%indices

       rotInd1 => intSlidingCell_1st(level,sps)%rotIndex
       rotInd2 => intSlidingCell_2nd(level,sps)%rotIndex

       ! Determine the offset in the rotation matrices for this color.

       offsetColor = 0
       do ii=1,(color-1)
         offsetColor = offsetColor + myInterfaces(ii)%nSlices1 - 1
       enddo

       ! Repeat the loops over the subfaces of the sliding interface,
       ! but now store the information.

       n1st    = n1stOld
       n2nd    = n2ndOld
       nSlices = myInterfaces(color)%nSlices1

       call updateRotViaSubfaces(nMySubfaces1, mySubfaces1)
       call updateRotViaSubfaces(nMySubfaces2, mySubfaces2)

       !=================================================================

       contains

         subroutine updateRotViaSubfaces(nMySubfaces, mySubfaces)
!
!        ****************************************************************
!        *                                                              *
!        * updateRotViaSubfaces updates the rotation info of            *
!        * intSlidingCell_1st and intSlidingCell_2nd using the data     *
!        * stored in mySubfaces. Some pointers have already been set    *
!        * to the members of intSlidingCell_1st and                     *
!        * intSlidingCell_2nd to make it more readable.                 *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nMySubfaces
         type(localSubfaceType), dimension(*), intent(inout) :: mySubfaces
!
!        Local variables.
!
         integer :: ierr
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         subfaceLoop: do nn=1,nMySubfaces
           quadLoop: do mm=1,mySubfaces(nn)%nQuad

             ! Check if the halo has been interpolated.

             if( mySubfaces(nn)%searchQuad(mm) ) then

               ! Determine the index of the corresponding rotation
               ! matrix for the velocities. Note that the number of
               ! rotations must be negated, as in nRotations the number
               ! of rotations from mySubfaces to the donor is stored and
               ! the inverse information is needed. Furthermore remember
               ! that nRotSliding = nSlices - 1, as the unity
               ! transformation is not stored.

               ii = -mySubfaces(nn)%nRotations(mm)
               if(ii < 0) ii = ii + nSlices
               if(ii > 0) ii = ii + offsetColor

               ! Update the information for the 1st level halo cells.

               n1st = n1st + 1
               block1(n1st)  = mySubfaces(nn)%blockId
               ind1(n1st,1)  = mySubfaces(nn)%indHalo1(mm,1)
               ind1(n1st,2)  = mySubfaces(nn)%indHalo1(mm,2)
               ind1(n1st,3)  = mySubfaces(nn)%indHalo1(mm,3)
               rotInd1(n1st) = ii

               ! Store the 1st halo cell also in the 2nd level halo
               ! communication pattern.

               n2nd = n2nd + 1
               block2(n2nd)  = mySubfaces(nn)%blockId
               ind2(n2nd,1)  = mySubfaces(nn)%indHalo1(mm,1)
               ind2(n2nd,2)  = mySubfaces(nn)%indHalo1(mm,2)
               ind2(n2nd,3)  = mySubfaces(nn)%indHalo1(mm,3)
               rotInd2(n2nd) = ii

               ! Only for direct halo's the 2nd level halo cell is
               ! needed, which happens here.

               if(mySubfaces(nn)%indHalo2(mm,1) > -1) then
                 n2nd = n2nd + 1
                 block2(n2nd)  = mySubfaces(nn)%blockId
                 ind2(n2nd,1)  = mySubfaces(nn)%indHalo2(mm,1)
                 ind2(n2nd,2)  = mySubfaces(nn)%indHalo2(mm,2)
                 ind2(n2nd,3)  = mySubfaces(nn)%indHalo2(mm,3)
                 rotInd2(n2nd) = ii
               endif

             endif
           enddo quadLoop

           ! Release the memory of nRotations inside mySubfaces.

           deallocate(mySubfaces(nn)%nRotations, stat=ierr)
           if(ierr /= 0)                          &
             call terminate("updateRotationInfo", &
                            "Deallocation error for nRotations of &
                            &mySubfaces.")
         enddo subfaceLoop

         end subroutine updateRotViaSubfaces

       end subroutine updateRotationInfo

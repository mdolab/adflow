!
!      ******************************************************************
!      *                                                                *
!      * File:          haloType.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-17-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine haloType(level)
!
!      ******************************************************************
!      *                                                                *
!      * haloType determines for the internal and 1st level halo cells  *
!      * what type of halo's they are. This information is needed to    *
!      * determine which halo cells/nodes should be interpolated and    *
!      * which not.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use communication
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: donorType

       integer(kind=intType) :: i, j, ii, jj, nn
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       integer(kind=intType), dimension(:,:), pointer :: haloInfo
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for donorDoms.

       allocate(donorDoms(nDom), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("statusDonorCells", &
                        "Memory allocation failure for donorDoms")

       ! Loop over the domains to allocate and initialize haloInfo.

       do ii=1,nDom

         ie = flowDoms(ii,level,1)%ie
         je = flowDoms(ii,level,1)%je
         ke = flowDoms(ii,level,1)%ke

         allocate(donorDoms(ii)%haloInfo(ie,je,ke), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("statusDonorCells", &
                          "Memory allocation failure for haloInfo")

         donorDoms(ii)%haloInfo = internalCell

       enddo

       ! Loop over the halo's of the 1st layer communication pattern to
       ! set haloInfo to either an owned or an unowned donor.
       ! First the internal communication pattern.

       localHalos: do ii=1,internalCell_1st(level)%ncopy

         ! Store the block and the indices of the donor a bit easier.

         d1 = internalCell_1st(level)%donorBlock(ii)
         i1 = internalCell_1st(level)%donorIndices(ii,1)
         j1 = internalCell_1st(level)%donorIndices(ii,2)
         k1 = internalCell_1st(level)%donorIndices(ii,3)

         ! Idem for the halo's.

         d2 = internalCell_1st(level)%haloBlock(ii)
         i2 = internalCell_1st(level)%haloIndices(ii,1)
         j2 = internalCell_1st(level)%haloIndices(ii,2)
         k2 = internalCell_1st(level)%haloIndices(ii,3)

         ! Set haloInfo to owned donor if the donor point is
         ! considered smaller than the halo point. Otherwise set
         ! it to an unowned donor.

         if(d1 < d2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = ownedDonor
           cycle
         else if(d1 > d2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = unownedDonor
           cycle
         endif

         if(k1 < k2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = ownedDonor
           cycle
         else if(k1 > k2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = unownedDonor
           cycle
         endif

         if(j1 < j2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = ownedDonor
           cycle
         else if(j1 > j2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = unownedDonor
           cycle
         endif

         if(i1 < i2) then
           donorDoms(d2)%haloInfo(i2,j2,k2) = ownedDonor
           cycle
         else
           donorDoms(d2)%haloInfo(i2,j2,k2) = unownedDonor
           cycle
         endif

       enddo localHalos

       ! The external communication pattern. As the donor is stored on
       ! a different processor by definition, only the processor id's
       ! need to be compared.

       recvProcLoop: do ii=1,commPatternCell_1st(level)%nProcRecv

         ! Set haloInfo to owned if the processor from which i
         ! receive the data has a smaller id than myID. Otherwise
         ! set haloInfo to unowned

         if(commPatternCell_1st(level)%recvProc(ii) < myID) then
           donorType = ownedDonor
         else
           donorType = unownedDonor
         endif

         ! Loop over the halo received from this processor.

         do jj=1,commPatternCell_1st(level)%nrecv(ii)

           ! Store the halo indices a bit easier.

           d2 = commPatternCell_1st(level)%recvList(ii)%block(jj)
           i2 = commPatternCell_1st(level)%recvList(ii)%indices(jj,1)
           j2 = commPatternCell_1st(level)%recvList(ii)%indices(jj,2)
           k2 = commPatternCell_1st(level)%recvList(ii)%indices(jj,3)

           ! Set haloInfo.

           donorDoms(d2)%haloInfo(i2,j2,k2) = donorType

         enddo
       enddo recvProcLoop

       ! Loop over the domains to set the sliding mesh and
       ! boundary halo's.

       domains: do ii=1,nDom

         ! Set the pointers for this block to make the code more
         ! readable.

         call setPointers(ii,level,1_intType)

         ! Loop over the sliding mesh interfaces and set the halo info
         ! to the highest sliding mesh interface id.

         sliding: do jj=1,nBocos
           if(BCType(jj) == slidingInterface) then

             ! Set the pointer for haloInfo for this subface and set the
             ! first level halo cells to the maximum of the id of the
             ! sliding interface and the currently stored value. The max
             ! function is present to be able to distinguish between two
             ! intersecting sliding interfaces.

             call haloTypePointers
             nn = abs(groupNum(jj))

             do j=BCData(jj)%jcBeg, BCData(jj)%jcEnd
               do i=BCData(jj)%icBeg, BCData(jj)%icEnd
                 haloInfo(i,j) = max(haloInfo(i,j),nn)
               enddo
             enddo

           endif
         enddo sliding

         ! Loop over the boundary subfaces and set the halo info
         ! to boundaryHalo.

         bocos: do jj=1,nBocos
           if(BCType(jj) /= slidingInterface) then

             ! Set the pointer for haloInfo for this subface and set
             ! the first level halo cells to boundaryHalo.

             call haloTypePointers

             do j=BCData(jj)%jcBeg, BCData(jj)%jcEnd
               do i=BCData(jj)%icBeg, BCData(jj)%icEnd
                 haloInfo(i,j) = boundaryHalo
               enddo
             enddo

           endif
         enddo bocos

       enddo domains

       !=================================================================

       contains

         !===============================================================

         subroutine haloTypePointers
!
!        ****************************************************************
!        *                                                              *
!        * HaloTypePointers set the pointer for haloInfo depending      *
!        * on the block face on which the subface is located.           *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         select case (BCFaceID(jj))

           case (iMin)
             haloInfo => donorDoms(ii)%haloInfo(1,:,:)

           case (iMax)
             haloInfo => donorDoms(ii)%haloInfo(ie,:,:)

           case (jMin)
             haloInfo => donorDoms(ii)%haloInfo(:,1,:)

           case (jMax)
             haloInfo => donorDoms(ii)%haloInfo(:,je,:)

           case (kMin)
             haloInfo => donorDoms(ii)%haloInfo(:,:,1)

           case (kMax)
             haloInfo => donorDoms(ii)%haloInfo(:,:,ke)

         end select

         end subroutine haloTypePointers

       end subroutine haloType

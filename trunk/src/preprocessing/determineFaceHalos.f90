!
!      ******************************************************************
!      *                                                                *
!      * File:          determineFaceHalos.F90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-28-2003                                      *
!      * Last modified: 11-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineFaceHalos(level)
!
!      ******************************************************************
!      *                                                                *
!      * determineFaceHalos determines the 1st level direct cell and    *
!      * node halo's. Direct halo means that at least one of the        *
!      * neighboring cell/nodes belongs is owned by the block.          *
!      * Consequently the halo can be found using the 1 to 1 block      *
!      * connectivity.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use haloList
       use periodicInfo
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: ii, jj, kk, ll, mm, nn
       integer(kind=intType) :: iH, jH, kH, iD, jD, kD
       integer(kind=intType) :: indexPeriodic

       integer(kind=intType), dimension(3) :: myOffset, donorOffset
       integer(kind=intType), dimension(3) :: step

       integer(kind=intType), dimension(3,3) :: tMat
       integer(kind=intType), dimension(3,2) :: myCellRange
       integer(kind=intType), dimension(3,2) :: myNodeRange
       integer(kind=intType), dimension(3,2) :: donorCellRange

       type(cgnsPeriodicType) :: key
!
!      Function definitions.
!
       integer(kind=intType) :: delta
       integer(kind=intType) :: bsearchCGNSPeriodicType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the counter variables for the 1st level halo's to 0.

       iicell1st = 0
       iinode1st = 0

       ! Determine the 1st level cell and node halo lists by looping over
       ! the boundary faces of the blocks stored on this processor.

       domains: do nn=1,nDom

         ! Set the pointers for this block on the given level.

         call setPointers(nn,level,1_intType)
!
!        ****************************************************************
!        *                                                              *
!        * Loop over the boundary halo's first. The reason is that a    *
!        * node could belong to both a boundary and an internal subface.*
!        * By looping first over the boundaries, the internal subfaces  *
!        * overwrite earlier set values by the boundary subface, which  *
!        * is desirable.                                                *
!        *                                                              *
!        ****************************************************************
!
         bocos: do mm=1,nBocos

           ! Determine the cell and nodal range for the halo's of this
           ! subface as well as the direction normal to the subface.

           call haloRanges(mm, myOffset, myCellRange, myNodeRange, step)
!
!          **************************************************************
!          *                                                            *
!          * First treat the nodes on the subface.                      *
!          *                                                            *
!          **************************************************************
!
           ! Loop over the nodes of the boundary subface and store
           ! the halo info. For the edges of the subface it is possible
           ! that the node is already stored. This must be checked;
           ! otherwise this node will occur more than once in the list,
           ! which is not correct.

           do k=myNodeRange(3,1),myNodeRange(3,2),step(3)
             do j=myNodeRange(2,1),myNodeRange(2,2),step(2)
               do i=myNodeRange(1,1),myNodeRange(1,2),step(1)

                 ! Determine the indices of my nodal halo node.

                 iH = i + myOffset(1)
                 jH = j + myOffset(2)
                 kH = k + myOffset(3)

                 ! Store the halo, if it has not been stored yet.

                 if(nodeIndex(nn)%entryList(iH,jH,kH) == 0) then

                   ! Update the counter iinode1st. Store it in ii for
                   ! convenience. Store the index in %entryList.

                   iinode1st = iinode1st +1
                   ii        = iinode1st

                   nodeIndex(nn)%entryList(iH,jH,kH) = ii

                   ! Store the info in nodeHalo1st. Note that
                   ! donorBlock contains the boundary condition, although
                   ! this is not really important for node halo's.
                   ! Furthermore donorProc is set to -1 to indicate a
                   ! boundary halo and the donor indices are set to i,j,k.
                   ! In this way the extrapolated coordinates can be
                   ! obtained easily later on.

                   nodeHalo1st(ii)%myBlock = nn
                   nodeHalo1st(ii)%myI     = iH
                   nodeHalo1st(ii)%myJ     = jH
                   nodeHalo1st(ii)%myK     = kH

                   nodeHalo1st(ii)%donorProc  = -1
                   nodeHalo1st(ii)%donorBlock = BCType(mm)

                   nodeHalo1st(ii)%dI = i
                   nodeHalo1st(ii)%dJ = j
                   nodeHalo1st(ii)%dK = k

                 endif

               enddo
             enddo
           enddo
!
!          **************************************************************
!          *                                                            *
!          * The cell halo's belonging to this subface. Direct cell     *
!          * halo's are unique and therefore info cannot already be     *
!          * written earlier.                                           *
!          *                                                            *
!          **************************************************************
!
           ! Loop over the halo cells located adjacent to the subface.

           do k=myCellRange(3,1),myCellRange(3,2),step(3)
             do j=myCellRange(2,1),myCellRange(2,2),step(2)
               do i=myCellRange(1,1),myCellRange(1,2),step(1)

                 ! Check in debug mode whether this halo is already
                 ! stored. This should not be the case.

                 if( debug ) then
                   if(cellIndex(nn)%entryList(i,j,k) /= 0) &
                     call terminate("determineFaceHalos",  &
                                    "boundary cell halo already stored")
                 endif

                 ! Update the counter iicell1st and store its value a
                 ! bit easier in ii and set entryList accordingly.

                 iicell1st = iicell1st +1
                 ii        = iicell1st

                 cellIndex(nn)%entryList(i,j,k) = ii

                 ! Store the info in cellHalo1st. Note that donorBlock
                 ! contains the boundary condition and donorProc is set
                 ! to -1 to indicate a boundary halo. The donor indices are
                 ! set to the owned cell on the other side of the subface.

                 cellHalo1st(ii)%myBlock = nn
                 cellHalo1st(ii)%myI     = i
                 cellHalo1st(ii)%myJ     = j
                 cellHalo1st(ii)%myK     = k

                 cellHalo1st(ii)%donorProc  = -1
                 cellHalo1st(ii)%donorBlock = BCType(mm)

                 cellHalo1st(ii)%dI = i - myOffset(1)
                 cellHalo1st(ii)%dJ = j - myOffset(2)
                 cellHalo1st(ii)%dK = k - myOffset(3)

               enddo
             enddo
           enddo

         enddo bocos
!
!        ****************************************************************
!        *                                                              *
!        * Loop over the 1 to 1 block to block boundaries.              *
!        *                                                              *
!        ****************************************************************
!
         n1to1Loop: do ll=1,n1to1

           ! Store the correct index for this subface, i.e. add the
           ! offset from the boundary subfaces.

           mm = nBocos + ll

           ! Check if the original subface is a periodic subface.
           ! Subfaces created by internal block splitting are indicated
           ! by 0 and are certainly not periodic. This must be tested
           ! first to avoid array overflow.

           indexPeriodic = 0

           kk = cgnsSubface(mm)
           if(kk > 0) then
             if(cgnsDoms(nbkGlobal)%conn1to1(kk)%periodic) then

               ! Determine the corresponding index in periodicGlobal.

               key%cgnsBlock   = nbkGlobal
               key%cgnsSubface = kk

               indexPeriodic = bsearchCGNSPeriodicType(key,            &
                                                       periodicGlobal, &
                                                       nPeriodicGlobal)
               if( debug ) then
                 if(indexPeriodic == 0)                 &
                   call terminate("determineFaceHalos", &
                                  "Entry not found in periodicGlobal")
               endif
             endif
           endif

           ! Determine the cell and nodal range for the halo's of this
           ! subface as well as the direction normal to the subface.

           call haloRanges(mm, myOffset, myCellRange, &
                           myNodeRange, step)

           ! Determine the complete transformation matrix from the
           ! given shorthand.

           tMat(1,1) = sign(1_intType,l1(mm)) * delta(l1(mm),1_intType)
           tMat(2,1) = sign(1_intType,l1(mm)) * delta(l1(mm),2_intType)
           tMat(3,1) = sign(1_intType,l1(mm)) * delta(l1(mm),3_intType)

           tMat(1,2) = sign(1_intType,l2(mm)) * delta(l2(mm),1_intType)
           tMat(2,2) = sign(1_intType,l2(mm)) * delta(l2(mm),2_intType)
           tMat(3,2) = sign(1_intType,l2(mm)) * delta(l2(mm),3_intType)

           tMat(1,3) = sign(1_intType,l3(mm)) * delta(l3(mm),1_intType)
           tMat(2,3) = sign(1_intType,l3(mm)) * delta(l3(mm),2_intType)
           tMat(3,3) = sign(1_intType,l3(mm)) * delta(l3(mm),3_intType)

           ! Determine the offset of the donor block.

           donorOffset(1) = tMat(1,1)*myOffset(1) &
                          + tMat(1,2)*myOffset(2) &
                          + tMat(1,3)*myOffset(3)
           donorOffset(2) = tMat(2,1)*myOffset(1) &
                          + tMat(2,2)*myOffset(2) &
                          + tMat(2,3)*myOffset(3)
           donorOffset(3) = tMat(3,1)*myOffset(1) &
                          + tMat(3,2)*myOffset(2) &
                          + tMat(3,3)*myOffset(3)
!
!          **************************************************************
!          *                                                            *
!          * First treat the nodes on the subface.                      *
!          *                                                            *
!          **************************************************************
!
           ! Loop over the nodal range for this subface.

           do k=myNodeRange(3,1),myNodeRange(3,2),step(3)
             do j=myNodeRange(2,1),myNodeRange(2,2),step(2)
               do i=myNodeRange(1,1),myNodeRange(1,2),step(1)

                 ! Determine the donor indices by applying the
                 ! transformation matrix to i,j,k and adding the
                 ! offset to obtain the halo.

                 ii = i - myNodeRange(1,1)
                 jj = j - myNodeRange(2,1)
                 kk = k - myNodeRange(3,1)

                 iD = donorOffset(1) + dinBeg(mm) &
                    + tMat(1,1)*ii + tMat(1,2)*jj + tMat(1,3)*kk
                 jD = donorOffset(2) + djnBeg(mm) &
                    + tMat(2,1)*ii + tMat(2,2)*jj + tMat(2,3)*kk
                 kD = donorOffset(3) + dknBeg(mm) &
                    + tMat(3,1)*ii + tMat(3,2)*jj + tMat(3,3)*kk

                 ! Determine the indices of my nodal halo node.

                 iH = i + myOffset(1)
                 jH = j + myOffset(2)
                 kH = k + myOffset(3)

                 ! It is possible that this halo is already stored,
                 ! either as a boundary or as an internal halo. In the
                 ! former case it should be overwritten; in the latter
                 ! this is not strictly necessary, but it does not hurt.
                 ! Therefore simply overwrite the old index. If the
                 ! halo has not been stored yet, update iinode1st.
                 ! The index to store the info will be ii.

                 if(nodeIndex(nn)%entryList(iH,jH,kH) == 0) then
                   iinode1st = iinode1st +1
                   ii        = iinode1st

                   nodeIndex(nn)%entryList(iH,jH,kH) = ii
                 else
                   ii = nodeIndex(nn)%entryList(iH,jH,kH)
                 endif

                 ! Store the info in the correct place in nodeHalo1st.

                 nodeHalo1st(ii)%myBlock = nn
                 nodeHalo1st(ii)%myI     = iH
                 nodeHalo1st(ii)%myJ     = jH
                 nodeHalo1st(ii)%myK     = kH

                 nodeHalo1st(ii)%donorProc  = neighProc(mm)
                 nodeHalo1st(ii)%donorBlock = neighBlock(mm)

                 nodeHalo1st(ii)%dI = iD
                 nodeHalo1st(ii)%dJ = jD
                 nodeHalo1st(ii)%dK = kD

                 ! Store the short hand of the transformation matrix
                 ! for this halo.

                 transformNode(ii,1) = l1(mm)
                 transformNode(ii,2) = l2(mm)
                 transformNode(ii,3) = l3(mm)

                 ! It is possible that ii is treated earlier and hence
                 ! periodic info may have been stored. Remove this.

                 if(nodeHalo1st(ii)%nPeriodicSubfaces > 0) then
                   deallocate(nodeHalo1st(ii)%periodicSubfaces, stat=ierr)
                   if(ierr /= 0)                              &
                     call terminate("determineFaceHalos",     &
                                    "Deallocation failure for &
                                    &periodicSubfaces")
                   nullify(nodeHalo1st(ii)%periodicSubfaces)
                   nodeHalo1st(ii)%nPeriodicSubfaces = 0
                 endif

                 ! If the subface is periodic store the periodic info.

                 if(indexPeriodic > 0) then
                   nodeHalo1st(ii)%nPeriodicSubfaces = 1
                   allocate(nodeHalo1st(ii)%periodicSubfaces(1), &
                            stat=ierr)
                   if(ierr /= 0)                                   &
                     call terminate("determineFaceHalos",          &
                                    "Memory allocation failure for &
                                    &periodicSubfaces")
                   nodeHalo1st(ii)%periodicSubfaces(1) = indexPeriodic
                 endif

               enddo
             enddo
           enddo
!
!          **************************************************************
!          *                                                            *
!          * The cell halo's belonging to this subface. Direct cell     *
!          * halo's are unique and therefore info cannot already be     *
!          * written earlier.                                           *
!          *                                                            *
!          **************************************************************
!
           ! First determine the cell range of the donor block on
           ! the subface. This equals the nodal range, except that 1 is
           ! added to the smallest index. As it is possible that the
           ! index is running negatively, this should be taken into account.

           donorCellRange(1,1) = dinBeg(mm)
           donorCellRange(2,1) = djnBeg(mm)
           donorCellRange(3,1) = dknBeg(mm)

           donorCellRange(1,2) = dinEnd(mm)
           donorCellRange(2,2) = djnEnd(mm)
           donorCellRange(3,2) = dknEnd(mm)

           ! The loop to add 1 to the lowest index and to correct the
           ! index corresponding to the face we are on.

           do i=1,3
             if(donorCellRange(i,1) == donorCellRange(i,2)) then

               ! If the face corresponds to a min face, indicated by
               ! donorCellRange(i,1) == 1 then 1 must be added;
               ! otherwise nothing needs to be done.

               if(donorCellRange(i,1) == 1) then
                 donorCellRange(i,1) = 2
                 donorCellRange(i,2) = 2
               endif

             else if(donorCellRange(i,1) > donorCellRange(i,2)) then
               donorCellRange(i,2) = donorCellRange(i,2) + 1
             else
               donorCellRange(i,1) = donorCellRange(i,1) + 1
             endif
           enddo

           ! Loop over the halo cells located adjacent to the subface.

           do k=myCellRange(3,1),myCellRange(3,2),step(3)
             do j=myCellRange(2,1),myCellRange(2,2),step(2)
               do i=myCellRange(1,1),myCellRange(1,2),step(1)

                 ! Check in debug mode whether this halo is already
                 ! stored. This should not be the case.

                 if( debug ) then
                   if(cellIndex(nn)%entryList(i,j,k) /= 0) &
                     call terminate("determineFaceHalos",  &
                                    "internal cell halo already stored")
                 endif

                 ! Determine the indices of the donor point by applying
                 ! the transformation matrix to i,j,k.

                 ii = i - myCellRange(1,1)
                 jj = j - myCellRange(2,1)
                 kk = k - myCellRange(3,1)

                 iD = donorCellRange(1,1) &
                    + tMat(1,1)*ii + tMat(1,2)*jj + tMat(1,3)*kk
                 jD = donorCellRange(2,1) &
                    + tMat(2,1)*ii + tMat(2,2)*jj + tMat(2,3)*kk
                 kD = donorCellRange(3,1) &
                    + tMat(3,1)*ii + tMat(3,2)*jj + tMat(3,3)*kk

                 ! Update the counter iicell1st and store its value a
                 ! bit easier in ii and set entryList accordingly.

                 iicell1st = iicell1st +1
                 ii        = iicell1st

                 cellIndex(nn)%entryList(i,j,k) = ii

                 ! Store the info in the correct place in cellHalo1st.

                 cellHalo1st(ii)%myBlock = nn
                 cellHalo1st(ii)%myI     = i
                 cellHalo1st(ii)%myJ     = j
                 cellHalo1st(ii)%myK     = k

                 cellHalo1st(ii)%donorProc  = neighProc(mm)
                 cellHalo1st(ii)%donorBlock = neighBlock(mm)

                 cellHalo1st(ii)%dI = iD
                 cellHalo1st(ii)%dJ = jD
                 cellHalo1st(ii)%dK = kD

                 ! Store the short hand of the transformation matrix
                 ! for this halo.

                 transformCell(ii,1) = l1(mm)
                 transformCell(ii,2) = l2(mm)
                 transformCell(ii,3) = l3(mm)

                 ! If the subface is periodic store the periodic info.
                 ! Note that for the cells it is not needed to check
                 ! if a previous transformation was already stored,
                 ! because cell halo's are unique.

                 if(indexPeriodic > 0) then
                   cellHalo1st(ii)%nPeriodicSubfaces = 1
                   allocate(cellHalo1st(ii)%periodicSubfaces(1), &
                            stat=ierr)
                   if(ierr /= 0)                                   &
                     call terminate("determineFaceHalos",          &
                                    "Memory allocation failure for &
                                    &periodicSubfaces")
                   cellHalo1st(ii)%periodicSubfaces(1) = indexPeriodic
                 endif

               enddo
             enddo
           enddo

         enddo n1to1Loop

       enddo domains

       end subroutine determineFaceHalos

!      ==================================================================

       subroutine haloRanges(mm, offset, cellRange, nodeRange, step)
!
!      ******************************************************************
!      *                                                                *
!      * haloRanges determines the cell and nodal ranges for the given  *
!      * subface as well as the direction normal to the subface,        *
!      * pointing outwards. In case of negative running indices of the  *
!      * subface, step is set to -1; otherwise it is 1.                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: mm
       integer(kind=intType), dimension(3), intent(out) :: offset, step
       integer(kind=intType), dimension(3,2), intent(out) :: cellRange
       integer(kind=intType), dimension(3,2), intent(out) :: nodeRange
!
!      Local variables.
!
       integer(kind=intType) :: i
       integer(kind=intType) :: cellHaloInd, cellHaloID
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the offset in i, j and k direction, depending on the
       ! faceID. This can be interpreted as the outward pointing normal
       ! in the index domain. Furthermore store cellHaloInd and
       ! cellHaloID. This info is needed to construct the cell halo
       ! range correctly.

       offset = 0

       select case (BCFaceID(mm))
         case (iMin)
           offset(1)   = -1
           cellHaloInd =  1
           cellHaloID  =  1
         case (iMax)
           offset(1)   =  1
           cellHaloInd =  1
           cellHaloID  =  ie
         case (jMin)
           offset(2)   = -1
           cellHaloInd =  2
           cellHaloID  =  1
         case (jMax)
           offset(2)   =  1
           cellHaloInd =  2
           cellHaloID  =  je
         case (kMin)
           offset(3)   = -1
           cellHaloInd =  3
           cellHaloID  =  1
         case (kMax)
           offset(3)   =  1
           cellHaloInd =  3
           cellHaloID  =  ke
       end select

       ! Copy the nodal range.

       nodeRange(1,1) = inBeg(mm)
       nodeRange(2,1) = jnBeg(mm)
       nodeRange(3,1) = knBeg(mm)

       nodeRange(1,2) = inEnd(mm)
       nodeRange(2,2) = jnEnd(mm)
       nodeRange(3,2) = knEnd(mm)

       ! Determine the cell range. The cell numbering of a block starts
       ! at index 2, i.e. 1 higher than the node numbering. Consequently
       ! 1 must be added to the smallest indices of the nodal range.
       ! Take negative running indices into account and set step
       ! accordingly.

       do i=1,3
         if(nodeRange(i,1) > nodeRange(i,2)) then

           ! Negative running index.

           step(i)        = -1
           cellRange(i,1) = nodeRange(i,1)
           cellRange(i,2) = nodeRange(i,2) + 1
         else

           ! Positive running index.

           step(i)        = 1
           cellRange(i,1) = nodeRange(i,1) + 1
           cellRange(i,2) = nodeRange(i,2)
         endif
       enddo

       ! Correct the cell range for the index corresponding to the face
       ! we are on.

       cellRange(cellHaloInd,1) = cellHaloID
       cellRange(cellHaloInd,2) = cellHaloID

       end subroutine haloRanges

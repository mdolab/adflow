!
!      ******************************************************************
!      *                                                                *
!      * File:          loadBalance.f90                                 *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 02-05-2003                                      *
!      * Last modified: 10-27-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine loadBalance
!
!      ******************************************************************
!      *                                                                *
!      * loadBalance determines the mapping of the blocks onto the      *
!      * processors. If the user allows so blocks my be split to obtain *
!      * a better load balance.                                         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use cgnsGrid
       use communication
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use iteration
       use partitionMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, nn, mm, ii, jj, kk
       integer(kind=intType) :: nViscBocos

       integer(kind=intType), dimension(0:nProc-1) :: nBlockPerProc

       integer(kind=intType), dimension(:), allocatable :: oldSubfaceID

       type(subblocksOfCGNSType), dimension(:), allocatable :: &
                                                          subblocksOfCGNS
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the block distribution over the processors.

       !nProc = 64
       !do nProc=250,400 !750,1000 !200,400
       !  write(*,*) 'nProc: ', nProc
       call blockDistribution
       !enddo
       !call terminate("loadBalance","Hack for distribution test")
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local block info.                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize nBlockPerProc to 0.

       nBlockPerProc = 0

       ! Determine the number of blocks the current processor will store
       ! and the local block number for every block.

       nDom = 0
       do i=1,nBlocks
         if(part(i) == myID) nDom = nDom +1

         nBlockPerProc(part(i)) = nBlockPerProc(part(i)) + 1
         blocks(i)%blockID = nBlockPerProc(part(i))
       enddo

       ! Allocate the memory for flowDoms and initialize its pointers
       ! to null pointers.

       call initFlowDoms

       ! Repeat the loop, but now store the info of the blocks
       ! in flowDoms. Store the number of time intervals for the spectral
       ! method a bit easier in mm. Note that this number is 1 for the
       ! steady and unsteady modes.

       nn = 0
       mm = nTimeIntervalsSpectral
       domains: do i=1,nBlocks
         myBlock: if(part(i) == myID) then

           ! Update the counter nn.

           nn = nn + 1

           ! Copy the dimensions of the block.

           flowDoms(nn,1,1:mm)%nx = blocks(i)%nx
           flowDoms(nn,1,1:mm)%ny = blocks(i)%ny
           flowDoms(nn,1,1:mm)%nz = blocks(i)%nz

           flowDoms(nn,1,1:mm)%il = blocks(i)%il
           flowDoms(nn,1,1:mm)%jl = blocks(i)%jl
           flowDoms(nn,1,1:mm)%kl = blocks(i)%kl

           ! The number of single halo quantities.

           flowDoms(nn,1,1:mm)%ie = blocks(i)%il + 1
           flowDoms(nn,1,1:mm)%je = blocks(i)%jl + 1
           flowDoms(nn,1,1:mm)%ke = blocks(i)%kl + 1

           ! The number of double halo quantities.

           flowDoms(nn,1,1:mm)%ib = blocks(i)%il + 2
           flowDoms(nn,1,1:mm)%jb = blocks(i)%jl + 2
           flowDoms(nn,1,1:mm)%kb = blocks(i)%kl + 2

           ! Relation to the original cgns grid.

           flowDoms(nn,1,1:mm)%cgnsBlockID = blocks(i)%cgnsBlockID

           flowDoms(nn,1,1:mm)%iBegOr = blocks(i)%iBegOr
           flowDoms(nn,1,1:mm)%jBegOr = blocks(i)%jBegOr
           flowDoms(nn,1,1:mm)%kBegOr = blocks(i)%kBegOr

           flowDoms(nn,1,1:mm)%iEndOr = blocks(i)%iEndOr
           flowDoms(nn,1,1:mm)%jEndOr = blocks(i)%jEndOr
           flowDoms(nn,1,1:mm)%kEndOr = blocks(i)%kEndOr

           ! Determine whether or not the block is moving.
           ! First initialize it to gridMotionSpecified. This is
           ! .true. if a rigid body motion was specified for the
           ! entire grid; otherwise it is .false.

           flowDoms(nn,1,1:mm)%blockIsMoving = gridMotionSpecified

           ! Check whether the corresponding cgns block is moving.
           ! Although it is possible that boundaries of a block rotate
           ! differently than the block itself, this should not be
           ! taken into account here; that's a matter of BC's.
           ! Here only the internal block structure is looked at.

           k = flowDoms(nn,1,1)%cgnsBlockID
           if( cgnsDoms(k)%rotatingFrameSpecified ) &
             flowDoms(nn,1,1:mm)%blockIsMoving = .true.

           ! For an unsteady computation on a deforming mesh
           ! blockIsMoving is always .true. Note that the time spectral
           ! method is also an unsteady computation.

           if(deforming_Grid .and.            &
              (equationMode == unsteady .or. &
               equationMode == timeSpectral)) &
             flowDoms(nn,1,1:mm)%blockIsMoving = .true.

           ! Set addGridVelocities to blockIsMoving. This could be
           ! overwritten later when the code is running in python mode.

           flowDoms(nn,1,1:mm)%addGridVelocities = &
                 flowDoms(nn,1,1:mm)%blockIsMoving

           ! Set the number of subfaces and allocate the memory for the
           ! subface info. Note that this memory is only allocated for
           ! the first spectral time value; the other ones are identical.

           flowDoms(nn,1,1:mm)%nBocos   = blocks(i)%nBocos
           flowDoms(nn,1,1:mm)%n1to1    = blocks(i)%n1to1
           flowDoms(nn,1,1:mm)%nSubface = blocks(i)%nSubface
           j                           = blocks(i)%nSubface

           allocate(flowDoms(nn,1,1)%BCType(j),      &
                    flowDoms(nn,1,1)%BCFaceID(j),    &
                    flowDoms(nn,1,1)%cgnsSubface(j), &
                    flowDoms(nn,1,1)%inBeg(j),       &
                    flowDoms(nn,1,1)%jnBeg(j),       &
                    flowDoms(nn,1,1)%knBeg(j),       &
                    flowDoms(nn,1,1)%inEnd(j),       &
                    flowDoms(nn,1,1)%jnEnd(j),       &
                    flowDoms(nn,1,1)%knEnd(j),       &
                    flowDoms(nn,1,1)%dinBeg(j),      &
                    flowDoms(nn,1,1)%djnBeg(j),      &
                    flowDoms(nn,1,1)%dknBeg(j),      &
                    flowDoms(nn,1,1)%dinEnd(j),      &
                    flowDoms(nn,1,1)%djnEnd(j),      &
                    flowDoms(nn,1,1)%dknEnd(j),      &
                    flowDoms(nn,1,1)%neighProc(j),   &
                    flowDoms(nn,1,1)%neighBlock(j),  &
                    flowDoms(nn,1,1)%l1(j),          &
                    flowDoms(nn,1,1)%l2(j),          &
                    flowDoms(nn,1,1)%l3(j),          &
                    flowDoms(nn,1,1)%groupNum(j),    stat=ierr)

           if(ierr /= 0)                   &
             call terminate("loadBalance", &
                            "Memory allocation failure for subface info")

           ! Determine the new numbering of the boundary subfaces, such
           ! that the viscous subfaces are numbered first, followed by
           ! the inViscid subfaces, etc.

           allocate(oldSubfaceID(blocks(i)%nBocos), stat=ierr)
           if(ierr /= 0)                   &
             call terminate("loadBalance", &
                            "Memory allocation failure for oldSubfaceID")

           call sortSubfaces(oldSubfaceID, blocks(i))

           ! Initialize the number of viscous boundary subfaces to 0.

           nViscBocos = 0

           ! Copy the info. Set the neighboring proc and block id to -1
           ! and 0 repectively in this loop. This is okay for boundary
           ! faces, but must be corrected for the internal block
           ! boundaries.

           do j=1,blocks(i)%nSubface

             ! Store the old subface id in k. For boundary faces the
             ! sorting is taken into account; for 1 to 1 subfaces the
             ! number is identical to the subface id in block.

             k = j
             if(j <= blocks(i)%nBocos) k = oldSubfaceID(j)

             ! Copy the info.

             flowDoms(nn,1,1)%BCType(j)      = blocks(i)%BCType(k)
             flowDoms(nn,1,1)%BCFaceID(j)    = blocks(i)%BCFaceID(k)
             flowDoms(nn,1,1)%cgnsSubface(j) = blocks(i)%cgnsSubface(k)

             flowDoms(nn,1,1)%inBeg(j) = blocks(i)%inBeg(k)
             flowDoms(nn,1,1)%jnBeg(j) = blocks(i)%jnBeg(k)
             flowDoms(nn,1,1)%knBeg(j) = blocks(i)%knBeg(k)
             flowDoms(nn,1,1)%inEnd(j) = blocks(i)%inEnd(k)
             flowDoms(nn,1,1)%jnEnd(j) = blocks(i)%jnEnd(k)
             flowDoms(nn,1,1)%knEnd(j) = blocks(i)%knEnd(k)

             flowDoms(nn,1,1)%dinBeg(j) = blocks(i)%dinBeg(k)
             flowDoms(nn,1,1)%djnBeg(j) = blocks(i)%djnBeg(k)
             flowDoms(nn,1,1)%dknBeg(j) = blocks(i)%dknBeg(k)
             flowDoms(nn,1,1)%dinEnd(j) = blocks(i)%dinEnd(k)
             flowDoms(nn,1,1)%djnEnd(j) = blocks(i)%djnEnd(k)
             flowDoms(nn,1,1)%dknEnd(j) = blocks(i)%dknEnd(k)

             flowDoms(nn,1,1)%neighProc(j)  = -1
             flowDoms(nn,1,1)%neighBlock(j) =  0

             flowDoms(nn,1,1)%l1(j) = blocks(i)%l1(k)
             flowDoms(nn,1,1)%l2(j) = blocks(i)%l2(k)
             flowDoms(nn,1,1)%l3(j) = blocks(i)%l3(k)

             flowDoms(nn,1,1)%groupNum(j) = blocks(i)%groupNum(k)

             ! Update the number of viscous boundaries if this
             ! is a viscous subface.

             if(flowDoms(nn,1,1)%BCType(j) == NSWallAdiabatic .or. &
                flowDoms(nn,1,1)%BCType(j) == NSWallIsothermal)    &
               nViscBocos = nViscBocos + 1

           enddo

           flowDoms(nn,1,1:mm)%nViscBocos = nViscBocos

           ! Correct the neighboring block and proc ID for internal
           ! block boundaries.

           do k=1,blocks(i)%n1to1
             j = blocks(i)%nBocos + k

             flowDoms(nn,1,1)%neighProc(j)  = part(blocks(i)%neighBlock(j))
             flowDoms(nn,1,1)%neighBlock(j) = &
                           blocks(blocks(i)%neighBlock(j))%blockID
           enddo

           ! Release the memory of oldSubfaceID.

           deallocate(oldSubfaceID, stat=ierr)
           if(ierr /= 0)                   &
             call terminate("loadBalance", &
                            "Deallocation error for oldSubfaceID")

           ! Allocate memory for the overset boundary info. Note the
           ! number of orphans is set to 0 because there is no support
           ! for input of orphans from CGNS or Plot3D yet.

           flowDoms(nn,1,1)%nCellsOverset = blocks(i)%nCellsOverset
           flowDoms(nn,1,1)%nOrphans      = 0

           j = blocks(i)%nCellsOverset
           k = nDonorWeights(oversetInterpType)

           allocate(flowDoms(nn,1,1)%ibndry(3,j),        &
                    flowDoms(nn,1,1)%idonor(3,j),        &
                    flowDoms(nn,1,1)%overint(k,j),       &
                    flowDoms(nn,1,1)%neighProcOver(j),   &
                    flowDoms(nn,1,1)%neighBlockOver(j),  &
                    stat=ierr)
           if(ierr /= 0)                   &
             call terminate("loadBalance", &
                           "Memory allocation failure for overset info")

           ! Find the neighboring processor and the local block id
           ! on that processor as in the internal subface case.

           do j = 1,blocks(i)%nCellsOverset
             flowDoms(nn,1,1)%neighProcOver(j) = &
                                            part(blocks(i)%neighOver(j))
             flowDoms(nn,1,1)%neighBlockOver(j) = &
                                 blocks(blocks(i)%neighOver(j))%blockId
           enddo

           ! Transfer the boundary and donor indices from the cgns doms
           ! to the flowDoms. Note the conversion to the local block
           ! indices, and the additional +1 to convert from a cgns cell
           ! index to the style used in this code.

           ii = blocks(i)%cgnsBlockId
           do j = 1,blocks(i)%nCellsOverset
             jj = blocks(i)%cgnsOver(j)
             kk = blocks(i)%ipntOver(j)

             flowDoms(nn,1,1)%ibndry(1,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%ibndry(1,kk) &
                             - blocks(i)%iBegOr + 2
             flowDoms(nn,1,1)%ibndry(2,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%ibndry(2,kk) &
                             - blocks(i)%jBegOr + 2
             flowDoms(nn,1,1)%ibndry(3,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%ibndry(3,kk) &
                             - blocks(i)%kBegOr + 2

             flowDoms(nn,1,1)%idonor(1,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%idonor(1,kk) &
                             - blocks(blocks(i)%neighOver(j))%iBegOr + 2
             flowDoms(nn,1,1)%idonor(2,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%idonor(2,kk) &
                             - blocks(blocks(i)%neighOver(j))%jBegOr + 2
             flowDoms(nn,1,1)%idonor(3,j) =                            &
                               cgnsDoms(ii)%connOver(jj)%idonor(3,kk) &
                             - blocks(blocks(i)%neighOver(j))%kBegOr + 2

             ! If the fine grid interpolants were input, check if they
             ! are the actual weights or the parametric coordinates. In
             ! the latter case, convert them to the weights.

             if (.not. oversetDonorsAreGuesses) then
               if (ubound(cgnsDoms(ii)%connOver(jj)%interp,1) == 3) then
                 call getWeights(cgnsDoms(ii)%connOver(jj)%interp(:,kk), &
                                 flowDoms(nn,1,1)%overint(:,j))
               else
                 flowDoms(nn,1,1)%overint(:,j) = &
                               cgnsDoms(ii)%connOver(jj)%interp(:,kk)
               end if
             end if
           end do

         endif myBlock

         ! Release the memory of the subface and overset info in
         ! this block.

         deallocate(blocks(i)%bcType,      blocks(i)%bcFaceid,    &
                    blocks(i)%cgnsSubface, blocks(i)%inBeg,       &
                    blocks(i)%jnBeg,       blocks(i)%knBeg,       &
                    blocks(i)%inEnd,       blocks(i)%jnEnd,       &
                    blocks(i)%knEnd,       blocks(i)%dinBeg,      &
                    blocks(i)%djnBeg,      blocks(i)%dknBeg,      &
                    blocks(i)%dinEnd,      blocks(i)%djnEnd,      &
                    blocks(i)%dknEnd,      blocks(i)%neighBlock,  &
                    blocks(i)%l1,          blocks(i)%l2,          &
                    blocks(i)%l3,          blocks(i)%groupNum,    &
                    blocks(i)%cgnsOver,    blocks(i)%ipntOver,    &
                    blocks(i)%neighOver,   blocks(i)%overComm,    &
                    stat=ierr)
         if(ierr /= 0)                   &
           call terminate("loadBalance", &
                          "Deallocation error for boundary info")
       enddo domains

!      ******************************************************************
!      *                                                                *
!      * Determine the number of processors, the processor ID's on      *
!      * which the original cgns blocks are stored, the local           *
!      * block ID's and the nodal ranges of the subblocks. As blocks    *
!      * can be split during run-time, multiple processors can store a  *
!      * part of original block.                                        *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for subblocksOfCGNS.

       allocate(subblocksOfCGNS(nBlocks), stat=ierr)
       if(ierr /= 0)                   &
         call terminate("loadBalance", &
                        "Memory allocation failure for subblocksOfCGNS")

       ! Copy the data into subblocksOfCGNS.

       do nn=1,nBlocks
         subblocksOfCGNS(nn)%cgnsBlockID = blocks(nn)%cgnsBlockID
         subblocksOfCGNS(nn)%procID      = part(nn)
         subblocksOfCGNS(nn)%blockID     = blocks(nn)%blockID

         subblocksOfCGNS(nn)%iBegOr = blocks(nn)%iBegOr
         subblocksOfCGNS(nn)%iEndOr = blocks(nn)%iEndOr
         subblocksOfCGNS(nn)%jBegOr = blocks(nn)%jBegOr
         subblocksOfCGNS(nn)%jEndOr = blocks(nn)%jEndOr
         subblocksOfCGNS(nn)%kBegOr = blocks(nn)%kBegOr
         subblocksOfCGNS(nn)%kEndOr = blocks(nn)%kEndOr
       enddo

       ! Sort subblocksOfCGNS in increasing order.

       call qsortSubblocksOfCGNSType(subblocksOfCGNS, nBlocks)

       ! Loop over the number of cgns blocks and find out the number of
       ! subblocks it contains.

       ii = 1
       subBlockLoop: do nn=1,cgnsNDom

         ! Determine the ending index jj in subblocksOfCGNS for this
         ! CGNS block. The starting index is ii.

         if(nn == cgnsNDom) then
           jj = nBlocks
         else
           jj = ii
           do
             if(subblocksOfCGNS(jj+1)%cgnsBlockID > nn) exit
             jj = jj + 1
           enddo
         endif

         ! Set nSubBlocks and allocate the memory for procStored,
         ! localBlockID, iBegOr, iEndOr, etc.

         cgnsDoms(nn)%nSubBlocks = jj - ii + 1
         k = cgnsDoms(nn)%nSubBlocks
         allocate(cgnsDoms(nn)%iBegOr(k),       cgnsDoms(nn)%iEndOr(k), &
                  cgnsDoms(nn)%jBegOr(k),       cgnsDoms(nn)%jEndOr(k), &
                  cgnsDoms(nn)%kBegOr(k),       cgnsDoms(nn)%kEndOr(k), &
                  cgnsDoms(nn)%procStored(k),                           &
                  cgnsDoms(nn)%localBlockID(k), stat=ierr)
         if(ierr /= 0)                   &
           call terminate("loadBalance", &
                          "Memory allocation failure for procStored, &
                          &localBlockID, iBegOr, iEndOr, etc.")

         ! Copy the processor ID's, the local block ID's
         ! and the subranges.

         do i=1,cgnsDoms(nn)%nSubBlocks
           j = i + ii - 1
           cgnsDoms(nn)%procStored(i)   = subblocksOfCGNS(j)%procID
           cgnsDoms(nn)%localBlockID(i) = subblocksOfCGNS(j)%blockID

           cgnsDoms(nn)%iBegOr(i) = subblocksOfCGNS(j)%iBegOr
           cgnsDoms(nn)%iEndOr(i) = subblocksOfCGNS(j)%iEndOr
           cgnsDoms(nn)%jBegOr(i) = subblocksOfCGNS(j)%jBegOr
           cgnsDoms(nn)%jEndOr(i) = subblocksOfCGNS(j)%jEndOr
           cgnsDoms(nn)%kBegOr(i) = subblocksOfCGNS(j)%kBegOr
           cgnsDoms(nn)%kEndOr(i) = subblocksOfCGNS(j)%kEndOr
         enddo

         ! Set ii for the next CGNS block.

         ii = jj + 1

       enddo subBlockLoop

       ! Release the memory of blocks, part and subblocksOfCGNS.

       deallocate(blocks, part, subblocksOfCGNS, stat=ierr)
       if(ierr /= 0)                   &
         call terminate("loadBalance", &
                        "Deallocation error for blocks, part and &
                        &subblocksOfCGNS")

      !j = 20+myID
      !do nn=1,ndom
      !  write(j,"(8I4)") nn, flowDoms(nn,1,1)%cgnsBlockID, &
      !                   flowDoms(nn,1,1)%iBegOr, flowDoms(nn,1,1)%iEndOr, &
      !                   flowDoms(nn,1,1)%jBegOr, flowDoms(nn,1,1)%jEndOr, &
      !                   flowDoms(nn,1,1)%kBegOr, flowDoms(nn,1,1)%kEndOr
      !enddo

       end subroutine loadBalance

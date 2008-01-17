!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSurfIndList.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-04-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdCreateSurfIndList(famID,startInd,endInd)
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateSurfIndList creates the list of cgns block ID's        *
!      * and indices of the surface for the given family ID. If the     *
!      * family ID == 0, the list contains all points on the solid      *
!      * surfaces. There is no need to specify the spectral solution,   *
!      * because the indices are the same for all of them.              *
!      * On return the variables startInd and endInd are set to the     *
!      * range this family occupies in the array of all families.       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use mdData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: famID
       integer(kind=intType), intent(out) :: startInd, endInd
!
!      Local variables.
!
       integer :: ierr, size
       integer, dimension(nProc) :: recvcounts, displs

       integer(kind=intType) :: ii, jj, mm, nn, i, j, k
       integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
       integer(kind=intType) :: nSurfNodesLoc, modFamID

       integer(kind=intType), dimension(:,:), allocatable :: indLoc

       logical :: storeSubface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Perform a check to see if this routine is called correctly.
       ! If not, terminate the program.

       if(famID == 0 .and. cgnsNfamilies > 0) then
         if(myID == 0)                           &
           call terminate("mdCreateSurfIndList", &
                          "Family ID 0 is only allowed when no family &
                          &info is present in the grid")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Determine the number of surface nodes per family if this
       ! information is not available.

       if(.not. allocated(mdNSurfNodes)) call mdCreateNSurfNodes

       ! Allocate the memory for the local surface ID's. ModFamID is
       ! introduced to take famID == 0 into account.

       modFamID = max(famID, 1_intType)
       nSurfNodesLoc = mdNSurfNodes(myID+1,modFamID) &
                     - mdNSurfNodes(myID,  modFamID)

       allocate(indLoc(4,nSurfNodesLoc), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("mdCreateSurfIndList", &
                        "Memory allocation failure for indLoc")

       ! Store the indices of the local surface nodes for the
       ! given family ID.

       ii = 0
       domains: do nn=1,nDom

         ! Have the pointers in blockPointers point to the 1st
         ! spectral solution of this block on the finest mg level.
         ! There is no need to distinguish between the different
         ! spectral solutions, because this info is the same
         ! for all of them.

         call setPointers(nn,1_intType,1_intType)

         ! Loop over the number of boundary subfaces of this block.

         bocos: do mm=1,nBocos

           ! Check if the data of this subface must be stored.

           storeSubface = .false.

           if(famID == 0) then

             ! No family info present; all solid wall points are stored.

             if(BCType(mm) == EulerWall       .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) storeSubface = .true.

           else

             ! Family info is present. Check if this subface belongs
             ! to the given familyID.

             jj = cgnsSubface(mm)
             jj = cgnsDoms(nbkGlobal)%bocoInfo(jj)%familyID
             if(jj == famID) storeSubface = .true.

           endif

           ! Store the data of this subface, if needed.

           storeSubfaceTest: if( storeSubface ) then

             ! Subface must be stored. Determine the face ID on which
             ! this subface is located and store its dimensions.

             select case (BCFaceID(mm))

               case (iMin)
                 kBeg = knBeg(mm); jBeg = jnBeg(mm); iBeg = 1
                 kEnd = knEnd(mm); jEnd = jnEnd(mm); iEnd = 1

               case (iMax)
                 kBeg = knBeg(mm); jBeg = jnBeg(mm); iBeg = il
                 kEnd = knEnd(mm); jEnd = jnEnd(mm); iEnd = il

               case (jMin)
                 kBeg = knBeg(mm); jBeg = 1; iBeg = inBeg(mm)
                 kEnd = knEnd(mm); jEnd = 1; iEnd = inEnd(mm)

               case (jMax)
                 kBeg = knBeg(mm); jBeg = jl; iBeg = inBeg(mm)
                 kEnd = knEnd(mm); jEnd = jl; iEnd = inEnd(mm)

               case (kMin)
                 kBeg = 1; jBeg = jnBeg(mm); iBeg = inBeg(mm)
                 kEnd = 1; jEnd = jnEnd(mm); iEnd = inEnd(mm)

               case (kMax)
                 kBeg = kl; jBeg = jnBeg(mm); iBeg = inBeg(mm)
                 kEnd = kl; jEnd = jnEnd(mm); iEnd = inEnd(mm)

             end select

             ! Loop over the nodes of the subface and store its indices
             ! in indLoc. Note that the original cgns block could have
             ! been split. This is taken into account via the offsets
             ! iBegor, jBegor, kBegor. The first three entries are used
             ! to store the i, j and k index of the node and the fourth
             ! is used for the original block ID.

             do k=kBeg,kEnd
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   ii = ii + 1

                   indLoc(1,ii) = i - 1 + iBegor
                   indLoc(2,ii) = j - 1 + jBegor
                   indLoc(3,ii) = k - 1 + kBegor
                   indLoc(4,ii) = nbkGlobal
                 enddo
               enddo
             enddo

           endif storeSubfaceTest
         enddo bocos
       enddo domains

       ! Test if the memory of mdSurfInd has already been allocated.
       ! If not, allocate it.

       if(.not. allocated(mdSurfInd) ) then

         jj = mdNSurfNodes(nProc,max(cgnsNfamilies,1_intType))
         allocate(mdSurfInd(4,jj), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("mdCreateSurfIndList", &
                          "Memory allocation failure for mdSurfInd")
       endif

       ! Construct the arrays recvcounts and displs needed for
       ! allgatherv.

       do nn=1,nProc
         recvcounts(nn) = 4*(mdNSurfNodes(nn,  modFamID) &
                        -    mdNSurfNodes(nn-1,modFamID))
         displs(nn)     = 4* mdNSurfNodes(nn-1,modFamID)
       enddo

       ! Call allgatherv to gather the data on processor 0.

       size = 4*nSurfNodesLoc
       call mpi_allgatherv(indLoc, size, sumb_integer, mdSurfInd, &
                           recvcounts, displs, sumb_integer,      &
                           SUmb_comm_world, ierr)

       ! Release the memory of indLoc.

       deallocate(indLoc, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("mdCreateSurfIndList", &
                        "Deallocation failure for indLoc")

       ! Set the values of startInd and endInd, which give the range
       ! in the array mdSurfInd where the info for the given family
       ! is stored.

       startInd = mdNSurfNodes(0,    modFamID) + 1
       endInd   = mdNSurfNodes(nProc,modFamID)

       end subroutine mdCreateSurfIndList

!      ==================================================================

       subroutine mdDeleteSurfIndList
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteSurfIndList deallocates the memory of mdSurfInd.       *
!      *                                                                *
!      ******************************************************************
!
       use mdData
       implicit none
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Deallocate the memory of mdSurfInd if it has been allocated.

       if( allocated(mdSurfInd) ) then

         deallocate(mdSurfInd, stat=ierr)
         if(ierr /= 0)                           &
           call terminate("mdDeleteSurfIndList", &
                          "Deallocation error for mdSurfInd")
       endif

       end subroutine mdDeleteSurfIndList

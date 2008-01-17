!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSurfVarList.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-26-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdCreateSurfVarList(sps,famID,startInd,endInd)
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateSurfVarList creates the list of values for the         *
!      * specified variable of the surface nodes for the given          *
!      * family ID. If the family ID == 0, the list contains all points *
!      * on the solid surfaces.                                         *
!      * On return the variables startInd and endInd are set to the     *
!      * range this family occupies in the array of all families.       *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use mdData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: sps, famID
       integer(kind=intType), intent(out) :: startInd, endInd
!
!      Local variables.
!
       integer :: ierr, size
       integer, dimension(nProc) :: recvcounts, displs

       integer(kind=intType) :: jj, nn, nSurfNodesLoc, modFamID

       real(kind=realType), dimension(:), allocatable :: valLoc
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
           call terminate("mdCreateSurfValList", &
                          "Family ID 0 is only allowed when no family &
                          &info is present in the grid")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Determine the number of surface nodes per family if this
       ! information is not available.

       if(.not. allocated(mdNSurfNodes)) call mdCreateNSurfNodes

       ! Allocate the memory for the local values and initialize it
       ! to zero.  ModFamID is introduced to take famID == 0
       ! into account.

       modFamID = max(famID, 1_intType)
       nSurfNodesLoc = mdNSurfNodes(myID+1,modFamID) &
                     - mdNSurfNodes(myID,  modFamID)

       allocate(valLoc(nSurfNodesLoc), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("mdCreateSurfValList", &
                        "Memory allocation failure for valLoc")

       ! Determine the variable to be stored and call the
       ! corresponding routine.

       valLoc = zero
       call mdStoreLocalCp(valLoc, sps, famID)

       ! Test if the memory of mdSurfVal has already been allocated.
       ! If not, allocate it.

       if(.not. allocated(mdSurfVal) ) then

         jj = mdNSurfNodes(nProc,max(cgnsNfamilies,1_intType))
         allocate(mdSurfVal(jj), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("mdCreateSurfValList", &
                          "Memory allocation failure for &
                          &mdSurfVal")
       endif

       ! Construct the arrays recvcounts and displs needed for gatherv.

       do nn=1,nProc
         recvcounts(nn) = mdNSurfNodes(nn,  modFamID) &
                        - mdNSurfNodes(nn-1,modFamID)
         displs(nn)     = mdNSurfNodes(nn-1,modFamID)
       enddo

       ! Call allgatherv to gather the data.

       size = nSurfNodesLoc
       call mpi_allgatherv(valLoc, size, sumb_real, mdSurfVal, &
                           recvcounts, displs, sumb_real,      &
                           SUmb_comm_world, ierr)

       ! Release the memory of valLoc.

       deallocate(valLoc, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("mdCreateSurfValList", &
                        "Deallocation failure for valLoc")

       ! Set the values of startInd and endInd, which give the range
       ! in the array mdSurfForce where the info for the given family
       ! is stored.

       startInd = mdNSurfNodes(0,    modFamID) + 1
       endInd   = mdNSurfNodes(nProc,modFamID)

       end subroutine mdCreateSurfVarList

!      ==================================================================

       subroutine mdDeleteSurfValList
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteSurfValList deallocates the memory of mdSurfVal.       *
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
       ! Deallocate the memory of mdSurfVal if it has been allocated.

       if( allocated(mdSurfVal) ) then

         deallocate(mdSurfVal, stat=ierr)
         if(ierr /= 0)                           &
           call terminate("mdDeleteSurfValList", &
                          "Deallocation error for mdSurfVal")
       endif

       end subroutine mdDeleteSurfValList

!      ==================================================================

       subroutine mdStoreLocalCp(valLoc, sps, famID)
!
!      ******************************************************************
!      *                                                                *
!      * mdStoreLocalCp stores the value of the pressure coefficient    *
!      * of the local surface nodes of the given family ID.             *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use flowVarRefState
       use inputPhysics
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps, famID

       real(kind=realType), dimension(*), intent(inout) :: valLoc
!
!      Local variables.
!
       integer(kind=intType) :: ii, jj, mm, nn, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

       real(kind=realType) :: fact, Cp

       real(kind=realType), dimension(:,:), pointer :: pp2, pp1

       logical :: storeSubface
       logical :: updateNode1, updateNode2
       logical :: updateNode3, updateNode4
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Factor multiplying p-pInf

       fact = two/(gammainf*pInf*MachCoef*MachCoef)

       ! Loop over the domains and accumulate the nodal Cp values.

       ii = 0
       domains: do nn=1,nDom

         ! Have the pointers in blockPointers point to this block
         ! on the finest mg level.

         call setPointers(nn,1_intType,sps)

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

             ! Subface must be stored. Set the pointers pp1 and pp2
             ! to make a generic treatment possible.

             select case (BCFaceID(mm))

               case (iMin)
                 pp2 => p(2 ,1:,1:); pp1 => p(1 ,1:,1:)

               case (iMax)
                 pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:)

               case (jMin)
                 pp2 => p(1:,2 ,1:); pp1 => p(1:,1 ,1:)

               case (jMax)
                 pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:)

               case (kMin)
                 pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1)

               case (kMax)
                 pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke)

             end select

             ! Store the cell range of the subfaces a bit easier.
             ! The data must be accumulated for all nodes and therefore
             ! the first layer of halo's must be included in this loop.

             jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd + 1
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd + 1

             ! Compute the Cp value on each of the faces and scatter
             ! it to the 4 nodes. Care must be taken at the boundaries.

             do j=jBeg, jEnd
               do i=iBeg, iEnd

                 ! Compute the Cp value in the center of the boundary
                 ! face, which is an average of pp1 and pp2. Take the
                 ! scaling of 1/4 already into account for scattering
                 ! it to the nodes.

                 Cp = half*(pp2(i,j) + pp1(i,j))
                 Cp = fourth*fact*(Cp - pInf)

                 ! Determine which nodes get an update.

                 updateNode1 = .true.
                 updateNode2 = .true.
                 updateNode3 = .true.
                 updateNode4 = .true.

                 if(j == jBeg) then
                   updateNode1 = .false.
                   updateNode2 = .false.
                 else if(j == jEnd) then
                   updateNode3 = .false.
                   updateNode4 = .false.
                 endif

                 if(i == iBeg) then
                   updateNode1 = .false.
                   updateNode3 = .false.
                 else if(i == iEnd) then
                   updateNode2 = .false.
                   updateNode4 = .false.
                 endif

                 ! Scatter the value to the 4 nodes.

                 jj = ii + (j-jBeg-1)*(iEnd-iBeg) + i-iBeg
                 if(updateNode1) valLoc(jj) = valLoc(jj) + Cp

                 jj = jj + 1
                 if(updateNode2) valLoc(jj) = valLoc(jj) + Cp

                 jj = jj + iEnd - iBeg - 1
                 if(updateNode3) valLoc(jj) = valLoc(jj) + Cp

                 jj = jj + 1
                 if(updateNode4) valLoc(jj) = valLoc(jj) + Cp

               enddo
             enddo

             ! Update the counter ii.

             ii = ii + (jEnd-jBeg)*(iEnd-iBeg)

           endif storeSubfaceTest
         enddo bocos
       enddo domains

       end subroutine mdStoreLocalCp

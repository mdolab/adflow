!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSurfCoorListLocal.f90                         *
!      * Author:        C.A.(Sandy) Mader, Edwin van der Weide          *
!      * Starting date: 10-24-2007                                      *
!      * Last modified: 10-24-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdCreateSurfCoorListLocal(sps,famID,startInd,endInd)
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateSurfCoorListLocal creates the list of coordinates of   *
!      * the surface points for the given spectral solution and         *
!      * family ID on the local processor. If the family ID == 0, the   *
!      * list contains all points                                       *
!      * on the solid surfaces.                                         *
!      * On return the variables startInd and endInd are set to the     *
!      * range this family occupies in the array of all families.       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use mdDataLocal
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

       integer(kind=intType) :: ii, jj, mm, nn, i, j, k
       integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
       integer(kind=intType) :: nSurfNodesLoc, modFamID

       real(kind=realType), dimension(:,:), allocatable :: xxLoc

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
!!$       print *,'ncgnsfamilies',cgnsNfamilies
!!$       if(famID == 0 .and. cgnsNfamilies > 0) then
!!$         if(myID == 0)                            &
!!$           call terminate("mdCreateSurfCoorListLocal", &
!!$                          "Family ID 0 is only allowed when no family &
!!$                          &info is present in the grid")
!!$         call mpi_barrier(SUmb_comm_world, ierr)
!!$       endif

       ! Determine the number of surface nodes per family if this
       ! information is not available.
       !print *,'creating Nsurfnodes local', mdNSurfNodesLocal
       if(.not. allocated(mdNSurfNodesLocal)) call mdCreateNSurfNodesLocal

       ! Allocate the memory for the local surface coordinates.
       ! ModFamID is introduced to take famID == 0 into account.

       modFamID = max(famID, 1_intType)
       !print *,'famid', modFamID
       nSurfNodesLoc = mdNSurfNodesLocal(modFamID) 
       !print *,'allocating xxloc',nsurfnodesloc
       allocate(xxLoc(3,nSurfNodesLoc), stat=ierr)
       if(ierr /= 0)                                &
         call terminate("mdCreateSurfCoorList", &
                        "Memory allocation failure for xxLoc")

       ! Store the coordinates of the local surface nodes.
       !print *,'entering domain loop'
       ii = 0
       domains: do nn=1,nDom

         ! Have the pointers in blockPointers point to this block
         ! on the finest mg level.

         call setPointers(nn,1_intType,sps)

         ! Loop over the number of boundary subfaces of this block.

         bocos: do mm=1,nBocos

           ! Check if the data of this subface must be stored.

           storeSubface = .false.
           !check to see whether family boundary conditions are present
           if (cgnsDoms(nbkGlobal)%BCFamilies .eqv. .true.) then
              !BC families are present

              ! Family info is present. Check if this subface belongs
              ! to the given familyID.
              
              jj = cgnsSubface(mm)
              jj = cgnsDoms(nbkGlobal)%bocoInfo(jj)%familyID
              if(jj == famID) storeSubface = .true.

           else
              
              !if(famID == 0) then
                 
              ! No family info present; all solid wall points are stored.
              
              if(BCType(mm) == EulerWall       .or. &
                   BCType(mm) == NSWallAdiabatic .or. &
                   BCType(mm) == NSWallIsothermal) storeSubface = .true.
                            
           endif
           !print *,'storing subface',mm,nBocos
           ! Store the data of this subface, if needed.

           storeSubfaceTest: if( storeSubface ) then
            
             ! Subface must be stored. 

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

             ! Loop over the nodes of the subface and store its
             ! coordinates in xxLoc.
            
             do k=kBeg,kEnd
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   ii = ii + 1
                   xxLoc(1,ii) = x(i,j,k,1)
                   xxLoc(2,ii) = x(i,j,k,2)
                   xxLoc(3,ii) = x(i,j,k,3)
                   
                 enddo
               enddo
             enddo
            
           endif storeSubfaceTest
         enddo bocos
       enddo domains
       !print *,'looping finished'
       ! Test if the memory of mdSurfxx has already been allocated.
       ! If not, allocate it.
       !print *,'aalocating mdSurfxxlocal'
       if(.not. allocated(mdSurfxxLocal) ) then

         jj = mdNSurfNodesLocal(max(cgnsNfamilies,1_intType))
         allocate(mdSurfxxLocal(3,jj), stat=ierr)
         if(ierr /= 0)                                &
           call terminate("mdCreateSurfCoorListLocal", &
                          "Memory allocation failure for mdSurfxxLocal")
       endif

       mdSurfxxLocal(:,:)=xxLoc(:,:)

       ! Release the memory of xxLoc.
       
       deallocate(xxLoc, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("mdCreateSurfCoorList", &
                        "Deallocation failure for xxLoc")

       ! Set the values of startInd and endInd, which give the range
       ! in the array mdSurfxx where the info for the given family
       ! is stored.

       startInd = mdNSurfNodesLocal(modFamID) + 1
       endInd   = mdNSurfNodesLocal(modFamID) + nSurfNodesLoc!mdNSurfNodesLocal(modFamID+1)

     end subroutine mdCreateSurfCoorListLocal

!      ==================================================================

     subroutine mdDeleteSurfCoorListLocal
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteSurfCoorList deallocates the memory of mdSurfxx.       *
!      *                                                                *
!      ******************************************************************
!
       use mdDataLocal
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
       ! Deallocate the memory of mdSurfxx if it has been allocated.

       if( allocated(mdSurfxxLocal) ) then

         deallocate(mdSurfxxLocal, stat=ierr)
         if(ierr /= 0) &
           call terminate("mdDeleteSurfCoorList", &
                          "Deallocation error for mdSurfxx")
       endif

     end subroutine mdDeleteSurfCoorListLocal

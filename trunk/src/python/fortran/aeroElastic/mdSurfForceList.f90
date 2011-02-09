!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSurfForceList.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdCreateSurfForceList(sps,famID,startInd,endInd)
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateSurfForceList creates the list of forces of the        *
!      * the surface nodes for the given spectral solution and          *
!      * family ID. If the family ID == 0, the list contains all points *
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
       use flowVarRefState
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

       integer(kind=intType) :: ii, jj, mm, nn, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType) :: nSurfNodesLoc, modFamID

       real(kind=realType) :: scaleDim, fact, pp, fx, fy, fz
       real(kind=realType) :: tauxx, tauyy, tauzz
       real(kind=realType) :: tauxy, tauxz, tauyz
       real(kind=realType), dimension(:,:), allocatable :: forceLoc

       real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
       real(kind=realType), dimension(:,:,:), pointer :: ss

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
         if(myID == 0)                             &
           call terminate("mdCreateSurfForceList", &
                          "Family ID 0 is only allowed when no family &
                          &info is present in the grid")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Determine the number of surface nodes per family if this
       ! information is not available.

       if(.not. allocated(mdNSurfNodes)) call mdCreateNSurfNodes

       ! Allocate the memory for the local forces and initialize it
       ! to zero.  ModFamID is introduced to take famID == 0
       ! into account.

       modFamID = max(famID, 1_intType)
       nSurfNodesLoc = mdNSurfNodes(myID+1,modFamID) &
                     - mdNSurfNodes(myID,  modFamID)

       allocate(forceLoc(3,nSurfNodesLoc), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("mdCreateSurfForceList", &
                        "Memory allocation failure for forceLoc")

       forceLoc = zero

       ! Compute the scaling factor to create the correct dimensional
       ! force in newton. As the coordinates are already in meters,
       ! this scaling factor is pRef.

       scaleDim = pRef

       ! Compute the local forces. Take the scaling factor into
       ! account to obtain the forces in SI-units, i.e. Newton.

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

             ! Subface must be stored. Set a bunch of variables depending
             ! on the face ID to make a generic treatment possible.

             select case (BCFaceID(mm))

               case (iMin)
                 pp2 => p( 2,1:,1:); pp1 => p( 1,1:,1:); ss => si( 1,:,:,:)
                 fact = -one

               case (iMax)
                 pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:); ss => si(il,:,:,:)
                 fact = one

               case (jMin)
                 pp2 => p(1:, 2,1:); pp1 => p(1:, 1,1:); ss => sj(:, 1,:,:)
                 fact = -one

               case (jMax)
                 pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:); ss => sj(:,jl,:,:)
                 fact = one

               case (kMin)
                 pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1); ss => sk(:,:, 1,:)
                 fact = -one

               case (kMax)
                 pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke); ss => sk(:,:,kl,:)
                 fact = one

             end select

             ! Store the cell range of the subfaces a bit easier.
             ! As only owned faces must be considered the nodal range
             ! in BCData must be used to obtain this data.

             jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

             ! Compute the inviscid force on each of the faces and
             ! scatter it to the 4 nodes, whose forces are updated.

             do j=jBeg, jEnd
               do i=iBeg, iEnd

                 ! Compute the pressure in the center of the boundary
                 ! face, which is an average between pp2 and pp1. The
                 ! value of pp is multiplied by 1/4 (the factor to
                 ! scatter to its 4 nodes, by scaleDim (to obtain the
                 ! correct dimensional value) and by fact (which takes
                 ! the possibility of inward or outward pointing normal
                 ! into account).

                 pp = half*(pp2(i,j) + pp1(i,j))
                 pp = fourth*fact*scaleDim*pp

                 ! Compute the corresponding force.

                 fx = pp*ss(i,j,1)
                 fy = pp*ss(i,j,2)
                 fz = pp*ss(i,j,3)

                 ! Distribute the force to the 4 nodes of the quad.
                 ! Note that the averaging factor 1/4 has already been
                 ! taken into account in pp.

                 jj = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + 1
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + iEnd - iBeg + 1
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + 1
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

               enddo
             enddo

             ! For a navier-stokes boundary also the viscous
             ! contribution is taken into account.

             viscTest: if(BCType(mm) == NSWallAdiabatic .or. &
                          BCType(mm) == NSWallIsothermal) then

               ! Compute the viscous force on each of the faces and
               ! scatter it to the 4 nodes, whose forces are updated.

               do j=jBeg,jEnd
                 do i=iBeg,iEnd

                   ! Store the components of the stress tensor
                   ! a bit easier.

                   tauxx = viscSubface(mm)%tau(i,j,1)
                   tauyy = viscSubface(mm)%tau(i,j,2)
                   tauzz = viscSubface(mm)%tau(i,j,3)
                   tauxy = viscSubface(mm)%tau(i,j,4)
                   tauxz = viscSubface(mm)%tau(i,j,5)
                   tauyz = viscSubface(mm)%tau(i,j,6)

                   ! Compute the viscous force on the face. A minus sign
                   ! is now present, due to the definition of this force.
                   ! Furthermore, the conversion to s.I. Units is made
                   ! and the averaging factor of 1/4 is taken into
                   ! account.

                   fx = -fourth*fact*scaleDim*(tauxx*ss(i,j,1) &
                      +                        tauxy*ss(i,j,2) &
                      +                        tauxz*ss(i,j,3))
                   fy = -fourth*fact*scaleDim*(tauxy*ss(i,j,1) &
                      +                        tauyy*ss(i,j,2) &
                      +                        tauyz*ss(i,j,3))
                   fy = -fourth*fact*scaleDim*(tauxz*ss(i,j,1) &
                      +                        tauyz*ss(i,j,2) &
                      +                        tauzz*ss(i,j,3))

                   ! Distribute the force to the 4 nodes of the quad.

                   jj = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                   forceLoc(1,jj) = forceLoc(1,jj) + fx
                   forceLoc(2,jj) = forceLoc(2,jj) + fy
                   forceLoc(3,jj) = forceLoc(3,jj) + fz

                   jj = jj + 1
                   forceLoc(1,jj) = forceLoc(1,jj) + fx
                   forceLoc(2,jj) = forceLoc(2,jj) + fy
                   forceLoc(3,jj) = forceLoc(3,jj) + fz

                   jj = jj + iEnd - iBeg + 1
                   forceLoc(1,jj) = forceLoc(1,jj) + fx
                   forceLoc(2,jj) = forceLoc(2,jj) + fy
                   forceLoc(3,jj) = forceLoc(3,jj) + fz

                   jj = jj + 1
                   forceLoc(1,jj) = forceLoc(1,jj) + fx
                   forceLoc(2,jj) = forceLoc(2,jj) + fy
                   forceLoc(3,jj) = forceLoc(3,jj) + fz

                 enddo
               enddo

             endif viscTest

             ! Update the counter ii.

             ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)

           endif storeSubfaceTest
         enddo bocos
       enddo domains

       ! Test if the memory of mdSurfForce has already been allocated.
       ! If not, allocate it.

       if(.not. allocated(mdSurfForce) ) then

         jj = mdNSurfNodes(nProc,max(cgnsNfamilies,1_intType))
         allocate(mdSurfForce(3,jj), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("mdCreateSurfForceList", &
                          "Memory allocation failure for &
                          &mdSurfForce")
       endif

       ! Construct the arrays recvcounts and displs needed for gatherv.

       do nn=1,nProc
         recvcounts(nn) = 3*(mdNSurfNodes(nn,  modFamID) &
                        -    mdNSurfNodes(nn-1,modFamID))
         displs(nn)     = 3* mdNSurfNodes(nn-1,modFamID)
       enddo

       ! Call allgatherv to gather the data.

       size = 3*nSurfNodesLoc
       call mpi_allgatherv(forceLoc, size, sumb_real, mdSurfForce, &
                           recvcounts, displs, sumb_real,          &
                           SUmb_comm_world, ierr)

       ! Release the memory of forceLoc.

       deallocate(forceLoc, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("mdCreateSurfForceList", &
                        "Deallocation failure for forceLoc")

       ! Set the values of startInd and endInd, which give the range
       ! in the array mdSurfForce where the info for the given family
       ! is stored.

       startInd = mdNSurfNodes(0,    modFamID) + 1
       endInd   = mdNSurfNodes(nProc,modFamID)

       end subroutine mdCreateSurfForceList

!      ==================================================================

       subroutine mdDeleteSurfForceList
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteSurfForceList deallocates the memory of                *
!      * mdSurfForce.                                                   *
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
       ! Deallocate the memory of mdSurfForce if it has been allocated.

       if( allocated(mdSurfForce) ) then

         deallocate(mdSurfForce, stat=ierr)
         if(ierr /= 0)                             &
           call terminate("mdDeleteSurfForceList", &
                          "Deallocation error for mdSurfForce")
       endif

       end subroutine mdDeleteSurfForceList

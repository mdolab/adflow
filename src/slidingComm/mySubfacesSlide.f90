!
!      ******************************************************************
!      *                                                                *
!      * File:          mySubfacesSlide.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-17-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mySubfacesSlide(level, sps, color)
!
!      ******************************************************************
!      *                                                                *
!      * mySubfacesSlide determines the local subfaces on the sliding   *
!      * mesh interface as well as its quadrilateral faces.             *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
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
       integer :: ierr

       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: slideID
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the slide id a bit easier.

       slideID = myInterfaces(color)%globalSlideID

       ! Determine the number of local subfaces for part 1 and part 2
       ! of the currently active sliding mesh interface. Note that for
       ! BCType, etc., the index 1 and not sps must be used, because
       ! this info is only allocated for the 1st spectral solution.
       ! It is the same for all spectral solutions.

       nMySubfaces1 = 0
       nMySubfaces2 = 0

       do nn=1,nDom
         do mm=1,flowDoms(nn,level,1)%nBocos
           if(flowDoms(nn,level,1)%BCType(mm) == slidingInterface) then

             if(flowDoms(nn,level,1)%groupNum(mm) == -slideID) &
               nMySubfaces1 = nMySubfaces1 + 1

             if(flowDoms(nn,level,1)%groupNum(mm) == slideID) &
               nMySubfaces2 = nMySubfaces2 + 1

           endif
         enddo
       enddo

       ! Allocate the memory for mySubfaces1 and mySubfaces2.

       allocate(mySubfaces1(nMySubfaces1), &
                mySubfaces2(nMySubfaces2), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mySubfacesSlide", &
                        "Memory allocation failure for mySubfaces1 &
                        &and mySubfaces2.")

       ! Determine and store the info of the sliding mesh subfaces.

       call storeMySubfaceInfoSlide(level, sps, -slideID, color, &
                                    mySubfaces1)
       call storeMySubfaceInfoSlide(level, sps,  slideID, color, &
                                    mySubfaces2)

       end subroutine mySubfacesSlide

!      ==================================================================

       subroutine storeMySubfaceInfoSlide(level, sps, slideID, &
                                          color, mySubfaces)
!
!      ******************************************************************
!      *                                                                *
!      * storeMySubfaceInfoSlide stores the sliding mesh surface        *
!      * grid info for all local subfaces involved for the given        *
!      * sliding mesh interface. Note that slideID can be both a        *
!      * negative and a positive number.                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, slideID, color
       type(localSubfaceType), dimension(*), intent(out) :: mySubfaces
!
!      Local parameters.
!
       integer(kind=porType), parameter :: izero = 0_porType
       integer(kind=porType), parameter :: ione  = 1_porType
       integer(kind=porType), parameter :: itwo  = 2_porType
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i, j, ii, m1
       integer(kind=intType) :: nQuad, nNode, nDual
       integer(kind=intType) :: n1, n2, n3, n4, slideIDAbs
       integer(kind=intType) :: nSub, nh1, nh2, nd, nnh
       integer(kind=intType) :: ind1, ind2, indc, nnzq

       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd

       integer(kind=intType), dimension(4) :: nk

       integer(kind=porType), dimension(:,:), allocatable :: nodeInfo

       integer(kind=intType), dimension(:), pointer :: donorDualQ

       integer(kind=intType), dimension(:,:), pointer :: haloInfo
       integer(kind=intType), dimension(:,:), pointer :: indHalo1
       integer(kind=intType), dimension(:,:), pointer :: indHalo2
       integer(kind=intType), dimension(:,:), pointer :: indHaloN
       integer(kind=intType), dimension(:,:), pointer :: connQuad
       integer(kind=intType), dimension(:,:), pointer :: connDual

       integer(kind=porType), dimension(:), pointer :: statusQuad
       integer(kind=porType), dimension(:), pointer :: statusDual

       real(kind=realType), dimension(:,:,:), pointer :: xFace, xFaceInt
       real(kind=realType), dimension(:,:),   pointer :: coorN, coorNInt
       real(kind=realType), dimension(:,:),   pointer :: coorQuad

       real(kind=realType), dimension(4) :: thetaN

       logical, dimension(:), pointer :: searchQuad, searchNode
       logical, dimension(:), pointer :: storeQuad,  storeDual
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the absolute value of the slide ID.

       slideIDAbs = abs(slideID)

       ! Initialize the counter nSub to 0 and perform the loop over the
       ! local boundary subfaces.

       nSub = 0

       domains: do nn=1,nDom

         ! Set the pointers for this block.

         call setPointers(nn, level, sps)

         ! Loop over the boundary conditions.

         bocos: do mm=1,nBocos

           ! Test if this is a sliding mesh interface and if so, check
           ! if it belongs to the part currently investigated.

           interfaceTest: if(BCType(mm)   == slidingInterface .and. &
                             groupNum(mm) == slideID) then

             ! Subface must be stored in mySubfaces.
             ! Determine the block face on which the subface is located
             ! and set a number of variables and pointers accordingly.
             ! Note that when subarrays are accessed, the lower boundary
             ! is automatically set to 1. Especially for the coordinates
             ! this should be remembered.

             select case (BCFaceID(mm))

               case (iMin)
                 nh1  = 1; nh2  = 0; nd   = 2; nnh = 0
                 ind1 = 2; ind2 = 3; indc = 1

                 haloInfo => donorDoms(nn)%haloInfo(2,:,:)
                 xFace    => x(1,:,:,:)
                 xFaceInt => x(2,:,:,:)

               case (iMax)
                 nh1  = ie; nh2  = ib; nd   = il; nnh = ie
                 ind1 = 2;  ind2 = 3;  indc = 1

                 haloInfo => donorDoms(nn)%haloInfo(il,:,:)
                 xFace    => x(il,:,:,:)
                 xFaceInt => x(nx,:,:,:)

               case (jMin)
                 nh1  = 1; nh2  = 0; nd   = 2; nnh = 0
                 ind1 = 1; ind2 = 3; indc = 2

                 haloInfo => donorDoms(nn)%haloInfo(:,2,:)
                 xFace    => x(:,1,:,:)
                 xFaceInt => x(:,2,:,:)

               case (jMax)
                 nh1  = je; nh2  = jb; nd   = jl; nnh = je
                 ind1 = 1;  ind2 = 3;  indc = 2

                 haloInfo => donorDoms(nn)%haloInfo(:,jl,:)
                 xFace    => x(:,jl,:,:)
                 xFaceInt => x(:,ny,:,:)

               case (kMin)
                 nh1  = 1; nh2  = 0; nd   = 2; nnh = 0
                 ind1 = 1; ind2 = 2; indc = 3

                 haloInfo => donorDoms(nn)%haloInfo(:,:,2)
                 xFace    => x(:,:,1,:)
                 xFaceInt => x(:,:,2,:)

               case (kMax)
                 nh1  = ke; nh2  = kb; nd   = kl; nnh = ke
                 ind1 = 1;  ind2 = 2;  indc = 3

                 haloInfo => donorDoms(nn)%haloInfo(:,:,kl)
                 xFace    => x(:,:,kl,:)
                 xFaceInt => x(:,:,nz,:)

             end select

             ! Store the cell range a bit easier. As the cell range of
             ! BCData may or may not contain the halo cells, it is safer
             ! to use the nodal range here. Note that the upper index
             ! needs to be adapted to obtain the cell range including
             ! the halo's.

             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd + 1
             jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd + 1

             ! Update the counter nSub, store the block and face ID and
             ! initialize zeroCrossed and piCrossed.

             nSub = nSub + 1

             mySubfaces(nSub)%blockID = nn
             mySubfaces(nSub)%faceID  = BCFaceID(mm)

             mySubfaces(nSub)%zeroCrossedQuad = .false.
             mySubfaces(nSub)%piCrossedQuad   = .false.

             mySubfaces(nSub)%zeroCrossedDual = .false.
             mySubfaces(nSub)%piCrossedDual   = .false.

             ! Determine the number of surface nodes, quadrilaterals and 
             ! dual quadrilaterals of the subface.

             m1 = (iEnd - iBeg + 1)
             nDual = (jEnd - jBeg)*(m1-1)
             nQuad = (jEnd - jBeg + 1)*m1
             nNode = (jEnd - jBeg + 2)*(m1+1)

             mySubfaces(nSub)%nQuad = nQuad
             mySubfaces(nSub)%nNode = nNode
             mySubfaces(nSub)%nDual = nDual

             ! Allocate the memory for the arrays of mySubfaces(nSub)
             ! which are used in this routine and some help arrays.

             allocate(mySubfaces(nSub)%searchQuad(nQuad),   &
                      mySubfaces(nSub)%searchNode(nNode),   &
                      mySubfaces(nSub)%statusQuad(nQuad),   &
                      mySubfaces(nSub)%statusDual(nDual),   &
                      mySubfaces(nSub)%indHalo1(nQuad,3),   &
                      mySubfaces(nSub)%indHalo2(nQuad,3),   &
                      mySubfaces(nSub)%indHaloN(nNode,3),   &
                      mySubfaces(nSub)%donorDualQ(nQuad),   &
                      mySubfaces(nSub)%connQuad(nQuad,4),   &
                      mySubfaces(nSub)%connDual(nDual,4),   &
                      mySubfaces(nSub)%storeQuad(nQuad),    &
                      mySubfaces(nSub)%storeDual(nDual),    &
                      mySubfaces(nSub)%coorQuad(nQuad,3),   &
                      mySubfaces(nSub)%coorN(nNode,3),      &
                      mySubfaces(nSub)%coorNInt(nNode,3),   &
                      nodeInfo(iBeg-1:iEnd,jBeg-1:jEnd),    &
                      stat = ierr)
             if(ierr /= 0)                                &
               call terminate("storeMySubfaceInfoSlide",  &
                              "Memory allocation failure for &
                              &mySubfaces and nodeInfo.")

             ! Set a couple of pointers to make the code more readable.

             searchQuad => mySubfaces(nSub)%searchQuad
             searchNode => mySubfaces(nSub)%searchNode
             statusQuad => mySubfaces(nSub)%statusQuad
             statusDual => mySubfaces(nSub)%statusDual
             indHalo1   => mySubfaces(nSub)%indHalo1
             indHalo2   => mySubfaces(nSub)%indHalo2
             indHaloN   => mySubfaces(nSub)%indHaloN
             donorDualQ => mySubfaces(nSub)%donorDualQ
             connQuad   => mySubfaces(nSub)%connQuad
             connDual   => mySubfaces(nSub)%connDual
             storeQuad  => mySubfaces(nSub)%storeQuad
             storeDual  => mySubfaces(nSub)%storeDual
             coorQuad   => mySubfaces(nSub)%coorQuad
             coorN      => mySubfaces(nSub)%coorN
             coorNInt   => mySubfaces(nSub)%coorNInt
!
!            ************************************************************
!            *                                                          *
!            * Determine the nodal coordinates.                         *
!            *                                                          *
!            ************************************************************
!
             indHaloN => mySubfaces(nSub)%indHaloN
             coorN    => mySubfaces(nSub)%coorN
             coorNInt => mySubfaces(nSub)%coorNInt
             call getNodalCylinderCoor

             ! Initialize nodeInfo to 0, which means that no cell
             ! to which this node belongs has been treated yet.

             nodeInfo = izero
!
!            ************************************************************
!            *                                                          *
!            * Determine all the face info, i.e. whether or not the     *
!            * face must be searched, the cylindrical coordinates of    *
!            * the face center and the connectivity. The quadrilateral  *
!            * surface grid is used to search the coordinates.          *
!            *                                                          *
!            ************************************************************
!
             nQuad = 0
             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 ! Update the local counter nQuad.

                 nQuad = nQuad + 1

                 ! Set the indices for the halo cells. If this is
                 ! not a direct halo, only the 1st level halo is
                 ! is needed even for second order schemes. This is
                 ! indicated by a -1 for the second halo index.

                 indHalo1(nQuad,ind1) = i
                 indHalo1(nQuad,ind2) = j
                 indHalo1(nQuad,indc) = nh1

                 indHalo2(nQuad,ind1) = i
                 indHalo2(nQuad,ind2) = j
                 indHalo2(nQuad,indc) = nh2

                 if(haloInfo(i,j) /= internalCell) indHalo2(nQuad,1) = -1

                 ! Determine whether or not the face must be
                 ! searched in the other part of the interface.

                 if(haloInfo(i,j) == internalCell .or. &
                    haloInfo(i,j) == ownedDonor   .or. &
                    haloInfo(i,j) == unownedDonor .or. &
                    haloInfo(i,j) == boundaryHalo .or. &
                    haloInfo(i,j) <  slideIDAbs) then

                   ! Internal cell, true donor halo, boundary halo or
                   ! a halo that belongs to a sliding mesh interface
                   ! with a lower ID. Search for a donor.

                   searchQuad(nQuad) = .true.

                   ! Set the nodeInfo for the 4 nodes such that
                   ! the nodes must be interpolated.

                   nodeInfo(i-1,j-1) = itwo
                   nodeInfo(i,  j-1) = itwo
                   nodeInfo(i-1,j)   = itwo
                   nodeInfo(i,  j)   = itwo

                 else

                   ! Halo belongs to a sliding mesh interface with a
                   ! higher ID. I am not allowed to construct the halo
                   ! here. Set donorDualQ to 0 just to be safe.

                   searchQuad(nQuad) = .false.
                   donorDualQ(nQuad) = 0

                   ! Set the nodeInfo for the 4 nodes, such that
                   ! interpolation info is not overwritten.

                   nodeInfo(i-1,j-1) = max(nodeInfo(i-1,j-1), ione)
                   nodeInfo(i,  j-1) = max(nodeInfo(i,  j-1), ione)
                   nodeInfo(i-1,j)   = max(nodeInfo(i-1,j),   ione)
                   nodeInfo(i,  j)   = max(nodeInfo(i,  j),   ione)

                 endif

                 ! Set the connectivity for the four nodes of this
                 ! quadrilateral.

                 n1 = (j-jBeg)*(m1+1) + i - iBeg + 1
                 n2 = n1 + 1
                 n4 = n2 + m1
                 n3 = n4 + 1

                 connQuad(nQuad,1) = n1
                 connQuad(nQuad,2) = n2
                 connQuad(nQuad,3) = n3
                 connQuad(nQuad,4) = n4

                 ! Determine whether or not the quad should be stored in
                 ! the surface mesh.

                 if(haloInfo(i,j) == unownedDonor) then
                   storeQuad(nQuad) = .false.
                 else
                   storeQuad(nQuad) = .true.
                 endif

                 ! Determine the cylindrical coordinates of the face
                 ! center on the interface.

                 call getFaceCenterCylinderCoor

               enddo
             enddo
!
!            ************************************************************
!            *                                                          *
!            * Loop again over the nodes to determine whether or not    *
!            * the node is to be searched.                              *
!            *                                                          *
!            ************************************************************
!
             nNode = 0
             do j=(jBeg-1),jEnd
               do i=(iBeg-1),iEnd

                 ! Update the counter nNode.

                 nNode = nNode + 1

                 ! Determine the case we are having here.

                 select case (nodeInfo(i,j))

                   case (itwo)

                     ! Node must be interpolated.

                     searchNode(nNode) = .true.

                   case (ione)

                     ! Sliding mesh node for which I am not
                     ! responsible. Node must not be interpolated.

                     searchNode(nNode) = .false.

                   case default

                     call terminate("storeMySubfaceInfoSlide", &
                                    "This should not happen")

                 end select

               enddo
             enddo
!
!            ************************************************************
!            *                                                          *
!            * Determine the connectivity of the dual mesh, which will  *
!            * be used for the interpolation of the flow variables.     *
!            *                                                          *
!            ************************************************************
!
             nDual = 0
             do j=(jBeg+1),jEnd
               do i=(iBeg+1),iEnd

                 ! Update nDual and determine the 4 indices of the cell
                 ! centers, which define the dual grid quadrilateral.

                 nDual = nDual + 1

                 n1 = (j-jBeg-1)*m1 + i - iBeg
                 n2 = n1 + 1
                 n3 = n2 + m1
                 n4 = n3 - 1

                 ! Store the connectivity.

                 connDual(nDual,1) = n1
                 connDual(nDual,2) = n2
                 connDual(nDual,3) = n3
                 connDual(nDual,4) = n4

                 ! Determine whether or not the dual quad should be
                 ! stored in the surface mesh. If only one halo type
                 ! of the 4 cell centers equals an unowned donor this
                 ! means that the dual quad does not have to be stored.

                 if(haloInfo(i-1,j-1) == unownedDonor .or. &
                    haloInfo(i  ,j-1) == unownedDonor .or. &
                    haloInfo(i-1,j)   == unownedDonor .or. &
                    haloInfo(i,  j)   == unownedDonor) then
                   storeDual(nDual) = .false.
                 else
                   storeDual(nDual) = .true.
                 endif

                 ! Determine the number of nodes in the 4 quadrants
                 ! of the polar plane.

                 thetaN(1) = coorQuad(n1,3)
                 thetaN(2) = coorQuad(n2,3)
                 thetaN(3) = coorQuad(n3,3)
                 thetaN(4) = coorQuad(n4,3)

                 nk(1) = 0; nk(2) = 0; nk(3) = 0; nk(4) = 0

                 do ii=1,4
                   if(thetaN(ii) > zero) then
                     if(thetaN(ii) > half*pi) then
                       nk(2) = nk(2) + 1
                     else
                       nk(1) = nk(1) + 1
                     endif
                   else
                     if(thetaN(ii) > -half*pi) then
                       nk(4) = nk(4) + 1
                     else
                       nk(3) = nk(3) + 1
                     endif
                   endif
                 enddo

                 ! Determine the number of nonzero quadrants.

                 nnzq = 0
                 if(nk(1) > 0) nnzq = nnzq + 1
                 if(nk(2) > 0) nnzq = nnzq + 1
                 if(nk(3) > 0) nnzq = nnzq + 1
                 if(nk(4) > 0) nnzq = nnzq + 1

                 ! Determine the status of the dual face.

                 statusDual(nDual) = normalQuad

                 select case (nnzq)

                   case (4_intType, 3_intType)

                     ! Rotation center is inside the quad.
                     ! Set the status to polar singular.

                     statusDual(nDual) = polarSingular

                   case (2_intType)

                     ! Dual cell has nodes in two quadrants.
                     ! Investigate the situation a bit further.

                     if(nk(2) > 0 .and. nk(3) > 0) then

                       ! The cell crosses the line theta == pi. Set the
                       ! status accordingly and set thetaPiCrossed to
                       ! .true. for this subface.

                       statusDual(nDual)              = piCrossed
                       mySubfaces(nSub)%piCrossedDual = .true.

                     else if(nk(1) > 0 .and. nk(4) > 0) then

                       ! The cell crosses the line theta == zero.
                       ! Nothing special for the quad, but set
                       ! zeroCrossedDual to .true.

                       mySubfaces(nSub)%zeroCrossedDual = .true.

                     else if((nk(1) > 0 .and. nk(3) > 0) .or. &
                             (nk(2) > 0 .and. nk(4) > 0)) then

                       ! Cell has nodes in two opposite quadrants.
                       ! This can only mean that the rotation center is
                       ! inside the quad. Set the status to polar
                       ! singular.

                       statusDual(nDual) = polarSingular

                     endif

                 end select

               enddo
             enddo

             ! Deallocate the memory of nodeInfo.

             deallocate(nodeInfo, stat=ierr)
             if(ierr /= 0)                             &
               call terminate("storeSubfaceInfoSlide", &
                              "Deallocation error for nodeInfo.")
           endif interfaceTest

         enddo bocos

       enddo domains

       !=================================================================

       contains

         !===============================================================

         subroutine getNodalCylinderCoor
!
!        ****************************************************************
!        *                                                              *
!        * getNodalCylinderCoor determines for the nodes of the         *
!        * currently active block face the cylinder coordinates         *
!        * relative to the rotation center of the sliding interface.    *
!        *                                                              *
!        ****************************************************************
!
         use interfaceGroups
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: ll

         real(kind=realType) :: xx, yy, zz, ax, r1, r2, r, theta
         real(kind=realType) :: scaleFact, axMin, rMin

         real(kind=realType), dimension(3) :: rotCenter, rotAxis
         real(kind=realType), dimension(3) :: radVec1, radVec2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Some abbreviations to make the code more readable.

         scaleFact = myInterfaces(color)%scale
         axMin     = myInterfaces(color)%axMin
         rMin      = myInterfaces(color)%rMin

         rotCenter = myInterfaces(color)%rotCenter
         rotAxis   = myInterfaces(color)%rotAxis
         radVec1   = myInterfaces(color)%radVec1
         radVec2   = myInterfaces(color)%radVec2

         ! Loop over the nodes of the subface.

         nNode = 0
         do j=(jBeg-1),jEnd
           do i=(iBeg-1),iEnd

             ! Update the counter, determine the corresponding halo
             ! indices and the cylindrical coordinates for the node.

             nNode = nNode + 1
             indHaloN(nNode,ind1) = i
             indHaloN(nNode,ind2) = j
             indHaloN(nNode,indc) = nnh

             ! Determine the Cartesian coordinates relative to the center
             ! of rotation. As xFace is a pointer to a slice in x, the
             ! lower indices of xFace are automatically set to 1. This in
             ! contrast to the original array x, which start at 0.
             ! Therefore i+1 and j+1 must be used.

             xx = xFace(i+1,j+1,1) - rotCenter(1)
             yy = xFace(i+1,j+1,2) - rotCenter(2)
             zz = xFace(i+1,j+1,3) - rotCenter(3)

             ! Determine the axial and two radial components for this point.

             ax = xx*rotAxis(1) + yy*rotAxis(2) + zz*rotAxis(3)
             r1 = xx*radVec1(1) + yy*radVec1(2) + zz*radVec1(3)
             r2 = xx*radVec2(1) + yy*radVec2(2) + zz*radVec2(3)

             ! Compute the polar coordinates r and theta from r1 and r2.
             ! Be careful when a polar singularity is present. In that
             ! case the angle of the nearest nonsingular neighbor is
             ! taken.

             r = sqrt(r1*r1 + r2*r2)

             if((abs(r1) < eps) .and. (abs(r2) < eps)) then

               ! Polar singularity. Find the nearest neighbor away
               ! from the singularity.

               ! First try the i-direction.
               ! Determine whether the node i+1 or i-1 must be taken.
               ! In the ll-index is already taken into account that it
               ! must be incremented by 1, because of this pointer stuff.

               ll = i
               if(i <= BCData(mm)%inBeg) ll = i+2

               ! Compute the corresponding values of r1 and r2.

               xx = xFace(ll,j+1,1) - rotCenter(1)
               yy = xFace(ll,j+1,2) - rotCenter(2)
               zz = xFace(ll,j+1,3) - rotCenter(3)

               r1 = xx*radVec1(1) + yy*radVec1(2) + zz*radVec1(3)
               r2 = xx*radVec2(1) + yy*radVec2(2) + zz*radVec2(3)

               ! Check if r1 and r2 still define a polar singularity.
               ! If so, it must be the j-direction.

               if((abs(r1) < eps) .and. (abs(r2) < eps)) then

                 ! Determine whether the node j+1 or j-1 must be taken.
                 ! Again the increment is already taken into account
                 ! in ll.

                 ll = j
                 if(j <= BCData(mm)%jnBeg) ll = j+2

                 ! Compute the corresponding values of r1 and r2.

                 xx = xFace(i+1,ll,1) - rotCenter(1)
                 yy = xFace(i+1,ll,2) - rotCenter(2)
                 zz = xFace(i+1,ll,3) - rotCenter(3)

                 r1 = xx*radVec1(1) + yy*radVec1(2) + zz*radVec1(3)
                 r2 = xx*radVec2(1) + yy*radVec2(2) + zz*radVec2(3)

               endif
             endif

             ! Compute the angle based on r1 and r2.

             theta = atan2(r2,r1)

             ! Store and scale the coordinates.

             coorN(nNode,1) = (ax - axMin)*scaleFact
             coorN(nNode,2) = (r  - rMin)* scaleFact
             coorN(nNode,3) = theta

             ! Compute the corresponding coordinates 1 layer into the
             ! block. First store the cartesian components relative to
             ! the center of rotation. Again the addition of 1 because
             ! of the pointers.

             xx = xFaceInt(i+1,j+1,1) - rotCenter(1)
             yy = xFaceInt(i+1,j+1,2) - rotCenter(2)
             zz = xFaceInt(i+1,j+1,3) - rotCenter(3)

             ! Determine the axial and two radial components for this point.

             ax = xx*rotAxis(1) + yy*rotAxis(2) + zz*rotAxis(3)
             r1 = xx*radVec1(1) + yy*radVec1(2) + zz*radVec1(3)
             r2 = xx*radVec2(1) + yy*radVec2(2) + zz*radVec2(3)

             ! Compute the polar coordinates r and theta from r1 and r2.
             ! Be careful when a polar singularity is present. In that
             ! case the angle of the coordinate on the face is taken.

             r = sqrt(r1*r1 + r2*r2)

             if((abs(r1) < eps) .and. (abs(r2) < eps)) then
               theta = coorN(nNode,3)
             else
               theta = atan2(r2,r1)
             endif

             ! Store and scale the axial and radial coordinates.

             coorNInt(nNode,1) =  (ax - axMin)*scaleFact
             coorNInt(nNode,2) =  (r  - rMin)* scaleFact

             ! Make sure that if the grid line theta == pi is crossed
             ! between the face and its internal counterpart to add or
             ! substract 2*pi to/from the angle of the internal point.

             ! Test for adding 2*pi, i.e. the node on the face is in the
             ! second quadrant and the internal node in the third.

             if(coorN(nNode,3) > half*pi .and. theta <= -half*pi) &
               theta = theta + two*pi

             ! Test for substracting 2*pi, i.e. the node on the face is
             ! in the third quadrant and the internal node in the second.

             if(coorN(nNode,3) <= -half*pi .and. theta > half*pi) &
               theta = theta - two*pi

             ! Store the angle.

             coorNInt(nNode,3) = theta

           enddo
         enddo

         end subroutine getNodalCylinderCoor

         !===============================================================

         subroutine getFaceCenterCylinderCoor
!
!        ****************************************************************
!        *                                                              *
!        * getFaceCenterCylinderCoor determines the cylindrical         *
!        * coordinates of the currently active face center and          *
!        * determines the status of the quad.                           *
!        *                                                              *
!        ****************************************************************
!
         use interfaceGroups
         implicit none
!
!        Local variables.
!
         integer(kind=intType), dimension(4) :: nk

         real(kind=realType) :: aAvg, rAvg, tAvg
         real(kind=realType) :: xx, yy, zz, ax, r1, r2, r, theta
         real(kind=realType) :: scaleFact, axMin, rMin

         real(kind=realType), dimension(3) :: rotCenter, rotAxis
         real(kind=realType), dimension(3) :: radVec1, radVec2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Some abbreviations to make the code more readable.

         scaleFact = myInterfaces(color)%scale
         axMin     = myInterfaces(color)%axMin
         rMin      = myInterfaces(color)%rMin

         rotCenter = myInterfaces(color)%rotCenter
         rotAxis   = myInterfaces(color)%rotAxis
         radVec1   = myInterfaces(color)%radVec1
         radVec2   = myInterfaces(color)%radVec2
 
         ! Initialize the status to a normal quad.

         statusQuad(nQuad)  = normalQuad

         ! Average the nodal coordinates in the local cylindrical
         ! coordinate system.

         aAvg = fourth*(coorN(n1,1) + coorN(n2,1) &
              +         coorN(n3,1) + coorN(n4,1))
         rAvg = fourth*(coorN(n1,2) + coorN(n2,2) &
              +         coorN(n3,2) + coorN(n4,2))
         tAvg = fourth*(coorN(n1,3) + coorN(n2,3) &
              +         coorN(n3,3) + coorN(n4,3))

         ! Determine the number of nodes in the 4 quadrants of the
         ! polar plane.

         thetaN(1) = coorN(n1,3)
         thetaN(2) = coorN(n2,3)
         thetaN(3) = coorN(n3,3)
         thetaN(4) = coorN(n4,3)

         nk(1) = 0; nk(2) = 0; nk(3) = 0; nk(4) = 0

         do ii=1,4
           if(thetaN(ii) > zero) then
             if(thetaN(ii) > half*pi) then
               nk(2) = nk(2) + 1
             else
               nk(1) = nk(1) + 1
             endif
           else
             if(thetaN(ii) > -half*pi) then
               nk(4) = nk(4) + 1
             else
               nk(3) = nk(3) + 1
             endif
           endif
         enddo

         ! Determine the number of nonzero quadrants.

         nnzq = 0
         if(nk(1) > 0) nnzq = nnzq + 1
         if(nk(2) > 0) nnzq = nnzq + 1
         if(nk(3) > 0) nnzq = nnzq + 1
         if(nk(4) > 0) nnzq = nnzq + 1

         ! Determine the situation we are having here.

         select case (nnzq)

           case (4_intType, 3_intType)

             ! Rotation center is inside the quad. Set the status to
             ! polar singular such that cartesian averaging is used to
             ! determine the cell center coordinates.

             statusQuad(nQuad) = polarSingular

           case (2_intType)

             ! Cell has nodes in two quadrants.
             ! Investigate the situation a bit further.

             if(nk(2) > 0 .and. nk(3) > 0) then

               ! The cell crosses the line theta = pi. Set the status
               ! accordingly and update the average angle first such that
               ! all angles correspond to a positive angle. Afterwards,
               ! restrict it to the domain -pi .. pi.

               statusQuad(nQuad)              = piCrossed
               mySubfaces(nSub)%piCrossedQuad = .true.

               tAvg = tAvg + nk(3)*half*pi
               if(tAvg > pi) tAvg = tAvg - two*pi

             else if(nk(1) > 0 .and. nk(4) > 0) then

               ! The cell crosses the line theta == zero. Nothing special
               ! for the quad, but set zeroCrossedQuad to .true.

               mySubfaces(nSub)%zeroCrossedQuad = .true.

             else if((nk(1) > 0 .and. nk(3) > 0) .or. &
                     (nk(2) > 0 .and. nk(4) > 0)) then

               ! Cell has nodes in two opposite quadrants. This can only
               ! mean that the rotation center is inside the quad.
               ! Set the status to polar singular such that cartesian
               ! averaging is used to determine the cell center coordinates.

               statusQuad(nQuad) = polarSingular

             endif

         end select

         ! Store the cylindrical coordinates of the face center.

         coorQuad(nQuad,1) = aAvg
         coorQuad(nQuad,2) = rAvg
         coorQuad(nQuad,3) = tAvg

         ! If the cell contains a polar singularity use Cartesian
         ! averaging to obtain the coordinates of the quadrilateral.

         if(statusQuad(nQuad) == polarSingular) then

           ! Compute the average of the cartesian coordinates. Again
           ! note the offset in the indices because of the usage of the
           ! pointer xFace.

           xx = fourth*(xFace(i,j,  1) + xFace(i+1,j,  1)  &
              +         xFace(i,j+1,1) + xFace(i+1,j+1,1)) &
              - rotCenter(1)
           yy = fourth*(xFace(i,j,  2) + xFace(i+1,j,  2)  &
              +         xFace(i,j+1,2) + xFace(i+1,j+1,2)) &
              - rotCenter(2)
           zz = fourth*(xFace(i,j,  3) + xFace(i+1,j,  3)  &
              +         xFace(i,j+1,3) + xFace(i+1,j+1,3)) &
              - rotCenter(3)

           ! Determine the axial and two radial components.

           ax = xx*rotAxis(1) + yy*rotAxis(2) + zz*rotAxis(3)
           r1 = xx*radVec1(1) + yy*radVec1(2) + zz*radVec1(3)
           r2 = xx*radVec2(1) + yy*radVec2(2) + zz*radVec2(3)

           ! Compute the polar coordinates r and theta from r1 and r2.
           ! Be careful when a polar singularity is present.

           r = sqrt(r1*r1 + r2*r2)
           if((abs(r1) < eps) .and. (abs(r2) < eps)) then
             theta = zero
           else
             theta = atan2(r2,r1)
           endif

           ! Store the coordinates in the local cylindrical frame.
           ! Apply the scaling.

           coorQuad(nQuad,1) = (ax - axMin)*scaleFact
           coorQuad(nQuad,2) = (r  - rMin)* scaleFact
           coorQuad(nQuad,3) = theta

         endif

         end subroutine getFaceCenterCylinderCoor

       end subroutine storeMySubfaceInfoSlide

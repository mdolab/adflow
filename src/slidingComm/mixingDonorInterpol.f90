!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingDonorInterpol.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-02-2005                                      *
!      * Last modified: 09-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mixingDonorInterpol(level, slideID, color, commPattern)
!
!      ******************************************************************
!      *                                                                *
!      * mixingDonorInterpol determines the donor interpolation data    *
!      * for the given slideID on the given multigrid level.            *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commMixing
       use interfaceGroups
       use localSubfacesMod
       use mixingData
       implicit none
!
!      Parameter definitions
! 
       integer(kind=intType), parameter :: below         = 1_intType
       integer(kind=intType), parameter :: contained     = 2_intType
       integer(kind=intType), parameter :: above         = 3_intType
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, slideID, color

       type(commMixingType), intent(inout) :: commPattern
!
!      Local variables.
!
       integer :: size, ierr

       integer(kind=intType) :: nn, mm, kk, ll, ii, jj, i, j, k
       integer(kind=intType) :: nInter, nDonor, nnzq, nSubSlide
       integer(kind=intType) :: iBeg,  iEnd,  jBeg,  jEnd
       integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
       integer(kind=intType) :: d1, d2, dind, din1, din2
       integer(kind=intType) :: start, end, interval, nNodesCV
       integer(kind=intType) :: newSize, nAlloc, n1, n2, n3

       integer(kind=intType), dimension(4) :: nk
       integer(kind=intType), dimension(4) :: statusNodes

       integer(kind=intType), dimension(:,:,:), allocatable :: subRange

       integer(kind=intType), dimension(:), pointer :: blockDonor
       integer(kind=intType), dimension(:), pointer :: nIntervalsDonor
       integer(kind=intType), dimension(:), pointer :: indListDonor

       integer(kind=intType), dimension(:,:),   pointer :: haloInfo
       integer(kind=intType), dimension(:,:,:), pointer :: indD

       integer(kind=porType) :: statusQuad

       real(kind=realType) :: thetaAvg, rMin, rMax, xMin, xMax, area
       real(kind=realType) :: cosTheta, sinTheta, r1, r2, ww1, ww2
       real(kind=realType) :: localWeight, xInter

       real(kind=realType), dimension(3)   :: axis, vec1, vec2
       real(kind=realType), dimension(4)   :: theta, rad, trueRad
       real(kind=realType), dimension(6,3) :: xCV

       real(kind=realType), dimension(commPattern%nInter) :: &
                                                    localSum, globalSum

       real(kind=realType), dimension(:,:,:), allocatable :: xFace

       real(kind=realType), dimension(:),     pointer :: weightDonor
       real(kind=realType), dimension(:,:,:), pointer :: rotMatDonor
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the vectors which define the local Cartesian
       ! coordinate system.

       axis = myInterfaces(color)%rotAxis
       vec1 = myInterfaces(color)%radVec1
       vec2 = myInterfaces(color)%radVec2

       ! Easier storage of the number of interpolation cells.

       nInter = commPattern%nInter
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local number of donors. Halo's are not included, *
!      * except if they correspond to a physical boundary or an         *
!      * intersecting sliding mesh interface.                           *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of contributing subfaces.

       nSubSlide = 0
       do nn=1,nDom
         do mm=1,flowDoms(nn,level,1)%nBocos
           if(flowDoms(nn,level,1)%BCType(mm)   == slidingInterface .and. &
              flowDoms(nn,level,1)%groupNum(mm) == slideID)               &
             nSubSlide = nSubSlide + 1
         enddo
       enddo

       ! Allocate the memory for the face ranges of these subfaces.

       allocate(subRange(nSubSlide,2,2), stat=ierr)
       if(ierr /= 0) &
         call terminate("mixingDonorInterpol", &
                        "Memory allocation failure for subRange.")

       ! Repeat the loop over the subfaces, but now determine the
       ! number of contributing faces. The sum over all subfaces is
       ! the number of local donors.

       nSubSlide = 0
       nDonor = 0

       domainLoop1: do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         bocoLoop1: do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           testInterface1: if(BCType(mm)   == slidingInterface .and. &
                              groupNum(mm) == slideID) then

             ! Determine the cell range, for the moment without halos.
             ! As the cell range in BCData may or may not contain the
             ! halo cells, it is easier to use the nodal range of
             ! BCData. This does not contain the halo. The +1 is
             ! present, because the cell numbering starts 1 higher
             ! than the nodes.

             jBeg = BCData(mm)%jnBeg + 1
             jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1
             iEnd = BCData(mm)%inEnd

             ! Set the pointer for haloInfo, depending on the block
             ! face on which the subface is located.

             select case (BCFaceID(mm))
               case (iMin)
                 haloInfo => donorDoms(nn)%haloInfo(2,:,:)
               case (iMax)
                 haloInfo => donorDoms(nn)%haloInfo(il,:,:)
               case (jMin)
                 haloInfo => donorDoms(nn)%haloInfo(:,2,:)
               case (jMax)
                 haloInfo => donorDoms(nn)%haloInfo(:,jl,:)
               case (kMin)
                 haloInfo => donorDoms(nn)%haloInfo(:,:,2)
               case (kMax)
                 haloInfo => donorDoms(nn)%haloInfo(:,:,kl)
             end select

             ! Check if the start and end indices of the cell range
             ! should be adapted because the halos either correspond to
             ! a physical boundary or a intersecting sliding interface.

             jjBeg = jBeg; jjEnd = jEnd; iiBeg = iBeg; iiEnd = iEnd

             ! Ibeg boundary.

             if(haloInfo(iiBeg-1,jjBeg) > 0 .or. &
                haloInfo(iiBeg-1,jjBeg) == boundaryHalo) then

               ! Check if all the halo's are of the same type.

               do j=jjBeg,jjEnd
                 if(haloInfo(iiBeg-1,j) /= haloInfo(iiBeg-1,jjBeg)) &
                   call terminate("mixingDonorInterpol",            &
                                  "Inconsistent halo info for iBeg &
                                  &boundary.")
               enddo
 
               ! Decrease iBeg.

               iBeg = iBeg - 1
             endif

             ! Iend boundary.

             if(haloInfo(iiEnd+1,jjBeg) > 0 .or. &
                haloInfo(iiEnd+1,jjBeg) == boundaryHalo) then

               ! Check if all the halo's are of the same type.

               do j=jjBeg,jjEnd
                 if(haloInfo(iiEnd+1,j) /= haloInfo(iiEnd+1,jjBeg)) &
                   call terminate("mixingDonorInterpol",            &
                                  "Inconsistent halo info for iEnd &
                                  &boundary.")
               enddo
 
               ! Increase iEnd.

               iEnd = iEnd + 1
             endif

             ! Jbeg boundary.

             if(haloInfo(iiBeg,jjBeg-1) > 0 .or. &
                haloInfo(iiBeg,jjBeg-1) == boundaryHalo) then

               ! Check if all the halo's are of the same type.

               do i=iiBeg,iiEnd
                 if(haloInfo(i,jjBeg-1) /= haloInfo(iiBeg,jjBeg-1)) &
                   call terminate("mixingDonorInterpol",            &
                                  "Inconsistent halo info for jBeg &
                                  &boundary.")
               enddo
 
               ! Decrease jBeg.

               jBeg = jBeg - 1
             endif

             ! Jend boundary.

             if(haloInfo(iiBeg,jjEnd+1) > 0 .or. &
                haloInfo(iiBeg,jjEnd+1) == boundaryHalo) then

               ! Check if all the halo's are of the same type.

               do i=iiBeg,iiEnd
                 if(haloInfo(i,jjEnd+1) /= haloInfo(iiBeg,jjEnd+1)) &
                   call terminate("mixingDonorInterpol",            &
                                  "Inconsistent halo info for jEnd &
                                  &boundary.")
               enddo
 
               ! Increase jEnd.

               jEnd = jEnd + 1
             endif

             ! Store the cell range in subRange and update the
             ! number of donors.

             nSubSlide = nSubSlide + 1

             subRange(nSubSlide,1,1) = iBeg
             subRange(nSubSlide,1,2) = iEnd
             subRange(nSubSlide,2,1) = jBeg
             subRange(nSubSlide,2,2) = jEnd

             nDonor = nDonor + (jEnd - jBeg + 1)*(iEnd - iBeg + 1)

           endif testInterface1
         enddo bocoLoop1
       enddo domainLoop1
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local interpolation information for the donors.  *
!      *                                                                *
!      ******************************************************************
!
       ! Store the number of donors and allocate some memory. Guess the
       ! total number of contributions and allocate the corresponding
       ! arrays accordingly.

       commPattern%nDonor = nDonor
       nAlloc             = 2*nDonor

       allocate(commPattern%blockDonor(nDonor),            &
                commPattern%indD(nDonor,3,2),              &
                commPattern%rotMatDonor(nDonor,3,3),       &
                commPattern%nIntervalsDonor(0:nDonor),     &
                commPattern%indListDonor(nAlloc),          &
                commPattern%weightDonor(nAlloc), stat=ierr)
       if(ierr /= 0) &
         call terminate("mixingDonorInterpol", &
                        "Memory allocation failure for donor data.")

       ! Set some pointers to make the code more readable.

       blockDonor      => commPattern%blockDonor
       indD            => commPattern%indD
       rotMatDonor     => commPattern%rotMatDonor
       nIntervalsDonor => commPattern%nIntervalsDonor
       indListDonor    => commPattern%indListDonor
       weightDonor     => commPattern%weightDonor

       ! Loop again over the local subfaces of this sliding interface;
       ! now the actual information is determined.

       nSubSlide = 0
       ii = 0
       nIntervalsDonor(0) = 0

       domainLoop2: do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         bocoLoop2: do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           testInterface2: if(BCType(mm)   == slidingInterface .and. &
                              groupNum(mm) == slideID) then

             ! Store the range for this subface a bit easier.

             nSubSlide = nSubSlide + 1

             iBeg = subRange(nSubSlide,1,1)
             iEnd = subRange(nSubSlide,1,2)
             jBeg = subRange(nSubSlide,2,1)
             jEnd = subRange(nSubSlide,2,2)

             ! Allocate the memory for the nodal cylindrical coordinates
             ! of the subface and determine them.

             allocate(xFace(iBeg-1:iEnd,jBeg-1:jEnd,3), stat=ierr)
             if(ierr /= 0) &
               call terminate("mixingDonorInterpol", &
                              "Memory allocation failure for xFace.")

             call cylCoorNodesOnBlockFace(xFace,  iBeg-1, iEnd,  &
                                          jBeg-1, jEnd,   color, &
                                          BCFaceID(mm))

             ! If this is a radial interface, store the radii in the
             ! first position such that no test is needed inside the
             ! loop. The axial coordinates are not needed in this case.

             if( radialInterface ) then
               do j=(jBeg-1),jEnd
                 do i=(iBeg-1),iEnd
                   xFace(i,j,1) = xFace(i,j,2)
                 enddo
               enddo
             endif

             ! Set some variables depending on the block face on which
             ! the subface is located.

             select case (BCFaceID(mm))
               case (iMin)
                 d1 =  2; d2 =  3; dind = 1; din1 = 2; din2 = 3
               case (iMax)
                 d1 = il; d2 = nx; dind = 1; din1 = 2; din2 = 3
               case (jMin)
                 d1 =  2; d2 =  3; dind = 2; din1 = 1; din2 = 3
               case (jMax)
                 d1 = jl; d2 = ny; dind = 2; din1 = 1; din2 = 3
               case (kMin)
                 d1 =  2; d2 =  3; dind = 3; din1 = 1; din2 = 2
               case (kMax)
                 d1 = kl; d2 = nz; dind = 3; din1 = 1; din2 = 2
             end select

             ! Loop over the faces of the subface and update the donor
             ! interpolation information.

             jloop: do j=jBeg,jEnd
               iloop: do i=iBeg,iEnd
!
!                ********************************************************
!                *                                                      *
!                * Determine the information of the donor itself, i.e.  *
!                * Indices and the rotation matrix.                     *
!                *                                                      *
!                ********************************************************
!
                 ! Store the donor information.

                 ii = ii + 1
                 blockDonor(ii) = nn

                 indD(ii,dind,1) = d1
                 indD(ii,din1,1) =  i
                 indD(ii,din2,1) =  j

                 indD(ii,dind,2) = d2
                 indD(ii,din1,2) =  i
                 indD(ii,din2,2) =  j

                 ! Store the angles and the radii of the 4 surrounding
                 ! nodes a bit easier. Note that for an axial interface
                 ! these are actually the axial coordinates.

                 rad(1) = xFace(i-1,j-1,1)
                 rad(2) = xFace(i,  j-1,1)
                 rad(3) = xFace(i,  j  ,1)
                 rad(4) = xFace(i-1,j,  1)

                 theta(1) = xFace(i-1,j-1,3)
                 theta(2) = xFace(i,  j-1,3)
                 theta(3) = xFace(i,  j  ,3)
                 theta(4) = xFace(i-1,j,  3)

                 ! Also store the true radii. Both for an axial and a
                 ! radial interface this info is stored in the second
                 ! coordinate of xFace.

                 trueRad(1) = xFace(i-1,j-1,2)
                 trueRad(2) = xFace(i,  j-1,2)
                 trueRad(3) = xFace(i,  j  ,2)
                 trueRad(4) = xFace(i-1,j,  2)

                 ! Compute the average value of theta. This may be
                 ! corrected later on.

                 thetaAvg = fourth*(theta(1) + theta(2) &
                           +        theta(3) + theta(4))

                 ! Compute the minimum and the maximum value of the
                 ! radius.

                 rMin = min(rad(1),rad(2),rad(3),rad(4))
                 rMax = max(rad(1),rad(2),rad(3),rad(4))

                 ! Determine the number of nodes in the 4 quadrants of
                 ! the polar plane.

                 nk(1) = 0; nk(2) = 0; nk(3) = 0; nk(4) = 0

                 do jj=1,4
                   if(theta(jj) > zero) then
                     if(theta(jj) > half*pi) then
                       nk(2) = nk(2) + 1
                     else
                       nk(1) = nk(1) + 1
                     endif
                   else
                     if(theta(jj) > -half*pi) then
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

                 statusQuad = normalQuad

                 select case (nnzq)

                   case (4_intType, 3_intType)

                     ! Rotation center is inside the quad, i.e. the
                     ! the polar singularity is inside the quad.
                     ! Set both the angle and the minimum radius to zero.

                     thetaAvg   = zero
                     rMin       = zero
                     statusQuad = polarSingular

                   !=====================================================

                   case (2_intType)

                     ! Cell has nodes in two quadrants.
                     ! Investigate the situation a bit further.

                     if(nk(2) > 0 .and. nk(3) > 0) then

                       ! The face crosses the line theta = pi. Update the
                       ! average angle such that all angles correspond to
                       ! a positive angle.

                       thetaAvg = thetaAvg + nk(3)*half*pi

                       ! Add two pi to the negative angles (the ones in
                       ! the third quadrant), such that they angles
                       ! change continuously.

                       if(theta(1) < zero) theta(1) = theta(1) + two*pi
                       if(theta(2) < zero) theta(2) = theta(2) + two*pi
                       if(theta(3) < zero) theta(3) = theta(3) + two*pi
                       if(theta(4) < zero) theta(4) = theta(4) + two*pi

                     else if((nk(1) > 0 .and. nk(3) > 0) .or. &
                             (nk(2) > 0 .and. nk(4) > 0)) then

                       ! Face has nodes in two opposite quadrants. Only
                       ! possible when the polar singularity is inside
                       ! the quad. Set both the angle and the minimum
                       ! radius to zero.

                       thetaAvg   = zero
                       rMin       = zero
                       statusQuad = polarSingular

                     endif

                 end select

                 ! Determine the transformation matrix to transform
                 ! the velocity componentes from the global Cartesian
                 ! frame to the local cylindrical frame.

                 cosTheta = cos(thetaAvg)
                 sinTheta = sin(thetaAvg)

                 rotMatDonor(ii,1,1) = axis(1)
                 rotMatDonor(ii,1,2) = axis(2)
                 rotMatDonor(ii,1,3) = axis(3)

                 rotMatDonor(ii,2,1) = cosTheta*vec1(1) &
                                     + sinTheta*vec2(1)
                 rotMatDonor(ii,2,2) = cosTheta*vec1(2) &
                                     + sinTheta*vec2(2)
                 rotMatDonor(ii,2,3) = cosTheta*vec1(3) &
                                     + sinTheta*vec2(3)

                 rotMatDonor(ii,3,1) = cosTheta*vec2(1) &
                                     - sinTheta*vec1(1)
                 rotMatDonor(ii,3,2) = cosTheta*vec2(2) &
                                     - sinTheta*vec1(2)
                 rotMatDonor(ii,3,3) = cosTheta*vec2(3) &
                                     - sinTheta*vec1(3)
!
!                ********************************************************
!                *                                                      *
!                * Determine the actual interpolation information.      *
!                *                                                      *
!                ********************************************************
!
                 ! Initialize the number of intervals to which this donor
                 ! contributes. As the array nIntervalsDonor is in
                 ! cumulative storage format, it is initialized to the
                 ! value of the previous donor.

                 jj = nIntervalsDonor(ii-1)

                 ! Binary algorithm to find the first interval from the
                 ! list that intersects with the current face.

                 start = 1
                 end   = nMixingPoints
                 do
                   ! Compare the median and get rid of either the
                   ! left or right neighbors.

                   interval = (start+end)/2

                   if(rmin < mixingPoints(interval+1) .and. &
                      rMax > mixingPoints(interval)) then
                     end = interval
                   else if(rmax <= mixingPoints(interval)) then
                     end = interval
                   else
                     start = interval
                   endif

                   ! Condition to exit the loop.

                   if((start == interval) .and. (end-start <= 1)) exit

                 enddo

                 ! End contains the correct interval. Copy the value
                 ! into interval, because this variable is used below.

                 interval = end

                 ! Loop to determine the contributions to the
                 ! different intervals.

                 interpolDataLoop: do

                   ! Conditions to exit the loop.

                   if(rmax <= mixingPoints(interval) .or. &
                      interval == nMixingPoints) exit interpolDataLoop

                   ! Check if enough memory has been allocated for the
                   ! donor interpolation info. If not, reallocate.
                   ! Update the counter jj afterwards.

                   if(jj == nAlloc) then
                     newSize = nAlloc + max(2*(nDonor-ii),10_intType)
                     call reallocDonorInterpolInfo(newSize)
                   endif

                   jj = jj + 1

                   ! Store the interval number.

                   indListDonor(jj) = interval

                   ! Store the two coordinates of the interval a bit
                   ! easier.

                   xMin = mixingPoints(interval)
                   xMax = mixingPoints(interval+1)

                   ! Check the status of the quad. A polar singular
                   ! quad requires special attention.

                   select case (statusQuad)

                     case (polarSingular)

                       ! Polar singularity is inside the quad. Set the
                       ! integral to the surface of a circle.

                       r1 = max(rmin, xMin)
                       r2 = min(rmax, xMax)

                       weightDonor(jj) = pi*(r2*r2 - r1*r1)

                     !===================================================

                     case default

                       ! Normal quad.

                       ! First determine the status of the nodes of the
                       ! quad. They could be below/above/contained inside
                       ! the lines xMin and xMax. This is done outside of
                       ! the loop determining the nodes of the CV so that
                       ! floating point comparisons do not give different
                       ! status for the same node when on the boundary.

                       nodesLoop: do kk=1,4

                         if(rad(kk) < xMin) then
                           statusNodes(kk) = below
                         else if(rad(kk) > xMax) then
                           statusNodes(kk) = above
                         else
                           statusNodes(kk) = contained
                         endif

                       enddo nodesLoop

                       ! Determine the nodes of the control volume
                       ! obtained by the current quad bounded by the
                       ! lines xMin and xMax. They are obtained by
                       ! looping over the 4 edges of the quadrilateral.

                       nNodesCV = 0
                       edgeCVLoop: do kk=1,4

                         ! Store the second node of the edge in ll.

                         ll = kk + 1
                         if(ll == 5) ll = 1

                         ! If this edge is either below or above the
                         ! interval, continue with the next one.

                         if ((statusNodes(kk) == below .and. &
                              statusNodes(ll) == below) .or. &
                             (statusNodes(kk) == above .and. &
                              statusNodes(ll) == above)) cycle

                         ! Edge is part of the control volume.
                         ! Determine the situation we are dealing with.

                         if (statusNodes(kk) == contained .and. &
                             statusNodes(ll) == contained) then

                           ! The entire edge belongs to the control
                           ! volume. Store the second node; the first
                           ! node is taken care of by the neighboring
                           ! edge.

                           nNodesCV = nNodesCV + 1
                           xCV(nNodesCV,1) = rad(ll)
                           xCV(nNodesCV,2) = theta(ll)
                           xCV(nNodesCV,3) = trueRad(ll)

                         else if((statusNodes(kk) == below .and. &
                                  statusNodes(ll) == above) .or. &
                                 (statusNodes(kk) == above .and. &
                                  statusNodes(ll) == below)) then

                           ! The edge is cut by both the lower and upper
                           ! bound of the interval. Both intersections
                           ! must be stored in the correct sequence.

                           ! The intersecion point closest to node kk
                           ! must be stored first. Determine the 1st
                           ! coordinate of this point as well as the
                           ! second intersection point.

                           if(rad(kk) > rad(ll)) then
                             xCV(nNodesCV+1,1) = xMax
                             xCV(nNodesCV+2,1) = xMin
                           else
                             xCV(nNodesCV+1,1) = xMin
                             xCV(nNodesCV+2,1) = xMax
                           endif

                           ! Interpolate the two other coordinates for
                           ! both intersection points.

                           do k=1,2
                             nNodesCV = nNodesCV + 1

                             ww1 = sign(one,rad(ll) - rad(kk))  &
                                 * max( abs(rad(ll) - rad(kk)), eps)
                             ww1 = (xCV(nNodesCV,1) - rad(kk)) &
                                 / ww1 
                             ww1 = min(one, max(zero, ww1))
                             ww2 = one - ww1

                             xCV(nNodesCV,2) = ww1*theta(ll) &
                                             + ww2*theta(kk)
                             xCV(nNodesCV,3) = ww1*trueRad(ll) &
                                             + ww2*trueRad(kk)
                           enddo

                         else if(statusNodes(kk) == contained .and. &
                                 statusNodes(ll) /= contained) then

                           ! Edge is only cut by xMin/xMax. Store the
                           ! intersection point.

                           ! Depending on if cut by xMin/xMax, situation
                           ! is different

                           if (statusNodes(ll) == above) then
                             xInter = xMax
                           else
                             xInter = xMin
                           endif 

                           nNodesCV = nNodesCV + 1

                           ww1 = sign(one, rad(ll) - rad(kk)) &
                               * max(  abs(rad(ll) - rad(kk)), eps)
                           ww1 = (xInter - rad(kk))         &
                               / ww1
                           ww1 = min(one, max(zero, ww1))
                           ww2 = one - ww1

                           xCV(nNodesCV,1) = xInter
                           xCV(nNodesCV,2) = ww1*theta(ll) &
                                           + ww2*theta(kk)
                           xCV(nNodesCV,3) = ww1*trueRad(ll) &
                                           + ww2*trueRad(kk)

                         else

                           ! kk is not contained and ll is contained.
                           ! Store the intersection

                           ! Depending on if cut by xMin/xMax, situation
                           ! is different

                           if (statusNodes(kk) == above) then
                             xInter = xMax
                           else
                             xInter = xMin
                           endif

                           nNodesCV = nNodesCV + 1

                           ww1 = sign(one, rad(ll) - rad(kk)) &
                               * max(  abs(rad(ll) - rad(kk)), eps)
                           ww1 = (xInter - rad(kk))         &
                               / ww1
                           ww1 = min(one, max(zero, ww1))
                           ww2 = one - ww1

                           xCV(nNodesCV,1) = xInter
                           xCV(nNodesCV,2) = ww1*theta(ll) &
                                           + ww2*theta(kk)
                           xCV(nNodesCV,3) = ww1*trueRad(ll) &
                                           + ww2*trueRad(kk)

                           ! store the coordinates of the second node as well. 

                           nNodesCV = nNodesCV + 1

                           xCV(nNodesCV,1) = rad(ll)
                           xCV(nNodesCV,2) = theta(ll)
                           xCV(nNodesCV,3) = trueRad(ll)

                         endif

                       enddo edgeCVLoop

                       ! To compute the integral over the control volume
                       ! the control volume is split into triangles and
                       ! the contributions from these triangles are
                       ! added.

                       localWeight = zero
                       triangleLoop: do kk=3,nNodesCV

                         ! Store the indices of the 3 contributing nodes
                         ! of the triangle in n1, n2 and n3.

                         n1 = 1
                         n2 = kk - 1
                         n3 = kk

                         ! Compute twice the area in the polar frame.

                         area = (xCV(n2,1) - xCV(n1,1)) &
                              * (xCV(n3,2) - xCV(n1,2)) &
                              - (xCV(n3,1) - xCV(n1,1)) &
                              * (xCV(n2,2) - xCV(n1,2))

                         ! Compute the contribution to the weight. The
                         ! abs is present, because area can be negative
                         ! due to the sequence of the cross product.

                         localWeight = localWeight                &
                                     + sixth*abs(area*(xCV(n1,3)  &
                                     +                 xCV(n2,3)  &
                                     +                 xCV(n3,3)))

                       enddo triangleLoop

                       weightDonor(jj) = localWeight

                   end select

                   ! Update interval such that the next interval can be
                   ! checked.

                   interval = interval + 1

                 enddo interpolDataLoop

                 ! Check if this face contributes to at least one
                 ! interval. If not print an error message and exit.

              !  if(jj == nIntervalsDonor(ii-1))         &
              !    call terminate("mixingDonorInterpol", &
              !                   "Face does not contribute to a &
              !                   &single interval.")

                 ! Store the nIntervalsDonor for this donor.

                 nIntervalsDonor(ii) = jj

               enddo iloop
             enddo jloop

             ! Deallocate the memory of xFace again.

             deallocate(xFace, stat=ierr)
             if(ierr /= 0)                           &
               call terminate("mixingDonorInterpol", &
                              "Deallocation failure for xFace.")

           endif testInterface2
         enddo bocoLoop2
       enddo domainLoop2

       ! Deallocate the memory for subRanges.

       deallocate(subRange, stat=ierr)
       if(ierr /= 0)                           &
         call terminate("mixingDonorInterpol", &
                        "Deallocation failure for subRange.")

       ! Reallocate the memory for the donor interpolation info, such
       ! that not too much memory is allocated. The waste of memory is
       ! probably not so much, but it is easier to find bugs when the
       ! correct amount has been allocated.

       call reallocDonorInterpolInfo(nIntervalsDonor(nDonor))
!
!      ******************************************************************
!      *                                                                *
!      * The currently stored interpolation weights of the donors are   *
!      * the  values of the local integrals. To obtain to true weights  *
!      * this must be divided by the global sum of these integrals.     *
!      *                                                                *
!      ******************************************************************
!
       ! First determine the local sum.

       do i=1,nInter
         localSum(i) = zero
       enddo

       mm = nIntervalsDonor(nDonor)
       do i=1,mm
         nn = indListDonor(i)
         localSum(nn) = localSum(nn) + weightDonor(i)
       enddo

       ! Determine the global sum via an allreduce.

       size = nInter
       call mpi_allreduce(localSum, globalSum, size, sumb_real, &
                          mpi_sum, myInterfaces(color)%commSlide, ierr)

       ! Invert the global sum to avoid some divisions later on.

       do i=1,nInter
         globalSum(i) = one/max(globalSum(i),eps)
       enddo

       ! Determine correct interpolation weights.

       mm = nIntervalsDonor(nDonor)
       do i=1,mm
         nn = indListDonor(i)
         weightDonor(i) = weightDonor(i)*globalSum(nn)
       enddo

       !=================================================================

       contains

         !===============================================================

         subroutine reallocDonorInterpolInfo(sizeNew)
!
!        ****************************************************************
!        *                                                              *
!        * ReallocDonorInterpolInfo reallocates the memory for the      *
!        * pointers indListDonor and weightDonor of the current         *
!        * mixing plane communication pattern.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: sizeNew
!
!        Local variables.
!
         integer :: ierr

         integer(kind=intType) :: i, sizeCopy
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the number of entities to be copied. This is the
         ! minimum of the currently allocated size and the new size.
         ! Set the new value of nAlloc afterwards.

         sizeCopy = min(nAlloc, sizeNew)
         nAlloc   = sizeNew

         ! Allocate the memory of the arrays indListDonor and
         ! weightDonor in commPattern. As pointers are already
         ! set to these arrays it is not needed to use a temporary
         ! pointer array here.

         allocate(commPattern%indListDonor(nAlloc), &
                  commPattern%weightDonor(nAlloc), stat=ierr)
         if(ierr /= 0) &
           call terminate("reallocDonorInterpolInfo", &
                          "Memory allocation failure for indListDonor &
                          &and weightDonor.")

         ! Copy the relevant data from the original arrays and
         ! deallocate them.

         do i=1,sizeCopy
           commPattern%indListDonor(i) = indListDonor(i)
           commPattern%weightDonor(i)  = weightDonor(i)
         enddo

         deallocate(indListDonor, weightDonor, stat=ierr)
         if(ierr /= 0)                                &
           call terminate("reallocDonorInterpolInfo", &
                          "Deallocation failure for indListDonor &
                          &and weightDonor.")

         ! Reset the pointers for indListDonor and weightDonor.

         indListDonor => commPattern%indListDonor
         weightDonor  => commPattern%weightDonor

         end subroutine reallocDonorInterpolInfo

       end subroutine mixingDonorInterpol

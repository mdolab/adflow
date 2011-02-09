!
!      ******************************************************************
!      *                                                                *
!      * File:          localViscousSurfaceMesh.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-12-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine localViscousSurfaceMesh(multSections, level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * localViscousSurfaceMesh stores the local viscous surface       *
!      * mesh (with possible periodic extensions in conn and coor.      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use communication
       use section
       use viscSurface
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
       integer(kind=intType), dimension(*), intent(in) :: multSections
!
!      Local variables.
!
       integer :: size, ierr

       integer(kind=intType) :: nn, mm, i, j, k
       integer(kind=intType) :: np, nq, npOld, np1, nqOld, mp, mq
       integer(kind=intType) :: nq1, nq2, nq3, nq4, sec, row, col
       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd

       integer(kind=intType), dimension(3) :: ind

       real(kind=realType) :: length, dot, xx, yy, zz, r1, r2, aa, bb
       real(kind=realType) :: theta, cosTheta, sinTheta

       real(kind=realType), dimension(3,3) :: a
       real(kind=realType), dimension(nSections) :: thetaNMin, &
                                                    thetaNMax, &
                                                    thetaPMin, &
                                                    thetaPMax, tmp
       real(kind=realType), dimension(nSections,3) :: rad1, rad2, axis

       real(kind=realType), dimension(:,:,:), pointer :: xface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the unit vectors of the local coordinate system      *
!      * aligned with the rotation axis of the possible rotational      *
!      * periodic section.                                              *
!      *                                                                *
!      ******************************************************************
!
       do nn=1,nSections
         if(sections(nn)%nSlices == 1) cycle

         ! Section is rotational periodic. First determine the rotation
         ! axis. This is the eigenvector which corresponds to the
         ! eigenvalue 1 of the transformation matrix.

         ! Store rot - i in a and initialize ind.

         ind(1) = 1; ind(2) = 2; ind(3) = 3

         a(1,1) = sections(nn)%rotMatrix(1,1) - one
         a(1,2) = sections(nn)%rotMatrix(1,2)
         a(1,3) = sections(nn)%rotMatrix(1,3)

         a(2,1) = sections(nn)%rotMatrix(2,1)
         a(2,2) = sections(nn)%rotMatrix(2,2) - one
         a(2,3) = sections(nn)%rotMatrix(2,3)

         a(3,1) = sections(nn)%rotMatrix(3,1)
         a(3,2) = sections(nn)%rotMatrix(3,2)
         a(3,3) = sections(nn)%rotMatrix(3,3) - one

         ! Loop over the two times that Gaussian elimination must be
         ! applied.

         loopGauss: do k=1,2

           ! Find the largest value in the sub-matrix.

           aa = abs(a(k,k)); row = k; col = k
           do j=k,3
             do i=k,3
               bb = abs(a(i,j))
               if(bb > aa) then
                 aa  = bb
                 row = i
                 col = j
               endif
             enddo
           enddo

           ! Swap the rows k and row.

           do j=1,3
             aa       = a(k,j)
             a(k,j)   = a(row,j)
             a(row,j) = aa
           enddo

           ! Swap the colums k and col; also swap ind(k) and ind(col).

           i        = ind(k)
           ind(k)   = ind(col)
           ind(col) = i
           do i=1,3
             aa       = a(i,k)
             a(i,k)   = a(i,col)
             a(i,col) = aa
           enddo

           ! Perform gaussian eliMination, because now it's sure that
           ! the element (k,k) is non-zero.

           aa = one/a(k,k)
           do i=(k+1),3
             bb = a(i,k)*aa
             do j=k,3
               a(i,j) = a(i,j) - bb*a(k,j)
             enddo
           enddo

         enddo loopGauss

         ! Due to the full pivoting it is now guaranteed that the elements
         ! a(ind(1),ind(1)) and a(ind(2),ind(2)) are nonzero and
         ! a(ind(3),ind(3)) == zero. Remember that the rotation matrix
         ! only has 1 eigenvalue of one. Set axis(ind(3)) to one and
         ! determine the other two elements of the eigen vector.

         axis(nn,ind(3)) =   one
         axis(nn,ind(2)) = -(a(2,3)*axis(nn,ind(3)))/a(2,2)
         axis(nn,ind(1)) = -(a(1,3)*axis(nn,ind(3)) &
                         +   a(1,2)*axis(nn,ind(2)))/a(1,1)

         ! Create a unit vector.

         length = one/sqrt(axis(nn,1)**2 + axis(nn,2)**2 + axis(nn,3)**2)
         axis(nn,1) = axis(nn,1)*length
         axis(nn,2) = axis(nn,2)*length
         axis(nn,3) = axis(nn,3)*length

         ! Make sure that the largest component of this vector is
         ! positive, such that a unique definition of the rotation
         ! axis is obtained. Use dot and length as a temporary
         ! storage.

         dot = axis(nn,1); length = abs(dot)
         if(abs(axis(nn,2)) > length) then
           dot = axis(nn,2); length = abs(dot)
         endif
         if(abs(axis(nn,3)) > length) then
           dot = axis(nn,3); length = abs(dot)
         endif

         if(dot < zero) then
           axis(nn,1) = -axis(nn,1)
           axis(nn,2) = -axis(nn,2)
           axis(nn,3) = -axis(nn,3)
         endif

         ! Determine the two vectors which determine the plane normal
         ! to the axis of rotation.

         ! Initial guess of rad1. First try the y-axis. If not good
         ! enough try the z-axis.

         if(abs(axis(nn,2)) < 0.707107_realType) then
           rad1(nn,1) = zero
           rad1(nn,2) = one
           rad1(nn,3) = zero
         else
           rad1(nn,1) = zero
           rad1(nn,2) = zero
           rad1(nn,3) = one
         endif

         ! Make sure that rad1 is normal to axis. Create a unit
         ! vector again.

         dot = rad1(nn,1)*axis(nn,1) + rad1(nn,2)*axis(nn,2) &
             + rad1(nn,3)*axis(nn,3)
         rad1(nn,1) = rad1(nn,1) - dot*axis(nn,1)
         rad1(nn,2) = rad1(nn,2) - dot*axis(nn,2)
         rad1(nn,3) = rad1(nn,3) - dot*axis(nn,3)

         length = one/(rad1(nn,1)**2 + rad1(nn,2)**2 + rad1(nn,3)**2)
         rad1(nn,1) = rad1(nn,1)*length
         rad1(nn,2) = rad1(nn,2)*length
         rad1(nn,3) = rad1(nn,3)*length

         ! Create the second vector which spans the radIal plane. This
         ! must be normal to both axis and rad1, i.e. the cross-product.

         rad2(nn,1) = axis(nn,2)*rad1(nn,3) - axis(nn,3)*rad1(nn,2)
         rad2(nn,2) = axis(nn,3)*rad1(nn,1) - axis(nn,1)*rad1(nn,3)
         rad2(nn,3) = axis(nn,1)*rad1(nn,2) - axis(nn,2)*rad1(nn,1)

       enddo

       ! Initialize the values of thetaNMin, etc.

       thetaNMin =  zero
       thetaNMax = -pi
       thetaPMin =  pi
       thetaPMax =  zero
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local values of thetaNMin, etc. for the          *
!      * different sections.                                            *
!      *                                                                *
!      ******************************************************************
!
       do nn=1,nDom

         ! Store the section id of the block a bit easier. Continue with
         ! the next block if this section consist of only one slice.

         sec = flowDoms(nn,level,sps)%sectionID
         if(sections(sec)%nSlices == 1) cycle

         ! Set the pointers for this block.

         call setPointers(nn, level,sps)

         ! Initialize nq1, nq2, nq3 and nq4 to 0. These integers store
         ! the number of nodes in the first, second, third and fourth
         ! quadrant respectivily.

         nq1 = 0; nq2 = 0; nq3 = 0; nq4 = 0

         ! Loop over the nodes of this block.

         do k=1,kl
           do j=1,jl
             do i=1,il

               ! Determine the coordinates relative to the
               ! center of rotation.

               xx = x(i,j,k,1) - sections(sec)%rotCenter(1)
               yy = x(i,j,k,2) - sections(sec)%rotCenter(2)
               zz = x(i,j,k,3) - sections(sec)%rotCenter(3)

               ! Determine the radIal components in the local
               ! cylindrical coordinate system of the section.

               r1 = xx*rad1(sec,1) + yy*rad1(sec,2) + zz*rad1(sec,3)
               r2 = xx*rad2(sec,1) + yy*rad2(sec,2) + zz*rad2(sec,3)

               ! Determine the angle if r1 or r2 is nonzero.

               if((abs(r1) >= eps) .or. (abs(r2) >= eps)) then

                 theta = atan2(r2,r1)

                 ! Update the minimum and maximum angle for this
                 ! section, depending on the sign of theta.

                 if(theta >= zero) then
                   thetaPMin(sec) = min(thetaPMin(sec),theta)
                   thetaPMax(sec) = max(thetaPMax(sec),theta)
                 endif

                 if(theta <= zero) then
                   thetaNMin(sec) = min(thetaNMin(sec),theta)
                   thetaNMax(sec) = max(thetaNMax(sec),theta)
                 endif

                 ! Determine the quadrant in which this node is located
                 ! and update the corresponding counter.

                 if(theta <= -half*pi) then
                   nq3 = nq3 + 1
                 else if(theta <= zero) then
                   nq4 = nq4 + 1
                 else if(theta <= half*pi) then
                   nq1 = nq1 + 1
                 else
                   nq2 = nq2 + 1
                 endif

               endif

             enddo
           enddo
         enddo

         ! Modify the minimum and maximum angles if nodes are present
         ! in multiple quadrants.

         if(nq1 > 0 .and. nq4 > 0) then

           ! Nodes in both the 1st and 4th quadrant. Update the
           ! corresponding minimum and maximum angle.

           thetaNMax(sec) = zero
           thetaPMin(sec) = zero

         endif

         if(nq2 > 0 .and. nq3 > 0) then

           ! Nodes in both the 2nd and 3rd quadrant. Update the
           ! corresponding minimum and maximum angle.

           thetaNMin(sec) = -pi
           thetaPMax(sec) =  pi

         endif

       enddo

       ! Determine the minimum of the minimum angles and the maximum of
       ! the maximum angles for all sections.

       size = nSections
       call mpi_allreduce(thetaNMax, tmp, size, sumb_real, mpi_max, &
                          SUmb_comm_world, ierr)
       thetaNMax = tmp

       call mpi_allreduce(thetaPMax, tmp, size, sumb_real, mpi_max, &
                          SUmb_comm_world, ierr)
       thetaPMax = tmp

       call mpi_allreduce(thetaNMin, tmp, size, sumb_real, mpi_min, &
                          SUmb_comm_world, ierr)
       thetaNMin = tmp

       call mpi_allreduce(thetaPMin, tmp, size, sumb_real, mpi_min, &
                          SUmb_comm_world, ierr)
       thetaPMin = tmp

       ! Allocate the memory for rotMatrixSections, the rotation
       ! matrices of the sections needed for the alignment.

       allocate(rotMatrixSections(nSections,3,3), stat=ierr)
       if(ierr /= 0)                                  &
         call terminate("localViscousSurfaceMesh", &
                        "Memory allocation failure for &
                        &rotMatrixSections")
!
!      ******************************************************************
!      *                                                                *
!      * Determine the rotation matrix for each section, which aligns   *
!      * the rotational periodic sections with other sections.          *
!      *                                                                *
!      ******************************************************************
!
       do nn=1,nSections

         ! Test if a rotation is actually needed.

         testRot: if(sections(nn)%nSlices == 1 .or. &
                      thetaPMin(nn) == zero) then

           ! Section consist out of 1 slice or the slice crosses the
           ! line theta == 0. For both cases the rotation matrix is
           ! the identity matrix.

           rotMatrixSections(nn,1,1) = one
           rotMatrixSections(nn,1,2) = zero
           rotMatrixSections(nn,1,3) = zero

           rotMatrixSections(nn,2,1) = zero
           rotMatrixSections(nn,2,2) = one
           rotMatrixSections(nn,2,3) = zero

           rotMatrixSections(nn,3,1) = zero
           rotMatrixSections(nn,3,2) = zero
           rotMatrixSections(nn,3,3) = one

         else testRot

           ! Section consist out of multiple slices and the current slice
           ! does not cross the line theta == 0. The rotation matrix
           ! for alignment must be computed.

           theta = two*pi/sections(nn)%nSlices

           ! Determine the number of rotations needed to align the mesh.

           if(thetaNMin(nn) < zero) then

             ! The section lies (at least partially) in the third and
             ! fourth quadrant. Determine the number of rotations for
             ! alignment; this is a positive number.

             mm = -thetaNMax(nn)/theta + 1

           else

             ! The section lies (completely) in the first and second
             ! quadrant. The number of rotations will be a negative
             ! number now.

             mm = -thetaPMin(nn)/theta - 1

           endif

           ! Compute the rotation angle in the local cylindrical frame
           ! and its sine and cosine.

           theta = mm*theta
           cosTheta = cos(theta)
           sinTheta = sin(theta)

           ! Apply the transformation to obtain the matrix in the
           ! original cartesian frame.

           rotMatrixSections(nn,1,1) = axis(nn,1)*axis(nn,1)             &
              + cosTheta*(rad1(nn,1)*rad1(nn,1) + rad2(nn,1)*rad2(nn,1))
           rotMatrixSections(nn,1,2) = axis(nn,1)*axis(nn,2)             &
              + cosTheta*(rad1(nn,1)*rad1(nn,2) + rad2(nn,1)*rad2(nn,2)) &
              + sinTheta*(rad1(nn,2)*rad2(nn,1) - rad1(nn,1)*rad2(nn,2))
           rotMatrixSections(nn,1,3) = axis(nn,1)*axis(nn,3)             &
              + cosTheta*(rad1(nn,1)*rad1(nn,3) + rad2(nn,1)*rad2(nn,3)) &
              + sinTheta*(rad1(nn,3)*rad2(nn,1) - rad1(nn,1)*rad2(nn,3))

           rotMatrixSections(nn,2,1) = axis(nn,1)*axis(nn,2)             &
              + cosTheta*(rad1(nn,1)*rad1(nn,2) + rad2(nn,1)*rad2(nn,2)) &
              - sinTheta*(rad1(nn,2)*rad2(nn,1) - rad1(nn,1)*rad2(nn,2))
           rotMatrixSections(nn,2,2) = axis(nn,2)*axis(nn,2)             &
              + cosTheta*(rad1(nn,2)*rad1(nn,2) + rad2(nn,2)*rad2(nn,2))
           rotMatrixSections(nn,2,3) = axis(nn,2)*axis(nn,3)             &
              + cosTheta*(rad1(nn,2)*rad1(nn,3) + rad2(nn,2)*rad2(nn,3)) &
              + sinTheta*(rad1(nn,3)*rad2(nn,2) - rad1(nn,2)*rad2(nn,3))

           rotMatrixSections(nn,3,1) = axis(nn,1)*axis(nn,3)             &
              + cosTheta*(rad1(nn,1)*rad1(nn,3) + rad2(nn,1)*rad2(nn,3)) &
              - sinTheta*(rad1(nn,3)*rad2(nn,1) - rad1(nn,1)*rad2(nn,3))
           rotMatrixSections(nn,3,2) = axis(nn,2)*axis(nn,3)             &
              + cosTheta*(rad1(nn,2)*rad1(nn,3) + rad2(nn,2)*rad2(nn,3)) &
              - sinTheta*(rad1(nn,3)*rad2(nn,2) - rad1(nn,2)*rad2(nn,3))
           rotMatrixSections(nn,3,3) = axis(nn,3)*axis(nn,3)             &
              + cosTheta*(rad1(nn,3)*rad1(nn,3) + rad2(nn,3)*rad2(nn,3))

         endif testRot

       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local viscous surface grid.                      *
!      *                                                                *
!      ******************************************************************
!
       np = 0
       nq = 0

       loopDomains: do nn=1,nDom

         ! Set the pointers for this block and store the section id
         ! a bit easier.

         call setPointers(nn, level, sps)
         sec = sectionID

         ! Loop over the subfaces of this block and test if this is
         ! a viscous subface.

         loopBocos: do mm=1,nBocos
           testViscous: if(BCType(mm) == NSWallAdiabatic .or. &
                           BCType(mm) == NSWallIsothermal) then

             ! Viscous subface. Set the pointer for the coordinates of
             ! the face.

             select case (BCFaceID(mm))

               case (iMin)
                 xface => x(1,1:,1:,:)

               case (iMax)
                 xface => x(il,1:,1:,:)

               case (jMin)
                 xface => x(1:,1,1:,:)

               case (jMax)
                 xface => x(1:,jl,1:,:)

               case (kMin)
                 xface => x(1:,1:,1,:)

               case (kMax)
                 xface => x(1:,1:,kl,:)

             end select

             ! Store the nodal range of this subface a bit easier.

             jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd

             ! Store the old value of the number of points stored and
             ! determine the new coordinates.

             npOld = np
             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 ! Determine the coordinates relative to the rotation
                 ! center of this section.

                 xx = xface(i,j,1) - sections(sec)%rotCenter(1)
                 yy = xface(i,j,2) - sections(sec)%rotCenter(2)
                 zz = xface(i,j,3) - sections(sec)%rotCenter(3)

                 ! Update the counter and determine the surface mesh
                 ! coordinates.

                 np = np + 1
                 coorVisc(1,np) = rotMatrixSections(sec,1,1)*xx &
                                + rotMatrixSections(sec,1,2)*yy &
                                + rotMatrixSections(sec,1,3)*zz &
                                + sections(sec)%rotCenter(1)

                 coorVisc(2,np) = rotMatrixSections(sec,2,1)*xx &
                                + rotMatrixSections(sec,2,2)*yy &
                                + rotMatrixSections(sec,2,3)*zz &
                                + sections(sec)%rotCenter(2)

                 coorVisc(3,np) = rotMatrixSections(sec,3,1)*xx &
                                + rotMatrixSections(sec,3,2)*yy &
                                + rotMatrixSections(sec,3,3)*zz &
                                + sections(sec)%rotCenter(3)
               enddo
             enddo

             ! Determine and store the connectivity of this subface.

             np1 = iEnd - iBeg + 1
             nqOld = nq

             do j=(jBeg+1),jEnd
               do i=(iBeg+1),iEnd

                 ! Update the counter nq and determine the 4 indices
                 ! of the surface quad.

                 nq = nq + 1

                 connVisc(1,nq) = npOld + (j-jBeg-1)*np1 + i - iBeg
                 connVisc(2,nq) = connVisc(1,nq) + 1
                 connVisc(3,nq) = connVisc(2,nq) + np1
                 connVisc(4,nq) = connVisc(3,nq) - 1

               enddo
             enddo

             ! Loop over the number of times the subface must be stored.
             ! This happens when the rotational periodicity differs from
             ! section to section. Note that this loop starts at k == 2.

             loopMultiplicity: do k=2,multSections(sec)

               ! Store the current number of nodes and quads
               ! in mp and mq respectivily.

               mp = np
               mq = nq

               ! Loop over the of points on this subface.

               do i=(npOld+1),mp

                 ! Determine the coordinates relative to the center
                 ! of rotation.

                 np = np + 1

                 xx = coorVisc(1,i) - sections(sec)%rotCenter(1)
                 yy = coorVisc(2,i) - sections(sec)%rotCenter(2)
                 zz = coorVisc(3,i) - sections(sec)%rotCenter(3)

                 ! Update the counter np and determine the new
                 ! coordinates after the transformation.

                 coorVisc(1,np) = sections(sec)%rotMatrix(1,1)*xx &
                                + sections(sec)%rotMatrix(1,2)*yy &
                                + sections(sec)%rotMatrix(1,3)*zz &
                                + sections(sec)%rotCenter(1)      &
                                + sections(sec)%translation(1)

                 coorVisc(2,np) = sections(sec)%rotMatrix(2,1)*xx &
                                + sections(sec)%rotMatrix(2,2)*yy &
                                + sections(sec)%rotMatrix(2,3)*zz &
                                + sections(sec)%rotCenter(2)      &
                                + sections(sec)%translation(2)
 
                 coorVisc(3,np) = sections(sec)%rotMatrix(3,1)*xx &
                                + sections(sec)%rotMatrix(3,2)*yy &
                                + sections(sec)%rotMatrix(3,3)*zz &
                                + sections(sec)%rotCenter(3)      &
                                + sections(sec)%translation(3)
               enddo

               ! Store the number of nodes in this subface in j
               ! and determine the connectivity of this rotated part.

               j = np - mp
               do i=(nqOld+1),mq

                 ! Update the counter nq and set the new connectivity,
                 ! which is the old connectivity plus an offset.

                 nq = nq + 1
                 connVisc(1,nq) = connVisc(1,i) + j
                 connVisc(2,nq) = connVisc(2,i) + j
                 connVisc(3,nq) = connVisc(3,i) + j
                 connVisc(4,nq) = connVisc(4,i) + j

               enddo

               ! Copy the values of mp and mq into npOld and nqOld for
               ! the next multiple of the slice. Idem for np in mp, etc.

               npOld = mp
               nqOld = mq

               mp = np
               mq = nq

             enddo loopMultiplicity

           endif testViscous
         enddo loopBocos
       enddo loopDomains

       end subroutine localViscousSurfaceMesh

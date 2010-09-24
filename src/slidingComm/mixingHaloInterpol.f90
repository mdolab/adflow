!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingHaloInterpol.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-01-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mixingHaloInterpol(level, slideID, color, commPattern)
!
!      ******************************************************************
!      *                                                                *
!      * mixingHaloInterpol determines the interpolation data for the   *
!      * halo cells for the given slideID on the given multigrid level. *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commMixing
       use interfaceGroups
       use mixingData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, slideID, color

       type(commMixingType), intent(inout) :: commPattern
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, nn, mm, i, j
       integer(kind=intType) :: nHalo, nFaces, nnzq, nInter
       integer(kind=intType) :: jBeg, jEnd, iBeg, iEnd
       integer(kind=intType) :: h1, h2, indh, ind1, ind2
       integer(kind=intType) :: start, end

       integer(kind=intType), dimension(4) :: nk

       real(kind=realType) :: thetaAvg, coorInt, weight1, weight2
       real(kind=realType) :: cosTheta, sinTheta

       real(kind=realType), dimension(3) :: axis, vec1, vec2
       real(kind=realType), dimension(4) :: theta

       real(kind=realType), dimension(:,:,:), allocatable :: xFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the vectors which define the local cartesian
       ! coordinate system.

       axis = myInterfaces(color)%rotAxis
       vec1 = myInterfaces(color)%radVec1
       vec2 = myInterfaces(color)%radVec2

       ! Easier storage of the number of interpolation cells.

       nInter = commPattern%nInter

       ! Determine the number of local halo cells to be interpolated.

       nHalo = 0
       do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           if(BCType(mm)   == slidingInterface .and. &
              groupNum(mm) == slideID) then

             ! Determine the number of faces, including the 1st layer
             ! halos. This number is then added to nHalo.
             ! The only way to guarantee that all halo cells are included
             ! is to use the nodal range. The cell range may or may not
             ! contain halo information, while the nodal range does not.
             ! This explains the +2 in the nodal range below.

             nFaces = (BCData(mm)%inEnd - BCData(mm)%inBeg + 2) &
                    * (BCData(mm)%jnEnd - BCData(mm)%jnBeg + 2)

             nHalo  = nHalo + nFaces

           endif
         enddo
       enddo

       ! Allocate the memory for the arrays to store the information
       ! to interpolate the halo data.

       commPattern%nHalo = nHalo

       allocate(commPattern%indListHalo(nHalo,2),  &
                commPattern%blockHalo(nHalo),      &
                commPattern%indH(nHalo,3,2),       &
                commPattern%weightHalo(nHalo,2),   &
                commPattern%rotMatHalo(nHalo,3,3), &
                stat=ierr)
       if(ierr /= 0)                          &
         call terminate("mixingHaloInterpol", &
                        "Memory allocation failure for halo data.")

       ! Repeat the loop over the subfaces of the sliding mesh interface
       ! currently treated and store the interpolation data.

       ii = 0
       domainLoop: do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         bocoLoop: do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           testInterface: if(BCType(mm)   == slidingInterface .and. &
                             groupNum(mm) == slideID) then

             ! Determine the cell range including the halo cells a
             ! bit easier. As the cell range may or may not contain
             ! halo cells, it is easier to use the nodal range.
             ! The latter contains the nodal range of the subface
             ! without any halos.

             jBeg = BCData(mm)%jnBeg
             jEnd = BCData(mm)%jnEnd + 1
             iBeg = BCData(mm)%inBeg
             iEnd = BCData(mm)%inEnd + 1

             ! Set some variables depending on the block face on which
             ! this subface is located.

             select case (BCFaceID(mm))
               case (iMin)
                 h1 =  1; h2 =  0; indh = 1; ind1 = 2; ind2 = 3
               case (iMax)
                 h1 = ie; h2 = ib; indh = 1; ind1 = 2; ind2 = 3
               case (jMin)
                 h1 =  1; h2 =  0; indh = 2; ind1 = 1; ind2 = 3
               case (jMax)
                 h1 = je; h2 = jb; indh = 2; ind1 = 1; ind2 = 3
               case (kMin)
                 h1 =  1; h2 =  0; indh = 3; ind1 = 1; ind2 = 2
               case (kMax)
                 h1 = ke; h2 = kb; indh = 3; ind1 = 1; ind2 = 2
             end select

             ! Allocate the memory for the nodal cylindrical coordinates
             ! of the subface, including the halo's, and determine them.

             allocate(xFace(iBeg-1:iEnd,jBeg-1:jEnd,3), stat=ierr)
             if(ierr /= 0)                          &
               call terminate("mixingHaloInterpol", &
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

             ! Loop over the faces of this subface to store the data
             ! to perform the interpolation.

             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 ! Store the indices of this halo in the communication
                 ! pattern for this mixing plane.

                 ii = ii + 1
                 commPattern%blockHalo(ii) = nn

                 commPattern%indH(ii,indh,1) = h1
                 commPattern%indH(ii,ind1,1) =  i
                 commPattern%indH(ii,ind2,1) =  j

                 commPattern%indH(ii,indh,2) = h2
                 commPattern%indH(ii,ind1,2) =  i
                 commPattern%indH(ii,ind2,2) =  j

                 ! Store the angles of the 4 surrounding nodes
                 ! a bit easier.

                 theta(1) = xFace(i-1,j-1,3)
                 theta(2) = xFace(i,  j-1,3)
                 theta(3) = xFace(i,  j  ,3)
                 theta(4) = xFace(i-1,j,  3)

                 ! Compute the average of theta and the coordinate
                 ! used for the interpolation. The latter is stored in
                 ! the first position of xFace.

                 thetaAvg = fourth*(theta(1) + theta(2) &
                          +         theta(3) + theta(4))
                 coorInt  = fourth*(xFace(i-1,j-1,1) + xFace(i,j-1,1) &
                          +         xFace(i-1,j,  1) + xFace(i,j,  1))

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

                 select case (nnzq)

                   case (4_intType, 3_intType)

                     ! Rotation center is inside the quad, i.e. the
                     ! the polar singularity is inside the quad.
                     ! Set both the angle and average interpolation
                     ! coordinate (radius) to zero.

                     thetaAvg = zero
                     coorInt  = zero

                   !=====================================================

                   case (2_intType)

                     ! Cell has nodes in two quadrants.
                     ! Investigate the situation a bit further.

                     if(nk(2) > 0 .and. nk(3) > 0) then

                       ! The face crosses the line theta = pi. Update the
                       ! average angle such that all angles correspond to
                       ! a positive angle.

                       thetaAvg = thetaAvg + nk(3)*half*pi

                     else if((nk(1) > 0 .and. nk(3) > 0) .or. &
                             (nk(2) > 0 .and. nk(4) > 0)) then

                       ! Face has nodes in two opposite quadrants. Only
                       ! possible when the polar singularity is inside
                       ! the quad. Set both the angle and average
                       ! interpolation coordinate (radius) to zero.

                       thetaAvg = zero
                       coorInt  = zero

                     endif

                 end select

                 ! Determine the two cell centers to be used in
                 ! the interpolation.

                 ! Take care of the exceptional cases.

                 if( coorInt <= mixingCells(1) ) then

                   commPattern%indListHalo(ii,1) = 1
                   commPattern%indListHalo(ii,2) = 2

                   commPattern%weightHalo(ii,1) = one
                   commPattern%weightHalo(ii,2) = zero

                 else if(coorInt >= mixingCells(nInter)) then

                   commPattern%indListHalo(ii,1) = nInter - 1
                   commPattern%indListHalo(ii,2) = nInter

                   commPattern%weightHalo(ii,1) = zero
                   commPattern%weightHalo(ii,2) = one

                 else

                   ! Coordinate is somewhere in the range. Use a
                   ! binary algorithm to find the correct one.

                   start = 1
                   end   = nInter
                   do
                     ! Condition to exit the loop.

                     if(end - start <= 1) exit

                     ! Compare the median and get rid of either the
                     ! left or right neighbors.

                     jj = (start+end)/2
                     if(mixingCells(jj) < coorInt) then
                       start = jj
                     else
                       end = jj
                     endif
                   enddo

                   ! Start and end contain the indices used for the
                   ! interpolation. Store them as well as the weights.

                   commPattern%indListHalo(ii,1) = start
                   commPattern%indListHalo(ii,2) = end

                   weight1 = (coorInt            - mixingCells(end)) &
                           / (mixingCells(start) - mixingCells(end))
                   weight2 = one - weight1

                   commPattern%weightHalo(ii,1) = weight1
                   commPattern%weightHalo(ii,2) = weight2

                 endif

                 ! Determine the transformation matrix to transform
                 ! the velocity componentes from the local cylindrical
                 ! frame to the global cartesian frame.

                 cosTheta = cos(thetaAvg)
                 sinTheta = sin(thetaAvg)

                 commPattern%rotMatHalo(ii,1,1) = axis(1)
                 commPattern%rotMatHalo(ii,2,1) = axis(2)
                 commPattern%rotMatHalo(ii,3,1) = axis(3)

                 commPattern%rotMatHalo(ii,1,2) = cosTheta*vec1(1) &
                                                + sinTheta*vec2(1)
                 commPattern%rotMatHalo(ii,2,2) = cosTheta*vec1(2) &
                                                + sinTheta*vec2(2)
                 commPattern%rotMatHalo(ii,3,2) = cosTheta*vec1(3) &
                                                + sinTheta*vec2(3)

                 commPattern%rotMatHalo(ii,1,3) = cosTheta*vec2(1) &
                                                - sinTheta*vec1(1)
                 commPattern%rotMatHalo(ii,2,3) = cosTheta*vec2(2) &
                                                - sinTheta*vec1(2)
                 commPattern%rotMatHalo(ii,3,3) = cosTheta*vec2(3) &
                                                - sinTheta*vec1(3)

               enddo
             enddo

             ! Deallocate the memory of xFace again.

             deallocate(xFace, stat=ierr)
             if(ierr /= 0)                          &
               call terminate("mixingHaloInterpol", &
                              "Deallocation failure for xFace.")

           endif testInterface
         enddo bocoLoop
       enddo domainLoop

       end subroutine mixingHaloInterpol

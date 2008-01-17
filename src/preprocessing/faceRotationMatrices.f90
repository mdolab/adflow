!
!      ******************************************************************
!      *                                                                *
!      * File:          faceRotationMatrices.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-20-2007                                      *
!      * Last modified: 12-01-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine faceRotationMatrices(level, allocMem)
!
!      ******************************************************************
!      *                                                                *
!      * faceRotationMatrices computes the rotation matrices on the     *
!      * faces, such that for a rotationally periodic the nonlinear     *
!      * reconstruction in the upwind schemes is consistent with its    *
!      * periodic neighbor.                                             *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputDiscretization
       use inputTimeSpectral
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
       logical,               intent(in) :: allocMem
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, sps, mm

       real(kind=realType), dimension(:,:,:),   pointer :: xFace
       real(kind=realType), dimension(:,:,:,:), pointer :: rotFace

       real(kind=realType), dimension(3) :: axis, vecR1, vecR2, rotCenter
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If this is not the finest level, return. The matrices are
       ! only needed on the finest grid.

       if(level /= 1) return

       ! Check if an upwind scheme is used. If not, return.

       if(spaceDiscr /= upwind) return

       ! Check if a limiter is used at all. The rotation matrices are
       ! only needed when a nonlinear construction is used.

       if(limiter == firstOrder .or. limiter == noLimiter) return

       ! Check if rotational periodicity occurs at all. If not, there
       ! is no need for the rotation matrices either.

       do nn=1,nSections
         if(sections(nn)%nSlices > 1) exit
       enddo
       if(nn > nSections) return

       ! Loop over the number of blocks..

       domains: do nn=1,nDom

         sectionID = flowDoms(nn,1,1)%sectionID

         ! Determine the two unit vectors in the plane normal to
         ! the rotation axis of this section.

         axis      = sections(sectionID)%rotAxis
         rotCenter = sections(sectionID)%rotCenter
         call unitVectorsInAxialPlane(axis, vecR1, vecR2)

         ! Loop over the number of time instances.

         spectral: do sps=1,nTimeIntervalsSpectral

           ! Check if the memory for the rotation matrices must be
           ! allocated and do so if needed.

           if( allocMem ) then

             il = flowDoms(nn,1,sps)%il
             jl = flowDoms(nn,1,sps)%jl
             kl = flowDoms(nn,1,sps)%kl

             allocate(flowDoms(nn,1,sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
                      flowDoms(nn,1,sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
                      flowDoms(nn,1,sps)%rotMatrixK(2:il,2:jl,kl,3,3), &
                      stat=ierr)
             if(ierr /= 0)                            &
               call terminate("faceRotationMatrices", &
                              "Memory allocation failure for the &
                              &rotation matrices.")
           endif

           ! Set the pointers to this block.

           call setPointers(nn,level,sps)

           ! The rotation matrices for the i-faces.

           do mm=1,il
             xFace   => x(mm,1:,1:,:);
             rotFace => rotMatrixI(mm,:,:,:,:)

             call computeRotMatrixFace(xFace, rotFace, jl, kl)
           enddo

           ! The rotation matrices for the j-faces.

           do mm=1,jl
             xFace   => x(1:,mm,1:,:);
             rotFace => rotMatrixJ(:,mm,:,:,:)

             call computeRotMatrixFace(xFace, rotFace, il, kl)
           enddo

           ! The rotation matrices for the k-faces.

           do mm=1,kl
             xFace   => x(1:,1:,mm,:);
             rotFace => rotMatrixK(:,:,mm,:,:)

             call computeRotMatrixFace(xFace, rotFace, il, jl)
           enddo

         enddo spectral
       enddo domains

       !=================================================================

       contains

         !===============================================================

         subroutine computeRotMatrixFace(xx, rotMat, iil, jjl)
!
!        ****************************************************************
!        *                                                              *
!        * computeRotMatrixFace is an internal subroutine, which        *
!        * computes the rotation matrix from Cartesian to local         *
!        * cylindrical velocity components for the face centers.        *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: iil, jjl

         real(kind=realType), dimension(:,:,:),     intent(in)  :: xx
         real(kind=realType), dimension(2:,2:,:,:), intent(out) :: rotMat
!
!        Local variables.
!
         integer(kind=intType) :: i, j

         real(kind=realType) :: r1, r2, rInv, cosTheta, sinTheta

         real(kind=realType), dimension(3) :: xF
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the face centers.

         do j=2,jjl
           do i=2,iil

             ! Compute the coordinates of the face center relative to
             ! the center of rotation.

             xF(1) = fourth*(xx(i-1,j-1,1) + xx(i-1,j,1) &
                   +         xx(i,  j-1,1) + xx(i,  j,1)) - rotCenter(1)
             xF(2) = fourth*(xx(i-1,j-1,2) + xx(i-1,j,2) &
                   +         xx(i,  j-1,2) + xx(i,  j,2)) - rotCenter(2)
             xF(3) = fourth*(xx(i-1,j-1,3) + xx(i-1,j,3) &
                   +         xx(i,  j-1,3) + xx(i,  j,3)) - rotCenter(3)

             ! Determine the two radial components for this point.

             r1 = xF(1)*vecR1(1) + xF(2)*vecR1(2) + xF(3)*vecR1(3)
             r2 = xF(1)*vecR2(1) + xF(2)*vecR2(2) + xF(3)*vecR2(3)

             ! Determine the sine and cosine of the polar angle.

             rInv     = one/sqrt(r1*r1 + r2*r2)
             cosTheta = r1*rInv
             sinTheta = r2*rInv

             ! Compute the transformation matrix.

             rotMat(i,j,1,1) = axis(1)
             rotMat(i,j,1,2) = axis(2)
             rotMat(i,j,1,3) = axis(3)

             rotMat(i,j,2,1) = cosTheta*vecR1(1) + sinTheta*vecR2(1)
             rotMat(i,j,2,2) = cosTheta*vecR1(2) + sinTheta*vecR2(2)
             rotMat(i,j,2,3) = cosTheta*vecR1(3) + sinTheta*vecR2(3)

             rotMat(i,j,3,1) = cosTheta*vecR2(1) - sinTheta*vecR1(1)
             rotMat(i,j,3,2) = cosTheta*vecR2(2) - sinTheta*vecR1(2)
             rotMat(i,j,3,3) = cosTheta*vecR2(3) - sinTheta*vecR1(3)

           enddo
         enddo

         end subroutine computeRotMatrixFace

       end subroutine faceRotationMatrices

!
!      ******************************************************************
!      *                                                                *
!      * File:          computeNormalSpacing.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-03-2005                                      *
!      * Last modified: 11-08-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeNormalSpacing(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * computeNormalSpacing computes the normal spacing of the first  *
!      * cell center from the viscous wall for the given multigrid      *
!      * level and spectral solution. This routine is called for        *
!      * turbulence models, which do not need the wall distance.        *
!      * However, they do need info of the first normal spacing for the *
!      * monitoring of y+ and possibly for the boundary conditions.     *
!      * This is computed in this routine.                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

       real(kind=realType) :: nnx, nny, nnz, vecx, vecy, vecz, dot
       
       real(kind=realType), dimension(:,:,:), pointer :: xFace, xInt
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall

!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the domains.

       domain: do nn=1, nDom

         ! Set the pointers for this block.

         call setPointers(nn, level, sps)

         ! Loop over the viscous subfaces of this block. Note that
         ! these are numbered first.

         bocos: do mm=1,nViscBocos

           ! Set the pointers for the plane on the surface, one plane
           ! into the computational domain and the wall distance.
           ! This depends on the block face on which the subface is
           ! located. Note that the starting index of d2Wall is 2 and
           ! therefore a pointer offset will be needed later on.

           select case (BCFaceID(mm))

             case (iMin)
               xFace   => x(1, 1:,1:,:); xInt => x(2, 1:,1:,:)
               dd2Wall => d2Wall(2, :,:)

             case (iMax)
               xFace   => x(il,1:,1:,:); xInt => x(nx,1:,1:,:)
               dd2Wall => d2Wall(il,:,:)

             case (jMin)
               xFace   => x(1:,1, 1:,:); xInt => x(1:,2, 1:,:)
               dd2Wall => d2Wall(:,2 ,:)

             case (jMax)
               xFace   => x(1:,jl,1:,:); xInt => x(1:,ny,1:,:)
               dd2Wall => d2Wall(:,jl,:)

             case (kMin)
               xFace   => x(1:,1:,1, :); xInt => x(1:,1:,2 ,:)
               dd2Wall => d2Wall(:,:,2 )

             case (kMax)
               xFace   => x(1:,1:,kl,:); xInt => x(1:,1:,nz,:)
               dd2Wall => d2Wall(:,:,kl)

           end select

           ! Store the face range of this subface a bit easier.

           jBeg = BCData(mm)%jnBeg+1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg+1; iEnd = BCData(mm)%inEnd

           ! Loop over the faces of the subfaces.

           do j=jBeg,jEnd
             do i=iBeg,iEnd

               ! Store the three components of the unit normal a
               ! bit easier.

               nnx = BCData(mm)%norm(i,j,1)
               nny = BCData(mm)%norm(i,j,2)
               nnz = BCData(mm)%norm(i,j,3)

               ! Compute the vector from centroid of the adjacent cell
               ! to the centroid of the face.

               vecx = eighth*(xFace(i-1,j-1,1) + xFace(i-1,j,1) &
                    +         xFace(i,  j-1,1) + xFace(i,  j,1) &
                    -          xInt(i-1,j-1,1) -  xInt(i-1,j,1) &
                    -          xInt(i,  j-1,1) -  xInt(i,  j,1))

               vecy = eighth*(xFace(i-1,j-1,2) + xFace(i-1,j,2) &
                    +         xFace(i,  j-1,2) + xFace(i,  j,2) &
                    -          xInt(i-1,j-1,2) -  xInt(i-1,j,2) &
                    -          xInt(i,  j-1,2) -  xInt(i,  j,2))

               vecz = eighth*(xFace(i-1,j-1,3) + xFace(i-1,j,3) &
                    +         xFace(i,  j-1,3) + xFace(i,  j,3) &
                    -          xInt(i-1,j-1,3) -  xInt(i-1,j,3) &
                    -          xInt(i,  j-1,3) -  xInt(i,  j,3))

               ! Compute the projection of this vector onto the normal
               ! vector of the face. For a decent mesh there will not be
               ! much of a difference between the projection and the
               ! original mesh, but it does not hurt to do it.

               dot = nnx*vecx + nny*vecy + nnz*vecz

               ! As (nnx,nny,nnz) is a unit vector the distance to the
               ! wall of the first cell center is given by the absolute
               ! value of dot. Due to the use of pointers and the fact
               ! that the original d2Wall array starts at 2 and offset
               ! of -1 must be used to store the data at the correct
               ! location.

               dd2Wall(i-1,j-1) = abs(dot)

             enddo
           enddo

         enddo bocos

       enddo domain

       end subroutine computeNormalSpacing

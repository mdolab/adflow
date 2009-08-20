!
!      ******************************************************************
!      *                                                                *
!      * File:          xhaloAdj.f90                                    *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 08-13-2009                                      *
!      * Last modified: 08-13-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine xhaloAdj(xAdj,xBlockCornerAdj,icell,jcell,kcell)
!
!      ******************************************************************
!      *                                                                *
!      * xhaloAdj determines the coordinates of the nodal halo's.       *
!      * First it sets all halo coordinates by simple extrapolation,    *
!      * then the symmetry planes are treated (also the unit normal of  *
!      * symmetry planes are determined) the internal halo exchange is  *
!      * ignored. The global indexing in the derivatives takes care of  *
!      * this.                                                          *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iCell, jCell, kCell
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3), &
            intent(out) :: xAdj
!integer(kind=intType) :: level
        real(kind=realType), dimension(2,2,2,3) :: xBlockCornerAdj
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, sps, i, j, k,ii,jj
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iiMax, jjMax
       !integer(kind=intType) :: iBeg2, iEnd2, jBeg2, jEnd2

       integer(kind=intType) :: iOffset, jOffset, kOffset
       logical :: iMinOverlap, jMinOverlap, kMinOverlap 
       logical :: iMaxOverlap, jMaxOverlap, kMaxOverlap

       real(kind=realType), dimension(-3:2,-3:2,1:3) :: x0, x1, x2
       !real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2

       real(kind=realType) :: length, dot

       real(kind=realType), dimension(3) :: v1, v2, norm
       real(kind=realType), dimension(2,2,3) :: xFaceCorner
       integer(kind=intType) ::istart,jstart!,iend,jend

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!!$       ! Loop over the number of spectral solutions and the local
!!$       ! number of blocks.
!!$
!!$       spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$         domains: do nn=1,nDom
!!$
!!$           ! Set the pointers to this block.
!!$
!!$           call setPointers(nn, level, sps)
!!$

       !Check to see if current cell contains a halo and if so which one...
       call checkXOverlapAdj(icell,jcell,kcell,ioffset,joffset,&
            koffset,iMinOverlap, jMinOverlap, kMinOverlap ,&
            iMaxOverlap, jMaxOverlap, kMaxOverlap)
!
!          **************************************************************
!          *                                                            *
!          * Extrapolation of the coordinates. First extrapolation in   *
!          * i-direction, without halo's, followed by extrapolation in  *
!          * j-direction, with i-halo's and finally extrapolation in    *
!          * k-direction, with both i- and j-halo's. In this way also   *
!          * the indirect halo's get a value, albeit a bit arbitrary.   *
!          *                                                            *
!          **************************************************************
!
           ! Extrapolation in i-direction.
       if (iMinOverlap)then
          if(ioffset>=-3)then
             do j = -3,2
                do k = -3,2
                   xAdj(ioffset,j,k,1) = two*xAdj(ioffset+1,j,k,1) - xAdj(ioffset+2,j,k,1)
                   xAdj(ioffset,j,k,2) = two*xAdj(ioffset+1,j,k,2) - xAdj(ioffset+2,j,k,2)
                   xAdj(ioffset,j,k,3) = two*xAdj(ioffset+1,j,k,3) - xAdj(ioffset+2,j,k,3)
                end do
             end do
          end if
       end if
       if (iMaxOverlap)then
          if(ioffset<=2)then
             do j = -3,2
                do k = -3,2
                   xAdj(ioffset,j,k,1) = two*xAdj(ioffset-1,j,k,1) - xAdj(ioffset-2,j,k,1)
                   xAdj(ioffset,j,k,2) = two*xAdj(ioffset-1,j,k,2) - xAdj(ioffset-2,j,k,2)
                   xAdj(ioffset,j,k,3) = two*xAdj(ioffset-1,j,k,3) - xAdj(ioffset-2,j,k,3)
                end do
             end do
          end if
       end if
!!$           do k=1,kl
!!$             do j=1,jl
!!$               x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
!!$               x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
!!$               x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)
!!$
!!$               x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
!!$               x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
!!$               x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
!!$             enddo
!!$           enddo

           ! Extrapolation in j-direction.
       if (jMinOverlap)then
          if(joffset>=-3)then
             do i = -3,2
                do k = -3,2
                   xAdj(i,joffset,k,1) = two*xAdj(i,joffset+1,k,1) - xAdj(i,joffset+2,k,1)
                   xAdj(i,joffset,k,2) = two*xAdj(i,joffset+1,k,2) - xAdj(i,joffset+2,k,2)
                   xAdj(i,joffset,k,3) = two*xAdj(i,joffset+1,k,3) - xAdj(i,joffset+2,k,3)
                end do
             end do
          end if
       end if
       if (jMaxOverlap)then
          if(joffset<=2)then
             do k = -3,2
                do i = -3,2
                   xAdj(i,joffset,k,1) = two*xAdj(i,joffset-1,k,1) - xAdj(i,joffset-2,k,1)
                   xAdj(i,joffset,k,2) = two*xAdj(i,joffset-1,k,2) - xAdj(i,joffset-2,k,2)
                   xAdj(i,joffset,k,3) = two*xAdj(i,joffset-1,k,3) - xAdj(i,joffset-2,k,3)
                end do
             end do
          end if
       end if
!!$           do k=1,kl
!!$             do i=0,ie
!!$               x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
!!$               x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
!!$               x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)
!!$
!!$               x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
!!$               x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
!!$               x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
!!$             enddo
!!$           enddo

           ! Extrapolation in k-direction.
      if (kMinOverlap)then
          if(koffset>=-3)then
             do j = -3,2
                do i = -3,2
                   xAdj(i,j,koffset,1) = two*xAdj(i,j,koffset+1,1) - xAdj(i,j,koffset+2,1)
                   xAdj(i,j,koffset,2) = two*xAdj(i,j,koffset+1,2) - xAdj(i,j,koffset+2,2)
                   xAdj(i,j,koffset,3) = two*xAdj(i,j,koffset+1,3) - xAdj(i,j,koffset+2,3)
                end do
             end do
          end if
       end if
       if (kMaxOverlap)then
          if(koffset<=2)then
             do j = -3,2
                do i = -3,2
                   xAdj(i,j,koffset,1) = two*xAdj(i,j,koffset-1,1) - xAdj(i,j,koffset-2,1)
                   xAdj(i,j,koffset,2) = two*xAdj(i,j,koffset-1,2) - xAdj(i,j,koffset-2,2)
                   xAdj(i,j,koffset,3) = two*xAdj(i,j,koffset-1,3) - xAdj(i,j,koffset-2,3)
                end do
             end do
          end if
       end if
!!$           do j=0,je
!!$             do i=0,ie
!!$               x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
!!$               x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
!!$               x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)
!!$
!!$               x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
!!$               x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
!!$               x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
!!$             enddo
!!$           enddo
!
!          **************************************************************
!          *                                                            *
!          * Mirror the halo coordinates adjacent to the symmetry       *
!          * planes                                                     *
!          *                                                            *
!          **************************************************************
!
           ! Loop over boundary subfaces.

           loopBocos: do mm=1,nBocos

             ! The actual correction of the coordinates only takes
             ! place for symmetry planes.

             testSymmetry: if(BCType(mm) == Symm) then

!!$               ! Set some variables, depending on the block face on
!!$               ! which the subface is located.
!!$
!!$               select case (BCFaceID(mm))
!!$                 case (iMin)
!!$                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
!!$                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
!!$                   x0 => x(0,:,:,:); x1 => x(1,:,:,:); x2 => x(2,:,:,:)
!!$
!!$                 case (iMax)
!!$                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
!!$                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
!!$                   x0 => x(ie,:,:,:); x1 => x(il,:,:,:); x2 => x(nx,:,:,:)
!!$
!!$                 case (jMin)
!!$                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
!!$                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
!!$                   x0 => x(:,0,:,:); x1 => x(:,1,:,:); x2 => x(:,2,:,:)
!!$
!!$                 case (jMax)
!!$                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
!!$                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
!!$                   x0 => x(:,je,:,:); x1 => x(:,jl,:,:); x2 => x(:,ny,:,:)
!!$
!!$                 case (kMin)
!!$                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
!!$                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
!!$                   x0 => x(:,:,0,:); x1 => x(:,:,1,:); x2 => x(:,:,2,:)
!!$
!!$                 case (kMax)
!!$                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
!!$                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
!!$                   x0 => x(:,:,ke,:); x1 => x(:,:,kl,:); x2 => x(:,:,nz,:)
!!$               end select

              ! Set some variables, depending on the block face on
               ! which the subface is located.
                istart = -3;iend = 2
                jstart = -3;jend = 2

               select case (BCFaceID(mm))
                 case (iMin)
                    iiMax = jl
                    jjMax = kl
                    !iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                    !jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                    !x0 => x(0,:,:,:); x1 => x(1,:,:,:); x2 => x(2,:,:,:)
                    istart = max(-3,(2-jcell)-2)
                    jstart = max(-3,(2-kcell)-2)
                    iend = min(2,(jl-jcell)+1)
                    jend = min(2,(kl-kcell)+1)
                    !print *,'bcfaceID,imin',BCFaceID(mm),imin,iMinOverlap,icell,jcell,kcell
                    if(iMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(1,:,:,:)
                       if (icell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-2,-3:2,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3,-3:2,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-2,-3:2,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling imin'
                       cycle loopBocos
                    endif

                 case (iMax)
                    iiMax = jl
                    jjMax = kl
                    !iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                    !jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                    !x0 => x(ie,:,:,:); x1 => x(il,:,:,:); x2 => x(nx,:,:,:)
                    istart = max(-3,(2-jcell)-2)
                    jstart = max(-3,(2-kcell)-2)
                    iend = min(2,(jl-jcell)+1)
                    jend = min(2,(kl-kcell)+1)
                    !print *,'bcfaceID,imax',BCFaceID(mm),imax,iMaxOverlap,icell,jcell,kcell
                    if(iMaxOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(2,:,:,:)
                       if(icell==il)then
                          x0(-3:2,-3:2,1:3) = xAdj(1,-3:2,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(2,-3:2,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(1,-3:2,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3)
                       endif
                    else
                       !skip remainder of loop
                       !print *,'cycling imax'
                       cycle loopBocos
                    endif

                 case (jMin)
                    iiMax = il
                    jjMax = kl
                   !iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   !jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   !x0 => x(:,0,:,:); x1 => x(:,1,:,:); x2 => x(:,2,:,:)
                    istart = max(-3,(2-icell)-2)
                    jstart = max(-3,(2-kcell)-2)
                    iend = min(2,(il-icell)+1)
                    jend = min(2,(kl-kcell)+1)
                    !print *,'bcfaceID,jmin',BCFaceID(mm),jmin,jMinOverlap,icell,jcell,kcell
                    if(jMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,1,:,:)
                       if (jcell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-2,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3,-3:2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-2,-3:2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling jmin'
                       cycle loopBocos
                    endif

                 case (jMax)
                    iiMax = il
                    jjMax = kl
                    !iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   !jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   !x0 => x(:,je,:,:); x1 => x(:,jl,:,:); x2 => x(:,ny,:,:)
                    istart = max(-3,(2-icell)-2)
                    jstart = max(-3,(2-kcell)-2)
                    iend = min(2,(il-icell)+1)
                    jend = min(2,(kl-kcell)+1)
                    !print *,'bcfaceID,jmax',BCFaceID(mm),jmax,jMaxOverlap,icell,jcell,kcell
                   if(jMaxOverlap)then
                      xFaceCorner(:,:,:) = xBlockCornerAdj(:,2,:,:)
                      if(jcell == jl)then
                         x0(-3:2,-3:2,1:3) = xAdj(-3:2,1,-3:2,1:3)
                         x1(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3)
                         x2(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3)
                      else
                         x0(-3:2,-3:2,1:3) = xAdj(-3:2,2,-3:2,1:3)
                         x1(-3:2,-3:2,1:3) = xAdj(-3:2,1,-3:2,1:3)
                         x2(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3)
                      endif
                   else
                      !skip remainder of loop
                      !print *,'cycling jmax'
                      cycle loopBocos
                   endif
                 case (kMin)
                    iiMax = il
                    jjMax = jl
                   !iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   !jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   !x0 => x(:,:,0,:); x1 => x(:,:,1,:); x2 => x(:,:,2,:)
                    istart = max(-3,(2-icell)-2)
                    jstart = max(-3,(2-jcell)-2)
                    iend = min(2,(il-icell)+1)
                    jend = min(2,(jl-jcell)+1)  
                    !print *,'bcfaceID,kmin',BCFaceID(mm),kmin,kMinOverlap,icell,jcell,kcell
                    if(kMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,:,1,:)
                       if (kcell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-3,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-2,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling kmin'
                       cycle loopBocos
                    endif
                 case (kMax)
                    iiMax = il
                    jjMax = jl
                   ! iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   ! jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   ! x0 => x(:,:,ke,:); x1 => x(:,:,kl,:); x2 => x(:,:,nz,:)
                    istart = max(-3,(2-icell)-2)
                    jstart = max(-3,(2-jcell)-2)
                    iend = min(2,(il-icell)+1)
                    jend = min(2,(jl-jcell)+1)            
                    !print *,'bcfaceID,kmax',BCFaceID(mm),kmax,kMaxOverlap,icell,jcell,kcell
                    if(kMaxOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,:,2,:)
                       if (kcell==kl)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,1,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,2,1:3)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,1,1:3)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3)
                       endif
                    else
                       !skip remainder of loop
                       !print *,'cycling kmax'
                       cycle loopBocos
                    endif
                    
               end select

               
!CAM commented out on Aug. 11, 2009 by C.A.Mader. Computation modified to fit 
!CAM into a single residual stencil...
!!$
!!$               ! Determine the sum of all face normals and store this
!!$               ! in norm. The sum of all faces is taken instead of the
!!$               ! cross product of the diagonals, because for some c-type
!!$               ! grids the diagonals of the subface are aligned.
!!$
!!$               norm = zero
!!$
!!$               do j=(jBeg+1),jEnd
!!$                 do i=(iBeg+1),iEnd
!!$
!!$                   ! Determine the vector from the lower left corner to
!!$                   ! the upper right corner. Due to the usage of pointers
!!$                   ! an offset of +1 must be used, because the original
!!$                   ! array x start at 0.
!!$
!!$                   v1(1) = x1(i+1,j+1,1) - x1(i,j,1)
!!$                   v1(2) = x1(i+1,j+1,2) - x1(i,j,2)
!!$                   v1(3) = x1(i+1,j+1,3) - x1(i,j,3)
!!$                   
!!$                   ! And the vector from the upper left corner to the
!!$                   ! lower right corner.
!!$
!!$                   v2(1) = x1(i+1,j,1) - x1(i,j+1,1)
!!$                   v2(2) = x1(i+1,j,2) - x1(i,j+1,2)
!!$                   v2(3) = x1(i+1,j,3) - x1(i,j+1,3)
!!$                   
!!$                   ! Determine the normal of the face by taking the cross
!!$                   ! product of v1 and v2 and add it to norm.
!!$
!!$                   norm(1) = norm(1) + v1(2)*v2(3) - v1(3)*v2(2)
!!$                   norm(2) = norm(2) + v1(3)*v2(1) - v1(1)*v2(3)
!!$                   norm(3) = norm(3) + v1(1)*v2(2) - v1(2)*v2(1)
!!$                   print *,'norm', x1(i+1,j+1,1) ,x1(i,j,1) ,norm(1),norm(2),norm(3),i,j,BCFaceID(mm)
!!$                 enddo
!!$               enddo
!!$
!!$               ! Compute the length of the normal and test if this is
!!$               ! larger than eps. If this is the case this means that
!!$               ! it is a nonsingular subface and the coordinates are
!!$               ! corrected.
!!$
!!$               length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
!!$
!!$               testSingular: if(length > eps) then
!!$
!!$                 ! Compute the unit normal of the subface.
!!$
!!$                 norm(1) = norm(1)/length
!!$                 norm(2) = norm(2)/length
!!$                 norm(3) = norm(3)/length
!!$
!!$                 ! Add an overlap to the symmetry subface if the
!!$                 ! boundaries coincide with the block boundaries.
!!$                 ! This way the indirect halo's are treated properly.
!!$
!!$                 if(iBeg == 1)     iBeg = 0
!!$                 if(iEnd == iiMax) iEnd = iiMax + 1
!!$
!!$                 if(jBeg == 1)     jBeg = 0
!!$                 if(jEnd == jjMax) jEnd = jjMax + 1
!!$
!!$                 ! Loop over the nodes of the subface and set the
!!$                 ! corresponding halo coordinates.
!!$
!!$                 do j=jBeg,jEnd
!!$                   do i=iBeg,iEnd
!!$
!!$                     ! Determine the vector from the internal node to the
!!$                     ! node on the face. Again an offset of +1 must be
!!$                     ! used, due to the usage of pointers.
!!$
!!$                     v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
!!$                     v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
!!$                     v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)
!!$
!!$                     ! Determine two times the normal component of this
!!$                     ! vector; this vector must be added to the
!!$                     ! coordinates of the internal node to obtain the
!!$                     ! halo coordinates. Again the offset of +1.
!!$
!!$                     dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
!!$                         +      v1(3)*norm(3))
!!$
!!$                     x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
!!$                     x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
!!$                     x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)
!!$                     print *,'xhalo', x0(i+1,j+1,1) ,x2(i+1,j+1,1) , dot,norm(1),i,j,BCFaceID(mm)
!!$
!!$                   enddo
!!$                 enddo
!!$
!!$                endif testSingular
               
               !compute a norm from the 4 corners of the subface. This
               ! will serve to check the subface for singularity and
               !reduce the number of points to be dealt with in the derivatives

               ! Determine the vector from the lower left corner to
               ! the upper right corner. Due to the usage of pointers
               ! an offset of +1 must be used, because the original
               ! array x start at 0.
               
               v1(1) = xFaceCorner(2,2,1) - xFaceCorner(1,1,1)
               v1(2) = xFaceCorner(2,2,2) - xFaceCorner(1,1,2)
               v1(3) = xFaceCorner(2,2,3) - xFaceCorner(1,1,3)
               !print *,'v1adj',v1,xFaceCorner(2,2,1), xFaceCorner(1,1,1)
               ! And the vector from the upper left corner to the
               ! lower right corner.
               
               v2(1) = xFaceCorner(2,1,1) - xFaceCorner(1,2,1)
               v2(2) = xFaceCorner(2,1,2) - xFaceCorner(1,2,2)
               v2(3) = xFaceCorner(2,1,3) - xFaceCorner(1,2,3)
!               print *,'v2adj',v2    
               ! Determine the normal of the face by taking the cross
               ! product of v1 and v2 and add it to norm.
                           
               norm(1) = v1(2)*v2(3) - v1(3)*v2(2)
               norm(2) = v1(3)*v2(1) - v1(1)*v2(3)
               norm(3) = v1(1)*v2(2) - v1(2)*v2(1)
               !print *,'normadj',norm
               ! Compute the length of the normal and test if this is
               ! larger than eps. If this is the case this means that
               ! it is a nonsingular subface and the coordinates are
               ! corrected.
               
               length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
               
               testSingular: if(length > eps) then
                  
                  ! Compute the unit normal of the subface.
                  
                  norm(1) = norm(1)/length
                  norm(2) = norm(2)/length
                  norm(3) = norm(3)/length
                  
                  
!!$                  ! Add an overlap to the symmetry subface if the
!!$                  ! boundaries coincide with the block boundaries.
!!$                  ! This way the indirect halo's are treated properly.
!!$                  
!!$                  if(iBeg == 1)     iBeg = 0
!!$                  if(iEnd == iiMax) iEnd = iiMax + 1
!!$                  
!!$                  if(jBeg == 1)     jBeg = 0
!!$                  if(jEnd == jjMax) jEnd = jjMax + 1
                  
                  ! Loop over the nodes of the subface and set the
                  ! corresponding halo coordinates.
                  
                  do j=jstart,jend!jBeg,jEnd
                     do i=istart,iend!iBeg,iEnd
                       
                        ! Determine the vector from the internal node to the
                        ! node on the face. Again an offset of +1 must be
                        ! used, due to the usage of pointers.
                        
!!$                        v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
!!$                        v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
!!$                        v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)
                        v1(1) = x1(i,j,1) - x2(i,j,1)
                        v1(2) = x1(i,j,2) - x2(i,j,2)
                        v1(3) = x1(i,j,3) - x2(i,j,3)
                        
                        ! Determine two times the normal component of this
                        ! vector; this vector must be added to the
                        ! coordinates of the internal node to obtain the
                        ! halo coordinates. Again the offset of +1.
                        
                        dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                             +      v1(3)*norm(3))
                        
                        !x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
                        !x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
                        !x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)
                        x0(i,j,1) = x2(i,j,1) + dot*norm(1)
                        x0(i,j,2) = x2(i,j,2) + dot*norm(2)
                        x0(i,j,3) = x2(i,j,3) + dot*norm(3)
                        !print *,'xhaloAdj', x0(i,j,1) ,x2(i,j,1) , dot,norm(1),icell+i,jcell+j,BCFaceID(mm)
                     enddo
                  enddo
               endif testSingular
               
               ! Reset modified xAdj, depending on the block face on
               ! which the subface is located.
                
               select case (BCFaceID(mm))
                 case (iMin)                   
                    if(iMinOverlap)then
                       if (icell==2)then
                          xAdj(-2,-3:2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3,-3:2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-2,-3:2,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                       end if
                    endif                  

                 case (iMax)                 
                    if(iMaxOverlap)then
                       if (icell==il)then
                          xAdj(1,-3:2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(2,-3:2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(1,-3:2,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                       endif
                    endif

                 case (jMin)
                    if(jMinOverlap)then
                       if (jcell==2)then
                          xAdj(-3:2,-2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-1,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,0,-3:2,1:3)=x2(-3:2,-3:2,1:3) 
                       else
                          xAdj(-3:2,-3,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-2,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-1,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                       end if
                    endif
                 case (jMax)
                   if(jMaxOverlap)then
                      if (jcell==jl)then
                         xAdj(-3:2,1,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                         xAdj(-3:2,0,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                         xAdj(-3:2,-1,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                      else
                         xAdj(-3:2,2,-3:2,1:3)=x0(-3:2,-3:2,1:3)
                         xAdj(-3:2,1,-3:2,1:3)=x1(-3:2,-3:2,1:3)
                         xAdj(-3:2,0,-3:2,1:3)=x2(-3:2,-3:2,1:3)
                      endif
                   endif
                 case (kMin)
                    if(kMinOverlap)then
                       if (kcell==2)then
                          xAdj(-3:2,-3:2,-2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3:2,-3:2,-3,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-2,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3)=x2(-3:2,-3:2,1:3)
                       end if
                    endif
                 case (kMax)
                    if(kMaxOverlap)then
                       if (kcell==kl)then
                          xAdj(-3:2,-3:2,1,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3:2,-3:2,2,1:3)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,1,1:3)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3)=x2(-3:2,-3:2,1:3)
                       endif
                    endif
               end select 
            endif testSymmetry
         enddo loopBocos
 !     enddo domains
 !  enddo spectralLoop
!!$!
!!$!      ******************************************************************
!!$!      *                                                                *
!!$!      * Exchange the coordinates for the internal halo's.              *
!!$!      *                                                                *
!!$!      ******************************************************************
!!$!
!!$   call exchangeCoor(level)
   
       end subroutine xhaloAdj

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
       subroutine xhaloAdjTS(xAdj,xBlockCornerAdj,icell,jcell,kcell,nn,level,sps)
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
       use inputTimeSpectral !nIntervalTimespectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iCell, jCell, kCell
       integer(kind=intType), intent(in) :: nn,level,sps
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), &
            intent(out) :: xAdj
!integer(kind=intType) :: level
        real(kind=realType), dimension(2,2,2,3,nTimeIntervalsSpectral) :: xBlockCornerAdj
!
!      Local variables.
!
       integer(kind=intType) ::  mm, i, j, k,ii,jj
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iiMax, jjMax
       !integer(kind=intType) :: iBeg2, iEnd2, jBeg2, jEnd2

       integer(kind=intType) :: iminOffset, jminOffset, kminOffset
       integer(kind=intType) :: imaxOffset, jmaxOffset, kmaxOffset

       logical :: iMinOverlap, jMinOverlap, kMinOverlap 
       logical :: iMaxOverlap, jMaxOverlap, kMaxOverlap
       logical :: iMinInternal, jMinInternal, kMinInternal 
       logical :: iMaxInternal, jMaxInternal, kMaxInternal

       real(kind=realType), dimension(-3:2,-3:2,1:3) :: x0, x1, x2
       !real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2

       real(kind=realType) :: length, dot!,dimension(nTimeIntervalsSpectral)

       real(kind=realType), dimension(3) :: v1, v2, norm
       real(kind=realType), dimension(2,2,3) :: xFaceCorner
       integer(kind=intType) ::istart,jstart!,iend,jend

       !FILEIO
       integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart2,jstart2,kstart2,iend2,jend2,kend2,n,nnn,itemp,jtemp
       integer ::unitxAD = 13

     
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
       call checkXOverlapAdjTS(icell,jcell,kcell,iminoffset,imaxoffset,&
            jminoffset,jmaxoffset,kminoffset,&
            kmaxoffset,iMinOverlap, jMinOverlap, kMinOverlap ,&
            iMaxOverlap, jMaxOverlap, kMaxOverlap)

!
!          **************************************************************
!          *                                                            *
!          * Because of exchangeCorr at the end of the original xhalo   *
!          * the internal face halos need to be left alone. Therefore   *
!          * loop over the Bocos to determine which block faces are     *
!          * internal.                                                  *
!          *                                                            *
!          **************************************************************
!
       iMinInternal=.false.
       jMinInternal=.false.
       kMinInternal=.false.
       iMaxInternal=.false.
       jMaxInternal=.false.
       kMaxInternal=.false.
       !print *,'internals',iMinInternal,jMinInternal,kMinInternal,&
        !    iMaxInternal,jMaxInternal,kMaxInternal
       loopSubface: do mm=1,nSubface!Bocos
          ! Check whether a given block face is internal
         ! print *,'bctype',BCType(mm),B2BMatch,mm
          testInternal: if(BCType(mm) ==B2BMatch)then
             select case (BCFaceID(mm))
             case (iMin)
                iMinInternal =.true.
             case (iMax)
                iMaxInternal =.true.
             case (jMin)
                jMinInternal =.true.
             case (jMax)
                jMaxInternal =.true. 
             case (kMin)
                kMinInternal =.true.
             case (kMax)
                kMaxInternal =.true.
             end select
          end if testInternal
       end do loopSubface
       !print *,'internals after',iMinInternal,jMinInternal,kMinInternal,&
       !     iMaxInternal,jMaxInternal,kMaxInternal
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
          if(.not.(iMinInternal ))then
             if(iminoffset>=-3)then
                do j = -3,2
                   do k = -3,2
                      jjj = jcell+j
                      kkk = kcell+k
                      if(kkk>=1 .and.jjj>=1)then
                         if(kkk<=kl.and.jjj<=jl)then
                            xAdj(iminoffset,j,k,1,sps) = two*xAdj(iminoffset+1,j,k,1,sps) - xAdj(iminoffset+2,j,k,1,sps)
                            xAdj(iminoffset,j,k,2,sps) = two*xAdj(iminoffset+1,j,k,2,sps) - xAdj(iminoffset+2,j,k,2,sps)
                            xAdj(iminoffset,j,k,3,sps) = two*xAdj(iminoffset+1,j,k,3,sps) - xAdj(iminoffset+2,j,k,3,sps)
                         endif
                      endif
                   end do
                end do
             end if
          end if
       endif
       if (iMaxOverlap)then
          if(.not.(iMaxInternal ))then
             if(imaxoffset<=2)then
                do j = -3,2
                   do k = -3,2
                      jjj = jcell+j
                      kkk = kcell+k
                      if(kkk>=1 .and.jjj>=1)then
                         if(kkk<=kl.and.jjj<=jl)then
                            xAdj(imaxoffset,j,k,1,sps) = two*xAdj(imaxoffset-1,j,k,1,sps) - xAdj(imaxoffset-2,j,k,1,sps)
                            xAdj(imaxoffset,j,k,2,sps) = two*xAdj(imaxoffset-1,j,k,2,sps) - xAdj(imaxoffset-2,j,k,2,sps)
                            xAdj(imaxoffset,j,k,3,sps) = two*xAdj(imaxoffset-1,j,k,3,sps) - xAdj(imaxoffset-2,j,k,3,sps)
                         endif
                      end if
                   end do
                end do
             end if
          end if
       endif
!!$!           do k=1,kl
!!$!             do j=1,jl
!!$!               x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
!!$!               x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
!!$!               x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)
!!$!
!!$!               x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
!!$!               x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
!!$!               x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
!!$!             enddo
!!$!           enddo

           ! Extrapolation in j-direction.
       if (jMinOverlap)then
          if(.not.(jMinInternal ))then
             if(jminoffset>=-3)then
                do i = -3,2
                   do k = -3,2 
                      iii = icell+i
                      kkk = kcell+k
                      if(kkk>=1 .and.iii>=0)then
                         if(kkk<=kl.and.iii<=ie)then
                            xAdj(i,jminoffset,k,1,sps) = two*xAdj(i,jminoffset+1,k,1,sps) - xAdj(i,jminoffset+2,k,1,sps)
                            xAdj(i,jminoffset,k,2,sps) = two*xAdj(i,jminoffset+1,k,2,sps) - xAdj(i,jminoffset+2,k,2,sps)
                            xAdj(i,jminoffset,k,3,sps) = two*xAdj(i,jminoffset+1,k,3,sps) - xAdj(i,jminoffset+2,k,3,sps)
                         endif
                      endif
                   end do
                end do
             end if
          end if
       endif
       if (jMaxOverlap)then
          if(.not.(jMaxInternal ))then
             if(jmaxoffset<=2)then
                do k = -3,2
                   do i = -3,2
                      iii = icell+i
                      kkk = kcell+k
                      if(kkk>=1 .and.iii>=0)then
                         if(kkk<=kl.and.iii<=ie)then
                            xAdj(i,jmaxoffset,k,1,sps) = two*xAdj(i,jmaxoffset-1,k,1,sps) - xAdj(i,jmaxoffset-2,k,1,sps)
                            xAdj(i,jmaxoffset,k,2,sps) = two*xAdj(i,jmaxoffset-1,k,2,sps) - xAdj(i,jmaxoffset-2,k,2,sps)
                            xAdj(i,jmaxoffset,k,3,sps) = two*xAdj(i,jmaxoffset-1,k,3,sps) - xAdj(i,jmaxoffset-2,k,3,sps)
                         endif
                      endif
                   end do
                end do
             end if
          end if
       endif
!           do k=1,kl
!             do i=0,ie
!               x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
!               x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
!               x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)
!
!               x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
!               x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
!               x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
!             enddo
!           enddo

           ! Extrapolation in k-direction.
       if (kMinOverlap)then
          if(.not.(kMinInternal ))then
             !print *,'kcellmin',kcell,kminoffset
             if(kminoffset>=-3)then
                do j = -3,2
                   do i = -3,2 
                      iii = icell+i
                      jjj = jcell+j
                      if(jjj>=0 .and.iii>=0)then
                         if(jjj<=je.and.iii<=ie)then
                            xAdj(i,j,kminoffset,1,sps) = two*xAdj(i,j,kminoffset+1,1,sps) - xAdj(i,j,kminoffset+2,1,sps)
                            xAdj(i,j,kminoffset,2,sps) = two*xAdj(i,j,kminoffset+1,2,sps) - xAdj(i,j,kminoffset+2,2,sps)
                            xAdj(i,j,kminoffset,3,sps) = two*xAdj(i,j,kminoffset+1,3,sps) - xAdj(i,j,kminoffset+2,3,sps)
                         endif
                      endif
                   end do
                end do
             end if
          end if
       endif
       if (kMaxOverlap)then
          if(.not.(kMaxInternal ))then
             !print *,'kcellmax',kcell,kmaxoffset
             if(kmaxoffset<=2)then
                do j = -3,2
                   do i = -3,2
                      iii = icell+i
                      jjj = jcell+j
                      if(jjj>=0 .and.iii>=0)then
                         if(jjj<=je.and.iii<=ie)then
                            xAdj(i,j,kmaxoffset,1,sps) = two*xAdj(i,j,kmaxoffset-1,1,sps) - xAdj(i,j,kmaxoffset-2,1,sps)
                            xAdj(i,j,kmaxoffset,2,sps) = two*xAdj(i,j,kmaxoffset-1,2,sps) - xAdj(i,j,kmaxoffset-2,2,sps)
                            xAdj(i,j,kmaxoffset,3,sps) = two*xAdj(i,j,kmaxoffset-1,3,sps) - xAdj(i,j,kmaxoffset-2,3,sps)     
                         endif
                      endif
                   end do
                end do
             end if
          end if
       endif
!           do j=0,je
!             do i=0,ie
!               x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
!               x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
!               x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)!!
!
!               x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
!               x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
!               x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
!             enddo
!           enddo
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
              !write(unitxAD,*)'loopBocos',mm,nbocos
             ! The actual correction of the coordinates only takes
             ! place for symmetry planes.

             testSymmetry: if(BCType(mm) == Symm) then
                !write(unitxAD,*)'testSymmetry',bcfaceID(mm)
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
                    itemp =jcell
                    jtemp =kcell
                    !print *,'bcfaceID,imin',BCFaceID(mm),imin,iMinOverlap,icell,jcell,kcell
                    if(iMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(1,:,:,:,sps)
                       if (icell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-2,-3:2,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3,sps)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3,-3:2,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-2,-3:2,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3,sps)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling imin'
                       !write(unitxAD,*)'cycling imin'
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
                    itemp =jcell
                    jtemp =kcell
                    if(iMaxOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(2,:,:,:,sps)
                       if(icell==il)then
                          x0(-3:2,-3:2,1:3) = xAdj(1,-3:2,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-1,-3:2,-3:2,1:3,sps)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(2,-3:2,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(1,-3:2,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(0,-3:2,-3:2,1:3,sps)
                       endif
                    else
                       !skip remainder of loop
                       !print *,'cycling imax'
                       !write(unitxAD,*)'cycling imax'
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
                    itemp =icell
                    jtemp =kcell
                    if(jMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,1,:,:,sps)
                       if (jcell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-2,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3,sps)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3,-3:2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-2,-3:2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3,sps)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling jmin'
                       !write(unitxAD,*)'cycling jmin'
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
                    itemp =icell
                    jtemp =kcell
                   if(jMaxOverlap)then
                      xFaceCorner(:,:,:) = xBlockCornerAdj(:,2,:,:,sps)
                      if(jcell == jl)then
                         x0(-3:2,-3:2,1:3) = xAdj(-3:2,1,-3:2,1:3,sps)
                         x1(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3,sps)
                         x2(-3:2,-3:2,1:3) = xAdj(-3:2,-1,-3:2,1:3,sps)
                      else
                         x0(-3:2,-3:2,1:3) = xAdj(-3:2,2,-3:2,1:3,sps)
                         x1(-3:2,-3:2,1:3) = xAdj(-3:2,1,-3:2,1:3,sps)
                         x2(-3:2,-3:2,1:3) = xAdj(-3:2,0,-3:2,1:3,sps)
                      endif
                   else
                      !skip remainder of loop
                      !print *,'cycling jmax'
                      !write(unitxAD,*)'cycling jmax'
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
                    itemp =icell
                    jtemp =jcell
                    if(kMinOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,:,1,:,sps)
                       if (kcell==2)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3,sps)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-3,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-2,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3,sps)
                       end if
                    else
                       !skip remainder of loop
                       !print *,'cycling kmin'
                       !write(unitxAD,*)'cycling kmin'
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
                    itemp =icell
                    jtemp =jcell
                    if(kMaxOverlap)then
                       xFaceCorner(:,:,:) = xBlockCornerAdj(:,:,2,:,sps)
                       if (kcell==kl)then
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,1,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,-1,1:3,sps)
                       else
                          x0(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,2,1:3,sps)
                          x1(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,1,1:3,sps)
                          x2(-3:2,-3:2,1:3) = xAdj(-3:2,-3:2,0,1:3,sps)
                       endif
                    else
                       !skip remainder of loop
                       !print *,'cycling kmax'
                       !write(unitxAD,*)'cycling kmax'
                       cycle loopBocos
                    endif
                    
               end select
!!$               write(unitxAD,*) 'kl',kcell,kl,BCFaceID(mm)
!!$               !print out xAdj
!!$               do j=jstart,jend!jBeg,jEnd
!!$                  do i=istart,iend!iBeg,iEnd
!!$                     do n = 1,3
!!$                        ii = itemp+i
!!$                        jj = jtemp+j
!!$                        
!                           write(unitxAD,10) i,j,k,n,nnn,v1(1),v1(2),v1(3)
!10                         format(1x,'res',5I8,3f20.14)
!!$                        write(unitxAD,10) ii,jj,n,itemp,jtemp,x0(i,j,n),x1(i,j,n),x2(i,j,n)!,dot*norm(n)
!!$10                      format(1x,'res',5I8,3f20.14)
!!$                     enddo
!!$                  enddo
!!$               enddo
          
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

                        !print out xAdj
                        
!!$               istart2 = -3
!!$               jstart2 = -3
!!$               kstart2 = -3
!!$               iend2 = 2
!!$               jend2 = 2
!!$               kend2 = 2
!!$               if(icell==2) istart2=-2
!!$               if(jcell==2) jstart2=-2
!!$               if(kcell==2) kstart2=-2
!!$               if(icell==il) iend2=1
!!$               if(jcell==jl) jend2=1
!!$               if(kcell==kl) kend2=1
               !do iiii = istart2,iend2
               !do jjjj = jstart2,jend2
                !     do kkkk = kstart2,kend2
!!$                        do n = 1,3
!!$                           ii = itemp+i
!!$                           jj = jtemp+j
!!$                           
!                           write(unitxAD,10) i,j,k,n,nnn,v1(1),v1(2),v1(3)
!10                         format(1x,'res',5I8,3f20.14)
!!$                           write(unitxAD,10) ii,jj,n,itemp,jtemp,x0(i,j,n),x2(i,j,n)!,dot*norm(n)
!!$10                         format(1x,'res',5I8,2f20.14)
!!$                        enddo
!
!!$                     enddo
!!$                  enddo
!!$               enddo
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
                          xAdj(-2,-3:2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3,-3:2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-2,-3:2,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                       end if
                    endif                  

                 case (iMax)                 
                    if(iMaxOverlap)then
                       if (icell==il)then
                          xAdj(1,-3:2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-1,-3:2,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(2,-3:2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(1,-3:2,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(0,-3:2,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                       endif
                    endif

                 case (jMin)
                    if(jMinOverlap)then
                       if (jcell==2)then
                          xAdj(-3:2,-2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-1,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,0,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3) 
                       else
                          xAdj(-3:2,-3,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-2,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-1,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                       end if
                    endif
                 case (jMax)
                   if(jMaxOverlap)then
                      if (jcell==jl)then
                         xAdj(-3:2,1,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                         xAdj(-3:2,0,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                         xAdj(-3:2,-1,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                      else
                         xAdj(-3:2,2,-3:2,1:3,sps)=x0(-3:2,-3:2,1:3)
                         xAdj(-3:2,1,-3:2,1:3,sps)=x1(-3:2,-3:2,1:3)
                         xAdj(-3:2,0,-3:2,1:3,sps)=x2(-3:2,-3:2,1:3)
                      endif
                   endif
                 case (kMin)
                    if(kMinOverlap)then
                       if (kcell==2)then
                          xAdj(-3:2,-3:2,-2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3,sps)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3:2,-3:2,-3,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-2,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3,sps)=x2(-3:2,-3:2,1:3)
                       end if
                    endif
                 case (kMax)
                    if(kMaxOverlap)then
                       if (kcell==kl)then
                          xAdj(-3:2,-3:2,1,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,-1,1:3,sps)=x2(-3:2,-3:2,1:3)
                       else
                          xAdj(-3:2,-3:2,2,1:3,sps)=x0(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,1,1:3,sps)=x1(-3:2,-3:2,1:3)
                          xAdj(-3:2,-3:2,0,1:3,sps)=x2(-3:2,-3:2,1:3)
                       endif
                    endif
               end select 
            endif testSymmetry
         enddo loopBocos
 !     enddo domains
 !  enddo spectralLoop

!!$         !print out xAdj
!!$         istart2 = -3
!!$         jstart2 = -3
!!$         kstart2 = -3
!!$         iend2 = 2
!!$         jend2 = 2
!!$         kend2 = 2
!!$         if(icell==2) istart2=-2
!!$         if(jcell==2) jstart2=-2
!!$         if(kcell==2) kstart2=-2
!!$         if(icell==il) iend2=1
!!$         if(jcell==jl) jend2=1
!!$         if(kcell==kl) kend2=1
!!$         do iiii = istart2,iend2
!!$            do jjjj = jstart2,jend2
!!$               do kkkk = kstart2,kend2
!!$                  do n = 1,3
!!$                     i = icell+iiii
!!$                     j = jcell+jjjj
!!$                     k = kcell+kkkk
!!$                     write(unitxAD,10) i,j,k,n,xAdj(iiii,jjjj,kkkk,n)
!!$10                   format(1x,'res',4I8,f20.14)
!!$                  enddo
!!$               enddo
!!$            enddo
!!$         enddo

!!$!
!!$!      ******************************************************************
!!$!      *                                                                *
!!$!      * Exchange the coordinates for the internal halo's.              *
!!$!      *                                                                *
!!$!      ******************************************************************
!!$!
!!$   call exchangeCoor(level)
   
       end subroutine xhaloAdjTS

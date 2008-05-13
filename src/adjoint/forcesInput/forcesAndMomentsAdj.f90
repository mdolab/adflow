!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesAndMomentsAdj.f90                            *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi                                    *
!      * Starting date: 12-18-2007                                      *
!      * Last modified: 12-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine forcesAndMomentsAdj(cFpAdj,cMpAdj, &
           cFpAdjOut,cMpAdjOut, &
           yplusMax,refPoint,siAdj,sjAdj,skAdj,normAdj,xAdj,pAdj,wAdj,&
           iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End, &
           level,mm,nn)
!
!      ******************************************************************
!      *                                                                *
!      * forcesAndMoments computes the contribution of the block        *
!      * given by the pointers in blockPointers to the force and        *
!      * moment coefficients of the geometry. A distinction is made     *
!      * between the inviscid and viscous parts. In case the maximum    *
!      * yplus value must be monitored (only possible for RANS), this   *
!      * value is also computed.                                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers !ie,je,ke
       use BCTypes
       use flowVarRefState
       use inputPhysics
       use communication ! myID for debug
       implicit none

!
!      Subroutine arguments
!
!s       integer(kind=intType),intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
!s
!s       real(kind=realType), dimension(3), intent(inout) :: cFpAdj, cFvAdj
!s       real(kind=realType), dimension(3), intent(inout) :: cMpAdj, cMvAdj
!s
!s       !add to allow for scaling!
!s       real(kind=realType), dimension(3), intent(out) :: cFpAdjOut, cFvAdjOut
!s       real(kind=realType), dimension(3), intent(out) :: cMpAdjOut, cMvAdjOut
!s       real(kind=realType), dimension(3), intent(in) :: refPoint
!s       real(kind=realType), intent(out), intent(in) :: yplusMax
!s       integer(kind=intType), intent(in) :: mm, level, nn

!s       real(kind=realType), dimension(2,1:je,1:ke,3), intent(in) :: siAdj
!s       real(kind=realType), dimension(1:ie,2,1:ke,3), intent(in) :: sjAdj
!s       real(kind=realType), dimension(1:ie,1:je,2,3), intent(in) :: skAdj
!s       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3), intent(in) :: normAdj
!s       ! Note that the range of xAdj should correspond to 0:ie,0:je,0:ke
!s       real(kind=realType), dimension(0:ie,0:je,0:ke,3), intent(in) :: xAdj


       integer(kind=intType) :: iiBeg,iiEnd,jjBeg,jjEnd
       integer(kind=intType) :: i2Beg,i2End,j2Beg,j2End
       integer(kind=intType) :: mm, level, nn

       real(kind=realType), dimension(3),intent(inout) :: cFpAdj !, cFvAdj
       real(kind=realType), dimension(3),intent(inout) :: cMpAdj !, cMvAdj
       !add to allow for scaling!
       real(kind=realType), dimension(3), intent(out) :: cFpAdjOut !, cFvAdjOut
       real(kind=realType), dimension(3),intent(out) :: cMpAdjOut !, cMvAdjOut
       real(kind=realType), dimension(3),intent(in) :: refPoint
       real(kind=realType),intent(in) :: yplusMax


       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3),intent(in) :: siAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3),intent(in) :: sjAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3),intent(in) :: skAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3),intent(in) :: normAdj
       ! Note that the range of xAdj should correspond to 0:ie,0:je,0:ke
       real(kind=realType), dimension(0:ie,0:je,0:ke,3),intent(in) :: xAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb,nw),intent(in) :: wAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb),intent(in) :: pAdj

!
!      Local variables.
!
       integer(kind=intType) ::  i, j, k, ii

       real(kind=realType) :: pm1, fx, fy, fz, fn
       real(kind=realType) :: xc, yc, zc
       real(kind=realType) :: fact,  dwall

!v       real(kind=realType) :: tauxx, tauyy, tauzz
!v       real(kind=realType) :: tauxy, tauxz, tauyz

       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd)   :: pp2, pp1
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd)   :: rho2, rho1
!       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: xx
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: xx
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: ss

!v       real(kind=realType), dimension(:,:),   pointer :: rlv2, rlv1
!v       real(kind=realType), dimension(:,:),   pointer :: dd2Wall
!v       real(kind=realType), dimension(:,:,:), pointer :: ss, xx


       logical :: viscousSubface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary subfaces of this block.
       
       invForce: if(BCType(mm) == EulerWall        .or. &
            BCType(mm) == NSWallAdiabatic .or.  &
            BCType(mm) == NSWallIsothermal) then
          
          ! Subface is a wall. Check if it is a viscous wall.
          
          
          viscousSubface = .true.
          if(BCType(mm) == EulerWall) viscousSubface = .false.

          ! Set a bunch of pointers depending on the face id to make
          ! a generic treatment possible. The routine setBcPointers
          ! is not used, because quite a few other ones are needed.
          
          select case (BCFaceID(mm))

             case (iMin)
!                print *,'imin'
               pp2(:,:)  = pAdj(2,iiBeg:iiEnd,jjBeg:jjEnd)
               pp1(:,:)  = pAdj(1,iiBeg:iiEnd,jjBeg:jjEnd)
               rho2(:,:) = wAdj(2,iiBeg:iiEnd,jjBeg:jjEnd,irho)
               rho1(:,:) = wAdj(1,iiBeg:iiEnd,jjBeg:jjEnd,irho)
               ss(:,:,:) = siAdj(1,iiBeg:iiEnd,jjBeg:jjEnd,:)
!               xx(:,:,:)   = xAdj(1,iiBeg-1:iiEnd,jjBeg-1:jjEnd,:)
!               print *,'shape',shape(xx),shape(xAdj)
               xx(:,:,:)   = xAdj(1,iiBeg-1:iiEnd-1,jjBeg-1:jjEnd-1,:)
!               xx(:,:,:)   = xAdj(1,:,:,:)
               fact = -one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(2,1:,1:); rlv1 => rlv(1,1:,1:)
!v               endif

             !===========================================================

             case (iMax)
!                print *,'imax'
               pp2(:,:)  = pAdj(il,iiBeg:iiEnd,jjBeg:jjEnd)     
               pp1(:,:)  = pAdj(ie,iiBeg:iiEnd,jjBeg:jjEnd)
               rho2(:,:) = wAdj(il,iiBeg:iiEnd,jjBeg:jjEnd,irho)
               rho1(:,:) = wAdj(ie,iiBeg:iiEnd,jjBeg:jjEnd,irho)
               ss(:,:,:)   = siAdj(2,iiBeg:iiEnd,jjBeg:jjEnd,:) 
!               xx(:,:,:)   = xAdj(il,iiBeg-1:iiEnd,jjBeg-1:jjEnd,:)
               xx(:,:,:)   = xAdj(il,iiBeg-1:iiEnd-1,jjBeg-1:jjEnd-1,:)
!               xx(:,:,:)   = xAdj(il,:,:,:)
               fact = one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(il,1:,1:); rlv1 => rlv(ie,1:,1:)
!v               endif

             !===========================================================

             case (jMin)
!                print *,'jmin'
               pp2(:,:)  = pAdj(iiBeg:iiEnd,2,jjBeg:jjEnd);
               pp1(:,:)  = pAdj(iiBeg:iiEnd,1,jjBeg:jjEnd)
!!$               pp4  => p(1:,2,1:);      pp3  => p(1:,1,1:)
!!$               ss2   => sj(:,1,:,:); xx2   => x(:,1,:,:)
!!$         
               !rho2 => w(1:,2,1:,irho); rho1 => w(1:,1,1:,irho)
               !ss   => sj(:,1,:,:);     xx   => x(:,1,:,:)
               !fact = -one

               rho2(:,:) = wAdj(iiBeg:iiEnd,2,jjBeg:jjEnd,irho);
               rho1(:,:) = wAdj(iiBeg:iiEnd,1,jjBeg:jjEnd,irho)
               ss(:,:,:)   = sjAdj(iiBeg:iiEnd,1,jjBeg:jjEnd,:);
!               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd,1,jjBeg-1:jjEnd,:)
!               print *,'indices',iibeg,i2beg,iiend,i2end
!               print *,'shape',shape(xx),shape(xAdj(iiBeg-1:iiEnd,1,jjBeg-1:jjEnd,:))
               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd-1,1,jjBeg-1:jjEnd-1,:)
               !xx(:,:,:)   = xAdj(:,1,:,:)
!               print *,'xxjmin',xx(1,1,1),xadj(1,1,1,1),xx(1,2,1),xadj(1,1,2,1)
!               print *,'shape2',shape(xx),shape(xAdj(iiBeg-1:iiEnd,1,jjBeg-1:jjEnd,:))
!               stop
!               print *,'xx',xx(1,1,1),xx2(1,1,1),xAdj(0,0,0,1),x(0,0,0,1),x(1,1,1,1)
!               stop
               fact = -one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(1:,2,1:); rlv1 => rlv(1:,1,1:)
!v               endif

             !===========================================================

             case (jMax)
!                print *,'jmax'
               pp2(:,:)  = pAdj(iiBeg:iiEnd,jl,jjBeg:jjEnd) 
               pp1(:,:)  = pAdj(iiBeg:iiEnd,je,jjBeg:jjEnd)

!!$               pp4  => p(1:,jl,1:);      pp3  => p(1:,je,1:)
!!$               ss2   => sj(:,jl,:,:);     xx2   => x(:,jl,:,:)
               
               rho2(:,:) = wAdj(iiBeg:iiEnd,jl,jjBeg:jjEnd,irho)
               rho1(:,:) = wAdj(iiBeg:iiEnd,je,jjBeg:jjEnd,irho)
               ss(:,:,:)   = sjAdj(iiBeg:iiEnd,2,jjBeg:jjEnd,:) 
!               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd,jl,jjBeg-1:jjEnd,:)
               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd-1,jl,jjBeg-1:jjEnd-1,:)
!               xx(:,:,:)   = xAdj(:,jl,:,:)

               fact = one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(1:,jl,1:); rlv1 => rlv(1:,je,1:)
!v               endif

             !===========================================================

             case (kMin)
!                print *,'kmin'
               pp2(:,:)  = pAdj(iiBeg:iiEnd,jjBeg:jjEnd,2) 
               pp1(:,:)  = pAdj(iiBeg:iiEnd,jjBeg:jjEnd,1)

!!$               pp4  => p(1:,1:,2);      pp3  => p(1:,1:,1)
!!$               ss2   => sk(:,:,1,:);     xx2   => x(:,:,1,:)

               rho2(:,:) = wAdj(iiBeg:iiEnd,jjBeg:jjEnd,2,irho)
               rho1(:,:) = wAdj(iiBeg:iiEnd,jjBeg:jjEnd,1,irho)
               ss(:,:,:)   = skAdj(iiBeg:iiEnd,jjBeg:jjEnd,1,:); 
!               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,1,:)
               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd-1,jjBeg-1:jjEnd-1,1,:)
!               xx(:,:,:)   = xAdj(:,:,1,:)
               fact = -one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(1:,1:,2); rlv1 => rlv(1:,1:,1)
!v               endif

             !===========================================================

             case (kMax)
!                print *,'kmax'
               pp2(:,:)  = pAdj(iiBeg:iiEnd,jjBeg:jjEnd,kl)
               pp1(:,:)  = pAdj(iiBeg:iiEnd,jjBeg:jjEnd,ke)

!!$               pp4  => p(1:,1:,kl);      pp3  => p(1:,1:,ke)
!!$               ss2   => sk(:,:,kl,:);     xx2   => x(:,:,kl,:)

               rho2(:,:) = wAdj(iiBeg:iiEnd,jjBeg:jjEnd,kl,irho) 
               rho1(:,:) = wAdj(iiBeg:iiEnd,jjBeg:jjEnd,ke,irho)
               ss(:,:,:)   = skAdj(iiBeg:iiEnd,jjBeg:jjEnd,2,:)
!               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd,jjBeg-1:jjEnd,kl,:)
               xx(:,:,:)   = xAdj(iiBeg-1:iiEnd-1,jjBeg-1:jjEnd-1,kl,:)
!               xx(:,:,:)   = xAdj(:,:,kl,:)
               fact = one

!v               if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)
!v               if( viscousSubface ) then
!v                 rlv2 => rlv(1:,1:,kl); rlv1 => rlv(1:,1:,ke)
!v               endif

           end select

           ! Loop over the quadrilateral faces of the subface. 
           ! Note that +1 to Beg, and -1 to End to have ranges for the owned cell.

           !print *,'indicies',i2beg,i2end,iiend,j2beg,j2end!,shape(xx)
           do j=j2Beg,j2End
             do i=i2Beg,i2End
                !print *,'indices',i,j,xx(i,j,1),xAdj(i-1,1,j-1,1),x(i-1,1,j-1,1)
               ! Compute the average pressure minus 1 and the coordinates
               ! of the centroid of the face relative from from the
               ! moment reference point. Due to the usage of pointers for
               ! the coordinates, whose original array starts at 0, an
               ! offset of 1 must be used. The pressure is multipled by
               ! fact to account for the possibility of an inward or
               ! outward pointing normal.
               
               pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)
               !pm2 = fact*(half*(pp4(i,j) + pp3(i,j)) - pInf)
               
               !print *,'pm comparison',pm1,pm2

               xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
                  +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
               yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
                  +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
               zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
                  +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

!!$               xc2 = fourth*(xx2(i,j,  1) + xx2(i+1,j,  1) &
!!$                  +         xx2(i,j+1,1) + xx2(i+1,j+1,1)) - refPoint(1)
!!$               yc2 = fourth*(xx2(i,j,  2) + xx2(i+1,j,  2) &
!!$                  +         xx2(i,j+1,2) + xx2(i+1,j+1,2)) - refPoint(2)
!!$               zc2 = fourth*(xx2(i,j,  3) + xx2(i+1,j,  3) &
!!$                  +         xx2(i,j+1,3) + xx2(i+1,j+1,3)) - refPoint(3)

               !print *,'shape',shape(xx),shape(xx2)
               !stop
!               print *,'xc Comp',xc,xc2,xx(i,j,1),xx2(i,j,1),refpoint(1)
!               print *,'yc Comp',yc,yc2,xx(i,j,2),xx2(i,j,2),refpoint(2)
!               print *,'zc Comp',zc,zc2,xx(i,j,3),xx2(i,j,3),refpoint(3)

!!$               pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)
!!$
!!$               xc = fourth*(xx(i,j,  1) + xx(i-1,j,  1) &
!!$                  +         xx(i,j-1,1) + xx(i-1,j-1,1)) - refPoint(1)
!!$               yc = fourth*(xx(i,j,  2) + xx(i-1,j,  2) &
!!$                  +         xx(i,j-1,2) + xx(i-1,j-1,2)) - refPoint(2)
!!$               zc = fourth*(xx(i,j,  3) + xx(i-1,j,  3) &
!!$                  +         xx(i,j-1,3) + xx(i-1,j-1,3)) - refPoint(3)

               ! Compute the force components.

               fx = pm1*ss(i,j,1)
               fy = pm1*ss(i,j,2)
               fz = pm1*ss(i,j,3)

!!$               fx2 = pm1*ss2(i,j,1)
!!$               fy2 = pm1*ss2(i,j,2)
!!$               fz2 = pm1*ss2(i,j,3)
!!$
!!$               !print *,'fx comparison',fx,fx2

               ! Update the inviscid force and moment coefficients.

               cFpAdj(1) = cFpAdj(1) + fx
               cFpAdj(2) = cFpAdj(2) + fy
               cFpAdj(3) = cFpAdj(3) + fz

               cMpAdj(1) = cMpAdj(1) + yc*fz - zc*fy
               cMpAdj(2) = cMpAdj(2) + zc*fx - xc*fz
               cMpAdj(3) = cMpAdj(3) + xc*fy - yc*fx
!!$               
!!$               print *,'forcemoment',cfpadj(1),cmpadj(1)

             enddo
           enddo

!
!          **************************************************************
!          *                                                            *
!          * Integration of the viscous forces.                         *
!          * Only for viscous boundaries.                               *
!          *                                                            *
!          **************************************************************
!
           visForce: if( viscousSubface ) then

              ! It hasn't been implemented yet. 

!!$               ! Initialize dwall for the laminar case and set the pointer
!!$             ! for the unit normals.
!!$
!!$             dwall = zero
!!$             norm => BCData(nn)%norm
!!$
!!$             ! Loop over the quadrilateral faces of the subface and
!!$             ! compute the viscous contribution to the force and
!!$             ! moment and update the maximum value of y+.
!!$
!!$             do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
!!$               do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd
!!$
!!$                 ! Store the viscous stress tensor a bit easier.
!!$
!!$                 tauXx = viscSubface(nn)%tau(i,j,1)
!!$                 tauYy = viscSubface(nn)%tau(i,j,2)
!!$                 tauZz = viscSubface(nn)%tau(i,j,3)
!!$                 tauXy = viscSubface(nn)%tau(i,j,4)
!!$                 tauXz = viscSubface(nn)%tau(i,j,5)
!!$                 tauYz = viscSubface(nn)%tau(i,j,6)
!!$
!!$                 ! Compute the viscous force on the face. A minus sign
!!$                 ! is now present, due to the definition of this force.
!!$
!!$                 fx = -fact*(tauXx*ss(i,j,1) + tauXy*ss(i,j,2) &
!!$                    +        tauXz*ss(i,j,3))
!!$                 fy = -fact*(tauXy*ss(i,j,1) + tauYy*ss(i,j,2) &
!!$                    +        tauYz*ss(i,j,3))
!!$                 fz = -fact*(tauXz*ss(i,j,1) + tauYz*ss(i,j,2) &
!!$                    +        tauZz*ss(i,j,3))
!!$
!!$                 ! Compute the coordinates of the centroid of the face
!!$                 ! relative from the moment reference point. Due to the
!!$                 ! usage of pointers for xx and offset of 1 is present,
!!$                 ! because x originally starts at 0.
!!$
!!$                 xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
!!$                    +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
!!$                 yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
!!$                    +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
!!$                 zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
!!$                    +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)
!!$
!!$                 ! Update the viscous force and moment coefficients.
!!$
!!$                 cFv(1) = cFv(1) + fx
!!$                 cFv(2) = cFv(2) + fy
!!$                 cFv(3) = cFv(3) + fz
!!$
!!$                 cMv(1) = cMv(1) + yc*fz - zc*fy
!!$                 cMv(2) = cMv(2) + zc*fx - xc*fz
!!$                 cMv(3) = cMv(3) + xc*fy - yc*fx
!!$
!!$                 ! Compute the tangential component of the stress tensor,
!!$                 ! which is needed to monitor y+. The result is stored
!!$                 ! in fx, fy, fz, although it is not really a force.
!!$                 ! As later on only the magnitude of the tangential
!!$                 ! component is important, there is no need to take the
!!$                 ! sign into account (it should be a minus sign).
!!$
!!$                 fx = tauXx*norm(i,j,1) + tauXy*norm(i,j,2) &
!!$                    + tauXz*norm(i,j,3)
!!$                 fy = tauXy*norm(i,j,1) + tauYy*norm(i,j,2) &
!!$                    + tauYz*norm(i,j,3)
!!$                 fz = tauXz*norm(i,j,1) + tauYz*norm(i,j,2) &
!!$                    + tauZz*norm(i,j,3)
!!$
!!$                 fn = fx*norm(i,j,1) + fy*norm(i,j,2) + fz*norm(i,j,3)
!!$
!!$                 fx = fx - fn*norm(i,j,1)
!!$                 fy = fy - fn*norm(i,j,2)
!!$                 fz = fz - fn*norm(i,j,3)
!!$
!!$                 ! Compute the local value of y+. Due to the usage
!!$                 ! of pointers there is on offset of -1 in dd2Wall..
!!$
!!$                 if(equations == RANSEquations) dwall = dd2Wall(i-1,j-1)
!!$
!!$                 rho   = half*(rho2(i,j) + rho1(i,j))
!!$                 mul   = half*(rlv2(i,j) + rlv1(i,j))
!!$                 yplus = sqrt(rho*sqrt(fx*fx + fy*fy + fz*fz))*dwall/mul
!!$
!!$                 ! Store this value if this value is larger than the
!!$                 ! currently stored value.
!!$
!!$                 yplusMax = max(yplusMax, yplus)
!!$
!!$               enddo
!!$             enddo

           endif visForce
        endif invForce
       ! Currently the coefficients only contain the surface integral
       ! of the pressure tensor. These values must be scaled to
       ! obtain the correct coefficients.

       fact = two/(gammaInf*pInf*MachCoef*MachCoef &
            *      surfaceRef*LRef*LRef)
       cFpAdjOut(1) = cFpAdj(1)*fact; 
       cFpAdjOut(2) = cFpAdj(2)*fact; 
       cFpAdjOut(3) = cFpAdj(3)*fact
!s       cFv(1) = cFv(1)*fact; cFv(2) = cFv(2)*fact; cFv(3) = cFv(3)*fact

       fact = fact/(lengthRef*LRef)
       cMpAdjOut(1) = cMpAdj(1)*fact; 
       cMpAdjOut(2) = cMpAdj(2)*fact; 
       cMpAdjOut(3) = cMpAdj(3)*fact
!s       cMv(1) = cMv(1)*fact; cMv(2) = cMv(2)*fact; cMv(3) = cMv(3)*fact

       end subroutine forcesAndMomentsAdj

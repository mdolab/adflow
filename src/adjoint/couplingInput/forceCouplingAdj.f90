!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesCouplingAdj.f90                           *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 08-17-2008                                      *
!      * Last modified: 08-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine forcesCouplingAdj(yplusMax,refPoint,siAdj,sjAdj,skAdj,&
            normAdj,xAdj,pAdj,wAdj,iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,&
            j2Beg,j2End,level,mm,nn,machCoefAdj,forceloc,nSurfNodesLoc,ii)
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
       integer(kind=intType) :: nSurfNodesLoc

       real(kind=realType), dimension(3),intent(in) :: refPoint
       real(kind=realType),intent(in) :: yplusMax


       real(kind=realType) :: scaleDim, pp

       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3),intent(in) :: siAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3),intent(in) :: sjAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3),intent(in) :: skAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3),intent(in) :: normAdj
       ! Note that the range of xAdj should correspond to 0:ie,0:je,0:ke
       real(kind=realType), dimension(0:ie,0:je,0:ke,3),intent(in) :: xAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb,nw),intent(in) :: wAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb),intent(in) :: pAdj
       real(kind=realType) :: machCoefAdj
!
!      Local variables.
!
       integer(kind=intType) ::  i, j, k, ii,jj

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

       real(kind=realType), dimension(3,nSurfNodesLoc) :: forceLoc

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

       !forceLoc = zero

       ! Compute the scaling factor to create the correct dimensional
       ! force in newton. As the coordinates are already in meters,
       ! this scaling factor is pRef.

       scaleDim = pRef

       ! Compute the local forces. Take the scaling factor into
       ! account to obtain the forces in SI-units, i.e. Newton.

       ! Loop over the boundary subfaces of this block.
       !print *,'bctype',BCType(mm),mm

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

                 !print *,'forces',fx,fy,fz,pp
                 
                 ! Distribute the force to the 4 nodes of the quad.
                 ! Note that the averaging factor 1/4 has already been
                 ! taken into account in pp.

                 jj = ii + (j-j2Beg)*(i2End-i2Beg+2) + i-i2Beg + 1
                 !print *,'jj1',jj
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + 1
                 !print *,'jj2',jj
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + i2End - i2Beg + 1
                 !print *,'jj3',jj
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 jj = jj + 1
                 !print *,'jj4',jj
                 forceLoc(1,jj) = forceLoc(1,jj) + fx
                 forceLoc(2,jj) = forceLoc(2,jj) + fy
                 forceLoc(3,jj) = forceLoc(3,jj) + fz

                 !print *,'forceLoc',forceLoc(1,jj),forceLoc(2,jj),forceLoc(3,jj)
             enddo
           enddo

           ! For a navier-stokes boundary also the viscous
           ! contribution is taken into account.
!!$
!!$           viscTest: if(BCType(mm) == NSWallAdiabatic .or. &
!!$                BCType(mm) == NSWallIsothermal) then
!!$
!!$              ! Compute the viscous force on each of the faces and
!!$              ! scatter it to the 4 nodes, whose forces are updated.
!!$              
!!$              do j=jBeg,jEnd
!!$                 do i=iBeg,iEnd
!!$                    
!!$                    ! Store the components of the stress tensor
!!$                    ! a bit easier.
!!$
!!$                   tauxx = viscSubface(mm)%tau(i,j,1)
!!$                   tauyy = viscSubface(mm)%tau(i,j,2)
!!$                   tauzz = viscSubface(mm)%tau(i,j,3)
!!$                   tauxy = viscSubface(mm)%tau(i,j,4)
!!$                   tauxz = viscSubface(mm)%tau(i,j,5)
!!$                   tauyz = viscSubface(mm)%tau(i,j,6)
!!$
!!$                   ! Compute the viscous force on the face. A minus sign
!!$                   ! is now present, due to the definition of this force.
!!$                   ! Furthermore, the conversion to s.I. Units is made
!!$                   ! and the averaging factor of 1/4 is taken into
!!$                   ! account.
!!$
!!$                   fx = -fourth*fact*scaleDim*(tauxx*ss(i,j,1) &
!!$                      +                        tauxy*ss(i,j,2) &
!!$                      +                        tauxz*ss(i,j,3))
!!$                   fy = -fourth*fact*scaleDim*(tauxy*ss(i,j,1) &
!!$                      +                        tauyy*ss(i,j,2) &
!!$                      +                        tauyz*ss(i,j,3))
!!$                   fy = -fourth*fact*scaleDim*(tauxz*ss(i,j,1) &
!!$                      +                        tauyz*ss(i,j,2) &
!!$                      +                        tauzz*ss(i,j,3))
!!$
!!$                   ! Distribute the force to the 4 nodes of the quad.
!!$
!!$                   jj = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
!!$                   forceLoc(1,jj) = forceLoc(1,jj) + fx
!!$                   forceLoc(2,jj) = forceLoc(2,jj) + fy
!!$                   forceLoc(3,jj) = forceLoc(3,jj) + fz
!!$
!!$                   jj = jj + 1
!!$                   forceLoc(1,jj) = forceLoc(1,jj) + fx
!!$                   forceLoc(2,jj) = forceLoc(2,jj) + fy
!!$                   forceLoc(3,jj) = forceLoc(3,jj) + fz
!!$
!!$                   jj = jj + iEnd - iBeg + 1
!!$                   forceLoc(1,jj) = forceLoc(1,jj) + fx
!!$                   forceLoc(2,jj) = forceLoc(2,jj) + fy
!!$                   forceLoc(3,jj) = forceLoc(3,jj) + fz
!!$
!!$                   jj = jj + 1
!!$                   forceLoc(1,jj) = forceLoc(1,jj) + fx
!!$                   forceLoc(2,jj) = forceLoc(2,jj) + fy
!!$                   forceLoc(3,jj) = forceLoc(3,jj) + fz
!!$
!!$                 enddo
!!$               enddo
!!$
!!$             endif viscTest

             ! Update the counter ii.

             ii = ii + (j2End-j2Beg+2)*(i2End-i2Beg+2)

          endif invForce


        end subroutine forcesCouplingAdj

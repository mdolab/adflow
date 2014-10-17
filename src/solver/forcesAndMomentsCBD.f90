!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesAndMomentsCBD.f90                         *
!      * Author:        Edwin van der Weide,Eran Arad (eran-CBD)        *
!      * Starting date: 04-01-2003                                      *
!      * Last modified: 09-09-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine forcesAndMomentsCBD(nd,cFp, cFv, cMp, cMv)
!
!      ******************************************************************
!      *                                                                *
!      * forcesAndMoments computes the contribution of the block        *
!      * given by the pointers in blockPointers to the force and        *
!      * moment coefficients of the geometry. A distinction is made     *
!      * between the inviscid and viscous parts. In case the maximum    *
!      * yplus value must be monitored (only possible for rans), this   *
!      * value is also computed.                                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use flowVarRefState
       use inputPhysics
       use cgnsGrid  ! eran-CBD 
       implicit none
!
!      Subroutine arguments
!
       real(kind=realType), dimension(3,0:*), intent(out) :: cFp, cFv
       real(kind=realType), dimension(3,0:*), intent(out) :: cMp, cMv

!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j
       integer(kind=intType) :: nd, idWBCB ! eran-CBD

       real(kind=realType) :: pm1, fx, fy, fz
       real(kind=realType) :: xc, yc, zc
       real(kind=realType) :: fact

       real(kind=realType) :: tauXx, tauYy, tauZz
       real(kind=realType) :: tauXy, tauXz, tauYz

       real(kind=realType), dimension(3) :: refPoint

       real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
       real(kind=realType), dimension(:,:),   pointer :: rho2, rho1
       real(kind=realType), dimension(:,:),   pointer :: rlv2, rlv1
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall
       real(kind=realType), dimension(:,:,:), pointer :: ss, xx
       real(kind=realType), dimension(:,:,:), pointer :: norm

       logical :: viscousSubface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the reference point for the moment computation in
       ! meters.

       refPoint(1) = LRef*pointRef(1)
       refPoint(2) = LRef*pointRef(2)
       refPoint(3) = LRef*pointRef(3)

       ! Initialize the force and moment coefficients to 0 as well as
       ! yplusMax.

       surfacesCounter1 : do i=0,cgnsNWallSurfaces 
          cFp(1,i) = zero; cFp(2,i) = zero; cFp(3,i) = zero
          cFv(1,i) = zero; cFv(2,i) = zero; cFv(3,i) = zero
          cMp(1,i) = zero; cMp(2,i) = zero; cMp(3,i) = zero
          cMv(1,i) = zero; cMv(2,i) = zero; cMv(3,i) = zero
       end do surfacesCounter1

       ! Loop over the boundary subfaces of this block.

       bocos: do nn=1,nBocos
!
!        ****************************************************************
!        *                                                              *
!        * Integrate the inviscid contribution over the solid walls,    *
!        * either inviscid or viscous. The integration is done with     *
!        * cp. For closed contours this is equal to the integration     *
!        * of p; for open contours this is not the case anymore.        *
!        * Question is whether a force for an open contour is           *
!        * meaningful anyway.                                           *
!        *                                                              *
!        ****************************************************************
!

         invForce: if(BCType(nn) == EulerWall        .or. &
                       BCType(nn) == NSWallAdiabatic .or. &
                       BCType(nn) == NSWallIsothermal) then
!
            idWBCB = idWBC(nn) !  Locate the relevant wall surface
            cFp(1,0) = zero ; cFp(2,0) = zero ;  cFp(3,0) = zero
            cFv(1,0) = zero ; cFv(2,0) = zero ;  cFv(3,0) = zero
            cMp(1,0) = zero ; cMp(2,0) = zero ;  cMp(3,0) = zero
            cMv(1,0) = zero ; cMv(2,0) = zero ;  cMv(3,0) = zero

           ! Subface is a wall. Check if it is a viscous wall.

           viscousSubface = .true.
           if(BCType(nn) == EulerWall) viscousSubface = .false.

           ! Set a bunch of pointers depending on the face id to make
           ! a generic treatment possible. The routine setBcPointers
           ! is not used, because quite a few other ones are needed.

           select case (BCFaceID(nn))

             case (iMin)
               pp2  => p(2,1:,1:);      pp1  => p(1,1:,1:)
               rho2 => w(2,1:,1:,irho); rho1 => w(1,1:,1:,irho)
               ss   => si(1,:,:,:);     xx   => x(1,:,:,:)
               fact = -one

               if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)
               if( viscousSubface ) then
                 rlv2 => rlv(2,1:,1:); rlv1 => rlv(1,1:,1:)
               endif

             !===========================================================

             case (iMax)
               pp2  => p(il,1:,1:);      pp1  => p(ie,1:,1:)
               rho2 => w(il,1:,1:,irho); rho1 => w(ie,1:,1:,irho)
               ss   => si(il,:,:,:);     xx   => x(il,:,:,:)
               fact = one

               if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)
               if( viscousSubface ) then
                 rlv2 => rlv(il,1:,1:); rlv1 => rlv(ie,1:,1:)
               endif

             !===========================================================

             case (jMin)
               pp2  => p(1:,2,1:);      pp1  => p(1:,1,1:)
               rho2 => w(1:,2,1:,irho); rho1 => w(1:,1,1:,irho)
               ss   => sj(:,1,:,:);     xx   => x(:,1,:,:)
               fact = -one

               if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)
               if( viscousSubface ) then
                 rlv2 => rlv(1:,2,1:); rlv1 => rlv(1:,1,1:)
               endif

             !===========================================================

             case (jMax)
               pp2  => p(1:,jl,1:);      pp1  => p(1:,je,1:)
               rho2 => w(1:,jl,1:,irho); rho1 => w(1:,je,1:,irho)
               ss   => sj(:,jl,:,:);     xx   => x(:,jl,:,:)
               fact = one

               if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)
               if( viscousSubface ) then
                 rlv2 => rlv(1:,jl,1:); rlv1 => rlv(1:,je,1:)
               endif

             !===========================================================

             case (kMin)
               pp2  => p(1:,1:,2);      pp1  => p(1:,1:,1)
               rho2 => w(1:,1:,2,irho); rho1 => w(1:,1:,1,irho)
               ss   => sk(:,:,1,:);     xx   => x(:,:,1,:)
               fact = -one

               if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)
               if( viscousSubface ) then
                 rlv2 => rlv(1:,1:,2); rlv1 => rlv(1:,1:,1)
               endif

             !===========================================================

             case (kMax)
               pp2  => p(1:,1:,kl);      pp1  => p(1:,1:,ke)
               rho2 => w(1:,1:,kl,irho); rho1 => w(1:,1:,ke,irho)
               ss   => sk(:,:,kl,:);     xx   => x(:,:,kl,:)
               fact = one

               if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)
               if( viscousSubface ) then
                 rlv2 => rlv(1:,1:,kl); rlv1 => rlv(1:,1:,ke)
               endif

           end select

           ! Loop over the quadrilateral faces of the subface. Note
           ! that the nodal range of BCData must be used and not the
           ! cell range, because the latter may include the halo's in i
           ! and j-direction. The offset +1 is there, because inBeg and
           ! jnBeg refer to nodal ranges and not to cell ranges.

           do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

               ! Compute the average pressure minus 1 and the coordinates
               ! of the centroid of the face relative from from the
               ! moment reference point. Due to the usage of pointers for
               ! the coordinates, whose original array starts at 0, an
               ! offset of 1 must be used. The pressure is multipled by
               ! fact to account for the possibility of an inward or
               ! outward pointing normal.

               pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)

               xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
                  +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
               yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
                  +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
               zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
                  +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

               ! Compute the force components.

               fx = pm1*ss(i,j,1)
               fy = pm1*ss(i,j,2)
               fz = pm1*ss(i,j,3)

               ! Update the inviscid force and moment coefficients.

               cFp(1,idWBCB) = cFp(1,idWBCB) + fx
               cFp(2,idWBCB) = cFp(2,idWBCB) + fy
               cFp(3,idWBCB) = cFp(3,idWBCB) + fz

               cMp(1,idWBCB) = cMp(1,idWBCB) + yc*fz - zc*fy
               cMp(2,idWBCB) = cMp(2,idWBCB) + zc*fx - xc*fz
               cMp(3,idWBCB) = cMp(3,idWBCB) + xc*fy - yc*fx

             enddo
           enddo

           direction1 : do j=1,3
              surfacesCounter : do i=1,cgnsNWallSurfaces
                 cFp(j,0) = cFp(j,0) + cFp(j,i) 
                 cMp(j,0) = cMp(j,0) + cMp(j,i) 
              end do surfacesCounter
           end do direction1

!
!          **************************************************************
!          *                                                            *
!          * Integration of the viscous forces.                         *
!          * Only for viscous boundaries.                               *
!          *                                                            *
!          **************************************************************
!
           visForce: if( viscousSubface ) then

             norm => BCData(nn)%norm

             ! Loop over the quadrilateral faces of the subface and
             ! compute the viscous contribution to the force and
             ! moment and update the maximum value of y+.

             do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
               do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

                 ! Store the viscous stress tensor a bit easier.

                 tauXx = viscSubface(nn)%tau(i,j,1)
                 tauYy = viscSubface(nn)%tau(i,j,2)
                 tauZz = viscSubface(nn)%tau(i,j,3)
                 tauXy = viscSubface(nn)%tau(i,j,4)
                 tauXz = viscSubface(nn)%tau(i,j,5)
                 tauYz = viscSubface(nn)%tau(i,j,6)

                 ! Compute the viscous force on the face. A minus sign
                 ! is now present, due to the definition of this force.

                 fx = -fact*(tauXx*ss(i,j,1) + tauXy*ss(i,j,2) &
                    +        tauXz*ss(i,j,3))
                 fy = -fact*(tauXy*ss(i,j,1) + tauYy*ss(i,j,2) &
                    +        tauYz*ss(i,j,3))
                 fz = -fact*(tauXz*ss(i,j,1) + tauYz*ss(i,j,2) &
                    +        tauZz*ss(i,j,3))

                 ! Compute the coordinates of the centroid of the face
                 ! relative from the moment reference point. Due to the
                 ! usage of pointers for xx and offset of 1 is present,
                 ! because x originally starts at 0.

                 xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
                    +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
                 yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
                    +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
                 zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
                    +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

                 ! Update the viscous force and moment coefficients.

                 cFv(1,idWBCB) = cFv(1,idWBCB) + fx
                 cFv(2,idWBCB) = cFv(2,idWBCB) + fy
                 cFv(3,idWBCB) = cFv(3,idWBCB) + fz

                 cMv(1,idWBCB) = cMv(1,idWBCB) + yc*fz - zc*fy
                 cMv(2,idWBCB) = cMv(2,idWBCB) + zc*fx - xc*fz
                 cMv(3,idWBCB) = cMv(3,idWBCB) + xc*fy - yc*fx

               enddo
             enddo

           directions2 : do j=1,3
              surfacesCounter2 : do i=1,cgnsNWallSurfaces
                 cFv(j,0) = cFv(j,0) + cFv(j,i)
                 cMv(j,0) = cMv(j,0) + cMv(j,i)
              end do surfacesCounter2
           end do directions2


           endif visForce
         endif invForce

       enddo bocos

       ! Currently the coefficients only contain the surface integral
       ! of the pressure tensor. These values must be scaled to
       ! obtain the correct coefficients.

       fact = two/(gammaInf*pInf*MachCoef*MachCoef &
            *      surfaceRef*LRef*LRef)
       surfacesCounter3 : do i=0,cgnsNWallSurfaces
          cFp(1,i) = cFp(1,i)*fact; cFp(2,i) = cFp(2,i)*fact; cFp(3,i) = cFp(3,i)*fact
          cFv(1,i) = cFv(1,i)*fact; cFv(2,i) = cFv(2,i)*fact; cFv(3,i) = cFv(3,i)*fact
       end do surfacesCounter3

       fact = fact/(lengthRef*LRef)
       surfacesCounter4 : do i=0,cgnsNWallSurfaces
          cMp(1,i) = cMp(1,i)*fact; cMp(2,i) = cMp(2,i)*fact; cMp(3,i) = cMp(3,i)*fact
          cMv(1,i) = cMv(1,i)*fact; cMv(2,i) = cMv(2,i)*fact; cMv(3,i) = cMv(3,i)*fact
       end do surfacesCounter4

     end subroutine forcesAndMomentsCBD

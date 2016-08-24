!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesAndMoments.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-01-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine forcesAndMomentsMask(cFp, cFv, cMp, cMv, yplusMax, sepSensor, &
     mask, maskCount, nMask)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * forcesAndMoments computes the contribution of the block        *
  !      * given by the pointers in blockPointers to the force and        *
  !      * moment coefficients of the geometry. A distinction is made     *
  !      * between the inviscid and viscous parts. In case the maximum    *
  !      * yplus value must be monitored (only possible for rans), this   *
  !      * value is also computed. The separation sensor is also computed *
  !      * here.                                                          *
  !      ******************************************************************
  !
  use blockPointers
  use flowVarRefState
  use inputPhysics
  use costFunctions
  implicit none
  !
  !      Subroutine arguments
  !
  real(kind=realType), dimension(3), intent(out) :: cFp, cFv
  real(kind=realType), dimension(3), intent(out) :: cMp, cMv
  integer(kind=intType), dimension(nMask), intent(in) :: mask
  integer(kind=intType) :: maskCount, nMask

  real(kind=realType), intent(out) :: yplusMax, sepSensor
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, i, j, isize, jsize

  real(kind=realType) :: pm1, fx, fy, fz, fn
  real(kind=realType) :: xc, yc, zc
  real(kind=realType) :: fact, rho, mul, yplus, dwall
  real(kind=realType) :: scaleDim, V(3), sensor
  real(kind=realType) :: tauXx, tauYy, tauZz
  real(kind=realType) :: tauXy, tauXz, tauYz

  real(kind=realType), dimension(3) :: refPoint

  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType), dimension(:,:),   pointer :: rho2, rho1
  real(kind=realType), dimension(:,:),   pointer :: rlv2, rlv1
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall
  real(kind=realType), dimension(:,:,:), pointer :: ss, xx
  real(kind=realType), dimension(:,:,:), pointer :: ww2
  real(kind=realType), dimension(:,:,:), pointer :: norm
  real(kind=realType) :: mx, my, mz, qa
  logical :: viscousSubface
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Set the actual scaling factor such that ACTUAL forces are computed
  scaleDim = pRef/pInf

  ! Determine the reference point for the moment computation in
  ! meters.

  refPoint(1) = LRef*pointRef(1)
  refPoint(2) = LRef*pointRef(2)
  refPoint(3) = LRef*pointRef(3)

  ! Initialize the force and moment coefficients to 0 as well as
  ! yplusMax.

  cFp(1) = zero; cFp(2) = zero; cFp(3) = zero
  cFv(1) = zero; cFv(2) = zero; cFv(3) = zero
  cMp(1) = zero; cMp(2) = zero; cMp(3) = zero
  cMv(1) = zero; cMv(2) = zero; cMv(3) = zero

  yplusMax = zero
  sepSensor = zero
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
           ww2 => w(2,1:,1:,:)
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
           ww2    => w(il,1:,1:,:)
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
           ww2    => w(1:,2,1:,:)
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
           ww2    => w(1:,jl,1:,:)
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
           ww2    => w(1:,1:,2,:)
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
           ww2    => w(1:,1:,kl,:)
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

        jsize = BCData(nn)%jnEnd - bcData(nn)%jnBeg
        isize = BCData(nn)%inEnd - bcData(nn)%inBeg
        do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
           do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

              ! Compute the average pressure minus 1 and the coordinates
              ! of the centroid of the face relative from from the
              ! moment reference point. Due to the usage of pointers for
              ! the coordinates, whose original array starts at 0, an
              ! offset of 1 must be used. The pressure is multipled by
              ! fact to account for the possibility of an inward or
              ! outward pointing normal.

              pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)*scaleDim

              ! Account for the mask:
              pm1 = pm1 * mask(maskCount)
              maskCount = maskCount + 1

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
              cFp(1) = cFp(1) + fx
              cFp(2) = cFp(2) + fy
              cFp(3) = cFp(3) + fz

              mx = yc*fz - zc*fy
              my = zc*fx - xc*fz
              mz = xc*fy - yc*fx

              cMp(1) = cMp(1) + mx
              cMp(2) = cMp(2) + my
              cMp(3) = cMp(3) + mz

              ! Divide by 4 so we can scatter
              fx = fourth*fx
              fy = fourth*fy
              fz = fourth*fz

              ! Scatter 1/4 of the force to each of the nodes:
              BCData(nn)%F(i-1,j-1,:) = (/fx, fy, fz/)
              BCData(nn)%F(i  ,j-1,:) = (/fx, fy, fz/)
              BCData(nn)%F(i-1,j  ,:) = (/fx, fy, fz/)
              BCData(nn)%F(i  ,j  ,:) = (/fx, fy, fz/)

              ! Scatter a quarter of the area to each node:
              qA = fourth*sqrt(ss(i, j, 1)**2 + ss(i, j, 2)**2 + ss(i, j, 3)**2)
              BCData(nn)%dualArea(i-1, j-1) = BCData(nn)%dualArea(i-1, j-1) + qA
              BCData(nn)%dualArea(i  , j-1) = BCData(nn)%dualArea(i  , j-1) + qA
              BCData(nn)%dualArea(i-1, j  ) = BCData(nn)%dualArea(i-1, j  ) + qA
              BCData(nn)%dualArea(i  , j  ) = BCData(nn)%dualArea(i  , j  ) + qA

              ! Get normalized surface velocity:
              v(1) = ww2(i, j, ivx)
              v(2) = ww2(i, j, ivy)
              v(3) = ww2(i, j, ivz)
              v = v / (sqrt(v(1)**2 + v(2)**2 + v(3)**2) + 1e-16)
              
              ! Dot product with free stream
              sensor = -(v(1)*velDirFreeStream(1) + &
                   v(2)*velDirFreeStream(2) + &
                   v(3)*velDirFreeStream(3))
              
              !Now run through a smooth heaviside function:
              sensor = one/(one + exp(-2*sepSensorSharpness*(sensor-sepSensorOffset)))

              ! And integrate over the area of this cell and save:
              sensor = sensor * four * qA
              sepSensor = sepSensor + sensor

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

           ! Rewind the maskCount
           maskCount = maskCount - isize*jsize

           ! Initialize dwall for the laminar case and set the pointer
           ! for the unit normals.

           dwall = zero
           norm => BCData(nn)%norm

           ! Loop over the quadrilateral faces of the subface and
           ! compute the viscous contribution to the force and
           ! moment and update the maximum value of y+.
           !DEC$ NOVECTOR
           do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
              !DEC$ NOVECTOR
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
                      +        tauXz*ss(i,j,3))*scaleDim
                 fy = -fact*(tauXy*ss(i,j,1) + tauYy*ss(i,j,2) &
                      +        tauYz*ss(i,j,3))*scaleDim
                 fz = -fact*(tauXz*ss(i,j,1) + tauYz*ss(i,j,2) &
                      +        tauZz*ss(i,j,3))*scaleDim

                 ! Account for the mask
                 fx = fx * mask(maskCount)
                 fy = fy * mask(maskCount)
                 fz = fz * mask(maskCount)
                 maskCount = maskCount + 1

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

                 cFv(1) = cFv(1) + fx
                 cFv(2) = cFv(2) + fy
                 cFv(3) = cFv(3) + fz

                 mx = yc*fz - zc*fy
                 my = zc*fx - xc*fz
                 mz = xc*fy - yc*fx

                 cMv(1) = cMv(1) + mx
                 cMv(2) = cMv(2) + my
                 cMv(3) = cMv(3) + mz

                 ! Compute the tangential component of the stress tensor,
                 ! which is needed to monitor y+. The result is stored
                 ! in fx, fy, fz, although it is not really a force.
                 ! As later on only the magnitude of the tangential
                 ! component is important, there is no need to take the
                 ! sign into account (it should be a minus sign).

                 fx = tauXx*BCData(nn)%norm(i,j,1) + tauXy*BCData(nn)%norm(i,j,2) &
                      + tauXz*BCData(nn)%norm(i,j,3)
                 fy = tauXy*BCData(nn)%norm(i,j,1) + tauYy*BCData(nn)%norm(i,j,2) &
                      + tauYz*BCData(nn)%norm(i,j,3)
                 fz = tauXz*BCData(nn)%norm(i,j,1) + tauYz*BCData(nn)%norm(i,j,2) &
                      + tauZz*BCData(nn)%norm(i,j,3)

                 fn = fx*BCData(nn)%norm(i,j,1) + fy*BCData(nn)%norm(i,j,2) + fz*BCData(nn)%norm(i,j,3)

                 fx = fx - fn*BCData(nn)%norm(i,j,1)
                 fy = fy - fn*BCData(nn)%norm(i,j,2)
                 fz = fz - fn*BCData(nn)%norm(i,j,3)

                 ! Compute the local value of y+. Due to the usage
                 ! of pointers there is on offset of -1 in dd2Wall..

                 if(equations == RANSEquations) then 
                    dwall = dd2Wall(i-1,j-1)
                    rho   = half*(rho2(i,j) + rho1(i,j))
                    mul   = half*(rlv2(i,j) + rlv1(i,j))
                    yplus = sqrt(rho*sqrt(fx*fx + fy*fy + fz*fz))*dwall/mul

                    ! Store this value if this value is larger than the
                    ! currently stored value.

                    yplusMax = max(yplusMax, yplus)
                 end if

                 ! Divide by 4 so we can scatter
                 fx = fourth*fx
                 fy = fourth*fy
                 fz = fourth*fz

                 ! Scatter 1/4 of the force to each of the nodes:
                 BCData(nn)%F(i-1,j-1,:) = BCData(nn)%F(i-1,j-1,:) + (/fx, fy, fz/)
                 BCData(nn)%F(i  ,j-1,:) = BCData(nn)%F(i  ,j-1,:) + (/fx, fy, fz/)
                 BCData(nn)%F(i-1,j  ,:) = BCData(nn)%F(i-1,j  ,:) + (/fx, fy, fz/)
                 BCData(nn)%F(i  ,j  ,:) = BCData(nn)%F(i  ,j  ,:) + (/fx, fy, fz/)
              enddo
           end do
           ! If forces are tractions we have to divide by the dual area:
           if (forcesAsTractions) then
              do j= BCData(nn)%jnBeg, BCData(nn)%jnEnd
                 do i=BCData(nn)%inBeg, BCData(nn)%inEnd
                    bcData(nn)%F(i, j, :) =  bcData(nn)%F(i, j, :) / bcData(nn)%dualArea(i, j)
                 end do
              end do
           end if
        endif visForce
     end if invForce
  enddo bocos

  ! Currently the coefficients only contain the surface integral
  ! of the pressure tensor. These values must be scaled to
  ! obtain the correct coefficients.

  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
       *surfaceRef*LRef*LRef*scaleDim)

  cFp(1) = cFp(1)*fact; cFp(2) = cFp(2)*fact; cFp(3) = cFp(3)*fact
  cFv(1) = cFv(1)*fact; cFv(2) = cFv(2)*fact; cFv(3) = cFv(3)*fact

  fact = fact/(lengthRef*LRef)
  cMp(1) = cMp(1)*fact; cMp(2) = cMp(2)*fact; cMp(3) = cMp(3)*fact
  cMv(1) = cMv(1)*fact; cMv(2) = cMv(2)*fact; cMv(3) = cMv(3)*fact

end subroutine forcesAndMomentsMask

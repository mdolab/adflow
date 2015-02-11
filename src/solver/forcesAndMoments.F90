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
subroutine forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, Cavitation)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * forcesAndMoments computes the contribution of the block        *
  !      * given by the pointers in blockPointers to the force and        *
  !      * moment coefficients of the geometry. A distinction is made     *
  !      * between the inviscid and viscous parts. In case the maximum    *
  !      * yplus value must be monitored (only possible for rans), this   *
  !      * value is also computed. The separation sensor and the cavita-  *
  !      * tion sensor is also computed                                   *
  !      * here.                                                          *
  !      ******************************************************************
  !
  use blockPointers
  use BCTypes
  use flowVarRefState
  use inputPhysics
  use bcroutines
  implicit none
  !
  !      Subroutine arguments
  !
  real(kind=realType), dimension(3), intent(out) :: cFp, cFv
  real(kind=realType), dimension(3), intent(out) :: cMp, cMv

  real(kind=realType), intent(out) :: yplusMax, sepSensor, Cavitation
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, i, j, ii

  real(kind=realType) :: pm1, fx, fy, fz, fn, sigma
  real(kind=realType) :: xc, yc, zc
  real(kind=realType) :: fact, rho, mul, yplus, dwall
  real(kind=realType) :: scaleDim, V(3), sensor, sensor1, Cp, tmp, plocal
  real(kind=realType) :: tauXx, tauYy, tauZz
  real(kind=realType) :: tauXy, tauXz, tauYz

  real(kind=realType), dimension(3) :: refPoint
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
  Cavitation = zero
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

     invForce: if(BCType(nn) == EulerWall .or. &
          BCType(nn) == NSWallAdiabatic .or. &
          BCType(nn) == NSWallIsothermal) then

        ! Subface is a wall. Check if it is a viscous wall.

        viscousSubface = .true.
        if(BCType(nn) == EulerWall) viscousSubface = .false.

        ! Set a bunch of pointers depending on the face id to make
        ! a generic treatment possible. The routine setBcPointers
        ! is not used, because quite a few other ones are needed.
        call setBCPointers(nn, .True.)

        select case (BCFaceID(nn))
        case (iMin)
           fact = -one
        case (iMax)
           fact = one
        case (jMin)
           fact = -one
        case (jMax)
           fact = one
        case (kMin)
           fact = -one
        case (kMax)
           fact = one
        end select

        ! Loop over the quadrilateral faces of the subface. Note
        ! that the nodal range of BCData must be used and not the
        ! cell range, because the latter may include the halo's in i
        ! and j-direction. The offset +1 is there, because inBeg and
        ! jnBeg refer to nodal ranges and not to cell ranges.

        bcData(nn)%oArea(:,:) = zero
        !$AD II-LOOP
        do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1 
           i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
           j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

              ! Compute the average pressure minus 1 and the coordinates
              ! of the centroid of the face relative from from the
              ! moment reference point. Due to the usage of pointers for
              ! the coordinates, whose original array starts at 0, an
              ! offset of 1 must be used. The pressure is multipled by
              ! fact to account for the possibility of an inward or
              ! outward pointing normal.

              pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)*scaleDim

              xc = fourth*(xx(i,j,  1) + xx(i+1,j,  1) &
                   +         xx(i,j+1,1) + xx(i+1,j+1,1)) - refPoint(1)
              yc = fourth*(xx(i,j,  2) + xx(i+1,j,  2) &
                   +         xx(i,j+1,2) + xx(i+1,j+1,2)) - refPoint(2)
              zc = fourth*(xx(i,j,  3) + xx(i+1,j,  3) &
                   +         xx(i,j+1,3) + xx(i+1,j+1,3)) - refPoint(3)

              ! Compute the force components.
              fx = pm1*ssi(i,j,1)
              fy = pm1*ssi(i,j,2)
              fz = pm1*ssi(i,j,3)

              ! Store Force data on face
              BCData(nn)%Fp(i,j,1) = fx
              BCData(nn)%Fp(i,j,2) = fy
              BCData(nn)%Fp(i,j,3) = fz

              ! Scatter a quarter of the area to each node:
              qA = fourth*sqrt(ssi(i, j, 1)**2 + ssi(i, j, 2)**2 + ssi(i, j, 3)**2)
              BCData(nn)%oArea(i-1, j-1) = BCData(nn)%oArea(i-1, j-1) + qA
              BCData(nn)%oArea(i  , j-1) = BCData(nn)%oArea(i  , j-1) + qA
              BCData(nn)%oArea(i-1, j  ) = BCData(nn)%oArea(i-1, j  ) + qA
              BCData(nn)%oArea(i  , j  ) = BCData(nn)%oArea(i  , j  ) + qA
              
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
              sensor = one/(one + exp(-2*10*sensor))
                 
              ! And integrate over the area of this cell and save:
              sensor = sensor * four * qA
              sepSensor = sepSensor + sensor
              bcData(nn)%sepSensor(i, j) = sensor
    
              plocal = pp2(i,j)
		
   	      tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
	      
	      Cp = tmp*(plocal-pinf)
	      Sigma = 1.4
	      Sensor1 = -Cp - Sigma
	      !IF (sense >= 0) THEN
	      !Sensor = 1
	      !ELSE 
	      !Sensor = 0
	      !END IF

	      Sensor1 = one/(one+exp(-2*10*Sensor1))

	      Sensor1 = Sensor1 * four * qA
	      Cavitation = Cavitation + Sensor1
	      bcData(nn)%Cavitation(i,j) = Sensor1

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

              ! Store Moment data on face
              BCData(nn)%M(i,j,1) = mx
              BCData(nn)%M(i,j,2) = my
              BCData(nn)%M(i,j,3) = mz
              
          
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

           ! Initialize dwall for the laminar case and set the pointer
           ! for the unit normals.

           dwall = zero
           ! Replace norm with BCData norm - Peter Lyu
           !norm => BCData(nn)%norm

           ! Loop over the quadrilateral faces of the subface and
           ! compute the viscous contribution to the force and
           ! moment and update the maximum value of y+.

           do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1 
              i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
              j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1
              
                 ! Store the viscous stress tensor a bit easier.

                 tauXx = viscSubface(nn)%tau(i,j,1)
                 tauYy = viscSubface(nn)%tau(i,j,2)
                 tauZz = viscSubface(nn)%tau(i,j,3)
                 tauXy = viscSubface(nn)%tau(i,j,4)
                 tauXz = viscSubface(nn)%tau(i,j,5)
                 tauYz = viscSubface(nn)%tau(i,j,6)

                 ! Compute the viscous force on the face. A minus sign
                 ! is now present, due to the definition of this force.

                 fx = -fact*(tauXx*ssi(i,j,1) + tauXy*ssi(i,j,2) &
                      +        tauXz*ssi(i,j,3))*scaleDim
                 fy = -fact*(tauXy*ssi(i,j,1) + tauYy*ssi(i,j,2) &
                      +        tauYz*ssi(i,j,3))*scaleDim
                 fz = -fact*(tauXz*ssi(i,j,1) + tauYz*ssi(i,j,2) &
                      +        tauZz*ssi(i,j,3))*scaleDim

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

                 ! Store Force data on face
                 BCData(nn)%Fv(i,j,1) = fx
                 BCData(nn)%Fv(i,j,2) = fy
                 BCData(nn)%Fv(i,j,3) = fz

                 mx = yc*fz - zc*fy
                 my = zc*fx - xc*fz
                 mz = xc*fy - yc*fx

                 cMv(1) = cMv(1) + mx
                 cMv(2) = cMv(2) + my
                 cMv(3) = cMv(3) + mz

                 ! Store Moment data on face
                 BCData(nn)%M(i,j,1) = BCData(nn)%M(i,j,1) + mx
                 BCData(nn)%M(i,j,2) = BCData(nn)%M(i,j,2) + my
                 BCData(nn)%M(i,j,3) = BCData(nn)%M(i,j,3) + mz

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
#ifndef TAPENADE_REVERSE                 
                 if(equations == RANSEquations) then 
                    dwall = dd2Wall(i-1,j-1)
                    rho   = half*(ww2(i,j,irho) + ww1(i,j,irho))
                    mul   = half*(rlv2(i,j) + rlv1(i,j))
                    yplus = sqrt(rho*sqrt(fx*fx + fy*fy + fz*fz))*dwall/mul
                    
                    ! Store this value if this value is larger than the
                    ! currently stored value.
                    
                    yplusMax = max(yplusMax, yplus)
                 end if
#endif
              enddo
      
        else
           ! Zero the viscous force contribution
           BCData(nn)%Fv = zero
        endif visForce

        ! We have to inverse the nodal areas
        do j=(BCData(nn)%jnBeg),BCData(nn)%jnEnd
           do i=(BCData(nn)%inBeg),BCData(nn)%inEnd
              bcData(nn)%oArea(i,j) = one/bcData(nn)%oArea(i,j)
           end do
        end do
        call resetBCPointers(nn, .True.)

     endif invForce
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

end subroutine forcesAndMoments

   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of forcesandmoments in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: gammainf pinf pref *w *x *(*bcdata.fp)
   !                *(*bcdata.fv) *(*bcdata.m) *(*bcdata.oarea) *(*bcdata.sepsensor)
   !                *(*bcdata.cavitation) lengthref machcoef pointref
   !                cfp cfv cmp cmv cavitation sepsensor
   !   with respect to varying inputs: gammainf pinf pref *p *w *x
   !                *si *sj *sk *(*viscsubface.tau) *(*bcdata.fp)
   !                *(*bcdata.fv) *(*bcdata.m) *(*bcdata.oarea) *(*bcdata.sepsensor)
   !                *(*bcdata.cavitation) veldirfreestream lengthref
   !                machcoef pointref
   !   Plus diff mem management of: p:in w:in x:in si:in sj:in sk:in
   !                viscsubface:in *viscsubface.tau:in bcdata:in *bcdata.fp:in
   !                *bcdata.fv:in *bcdata.m:in *bcdata.oarea:in *bcdata.sepsensor:in
   !                *bcdata.cavitation:in
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
   SUBROUTINE FORCESANDMOMENTS_B(cfp, cfpd, cfv, cfvd, cmp, cmpd, cmv, cmvd&
   & , yplusmax, sepsensor, sepsensord, cavitation, cavitationd)
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
   USE BLOCKPOINTERS
   USE BCTYPES
   USE FLOWVARREFSTATE
   USE INPUTPHYSICS
   IMPLICIT NONE
   !
   !      Subroutine arguments
   !
   REAL(kind=realtype), DIMENSION(3) :: cfp, cfv
   REAL(kind=realtype), DIMENSION(3) :: cfpd, cfvd
   REAL(kind=realtype), DIMENSION(3) :: cmp, cmv
   REAL(kind=realtype), DIMENSION(3) :: cmpd, cmvd
   REAL(kind=realtype) :: yplusmax, sepsensor, cavitation
   REAL(kind=realtype) :: sepsensord, cavitationd
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, i, j
   REAL(kind=realtype) :: pm1, fx, fy, fz, fn, sigma
   REAL(kind=realtype) :: pm1d, fxd, fyd, fzd
   REAL(kind=realtype) :: xc, yc, zc
   REAL(kind=realtype) :: xcd, ycd, zcd
   REAL(kind=realtype) :: fact, rho, mul, yplus, dwall
   REAL(kind=realtype) :: factd
   REAL(kind=realtype) :: scaledim, v(3), sensor, sensor1, cp, tmp, &
   & plocal
   REAL(kind=realtype) :: scaledimd, vd(3), sensord, sensor1d, cpd, tmpd&
   & , plocald
   REAL(kind=realtype) :: tauxx, tauyy, tauzz
   REAL(kind=realtype) :: tauxxd, tauyyd, tauzzd
   REAL(kind=realtype) :: tauxy, tauxz, tauyz
   REAL(kind=realtype) :: tauxyd, tauxzd, tauyzd
   REAL(kind=realtype), DIMENSION(3) :: refpoint
   REAL(kind=realtype), DIMENSION(3) :: refpointd
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, nw) :: ww1, ww2
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, nw) :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: pp1, pp2
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rho2, rho1
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rho2d, rho1d
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rev1, rev2
   REAL(kind=realtype), DIMENSION(imaxdim - 2, jmaxdim - 2) :: dd2wall
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ss
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ssd
   REAL(kind=realtype), DIMENSION(imaxdim + 1, jmaxdim + 1, 3) :: xx
   REAL(kind=realtype), DIMENSION(imaxdim+1, jmaxdim+1, 3) :: xxd
   REAL(kind=realtype) :: mx, my, mz, qa
   REAL(kind=realtype) :: mxd, myd, mzd, qad
   LOGICAL :: viscoussubface
   INTRINSIC SQRT
   INTRINSIC EXP
   INTRINSIC MAX
   REAL(kind=realtype), DIMENSION(3) :: tmp0
   INTEGER :: ad_from
   INTEGER :: ad_to
   INTEGER :: ad_from0
   INTEGER :: ad_to0
   INTEGER :: ad_from1
   INTEGER :: ad_to1
   INTEGER :: ad_from2
   INTEGER :: ad_to2
   INTEGER :: ad_from3
   INTEGER :: ad_to3
   INTEGER :: ad_from4
   INTEGER :: ad_to4
   INTEGER :: branch
   REAL(kind=realtype) :: temp3
   REAL(kind=realtype) :: tempd14
   REAL(kind=realtype) :: temp2
   REAL(kind=realtype) :: tempd13
   REAL(kind=realtype) :: temp1
   REAL(kind=realtype) :: tempd12
   REAL(kind=realtype) :: temp0
   REAL(kind=realtype) :: tempd11
   REAL(kind=realtype) :: tempd10
   REAL(kind=realtype) :: tempd9
   REAL(kind=realtype) :: tempd
   REAL(kind=realtype) :: tempd8
   REAL(kind=realtype) :: tempd7
   REAL(kind=realtype) :: tempd6
   REAL(kind=realtype) :: tempd5
   REAL(kind=realtype) :: tempd4
   REAL(kind=realtype) :: tempd3
   REAL(kind=realtype) :: tempd2(3)
   REAL(kind=realtype) :: tempd1
   REAL(kind=realtype) :: tempd0
   REAL(kind=realtype) :: tmpd0(3)
   INTEGER :: ii1
   REAL(kind=realtype) :: temp
   REAL(kind=realtype) :: temp9
   REAL(kind=realtype) :: temp8
   REAL(kind=realtype) :: temp7
   REAL(kind=realtype) :: temp6
   REAL(kind=realtype) :: temp5
   REAL(kind=realtype) :: tempd16
   REAL(kind=realtype) :: temp4
   REAL(kind=realtype) :: tempd15
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Set the actual scaling factor such that ACTUAL forces are computed
   scaledim = pref/pinf
   ! Determine the reference point for the moment computation in
   ! meters.
   refpoint(1) = lref*pointref(1)
   refpoint(2) = lref*pointref(2)
   refpoint(3) = lref*pointref(3)
   ! Initialize the force and moment coefficients to 0 as well as
   ! yplusMax.
   cfp(1) = zero
   cfp(2) = zero
   cfp(3) = zero
   cfv(1) = zero
   cfv(2) = zero
   cfv(3) = zero
   cmp(1) = zero
   cmp(2) = zero
   cmp(3) = zero
   cmv(1) = zero
   cmv(2) = zero
   cmv(3) = zero
   ! Loop over the boundary subfaces of this block.
   bocos:DO nn=1,nbocos
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
   IF ((bctype(nn) .EQ. eulerwall .OR. bctype(nn) .EQ. nswalladiabatic)&
   &       .OR. bctype(nn) .EQ. nswallisothermal) THEN
   ! Subface is a wall. Check if it is a viscous wall.
   viscoussubface = .true.
   IF (bctype(nn) .EQ. eulerwall) viscoussubface = .false.
   ! Set a bunch of pointers depending on the face id to make
   ! a generic treatment possible. The routine setBcPointers
   ! is not used, because quite a few other ones are needed.
   CALL PUSHREAL8ARRAY(pp2, imaxdim*jmaxdim)
   CALL PUSHREAL8ARRAY(pp1, imaxdim*jmaxdim)
   CALL SETBCPOINTERSBWD(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, rev1&
   &                        , rev2, 0)
   CALL PUSHREAL8ARRAY(ss, imaxdim*jmaxdim*3)
   CALL SETXXSSRHODD2WALLBWD(nn, xx, ss, rho1, rho2, dd2wall)
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   CALL PUSHREAL8(fact)
   fact = -one
   CALL PUSHCONTROL3B(1)
   CASE (imax) 
   !===========================================================
   CALL PUSHREAL8(fact)
   fact = one
   CALL PUSHCONTROL3B(2)
   CASE (jmin) 
   !===========================================================
   CALL PUSHREAL8(fact)
   fact = -one
   CALL PUSHCONTROL3B(3)
   CASE (jmax) 
   !===========================================================
   CALL PUSHREAL8(fact)
   fact = one
   CALL PUSHCONTROL3B(4)
   CASE (kmin) 
   !===========================================================
   CALL PUSHREAL8(fact)
   fact = -one
   CALL PUSHCONTROL3B(5)
   CASE (kmax) 
   !===========================================================
   CALL PUSHREAL8(fact)
   fact = one
   CALL PUSHCONTROL3B(6)
   CASE DEFAULT
   CALL PUSHCONTROL3B(0)
   END SELECT
   ! Loop over the quadrilateral faces of the subface. Note
   ! that the nodal range of BCData must be used and not the
   ! cell range, because the latter may include the halo's in i
   ! and j-direction. The offset +1 is there, because inBeg and
   ! jnBeg refer to nodal ranges and not to cell ranges.
   CALL PUSHREAL8ARRAY(bcdata(nn)%oarea, SIZE(bcdata(nn)%oarea, 1)*&
   &                   SIZE(bcdata(nn)%oarea, 2))
   bcdata(nn)%oarea(:, :) = zero
   ad_from0 = bcdata(nn)%jnbeg + 1
   DO j=ad_from0,bcdata(nn)%jnend
   ad_from = bcdata(nn)%inbeg + 1
   DO i=ad_from,bcdata(nn)%inend
   ! Compute the average pressure minus 1 and the coordinates
   ! of the centroid of the face relative from from the
   ! moment reference point. Due to the usage of pointers for
   ! the coordinates, whose original array starts at 0, an
   ! offset of 1 must be used. The pressure is multipled by
   ! fact to account for the possibility of an inward or
   ! outward pointing normal.
   CALL PUSHREAL8(pm1)
   pm1 = fact*(half*(pp2(i, j)+pp1(i, j))-pinf)*scaledim
   CALL PUSHREAL8(xc)
   xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j&
   &           +1, 1)) - refpoint(1)
   CALL PUSHREAL8(yc)
   yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j&
   &           +1, 2)) - refpoint(2)
   CALL PUSHREAL8(zc)
   zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j&
   &           +1, 3)) - refpoint(3)
   ! Compute the force components.
   CALL PUSHREAL8(fx)
   fx = pm1*ss(i, j, 1)
   CALL PUSHREAL8(fy)
   fy = pm1*ss(i, j, 2)
   CALL PUSHREAL8(fz)
   fz = pm1*ss(i, j, 3)
   ! Store Force data on face
   ! Scatter a quarter of the area to each node:
   CALL PUSHREAL8(qa)
   qa = fourth*SQRT(ss(i, j, 1)**2+ss(i, j, 2)**2+ss(i, j, 3)**2)
   CALL PUSHREAL8(bcdata(nn)%oarea(i-1, j-1))
   bcdata(nn)%oarea(i-1, j-1) = bcdata(nn)%oarea(i-1, j-1) + qa
   CALL PUSHREAL8(bcdata(nn)%oarea(i, j-1))
   bcdata(nn)%oarea(i, j-1) = bcdata(nn)%oarea(i, j-1) + qa
   CALL PUSHREAL8(bcdata(nn)%oarea(i-1, j))
   bcdata(nn)%oarea(i-1, j) = bcdata(nn)%oarea(i-1, j) + qa
   CALL PUSHREAL8(bcdata(nn)%oarea(i, j))
   bcdata(nn)%oarea(i, j) = bcdata(nn)%oarea(i, j) + qa
   ! Get normalized surface velocity:
   CALL PUSHREAL8(v(1))
   v(1) = ww2(i, j, ivx)
   CALL PUSHREAL8(v(2))
   v(2) = ww2(i, j, ivy)
   CALL PUSHREAL8(v(3))
   v(3) = ww2(i, j, ivz)
   tmp0 = v/(SQRT(v(1)**2+v(2)**2+v(3)**2)+1e-16)
   CALL PUSHREAL8ARRAY(v, 3)
   v = tmp0
   ! Dot product with free stream
   CALL PUSHREAL8(sensor)
   sensor = -(v(1)*veldirfreestream(1)+v(2)*veldirfreestream(2)+v&
   &           (3)*veldirfreestream(3))
   !Now run through a smooth heaviside function:
   CALL PUSHREAL8(sensor)
   sensor = one/(one+EXP(-(2*10*sensor)))
   ! And integrate over the area of this cell and save:
   CALL PUSHREAL8(plocal)
   plocal = pp2(i, j)
   CALL PUSHREAL8(tmp)
   tmp = two/(gammainf*pinf*machcoef*machcoef)
   cp = tmp*(plocal-pinf)
   sigma = 1.4
   CALL PUSHREAL8(sensor1)
   sensor1 = -cp - sigma
   !IF (sense >= 0) THEN
   !Sensor = 1
   !ELSE 
   !Sensor = 0
   !END IF
   CALL PUSHREAL8(sensor1)
   sensor1 = one/(one+EXP(-(2*10*sensor1)))
   ! Update the inviscid force and moment coefficients.
   cfp(1) = cfp(1) + fx
   cfp(2) = cfp(2) + fy
   cfp(3) = cfp(3) + fz
   mx = yc*fz - zc*fy
   my = zc*fx - xc*fz
   mz = xc*fy - yc*fx
   cmp(1) = cmp(1) + mx
   cmp(2) = cmp(2) + my
   cmp(3) = cmp(3) + mz
   ! Store Moment data on face
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from0)
   !
   !          **************************************************************
   !          *                                                            *
   !          * Integration of the viscous forces.                         *
   !          * Only for viscous boundaries.                               *
   !          *                                                            *
   !          **************************************************************
   !
   IF (viscoussubface) THEN
   ad_from2 = bcdata(nn)%jnbeg + 1
   ! Replace norm with BCData norm - Peter Lyu
   !norm => BCData(nn)%norm
   ! Loop over the quadrilateral faces of the subface and
   ! compute the viscous contribution to the force and
   ! moment and update the maximum value of y+.
   !DEC$ NOVECTOR
   DO j=ad_from2,bcdata(nn)%jnend
   ad_from1 = bcdata(nn)%inbeg + 1
   !DEC$ NOVECTOR
   DO i=ad_from1,bcdata(nn)%inend
   ! Store the viscous stress tensor a bit easier.
   tauxx = viscsubface(nn)%tau(i, j, 1)
   tauyy = viscsubface(nn)%tau(i, j, 2)
   tauzz = viscsubface(nn)%tau(i, j, 3)
   tauxy = viscsubface(nn)%tau(i, j, 4)
   tauxz = viscsubface(nn)%tau(i, j, 5)
   tauyz = viscsubface(nn)%tau(i, j, 6)
   ! Compute the viscous force on the face. A minus sign
   ! is now present, due to the definition of this force.
   CALL PUSHREAL8(fx)
   fx = -(fact*(tauxx*ss(i, j, 1)+tauxy*ss(i, j, 2)+tauxz*ss(i&
   &             , j, 3))*scaledim)
   CALL PUSHREAL8(fy)
   fy = -(fact*(tauxy*ss(i, j, 1)+tauyy*ss(i, j, 2)+tauyz*ss(i&
   &             , j, 3))*scaledim)
   CALL PUSHREAL8(fz)
   fz = -(fact*(tauxz*ss(i, j, 1)+tauyz*ss(i, j, 2)+tauzz*ss(i&
   &             , j, 3))*scaledim)
   ! Compute the coordinates of the centroid of the face
   ! relative from the moment reference point. Due to the
   ! usage of pointers for xx and offset of 1 is present,
   ! because x originally starts at 0.
   CALL PUSHREAL8(xc)
   xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1&
   &             , j+1, 1)) - refpoint(1)
   CALL PUSHREAL8(yc)
   yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1&
   &             , j+1, 2)) - refpoint(2)
   CALL PUSHREAL8(zc)
   zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1&
   &             , j+1, 3)) - refpoint(3)
   ! Update the viscous force and moment coefficients.
   cfv(1) = cfv(1) + fx
   cfv(2) = cfv(2) + fy
   cfv(3) = cfv(3) + fz
   ! Store Force data on face
   mx = yc*fz - zc*fy
   my = zc*fx - xc*fz
   mz = xc*fy - yc*fx
   cmv(1) = cmv(1) + mx
   cmv(2) = cmv(2) + my
   cmv(3) = cmv(3) + mz
   ! Store Moment data on face
   ! Compute the tangential component of the stress tensor,
   ! which is needed to monitor y+. The result is stored
   ! in fx, fy, fz, although it is not really a force.
   ! As later on only the magnitude of the tangential
   ! component is important, there is no need to take the
   ! sign into account (it should be a minus sign).
   ! Compute the local value of y+. Due to the usage
   ! of pointers there is on offset of -1 in dd2Wall..
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from1)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from2)
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   ad_from4 = bcdata(nn)%jnbeg
   ! We have to inverse the nodal areas
   DO j=ad_from4,bcdata(nn)%jnend
   ad_from3 = bcdata(nn)%inbeg
   DO i=ad_from3,bcdata(nn)%inend
   CALL PUSHREAL8(bcdata(nn)%oarea(i, j))
   bcdata(nn)%oarea(i, j) = one/bcdata(nn)%oarea(i, j)
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from3)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from4)
   CALL PUSHREAL8ARRAY(rlv, SIZE(rlv, 1)*SIZE(rlv, 2)*SIZE(rlv, 3))
   CALL PUSHREAL8ARRAY(rev, SIZE(rev, 1)*SIZE(rev, 2)*SIZE(rev, 3))
   CALL RESETBCPOINTERSBWD(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &                          rev1, rev2, 0)
   CALL PUSHREAL8ARRAY(sk, SIZE(sk, 1)*SIZE(sk, 2)*SIZE(sk, 3)*SIZE(&
   &                   sk, 4))
   CALL PUSHREAL8ARRAY(sj, SIZE(sj, 1)*SIZE(sj, 2)*SIZE(sj, 3)*SIZE(&
   &                   sj, 4))
   CALL PUSHREAL8ARRAY(si, SIZE(si, 1)*SIZE(si, 2)*SIZE(si, 3)*SIZE(&
   &                   si, 4))
   CALL PUSHREAL8ARRAY(d2wall, SIZE(d2wall, 1)*SIZE(d2wall, 2)*SIZE(&
   &                   d2wall, 3))
   CALL PUSHREAL8ARRAY(x, SIZE(x, 1)*SIZE(x, 2)*SIZE(x, 3)*SIZE(x, 4)&
   &                  )
   CALL RESETXXSSRHODD2WALLBWD(nn, xx, ss, rho1, rho2, dd2wall)
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO bocos
   ! Currently the coefficients only contain the surface integral
   ! of the pressure tensor. These values must be scaled to
   ! obtain the correct coefficients.
   CALL PUSHREAL8(fact)
   fact = two/(gammainf*pinf*machcoef*machcoef*surfaceref*lref*lref*&
   &   scaledim)
   CALL PUSHREAL8(cfp(1))
   cfp(1) = cfp(1)*fact
   CALL PUSHREAL8(cfp(2))
   cfp(2) = cfp(2)*fact
   CALL PUSHREAL8(cfv(1))
   cfv(1) = cfv(1)*fact
   CALL PUSHREAL8(cfv(2))
   cfv(2) = cfv(2)*fact
   CALL PUSHREAL8(fact)
   fact = fact/(lengthref*lref)
   CALL PUSHREAL8(cmp(1))
   cmp(1) = cmp(1)*fact
   CALL PUSHREAL8(cmp(2))
   cmp(2) = cmp(2)*fact
   CALL PUSHREAL8(cmv(1))
   cmv(1) = cmv(1)*fact
   CALL PUSHREAL8(cmv(2))
   cmv(2) = cmv(2)*fact
   factd = cmv(3)*cmvd(3)
   cmvd(3) = fact*cmvd(3)
   CALL POPREAL8(cmv(2))
   factd = factd + cmv(2)*cmvd(2)
   cmvd(2) = fact*cmvd(2)
   CALL POPREAL8(cmv(1))
   factd = factd + cmp(3)*cmpd(3) + cmv(1)*cmvd(1)
   cmvd(1) = fact*cmvd(1)
   cmpd(3) = fact*cmpd(3)
   CALL POPREAL8(cmp(2))
   factd = factd + cmp(2)*cmpd(2)
   cmpd(2) = fact*cmpd(2)
   CALL POPREAL8(cmp(1))
   factd = factd + cmp(1)*cmpd(1)
   cmpd(1) = fact*cmpd(1)
   CALL POPREAL8(fact)
   tempd14 = factd/(lref*lengthref)
   lengthrefd = lengthrefd - fact*tempd14/lengthref
   factd = cfv(3)*cfvd(3) + tempd14
   cfvd(3) = fact*cfvd(3)
   CALL POPREAL8(cfv(2))
   factd = factd + cfv(2)*cfvd(2)
   cfvd(2) = fact*cfvd(2)
   CALL POPREAL8(cfv(1))
   factd = factd + cfp(3)*cfpd(3) + cfv(1)*cfvd(1)
   cfvd(1) = fact*cfvd(1)
   cfpd(3) = fact*cfpd(3)
   CALL POPREAL8(cfp(2))
   factd = factd + cfp(2)*cfpd(2)
   cfpd(2) = fact*cfpd(2)
   CALL POPREAL8(cfp(1))
   factd = factd + cfp(1)*cfpd(1)
   cfpd(1) = fact*cfpd(1)
   CALL POPREAL8(fact)
   temp9 = machcoef**2*scaledim
   temp8 = surfaceref*lref**2
   temp7 = temp8*gammainf*pinf
   tempd15 = -(two*factd/(temp7**2*temp9**2))
   tempd16 = temp9*temp8*tempd15
   gammainfd = gammainfd + pinf*tempd16
   pinfd = pinfd + gammainf*tempd16
   machcoefd = machcoefd + scaledim*temp7*2*machcoef*tempd15
   scaledimd = temp7*machcoef**2*tempd15
   pd = 0.0_8
   sid = 0.0_8
   sjd = 0.0_8
   skd = 0.0_8
   DO ii1=1,SIZE(viscsubfaced)
   viscsubfaced(ii1)%tau = 0.0_8
   END DO
   veldirfreestreamd = 0.0_8
   rho1d = 0.0_8
   rho2d = 0.0_8
   vd = 0.0_8
   refpointd = 0.0_8
   xxd = 0.0_8
   ssd = 0.0_8
   DO nn=nbocos,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   CALL POPREAL8ARRAY(x, SIZE(x, 1)*SIZE(x, 2)*SIZE(x, 3)*SIZE(x, 4))
   CALL POPREAL8ARRAY(d2wall, SIZE(d2wall, 1)*SIZE(d2wall, 2)*SIZE(&
   &                  d2wall, 3))
   CALL POPREAL8ARRAY(si, SIZE(si, 1)*SIZE(si, 2)*SIZE(si, 3)*SIZE(si&
   &                  , 4))
   CALL POPREAL8ARRAY(sj, SIZE(sj, 1)*SIZE(sj, 2)*SIZE(sj, 3)*SIZE(sj&
   &                  , 4))
   CALL POPREAL8ARRAY(sk, SIZE(sk, 1)*SIZE(sk, 2)*SIZE(sk, 3)*SIZE(sk&
   &                  , 4))
   CALL RESETXXSSRHODD2WALLBWD_B(nn, xx, xxd, ss, ssd, rho1, rho1d, &
   &                             rho2, rho2d, dd2wall)
   CALL POPREAL8ARRAY(rev, SIZE(rev, 1)*SIZE(rev, 2)*SIZE(rev, 3))
   CALL POPREAL8ARRAY(rlv, SIZE(rlv, 1)*SIZE(rlv, 2)*SIZE(rlv, 3))
   CALL RESETBCPOINTERSBWD_B(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2&
   &                         , pp2d, rlv1, rlv2, rev1, rev2, 0)
   CALL POPINTEGER4(ad_from4)
   CALL POPINTEGER4(ad_to4)
   DO j=ad_to4,ad_from4,-1
   CALL POPINTEGER4(ad_from3)
   CALL POPINTEGER4(ad_to3)
   DO i=ad_to3,ad_from3,-1
   CALL POPREAL8(bcdata(nn)%oarea(i, j))
   temp6 = bcdata(nn)%oarea(i, j)
   bcdatad(nn)%oarea(i, j) = -(one*bcdatad(nn)%oarea(i, j)/temp6&
   &           **2)
   END DO
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   bcdatad(nn)%fv = 0.0_8
   ELSE
   CALL POPINTEGER4(ad_from2)
   CALL POPINTEGER4(ad_to2)
   DO j=ad_to2,ad_from2,-1
   CALL POPINTEGER4(ad_from1)
   CALL POPINTEGER4(ad_to1)
   DO i=ad_to1,ad_from1,-1
   tauxx = viscsubface(nn)%tau(i, j, 1)
   tauyy = viscsubface(nn)%tau(i, j, 2)
   tauxy = viscsubface(nn)%tau(i, j, 4)
   mzd = cmvd(3) + bcdatad(nn)%m(i, j, 3)
   myd = cmvd(2) + bcdatad(nn)%m(i, j, 2)
   mxd = cmvd(1) + bcdatad(nn)%m(i, j, 1)
   xcd = fy*mzd - fz*myd
   ycd = fz*mxd - fx*mzd
   zcd = fx*myd - fy*mxd
   fzd = yc*mxd + bcdatad(nn)%fv(i, j, 3) - xc*myd
   bcdatad(nn)%fv(i, j, 3) = 0.0_8
   fyd = bcdatad(nn)%fv(i, j, 2) - zc*mxd + xc*mzd
   bcdatad(nn)%fv(i, j, 2) = 0.0_8
   fxd = zc*myd + bcdatad(nn)%fv(i, j, 1) - yc*mzd
   bcdatad(nn)%fv(i, j, 1) = 0.0_8
   fzd = fzd + cfvd(3)
   fyd = fyd + cfvd(2)
   fxd = fxd + cfvd(1)
   CALL POPREAL8(zc)
   tempd8 = fourth*zcd
   xxd(i, j, 3) = xxd(i, j, 3) + tempd8
   xxd(i+1, j, 3) = xxd(i+1, j, 3) + tempd8
   xxd(i, j+1, 3) = xxd(i, j+1, 3) + tempd8
   xxd(i+1, j+1, 3) = xxd(i+1, j+1, 3) + tempd8
   refpointd(3) = refpointd(3) - zcd
   CALL POPREAL8(yc)
   tempd9 = fourth*ycd
   xxd(i, j, 2) = xxd(i, j, 2) + tempd9
   xxd(i+1, j, 2) = xxd(i+1, j, 2) + tempd9
   xxd(i, j+1, 2) = xxd(i, j+1, 2) + tempd9
   xxd(i+1, j+1, 2) = xxd(i+1, j+1, 2) + tempd9
   refpointd(2) = refpointd(2) - ycd
   CALL POPREAL8(xc)
   tempd10 = fourth*xcd
   xxd(i, j, 1) = xxd(i, j, 1) + tempd10
   xxd(i+1, j, 1) = xxd(i+1, j, 1) + tempd10
   xxd(i, j+1, 1) = xxd(i, j+1, 1) + tempd10
   xxd(i+1, j+1, 1) = xxd(i+1, j+1, 1) + tempd10
   refpointd(1) = refpointd(1) - xcd
   tauzz = viscsubface(nn)%tau(i, j, 3)
   tauxz = viscsubface(nn)%tau(i, j, 5)
   tauyz = viscsubface(nn)%tau(i, j, 6)
   CALL POPREAL8(fz)
   tempd11 = -(fact*scaledim*fzd)
   ssd(i, j, 1) = ssd(i, j, 1) + tauxz*tempd11
   ssd(i, j, 2) = ssd(i, j, 2) + tauyz*tempd11
   tauzzd = ss(i, j, 3)*tempd11
   ssd(i, j, 3) = ssd(i, j, 3) + tauzz*tempd11
   scaledimd = scaledimd - fact*(tauxy*ss(i, j, 1)+tauyy*ss(i, &
   &             j, 2)+tauyz*ss(i, j, 3))*fyd - fact*(tauxx*ss(i, j, 1)+&
   &             tauxy*ss(i, j, 2)+tauxz*ss(i, j, 3))*fxd - fact*(tauxz*ss(&
   &             i, j, 1)+tauyz*ss(i, j, 2)+tauzz*ss(i, j, 3))*fzd
   CALL POPREAL8(fy)
   tempd13 = -(fact*scaledim*fyd)
   tauyzd = ss(i, j, 3)*tempd13 + ss(i, j, 2)*tempd11
   ssd(i, j, 1) = ssd(i, j, 1) + tauxy*tempd13
   tauyyd = ss(i, j, 2)*tempd13
   ssd(i, j, 2) = ssd(i, j, 2) + tauyy*tempd13
   ssd(i, j, 3) = ssd(i, j, 3) + tauyz*tempd13
   CALL POPREAL8(fx)
   tempd12 = -(fact*scaledim*fxd)
   tauxzd = ss(i, j, 3)*tempd12 + ss(i, j, 1)*tempd11
   tauxyd = ss(i, j, 2)*tempd12 + ss(i, j, 1)*tempd13
   tauxxd = ss(i, j, 1)*tempd12
   ssd(i, j, 1) = ssd(i, j, 1) + tauxx*tempd12
   ssd(i, j, 2) = ssd(i, j, 2) + tauxy*tempd12
   ssd(i, j, 3) = ssd(i, j, 3) + tauxz*tempd12
   viscsubfaced(nn)%tau(i, j, 6) = viscsubfaced(nn)%tau(i, j, 6&
   &             ) + tauyzd
   viscsubfaced(nn)%tau(i, j, 5) = viscsubfaced(nn)%tau(i, j, 5&
   &             ) + tauxzd
   viscsubfaced(nn)%tau(i, j, 4) = viscsubfaced(nn)%tau(i, j, 4&
   &             ) + tauxyd
   viscsubfaced(nn)%tau(i, j, 3) = viscsubfaced(nn)%tau(i, j, 3&
   &             ) + tauzzd
   viscsubfaced(nn)%tau(i, j, 2) = viscsubfaced(nn)%tau(i, j, 2&
   &             ) + tauyyd
   viscsubfaced(nn)%tau(i, j, 1) = viscsubfaced(nn)%tau(i, j, 1&
   &             ) + tauxxd
   END DO
   END DO
   END IF
   CALL POPINTEGER4(ad_from0)
   CALL POPINTEGER4(ad_to0)
   DO j=ad_to0,ad_from0,-1
   CALL POPINTEGER4(ad_from)
   CALL POPINTEGER4(ad_to)
   DO i=ad_to,ad_from,-1
   mzd = bcdatad(nn)%m(i, j, 3)
   bcdatad(nn)%m(i, j, 3) = 0.0_8
   myd = bcdatad(nn)%m(i, j, 2)
   bcdatad(nn)%m(i, j, 2) = 0.0_8
   mxd = bcdatad(nn)%m(i, j, 1)
   bcdatad(nn)%m(i, j, 1) = 0.0_8
   mzd = mzd + cmpd(3)
   myd = myd + cmpd(2)
   mxd = mxd + cmpd(1)
   xcd = fy*mzd - fz*myd
   fyd = cfpd(2) - zc*mxd + xc*mzd
   ycd = fz*mxd - fx*mzd
   fxd = zc*myd + cfpd(1) - yc*mzd
   zcd = fx*myd - fy*mxd
   fzd = yc*mxd + cfpd(3) - xc*myd
   sensor1d = bcdatad(nn)%cavitation(i, j)
   bcdatad(nn)%cavitation(i, j) = 0.0_8
   sensor1d = sensor1d + cavitationd
   qad = four*sensor1*sensor1d
   sensor1d = four*qa*sensor1d
   CALL POPREAL8(sensor1)
   temp5 = -(10*2*sensor1)
   temp4 = one + EXP(temp5)
   sensor1d = EXP(temp5)*one*10*2*sensor1d/temp4**2
   CALL POPREAL8(sensor1)
   cpd = -sensor1d
   tmpd = (plocal-pinf)*cpd
   plocald = tmp*cpd
   temp3 = gammainf*pinf*machcoef**2
   tempd0 = -(two*tmpd/temp3**2)
   tempd = machcoef**2*tempd0
   pinfd = pinfd + gammainf*tempd - tmp*cpd
   CALL POPREAL8(tmp)
   gammainfd = gammainfd + pinf*tempd
   machcoefd = machcoefd + gammainf*pinf*2*machcoef*tempd0
   CALL POPREAL8(plocal)
   sensord = bcdatad(nn)%sepsensor(i, j)
   bcdatad(nn)%sepsensor(i, j) = 0.0_8
   sensord = sensord + sepsensord
   qad = qad + four*sensor*sensord
   sensord = four*qa*sensord
   CALL POPREAL8(sensor)
   temp2 = -(10*2*sensor)
   temp1 = one + EXP(temp2)
   sensord = EXP(temp2)*one*10*2*sensord/temp1**2
   CALL POPREAL8(sensor)
   vd(1) = vd(1) - veldirfreestream(1)*sensord
   veldirfreestreamd(1) = veldirfreestreamd(1) - v(1)*sensord
   vd(2) = vd(2) - veldirfreestream(2)*sensord
   veldirfreestreamd(2) = veldirfreestreamd(2) - v(2)*sensord
   vd(3) = vd(3) - veldirfreestream(3)*sensord
   veldirfreestreamd(3) = veldirfreestreamd(3) - v(3)*sensord
   CALL POPREAL8ARRAY(v, 3)
   tmpd0 = vd
   temp = v(1)**2 + v(2)**2 + v(3)**2
   temp0 = SQRT(temp)
   tempd2 = tmpd0/(temp0+1e-16)
   vd = tempd2
   IF (temp .EQ. 0.0_8) THEN
   tempd3 = 0.0
   ELSE
   tempd3 = SUM(-(v*tempd2/(temp0+1e-16)))/(2.0*temp0)
   END IF
   vd(1) = vd(1) + 2*v(1)*tempd3
   vd(2) = vd(2) + 2*v(2)*tempd3
   vd(3) = vd(3) + 2*v(3)*tempd3
   CALL POPREAL8(v(3))
   ww2d(i, j, ivz) = ww2d(i, j, ivz) + vd(3)
   vd(3) = 0.0_8
   CALL POPREAL8(v(2))
   ww2d(i, j, ivy) = ww2d(i, j, ivy) + vd(2)
   vd(2) = 0.0_8
   CALL POPREAL8(v(1))
   ww2d(i, j, ivx) = ww2d(i, j, ivx) + vd(1)
   vd(1) = 0.0_8
   CALL POPREAL8(bcdata(nn)%oarea(i, j))
   qad = qad + bcdatad(nn)%oarea(i-1, j) + bcdatad(nn)%oarea(i-1&
   &           , j-1) + bcdatad(nn)%oarea(i, j-1) + bcdatad(nn)%oarea(i, j)
   CALL POPREAL8(bcdata(nn)%oarea(i-1, j))
   CALL POPREAL8(bcdata(nn)%oarea(i, j-1))
   CALL POPREAL8(bcdata(nn)%oarea(i-1, j-1))
   CALL POPREAL8(qa)
   IF (ss(i, j, 1)**2 + ss(i, j, 2)**2 + ss(i, j, 3)**2 .EQ. &
   &             0.0_8) THEN
   tempd4 = 0.0
   ELSE
   tempd4 = fourth*qad/(2.0*SQRT(ss(i, j, 1)**2+ss(i, j, 2)**2+&
   &             ss(i, j, 3)**2))
   END IF
   ssd(i, j, 1) = ssd(i, j, 1) + 2*ss(i, j, 1)*tempd4
   ssd(i, j, 2) = ssd(i, j, 2) + 2*ss(i, j, 2)*tempd4
   fzd = fzd + bcdatad(nn)%fp(i, j, 3)
   ssd(i, j, 3) = ssd(i, j, 3) + pm1*fzd + 2*ss(i, j, 3)*tempd4
   bcdatad(nn)%fp(i, j, 3) = 0.0_8
   fyd = fyd + bcdatad(nn)%fp(i, j, 2)
   bcdatad(nn)%fp(i, j, 2) = 0.0_8
   fxd = fxd + bcdatad(nn)%fp(i, j, 1)
   bcdatad(nn)%fp(i, j, 1) = 0.0_8
   CALL POPREAL8(fz)
   pm1d = ss(i, j, 2)*fyd + ss(i, j, 1)*fxd + ss(i, j, 3)*fzd
   CALL POPREAL8(fy)
   ssd(i, j, 2) = ssd(i, j, 2) + pm1*fyd
   CALL POPREAL8(fx)
   ssd(i, j, 1) = ssd(i, j, 1) + pm1*fxd
   CALL POPREAL8(zc)
   tempd5 = fourth*zcd
   xxd(i, j, 3) = xxd(i, j, 3) + tempd5
   xxd(i+1, j, 3) = xxd(i+1, j, 3) + tempd5
   xxd(i, j+1, 3) = xxd(i, j+1, 3) + tempd5
   xxd(i+1, j+1, 3) = xxd(i+1, j+1, 3) + tempd5
   refpointd(3) = refpointd(3) - zcd
   CALL POPREAL8(yc)
   tempd6 = fourth*ycd
   xxd(i, j, 2) = xxd(i, j, 2) + tempd6
   xxd(i+1, j, 2) = xxd(i+1, j, 2) + tempd6
   xxd(i, j+1, 2) = xxd(i, j+1, 2) + tempd6
   xxd(i+1, j+1, 2) = xxd(i+1, j+1, 2) + tempd6
   refpointd(2) = refpointd(2) - ycd
   CALL POPREAL8(xc)
   tempd7 = fourth*xcd
   xxd(i, j, 1) = xxd(i, j, 1) + tempd7
   xxd(i+1, j, 1) = xxd(i+1, j, 1) + tempd7
   xxd(i, j+1, 1) = xxd(i, j+1, 1) + tempd7
   xxd(i+1, j+1, 1) = xxd(i+1, j+1, 1) + tempd7
   refpointd(1) = refpointd(1) - xcd
   CALL POPREAL8(pm1)
   tempd1 = fact*scaledim*pm1d
   pp2d(i, j) = pp2d(i, j) + half*tempd1 + plocald
   pp1d(i, j) = pp1d(i, j) + half*tempd1
   pinfd = pinfd - tempd1
   scaledimd = scaledimd + fact*(half*(pp2(i, j)+pp1(i, j))-pinf)&
   &           *pm1d
   END DO
   END DO
   CALL POPREAL8ARRAY(bcdata(nn)%oarea, SIZE(bcdata(nn)%oarea, 1)*&
   &                  SIZE(bcdata(nn)%oarea, 2))
   bcdatad(nn)%oarea = 0.0_8
   CALL POPCONTROL3B(branch)
   IF (branch .LT. 3) THEN
   IF (branch .NE. 0) THEN
   IF (branch .EQ. 1) THEN
   CALL POPREAL8(fact)
   ELSE
   CALL POPREAL8(fact)
   END IF
   END IF
   ELSE IF (branch .LT. 5) THEN
   IF (branch .EQ. 3) THEN
   CALL POPREAL8(fact)
   ELSE
   CALL POPREAL8(fact)
   END IF
   ELSE IF (branch .EQ. 5) THEN
   CALL POPREAL8(fact)
   ELSE
   CALL POPREAL8(fact)
   END IF
   CALL POPREAL8ARRAY(ss, imaxdim*jmaxdim*3)
   CALL SETXXSSRHODD2WALLBWD_B(nn, xx, xxd, ss, ssd, rho1, rho1d, &
   &                           rho2, rho2d, dd2wall)
   CALL POPREAL8ARRAY(pp1, imaxdim*jmaxdim)
   CALL POPREAL8ARRAY(pp2, imaxdim*jmaxdim)
   CALL SETBCPOINTERSBWD_B(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                       pp2d, rlv1, rlv2, rev1, rev2, 0)
   END IF
   END DO
   pointrefd(3) = pointrefd(3) + lref*refpointd(3)
   refpointd(3) = 0.0_8
   pointrefd(2) = pointrefd(2) + lref*refpointd(2)
   refpointd(2) = 0.0_8
   pointrefd(1) = pointrefd(1) + lref*refpointd(1)
   prefd = prefd + scaledimd/pinf
   pinfd = pinfd - pref*scaledimd/pinf**2
   END SUBROUTINE FORCESANDMOMENTS_B

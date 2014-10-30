   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of bceulerwall in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *rev *p *gamma *w *rlv
   !   with respect to varying inputs: tref rgas *rev *p *s *gamma
   !                *w *rlv *si *sj *sk *(*bcdata.norm) *(*bcdata.rface)
   !   Plus diff mem management of: rev:in p:in gamma:in w:in rlv:in
   !                si:in sj:in sk:in bcdata:in *bcdata.norm:in *bcdata.rface:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcEulerWall.f90                                 *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-07-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCEULERWALL_D(secondhalo, correctfork)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcEulerWall applies the inviscid wall boundary condition to    *
   !      * a block. It is assumed that the pointers in blockPointers are  *
   !      * already set to the correct block on the correct grid level.    *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_D
   USE BCTYPES
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE INPUTDISCRETIZATION
   USE INPUTPHYSICS
   USE ITERATION
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   LOGICAL, INTENT(IN) :: secondhalo, correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, j, k, l
   INTEGER(kind=inttype) :: jm1, jp1, km1, kp1
   INTEGER(kind=inttype) :: walltreatment
   REAL(kind=realtype) :: sixa, siya, siza, sjxa, sjya, sjza
   REAL(kind=realtype) :: sixad, siyad, sizad, sjxad, sjyad, sjzad
   REAL(kind=realtype) :: skxa, skya, skza, a1, b1
   REAL(kind=realtype) :: skxad, skyad, skzad
   REAL(kind=realtype) :: rxj, ryj, rzj, rxk, ryk, rzk
   REAL(kind=realtype) :: rxjd, ryjd, rzjd, rxkd, rykd, rzkd
   REAL(kind=realtype) :: dpj, dpk, ri, rj, rk, qj, qk, vn
   REAL(kind=realtype) :: dpjd, dpkd, rid, rjd, rkd, qjd, qkd, vnd
   REAL(kind=realtype) :: ux, uy, uz
   REAL(kind=realtype) :: uxd, uyd, uzd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3d, pp4d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1d, rev2d
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssi, ssj, ssk
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssid, ssjd, sskd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssd
   INTERFACE 
   SUBROUTINE SETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &       rev1, rev2, offset)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS
   SUBROUTINE RESETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &       rev1, rev2, offset)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE RESETBCPOINTERS
   SUBROUTINE SETPP3PP4(nn, pp3, pp4)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   END SUBROUTINE SETPP3PP4
   SUBROUTINE RESETPP3PP4(nn, pp3, pp4)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   END SUBROUTINE RESETPP3PP4
   SUBROUTINE SETSS(nn, ssi, ssj, ssk, ss)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssi, ssj, &
   &       ssk
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss
   END SUBROUTINE SETSS
   SUBROUTINE RESETSS(nn, ssi, ssj, ssk, ss)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssi, ssj, &
   &       ssk
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss
   END SUBROUTINE RESETSS
   END INTERFACE
      INTERFACE 
   SUBROUTINE SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
   &       pp2, pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev1d, rev2, rev2d, &
   &       offset)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1d, rev2d
   END SUBROUTINE SETBCPOINTERS_D
   SUBROUTINE SETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3d, pp4d
   END SUBROUTINE SETPP3PP4_D
   SUBROUTINE SETSS_D(nn, ssi, ssid, ssj, ssjd, ssk, sskd, ss, ssd)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssi, ssj, &
   &       ssk
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssid, ssjd, &
   &       sskd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssd
   END SUBROUTINE SETSS_D
   END INTERFACE
      INTRINSIC MAX
   INTRINSIC MIN
   REAL(kind=realtype) :: DIM
   REAL(kind=realtype) :: DIM_D
   REAL(kind=realtype) :: tmpresult
   INTEGER(kind=inttype) :: max2
   INTEGER(kind=inttype) :: max1
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Make sure that on the coarser grids the constant pressure
   ! boundary condition is used.
   walltreatment = wallbctreatment
   IF (currentlevel .GT. groundlevel) walltreatment = constantpressure
   ! Loop over the boundary condition subfaces of this block.
   bocos:DO nn=1,nbocos
   ! Check for Euler wall boundary condition.
   IF (bctype(nn) .EQ. eulerwall) THEN
   ! Set the pointers for the unit normal and the normal
   ! velocity to make the code more readable.
   ! Modify to use actual pointer - Peter Lyu
   !norm  => BCData(nn)%norm
   !rface => BCData(nn)%rface
   ! Nullify the pointers and set them to the correct subface.
   ! They are nullified first, because some compilers require
   ! that.
   !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
   CALL SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                    pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev1d, rev2, &
   &                    rev2d, 0)
   !
   !          **************************************************************
   !          *                                                            *
   !          * Determine the boundary condition treatment and compute the *
   !          * undivided pressure gradient accordingly. This gradient is  *
   !          * temporarily stored in the halo pressure.                   *
   !          *                                                            *
   !          **************************************************************
   !
   SELECT CASE  (walltreatment) 
   CASE (constantpressure) 
   ! Constant pressure. Set the gradient to zero.
   DO k=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO j=bcdata(nn)%icbeg,bcdata(nn)%icend
   pp1d(j, k) = 0.0_8
   pp1(j, k) = zero
   END DO
   END DO
   CASE (linextrapolpressure) 
   !===========================================================
   ! Linear extrapolation. First set the additional pointer
   ! for pp3, depending on the block face.
   CALL SETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
   ! Compute the gradient.
   DO k=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO j=bcdata(nn)%icbeg,bcdata(nn)%icend
   pp1d(j, k) = pp3d(j, k) - pp2d(j, k)
   pp1(j, k) = pp3(j, k) - pp2(j, k)
   END DO
   END DO
   CALL RESETPP3PP4(nn, pp3, pp4)
   CASE (quadextrapolpressure) 
   !===========================================================
   ! Quadratic extrapolation. First set the additional
   ! pointers for pp3 and pp4, depending on the block face.
   CALL SETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
   ! Compute the gradient.
   DO k=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO j=bcdata(nn)%icbeg,bcdata(nn)%icend
   pp1d(j, k) = two*pp3d(j, k) - 1.5_realType*pp2d(j, k) - half&
   &             *pp4d(j, k)
   pp1(j, k) = two*pp3(j, k) - 1.5_realType*pp2(j, k) - half*&
   &             pp4(j, k)
   END DO
   END DO
   CALL RESETPP3PP4(nn, pp3, pp4)
   CASE (normalmomentum) 
   !===========================================================
   ! Pressure gradient is computed using the normal momentum
   ! equation. First set a couple of additional variables for
   ! the normals, depending on the block face. Note that the
   ! construction 1: should not be used in these pointers,
   ! because element 0 is needed. Consequently there will be
   ! an offset of 1 for these normals. This is commented in
   ! the code. For moving faces also the grid velocity of
   ! the 1st cell center from the wall is needed.
   CALL SETSS_D(nn, ssi, ssid, ssj, ssjd, ssk, sskd, ss, ssd)
   ! Loop over the faces of the generic subface.
   DO k=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   ! Store the indices k+1, k-1 a bit easier and make
   ! sure that they do not exceed the range of the arrays.
   km1 = k - 1
   IF (bcdata(nn)%jcbeg .LT. km1) THEN
   km1 = km1
   ELSE
   km1 = bcdata(nn)%jcbeg
   END IF
   kp1 = k + 1
   IF (bcdata(nn)%jcend .GT. kp1) THEN
   kp1 = kp1
   ELSE
   kp1 = bcdata(nn)%jcend
   END IF
   IF (1_intType .LT. kp1 - km1) THEN
   max1 = kp1 - km1
   ELSE
   max1 = 1_intType
   END IF
   ! Compute the scaling factor for the central difference
   ! in the k-direction.
   b1 = one/max1
   ! The j-loop.
   DO j=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! The indices j+1 and j-1. Make sure that they
   ! do not exceed the range of the arrays.
   jm1 = j - 1
   IF (bcdata(nn)%icbeg .LT. jm1) THEN
   jm1 = jm1
   ELSE
   jm1 = bcdata(nn)%icbeg
   END IF
   jp1 = j + 1
   IF (bcdata(nn)%icend .GT. jp1) THEN
   jp1 = jp1
   ELSE
   jp1 = bcdata(nn)%icend
   END IF
   IF (1_intType .LT. jp1 - jm1) THEN
   max2 = jp1 - jm1
   ELSE
   max2 = 1_intType
   END IF
   ! Compute the scaling factor for the central
   ! difference in the j-direction.
   a1 = one/max2
   ! Compute (twice) the average normal in the generic i,
   ! j and k-direction. Note that in j and k-direction
   ! the average in the original indices should be taken
   ! using j-1 and j (and k-1 and k). However due to the
   ! usage of pointers ssj and ssk there is an offset in
   ! the indices of 1 and therefore now the correct
   ! average is obtained with the indices j and j+1
   ! (k and k+1).
   sixad = two*ssid(j, k, 1)
   sixa = two*ssi(j, k, 1)
   siyad = two*ssid(j, k, 2)
   siya = two*ssi(j, k, 2)
   sizad = two*ssid(j, k, 3)
   siza = two*ssi(j, k, 3)
   sjxad = ssjd(j, k, 1) + ssjd(j+1, k, 1)
   sjxa = ssj(j, k, 1) + ssj(j+1, k, 1)
   sjyad = ssjd(j, k, 2) + ssjd(j+1, k, 2)
   sjya = ssj(j, k, 2) + ssj(j+1, k, 2)
   sjzad = ssjd(j, k, 3) + ssjd(j+1, k, 3)
   sjza = ssj(j, k, 3) + ssj(j+1, k, 3)
   skxad = sskd(j, k, 1) + sskd(j, k+1, 1)
   skxa = ssk(j, k, 1) + ssk(j, k+1, 1)
   skyad = sskd(j, k, 2) + sskd(j, k+1, 2)
   skya = ssk(j, k, 2) + ssk(j, k+1, 2)
   skzad = sskd(j, k, 3) + sskd(j, k+1, 3)
   skza = ssk(j, k, 3) + ssk(j, k+1, 3)
   ! Compute the difference of the normal vector and
   ! pressure in j and k-direction. As the indices are
   ! restricted to the 1st halo-layer, the computation
   ! of the internal halo values is not consistent;
   ! however this is not really a problem, because these
   ! values are overwritten in the communication pattern.
   rxjd = a1*(bcdatad(nn)%norm(jp1, k, 1)-bcdatad(nn)%norm(jm1&
   &             , k, 1))
   rxj = a1*(bcdata(nn)%norm(jp1, k, 1)-bcdata(nn)%norm(jm1, k&
   &             , 1))
   ryjd = a1*(bcdatad(nn)%norm(jp1, k, 2)-bcdatad(nn)%norm(jm1&
   &             , k, 2))
   ryj = a1*(bcdata(nn)%norm(jp1, k, 2)-bcdata(nn)%norm(jm1, k&
   &             , 2))
   rzjd = a1*(bcdatad(nn)%norm(jp1, k, 3)-bcdatad(nn)%norm(jm1&
   &             , k, 3))
   rzj = a1*(bcdata(nn)%norm(jp1, k, 3)-bcdata(nn)%norm(jm1, k&
   &             , 3))
   dpjd = a1*(pp2d(jp1, k)-pp2d(jm1, k))
   dpj = a1*(pp2(jp1, k)-pp2(jm1, k))
   rxkd = b1*(bcdatad(nn)%norm(j, kp1, 1)-bcdatad(nn)%norm(j, &
   &             km1, 1))
   rxk = b1*(bcdata(nn)%norm(j, kp1, 1)-bcdata(nn)%norm(j, km1&
   &             , 1))
   rykd = b1*(bcdatad(nn)%norm(j, kp1, 2)-bcdatad(nn)%norm(j, &
   &             km1, 2))
   ryk = b1*(bcdata(nn)%norm(j, kp1, 2)-bcdata(nn)%norm(j, km1&
   &             , 2))
   rzkd = b1*(bcdatad(nn)%norm(j, kp1, 3)-bcdatad(nn)%norm(j, &
   &             km1, 3))
   rzk = b1*(bcdata(nn)%norm(j, kp1, 3)-bcdata(nn)%norm(j, km1&
   &             , 3))
   dpkd = b1*(pp2d(j, kp1)-pp2d(j, km1))
   dpk = b1*(pp2(j, kp1)-pp2(j, km1))
   ! Compute the dot product between the unit vector
   ! and the normal vectors in i, j and k-direction.
   rid = bcdatad(nn)%norm(j, k, 1)*sixa + bcdata(nn)%norm(j, k&
   &             , 1)*sixad + bcdatad(nn)%norm(j, k, 2)*siya + bcdata(nn)%&
   &             norm(j, k, 2)*siyad + bcdatad(nn)%norm(j, k, 3)*siza + &
   &             bcdata(nn)%norm(j, k, 3)*sizad
   ri = bcdata(nn)%norm(j, k, 1)*sixa + bcdata(nn)%norm(j, k, 2&
   &             )*siya + bcdata(nn)%norm(j, k, 3)*siza
   rjd = bcdatad(nn)%norm(j, k, 1)*sjxa + bcdata(nn)%norm(j, k&
   &             , 1)*sjxad + bcdatad(nn)%norm(j, k, 2)*sjya + bcdata(nn)%&
   &             norm(j, k, 2)*sjyad + bcdatad(nn)%norm(j, k, 3)*sjza + &
   &             bcdata(nn)%norm(j, k, 3)*sjzad
   rj = bcdata(nn)%norm(j, k, 1)*sjxa + bcdata(nn)%norm(j, k, 2&
   &             )*sjya + bcdata(nn)%norm(j, k, 3)*sjza
   rkd = bcdatad(nn)%norm(j, k, 1)*skxa + bcdata(nn)%norm(j, k&
   &             , 1)*skxad + bcdatad(nn)%norm(j, k, 2)*skya + bcdata(nn)%&
   &             norm(j, k, 2)*skyad + bcdatad(nn)%norm(j, k, 3)*skza + &
   &             bcdata(nn)%norm(j, k, 3)*skzad
   rk = bcdata(nn)%norm(j, k, 1)*skxa + bcdata(nn)%norm(j, k, 2&
   &             )*skya + bcdata(nn)%norm(j, k, 3)*skza
   ! Store the velocity components in ux, uy and uz and
   ! subtract the mesh velocity if the face is moving.
   uxd = ww2d(j, k, ivx)
   ux = ww2(j, k, ivx)
   uyd = ww2d(j, k, ivy)
   uy = ww2(j, k, ivy)
   uzd = ww2d(j, k, ivz)
   uz = ww2(j, k, ivz)
   IF (addgridvelocities) THEN
   uxd = uxd - ssd(j, k, 1)
   ux = ux - ss(j, k, 1)
   uyd = uyd - ssd(j, k, 2)
   uy = uy - ss(j, k, 2)
   uzd = uzd - ssd(j, k, 3)
   uz = uz - ss(j, k, 3)
   END IF
   ! Compute the velocity components in j and
   ! k-direction.
   qjd = uxd*sjxa + ux*sjxad + uyd*sjya + uy*sjyad + uzd*sjza +&
   &             uz*sjzad
   qj = ux*sjxa + uy*sjya + uz*sjza
   qkd = uxd*skxa + ux*skxad + uyd*skya + uy*skyad + uzd*skza +&
   &             uz*skzad
   qk = ux*skxa + uy*skya + uz*skza
   ! Compute the pressure gradient, which is stored
   ! in pp1. I'm not entirely sure whether this
   ! formulation is correct for moving meshes. It could
   ! be that an additional term is needed there.
   pp1d(j, k) = (((qjd*(ux*rxj+uy*ryj+uz*rzj)+qj*(uxd*rxj+ux*&
   &             rxjd+uyd*ryj+uy*ryjd+uzd*rzj+uz*rzjd)+qkd*(ux*rxk+uy*ryk+&
   &             uz*rzk)+qk*(uxd*rxk+ux*rxkd+uyd*ryk+uy*rykd+uzd*rzk+uz*&
   &             rzkd))*ww2(j, k, irho)+(qj*(ux*rxj+uy*ryj+uz*rzj)+qk*(ux*&
   &             rxk+uy*ryk+uz*rzk))*ww2d(j, k, irho)-rjd*dpj-rj*dpjd-rkd*&
   &             dpk-rk*dpkd)*ri-((qj*(ux*rxj+uy*ryj+uz*rzj)+qk*(ux*rxk+uy*&
   &             ryk+uz*rzk))*ww2(j, k, irho)-rj*dpj-rk*dpk)*rid)/ri**2
   pp1(j, k) = ((qj*(ux*rxj+uy*ryj+uz*rzj)+qk*(ux*rxk+uy*ryk+uz&
   &             *rzk))*ww2(j, k, irho)-rj*dpj-rk*dpk)/ri
   END DO
   END DO
   CALL RESETSS(nn, ssi, ssj, ssk, ss)
   END SELECT
   ! Determine the state in the halo cell. Again loop over
   ! the cell range for this subface.
   DO k=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO j=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Compute the pressure density and velocity in the
   ! halo cell. Note that rface is the grid velocity
   ! component in the direction of norm, i.e. outward
   ! pointing.
   pp1d(j, k) = DIM_D(pp2(j, k), pp2d(j, k), pp1(j, k), pp1d(j, k&
   &           ), tmpresult)
   pp1(j, k) = tmpresult
   vnd = two*(bcdatad(nn)%rface(j, k)-ww2d(j, k, ivx)*bcdata(nn)%&
   &           norm(j, k, 1)-ww2(j, k, ivx)*bcdatad(nn)%norm(j, k, 1)-ww2d(&
   &           j, k, ivy)*bcdata(nn)%norm(j, k, 2)-ww2(j, k, ivy)*bcdatad(&
   &           nn)%norm(j, k, 2)-ww2d(j, k, ivz)*bcdata(nn)%norm(j, k, 3)-&
   &           ww2(j, k, ivz)*bcdatad(nn)%norm(j, k, 3))
   vn = two*(bcdata(nn)%rface(j, k)-ww2(j, k, ivx)*bcdata(nn)%&
   &           norm(j, k, 1)-ww2(j, k, ivy)*bcdata(nn)%norm(j, k, 2)-ww2(j&
   &           , k, ivz)*bcdata(nn)%norm(j, k, 3))
   ww1d(j, k, irho) = ww2d(j, k, irho)
   ww1(j, k, irho) = ww2(j, k, irho)
   ww1d(j, k, ivx) = ww2d(j, k, ivx) + vnd*bcdata(nn)%norm(j, k, &
   &           1) + vn*bcdatad(nn)%norm(j, k, 1)
   ww1(j, k, ivx) = ww2(j, k, ivx) + vn*bcdata(nn)%norm(j, k, 1)
   ww1d(j, k, ivy) = ww2d(j, k, ivy) + vnd*bcdata(nn)%norm(j, k, &
   &           2) + vn*bcdatad(nn)%norm(j, k, 2)
   ww1(j, k, ivy) = ww2(j, k, ivy) + vn*bcdata(nn)%norm(j, k, 2)
   ww1d(j, k, ivz) = ww2d(j, k, ivz) + vnd*bcdata(nn)%norm(j, k, &
   &           3) + vn*bcdatad(nn)%norm(j, k, 3)
   ww1(j, k, ivz) = ww2(j, k, ivz) + vn*bcdata(nn)%norm(j, k, 3)
   ! Just copy the turbulent variables.
   DO l=nt1mg,nt2mg
   ww1d(j, k, l) = ww2d(j, k, l)
   ww1(j, k, l) = ww2(j, k, l)
   END DO
   ! The laminar and eddy viscosity, if present.
   IF (viscous) THEN
   rlv1d(j, k) = rlv2d(j, k)
   rlv1(j, k) = rlv2(j, k)
   END IF
   IF (eddymodel) THEN
   rev1d(j, k) = rev2d(j, k)
   rev1(j, k) = rev2(j, k)
   END IF
   END DO
   END DO
   ! deallocation all pointer
   CALL RESETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, &
   &                       rev2, 0)
   ! Compute the energy for these halo's.
   CALL COMPUTEETOT_D(icbeg(nn), icend(nn), jcbeg(nn), jcend(nn), &
   &                  kcbeg(nn), kcend(nn), correctfork)
   ! Extrapolate the state vectors in case a second halo
   ! is needed.
   IF (secondhalo) CALL EXTRAPOLATE2NDHALO_D(nn, correctfork)
   END IF
   END DO bocos
   END SUBROUTINE BCEULERWALL_D

   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of bceulerwall in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *p *w *rlv
   !   with respect to varying inputs: *p *w *rlv
   !   Plus diff mem management of: p:in w:in rlv:in
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
   REAL(kind=realtype) :: skxa, skya, skza, a1, b1
   REAL(kind=realtype) :: rxj, ryj, rzj, rxk, ryk, rzk
   REAL(kind=realtype) :: dpj, dpk, ri, rj, rk, qj, qk, vn
   REAL(kind=realtype) :: vnd
   REAL(kind=realtype) :: ux, uy, uz
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3d, pp4d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ssi, ssj, ssk
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss
   INTERFACE 
   SUBROUTINE SETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &       rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS
   SUBROUTINE RESETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &       rev1, rev2, offset)
   USE BLOCKPOINTERS_D
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
   &       pp2, pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS_D
   SUBROUTINE RESETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
   &       pp2, pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE RESETBCPOINTERS_D
   SUBROUTINE SETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3d, pp4d
   END SUBROUTINE SETPP3PP4_D
   SUBROUTINE RESETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
   USE BCTYPES
   USE BLOCKPOINTERS_D
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3, pp4
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp3d, pp4d
   END SUBROUTINE RESETPP3PP4_D
   END INTERFACE
      REAL(kind=realtype) :: DIM
   REAL(kind=realtype) :: DIM_D
   REAL(kind=realtype) :: tmpresult
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
   !norm  => BCData(nn)%norm
   !rface => BCData(nn)%rface
   ! Nullify the pointers and set them to the correct subface.
   ! They are nullified first, because some compilers require
   ! that.
   !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
   CALL SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                    pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev2, 0)
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
   CALL RESETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
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
   CALL RESETPP3PP4_D(nn, pp3, pp3d, pp4, pp4d)
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
   CALL SETSS(nn, ssi, ssj, ssk, ss)
   !!$
   !!$               ! Loop over the faces of the generic subface.
   !!$
   !!$               do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd
   !!$
   !!$                 ! Store the indices k+1, k-1 a bit easier and make
   !!$                 ! sure that they do not exceed the range of the arrays.
   !!$
   !!$                 km1 = k-1; km1 = max(BCData(nn)%jcBeg,km1)
   !!$                 kp1 = k+1; kp1 = min(BCData(nn)%jcEnd,kp1)
   !!$
   !!$                 ! Compute the scaling factor for the central difference
   !!$                 ! in the k-direction.
   !!$
   !!$                 b1 = one/max(1_intType,(kp1-km1))
   !!$
   !!$                 ! The j-loop.
   !!$
   !!$                 do j=BCData(nn)%icBeg, BCData(nn)%icEnd
   !!$
   !!$                   ! The indices j+1 and j-1. Make sure that they
   !!$                   ! do not exceed the range of the arrays.
   !!$
   !!$                   jm1 = j-1; jm1 = max(BCData(nn)%icBeg,jm1)
   !!$                   jp1 = j+1; jp1 = min(BCData(nn)%icEnd,jp1)
   !!$
   !!$                   ! Compute the scaling factor for the central
   !!$                   ! difference in the j-direction.
   !!$
   !!$                   a1 = one/max(1_intType,(jp1-jm1))
   !!$
   !!$                   ! Compute (twice) the average normal in the generic i,
   !!$                   ! j and k-direction. Note that in j and k-direction
   !!$                   ! the average in the original indices should be taken
   !!$                   ! using j-1 and j (and k-1 and k). However due to the
   !!$                   ! usage of pointers ssj and ssk there is an offset in
   !!$                   ! the indices of 1 and therefore now the correct
   !!$                   ! average is obtained with the indices j and j+1
   !!$                   ! (k and k+1).
   !!$
   !!$                   sixa = two*ssi(j,k,1)
   !!$                   siya = two*ssi(j,k,2)
   !!$                   siza = two*ssi(j,k,3)
   !!$
   !!$                   sjxa = ssj(j,k,1) + ssj(j+1,k,1)
   !!$                   sjya = ssj(j,k,2) + ssj(j+1,k,2)
   !!$                   sjza = ssj(j,k,3) + ssj(j+1,k,3)
   !!$
   !!$                   skxa = ssk(j,k,1) + ssk(j,k+1,1)
   !!$                   skya = ssk(j,k,2) + ssk(j,k+1,2)
   !!$                   skza = ssk(j,k,3) + ssk(j,k+1,3)
   !!$
   !!$                   ! Compute the difference of the normal vector and
   !!$                   ! pressure in j and k-direction. As the indices are
   !!$                   ! restricted to the 1st halo-layer, the computation
   !!$                   ! of the internal halo values is not consistent;
   !!$                   ! however this is not really a problem, because these
   !!$                   ! values are overwritten in the communication pattern.
   !!$
   !!$                   rxj = a1*(norm(jp1,k,1) - norm(jm1,k,1))
   !!$                   ryj = a1*(norm(jp1,k,2) - norm(jm1,k,2))
   !!$                   rzj = a1*(norm(jp1,k,3) - norm(jm1,k,3))
   !!$                   dpj = a1*(pp2(jp1,k)    - pp2(jm1,k))
   !!$
   !!$                   rxk = b1*(norm(j,kp1,1) - norm(j,km1,1))
   !!$                   ryk = b1*(norm(j,kp1,2) - norm(j,km1,2))
   !!$                   rzk = b1*(norm(j,kp1,3) - norm(j,km1,3))
   !!$                   dpk = b1*(pp2(j,kp1)    - pp2(j,km1))
   !!$
   !!$                   ! Compute the dot product between the unit vector
   !!$                   ! and the normal vectors in i, j and k-direction.
   !!$
   !!$                   ri = norm(j,k,1)*sixa + norm(j,k,2)*siya &
   !!$                      + norm(j,k,3)*siza
   !!$                   rj = norm(j,k,1)*sjxa + norm(j,k,2)*sjya &
   !!$                      + norm(j,k,3)*sjza
   !!$                   rk = norm(j,k,1)*skxa + norm(j,k,2)*skya &
   !!$                      + norm(j,k,3)*skza
   !!$
   !!$                   ! Store the velocity components in ux, uy and uz and
   !!$                   ! subtract the mesh velocity if the face is moving.
   !!$
   !!$                   ux = ww2(j,k,ivx)
   !!$                   uy = ww2(j,k,ivy)
   !!$                   uz = ww2(j,k,ivz)
   !!$
   !!$                   if( addGridVelocities ) then
   !!$                     ux = ux - ss(j,k,1)
   !!$                     uy = uy - ss(j,k,2)
   !!$                     uz = uz - ss(j,k,3)
   !!$                   endif
   !!$
   !!$                   ! Compute the velocity components in j and
   !!$                   ! k-direction.
   !!$
   !!$                   qj = ux*sjxa + uy*sjya + uz*sjza
   !!$                   qk = ux*skxa + uy*skya + uz*skza
   !!$
   !!$                   ! Compute the pressure gradient, which is stored
   !!$                   ! in pp1. I'm not entirely sure whether this
   !!$                   ! formulation is correct for moving meshes. It could
   !!$                   ! be that an additional term is needed there.
   !!$
   !!$                   pp1(j,k) = ((qj*(ux*rxj + uy*ryj + uz*rzj)      &
   !!$                            +   qk*(ux*rxk + uy*ryk + uz*rzk))     &
   !!$                            *  ww2(j,k,irho) - rj*dpj - rk*dpk)/ri
   !!$                 enddo
   !!$               enddo
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
   vnd = two*(-(bcdata(nn)%norm(j, k, 1)*ww2d(j, k, ivx))-bcdata(&
   &           nn)%norm(j, k, 2)*ww2d(j, k, ivy)-bcdata(nn)%norm(j, k, 3)*&
   &           ww2d(j, k, ivz))
   vn = two*(bcdata(nn)%rface(j, k)-ww2(j, k, ivx)*bcdata(nn)%&
   &           norm(j, k, 1)-ww2(j, k, ivy)*bcdata(nn)%norm(j, k, 2)-ww2(j&
   &           , k, ivz)*bcdata(nn)%norm(j, k, 3))
   ww1d(j, k, irho) = ww2d(j, k, irho)
   ww1(j, k, irho) = ww2(j, k, irho)
   ww1d(j, k, ivx) = ww2d(j, k, ivx) + bcdata(nn)%norm(j, k, 1)*&
   &           vnd
   ww1(j, k, ivx) = ww2(j, k, ivx) + vn*bcdata(nn)%norm(j, k, 1)
   ww1d(j, k, ivy) = ww2d(j, k, ivy) + bcdata(nn)%norm(j, k, 2)*&
   &           vnd
   ww1(j, k, ivy) = ww2(j, k, ivy) + vn*bcdata(nn)%norm(j, k, 2)
   ww1d(j, k, ivz) = ww2d(j, k, ivz) + bcdata(nn)%norm(j, k, 3)*&
   &           vnd
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
   IF (eddymodel) rev1(j, k) = rev2(j, k)
   END DO
   END DO
   ! Compute the energy for these halo's.
   CALL COMPUTEETOT_D(icbeg(nn), icend(nn), jcbeg(nn), jcend(nn), &
   &                  kcbeg(nn), kcend(nn), correctfork)
   ! Extrapolate the state vectors in case a second halo
   ! is needed.
   IF (secondhalo) CALL EXTRAPOLATE2NDHALO_D(nn, correctfork)
   ! deallocation all pointer
   CALL RESETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                      pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev2, 0)
   END IF
   END DO bocos
   END SUBROUTINE BCEULERWALL_D

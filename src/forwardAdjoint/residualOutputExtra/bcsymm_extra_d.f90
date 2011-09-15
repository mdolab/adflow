   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of bcsymm in forward (tangent) mode:
   !   variations   of useful results: *p
   !   with respect to varying inputs: *p
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcSymm.f90                                      *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-07-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCSYMM_EXTRA_D(secondhalo)
   USE FLOWVARREFSTATE
   USE BCTYPES
   USE BLOCKPOINTERS_D
   USE ITERATION
   USE CONSTANTS
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcSymm applies the symmetry boundary conditions to a block.    *
   !      * It is assumed that the pointers in blockPointers are already   *
   !      * set to the correct block on the correct grid level.            *
   !      *                                                                *
   !      * In case also the second halo must be set the loop over the     *
   !      * boundary subfaces is executed twice. This is the only correct  *
   !      * way in case the block contains only 1 cell between two         *
   !      * symmetry planes, i.e. a 2D problem.                            *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   LOGICAL, INTENT(IN) :: secondhalo
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: kk, mm, nn, i, j, l
   REAL(kind=realtype) :: vn, nnx, nny, nnz
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: gamma1, gamma2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww2d
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d
   INTERFACE 
   SUBROUTINE SETBCPOINTERS_EXTRA_D(nn, ww1, ww1d, ww2, ww2d, pp1, &
   &        pp1d, pp2, pp2d, rlv1, rlv2, rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS_EXTRA_D
   END INTERFACE
      INTERFACE 
   SUBROUTINE SETBCPOINTERS(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &        rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS
   END INTERFACE
      !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Set the value of kk; kk == 0 means only single halo, kk == 1
   ! double halo.
   kk = 0
   IF (secondhalo) kk = 1
   ! Loop over the number of times the halo computation must be done.
   nhalo:DO mm=0,kk
   ! Loop over the boundary condition subfaces of this block.
   bocos:DO nn=1,nbocos
   ! Check for symmetry boundary condition.
   IF (bctype(nn) .EQ. symm) THEN
   ! Nullify the pointers, because some compilers require that.
   !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
   ! Set the pointers to the correct subface.
   CALL SETBCPOINTERS_EXTRA_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
   &                             pp2, pp2d, rlv1, rlv2, rev1, rev2, mm)
   ! Set the additional pointers for gamma1 and gamma2.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   gamma1 => gamma(1, 1:, 1:)
   gamma2 => gamma(2, 1:, 1:)
   CASE (imax) 
   gamma1 => gamma(ie, 1:, 1:)
   gamma2 => gamma(il, 1:, 1:)
   CASE (jmin) 
   gamma1 => gamma(1:, 1, 1:)
   gamma2 => gamma(1:, 2, 1:)
   CASE (jmax) 
   gamma1 => gamma(1:, je, 1:)
   gamma2 => gamma(1:, jl, 1:)
   CASE (kmin) 
   gamma1 => gamma(1:, 1:, 1)
   gamma2 => gamma(1:, 1:, 2)
   CASE (kmax) 
   gamma1 => gamma(1:, 1:, ke)
   gamma2 => gamma(1:, 1:, kl)
   END SELECT
   ! Loop over the generic subface to set the state in the
   ! halo cells.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Store the three components of the unit normal a
   ! bit easier.
   nnx = bcdata(nn)%norm(i, j, 1)
   nny = bcdata(nn)%norm(i, j, 2)
   nnz = bcdata(nn)%norm(i, j, 3)
   ! Determine twice the normal velocity component,
   ! which must be substracted from the donor velocity
   ! to obtain the halo velocity.
   vn = two*(ww2(i, j, ivx)*nnx+ww2(i, j, ivy)*nny+ww2(i, j, &
   &              ivz)*nnz)
   ! Determine the flow variables in the halo cell.
   ww1(i, j, irho) = ww2(i, j, irho)
   ww1(i, j, ivx) = ww2(i, j, ivx) - vn*nnx
   ww1(i, j, ivy) = ww2(i, j, ivy) - vn*nny
   ww1(i, j, ivz) = ww2(i, j, ivz) - vn*nnz
   ww1(i, j, irhoe) = ww2(i, j, irhoe)
   ! Simply copy the turbulent variables.
   DO l=nt1mg,nt2mg
   ww1(i, j, l) = ww2(i, j, l)
   END DO
   ! Set the pressure and gamma and possibly the
   ! laminar and eddy viscosity in the halo.
   gamma1(i, j) = gamma2(i, j)
   pp1d(i, j) = pp2d(i, j)
   pp1(i, j) = pp2(i, j)
   IF (viscous) rlv1(i, j) = rlv2(i, j)
   IF (eddymodel) rev1(i, j) = rev2(i, j)
   END DO
   END DO
   END IF
   END DO bocos
   END DO nhalo
   END SUBROUTINE BCSYMM_EXTRA_D

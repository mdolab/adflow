!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of bcsymmpolar in forward (tangent) mode:
!   variations   of useful results: *p *w
!   with respect to varying inputs: *p *w
!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymmPolar.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-02-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE BCSYMMPOLAR_D(secondhalo)
  USE FLOWVARREFSTATE
  USE BLOCKPOINTERS
  USE BCTYPES
  USE CONSTANTS
  USE ITERATION
  IMPLICIT NONE
!
!      ******************************************************************
!      *                                                                *
!      * bcSymmPolar applies the polar symmetry boundary conditions     *
!      * to a singular line of a block. It is assumed that the pointers *
!      * in blockPointers are already set to the correct block on the   *
!      * correct grid level.                                            *
!      * The polar symmetry condition is a special case of a degenerate *
!      * line, as this line is the axi-symmetric centerline.            *
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
  INTEGER(kind=inttype) :: i, j, l, kk, mm, nn
  REAL(kind=realtype) :: nnx, nny, nnz, tmp, vtx, vty, vtz
  REAL(kind=realtype) :: tmpd, vtxd, vtyd, vtzd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: xline
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
  REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
  REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
  REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
  REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
  REAL(kind=realtype) :: arg1
  REAL(kind=realtype) :: result1
  INTRINSIC SQRT
  INTERFACE 
      SUBROUTINE SETBCPOINTERS_CD1(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
&        rev1, rev2, offset)
        USE BLOCKPOINTERS_D
        INTEGER(kind=inttype), INTENT(IN) :: nn, offset
        REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
      END SUBROUTINE SETBCPOINTERS_CD1
  END INTERFACE

  INTERFACE 
      SUBROUTINE SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
&        pp2, pp2d, rlv1, rlv2, rev1, rev2, offset)
        USE BLOCKPOINTERS_D
        INTEGER(kind=inttype), INTENT(IN) :: nn, offset
        REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
        REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
        REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
      END SUBROUTINE SETBCPOINTERS_D
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
! Loop over the boundary condition subfaces of this block.
bocos:DO nn=1,nbocos
! Check for the polar symmetry boundary condition.
    IF (bctype(nn) .EQ. symmpolar) THEN
! Nullify the pointers, because some compilers require that.
!nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
! Set the pointers for for the coordinates of the centerline.
! Needed to determine the direction of the velocity.
! This depends on the block face on which this subface is
! located.
      SELECT CASE  (bcfaceid(nn)) 
      CASE (imin) 
        xline => x(1, :, :, :)
      CASE (imax) 
        xline => x(il, :, :, :)
      CASE (jmin) 
        xline => x(:, 1, :, :)
      CASE (jmax) 
        xline => x(:, jl, :, :)
      CASE (kmin) 
        xline => x(:, :, 1, :)
      CASE (kmax) 
        xline => x(:, :, kl, :)
      END SELECT
! Loop over the number of times the halo computation must
! be done.
nhalo:DO mm=0,kk
! Set the pointers to the correct subface.
        CALL SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
&                       pp2d, rlv1, rlv2, rev1, rev2, mm)
! Loop over the generic subface to set the state in the
! halo cells.
        DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
! Determine the unit vector along the degenerated face.
! However it is not known which is the singular
! direction and therefore determine the direction along
! the diagonal (i,j) -- (i-1,j-1), which is correct for
! both singular i and j-direction. Note that due to the
! usage of the pointer xline there is an offset of +1
! in the indices and therefore (i+1,j+1) - (i,j) must
! be used to determine this vector.
            nnx = xline(i+1, j+1, 1) - xline(i, j, 1)
            nny = xline(i+1, j+1, 2) - xline(i, j, 2)
            nnz = xline(i+1, j+1, 3) - xline(i, j, 3)
! Determine the unit vector in this direction.
            arg1 = nnx*nnx + nny*nny + nnz*nnz
            result1 = SQRT(arg1)
            tmp = one/result1
            nnx = nnx*tmp
            nny = nny*tmp
            nnz = nnz*tmp
! Determine twice the tangential velocity vector of the
! internal cell.
            tmpd = two*(nnx*ww2d(i, j, ivx)+nny*ww2d(i, j, ivy)+nnz*ww2d&
&              (i, j, ivz))
            tmp = two*(ww2(i, j, ivx)*nnx+ww2(i, j, ivy)*nny+ww2(i, j, &
&              ivz)*nnz)
            vtxd = nnx*tmpd
            vtx = tmp*nnx
            vtyd = nny*tmpd
            vty = tmp*nny
            vtzd = nnz*tmpd
            vtz = tmp*nnz
! Determine the flow variables in the halo cell. The
! velocity is constructed such that the average of the
! internal and the halo cell is along the centerline.
! Note that the magnitude of the velocity does not
! change and thus the energy is identical.
            ww1d(i, j, irho) = ww2d(i, j, irho)
            ww1(i, j, irho) = ww2(i, j, irho)
            ww1d(i, j, ivx) = vtxd - ww2d(i, j, ivx)
            ww1(i, j, ivx) = vtx - ww2(i, j, ivx)
            ww1d(i, j, ivy) = vtyd - ww2d(i, j, ivy)
            ww1(i, j, ivy) = vty - ww2(i, j, ivy)
            ww1d(i, j, ivz) = vtzd - ww2d(i, j, ivz)
            ww1(i, j, ivz) = vtz - ww2(i, j, ivz)
            ww1d(i, j, irhoe) = ww2d(i, j, irhoe)
            ww1(i, j, irhoe) = ww2(i, j, irhoe)
! Simply copy the turbulent variables.
            DO l=nt1mg,nt2mg
              ww1d(i, j, l) = ww2d(i, j, l)
              ww1(i, j, l) = ww2(i, j, l)
            END DO
! Set the pressure and possibly the laminar and
! eddy viscosity in the halo.
            pp1d(i, j) = pp2d(i, j)
            pp1(i, j) = pp2(i, j)
            IF (viscous) rlv1(i, j) = rlv2(i, j)
            IF (eddymodel) rev1(i, j) = rev2(i, j)
          END DO
        END DO
      END DO nhalo
    END IF
  END DO bocos
END SUBROUTINE BCSYMMPOLAR_D

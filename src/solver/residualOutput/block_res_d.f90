!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of block_res in forward (tangent) mode:
!   variations   of useful results: *dw *w
!   with respect to varying inputs: *w
!   RW status of diff variables: *dw:out *w:in-out
! This is a super-combined function that combines the original
! functionality of: 
! Pressure Computation
! timeStep
! applyAllBCs
! initRes
! residual 
! The real difference between this and the original modules is that it
! it only operates on a single block at a time and as such the
! block/sps loop is outside the calculation. This routine is suitable
! for forward mode AD with Tapenade
SUBROUTINE BLOCK_RES_D(nn, sps)
  USE FLOWVARREFSTATE
  USE BLOCKPOINTERS
  USE INPUTPHYSICS
  IMPLICIT NONE
! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
! nw
  REAL(kind=realtype) :: gm1
  INTEGER(kind=inttype) :: nn, sps
  INTEGER :: k
  INTEGER :: j
  INTEGER :: i
  REAL :: v2
  REAL :: v2d
  INTRINSIC MAX
! Compute the pressures
  gm1 = gammaconstant - one
  pd = 0.0
! Compute P 
  DO k=2,kl
    DO j=2,jl
      DO i=2,il
        v2d = 2*w(i, j, k, ivx)*wd(i, j, k, ivx) + 2*w(i, j, k, ivy)*wd(&
&          i, j, k, ivy) + 2*w(i, j, k, ivz)*wd(i, j, k, ivz)
        v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**&
&          2
        pd(i, j, k) = gm1*(wd(i, j, k, irhoe)-half*(wd(i, j, k, irho)*v2&
&          +w(i, j, k, irho)*v2d))
        p(i, j, k) = gm1*(w(i, j, k, irhoe)-half*w(i, j, k, irho)*v2)
        IF (p(i, j, k) .LT. 1.e-4_realType*pinfcorr) THEN
          pd(i, j, k) = 0.0
          p(i, j, k) = 1.e-4_realType*pinfcorr
        ELSE
          p(i, j, k) = p(i, j, k)
        END IF
      END DO
    END DO
  END DO
  CALL TIMESTEP_BLOCK_D(.false.)
  CALL APPLYALLBC_BLOCK_D(.true.)
  CALL INITRES_BLOCK_D(1, nwf, nn, sps)
  CALL RESIDUAL_BLOCK_D()
END SUBROUTINE BLOCK_RES_D

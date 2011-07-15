!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!
!      ******************************************************************
!      *                                                                *
!      * File:          residual.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher (blanking)   *
!      * Starting date: 03-15-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE RESIDUAL_BLOCK_CD()
  USE INPUTITERATION_SPATIAL_D
  USE ITERATION_SPATIAL_D
  USE CGNSGRID_SPATIAL_D
  USE INPUTDISCRETIZATION_SPATIAL_D
  USE INPUTTIMESPECTRAL_SPATIAL_D
  USE BLOCKPOINTERS_SPATIAL_D
  USE FLOWVARREFSTATE_SPATIAL_D
  IMPLICIT NONE
!write(14,40),i,j,k,dw(i,j,k,l),fw(i,j,k,l)
!
!      ******************************************************************
!      *                                                                *
!      * residual computes the residual of the mean flow equations on   *
!      * the current MG level.                                          *
!      *                                                                *
!      ******************************************************************
!
!
!      Local variables.
!
  INTEGER(kind=inttype) :: sps, nn, discr
  INTEGER(kind=inttype) :: i, j, k, l
  LOGICAL :: finegrid
  REAL :: result1
  INTRINSIC REAL
  IF (smoother .EQ. rungekutta) THEN
    rfil = cdisrk(rkstage+1)
  ELSE
    rfil = one
  END IF
! Initialize the local arrays to monitor the massflows to zero.
  massflowfamilyinv = zero
  massflowfamilydiss = zero
! Set the value of the discretization, depending on the grid level,
! and the logical fineGrid, which indicates whether or not this
! is the finest grid level of the current mg cycle.
  discr = spacediscrcoarse
  IF (currentlevel .EQ. 1) discr = spacediscr
  finegrid = .false.
  IF (currentlevel .EQ. groundlevel) finegrid = .true.
  IF (finegrid .EQV. .false.) THEN
    PRINT*, 'Fine Grid should not be false here'
    STOP
  ELSE
    CALL INVISCIDCENTRALFLUX_CD()
! Compute the artificial dissipation fluxes.
! This depends on the parameter discr.
    SELECT CASE  (discr) 
    CASE (dissscalar) 
! Standard scalar dissipation scheme.
      IF (finegrid) THEN
        CALL INVISCIDDISSFLUXSCALAR_CD()
      ELSE
        CALL INVISCIDDISSFLUXSCALARCOARSE_CD()
      END IF
    CASE (dissmatrix) 
!===========================================================
! Matrix dissipation scheme.
      IF (finegrid) THEN
        CALL INVISCIDDISSFLUXMATRIX_CD()
      ELSE
        CALL INVISCIDDISSFLUXMATRIXCOARSE_CD()
      END IF
    CASE (disscusp) 
!===========================================================
! Cusp dissipation scheme.
      IF (finegrid) THEN
        CALL INVISCIDDISSFLUXCUSP_CD()
      ELSE
        CALL INVISCIDDISSFLUXCUSPCOARSE_CD()
      END IF
    CASE (upwind) 
!===========================================================
! Dissipation via an upwind scheme.
      CALL INVISCIDUPWINDFLUX_CD(finegrid)
    END SELECT
! Compute the viscous flux in case of a viscous computation.
    IF (viscous) CALL VISCOUSFLUX_CD()
! add the dissipative and possibly viscous fluxes to the
! Euler fluxes. Loop over the owned cells and add fw to dw.
! Also multiply by iblank so that no updates occur in holes
! or on the overset boundary.
    DO l=1,nwf
      DO k=2,kl
        DO j=2,jl
          DO i=2,il
            result1 = REAL(iblank(i, j, k), realtype)
            dw(i, j, k, l) = (dw(i, j, k, l)+fw(i, j, k, l))*result1
          END DO
        END DO
      END DO
    END DO
  END IF
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! Add the source terms from the level 0 cooling model.
!call level0CoolingModel
! Set the value of rFil, which controls the fraction of the old
! dissipation residual to be used. This is only for the runge-kutta
! schemes; for other smoothers rFil is simply set to 1.0.
! Note the index rkStage+1 for cdisRK. The reason is that the
! residual computation is performed before rkStage is incremented.
 40 FORMAT(1x,i4,i4,i4,e20.6,e20.6)
END SUBROUTINE RESIDUAL_BLOCK_CD

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of initres_block_ts in forward (tangent) mode:
!   variations   of useful results: *dw
!   with respect to varying inputs: *dw *w_offtimeinstance
!
!      ******************************************************************
!      *                                                                *
!      * File:          initres.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-18-2003                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE INITRES_BLOCK_TS_SPATIAL_D(varstart, varend, nn, sps, mm)
  USE INPUTITERATION_SPATIAL_D
  USE INPUTPHYSICS_SPATIAL_D
  USE ITERATION_SPATIAL_D
  USE INPUTTIMESPECTRAL_SPATIAL_D
  USE BLOCKPOINTERS_SPATIAL_D
  USE INPUTUNSTEADY_SPATIAL_D
  USE FLOWVARREFSTATE_SPATIAL_D
  IMPLICIT NONE
!
!      ******************************************************************
!      *                                                                *
!      * initres initializes the given range of the residual. Either to *
!      * zero, steady computation, or to an unsteady term for the time  *
!      * spectral and unsteady modes. For the coarser grid levels the   *
!      * residual forcing term is taken into account.                   *
!      *                                                                *
!      ******************************************************************
!
!
!      Subroutine arguments.
!
  INTEGER(kind=inttype), INTENT(IN) :: varstart, varend
!
!      Local variables.
!
  INTEGER(kind=inttype) :: sps, nn, mm, ll, ii, jj, i, j, k, l, m
  REAL(kind=realtype) :: oneoverdt, tmp
  REAL(kind=realtype) :: tmpd
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: ww, wsp1, wsp
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: wspd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: volsp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! Return immediately of no variables are in the range.
  IF (varend .LT. varstart) THEN
    RETURN
  ELSE
! Determine the equation mode and act accordingly.
    SELECT CASE  (equationmode) 
    CASE (timespectral) 
! Time spectral computation. The time derivative of the
! current solution is given by a linear combination of
! all other solutions, i.e. a matrix vector product.
! First store the section to which this block belongs
! in jj.
      jj = sectionid
! Determine the currently active multigrid level.
! Finest multigrid level. The residual must be
! initialized to the time derivative.
! Initialize it to zero.
! Loop over the number of terms which contribute
! to the time derivative.
! Store the pointer for the variable to be used to
! compute the unsteady source term and the volume.
! Also store in ii the offset needed for vector
! quantities.
!wsp   => flowDoms(nn,currentLevel,mm)%w
!volsp => flowDoms(nn,currentLevel,mm)%vol
      wspd => w_offtimeinstanced
      wsp => w_offtimeinstance
      volsp => vol_offtimeinstance
      ii = 3*(mm-1)
! Loop over the number of variables to be set.
varloopfine:DO l=varstart,varend
! Test for a momentum variable.
        IF ((l .EQ. ivx .OR. l .EQ. ivy) .OR. l .EQ. ivz) THEN
! Momentum variable. A special treatment is
! needed because it is a vector and the velocities
! are stored instead of the momentum. Set the
! coefficient ll, which defines the row of the
! matrix used later on.
          IF (l .EQ. ivx) ll = 3*sps - 2
          IF (l .EQ. ivy) ll = 3*sps - 1
          IF (l .EQ. ivz) ll = 3*sps
! Loop over the owned cell centers to add the
! contribution from wsp.
          DO k=2,kl
            DO j=2,jl
              DO i=2,il
! Store the matrix vector product with the
! velocity in tmp.
                tmpd = dvector(jj, ll, ii+1)*wspd(i, j, k, ivx) + &
&                  dvector(jj, ll, ii+2)*wspd(i, j, k, ivy) + dvector(jj&
&                  , ll, ii+3)*wspd(i, j, k, ivz)
                tmp = dvector(jj, ll, ii+1)*wsp(i, j, k, ivx) + dvector(&
&                  jj, ll, ii+2)*wsp(i, j, k, ivy) + dvector(jj, ll, ii+3&
&                  )*wsp(i, j, k, ivz)
! Update the residual. Note the
! multiplication with the density to obtain
! the correct time derivative for the
! momentum variable.
                dwd(i, j, k, l) = dwd(i, j, k, l) + volsp(i, j, k)*(tmpd&
&                  *wsp(i, j, k, irho)+tmp*wspd(i, j, k, irho))
                dw(i, j, k, l) = dw(i, j, k, l) + tmp*volsp(i, j, k)*wsp&
&                  (i, j, k, irho)
              END DO
            END DO
          END DO
        ELSE
! Scalar variable.  Loop over the owned cells to
! add the contribution of wsp to the time
! derivative.
          DO k=2,kl
            DO j=2,jl
              DO i=2,il
                dwd(i, j, k, l) = dwd(i, j, k, l) + dscalar(jj, sps, mm)&
&                  *volsp(i, j, k)*wspd(i, j, k, l)
                dw(i, j, k, l) = dw(i, j, k, l) + dscalar(jj, sps, mm)*&
&                  volsp(i, j, k)*wsp(i, j, k, l)
              END DO
            END DO
          END DO
        END IF
      END DO varloopfine
    END SELECT
  END IF
END SUBROUTINE INITRES_BLOCK_TS_SPATIAL_D

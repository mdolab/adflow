   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of residual_block in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *dw *w
   !   with respect to varying inputs: *rev *p *gamma *dw *w *rlv
   !                *radi *radj *radk
   !   Plus diff mem management of: rev:in p:in gamma:in dw:in w:in
   !                rlv:in fw:in radi:in radj:in radk:in
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
   SUBROUTINE RESIDUAL_BLOCK_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * residual computes the residual of the mean flow equations on   *
   !      * the current MG level.                                          *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_D
   USE CGNSGRID
   USE FLOWVARREFSTATE
   USE INPUTITERATION
   USE INPUTDISCRETIZATION
   USE INPUTTIMESPECTRAL
   USE ITERATION
   IMPLICIT NONE
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: sps, nn, discr
   INTEGER(kind=inttype) :: i, j, k, l
   LOGICAL :: finegrid
   INTRINSIC REAL
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Add the source terms from the level 0 cooling model.
   ! Set the value of rFil, which controls the fraction of the old
   ! dissipation residual to be used. This is only for the runge-kutta
   ! schemes; for other smoothers rFil is simply set to 1.0.
   ! Note the index rkStage+1 for cdisRK. The reason is that the
   ! residual computation is performed before rkStage is incremented.
   IF (smoother .EQ. rungekutta) THEN
   rfil = cdisrk(rkstage+1)
   ELSE
   rfil = one
   END IF
   ! Initialize the local arrays to monitor the massflows to zero.
   ! Set the value of the discretization, depending on the grid level,
   ! and the logical fineGrid, which indicates whether or not this
   ! is the finest grid level of the current mg cycle.
   discr = spacediscrcoarse
   IF (currentlevel .EQ. 1) discr = spacediscr
   finegrid = .false.
   IF (currentlevel .EQ. groundlevel) finegrid = .true.
   CALL INVISCIDCENTRALFLUX_D()
   ! Compute the artificial dissipation fluxes.
   ! This depends on the parameter discr.
   SELECT CASE  (discr) 
   CASE (dissscalar) 
   ! Standard scalar dissipation scheme.
   IF (finegrid) THEN
   CALL INVISCIDDISSFLUXSCALAR_D()
   ELSE
   CALL INVISCIDDISSFLUXSCALARCOARSE_D()
   END IF
   CASE (dissmatrix) 
   ! Reverse adjoint currently only work with invisciddissscalar
   !===========================================================
   ! Matrix dissipation scheme.
   IF (finegrid) THEN
   CALL INVISCIDDISSFLUXMATRIX_D()
   ELSE
   CALL INVISCIDDISSFLUXMATRIXCOARSE_D()
   END IF
   CASE (disscusp) 
   !===========================================================
   ! Cusp dissipation scheme.
   IF (finegrid) THEN
   CALL INVISCIDDISSFLUXCUSP()
   fwd = 0.0_8
   ELSE
   CALL INVISCIDDISSFLUXCUSPCOARSE()
   fwd = 0.0_8
   END IF
   CASE (upwind) 
   !===========================================================
   ! Dissipation via an upwind scheme.
   CALL INVISCIDUPWINDFLUX_D(finegrid)
   CASE DEFAULT
   fwd = 0.0_8
   END SELECT
   ! Compute the viscous flux in case of a viscous computation.
   IF (viscous) CALL VISCOUSFLUX_D()
   ! Add the dissipative and possibly viscous fluxes to the
   ! Euler fluxes. Loop over the owned cells and add fw to dw.
   ! Also multiply by iblank so that no updates occur in holes
   ! or on the overset boundary.
   DO l=1,nwf
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   dwd(i, j, k, l) = REAL(iblank(i, j, k), realtype)*(dwd(i, j, k&
   &           , l)+fwd(i, j, k, l))
   dw(i, j, k, l) = (dw(i, j, k, l)+fw(i, j, k, l))*REAL(iblank(i&
   &           , j, k), realtype)
   END DO
   END DO
   END DO
   END DO
   END SUBROUTINE RESIDUAL_BLOCK_D

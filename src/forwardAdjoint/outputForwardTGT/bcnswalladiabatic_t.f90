   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of bcnswalladiabatic in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: *rev *bvtj1 *bvtj2 *p *gamma
   !                *bmtk1 *w *bmtk2 *rlv *bvtk1 *bvtk2 *bmti1 *bmti2
   !                *bvti1 *bvti2 *bmtj1 *bmtj2
   !   with respect to varying inputs: *rev *p *w *rlv rgas
   !   Plus diff mem management of: rev:in bvtj1:in bvtj2:in p:in
   !                gamma:in bmtk1:in w:in bmtk2:in rlv:in bvtk1:in
   !                bvtk2:in bmti1:in bmti2:in bvti1:in bvti2:in bmtj1:in
   !                bmtj2:in bcdata:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcNsWallAdiabatic.f90                           *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-10-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCNSWALLADIABATIC_T(secondhalo, correctfork)
   USE CONSTANTS
   USE ITERATION
   USE FLOWVARREFSTATE
   USE BLOCKPOINTERS_D
   USE BCTYPES
   USE DIFFSIZES
   !  Hint: ISIZE3OFDrfrlv should be the size of dimension 3 of array *rlv
   !  Hint: ISIZE2OFDrfrlv should be the size of dimension 2 of array *rlv
   !  Hint: ISIZE1OFDrfrlv should be the size of dimension 1 of array *rlv
   !  Hint: ISIZE4OFDrfw should be the size of dimension 4 of array *w
   !  Hint: ISIZE3OFDrfw should be the size of dimension 3 of array *w
   !  Hint: ISIZE2OFDrfw should be the size of dimension 2 of array *w
   !  Hint: ISIZE1OFDrfw should be the size of dimension 1 of array *w
   !  Hint: ISIZE3OFDrfp should be the size of dimension 3 of array *p
   !  Hint: ISIZE2OFDrfp should be the size of dimension 2 of array *p
   !  Hint: ISIZE1OFDrfp should be the size of dimension 1 of array *p
   !  Hint: ISIZE3OFDrfrev should be the size of dimension 3 of array *rev
   !  Hint: ISIZE2OFDrfrev should be the size of dimension 2 of array *rev
   !  Hint: ISIZE1OFDrfrev should be the size of dimension 1 of array *rev
   !  Hint: ISIZE4OFDrfbmtj2 should be the size of dimension 4 of array *bmtj2
   !  Hint: ISIZE3OFDrfbmtj2 should be the size of dimension 3 of array *bmtj2
   !  Hint: ISIZE2OFDrfbmtj2 should be the size of dimension 2 of array *bmtj2
   !  Hint: ISIZE1OFDrfbmtj2 should be the size of dimension 1 of array *bmtj2
   !  Hint: ISIZE4OFDrfbmtj1 should be the size of dimension 4 of array *bmtj1
   !  Hint: ISIZE3OFDrfbmtj1 should be the size of dimension 3 of array *bmtj1
   !  Hint: ISIZE2OFDrfbmtj1 should be the size of dimension 2 of array *bmtj1
   !  Hint: ISIZE1OFDrfbmtj1 should be the size of dimension 1 of array *bmtj1
   !  Hint: ISIZE3OFDrfbvti2 should be the size of dimension 3 of array *bvti2
   !  Hint: ISIZE2OFDrfbvti2 should be the size of dimension 2 of array *bvti2
   !  Hint: ISIZE1OFDrfbvti2 should be the size of dimension 1 of array *bvti2
   !  Hint: ISIZE3OFDrfbvti1 should be the size of dimension 3 of array *bvti1
   !  Hint: ISIZE2OFDrfbvti1 should be the size of dimension 2 of array *bvti1
   !  Hint: ISIZE1OFDrfbvti1 should be the size of dimension 1 of array *bvti1
   !  Hint: ISIZE4OFDrfbmti2 should be the size of dimension 4 of array *bmti2
   !  Hint: ISIZE3OFDrfbmti2 should be the size of dimension 3 of array *bmti2
   !  Hint: ISIZE2OFDrfbmti2 should be the size of dimension 2 of array *bmti2
   !  Hint: ISIZE1OFDrfbmti2 should be the size of dimension 1 of array *bmti2
   !  Hint: ISIZE4OFDrfbmti1 should be the size of dimension 4 of array *bmti1
   !  Hint: ISIZE3OFDrfbmti1 should be the size of dimension 3 of array *bmti1
   !  Hint: ISIZE2OFDrfbmti1 should be the size of dimension 2 of array *bmti1
   !  Hint: ISIZE1OFDrfbmti1 should be the size of dimension 1 of array *bmti1
   !  Hint: ISIZE3OFDrfbvtk2 should be the size of dimension 3 of array *bvtk2
   !  Hint: ISIZE2OFDrfbvtk2 should be the size of dimension 2 of array *bvtk2
   !  Hint: ISIZE1OFDrfbvtk2 should be the size of dimension 1 of array *bvtk2
   !  Hint: ISIZE3OFDrfbvtk1 should be the size of dimension 3 of array *bvtk1
   !  Hint: ISIZE2OFDrfbvtk1 should be the size of dimension 2 of array *bvtk1
   !  Hint: ISIZE1OFDrfbvtk1 should be the size of dimension 1 of array *bvtk1
   !  Hint: ISIZE4OFDrfbmtk2 should be the size of dimension 4 of array *bmtk2
   !  Hint: ISIZE3OFDrfbmtk2 should be the size of dimension 3 of array *bmtk2
   !  Hint: ISIZE2OFDrfbmtk2 should be the size of dimension 2 of array *bmtk2
   !  Hint: ISIZE1OFDrfbmtk2 should be the size of dimension 1 of array *bmtk2
   !  Hint: ISIZE4OFDrfbmtk1 should be the size of dimension 4 of array *bmtk1
   !  Hint: ISIZE3OFDrfbmtk1 should be the size of dimension 3 of array *bmtk1
   !  Hint: ISIZE2OFDrfbmtk1 should be the size of dimension 2 of array *bmtk1
   !  Hint: ISIZE1OFDrfbmtk1 should be the size of dimension 1 of array *bmtk1
   !  Hint: ISIZE3OFDrfgamma should be the size of dimension 3 of array *gamma
   !  Hint: ISIZE2OFDrfgamma should be the size of dimension 2 of array *gamma
   !  Hint: ISIZE1OFDrfgamma should be the size of dimension 1 of array *gamma
   !  Hint: ISIZE3OFDrfbvtj2 should be the size of dimension 3 of array *bvtj2
   !  Hint: ISIZE2OFDrfbvtj2 should be the size of dimension 2 of array *bvtj2
   !  Hint: ISIZE1OFDrfbvtj2 should be the size of dimension 1 of array *bvtj2
   !  Hint: ISIZE3OFDrfbvtj1 should be the size of dimension 3 of array *bvtj1
   !  Hint: ISIZE2OFDrfbvtj1 should be the size of dimension 2 of array *bvtj1
   !  Hint: ISIZE1OFDrfbvtj1 should be the size of dimension 1 of array *bvtj1
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcNSWallAdiabatic applies the viscous adiabatic wall           *
   !      * boundary condition to a block. It is assumed that the pointers *
   !      * in blockPointers are already set to the correct block on the   *
   !      * correct grid level.                                            *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   LOGICAL, INTENT(IN) :: secondhalo, correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, i, j
   REAL(kind=realtype) :: rhok
   REAL(kind=realtype) :: rhokd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: uslip
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1d, rev2d
   INTERFACE 
   SUBROUTINE SETBCPOINTERS_T(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
   &        pp2, pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev1d, rev2, rev2d, &
   &        offset)
   USE BLOCKPOINTERS_D
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1d, rlv2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1d, rev2d
   END SUBROUTINE SETBCPOINTERS_T
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
      EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('rev', rev, revd, ISIZE1OFDrfrev*&
   &                        ISIZE2OFDrfrev*ISIZE3OFDrfrev)
   CALL DEBUG_TGT_REAL8ARRAY('p', p, pd, ISIZE1OFDrfp*ISIZE2OFDrfp*&
   &                        ISIZE3OFDrfp)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8ARRAY('rlv', rlv, rlvd, ISIZE1OFDrfrlv*&
   &                        ISIZE2OFDrfrlv*ISIZE3OFDrfrlv)
   CALL DEBUG_TGT_REAL8('rgas', rgas, rgasd)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! In case the turbulent transport equations are solved
   ! together with the mean flow equations, aplly the viscous
   ! wall boundary conditions for the turbulent variables.
   ! No need to extrapolate the secondary halo's, because this
   ! is done in extrapolate2ndHalo.
   IF (turbcoupled) THEN
   bmtj2d = 0.0_8
   bmtj1d = 0.0_8
   bvti2d = 0.0_8
   bvti1d = 0.0_8
   bmti2d = 0.0_8
   bmti1d = 0.0_8
   bvtk2d = 0.0_8
   bvtk1d = 0.0_8
   bmtk2d = 0.0_8
   bmtk1d = 0.0_8
   bvtj2d = 0.0_8
   bvtj1d = 0.0_8
   CALL DEBUG_TGT_CALL('TURBBCNSWALL', .TRUE., .FALSE.)
   CALL TURBBCNSWALL_T(.false.)
   CALL DEBUG_TGT_EXIT()
   gammad = 0.0_8
   ELSE
   bvtj1d = 0.0_8
   bvtj2d = 0.0_8
   gammad = 0.0_8
   bmtk1d = 0.0_8
   bmtk2d = 0.0_8
   bvtk1d = 0.0_8
   bvtk2d = 0.0_8
   bmti1d = 0.0_8
   bmti2d = 0.0_8
   bvti1d = 0.0_8
   bvti2d = 0.0_8
   bmtj1d = 0.0_8
   bmtj2d = 0.0_8
   END IF
   ! Loop over the viscous subfaces of this block. Note that
   ! these are numbered first.
   bocos:DO nn=1,nviscbocos
   ! Check for adiabatic viscous wall boundary conditions.
   IF (bctype(nn) .EQ. nswalladiabatic) THEN
   ! Set the pointer for uSlip to make the code more readable.
   uslip => bcdata(nn)%uslip
   CALL DEBUG_TGT_CALL('SETBCPOINTERS', .TRUE., .FALSE.)
   ! Nullify the pointers and set them to the correct subface.
   ! They are nullified first, because some compilers require
   ! that.
   !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
   CALL SETBCPOINTERS_T(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                     pp2d, rlv1, rlv1d, rlv2, rlv2d, rev1, rev1d, rev2, &
   &                     rev2d, 0)
   CALL DEBUG_TGT_EXIT()
   ! Initialize rhok to zero. This will be overwritten if a
   ! correction for k must be applied.
   rhok = zero
   rhokd = 0.0_8
   ! Loop over the generic subface to set the state in the
   ! halo cells.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('rev', rev, revd, ISIZE1OFDrfrev*&
   &                                ISIZE2OFDrfrev*ISIZE3OFDrfrev)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj1', bvtj1, bvtj1d, &
   &                                ISIZE1OFDrfbvtj1*ISIZE2OFDrfbvtj1*&
   &                                ISIZE3OFDrfbvtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj2', bvtj2, bvtj2d, &
   &                                ISIZE1OFDrfbvtj2*ISIZE2OFDrfbvtj2*&
   &                                ISIZE3OFDrfbvtj2)
   CALL DEBUG_TGT_REAL8ARRAY('p', p, pd, ISIZE1OFDrfp*&
   &                                ISIZE2OFDrfp*ISIZE3OFDrfp)
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, &
   &                                ISIZE1OFDrfgamma*ISIZE2OFDrfgamma*&
   &                                ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk1', bmtk1, bmtk1d, &
   &                                ISIZE1OFDrfbmtk1*ISIZE2OFDrfbmtk1*&
   &                                ISIZE3OFDrfbmtk1*ISIZE4OFDrfbmtk1)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*&
   &                                ISIZE2OFDrfw*ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk2', bmtk2, bmtk2d, &
   &                                ISIZE1OFDrfbmtk2*ISIZE2OFDrfbmtk2*&
   &                                ISIZE3OFDrfbmtk2*ISIZE4OFDrfbmtk2)
   CALL DEBUG_TGT_REAL8ARRAY('rlv', rlv, rlvd, ISIZE1OFDrfrlv*&
   &                                ISIZE2OFDrfrlv*ISIZE3OFDrfrlv)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk1', bvtk1, bvtk1d, &
   &                                ISIZE1OFDrfbvtk1*ISIZE2OFDrfbvtk1*&
   &                                ISIZE3OFDrfbvtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk2', bvtk2, bvtk2d, &
   &                                ISIZE1OFDrfbvtk2*ISIZE2OFDrfbvtk2*&
   &                                ISIZE3OFDrfbvtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bmti1', bmti1, bmti1d, &
   &                                ISIZE1OFDrfbmti1*ISIZE2OFDrfbmti1*&
   &                                ISIZE3OFDrfbmti1*ISIZE4OFDrfbmti1)
   CALL DEBUG_TGT_REAL8ARRAY('bmti2', bmti2, bmti2d, &
   &                                ISIZE1OFDrfbmti2*ISIZE2OFDrfbmti2*&
   &                                ISIZE3OFDrfbmti2*ISIZE4OFDrfbmti2)
   CALL DEBUG_TGT_REAL8ARRAY('bvti1', bvti1, bvti1d, &
   &                                ISIZE1OFDrfbvti1*ISIZE2OFDrfbvti1*&
   &                                ISIZE3OFDrfbvti1)
   CALL DEBUG_TGT_REAL8ARRAY('bvti2', bvti2, bvti2d, &
   &                                ISIZE1OFDrfbvti2*ISIZE2OFDrfbvti2*&
   &                                ISIZE3OFDrfbvti2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj1', bmtj1, bmtj1d, &
   &                                ISIZE1OFDrfbmtj1*ISIZE2OFDrfbmtj1*&
   &                                ISIZE3OFDrfbmtj1*ISIZE4OFDrfbmtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj2', bmtj2, bmtj2d, &
   &                                ISIZE1OFDrfbmtj2*ISIZE2OFDrfbmtj2*&
   &                                ISIZE3OFDrfbmtj2*ISIZE4OFDrfbmtj2)
   CALL DEBUG_TGT_REAL8('rgas', rgas, rgasd)
   CALL DEBUG_TGT_REAL8('rhok', rhok, rhokd)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   ! Set the value of rhok if a correcton must be applied.
   ! It probably does not matter too much, because k is very
   ! small near the wall.
   IF (correctfork) THEN
   rhokd = ww2d(i, j, irho)*ww2(i, j, itu1) + ww2(i, j, irho)*&
   &              ww2d(i, j, itu1)
   rhok = ww2(i, j, irho)*ww2(i, j, itu1)
   END IF
   ! Determine the variables in the halo. As the spacing
   ! is very small a constant pressure boundary condition
   ! (except for the k correction) is okay. Take the slip
   ! velocity into account.
   ww1d(i, j, irho) = ww2d(i, j, irho)
   ww1(i, j, irho) = ww2(i, j, irho)
   ww1d(i, j, ivx) = -ww2d(i, j, ivx)
   ww1(i, j, ivx) = -ww2(i, j, ivx) + two*uslip(i, j, 1)
   ww1d(i, j, ivy) = -ww2d(i, j, ivy)
   ww1(i, j, ivy) = -ww2(i, j, ivy) + two*uslip(i, j, 2)
   ww1d(i, j, ivz) = -ww2d(i, j, ivz)
   ww1(i, j, ivz) = -ww2(i, j, ivz) + two*uslip(i, j, 3)
   pp1d(i, j) = pp2d(i, j) - four*third*rhokd
   pp1(i, j) = pp2(i, j) - four*third*rhok
   ! Set the viscosities. There is no need to test for a
   ! viscous problem of course. The eddy viscosity is
   ! set to the negative value, as it should be zero on
   ! the wall.
   rlv1d(i, j) = rlv2d(i, j)
   rlv1(i, j) = rlv2(i, j)
   IF (eddymodel) THEN
   rev1d(i, j) = -rev2d(i, j)
   rev1(i, j) = -rev2(i, j)
   END IF
   END DO
   END DO
   CALL DEBUG_TGT_CALL('COMPUTEETOT', .TRUE., .FALSE.)
   ! Compute the energy for these halo's.
   CALL COMPUTEETOT_T(icbeg(nn), icend(nn), jcbeg(nn), jcend(nn), &
   &                   kcbeg(nn), kcend(nn), correctfork)
   CALL DEBUG_TGT_EXIT()
   ! Extrapolate the state vectors in case a second halo
   ! is needed.
   IF (secondhalo) THEN
   CALL DEBUG_TGT_CALL('EXTRAPOLATE2NDHALO', .TRUE., .FALSE.)
   CALL EXTRAPOLATE2NDHALO_T(nn, correctfork)
   CALL DEBUG_TGT_EXIT()
   END IF
   END IF
   END DO bocos
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('rev', rev, revd, ISIZE1OFDrfrev*&
   &                        ISIZE2OFDrfrev*ISIZE3OFDrfrev)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj1', bvtj1, bvtj1d, ISIZE1OFDrfbvtj1*&
   &                        ISIZE2OFDrfbvtj1*ISIZE3OFDrfbvtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj2', bvtj2, bvtj2d, ISIZE1OFDrfbvtj2*&
   &                        ISIZE2OFDrfbvtj2*ISIZE3OFDrfbvtj2)
   CALL DEBUG_TGT_REAL8ARRAY('p', p, pd, ISIZE1OFDrfp*ISIZE2OFDrfp*&
   &                        ISIZE3OFDrfp)
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, ISIZE1OFDrfgamma*&
   &                        ISIZE2OFDrfgamma*ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk1', bmtk1, bmtk1d, ISIZE1OFDrfbmtk1*&
   &                        ISIZE2OFDrfbmtk1*ISIZE3OFDrfbmtk1*&
   &                        ISIZE4OFDrfbmtk1)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk2', bmtk2, bmtk2d, ISIZE1OFDrfbmtk2*&
   &                        ISIZE2OFDrfbmtk2*ISIZE3OFDrfbmtk2*&
   &                        ISIZE4OFDrfbmtk2)
   CALL DEBUG_TGT_REAL8ARRAY('rlv', rlv, rlvd, ISIZE1OFDrfrlv*&
   &                        ISIZE2OFDrfrlv*ISIZE3OFDrfrlv)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk1', bvtk1, bvtk1d, ISIZE1OFDrfbvtk1*&
   &                        ISIZE2OFDrfbvtk1*ISIZE3OFDrfbvtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk2', bvtk2, bvtk2d, ISIZE1OFDrfbvtk2*&
   &                        ISIZE2OFDrfbvtk2*ISIZE3OFDrfbvtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bmti1', bmti1, bmti1d, ISIZE1OFDrfbmti1*&
   &                        ISIZE2OFDrfbmti1*ISIZE3OFDrfbmti1*&
   &                        ISIZE4OFDrfbmti1)
   CALL DEBUG_TGT_REAL8ARRAY('bmti2', bmti2, bmti2d, ISIZE1OFDrfbmti2*&
   &                        ISIZE2OFDrfbmti2*ISIZE3OFDrfbmti2*&
   &                        ISIZE4OFDrfbmti2)
   CALL DEBUG_TGT_REAL8ARRAY('bvti1', bvti1, bvti1d, ISIZE1OFDrfbvti1*&
   &                        ISIZE2OFDrfbvti1*ISIZE3OFDrfbvti1)
   CALL DEBUG_TGT_REAL8ARRAY('bvti2', bvti2, bvti2d, ISIZE1OFDrfbvti2*&
   &                        ISIZE2OFDrfbvti2*ISIZE3OFDrfbvti2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj1', bmtj1, bmtj1d, ISIZE1OFDrfbmtj1*&
   &                        ISIZE2OFDrfbmtj1*ISIZE3OFDrfbmtj1*&
   &                        ISIZE4OFDrfbmtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj2', bmtj2, bmtj2d, ISIZE1OFDrfbmtj2*&
   &                        ISIZE2OFDrfbmtj2*ISIZE3OFDrfbmtj2*&
   &                        ISIZE4OFDrfbmtj2)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   END SUBROUTINE BCNSWALLADIABATIC_T

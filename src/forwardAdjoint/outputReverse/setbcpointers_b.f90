   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of setbcpointers in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   Plus diff mem management of: rev:in p:in w:in rlv:in rev1:in-out
   !                rev2:in-out pp1:in-out pp2:in-out rlv1:in-out
   !                rlv2:in-out ww1:in-out ww2:in-out
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          setBcPointers.f90                               *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 02-17-2004                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE SETBCPOINTERS_B(nn, ww1, ww1b, ww2, ww2b, pp1, pp1b, pp2, &
   & pp2b, rlv1, rlv1b, rlv2, rlv2b, rev1, rev1b, rev2, rev2b, offset)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * setBCPointers sets the pointers needed for the boundary        *
   !      * condition treatment on a general face, such that the boundary  *
   !      * routines are only implemented once instead of 6 times.         *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BCTYPES
   USE BLOCKPOINTERS_B
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1b, ww2b
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1b, pp2b
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1b, rlv2b
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1b, rev2b
   !
   !      Local variables
   !
   INTEGER(kind=inttype) :: id, ih
   INTEGER :: branch
   INTERFACE 
   SUBROUTINE PUSHPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE PUSHPOINTER4
   SUBROUTINE LOOKPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE LOOKPOINTER4
   SUBROUTINE POPPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE POPPOINTER4
   END INTERFACE
      !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the face id on which the subface is located and set
   ! the pointers accordinly.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(ih, 1:, 1:, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(id, 1:, 1:, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(ih, 1:, 1:)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(id, 1:, 1:)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(ih, 1:, 1:)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(id, 1:, 1:)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(ih, 1:, 1:)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(id, 1:, 1:)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   CASE (imax) 
   !===============================================================
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(ih, 1:, 1:, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(id, 1:, 1:, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(ih, 1:, 1:)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(id, 1:, 1:)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(ih, 1:, 1:)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(id, 1:, 1:)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(ih, 1:, 1:)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(id, 1:, 1:)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   CASE (jmin) 
   !===============================================================
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(1:, ih, 1:, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(1:, id, 1:, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(1:, ih, 1:)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(1:, id, 1:)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(1:, ih, 1:)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(1:, id, 1:)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(1:, ih, 1:)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(1:, id, 1:)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   CASE (jmax) 
   !===============================================================
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(1:, ih, 1:, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(1:, id, 1:, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(1:, ih, 1:)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(1:, id, 1:)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(1:, ih, 1:)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(1:, id, 1:)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(1:, ih, 1:)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(1:, id, 1:)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   CASE (kmin) 
   !===============================================================
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(1:, 1:, ih, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(1:, 1:, id, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(1:, 1:, ih)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(1:, 1:, id)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(1:, 1:, ih)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(1:, 1:, id)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(1:, 1:, ih)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(1:, 1:, id)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   CASE (kmax) 
   !===============================================================
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(1:, 1:, ih, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(1:, 1:, id, :)
   CALL PUSHPOINTER4(pp1b)
   pp1b => pb(1:, 1:, ih)
   CALL PUSHPOINTER4(pp2b)
   pp2b => pb(1:, 1:, id)
   IF (viscous) THEN
   CALL PUSHPOINTER4(rlv1b)
   rlv1b => rlvb(1:, 1:, ih)
   CALL PUSHPOINTER4(rlv2b)
   rlv2b => rlvb(1:, 1:, id)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(1:, 1:, ih)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(1:, 1:, id)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1b)
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(rlv2b)
   CALL POPPOINTER4(rlv1b)
   END IF
   CALL POPPOINTER4(pp2b)
   CALL POPPOINTER4(pp1b)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1b)
   END SELECT
   END SUBROUTINE SETBCPOINTERS_B

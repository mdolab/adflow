   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of saeddyviscosity in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: *rev
   !   with respect to varying inputs: *rev *w *rlv
   !   Plus diff mem management of: rev:in w:in rlv:in
   !      ==================================================================
   !      ==================================================================
   SUBROUTINE SAEDDYVISCOSITY_T()
   USE CONSTANTS
   USE BLOCKPOINTERS_D
   USE PARAMTURB
   USE DIFFSIZES
   !  Hint: ISIZE3OFDrfrlv should be the size of dimension 3 of array *rlv
   !  Hint: ISIZE2OFDrfrlv should be the size of dimension 2 of array *rlv
   !  Hint: ISIZE1OFDrfrlv should be the size of dimension 1 of array *rlv
   !  Hint: ISIZE4OFDrfw should be the size of dimension 4 of array *w
   !  Hint: ISIZE3OFDrfw should be the size of dimension 3 of array *w
   !  Hint: ISIZE2OFDrfw should be the size of dimension 2 of array *w
   !  Hint: ISIZE1OFDrfw should be the size of dimension 1 of array *w
   !  Hint: ISIZE3OFDrfrev should be the size of dimension 3 of array *rev
   !  Hint: ISIZE2OFDrfrev should be the size of dimension 2 of array *rev
   !  Hint: ISIZE1OFDrfrev should be the size of dimension 1 of array *rev
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * saEddyViscosity computes the eddy-viscosity according to the   *
   !      * Spalart-Allmaras model for the block given in blockPointers.   *
   !      * This routine for both the original version as well as the      *
   !      * modified version according to Edwards.                         *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k
   REAL(kind=realtype) :: chi, chi3, fv1, rnusa, cv13
   REAL(kind=realtype) :: chid, chi3d, fv1d, rnusad
   EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8ARRAY('rlv', rlv, rlvd, ISIZE1OFDrfrlv*&
   &                        ISIZE2OFDrfrlv*ISIZE3OFDrfrlv)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Store the cv1^3; cv1 is a constant of the Spalart-Allmaras model.
   cv13 = rsacv1**3
   ! Loop over the cells of this block and compute the eddy viscosity.
   ! Do not include halo's.
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('rev', rev, revd, ISIZE1OFDrfrev*&
   &                              ISIZE2OFDrfrev*ISIZE3OFDrfrev)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*&
   &                              ISIZE2OFDrfw*ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8ARRAY('rlv', rlv, rlvd, ISIZE1OFDrfrlv*&
   &                              ISIZE2OFDrfrlv*ISIZE3OFDrfrlv)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   rnusad = wd(i, j, k, itu1)*w(i, j, k, irho) + w(i, j, k, itu1)*&
   &          wd(i, j, k, irho)
   rnusa = w(i, j, k, itu1)*w(i, j, k, irho)
   chid = (rnusad*rlv(i, j, k)-rnusa*rlvd(i, j, k))/rlv(i, j, k)**2
   chi = rnusa/rlv(i, j, k)
   chi3d = 3*chi**2*chid
   chi3 = chi**3
   fv1d = (chi3d*(chi3+cv13)-chi3*chi3d)/(chi3+cv13)**2
   fv1 = chi3/(chi3+cv13)
   revd(i, j, k) = fv1d*rnusa + fv1*rnusad
   rev(i, j, k) = fv1*rnusa
   END DO
   END DO
   END DO
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('rev', rev, revd, ISIZE1OFDrfrev*&
   &                        ISIZE2OFDrfrev*ISIZE3OFDrfrev)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   END SUBROUTINE SAEDDYVISCOSITY_T

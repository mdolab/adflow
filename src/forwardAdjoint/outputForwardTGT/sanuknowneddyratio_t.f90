   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of sanuknowneddyratio in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: sanuknowneddyratio
   !   with respect to varying inputs: nulam
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          saNuFromEddyRatio.f90                           *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 06-22-2003                                      *
   !      * Last modified: 04-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   FUNCTION SANUKNOWNEDDYRATIO_T(eddyratio, nulam, nulamd, &
   &  sanuknowneddyratio)
   USE CONSTANTS
   USE PARAMTURB
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * saNuKnownEddyRatio computes the Spalart-Allmaras transport     *
   !      * variable nu for the given eddy viscosity ratio.                *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Function type.
   !
   REAL(kind=realtype) :: sanuknowneddyratio
   REAL(kind=realtype) :: sanuknowneddyratio_t
   !
   !      Function arguments.
   !
   REAL(kind=realtype), INTENT(IN) :: eddyratio, nulam
   REAL(kind=realtype), INTENT(IN) :: nulamd
   !
   !      Local variables.
   !
   REAL(kind=realtype) :: cv13, chi, chi2, chi3, chi4, f, df, dchi
   INTRINSIC ABS
   EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   REAL(kind=realtype) :: abs0
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8('nulam', nulam, nulamd)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Take care of the exceptional cases.
   IF (eddyratio .LE. zero) THEN
   sanuknowneddyratio = zero
   sanuknowneddyratio_t = 0.0_8
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8('sanuknowneddyratio', sanuknowneddyratio, &
   &                     sanuknowneddyratio_t)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   RETURN
   ELSE
   ! Set the value of cv1^3, which is the constant appearing in the
   ! sa function fv1 to compute the eddy viscosity
   cv13 = rsacv1**3
   ! Determine the value of chi, which is given by the quartic
   ! polynomial chi^4 - ratio*(chi^3 + cv1^3) = 0.
   ! First determine the start value, depending on the eddyRatio.
   IF (eddyratio .LT. 1.e-4_realType) THEN
   chi = 0.5_realType
   ELSE IF (eddyratio .LT. 1.0_realType) THEN
   chi = 5.0_realType
   ELSE IF (eddyratio .LT. 10.0_realType) THEN
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8('nulam', nulam, nulamd)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   chi = 10.0_realType
   ELSE
   chi = eddyratio
   END IF
   ! The actual newton algorithm.
   DO 
   ! Compute the function value and the derivative.
   chi2 = chi*chi
   chi3 = chi*chi2
   chi4 = chi*chi3
   f = chi4 - eddyratio*(chi3+cv13)
   df = four*chi3 - three*eddyratio*chi2
   ! Compute the negative update and the new value of chi.
   dchi = f/df
   chi = chi - dchi
   IF (dchi/chi .GE. 0.) THEN
   abs0 = dchi/chi
   ELSE
   abs0 = -(dchi/chi)
   END IF
   ! Condition to exit the loop.
   IF (abs0 .LE. thresholdreal) GOTO 100
   END DO
   ! Chi is the ratio of the spalart allmaras transport variable and
   ! the laminar viscosity. So multiply chi with the laminar viscosity
   ! to obtain the correct value.
   100 sanuknowneddyratio_t = chi*nulamd
   sanuknowneddyratio = nulam*chi
   END IF
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8('sanuknowneddyratio', sanuknowneddyratio, &
   &                   sanuknowneddyratio_t)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   END FUNCTION SANUKNOWNEDDYRATIO_T

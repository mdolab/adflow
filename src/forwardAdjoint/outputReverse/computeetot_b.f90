   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of computeetot in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *p *w
   !   with respect to varying inputs: *p *w
   !   Plus diff mem management of: p:in w:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          computeEtot.F90                                 *
   !      * Author:        Edwin van der Weide, Steve Repsher              *
   !      * Starting date: 08-13-2003                                      *
   !      * Last modified: 10-14-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE COMPUTEETOT_B(istart, iend, jstart, jend, kstart, kend, &
   & correctfork)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * ComputeEtot computes the total energy from the given density,  *
   !      * velocity and presssure. For a calorically and thermally        *
   !      * perfect gas the well-known expression is used; for only a      *
   !      * thermally perfect gas, cp is a function of temperature, curve  *
   !      * fits are used and a more complex expression is obtained.       *
   !      * It is assumed that the pointers in blockPointers already       *
   !      * point to the correct block.                                    *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_B
   USE FLOWVARREFSTATE
   USE INPUTPHYSICS
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: istart, iend, jstart, jend
   INTEGER(kind=inttype), INTENT(IN) :: kstart, kend
   LOGICAL, INTENT(IN) :: correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k
   REAL(kind=realtype) :: ovgm1, factk, scale
   REAL(kind=realtype) :: tmp
   REAL(kind=realtype) :: tmp0
   REAL(kind=realtype) :: temp1
   REAL(kind=realtype) :: temp0
   REAL(kind=realtype) :: tmpb
   REAL(kind=realtype) :: tmpb0
   REAL(kind=realtype) :: tempb
   REAL(kind=realtype) :: temp
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the cp model used in the computation.
   SELECT CASE  (cpmodel) 
   CASE (cpconstant) 
   ! Constant cp and thus constant gamma.
   ! Abbreviate 1/(gamma -1) a bit easier.
   ovgm1 = one/(gammaconstant-one)
   ! Loop over the given range of the block and compute the first
   ! step of the energy.
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   tmp = ovgm1*p(i, j, k) + half*w(i, j, k, irho)*(w(i, j, k, ivx&
   &           )**2+w(i, j, k, ivy)**2+w(i, j, k, ivz)**2)
   CALL PUSHREAL8(w(i, j, k, irhoe))
   w(i, j, k, irhoe) = tmp
   END DO
   END DO
   END DO
   ! Second step. Correct the energy in case a turbulent kinetic
   ! energy is present.
   IF (correctfork) THEN
   factk = ovgm1*(five*third-gammaconstant)
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   tmp0 = w(i, j, k, irhoe) - factk*w(i, j, k, irho)*w(i, j, k&
   &             , itu1)
   CALL PUSHREAL8(w(i, j, k, irhoe))
   w(i, j, k, irhoe) = tmp0
   END DO
   END DO
   END DO
   DO k=kend,kstart,-1
   DO j=jend,jstart,-1
   DO i=iend,istart,-1
   CALL POPREAL8(w(i, j, k, irhoe))
   tmpb0 = wb(i, j, k, irhoe)
   wb(i, j, k, irhoe) = tmpb0
   wb(i, j, k, irho) = wb(i, j, k, irho) - factk*w(i, j, k, &
   &             itu1)*tmpb0
   wb(i, j, k, itu1) = wb(i, j, k, itu1) - factk*w(i, j, k, &
   &             irho)*tmpb0
   END DO
   END DO
   END DO
   END IF
   DO k=kend,kstart,-1
   DO j=jend,jstart,-1
   DO i=iend,istart,-1
   CALL POPREAL8(w(i, j, k, irhoe))
   tmpb = wb(i, j, k, irhoe)
   wb(i, j, k, irhoe) = 0.0_8
   temp1 = w(i, j, k, ivz)
   temp0 = w(i, j, k, ivy)
   temp = w(i, j, k, ivx)
   tempb = half*w(i, j, k, irho)*tmpb
   pb(i, j, k) = pb(i, j, k) + ovgm1*tmpb
   wb(i, j, k, irho) = wb(i, j, k, irho) + half*(temp**2+temp0**2&
   &           +temp1**2)*tmpb
   wb(i, j, k, ivx) = wb(i, j, k, ivx) + 2*temp*tempb
   wb(i, j, k, ivy) = wb(i, j, k, ivy) + 2*temp0*tempb
   wb(i, j, k, ivz) = wb(i, j, k, ivz) + 2*temp1*tempb
   END DO
   END DO
   END DO
   CASE (cptempcurvefits) 
   !        ================================================================
   ! Cp as function of the temperature is given via curve fits.
   ! Store a scale factor to compute the nonDimensional
   ! internal energy.
   scale = rgas/tref
   ! Loop over the given range of the block.
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   CALL PUSHREAL8ARRAY(w, SIZE(w, 1)*SIZE(w, 2)*SIZE(w, 3)*SIZE(w&
   &                       , 4))
   CALL COMPUTEETOTCELLCPFIT(i, j, k, scale, correctfork)
   END DO
   END DO
   END DO
   DO k=kend,kstart,-1
   DO j=jend,jstart,-1
   DO i=iend,istart,-1
   CALL POPREAL8ARRAY(w, SIZE(w, 1)*SIZE(w, 2)*SIZE(w, 3)*SIZE(w&
   &                      , 4))
   CALL COMPUTEETOTCELLCPFIT_B(i, j, k, scale, correctfork)
   END DO
   END DO
   END DO
   END SELECT
   !
   40 FORMAT(1x,i4,i4,i4,e20.6)
   END SUBROUTINE COMPUTEETOT_B

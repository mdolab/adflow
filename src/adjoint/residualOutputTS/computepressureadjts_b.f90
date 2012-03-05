   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
   !  
   !  Differentiation of computepressureadjts in reverse (adjoint) mode:
   !   gradient, with respect to input variables: wadj
   !   of linear combination of output variables: padj wadj
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          computePressureAdj.f90                          *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-19-2006                                      *
   !      * Last modified: 03-20-2006                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE COMPUTEPRESSUREADJTS_B(wadj, wadjb, padj, padjb, nn, level, &
   &  sps)
   USE flowvarrefstate
   USE inputphysics
   USE inputtimespectral
   IMPLICIT NONE
   INTEGER(KIND=INTTYPE) :: level, nn, sps
   REAL(KIND=REALTYPE) :: padj(-2:2, -2:2, -2:2, ntimeintervalsspectral)&
   &  , padjb(-2:2, -2:2, -2:2, ntimeintervalsspectral)
   REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral), INTENT(IN) :: wadj
   REAL(KIND=REALTYPE) :: wadjb(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral)
   INTEGER(KIND=INTTYPE) :: i, j, k
   REAL(KIND=REALTYPE) :: tempb, tempb0
   REAL(KIND=REALTYPE) :: factk, gm1, v2, v2b
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Simple routine to compute the pressure from the variables w.   *
   !      * A calorically perfect gas, i.e. constant gamma, is assumed.    *
   !      *                                                                *
   !      ******************************************************************
   !
   !nIntervalTimespectral
   !
   !      Subroutine arguments
   !
   !
   !      Local variables
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   gm1 = gammaconstant - one
   ! Check the situation.
   IF (kpresent) THEN
   ! A separate equation for the turbulent kinetic energy is
   ! present. This variable must be taken into account.
   factk = five*third - gammaconstant
   DO k=-2,2
   DO j=-2,2
   DO i=-2,2
   CALL PUSHREAL8(v2)
   v2 = wadj(i, j, k, ivx, sps)**2 + wadj(i, j, k, ivy, sps)**2 +&
   &           wadj(i, j, k, ivz, sps)**2
   END DO
   END DO
   END DO
   DO k=2,-2,-1
   DO j=2,-2,-1
   DO i=2,-2,-1
   tempb = gm1*padjb(i, j, k, sps)
   wadjb(i, j, k, irhoe, sps) = wadjb(i, j, k, irhoe, sps) + &
   &            tempb
   wadjb(i, j, k, irho, sps) = wadjb(i, j, k, irho, sps) + factk*&
   &            wadj(i, j, k, itu1, sps)*padjb(i, j, k, sps) - half*v2*&
   &            tempb
   v2b = -(half*wadj(i, j, k, irho, sps)*tempb)
   wadjb(i, j, k, itu1, sps) = wadjb(i, j, k, itu1, sps) + factk*&
   &            wadj(i, j, k, irho, sps)*padjb(i, j, k, sps)
   padjb(i, j, k, sps) = 0.0
   CALL POPREAL8(v2)
   wadjb(i, j, k, ivx, sps) = wadjb(i, j, k, ivx, sps) + 2*wadj(i&
   &            , j, k, ivx, sps)*v2b
   wadjb(i, j, k, ivy, sps) = wadjb(i, j, k, ivy, sps) + 2*wadj(i&
   &            , j, k, ivy, sps)*v2b
   wadjb(i, j, k, ivz, sps) = wadjb(i, j, k, ivz, sps) + 2*wadj(i&
   &            , j, k, ivz, sps)*v2b
   END DO
   END DO
   END DO
   ELSE
   ! No separate equation for the turbulent kinetic enery.
   ! Use the standard formula.
   DO k=-2,2
   DO j=-2,2
   DO i=-2,2
   CALL PUSHREAL8(v2)
   v2 = wadj(i, j, k, ivx, sps)**2 + wadj(i, j, k, ivy, sps)**2 +&
   &           wadj(i, j, k, ivz, sps)**2
   END DO
   END DO
   END DO
   DO k=2,-2,-1
   DO j=2,-2,-1
   DO i=2,-2,-1
   tempb0 = gm1*padjb(i, j, k, sps)
   wadjb(i, j, k, irhoe, sps) = wadjb(i, j, k, irhoe, sps) + &
   &            tempb0
   wadjb(i, j, k, irho, sps) = wadjb(i, j, k, irho, sps) - half*&
   &            v2*tempb0
   v2b = -(half*wadj(i, j, k, irho, sps)*tempb0)
   padjb(i, j, k, sps) = 0.0
   CALL POPREAL8(v2)
   wadjb(i, j, k, ivx, sps) = wadjb(i, j, k, ivx, sps) + 2*wadj(i&
   &            , j, k, ivx, sps)*v2b
   wadjb(i, j, k, ivy, sps) = wadjb(i, j, k, ivy, sps) + 2*wadj(i&
   &            , j, k, ivy, sps)*v2b
   wadjb(i, j, k, ivz, sps) = wadjb(i, j, k, ivz, sps) + 2*wadj(i&
   &            , j, k, ivz, sps)*v2b
   END DO
   END DO
   END DO
   END IF
   END SUBROUTINE COMPUTEPRESSUREADJTS_B

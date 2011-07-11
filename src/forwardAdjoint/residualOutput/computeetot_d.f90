!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of computeetot in forward (tangent) mode:
!   variations   of useful results: *gamma *w
!   with respect to varying inputs: *p *gamma *w
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
SUBROUTINE COMPUTEETOT_D(istart, iend, jstart, jend, kstart, kend, &
     &  correctfork)
  USE FLOWVARREFSTATE
  USE BLOCKPOINTERS_D
  USE INPUTPHYSICS
  IMPLICIT NONE
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
  !
  40 format (1x,I4,I4,I4,E20.6)
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
              wd(i, j, k, irhoe) = ovgm1*pd(i, j, k) + half*(wd(i, j, k, &
                   &            irho)*(w(i, j, k, ivx)**2+w(i, j, k, ivy)**2+w(i, j, k, ivz)&
                   &            **2)+w(i, j, k, irho)*(2*w(i, j, k, ivx)*wd(i, j, k, ivx)+2*&
                   &            w(i, j, k, ivy)*wd(i, j, k, ivy)+2*w(i, j, k, ivz)*wd(i, j, &
                   &            k, ivz)))
              w(i, j, k, irhoe) = ovgm1*p(i, j, k) + half*w(i, j, k, irho)*(&
                   &            w(i, j, k, ivx)**2+w(i, j, k, ivy)**2+w(i, j, k, ivz)**2)
              !write(13,40),i,j,k,w(i,j,k,irhoe)
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
                 wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - factk*(wd(i, j, k&
                      &              , irho)*w(i, j, k, itu1)+w(i, j, k, irho)*wd(i, j, k, itu1&
                      &              ))
                 w(i, j, k, irhoe) = w(i, j, k, irhoe) - factk*w(i, j, k, &
                      &              irho)*w(i, j, k, itu1)
                 !write(13,40),i,j,k,w(i,j,k,irhoe)
              END DO
           END DO
        END DO
     END IF
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
              CALL COMPUTEETOTCELLCPFIT_D(i, j, k, scale, correctfork)
           END DO
        END DO
     END DO
  END SELECT
END SUBROUTINE COMPUTEETOT_D

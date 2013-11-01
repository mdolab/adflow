   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of etotarray in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: etot
   !   with respect to varying inputs: k p u v w etot rho (global)tref
   !                (global)rgas
   !      ==================================================================
   SUBROUTINE ETOTARRAY_D(rho, rhod, u, ud, v, vd, w, wd, p, pd, k, kd, &
   &  etot, etotd, correctfork, kk)
   USE CONSTANTS
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * EtotArray computes the total energy from the given density,    *
   !      * velocity and presssure for the given kk elements of the arrays.*
   !      * First the internal energy per unit mass is computed and after  *
   !      * that the kinetic energy is added as well the conversion to     *
   !      * energy per unit volume.                                        *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   REAL(kind=realtype), DIMENSION(*), INTENT(IN) :: rho, p, k
   REAL(kind=realtype), DIMENSION(*), INTENT(IN) :: rhod, pd, kd
   REAL(kind=realtype), DIMENSION(*), INTENT(IN) :: u, v, w
   REAL(kind=realtype), DIMENSION(*), INTENT(IN) :: ud, vd, wd
   REAL(kind=realtype), DIMENSION(*), INTENT(OUT) :: etot
   REAL(kind=realtype), DIMENSION(*), INTENT(OUT) :: etotd
   LOGICAL, INTENT(IN) :: correctfork
   INTEGER(kind=inttype), INTENT(IN) :: kk
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Compute the internal energy for unit mass.
   CALL EINTARRAY_D(rho, rhod, p, pd, k, kd, etot, etotd, correctfork, kk&
   &            )
   ! Add the kinetic energy.
   DO i=1,kk
   etotd(i) = rhod(i)*(etot(i)+half*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i))) + &
   &      rho(i)*(etotd(i)+half*(ud(i)*u(i)+u(i)*ud(i)+vd(i)*v(i)+v(i)*vd(i)&
   &      +wd(i)*w(i)+w(i)*wd(i)))
   etot(i) = rho(i)*(etot(i)+half*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)))
   END DO
   END SUBROUTINE ETOTARRAY_D

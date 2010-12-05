   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of eintarraynkpc in forward (tangent) mode:
   !   variations   of useful results: eint
   !   with respect to varying inputs: gammaconstant k p rho
   !      ==================================================================
   SUBROUTINE EINTARRAYNKPC_D(rho, rhod, p, pd, k, kd, eint, eintd, &
   &  correctfork, kk)
   USE FLOWVARREFSTATE
   USE CPCURVEFITS
   USE INPUTPHYSICS
   USE CONSTANTS
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: kk
   !
   !      ******************************************************************
   !      *                                                                *
   !      * EintArray computes the internal energy per unit mass from the  *
   !      * given density and pressure (and possibly turbulent energy) for *
   !      * the given kk elements of the arrays.                           *
   !      * For a calorically and thermally perfect gas the well-known     *
   !      * expression is used; for only a thermally perfect gas, cp is a  *
   !      * function of temperature, curve fits are used and a more        *
   !      * complex expression is obtained.                                *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   REAL(kind=realtype), DIMENSION(kk), INTENT(IN) :: rho, p, k
   REAL(kind=realtype), DIMENSION(kk), INTENT(IN) :: rhod, pd, kd
   REAL(kind=realtype), DIMENSION(kk), INTENT(OUT) :: eint
   REAL(kind=realtype), DIMENSION(kk), INTENT(OUT) :: eintd
   !!$       real(kind=realType), dimension(*), intent(in)  :: rho, p, k
   !!$       real(kind=realType), dimension(*), intent(out) :: eint
   LOGICAL, INTENT(IN) :: correctfork
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: twothird=two*third
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, nn, mm, ii, start
   REAL(kind=realtype) :: ovgm1, factk, pp, t, t2, scale
   !!$           ! Loop over the number of elements of the array
   !!$
   !!$           do i=1,kk
   !!$
   !!$             ! Compute the dimensional temperature.
   !!$
   !!$             pp = p(i)
   !!$             if( correctForK ) pp = pp - twoThird*rho(i)*k(i)
   !!$             t = Tref*pp/(RGas*rho(i))
   !!$
   !!$             ! Determine the case we are having here.
   !!$
   !!$             if(t <= cpTrange(0)) then
   !!$
   !!$               ! Temperature is less than the smallest temperature
   !!$               ! in the curve fits. Use extrapolation using
   !!$               ! constant cv.
   !!$
   !!$               eint(i) = scale*(cpEint(0) + cv0*(t - cpTrange(0)))
   !!$
   !!$             else if(t >= cpTrange(cpNparts)) then
   !!$
   !!$               ! Temperature is larger than the largest temperature
   !!$               ! in the curve fits. Use extrapolation using
   !!$               ! constant cv.
   !!$
   !!$               eint(i) = scale*(cpEint(cpNparts) &
   !!$                       +        cvn*(t - cpTrange(cpNparts)))
   !!$
   !!$             else
   !!$
   !!$               ! Temperature is in the curve fit range.
   !!$               ! First find the valid range.
   !!$
   !!$               ii    = cpNparts
   !!$               start = 1
   !!$               interval: do
   !!$
   !!$                 ! Next guess for the interval.
   !!$
   !!$                 nn = start + ii/2
   !!$
   !!$                 ! Determine the situation we are having here.
   !!$
   !!$                 if(t > cpTrange(nn)) then
   !!$
   !!$                   ! Temperature is larger than the upper boundary of
   !!$                   ! the current interval. Update the lower boundary.
   !!$
   !!$                   start = nn + 1
   !!$                   ii    = ii - 1
   !!$
   !!$                 else if(t >= cpTrange(nn-1)) then
   !!$
   !!$                   ! This is the correct range. Exit the do-loop.
   !!$
   !!$                   exit
   !!$
   !!$                 endif
   !!$
   !!$                 ! Modify ii for the next branch to search.
   !!$
   !!$                 ii = ii/2
   !!$
   !!$               enddo interval
   !!$
   !!$               ! Nn contains the correct curve fit interval.
   !!$               ! Integrate cv to compute eint.
   !!$
   !!$               eint(i) = cpTempFit(nn)%eint0 - t
   !!$               do ii=1,cpTempFit(nn)%nterm
   !!$                 if(cpTempFit(nn)%exponents(ii) == -1_intType) then
   !!$                   eint(i) = eint(i) &
   !!$                           + cpTempFit(nn)%constants(ii)*log(t)
   !!$                 else
   !!$                   mm   = cpTempFit(nn)%exponents(ii) + 1
   !!$                   t2   = t**mm
   !!$                   eint(i) = eint(i) &
   !!$                           + cpTempFit(nn)%constants(ii)*t2/mm
   !!$                 endif
   !!$               enddo
   !!$
   !!$               eint(i) = scale*eint(i)
   !!$
   !!$             endif
   !!$
   !!$             ! Add the turbulent energy if needed.
   !!$
   !!$             if( correctForK ) eint(i) = eint(i) + k(i)
   !!$
   !!$           enddo
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the cp model used in the computation.
   SELECT CASE  (cpmodel) 
   CASE (cpconstant) 
   ! Abbreviate 1/(gamma -1) a bit easier.
   ovgm1 = one/(gammaconstant-one)
   eintd = 0.0
   ! Loop over the number of elements of the array and compute
   ! the total energy.
   DO i=1,kk
   eintd(i) = (ovgm1*pd(i)*rho(i)-ovgm1*p(i)*rhod(i))/rho(i)**2
   eint(i) = ovgm1*p(i)/rho(i)
   END DO
   ! Second step. Correct the energy in case a turbulent kinetic
   ! energy is present.
   IF (correctfork) THEN
   factk = ovgm1*(five*third-gammaconstant)
   DO i=1,kk
   eintd(i) = eintd(i) - factk*kd(i)
   eint(i) = eint(i) - factk*k(i)
   END DO
   END IF
   CASE (cptempcurvefits) 
   !        ================================================================
   ! Cp as function of the temperature is given via curve fits.
   ! Store a scale factor to compute the nonDimensional
   ! internal energy.
   scale = rgas/tref
   eintd = 0.0
   CASE DEFAULT
   eintd = 0.0
   END SELECT
   END SUBROUTINE EINTARRAYNKPC_D

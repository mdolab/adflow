!
!      ******************************************************************
!      *                                                                *
!      * File:          computePtot.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-16-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computePtot(rho, u, v, w, p, ptot, kk)
!
!      ******************************************************************
!      *                                                                *
!      * ComputePtot computes the total pressure for the given          *
!      * pressures, densities and velocities.                           *
!      *                                                                *
!      ******************************************************************
!
       use cpCurveFits
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: kk

       real(kind=realType), dimension(*), intent(in)  :: rho, p, u, v, w
       real(kind=realType), dimension(*), intent(out) :: ptot
!
!      Local parameters.
!
       real(kind=realType), parameter :: dtStop  = 0.01_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, ii, mm, nn, nnt, start

       real(kind=realType) :: govgm1, kin
       real(kind=realType) :: T, T2, tt, dt, h, htot, cp, scale, alp
       real(kind=realType) :: intCport, intCportT, intCportTt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the cp model used.

       select case (cpModel)

         case (cpConstant)

           ! Constant cp and thus constant gamma. The well-known
           ! formula is valid.

           govgm1 = gammaInf/(gammaInf-one)

           do i=1,kk
             kin     = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
             ptot(i) = p(i)*((one + rho(i)*kin/(govgm1*p(i)))**govgm1)
           enddo

         !===============================================================

         case (cpTempCurveFits)

           ! Cp is a function of the temperature. The formula used for
           ! constant cp is not valid anymore and a more complicated
           ! procedure must be followed.

           do i=1,kk

             ! Compute the dimensional temperature and the scale
             ! factor to convert the integral of cp to the correct
             ! nonDimensional value.

             T     = Tref*p(i)/(RGas*rho(i))
             scale = RGas/Tref

             ! Compute the enthalpy and the integrand of cp/(r*t) at the
             ! given temperature. Take care of the exceptional situations.

             if(T <= cpTrange(0)) then

               ! Temperature is smaller than the smallest temperature in
               ! the curve fit range. Use extrapolation using constant cp.

               nn = 0
               cp = cv0+one
               h  = scale*(cpHint(0) + cp*(T-cpTrange(0)))

               intCportT = cp*log(T)

             else if(T >= cpTrange(cpNparts)) then

               ! Temperature is larger than the largest temperature in the
               ! curve fit range. Use extrapolation using constant cp.

               nn = cpNparts + 1
               cp = cvn+one
               h  = scale*(cpHint(cpNparts) + cp*(T-cpTrange(cpNparts)))

               intCportT = cp*log(T)

             else

               ! Temperature lies in the curve fit range. Find the correct
               ! interval.

               ii    = cpNparts
               start = 1
               interval: do

                 ! Next guess for the interval.

                 nn = start + ii/2

                 ! Determine the situation we are having here.

                 if(T > cpTrange(nn)) then

                   ! Temperature is larger than the upper boundary of
                   ! the current interval. Update the lower boundary.

                   start = nn + 1
                   ii    = ii - 1

                 else if(T >= cpTrange(nn-1)) then

                   ! This is the correct range. Exit the do-loop.

                   exit

                 endif

                 ! Modify ii for the next branch to search.

                 ii = ii/2

               enddo interval

               ! nn contains the correct curve fit interval.
               ! Integrate cp to compute h and the integrand of cp/(r*t)

               h           = cpTempFit(nn)%eint0
               intCportT = zero

               do ii=1,cpTempFit(nn)%nterm
                 if(cpTempFit(nn)%exponents(ii) == -1_intType) then
                   h = h + cpTempFit(nn)%constants(ii)*log(T)
                 else
                   mm = cpTempFit(nn)%exponents(ii) + 1
                   T2 = T**mm
                   h  = h + cpTempFit(nn)%constants(ii)*T2/mm
                 endif

                 if(cpTempFit(nn)%exponents(ii) == 0_intType) then
                   intCportT = intCportT &
                             + cpTempFit(nn)%constants(ii)*log(T)
                 else
                   mm = cpTempFit(nn)%exponents(ii)
                   T2 = T**mm
                   intCportT = intCportT &
                             + cpTempFit(nn)%constants(ii)*T2/mm
                 endif
               enddo

               h = scale*h

             endif

             ! Compute the total enthalpy. Divide by scale to get the same
             ! dimensions as for the integral of cp/r.

             htot = (h + half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)))/scale

             ! Compute the corresponding total temperature. First determine
             ! the situation we are having here.

             if(htot <= cpHint(0)) then

               ! Total enthalpy is smaller than the lowest value of the
               ! curve fit. Use extrapolation using constant cp.

               nnt = 0
               Tt  = cpTrange(0) + (htot - cpHint(0))/(cv0+one)

             else if(htot >= cpHint(cpNparts)) then

               ! Total enthalpy is larger than the largest value of the
               ! curve fit. Use extrapolation using constant cp.

               nnt = cpNparts + 1
               Tt  = cpTrange(cpNparts) &
                   + (htot - cpHint(cpNparts))/(cvn+one)

             else

               ! Total temperature is in the range of the curve fits.
               ! Use a newton algorithm to find the correct temperature.
               ! First find the correct interval.

               ii    = cpNparts
               start = 1
               intervalTt: do

                 ! Next guess for the interval.

                 nnt = start + ii/2

                 ! Determine the situation we are having here.

                 if(htot > cpHint(nnt)) then

                   ! Enthalpy is larger than the upper boundary of
                   ! the current interval. Update the lower boundary.

                   start = nnt + 1
                   ii    = ii - 1

                 else if(htot >= cpHint(nnt-1)) then

                   ! This is the correct range. Exit the do-loop.

                   exit

                 endif

                 ! Modify ii for the next branch to search.

                 ii = ii/2

               enddo intervalTt

               ! Nnt contains the range in which the newton algorithm must
               ! be applied. Initial guess of the total temperature.

               alp = (cpHint(nnt) - htot)/(cpHint(nnt) - cpHint(nnt-1))
               Tt  = alp*cpTrange(nnt-1) + (one-alp)*cpTrange(nnt)

               ! The actual newton algorithm to compute the total
               ! temperature.

               newton: do

                 ! Compute the energy as well as the value of cv/r for the
                 ! given temperature.

                 cp = zero
                 h  = cpTempFit(nnt)%eint0

                   do ii=1,cpTempFit(nnt)%nterm

                   ! Update cp.

                   T2 = Tt**(cpTempFit(nnt)%exponents(ii))
                   cp = cp + cpTempFit(nnt)%constants(ii)*t2

                   ! Update h, for which this contribution must be
                   ! integrated. Take the exceptional case that the
                   ! exponent == -1 into account.

                   if(cpTempFit(nnt)%exponents(ii) == -1_intType) then
                     h = h + cpTempFit(nnt)%constants(ii)*log(Tt)
                   else
                     h = h + cpTempFit(nnt)%constants(ii)*t2*Tt &
                       / (cpTempFit(nnt)%exponents(ii) + 1)
                   endif

                 enddo

                 ! Compute the update and the new total temperature.

                 dT = (htot - h)/cp
                 Tt = Tt + dT

                 ! Exit the newton loop if the update is smaller than the
                 ! threshold value.

                 if(abs(dT) < dTStop) exit

               enddo newton

             endif

             ! To compute the total pressure, the integral of cp/(r*T)
             ! must be computed from T = T to T = Tt. Compute the integrand
             ! at T = Tt; take care of the exceptional situations.

             if(Tt <= cpTrange(0)) then
               intCportTt = (cv0+one)*log(Tt)
             else if( Tt >= cpTrange(cpNparts)) then
               intCportTt = (cvn+one)*log(Tt)
             else

               intCportTt = zero
               do ii=1,cpTempFit(nnt)%nterm
                 if(cpTempFit(nnt)%exponents(ii) == 0_intType) then
                   intCportTt = intCportTt &
                                + cpTempFit(nnt)%constants(ii)*log(Tt)
                 else
                   mm = cpTempFit(nnt)%exponents(ii)
                   T2 = Tt**mm
                   intCportTt = intCportTt &
                              + cpTempFit(nnt)%constants(ii)*T2/mm
                 endif
               enddo

             endif

             ! Compute the integral of cp/(r*T) from T to Tt. First
             ! substract the lower boundary from the upper boundary.

             intCport = intCportTt - intCportT

             ! Add the contributions from the possible internal curve fit
             ! boundaries if Tt and T are in different curve fit intervals.

             do mm=(nn+1),nnt
               ii = mm - 1

               if(ii == 0_intType) then
                 intCport = intCport + (cv0+one)*log(cpTrange(0))
               else
                 intCport = intCport + cpTempFit(ii)%intCpovrt_2
               endif

               if(mm > cpNparts) then
                 intCport = intCport - (cvn+one)*log(cpTrange(cpNparts))
               else
                 intCport = intCport - cpTempFit(mm)%intCpovrt_1
               endif
             enddo

             ! And finally, compute the total pressure.

             ptot(i) = p(i)*exp(intCport)

           enddo

       end select

       end subroutine computePtot

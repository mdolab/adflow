!
!      ******************************************************************
!      *                                                                *
!      * File:          computeHtot.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-07-2004                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeHtot(tt, ht)
!
!      ******************************************************************
!      *                                                                *
!      * computeHtot computes the total enthalpy from the given total   *
!      * temperature. The total enthalpy is the integral of cp, which   *
!      * is a very simple expression for constant cp. For a variable cp *
!      * it is a bit more work.                                         *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use cpCurveFits
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(in)  :: tt
       real(kind=realType), intent(out) :: ht
!
!      Local variables.
!
       integer(kind=intType) :: ii, nn, mm, start

       real(kind=realType) :: t2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the cp model used in the computation.
 
       select case (cpModel)

         case (cpConstant)
 
           ! Constant cp. The total enthalpy is simply cp*tt.

           ht = gammaConstant*RGasDim*tt/(gammaConstant - one)

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.
           ! The actual integral must be computed.

           ! Determine the case we are having here.

           if(tt < cpTrange(0)) then

             ! Temperature is less than the smallest value in the
             ! curve fits. Print a warning and use extrapolation using
             ! constant cp.

             if(myId == 0) then
               print "(a)", "#"
               print "(a)", "#                    Warning"
               print 100, tt, cpTrange(0)
               print "(a)", "# Extrapolation with constant cp is used."
               print "(a)", "#"
 100           format("# Prescribed total temperature ",e12.5,          &
                      " is less than smallest curve fit value, ",e12.5, &
                      ".")
             endif

             ht = RGasDim*(cpEint(0) + tt + cv0*(tt - cpTrange(0)))

           else if(tt > cpTrange(cpNparts)) then

             ! Temperature is larger than the largest value in the
             ! curve fits. Print a warning and use extrapolation using
             ! constant cp.

             if(myId == 0) then
               print "(a)", "#"
               print "(a)", "#                    Warning"
               print 101, tt, cpTrange(cpNparts)
               print "(a)", "# Extrapolation with constant cp is used."
               print "(a)", "#"
 101           format("# Prescribed total temperature ",e12.5,     &
                      " is larger than largest curve fit value, ", &
                      e12.5, ".")
             endif

             ht = RGasDim*(cpEint(cpNparts) + tt &
                +           cvn*(tt - cpTrange(cpNparts)))

           else

             ! Temperature is in the curve fit range.
             ! First find the correct range for this temperature.

             ii    = cpNparts
             start = 1
             interval: do

               ! Next guess for the interval.

               nn = start + ii/2

               ! Determine the situation we are having here.

               if(tt > cpTrange(nn)) then

                 ! Temperature is larger than the upper boundary of
                 ! the current interval. Update the lower boundary.

                 start = nn + 1
                 ii    = ii - 1

               else if(tt >= cpTrange(nn-1)) then

                 ! This is the correct range. Exit the do-loop.

                 exit

               endif

               ! Modify ii for the next branch to search.

               ii = ii/2

             enddo interval

             ! nn contains the correct curve fit interval.
             ! Integrate cp to get ht.

             ht = cpTempFit(nn)%eint0
             do ii=1,cpTempFit(nn)%nterm
               if(cpTempFit(nn)%exponents(ii) == -1_intType) then
                 ht = ht + cpTempFit(nn)%constants(ii)*log(tt)
               else
                 mm = cpTempFit(nn)%exponents(ii) + 1
                 t2 = tt**mm
                 ht = ht + cpTempFit(nn)%constants(ii)*t2/mm
               endif
             enddo

             ! Multiply ht by RGasDim to obtain the correct
             ! dimensional value.

             ht = RGasDim*ht

           endif

       end select

       end subroutine computeHtot

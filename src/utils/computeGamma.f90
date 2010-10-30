!
!      ******************************************************************
!      *                                                                *
!      * File:          computeGamma.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-16-2003                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeGamma(T, gamma, mm)
!
!      ******************************************************************
!      *                                                                *
!      * computeGamma computes the corresponding values of gamma for    *
!      * the given dimensional temperatures.                            *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use cpCurveFits
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(*), intent(in)  :: T
       real(kind=realType), dimension(*), intent(out) :: gamma
       integer(kind=intType), intent(in)              :: mm
!
!      Local variables.
!
       integer(kind=intType) :: i, ii, nn, start
       real(kind=realType)   :: cp, T2
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

           ! Constant cp and thus constant gamma. Set the values.

           do i=1,mm
             gamma(i) = gammaConstant
           enddo

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           do i=1,mm

             ! Determine the case we are having here.

             if(T(i) <= cpTrange(0)) then

               ! Temperature is less than the smallest temperature
               ! in the curve fits. Use cv0 to compute gamma.

               gamma(i) = (cv0 + one)/cv0

             else if(T(i) >= cpTrange(cpNparts)) then

               ! Temperature is larger than the largest temperature
               ! in the curve fits. Use cvn to compute gamma.

               gamma(i) = (cvn + one)/cvn

             else

               ! Temperature is in the curve fit range.
               ! First find the valid range.

               ii    = cpNparts
               start = 1
               interval: do

                 ! Next guess for the interval.

                 nn = start + ii/2

                 ! Determine the situation we are having here.

                 if(T(i) > cpTrange(nn)) then

                   ! Temperature is larger than the upper boundary of
                   ! the current interval. Update the lower boundary.

                   start = nn + 1
                   ii    = ii - 1

                 else if(T(i) >= cpTrange(nn-1)) then

                   ! This is the correct range. Exit the do-loop.

                   exit

                 endif

                 ! Modify ii for the next branch to search.

                 ii = ii/2

               enddo interval

               ! Nn contains the correct curve fit interval.
               ! Compute the value of cp.

               cp = zero
               do ii=1,cpTempFit(nn)%nterm
                 T2 = T(i)**(cpTempFit(nn)%exponents(ii))
                 cp = cp + cpTempFit(nn)%constants(ii)*T2
               enddo

               ! Compute the corresponding value of gamma.

               gamma(i) = cp/(cp-one)

             endif

           enddo

       end select

       end subroutine computeGamma

!
!      ******************************************************************
!      *                                                                *
!      * File:          curveFit.f90                                    *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 03-02-2004                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function curveUpRe(Re)
!
!      ******************************************************************
!      *                                                                *
!      * curveUpRe determines the value of the nonDimensional           *
!      * tangential velocity (made nonDimensional with the skin         *
!      * friction velocity) for the given Reynolds number.              *
!      * This data has been curve fitted with cubic splines.            *
!      *                                                                *
!      ******************************************************************
!
       use paramTurb
       implicit none
!
!      Function type.
!
       real(kind=realType) :: curveUpRe
!
!      Function arguments.
!
       real(kind=realType), intent(in) :: Re
!
!      Local variables.
!
       integer(kind=intType) :: ii, nn, start
       real(kind=realType)   :: x, x2, x3, upRe
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are dealing with.

       if(Re <= reT(0)) then

         ! Reynolds number is less than the smallest number in the curve
         ! fit. Use extrapolation.

         x     = sqrt(Re/reT(0))
         upRe = x*up0(1)

       else if(Re >= reT(nFit)) then

         ! Reynolds number is larger than the largest number in the curve
         ! fit. Set upRe to the largest value available.

         nn = nFit
         x  = reT(nn) - reT(nn-1)
         x2 = x*x
         x3 = x*x2

         upRe = up0(nn) + up1(nn)*x + up2(nn)*x2 + up3(nn)*x3

       else

         ! Reynolds number is in the range of the curve fits.
         ! First find the correct interval.

         ii    = nFit
         start = 1
         interval: do

           ! Next guess for the interval.

           nn = start + ii/2

           ! Determine the situation we are having here.

           if(Re > reT(nn)) then

             ! Reynoldls number is larger than the upper boundary of
             ! the current interval. Update the lower boundary.

             start = nn + 1
             ii    = ii - 1

           else if(Re >= reT(nn-1)) then

             ! This is the correct range. Exit the do-loop.

             exit

           endif

           ! Modify ii for the next branch to search.

           ii = ii/2

         enddo interval

         ! Compute upRe using the cubic polynomial for this interval.

         x  = Re - reT(nn-1)
         x2 = x*x
         x3 = x*x2

         upRe = up0(nn) + up1(nn)*x + up2(nn)*x2 + up3(nn)*x3

       endif

       ! And set the function value.

       curveUpRe = upRe

       end function curveUpRe
!
!      ==================================================================
!
       subroutine curveTupYp(tup, yp, ntu1, ntu2)
!
!      ******************************************************************
!      *                                                                *
!      * CurveTupYp determines the value of the turbulent variables     *
!      * ntu1 to ntu2 for the given yplus.                              *
!      * This data has been curve fitted with cubic splines.            *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputPhysics
       use paramTurb
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: ntu1, ntu2
       real(kind=realType),   intent(in) :: yp

       real(kind=realType), dimension(ntu1:ntu2), intent(out) :: tup
!
!      Local variables.
!
       integer(kind=intType) :: ii, nn, start, mm
       real(kind=realType)   :: x, x2, x3, epsWall, fWall
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are dealing with.

       if(yp <= ypT(0)) then

         ! Yplus is less than the smallest number in the curve
         ! fit. The treatment is turbulence model dependent.

         select case(turbModel)

           case (spalartAllmaras, spalartAllmarasEdwards)

             ! Transport variable is zero on the wall. Use linear
             ! interpolation.

             x = yp/ypT(0)
             do mm=ntu1,ntu2
               tup(mm) = x*tup0(1,mm)
             enddo

           !=============================================================

           case (komegaWilcox, komegaModified, menterSST)

             ! Use the near wall expressions for k and omega.

             x = yp/ypT(0)
             do mm=ntu1,ntu2
               select case(mm)
                 case (itu1)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))*(x**3.23_realType)
                   else
                     tup(mm) = tup0(1,mm)*(x**3.23_realType)
                   endif

                 case (itu2)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))/(max(x,eps)**2)
                   else
                     tup(mm) = tup0(1,mm)/(max(x,eps)**2)
                   endif
               end select
             enddo

           !=============================================================

           case (ktau)

             ! Use the near wall expressions for k and tau.

             x = yp/ypT(0)
             do mm=ntu1,ntu2
               select case(mm)
                 case (itu1)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))*(x**3.23_realType)
                   else
                     tup(mm) = tup0(1,mm)*(x**3.23_realType)
                   endif

                 case (itu2)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))*x*x
                   else
                     tup(mm) = tup0(1,mm)*x*x
                   endif
               end select
             enddo

           !=============================================================

           case (v2f)

             ! Use the near wall expressions for k, epsilon, v2 and f.

             x = yp/ypT(0)
             do mm=ntu1,ntu2
               select case(mm)
                 case (itu1)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))*x**2
                   else
                     tup(mm) = tup0(1,mm)*x**2
                   endif

                 case (itu2)  ! epsilon cannot be fitted logarithmically.
                   if( tuLogFit(mm) ) then
                      call terminate("curveTupYp", &
                                     "Check curveFit, &
                                     &epsilon cannot be fitted with log")
                   else
                      if(rvfN == 1) epsWall = 0.33_realType
                      if(rvfN == 6) epsWall = 0.27_realType
                      tup(mm) = epsWall + (tup0(1,mm)-epsWall)*x
                   endif

                 case (itu3)
                   if( tuLogFit(mm) ) then
                     tup(mm) = exp(tup0(1,mm))*x**4
                   else
                     tup(mm) = tup0(1,mm)*x**4
                   endif

                 case (itu4)
                   if( tuLogFit(mm) ) then
                      if(rvfN == 1) &
                        call terminate("curveTupYp", &
                                     "Check curveFit, &
                                     &f cannot be fitted with log")
                      if(rvfN == 6) tup(mm) = exp(tup(mm))*x
                   else
                      if(rvfN == 1) fWall =-0.0035_realType
                      if(rvfN == 6) fWall = zero
                      tup(mm) = fWall + (tup0(1,mm)-fWall)*x
                   endif

                 case (itu5)
                   if( tuLogFit(mm) ) then
                      tup(mm) = exp(tup(mm))*x**4
                   else
                      tup(mm) = tup0(1,mm)*x**4
                   endif
               end select
             enddo

         end select

       !=================================================================

       else if(yp >= ypT(nFit)) then

         ! Yplus is larger than the largest number in the curve
         ! fit. Set tup to the largest value available.

         nn = nFit
         x  = ypT(nn) - ypT(nn-1)
         x2 = x*x
         x3 = x*x2

         do mm=ntu1,ntu2
           tup(mm) = tup0(nn,mm)    + tup1(nn,mm)*x &
                   + tup2(nn,mm)*x2 + tup3(nn,mm)*x3
           if( tuLogFit(mm) ) tup(mm) = exp(tup(mm))
         enddo

       !=================================================================

       else

         ! y-plus is in the range of the curve fits.
         ! First find the correct interval.

         ii    = nFit
         start = 1
         interval: do

           ! Next guess for the interval.

           nn = start + ii/2

           ! Determine the situation we are having here.

           if(yp > ypT(nn)) then

             ! Yplus is larger than the upper boundary of
             ! the current interval. Update the lower boundary.

             start = nn + 1
             ii    = ii - 1

           else if(yp >= ypT(nn-1)) then

             ! This is the correct range. Exit the do-loop.

             exit

           endif

           ! Modify ii for the next branch to search.

           ii = ii/2

         enddo interval

         ! Compute tup using the cubic polynomial for this interval.

         x  = yp - ypT(nn-1)
         x2 = x*x
         x3 = x*x2

         do mm=ntu1,ntu2
           tup(mm) = tup0(nn,mm)    + tup1(nn,mm)*x &
                   + tup2(nn,mm)*x2 + tup3(nn,mm)*x3
           if( tuLogFit(mm) ) tup(mm) = exp(tup(mm))
         enddo

       endif

       end subroutine curveTupYp

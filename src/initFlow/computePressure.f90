!
!      ******************************************************************
!      *                                                                *
!      * File:          computePressure.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-14-2003                                      *
!      * Last modified: 11-19-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computePressure(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, &
                                  pointerOffset)
!
!      ******************************************************************
!      *                                                                *
!      * computePressure computes the pressure from the total energy,   *
!      * density and velocities in the given cell range of the block to *
!      * which the pointers in blockPointers currently point.           *
!      * It is possible to specify a possible pointer offset, because   *
!      * this routine is also used when reading a restart file.         *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use blockPointers
       use flowVarRefState
       use constants
       use iteration
       use cpCurveFits
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBeg,iEnd,jBeg,jEnd,kBeg,kEnd
       integer(kind=intType), intent(in) :: pointerOffset
!
!      Local parameters.
!
       real(kind=realType), parameter :: dTStop  = 0.01_realType
       real(kind=realType), parameter :: twothird = two*third
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ip, jp, kp, nn, ii, start

       real(kind=realType) :: gm1, factK, v2, scale, e0, e
       real(kind=realType) :: TRefInv, T, dT, T2, alp, cv
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

           ! Constant cp and thus constant gamma. The relation
           ! eint = cv*T can be used and consequently the standard
           ! relation between pressure and internal energy is valid.

           ! Abbreviate some constants that occur in the pressure
           ! computation.

           gm1   = gammaConstant - one
           factK = five*third - gammaConstant

           ! Loop over the cells. Take the possible pointer
           ! offset into account and store the pressure in the
           ! position w(:,:,:, irhoE).

           do k=kBeg,kEnd
             kp = k + pointerOffset
             do j=jBeg,jEnd
               jp = j + pointerOffset
               do i=iBeg,iEnd
                 ip = i + pointerOffset

                 v2 = w(ip,jp,kp,ivx)**2 + w(ip,jp,kp,ivy)**2 &
                    + w(ip,jp,kp,ivz)**2
                 w(i,j,k,irhoE) = gm1*(w(ip,jp,kp,irhoE) &
                                - half*w(ip,jp,kp,irho)*v2)
                 w(i,j,k,irhoE) = max(w(i,j,k,irhoE), &
                                      1.e-5_realType*pInf)
               enddo
             enddo
           enddo

           ! Correct p if a k-equation is present.

           if( kPresent ) then
             do k=kBeg,kEnd
               kp = k + pointerOffset
               do j=jBeg,jEnd
                 jp = j + pointerOffset
                 do i=iBeg,iEnd
                   ip = i + pointerOffset

                   w(i,j,k,irhoE) = w(i,j,k,irhoE) &
                                   + factK*w(ip,jp,kp,irho) &
                                   * w(ip,jp,kp,itu1)
                 enddo
               enddo
             enddo
           endif

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           ! Store a scale factor when converting the nonDimensional
           ! energy to the units of cpEint

           TRefInv = one/TRef
           scale   = TRef/RGas

           ! Loop over the cells to compute the internal energy per
           ! unit mass. This is stored in w(:,:,:,irhoE) for the moment.

           do k=kBeg,kEnd
             kp = k + pointerOffset
             do j=jBeg,jEnd
               jp = j + pointerOffset
               do i=iBeg,iEnd
                 ip = i + pointerOffset

                 w(i,j,k,irhoE) = w(ip,jp,kp,irhoE)/w(ip,jp,kp,irho)
                 if( kPresent ) &
                   w(i,j,k,irhoE) = w(i,j,k,irhoE) - w(ip,jp,kp,itu1)

                 v2 = w(ip,jp,kp,ivx)**2 + w(ip,jp,kp,ivy)**2 &
                    + w(ip,jp,kp,ivz)**2
                 w(i,j,k,irhoE) = w(i,j,k,irhoE) - half*v2
               enddo
             enddo
           enddo

           ! Newton algorithm to compute the temperature from the known
           ! value of the internal energy.

           do k=kBeg,kEnd
             kp = k + pointerOffset
             do j=jBeg,jEnd
               jp = j + pointerOffset
               do i=iBeg,iEnd
                 ip = i + pointerOffset

                 ! Store the internal energy in the same dimensional
                 ! units as cpEint.

                 e0 = scale*w(i,j,k,irhoE)

                 ! Take care of the exceptional cases.

                 if(e0 <= cpEint(0)) then

                   ! Energy smaller than the lowest value of the curve
                   ! fit. Use extrapolation using constant cv.

                   T = TRefInv*(cpTrange(0) + (e0 - cpEint(0))/cv0)

                 else if(e0 >= cpEint(cpNparts)) then

                   ! Energy larger than the largest value of the curve
                   ! fit. Use extrapolation using constant cv.

                   T = TRefInv*(cpTrange(cpNparts) &
                     +         (e0 - cpEint(cpNparts))/cvn)

                 else

                   ! The value is in the range of the curve fits.
                   ! A Newton algorithm is used to find the temperature.

                   ! First find the curve fit interval to be searched.

                   ii    = cpNparts
                   start = 1
                   interval: do

                     ! Next guess for the interval.

                     nn = start + ii/2

                     ! Determine the situation we are having here.

                     if(e0 > cpEint(nn)) then

                       ! Energy is larger than the upper boundary of
                       ! the current interval. Update the lower
                       ! boundary.

                       start = nn + 1
                       ii    = ii - 1

                     else if(e0 >= cpEint(nn-1)) then

                       ! This is the correct range. Exit the do-loop.

                       exit

                     endif

                     ! Modify ii for the next branch to search.

                     ii = ii/2

                   enddo interval

                   ! nn contains the range in which the Newton algorithm
                   ! must be applied.

                   ! Initial guess of the dimensional temperature.

                   alp = (cpEint(nn) - e0)/(cpEint(nn) - cpEint(nn-1))
                   T   = alp*cpTrange(nn-1) + (one-alp)*cpTrange(nn)

                   ! The actual Newton algorithm to compute the
                   ! temperature.

                   Newton: do

                     ! Compute the internal energy as well as the
                     ! value of cv/r for the given temperature.

                     cv = -one                      ! cv/r = cp/r - 1.0
                     e  =  cpTempFit(nn)%eint0 - T  ! e = integral of cv,
                                                    ! Not of cp.
                     do ii=1,cpTempFit(nn)%nterm

                       ! Update cv.

                       T2 = T**(cpTempFit(nn)%exponents(ii))
                       cv = cv + cpTempFit(nn)%constants(ii)*T2

                       ! Update e, for which this contribution must be
                       ! integrated. Take the exceptional case that the
                       ! exponent == -1 into account.

                       if(cpTempFit(nn)%exponents(ii) == -1_intType) then
                         e = e + cpTempFit(nn)%constants(ii)*log(T)
                       else
                         e = e + cpTempFit(nn)%constants(ii)*T2*T &
                           / (cpTempFit(nn)%exponents(ii) + 1)
                       endif

                     enddo

                     ! Compute the update and the new temperature.

                     dT = (e0 - e)/cv
                     T  = T + dT

                     ! Exit the Newton loop if the update is smaller
                     ! than the threshold value.

                     if(abs(dT) < dTStop) exit

                   enddo Newton

                   ! Create the nonDimensional temperature.

                   T = T*TRefInv

                 endif

                 ! Compute the pressure from the known temperature
                 ! and density. Include the correction if a k-equation
                 ! is present.

                 w(i,j,k,irhoE) = w(ip,jp,kp,irho)*RGas*T
                 w(i,j,k,irhoE) = max(w(i,j,k,irhoE), &
                                      1.e-5_realType*pInf)
                 if( kPresent ) &
                   w(i,j,k,irhoE) = w(i,j,k,irhoE) &
                                  + twothird*w(ip,jp,kp,irho) &
                                  * w(ip,jp,kp,itu1)
               enddo
             enddo
           enddo

       end select

       end subroutine computePressure

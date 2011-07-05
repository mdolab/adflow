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
       subroutine computeEtot(iStart,iEnd, jStart,jEnd, kStart, kEnd, &
                              correctForK)
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
       use blockPointers
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iStart,iEnd, jStart,jEnd
       integer(kind=intType), intent(in) :: kStart, kEnd
       logical,               intent(in) :: correctForK
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       real(kind=realType)   :: ovgm1, factK, scale
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

           ! Constant cp and thus constant gamma.
           ! Abbreviate 1/(gamma -1) a bit easier.

           ovgm1 = one/(gammaConstant - one)

           ! Loop over the given range of the block and compute the first
           ! step of the energy.

           do k=kStart,kEnd
             do j=jStart,jEnd
               do i=iStart,iEnd
                 w(i,j,k,irhoE) = ovgm1*p(i,j,k) &
                                + half*w(i,j,k,irho)*(w(i,j,k,ivx)**2 &
                                +                     w(i,j,k,ivy)**2 &
                                +                     w(i,j,k,ivz)**2)
               enddo
             enddo
           enddo

           ! Second step. Correct the energy in case a turbulent kinetic
           ! energy is present.

           if( correctForK ) then

             factK = ovgm1*(five*third - gammaConstant)

             do k=kStart,kEnd
               do j=jStart,jEnd
                 do i=iStart,iEnd
                   w(i,j,k,irhoE) = w(i,j,k,irhoE) &
                                  - factK*w(i,j,k,irho)*w(i,j,k,itu1)
                 enddo
               enddo
             enddo
           endif

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           ! Store a scale factor to compute the nonDimensional
           ! internal energy.

           scale = RGas/Tref

           ! Loop over the given range of the block.

           do k=kStart,kEnd
             do j=jStart,jEnd
               do i=iStart,iEnd
                 call computeEtotCellCpfit(i, j, k, scale, &
                                           correctForK)
               enddo
             enddo
           enddo

       end select
 
       end subroutine computeEtot
 
!      ==================================================================

       subroutine computeEtotBndryList(correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * ComputeEtotBndryList is identical to computeEtot except        *
!      * it computes total energy for the overset bndry index list.     *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: correctForK
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, mm
       real(kind=realType)   :: ovgm1, factK, scale
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

           ! Constant cp and thus constant gamma.
           ! Abbreviate 1/(gamma -1) a bit easier.

           ovgm1 = one/(gammaConstant - one)

           ! Loop over the given indices of the block and compute the first
           ! step of the energy.

           do mm = 1,nCellsOversetAll
             i = ibndry(1,mm)
             j = ibndry(2,mm)
             k = ibndry(3,mm)
             w(i,j,k,irhoE) = ovgm1*p(i,j,k) &
                            + half*w(i,j,k,irho)*(w(i,j,k,ivx)**2 &
                            +                     w(i,j,k,ivy)**2 &
                            +                     w(i,j,k,ivz)**2)
           enddo

           ! Second step. Correct the energy in case a turbulent kinetic
           ! energy is present.

           if( correctForK ) then

             factK = ovgm1*(five*third - gammaConstant)

             do mm = 1,nCellsOversetAll
               i = ibndry(1,mm)
               j = ibndry(2,mm)
               k = ibndry(3,mm)
               w(i,j,k,irhoE) = w(i,j,k,irhoE) &
                              - factK*w(i,j,k,irho)*w(i,j,k,itu1)
             enddo
           endif

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           ! Store a scale factor to compute the nonDimensional
           ! internal energy.

           scale = RGas/Tref

           ! Loop over the given range of the block.

           do mm = 1,nCellsOversetAll
             i = ibndry(1,mm)
             j = ibndry(2,mm)
             k = ibndry(3,mm)
             call computeEtotCellCpfit(i, j, k, scale, &
                                       correctForK)
           enddo

       end select
 
       end subroutine computeEtotBndryList
 
!      ==================================================================

       subroutine computeEtotCellCpfit(i, j, k, scale, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * ComputeEtotCellCpfit will compute the total energy for the     *
!      * given cell of the block given by the current pointers with the *
!      * cp temperature curve fit model.                                *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cpCurveFits
       use flowVarRefState
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: twoThird = two*third
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: i, j, k
       real(kind=realType),   intent(in) :: scale
       logical,                intent(in) :: correctForK
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ii, start

       real(kind=realType) :: pp, t, t2, cv, eint
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the dimensional temperature.

       pp = p(i,j,k)
       if( correctForK ) pp = pp - twoThird &
                                 * w(i,j,k,irho)*w(i,j,k,itu1)
       t = Tref*pp/(RGas*w(i,j,k,irho))

       ! Determine the case we are having here.

       if(t <= cpTrange(0)) then

         ! Temperature is less than the smallest temperature
         ! in the curve fits. Use extrapolation using
         ! constant cv.

         eint = scale*(cpEint(0) + cv0*(t - cpTrange(0)))
         gamma(i,j,k) = (cv0 + one)/cv0

       else if(t >= cpTrange(cpNparts)) then

         ! Temperature is larger than the largest temperature
         ! in the curve fits. Use extrapolation using
         ! constant cv.

         eint = scale*(cpEint(cpNparts) &
              +        cvn*(t - cpTrange(cpNparts)))

         gamma(i,j,k) = (cvn + one)/cvn

       else

         ! Temperature is in the curve fit range.
         ! First find the valid range.

         ii    = cpNparts
         start = 1
         interval: do

           ! Next guess for the interval.

           nn = start + ii/2

           ! Determine the situation we are having here.

           if(t > cpTrange(nn)) then

             ! Temperature is larger than the upper boundary of
             ! the current interval. Update the lower boundary.

             start = nn + 1
             ii    = ii - 1

           else if(t >= cpTrange(nn-1)) then

             ! This is the correct range. Exit the do-loop.

             exit

           endif

           ! Modify ii for the next branch to search.

           ii = ii/2

         enddo interval

         ! Nn contains the correct curve fit interval.
         ! Integrate cv to compute eint.

         eint = cpTempFit(nn)%eint0 - t
         cv   = -one

         do ii=1,cpTempFit(nn)%nterm
           t2 = t**cpTempFit(nn)%exponents(ii)
           cv = cv + cpTempFit(nn)%constants(ii)*t2

           if(cpTempFit(nn)%exponents(ii) == -1_intType) then
             eint = eint + cpTempFit(nn)%constants(ii)*log(t)
           else
             mm   = cpTempFit(nn)%exponents(ii) + 1
             t2   = t*t2
             eint = eint + cpTempFit(nn)%constants(ii)*t2/mm
           endif
         enddo

         eint = scale*eint
         gamma(i,j,k) = (cv + one)/cv

       endif

       ! Compute the total energy per unit volume.

       w(i,j,k,irhoE) = w(i,j,k,irho)*(eint &
                      +                half*(w(i,j,k,ivx)**2 &
                      +                      w(i,j,k,ivy)**2 &
                      +                      w(i,j,k,ivz)**2))

       if( correctForK )                 &
         w(i,j,k,irhoE) = w(i,j,k,irhoE) &
                        + w(i,j,k,irho)*w(i,j,k,itu1)

       end subroutine computeEtotCellCpfit

!      ==================================================================

       subroutine etotArray(rho, u, v, w, p, k, etot, correctForK, kk)
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
       use constants
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(*), intent(in)  :: rho, p, k
       real(kind=realType), dimension(*), intent(in)  :: u, v, w
       real(kind=realType), dimension(*), intent(out) :: etot
       logical, intent(in)                             :: correctForK
       integer(kind=intType), intent(in)              :: kk
!
!      Local variables.
!
       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the internal energy for unit mass.

       call eintArray(rho, p, k, etot, correctForK, kk)

       ! Add the kinetic energy.

       do i=1,kk
         etot(i) = rho(i)*(etot(i) &
                 +         half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)))
       enddo

       end subroutine etotArray

!      ==================================================================

       subroutine eintArray(rho, p, k, eint, correctForK, kk)
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
       use constants
       use cpCurveFits
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(*), intent(in)  :: rho, p, k
       real(kind=realType), dimension(*), intent(out) :: eint
       logical, intent(in)                             :: correctForK
       integer(kind=intType), intent(in)              :: kk
!
!      Local parameter.
!
       real(kind=realType), parameter :: twoThird = two*third
!
!      Local variables.
!
       integer(kind=intType) :: i, nn, mm, ii, start

       real(kind=realType) :: ovgm1, factK, pp, t, t2, scale
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

           ! Abbreviate 1/(gamma -1) a bit easier.

           ovgm1 = one/(gammaConstant - one)

           ! Loop over the number of elements of the array and compute
           ! the total energy.

           do i=1,kk
             eint(i) = ovgm1*p(i)/rho(i)
           enddo

           ! Second step. Correct the energy in case a turbulent kinetic
           ! energy is present.

           if( correctForK ) then

             factK = ovgm1*(five*third - gammaConstant)

             do i=1,kk
               eint(i) = eint(i) - factK*k(i)
             enddo

           endif

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           ! Store a scale factor to compute the nonDimensional
           ! internal energy.

           scale = RGas/Tref

           ! Loop over the number of elements of the array

           do i=1,kk

             ! Compute the dimensional temperature.

             pp = p(i)
             if( correctForK ) pp = pp - twoThird*rho(i)*k(i)
             t = Tref*pp/(RGas*rho(i))

             ! Determine the case we are having here.

             if(t <= cpTrange(0)) then

               ! Temperature is less than the smallest temperature
               ! in the curve fits. Use extrapolation using
               ! constant cv.

               eint(i) = scale*(cpEint(0) + cv0*(t - cpTrange(0)))

             else if(t >= cpTrange(cpNparts)) then

               ! Temperature is larger than the largest temperature
               ! in the curve fits. Use extrapolation using
               ! constant cv.

               eint(i) = scale*(cpEint(cpNparts) &
                       +        cvn*(t - cpTrange(cpNparts)))

             else

               ! Temperature is in the curve fit range.
               ! First find the valid range.

               ii    = cpNparts
               start = 1
               interval: do

                 ! Next guess for the interval.

                 nn = start + ii/2

                 ! Determine the situation we are having here.

                 if(t > cpTrange(nn)) then

                   ! Temperature is larger than the upper boundary of
                   ! the current interval. Update the lower boundary.

                   start = nn + 1
                   ii    = ii - 1

                 else if(t >= cpTrange(nn-1)) then

                   ! This is the correct range. Exit the do-loop.

                   exit

                 endif

                 ! Modify ii for the next branch to search.

                 ii = ii/2

               enddo interval

               ! Nn contains the correct curve fit interval.
               ! Integrate cv to compute eint.

               eint(i) = cpTempFit(nn)%eint0 - t
               do ii=1,cpTempFit(nn)%nterm
                 if(cpTempFit(nn)%exponents(ii) == -1_intType) then
                   eint(i) = eint(i) &
                           + cpTempFit(nn)%constants(ii)*log(t)
                 else
                   mm   = cpTempFit(nn)%exponents(ii) + 1
                   t2   = t**mm
                   eint(i) = eint(i) &
                           + cpTempFit(nn)%constants(ii)*t2/mm
                 endif
               enddo

               eint(i) = scale*eint(i)

             endif

             ! Add the turbulent energy if needed.

             if( correctForK ) eint(i) = eint(i) + k(i)

           enddo

       end select

       end subroutine eintArray

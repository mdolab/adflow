!
!      ******************************************************************
!      *                                                                *
!      * File:          etotArrayAdj.f90                              *
!      * Author:        Edwin van der Weide, C.A.(Sandy) Mader          *
!      * Starting date: 04-25-2008                                      *
!      * Last modified: 04-25-2008                                      *
!      *                                                                *
!      ******************************************************************
!

       subroutine etotArrayForcesAdj(rho, u, v, w, p, k, etot, correctForK, kk)
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
!!$       real(kind=realType), dimension(*), intent(in)  :: rho, p, k
!!$       real(kind=realType), dimension(*), intent(in)  :: u, v, w
!!$       real(kind=realType), dimension(*), intent(out) :: etot
       real(kind=realType), dimension(kk), intent(in)  :: rho, p, k
       real(kind=realType), dimension(kk), intent(in)  :: u, v, w
       real(kind=realType), dimension(kk), intent(out) :: etot
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

       call eintArrayForcesAdj(rho, p, k, etot, correctForK, kk)

       ! Add the kinetic energy.

       do i=1,kk
         etot(i) = rho(i)*(etot(i) &
                 +         half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)))
       enddo

     end subroutine etotArrayForcesAdj

!      ==================================================================

       subroutine eintArrayForcesAdj(rho, p, k, eint, correctForK, kk)
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
       real(kind=realType), dimension(kk), intent(in)  :: rho, p, k
       real(kind=realType), dimension(kk), intent(out) :: eint
!!$       real(kind=realType), dimension(*), intent(in)  :: rho, p, k
!!$       real(kind=realType), dimension(*), intent(out) :: eint
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

       end select
       
     end subroutine eintArrayForcesAdj

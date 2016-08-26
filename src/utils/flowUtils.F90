module flowUtils
contains

  subroutine computeTtot(rho, u, v, w, p, Ttot)
    !
    !       computeTtot computes the total temperature for the given       
    !       pressures, densities and velocities.                           
    !
    use constants
    use inputPhysics, only : cpModel
    use flowVarRefState, only : RGas, gammaInf
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), intent(in)  :: rho, p, u, v, w
    real(kind=realType), intent(out) :: Ttot
    !
    !      Local variables.
    !
    integer(kind=intType) :: i

    real(kind=realType) :: govgm1, T, kin

    ! Determine the cp model used.

    select case (cpModel)

    case (cpConstant)

       ! Constant cp and thus constant gamma. The well-known
       ! formula is valid.

       govgm1 = gammainf/(gammainf-one)
       T       = p/(rho*RGas)
       kin     = half*(u*u + v*v + w*w)
       Ttot = T*(one + rho*kin/(govgm1*p))

       !===============================================================

    case (cpTempCurveFits)

       ! Cp is a function of the temperature. The formula used for
       ! constant cp is not valid anymore and a more complicated
       ! procedure must be followed.

       call terminate("computeTtot", &
            "Variable cp formulation not implemented yet")

    end select

  end subroutine computeTtot

  subroutine computeGamma(T, gamma, mm)
    !
    !       computeGamma computes the corresponding values of gamma for    
    !       the given dimensional temperatures.                            
    !
    use constants
    use cpCurveFits
    use inputPhysics, only : cpModel, gammaConstant
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), dimension(mm), intent(in)  :: T
    real(kind=realType), dimension(mm), intent(out) :: gamma
    integer(kind=intType), intent(in)              :: mm
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, ii, nn, start
    real(kind=realType)   :: cp, T2

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

  subroutine computePtot(rho, u, v, w, p, ptot)
    !
    !       ComputePtot computes the total pressure for the given          
    !       pressures, densities and velocities.                           
    !
    use constants
    use cpCurveFits
    use flowVarRefState, only : tref, RGas, gammaInf
    use inputPhysics, only : cpModel
    implicit none

    real(kind=realType), intent(in)  :: rho, p, u, v, w
    real(kind=realType), intent(out) :: ptot
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
    ! Determine the cp model used.

    select case (cpModel)

    case (cpConstant)

       ! Constant cp and thus constant gamma. The well-known
       ! formula is valid.

       govgm1 = gammaInf/(gammaInf-one)

       kin     = half*(u*u + v*v + w*w)
       ptot = p*((one + rho*kin/(govgm1*p))**govgm1)

       !===============================================================

    case (cpTempCurveFits)

       ! Cp is a function of the temperature. The formula used for
       ! constant cp is not valid anymore and a more complicated
       ! procedure must be followed.


       ! Compute the dimensional temperature and the scale
       ! factor to convert the integral of cp to the correct
       ! nonDimensional value.

       T     = Tref*p/(RGas*rho)
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

       htot = (h + half*(u*u + v*v + w*w))/scale

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

       ptot = p*exp(intCport)

    end select

  end subroutine computePtot

  subroutine computeSpeedOfSoundSquared

    !
    !       computeSpeedOfSoundSquared does what it says.                  
    !
    use constants
    use blockPointers, only : ie, je, ke, w, p, aa, gamma
    use utils, only : getCorrectForK
    implicit none
    !
    !      Local variables.
    !
    real(kind=realType), PARAMETER :: twothird=two*third
    integer(kind=intType) :: i, j, k, ii
    real(kind=realType) :: pp
    logical :: correctForK

    ! Determine if we need to correct for K
    correctForK = getCorrectForK()

    if (correctForK) then 
#ifdef TAPENADE_REVERSE
       !$AD II-LOOP
       do ii=0,ie*je*ke - 1
          i = mod(ii, ie) + 1
          j = mod(ii/ie, je) + 1
          k = ii/(ie*je) + 1
#else
          do k=1,ke
             do j=1,je
                do i=1,ie
#endif             
                   pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
                   aa(i,j,k) = gamma(i,j,k)*pp/w(i,j,k,irho)
#ifdef TAPENADE_REVERSE
                end do
#else
             enddo
          enddo
       enddo
#endif   
    else
#ifdef TAPENADE_REVERSE
       !$AD II-LOOP
       do ii=0,ie*je*ke - 1
          i = mod(ii, ie) + 1
          j = mod(ii/ie, je) + 1
          k = ii/(ie*je) + 1
#else
          do k=1,ke
             do j=1,je
                do i=1,ie
#endif             
                   aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
#ifdef TAPENADE_REVERSE
                end do
#else
             enddo
          enddo
       enddo
#endif   
    end if
  end subroutine computeSpeedOfSoundSquared
 subroutine computeEtotBlock(iStart,iEnd, jStart,jEnd, kStart, kEnd, &
       correctForK)
    !
    !       ComputeEtot computes the total energy from the given density,  
    !       velocity and presssure. For a calorically and thermally        
    !       perfect gas the well-known expression is used; for only a      
    !       thermally perfect gas, cp is a function of temperature, curve  
    !       fits are used and a more complex expression is obtained.       
    !       It is assumed that the pointers in blockPointers already       
    !       point to the correct block.                                    
    !
    use constants
    use blockPointers, only : w, p
    use flowVarRefState, only : RGas, Tref
    use inputPhysics, only : cpModel, gammaConstant
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
40  format (1x,I4,I4,I4,E20.6)

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
#ifndef USE_TAPENADE
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
#endif

    end select
  end subroutine computeEtotBlock

  subroutine etot(rho, u, v, w, p, k, etotal, correctForK)
    !
    !       EtotArray computes the total energy from the given density,    
    !       velocity and presssure.                                        
    !       First the internal energy per unit mass is computed and after  
    !       that the kinetic energy is added as well the conversion to     
    !       energy per unit volume.                                        
    !
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), intent(in)  :: rho, p, k
    real(kind=realType), intent(in)  :: u, v, w
    real(kind=realType), intent(out) :: etotal
    logical, intent(in) :: correctForK

    !
    !      Local variables.
    !
    integer(kind=intType) :: i

    ! Compute the internal energy for unit mass.

    call eint(rho, p, k, etotal, correctForK)
    etotal = rho*(etotal &
         +         half*(u*u + v*v + w*w))

  end subroutine etot

  !      ==================================================================

  subroutine eint(rho, p, k, einternal, correctForK)
    !
    !       EintArray computes the internal energy per unit mass from the  
    !       given density and pressure (and possibly turbulent energy)     
    !       For a calorically and thermally perfect gas the well-known     
    !       expression is used; for only a thermally perfect gas, cp is a  
    !       function of temperature, curve fits are used and a more        
    !       complex expression is obtained.                                
    !
    use constants
    use cpCurveFits
    use flowVarRefState, only : RGas, Tref
    use inputPhysics, only : cpModel, gammaConstant
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), intent(in)  :: rho, p, k
    real(kind=realType), intent(out) :: eInternal
    logical, intent(in)                             :: correctForK
    !
    !      Local parameter.
    !
    real(kind=realType), parameter :: twoThird = two*third
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, nn, mm, ii, start

    real(kind=realType) :: ovgm1, factK, pp, t, t2, scale

    ! Determine the cp model used in the computation.

    select case (cpModel)

    case (cpConstant)

       ! Abbreviate 1/(gamma -1) a bit easier.

       ovgm1 = one/(gammaConstant - one)

       ! Loop over the number of elements of the array and compute
       ! the total energy.

       eInternal = ovgm1*p/rho

       ! Second step. Correct the energy in case a turbulent kinetic
       ! energy is present.

       if( correctForK ) then

          factK = ovgm1*(five*third - gammaConstant)

          eInternal = eInternal - factK*k

       endif

#ifndef USE_TAPENADE       
    case (cpTempCurveFits)

       ! Cp as function of the temperature is given via curve fits.

       ! Store a scale factor to compute the nonDimensional
       ! internal energy.

       scale = RGas/Tref

       ! Loop over the number of elements of the array


       ! Compute the dimensional temperature.

       pp = p
       if( correctForK ) pp = pp - twoThird*rho*k
       t = Tref*pp/(RGas*rho)

       ! Determine the case we are having here.

       if(t <= cpTrange(0)) then

          ! Temperature is less than the smallest temperature
          ! in the curve fits. Use extrapolation using
          ! constant cv.

          eInternal = scale*(cpEint(0) + cv0*(t - cpTrange(0)))

       else if(t >= cpTrange(cpNparts)) then

          ! Temperature is larger than the largest temperature
          ! in the curve fits. Use extrapolation using
          ! constant cv.

          eInternal = scale*(cpEint(cpNparts) &
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
          ! Integrate cv to compute eInternal.

          eInternal = cpTempFit(nn)%eint0 - t
          do ii=1,cpTempFit(nn)%nterm
             if(cpTempFit(nn)%exponents(ii) == -1) then
                eInternal = eInternal &
                     + cpTempFit(nn)%constants(ii)*log(t)
             else
                mm   = cpTempFit(nn)%exponents(ii) + 1
                t2   = t**mm
                eInternal = eInternal &
                     + cpTempFit(nn)%constants(ii)*t2/mm
             endif
          enddo

          eInternal = scale*eInternal

       endif

       ! Add the turbulent energy if needed.

       if( correctForK ) eInternal = eInternal + k

#endif
    end select

  end subroutine eint

  subroutine computePressure(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, &
       pointerOffset)
    !
    !       computePressure computes the pressure from the total energy,   
    !       density and velocities in the given cell range of the block to 
    !       which the pointers in blockPointers currently point.           
    !       It is possible to specify a possible pointer offset, because   
    !       this routine is also used when reading a restart file.         
    !
    use constants
    use inputPhysics, only : cpModel, gammaConstant
    use blockPointers, only : w, p
    use flowVarRefState, only : kPresent, Pinf, RGas, TRef
    use cpCurveFits
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

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef  USE_TAPENADE


  subroutine computeEtotCellCpfit(i, j, k, scale, correctForK)
    !
    !       ComputeEtotCellCpfit will compute the total energy for the     
    !       given cell of the block given by the current pointers with the 
    !       cp temperature curve fit model.                                
    !
    use constants
    use cpCurveFits
    use blockPointers, only : gamma, w, P
    use flowVarRefState, only : Rgas, TRef
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

          if(cpTempFit(nn)%exponents(ii) == -1) then
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


  subroutine updateGamma
    !
    !       This is a utility routine to update the gamma variable from    
    !       from gammaConstant if gammaConstant has changed.               
    !
    use constants
    use blockPointers, only : nDom, gamma
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputPhysics, only : gammaConstant
    use utils, only : setPointers
    implicit none

    integer :: nn, sps

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn, 1_intType, sps)
          gamma = gammaConstant
       end do
    end do
  end subroutine updateGamma
#endif
end module flowUtils

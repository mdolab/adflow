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

  subroutine computePressureSimple(includeHalos) 

    ! Compute the pressure on a block with the pointers already set. This
    ! routine is used by the forward mode AD code only. 

    use constants
    use blockPointers
    use flowVarRefState
    use inputPhysics
    implicit none

    ! Input parameter
    logical, intent(in) :: includeHalos

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii
    real(kind=realType) :: gm1, v2
    integer(kind=intType) :: iBeg, iEnd, iSize, jBeg, jEnd, jSize, kBeg, kEnd, kSize
    ! Compute the pressures
    gm1 = gammaConstant - one

    if (includeHalos) then 
       iBeg = 0
       jBeg = 0
       kBeg = 0
       iEnd = ib
       jEnd = jb
       kEnd = kb
    else
       iBeg = 2
       jBeg = 2
       kBeg = 2
       iEnd = il
       jEnd = jl
       kEnd = kl
    end if

#ifdef TAPENADE_FAST
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/(iSize), jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif             
                v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
                p(i, j, k) = gm1*(w(i, j, k, irhoE) - half*w( i, j, k, irho)*v2)
                p(i, j, k) = max(p(i, j, k), 1.e-4_realType*pInfCorr)

#ifdef TAPENADE_FAST
             end do
#else
          end do
       end do
    end do
#endif
  end subroutine computePressureSimple

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

subroutine computeLamViscosity(includeHalos) 
  !
  !       computeLamViscosity computes the laminar viscosity ratio in    
  !       the owned cell centers of the given block. Sutherland's law is 
  !       used. It is assumed that the pointes already point to the      
  !       correct block before entering this subroutine.                 
  !
  use blockPointers
  use constants
  use flowVarRefState
  use inputPhysics
  use iteration
  use utils, only : getCorrectForK
  implicit none

  ! input variables
  logical, intent(in) :: includeHalos
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: twoThird = two*third
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, ii
  real(kind=realType)   :: muSuth, TSuth, SSuth, T, pp
  logical               :: correctForK
  integer(kind=intType) :: iBeg, iEnd, iSize, jBeg, jEnd, jSize, kBeg, kEnd, kSize

  ! Return immediately if no laminar viscosity needs to be computed.

  if(.not. viscous ) return

  ! Determine whether or not the pressure must be corrected
  ! for the presence of the turbulent kinetic energy.

  correctForK = getCorrectForK()

  ! Compute the nonDimensional constants in sutherland's law.

  muSuth = muSuthDim/muRef
  TSuth  = TSuthDim/Tref
  SSuth  = SSuthDim/Tref

    if (includeHalos) then 
       iBeg = 1
       jBeg = 1
       kBeg = 1
       iEnd = ie
       jEnd = je
       kEnd = ke
    else
       iBeg = 2
       jBeg = 2
       kBeg = 2
       iEnd = il
       jEnd = jl
       kEnd = kl
    end if

  ! Substract 2/3 rho k, which is a part of the normal turbulent
  ! stresses, in case the pressure must be corrected.

  if( correctForK ) then
#ifdef TAPENADE_FAST
     iSize = (iEnd-iBeg)+1
     jSize = (jEnd-jBeg)+1
     kSize = (kEnd-kBeg)+1

     !$AD II-LOOP
     do ii=0, iSize*jSize*kSize-1
        i = mod(ii, iSize) + iBeg
        j = mod(ii/(iSize), jSize) + jBeg
        k = ii/((iSize*jSize)) + kBeg
#else
        do k=kBeg, kEnd
           do j=jBeg, jEnd
              do i=iBeg, iEnd
#endif             
                 pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
                 T = pp/(RGas*w(i,j,k,irho))
                 rlv(i,j,k) = muSuth*((TSuth + SSuth)/(T + SSuth)) &
                      * ((T/TSuth)**1.5_realType)
#ifdef TAPENADE_FAST
              end do
#else
           enddo
        enddo
     enddo
#endif   
  else
     ! Loop over the owned cells *AND* first level halos of this
     ! block and compute the laminar viscosity ratio.
#ifdef TAPENADE_FAST
     iSize = (iEnd-iBeg)+1
     jSize = (jEnd-jBeg)+1
     kSize = (kEnd-kBeg)+1

     !$AD II-LOOP
     do ii=0, iSize*jSize*kSize-1
        i = mod(ii, iSize) + iBeg
        j = mod(ii/(iSize), jSize) + jBeg
        k = ii/((iSize*jSize)) + kBeg
#else
        do k=kBeg, kEnd
           do j=jBeg, jEnd
              do i=iBeg, iEnd
#endif             

                 ! Compute the nonDimensional temperature and the
                 ! nonDimensional laminar viscosity.
                 T = p(i,j,k)/(RGas*w(i,j,k,irho))
                 rlv(i,j,k) = muSuth*((TSuth + SSuth)/(T + SSuth)) &
                      * ((T/TSuth)**1.5_realType)
#ifdef TAPENADE_FAST
              end do
#else
           enddo
        enddo
     enddo
#endif   
  end if
end subroutine computeLamViscosity

subroutine adjustInflowAngle(alpha, beta, liftIndex)

  use constants
  use inputPhysics

  implicit none

  !Subroutine Vars
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) :: liftIndex

  !Local Vars
  real(kind=realType), dimension(3) :: refDirection

  ! Velocity direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alpha, beta, velDirFreestream,&
       liftIndex)

  ! Drag direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alpha, beta, dragDirection, &
       liftIndex)

  ! Lift direction given by the rotation of a unit vector
  ! initially aligned along the positive z-direction (0,0,1)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(liftIndex) = one

  call getDirVector(refDirection, alpha, beta,liftDirection, liftIndex)

end subroutine adjustInflowAngle

       subroutine derivativeRotMatrixRigid(rotationMatrix, &
                                           rotationPoint, t)
!
!       derivativeRotMatrixRigid determines the derivative of the      
!       rotation matrix at the given time for the rigid body rotation, 
!       such that the grid velocities can be determined analytically.  
!       Also the rotation point of the current time level is           
!       determined. This value can change due to translation of the    
!       entire grid.                                                   
!
       use constants
       use flowVarRefState
       use inputMotion
       use monitor
       use utils, only : rigidRotAngle, derivativeRigidRotAngle
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(in) :: t

       real(kind=realType), dimension(3),   intent(out) :: rotationPoint
       real(kind=realType), dimension(3,3), intent(out) :: rotationMatrix
!
!      Local variables.
!
       integer(kind=intType) :: i, j

       real(kind=realType) :: phi, dphiX, dphiY, dphiZ
       real(kind=realType) :: cosX, cosY, cosZ, sinX, sinY, sinZ

       real(kind=realType), dimension(3,3) :: dm, m

       ! Determine the rotation angle around the x-axis for the new
       ! time level and the corresponding values of the sine and cosine.

       phi = rigidRotAngle(degreePolXRot,   coefPolXRot,      &
                           degreeFourXRot,  omegaFourXRot,    &
                           cosCoefFourXRot, sinCoefFourXRot, t)
       sinX = sin(phi)
       cosX = cos(phi)
       ! Idem for the y-axis.

       phi = rigidRotAngle(degreePolYRot,   coefPolYRot,      &
                           degreeFourYRot,  omegaFourYRot,    &
                           cosCoefFourYRot, sinCoefFourYRot, t)
       sinY = sin(phi)
       cosY = cos(phi)
       ! Idem for the z-axis.

       phi = rigidRotAngle(degreePolZRot,   coefPolZRot,      &
                           degreeFourZRot,  omegaFourZRot,    &
                           cosCoefFourZRot, sinCoefFourZRot, t)
       sinZ = sin(phi)
       cosZ = cos(phi)
       ! Compute the time derivative of the rotation angles around the
       ! x-axis, y-axis and z-axis.

       dphiX = derivativeRigidRotAngle(degreePolXRot,   &
                                       coefPolXRot,     &
                                       degreeFourXRot,  &
                                       omegaFourXRot,   &
                                       cosCoefFourXRot, &
                                       sinCoefFourXRot, t)

       dphiY = derivativeRigidRotAngle(degreePolYRot,   &
                                       coefPolYRot,     &
                                       degreeFourYRot,  &
                                       omegaFourYRot,   &
                                       cosCoefFourYRot, &
                                       sinCoefFourYRot, t)

       dphiZ = derivativeRigidRotAngle(degreePolZRot,   &
                                       coefPolZRot,     &
                                       degreeFourZRot,  &
                                       omegaFourZRot,   &
                                       cosCoefFourZRot, &
                                       sinCoefFourZRot, t)

       ! Compute the time derivative of the rotation matrix applied to
       ! the coordinates at t == 0.

       ! Part 1. Derivative of the z-rotation matrix multiplied by the
       ! x and y rotation matrix, i.e. dmz * my * mx

       dm(1,1) = -cosY*sinZ*dphiZ
       dm(1,2) = (-cosX*cosZ - sinX*sinY*sinZ)*dphiZ
       dm(1,3) = ( sinX*cosZ - cosX*sinY*sinZ)*dphiZ

       dm(2,1) = cosY*cosZ*dphiZ
       dm(2,2) = (sinX*sinY*cosZ - cosX*sinZ)*dphiZ
       dm(2,3) = (cosX*sinY*cosZ + sinX*sinZ)*dphiZ

       dm(3,1) = zero
       dm(3,2) = zero
       dm(3,3) = zero

       ! Part 2: mz * dmy * mx.

       dm(1,1) = dm(1,1) - sinY*cosZ*dphiY
       dm(1,2) = dm(1,2) + sinX*cosY*cosZ*dphiY
       dm(1,3) = dm(1,3) + cosX*cosY*cosZ*dphiY

       dm(2,1) = dm(2,1) - sinY*sinZ*dphiY
       dm(2,2) = dm(2,2) + sinX*cosY*sinZ*dphiY
       dm(2,3) = dm(2,3) + cosX*cosY*sinZ*dphiY

       dm(3,1) = dm(3,1) - cosY*dphiY
       dm(3,2) = dm(3,2) - sinX*sinY*dphiY
       dm(3,3) = dm(3,3) - cosX*sinY*dphiY

       ! Part 3: mz * my * dmx

       dm(1,2) = dm(1,2) + (sinX*sinZ + cosX*sinY*cosZ)*dphiX
       dm(1,3) = dm(1,3) + (cosX*sinZ - sinX*sinY*cosZ)*dphiX

       dm(2,2) = dm(2,2) + (cosX*sinY*sinZ - sinX*cosZ)*dphiX
       dm(2,3) = dm(2,3) - (sinX*sinY*sinZ + cosX*cosZ)*dphiX

       dm(3,2) = dm(3,2) + cosX*cosY*dphiX
       dm(3,3) = dm(3,3) - sinX*cosY*dphiX

       ! Determine the rotation matrix at t == t.

       m(1,1) =  cosY*cosZ
       m(2,1) =  cosY*sinZ
       m(3,1) = -sinY

       m(1,2) = sinX*sinY*cosZ - cosX*sinZ
       m(2,2) = sinX*sinY*sinZ + cosX*cosZ
       m(3,2) = sinX*cosY

       m(1,3) = cosX*sinY*cosZ + sinX*sinZ
       m(2,3) = cosX*sinY*sinZ - sinX*cosZ
       m(3,3) = cosX*cosY

       ! Determine the matrix product dm * inverse(m), which will give
       ! the derivative of the rotation matrix when applied to the
       ! current coordinates. Note that inverse(m) == transpose(m).

       do j=1,3
         do i=1,3
           rotationMatrix(i,j) = dm(i,1)*m(j,1) + dm(i,2)*m(j,2) &
                               + dm(i,3)*m(j,3)
         enddo
       enddo

       ! Determine the rotation point at the new time level; it is
       ! possible that this value changes due to translation of the grid.

    !  aInf = sqrt(gammaInf*pInf/rhoInf)

    !  RotationPoint(1) = LRef*rotPoint(1) &
    !                    + MachGrid(1)*aInf*t/timeRef
    !  rotationPoint(2) = LRef*rotPoint(2) &
    !                    + MachGrid(2)*aInf*t/timeRef
    !  rotationPoint(3) = LRef*rotPoint(3) &
    !                    + MachGrid(3)*aInf*t/timeRef

       rotationPoint(1) = LRef*rotPoint(1)
       rotationPoint(2) = LRef*rotPoint(2)
       rotationPoint(3) = LRef*rotPoint(3)

       end subroutine derivativeRotMatrixRigid


      subroutine getDirVector(refDirection, alpha, beta,&
           windDirection,liftIndex)
        !(xb,yb,zb,alpha,beta,xw,yw,zw)
!
!      Convert the angle of attack and side slip angle to wind axes.  
!      The components of the wind direction vector (xw,yw,zw) are     
!      computed given the direction angles in radians and the body    
!      direction by performing two rotations on the original          
!      direction vector:                                              
!        1) Rotation about the zb or yb-axis: alpha clockwise (CW)    
!           (xb,yb,zb) -> (x1,y1,z1)                                  
!        2) Rotation about the yl or z1-axis: beta counter-clockwise  
!           (CCW)  (x1,y1,z1) -> (xw,yw,zw)                           
!         input arguments:                                            
!            alpha    = angle of attack in radians                    
!            beta     = side slip angle in radians                    
!            refDirection = reference direction vector                
!         output arguments:                                           
!            windDirection = unit wind vector in body axes            
!
      use constants
      use utils, only : terminate
      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType),dimension(3), intent(in)  :: refDirection
      real(kind=realType)  :: alpha, beta
      real(kind=realType),dimension (3), intent(out) :: windDirection
      integer(kind=intType)::liftIndex
!
!     Local variables.
!
      real(kind=realType) :: rnorm,x1,y1,z1,xbn,ybn,zbn,xw,yw,zw
      real(kind=realType) :: tmp

      ! Normalize the input vector.

      rnorm = sqrt( refDirection(1)**2 + refDirection(2)**2 + refDirection(3)**2 )
      xbn = refDirection(1) / rnorm
      ybn = refDirection(2) / rnorm
      zbn = refDirection(3) / rnorm

!!$      ! Compute the wind direction vector.
!!$
!!$      ! 1) rotate alpha radians cw about y-axis
!!$      !    ( <=> rotate y-axis alpha radians ccw)
!!$
!!$      call vectorRotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
!!$
!!$      ! 2) rotate beta radians ccw about z-axis
!!$      !    ( <=> rotate z-axis -beta radians ccw)
!!$
!!$      call vectorRotation(xw, yw, zw, 3, -beta, x1, y1, z1)

      if (liftIndex==2)then
         ! Compute the wind direction vector.Aerosurf axes different!!
         
         ! 1) rotate alpha radians cw about z-axis
         !    ( <=> rotate z-axis alpha radians ccw)
         
         tmp = -alpha
         call vectorRotation(x1, y1, z1, 3, tmp, xbn, ybn, zbn)
         
         ! 2) rotate beta radians ccw about y-axis
         !    ( <=> rotate z-axis -beta radians ccw)
         tmp = -beta
         call vectorRotation(xw, yw, zw, 2, tmp, x1, y1, z1)

      elseif(liftIndex==3)then
         ! Compute the wind direction vector.Aerosurf axes different!!
         
         ! 1) rotate alpha radians cw about z-axis
         !    ( <=> rotate z-axis alpha radians ccw)
         
         call vectorRotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
         
         ! 2) rotate beta radians ccw about y-axis
         !    ( <=> rotate z-axis -beta radians ccw)
         
         call vectorRotation(xw, yw, zw, 3, beta, x1, y1, z1)
         
         
      else
         call terminate('getDirVector', 'Invalid Lift Direction')
         
      endif
      
      windDirection(1) = xw
      windDirection(2) = yw
      windDirection(3) = zw
      
      end subroutine getDirVector
      subroutine vectorRotation(xp, yp, zp, iaxis, angle, x, y, z)
!
!      vectorRotation rotates a given vector with respect to a      
!      specified axis by a given angle.                             
!         input arguments:                                          
!            iaxis      = rotation axis (1-x, 2-y, 3-z)             
!            angle      = rotation angle (measured ccw in radians)  
!            x, y, z    = coordinates in original system            
!         output arguments:                                         
!            xp, yp, zp = coordinates in rotated system             
!
      use precision
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: iaxis
      real(kind=realType), intent(in)   :: angle, x, y, z
      real(kind=realType), intent(out)  :: xp, yp, zp

      ! rotation about specified axis by specified angle

      select case(iaxis)

        ! rotation about the x-axis

        case(1)
          xp =        1.    * x +     0.     * y +     0.     * z
          yp =        0.    * x + cos(angle) * y + sin(angle) * z
          zp =        0.    * x - sin(angle) * y + cos(angle) * z

        ! rotation about the y-axis

        case(2)
          xp =   cos(angle) * x +     0.     * y - sin(angle) * z
          yp =        0.    * x +     1.     * y +     0.     * z
          zp =   sin(angle) * x +     0.     * y + cos(angle) * z

        ! rotation about the z-axis

        case(3)
          xp =   cos(angle) * x + sin(angle) * y +     0.     * z
          yp = - sin(angle) * x + cos(angle) * y +     0.     * z
          zp =       0.     * x +     0.     * y +     1.     * z

        case default
          write(*,*) "vectorRotation called with invalid arguments"
          stop

      end select

      end subroutine vectorRotation

subroutine allNodalGradients
  !
  !         nodalGradients computes the nodal velocity gradients and     
  !         minus the gradient of the speed of sound squared. The minus  
  !         sign is present, because this is the definition of the heat  
  !         flux. These gradients are computed for all nodes.            
  !
  use constants
  use blockPointers
  implicit none
  !        Local variables.
  integer(kind=intType) :: i, j, k
  integer(kind=intType) :: k1, kk
  integer(kind=intType) :: istart, iend, isize, ii
  integer(kind=intType) :: jstart, jend, jsize
  integer(kind=intType) :: kstart, kend, ksize

  real(kind=realType) :: oneOverV, ubar, vbar, wbar, a2
  real(kind=realType) :: sx, sx1, sy, sy1, sz, sz1


  ! Zero all nodeal gradients:
  ux = zero; uy = zero; uz = zero;
  vx = zero; vy = zero; vz = zero;
  wx = zero; wy = zero; wz = zero;
  qx = zero; qy = zero; qz = zero;

  ! First part. Contribution in the k-direction.
  ! The contribution is scattered to both the left and right node
  ! in k-direction.

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,il*jl*ke-1
     i = mod(ii, il) + 1
     j = mod(ii/il, jl) + 1
     k = ii/(il*jl) + 1
#else
     do k=1, ke
        do j=1, jl
           do i=1, il
#endif      
              ! Compute 8 times the average normal for this part of
              ! the control volume. The factor 8 is taken care of later
              ! on when the division by the volume takes place.

              sx = sk(i,j,k-1,  1) + sk(i+1,j,k-1,  1) &
                   + sk(i,j+1,k-1,1) + sk(i+1,j+1,k-1,1) &
                   + sk(i,j,  k,  1) + sk(i+1,j,  k,  1) &
                   + sk(i,j+1,k  ,1) + sk(i+1,j+1,k  ,1)
              sy = sk(i,j,k-1,  2) + sk(i+1,j,k-1,  2) &
                   + sk(i,j+1,k-1,2) + sk(i+1,j+1,k-1,2) &
                   + sk(i,j,  k,  2) + sk(i+1,j,  k,  2) &
                   + sk(i,j+1,k  ,2) + sk(i+1,j+1,k  ,2)
              sz = sk(i,j,k-1,  3) + sk(i+1,j,k-1,  3) &
                   + sk(i,j+1,k-1,3) + sk(i+1,j+1,k-1,3) &
                   + sk(i,j,  k,  3) + sk(i+1,j,  k,  3) &
                   + sk(i,j+1,k  ,3) + sk(i+1,j+1,k  ,3)

              ! Compute the average velocities and speed of sound squared
              ! for this integration point. Node that these variables are
              ! stored in w(ivx), w(ivy), w(ivz) and p.

              ubar = fourth*(w(i,j,  k,ivx) + w(i+1,j,  k,ivx) &
                   +         w(i,j+1,k,ivx) + w(i+1,j+1,k,ivx))
              vbar = fourth*(w(i,j,  k,ivy) + w(i+1,j,  k,ivy) &
                   +         w(i,j+1,k,ivy) + w(i+1,j+1,k,ivy))
              wbar = fourth*(w(i,j,  k,ivz) + w(i+1,j,  k,ivz) &
                   +         w(i,j+1,k,ivz) + w(i+1,j+1,k,ivz))

              a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j+1,k) + aa(i+1,j+1,k))


              ! Add the contributions to the surface integral to the node
              ! j-1 and substract it from the node j. For the heat flux it
              ! is reversed, because the negative of the gradient of the
              ! speed of sound must be computed.

              if(k > 1) then
                 ux(i,j,k-1) = ux(i,j,k-1) + ubar*sx
                 uy(i,j,k-1) = uy(i,j,k-1) + ubar*sy
                 uz(i,j,k-1) = uz(i,j,k-1) + ubar*sz

                 vx(i,j,k-1) = vx(i,j,k-1) + vbar*sx
                 vy(i,j,k-1) = vy(i,j,k-1) + vbar*sy
                 vz(i,j,k-1) = vz(i,j,k-1) + vbar*sz

                 wx(i,j,k-1) = wx(i,j,k-1) + wbar*sx
                 wy(i,j,k-1) = wy(i,j,k-1) + wbar*sy
                 wz(i,j,k-1) = wz(i,j,k-1) + wbar*sz

                 qx(i,j,k-1) = qx(i,j,k-1) - a2*sx
                 qy(i,j,k-1) = qy(i,j,k-1) - a2*sy
                 qz(i,j,k-1) = qz(i,j,k-1) - a2*sz
              endif

              if(k < ke) then
                 ux(i,j,k) = ux(i,j,k) - ubar*sx
                 uy(i,j,k) = uy(i,j,k) - ubar*sy
                 uz(i,j,k) = uz(i,j,k) - ubar*sz

                 vx(i,j,k) = vx(i,j,k) - vbar*sx
                 vy(i,j,k) = vy(i,j,k) - vbar*sy
                 vz(i,j,k) = vz(i,j,k) - vbar*sz

                 wx(i,j,k) = wx(i,j,k) - wbar*sx
                 wy(i,j,k) = wy(i,j,k) - wbar*sy
                 wz(i,j,k) = wz(i,j,k) - wbar*sz

                 qx(i,j,k) = qx(i,j,k) + a2*sx
                 qy(i,j,k) = qy(i,j,k) + a2*sy
                 qz(i,j,k) = qz(i,j,k) + a2*sz
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif   

  ! Second part. Contribution in the j-direction.
  ! The contribution is scattered to both the left and right node
  ! in j-direction.

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,il*je*kl-1
     i = mod(ii, il) + 1
     j = mod(ii/il, je) + 1
     k = ii/(il*je) + 1
#else
     do k=1, kl
        do j=1, je
           do i=1, il
#endif   

              ! Compute 8 times the average normal for this part of
              ! the control volume. The factor 8 is taken care of later
              ! on when the division by the volume takes place.

              sx = sj(i,j-1,k,  1) + sj(i+1,j-1,k,  1) &
                   + sj(i,j-1,k+1,1) + sj(i+1,j-1,k+1,1) &
                   + sj(i,j,  k,  1) + sj(i+1,j,  k,  1) &
                   + sj(i,j,  k+1,1) + sj(i+1,j,  k+1,1)
              sy = sj(i,j-1,k,  2) + sj(i+1,j-1,k,  2) &
                   + sj(i,j-1,k+1,2) + sj(i+1,j-1,k+1,2) &
                   + sj(i,j,  k,  2) + sj(i+1,j,  k,  2) &
                   + sj(i,j,  k+1,2) + sj(i+1,j,  k+1,2)
              sz = sj(i,j-1,k,  3) + sj(i+1,j-1,k,  3) &
                   + sj(i,j-1,k+1,3) + sj(i+1,j-1,k+1,3) &
                   + sj(i,j,  k,  3) + sj(i+1,j,  k,  3) &
                   + sj(i,j,  k+1,3) + sj(i+1,j,  k+1,3)

              ! Compute the average velocities and speed of sound squared
              ! for this integration point. Node that these variables are
              ! stored in w(ivx), w(ivy), w(ivz) and p.

              ubar = fourth*(w(i,j,k,  ivx) + w(i+1,j,k,  ivx) &
                   +         w(i,j,k+1,ivx) + w(i+1,j,k+1,ivx))
              vbar = fourth*(w(i,j,k,  ivy) + w(i+1,j,k,  ivy) &
                   +         w(i,j,k+1,ivy) + w(i+1,j,k+1,ivy))
              wbar = fourth*(w(i,j,k,  ivz) + w(i+1,j,k,  ivz) &
                   +         w(i,j,k+1,ivz) + w(i+1,j,k+1,ivz))

              a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j,k+1) + aa(i+1,j,k+1))

              ! Add the contributions to the surface integral to the node
              ! j-1 and substract it from the node j. For the heat flux it
              ! is reversed, because the negative of the gradient of the
              ! speed of sound must be computed.

              if(j > 1) then
                 ux(i,j-1,k) = ux(i,j-1,k) + ubar*sx
                 uy(i,j-1,k) = uy(i,j-1,k) + ubar*sy
                 uz(i,j-1,k) = uz(i,j-1,k) + ubar*sz

                 vx(i,j-1,k) = vx(i,j-1,k) + vbar*sx
                 vy(i,j-1,k) = vy(i,j-1,k) + vbar*sy
                 vz(i,j-1,k) = vz(i,j-1,k) + vbar*sz

                 wx(i,j-1,k) = wx(i,j-1,k) + wbar*sx
                 wy(i,j-1,k) = wy(i,j-1,k) + wbar*sy
                 wz(i,j-1,k) = wz(i,j-1,k) + wbar*sz

                 qx(i,j-1,k) = qx(i,j-1,k) - a2*sx
                 qy(i,j-1,k) = qy(i,j-1,k) - a2*sy
                 qz(i,j-1,k) = qz(i,j-1,k) - a2*sz
              endif

              if(j < je) then
                 ux(i,j,k) = ux(i,j,k) - ubar*sx
                 uy(i,j,k) = uy(i,j,k) - ubar*sy
                 uz(i,j,k) = uz(i,j,k) - ubar*sz

                 vx(i,j,k) = vx(i,j,k) - vbar*sx
                 vy(i,j,k) = vy(i,j,k) - vbar*sy
                 vz(i,j,k) = vz(i,j,k) - vbar*sz

                 wx(i,j,k) = wx(i,j,k) - wbar*sx
                 wy(i,j,k) = wy(i,j,k) - wbar*sy
                 wz(i,j,k) = wz(i,j,k) - wbar*sz

                 qx(i,j,k) = qx(i,j,k) + a2*sx
                 qy(i,j,k) = qy(i,j,k) + a2*sy
                 qz(i,j,k) = qz(i,j,k) + a2*sz
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif 

  ! Third part. Contribution in the i-direction.
  ! The contribution is scattered to both the left and right node
  ! in i-direction.

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,ie*jl*kl-1
     i = mod(ii, ie) + 1
     j = mod(ii/ie, jl) + 1
     k = ii/(ie*jl) + 1
#else
     do k=1,kl
        do j=1,jl
           do i=1,ie
#endif   

              ! Compute 8 times the average normal for this part of
              ! the control volume. The factor 8 is taken care of later
              ! on when the division by the volume takes place.

              sx = si(i-1,j,k,  1) + si(i-1,j+1,k,  1) &
                   + si(i-1,j,k+1,1) + si(i-1,j+1,k+1,1) &
                   + si(i,  j,k,  1) + si(i,  j+1,k,  1) &
                   + si(i,  j,k+1,1) + si(i,  j+1,k+1,1)
              sy = si(i-1,j,k,  2) + si(i-1,j+1,k,  2) &
                   + si(i-1,j,k+1,2) + si(i-1,j+1,k+1,2) &
                   + si(i,  j,k,  2) + si(i,  j+1,k,  2) &
                   + si(i,  j,k+1,2) + si(i,  j+1,k+1,2)
              sz = si(i-1,j,k,  3) + si(i-1,j+1,k,  3) &
                   + si(i-1,j,k+1,3) + si(i-1,j+1,k+1,3) &
                   + si(i,  j,k,  3) + si(i,  j+1,k,  3) &
                   + si(i,  j,k+1,3) + si(i,  j+1,k+1,3)

              ! Compute the average velocities and speed of sound squared
              ! for this integration point. Node that these variables are
              ! stored in w(ivx), w(ivy), w(ivz) and p.

              ubar = fourth*(w(i,j,k,  ivx) + w(i,j+1,k,  ivx) &
                   +         w(i,j,k+1,ivx) + w(i,j+1,k+1,ivx))
              vbar = fourth*(w(i,j,k,  ivy) + w(i,j+1,k,  ivy) &
                   +         w(i,j,k+1,ivy) + w(i,j+1,k+1,ivy))
              wbar = fourth*(w(i,j,k,  ivz) + w(i,j+1,k,  ivz) &
                   +         w(i,j,k+1,ivz) + w(i,j+1,k+1,ivz))

              a2 = fourth*(aa(i,j,k) + aa(i,j+1,k) + aa(i,j,k+1) + aa(i,j+1,k+1))

              ! Add the contributions to the surface integral to the node
              ! j-1 and substract it from the node j. For the heat flux it
              ! is reversed, because the negative of the gradient of the
              ! speed of sound must be computed.

              if(i > 1) then
                 ux(i-1,j,k) = ux(i-1,j,k) + ubar*sx
                 uy(i-1,j,k) = uy(i-1,j,k) + ubar*sy
                 uz(i-1,j,k) = uz(i-1,j,k) + ubar*sz

                 vx(i-1,j,k) = vx(i-1,j,k) + vbar*sx
                 vy(i-1,j,k) = vy(i-1,j,k) + vbar*sy
                 vz(i-1,j,k) = vz(i-1,j,k) + vbar*sz

                 wx(i-1,j,k) = wx(i-1,j,k) + wbar*sx
                 wy(i-1,j,k) = wy(i-1,j,k) + wbar*sy
                 wz(i-1,j,k) = wz(i-1,j,k) + wbar*sz

                 qx(i-1,j,k) = qx(i-1,j,k) - a2*sx
                 qy(i-1,j,k) = qy(i-1,j,k) - a2*sy
                 qz(i-1,j,k) = qz(i-1,j,k) - a2*sz
              endif

              if(i < ie) then
                 ux(i,j,k) = ux(i,j,k) - ubar*sx
                 uy(i,j,k) = uy(i,j,k) - ubar*sy
                 uz(i,j,k) = uz(i,j,k) - ubar*sz

                 vx(i,j,k) = vx(i,j,k) - vbar*sx
                 vy(i,j,k) = vy(i,j,k) - vbar*sy
                 vz(i,j,k) = vz(i,j,k) - vbar*sz

                 wx(i,j,k) = wx(i,j,k) - wbar*sx
                 wy(i,j,k) = wy(i,j,k) - wbar*sy
                 wz(i,j,k) = wz(i,j,k) - wbar*sz

                 qx(i,j,k) = qx(i,j,k) + a2*sx
                 qy(i,j,k) = qy(i,j,k) + a2*sy
                 qz(i,j,k) = qz(i,j,k) + a2*sz
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif   
  ! Divide by 8 times the volume to obtain the correct gradients.

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,il*jl*kl-1
     i = mod(ii, il) + 1
     j = mod(ii/il, jl) + 1
     k = ii/(il*jl) + 1
#else
     do k=1,kl
        do j=1,jl
           do i=1,il
#endif  
              ! Compute the inverse of 8 times the volume for this node.

              oneOverV = one/(vol(i,  j,  k) + vol(i,  j,  k+1) &
                   +      vol(i+1,j,  k) + vol(i+1,j,  k+1) &
                   +      vol(i,  j+1,k) + vol(i,  j+1,k+1) &
                   +      vol(i+1,j+1,k) + vol(i+1,j+1,k+1))

              ! Compute the correct velocity gradients and "unit" heat
              ! fluxes. The velocity gradients are stored in ux, etc.

              ux(i,j,k) = ux(i,j,k)*oneOverV
              uy(i,j,k) = uy(i,j,k)*oneOverV
              uz(i,j,k) = uz(i,j,k)*oneOverV

              vx(i,j,k) = vx(i,j,k)*oneOverV
              vy(i,j,k) = vy(i,j,k)*oneOverV
              vz(i,j,k) = vz(i,j,k)*oneOverV

              wx(i,j,k) = wx(i,j,k)*oneOverV
              wy(i,j,k) = wy(i,j,k)*oneOverV
              wz(i,j,k) = wz(i,j,k)*oneOverV

              qx(i,j,k) = qx(i,j,k)*oneOverV
              qy(i,j,k) = qy(i,j,k)*oneOverV
              qz(i,j,k) = qz(i,j,k)*oneOverV
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif
end subroutine allNodalGradients

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

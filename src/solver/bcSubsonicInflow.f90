!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSubsonicInflow.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-03-2003                                      *
!      * Last modified: 09-27-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcSubsonicInflow(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcSubsonicInflow applies the subsonic outflow boundary         *
!      * condition, total pressure, total density and flow direction    *
!      * prescribed,  to a block. It is assumed that the pointers in    *
!      * blockPointers are already set to the correct block on the      *
!      * correct grid level.                                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
       use flowVarRefState
       use inputDiscretization
       use inputPhysics
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo, correctForK
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, nn

       real(kind=realType) :: gm1, ovgm1
       real(kind=realType) :: ptot, ttot, htot, a2tot, r, alpha, beta
       real(kind=realType) :: aa, bb, cc, dd, q, q2, a2, m2, scaleFact
       real(kind=realType) :: ssx, ssy, ssz, nnx, nny, nnz
       real(kind=realType) :: rho, velx, vely, velz

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: gamma2
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!
!      Interfaces
!
       interface
         subroutine setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                  rev1, rev2, offset)
           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn, offset
           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
         end subroutine setBCPointers
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for the subsonic inflow boundary condition.

         inflowSubsonic: if(BCType(nn) == SubsonicInflow) then

           ! Nullify the pointers and set them to the correct subface.
           ! They are nullified first, because some compilers require
           ! that.

           nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
           call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0_intType)

           ! Set the additional pointer for gamma2.

           select case (BCFaceID(nn))
             case (iMin)
               gamma2 => gamma(2,1:,1:)
             case (iMax)
               gamma2 => gamma(il,1:,1:)
             case (jMin)
               gamma2 => gamma(1:,2,1:)
             case (jMax)
               gamma2 => gamma(1:,jl,1:)
             case (kMin)
               gamma2 => gamma(1:,1:,2)
             case (kMax)
               gamma2 => gamma(1:,1:,kl)
           end select

           ! Determine the boundary treatment to be used.

           select case (BCData(nn)%subsonicInletTreatment)

             case (totalConditions)

               ! The total conditions have been prescribed.

               ! Loop over the generic subface to set the state in the
               ! halo cells.

               do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                   ! Store a couple of variables, such as the total
                   ! pressure, total temperature, total enthalpy, flow
                   ! direction and grid unit outward normal, a bit easier.

                   ptot = BCData(nn)%ptInlet(i,j)
                   ttot = BCData(nn)%ttInlet(i,j)
                   htot = BCData(nn)%htInlet(i,j)

                   ssx  = BCData(nn)%flowXdirInlet(i,j)
                   ssy  = BCData(nn)%flowYdirInlet(i,j)
                   ssz  = BCData(nn)%flowZdirInlet(i,j)

                   nnx = BCData(nn)%norm(i,j,1)
                   nny = BCData(nn)%norm(i,j,2)
                   nnz = BCData(nn)%norm(i,j,3)

                   ! Some abbreviations in which gamma occurs.

                   gm1   = gamma2(i,j) - one
                   ovgm1 = one/gm1

                   ! Determine the acoustic Riemann variable that must be
                   ! extrapolated from the domain.

                   r    = one/ww2(i,j,irho)
                   a2   = gamma2(i,j)*pp2(i,j)*r
                   beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
                        + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

                   ! Correct the value of the Riemann invariant if total
                   ! enthalpy scaling must be applied. This scaling may
                   ! be needed for stability if large gradients of the
                   ! total temperature are prescribed.

                   scaleFact = one
                   if( hScalingInlet ) &
                    scaleFact = sqrt(htot/(r*(ww2(i,j,irhoE) + pp2(i,j))))

                   beta = beta*scaleFact

                   ! Compute the value of a2 + 0.5*gm1*q2, which is the
                   ! total speed of sound for constant cp. However, the
                   ! expression below is also valid for variable cp,
                   ! although a linearization around the value of the
                   ! internal cell is performed.

                   q2    = ww2(i,j,ivx)**2 + ww2(i,j,ivy)**2 &
                         + ww2(i,j,ivz)**2
                   a2tot = gm1*(htot - r*(ww2(i,j,irhoE) + pp2(i,j)) &
                         +      half*q2) + a2

                   ! Compute the dot product between the normal and the
                   ! velocity direction. This value should be negative.

                   alpha = nnx*ssx + nny*ssy + nnz*ssz

                   ! Compute the coefficients in the quadratic equation
                   ! for the magnitude of the velocity.

                   aa =  half*gm1*alpha*alpha + one
                   bb = -gm1*alpha*beta
                   cc =  half*gm1*beta*beta - two*ovgm1*a2tot

                   ! Solve the equation for the magnitude of the
                   ! velocity. As this value must be positive and both aa
                   ! and bb are positive (alpha is negative and beta is
                   ! positive up till Mach = 5.0 or so, which is not
                   ! really subsonic anymore), it is clear which of the
                   ! two possible solutions must be taken. Some clipping
                   ! is present, but this is normally not active.

                   dd = bb*bb - four*aa*cc
                   dd = sqrt(max(zero,dd))
                   q  = (-bb + dd)/(two*aa)
                   q  = max(zero,q)
                   q2 = q*q

                   ! Compute the speed of sound squared from the total
                   ! speed of sound equation (== total enthalpy equation
                   ! for constant cp).

                   a2 = a2tot - half*gm1*q2

                   ! Compute the Mach number squared and cut it between
                   ! 0.0 and 1.0. Adapt the velocity and speed of sound
                   ! squared accordingly.

                   m2 = q2/a2
                   m2 = min(one,m2)
                   q2 = m2*a2
                   q  = sqrt(q2)
                   a2 = a2tot - half*gm1*q2

                   ! Compute the velocities in the halo cell and use rho,
                   ! rhoe and p as temporary buffers to store the total
                   ! temperature, total pressure and static temperature.

                   ww1(i,j,ivx)  = q*ssx
                   ww1(i,j,ivy)  = q*ssy
                   ww1(i,j,ivz)  = q*ssz

                   ww1(i,j,irho)  = ttot
                   pp1(i,j)       = ptot
                   ww1(i,j,irhoE) = a2/(gamma2(i,j)*RGas)

                   ! Compute the turbulent variables, which are
                   ! prescribed.

                   do l=nt1MG,nt2MG
                     ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
                   enddo

                   ! Set the viscosities in the halo to the viscosities
                   ! in the donor cell.

                   if( viscous )   rlv1(i,j) = rlv2(i,j)
                   if( eddyModel ) rev1(i,j) = rev2(i,j)

                 enddo
               enddo

               ! Compute the pressure and density for these halo's.

               call pRhoSubsonicInlet(icBeg(nn),icEnd(nn), &
                                      jcBeg(nn),jcEnd(nn), &
                                      kcBeg(nn),kcEnd(nn), &
                                      correctForK)

             !===========================================================

             case (massFlow)

               ! Density and velocity vector prescribed.

               ! Loop over the generic subface to set the state in the
               ! halo cells.

               do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                   ! Store a couple of variables, such as the density,
                   ! velocity and grid unit outward normal, a bit easier.

                   rho  = BCData(nn)%rho(i,j)
                   velx = BCData(nn)%velx(i,j)
                   vely = BCData(nn)%vely(i,j)
                   velz = BCData(nn)%velz(i,j)

                   nnx = BCData(nn)%norm(i,j,1)
                   nny = BCData(nn)%norm(i,j,2)
                   nnz = BCData(nn)%norm(i,j,3)

                   ! Some abbreviations in which gamma occurs.

                   gm1   = gamma2(i,j) - one
                   ovgm1 = one/gm1

                   ! Determine the acoustic Riemann variable that must be
                   ! extrapolated from the domain.

                   r    = one/ww2(i,j,irho)
                   a2   = gamma2(i,j)*pp2(i,j)*r
                   beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
                        + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

                   ! Compute the speed of sound squared in the halo.

                   a2 = half*gm1*(beta - velx*nnx - vely*nny - velz*nnz)
                   a2 = max(zero,a2)
                   a2 = a2*a2

                   ! Compute the pressure in the halo, assuming a
                   ! constant value of gamma.

                   pp1(i,j) = rho*a2/gamma2(i,j)

                   ! Simply copy the density and velocities.

                   ww1(i,j,irho) = rho
                   ww1(i,j,ivx)  = velx
                   ww1(i,j,ivy)  = vely
                   ww1(i,j,ivz)  = velz

                   ! Compute the turbulent variables, which are
                   ! prescribed.

                   do l=nt1MG,nt2MG
                     ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
                   enddo

                   ! Set the viscosities in the halo to the viscosities
                   ! in the donor cell.

                   if( viscous )   rlv1(i,j) = rlv2(i,j)
                   if( eddyModel ) rev1(i,j) = rev2(i,j)

                 enddo
               enddo

           end select

           ! Compute the total energy for these halo cells.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                            kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)

         endif inflowSubsonic
       enddo bocos

       end subroutine bcSubsonicInflow

!      ==================================================================

       subroutine pRhoSubsonicInlet(iStart,iEnd, jStart,jEnd, &
                                    kStart,kEnd, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * pRhoSubsonicInlet computes the pressure and density for the    *
!      * given range of the block to which the pointers in              *
!      * blockPointers currently point.                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cpCurveFits
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: twoThird = two*third
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iStart,iEnd, jStart,jEnd
       integer(kind=intType), intent(in) :: kStart, kEnd
       logical, intent(in)               :: correctForK
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, mm, nns, nnt
       real(kind=realType)   :: govgm1, tt, ts, pt, ratio
       real(kind=realType)   :: intTs, intTt, val
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

           ! Constant cp and thus constant gamma. Compute the coefficient
           ! gamma/(gamma-1), which occurs in the isentropic expression
           ! for the total pressure.

           govgm1 = gammaConstant/(gammaConstant - one)

           ! Loop over the given range of cells.

           do k=kStart,kEnd
             do j=jStart,jEnd
               do i=iStart,iEnd

                 ! Store the total temperature, total pressure and
                 ! static temperature a bit easier.

                 tt = w(i,j,k,irho)
                 pt = p(i,j,k)
                 ts = w(i,j,k,irhoE)

                 ! Compute the static pressure from the total pressure
                 ! and the temperature ratio. Compute the density using
                 ! the gas law.

                 ratio         = (ts/tt)**govgm1
                 p(i,j,k)      = pt*ratio
                 w(i,j,k,irho) = p(i,j,k)/(RGas*ts)

               enddo
             enddo
           enddo

!        ================================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.
           ! The ratio pt/ps is given by exp(a), where a is the integral
           ! from ts to tt of cp/(r*t).

           ! Loop over the given range of cells.

           do k=kStart,kEnd
             do j=jStart,jEnd
               do i=iStart,iEnd

                 ! Store the total temperature, total pressure and
                 ! static temperature a bit easier. Note that the
                 ! temperatures get their dimensional value.

                 tt = Tref*w(i,j,k,irho)
                 pt = p(i,j,k)
                 ts = Tref*w(i,j,k,irhoE)

                 ! Determine the integrant of cp/(r*t) for the static
                 ! temperature ts and the total temperature tt.

                 call cportIntegrant(ts, nns, intTs)
                 call cportIntegrant(tt, nnt, intTt)

                 ! Compute the value of the integral of cp/(r*t) from
                 ! ts to tt. First part is the initialization where it
                 ! is assumed that both ts and tt lie in the same
                 ! interval of the curve fits.

                 val = intTt - intTs

                 ! Correct this value if ts and tt belong to different
                 ! curve fit intervals.

                 do mm=(nns+1),nnt

                   ! The contribution from the interval mm-1. Add the
                   ! value of the integrant at the upper boundary.

                   ii = mm - 1
                   if(ii == 0_intType) then
                     val = val + (cv0+one)*log(cpTrange(0))
                   else
                     val = val + cpTempFit(ii)%intCpovrt_2
                   endif

                   ! The contribution from the interval mm. Substract
                   ! the value of integrant at the lower boundary.

                   if(mm > cpNparts) then
                     val = val - (cvn+one)*log(cpTrange(cpNparts))
                   else
                     val = val - cpTempFit(mm)%intCpovrt_1
                   endif

                 enddo

                 ! Compute the static pressure from the known
                 ! total pressure.

                 ratio    = exp(val)
                 p(i,j,k) = pt/ratio

                 ! Compute the density using the gas law.

                 ts = w(i,j,k,irhoE)
                 w(i,j,k,irho) = p(i,j,k)/(RGas*ts)

               enddo
             enddo
           enddo

       end select

       ! Add 2*rho*k/3 to the pressure if a k-equation is present.

       if( correctForK ) then
         do k=kStart,kEnd
           do j=jStart,jEnd
             do i=iStart,iEnd
               p(i,j,k) = p(i,j,k) &
                        + twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
             enddo
           enddo
         enddo
       endif

       contains

         subroutine cportIntegrant(T, nn, int)
!
!        ****************************************************************
!        *                                                              *
!        * cportIntegrant computes the integrant of the function        *
!        * cp/(r*t) for the given temperature. It also stores the       *
!        * correct curve fit interval, which is needed to determine the *
!        * entire integral in the main subroutine.                      *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(out) :: nn
         real(kind=realType),   intent(in)  :: t
         real(kind=realType),   intent(out) :: int
!
!        Local variables.
!
         integer(kind=intType) :: mm, ii, start
         real(kind=realType)   :: T2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the situation we are having here for the temperature.

         if(T <= cpTrange(0)) then

           ! Temperature is less than the smallest temperature of the
           ! curve fits. Use extrapolation using constant cp.
           ! Set nn to 0 to indicate this.

           nn  = 0
           int = (cv0+one)*log(T)

         else if(T >= cpTrange(cpNparts)) then

           ! Temperature is larger than the largest temperature of the
           ! curve fits. Use extrapolation using constant cp.
           ! Set nn to cpNparts+1 to indicate this.

           nn  = cpNparts + 1
           int = (cvn+one)*log(T)

         else

           ! Temperature is within the curve fit range. Determine
           ! the correct interval.

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

           ! Nn contains the correct curve fit interval.
           ! Compute the value of the integrant.

           int = zero
           do ii=1,cpTempFit(nn)%nterm

             mm = cpTempFit(nn)%exponents(ii)
             if(mm == 0_intType) then
               int = int + cpTempFit(nn)%constants(ii)*log(T)
             else
               T2  = T**mm
               int = int + cpTempFit(nn)%constants(ii)*T2/mm
             endif

           enddo

         endif

         end subroutine cportIntegrant

       end subroutine pRhoSubsonicInlet

!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxMatrixCoarse.f90                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidDissFluxMatrixCoarse
!
!      ******************************************************************
!      *                                                                *
!      * inviscidDissFluxMatrixCoarse computes the matrix artificial    *
!      * dissipation term. Instead of the spectral radius, as used in   *
!      * the scalar dissipation scheme, the absolute value of the flux  *
!      * jacobian is used. This routine is used on the coarser grids in *
!      * the multigrid cycle and only computes the first order          *
!      * dissipation term. It is assumed that the pointers in           *
!      * blockPointers already point to the correct block.              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use flowVarRefState
       use inputDiscretization
       use inputPhysics
       use iteration
       implicit none
!
!      Local parameters.
!
       real(kind=realType), parameter :: epsAcoustic = 0.25_realType
       real(kind=realType), parameter :: epsShear    = 0.025_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: sfil, fis0, dis0, ppor, rrad, sface
       real(kind=realType) :: gammaAvg, gm1, ovgm1, gm53, tmp, fs
       real(kind=realType) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
       real(kind=realType) :: uAvg, vAvg, wAvg, a2Avg, aAvg, hAvg
       real(kind=realType) :: alphaAvg, unAvg, ovaAvg, ova2Avg
       real(kind=realType) :: kAvg, lam1, lam2, lam3, area
       real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if rFil == 0. If so, the dissipative flux needs not to
       ! be computed.

       if(rFil == zero) return

       ! Determine whether or not the total energy must be corrected
       ! for the presence of the turbulent kinetic energy.

       if( kPresent ) then
         if((currentLevel == groundLevel) .or. turbCoupled) then
           correctForK = .true.
         else
           correctForK = .false.
         endif
       else
         correctForK = .false.
       endif

       ! Initialize sface to zero. This value will be used if the
       ! block is not moving.

       sface = zero

       ! Set a couple of constants for the scheme.

       fis0 = rFil*vis2Coarse
       sfil = one - rFil

       ! Initialize the dissipative residual to a certain times,
       ! possibly zero, the previously stored value. Owned cells
       ! only, because the halo values do not matter.

       do k=2,kl
         do j=2,jl
           do i=2,il
             fw(i,j,k,irho)  = sfil*fw(i,j,k,irho)
             fw(i,j,k,imx)   = sfil*fw(i,j,k,imx)
             fw(i,j,k,imy)   = sfil*fw(i,j,k,imy)
             fw(i,j,k,imz)   = sfil*fw(i,j,k,imz)
             fw(i,j,k,irhoE) = sfil*fw(i,j,k,irhoE)
           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the i-direction.                         *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=1,il

             ! Compute the dissipation coefficient for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = one

             dis0 = fis0*ppor

             ! Construct the vector of the first differences multiplied
             ! by dis0.

             dr  = dis0*(w(i+1,j,k,irho) - w(i,j,k,irho))
             dru = dis0*(w(i+1,j,k,irho)*w(i+1,j,k,ivx) &
                 -       w(i,  j,k,irho)*w(i,  j,k,ivx))
             drv = dis0*(w(i+1,j,k,irho)*w(i+1,j,k,ivy) &
                 -       w(i,  j,k,irho)*w(i,j,k,ivy))
             drw = dis0*(w(i+1,j,k,irho)*w(i+1,j,k,ivz) &
                 -       w(i,  j,k,irho)*w(i,j,k,ivz))
             dre = dis0*(w(i+1,j,k,irhoE) - w(i,j,k,irhoE))

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
               drk  = dis0*(w(i+1,j,k,irho)*w(i+1,j,k,itu1) &
                    -       w(i,  j,k,irho)*w(i,  j,k,itu1))
               kAvg = half*(w(i+1,j,k,itu1) + w(i,j,k,itu1))
             else
               drk  = zero
               kAvg = zero
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i+1,j,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i+1,j,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i+1,j,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i+1,j,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i+1,j,k)*p(i+1,j,k)/w(i+1,j,k,irho) &
                   +       gamma(i,  j,k)*p(i,  j,k)/w(i,  j,k,irho))

             sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
             area = sqrt(sx**2 + sy**2 + sz**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sx*tmp
             sy   = sy*tmp
             sz   = sz*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             if( addGridVelocities ) sface = sFaceI(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             lam1 = max(lam1,epsAcoustic*rrad)
             lam2 = max(lam2,epsAcoustic*rrad)
             lam3 = max(lam3,epsShear*rrad)

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = lam1*area
             lam2 = lam2*area
             lam3 = lam3*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the j-direction.                         *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=1,jl
           do i=2,il

             ! Compute the dissipation coefficient for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = one

             dis0 = fis0*ppor

             ! Construct the vector of the first differences multiplied
             ! by dis0.

             dr  = dis0*(w(i,j+1,k,irho) - w(i,j,k,irho))
             dru = dis0*(w(i,j+1,k,irho)*w(i,j+1,k,ivx) &
                 -       w(i,j,  k,irho)*w(i,j,  k,ivx))
             drv = dis0*(w(i,j+1,k,irho)*w(i,j+1,k,ivy) &
                 -       w(i,j,  k,irho)*w(i,j,  k,ivy))
             drw = dis0*(w(i,j+1,k,irho)*w(i,j+1,k,ivz) &
                 -       w(i,j,  k,irho)*w(i,j,  k,ivz))
             dre = dis0*(w(i,j+1,k,irhoE) - w(i,j,k,irhoE))

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
               drk  = dis0*(w(i,j+1,k,irho)*w(i,j+1,k,itu1) &
                    -       w(i,j,  k,irho)*w(i,j,  k,itu1))
               kAvg = half*(w(i,j+1,k,itu1) + w(i,j,k,itu1))
             else
               drk  = zero
               kAvg = zero
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j+1,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j+1,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j+1,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j+1,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j+1,k)*p(i,j+1,k)/w(i,j+1,k,irho) &
                    +      gamma(i,j,  k)*p(i,j,  k)/w(i,j,  k,irho))

             sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
             area = sqrt(sx**2 + sy**2 + sz**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sx*tmp
             sy   = sy*tmp
             sz   = sz*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             if( addGridVelocities ) sface = sFaceJ(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             lam1 = max(lam1,epsAcoustic*rrad)
             lam2 = max(lam2,epsAcoustic*rrad)
             lam3 = max(lam3,epsShear*rrad)

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = lam1*area
             lam2 = lam2*area
             lam3 = lam3*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the k-direction.                         *
!      *                                                                *
!      ******************************************************************
!
       do k=1,kl
         do j=2,jl
           do i=2,il

             ! Compute the dissipation coefficient for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = one

             dis0 = fis0*ppor

             ! Construct the vector of the first differences multiplied
             ! by dis0.

             dr  = dis0*(w(i,j,k+1,irho) - w(i,j,k,irho))
             dru = dis0*(w(i,j,k+1,irho)*w(i,j,k+1,ivx) &
                 -       w(i,j,k,  irho)*w(i,j,k,  ivx))
             drv = dis0*(w(i,j,k+1,irho)*w(i,j,k+1,ivy) &
                 -       w(i,j,k,  irho)*w(i,j,k,  ivy))
             drw = dis0*(w(i,j,k+1,irho)*w(i,j,k+1,ivz) &
                 -       w(i,j,k,  irho)*w(i,j,k,  ivz))
             dre = dis0*(w(i,j,k+1,irhoE) - w(i,j,k,irhoE))

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
               drk  = dis0*(w(i,j,k+1,irho)*w(i,j,k+1,itu1) &
                    -       w(i,j,k,  irho)*w(i,j,k,  itu1))
               kAvg = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
             else
               drk  = zero
               kAvg = zero
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j,k+1) + gamma(i,j,k))
             gm1       = gammaAvg - one
             ovgm1     = one/gm1
             gm53      = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j,k+1,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j,k+1,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j,k+1,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j,k+1)*p(i,j,k+1)/w(i,j,k+1,irho) &
                   +       gamma(i,j,k)  *p(i,j,k)  /w(i,j,k,  irho))

             sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
             area = sqrt(sx**2 + sy**2 + sz**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sx*tmp
             sy   = sy*tmp
             sz   = sz*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             if( addGridVelocities ) sface = sFaceK(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             lam1 = max(lam1,epsAcoustic*rrad)
             lam2 = max(lam2,epsAcoustic*rrad)
             lam3 = max(lam3,epsShear*rrad)

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = lam1*area
             lam2 = lam2*area
             lam3 = lam3*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

           enddo
         enddo
       enddo

       end subroutine inviscidDissFluxMatrixCoarse

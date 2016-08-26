!
!       File:          inviscidDissFluxMatrix.f90                      
!       Author:        Edwin van der Weide                             
!       Starting date: 03-25-2003                                      
!       Last modified: 10-29-2007                                      
!
subroutine inviscidDissFluxMatrix
  !
  !       inviscidDissFluxMatrix computes the matrix artificial          
  !       dissipation term. Instead of the spectral radius, as used in   
  !       the scalar dissipation scheme, the absolute value of the flux  
  !       jacobian is used. This leads to a less diffusive and           
  !       consequently more accurate scheme. It is assumed that the      
  !       pointers in blockPointers already point to the correct block.  
  !
  use constants
  use blockPointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb, &
       w, p, porI, porJ, porK, fw, gamma, si, sj, sk, &
       indFamilyI, indFamilyJ, indFamilyK, spectralSol, addGridVelocities, &
       sFaceI, sfaceJ, sFacek, factFamilyI, factFamilyJ, factFamilyK
  use flowVarRefState, only : pInfCorr
  use inputDiscretization, only: vis2, vis4
  use inputPhysics, only : equations
  use iteration, only : rFil
  use cgnsGrid, only: massFlowFamilyDiss
  use utils, only : getCorrectForK, myDim
  implicit none
  !
  !      Local parameters.
  !
  real(kind=realType), parameter :: dpMax        = 0.25_realType
  real(kind=realType), parameter :: epsAcoustic  = 0.25_realType
  real(kind=realType), parameter :: epsShear     = 0.025_realType
  real(kind=realType), parameter :: omega        = 0.5_realType
  real(kind=realType), parameter :: oneMinOmega = one - omega
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, ind, ii

  real(kind=realType) :: plim, sface
  real(kind=realType) :: sfil, fis2, fis4
  real(kind=realType) :: gammaAvg, gm1, ovgm1, gm53
  real(kind=realType) :: ppor, rrad, dis2, dis4
  real(kind=realType) :: dp1, dp2, tmp, fs
  real(kind=realType) :: ddw1, ddw2, ddw3, ddw4, ddw5, ddw6
  real(kind=realType) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
  real(kind=realType) :: uAvg, vAvg, wAvg, a2Avg, aAvg, hAvg
  real(kind=realType) :: alphaAvg, unAvg, ovaAvg, ova2Avg
  real(kind=realType) :: kAvg, lam1, lam2, lam3, area
  real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
  real(kind=realType),dimension(1:ie,1:je,1:ke,3) :: dss
  logical :: correctForK

  ! Check if rFil == 0. If so, the dissipative flux needs not to
  ! be computed.

  if(abs(rFil) < thresholdReal) return

  ! Set the value of plim. To be fully consistent this must have
  ! the dimension of a pressure. Therefore a fraction of pInfCorr
  ! is used.

  plim = 0.001_realType*pInfCorr

  ! Determine whether or not the total energy must be corrected
  ! for the presence of the turbulent kinetic energy.

  correctForK = getCorrectForK()

  ! Initialize sface to zero. This value will be used if the
  ! block is not moving.

  sface = zero

  ! Set a couple of constants for the scheme.

  fis2 = rFil*vis2
  fis4 = rFil*vis4
  sfil = one - rFil

  ! Initialize the dissipative residual to a certain times,
  ! possibly zero, the previously stored value. 

  fw = sfil*fw

  ! Compute the pressure sensor for each cell, in each direction:
#ifdef TAPENADE_REVERSE
  !$AD II-LOOP
  do ii=0,ie*je*ke-1
     i = mod(ii, ie) + 1
     j = mod(ii/ie, je) + 1
     k = ii/(ie*je) + 1
#else
     do k=1,ke
        do j=1,je
           do i=1,ie
#endif             
              dss(i,j,k,1) =abs((p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k))        &
                   /     (omega*(p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k)) &
                   +      oneMinOmega*(abs(p(i+1,j,k) - p(i,j,k))      &
                   +                   abs(p(i,j,k) - p(i-1,j,k))) + plim))


              dss(i,j,k,2) =abs((p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k))        &
                   /     (omega*(p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k)) &
                   +      oneMinOmega*(abs(p(i,j+1,k) - p(i,j,k))      &
                   +                   abs(p(i,j,k) - p(i,j-1,k))) + plim))

              dss(i,j,k,3) =  abs((p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1))        &
                   /     (omega*(p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1)) &
                   +      oneMinOmega*(abs(p(i,j,k+1) - p(i,j,k))      &
                   +                   abs(p(i,j,k) - p(i,j,k-1))) + plim))
#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif
  !
  !       Dissipative fluxes in the i-direction.                         
  !
#ifdef TAPENADE_REVERSE
  !$AD II-LOOP
  do ii=0,il*ny*nz-1
     i = mod(ii, il) + 1
     j = mod(ii/il, ny) + 2
     k = ii/(il*ny) + 2
#else
     do k=2,kl
        do j=2,jl
           do i=1,il
#endif 
              ! Compute the dissipation coefficients for this face.

              ppor = zero
              if(porI(i,j,k) == normalFlux) ppor = one
              dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
              dis4 = myDim(ppor*fis4, dis2)

              ! Construct the vector of the first and third differences
              ! multiplied by the appropriate constants.

              ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
              dr  = dis2*ddw1 &
                   - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

              ddw2 = w(i+1,j,k,irho)*w(i+1,j,k,ivx) &
                   - w(i,j,k,irho)*w(i,j,k,ivx)
              dru = dis2*ddw2                             &
                   - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivx) &
                   -       w(i-1,j,k,irho)*w(i-1,j,k,ivx) - three*ddw2)

              ddw3 = w(i+1,j,k,irho)*w(i+1,j,k,ivy) &
                   - w(i,j,k,irho)*w(i,j,k,ivy)
              drv = dis2*ddw3                             &
                   - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivy) &
                   -       w(i-1,j,k,irho)*w(i-1,j,k,ivy) - three*ddw3)

              ddw4 = w(i+1,j,k,irho)*w(i+1,j,k,ivz) &
                   - w(i,j,k,irho)*w(i,j,k,ivz)
              drw = dis2*ddw4                             &
                   - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivz) &
                   -       w(i-1,j,k,irho)*w(i-1,j,k,ivz) - three*ddw4)

              ddw5 = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
              dre = dis2*ddw5 &
                   - dis4*(w(i+2,j,k,irhoE) - w(i-1,j,k,irhoE) - three*ddw5)

              ! In case a k-equation is present, compute the difference
              ! of rhok and store the average value of k. If not present,
              ! set both these values to zero, such that later on no
              ! decision needs to be made anymore.

              if( correctForK ) then
                 ddw6 = w(i+1,j,k,irho)*w(i+1,j,k,itu1) &
                      - w(i,j,k,irho)*w(i,j,k,itu1)
                 drk = dis2*ddw6                              &
                      - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,itu1) &
                      -       w(i-1,j,k,irho)*w(i-1,j,k,itu1) - three*ddw6)

                 kAvg = half*(w(i,j,k,itu1) + w(i+1,j,k,itu1))
              else
                 drk   = zero
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

              area = sqrt(si(i,j,k,1)**2 + si(i,j,k,2)**2 + si(i,j,k,3)**2)
              tmp  = one/max(1.e-25_realType,area)
              sx   = si(i,j,k,1)*tmp
              sy   = si(i,j,k,2)*tmp
              sz   = si(i,j,k,3)*tmp

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

              ! Multiply the eigenvalues by the area to obtain
              ! the correct values for the dissipation term.

              lam1 = max(lam1,epsAcoustic*rrad)*area
              lam2 = max(lam2,epsAcoustic*rrad)*area
              lam3 = max(lam3,epsShear*rrad)*area

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
#ifndef USE_TAPENADE
              ind = indFamilyI(i,j,k)
              massFlowFamilyDiss(ind,spectralSol) =       &
                   massFlowFamilyDiss(ind,spectralSol) &
                   - factFamilyI(i,j,k)*fs
#endif
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

#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif
  !
  !       Dissipative fluxes in the j-direction.                         
  !
#ifdef TAPENADE_REVERSE
  !$AD II-LOOP
  do ii=0,nx*jl*nz-1
     i = mod(ii, nx) + 2
     j = mod(ii/nx, jl) + 1
     k = ii/(nx*jl) + 2
#else
     do k=2,kl
        do j=1,jl
           do i=2,il
#endif

              ! Compute the dissipation coefficients for this face.

              ppor = zero
              if(porJ(i,j,k) == normalFlux) ppor = one

              dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,2), dss(i,j+1,k,2)))
              dis4 = myDim(ppor*fis4, dis2)

              ! Construct the vector of the first and third differences
              ! multiplied by the appropriate constants.

              ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
              dr  = dis2*ddw1 &
                   - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

              ddw2 = w(i,j+1,k,irho)*w(i,j+1,k,ivx) &
                   - w(i,j,k,irho)*w(i,j,k,ivx)
              dru = dis2*ddw2                             &
                   - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivx) &
                   -       w(i,j-1,k,irho)*w(i,j-1,k,ivx) - three*ddw2)

              ddw3 = w(i,j+1,k,irho)*w(i,j+1,k,ivy) &
                   - w(i,j,k,irho)*w(i,j,k,ivy)
              drv = dis2*ddw3 &
                   - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivy) &
                   -       w(i,j-1,k,irho)*w(i,j-1,k,ivy) - three*ddw3)

              ddw4 = w(i,j+1,k,irho)*w(i,j+1,k,ivz) &
                   - w(i,j,k,irho)*w(i,j,k,ivz)
              drw = dis2*ddw4 &
                   - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivz) &
                   -       w(i,j-1,k,irho)*w(i,j-1,k,ivz) - three*ddw4)

              ddw5 = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
              dre = dis2*ddw5 &
                   - dis4*(w(i,j+2,k,irhoE) - w(i,j-1,k,irhoE) - three*ddw5)

              ! In case a k-equation is present, compute the difference
              ! of rhok and store the average value of k. If not present,
              ! set both these values to zero, such that later on no
              ! decision needs to be made anymore.

              if( correctForK ) then
                 ddw6 = w(i,j+1,k,irho)*w(i,j+1,k,itu1) &
                      - w(i,j,k,irho)*w(i,j,k,itu1)
                 drk = dis2*ddw6                              &
                      - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,itu1) &
                      -       w(i,j-1,k,irho)*w(i,j-1,k,itu1) - three*ddw6)

                 kAvg = half*(w(i,j,k,itu1) + w(i,j+1,k,itu1))
              else
                 drk   = zero
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
                   +       gamma(i,j,  k)*p(i,j,  k)/w(i,j,  k,irho))

              area = sqrt(sj(i,j,k,1)**2 + sj(i,j,k,2)**2 + sj(i,j,k,3)**2)
              tmp  = one/max(1.e-25_realType,area)
              sx   = sj(i,j,k,1)*tmp
              sy   = sj(i,j,k,2)*tmp
              sz   = sj(i,j,k,3)*tmp

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

              ! Multiply the eigenvalues by the area to obtain
              ! the correct values for the dissipation term.

              lam1 = max(lam1,epsAcoustic*rrad)*area
              lam2 = max(lam2,epsAcoustic*rrad)*area
              lam3 = max(lam3,epsShear*rrad)*area

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
#ifndef USE_TAPENADE
              ind = indFamilyJ(i,j,k)
              massFlowFamilyDiss(ind,spectralSol) =       &
                   massFlowFamilyDiss(ind,spectralSol) &
                   - factFamilyJ(i,j,k)*fs
#endif
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

#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif
  !
  !       Dissipative fluxes in the k-direction.                         
  !
#ifdef TAPENADE_REVERSE
  !$AD II-LOOP
  do ii=0,nx*ny*kl-1
     i = mod(ii, nx) + 2
     j = mod(ii/nx, ny) + 2
     k = ii/(nx*ny) + 1
#else
     do k=1,kl
        do j=2,jl
           do i=2,il
#endif     
              ! Compute the dissipation coefficients for this face.

              ppor = zero
              if(porK(i,j,k) == normalFlux) ppor = one

              dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
              dis4 = myDim(ppor*fis4, dis2)

              ! Construct the vector of the first and third differences
              ! multiplied by the appropriate constants.

              ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
              dr  = dis2*ddw1 &
                   - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

              ddw2 = w(i,j,k+1,irho)*w(i,j,k+1,ivx) &
                   - w(i,j,k,irho)*w(i,j,k,ivx)
              dru = dis2*ddw2                             &
                   - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivx) &
                   -       w(i,j,k-1,irho)*w(i,j,k-1,ivx) - three*ddw2)

              ddw3 = w(i,j,k+1,irho)*w(i,j,k+1,ivy) &
                   - w(i,j,k,irho)*w(i,j,k,ivy)
              drv = dis2*ddw3 &
                   - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivy) &
                   -       w(i,j,k-1,irho)*w(i,j,k-1,ivy) - three*ddw3)

              ddw4 = w(i,j,k+1,irho)*w(i,j,k+1,ivz) &
                   - w(i,j,k,irho)*w(i,j,k,ivz)
              drw = dis2*ddw4 &
                   - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivz) &
                   -       w(i,j,k-1,irho)*w(i,j,k-1,ivz) - three*ddw4)

              ddw5 = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
              dre = dis2*ddw5 &
                   - dis4*(w(i,j,k+2,irhoE) - w(i,j,k-1,irhoE) - three*ddw5)

              ! In case a k-equation is present, compute the difference
              ! of rhok and store the average value of k. If not present,
              ! set both these values to zero, such that later on no
              ! decision needs to be made anymore.

              if( correctForK ) then
                 ddw6 = w(i,j,k+1,irho)*w(i,j,k+1,itu1) &
                      - w(i,j,k,irho)*w(i,j,k,itu1)
                 drk = dis2*ddw6                              &
                      - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,itu1) &
                      -       w(i,j,k-1,irho)*w(i,j,k-1,itu1) - three*ddw6)

                 kAvg = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
              else
                 drk   = zero
                 kAvg = zero
              endif

              ! Compute the average value of gamma and compute some
              ! expressions in which it occurs.

              gammaAvg = half*(gamma(i,j,k+1) + gamma(i,j,k))
              gm1      = gammaAvg - one
              ovgm1    = one/gm1
              gm53     = gammaAvg - five*third

              ! Compute the average state at the interface.

              uAvg  = half*(w(i,j,k+1,ivx) + w(i,j,k,ivx))
              vAvg  = half*(w(i,j,k+1,ivy) + w(i,j,k,ivy))
              wAvg  = half*(w(i,j,k+1,ivz) + w(i,j,k,ivz))
              a2Avg = half*(gamma(i,j,k+1)*p(i,j,k+1)/w(i,j,k+1,irho) &
                   +       gamma(i,j,k)  *p(i,j,k)  /w(i,j,k,  irho))

              area = sqrt(sk(i,j,k,1)**2 + sk(i,j,k,2)**2 + sk(i,j,k,3)**2)
              tmp  = one/max(1.e-25_realType,area)
              sx   = sk(i,j,k,1)*tmp
              sy   = sk(i,j,k,2)*tmp
              sz   = sk(i,j,k,3)*tmp

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

              ! Multiply the eigenvalues by the area to obtain
              ! the correct values for the dissipation term.

              lam1 = max(lam1,epsAcoustic*rrad)*area
              lam2 = max(lam2,epsAcoustic*rrad)*area
              lam3 = max(lam3,epsShear*rrad)*area

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
#ifndef USE_TAPENADE
              ind = indFamilyK(i,j,k)
              massFlowFamilyDiss(ind,spectralSol) =       &
                   massFlowFamilyDiss(ind,spectralSol) &
                   - factFamilyK(i,j,k)*fs
#endif
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

#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif
end subroutine inviscidDissFluxMatrix

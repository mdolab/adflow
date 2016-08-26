subroutine inviscidDissFluxScalar
  !
  !       inviscidDissFluxScalar computes the scalar artificial          
  !       dissipation, see AIAA paper 81-1259, for a given block.        
  !       Therefore it is assumed that the pointers in  blockPointers    
  !       already point to the correct block.                            
  !
  use constants
  use blockPointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb, &
       w, p, porI, porJ, porK, fw, radI, radJ, radK, gamma
  use flowVarRefState, only : gammaInf, pInfCorr, rhoInf
  use inputDiscretization, only: vis2, vis4
  use inputPhysics, only : equations
  use iteration, only : rFil
  use utils, only : myDim
  implicit none
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: dssMax = 0.25_realType
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, ind, ii

  real(kind=realType) :: sslim, rhoi
  real(kind=realType) :: sfil, fis2, fis4
  real(kind=realType) :: ppor, rrad, dis2, dis4
  real(kind=realType) :: ddw1,ddw2,ddw3,ddw4,ddw5,fs
  real(kind=realType),dimension(1:ie,1:je,1:ke,3) :: dss
  real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ss

  ! Check if rFil == 0. If so, the dissipative flux needs not to
  ! be computed.

  if(abs(rFil) < thresholdReal) return

  ! Determine the variables used to compute the switch.
  ! For the inviscid case this is the pressure; for the viscous
  ! case it is the entropy.

  select case (equations)
  case (EulerEquations)

     ! Inviscid case. Pressure switch is based on the pressure.
     ! Also set the value of sslim. To be fully consistent this
     ! must have the dimension of pressure and it is therefore
     ! set to a fraction of the free stream value.

     sslim = 0.001_realType*pInfCorr

     ! Copy the pressure in ss. Only need the entries used in the
     ! discretization, i.e. not including the corner halo's, but we'll
     ! just copy all anyway. 

     ss = P
     !===============================================================

  case (NSEquations, RANSEquations)

     ! Viscous case. Pressure switch is based on the entropy.
     ! Also set the value of sslim. To be fully consistent this
     ! must have the dimension of entropy and it is therefore
     ! set to a fraction of the free stream value.

     sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

     ! Store the entropy in ss. See above. 

#ifdef TAPENADE_REVERSE
     !$AD II-LOOP
     do ii=0,(ib+1)*(jb+1)*(kb+1)-1
        i = mod(ii, ib+1)
        j = mod(ii/(ib+1), jb+1) 
        k = ii/((ib+1)*(jb+1))
#else
        do k=0,kb
           do j=0,jb
              do i=0,ib
#endif      
                 ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
#ifdef TAPENADE_REVERSE
              end do
#else
           end do
        end do
     end do
#endif
  end select

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
              dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
                   /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

              dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
                   /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

              dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
                   /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif

  ! Set a couple of constants for the scheme.

  fis2 = rFil*vis2
  fis4 = rFil*vis4
  sfil = one - rFil

  ! Initialize the dissipative residual to a certain times,
  ! possibly zero, the previously stored value. Owned cells
  ! only, because the halo values do not matter.

  fw = sfil*fw
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
              if(porI(i,j,k) == normalFlux) ppor = half
              rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

              dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
              dis4 = myDim(fis4*rrad, dis2)

              ! Compute and scatter the dissipative flux.
              ! Density. Store it in the mass flow of the
              ! appropriate sliding mesh interface.

              ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
              fs  = dis2*ddw1 &
                   - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

              fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
              fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

              ! X-momentum.

              ddw2 = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
              fs  = dis2*ddw2 &
                   - dis4*(w(i+2,j,k,ivx)*w(i+2,j,k,irho) - w(i-1,j,k,ivx)*w(i-1,j,k,irho) - three*ddw2)

              fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
              fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

              ! Y-momentum.

              ddw3 = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
              fs  = dis2*ddw3 &
                   - dis4*(w(i+2,j,k,ivy)*w(i+2,j,k,irho) - w(i-1,j,k,ivy)*w(i-1,j,k,irho) - three*ddw3)

              fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
              fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

              ! Z-momentum.

              ddw4 = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
              fs  = dis2*ddw4 &
                   - dis4*(w(i+2,j,k,ivz)*w(i+2,j,k,irho) - w(i-1,j,k,ivz)*w(i-1,j,k,irho) - three*ddw4)

              fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
              fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

              ! Energy.

              ddw5 = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
              fs  = dis2*ddw5 &
                   - dis4*((w(i+2,j,k,irhoE) + P(i+2,j,k)) - (w(i-1,j,k,irhoE) + P(i-1,j,k)) - three*ddw5)

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
              if(porJ(i,j,k) == normalFlux) ppor = half
              rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

              dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2)))
              dis4 = myDim(fis4*rrad, dis2)

              ! Compute and scatter the dissipative flux.
              ! Density. Store it in the mass flow of the
              ! appropriate sliding mesh interface.

              ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
              fs  = dis2*ddw1 &
                   - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

              fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
              fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

              ! X-momentum.

              ddw2 = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
              fs  = dis2*ddw2 &
                   - dis4*(w(i,j+2,k,ivx)*w(i,j+2,k,irho) - w(i,j-1,k,ivx)*w(i,j-1,k,irho) - three*ddw2)

              fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
              fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

              ! Y-momentum.

              ddw3 = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
              fs  = dis2*ddw3 &
                   - dis4*(w(i,j+2,k,ivy)*w(i,j+2,k,irho) - w(i,j-1,k,ivy)*w(i,j-1,k,irho) - three*ddw3)

              fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
              fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

              ! Z-momentum.

              ddw4 = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
              fs  = dis2*ddw4 &
                   - dis4*(w(i,j+2,k,ivz)*w(i,j+2,k,irho) - w(i,j-1,k,ivz)*w(i,j-1,k,irho) - three*ddw4)

              fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
              fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

              ! Energy.

              ddw5 = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
              fs  = dis2*ddw5 &
                   - dis4*((w(i,j+2,k,irhoE) + P(i,j+2,k)) - (w(i,j-1,k,irhoE) + P(i,j-1,k)) - three*ddw5)

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
              if(porK(i,j,k) == normalFlux) ppor = half
              rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

              dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
              dis4 = myDim(fis4*rrad, dis2)

              ! Compute and scatter the dissipative flux.
              ! Density. Store it in the mass flow of the
              ! appropriate sliding mesh interface.

              ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
              fs  = dis2*ddw1 &
                   - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

              fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
              fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

              ! X-momentum.

              ddw2 = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
              fs  = dis2*ddw2 &
                   - dis4*(w(i,j,k+2,ivx)*w(i,j,k+2,irho) - w(i,j,k-1,ivx)*w(i,j,k-1,irho) - three*ddw2)

              fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
              fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

              ! Y-momentum.

              ddw3 = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
              fs  = dis2*ddw3 &
                   - dis4*(w(i,j,k+2,ivy)*w(i,j,k+2,irho) - w(i,j,k-1,ivy)*w(i,j,k-1,irho) - three*ddw3)

              fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
              fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

              ! Z-momentum.

              ddw4 = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
              fs  = dis2*ddw4 &
                   - dis4*(w(i,j,k+2,ivz)*w(i,j,k+2,irho) - w(i,j,k-1,ivz)*w(i,j,k-1,irho) - three*ddw4)

              fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
              fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

              ! Energy.

              ddw5 = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
              fs  = dis2*ddw5 &
                   - dis4*((w(i,j,k+2,irhoE) + P(i,j,k+2)) - (w(i,j,k-1,irhoE) + P(i,j,k-1)) - three*ddw5)

              fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
              fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

#ifdef TAPENADE_REVERSE
           end do
#else
        end do
     end do
  end do
#endif
end subroutine inviscidDissFluxScalar

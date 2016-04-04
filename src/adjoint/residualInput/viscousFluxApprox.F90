subroutine viscousFluxApprox

  use blockPointers
  use flowVarRefState
  use inputPhysics
  use iteration
  implicit none
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: twoThird = two*third
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k
  integer(kind=intType) :: ii, jj, kk

  real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
  real(kind=realType) :: gm1, factLamHeat, factTurbHeat
  real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
  real(kind=realType) :: q_x, q_y, q_z, ubar, vbar, wbar
  real(kind=realType) :: corr, ssx, ssy, ssz, ss, fracDiv
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz
  real(kind=realType) :: fmx, fmy, fmz, frhoE
  real(kind=realType) :: dd
  logical :: correctForK

  mue = zero
  rFilv = rFil

  ! Viscous fluxes in the I-direction

  do k=2,kl
     do j=2,jl
        do i=1,il

           ! Compute the vector from the center of cell i to cell i+1           
           ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
                +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
                +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
                +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
           ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
                +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
                +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
                +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
           ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
                +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
                +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
                +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

           ! And determine one/ length of vector squared
           ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
           ssx = ss*ssx
           ssy = ss*ssy
           ssz = ss*ssz

           ! Now compute each gradient
           dd = w(i+1, j, k, ivx)-w(i, j, k, ivx)
           u_x = dd*ssx
           u_y = dd*ssy
           u_z = dd*ssz

           dd = w(i+1, j, k, ivy)-w(i, j, k, ivy)
           v_x = dd*ssx
           v_y = dd*ssy
           v_z = dd*ssz

           dd = w(i+1, j, k, ivz)-w(i, j, k, ivz)
           w_x = dd*ssx
           w_y = dd*ssy
           w_z = dd*ssz

           dd = aa(i+1, j, k)-aa(i, j, k)
           q_x = -dd*ssx
           q_y = -dd*ssy
           q_z = -dd*ssz

           por = half*rFilv
           if(porI(i,j,k) == noFlux) por = zero

           ! Compute the laminar and (if present) the eddy viscosities
           ! multiplied by the porosity. Compute the factor in front of
           ! the gradients of the speed of sound squared for the heat
           ! flux.

           mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
           if( eddyModel ) mue = por*(rev(i,j,k) + rev(i+1,j,k))
           mut = mul + mue

           gm1          = half*(gamma(i,j,k) +gamma(i+1,j,k))- one
           factLamHeat  = one/(prandtl*gm1)
           factTurbHeat = one/(prandtlTurb*gm1)

           heatCoef = mul*factLamHeat + mue*factTurbHeat

           ! Compute the stress tensor and the heat flux vector.

           fracDiv = twoThird*(u_x + v_y + w_z)

           tauxx = mut*(two*u_x - fracDiv)
           tauyy = mut*(two*v_y - fracDiv)
           tauzz = mut*(two*w_z - fracDiv)

           tauxy = mut*(u_y + v_x)
           tauxz = mut*(u_z + w_x)
           tauyz = mut*(v_z + w_y)

           q_x = heatCoef*q_x
           q_y = heatCoef*q_y
           q_z = heatCoef*q_z

           ! Compute the average velocities for the face. Remember that
           ! the velocities are stored and not the momentum.

           ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
           vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
           wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

           ! Compute the viscous fluxes for this i-face.

           fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) + tauxz*si(i,j,k,3)
           fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) + tauyz*si(i,j,k,3)
           fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) + tauzz*si(i,j,k,3)
           frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
                + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
                + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
                - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

           ! Update the residuals of cell i and i+1.

           fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
           fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
           fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
           fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

           fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
           fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
           fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
           fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

        end do
     end do
  end do

  ! Viscous fluxes in the J-direction

  do k=2,kl
     do j=1,jl
        do i=2,il

           ! Compute the vector from the center of cell j to cell j+1           
           ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
                +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
                +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
                +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
           ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
                +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
                +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
                +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
           ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
                +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
                +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
                +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

           ! And determine one/ length of vector squared
           ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
           ssx = ss*ssx
           ssy = ss*ssy
           ssz = ss*ssz

           ! Now compute each gradient
           dd = w(i, j+1, k, ivx)-w(i, j, k, ivx)
           u_x = dd*ssx
           u_y = dd*ssy
           u_z = dd*ssz

           dd = w(i, j+1, k, ivy)-w(i, j, k, ivy)
           v_x = dd*ssx
           v_y = dd*ssy
           v_z = dd*ssz

           dd = w(i, j+1, k, ivz)-w(i, j, k, ivz)
           w_x = dd*ssx
           w_y = dd*ssy
           w_z = dd*ssz

           dd = aa(i, j+1, k)-aa(i, j, k)
           q_x = -dd*ssx
           q_y = -dd*ssy
           q_z = -dd*ssz

           por = half*rFilv
           if(porJ(i,j,k) == noFlux) por = zero

           ! Compute the laminar and (if present) the eddy viscosities
           ! multiplied by the porosity. Compute the factor in front of
           ! the gradients of the speed of sound squared for the heat
           ! flux.

           mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
           if( eddyModel ) mue = por*(rev(i,j,k) + rev(i,j+1,k))
           mut = mul + mue

           gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
           factLamHeat  = one/(prandtl*gm1)
           factTurbHeat = one/(prandtlTurb*gm1)

           heatCoef = mul*factLamHeat + mue*factTurbHeat

           ! Compute the stress tensor and the heat flux vector.

           fracDiv = twoThird*(u_x + v_y + w_z)

           tauxx = mut*(two*u_x - fracDiv)
           tauyy = mut*(two*v_y - fracDiv)
           tauzz = mut*(two*w_z - fracDiv)

           tauxy = mut*(u_y + v_x)
           tauxz = mut*(u_z + w_x)
           tauyz = mut*(v_z + w_y)

           q_x = heatCoef*q_x
           q_y = heatCoef*q_y
           q_z = heatCoef*q_z

           ! Compute the average velocities for the face. Remember that
           ! the velocities are stored and not the momentum.

           ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
           vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
           wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

           ! Compute the viscous fluxes for this j-face.

           fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) + tauxz*sj(i,j,k,3)
           fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) + tauyz*sj(i,j,k,3)
           fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) + tauzz*sj(i,j,k,3)
           frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
                + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
                + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
                - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

           ! Update the residuals of cell j and j+1.

           fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
           fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
           fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
           fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

           fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
           fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
           fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
           fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE

        end do
     end do
  end do

  ! Viscous fluxes in the K-direction

  do k=1,kl
     do j=2,jl
        do i=2,il

           ! Compute the vector from the center of cell k to cell k+1           
           ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
                +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
                +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
                +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
           ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
                +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
                +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
                +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
           ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
                +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
                +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
                +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))
           ! And determine one/ length of vector squared
           ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
           ssx = ss*ssx
           ssy = ss*ssy
           ssz = ss*ssz

           ! Now compute each gradient
           dd = w(i, j, k+1, ivx)-w(i, j, k, ivx)
           u_x = dd*ssx
           u_y = dd*ssy
           u_z = dd*ssz

           dd = w(i, j, k+1, ivy)-w(i, j, k, ivy)
           v_x = dd*ssx
           v_y = dd*ssy
           v_z = dd*ssz

           dd = w(i, j, k+1, ivz)-w(i, j, k, ivz)
           w_x = dd*ssx
           w_y = dd*ssy
           w_z = dd*ssz

           dd = aa(i, j, k+1)-aa(i, j, k)
           q_x = -dd*ssx
           q_y = -dd*ssy
           q_z = -dd*ssz

           por = half*rFilv
           if(porK(i,j,k) == noFlux) por = zero

           ! Compute the laminar and (if present) the eddy viscosities
           ! multiplied by the porosity. Compute the factor in front of
           ! the gradients of the speed of sound squared for the heat
           ! flux.

           mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
           if( eddyModel ) mue = por*(rev(i,j,k) + rev(i,j,k+1))
           mut = mul + mue

           gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
           factLamHeat  = one/(prandtl*gm1)
           factTurbHeat = one/(prandtlTurb*gm1)

           heatCoef = mul*factLamHeat + mue*factTurbHeat

           ! Compute the stress tensor and the heat flux vector.

           fracDiv = twoThird*(u_x + v_y + w_z)

           tauxx = mut*(two*u_x - fracDiv)
           tauyy = mut*(two*v_y - fracDiv)
           tauzz = mut*(two*w_z - fracDiv)

           tauxy = mut*(u_y + v_x)
           tauxz = mut*(u_z + w_x)
           tauyz = mut*(v_z + w_y)

           q_x = heatCoef*q_x
           q_y = heatCoef*q_y
           q_z = heatCoef*q_z

           ! Compute the average velocities for the face. Remember that
           ! the velocities are stored and not the momentum.

           ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
           vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
           wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

           ! Compute the viscous fluxes for this j-face.

           fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) + tauxz*sk(i,j,k,3)
           fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) + tauyz*sk(i,j,k,3)
           fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) + tauzz*sk(i,j,k,3)
           frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1) &
                + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2) &
                + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3) &
                - q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

           ! Update the residuals of cell j and j+1.

           fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
           fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
           fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
           fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

           fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
           fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
           fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
           fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

        end do
     end do
  end do

end subroutine viscousFluxApprox

subroutine viscousFlux
  !
  !       viscousFlux computes the viscous fluxes using a central        
  !       difference scheme for a block.                                 
  !       It is assumed that the pointers in block pointer already point 
  !       to the correct block.                                          
  !
  use constants
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
  integer(kind=intType) :: i, j, k, ii

  real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
  real(kind=realType) :: gm1, factLamHeat, factTurbHeat
  real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
  real(kind=realType) :: q_x, q_y, q_z, ubar, vbar, wbar
  real(kind=realType) :: corr, ssx, ssy, ssz, ss, fracDiv
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz
  real(kind=realType) :: fmx, fmy, fmz, frhoE
  logical :: correctForK, storeWallTensor


  ! Set rFilv to rFil to indicate that this is the viscous part.
  ! If rFilv == 0 the viscous residuals need not to be computed
  ! and a return can be made.

  rFilv = rFil

  if(abs(rFilv) < thresholdReal) return

  ! Determine whether or not the wall stress tensor and wall heat
  ! flux must be stored for viscous walls.

  storeWallTensor = .false.
  if( wallFunctions ) then
     storeWallTensor = .true.
  else if(rkStage == 0 .and. currentLevel == groundLevel) then
     storeWallTensor = .true.
  endif

  !
  !         viscous fluxes in the k-direction.                           
  !
  continue
  !$AD CHECKPOINT-START
  mue = zero
#ifdef TAPENADE_FAST
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

              ! Set the value of the porosity. If not zero, it is set
              ! to average the eddy-viscosity and to take the factor
              ! rFilv into account.

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

              ! Compute the gradients at the face by averaging the four
              ! nodal values.

              u_x = fourth*(ux(i-1,j-1,k) + ux(i,j-1,k) &
                   +         ux(i-1,j,  k) + ux(i,j,  k))
              u_y = fourth*(uy(i-1,j-1,k) + uy(i,j-1,k) &
                   +         uy(i-1,j,  k) + uy(i,j,  k))
              u_z = fourth*(uz(i-1,j-1,k) + uz(i,j-1,k) &
                   +         uz(i-1,j,  k) + uz(i,j,  k))

              v_x = fourth*(vx(i-1,j-1,k) + vx(i,j-1,k) &
                   +         vx(i-1,j,  k) + vx(i,j,  k))
              v_y = fourth*(vy(i-1,j-1,k) + vy(i,j-1,k) &
                   +         vy(i-1,j,  k) + vy(i,j,  k))
              v_z = fourth*(vz(i-1,j-1,k) + vz(i,j-1,k) &
                   +         vz(i-1,j,  k) + vz(i,j,  k))

              w_x = fourth*(wx(i-1,j-1,k) + wx(i,j-1,k) &
                   +         wx(i-1,j,  k) + wx(i,j,  k))
              w_y = fourth*(wy(i-1,j-1,k) + wy(i,j-1,k) &
                   +         wy(i-1,j,  k) + wy(i,j,  k))
              w_z = fourth*(wz(i-1,j-1,k) + wz(i,j-1,k) &
                   +         wz(i-1,j,  k) + wz(i,j,  k))

              q_x = fourth*(qx(i-1,j-1,k) + qx(i,j-1,k) &
                   +         qx(i-1,j,  k) + qx(i,j,  k))
              q_y = fourth*(qy(i-1,j-1,k) + qy(i,j-1,k) &
                   +         qy(i-1,j,  k) + qy(i,j,  k))
              q_z = fourth*(qz(i-1,j-1,k) + qz(i,j-1,k) &
                   +         qz(i-1,j,  k) + qz(i,j,  k))


              ! The gradients in the normal direction are corrected, such
              ! that no averaging takes places here.
              ! First determine the vector in the direction from the
              ! cell center k to cell center k+1.

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

              ! Determine the length of this vector and create the
              ! unit normal.

              ss  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
              ssx = ss*ssx
              ssy = ss*ssy
              ssz = ss*ssz

              ! Correct the gradients.

              corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                   - (w(i,j,k+1,ivx) - w(i,j,k,ivx))*ss
              u_x  = u_x - corr*ssx
              u_y  = u_y - corr*ssy
              u_z  = u_z - corr*ssz

              corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                   - (w(i,j,k+1,ivy) - w(i,j,k,ivy))*ss
              v_x  = v_x - corr*ssx
              v_y  = v_y - corr*ssy
              v_z  = v_z - corr*ssz

              corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                   - (w(i,j,k+1,ivz) - w(i,j,k,ivz))*ss
              w_x  = w_x - corr*ssx
              w_y  = w_y - corr*ssy
              w_z  = w_z - corr*ssz

              corr = q_x*ssx + q_y*ssy + q_z*ssz &
                   + (aa(i,j,k+1) - aa(i,j,k))*ss
              q_x  = q_x - corr*ssx
              q_y  = q_y - corr*ssy
              q_z  = q_z - corr*ssz

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

              ! Compute the viscous fluxes for this k-face.

              fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) &
                   + tauxz*sk(i,j,k,3)
              fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) &
                   + tauyz*sk(i,j,k,3)
              fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) &
                   + tauzz*sk(i,j,k,3)
              frhoE =         (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1)
              frhoE = frhoE + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2) 
              frhoE = frhoE + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3) 
              frhoE = frhoE -  q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

              ! Update the residuals of cell k and k+1.

              fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
              fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
              fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
              fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

              fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
              fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
              fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
              fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

              ! Store the stress tensor and the heat flux vector if this
              ! face is part of a viscous subface. Both the cases k == 1
              ! and k == kl must be tested.

              if(k == 1 .and. storeWallTensor .and. &
                   viscKminPointer(i,j) > 0) then
                 ! We need to index viscSubface with viscKminPointer(i,j) 
                 ! since Tapenade does not like temporary indexes 

                 viscSubface(viscKminPointer(i,j))%tau(i,j,1) = tauxx
                 viscSubface(viscKminPointer(i,j))%tau(i,j,2) = tauyy
                 viscSubface(viscKminPointer(i,j))%tau(i,j,3) = tauzz
                 viscSubface(viscKminPointer(i,j))%tau(i,j,4) = tauxy
                 viscSubface(viscKminPointer(i,j))%tau(i,j,5) = tauxz
                 viscSubface(viscKminPointer(i,j))%tau(i,j,6) = tauyz

                 viscSubface(viscKminPointer(i,j))%q(i,j,1) = q_x
                 viscSubface(viscKminPointer(i,j))%q(i,j,2) = q_y
                 viscSubface(viscKminPointer(i,j))%q(i,j,3) = q_z
              endif

              ! And the k == kl case.
              if(k == kl .and. storeWallTensor .and. &
                   viscKmaxPointer(i,j) > 0) then
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,1) = tauxx
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,2) = tauyy
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,3) = tauzz
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,4) = tauxy
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,5) = tauxz
                 viscSubface(viscKmaxPointer(i,j))%tau(i,j,6) = tauyz

                 viscSubface(viscKmaxPointer(i,j))%q(i,j,1) = q_x
                 viscSubface(viscKmaxPointer(i,j))%q(i,j,2) = q_y
                 viscSubface(viscKmaxPointer(i,j))%q(i,j,3) = q_z
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif   
  continue
  !$AD CHECKPOINT-END

  !
  !         Viscous fluxes in the j-direction.                           
  !
  continue
  !$AD CHECKPOINT-START
  mue = zero
#ifdef TAPENADE_FAST
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

              ! Set the value of the porosity. If not zero, it is set
              ! to average the eddy-viscosity and to take the factor
              ! rFilv into account.

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

              ! Compute the gradients at the face by averaging the four
              ! nodal values.

              u_x = fourth*(ux(i-1,j,k-1) + ux(i,j,k-1) &
                   +        ux(i-1,j,k  ) + ux(i,j,k  ))
              u_y = fourth*(uy(i-1,j,k-1) + uy(i,j,k-1) &
                   +        uy(i-1,j,k  ) + uy(i,j,k  ))
              u_z = fourth*(uz(i-1,j,k-1) + uz(i,j,k-1) &
                   +        uz(i-1,j,k  ) + uz(i,j,k  ))

              v_x = fourth*(vx(i-1,j,k-1) + vx(i,j,k-1) &
                   +        vx(i-1,j,k  ) + vx(i,j,k  ))
              v_y = fourth*(vy(i-1,j,k-1) + vy(i,j,k-1) &
                   +        vy(i-1,j,k  ) + vy(i,j,k  ))
              v_z = fourth*(vz(i-1,j,k-1) + vz(i,j,k-1) &
                   +        vz(i-1,j,k  ) + vz(i,j,k  ))

              w_x = fourth*(wx(i-1,j,k-1) + wx(i,j,k-1) &
                   +        wx(i-1,j,k  ) + wx(i,j,k  ))
              w_y = fourth*(wy(i-1,j,k-1) + wy(i,j,k-1) &
                   +        wy(i-1,j,k  ) + wy(i,j,k  ))
              w_z = fourth*(wz(i-1,j,k-1) + wz(i,j,k-1) &
                   +        wz(i-1,j,k  ) + wz(i,j,k  ))

              q_x = fourth*(qx(i-1,j,k-1) + qx(i,j,k-1) &
                   +        qx(i-1,j,k  ) + qx(i,j,k  ))
              q_y = fourth*(qy(i-1,j,k-1) + qy(i,j,k-1) &
                   +        qy(i-1,j,k  ) + qy(i,j,k  ))
              q_z = fourth*(qz(i-1,j,k-1) + qz(i,j,k-1) &
                   +        qz(i-1,j,k  ) + qz(i,j,k  ))

              ! The gradients in the normal direction are corrected, such
              ! that no averaging takes places here.
              ! First determine the vector in the direction from the
              ! cell center j to cell center j+1.

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

              ! Determine the length of this vector and create the
              ! unit normal.

              ss  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
              ssx = ss*ssx
              ssy = ss*ssy
              ssz = ss*ssz

              ! Correct the gradients.

              corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                   - (w(i,j+1,k,ivx) - w(i,j,k,ivx))*ss
              u_x  = u_x - corr*ssx
              u_y  = u_y - corr*ssy
              u_z  = u_z - corr*ssz

              corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                   - (w(i,j+1,k,ivy) - w(i,j,k,ivy))*ss
              v_x  = v_x - corr*ssx
              v_y  = v_y - corr*ssy
              v_z  = v_z - corr*ssz

              corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                   - (w(i,j+1,k,ivz) - w(i,j,k,ivz))*ss
              w_x  = w_x - corr*ssx
              w_y  = w_y - corr*ssy
              w_z  = w_z - corr*ssz

              corr = q_x*ssx + q_y*ssy + q_z*ssz &
                   + (aa(i,j+1,k) - aa(i,j,k))*ss
              q_x  = q_x - corr*ssx
              q_y  = q_y - corr*ssy
              q_z  = q_z - corr*ssz

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

              fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) &
                   + tauxz*sj(i,j,k,3)
              fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) &
                   + tauyz*sj(i,j,k,3)
              fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) &
                   + tauzz*sj(i,j,k,3)
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

              ! Store the stress tensor and the heat flux vector if this
              ! face is part of a viscous subface. Both the cases j == 1
              ! and j == jl must be tested.

              if(j == 1 .and. storeWallTensor .and. &
                   viscJminPointer(i,k) > 0) then
                 ! We need to index viscSubface with viscJminPointer(i,k) 
                 ! since Tapenade does not like temporary indexes 

                 viscSubface(viscJminPointer(i,k))%tau(i,k,1) = tauxx
                 viscSubface(viscJminPointer(i,k))%tau(i,k,2) = tauyy
                 viscSubface(viscJminPointer(i,k))%tau(i,k,3) = tauzz
                 viscSubface(viscJminPointer(i,k))%tau(i,k,4) = tauxy
                 viscSubface(viscJminPointer(i,k))%tau(i,k,5) = tauxz
                 viscSubface(viscJminPointer(i,k))%tau(i,k,6) = tauyz

                 viscSubface(viscJminPointer(i,k))%q(i,k,1) = q_x
                 viscSubface(viscJminPointer(i,k))%q(i,k,2) = q_y
                 viscSubface(viscJminPointer(i,k))%q(i,k,3) = q_z
              endif

              ! And the j == jl case.

              if(j == jl .and. storeWallTensor .and. &
                   viscJmaxPointer(i,k) > 0) then
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,1) = tauxx
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,2) = tauyy
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,3) = tauzz
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,4) = tauxy
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,5) = tauxz
                 viscSubface(viscJmaxPointer(i,k))%tau(i,k,6) = tauyz

                 viscSubface(viscJmaxPointer(i,k))%q(i,k,1) = q_x
                 viscSubface(viscJmaxPointer(i,k))%q(i,k,2) = q_y
                 viscSubface(viscJmaxPointer(i,k))%q(i,k,3) = q_z
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif 
  continue
  !$AD CHECKPOINT-END



  !
  !         Viscous fluxes in the i-direction.                           
  !
  continue
  !$AD CHECKPOINT-START
  mue= zero
#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,il*ny*nz-1
     i = mod(ii, il) + 1
     j = mod(ii/il, ny) + 2
     k = ii/(il*ny) + 2
#else
     do k=2, kl
        do j=2, jl
           do i=1, il
#endif 

              ! Set the value of the porosity. If not zero, it is set
              ! to average the eddy-viscosity and to take the factor
              ! rFilv into account.

              por = half*rFilv
              if(porI(i,j,k) == noFlux) por = zero

              ! Compute the laminar and (if present) the eddy viscosities
              ! multiplied the porosity. Compute the factor in front of
              ! the gradients of the speed of sound squared for the heat
              ! flux.

              mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
              if( eddyModel ) mue = por*(rev(i,j,k) + rev(i+1,j,k))
              mut = mul + mue

              gm1          = half*(gamma(i,j,k) + gamma(i+1,j,k)) - one
              factLamHeat  = one/(prandtl*gm1)
              factTurbHeat = one/(prandtlTurb*gm1)

              heatCoef = mul*factLamHeat + mue*factTurbHeat

              ! Compute the gradients at the face by averaging the four
              ! nodal values.

              u_x = fourth*(ux(i,j-1,k-1) + ux(i,j,k-1) &
                   +        ux(i,j-1,k  ) + ux(i,j,k  ))
              u_y = fourth*(uy(i,j-1,k-1) + uy(i,j,k-1) &
                   +        uy(i,j-1,k  ) + uy(i,j,k  ))
              u_z = fourth*(uz(i,j-1,k-1) + uz(i,j,k-1) &
                   +        uz(i,j-1,k  ) + uz(i,j,k  ))

              v_x = fourth*(vx(i,j-1,k-1) + vx(i,j,k-1) &
                   +        vx(i,j-1,k  ) + vx(i,j,k  ))
              v_y = fourth*(vy(i,j-1,k-1) + vy(i,j,k-1) &
                   +        vy(i,j-1,k  ) + vy(i,j,k  ))
              v_z = fourth*(vz(i,j-1,k-1) + vz(i,j,k-1) &
                   +        vz(i,j-1,k  ) + vz(i,j,k  ))

              w_x = fourth*(wx(i,j-1,k-1) + wx(i,j,k-1) &
                   +        wx(i,j-1,k  ) + wx(i,j,k  ))
              w_y = fourth*(wy(i,j-1,k-1) + wy(i,j,k-1) &
                   +        wy(i,j-1,k  ) + wy(i,j,k  ))
              w_z = fourth*(wz(i,j-1,k-1) + wz(i,j,k-1) &
                   +        wz(i,j-1,k  ) + wz(i,j,k  ))

              q_x = fourth*(qx(i,j-1,k-1) + qx(i,j,k-1) &
                   +        qx(i,j-1,k  ) + qx(i,j,k  ))
              q_y = fourth*(qy(i,j-1,k-1) + qy(i,j,k-1) &
                   +        qy(i,j-1,k  ) + qy(i,j,k  ))
              q_z = fourth*(qz(i,j-1,k-1) + qz(i,j,k-1) &
                   +        qz(i,j-1,k  ) + qz(i,j,k  ))

              ! The gradients in the normal direction are corrected, such
              ! that no averaging takes places here.
              ! First determine the vector in the direction from the
              ! cell center i to cell center i+1.

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

              ! Determine the length of this vector and create the
              ! unit normal.

              ss  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
              ssx = ss*ssx
              ssy = ss*ssy
              ssz = ss*ssz

              ! Correct the gradients.

              corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                   - (w(i+1,j,k,ivx) - w(i,j,k,ivx))*ss
              u_x  = u_x - corr*ssx
              u_y  = u_y - corr*ssy
              u_z  = u_z - corr*ssz

              corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                   - (w(i+1,j,k,ivy) - w(i,j,k,ivy))*ss
              v_x  = v_x - corr*ssx
              v_y  = v_y - corr*ssy
              v_z  = v_z - corr*ssz

              corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                   - (w(i+1,j,k,ivz) - w(i,j,k,ivz))*ss
              w_x  = w_x - corr*ssx
              w_y  = w_y - corr*ssy
              w_z  = w_z - corr*ssz

              corr = q_x*ssx + q_y*ssy + q_z*ssz &
                   + (aa(i+1,j,k) - aa(i,j,k))*ss
              q_x  = q_x - corr*ssx
              q_y  = q_y - corr*ssy
              q_z  = q_z - corr*ssz

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

              fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) &
                   + tauxz*si(i,j,k,3)
              fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) &
                   + tauyz*si(i,j,k,3)
              fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) &
                   + tauzz*si(i,j,k,3)
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

              ! Store the stress tensor and the heat flux vector if this
              ! face is part of a viscous subface. Both the cases i == 1
              ! and i == il must be tested.

              if(i == 1 .and. storeWallTensor .and. &
                   viscIminPointer(j,k) > 0) then
                 ! We need to index viscSubface with viscIminPointer(j,k) 
                 ! since Tapenade does not like temporary indexes 

                 viscSubface(viscIminPointer(j,k))%tau(j,k,1) = tauxx
                 viscSubface(viscIminPointer(j,k))%tau(j,k,2) = tauyy
                 viscSubface(viscIminPointer(j,k))%tau(j,k,3) = tauzz
                 viscSubface(viscIminPointer(j,k))%tau(j,k,4) = tauxy
                 viscSubface(viscIminPointer(j,k))%tau(j,k,5) = tauxz
                 viscSubface(viscIminPointer(j,k))%tau(j,k,6) = tauyz

                 viscSubface(viscIminPointer(j,k))%q(j,k,1) = q_x
                 viscSubface(viscIminPointer(j,k))%q(j,k,2) = q_y
                 viscSubface(viscIminPointer(j,k))%q(j,k,3) = q_z
              endif

              ! And the i == il case.

              if(i == il .and. storeWallTensor .and. &
                   viscImaxPointer(j,k) > 0) then
                 ! We need to index viscSubface with viscImaxPointer(j,k) 
                 ! since Tapenade does not like temporary indexes 

                 viscSubface(viscImaxPointer(j,k))%tau(j,k,1) = tauxx
                 viscSubface(viscImaxPointer(j,k))%tau(j,k,2) = tauyy
                 viscSubface(viscImaxPointer(j,k))%tau(j,k,3) = tauzz
                 viscSubface(viscImaxPointer(j,k))%tau(j,k,4) = tauxy
                 viscSubface(viscImaxPointer(j,k))%tau(j,k,5) = tauxz
                 viscSubface(viscImaxPointer(j,k))%tau(j,k,6) = tauyz

                 viscSubface(viscImaxPointer(j,k))%q(j,k,1) = q_x
                 viscSubface(viscImaxPointer(j,k))%q(j,k,2) = q_y
                 viscSubface(viscImaxPointer(j,k))%q(j,k,3) = q_z
              endif
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif 
  !$AD CHECKPOINT-END
  continue
  ! Possibly correct the wall shear stress.
  ! Wall function is not ADed
#ifndef USE_TAPENADE       
  call utauWF(rFilv)
#endif
end subroutine viscousFlux

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
  !
  !         Begin execution                                              
  !

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

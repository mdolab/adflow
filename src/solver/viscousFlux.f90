!
!      ******************************************************************
!      *                                                                *
!      * File:          viscousFlux.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-21-2003                                      *
!      * Last modified: 04-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine viscousFlux
!
!      ******************************************************************
!      *                                                                *
!      * viscousFlux computes the viscous fluxes using a central        *
!      * difference scheme for a block.                                 *
!      * It is assumed that the pointers in block pointer already point *
!      * to the correct block.                                          *
!      *                                                                *
!      ******************************************************************
!
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
       integer(kind=intType) :: i, j, k, nn
       integer(kind=intType) :: k1, k2, kk

       real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
       real(kind=realType) :: gm1, factLamHeat, factTurbHeat
       real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
       real(kind=realType) :: q_x, q_y, q_z, ubar, vbar, wbar
       real(kind=realType) :: corr, ssx, ssy, ssz, ss, fracDiv
       real(kind=realType) :: tauxx, tauyy, tauzz
       real(kind=realType) :: tauxy, tauxz, tauyz
       real(kind=realType) :: fmx, fmy, fmz, frhoE

       real(kind=realType), dimension(il,jl,2) :: ux, uy, uz
       real(kind=realType), dimension(il,jl,2) :: vx, vy, vz
       real(kind=realType), dimension(il,jl,2) :: wx, wy, wz
       real(kind=realType), dimension(il,jl,2) :: qx, qy, qz

       logical :: correctForK, storeWallTensor
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set rFilv to rFil to indicate that this is the viscous part.
       ! If rFilv == 0 the viscous residuals need not to be computed
       ! and a return can be made.

       rFilv = rFil
       if(rFilv == zero) return

       ! Determine whether or not the pressure must be corrected
       ! for the presence of the turbulent kinetic energy.

       if( kPresent ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
           correctForK = .true.
         else
           correctForK = .false.
         endif
       else
         correctForK = .false.
       endif

       ! Determine whether or not the wall stress tensor and wall heat
       ! flux must be stored for viscous walls.

       storeWallTensor = .false.
       if( wallFunctions ) then
         storeWallTensor = .true.
       else if(rkStage == 0 .and. currentLevel == groundLevel) then
         storeWallTensor = .true.
       endif

       ! Store the speed of sound squared instead of the pressure.
       ! To be 100 percent correct, substract 2/3*rho*k (if present)
       ! from the pressure to obtain the true presssure. First layer of
       ! halo's, because that's what is needed by the viscous stencil.

       do k=1,ke
         do j=1,je
           do i=1,ie
             if( correctForK ) &
               p(i,j,k) = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
             p(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
           enddo
         enddo
       enddo

       ! Compute the nodal gradients for the nodes in the plane k = 1.

       k  = 1
       k1 = 1
       k2 = 2
       call nodalGradients(ux, uy, uz, vx, vy, vz, wx, wy, wz, &
                           qx, qy, qz)

       ! Compute the viscous fluxes for the faces k == 1 and update
       ! the residuals for the cells k == 2.

       mue = zero
       do j=2,jl
         do i=2,il

           ! Set the value of the porosity. If not zero, it is set
           ! to average the eddy-viscosity and to take the factor
           ! rFilv into account.

           por = half*rFilv
           if(porK(i,j,1) == noFlux) por = zero

           ! Compute the laminar and (if present) the eddy viscosities
           ! multiplied by the porosity. Compute the factor in front of
           ! the gradients of the speed of sound squared for the heat
           ! flux.

           mul = por*(rlv(i,j,1) + rlv(i,j,2))
           if( eddyModel ) mue = por*(rev(i,j,1) + rev(i,j,2))
           mut = mul + mue

           gm1          = half*(gamma(i,j,1) + gamma(i,j,2)) - one
           factLamHeat  = one/(prandtl*gm1)
           factTurbHeat = one/(prandtlTurb*gm1)

           heatCoef = mul*factLamHeat + mue*factTurbHeat

           ! Compute the gradients at the face by averaging the four
           ! nodal values.

           u_x = fourth*(ux(i-1,j-1,k1) + ux(i,j-1,k1) &
               +         ux(i-1,j,  k1) + ux(i,j,  k1))
           u_y = fourth*(uy(i-1,j-1,k1) + uy(i,j-1,k1) &
               +         uy(i-1,j,  k1) + uy(i,j,  k1))
           u_z = fourth*(uz(i-1,j-1,k1) + uz(i,j-1,k1) &
               +         uz(i-1,j,  k1) + uz(i,j,  k1))

           v_x = fourth*(vx(i-1,j-1,k1) + vx(i,j-1,k1) &
               +         vx(i-1,j,  k1) + vx(i,j,  k1))
           v_y = fourth*(vy(i-1,j-1,k1) + vy(i,j-1,k1) &
               +         vy(i-1,j,  k1) + vy(i,j,  k1))
           v_z = fourth*(vz(i-1,j-1,k1) + vz(i,j-1,k1) &
               +         vz(i-1,j,  k1) + vz(i,j,  k1))

           w_x = fourth*(wx(i-1,j-1,k1) + wx(i,j-1,k1) &
               +         wx(i-1,j,  k1) + wx(i,j,  k1))
           w_y = fourth*(wy(i-1,j-1,k1) + wy(i,j-1,k1) &
               +         wy(i-1,j,  k1) + wy(i,j,  k1))
           w_z = fourth*(wz(i-1,j-1,k1) + wz(i,j-1,k1) &
               +         wz(i-1,j,  k1) + wz(i,j,  k1))

           q_x = fourth*(qx(i-1,j-1,k1) + qx(i,j-1,k1) &
               +         qx(i-1,j,  k1) + qx(i,j,  k1))
           q_y = fourth*(qy(i-1,j-1,k1) + qy(i,j-1,k1) &
               +         qy(i-1,j,  k1) + qy(i,j,  k1))
           q_z = fourth*(qz(i-1,j-1,k1) + qz(i,j-1,k1) &
               +         qz(i-1,j,  k1) + qz(i,j,  k1))

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
                + (p(i,j,k+1) - p(i,j,k))*ss
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

           ubar = half*(w(i,j,1,ivx) + w(i,j,2,ivx))
           vbar = half*(w(i,j,1,ivy) + w(i,j,2,ivy))
           wbar = half*(w(i,j,1,ivz) + w(i,j,2,ivz))

           ! Compute the viscous fluxes for this k-face.

           fmx   = tauxx*sk(i,j,1,1) + tauxy*sk(i,j,1,2) &
                 + tauxz*sk(i,j,1,3)
           fmy   = tauxy*sk(i,j,1,1) + tauyy*sk(i,j,1,2) &
                 + tauyz*sk(i,j,1,3)
           fmz   = tauxz*sk(i,j,1,1) + tauyz*sk(i,j,1,2) &
                 + tauzz*sk(i,j,1,3)
           frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,1,1) &
                 + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,1,2) &
                 + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,1,3) &
                 - q_x*sk(i,j,1,1) - q_y*sk(i,j,1,2) - q_z*sk(i,j,1,3)

           ! Update the residuals of cell k == 2.

           fw(i,j,2,imx)   = fw(i,j,2,imx)   + fmx
           fw(i,j,2,imy)   = fw(i,j,2,imy)   + fmy
           fw(i,j,2,imz)   = fw(i,j,2,imz)   + fmz
           fw(i,j,2,irhoE) = fw(i,j,2,irhoE) + frhoE

           ! Store the stress tensor and the heat flux vector if this
           ! face is part of a viscous subface.

           if(storeWallTensor .and. viscKminPointer(i,j) > 0) then
             nn = viscKminPointer(i,j)

             viscSubface(nn)%tau(i,j,1) = tauxx
             viscSubface(nn)%tau(i,j,2) = tauyy
             viscSubface(nn)%tau(i,j,3) = tauzz
             viscSubface(nn)%tau(i,j,4) = tauxy
             viscSubface(nn)%tau(i,j,5) = tauxz
             viscSubface(nn)%tau(i,j,6) = tauyz

             viscSubface(nn)%q(i,j,1) = q_x
             viscSubface(nn)%q(i,j,2) = q_y
             viscSubface(nn)%q(i,j,3) = q_z
           endif

         enddo
       enddo

       ! Loop over the k-planes.

       kLoop: do k=2,kl

         ! Switch the indices k1 and k2.

         kk = k1
         k1 = k2
         k2 = kk

         ! Compute the nodal gradients for the nodes in this k-plane.
         ! The results are stored in ux(:,:,k1), etc.

         call nodalGradients(ux, uy, uz, vx, vy, vz, wx, wy, wz, &
                             qx, qy, qz)
!
!        ****************************************************************
!        *                                                              *
!        * Viscous fluxes in the k-direction.                           *
!        *                                                              *
!        ****************************************************************
!
         do j=2,jl
           do i=2,il

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

             u_x = fourth*(ux(i-1,j-1,k1) + ux(i,j-1,k1) &
                 +         ux(i-1,j,  k1) + ux(i,j,  k1))
             u_y = fourth*(uy(i-1,j-1,k1) + uy(i,j-1,k1) &
                 +         uy(i-1,j,  k1) + uy(i,j,  k1))
             u_z = fourth*(uz(i-1,j-1,k1) + uz(i,j-1,k1) &
                 +         uz(i-1,j,  k1) + uz(i,j,  k1))

             v_x = fourth*(vx(i-1,j-1,k1) + vx(i,j-1,k1) &
                 +         vx(i-1,j,  k1) + vx(i,j,  k1))
             v_y = fourth*(vy(i-1,j-1,k1) + vy(i,j-1,k1) &
                 +         vy(i-1,j,  k1) + vy(i,j,  k1))
             v_z = fourth*(vz(i-1,j-1,k1) + vz(i,j-1,k1) &
                 +         vz(i-1,j,  k1) + vz(i,j,  k1))

             w_x = fourth*(wx(i-1,j-1,k1) + wx(i,j-1,k1) &
                 +         wx(i-1,j,  k1) + wx(i,j,  k1))
             w_y = fourth*(wy(i-1,j-1,k1) + wy(i,j-1,k1) &
                 +         wy(i-1,j,  k1) + wy(i,j,  k1))
             w_z = fourth*(wz(i-1,j-1,k1) + wz(i,j-1,k1) &
                 +         wz(i-1,j,  k1) + wz(i,j,  k1))

             q_x = fourth*(qx(i-1,j-1,k1) + qx(i,j-1,k1) &
                 +         qx(i-1,j,  k1) + qx(i,j,  k1))
             q_y = fourth*(qy(i-1,j-1,k1) + qy(i,j-1,k1) &
                 +         qy(i-1,j,  k1) + qy(i,j,  k1))
             q_z = fourth*(qz(i-1,j-1,k1) + qz(i,j-1,k1) &
                 +         qz(i-1,j,  k1) + qz(i,j,  k1))

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
                  + (p(i,j,k+1) - p(i,j,k))*ss
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
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1) &
                   + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2) &
                   + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3) &
                   - q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

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
             ! face is part of a viscous subface.

             if(k == kl .and. storeWallTensor .and. &
                viscKmaxPointer(i,j) > 0) then
               nn = viscKmaxPointer(i,j)

               viscSubface(nn)%tau(i,j,1) = tauxx
               viscSubface(nn)%tau(i,j,2) = tauyy
               viscSubface(nn)%tau(i,j,3) = tauzz
               viscSubface(nn)%tau(i,j,4) = tauxy
               viscSubface(nn)%tau(i,j,5) = tauxz
               viscSubface(nn)%tau(i,j,6) = tauyz

               viscSubface(nn)%q(i,j,1) = q_x
               viscSubface(nn)%q(i,j,2) = q_y
               viscSubface(nn)%q(i,j,3) = q_z
             endif

           enddo
         enddo
!
!        ****************************************************************
!        *                                                              *
!        * Viscous fluxes in the j-direction.                           *
!        *                                                              *
!        ****************************************************************
!
         do j=1,jl
           do i=2,il

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

             u_x = fourth*(ux(i-1,j,k1) + ux(i,j,k1) &
                 +         ux(i-1,j,k2) + ux(i,j,k2))
             u_y = fourth*(uy(i-1,j,k1) + uy(i,j,k1) &
                 +         uy(i-1,j,k2) + uy(i,j,k2))
             u_z = fourth*(uz(i-1,j,k1) + uz(i,j,k1) &
                 +         uz(i-1,j,k2) + uz(i,j,k2))

             v_x = fourth*(vx(i-1,j,k1) + vx(i,j,k1) &
                 +         vx(i-1,j,k2) + vx(i,j,k2))
             v_y = fourth*(vy(i-1,j,k1) + vy(i,j,k1) &
                 +         vy(i-1,j,k2) + vy(i,j,k2))
             v_z = fourth*(vz(i-1,j,k1) + vz(i,j,k1) &
                 +         vz(i-1,j,k2) + vz(i,j,k2))

             w_x = fourth*(wx(i-1,j,k1) + wx(i,j,k1) &
                 +         wx(i-1,j,k2) + wx(i,j,k2))
             w_y = fourth*(wy(i-1,j,k1) + wy(i,j,k1) &
                 +         wy(i-1,j,k2) + wy(i,j,k2))
             w_z = fourth*(wz(i-1,j,k1) + wz(i,j,k1) &
                 +         wz(i-1,j,k2) + wz(i,j,k2))

             q_x = fourth*(qx(i-1,j,k1) + qx(i,j,k1) &
                 +         qx(i-1,j,k2) + qx(i,j,k2))
             q_y = fourth*(qy(i-1,j,k1) + qy(i,j,k1) &
                 +         qy(i-1,j,k2) + qy(i,j,k2))
             q_z = fourth*(qz(i-1,j,k1) + qz(i,j,k1) &
                 +         qz(i-1,j,k2) + qz(i,j,k2))

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
                  + (p(i,j+1,k) - p(i,j,k))*ss
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
               nn = viscJminPointer(i,k)

               viscSubface(nn)%tau(i,k,1) = tauxx
               viscSubface(nn)%tau(i,k,2) = tauyy
               viscSubface(nn)%tau(i,k,3) = tauzz
               viscSubface(nn)%tau(i,k,4) = tauxy
               viscSubface(nn)%tau(i,k,5) = tauxz
               viscSubface(nn)%tau(i,k,6) = tauyz

               viscSubface(nn)%q(i,k,1) = q_x
               viscSubface(nn)%q(i,k,2) = q_y
               viscSubface(nn)%q(i,k,3) = q_z
             endif

             ! And the j == jl case.

             if(j == jl .and. storeWallTensor .and. &
                viscJmaxPointer(i,k) > 0) then
               nn = viscJmaxPointer(i,k)

               viscSubface(nn)%tau(i,k,1) = tauxx
               viscSubface(nn)%tau(i,k,2) = tauyy
               viscSubface(nn)%tau(i,k,3) = tauzz
               viscSubface(nn)%tau(i,k,4) = tauxy
               viscSubface(nn)%tau(i,k,5) = tauxz
               viscSubface(nn)%tau(i,k,6) = tauyz

               viscSubface(nn)%q(i,k,1) = q_x
               viscSubface(nn)%q(i,k,2) = q_y
               viscSubface(nn)%q(i,k,3) = q_z
             endif

           enddo
         enddo
!
!        ****************************************************************
!        *                                                              *
!        * Viscous fluxes in the i-direction.                           *
!        *                                                              *
!        ****************************************************************
!
         do j=2,jl
           do i=1,il

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

             u_x = fourth*(ux(i,j-1,k1) + ux(i,j,k1) &
                 +         ux(i,j-1,k2) + ux(i,j,k2))
             u_y = fourth*(uy(i,j-1,k1) + uy(i,j,k1) &
                 +         uy(i,j-1,k2) + uy(i,j,k2))
             u_z = fourth*(uz(i,j-1,k1) + uz(i,j,k1) &
                 +         uz(i,j-1,k2) + uz(i,j,k2))

             v_x = fourth*(vx(i,j-1,k1) + vx(i,j,k1) &
                 +         vx(i,j-1,k2) + vx(i,j,k2))
             v_y = fourth*(vy(i,j-1,k1) + vy(i,j,k1) &
                 +         vy(i,j-1,k2) + vy(i,j,k2))
             v_z = fourth*(vz(i,j-1,k1) + vz(i,j,k1) &
                 +         vz(i,j-1,k2) + vz(i,j,k2))

             w_x = fourth*(wx(i,j-1,k1) + wx(i,j,k1) &
                 +         wx(i,j-1,k2) + wx(i,j,k2))
             w_y = fourth*(wy(i,j-1,k1) + wy(i,j,k1) &
                 +         wy(i,j-1,k2) + wy(i,j,k2))
             w_z = fourth*(wz(i,j-1,k1) + wz(i,j,k1) &
                 +         wz(i,j-1,k2) + wz(i,j,k2))

             q_x = fourth*(qx(i,j-1,k1) + qx(i,j,k1) &
                 +         qx(i,j-1,k2) + qx(i,j,k2))
             q_y = fourth*(qy(i,j-1,k1) + qy(i,j,k1) &
                 +         qy(i,j-1,k2) + qy(i,j,k2))
             q_z = fourth*(qz(i,j-1,k1) + qz(i,j,k1) &
                 +         qz(i,j-1,k2) + qz(i,j,k2))

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
                  + (p(i+1,j,k) - p(i,j,k))*ss
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
               nn = viscIminPointer(j,k)

               viscSubface(nn)%tau(j,k,1) = tauxx
               viscSubface(nn)%tau(j,k,2) = tauyy
               viscSubface(nn)%tau(j,k,3) = tauzz
               viscSubface(nn)%tau(j,k,4) = tauxy
               viscSubface(nn)%tau(j,k,5) = tauxz
               viscSubface(nn)%tau(j,k,6) = tauyz

               viscSubface(nn)%q(j,k,1) = q_x
               viscSubface(nn)%q(j,k,2) = q_y
               viscSubface(nn)%q(j,k,3) = q_z
             endif

             ! And the i == il case.

             if(i == il .and. storeWallTensor .and. &
                viscImaxPointer(j,k) > 0) then
               nn = viscImaxPointer(j,k)

               viscSubface(nn)%tau(j,k,1) = tauxx
               viscSubface(nn)%tau(j,k,2) = tauyy
               viscSubface(nn)%tau(j,k,3) = tauzz
               viscSubface(nn)%tau(j,k,4) = tauxy
               viscSubface(nn)%tau(j,k,5) = tauxz
               viscSubface(nn)%tau(j,k,6) = tauyz

               viscSubface(nn)%q(j,k,1) = q_x
               viscSubface(nn)%q(j,k,2) = q_y
               viscSubface(nn)%q(j,k,3) = q_z
             endif

           enddo
         enddo

       enddo kLoop

       ! Restore the pressure in p. Again only the first layer of
       ! halo cells.

       do k=1,ke
         do j=1,je
           do i=1,ie
             p(i,j,k) = w(i,j,k,irho)*p(i,j,k)/gamma(i,j,k)
           enddo
         enddo
       enddo

       if( correctForK ) then
         do k=1,ke
           do j=1,je
             do i=1,ie
               p(i,j,k) = p(i,j,k) + twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
             enddo
           enddo
         enddo
       endif

       ! Possibly correct the wall shear stress.

       call utauWF(rFilv)

       contains

         subroutine nodalGradients(ux, uy, uz, vx, vy, vz, wx, wy, wz, &
                                   qx, qy, qz)
!
!        ****************************************************************
!        *                                                              *
!        * nodalGradients computes the nodal velocity gradients and     *
!        * minus the gradient of the speed of sound squared. The minus  *
!        * sign is present, because this is the definition of the heat  *
!        * flux. These gradients are computed for the k-plane k. The    *
!        * results are stored in ux(:,:,k1), etc. Here k1 is either 1   *
!        * or 2. The gradients have the intent(inout) label, because    *
!        * only the k1 elements are changed; the others remain the same.*
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments
!
         real(kind=realType), dimension(il,jl,2), intent(inout) :: &
                        ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz
!
!        Local variables.
!
         real(kind=realType) :: oneOverV, ubar, vbar, wbar, a2
         real(kind=realType) :: sx, sx1, sy, sy1, sz, sz1
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! First part. Contribution in the k-direction. As the viscous
         ! fluxes are computed per k-plane, the values of the surface
         ! integrals in gauss's theorem are only scattered to one node.
         ! Consequently this part is an initialization as well.

         ! Loop over the j and i direction.

         do j=1,jl
           do i=1,il

             ! Store 8 times the average normal for the contribution from
             ! the k-layer of cells. The factor 8 drops out later when
             ! dividing by the volume.

             sx1 = sk(i,j,  k,1) + sk(i+1,j,  k,1) &
                 + sk(i,j+1,k,1) + sk(i+1,j+1,k,1)
             sy1 = sk(i,j,  k,2) + sk(i+1,j,  k,2) &
                 + sk(i,j+1,k,2) + sk(i+1,j+1,k,2)
             sz1 = sk(i,j,  k,3) + sk(i+1,j,  k,3) &
                 + sk(i,j+1,k,3) + sk(i+1,j+1,k,3)

             sx = sx1 + sk(i,j,  k-1,1) + sk(i+1,j,  k-1,1) &
                +       sk(i,j+1,k-1,1) + sk(i+1,j+1,k-1,1)
             sy = sy1 + sk(i,j,  k-1,2) + sk(i+1,j,  k-1,2) &
                +       sk(i,j+1,k-1,2) + sk(i+1,j+1,k-1,2)
             sz = sz1 + sk(i,j,  k-1,3) + sk(i+1,j,  k-1,3) &
                +       sk(i,j+1,k-1,3) + sk(i+1,j+1,k-1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,  k,ivx) + w(i+1,j,  k,ivx) &
                  +         w(i,j+1,k,ivx) + w(i+1,j+1,k,ivx))
             vbar = fourth*(w(i,j,  k,ivy) + w(i+1,j,  k,ivy) &
                  +         w(i,j+1,k,ivy) + w(i+1,j+1,k,ivy))
             wbar = fourth*(w(i,j,  k,ivz) + w(i+1,j,  k,ivz) &
                  +         w(i,j+1,k,ivz) + w(i+1,j+1,k,ivz))

             a2 = fourth*(p(i,j,k) + p(i+1,j,k) + p(i,j+1,k) + p(i+1,j+1,k))

             ! Set the velocity and speed of sound square gradients.
             ! The minus sign is there, because these normals are inward
             ! pointing for the nodal k-layer. The exception is a2,
             ! because the gradient of -a2 is stored, as this is needed
             ! in the heat fluxes.

             ux(i,j,k1) = -ubar*sx
             uy(i,j,k1) = -ubar*sy
             uz(i,j,k1) = -ubar*sz

             vx(i,j,k1) = -vbar*sx
             vy(i,j,k1) = -vbar*sy
             vz(i,j,k1) = -vbar*sz

             wx(i,j,k1) = -wbar*sx
             wy(i,j,k1) = -wbar*sy
             wz(i,j,k1) = -wbar*sz

             qx(i,j,k1) = a2*sx
             qy(i,j,k1) = a2*sy
             qz(i,j,k1) = a2*sz

             ! Store 8 times the average normal for the contribution from
             ! the k+1 layer of cells. The factor 8 drops out later when
             ! dividing by the volume.

             sx = sx1 + sk(i,j,  k+1,1) + sk(i+1,j,  k+1,1) &
                +       sk(i,j+1,k+1,1) + sk(i+1,j+1,k+1,1)
             sy = sy1 + sk(i,j,  k+1,2) + sk(i+1,j,  k+1,2) &
                +       sk(i,j+1,k+1,2) + sk(i+1,j+1,k+1,2)
             sz = sz1 + sk(i,j,  k+1,3) + sk(i+1,j,  k+1,3) &
                +       sk(i,j+1,k+1,3) + sk(i+1,j+1,k+1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,  k+1,ivx) + w(i+1,j,  k+1,ivx) &
                  +         w(i,j+1,k+1,ivx) + w(i+1,j+1,k+1,ivx))
             vbar = fourth*(w(i,j,  k+1,ivy) + w(i+1,j,  k+1,ivy) &
                  +         w(i,j+1,k+1,ivy) + w(i+1,j+1,k+1,ivy))
             wbar = fourth*(w(i,j,  k+1,ivz) + w(i+1,j,  k+1,ivz) &
                  +         w(i,j+1,k+1,ivz) + w(i+1,j+1,k+1,ivz))

             a2 = fourth*(p(i,j,  k+1) + p(i+1,j,  k+1) &
                +         p(i,j+1,k+1) + p(i+1,j+1,k+1))

             ! Update the velocity and (minus) the speed of sound
             ! gradients. As the normals are now outward pointing for
             ! the nodes in the k-layer, there is a plus sign for the
             ! velocity gradients and a minus sign for minus the speed
             ! of sound squared.

             ux(i,j,k1) = ux(i,j,k1) + ubar*sx
             uy(i,j,k1) = uy(i,j,k1) + ubar*sy
             uz(i,j,k1) = uz(i,j,k1) + ubar*sz

             vx(i,j,k1) = vx(i,j,k1) + vbar*sx
             vy(i,j,k1) = vy(i,j,k1) + vbar*sy
             vz(i,j,k1) = vz(i,j,k1) + vbar*sz

             wx(i,j,k1) = wx(i,j,k1) + wbar*sx
             wy(i,j,k1) = wy(i,j,k1) + wbar*sy
             wz(i,j,k1) = wz(i,j,k1) + wbar*sz

             qx(i,j,k1) = qx(i,j,k1) - a2*sx
             qy(i,j,k1) = qy(i,j,k1) - a2*sy
             qz(i,j,k1) = qz(i,j,k1) - a2*sz

           enddo
         enddo

         ! Second part. Contribution in the j-direction.
         ! The contribution is scattered to both the left and right node
         ! in j-direction.

         do j=1,je
           do i=1,il

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

             a2 = fourth*(p(i,j,k) + p(i+1,j,k) + p(i,j,k+1) + p(i+1,j,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(j > 1) then
               ux(i,j-1,k1) = ux(i,j-1,k1) + ubar*sx
               uy(i,j-1,k1) = uy(i,j-1,k1) + ubar*sy
               uz(i,j-1,k1) = uz(i,j-1,k1) + ubar*sz

               vx(i,j-1,k1) = vx(i,j-1,k1) + vbar*sx
               vy(i,j-1,k1) = vy(i,j-1,k1) + vbar*sy
               vz(i,j-1,k1) = vz(i,j-1,k1) + vbar*sz

               wx(i,j-1,k1) = wx(i,j-1,k1) + wbar*sx
               wy(i,j-1,k1) = wy(i,j-1,k1) + wbar*sy
               wz(i,j-1,k1) = wz(i,j-1,k1) + wbar*sz

               qx(i,j-1,k1) = qx(i,j-1,k1) - a2*sx
               qy(i,j-1,k1) = qy(i,j-1,k1) - a2*sy
               qz(i,j-1,k1) = qz(i,j-1,k1) - a2*sz
             endif

             if(j < je) then
               ux(i,j,k1) = ux(i,j,k1) - ubar*sx
               uy(i,j,k1) = uy(i,j,k1) - ubar*sy
               uz(i,j,k1) = uz(i,j,k1) - ubar*sz

               vx(i,j,k1) = vx(i,j,k1) - vbar*sx
               vy(i,j,k1) = vy(i,j,k1) - vbar*sy
               vz(i,j,k1) = vz(i,j,k1) - vbar*sz

               wx(i,j,k1) = wx(i,j,k1) - wbar*sx
               wy(i,j,k1) = wy(i,j,k1) - wbar*sy
               wz(i,j,k1) = wz(i,j,k1) - wbar*sz

               qx(i,j,k1) = qx(i,j,k1) + a2*sx
               qy(i,j,k1) = qy(i,j,k1) + a2*sy
               qz(i,j,k1) = qz(i,j,k1) + a2*sz
             endif

           enddo
         enddo

         ! Third part. Contribution in the i-direction.
         ! The contribution is scattered to both the left and right node
         ! in i-direction.

         do j=1,jl
           do i=1,ie

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

             a2 = fourth*(p(i,j,k) + p(i,j+1,k) + p(i,j,k+1) + p(i,j+1,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(i > 1) then
               ux(i-1,j,k1) = ux(i-1,j,k1) + ubar*sx
               uy(i-1,j,k1) = uy(i-1,j,k1) + ubar*sy
               uz(i-1,j,k1) = uz(i-1,j,k1) + ubar*sz

               vx(i-1,j,k1) = vx(i-1,j,k1) + vbar*sx
               vy(i-1,j,k1) = vy(i-1,j,k1) + vbar*sy
               vz(i-1,j,k1) = vz(i-1,j,k1) + vbar*sz

               wx(i-1,j,k1) = wx(i-1,j,k1) + wbar*sx
               wy(i-1,j,k1) = wy(i-1,j,k1) + wbar*sy
               wz(i-1,j,k1) = wz(i-1,j,k1) + wbar*sz

               qx(i-1,j,k1) = qx(i-1,j,k1) - a2*sx
               qy(i-1,j,k1) = qy(i-1,j,k1) - a2*sy
               qz(i-1,j,k1) = qz(i-1,j,k1) - a2*sz
             endif

             if(i < ie) then
               ux(i,j,k1) = ux(i,j,k1) - ubar*sx
               uy(i,j,k1) = uy(i,j,k1) - ubar*sy
               uz(i,j,k1) = uz(i,j,k1) - ubar*sz

               vx(i,j,k1) = vx(i,j,k1) - vbar*sx
               vy(i,j,k1) = vy(i,j,k1) - vbar*sy
               vz(i,j,k1) = vz(i,j,k1) - vbar*sz

               wx(i,j,k1) = wx(i,j,k1) - wbar*sx
               wy(i,j,k1) = wy(i,j,k1) - wbar*sy
               wz(i,j,k1) = wz(i,j,k1) - wbar*sz

               qx(i,j,k1) = qx(i,j,k1) + a2*sx
               qy(i,j,k1) = qy(i,j,k1) + a2*sy
               qz(i,j,k1) = qz(i,j,k1) + a2*sz
             endif

           enddo
         enddo

         ! Divide by 8 times the volume to obtain the correct gradients.

         do j=1,jl
           do i=1,il

             ! Compute the inverse of 8 times the volume for this node.

             oneOverV = one/(vol(i,  j,  k) + vol(i,  j,  k+1) &
                      +      vol(i+1,j,  k) + vol(i+1,j,  k+1) &
                      +      vol(i,  j+1,k) + vol(i,  j+1,k+1) &
                      +      vol(i+1,j+1,k) + vol(i+1,j+1,k+1))

             ! Compute the correct velocity gradients and "unit" heat
             ! fluxes. The velocity gradients are stored in ux, etc.

             ux(i,j,k1) = ux(i,j,k1)*oneOverV
             uy(i,j,k1) = uy(i,j,k1)*oneOverV
             uz(i,j,k1) = uz(i,j,k1)*oneOverV

             vx(i,j,k1) = vx(i,j,k1)*oneOverV
             vy(i,j,k1) = vy(i,j,k1)*oneOverV
             vz(i,j,k1) = vz(i,j,k1)*oneOverV

             wx(i,j,k1) = wx(i,j,k1)*oneOverV
             wy(i,j,k1) = wy(i,j,k1)*oneOverV
             wz(i,j,k1) = wz(i,j,k1)*oneOverV

             qx(i,j,k1) = qx(i,j,k1)*oneOverV
             qy(i,j,k1) = qy(i,j,k1)*oneOverV
             qz(i,j,k1) = qz(i,j,k1)*oneOverV

           enddo
         enddo

         end subroutine nodalGradients

       end subroutine viscousFlux

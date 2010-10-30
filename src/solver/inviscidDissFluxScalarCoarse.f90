!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxScalarCoarse.f90                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidDissFluxScalarCoarse
!
!      ******************************************************************
!      *                                                                *
!      * inviscidDissFluxScalarCoarse computes the coarse grid, i.e.    *
!      * 1st order, artificial dissipation flux for the scalar          *
!      * dissipation scheme for a given block. Therefore it is assumed  *
!      * that the pointers in blockPointers already point to the        *
!      * correct block.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use inputDiscretization
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: sfil, fis0, dis0, ppor, fs, rhoi
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

       ! Set a couple of constants for the scheme.

       fis0 = rFil*vis2Coarse
       sfil = one - rFil

       ! Replace the total energy by rho times the total enthalpy.
       ! In this way the numerical solution is total enthalpy preserving
       ! for the steady Euler equations. Also replace the velocities by
       ! the momentum. As only first order halo's are needed, only include
       ! the first order halo's.

       do k=1,ke
         do j=1,je
           do i=1,ie
             w(i,j,k,ivx)   = w(i,j,k,irho)*w(i,j,k,ivx)
             w(i,j,k,ivy)   = w(i,j,k,irho)*w(i,j,k,ivy)
             w(i,j,k,ivz)   = w(i,j,k,irho)*w(i,j,k,ivz)
             w(i,j,k,irhoE) = w(i,j,k,irhoE) + p(i,j,k)
           enddo
         enddo
       enddo

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

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half

             dis0 = fis0*ppor*(radI(i,j,k) + radI(i+1,j,k))

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs = dis0*(w(i+1,j,k,irho) - w(i,j,k,irho))
             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs = dis0*(w(i+1,j,k,ivx) - w(i,j,k,ivx))
             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs = dis0*(w(i+1,j,k,ivy) - w(i,j,k,ivy))
             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs = dis0*(w(i+1,j,k,ivz) - w(i,j,k,ivz))
             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs = dis0*(w(i+1,j,k,irhoE) - w(i,j,k,irhoE))
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

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half

             dis0 = fis0*ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs = dis0*(w(i,j+1,k,irho) - w(i,j,k,irho))
             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs = dis0*(w(i,j+1,k,ivx) - w(i,j,k,ivx))
             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs = dis0*(w(i,j+1,k,ivy) - w(i,j,k,ivy))
             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs = dis0*(w(i,j+1,k,ivz) - w(i,j,k,ivz))
             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy

             fs = dis0*(w(i,j+1,k,irhoE) - w(i,j,k,irhoE))
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

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half

             dis0 = fis0*ppor*(radK(i,j,k) + radK(i,j,k+1))

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs = dis0*(w(i,j,k+1,irho) - w(i,j,k,irho))
             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs = dis0*(w(i,j,k+1,ivx) - w(i,j,k,ivx))
             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs = dis0*(w(i,j,k+1,ivy) - w(i,j,k,ivy))
             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs = dis0*(w(i,j,k+1,ivz) - w(i,j,k,ivz))
             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy

             fs = dis0*(w(i,j,k+1,irhoE) - w(i,j,k,irhoE))
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

           enddo
         enddo
       enddo

       ! Replace rho times the total enthalpy by the total energy and
       ! store the velocities again instead of the momentum. As only
       ! the first halo cells are included, this must be done here again.

       do k=1,ke
         do j=1,je
           do i=1,ie
             rhoi = one/w(i,j,k,irho)
             w(i,j,k,ivx) = w(i,j,k,ivx)*rhoi
             w(i,j,k,ivy) = w(i,j,k,ivy)*rhoi
             w(i,j,k,ivz) = w(i,j,k,ivz)*rhoi

             w(i,j,k,irhoE) = w(i,j,k,irhoE) - p(i,j,k)
           enddo
         enddo
       enddo

       end subroutine inviscidDissFluxScalarCoarse

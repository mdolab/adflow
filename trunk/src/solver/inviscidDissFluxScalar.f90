!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxScalar.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-24-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidDissFluxScalar
!
!      ******************************************************************
!      *                                                                *
!      * inviscidDissFluxScalar computes the scalar artificial          *
!      * dissipation, see AIAA paper 81-1259, for a given block.        *
!      * Therefore it is assumed that the pointers in  blockPointers    *
!      * already point to the correct block.                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use constants
       use flowVarRefState
       use inputDiscretization
       use inputPhysics
       use iteration
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: dssMax = 0.25_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ind

       real(kind=realType) :: sslim, rhoi
       real(kind=realType) :: sfil, fis2, fis4
       real(kind=realType) :: ppor, rrad, dis2, dis4
       real(kind=realType) :: dss1, dss2, ddw, fs

       real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ss
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

           ! Copy the pressure in ss. Only fill the entries used in
           ! the discretization, i.e. ignore the corner halo's.

           do k=0,kb
             do j=2,jl
               do i=2,il
                 ss(i,j,k) = p(i,j,k)
               enddo
             enddo
           enddo

           do k=2,kl
             do j=2,jl
               ss(0, j,k) = p(0, j,k); ss(1, j,k) = p(1, j,k)
               ss(ie,j,k) = p(ie,j,k); ss(ib,j,k) = p(ib,j,k)
             enddo
           enddo

           do k=2,kl
             do i=2,il
               ss(i,0, k) = p(i,0, k); ss(i,1, k) = p(i,1, k)
               ss(i,je,k) = p(i,je,k); ss(i,jb,k) = p(i,jb,k)
             enddo
           enddo

         !===============================================================

         case (NSEquations, RANSEquations)

           ! Viscous case. Pressure switch is based on the entropy.
           ! Also set the value of sslim. To be fully consistent this
           ! must have the dimension of entropy and it is therefore
           ! set to a fraction of the free stream value.

           sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

           ! Store the entropy in ss. Only fill the entries used in
           ! the discretization, i.e. ignore the corner halo's.

           do k=0,kb
             do j=2,jl
               do i=2,il
                 ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
               enddo
             enddo
           enddo

           do k=2,kl
             do j=2,jl
               ss(0, j,k) = p(0, j,k)/(w(0, j,k,irho)**gamma(0, j,k))
               ss(1, j,k) = p(1, j,k)/(w(1, j,k,irho)**gamma(1, j,k))
               ss(ie,j,k) = p(ie,j,k)/(w(ie,j,k,irho)**gamma(ie,j,k))
               ss(ib,j,k) = p(ib,j,k)/(w(ib,j,k,irho)**gamma(ib,j,k))
             enddo
           enddo

           do k=2,kl
             do i=2,il
               ss(i,0, k) = p(i,0, k)/(w(i,0, k,irho)**gamma(i,0, k))
               ss(i,1, k) = p(i,1, k)/(w(i,1, k,irho)**gamma(i,1, k))
               ss(i,je,k) = p(i,je,k)/(w(i,je,k,irho)**gamma(i,je,k))
               ss(i,jb,k) = p(i,jb,k)/(w(i,jb,k,irho)**gamma(i,jb,k))
             enddo
           enddo

       end select

       ! Set a couple of constants for the scheme.

       fis2 = rFil*vis2
       fis4 = rFil*vis4
       sfil = one - rFil

       ! Replace the total energy by rho times the total enthalpy.
       ! In this way the numerical solution is total enthalpy preserving
       ! for the steady Euler equations. Also replace the velocities by
       ! the momentum. Only done for the entries used in the
       ! discretization, i.e. ignore the corner halo's.

       do k=0,kb
         do j=2,jl
           do i=2,il
             w(i,j,k,ivx)   = w(i,j,k,irho)*w(i,j,k,ivx)
             w(i,j,k,ivy)   = w(i,j,k,irho)*w(i,j,k,ivy)
             w(i,j,k,ivz)   = w(i,j,k,irho)*w(i,j,k,ivz)
             w(i,j,k,irhoE) = w(i,j,k,irhoE) + p(i,j,k)
           enddo
         enddo
       enddo

       do k=2,kl
         do j=2,jl
           w(0,j,k,ivx)   = w(0,j,k,irho)*w(0,j,k,ivx)
           w(0,j,k,ivy)   = w(0,j,k,irho)*w(0,j,k,ivy)
           w(0,j,k,ivz)   = w(0,j,k,irho)*w(0,j,k,ivz)
           w(0,j,k,irhoE) = w(0,j,k,irhoE) + p(0,j,k)

           w(1,j,k,ivx)   = w(1,j,k,irho)*w(1,j,k,ivx)
           w(1,j,k,ivy)   = w(1,j,k,irho)*w(1,j,k,ivy)
           w(1,j,k,ivz)   = w(1,j,k,irho)*w(1,j,k,ivz)
           w(1,j,k,irhoE) = w(1,j,k,irhoE) + p(1,j,k)

           w(ie,j,k,ivx)   = w(ie,j,k,irho)*w(ie,j,k,ivx)
           w(ie,j,k,ivy)   = w(ie,j,k,irho)*w(ie,j,k,ivy)
           w(ie,j,k,ivz)   = w(ie,j,k,irho)*w(ie,j,k,ivz)
           w(ie,j,k,irhoE) = w(ie,j,k,irhoE) + p(ie,j,k)

           w(ib,j,k,ivx)   = w(ib,j,k,irho)*w(ib,j,k,ivx)
           w(ib,j,k,ivy)   = w(ib,j,k,irho)*w(ib,j,k,ivy)
           w(ib,j,k,ivz)   = w(ib,j,k,irho)*w(ib,j,k,ivz)
           w(ib,j,k,irhoE) = w(ib,j,k,irhoE) + p(ib,j,k)
         enddo
       enddo

       do k=2,kl
         do i=2,il
           w(i,0,k,ivx)   = w(i,0,k,irho)*w(i,0,k,ivx)
           w(i,0,k,ivy)   = w(i,0,k,irho)*w(i,0,k,ivy)
           w(i,0,k,ivz)   = w(i,0,k,irho)*w(i,0,k,ivz)
           w(i,0,k,irhoE) = w(i,0,k,irhoE) + p(i,0,k)

           w(i,1,k,ivx)   = w(i,1,k,irho)*w(i,1,k,ivx)
           w(i,1,k,ivy)   = w(i,1,k,irho)*w(i,1,k,ivy)
           w(i,1,k,ivz)   = w(i,1,k,irho)*w(i,1,k,ivz)
           w(i,1,k,irhoE) = w(i,1,k,irhoE) + p(i,1,k)

           w(i,je,k,ivx)   = w(i,je,k,irho)*w(i,je,k,ivx)
           w(i,je,k,ivy)   = w(i,je,k,irho)*w(i,je,k,ivy)
           w(i,je,k,ivz)   = w(i,je,k,irho)*w(i,je,k,ivz)
           w(i,je,k,irhoE) = w(i,je,k,irhoE) + p(i,je,k)

           w(i,jb,k,ivx)   = w(i,jb,k,irho)*w(i,jb,k,ivx)
           w(i,jb,k,ivy)   = w(i,jb,k,irho)*w(i,jb,k,ivy)
           w(i,jb,k,ivz)   = w(i,jb,k,irho)*w(i,jb,k,ivz)
           w(i,jb,k,irhoE) = w(i,jb,k,irhoE) + p(i,jb,k)
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

           ! Compute the pressure sensor in the first cell, which
           ! is a halo cell.

           dss1 = abs((ss(2,j,k) - two*ss(1,j,k) + ss(0,j,k)) &
                /     (ss(2,j,k) + two*ss(1,j,k) + ss(0,j,k) + sslim))

           ! Loop in i-direction.

           do i=1,il

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((ss(i+2,j,k) - two*ss(i+1,j,k) + ss(i,j,k)) &
                  /     (ss(i+2,j,k) + two*ss(i+1,j,k) + ss(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw &
                 - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw)

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ind = indFamilyI(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyI(i,j,k)*fs

             ! X-momentum.

             ddw = w(i+1,j,k,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw &
                 - dis4*(w(i+2,j,k,ivx) - w(i-1,j,k,ivx) - three*ddw)

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i+1,j,k,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw &
                 - dis4*(w(i+2,j,k,ivy) - w(i-1,j,k,ivy) - three*ddw)

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i+1,j,k,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw &
                 - dis4*(w(i+2,j,k,ivz) - w(i-1,j,k,ivz) - three*ddw)

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw &
                 - dis4*(w(i+2,j,k,irhoE) - w(i-1,j,k,irhoE) - three*ddw)

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

             ! Set dss1 to dss2 for the next face.

             dss1 = dss2

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
         do i=2,il

           ! Compute the pressure sensor in the first cell, which
           ! is a halo cell.

           dss1 = abs((ss(i,2,k) - two*ss(i,1,k) + ss(i,0,k)) &
                /     (ss(i,2,k) + two*ss(i,1,k) + ss(i,0,k) + sslim))

           ! Loop in j-direction.

           do j=1,jl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((ss(i,j+2,k) - two*ss(i,j+1,k) + ss(i,j,k)) &
                  /     (ss(i,j+2,k) + two*ss(i,j+1,k) + ss(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw &
                 - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw)

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ind = indFamilyJ(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyJ(i,j,k)*fs

             ! X-momentum.

             ddw = w(i,j+1,k,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw &
                 - dis4*(w(i,j+2,k,ivx) - w(i,j-1,k,ivx) - three*ddw)

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j+1,k,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw &
                 - dis4*(w(i,j+2,k,ivy) - w(i,j-1,k,ivy) - three*ddw)

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j+1,k,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw &
                 - dis4*(w(i,j+2,k,ivz) - w(i,j-1,k,ivz) - three*ddw)

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw &
                 - dis4*(w(i,j+2,k,irhoE) - w(i,j-1,k,irhoE) - three*ddw)

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

             ! Set dss1 to dss2 for the next face.

             dss1 = dss2

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
       do j=2,jl
         do i=2,il

           ! Compute the pressure sensor in the first cell, which
           ! is a halo cell.

           dss1 = abs((ss(i,j,2) - two*ss(i,j,1) + ss(i,j,0)) &
                /     (ss(i,j,2) + two*ss(i,j,1) + ss(i,j,0) + sslim))

           ! Loop in k-direction.

           do k=1,kl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((ss(i,j,k+2) - two*ss(i,j,k+1) + ss(i,j,k)) &
                  /     (ss(i,j,k+2) + two*ss(i,j,k+1) + ss(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw &
                 - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw)

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ind = indFamilyK(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyK(i,j,k)*fs

             ! X-momentum.

             ddw = w(i,j,k+1,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw &
                 - dis4*(w(i,j,k+2,ivx) - w(i,j,k-1,ivx) - three*ddw)

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j,k+1,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw &
                 - dis4*(w(i,j,k+2,ivy) - w(i,j,k-1,ivy) - three*ddw)

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j,k+1,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw &
                 - dis4*(w(i,j,k+2,ivz) - w(i,j,k-1,ivz) - three*ddw)

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw &
                 - dis4*(w(i,j,k+2,irhoE) - w(i,j,k-1,irhoE) - three*ddw)

             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

             ! Set dss1 to dss2 for the next face.

             dss1 = dss2

           enddo
         enddo
       enddo

       ! Replace rho times the total enthalpy by the total energy and
       ! store the velocities again instead of the momentum. Only for
       ! those entries that have been altered, i.e. ignore the
       ! corner halo's.

       do k=0,kb
         do j=2,jl
           do i=2,il
             rhoi           = one/w(i,j,k,irho)
             w(i,j,k,ivx)   = w(i,j,k,ivx)*rhoi
             w(i,j,k,ivy)   = w(i,j,k,ivy)*rhoi
             w(i,j,k,ivz)   = w(i,j,k,ivz)*rhoi
             w(i,j,k,irhoE) = w(i,j,k,irhoE) - p(i,j,k)
           enddo
         enddo
       enddo

       do k=2,kl
         do j=2,jl
           rhoi           = one/w(0,j,k,irho)
           w(0,j,k,ivx)   = w(0,j,k,ivx)*rhoi
           w(0,j,k,ivy)   = w(0,j,k,ivy)*rhoi
           w(0,j,k,ivz)   = w(0,j,k,ivz)*rhoi
           w(0,j,k,irhoE) = w(0,j,k,irhoE) - p(0,j,k)

           rhoi           = one/w(1,j,k,irho)
           w(1,j,k,ivx)   = w(1,j,k,ivx)*rhoi
           w(1,j,k,ivy)   = w(1,j,k,ivy)*rhoi
           w(1,j,k,ivz)   = w(1,j,k,ivz)*rhoi
           w(1,j,k,irhoE) = w(1,j,k,irhoE) - p(1,j,k)

           rhoi            = one/w(ie,j,k,irho)
           w(ie,j,k,ivx)   = w(ie,j,k,ivx)*rhoi
           w(ie,j,k,ivy)   = w(ie,j,k,ivy)*rhoi
           w(ie,j,k,ivz)   = w(ie,j,k,ivz)*rhoi
           w(ie,j,k,irhoE) = w(ie,j,k,irhoE) - p(ie,j,k)

           rhoi            = one/w(ib,j,k,irho)
           w(ib,j,k,ivx)   = w(ib,j,k,ivx)*rhoi
           w(ib,j,k,ivy)   = w(ib,j,k,ivy)*rhoi
           w(ib,j,k,ivz)   = w(ib,j,k,ivz)*rhoi
           w(ib,j,k,irhoE) = w(ib,j,k,irhoE) - p(ib,j,k)
         enddo
       enddo

       do k=2,kl
         do i=2,il
           rhoi           = one/w(i,0,k,irho)
           w(i,0,k,ivx)   = w(i,0,k,ivx)*rhoi
           w(i,0,k,ivy)   = w(i,0,k,ivy)*rhoi
           w(i,0,k,ivz)   = w(i,0,k,ivz)*rhoi
           w(i,0,k,irhoE) = w(i,0,k,irhoE) - p(i,0,k)

           rhoi           = one/w(i,1,k,irho)
           w(i,1,k,ivx)   = w(i,1,k,ivx)*rhoi
           w(i,1,k,ivy)   = w(i,1,k,ivy)*rhoi
           w(i,1,k,ivz)   = w(i,1,k,ivz)*rhoi
           w(i,1,k,irhoE) = w(i,1,k,irhoE) - p(i,1,k)

           rhoi            = one/w(i,je,k,irho)
           w(i,je,k,ivx)   = w(i,je,k,ivx)*rhoi
           w(i,je,k,ivy)   = w(i,je,k,ivy)*rhoi
           w(i,je,k,ivz)   = w(i,je,k,ivz)*rhoi
           w(i,je,k,irhoE) = w(i,je,k,irhoE) - p(i,je,k)

           rhoi            = one/w(i,jb,k,irho)
           w(i,jb,k,ivx)   = w(i,jb,k,ivx)*rhoi
           w(i,jb,k,ivy)   = w(i,jb,k,ivy)*rhoi
           w(i,jb,k,ivz)   = w(i,jb,k,ivz)*rhoi
           w(i,jb,k,irhoE) = w(i,jb,k,irhoE) - p(i,jb,k)
         enddo
       enddo

       end subroutine inviscidDissFluxScalar

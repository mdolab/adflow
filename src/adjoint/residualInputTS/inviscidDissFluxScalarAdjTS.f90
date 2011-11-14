!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxScalarAdj.f90                   *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 06-10-2009                                      *
!      * Last modified: 06-10-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidDissFluxScalarAdjTS(wAdj,  pAdj,  dwadj,   &
                                            radIAdj,radJAdj,radKAdj, &
                                            iCell, jCell, kCell,nn,level,sps)
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
       use inputTimeSpectral !nTimeIntervalsSpectral
       use inputDiscretization
       implicit none

!
!      Subroutine arguments
!
       integer(kind=intType) :: iCell, jCell, kCell,nn,level,sps

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
                                                      intent(inout) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),    &
                                                      intent(in) :: pAdj
       real(kind=realType), dimension(-1:1,-1:1,-1:1,nTimeIntervalsSpectral) :: radIAdj,radJAdj,radKAdj
       real(kind=realType), dimension(nw,nTimeIntervalsSpectral), intent(inout) :: dwadj

!
!      Local parameter.
!
       real(kind=realType), parameter :: dssMax = 0.25_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ind
       integer(kind=intType) :: ii, jj, kk

       real(kind=realType) :: sslim, rhoi
       real(kind=realType) :: sfil, fis2, fis4
       real(kind=realType) :: ppor, rrad, dis2, dis4
       real(kind=realType) :: dss1, dss2, ddw, fs
       real(kind=realType) :: fact

       !real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ss
       real(kind=realType), dimension(-2:2,-2:2,-2:2) :: ss
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
           !do we need to ignore the corners in the ADjoint?... leave in for now...
           do k=-2,2!0,kb
             do j=-2,2!2,jl
               do i=-2,2!2,il
                 ss(i,j,k) = pAdj(i,j,k,sps)
               enddo
             enddo
           enddo

!!$           do k=2,kl
!!$             do j=2,jl
!!$               ss(0, j,k) = p(0, j,k); ss(1, j,k) = p(1, j,k)
!!$               ss(ie,j,k) = p(ie,j,k); ss(ib,j,k) = p(ib,j,k)
!!$             enddo
!!$           enddo
!!$
!!$           do k=2,kl
!!$             do i=2,il
!!$               ss(i,0, k) = p(i,0, k); ss(i,1, k) = p(i,1, k)
!!$               ss(i,je,k) = p(i,je,k); ss(i,jb,k) = p(i,jb,k)
!!$             enddo
!!$           enddo

         !===============================================================

         case (NSEquations, RANSEquations)

            print *,'NSEquations and RANSEquations not yet supported'
            stop
!!$
           ! Viscous case. Pressure switch is based on the entropy.
           ! Also set the value of sslim. To be fully consistent this
           ! must have the dimension of entropy and it is therefore
           ! set to a fraction of the free stream value.

 !!          sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

           ! Store the entropy in ss. Only fill the entries used in
           ! the discretization, i.e. ignore the corner halo's.
!!$
!!$           do k=0,kb
!!$             do j=2,jl
!!$               do i=2,il
!!$                 ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
!!$               enddo
!!$             enddo
!!$           enddo

   !        do k = -2,2
   !           do j = -2,2
   !              do i = -2,2
   !                 ss(i,j,k) = pAdj(i,j,k,sps)/(wAdj(i,j,k,irho,sps)**gamma(iCell+i,jCell+j,kCell+k))
   !              enddo
   !           enddo
   !        enddo
!!$
!!$           do k=2,kl
!!$             do j=2,jl
!!$               ss(0, j,k) = p(0, j,k)/(w(0, j,k,irho)**gamma(0, j,k))
!!$               ss(1, j,k) = p(1, j,k)/(w(1, j,k,irho)**gamma(1, j,k))
!!$               ss(ie,j,k) = p(ie,j,k)/(w(ie,j,k,irho)**gamma(ie,j,k))
!!$               ss(ib,j,k) = p(ib,j,k)/(w(ib,j,k,irho)**gamma(ib,j,k))
!!$             enddo
!!$           enddo
!!$
!!$           do k=2,kl
!!$             do i=2,il
!!$               ss(i,0, k) = p(i,0, k)/(w(i,0, k,irho)**gamma(i,0, k))
!!$               ss(i,1, k) = p(i,1, k)/(w(i,1, k,irho)**gamma(i,1, k))
!!$               ss(i,je,k) = p(i,je,k)/(w(i,je,k,irho)**gamma(i,je,k))
!!$               ss(i,jb,k) = p(i,jb,k)/(w(i,jb,k,irho)**gamma(i,jb,k))
!!$             enddo
!!$           enddo

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

       do k=-2,2!0,kb
         do j=-2,2!2,jl
           do i=-2,2!2,il
             wAdj(i,j,k,ivx,sps)   = wAdj(i,j,k,irho,sps)*wAdj(i,j,k,ivx,sps)
             wAdj(i,j,k,ivy,sps)   = wAdj(i,j,k,irho,sps)*wAdj(i,j,k,ivy,sps)
             wAdj(i,j,k,ivz,sps)   = wAdj(i,j,k,irho,sps)*wAdj(i,j,k,ivz,sps)
             wAdj(i,j,k,irhoE,sps) = wAdj(i,j,k,irhoE,sps) + pAdj(i,j,k,sps)
           enddo
         enddo
       enddo

!!$       do k=2,kl
!!$         do j=2,jl
!!$           w(0,j,k,ivx)   = w(0,j,k,irho)*w(0,j,k,ivx)
!!$           w(0,j,k,ivy)   = w(0,j,k,irho)*w(0,j,k,ivy)
!!$           w(0,j,k,ivz)   = w(0,j,k,irho)*w(0,j,k,ivz)
!!$           w(0,j,k,irhoE) = w(0,j,k,irhoE) + p(0,j,k)
!!$
!!$           w(1,j,k,ivx)   = w(1,j,k,irho)*w(1,j,k,ivx)
!!$           w(1,j,k,ivy)   = w(1,j,k,irho)*w(1,j,k,ivy)
!!$           w(1,j,k,ivz)   = w(1,j,k,irho)*w(1,j,k,ivz)
!!$           w(1,j,k,irhoE) = w(1,j,k,irhoE) + p(1,j,k)
!!$
!!$           w(ie,j,k,ivx)   = w(ie,j,k,irho)*w(ie,j,k,ivx)
!!$           w(ie,j,k,ivy)   = w(ie,j,k,irho)*w(ie,j,k,ivy)
!!$           w(ie,j,k,ivz)   = w(ie,j,k,irho)*w(ie,j,k,ivz)
!!$           w(ie,j,k,irhoE) = w(ie,j,k,irhoE) + p(ie,j,k)
!!$
!!$           w(ib,j,k,ivx)   = w(ib,j,k,irho)*w(ib,j,k,ivx)
!!$           w(ib,j,k,ivy)   = w(ib,j,k,irho)*w(ib,j,k,ivy)
!!$           w(ib,j,k,ivz)   = w(ib,j,k,irho)*w(ib,j,k,ivz)
!!$           w(ib,j,k,irhoE) = w(ib,j,k,irhoE) + p(ib,j,k)
!!$         enddo
!!$       enddo

!!$       do k=2,kl
!!$         do i=2,il
!!$           w(i,0,k,ivx)   = w(i,0,k,irho)*w(i,0,k,ivx)
!!$           w(i,0,k,ivy)   = w(i,0,k,irho)*w(i,0,k,ivy)
!!$           w(i,0,k,ivz)   = w(i,0,k,irho)*w(i,0,k,ivz)
!!$           w(i,0,k,irhoE) = w(i,0,k,irhoE) + p(i,0,k)
!!$
!!$           w(i,1,k,ivx)   = w(i,1,k,irho)*w(i,1,k,ivx)
!!$           w(i,1,k,ivy)   = w(i,1,k,irho)*w(i,1,k,ivy)
!!$           w(i,1,k,ivz)   = w(i,1,k,irho)*w(i,1,k,ivz)
!!$           w(i,1,k,irhoE) = w(i,1,k,irhoE) + p(i,1,k)
!!$
!!$           w(i,je,k,ivx)   = w(i,je,k,irho)*w(i,je,k,ivx)
!!$           w(i,je,k,ivy)   = w(i,je,k,irho)*w(i,je,k,ivy)
!!$           w(i,je,k,ivz)   = w(i,je,k,irho)*w(i,je,k,ivz)
!!$           w(i,je,k,irhoE) = w(i,je,k,irhoE) + p(i,je,k)
!!$
!!$           w(i,jb,k,ivx)   = w(i,jb,k,irho)*w(i,jb,k,ivx)
!!$           w(i,jb,k,ivy)   = w(i,jb,k,irho)*w(i,jb,k,ivy)
!!$           w(i,jb,k,ivz)   = w(i,jb,k,irho)*w(i,jb,k,ivz)
!!$           w(i,jb,k,irhoE) = w(i,jb,k,irhoE) + p(i,jb,k)
!!$         enddo
!!$       enddo


!Following method in the upwind scheme, take the residual onto dwAdj instead 
!of a separate fw. If it needs to be switched back, fw is dissiptive, dw is
!inviscid...
!!$       ! Initialize the dissipative residual to a certain times,
!!$       ! possibly zero, the previously stored value. Owned cells
!!$       ! only, because the halo values do not matter.
!!$
!!$!       do k=2,kl
!!$!         do j=2,jl
!!$!           do i=2,il
!!$             fw(irho)  = sfil*fw(irho)
!!$             fw(imx)   = sfil*fw(imx)
!!$             fw(imy)   = sfil*fw(imy)
!!$             fw(imz)   = sfil*fw(imz)
!!$             fw(irhoE) = sfil*fw(irhoE)
!!$ !          enddo
!!$ !        enddo
!!$ !      enddo
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the i-direction.                         *
!      *                                                                *
!      ******************************************************************
! 
             !set some indices for use later
       i    = iCell-1; j = jCell; k = kCell

       fact = one
       
       !do k=2,kl
       !  do j=2,jl

           ! Compute the pressure sensor in the first cell, which
           ! is a halo cell.

           !!dss1 = abs((ss(2,j,k) - two*ss(1,j,k) + ss(0,j,k)) &
           !!     /     (ss(2,j,k) + two*ss(1,j,k) + ss(0,j,k) + sslim))

           ! Loop in i-direction.
           do ii=-1,0
              !do i=1,il
              dss1 = abs((ss(ii+1,0,0) - two*ss(ii,0,0) + ss(ii-1,0,0)) &
                   /     (ss(ii+1,0,0) + two*ss(ii,0,0) + ss(ii-1,0,0) + sslim))
             ! print *,'dss1',dss1
             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((ss(ii+2,0,0) - two*ss(ii+1,0,0) + ss(ii,0,0)) &
                  /     (ss(ii+2,0,0) + two*ss(ii+1,0,0) + ss(ii,0,0) + sslim))
             !print *,'dss2',dss2
             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             !rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))
             rrad = ppor*(radIAdj(ii,0,0,sps) + radIAdj(ii+1,0,0,sps))
             !print *,'radI',radIAdj(ii,0,0),radI(icell+ii,jcell,kcell),icell,jcell,kcell,radIAdj(ii+1,0,0),radI(icell+ii+1,jcell,kcell)
               
             !lumped Dissipation for preconditioner
             if(lumpedDiss)then
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad
                dis4 = 0.0
             else
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
                !dis4 = dim(fis4*rrad, dis2)
                if ((fis4*rrad- dis2)>0.0)then
                   dis4 =fis4*rrad- dis2
                else
                   dis4 =0.0
                endif
             endif
             !print *,'dis2,4',dis2,dis4

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = wAdj(ii+1,0,0,irho,sps) - wAdj(ii,0,0,irho,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(ii+2,0,0,irho,sps) - wAdj(ii-1,0,0,irho,sps) - three*ddw)

             !fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             !fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs
             dwadj(irho,sps) = dwadj(irho,sps) + fact*fs

             ind = indFamilyI(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyI(i,j,k)*fs

             ! X-momentum.

             ddw = wAdj(ii+1,0,0,ivx,sps) - wAdj(ii,0,0,ivx,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(ii+2,0,0,ivx,sps) - wAdj(ii-1,0,0,ivx,sps) - three*ddw)

             !fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             !fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs
             dwadj(imx,sps) = dwadj(imx,sps) + fact*fs

             ! Y-momentum.

             ddw = wAdj(ii+1,0,0,ivy,sps) - wAdj(ii,0,0,ivy,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(ii+2,0,0,ivy,sps) - wAdj(ii-1,0,0,ivy,sps) - three*ddw)

             !fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             !fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs
             dwadj(imy,sps) = dwadj(imy,sps) + fact*fs

             ! Z-momentum.

             ddw = wAdj(ii+1,0,0,ivz,sps) - wAdj(ii,0,0,ivz,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(ii+2,0,0,ivz,sps) - wAdj(ii-1,0,0,ivz,sps) - three*ddw)

             !fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             !fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs
             dwadj(imz,sps) = dwadj(imz,sps) + fact*fs

             ! Energy.

             ddw = wAdj(ii+1,0,0,irhoE,sps) - wAdj(ii,0,0,irhoE,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(ii+2,0,0,irhoE,sps) - wAdj(ii-1,0,0,irhoE,sps) - three*ddw)

             !fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             !fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
             dwadj(irhoE,sps) = dwadj(irhoE,sps) + fact*fs

             ! Update i and set fact to 1 for the second face.

             i    = i + 1
             fact = -one
             !!! Set dss1 to dss2 for the next face.
             !!
             !!dss1 = dss2

           enddo
!         enddo
!       enddo

!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the j-direction.                         *
!      *                                                                *
!      ******************************************************************
!
      ! do k=2,kl
      !   do i=2,il

         i    = iCell; j = jCell-1; k = kCell
         fact = one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do jj=-1,0
!!           ! Compute the pressure sensor in the first cell, which
!!           ! is a halo cell.!
!!
           dss1 = abs((ss(0,jj+1,0) - two*ss(0,jj,0) + ss(0,jj-1,0)) &
                /     (ss(0,jj+1,0) + two*ss(0,jj,0) + ss(0,jj-1,0) + sslim))

!           ! Loop in j-direction.
!
!           do j=1,jl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((ss(0,jj+2,0) - two*ss(0,jj+1,0) + ss(0,jj,0)) &
                  /     (ss(0,jj+2,0) + two*ss(0,jj+1,0) + ss(0,jj,0) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             !rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))
             rrad = ppor*(radJAdj(0,jj,0,sps) + radJAdj(0,jj+1,0,sps))
   
             !lumped Dissipation for preconditioner
             if(lumpedDiss)then
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad
                dis4 = 0.0
             else
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
                !dis4 = dim(fis4*rrad, dis2)
                if ((fis4*rrad- dis2)>0.0)then
                   dis4 =fis4*rrad- dis2
                else
                   dis4 =0.0
                endif
             endif
             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = wAdj(0,jj+1,0,irho,sps) - wAdj(0,jj,0,irho,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,jj+2,0,irho,sps) - wAdj(0,jj-1,0,irho,sps) - three*ddw)

             !fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             !fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs
             dwadj(irho,sps) = dwadj(irho,sps) + fact*fs
             
             ind = indFamilyJ(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyJ(i,j,k)*fs

             ! X-momentum.

             ddw = wAdj(0,jj+1,0,ivx,sps) - wAdj(0,jj,0,ivx,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,jj+2,0,ivx,sps) - wAdj(0,jj-1,0,ivx,sps) - three*ddw)

             !fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             !fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs
             dwadj(imx,sps) = dwadj(imx,sps) + fact*fs
             
             ! Y-momentum.

             ddw = wAdj(0,jj+1,0,ivy,sps) - wAdj(0,jj,0,ivy,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,jj+2,0,ivy,sps) - wAdj(0,jj-1,0,ivy,sps) - three*ddw)

             !fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             !fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs
             dwadj(imy,sps) = dwadj(imy,sps) + fact*fs
             
             ! Z-momentum.

             ddw = wAdj(0,jj+1,0,ivz,sps) - wAdj(0,jj,0,ivz,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,jj+2,0,ivz,sps) - wAdj(0,jj-1,0,ivz,sps) - three*ddw)

             !fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             !fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs
             dwadj(imz,sps) = dwadj(imz,sps) + fact*fs
             
             ! Energy.

             ddw = wAdj(0,jj+1,0,irhoE,sps) - wAdj(0,jj,0,irhoE,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,jj+2,0,irhoE,sps) - wAdj(0,jj-1,0,irhoE,sps) - three*ddw)

             !fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             !fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
             dwadj(irhoE,sps) = dwadj(irhoE,sps) + fact*fs
             
             ! Update j and set fact to 1 for the second face.
             
             j    = j + 1
             fact = -one
             !!! Set dss1 to dss2 for the next face.
             !!
             !!dss1 = dss2

           enddo
        !enddo
       !enddo
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the k-direction.                         *
!      *                                                                *
!      ******************************************************************
!    

       ! Fluxes in k-direction.

       i    = iCell; j = jCell; k = kCell-1
       fact = one
!       do j=2,jl
!         do i=2,il
       ! Loop over the two faces which contribute to the residual of
       ! the cell considered.

       do kk=-1,0
           ! Compute the pressure sensor in the first cell, which
           ! is a halo cell.

           dss1 = abs((ss(0,0,kk+1) - two*ss(0,0,kk) + ss(0,0,kk-1)) &
                /     (ss(0,0,kk+1) + two*ss(0,0,kk) + ss(0,0,kk-1) + sslim))

           !!! Loop in k-direction.
           !!
           !!do k=1,kl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

           dss2 = abs((ss(0,0,kk+2) - two*ss(0,0,kk+1) + ss(0,0,kk)) &
                  /     (ss(0,0,kk+2) + two*ss(0,0,kk+1) + ss(0,0,kk) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             !rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))
             rrad = ppor*(radKAdj(0,0,kk,sps) + radKAdj(0,0,kk+1,sps))
             
             !lumped Dissipation for preconditioner
             if(lumpedDiss)then
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad
                dis4 = 0.0
             else
                dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))
                !dis4 = dim(fis4*rrad, dis2)
                if ((fis4*rrad- dis2)>0.0)then
                   dis4 =fis4*rrad- dis2
                else
                   dis4 =0.0
                endif
             endif

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = wAdj(0,0,kk+1,irho,sps) - wAdj(0,0,kk,irho,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,0,kk+2,irho,sps) - wAdj(0,0,kk-1,irho,sps) - three*ddw)

             !fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             !fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs
             dwadj(irho,sps) = dwadj(irho,sps) + fact*fs
             
             ind = indFamilyK(i,j,k)
             massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                                                  - factFamilyK(i,j,k)*fs

             ! X-momentum.

             ddw = wAdj(0,0,kk+1,ivx,sps) - wAdj(0,0,kk,ivx,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,0,kk+2,ivx,sps) - wAdj(0,0,kk-1,ivx,sps) - three*ddw)

             !fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             !fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs
             dwadj(imx,sps) = dwadj(imx,sps) + fact*fs

             ! Y-momentum.

             ddw = wAdj(0,0,kk+1,ivy,sps) - wAdj(0,0,kk,ivy,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,0,kk+2,ivy,sps) - wAdj(0,0,kk-1,ivy,sps) - three*ddw)

             !fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             !fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs
             dwadj(imy,sps) = dwadj(imy,sps) + fact*fs

             ! Z-momentum.

             ddw = wAdj(0,0,kk+1,ivz,sps) - wAdj(0,0,kk,ivz,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,0,kk+2,ivz,sps) - wAdj(0,0,kk-1,ivz,sps) - three*ddw)

             !fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             !fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs
             dwadj(imz,sps) = dwadj(imz,sps) + fact*fs
             
             ! Energy.

             ddw = wAdj(0,0,kk+1,irhoE,sps) - wAdj(0,0,kk,irhoE,sps)
             fs  = dis2*ddw &
                 - dis4*(wAdj(0,0,kk+2,irhoE,sps) - wAdj(0,0,kk-1,irhoE,sps) - three*ddw)

             !fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             !fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
             dwadj(irhoE,sps) = dwadj(irhoE,sps) + fact*fs

             ! Update k and set fact to 1 for the second face.
             
             k    = k + 1
             fact = -one
             !!! Set dss1 to dss2 for the next face.
             !!
             !dss1 = dss2

           enddo
         !enddo
       !enddo

       ! Replace rho times the total enthalpy by the total energy and
       ! store the velocities again instead of the momentum. Only for
       ! those entries that have been altered, i.e. ignore the
       ! corner halo's.

       ! Again, do we need to ignore the halo's? Simpler if we just 
       ! copy everything....
       do k=-2,2!0,kb
         do j=-2,2!2,jl
           do i=-2,2!2,il
             rhoi           = one/wAdj(i,j,k,irho,sps)
             wAdj(i,j,k,ivx,sps)   = wAdj(i,j,k,ivx,sps)*rhoi
             wAdj(i,j,k,ivy,sps)   = wAdj(i,j,k,ivy,sps)*rhoi
             wAdj(i,j,k,ivz,sps)   = wAdj(i,j,k,ivz,sps)*rhoi
             wAdj(i,j,k,irhoE,sps) = wAdj(i,j,k,irhoE,sps) - pAdj(i,j,k,sps)
           enddo
         enddo
       enddo

!!$       do k=2,kl
!!$         do j=2,jl
!!$           rhoi           = one/w(0,j,k,irho)
!!$           w(0,j,k,ivx)   = w(0,j,k,ivx)*rhoi
!!$           w(0,j,k,ivy)   = w(0,j,k,ivy)*rhoi
!!$           w(0,j,k,ivz)   = w(0,j,k,ivz)*rhoi
!!$           w(0,j,k,irhoE) = w(0,j,k,irhoE) - p(0,j,k)
!!$
!!$           rhoi           = one/w(1,j,k,irho)
!!$           w(1,j,k,ivx)   = w(1,j,k,ivx)*rhoi
!!$           w(1,j,k,ivy)   = w(1,j,k,ivy)*rhoi
!!$           w(1,j,k,ivz)   = w(1,j,k,ivz)*rhoi
!!$           w(1,j,k,irhoE) = w(1,j,k,irhoE) - p(1,j,k)
!!$
!!$           rhoi            = one/w(ie,j,k,irho)
!!$           w(ie,j,k,ivx)   = w(ie,j,k,ivx)*rhoi
!!$           w(ie,j,k,ivy)   = w(ie,j,k,ivy)*rhoi
!!$           w(ie,j,k,ivz)   = w(ie,j,k,ivz)*rhoi
!!$           w(ie,j,k,irhoE) = w(ie,j,k,irhoE) - p(ie,j,k)
!!$
!!$           rhoi            = one/w(ib,j,k,irho)
!!$           w(ib,j,k,ivx)   = w(ib,j,k,ivx)*rhoi
!!$           w(ib,j,k,ivy)   = w(ib,j,k,ivy)*rhoi
!!$           w(ib,j,k,ivz)   = w(ib,j,k,ivz)*rhoi
!!$           w(ib,j,k,irhoE) = w(ib,j,k,irhoE) - p(ib,j,k)
!!$         enddo
!!$       enddo
!!$
!!$       do k=2,kl
!!$         do i=2,il
!!$           rhoi           = one/w(i,0,k,irho)
!!$           w(i,0,k,ivx)   = w(i,0,k,ivx)*rhoi
!!$           w(i,0,k,ivy)   = w(i,0,k,ivy)*rhoi
!!$           w(i,0,k,ivz)   = w(i,0,k,ivz)*rhoi
!!$           w(i,0,k,irhoE) = w(i,0,k,irhoE) - p(i,0,k)
!!$
!!$           rhoi           = one/w(i,1,k,irho)
!!$           w(i,1,k,ivx)   = w(i,1,k,ivx)*rhoi
!!$           w(i,1,k,ivy)   = w(i,1,k,ivy)*rhoi
!!$           w(i,1,k,ivz)   = w(i,1,k,ivz)*rhoi
!!$           w(i,1,k,irhoE) = w(i,1,k,irhoE) - p(i,1,k)
!!$
!!$           rhoi            = one/w(i,je,k,irho)
!!$           w(i,je,k,ivx)   = w(i,je,k,ivx)*rhoi
!!$           w(i,je,k,ivy)   = w(i,je,k,ivy)*rhoi
!!$           w(i,je,k,ivz)   = w(i,je,k,ivz)*rhoi
!!$           w(i,je,k,irhoE) = w(i,je,k,irhoE) - p(i,je,k)
!!$
!!$           rhoi            = one/w(i,jb,k,irho)
!!$           w(i,jb,k,ivx)   = w(i,jb,k,ivx)*rhoi
!!$           w(i,jb,k,ivy)   = w(i,jb,k,ivy)*rhoi
!!$           w(i,jb,k,ivz)   = w(i,jb,k,ivz)*rhoi
!!$           w(i,jb,k,irhoE) = w(i,jb,k,irhoE) - p(i,jb,k)
!!$         enddo
!!$       enddo

     end subroutine inviscidDissFluxScalarAdjTS

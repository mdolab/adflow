module fluxes

contains
  subroutine inviscidCentralFlux
    !
    !       inviscidCentralFlux computes the Euler fluxes using a central  
    !       discretization for a given block. Therefore it is assumed that 
    !       the pointers in block pointer already point to the correct     
    !       block on the correct multigrid level.                          
    !
    use blockPointers
    use cgnsGrid
    use constants
    use flowVarRefState
    use inputPhysics
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ind, ii
    real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
    real(kind=realType) :: pa, fs, sFace, vnp, vnm
    real(kind=realType) :: wwx, wwy, wwz, rvol

    continue
    !$AD CHECKPOINT-START
    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.
    sFace = zero
    !
    !       Advective fluxes in the i-direction.                           
    !
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
                ! Set the dot product of the grid velocity and the
                ! normal in i-direction for a moving face.

                if( addGridVelocities ) sFace = sFaceI(i,j,k)

                ! Compute the normal velocities of the left and right state.

                vnp = w(i+1,j,k,ivx)*sI(i,j,k,1) &
                     + w(i+1,j,k,ivy)*sI(i,j,k,2) &
                     + w(i+1,j,k,ivz)*sI(i,j,k,3)
                vnm = w(i,  j,k,ivx)*sI(i,j,k,1) &
                     + w(i,  j,k,ivy)*sI(i,j,k,2) &
                     + w(i,  j,k,ivz)*sI(i,j,k,3)
                ! Set the values of the porosities for this face.
                ! porVel defines the porosity w.r.t. velocity;
                ! porFlux defines the porosity w.r.t. the entire flux.
                ! The latter is only zero for a discontinuous block
                ! boundary that must be treated conservatively.
                ! The default value of porFlux is 0.5, such that the
                ! correct central flux is scattered to both cells.
                ! In case of a boundFlux the normal velocity is set
                ! to sFace.

                porVel  = one
                porFlux = half
                if(porI(i,j,k) == noFlux)    porFlux = zero
                if(porI(i,j,k) == boundFlux) then
                   porVel = zero
                   vnp    = sFace
                   vnm    = sFace
                endif

                ! Incorporate porFlux in porVel.

                porVel = porVel*porFlux

                ! Compute the normal velocities relative to the grid for
                ! the face as well as the mass fluxes.

                qsp = (vnp -sFace)*porVel
                qsm = (vnm -sFace)*porVel

                rqsp = qsp*w(i+1,j,k,irho)
                rqsm = qsm*w(i,  j,k,irho)

                ! Compute the sum of the pressure multiplied by porFlux.
                ! For the default value of porFlux, 0.5, this leads to
                ! the average pressure.

                pa = porFlux*(p(i+1,j,k) + p(i,j,k))

                ! Compute the fluxes and scatter them to the cells
                ! i,j,k and i+1,j,k. Store the density flux in the
                ! mass flow of the appropriate sliding mesh interface.

                fs = rqsp + rqsm
                dw(i+1,j,k,irho) = dw(i+1,j,k,irho) - fs
                dw(i,  j,k,irho) = dw(i,  j,k,irho) + fs
#ifndef USE_TAPENADE
                ind = indFamilyI(i,j,k)
                massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                     + factFamilyI(i,j,k)*fs
#endif
                fs = rqsp*w(i+1,j,k,ivx) + rqsm*w(i,j,k,ivx) &
                     + pa*sI(i,j,k,1)
                dw(i+1,j,k,imx) = dw(i+1,j,k,imx) - fs
                dw(i,  j,k,imx) = dw(i,  j,k,imx) + fs

                fs = rqsp*w(i+1,j,k,ivy) + rqsm*w(i,j,k,ivy) &
                     + pa*sI(i,j,k,2)
                dw(i+1,j,k,imy) = dw(i+1,j,k,imy) - fs
                dw(i,  j,k,imy) = dw(i,  j,k,imy) + fs

                fs = rqsp*w(i+1,j,k,ivz) + rqsm*w(i,j,k,ivz) &
                     + pa*sI(i,j,k,3)
                dw(i+1,j,k,imz) = dw(i+1,j,k,imz) - fs
                dw(i,  j,k,imz) = dw(i,  j,k,imz) + fs

                fs = qsp*w(i+1,j,k,irhoE) + qsm*w(i,j,k,irhoE) &
                     + porFlux*(vnp*p(i+1,j,k) + vnm*p(i,j,k))
                dw(i+1,j,k,irhoE) = dw(i+1,j,k,irhoE) - fs
                dw(i,  j,k,irhoE) = dw(i,  j,k,irhoE) + fs
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
    !       Advective fluxes in the j-direction.                           
    !
    continue
    !$AD CHECKPOINT-START
    sface = zero
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
                ! Set the dot product of the grid velocity and the
                ! normal in j-direction for a moving face.

                if( addGridVelocities ) sFace = sFaceJ(i,j,k)

                ! Compute the normal velocities of the left and right state.

                vnp = w(i,j+1,k,ivx)*sJ(i,j,k,1) &
                     + w(i,j+1,k,ivy)*sJ(i,j,k,2) &
                     + w(i,j+1,k,ivz)*sJ(i,j,k,3)
                vnm = w(i,j,  k,ivx)*sJ(i,j,k,1) &
                     + w(i,j,  k,ivy)*sJ(i,j,k,2) &
                     + w(i,j,  k,ivz)*sJ(i,j,k,3)

                ! Set the values of the porosities for this face.
                ! porVel defines the porosity w.r.t. velocity;
                ! porFlux defines the porosity w.r.t. the entire flux.
                ! The latter is only zero for a discontinuous block
                ! boundary that must be treated conservatively.
                ! The default value of porFlux is 0.5, such that the
                ! correct central flux is scattered to both cells.
                ! In case of a boundFlux the normal velocity is set
                ! to sFace.

                porVel  = one
                porFlux = half
                if(porJ(i,j,k) == noFlux)    porFlux = zero
                if(porJ(i,j,k) == boundFlux) then
                   porVel = zero
                   vnp    = sFace
                   vnm    = sFace
                endif

                ! Incorporate porFlux in porVel.

                porVel = porVel*porFlux

                ! Compute the normal velocities for the face as well as the
                ! mass fluxes.

                qsp = (vnp - sFace)*porVel
                qsm = (vnm - sFace)*porVel

                rqsp = qsp*w(i,j+1,k,irho)
                rqsm = qsm*w(i,j,  k,irho)

                ! Compute the sum of the pressure multiplied by porFlux.
                ! For the default value of porFlux, 0.5, this leads to
                ! the average pressure.

                pa = porFlux*(p(i,j+1,k) + p(i,j,k))

                ! Compute the fluxes and scatter them to the cells
                ! i,j,k and i,j+1,k. Store the density flux in the
                ! mass flow of the appropriate sliding mesh interface.

                fs = rqsp + rqsm
                dw(i,j+1,k,irho) = dw(i,j+1,k,irho) - fs
                dw(i,j,  k,irho) = dw(i,j,  k,irho) + fs
#ifndef USE_TAPENADE
                ind = indFamilyJ(i,j,k)
                massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                     + factFamilyJ(i,j,k)*fs
#endif
                fs = rqsp*w(i,j+1,k,ivx) + rqsm*w(i,j,k,ivx) &
                     + pa*sJ(i,j,k,1)
                dw(i,j+1,k,imx) = dw(i,j+1,k,imx) - fs
                dw(i,j,  k,imx) = dw(i,j,  k,imx) + fs

                fs = rqsp*w(i,j+1,k,ivy) + rqsm*w(i,j,k,ivy) &
                     + pa*sJ(i,j,k,2)
                dw(i,j+1,k,imy) = dw(i,j+1,k,imy) - fs
                dw(i,j,  k,imy) = dw(i,j,  k,imy) + fs

                fs = rqsp*w(i,j+1,k,ivz) + rqsm*w(i,j,k,ivz) &
                     + pa*sJ(i,j,k,3)
                dw(i,j+1,k,imz) = dw(i,j+1,k,imz) - fs
                dw(i,j,  k,imz) = dw(i,j,  k,imz) + fs

                fs = qsp*w(i,j+1,k,irhoE) + qsm*w(i,j,k,irhoE) &
                     + porFlux*(vnp*p(i,j+1,k) + vnm*p(i,j,k))
                dw(i,j+1,k,irhoE) = dw(i,j+1,k,irhoE) - fs
                dw(i,j,  k,irhoE) = dw(i,j,  k,irhoE) + fs

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
    !       Advective fluxes in the k-direction.                           
    continue
    !$AD CHECKPOINT-START
    sface = zero
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
                ! Set the dot product of the grid velocity and the
                ! normal in k-direction for a moving face.

                if( addGridVelocities ) sFace = sFaceK(i,j,k)

                ! Compute the normal velocities of the left and right state.

                vnp = w(i,j,k+1,ivx)*sK(i,j,k,1) &
                     + w(i,j,k+1,ivy)*sK(i,j,k,2) &
                     + w(i,j,k+1,ivz)*sK(i,j,k,3)
                vnm = w(i,j,k,  ivx)*sK(i,j,k,1) &
                     + w(i,j,k,  ivy)*sK(i,j,k,2) &
                     + w(i,j,k,  ivz)*sK(i,j,k,3)

                ! Set the values of the porosities for this face.
                ! porVel defines the porosity w.r.t. velocity;
                ! porFlux defines the porosity w.r.t. the entire flux.
                ! The latter is only zero for a discontinuous block
                ! block boundary that must be treated conservatively.
                ! The default value of porFlux is 0.5, such that the
                ! correct central flux is scattered to both cells.
                ! In case of a boundFlux the normal velocity is set
                ! to sFace.

                porVel  = one
                porFlux = half

                if(porK(i,j,k) == noFlux)    porFlux = zero
                if(porK(i,j,k) == boundFlux) then
                   porVel = zero
                   vnp    = sFace
                   vnm    = sFace
                endif

                ! Incorporate porFlux in porVel.

                porVel = porVel*porFlux

                ! Compute the normal velocities for the face as well as the
                ! mass fluxes.

                qsp = (vnp - sFace)*porVel
                qsm = (vnm - sFace)*porVel

                rqsp = qsp*w(i,j,k+1,irho)
                rqsm = qsm*w(i,j,k,  irho)

                ! Compute the sum of the pressure multiplied by porFlux.
                ! For the default value of porFlux, 0.5, this leads to
                ! the average pressure.

                pa = porFlux*(p(i,j,k+1) + p(i,j,k))

                ! Compute the fluxes and scatter them to the cells
                ! i,j,k and i,j,k+1. Store the density flux in the
                ! mass flow of the appropriate sliding mesh interface.

                fs = rqsp + rqsm
                dw(i,j,k+1,irho) = dw(i,j,k+1,irho) - fs
                dw(i,j,k,  irho) = dw(i,j,k,  irho) + fs
#ifndef USE_TAPENADE
                ind = indFamilyK(i,j,k)
                massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                     + factFamilyK(i,j,k)*fs
#endif
                fs = rqsp*w(i,j,k+1,ivx) + rqsm*w(i,j,k,ivx) &
                     + pa*sK(i,j,k,1)
                dw(i,j,k+1,imx) = dw(i,j,k+1,imx) - fs
                dw(i,j,k,  imx) = dw(i,j,k,  imx) + fs

                fs = rqsp*w(i,j,k+1,ivy) + rqsm*w(i,j,k,ivy) &
                     + pa*sK(i,j,k,2)
                dw(i,j,k+1,imy) = dw(i,j,k+1,imy) - fs
                dw(i,j,k,  imy) = dw(i,j,k,  imy) + fs

                fs = rqsp*w(i,j,k+1,ivz) + rqsm*w(i,j,k,ivz) &
                     + pa*sK(i,j,k,3)
                dw(i,j,k+1,imz) = dw(i,j,k+1,imz) - fs
                dw(i,j,k,  imz) = dw(i,j,k,  imz) + fs

                fs = qsp*w(i,j,k+1,irhoE) + qsm*w(i,j,k,irhoE) &
                     + porFlux*(vnp*p(i,j,k+1) + vnm*p(i,j,k))
                dw(i,j,k+1,irhoE) = dw(i,j,k+1,irhoE) - fs
                dw(i,j,k,  irhoE) = dw(i,j,k,  irhoE) + fs
#ifdef TAPENADE_FAST
             end do
#else
          enddo
       enddo
    enddo
#endif   
    continue
    !$AD CHECKPOINT-END

    ! Add the rotational source terms for a moving block in a
    ! steady state computation. These source terms account for the
    ! centrifugal acceleration and the coriolis term. However, as
    ! the the equations are solved in the inertial frame and not
    ! in the moving frame, the form is different than what you
    ! normally find in a text book.

    continue
    !$AD CHECKPOINT-START
    rotation: if(blockIsMoving .and. equationMode == steady) then

       ! Compute the three nonDimensional angular velocities.

       wwx = timeRef*cgnsDoms(nbkGlobal)%rotRate(1)
       wwy = timeRef*cgnsDoms(nbkGlobal)%rotRate(2)
       wwz = timeRef*cgnsDoms(nbkGlobal)%rotRate(3)

       ! Loop over the internal cells of this block to compute the
       ! rotational terms for the momentum equations.
       !$AD II-LOOP
       do ii=0,nx*ny*nz-1
          i = mod(ii, nx) + 2
          j = mod(ii/nx, ny) + 2
          k = ii/(nx*ny) + 2
          rvol = w(i,j,k,irho)*vol(i,j,k)

          dw(i,j,k,imx) = dw(i,j,k,imx) &
               + rvol*(wwy*w(i,j,k,ivz) - wwz*w(i,j,k,ivy))
          dw(i,j,k,imy) = dw(i,j,k,imy) &
               + rvol*(wwz*w(i,j,k,ivx) - wwx*w(i,j,k,ivz))
          dw(i,j,k,imz) = dw(i,j,k,imz) &
               + rvol*(wwx*w(i,j,k,ivy) - wwy*w(i,j,k,ivx))
       enddo

    endif rotation

    !$AD CHECKPOINT-END
    continue
  end subroutine inviscidCentralFlux

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

  subroutine inviscidUpwindFlux(fineGrid)
    !
    !       inviscidUpwindFlux computes the artificial dissipation part of 
    !       the Euler fluxes by means of an approximate solution of the 1D 
    !       Riemann problem on the face. For first order schemes,          
    !       fineGrid == .false., the states in the cells are assumed to    
    !       be constant; for the second order schemes on the fine grid a   
    !       nonlinear reconstruction of the left and right state is done   
    !       for which several options exist.                               
    !       It is assumed that the pointers in blockPointers already       
    !       point to the correct block.                                    
    !
    use constants
    use blockPointers, only : il, jl, kl, ie, je, ke, ib, jb, kb, w, p, &
         porI, porJ, porK, fw, gamma, si, sj, sk, &
         indFamilyI, indFamilyJ, indFamilyK, spectralSol, addGridVelocities, &
         sFaceI, sfaceJ, sFacek, rotMatrixI, rotMatrixJ, rotMatrixK, &
         factFamilyI, factFamilyJ, factFamilyK 
    use flowVarRefState, only : kPresent, nw, nwf, rgas, tref
    use inputDiscretization, only: limiter, lumpedDiss, precond, riemann, &
         riemannCoarse, orderTurb, kappaCoef
    use inputPhysics, only : equations
    use iteration, only : rFil, currentLevel, groundLevel
    use cgnsGrid, only: massFlowFamilyDiss
    use utils, only : getCorrectForK, terminate
    use flowUtils, only : eTot
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: fineGrid
    !
    !      Local variables.
    !
    integer(kind=porType) :: por

    integer(kind=intType) :: nwInt
    integer(kind=intType) :: i, j, k, ind
    integer(kind=intType) :: limUsed, riemannUsed

    real(kind=realType) :: sx, sy, sz, omk, opk, sFil, gammaFace
    real(kind=realType) :: factMinmod, sFace

    real(kind=realType), dimension(nw)  :: left, right
    real(kind=realType), dimension(nw)  :: du1, du2, du3
    real(kind=realType), dimension(nwf) :: flux

    logical :: firstOrderK, correctForK, rotationalPeriodic
    !
    ! Check if rFil == 0. If so, the dissipative flux needs not to
    ! be computed.

    if(abs(rFil) < thresholdReal) return

    ! Check if the formulation for rotational periodic problems
    ! must be used.

    if( associated(rotMatrixI) ) then
       rotationalPeriodic = .true.
    else
       rotationalPeriodic = .false.
    endif

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value. Owned cells
    ! only, because the halo values do not matter.

    sFil = one - rFil

    do k=2,kl
       do j=2,jl
          do i=2,il
             fw(i,j,k,irho)  = sFil*fw(i,j,k,irho)
             fw(i,j,k,imx)   = sFil*fw(i,j,k,imx)
             fw(i,j,k,imy)   = sFil*fw(i,j,k,imy)
             fw(i,j,k,imz)   = sFil*fw(i,j,k,imz)
             fw(i,j,k,irhoE) = sFil*fw(i,j,k,irhoE)
          enddo
       enddo
    enddo

    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    correctForK = getCorrectForK()

    ! Compute the factor used in the minmod limiter.

    factMinmod = (three-kappaCoef) &
         / max(1.e-10_realType, one-kappaCoef)

    ! Determine the limiter scheme to be used. On the fine grid the
    ! user specified scheme is used; on the coarse grid a first order
    ! scheme is computed.

    limUsed = firstOrder
    if( fineGrid ) limUsed = limiter

    ! Lumped diss is true for doing approx PC
    if(lumpedDiss) then
       limUsed = firstOrder
    end if

    ! Determine the riemann solver which must be used.

    riemannUsed = riemannCoarse
    if( fineGrid ) riemannUsed = riemann

    ! Store 1-kappa and 1+kappa a bit easier and multiply it by 0.25.

    omk = fourth*(one - kappaCoef)
    opk = fourth*(one + kappaCoef)

    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.

    sFace = zero

    ! Set the number of variables to be interpolated depending
    ! whether or not a k-equation is present. If a k-equation is
    ! present also set the logical firstOrderK. This indicates
    ! whether or not only a first order approximation is to be used
    ! for the turbulent kinetic energy.

    if( correctForK ) then
       if(orderTurb == firstOrder) then
          nwInt       = nwf
          firstOrderK = .true.
       else
          nwInt       = itu1
          firstOrderK = .false.
       endif
    else
       nwInt       = nwf
       firstOrderK = .false.
    endif
    !
    !       Flux computation. A distinction is made between first and      
    !       second order schemes to avoid the overhead for the first order 
    !       scheme.                                                        
    !
    orderTest: if(limUsed == firstOrder) then
       !
       !         First order reconstruction. The states in the cells are      
       !         constant. The left and right states are constructed easily.  
       !
       ! Fluxes in the i-direction.

       do k=2,kl
          do j=2,jl
             do i=1,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
                por = porI(i,j,k)
                if( addGridVelocities ) sFace = sFaceI(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i+1,j,k,irho)
                right(ivx)   = w(i+1,j,k,ivx)
                right(ivy)   = w(i+1,j,k,ivy)
                right(ivz)   = w(i+1,j,k,ivz)
                right(irhoE) = p(i+1,j,k)
                if( correctForK ) right(itu1) = w(i+1,j,k,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
                fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
                fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
                fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
                fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyI(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyI(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

       ! Fluxes in j-direction.

       do k=2,kl
          do j=1,jl
             do i=2,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
                por = porJ(i,j,k)
                if( addGridVelocities ) sFace = sFaceJ(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i,j+1,k,irho)
                right(ivx)   = w(i,j+1,k,ivx)
                right(ivy)   = w(i,j+1,k,ivy)
                right(ivz)   = w(i,j+1,k,ivz)
                right(irhoE) = p(i,j+1,k)
                if( correctForK ) right(itu1) = w(i,j+1,k,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
                fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
                fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
                fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
                fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyJ(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyJ(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

       ! Fluxes in k-direction.

       do k=1,kl
          do j=2,jl
             do i=2,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
                por = porK(i,j,k)
                if( addGridVelocities ) sFace = sFaceK(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i,j,k+1,irho)
                right(ivx)   = w(i,j,k+1,ivx)
                right(ivy)   = w(i,j,k+1,ivy)
                right(ivz)   = w(i,j,k+1,ivz)
                right(irhoE) = p(i,j,k+1)
                if( correctForK ) right(itu1) = w(i,j,k+1,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
                fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
                fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
                fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
                fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyK(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyK(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

       !      ==================================================================

    else orderTest

       !      ==================================================================
       !
       !         Second order reconstruction of the left and right state.     
       !         The three differences used in the, possibly nonlinear,       
       !         interpolation are constructed here; the actual left and      
       !         right states, or at least the differences from the first     
       !         order interpolation, are computed in the subroutine          
       !         leftRightState.                                              
       !
       ! Fluxes in the i-direction.

       do k=2,kl
          do j=2,jl
             do i=1,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i-1,j,k,irho)
                du2(irho) = w(i+1,j,k,irho) - w(i,  j,k,irho)
                du3(irho) = w(i+2,j,k,irho) - w(i+1,j,k,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i-1,j,k,ivx)
                du2(ivx) = w(i+1,j,k,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i+2,j,k,ivx) - w(i+1,j,k,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i-1,j,k,ivy)
                du2(ivy) = w(i+1,j,k,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i+2,j,k,ivy) - w(i+1,j,k,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i-1,j,k,ivz)
                du2(ivz) = w(i+1,j,k,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i+2,j,k,ivz) - w(i+1,j,k,ivz)

                du1(irhoE) = p(i,  j,k) - p(i-1,j,k)
                du2(irhoE) = p(i+1,j,k) - p(i,  j,k)
                du3(irhoE) = p(i+2,j,k) - p(i+1,j,k)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i-1,j,k,itu1)
                   du2(itu1) = w(i+1,j,k,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i+2,j,k,itu1) - w(i+1,j,k,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixI, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i+1,j,k,irho)
                right(ivx)   = right(ivx)   + w(i+1,j,k,ivx)
                right(ivy)   = right(ivy)   + w(i+1,j,k,ivy)
                right(ivz)   = right(ivz)   + w(i+1,j,k,ivz)
                right(irhoE) = right(irhoE) + p(i+1,j,k)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i+1,j,k,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
                por = porI(i,j,k)
                if( addGridVelocities ) sFace = sFaceI(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
                fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
                fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
                fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
                fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyI(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyI(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

       ! Fluxes in the j-direction.

       do k=2,kl
          do j=1,jl
             do i=2,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i,j-1,k,irho)
                du2(irho) = w(i,j+1,k,irho) - w(i,  j,k,irho)
                du3(irho) = w(i,j+2,k,irho) - w(i,j+1,k,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i,j-1,k,ivx)
                du2(ivx) = w(i,j+1,k,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i,j+2,k,ivx) - w(i,j+1,k,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i,j-1,k,ivy)
                du2(ivy) = w(i,j+1,k,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i,j+2,k,ivy) - w(i,j+1,k,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i,j-1,k,ivz)
                du2(ivz) = w(i,j+1,k,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i,j+2,k,ivz) - w(i,j+1,k,ivz)

                du1(irhoE) = p(i,  j,k) - p(i,j-1,k)
                du2(irhoE) = p(i,j+1,k) - p(i,  j,k)
                du3(irhoE) = p(i,j+2,k) - p(i,j+1,k)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i,j-1,k,itu1)
                   du2(itu1) = w(i,j+1,k,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i,j+2,k,itu1) - w(i,j+1,k,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixJ, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i,j+1,k,irho)
                right(ivx)   = right(ivx)   + w(i,j+1,k,ivx)
                right(ivy)   = right(ivy)   + w(i,j+1,k,ivy)
                right(ivz)   = right(ivz)   + w(i,j+1,k,ivz)
                right(irhoE) = right(irhoE) + p(i,j+1,k)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i,j+1,k,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
                por = porJ(i,j,k)
                if( addGridVelocities ) sFace = sFaceJ(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
                fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
                fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
                fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
                fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyJ(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyJ(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

       ! Fluxes in the k-direction.

       do k=1,kl
          do j=2,jl
             do i=2,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i,j,k-1,irho)
                du2(irho) = w(i,j,k+1,irho) - w(i,  j,k,irho)
                du3(irho) = w(i,j,k+2,irho) - w(i,j,k+1,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i,j,k-1,ivx)
                du2(ivx) = w(i,j,k+1,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i,j,k+2,ivx) - w(i,j,k+1,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i,j,k-1,ivy)
                du2(ivy) = w(i,j,k+1,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i,j,k+2,ivy) - w(i,j,k+1,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i,j,k-1,ivz)
                du2(ivz) = w(i,j,k+1,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i,j,k+2,ivz) - w(i,j,k+1,ivz)

                du1(irhoE) = p(i,  j,k) - p(i,j,k-1)
                du2(irhoE) = p(i,j,k+1) - p(i,  j,k)
                du3(irhoE) = p(i,j,k+2) - p(i,j,k+1)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i,j,k-1,itu1)
                   du2(itu1) = w(i,j,k+1,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i,j,k+2,itu1) - w(i,j,k+1,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixK, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i,j,k+1,irho)
                right(ivx)   = right(ivx)   + w(i,j,k+1,ivx)
                right(ivy)   = right(ivy)   + w(i,j,k+1,ivy)
                right(ivz)   = right(ivz)   + w(i,j,k+1,ivz)
                right(irhoE) = right(irhoE) + p(i,j,k+1)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i,j,k+1,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
                por = porK(i,j,k)
                if( addGridVelocities ) sFace = sFaceK(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
                fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
                fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
                fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
                fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

                ! Store the density flux in the mass flow of the
                ! appropriate sliding mesh interface.
#ifndef USE_TAPENADE
                ind = indFamilyK(i,j,k)
                massFlowFamilyDiss(ind,spectralSol) =       &
                     massFlowFamilyDiss(ind,spectralSol) &
                     + factFamilyK(i,j,k)*flux(irho)
#endif
             enddo
          enddo
       enddo

    endif orderTest

    !      ==================================================================

  contains

    subroutine leftRightState(du1, du2, du3, rotMatrix, left, right)
      !
      !         leftRightState computes the differences in the left and      
      !         right state compared to the first order interpolation. For a 
      !         monotonic second order discretization the interpolations     
      !         need to be nonlinear. The linear second order scheme can be  
      !         stable (depending on the value of kappa), but it will have   
      !         oscillations near discontinuities.                           
      !
      implicit none
      !
      !        Local parameter.
      !
      real(kind=realType), parameter :: epsLim = 1.e-10_realType
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(:), intent(inout) :: du1, du2, du3
      real(kind=realType), dimension(:), intent(out)   :: left, right

      real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrix
      !
      !        Local variables.
      !
      integer(kind=intType) :: l

      real(kind=realType) :: rl1, rl2, rr1, rr2, tmp, dvx, dvy, dvz

      real(kind=realType), dimension(3,3) :: rot

      ! Check if the velocity components should be transformed to
      ! the cylindrical frame.

      if( rotationalPeriodic ) then

         ! Store the rotation matrix a bit easier. Note that the i,j,k
         ! come from the main subroutine.

         rot(1,1) = rotMatrix(i,j,k,1,1)
         rot(1,2) = rotMatrix(i,j,k,1,2)
         rot(1,3) = rotMatrix(i,j,k,1,3)

         rot(2,1) = rotMatrix(i,j,k,2,1)
         rot(2,2) = rotMatrix(i,j,k,2,2)
         rot(2,3) = rotMatrix(i,j,k,2,3)

         rot(3,1) = rotMatrix(i,j,k,3,1)
         rot(3,2) = rotMatrix(i,j,k,3,2)
         rot(3,3) = rotMatrix(i,j,k,3,3)

         ! Apply the transformation to the velocity components
         ! of du1, du2 and du3.

         dvx = du1(ivx); dvy = du1(ivy); dvz = du1(ivz)
         du1(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du1(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du1(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

         dvx = du2(ivx); dvy = du2(ivy); dvz = du2(ivz)
         du2(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du2(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du2(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

         dvx = du3(ivx); dvy = du3(ivy); dvz = du3(ivz)
         du3(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du3(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du3(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

      endif

      ! Determine the limiter used.

      select case (limUsed)

      case (noLimiter)

         ! Linear interpolation; no limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt
            left(l)  =  omk*du1(l) + opk*du2(l)
            right(l) = -omk*du3(l) - opk*du2(l)
         enddo

         !          ==============================================================

      case (vanAlbeda)

         ! Nonlinear interpolation using the van albeda limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt

            ! Compute the limiter argument rl1, rl2, rr1 and rr2.
            ! Note the cut off to 0.0.

            tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
            rl1 = max(zero, &
                 du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
            rl2 = max(zero,du1(l)*tmp)

            rr1 = max(zero,du3(l)*tmp)
            rr2 = max(zero, &
                 du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

            ! Compute the corresponding limiter values.

            rl1 = rl1*(rl1 + one)/(rl1*rl1 + one)
            rl2 = rl2*(rl2 + one)/(rl2*rl2 + one)
            rr1 = rr1*(rr1 + one)/(rr1*rr1 + one)
            rr2 = rr2*(rr2 + one)/(rr2*rr2 + one)

            ! Compute the nonlinear corrections to the first order
            ! scheme.

            left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
            right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

         enddo

         !          ==============================================================

      case (minmod)

         ! Nonlinear interpolation using the minmod limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt

            ! Compute the limiter argument rl1, rl2, rr1 and rr2.
            ! Note the cut off to 0.0.

            tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
            rl1 = max(zero, &
                 du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
            rl2 = max(zero,du1(l)*tmp)

            rr1 = max(zero,du3(l)*tmp)
            rr2 = max(zero, &
                 du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

            ! Compute the corresponding limiter values.

            rl1 = min(one, factMinmod*rl1)
            rl2 = min(one, factMinmod*rl2)
            rr1 = min(one, factMinmod*rr1)
            rr2 = min(one, factMinmod*rr2)

            ! Compute the nonlinear corrections to the first order
            ! scheme.

            left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
            right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

         enddo

      end select

      ! In case only a first order scheme must be used for the
      ! turbulent transport equations, set the correction for the
      ! turbulent kinetic energy to 0.

      if( firstOrderK ) then
         left(itu1)  = zero
         right(itu1) = zero
      endif

      ! For rotational periodic problems transform the velocity
      ! differences back to Cartesian again. Note that now the
      ! transpose of the rotation matrix must be used.

      if( rotationalPeriodic ) then

         ! Left state.

         dvx = left(ivx); dvy = left(ivy); dvz = left(ivz)
         left(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
         left(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
         left(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

         ! Right state.

         dvx = right(ivx); dvy = right(ivy); dvz = right(ivz)
         right(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
         right(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
         right(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

      endif

    end subroutine leftRightState

    !        ================================================================

    subroutine riemannFlux(left, right, flux)
      !
      !         riemannFlux computes the flux for the given face and left    
      !         and right states.                                            
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(*), intent(in)  :: left, right
      real(kind=realType), dimension(*), intent(out) :: flux
      !
      !        Local variables.
      !
      real(kind=realType) :: porFlux, rFace
      real(kind=realType) :: Etl, Etr, z1l, z1r, tmp
      real(kind=realType) :: dr, dru, drv, drw, drE, drk
      real(kind=realType) :: rAvg, uAvg, vAvg, wAvg, hAvg, kAvg
      real(kind=realType) :: alphaAvg, a2Avg, aAvg, unAvg
      real(kind=realType) :: ovaAvg, ova2Avg, area, Eta
      real(kind=realType) :: gm1, gm53
      real(kind=realType) :: lam1, lam2, lam3
      real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
      real(kind=realType), dimension(2) :: ktmp

      ! Set the porosity for the flux. The default value, 0.5*rFil, is
      ! a scaling factor where an rFil != 1 is taken into account.

      porFlux = half*rFil
      if(por == noFlux .or. por == boundFlux)    porFlux = zero

      ! Abbreviate some expressions in which gamma occurs.

      gm1  = gammaFace - one
      gm53 = gammaFace - five*third

      ! Determine which riemann solver must be solved.

      select case (riemannUsed)

      case (Roe)

         ! Determine the preconditioner used.

         select case (precond)

         case (noPrecond)

            ! No preconditioner used. Use the Roe scheme of the
            ! standard equations.

            ! Compute the square root of the left and right densities
            ! and the inverse of the sum.

            z1l = sqrt(left(irho))
            z1r = sqrt(right(irho))
            tmp = one/(z1l + z1r)

            ! Compute some variables depending whether or not a
            ! k-equation is present.

            if( correctForK ) then

               ! Store the left and right kinetic energy in ktmp,
               ! which is needed to compute the total energy.

               ktmp(1) = left(itu1)
               ktmp(2) = right(itu1)

               ! Store the difference of the turbulent kinetic energy
               ! per unit volume, i.e. the conserved variable.

               drk = right(irho)*right(itu1) - left(irho)*left(itu1)

               ! Compute the average turbulent energy per unit mass
               ! using Roe averages.

               kAvg = tmp*(z1l*left(itu1) + z1r*right(itu1))

            else

               ! Set the difference of the turbulent kinetic energy
               ! per unit volume and the averaged kinetic energy per
               ! unit mass to zero.

               drk  = 0.0
               kAvg = 0.0

            endif

            ! Compute the total energy of the left and right state.
            call etot(left(irho), left(ivx), left(ivy), left(ivz), &
                 left(irhoe), ktmp(1), Etl, correctForK)

            call etot(right(irho), right(ivx), right(ivy), right(ivz), &
                 right(irhoe), ktmp(2), Etr, correctForK)

            ! Compute the difference of the conservative mean
            ! flow variables.

            dr  = right(irho) - left(irho)
            dru = right(irho)*right(ivx) - left(irho)*left(ivx)
            drv = right(irho)*right(ivy) - left(irho)*left(ivy)
            drw = right(irho)*right(ivz) - left(irho)*left(ivz)
            drE = Etr - Etl

            ! Compute the Roe average variables, which can be
            ! computed directly from the average Roe vector.

            rAvg = fourth*(z1r + z1l)**2
            uAvg = tmp*(z1l*left(ivx) + z1r*right(ivx))
            vAvg = tmp*(z1l*left(ivy) + z1r*right(ivy))
            wAvg = tmp*(z1l*left(ivz) + z1r*right(ivz))
            hAvg = tmp*((Etl+left(irhoE)) /z1l &
                 +      (Etr+right(irhoE))/z1r)

            ! Compute the unit vector and store the area of the
            ! normal. Also compute the unit normal velocity of the face.

            area  = sqrt(sx**2 + sy**2 + sz**2)
            tmp   = one/max(1.e-25_realType,area)
            sx    = sx*tmp
            sy    = sy*tmp
            sz    = sz*tmp
            rFace = sFace*tmp

            ! Compute some dependent variables at the Roe
            ! average state.

            alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
            a2Avg    = abs(gm1*(hAvg - alphaAvg) - gm53*kAvg)
            aAvg     = sqrt(a2Avg)
            unAvg    = uAvg*sx + vAvg*sy + wAvg*sz

            ovaAvg  = one/aAvg
            ova2Avg = one/a2Avg

            ! Set for a boundary the normal velocity to rFace, the
            ! normal velocity of the boundary.

            if(por == boundFlux) unAvg = rFace

            ! Compute the coefficient eta for the entropy correction.
            ! At the moment a 1D entropy correction is used, which
            ! removes expansion shocks. Although it also reduces the
            ! carbuncle phenomenon, it does not remove it completely.
            ! In other to do that a multi-dimensional entropy fix is
            ! needed, see Sanders et. al, JCP, vol. 145, 1998,
            ! pp. 511 - 537. Although relatively easy to implement,
            ! an efficient implementation requires the storage of
            ! all the left and right states, which is rather
            ! expensive in terms of memory.

            eta = half*(abs((left(ivx) - right(ivx))*sx        &
                 +           (left(ivy) - right(ivy))*sy        &
                 +           (left(ivz) - right(ivz))*sz)       &
                 +       abs(sqrt(gammaFace*left(irhoE)/left(irho)) &
                 -           sqrt(gammaFace*right(irhoE)/right(irho))))

            ! Compute the absolute values of the three eigenvalues.

            lam1 = abs(unAvg - rFace + aAvg)
            lam2 = abs(unAvg - rFace - aAvg)
            lam3 = abs(unAvg - rFace)

            ! Apply the entropy correction to the eigenvalues.

            tmp = two*eta
            if(lam1 < tmp) lam1 = eta + fourth*lam1*lam1/eta
            if(lam2 < tmp) lam2 = eta + fourth*lam2*lam2/eta
            if(lam3 < tmp) lam3 = eta + fourth*lam3*lam3/eta

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
                 -      wAvg*drw + drE) - gm53*drk
            abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

            abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
            abv7 = abv2*abv4*ovaAvg  + abv3*abv5

            ! Compute the dissipation term, -|a| (wr - wl), which is
            ! multiplied by porFlux. Note that porFlux is either
            ! 0.0 or 0.5*rFil.

            flux(irho)  = -porFlux*(lam3*dr  + abv6)
            flux(imx)   = -porFlux*(lam3*dru + uAvg*abv6 &
                 + sx*abv7)
            flux(imy)   = -porFlux*(lam3*drv + vAvg*abv6 &
                 + sy*abv7)
            flux(imz)   = -porFlux*(lam3*drw + wAvg*abv6 &
                 + sz*abv7)
            flux(irhoE) = -porFlux*(lam3*drE + hAvg*abv6 &
                 + unAvg*abv7)

            !          tmp = max(lam1,lam2,lam3)

            !          flux(irho)  = -porFlux*(tmp*dr)
            !          flux(imx)   = -porFlux*(tmp*dru)
            !          flux(imy)   = -porFlux*(tmp*drv)
            !          flux(imz)   = -porFlux*(tmp*drw)
            !          flux(irhoE) = -porFlux*(tmp*drE)

         case (Turkel)
            call terminate(&
                 "riemannFlux",&
                 "Turkel preconditioner not implemented yet")

         case (ChoiMerkle)
            call terminate("riemannFlux",&
                 "choi merkle preconditioner not implemented yet")

         end select

      case (vanLeer)
         call terminate("riemannFlux", "van leer fvs not implemented yet")

      case (ausmdv)
         call terminate("riemannFlux","ausmdv fvs not implemented yet")

      end select

    end subroutine riemannFlux

  end subroutine inviscidUpwindFlux

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
#ifndef USE_TAPENADE
    use solverUtils, only : utauWf
#endif
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

  subroutine viscousFluxApprox
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

  subroutine inviscidDissFluxScalarApprox
    !
    !       inviscidDissFluxScalar computes the scalar artificial          
    !       dissipation, see AIAA paper 81-1259, for a given block.        
    !       Therefore it is assumed that the pointers in  blockPointers    
    !       already point to the correct block.                            
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
    real(kind=realType) :: ppor, rrad, dis2
    real(kind=realType) :: dss1, dss2, ddw, fs


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


       !===============================================================

    case (NSEquations, RANSEquations)

       ! Viscous case. Pressure switch is based on the entropy.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of entropy and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

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
    !       Dissipative fluxes in the i-direction.                         
    !
    do k=2,kl
       do j=2,jl

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dss1 = abs((shockSensor(2,j,k) - two*shockSensor(1,j,k) + shockSensor(0,j,k)) &
               /     (shockSensor(2,j,k) + two*shockSensor(1,j,k) + shockSensor(0,j,k) + sslim))

          ! Loop in i-direction.

          do i=1,il

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((shockSensor(i+2,j,k) - two*shockSensor(i+1,j,k) + shockSensor(i,j,k)) &
                  /     (shockSensor(i+2,j,k) + two*shockSensor(i+1,j,k) + shockSensor(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             ! Modification for FD Preconditioner Note: This lumping
             ! actually still results in a greater than 3 cell stencil
             ! in any direction. Since this seems to work slightly
             ! better than the dis2=sigma*fis4*rrad, we will just use
             ! a 5-cell stencil for doing the PC

             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i+1,j,k,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i+1,j,k,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i+1,j,k,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

             ! Set dss1 to dss2 for the next face.

             dss1 = dss2

          enddo
       enddo
    enddo
    !
    !       Dissipative fluxes in the j-direction.                         
    !
    do k=2,kl
       do i=2,il

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dss1 = abs((shockSensor(i,2,k) - two*shockSensor(i,1,k) + shockSensor(i,0,k)) &
               /     (shockSensor(i,2,k) + two*shockSensor(i,1,k) + shockSensor(i,0,k) + sslim))

          ! Loop in j-direction.

          do j=1,jl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((shockSensor(i,j+2,k) - two*shockSensor(i,j+1,k) + shockSensor(i,j,k)) &
                  /     (shockSensor(i,j+2,k) + two*shockSensor(i,j+1,k) + shockSensor(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             ! Modification for FD Preconditioner
             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j+1,k,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j+1,k,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j+1,k,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

             ! Set dss1 to dss2 for the next face.

             dss1 = dss2

          enddo
       enddo
    enddo
    !
    !       Dissipative fluxes in the k-direction.                         
    !
    do j=2,jl
       do i=2,il

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dss1 = abs((shockSensor(i,j,2) - two*shockSensor(i,j,1) + shockSensor(i,j,0)) &
               /     (shockSensor(i,j,2) + two*shockSensor(i,j,1) + shockSensor(i,j,0) + sslim))

          ! Loop in k-direction.

          do k=1,kl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dss2 = abs((shockSensor(i,j,k+2) - two*shockSensor(i,j,k+1) + shockSensor(i,j,k)) &
                  /     (shockSensor(i,j,k+2) + two*shockSensor(i,j,k+1) + shockSensor(i,j,k) + sslim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             ! Modification for FD Preconditioner
             dis2 = fis2*rrad*min(dssMax, max(dss1,dss2))+sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j,k+1,ivx) - w(i,j,k,ivx)
             fs  = dis2*ddw

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j,k+1,ivy) - w(i,j,k,ivy)
             fs  = dis2*ddw

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j,k+1,ivz) - w(i,j,k,ivz)
             fs  = dis2*ddw

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
             fs  = dis2*ddw

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

  end subroutine inviscidDissFluxScalarApprox

  subroutine inviscidDissFluxMatrixApprox
    !
    !       inviscidDissFluxMatrix computes the matrix artificial          
    !       dissipation term. Instead of the spectral radius, as used in   
    !       the scalar dissipation scheme, the absolute value of the flux  
    !       jacobian is used. This leads to a less diffusive and           
    !       consequently more accurate scheme. It is assumed that the      
    !       pointers in blockPointers already point to the correct block.  
    !
    use blockPointers
    use cgnsGrid
    use constants
    use flowVarRefState
    use inputDiscretization
    use inputPhysics
    use iteration
    use utils, only : getCorrectForK
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
    integer(kind=intType) :: i, j, k, ind

    real(kind=realType) :: plim, sface
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: gammaAvg, gm1, ovgm1, gm53
    real(kind=realType) :: ppor, rrad, dis2
    real(kind=realType) :: dp1, dp2, ddw, tmp, fs
    real(kind=realType) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
    real(kind=realType) :: uAvg, vAvg, wAvg, a2Avg, aAvg, hAvg
    real(kind=realType) :: alphaAvg, unAvg, ovaAvg, ova2Avg
    real(kind=realType) :: kAvg, lam1, lam2, lam3, area
    real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
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
    !       Dissipative fluxes in the i-direction.                         
    !
    do k=2,kl
       do j=2,jl

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dp1 = abs((shockSensor(2,j,k) - two*shockSensor(1,j,k) + shockSensor(0,j,k))        &
               /     (omega*(shockSensor(2,j,k) + two*shockSensor(1,j,k) + shockSensor(0,j,k)) &
               +      oneMinOmega*(abs(shockSensor(2,j,k)-shockSensor(1,j,k))        &
               +                   abs(shockSensor(1,j,k)-shockSensor(0,j,k))) + plim))

          ! Loop in i-direction.

          do i=1,il

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dp2 = abs((shockSensor(i+2,j,k) - two*shockSensor(i+1,j,k) + shockSensor(i,j,k))        &
                  /     (omega*(shockSensor(i+2,j,k) + two*shockSensor(i+1,j,k) + shockSensor(i,j,k)) &
                  +      oneMinOmega*(abs(shockSensor(i+2,j,k) - shockSensor(i+1,j,k))      &
                  +                   abs(shockSensor(i+1,j,k) - shockSensor(i,j,k))) + plim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dp1,dp2))+sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
                ddw = w(i+1,j,k,irho)*w(i+1,j,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw

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

             ! Set dp1 to dp2 for the next face.

             dp1 = dp2

          enddo
       enddo
    enddo
    !
    !       Dissipative fluxes in the j-direction.                         
    !
    do k=2,kl
       do i=2,il

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dp1 = abs((shockSensor(i,2,k) - two*shockSensor(i,1,k) + shockSensor(i,0,k))        &
               /     (omega*(shockSensor(i,2,k) + two*shockSensor(i,1,k) + shockSensor(i,0,k)) &
               +      oneMinOmega*(abs(shockSensor(i,2,k)-shockSensor(i,1,k))        &
               +                    abs(shockSensor(i,1,k)-shockSensor(i,0,k))) + plim))

          ! Loop in j-direction.

          do j=1,jl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dp2 = abs((shockSensor(i,j+2,k) - two*shockSensor(i,j+1,k) + shockSensor(i,j,k))        &
                  /     (omega*(shockSensor(i,j+2,k) + two*shockSensor(i,j+1,k) + shockSensor(i,j,k)) &
                  +      oneMinOmega*(abs(shockSensor(i,j+2,k) - shockSensor(i,j+1,k))      &
                  +                   abs(shockSensor(i,j+1,k) - shockSensor(i,j,k))) + plim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dp1,dp2))+sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
                ddw = w(i,j+1,k,irho)*w(i,j+1,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw

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

             ! Set dp1 to dp2 for the next face.

             dp1 = dp2

          enddo
       enddo
    enddo
    !
    !       Dissipative fluxes in the k-direction.                         
    !
    do j=2,jl
       do i=2,il

          ! Compute the pressure sensor in the first cell, which
          ! is a halo cell.

          dp1 = abs((shockSensor(i,j,2) - two*shockSensor(i,j,1) + shockSensor(i,j,0))        &
               /     (omega*(shockSensor(i,j,2) + two*shockSensor(i,j,1) + shockSensor(i,j,0)) &
               +      oneMinOmega*(abs(shockSensor(i,j,2)-shockSensor(i,j,1))        &
               +                   abs(shockSensor(i,j,1)-shockSensor(i,j,0))) + plim))

          ! Loop in k-direction.

          do k=1,kl

             ! Compute the pressure sensor in the cell to the right
             ! of the face.

             dp2 = abs((shockSensor(i,j,k+2) - two*shockSensor(i,j,k+1) + shockSensor(i,j,k))        &
                  /     (omega*(shockSensor(i,j,k+2) + two*shockSensor(i,j,k+1) + shockSensor(i,j,k)) &
                  +      oneMinOmega*(abs(shockSensor(i,j,k+2) - shockSensor(i,j,k+1))      &
                  +                   abs(shockSensor(i,j,k+1) - shockSensor(i,j,k))) + plim))

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dp1,dp2))+sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.

             if( correctForK ) then
                ddw = w(i,j,k+1,irho)*w(i,j,k+1,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw
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

             ! Set dp1 to dp2 for the next face.

             dp1 = dp2

          enddo
       enddo
    enddo

  end subroutine inviscidDissFluxMatrixApprox

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef  USE_TAPENADE

  subroutine inviscidDissFluxScalarCoarse
    !
    !       inviscidDissFluxScalarCoarse computes the coarse grid, i.e.    
    !       1st order, artificial dissipation flux for the scalar          
    !       dissipation scheme for a given block. Therefore it is assumed  
    !       that the pointers in blockPointers already point to the        
    !       correct block.                                                 
    !
    use constants
    use blockPointers, only : il, jl, kl, ie, je, ke, w, p, &
         porI, porJ, porK, fw, radI, radJ, radK, gamma
    use inputDiscretization, only: vis2Coarse
    use iteration, only : rFil

    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k

    real(kind=realType) :: sfil, fis0, dis0, ppor, fs, rhoi

    ! Check if rFil == 0. If so, the dissipative flux needs not to
    ! be computed.

    if(abs(rFil) < thresholdReal) return

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
    !       Dissipative fluxes in the i-direction.                         
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
    !       Dissipative fluxes in the j-direction.                         
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
    !       Dissipative fluxes in the k-direction.                         
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

  subroutine inviscidDissFluxMatrixCoarse
    !
    !       inviscidDissFluxMatrixCoarse computes the matrix artificial    
    !       dissipation term. Instead of the spectral radius, as used in   
    !       the scalar dissipation scheme, the absolute value of the flux  
    !       jacobian is used. This routine is used on the coarser grids in 
    !       the multigrid cycle and only computes the first order          
    !       dissipation term. It is assumed that the pointers in           
    !       blockPointers already point to the correct block.              
    !
    use constants
    use blockPointers, only : il, jl, kl, ie, je, ke, ib, jb, kb, w, p, &
         porI, porJ, porK, fw, gamma, si, sj, sk, &
         addGridVelocities, sFaceI, sfaceJ, sFacek
    use inputDiscretization, only: vis2Coarse
    use inputPhysics, only : equations
    use iteration, only : rFil
    use utils, only : getCorrectForK
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

    ! Check if rFil == 0. If so, the dissipative flux needs not to
    ! be computed.

    if(abs(rFil) < thresholdReal) return

    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.

    correctForK = getCorrectForK()

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
    !       Dissipative fluxes in the i-direction.                         
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
    !       Dissipative fluxes in the j-direction.                         
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
    !       Dissipative fluxes in the k-direction.                         
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
#endif
end module fluxes

module solverUtils
contains
#ifndef USE_TAPENADE
  subroutine timeStep(onlyRadii)
    !
    ! Shell function to call timeStep_block on all blocks
    !
    use constants
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine argument.
    !
    logical, intent(in) :: onlyRadii
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, nn

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, currentLevel, sps)

          call timeStep_Block(onlyRadii)

       end do domains

    end do spectralLoop

  end subroutine timeStep
#endif
  subroutine timeStep_block(onlyRadii)
    !
    !       timeStep computes the time step, or more precisely the time
    !       step divided by the volume per unit CFL, in the owned cells.
    !       However, for the artificial dissipation schemes, the spectral
    !       radIi in the halo's are needed. Therefore the loop is taken
    !       over the the first level of halo cells. The spectral radIi are
    !       stored and possibly modified for high aspect ratio cells.
    !
    use constants
    use blockPointers, only : ie, je, ke, il, jl, kl, w, p, rlv, rev, &
         radi, radj, radk, si, sj, sk, sFaceI, sfaceJ, sfaceK, dtl, gamma, vol, &
         addGridVelocities, sectionID
    use flowVarRefState, only : timeRef, eddyModel, gammaInf, pInfCorr, &
         viscous, rhoInf
    use inputDiscretization, only : adis, dirScaling, radiiNeededCoarse, &
         radiiNeededFine, precond
    use inputPhysics, only : equationMode
    use iteration, only : groundLevel, currentLevel
    use section, only : sections
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : terminate
    implicit none
    !
    !      Subroutine argument.
    !
    logical, intent(in) :: onlyRadii
    !
    !      Local parameters.
    !
    real(kind=realType), parameter :: b = 2.0_realType
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii

    real(kind=realType) :: plim, rlim, clim2
    real(kind=realType) :: uux, uuy, uuz, cc2, qsi, qsj, qsk, sx, sy, sz, rmu
    real(kind=realType) :: ri, rj, rk, rij, rjk, rki
    real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
    real(kind=realType) :: sFace, tmp

    logical :: radiiNeeded, doScaling

    ! Determine whether or not the spectral radii are needed for the
    ! flux computation.

    radiiNeeded = radiiNeededCoarse
    if(currentLevel <= groundLevel) radiiNeeded = radiiNeededFine

    ! Return immediately if only the spectral radii must be computed
    ! and these are not needed for the flux computation.

    if(onlyRadii .and. (.not. radiiNeeded)) return

    ! Set the value of plim. To be fully consistent this must have
    ! the dimension of a pressure. Therefore a fraction of pInfCorr
    ! is used. Idem for rlim; compute clim2 as well.

    plim  = 0.001_realType*pInfCorr
    rlim  = 0.001_realType*rhoInf
    clim2 = 0.000001_realType*gammaInf*pInfCorr/rhoInf

    doScaling = (dirScaling .and. currentLevel <= groundLevel)

    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.

    sFace = zero
    !
    !           Inviscid contribution, depending on the preconditioner.
    !           Compute the cell centered values of the spectral radii.
    !
    select case (precond)

    case (noPrecond)

       ! No preconditioner. Simply the standard spectral radius.
       ! Loop over the cells, including the first level halo.

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
                   ! Compute the velocities and speed of sound squared.

                   uux  = w(i,j,k,ivx)
                   uuy  = w(i,j,k,ivy)
                   uuz  = w(i,j,k,ivz)
                   cc2 = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
                   cc2 = max(cc2,clim2)

                   ! Set the dot product of the grid velocity and the
                   ! normal in i-direction for a moving face. To avoid
                   ! a number of multiplications by 0.5 simply the sum
                   ! is taken.

                   if( addGridVelocities ) &
                        sFace = sFaceI(i-1,j,k) + sFaceI(i,j,k)

                   ! Spectral radius in i-direction.

                   sx = si(i-1,j,k,1) + si(i,j,k,1)
                   sy = si(i-1,j,k,2) + si(i,j,k,2)
                   sz = si(i-1,j,k,3) + si(i,j,k,3)

                   qsi = uux*sx + uuy*sy + uuz*sz - sFace

                   ri = half*(abs(qsi) &
                        +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                   ! The grid velocity in j-direction.

                   if( addGridVelocities ) &
                        sFace = sFaceJ(i,j-1,k) + sFaceJ(i,j,k)

                   ! Spectral radius in j-direction.

                   sx = sj(i,j-1,k,1) + sj(i,j,k,1)
                   sy = sj(i,j-1,k,2) + sj(i,j,k,2)
                   sz = sj(i,j-1,k,3) + sj(i,j,k,3)

                   qsj = uux*sx + uuy*sy + uuz*sz - sFace

                   rj = half*(abs(qsj) &
                        +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                   ! The grid velocity in k-direction.

                   if( addGridVelocities ) &
                        sFace = sFaceK(i,j,k-1) + sFaceK(i,j,k)

                   ! Spectral radius in k-direction.

                   sx = sk(i,j,k-1,1) + sk(i,j,k,1)
                   sy = sk(i,j,k-1,2) + sk(i,j,k,2)
                   sz = sk(i,j,k-1,3) + sk(i,j,k,3)

                   qsk = uux*sx + uuy*sy + uuz*sz - sFace

                   rk = half*(abs(qsk) &
                        +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                   ! Compute the inviscid contribution to the time step.

                   if (.not. onlyRadii) dtl(i,j,k) = ri + rj + rk

                   !
                   !           Adapt the spectral radii if directional scaling must be
                   !           applied.
                   !
                   if(doScaling) then

                      ! Avoid division by zero by clipping radi, radJ and
                      ! radK.

                      ri = max(ri, eps)
                      rj = max(rj, eps)
                      rk = max(rk, eps)

                      ! Compute the scaling in the three coordinate
                      ! directions.

                      rij = (ri/rj)**adis
                      rjk = (rj/rk)**adis
                      rki = (rk/ri)**adis

                      ! Create the scaled versions of the aspect ratios.
                      ! Note that the multiplication is done with radi, radJ
                      ! and radK, such that the influence of the clipping
                      ! is negligible.

                      radi(i,j,k) = ri*(one + one/rij + rki)
                      radJ(i,j,k) = rj*(one + one/rjk + rij)
                      radK(i,j,k) = rk*(one + one/rki + rjk)
                   else
                      radi(i,j,k) = ri
                      radj(i,j,k) = rj
                      radk(i,j,k) = rk
                   end if
#ifdef TAPENADE_REVERSE
                end do
#else
             enddo
          enddo
       enddo
#endif

    case (Turkel)
       call terminate("timeStep","Turkel preconditioner not implemented yet")


    case (ChoiMerkle)
       call terminate("timeStep", &
            "choi merkle preconditioner not implemented yet")
    end select


    ! The rest of this file can be skipped if only the spectral
    ! radii need to be computed.
#ifndef USE_TAPENADE
    testRadiiOnly: if(.not. onlyRadii) then

       ! The viscous contribution, if needed.

       viscousTerm: if( viscous ) then

          ! Loop over the owned cell centers.

          do k=2,kl
             do j=2,jl
                do i=2,il

                   ! Compute the effective viscosity coefficient. The
                   ! factor 0.5 is a combination of two things. In the
                   ! standard central discretization of a second
                   ! derivative there is a factor 2 multiplying the
                   ! central node. However in the code below not the
                   ! average but the sum of the left and the right face
                   ! is taken and squared. This leads to a factor 4.
                   ! Combining both effects leads to 0.5. Furthermore,
                   ! it is divided by the volume and density to obtain
                   ! the correct dimensions and multiplied by the
                   ! non-dimensional factor factVis.

                   rmu = rlv(i,j,k)
                   if( eddyModel ) rmu = rmu + rev(i,j,k)
                   rmu = half*rmu/(w(i,j,k,irho)*vol(i,j,k))

                   ! Add the viscous contribution in i-direction to the
                   ! (inverse) of the time step.

                   sx = si(i,j,k,1) + si(i-1,j,k,1)
                   sy = si(i,j,k,2) + si(i-1,j,k,2)
                   sz = si(i,j,k,3) + si(i-1,j,k,3)

                   vsi        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsi

                   ! Add the viscous contribution in j-direction to the
                   ! (inverse) of the time step.

                   sx = sj(i,j,k,1) + sj(i,j-1,k,1)
                   sy = sj(i,j,k,2) + sj(i,j-1,k,2)
                   sz = sj(i,j,k,3) + sj(i,j-1,k,3)

                   vsj        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsj

                   ! Add the viscous contribution in k-direction to the
                   ! (inverse) of the time step.

                   sx = sk(i,j,k,1) + sk(i,j,k-1,1)
                   sy = sk(i,j,k,2) + sk(i,j,k-1,2)
                   sz = sk(i,j,k,3) + sk(i,j,k-1,3)

                   vsk        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsk

                enddo
             enddo
          enddo

       endif viscousTerm

       ! For the spectral mode an additional term term must be
       ! taken into account, which corresponds to the contribution
       ! of the highest frequency.

       if(equationMode == timeSpectral) then

          tmp = nTimeIntervalsSpectral*pi*timeRef &
               / sections(sectionID)%timePeriod

          ! Loop over the owned cell centers and add the term.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   dtl(i,j,k) = dtl(i,j,k) + tmp*vol(i,j,k)
                enddo
             enddo
          enddo

       endif

       ! Currently the inverse of dt/vol is stored in dtl. Invert
       ! this value such that the time step per unit cfl number is
       ! stored and correct in cases of high gradients.

       do k=2,kl
          do j=2,jl
             do i=2,il
                dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                     /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                     /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                     /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                rfl = one/(one + b*(dpi  +dpj  +dpk))

                dtl(i,j,k) = rfl/dtl(i,j,k)
             enddo
          enddo
       enddo

    endif testRadiiOnly
#endif
  end subroutine timeStep_block

#ifndef USE_TAPENADE
  subroutine gridVelocitiesFineLevel(useOldCoor, t, sps)
    !
    ! Shell function to call gridVelocitiesFineLevel on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor
    real(kind=realType), dimension(*), intent(in) :: t  !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Loop over the number of blocks.

    domains: do nn=1,nDom

       ! Set the pointers for this block.

       call setPointers(nn, groundLevel, sps)
       call gridVelocitiesFineLevel_block(useOldCoor, t, sps)

    end do domains

  end subroutine gridVelocitiesFineLevel
#endif


#ifndef USE_TAPENADE
  subroutine gridVelocitiesFineLevel_TS(sps)
    !
    ! Shell function to call gridVelocitiesFineLevel on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    integer(kind=intType) :: nn

    ! Loop over the number of blocks.

    domains: do nn=1,nDom

       ! Set the pointers for this block.

       call setPointers(nn, groundLevel, sps)
       call gridVelocitiesFineLevel_TS_block(nn, sps)

    end do domains

  end subroutine gridVelocitiesFineLevel_TS
#endif

  subroutine gridVelocitiesFineLevel_TS_block(nn, sps)

    use precision
    use constants
    use blockPointers!, only: nDom, ie, je, ke, x, s, sFaceI, sI, sFaceJ, sJ, sFaceK, sK
    use inputPhysics, only: machgrid, velDirFreestream
    use flowVarRefState, only: gammaInf, pInf, rhoInf
    !use partitioning, only: timePeriodSpectral
    use inputTimeSpectral, only: dscalar, nTimeIntervalsSpectral

    !      local variables
    integer(kind=intType), intent(in) :: nn, sps
    integer :: i, j, k, mm, ii, ie_l, je_l, ke_l                                                   ! index variables
    real(kind=realType) :: x_vc, y_vc, z_vc                                        ! cell volume center coords
    real(kind=realType) :: x_fc, y_fc, z_fc                                        ! cell face center coords
    real(kind=realType) :: aInf                                                              ! sound speed
    real(kind=realType) :: velxGrid, velyGrid, velzGrid                                      ! infinite speed

     


    ! get the grid free stream velocity
    aInf = sqrt(gammaInf*pInf/rhoInf)
    velxGrid = (aInf*machgrid)*(-velDirFreestream(1))
    velyGrid = (aInf*machgrid)*(-velDirFreestream(2))
    velzGrid = (aInf*machgrid)*(-velDirFreestream(3))

    ! get the temporal info (T ect.)
    !call timePeriodSpectral


!
!            ************************************************************
!            *                                                          *
!            * Grid velocities of the cell centers, including the       *
!            * 1st level halo cells.                                    *
!            *                                                          *
!            ************************************************************
!

    ! initialize with free stream velocity
    ie_l = flowDoms(nn, 1, sps)%ie
    je_l = flowDoms(nn, 1, sps)%je
    ke_l = flowDoms(nn, 1, sps)%ke  

    do k=1, ke_l
       do j=1, je_l
          do i=1, ie_l

             s(i, j, k, 1) = velxGrid
             s(i, j, k, 2) = velyGrid             
             s(i, j, k, 3) = velzGrid

          end do
       end do
    end do

    ! the velocity contributed from mesh deformation 
    do mm=1, nTimeIntervalsSpectral

       ie_l = flowDoms(nn, 1, mm)%ie
       je_l = flowDoms(nn, 1, mm)%je
       ke_l = flowDoms(nn, 1, mm)%ke
       
       do k=1,ke_l
          do j=1,je_l
             do i=1,ie_l

                x_vc = eighth*(flowDoms(nn, 1, mm)%x(i-1,j-1,k-1,1) + flowDoms(nn, 1, mm)%x(i,j-1,k-1,1) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k-1,1) + flowDoms(nn, 1, mm)%x(i,j,  k-1,1) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j-1,k,  1) + flowDoms(nn, 1, mm)%x(i,j-1,k,  1) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k,  1) + flowDoms(nn, 1, mm)%x(i,j,  k,  1))

                y_vc = eighth*(flowDoms(nn, 1, mm)%x(i-1,j-1,k-1,2) + flowDoms(nn, 1, mm)%x(i,j-1,k-1,2) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k-1,2) + flowDoms(nn, 1, mm)%x(i,j,  k-1,2) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j-1,k,  2) + flowDoms(nn, 1, mm)%x(i,j-1,k,  2) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k,  2) + flowDoms(nn, 1, mm)%x(i,j,  k,  2))

	              z_vc = eighth*(flowDoms(nn, 1, mm)%x(i-1,j-1,k-1,3) + flowDoms(nn, 1, mm)%x(i,j-1,k-1,3) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k-1,3) + flowDoms(nn, 1, mm)%x(i,j,  k-1,3) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j-1,k,  3) + flowDoms(nn, 1, mm)%x(i,j-1,k,  3) &
                     +         flowDoms(nn, 1, mm)%x(i-1,j,  k,  3) + flowDoms(nn, 1, mm)%x(i,j,  k,  3))


                s(i,j,k,1) = s(i,j,k,1) + dscalar(1, sps, mm)*x_vc
                s(i,j,k,2) = s(i,j,k,2) + dscalar(1, sps, mm)*y_vc
                s(i,j,k,3) = s(i,j,k,3) + dscalar(1, sps, mm)*z_vc

             end do
          end do
       end do

    end do
     
!
!            ************************************************************
!            *                                                          *
!            * Normal grid velocities of the faces.                     *
!            *                                                          *
!            ************************************************************
!
    ! sFaceI=	dot(sI, v)=dot(sI, v_freestream + v_grid)=dot(sI, v_freestream) + dot(sI, v_grid)
    ! sFaceJ, sFaceK same rule!

    ! dot(sI, v_freestream)
    ie_l = flowDoms(nn, 1, sps)%ie
    je_l = flowDoms(nn, 1, sps)%je
    ke_l = flowDoms(nn, 1, sps)%ke  


    ! i
    do k=1, ke_l
       do j=1, je_l
          do i=0, ie_l

             sFaceI(i, j, k) = velxGrid*sI(i, j, k, 1) + velyGrid*sI(i, j, k, 2) &
                             + velzGrid*sI(i, j, k, 3)


          end do
       end do
    end do

    ! j
    do k=1, ke_l
       do j=0, je_l
          do i=1, ie_l

             sFaceJ(i, j, k) = velxGrid*sJ(i, j, k, 1) + velyGrid*sJ(i, j, k, 2) &
                             + velzGrid*sJ(i, j, k, 3)


          end do
       end do
    end do

    ! k
    do k=0, ke_l
       do j=1, je_l
          do i=1, ie_l

             sFaceK(i, j, k) = velxGrid*sK(i, j, k, 1) + velyGrid*sK(i, j, k, 2) &
                             + velzGrid*sK(i, j, k, 3)


          end do
       end do
    end do


    !  dot(sI, v_grid)
    ! Loop over inner cells (also 1st layer halo of 3 surfs sI, sJ and sK are handled here, the left will be handled in the next section)

    do mm=1,nTimeIntervalsSpectral

       ie_l = flowDoms(nn, 1, mm)%ie
       je_l = flowDoms(nn, 1, mm)%je
       ke_l = flowDoms(nn, 1, mm)%ke


       ! i
       do k=1, ke_l
          do j=1, je_l
             do i=0, ie_l
   
                x_fc = fourth*(flowDoms(nn, 1, mm)%x(  i,j-1,k-1, 1) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 1)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 1) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 1))
                y_fc = fourth*(flowDoms(nn, 1, mm)%x(  i,j-1,k-1, 2) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 2)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 2) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 2))
                z_fc = fourth*(flowDoms(nn, 1, mm)%x(  i,j-1,k-1, 3) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 3)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 3) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 3))

                sFaceI(i, j, k) = sFaceI(i, j, k) & 
                                + dscalar(1, sps, mm)*x_fc*sI(i, j, k, 1) &
                                + dscalar(1, sps, mm)*y_fc*sI(i, j, k, 2) &
                                + dscalar(1, sps, mm)*z_fc*sI(i, j, k, 3)                

             end do
          end do
       end do       

       ! j
       do k=1, ke_l
          do j=0, je_l
             do i=1, ie_l
   
                x_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,  j,k-1, 1) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 1)&
                             + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 1) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 1))
                y_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,  j,k-1, 2) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 2)&
                             + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 2) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 2))
                z_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,  j,k-1, 3) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 3)&
                             + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 3) + flowDoms(nn, 1, mm)%x(  i,  j,k-1, 3))

                sFaceJ(i, j, k) = sFaceJ(i, j, k)                          & 
                                + dscalar(1, sps, mm)*x_fc*sJ(i, j, k, 1) &
                                + dscalar(1, sps, mm)*y_fc*sJ(i, j, k, 2) &
                                + dscalar(1, sps, mm)*z_fc*sJ(i, j, k, 3)              

             end do
          end do
       end do

       ! k
       do k=0, ke_l
          do j=1, je_l
             do i=1, ie_l
   
                x_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,j-1,  k, 1) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 1)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 1) + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 1))
                y_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,j-1,  k, 2) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 2)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 2) + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 2))
                z_fc = fourth*(flowDoms(nn, 1, mm)%x(i-1,j-1,  k, 3) + flowDoms(nn, 1, mm)%x(  i,  j,  k, 3)&
                             + flowDoms(nn, 1, mm)%x(  i,j-1,  k, 3) + flowDoms(nn, 1, mm)%x(i-1,  j,  k, 3))

                sFaceK(i, j, k) = sFaceK(i, j, k)                          & 
                                + dscalar(1, sps, mm)*x_fc*sK(i, j, k, 1) &
                                + dscalar(1, sps, mm)*y_fc*sK(i, j, k, 2) &
                                + dscalar(1, sps, mm)*z_fc*sK(i, j, k, 3)               

             end do
          end do
       end do  
     
    end do


   
  end subroutine gridVelocitiesFineLevel_TS_block






  subroutine gridVelocitiesFineLevel_block(useOldCoor, t, sps)
    !
    !       gridVelocitiesFineLevel computes the grid velocities for
    !       the cell centers and the normal grid velocities for the faces
    !       of moving blocks for the currently finest grid, i.e.
    !       groundLevel. The velocities are computed at time t for
    !       spectral mode sps. If useOldCoor is .true. the velocities
    !       are determined using the unsteady time integrator in
    !       combination with the old coordinates; otherwise the analytic
    !       form is used.
    !
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication
    use flowUtils, only :  derivativeRotMatrixRigid, getDirVector
    use utils, only : setCoefTimeIntegrator,tsAlpha, tsBeta, tsMach, terminate, &
         rotMatrixRigidBody, getDirAngle

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor

    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: i, j, k, ii, iie, jje, kke

    real(kind=realType) :: oneOver4dt, oneOver8dt
    real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
    real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

    real(kind=realType), dimension(3) :: sc, xc, xxc
    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix,&
         derivRotationMatrix

    real(kind=realType) :: tNew, tOld
    real(kind=realType), dimension(:,:), pointer :: sFace

    real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
    real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

    real(kind=realType) :: intervalMach,alphaTS,alphaIncrement,&
         betaTS,betaIncrement
    real(kind=realType), dimension(3) ::velDir
    real(kind=realType), dimension(3) :: refDirection


    ! Compute the mesh velocity from the given mesh Mach number.

    ! vel{x,y,z}Grid0 is the ACTUAL velocity you want at the
    ! geometry.
    aInf = sqrt(gammaInf*pInf/rhoInf)
    velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
    velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
    velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

    ! Compute the derivative of the rotation matrix and the rotation
    ! point; needed for velocity due to the rigid body rotation of
    ! the entire grid. It is assumed that the rigid body motion of
    ! the grid is only specified if there is only 1 section present.

    call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, t(1))

    !compute the rotation matrix to update the velocities for the time
    !spectral stability derivative case...

    if(TSStability)then
       ! Determine the time values of the old and new time level.
       ! It is assumed that the rigid body rotation of the mesh is only
       ! used when only 1 section is present.

       tNew = timeUnsteady + timeUnsteadyRestart
       tOld = tNew - t(1)

       if(TSpMode.or. TSqMode .or.TSrMode) then
          ! Compute the rotation matrix of the rigid body rotation as
          ! well as the rotation point; the latter may vary in time due
          ! to rigid body translation.

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

          if(TSAlphaFollowing) then

             velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
                  + rotationMatrix(1,2)*velygrid0 &
                  + rotationMatrix(1,3)*velzgrid0
             velygrid0 = rotationMatrix(2,1)*velxgrid0 &
                  + rotationMatrix(2,2)*velygrid0 &
                  + rotationMatrix(2,3)*velzgrid0
             velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
                  + rotationMatrix(3,2)*velygrid0 &
                  + rotationMatrix(3,3)*velzgrid0

          end if

       elseif(tsAlphaMode)then

          !Determine the alpha for this time instance
          alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
               degreeFourAlpha,  omegaFourAlpha,     &
               cosCoefFourAlpha, sinCoefFourAlpha, t(1))

          alphaTS = alpha+alphaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))

       elseif(tsBetaMode)then

          !Determine the alpha for this time instance
          betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
               degreeFourBeta,  omegaFourBeta,     &
               cosCoefFourBeta, sinCoefFourBeta, t(1))

          betaTS = beta+betaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))
       elseif(TSMachMode)then
          !determine the mach number at this time interval
          IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
               degreeFourMach,  omegaFourMach,     &
               cosCoefFourMach, sinCoefFourMach, t(1))
          !set the effective grid velocity
          velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
          velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
          velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

       elseif(TSAltitudeMode)then
          call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
       else
          call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
       end if
    endif

    testMoving: if( blockIsMoving ) then
       ! Determine the situation we are having here.

       testUseOldCoor: if( useOldCoor ) then
          !
          !             The velocities must be determined via a finite
          !             difference formula using the coordinates of the old
          !             levels.
          !
          ! Set the coefficients for the time integrator and store
          ! the inverse of the physical nonDimensional time step,
          ! divided by 4 and 8, a bit easier.

          call setCoefTimeIntegrator
          oneOver4dt = fourth*timeRef/deltaT
          oneOver8dt = half*oneOver4dt
          !
          !             Grid velocities of the cell centers, including the
          !             1st level halo cells.
          !
          ! Loop over the cells, including the 1st level halo's.

          do k=1,ke
             do j=1,je
                do i=1,ie

                   ! The velocity of the cell center is determined
                   ! by a finite difference formula. First store
                   ! the current coordinate, multiplied by 8 and
                   ! coefTime(0) in sc.

                   sc(1) = (x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                        +  x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                        +  x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                        +  x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
                        * coefTime(0)
                   sc(2) = (x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                        +  x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                        +  x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                        +  x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
                        * coefTime(0)
                   sc(3) = (x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                        +  x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                        +  x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                        +  x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
                        * coefTime(0)

                   ! Loop over the older levels to complete the
                   ! finite difference formula.

                   do ii=1,nOldLevels
                      sc(1) = sc(1) + (xOld(ii,i-1,j-1,k-1,1)  &
                           +          xOld(ii,i,  j-1,k-1,1)  &
                           +          xOld(ii,i-1,j,  k-1,1)  &
                           +          xOld(ii,i,  j,  k-1,1)  &
                           +          xOld(ii,i-1,j-1,k,  1)  &
                           +          xOld(ii,i,  j-1,k,  1)  &
                           +          xOld(ii,i-1,j,  k,  1)  &
                           +          xOld(ii,i,  j,  k,  1)) &
                           * coefTime(ii)
                      sc(2) = sc(2) + (xOld(ii,i-1,j-1,k-1,2)  &
                           +          xOld(ii,i,  j-1,k-1,2)  &
                           +          xOld(ii,i-1,j,  k-1,2)  &
                           +          xOld(ii,i,  j,  k-1,2)  &
                           +          xOld(ii,i-1,j-1,k,  2)  &
                           +          xOld(ii,i,  j-1,k,  2)  &
                           +          xOld(ii,i-1,j,  k,  2)  &
                           +          xOld(ii,i,  j,  k,  2)) &
                           * coefTime(ii)
                      sc(3) = sc(3) + (xOld(ii,i-1,j-1,k-1,3)  &
                           +          xOld(ii,i,  j-1,k-1,3)  &
                           +          xOld(ii,i-1,j,  k-1,3)  &
                           +          xOld(ii,i,  j,  k-1,3)  &
                           +          xOld(ii,i-1,j-1,k,  3)  &
                           +          xOld(ii,i,  j-1,k,  3)  &
                           +          xOld(ii,i-1,j,  k,  3)  &
                           +          xOld(ii,i,  j,  k,  3)) &
                           * coefTime(ii)
                   enddo

                   ! Divide by 8 delta t to obtain the correct
                   ! velocities.

                   s(i,j,k,1) = sc(1)*oneOver8dt
                   s(i,j,k,2) = sc(2)*oneOver8dt
                   s(i,j,k,3) = sc(3)*oneOver8dt
                enddo
             enddo
          enddo
          !
          !             Normal grid velocities of the faces.
          !
          ! Loop over the three directions.

          loopDir: do mm=1,3

             ! Set the upper boundaries depending on the direction.

             select case (mm)
             case (1_intType)       ! normals in i-direction
                iie = ie; jje = je; kke = ke

             case (2_intType)       ! normals in j-direction
                iie = je; jje = ie; kke = ke

             case (3_intType)       ! normals in k-direction
                iie = ke; jje = ie; kke = je
             end select
             !
             !               Normal grid velocities in generalized i-direction.
             !               Mm == 1: i-direction
             !               mm == 2: j-direction
             !               mm == 3: k-direction
             !
             do i=0,iie

                ! Set the pointers for the coordinates, normals and
                ! normal velocities for this generalized i-plane.
                ! This depends on the value of mm.

                select case (mm)
                case (1_intType)       ! normals in i-direction
                   xx =>  x(i,:,:,:);  xxOld => xOld(:,i,:,:,:)
                   ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)

                case (2_intType)       ! normals in j-direction
                   xx =>  x(:,i,:,:);  xxOld => xOld(:,:,i,:,:)
                   ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)

                case (3_intType)       ! normals in k-direction
                   xx =>  x(:,:,i,:);  xxOld => xOld(:,:,:,i,:)
                   ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
                end select

                ! Loop over the k and j-direction of this
                ! generalized i-face. Note that due to the usage of
                ! the pointers xx and xxOld an offset of +1 must be
                ! used in the coordinate arrays, because x and xOld
                ! originally start at 0 for the i, j and k indices.

                do k=1,kke
                   do j=1,jje

                      ! The velocity of the face center is determined
                      ! by a finite difference formula. First store
                      ! the current coordinate, multiplied by 4 and
                      ! coefTime(0) in sc.

                      sc(1) = coefTime(0)*(xx(j+1,k+1,1) + xx(j,k+1,1) &
                           +              xx(j+1,k,  1) + xx(j,k,  1))
                      sc(2) = coefTime(0)*(xx(j+1,k+1,2) + xx(j,k+1,2) &
                           +              xx(j+1,k,  2) + xx(j,k,  2))
                      sc(3) = coefTime(0)*(xx(j+1,k+1,3) + xx(j,k+1,3) &
                           +              xx(j+1,k,  3) + xx(j,k,  3))

                      ! Loop over the older levels to complete the
                      ! finite difference.

                      do ii=1,nOldLevels

                         sc(1) = sc(1) + coefTime(ii)         &
                              *         (xxOld(ii,j+1,k+1,1) &
                              +          xxOld(ii,j,  k+1,1) &
                              +          xxOld(ii,j+1,k,  1) &
                              +          xxOld(ii,j,  k,  1))
                         sc(2) = sc(2) + coefTime(ii)         &
                              *         (xxOld(ii,j+1,k+1,2) &
                              +          xxOld(ii,j,  k+1,2) &
                              +          xxOld(ii,j+1,k,  2) &
                              +          xxOld(ii,j,  k,  2))
                         sc(3) = sc(3) + coefTime(ii)         &
                              *         (xxOld(ii,j+1,k+1,3) &
                              +          xxOld(ii,j,  k+1,3) &
                              +          xxOld(ii,j+1,k,  3) &
                              +          xxOld(ii,j,  k,  3))
                      enddo

                      ! Determine the dot product of sc and the normal
                      ! and divide by 4 deltaT to obtain the correct
                      ! value of the normal velocity.

                      sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
                           + sc(3)*ss(j,k,3)
                      sFace(j,k) = sFace(j,k)*oneOver4dt

                   enddo
                enddo
             enddo

          enddo loopDir

       else testUseOldCoor
          !
          !             The velocities must be determined analytically.
          !
          ! Store the rotation center and determine the
          ! nonDimensional rotation rate of this block. As the
          ! reference length is 1 timeRef == 1/uRef and at the end
          ! the nonDimensional velocity is computed.

          j = nbkGlobal

          rotCenter = cgnsDoms(j)%rotCenter
          rotRate   = timeRef*cgnsDoms(j)%rotRate

          velXgrid = velXGrid0
          velYgrid = velYGrid0
          velZgrid = velZGrid0
          !
          !             Grid velocities of the cell centers, including the
          !             1st level halo cells.
          !
          ! Loop over the cells, including the 1st level halo's.

          do k=1,ke
             do j=1,je
                do i=1,ie

                   ! Determine the coordinates of the cell center,
                   ! which are stored in xc.

                   xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1) &
                        +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1) &
                        +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1) &
                        +         x(i-1,j,  k,  1) + x(i,j,  k,  1))
                   xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2) &
                        +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2) &
                        +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2) &
                        +         x(i-1,j,  k,  2) + x(i,j,  k,  2))
                   xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3) &
                        +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3) &
                        +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3) &
                        +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

                   ! Determine the coordinates relative to the
                   ! center of rotation.

                   xxc(1) = xc(1) - rotCenter(1)
                   xxc(2) = xc(2) - rotCenter(2)
                   xxc(3) = xc(3) - rotCenter(3)

                   ! Determine the rotation speed of the cell center,
                   ! which is omega*r.

                   sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                   sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                   sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                   ! Determine the coordinates relative to the
                   ! rigid body rotation point.

                   xxc(1) = xc(1) - rotationPoint(1)
                   xxc(2) = xc(2) - rotationPoint(2)
                   xxc(3) = xc(3) - rotationPoint(3)

                   ! Determine the total velocity of the cell center.
                   ! This is a combination of rotation speed of this
                   ! block and the entire rigid body rotation.

                   s(i,j,k,1) = sc(1) + velxGrid           &
                        + derivRotationMatrix(1,1)*xxc(1) &
                        + derivRotationMatrix(1,2)*xxc(2) &
                        + derivRotationMatrix(1,3)*xxc(3)
                   s(i,j,k,2) = sc(2) + velyGrid           &
                        + derivRotationMatrix(2,1)*xxc(1) &
                        + derivRotationMatrix(2,2)*xxc(2) &
                        + derivRotationMatrix(2,3)*xxc(3)
                   s(i,j,k,3) = sc(3) + velzGrid           &
                        + derivRotationMatrix(3,1)*xxc(1) &
                        + derivRotationMatrix(3,2)*xxc(2) &
                        + derivRotationMatrix(3,3)*xxc(3)
                enddo
             enddo
          enddo
          !
          !             Normal grid velocities of the faces.
          !
          ! Loop over the three directions.

          loopDirection: do mm=1,3

             ! Set the upper boundaries depending on the direction.

             select case (mm)
             case (1_intType)       ! Normals in i-direction
                iie = ie; jje = je; kke = ke

             case (2_intType)       ! Normals in j-direction
                iie = je; jje = ie; kke = ke

             case (3_intType)       ! Normals in k-direction
                iie = ke; jje = ie; kke = je
             end select
             !
             !               Normal grid velocities in generalized i-direction.
             !               mm == 1: i-direction
             !               mm == 2: j-direction
             !               mm == 3: k-direction
             !
             do i=0,iie

                ! Set the pointers for the coordinates, normals and
                ! normal velocities for this generalized i-plane.
                ! This depends on the value of mm.

                select case (mm)
                case (1_intType)       ! normals in i-direction
                   xx =>  x(i,:,:,:)
                   ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)

                case (2_intType)       ! normals in j-direction
                   xx =>  x(:,i,:,:)
                   ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)

                case (3_intType)       ! normals in k-direction
                   xx =>  x(:,:,i,:)
                   ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
                end select

                ! Loop over the k and j-direction of this generalized
                ! i-face. Note that due to the usage of the pointer
                ! xx an offset of +1 must be used in the coordinate
                ! array, because x originally starts at 0 for the
                ! i, j and k indices.

                do k=1,kke
                   do j=1,jje

                      ! Determine the coordinates of the face center,
                      ! which are stored in xc.

                      xc(1) = fourth*(xx(j+1,k+1,1) + xx(j,k+1,1) &
                           +         xx(j+1,k,  1) + xx(j,k,  1))
                      xc(2) = fourth*(xx(j+1,k+1,2) + xx(j,k+1,2) &
                           +         xx(j+1,k,  2) + xx(j,k,  2))
                      xc(3) = fourth*(xx(j+1,k+1,3) + xx(j,k+1,3) &
                           +         xx(j+1,k,  3) + xx(j,k,  3))

                      ! Determine the coordinates relative to the
                      ! center of rotation.

                      xxc(1) = xc(1) - rotCenter(1)
                      xxc(2) = xc(2) - rotCenter(2)
                      xxc(3) = xc(3) - rotCenter(3)

                      ! Determine the rotation speed of the face center,
                      ! which is omega*r.

                      sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                      sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                      sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                      ! Determine the coordinates relative to the
                      ! rigid body rotation point.

                      xxc(1) = xc(1) - rotationPoint(1)
                      xxc(2) = xc(2) - rotationPoint(2)
                      xxc(3) = xc(3) - rotationPoint(3)

                      ! Determine the total velocity of the cell face.
                      ! This is a combination of rotation speed of this
                      ! block and the entire rigid body rotation.

                      sc(1) = sc(1) + velxGrid           &
                           + derivRotationMatrix(1,1)*xxc(1) &
                           + derivRotationMatrix(1,2)*xxc(2) &
                           + derivRotationMatrix(1,3)*xxc(3)
                      sc(2) = sc(2) + velyGrid           &
                           + derivRotationMatrix(2,1)*xxc(1) &
                           + derivRotationMatrix(2,2)*xxc(2) &
                           + derivRotationMatrix(2,3)*xxc(3)
                      sc(3) = sc(3) + velzGrid           &
                           + derivRotationMatrix(3,1)*xxc(1) &
                           + derivRotationMatrix(3,2)*xxc(2) &
                           + derivRotationMatrix(3,3)*xxc(3)

                      ! Store the dot product of grid velocity sc and
                      ! the normal ss in sFace.

                      sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
                           + sc(3)*ss(j,k,3)

                   enddo
                enddo
             enddo

          enddo loopDirection
       endif testUseOldCoor
    endif testMoving

  end subroutine gridVelocitiesFineLevel_block

#ifndef USE_TAPENADE
  subroutine slipVelocitiesFineLevel(useOldCoor, t, sps)
    !
    ! Shell function to call slipVelocitiesFineLevel on all blocks
    !
    use constants
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : groundLevel
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor
    real(kind=realType), dimension(*), intent(in) :: t  !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Loop over the number of blocks.

    domains: do nn=1,nDom

       ! Set the pointers for this block.

       call setPointers(nn, groundLevel, sps)

       call slipVelocitiesFineLevel_block(useOldCoor, t, sps)

    end do domains

  end subroutine slipVelocitiesFineLevel
#endif

  subroutine slipVelocitiesFineLevel_block(useOldCoor, t, sps)
    !
    !       slipVelocitiesFineLevel computes the slip velocities for
    !       viscous subfaces on all viscous boundaries on groundLevel for
    !       the given spectral solution. If useOldCoor is .true. the
    !       velocities are determined using the unsteady time integrator;
    !       otherwise the analytic form is used.
    !
    use constants
    use inputTimeSpectral
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication
   use flowUtils, only :  derivativeRotMatrixRigid, getDirVector
    use utils, only : tsAlpha, tsBeta, tsMach, terminate, rotMatrixRigidBody, &
         setCoefTimeIntegrator, getDirAngle
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor

    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, i, j, level

    real(kind=realType) :: oneOver4dt
    real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
    real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

    real(kind=realType), dimension(3) :: xc, xxc
    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix,&
         derivRotationMatrix

    real(kind=realType) :: tNew, tOld

    real(kind=realType), dimension(:,:,:),   pointer :: uSlip
    real(kind=realType), dimension(:,:,:),   pointer :: xFace
    real(kind=realType), dimension(:,:,:,:), pointer :: xFaceOld

    real(kind=realType) :: intervalMach,alphaTS,alphaIncrement,&
         betaTS,betaIncrement
    real(kind=realType), dimension(3) ::velDir
    real(kind=realType), dimension(3) :: refDirection

    ! Determine the situation we are having here.

    testUseOldCoor: if( useOldCoor ) then

       ! The velocities must be determined via a finite difference
       ! formula using the coordinates of the old levels.

       ! Set the coefficients for the time integrator and store the
       ! inverse of the physical nonDimensional time step, divided
       ! by 4, a bit easier.

       call setCoefTimeIntegrator
       oneOver4dt = fourth*timeRef/deltaT

       ! Loop over the number of viscous subfaces.

       bocoLoop1: do mm=1,nViscBocos

          ! Set the pointer for uSlip to make the code more
          ! readable.

          uSlip => BCData(mm)%uSlip

          ! Determine the grid face on which the subface is located
          ! and set some variables accordingly.

          select case (BCFaceID(mm))

          case (iMin)
             xFace => x(1,:,:,:);  xFaceOld => xOld(:,1,:,:,:)

          case (iMax)
             xFace => x(il,:,:,:); xFaceOld => xOld(:,il,:,:,:)

          case (jMin)
             xFace => x(:,1,:,:);  xFaceOld => xOld(:,:,1,:,:)

          case (jMax)
             xFace => x(:,jl,:,:); xFaceOld => xOld(:,:,jl,:,:)

          case (kMin)
             xFace => x(:,:,1,:);  xFaceOld => xOld(:,:,:,1,:)

          case (kMax)
             xFace => x(:,:,kl,:); xFaceOld => xOld(:,:,:,kl,:)

          end select

          ! Some boundary faces have a different rotation speed than
          ! the corresponding block. This happens e.g. in the tip gap
          ! region of turboMachinary problems where the casing does
          ! not rotate. As the coordinate difference corresponds to
          ! the rotation rate of the block, a correction must be
          ! computed. Therefore compute the difference in rotation
          ! rate and store the rotation center a bit easier. Note that
          ! the rotation center of subface is taken, because if there
          ! is a difference in rotation rate this info for the subface
          ! must always be specified.

          j = nbkGlobal
          i = cgnsSubface(mm)

          rotCenter = cgnsDoms(j)%bocoInfo(i)%rotCenter
          rotRate   = timeRef*(cgnsDoms(j)%bocoInfo(i)%rotRate &
               -          cgnsDoms(j)%rotRate)

          ! Loop over the quadrilateral faces of the viscous subface.
          ! Note that due to the usage of the pointers xFace and
          ! xFaceOld an offset of +1 must be used in the coordinate
          ! arrays, because x and xOld originally start at 0 for the
          ! i, j and k indices.

          do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                ! Determine the coordinates of the centroid of the
                ! face, multiplied by 4.

                xc(1) = xFace(i+1,j+1,1) + xFace(i+1,j,1) &
                     + xFace(i,  j+1,1) + xFace(i,  j,1)
                xc(2) = xFace(i+1,j+1,2) + xFace(i+1,j,2) &
                     + xFace(i,  j+1,2) + xFace(i,  j,2)
                xc(3) = xFace(i+1,j+1,3) + xFace(i+1,j,3) &
                     + xFace(i,  j+1,3) + xFace(i,  j,3)

                ! Multiply the sum of the 4 vertex coordinates with
                ! coefTime(0) to obtain the contribution for the
                ! current time level. The division by 4*deltaT will
                ! take place later. This is both more efficient and
                ! more accurate for extremely small time steps.

                uSlip(i,j,1) = coefTime(0)*xc(1)
                uSlip(i,j,2) = coefTime(0)*xc(2)
                uSlip(i,j,3) = coefTime(0)*xc(3)

                ! Loop over the older time levels and take their
                ! contribution into account.

                do level=1,nOldLevels

                   uSlip(i,j,1) = uSlip(i,j,1) + coefTime(level)      &
                        *          (xFaceOld(level,i+1,j+1,1) &
                        +           xFaceOld(level,i+1,j,  1) &
                        +           xFaceOld(level,i,  j+1,1) &
                        +           xFaceOld(level,i,  j,  1))

                   uSlip(i,j,2) = uSlip(i,j,2) + coefTime(level)      &
                        *          (xFaceOld(level,i+1,j+1,2) &
                        +           xFaceOld(level,i+1,j,  2) &
                        +           xFaceOld(level,i,  j+1,2) &
                        +           xFaceOld(level,i,  j,  2))

                   uSlip(i,j,3) = uSlip(i,j,3) + coefTime(level)      &
                        *          (xFaceOld(level,i+1,j+1,3) &
                        +           xFaceOld(level,i+1,j,  3) &
                        +           xFaceOld(level,i,  j+1,3) &
                        +           xFaceOld(level,i,  j,  3))
                enddo

                ! Divide by 4 times the time step to obtain the
                ! correct velocity.

                uSlip(i,j,1) = uSlip(i,j,1)*oneOver4dt
                uSlip(i,j,2) = uSlip(i,j,2)*oneOver4dt
                uSlip(i,j,3) = uSlip(i,j,3)*oneOver4dt

                ! Determine the correction due to the difference
                ! in rotation rate between the block and subface.

                ! First determine the coordinates relative to the
                ! rotation center. Remember that 4 times this value
                ! is currently stored in xc.

                xc(1) = fourth*xc(1) - rotCenter(1)
                xc(2) = fourth*xc(2) - rotCenter(2)
                xc(3) = fourth*xc(3) - rotCenter(3)

                ! Compute the velocity, which is the cross product
                ! of rotRate and xc and add it to uSlip.

                uSlip(i,j,1) = uSlip(i,j,1) &
                     + rotRate(2)*xc(3) - rotRate(3)*xc(2)
                uSlip(i,j,2) = uSlip(i,j,2) &
                     + rotRate(3)*xc(1) - rotRate(1)*xc(3)
                uSlip(i,j,3) = uSlip(i,j,3) &
                     + rotRate(1)*xc(2) - rotRate(2)*xc(1)

             enddo
          enddo

       enddo bocoLoop1

    else

       ! The velocities must be determined analytically.

       ! Compute the mesh velocity from the given mesh Mach number.

       !  aInf = sqrt(gammaInf*pInf/rhoInf)
       !  velxGrid = aInf*MachGrid(1)
       !  velyGrid = aInf*MachGrid(2)
       !  velzGrid = aInf*MachGrid(3)

       aInf = sqrt(gammaInf*pInf/rhoInf)
       velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
       velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
       velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

       ! Compute the derivative of the rotation matrix and the rotation
       ! point; needed for velocity due to the rigid body rotation of
       ! the entire grid. It is assumed that the rigid body motion of
       ! the grid is only specified if there is only 1 section present.

       call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, &
            t(1))

       !compute the rotation matrix to update the velocities for the time
       !spectral stability derivative case...

       if(TSStability)then
          ! Determine the time values of the old and new time level.
          ! It is assumed that the rigid body rotation of the mesh is only
          ! used when only 1 section is present.

          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)

          if(TSpMode.or. TSqMode .or.TSrMode) then
             ! Compute the rotation matrix of the rigid body rotation as
             ! well as the rotation point; the latter may vary in time due
             ! to rigid body translation.

             call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

             if(TSAlphaFollowing) then

                velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
                     + rotationMatrix(1,2)*velygrid0 &
                     + rotationMatrix(1,3)*velzgrid0
                velygrid0 = rotationMatrix(2,1)*velxgrid0 &
                     + rotationMatrix(2,2)*velygrid0 &
                     + rotationMatrix(2,3)*velzgrid0
                velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
                     + rotationMatrix(3,2)*velygrid0 &
                     + rotationMatrix(3,3)*velzgrid0

             endif
          elseif(tsAlphaMode)then
             !Determine the alpha for this time instance
             alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
                  degreeFourAlpha,  omegaFourAlpha,     &
                  cosCoefFourAlpha, sinCoefFourAlpha, t(1))

             alphaTS = alpha+alphaIncrement
             !Determine the grid velocity for this alpha
             refDirection(:) = zero
             refDirection(1) = one
             call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgrid)*(-velDir(1))
             velyGrid0 = (aInf*machgrid)*(-velDir(2))
             velzGrid0 = (aInf*machgrid)*(-velDir(3))

          elseif(tsBetaMode)then

             !Determine the alpha for this time instance
             betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
                  degreeFourBeta,  omegaFourBeta,     &
                  cosCoefFourBeta, sinCoefFourBeta, t(1))

             betaTS = beta+betaIncrement
             !Determine the grid velocity for this alpha
             refDirection(:) = zero
             refDirection(1) = one
             call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgrid)*(-velDir(1))
             velyGrid0 = (aInf*machgrid)*(-velDir(2))
             velzGrid0 = (aInf*machgrid)*(-velDir(3))
          elseif(TSMachMode)then
             !determine the mach number at this time interval
             IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
                  degreeFourMach,  omegaFourMach,     &
                  cosCoefFourMach, sinCoefFourMach, t(1))
             !set the effective grid velocity
             velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
             velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
             velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

          elseif(TSAltitudeMode)then
             call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
          else
             call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
          end if
       endif

       ! Loop over the number of viscous subfaces.

       bocoLoop2: do mm=1,nViscBocos

          ! Determine the grid face on which the subface is located
          ! and set some variables accordingly.

          select case (BCFaceID(mm))

          case (iMin)
             xFace => x(1,:,:,:)

          case (iMax)
             xFace => x(il,:,:,:)

          case (jMin)
             xFace => x(:,1,:,:)

          case (jMax)
             xFace => x(:,jl,:,:)

          case (kMin)
             xFace => x(:,:,1,:)

          case (kMax)
             xFace => x(:,:,kl,:)

          end select

          ! Store the rotation center and the rotation rate
          ! for this subface.

          j = nbkGlobal
          i = cgnsSubface(mm)

          rotCenter = cgnsDoms(j)%bocoInfo(i)%rotCenter
          rotRate   = timeRef*cgnsDoms(j)%bocoInfo(i)%rotRate

          ! useWindAxis should go back here!
          velXgrid = velXGrid0
          velYgrid = velYGrid0
          velZgrid = velZGrid0

          ! Loop over the quadrilateral faces of the viscous
          ! subface.

          do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                ! Compute the coordinates of the centroid of the face.
                ! Normally this would be an average of i-1 and i, but
                ! due to the usage of the pointer xFace and the fact
                ! that x starts at index 0 this is shifted 1 index.

                xc(1) = fourth*(xFace(i+1,j+1,1) + xFace(i+1,j,1) &
                     +         xFace(i,  j+1,1) + xFace(i,  j,1))
                xc(2) = fourth*(xFace(i+1,j+1,2) + xFace(i+1,j,2) &
                     +         xFace(i,  j+1,2) + xFace(i,  j,2))
                xc(3) = fourth*(xFace(i+1,j+1,3) + xFace(i+1,j,3) &
                     +         xFace(i,  j+1,3) + xFace(i,  j,3))

                ! Determine the coordinates relative to the center
                ! of rotation.

                xxc(1) = xc(1) - rotCenter(1)
                xxc(2) = xc(2) - rotCenter(2)
                xxc(3) = xc(3) - rotCenter(3)

                ! Compute the velocity, which is the cross product
                ! of rotRate and xc.

                BCData(mm)%uSlip(i,j,1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                BCData(mm)%uSlip(i,j,2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                BCData(mm)%uSlip(i,j,3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                ! Determine the coordinates relative to the
                ! rigid body rotation point.

                xxc(1) = xc(1) - rotationPoint(1)
                xxc(2) = xc(2) - rotationPoint(2)
                xxc(3) = xc(3) - rotationPoint(3)

                ! Determine the total velocity of the cell center.
                ! This is a combination of rotation speed of this
                ! block and the entire rigid body rotation.

                BCData(mm)% uSlip(i,j,1) = BCData(mm)%uSlip(i,j,1) + velxGrid    &
                     + derivRotationMatrix(1,1)*xxc(1) &
                     + derivRotationMatrix(1,2)*xxc(2) &
                     + derivRotationMatrix(1,3)*xxc(3)
                BCData(mm)%uSlip(i,j,2) = BCData(mm)%uSlip(i,j,2) + velyGrid    &
                     + derivRotationMatrix(2,1)*xxc(1) &
                     + derivRotationMatrix(2,2)*xxc(2) &
                     + derivRotationMatrix(2,3)*xxc(3)
                BCData(mm)%uSlip(i,j,3) = BCData(mm)%uSlip(i,j,3) + velzGrid    &
                     + derivRotationMatrix(3,1)*xxc(1) &
                     + derivRotationMatrix(3,2)*xxc(2) &
                     + derivRotationMatrix(3,3)*xxc(3)
             enddo
          enddo

       enddo bocoLoop2


    endif testUseOldCoor

  end subroutine slipVelocitiesFineLevel_block




#ifndef USE_TAPENADE
  subroutine normalVelocitiesAllLevels(sps)
    !
    ! Shell function to call normalVelocities_block on all blocks/levels
    !
    use constants
    use blockPointers, only : nDom, flowDoms
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : groundLevel
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    !      Local variables.
    !
    integer(kind=intType) :: nn, level, nLevels

    nLevels = ubound(flowDoms,2)
    levelLoop: do level=groundLevel, nLevels
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, level, sps)

          call normalVelocities_block(sps)

       end do domains
    end do levelLoop
  end subroutine normalVelocitiesAllLevels
#endif

  subroutine normalVelocities_block(sps)
    !
    !       normalVelocitiesAllLevels computes the normal grid
    !       velocities of some boundary faces of the moving blocks for
    !       spectral mode sps. All grid levels from ground level to the
    !       coarsest level are considered.
    !
    use constants
    use blockPointers, only : il, jl, kl, addGridVelocities, nBocos, BCData, &
         sfaceI, sfaceJ, sfaceK, bcFaceID, si, sj, sk
    !use iteration
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm
    integer(kind=intType) :: i, j

    real(kind=realType) :: weight, mult

    real(kind=realType), dimension(:,:),   pointer :: sFace
    real(kind=realType), dimension(:,:,:), pointer :: ss

    ! Check for a moving block. As it is possible that in a
    ! multidisicplinary environment additional grid velocities
    ! are set, the test should be done on addGridVelocities
    ! and not on blockIsMoving.

    testMoving: if( addGridVelocities ) then
       !
       !             Determine the normal grid velocities of the boundaries.
       !             As these values are based on the unit normal. A division
       !             by the length of the normal is needed.
       !             Furthermore the boundary unit normals are per definition
       !             outward pointing, while on the iMin, jMin and kMin
       !             boundaries the face normals are inward pointing. This
       !             is taken into account by the factor mult.
       !
       ! Loop over the boundary subfaces.

       bocoLoop: do mm=1,nBocos

          ! Check whether rFace is allocated.

          testAssoc: if( associated(BCData(mm)%rFace) ) then

             ! Determine the block face on which the subface is
             ! located and set some variables accordingly.

             select case (BCFaceID(mm))

             case (iMin)
                mult = -one
                ss => si(1,:,:,:);  sFace => sFaceI(1,:,:)
             case (iMax)
                mult = one
                ss => si(il,:,:,:); sFace => sFaceI(il,:,:)
             case (jMin)
                mult = -one
                ss => sj(:,1,:,:);  sFace => sFaceJ(:,1,:)
             case (jMax)
                mult = one
                ss => sj(:,jl,:,:); sFace => sFaceJ(:,jl,:)
             case (kMin)
                mult = -one
                ss => sk(:,:,1,:);  sFace => sFaceK(:,:,1)
             case (kMax)
                mult = one
                ss => sk(:,:,kl,:); sFace => sFaceK(:,:,kl)

             end select

             ! Loop over the faces of the subface.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   ! Compute the inverse of the length of the normal
                   ! vector and possibly correct for inward pointing.

                   weight = sqrt(ss(i,j,1)**2 + ss(i,j,2)**2 &
                        +      ss(i,j,3)**2)
                   if(weight > zero) weight = mult/weight

                   ! Compute the normal velocity based on the outward
                   ! pointing unit normal.

                   BCData(mm)%rFace(i,j) = weight*sFace(i,j)
                enddo
             enddo

          endif testAssoc
       enddo bocoLoop

    else testMoving

       ! Block is not moving. Loop over the boundary faces and set
       ! the normal grid velocity to zero if allocated.

       do mm=1,nBocos
          if( associated(BCData(mm)%rFace) ) &
               BCData(mm)%rFace = zero
       enddo

    endif testMoving

  end subroutine normalVelocities_block

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE

  subroutine shiftSolution
    !
    !       shiftSolution shifts the solution of the older time levels,
    !       such that a new time step can be started.
    !
    use constants
    use blockPointers, only: il, jl, kl, nbkglobal, wOld, w, nDom
    use cgnsGrid, only : cgnsDoms
    use flowvarrefstate, only : nw
    use iteration, only : groundLevel, nOldLevels
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputUnsteady, only : deltaT
    use monitor, only : timeUnsteadyRestart, timeUnsteady
    use utils, only : setPointers, rotMatrixRigidBody
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l, sps, nn, mm, ll

    real(kind=realType) :: tOld, tNew
    real(kind=realType) :: vvX, vvY, vvZ, vXi, vEta, vZeta
    real(kind=realType) :: t, angleX, angleY, angleZ
    real(kind=realType) :: phi, cosPhi, sinPhi
    real(kind=realType) :: xiX, xiY, xiZ, etaX, etaY, etaZ
    real(kind=realType) :: zetaX, zetaY, zetaZ

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix

    ! Compute the rotation matrix of the rigid body rotation as well
    ! as the rotation point; the latter is not needed to correct the
    ! velocities, but the routine rotMatrixRigidBody is also used
    ! for coordinates.

    tNew = timeUnsteady + timeUnsteadyRestart
    tOld = tNew         - deltaT

    call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

    ! Loop over the number of spectral solutions and local blocks.
    ! Although this routine is only called in unsteady mode where the
    ! number of spectral solutions is 1, this loop is there just for
    ! consistency.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block on the ground level.

          call setPointers(nn, groundLevel, sps)

          ! Shift the solution already stored in wOld.

          loopOldLevels: do mm=nOldLevels,2,-1

             ! Shift the owned solution variables from level mm-1 to mm.

             ll = mm - 1

             do l=1,nw
                do k=2,kl
                   do j=2,jl
                      do i=2,il
                         wOld(mm,i,j,k,l) = wOld(ll,i,j,k,l)
                      enddo
                   enddo
                enddo
             enddo

          enddo loopOldLevels

          ! Shift the current solution into the 1st level of wOld.
          ! Note that in wOld the conservative flow variables are stored,
          ! while in w the velocity components are stored and not
          ! the momentum. Therefore this must be corrected.
          ! Also the turbulent primitive variables are stored, but this
          ! is okay, because the quasi-linear form of the turbulent
          ! transport equations is solved and not the conservative one.

          do l=1,nw
             do k=2,kl
                do j=2,jl
                   do i=2,il
                      wOld(1,i,j,k,l) = w(i,j,k,l)
                   enddo
                enddo
             enddo
          enddo

          ! Make sure that the momentum variables are stored in wOld.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   wOld(1,i,j,k,ivx) = wOld(1,i,j,k,ivx)*wOld(1,i,j,k,irho)
                   wOld(1,i,j,k,ivy) = wOld(1,i,j,k,ivy)*wOld(1,i,j,k,irho)
                   wOld(1,i,j,k,ivz) = wOld(1,i,j,k,ivz)*wOld(1,i,j,k,irho)
                enddo
             enddo
          enddo

          ! To improve the initial guess of the velocity field the
          ! velocity of rotating parts is rotated. First the rigid
          ! body motion.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   vvX = w(i,j,k,ivx)
                   vvY = w(i,j,k,ivy)
                   vvZ = w(i,j,k,ivz)

                   w(i,j,k,ivx) = rotationMatrix(1,1)*vvX &
                        + rotationMatrix(1,2)*vvY &
                        + rotationMatrix(1,3)*vvZ
                   w(i,j,k,ivy) = rotationMatrix(2,1)*vvX &
                        + rotationMatrix(2,2)*vvY &
                        + rotationMatrix(2,3)*vvZ
                   w(i,j,k,ivz) = rotationMatrix(3,1)*vvX &
                        + rotationMatrix(3,2)*vvY &
                        + rotationMatrix(3,3)*vvZ
                enddo
             enddo
          enddo

          ! Apply an additional correction for the velocity components
          ! if a rotation rate is prescribed for this block.

          rotTest: if( cgnsDoms(nbkGlobal)%rotatingFrameSpecified ) then

             ! Compute the rotation angles.

             angleX = deltaT*cgnsDoms(nbkGlobal)%rotRate(1)
             angleY = deltaT*cgnsDoms(nbkGlobal)%rotRate(2)
             angleZ = deltaT*cgnsDoms(nbkGlobal)%rotRate(3)

             ! Compute the unit vector in the direction of the rotation
             ! axis, which will be called the xi-direction.

             t   = one/max(eps,sqrt(angleX**2 + angleY**2 + angleZ**2))
             xiX = t*angleX
             xiY = t*angleY
             xiZ = t*angleZ

             ! Determine the rotation angle in xi-direction and its sine
             ! and cosine. Due to the definition of the xi-direction this
             ! angle will always be positive.

             phi    = xiX*angleX + xiY*angleY + xiZ*angleZ
             cosPhi = cos(phi)
             sinPhi = sin(phi)

             ! Loop over the cell centers.

             do k=2,kl
                do j=2,jl
                   do i=2,il

                      ! Abbreviate the velocity components a bit easier.

                      vvX = w(i,j,k,ivx)
                      vvY = w(i,j,k,ivy)
                      vvZ = w(i,j,k,ivz)

                      ! Determine the component of the velocity vector
                      ! in xi direction and determine the direction eta,
                      ! the direction of the velocity when the xi component
                      ! is substracted.

                      vXi = vvX*xiX + vvY*xiY + vvZ*xiZ

                      etaX = vvX - vXi*xiX
                      etaY = vvY - vXi*xiY
                      etaZ = vvZ - vXi*xiZ

                      t    = one/max(eps,sqrt(etaX**2 + etaY**2 + etaZ**2))
                      etaX = t*etaX
                      etaY = t*etaY
                      etaZ = t*etaZ

                      ! Determine the velocity component in eta direction.

                      vEta = vvX*etaX + vvY*etaY + vvZ*etaZ

                      ! Determine the unit vector in zeta-direction. This is
                      ! the cross product of the unit vectors in xi and in
                      ! eta-direction.

                      zetaX = xiY*etaZ - xiZ*etaY
                      zetaY = xiZ*etaX - xiX*etaZ
                      zetaZ = xiX*etaY - xiY*etaX

                      ! Determine the velocity components in eta and zeta
                      ! direction after the rotation.

                      vZeta = vEta*sinPhi
                      vEta  = vEta*cosPhi

                      ! Compute the new Cartesian velocity components.

                      w(i,j,k,ivx) = vXi*xiX + vEta*etaX + vZeta*zetaX
                      w(i,j,k,ivy) = vXi*xiY + vEta*etaY + vZeta*zetaY
                      w(i,j,k,ivz) = vXi*xiZ + vEta*etaZ + vZeta*zetaZ

                   enddo
                enddo
             enddo

          endif rotTest

       enddo domains
    enddo spectralLoop

  end subroutine shiftSolution

  subroutine computeUtau
    !
    ! Shell function to call computUTau  on all blocks
    !
    use constants
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use Iteration, only : groundLevel
    use utils, only : setPointers
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, nn

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, groundLevel, sps)

          call computeUtau_block

       end do domains

    end do spectralLoop

  end subroutine computeUtau

  subroutine computeUtau_block
    !
    !       computeUtau computes the skin friction velocity for the
    !       viscous subfaces. This data is only needed if wall functions
    !       are used.
    !
    use constants
    use blockPointers
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use turbCurveFits, only : curveUpRe

    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm, i, j

    real(kind=realType) :: re, vvx, vvy, vvz, veln, veltmag

    real(kind=realType), dimension(:,:,:), pointer :: ww, norm, uSlip
    real(kind=realType), dimension(:,:),   pointer :: dd2Wall, rrlv
    real(kind=realType), dimension(:,:),   pointer :: utau

    ! Return immediately if no wall functions must be used.

    if(.not. wallFunctions) return

    ! Loop over the viscous subfaces of this block.

    viscSubfaces: do mm=1,nViscBocos

       ! Set a bunch of pointers depending on the face id to make
       ! a generic treatment possible.

       select case (BCFaceID(mm))

       case (iMin)
          ww => w(2,1:,1:,:);
          dd2Wall => d2Wall(2,:,:); rrlv => rlv(2,1:,1:)

          !=========================================================

       case (iMax)
          ww => w(il,1:,1:,:)
          dd2Wall => d2Wall(il,:,:); rrlv => rlv(il,1:,1:)

          !=========================================================

       case (jMin)
          ww => w(1:,2,1:,:)
          dd2Wall => d2Wall(:,2,:); rrlv => rlv(1:,2,1:)

          !=========================================================

       case (jMax)
          ww => w(1:,jl,1:,:)
          dd2Wall => d2Wall(:,jl,:); rrlv => rlv(1:,jl,1:)

          !=========================================================

       case (kMin)
          ww => w(1:,1:,2,:)
          dd2Wall => d2Wall(:,:,2); rrlv => rlv(1:,1:,2)

          !=========================================================

       case (kMax)
          ww => w(1:,1:,kl,:)
          dd2Wall => d2Wall(:,:,kl); rrlv => rlv(1:,1:,kl)

       end select

       ! Set the pointers for the unit outward normals, uSlip
       ! and utau to make the code more readable.

       norm  => BCData(mm)%norm
       uSlip => BCData(mm)%uSlip
       utau  => viscSubface(mm)%utau

       ! Loop over the quadrilateral faces of the subface. Note
       ! that the nodal range of BCData must be used and not the
       ! cell range, because the latter may include the halo's in i
       ! and j-direction. The offset +1 is there, because inBeg and
       ! jnBeg refer to nodal ranges and not to cell ranges.
       ! Note that an offset of -1 must be used in dd2Wall, because
       ! the original array d2Wall starts at 2.

       do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
          do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

             ! Compute the velocity difference between the internal
             ! cell and the wall.

             vvx = ww(i,j,ivx) - uSlip(i,j,1)
             vvy = ww(i,j,ivy) - uSlip(i,j,2)
             vvz = ww(i,j,ivz) - uSlip(i,j,3)

             ! Compute the normal velocity of the internal cell.

             veln  = vvx*norm(i,j,1) + vvy*norm(i,j,2) + vvz*norm(i,j,3)

             ! Compute the magnitude of the tangential velocity.

             veltmag = max(eps,sqrt(vvx*vvx + vvy*vvy + vvz*vvz - veln*veln))

             ! Compute the Reynolds number. Note that an offset of -1
             ! must be used in dd2Wall, because the original array
             ! d2Wall starts at 2.
             ! Afterwards compute utau.

             re = ww(i,j,irho)*veltmag*dd2Wall(i-1,j-1)/rrlv(i,j)
             utau(i,j) = veltmag/max(curveUpRe(re),eps)

          enddo
       enddo

    enddo viscSubfaces

  end subroutine computeUtau_block


  subroutine gridVelocitiesFineLevelPart1(useOldCoor, t, sps)
    !
    ! Shell function to call gridVelocitiesFineLevel on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor
    real(kind=realType), dimension(*), intent(in) :: t  !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Loop over the number of blocks.

    domains: do nn=1,nDom

       ! Set the pointers for this block.

       call setPointers(nn, groundLevel, sps)
       call gridVelocitiesFineLevelPart1_block(useOldCoor, t, sps)

    end do domains

  end subroutine gridVelocitiesFineLevelPart1

  subroutine gridVelocitiesFineLevelPart1_block(useOldCoor, t, sps)
    !
    !       gridVelocitiesFineLevel computes the grid velocities for
    !       the cell centers and the normal grid velocities for the faces
    !       of moving blocks for the currently finest grid, i.e.
    !       groundLevel. The velocities are computed at time t for
    !       spectral mode sps. If useOldCoor is .true. the velocities
    !       are determined using the unsteady time integrator in
    !       combination with the old coordinates; otherwise the analytic
    !       form is used.
    !       Now it is split up into two parts.
    !       First part calculate the grid velocity using FIRST order BDF.
    !       Second part calculate the surface normal and normal velocity.
    !
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication
    use utils, only : TSAlpha, TSBeta, TSMach, terminate, rotMatrixRigidBody, &
         setCoefTimeIntegrator, getDirAngle
    use flowUtils, only : derivativeRotMatrixRigid, getDirVector
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor

    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: i, j, k, ii, iie, jje, kke

    real(kind=realType) :: oneOver4dt, oneOver8dt
    real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
    real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

    real(kind=realType), dimension(3) :: sc, xc, xxc
    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix,&
         derivRotationMatrix

    real(kind=realType) :: tNew, tOld
    real(kind=realType), dimension(:,:), pointer :: sFace
    real(kind=realType), dimension(:,:,:),   pointer :: sVelo

    real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
    real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

    real(kind=realType) :: intervalMach,alphaTS,alphaIncrement,&
         betaTS,betaIncrement
    real(kind=realType), dimension(3) ::velDir
    real(kind=realType), dimension(3) :: refDirection

    ! Compute the mesh velocity from the given mesh Mach number.

    ! vel{x,y,z}Grid0 is the ACTUAL velocity you want at the
    ! geometry.
    aInf = sqrt(gammaInf*pInf/rhoInf)
    velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
    velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
    velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

    ! Compute the derivative of the rotation matrix and the rotation
    ! point; needed for velocity due to the rigid body rotation of
    ! the entire grid. It is assumed that the rigid body motion of
    ! the grid is only specified if there is only 1 section present.

    call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, t(1))

    !compute the rotation matrix to update the velocities for the time
    !spectral stability derivative case...

    if(TSStability)then
       ! Determine the time values of the old and new time level.
       ! It is assumed that the rigid body rotation of the mesh is only
       ! used when only 1 section is present.

       tNew = timeUnsteady + timeUnsteadyRestart
       tOld = tNew - t(1)

       if(TSpMode.or. TSqMode .or.TSrMode) then
          ! Compute the rotation matrix of the rigid body rotation as
          ! well as the rotation point; the latter may vary in time due
          ! to rigid body translation.

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

          velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
               + rotationMatrix(1,2)*velygrid0 &
               + rotationMatrix(1,3)*velzgrid0
          velygrid0 = rotationMatrix(2,1)*velxgrid0 &
               + rotationMatrix(2,2)*velygrid0 &
               + rotationMatrix(2,3)*velzgrid0
          velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
               + rotationMatrix(3,2)*velygrid0 &
               + rotationMatrix(3,3)*velzgrid0

       elseif(tsAlphaMode)then
          !Determine the alpha for this time instance
          alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
               degreeFourAlpha,  omegaFourAlpha,     &
               cosCoefFourAlpha, sinCoefFourAlpha, t(1))

          alphaTS = alpha+alphaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))

       elseif(tsBetaMode)then

          !Determine the alpha for this time instance
          betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
               degreeFourBeta,  omegaFourBeta,     &
               cosCoefFourBeta, sinCoefFourBeta, t(1))

          betaTS = beta+betaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))
       elseif(TSMachMode)then
          !determine the mach number at this time interval
          IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
               degreeFourMach,  omegaFourMach,     &
               cosCoefFourMach, sinCoefFourMach, t(1))
          !set the effective grid velocity
          velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
          velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
          velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

       elseif(TSAltitudeMode)then
          call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
       else
          call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
       end if
    endif

    testMoving: if( blockIsMoving ) then
       ! REMOVED the rigid body rotation part for simplicity

       !
       !             The velocities must be determined via a finite
       !             difference formula using the coordinates of the old
       !             levels.
       !
       ! Set the coefficients for the time integrator and store
       ! the inverse of the physical nonDimensional time step,
       ! divided by 4 and 8, a bit easier.

       call setCoefTimeIntegrator
       oneOver4dt = fourth*timeRef/deltaT
       oneOver8dt = half*oneOver4dt
       !
       !             Grid velocities of the cell centers, including the
       !             1st level halo cells.
       !
       ! Loop over the cells, including the 1st level halo's.

       do k=1,ke
          do j=1,je
             do i=1,ie

                ! Using FIRST order BDF for all cases
                ! Refer to eq. 11b, found paper by C.Farhat http://dx.doi.org/10.1016/S0021-9991(03)00311-5
                ! Same applies for the velocities of the faces below. Theta(n+1) = 1, Theta(n) = -1 therfore
                ! it becoms a first order scheme.

                ! The velocity of the cell center is determined
                ! by a finite difference formula. First store
                ! the current coordinate, multiplied by 8 and
                ! coefTime(0) in sc.

                sc(1) = (x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                     +  x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                     +  x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                     +  x(i-1,j,  k,  1) + x(i,j,  k,  1))
                sc(2) = (x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                     +  x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                     +  x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                     +  x(i-1,j,  k,  2) + x(i,j,  k,  2))
                sc(3) = (x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                     +  x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                     +  x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                     +  x(i-1,j,  k,  3) + x(i,j,  k,  3))

                ! Loop over the older levels to complete the
                ! finite difference formula.

                ii = 1 ! There was a loop over all old levels
                sc(1) = sc(1) + (xOld(ii,i-1,j-1,k-1,1)  &
                     +          xOld(ii,i,  j-1,k-1,1)  &
                     +          xOld(ii,i-1,j,  k-1,1)  &
                     +          xOld(ii,i,  j,  k-1,1)  &
                     +          xOld(ii,i-1,j-1,k,  1)  &
                     +          xOld(ii,i,  j-1,k,  1)  &
                     +          xOld(ii,i-1,j,  k,  1)  &
                     +          xOld(ii,i,  j,  k,  1)) &
                     * (-1.0_realType)
                sc(2) = sc(2) + (xOld(ii,i-1,j-1,k-1,2)  &
                     +          xOld(ii,i,  j-1,k-1,2)  &
                     +          xOld(ii,i-1,j,  k-1,2)  &
                     +          xOld(ii,i,  j,  k-1,2)  &
                     +          xOld(ii,i-1,j-1,k,  2)  &
                     +          xOld(ii,i,  j-1,k,  2)  &
                     +          xOld(ii,i-1,j,  k,  2)  &
                     +          xOld(ii,i,  j,  k,  2)) &
                     * (-1.0_realType)
                sc(3) = sc(3) + (xOld(ii,i-1,j-1,k-1,3)  &
                     +          xOld(ii,i,  j-1,k-1,3)  &
                     +          xOld(ii,i-1,j,  k-1,3)  &
                     +          xOld(ii,i,  j,  k-1,3)  &
                     +          xOld(ii,i-1,j-1,k,  3)  &
                     +          xOld(ii,i,  j-1,k,  3)  &
                     +          xOld(ii,i-1,j,  k,  3)  &
                     +          xOld(ii,i,  j,  k,  3)) &
                     * (-1.0_realType)

                ! Divide by 8 delta t to obtain the correct
                ! velocities.

                s(i,j,k,1) = sc(1)*oneOver8dt
                s(i,j,k,2) = sc(2)*oneOver8dt
                s(i,j,k,3) = sc(3)*oneOver8dt
             enddo
          enddo
       enddo

       !
       !             Velocities of the faces, vector.
       !
       ! Loop over the three directions.

       loopDir: do mm=1,3

          ! Set the upper boundaries depending on the direction.

          select case (mm)
          case (1_intType)       ! normals in i-direction
             iie = ie; jje = je; kke = ke

          case (2_intType)       ! normals in j-direction
             iie = je; jje = ie; kke = ke

          case (3_intType)       ! normals in k-direction
             iie = ke; jje = ie; kke = je
          end select
          !
          !               Face velocities in generalized i-direction.
          !               mm == 1: i-direction
          !               mm == 2: j-direction
          !               mm == 3: k-direction
          !
          do i=0,iie

             ! Set the pointers for the coordinates, normals and
             ! normal velocities for this generalized i-plane.
             ! This depends on the value of mm.

             select case (mm)
             case (1_intType)       ! normals in i-direction
                xx =>  x(i,:,:,:);  xxOld => xOld(:,i,:,:,:)
                sVelo => sVeloIALE(i,:,:,:)

             case (2_intType)       ! normals in j-direction
                xx =>  x(:,i,:,:);  xxOld => xOld(:,:,i,:,:)
                sVelo => sVeloJALE(:,i,:,:)

             case (3_intType)       ! normals in k-direction
                xx =>  x(:,:,i,:);  xxOld => xOld(:,:,:,i,:)
                sVelo => sVeloKALE(:,:,i,:)
             end select

             ! Loop over the k and j-direction of this
             ! generalized i-face. Note that due to the usage of
             ! the pointers xx and xxOld an offset of +1 must be
             ! used in the coordinate arrays, because x and xOld
             ! originally start at 0 for the i, j and k indices.
             ! print *, mm
             do k=1,kke
                do j=1,jje

                   ! The velocity of the face center is determined
                   ! by a finite difference formula. First store
                   ! the current coordinate, multiplied by 4 and
                   ! coefTime(0) in sc.

                   sc(1) = (xx(j+1,k+1,1) + xx(j,k+1,1) &
                        +  xx(j+1,k,  1) + xx(j,k,  1))
                   sc(2) = (xx(j+1,k+1,2) + xx(j,k+1,2) &
                        +  xx(j+1,k,  2) + xx(j,k,  2))
                   sc(3) = (xx(j+1,k+1,3) + xx(j,k+1,3) &
                        +  xx(j+1,k,  3) + xx(j,k,  3))

                   ii = 1 ! There was a loop who looped over nOldLevels
                   sc(1) = sc(1) + (xxOld(ii,j+1,k+1,1) &
                        +          xxOld(ii,j,  k+1,1) &
                        +          xxOld(ii,j+1,k,  1) &
                        +          xxOld(ii,j,  k,  1)) &
                        * (-1.0_realType)
                   sc(2) = sc(2) + (xxOld(ii,j+1,k+1,2) &
                        +          xxOld(ii,j,  k+1,2) &
                        +          xxOld(ii,j+1,k,  2) &
                        +          xxOld(ii,j,  k,  2)) &
                        * (-1.0_realType)
                   sc(3) = sc(3) + (xxOld(ii,j+1,k+1,3) &
                        +          xxOld(ii,j,  k+1,3) &
                        +          xxOld(ii,j+1,k,  3) &
                        +          xxOld(ii,j,  k,  3)) &
                        * (-1.0_realType)

                   ! Determine the dot product of sc and the normal
                   ! and divide by 4 deltaT to obtain the correct
                   ! value of the normal velocity.

                   sVelo(j,k,1) = sc(1)*oneOver4dt
                   sVelo(j,k,2) = sc(2)*oneOver4dt
                   sVelo(j,k,3) = sc(3)*oneOver4dt

                   ! if ((i.ge.2) .and. (i.le.3) .and. (j.ge.2) .and. (j.le.3) .and. (k.ge.2) .and. (k.le.3)) then
                   !    print *, i,j,k, sVelo(j,k,:)
                   !    print *, '                                   ', xx(j,k,:)
                   !    print *, '                                   ', xxOld(1,j,k,:)
                   ! end if

                enddo
             enddo
          enddo

       enddo loopDir

    endif testMoving

  end subroutine gridVelocitiesFineLevelPart1_block

  !
  !      Here begins the second part
  !

  subroutine gridVelocitiesFineLevelPart2(useOldCoor, t, sps)
    !
    ! Shell function to call gridVelocitiesFineLevel on all blocks
    !
    use blockPointers
    use constants
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor
    real(kind=realType), dimension(*), intent(in) :: t  !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Loop over the number of blocks.

    domains: do nn=1,nDom

       ! Set the pointers for this block.

       call setPointers(nn, groundLevel, sps)
       call gridVelocitiesFineLevelPart2_block(useOldCoor, t, sps)

    end do domains

  end subroutine gridVelocitiesFineLevelPart2

  subroutine gridVelocitiesFineLevelPart2_block(useOldCoor, t, sps)
    !
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor
    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: i, j, k, ii, iie, jje, kke

    real(kind=realType) :: oneOver4dt, oneOver8dt
    real(kind=realType), dimension(3) :: sc, xc, xxc
    real(kind=realType), dimension(:,:), pointer :: sFace
    real(kind=realType), dimension(:,:,:),   pointer :: sVelo
    real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
    real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

    testMoving: if( blockIsMoving ) then
       !
       !             Normal grid velocities of the faces.
       !
       ! Loop over the three directions.

       loopDir: do mm=1,3

          ! Set the upper boundaries depending on the direction.

          select case (mm)
          case (1_intType)       ! normals in i-direction
             iie = ie; jje = je; kke = ke

          case (2_intType)       ! normals in j-direction
             iie = je; jje = ie; kke = ke

          case (3_intType)       ! normals in k-direction
             iie = ke; jje = ie; kke = je
          end select
          !
          !               Normal grid velocities in generalized i-direction.
          !               Mm == 1: i-direction
          !               mm == 2: j-direction
          !               mm == 3: k-direction
          !
          do i=0,iie

             ! Set the pointers for the coordinates, normals and
             ! normal velocities for this generalized i-plane.
             ! This depends on the value of mm.

             select case (mm)
             case (1_intType)       ! normals in i-direction
                ss => si(i,:,:,:);  sFace => sFaceI(i,:,:)
                sVelo => sVeloIALE(i,:,:,:)

             case (2_intType)       ! normals in j-direction
                ss => sj(:,i,:,:);  sFace => sFaceJ(:,i,:)
                sVelo => sVeloJALE(:,i,:,:)

             case (3_intType)       ! normals in k-direction
                ss => sk(:,:,i,:);  sFace => sFaceK(:,:,i)
                sVelo => sVeloKALE(:,:,i,:)
             end select

             ! Loop over the k and j-direction of this
             ! generalized i-face. Note that due to the usage of
             ! the pointers xx and xxOld an offset of +1 must be
             ! used in the coordinate arrays, because x and xOld
             ! originally start at 0 for the i, j and k indices.

             do k=1,kke
                do j=1,jje

                   ! Determine the dot product of sc and the normal
                   ! and divide by 4 deltaT to obtain the correct
                   ! value of the normal velocity.

                   sFace(j,k) = sVelo(j,k,1)*ss(j,k,1) &
                        + sVelo(j,k,2)*ss(j,k,2) &
                        + sVelo(j,k,3)*ss(j,k,3)

                enddo
             enddo
          enddo

       enddo loopDir
    endif testMoving

  end subroutine gridVelocitiesFineLevelPart2_block

  subroutine utauWF(rFilv)
    !
    !       utauWF substitutes the wall shear stress with values from a
    !       look-up table, if desired.
    !
    use constants
    use blockPointers, only : si, sj, sk, fw, rlv, d2wall, w, BCData, viscSubFace, &
         ie, je, ke, il, jl, kl, nViscBocos, BCFaceID
    use inputPhysics, only : wallFunctions
    use turbCurveFits, only : curveUpRe
    implicit none
    !
    !      Subroutine argument.
    !
    real(kind=realType), intent(in) :: rFilv
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, nn

    real(kind=realType) :: fact
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: rbar, ubar, vbar, wbar, vvx, vvy, vvz
    real(kind=realType) :: fmx, fmy, fmz, frhoe
    real(kind=realType) :: veln, velnx, velny, velnz, tx, ty, tz
    real(kind=realType) :: veltx, velty, veltz, veltmag
    real(kind=realType) :: txnx, txny, txnz, tynx, tyny, tynz
    real(kind=realType) :: tznx, tzny, tznz
    real(kind=realType) :: tautn, tauWall, utau, re

    real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
    real(kind=realType), dimension(:,:,:), pointer :: ss, rres
    real(kind=realType), dimension(:,:,:), pointer :: norm
    real(kind=realType), dimension(:,:),   pointer :: rrlv2, dd2Wall2

    ! Return immediately if no wall functions must be used.

    if(.not. wallFunctions) return

    ! Loop over the viscous subfaces of this block.

    viscSubfaces: do nn=1,nViscBocos

       ! Set a bunch of variables depending on the face id to make
       ! a generic treatment possible.

       select case (BCFaceID(nn))

       case (iMin)
          fact = -one

          ss  => si(1,:,:,:);  rres => fw(2,1:,1:,:)
          ww2 => w(2,1:,1:,:); ww1  => w(1,1:,1:,:)
          dd2Wall2 => d2Wall(2,:,:); rrlv2 => rlv(2,1:,1:)

          !===========================================================

       case (iMax)
          fact = one

          ss  => si(il,:,:,:);  rres => fw(il,1:,1:,:)
          ww2 => w(il,1:,1:,:); ww1  => w(ie,1:,1:,:)
          dd2Wall2 => d2Wall(il,:,:); rrlv2 => rlv(il,1:,1:)

          !===========================================================

       case (jMin)
          fact = -one

          ss  => sj(:,1,:,:);  rres => fw(1:,2,1:,:)
          ww2 => w(1:,2,1:,:); ww1  => w(1:,1,1:,:)
          dd2Wall2 => d2Wall(:,2,:); rrlv2 => rlv(1:,2,1:)

          !===========================================================

       case (jMax)
          fact = one

          ss  => sj(:,jl,:,:);  rres => fw(1:,jl,1:,:)
          ww2 => w(1:,jl,1:,:); ww1  => w(1:,je,1:,:)
          dd2Wall2 => d2Wall(:,jl,:); rrlv2 => rlv(1:,jl,1:)

          !===========================================================

       case (kMin)
          fact = -one

          ss  => sk(:,:,1,:);  rres => fw(1:,1:,2,:)
          ww2 => w(1:,1:,2,:); ww1  => w(1:,1:,1,:)
          dd2Wall2 => d2Wall(:,:,2); rrlv2 => rlv(1:,1:,2)

          !===========================================================

       case (kMax)
          fact = one

          ss  => sk(:,:,kl,:);  rres => fw(1:,1:,kl,:)
          ww2 => w(1:,1:,kl,:); ww1  => w(1:,1:,ke,:)
          dd2Wall2 => d2Wall(:,:,kl); rrlv2 => rlv(1:,1:,kl)

       end select

       ! Set the pointer for the unit outward normals.

       norm => BCData(nn)%norm

       ! Loop over the quadrilateral faces of the subface. Note
       ! that the nodal range of BCData must be used and not the
       ! cell range, because the latter may include the halo's in i
       ! and j-direction. The offset +1 is there, because inBeg and
       ! jnBeg refer to nodal ranges and not to cell ranges.

       do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
          do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

             ! Store the viscous stress tensor a bit easier.

             tauxx = viscSubface(nn)%tau(i,j,1)
             tauyy = viscSubface(nn)%tau(i,j,2)
             tauzz = viscSubface(nn)%tau(i,j,3)
             tauxy = viscSubface(nn)%tau(i,j,4)
             tauxz = viscSubface(nn)%tau(i,j,5)
             tauyz = viscSubface(nn)%tau(i,j,6)

             ! Compute the velocities at the wall face; these are only
             ! non-zero for moving a block. Also compute the density,
             ! which is needed to compute the wall shear stress via
             ! wall functions.

             rbar = half*(ww2(i,j,irho) + ww1(i,j,irho))
             ubar = half*(ww2(i,j,ivx)  + ww1(i,j,ivx))
             vbar = half*(ww2(i,j,ivy)  + ww1(i,j,ivy))
             wbar = half*(ww2(i,j,ivz)  + ww1(i,j,ivz))

             ! Compute the velocity difference between the internal cell
             ! and the wall.

             vvx = ww2(i,j,ivx) - ubar
             vvy = ww2(i,j,ivy) - vbar
             vvz = ww2(i,j,ivz) - wbar

             ! Compute the normal velocity of the internal cell.

             veln  = vvx*norm(i,j,1) + vvy*norm(i,j,2) + vvz*norm(i,j,3)
             velnx = veln*norm(i,j,1)
             velny = veln*norm(i,j,2)
             velnz = veln*norm(i,j,3)

             ! Compute the tangential velocity, its magnitude and its
             ! unit vector of the internal cell.

             veltx = vvx - velnx
             velty = vvy - velny
             veltz = vvz - velnz

             veltmag = max(eps,sqrt(veltx**2 + velty**2 + veltz**2))

             tx = veltx/veltmag
             ty = velty/veltmag
             tz = veltz/veltmag

             ! Compute some coefficients needed for the transformation
             ! between the cartesian frame and the frame defined by the
             ! tangential direction (tx,ty,tz) and the normal direction.
             ! The minus sign is present, because for this transformation
             ! the normal direction should be inward pointing and norm
             ! is outward pointing.

             txnx = -tx*norm(i,j,1)
             txny = -tx*norm(i,j,2)
             txnz = -tx*norm(i,j,3)

             tynx = -ty*norm(i,j,1)
             tyny = -ty*norm(i,j,2)
             tynz = -ty*norm(i,j,3)

             tznx = -tz*norm(i,j,1)
             tzny = -tz*norm(i,j,2)
             tznz = -tz*norm(i,j,3)

             ! Compute the tn component of the wall shear stress
             ! tensor. Normally this is the only nonzero shear
             ! stress component in the t-n frame.

             tautn = tauxx*txnx + tauyy*tyny + tauzz*tznz &
                  + tauxy*(txny + tynx)                  &
                  + tauxz*(txnz + tznx)                  &
                  + tauyz*(tynz + tzny)

             ! Compute the Reynolds number using the velocity, density,
             ! laminar viscosity and wall distance. Note that an offset
             ! of -1 must be used in dd2Wall2, because the original array
             ! d2Wall starts at 2.

             re = ww2(i,j,irho)*veltmag*dd2Wall2(i-1,j-1)/rrlv2(i,j)

             ! Determine the friction velocity from the table and
             ! compute the wall shear stress from it.

             utau    = veltmag/max(curveUpRe(re),eps)
             tauWall = rbar*utau*utau

             ! Compute the correction to the wall shear stress tautn and
             ! transform this correction back to the cartesian frame.
             ! Take rFilv into account, such that the correction to the
             ! stress tensor is computed correctly.

             tautn = rFilv*tauWall - tautn

             tauxx = two*tautn*txnx
             tauyy = two*tautn*tyny
             tauzz = two*tautn*tznz

             tauxy = tautn*(txny + tynx)
             tauxz = tautn*(txnz + tznx)
             tauyz = tautn*(tynz + tzny)

             ! Compute the correction to the viscous flux at the wall.

             fmx   = tauxx*ss(i,j,1) + tauxy*ss(i,j,2) &
                  + tauxz*ss(i,j,3)
             fmy   = tauxy*ss(i,j,1) + tauyy*ss(i,j,2) &
                  + tauyz*ss(i,j,3)
             fmz   = tauxz*ss(i,j,1) + tauyz*ss(i,j,2) &
                  + tauzz*ss(i,j,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*ss(i,j,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*ss(i,j,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*ss(i,j,3)

             ! Add them to the residual. Note that now the factor rFilv
             ! is already taken into account via tau. Fact is present to
             ! take inward/outward pointing normals into account

             rres(i,j,imx)   = rres(i,j,imx)   - fact*fmx
             rres(i,j,imy)   = rres(i,j,imy)   - fact*fmy
             rres(i,j,imz)   = rres(i,j,imz)   - fact*fmz
             rres(i,j,irhoE) = rres(i,j,irhoE) - fact*frhoE

             ! Store the friction velocity for later use.

             viscSubface(nn)%utau(i,j) = utau

             ! Also add the correction to the wall stress tensor.

             viscSubface(nn)%tau(i,j,1) = &
                  viscSubface(nn)%tau(i,j,1) + tauxx
             viscSubface(nn)%tau(i,j,2) = &
                  viscSubface(nn)%tau(i,j,2) + tauyy
             viscSubface(nn)%tau(i,j,3) = &
                  viscSubface(nn)%tau(i,j,3) + tauzz
             viscSubface(nn)%tau(i,j,4) = &
                  viscSubface(nn)%tau(i,j,4) + tauxy
             viscSubface(nn)%tau(i,j,5) = &
                  viscSubface(nn)%tau(i,j,5) + tauxz
             viscSubface(nn)%tau(i,j,6) = &
                  viscSubface(nn)%tau(i,j,6) + tauyz
          enddo
       enddo

    enddo viscSubfaces

  end subroutine utauWF



#ifndef USE_TAPENADE
  subroutine slipVelocitiesCoarseLevels(sps)
    !
    !       slipVelocitiesCoarseLevels determines the slip velocities
    !       for the given spectral solution starting from the known
    !       velocities on the finer level.
    !
    use constants
    use blockPointers
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, iiMax, jjMax
    integer(kind=intType) :: if1, if2, jf1, jf2
    integer(kind=intType) :: nLevels, level, levm1, nn, mm

    integer(kind=intType), dimension(:,:), pointer :: iFine, jFine

    real(kind=realType), dimension(:,:,:), pointer :: uSlip
    real(kind=realType), dimension(:,:,:), pointer :: uSlipFine

    ! Determine the number of multigrid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over coarser grid levels, where ground level is considered
    ! as the finest grid.

    levelLoop: do level=(groundLevel+1),nLevels

       ! Set levm1 for the finer level.

       levm1 = level - 1

       ! Loop over the number of local blocks.

       domains: do nn=1,nDom

          ! Set the pointers to the coarse block.

          call setPointers(nn, level, sps)

          ! Loop over the number of viscous subfaces.

          bocoLoop: do mm=1,nViscBocos

             ! Set the pointers for uSlip and uSlipFine to make the
             ! code more readable.

             uSlip     => BCData(mm)%uSlip
             uSlipFine => flowDoms(nn,levm1,sps)%BCData(mm)%uSlip

             ! Determine the grid face on which the subface is located
             ! and set some variables accordingly.

             select case (BCFaceID(mm))

             case (iMin,iMax)
                iiMax = jl; jjMax = kl
                iFine => mgJFine; jFine => mgKFine

             case (jMin,jMax)
                iiMax = il; jjMax = kl
                iFine => mgIFine; jFine => mgKFine

             case (kMin,kMax)
                iiMax = il; jjMax = jl
                iFine => mgIFine; jFine => mgJFine

             end select

             ! Loop over the number of faces of the viscous subface.
             ! First in the generalized j-direction.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd

                ! Determine the two children in this direction.
                ! Take care of the halo's, as this info is only
                ! available for owned cells.

                if(j < 2) then
                   jf1 = 1; jf2 = 1
                else if(j > jjMax) then
                   jf1 = jFine(jjMax,2) +1; jf2 = jf1
                else
                   jf1 = jFine(j,1); jf2 = jFine(j,2)
                endif

                ! Loop in the generalized i-direction.

                do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   ! Determine the two children in this direction.
                   ! Same story as in j-direction.

                   if(i < 2) then
                      if1 = 1; if2 = 1
                   else if(i > iiMax) then
                      if1 = iFine(iiMax,2) +1; if2 = if1
                   else
                      if1 = iFine(i,1); if2 = iFine(i,2)
                   endif

                   ! Average the fine grid velocities to the
                   ! coarse grid velocities.

                   uSlip(i,j,1) = fourth*(uSlipFine(if1,jf1,1) &
                        +         uSlipFine(if2,jf1,1) &
                        +         uSlipFine(if1,jf2,1) &
                        +         uSlipFine(if2,jf2,1))

                   uSlip(i,j,2) = fourth*(uSlipFine(if1,jf1,2) &
                        +         uSlipFine(if2,jf1,2) &
                        +         uSlipFine(if1,jf2,2) &
                        +         uSlipFine(if2,jf2,2))

                   uSlip(i,j,3) = fourth*(uSlipFine(if1,jf1,3) &
                        +         uSlipFine(if2,jf1,3) &
                        +         uSlipFine(if1,jf2,3) &
                        +         uSlipFine(if2,jf2,3))
                enddo
             enddo

          enddo bocoLoop
       enddo domains
    enddo levelLoop

  end subroutine slipVelocitiesCoarseLevels
#endif

  subroutine gridVelocitiesCoarseLevels(sps)
    !
    !       gridVelocitiesCoarseLevels computes the grid velocities for
    !       the cell centers and the normal grid velocities for the faces
    !       of moving blocks on the coarser grid levels. GroundLevel is
    !       considered the fine grid level.
    !
    use constants
    use blockPointers
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, level, levm1, nn, mm
    integer(kind=intType) :: i, j, k, iie, jje, kke
    integer(kind=intType) :: ii, ii1, jj, jj1, kk, kk1

    integer(kind=intType), dimension(:,:), pointer :: iFine, jFine
    integer(kind=intType), dimension(:,:), pointer :: kFine

    real(kind=realType) :: jjWeight, kkWeight, weight

    real(kind=realType), dimension(:), pointer :: jWeight, kWeight

    real(kind=realType), dimension(:,:),     pointer :: sFine, sFace
    real(kind=realType), dimension(:,:,:,:), pointer :: sf

    ! Loop over the number of coarse grid levels, starting at
    ! groundLevel+1,

    nLevels = ubound(flowDoms,2)
    levelLoop: do level=groundLevel+1,nLevels

       ! Loop over the number of local blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, level, sps)

          ! Check for a moving block. As it is possible that in a
          ! multidisicplinary environment additional grid velocities
          ! are set, the test should be done on addGridVelocities
          ! and not on blockIsMoving.

          testMoving: if( addGridVelocities ) then
             !
             !             Grid velocities of the cell centers, including the 1st
             !             level halo cells. These are determined by accumulating
             !             the fine grid values. At the end the internal halo's are
             !             communicated to obtain the correct values.
             !
             levm1 = level - 1

             ! Set the pointer sf to the cell velocities on the fine mesh.

             sf => flowDoms(nn,levm1,sps)%s

             ! Loop over the cells, including the 1st level halo's.
             ! The indices kk, kk1 contain the corresponding fine
             ! grid indices in k-direction. Idem for jj, jj1, ii
             ! and ii1.

             do k=1,ke
                if(k == 1) then
                   kk = 1; kk1 = 1
                else if(k == ke) then
                   kk = flowDoms(nn,levm1,sps)%ke; kk1 = kk
                else
                   kk = mgKFine(k,1); kk1 = mgKFine(k,2)
                endif

                do j=1,je
                   if(j == 1) then
                      jj = 1; jj1 = 1
                   else if(j == je) then
                      jj = flowDoms(nn,levm1,sps)%je; jj1 = jj
                   else
                      jj = mgJFine(j,1); jj1 = mgJFine(j,2)
                   endif

                   do i=1,ie
                      if(i == 1) then
                         ii = 1; ii1 = 1
                      else if(i == ie) then
                         ii = flowDoms(nn,levm1,sps)%ie; ii1 = ii
                      else
                         ii = mgIFine(i,1); ii1 = mgIFine(i,2)
                      endif

                      ! Determine the coarse grid velocity by
                      ! averaging the fine grid values.

                      s(i,j,k,1) = (sf(ii1,jj1,kk1,1) + sf(ii,jj1,kk1,1)  &
                           +  sf(ii1,jj, kk1,1) + sf(ii,jj, kk1,1)  &
                           +  sf(ii1,jj1,kk, 1) + sf(ii,jj1,kk, 1)  &
                           +  sf(ii1,jj, kk, 1) + sf(ii,jj, kk, 1)) &
                           * eighth
                      s(i,j,k,2) = (sf(ii1,jj1,kk1,2) + sf(ii,jj1,kk1,2)  &
                           +  sf(ii1,jj, kk1,2) + sf(ii,jj, kk1,2)  &
                           +  sf(ii1,jj1,kk, 2) + sf(ii,jj1,kk, 2)  &
                           +  sf(ii1,jj, kk, 2) + sf(ii,jj, kk, 2)) &
                           * eighth
                      s(i,j,k,3) = (sf(ii1,jj1,kk1,3) + sf(ii,jj1,kk1,3)  &
                           +  sf(ii1,jj, kk1,3) + sf(ii,jj, kk1,3)  &
                           +  sf(ii1,jj1,kk, 3) + sf(ii,jj1,kk, 3)  &
                           +  sf(ii1,jj, kk, 3) + sf(ii,jj, kk, 3)) &
                           * eighth
                   enddo
                enddo
             enddo
             !
             !             Normal grid velocities of the faces.
             !
             ! Loop over the three directions.

             loopCoarseDir: do mm=1,3

                ! Set some values depending on the situation.

                select case (mm)

                case (1_intType)       ! Normals in i-direction
                   iie = ie; jje = je; kke = ke
                   iFine => mgIFine; jFine => mgJFine; kFine => mgKFine
                   jWeight => mgJWeight; kWeight => mgKWeight

                case (2_intType)       ! Normals in j-direction
                   iie = je; jje = ie; kke = ke
                   iFine => mgJFine; jFine => mgIFine; kFine => mgKFine
                   jWeight => mgIWeight; kWeight => mgKWeight

                case (3_intType)       ! Normals in k-direction
                   iie = ke; jje = ie; kke = je
                   iFine => mgKFine; jFine => mgIFine; kFine => mgJFine
                   jWeight => mgIWeight; kWeight => mgJWeight

                end select
                !
                !               Normal grid velocities in generalized i-direction.
                !               mm == 1: i-direction
                !               mm == 2: j-direction
                !               mm == 3: k-direction
                !
                do i=0,iie

                   ! Determine the i-index of the corresponding plane on
                   ! the fine grid. Note that halo planes are not entirely
                   ! correct. This is not really a problem.

                   if(i < 2) then
                      ii = i
                   else if(i < iie) then
                      ii = iFine(i,2)
                   else
                      ii = iFine(iie-1,2) + 1
                   endif

                   ! Set the pointers for sFine and sFace, which will
                   ! contain the mesh velocities for this particular
                   ! plane. The index depends on the value of mm.

                   select case (mm)
                   case (1_intType)
                      sFine => flowDoms(nn,levm1,sps)%sFaceI(ii,:,:)
                      sFace => sFaceI(i,:,:)
                   case (2_intType)
                      sFine => flowDoms(nn,levm1,sps)%sFaceJ(:,ii,:)
                      sFace => sFaceJ(:,i,:)
                   case (3_intType)
                      sFine => flowDoms(nn,levm1,sps)%sFaceK(:,:,ii)
                      sFace => sFaceK(:,:,i)
                   end select

                   ! Loop over the k and j faces for this general i-plane.
                   ! Again the halo's are not entirely correct. Kk, kk1,
                   ! jj and jj1 are the children in k and j-direction
                   ! respectively.

                   do k=1,kke
                      if(k == 1) then
                         kk = 1; kk1 = 1
                         kkWeight = kWeight(2)
                      else if(k == kke) then
                         kk = kFine(kke-1,2) + 1; kk1 = kk
                         kkWeight = kWeight(kke-1)
                      else
                         kk = kFine(k,1); kk1 = kFine(k,2)
                         kWeight = kWeight(k)
                      endif

                      do j=1,jje
                         if(j == 1) then
                            jj = 1; jj1 = 1
                            jjWeight = jWeight(2)
                         else if(j == jje) then
                            jj = jFine(jje-1,2) + 1; jj1 = jj
                            jWeight = jWeight(jje-1)
                         else
                            jj = jFine(j,1); jj1 = jFine(j,2)
                            jjWeight = jWeight(j)
                         endif

                         ! Determine the coarse grid normal velocity.
                         ! Take the averaging weight into account; for
                         ! a normal coarsening this weight is 1.0.

                         weight = kkWeight*jjWeight
                         sFace(j,k) = weight*(sFine(jj1,kk1) &
                              +         sFine(jj ,kk1) &
                              +         sFine(jj1,kk)  &
                              +         sFine(jj ,kk))
                      enddo
                   enddo
                enddo

             enddo loopCoarseDir
          endif testMoving
       enddo domains

       ! Exchange the cell centered velocities.

       call exchangeCellGridVelocities(level,sps)

    enddo levelLoop

  end subroutine gridVelocitiesCoarseLevels

  !      ==================================================================

  subroutine exchangeCellGridVelocities(level,sps)
    !
    !       exchangeCellGridVelocities exchanges the grid velocities in
    !       the cell centers for the given grid level and spectral
    !       solution.
    !
    use constants
    use block
    use communication
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

    real(kind=realType) :: alp
    real(kind=realType), dimension(3) :: vv

    !       The 1 to 1 communication.
    !
    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPatternCell_1st(level)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPatternCell_1st(level)%sendProc(i)
       size   = 3*commPatternCell_1st(level)%nSend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPatternCell_1st(level)%nSend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPatternCell_1st(level)%sendList(i)%block(j)
          i1 = commPatternCell_1st(level)%sendList(i)%indices(j,1)
          j1 = commPatternCell_1st(level)%sendList(i)%indices(j,2)
          k1 = commPatternCell_1st(level)%sendList(i)%indices(j,3)

          ! Store the grid velocities in sendBuffer, if they
          ! are allocated. Otherwise they are simply zero.

          if( flowDoms(d1,level,sps)%addGridVelocities ) then
             sendBuffer(jj)   = flowDoms(d1,level,sps)%s(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDoms(d1,level,sps)%s(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDoms(d1,level,sps)%s(i1,j1,k1,3)
          else
             sendBuffer(jj)   = zero
             sendBuffer(jj+1) = zero
             sendBuffer(jj+2) = zero
          endif

          jj = jj + 3
       enddo

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj
    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPatternCell_1st(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPatternCell_1st(level)%recvProc(i)
       size   = 3*commPatternCell_1st(level)%nRecv(i)

       ! Post the receive.

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size
    enddo receives

    ! Copy the local data.

    localCopy: do i=1,internalCell_1st(level)%nCopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internalCell_1st(level)%donorBlock(i)
       i1 = internalCell_1st(level)%donorIndices(i,1)
       j1 = internalCell_1st(level)%donorIndices(i,2)
       k1 = internalCell_1st(level)%donorIndices(i,3)

       ! Idem for the halo's.

       d2 = internalCell_1st(level)%haloBlock(i)
       i2 = internalCell_1st(level)%haloIndices(i,1)
       j2 = internalCell_1st(level)%haloIndices(i,2)
       k2 = internalCell_1st(level)%haloIndices(i,3)

       ! Copy the grid velocities, if they are both allocated.
       ! Otherwise they are either set to zero or nothing is done.

       if( flowDoms(d2,level,sps)%addGridVelocities ) then
          if( flowDoms(d1,level,sps)%addGridVelocities ) then
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = &
                  flowDoms(d1,level,sps)%s(i1,j1,k1,1)
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = &
                  flowDoms(d1,level,sps)%s(i1,j1,k1,2)
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = &
                  flowDoms(d1,level,sps)%s(i1,j1,k1,3)
          else
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = zero
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = zero
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = zero
          endif
       endif

    enddo localCopy

    ! Correct the periodic halo's of the internal communication
    ! pattern, if present.

    if(internalCell_1st(level)%nPeriodic > 0) &
         call correctPeriodicGridVel(level, sps,                        &
         internalCell_1st(level)%nPeriodic, &
         internalCell_1st(level)%periodicData)

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPatternCell_1st(level)%nProcRecv
    completeRecvs: do i=1,commPatternCell_1st(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = 3*commPatternCell_1st(level)%nRecvCum(ii-1)
       do j=1,commPatternCell_1st(level)%nRecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPatternCell_1st(level)%recvList(ii)%block(j)
          i2 = commPatternCell_1st(level)%recvList(ii)%indices(j,1)
          j2 = commPatternCell_1st(level)%recvList(ii)%indices(j,2)
          k2 = commPatternCell_1st(level)%recvList(ii)%indices(j,3)

          ! Copy the grid velocities from recvBuffer if they are
          ! both allocated.

          if( flowDoms(d2,level,sps)%addGridVelocities ) then
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = recvBuffer(jj+1)
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = recvBuffer(jj+2)
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = recvBuffer(jj+3)
          endif

          jj = jj + 3
       enddo

    enddo completeRecvs

    ! Correct the periodic halo's of the external communication
    ! pattern, if present.

    if(commPatternCell_1st(level)%nPeriodic > 0) &
         call correctPeriodicGridVel(level, sps,                           &
         commPatternCell_1st(level)%nPeriodic, &
         commPatternCell_1st(level)%periodicData)

    ! Complete the nonblocking sends.

    size = commPatternCell_1st(level)%nProcSend
    do i=1,commPatternCell_1st(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

  end subroutine exchangeCellGridVelocities

  !      ==================================================================

  subroutine correctPeriodicGridVel(level, sps, nPeriodic, &
       periodicData)
    !
    !       correctPeriodicGridVel applies the periodic transformation
    !       to the grid velocities of the cell halo's in periodicData.
    !
    use block
    use communication
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps, nPeriodic
    type(periodicDataType), dimension(:), pointer :: periodicData
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, ii, i, j, k
    real(kind=realType)   :: vx, vy, vz

    real(kind=realType), dimension(3,3) :: rotMatrix

    ! Loop over the number of periodic transformations.

    do nn=1,nPeriodic

       ! Store the rotation matrix a bit easier.

       rotMatrix = periodicData(nn)%rotMatrix

       ! Loop over the number of halo cells for this transformation.

       do ii=1,periodicData(nn)%nHalos

          ! Store the block and the indices a bit easier.

          mm = periodicData(nn)%block(ii)
          i  = periodicData(nn)%indices(ii,1)
          j  = periodicData(nn)%indices(ii,2)
          k  = periodicData(nn)%indices(ii,3)

          ! Check if the grid velocities have been allocated.

          if( flowDoms(mm,level,sps)%addGridVelocities ) then

             ! Store the original velocities in vx, vy, vz.

             vx = flowDoms(mm,level,sps)%s(i,j,k,1)
             vy = flowDoms(mm,level,sps)%s(i,j,k,2)
             vz = flowDoms(mm,level,sps)%s(i,j,k,3)

             ! Compute the new velocity vector.

             flowDoms(mm,level,sps)%s(i,j,k,1) = rotMatrix(1,1)*vx &
                  + rotMatrix(1,2)*vy &
                  + rotMatrix(1,3)*vz
             flowDoms(mm,level,sps)%s(i,j,k,2) = rotMatrix(2,1)*vx &
                  + rotMatrix(2,2)*vy &
                  + rotMatrix(2,3)*vz
             flowDoms(mm,level,sps)%s(i,j,k,3) = rotMatrix(3,1)*vx &
                  + rotMatrix(3,2)*vy &
                  + rotMatrix(3,3)*vz

          endif
       enddo

    enddo

  end subroutine correctPeriodicGridVel



  !      ==================================================================


#endif
end module solverUtils

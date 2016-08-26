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
!
!       File:          timeStep.f90                                    
!       Author:        Edwin van der Weide                             
!       Starting date: 03-17-2003                                      
!       Last modified: 06-28-2005                                      
!
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
                 
                 dtl(i,j,k) = ri + rj + rk
                 
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

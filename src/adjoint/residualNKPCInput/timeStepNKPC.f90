!
!      ******************************************************************
!      *                                                                *
!      * File:          timeStepAdj.f90                                 *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 08-05-2009                                      *
!      * Last modified: 08-05-2009                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine timeStepNKPC(onlyRadii,wAdj,pAdj,siAdj, sjAdj, skAdj,&
     sFaceIAdj,sFaceJAdj,sFaceKAdj,volAdj,radIAdj,radJAdj,radKAdj,&
     iCell, jCell, kCell,pInfCorrAdj,rhoInfAdj,nn,level,sps,sps2)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * timeStep computes the time step, or more precisely the time    *
  !      * step divided by the volume per unit CFL, in the owned cells.   *
  !      * However, for the artificial dissipation schemes, the spectral  *
  !      * radIi in the halo's are needed. Therefore the loop is taken    *
  !      * over the the first level of halo cells. The spectral radIi are *
  !      * stored and possibly modified for high aspect ratio cells.      *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use constants
  use flowVarRefState
  use inputDiscretization
  use inputIteration
  use inputPhysics
  use inputTimeSpectral
  use iteration
  use section
  implicit none
  !
  !      Subroutine argument.
  !
  logical, intent(in) :: onlyRadii
  integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,level,sps,sps2
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: pAdj

  real(kind=realType), dimension(-1:1,-1:1,-1:1,nTimeIntervalsSpectral) :: radIAdj,radJAdj,radKAdj


  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType) ::pInfCorrAdj,rhoinfAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral):: volAdj
  !
  !      Local parameters.
  !
  real(kind=realType), parameter :: b = 2.0_realType
  !
  !      Local variables.
  !
  !integer(kind=intType) :: sps, nn, i, j, k
  integer(kind=intType) :: i, j, k,ii,jj,kk

  real(kind=realType) :: plim, rlim, clim2
  real(kind=realType) :: ux, uy, uz, cc2, qs, sx, sy, sz, rmu
  real(kind=realType) :: ri, rj, rk, rij, rjk, rki
  real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
  real(kind=realType) :: sFace, tmp
  real(kind=realType) :: dtladj

  logical :: radiiNeeded
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
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

  plim  = 0.001_realType*pInfCorrAdj
  rlim  = 0.001_realType*rhoInfAdj
  clim2 = 0.000001_realType*gammaInf*pInfCorrAdj/rhoInfAdj

  ! Initialize sFace to zero. This value will be used if the
  ! block is not moving.

  sFace = zero
  !
  !          **************************************************************
  !          *                                                            *
  !          * Inviscid contribution, depending on the preconditioner.    *
  !          * Compute the cell centered values of the spectral radii.    *
  !          *                                                            *
  !          **************************************************************
  !
  !print *,'precond:',precond
  select case (precond)

  case (noPrecond)

     ! No preconditioner. Simply the standard spectral radius.
     ! Loop over the cells, including the first level halo.
     do k=-1,1
        do j=-1,1
           do i=-1,1
              ii = icell+i
              jj = jcell+j
              kk = kcell+k
              ! Compute the velocities and speed of sound squared.
           
              ux  = wAdj(i,j,k,ivx,sps2)
              uy  = wAdj(i,j,k,ivy,sps2)
              uz  = wAdj(i,j,k,ivz,sps2)
              cc2 = gamma(ii,jj,kk)*pAdj(i,j,k,sps2)/wAdj(i,j,k,irho,sps2)
              cc2 = max(cc2,clim2)

              ! Set the dot product of the grid velocity and the
              ! normal in i-direction for a moving face. To avoid
              ! a number of multiplications by 0.5 simply the sum
              ! is taken.

              if( addGridVelocities ) &
                   sFace = sFaceIAdj(i-1,j,k,sps2) + sFaceIAdj(i,j,k,sps2)

              ! Spectral radius in i-direction.

              sx = siAdj(i-1,j,k,1,sps2) + siAdj(i,j,k,1,sps2)
              sy = siAdj(i-1,j,k,2,sps2) + siAdj(i,j,k,2,sps2)
              sz = siAdj(i-1,j,k,3,sps2) + siAdj(i,j,k,3,sps2)

              qs = ux*sx + uy*sy + uz*sz - sFace

              if( sx**2>zero .or. sy**2>zero .or. sz**2>zero)then
                 radiAdj(i,j,k,sps2) = half*(abs(qs) &
                      +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))
              else
                 radiAdj(i,j,k,sps2) = half*(abs(qs))
              endif

              ! The grid velocity in j-direction.

              if( addGridVelocities ) &
                   sFace = sFaceJAdj(i,j-1,k,sps2) + sFaceJAdj(i,j,k,sps2)

              ! Spectral radius in j-direction.

              sx = sjAdj(i,j-1,k,1,sps2) + sjAdj(i,j,k,1,sps2)
              sy = sjAdj(i,j-1,k,2,sps2) + sjAdj(i,j,k,2,sps2)
              sz = sjAdj(i,j-1,k,3,sps2) + sjAdj(i,j,k,3,sps2)

              qs = ux*sx + uy*sy + uz*sz - sFace
           
              if( sx**2>zero .or. sy**2>zero .or. sz**2>zero)then
                 radJAdj(i,j,k,sps2) = half*(abs(qs) &
                      +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))
              else
                 radJAdj(i,j,k,sps2) = half*(abs(qs))
              endif
              ! The grid velocity in k-direction.

              if( addGridVelocities ) &
                   sFace = sFaceKAdj(i,j,k-1,sps2) + sFaceKAdj(i,j,k,sps2)

              ! Spectral radius in k-direction.

              sx = skAdj(i,j,k-1,1,sps2) + skAdj(i,j,k,1,sps2)
              sy = skAdj(i,j,k-1,2,sps2) + skAdj(i,j,k,2,sps2)
              sz = skAdj(i,j,k-1,3,sps2) + skAdj(i,j,k,3,sps2)

              qs = ux*sx + uy*sy + uz*sz - sFace

              if( sx**2>zero .or. sy**2>zero .or. sz**2>zero)then

                 radKAdj(i,j,k,sps2) = half*(abs(qs) &
                      +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))
              else
                 radKAdj(i,j,k,sps2) = half*(abs(qs) )
              endif

              ! Compute the inviscid contribution to the time step.

              dtlAdj = radiAdj(i,j,k,sps2) + radJAdj(i,j,k,sps2) + radKAdj(i,j,k,sps2)

           enddo
        enddo
     enddo
  case (Turkel)
     call terminate("timeStep","Turkel preconditioner not implemented yet")
  case (ChoiMerkle)
     call terminate("timeStep", "choi merkle preconditioner not implemented yet")
  end select
  !
  !          **************************************************************
  !          *                                                            *
  !          * Adapt the spectral radii if directional scaling must be    *
  !          * applied.                                                   *
  !          *                                                            *
  !          **************************************************************
  !
  if(dirScaling .and. currentLevel <= groundLevel) then
     do k=-1,1
        do j=-1,1
           do i=-1,1

              ! Avoid division by zero by clipping radi, radJ and
              ! radK.
              ri = max(radiAdj(i,j,k,sps2),eps)
              rj = max(radJAdj(i,j,k,sps2),eps)
              rk = max(radKAdj(i,j,k,sps2),eps)

              ! Compute the scaling in the three coordinate
              ! directions.

              rij = (ri/rj)**adis
              rjk = (rj/rk)**adis
              rki = (rk/ri)**adis

              ! Create the scaled versions of the aspect ratios.
              ! Note that the multiplication is done with radi, radJ
              ! and radK, such that the influence of the clipping
              ! is negligible.

              !   radi(i,j,k) = third*radi(i,j,k)*(one + one/rij + rki)
              !   radJ(i,j,k) = third*radJ(i,j,k)*(one + one/rjk + rij)
              !   radK(i,j,k) = third*radK(i,j,k)*(one + one/rki + rjk)

              radiAdj(i,j,k,sps2) = radiAdj(i,j,k,sps2)*(one + one/rij + rki)
              radJAdj(i,j,k,sps2) = radJAdj(i,j,k,sps2)*(one + one/rjk + rij)
              radKAdj(i,j,k,sps2) = radKAdj(i,j,k,sps2)*(one + one/rki + rjk)

           enddo
        enddo
     enddo

  endif

  ! The rest of this file can be skipped if only the spectral
  ! radii need to be computed.

  testRadiiOnly: if(.not. onlyRadii) then

     ! The viscous contribution, if needed.

     viscousTerm: if( viscous ) then

        print *,'Viscous not yet implemented'
        stop

     endif viscousTerm

     ! For the spectral mode an additional term term must be
     ! taken into account, which corresponds to the contribution
     ! of the highest frequency.

     if(equationMode == timeSpectral) then

        print *,'time spectral not yet implemented'
        tmp = nTimeIntervalsSpectral*pi*timeRef &
             / sections(sectionID)%timePeriod

        ! Loop over the owned cell centers and add the term.

        dtlAdj = dtlAdj + tmp*volAdj(sps2)

     endif

     ! Currently the inverse of dt/vol is stored in dtl. Invert
     ! this value such that the time step per unit cfl number is
     ! stored and correct in cases of high gradients.
     do k=-1,1
        do j=-1,1
           do i=-1,1
              dpi = abs(pAdj(i+1,j,k,sps2) - two*pAdj(i,j,k,sps2) + pAdj(i-1,j,k,sps2)) &
                   /    (pAdj(i+1,j,k,sps2) + two*pAdj(i,j,k,sps2) + pAdj(i-1,j,k,sps2) + plim)
              dpj = abs(pAdj(i,j+1,k,sps2) - two*pAdj(i,j,k,sps2) + pAdj(i,j-1,k,sps2)) &
                   /    (pAdj(i,j+1,k,sps2) + two*pAdj(i,j,k,sps2) + pAdj(i,j-1,k,sps2) + plim)
              dpk = abs(pAdj(i,j,k+1,sps2) - two*pAdj(i,j,k,sps2) + pAdj(i,j,k-1,sps2)) &
                   /    (pAdj(i,j,k+1,sps2) + two*pAdj(i,j,k,sps2) + pAdj(i,j,k-1,sps2) + plim)
              rfl = one/(one + b*(dpi  +dpj  +dpk))

              dtlAdj = rfl/dtlAdj
           enddo
        enddo
     enddo

  endif testRadiiOnly
end subroutine timeStepNKPC

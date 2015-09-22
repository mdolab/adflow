subroutine residual
  !
  ! Shell function to call residual_block on all blocks
  !
  use blockPointers
  use constants
  use inputTimeSpectral
  use Iteration

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

        call setPointers(nn, currentLevel, sps)

        call residual_block

     end do domains

  end do spectralLoop

end subroutine residual

!
!      ******************************************************************
!      *                                                                *
!      * File:          residual.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher (blanking)   *
!      * Starting date: 03-15-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine residual_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * residual computes the residual of the mean flow equations on   *
  !      * the current MG level.                                          *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use cgnsGrid
  use flowVarRefState
  use inputIteration
  use inputDiscretization
  use inputTimeSpectral
  use inputUnsteady ! Added by HDN
  use iteration
  use inputAdjoint
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: discr
  integer(kind=intType) :: i, j, k, l
  integer(kind=intType) :: iale, jale, kale, lale, male ! For loops of ALE
  real(kind=realType), parameter :: K1 = 1.05_realType
  real(kind=realType), parameter :: K2 = 0.6_realType ! Random given number
  real(kind=realType), parameter :: M0 = 0.2_realType ! Mach number preconditioner activation
  real(kind=realType), parameter :: alpha = 0_realType
  real(kind=realType), parameter :: delta = 0_realType
  !real(kind=realType), parameter :: hinf = 2_realType ! Test phase 
  real(kind=realType), parameter :: Cpres = 4.18_realType ! Test phase
  real(kind=realType), parameter :: TEMP= 297.15_realType

  !
  !     Local variables
  !
  real(kind=realType) :: K3, h, velXrho, velYrho, velZrho,SoS,hinf
  real(kind=realType) :: resM,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35
  real(kind=realType) :: A41,A42,A43,A44,A45,A51,A52,A53,A54,A55,B11,B12,B13,B14,B15
  real(kind=realType) :: B21,B22,B23,B24,B25,B31,B32,B33,B34,B35
  real(kind=realType) :: B41,B42,B43,B44,B45,B51,B52,B53,B54,B55
  real(kind=realType) :: rhoHdash, betaMr2
  real(kind=realType) :: G, q
  real(kind=realType) :: b1, b2, b3, b4, b5
  real(kind=realType) :: dwo(nwf)
  logical :: fineGrid
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Add the source terms from the level 0 cooling model.
#ifndef USE_TAPENADE
  call level0CoolingModel
#endif
  ! Set the value of rFil, which controls the fraction of the old
  ! dissipation residual to be used. This is only for the runge-kutta
  ! schemes; for other smoothers rFil is simply set to 1.0.
  ! Note the index rkStage+1 for cdisRK. The reason is that the
  ! residual computation is performed before rkStage is incremented.

  if(smoother == RungeKutta) then
     rFil = cdisRK(rkStage+1)
  else
     rFil = one
  endif

  ! Initialize the local arrays to monitor the massflows to zero.
#ifndef USE_TAPENADE
  massFlowFamilyInv  = zero
  massFlowFamilyDiss = zero
#endif
  ! Set the value of the discretization, depending on the grid level,
  ! and the logical fineGrid, which indicates whether or not this
  ! is the finest grid level of the current mg cycle.

  discr = spaceDiscrCoarse
  if(currentLevel == 1) discr = spaceDiscr

  fineGrid = .false.
  if(currentLevel == groundLevel) fineGrid = .true.


  ! ===========================================================
  !
  ! Assuming ALE has nothing to do with MG
  ! The geometric data will be interpolated if in MD mode
  !
  ! ===========================================================
#ifndef USE_TAPENADE
  call interpLevelALE_block
#endif
  ! ===========================================================
  !
  ! The fluxes are calculated as usual
  !
  ! ===========================================================

  call inviscidCentralFlux

  select case (discr)

  case (dissScalar) ! Standard scalar dissipation scheme.

     if( fineGrid) then 
        if (.not. lumpedDiss) then
           call inviscidDissFluxScalar
        else
           call inviscidDissFluxScalarApprox
        end if
     else
#ifndef USE_TAPENADE
        call inviscidDissFluxScalarCoarse
#endif
     endif

     !===========================================================

  case (dissMatrix) ! Matrix dissipation scheme.

     if( fineGrid ) then
        if (.not. lumpedDiss) then 
           call inviscidDissFluxMatrix
        else
           call inviscidDissFluxMatrixApprox
        end if
     else
#ifndef USE_TAPENADE
        call inviscidDissFluxMatrixCoarse
#endif
     endif

     !===========================================================

  case (dissCusp) ! Cusp dissipation scheme.
#ifndef USE_TAPENADE
     if( fineGrid ) then
        call inviscidDissFluxCusp
     else
        call inviscidDissFluxCuspCoarse
     endif
#endif
     !===========================================================

  case (upwind) ! Dissipation via an upwind scheme.

     call inviscidUpwindFlux(fineGrid)

  end select

  !-------------------------------------------------------
  ! Lastly, recover the old s[I,J,K], sFace[I,J,K]
  ! This shall be done before difussive and source terms
  ! are computed.
  !-------------------------------------------------------
#ifndef USE_TAPENADE
  call recoverLevelALE_block
#endif


  if( viscous ) then 
     ! Only compute viscous fluxes if rFil > 0
     if(abs(rFil) > thresholdReal) then 
        ! not lumpedDiss means it isn't the PC...call the vicousFlux
        if (.not. lumpedDiss) then 
           call computeSpeedOfSoundSquared
           call allNodalGradients
           call viscousFlux
        else
           ! This is a PC calc...only include viscous fluxes if viscPC
           ! is used
           if (viscPC) then
              call computeSpeedOfSoundSquared
              call allNodalGradients
              call viscousFlux 
           else
              call viscousFluxApprox
           end if
        end if
     end if
  end if

  !===========================================================

  ! Add the dissipative and possibly viscous fluxes to the
  ! Euler fluxes. Loop over the owned cells and add fw to dw.
  ! Also multiply by iblank so that no updates occur in holes
  if ( lowspeedpreconditioner ) then
     do k=2,kl
        do j=2,jl
           do i=2,il
              !    Compute speed of sound
              SoS = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))

              ! Coompute velocities without rho from state vector
              velXrho = w(i,j,k,ivx)
              velYrho = w(i,j,k,ivy)
              velZrho = w(i,j,k,ivz)

              q = (velXrho**2 + velYrho**2 + velZrho**2)
              resM=sqrt(q)/SoS

              !
              !    Compute K3
              K3 = K1 * ( 1 + ((1-K1*M0**2)*resM**2)/(K1*M0**4) )
              !    Compute BetaMr2
              betaMr2 = min( max( K3*(velXrho**2 + velYrho**2  &
                   + velZrho**2), ((K2)*(wInf(ivx)**2 &
                   + wInf(ivy)**2 + wInf(ivz)**2))) , SoS**2 )

              A11=  (betaMr2)*(1/SoS**4)
              A12 = zero
              A13 = zero
              A14 = zero
              A15 = (-betaMr2)/SoS**4 

              A21 =one*velXrho/SoS**2 
              A22 = one*w(i,j,k,irho)
              A23 = zero
              A24 = zero
              A25 = one*(-velXrho)/SoS**2

              A31 = one*velYrho/SoS**2
              A32 = zero
              A33 = one*w(i,j,k,irho)
              A34 = zero 
              A35 = one*(-velYrho)/SoS**2

              A41 = one*velZrho/SoS**2 
              A42 = zero
              A43 = zero
              A44 = one*w(i,j,k,irho) 
              A45 = zero + one *(-velZrho)/SoS**2 

              A51=  one*((1/(gamma(i,j,k)-1))+(resM**2)/2)
              A52 = one * w(i,j,k,irho)*velXrho
              A53 = one * w(i,j,k,irho)*velYrho
              A54 = one * w(i,j,k,irho)*velzrho
              A55 = one * ((-(resM**2))/2)

              B11 = A11*(gamma(i,j,k)-1)*q/2+A12*(-velXrho)&
                   /w(i,j,k,irho)+A13*(-velYrho)/w(i,j,k,irho)+A14*(-velZrho)/w(i,j,k,irho)&
                   + A15* (((gamma(i,j,k)-1)*q/2)-SoS**2)
              B12 = A11*(1-gamma(i,j,k))*velXrho+A12*1/w(i,j,k,irho)&
                   +A15*(1-gamma(i,j,k))*velXrho
              B13 = A11*(1-gamma(i,j,k))*velYrho+A13&
                   /w(i,j,k,irho)+A15*(1-gamma(i,j,k))*velYrho
              B14 = A11*(1-gamma(i,j,k))*velZrho&
                   +A14/w(i,j,k,irho)+A15*(1-gamma(i,j,k))*velZrho
              B15 = A11*(gamma(i,j,k)-1)+A15*(gamma(i,j,k)-1)

              B21 = A21*(gamma(i,j,k)-1)*q/2+A22*(-velXrho)&
                   /w(i,j,k,irho)+A23*(-velYrho)/w(i,j,k,irho)+A24*(-velZrho)&
                   /w(i,j,k,irho)+ A25* (((gamma(i,j,k)-1)*q/2)-SoS**2)
              B22 = A21*(1-gamma(i,j,k))*velXrho+A22&
                   /w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velXrho
              B23 = A21*(1-gamma(i,j,k))*velYrho&
                   +A23*1/w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velYrho
              B24 = A21*(1-gamma(i,j,k))*velZrho&
                   +A24*1/w(i,j,k,irho)+A25*(1-gamma(i,j,k))*velZrho
              B25 = A21*(gamma(i,j,k)-1)+A25*(gamma(i,j,k)-1)

              B31 = A31*(gamma(i,j,k)-1)*q/2+A32*(-velXrho)&
                   /w(i,j,k,irho)+A33*(-velYrho)/w(i,j,k,irho)+A34*(-velZrho)/w(i,j,k,irho)&
                   + A35* (((gamma(i,j,k)-1)*q/2)-SoS**2)
              B32 = A31*(1-gamma(i,j,k))*velXrho+A32&
                   /w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velXrho
              B33 = A31*(1-gamma(i,j,k))*velYrho&
                   +A33*1/w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velYrho
              B34 = A31*(1-gamma(i,j,k))*velZrho&
                   +A34*1/w(i,j,k,irho)+A35*(1-gamma(i,j,k))*velZrho
              B35 = A31*(gamma(i,j,k)-1)+A35*(gamma(i,j,k)-1)

              B41 = A41*(gamma(i,j,k)-1)*q/2+A42*(-velXrho)&
                   /w(i,j,k,irho)+A43*(-velYrho)/w(i,j,k,irho)+A44*(-velZrho) &
                   /w(i,j,k,irho)+ A45* (((gamma(i,j,k)-1)*q/2)-SoS**2)
              B42 = A41*(1-gamma(i,j,k))*velXrho+A42&
                   /w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velXrho
              B43 = A41*(1-gamma(i,j,k))*velYrho&
                   +A43*1/w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velYrho
              B44 = A41*(1-gamma(i,j,k))*velZrho&
                   +A44*1/w(i,j,k,irho)+A45*(1-gamma(i,j,k))*velZrho
              B45 = A41*(gamma(i,j,k)-1)+A45*(gamma(i,j,k)-1)

              B51 = A51*(gamma(i,j,k)-1)*q/2+A52*(-velXrho)&
                   /w(i,j,k,irho)+A53*(-velYrho)/w(i,j,k,irho)+A54*(-velZrho) &
                   /w(i,j,k,irho)+ A55* (((gamma(i,j,k)-1)*q/2)-SoS**2)
              B52 = A51*(1-gamma(i,j,k))*velXrho+A52&
                   /w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velXrho
              B53 = A51*(1-gamma(i,j,k))*velYrho&
                   +A53*1/w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velYrho
              B54 = A51*(1-gamma(i,j,k))*velZrho&
                   +A54*1/w(i,j,k,irho)+A55*(1-gamma(i,j,k))*velZrho
              B55 = A51*(gamma(i,j,k)-1)+A55*(gamma(i,j,k)-1)

              ! dwo is the orginal redisual
              do l=1,nwf
                 dwo(l) = (dw(i,j,k,l) + fw(i,j,k,l))* max(real(iblank(i,j,k), realType), zero)
              end do

              dw(i,j,k,1)=B11*dwo(1) + B12*dwo(2)+ B13*dwo(3) + B14*dwo(4) + B15*dwo(5)
              dw(i,j,k,2)=B21*dwo(1) + B22*dwo(2)+ B23*dwo(3) + B24*dwo(4) + B25*dwo(5)
              dw(i,j,k,3)=B31*dwo(1) + B32*dwo(2)+ B33*dwo(3) + B34*dwo(4) + B35*dwo(5)
              dw(i,j,k,4)=B41*dwo(1) + B42*dwo(2)+ B43*dwo(3) + B44*dwo(4) + B45*dwo(5)
              dw(i,j,k,5)=B51*dwo(1) + B52*dwo(2)+ B53*dwo(3) + B54*dwo(4) + B55*dwo(5)

           enddo
        enddo
     enddo
  else
     do l=1,nwf
        do k=2,kl
           do j=2,jl
              do i=2,il
                 dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
                      * max(real(iblank(i,j,k), realType), zero)
              enddo
           enddo
        enddo
     enddo
  endif

end subroutine residual_block

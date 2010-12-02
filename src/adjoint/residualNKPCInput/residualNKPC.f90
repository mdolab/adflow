!
!      ******************************************************************
!      *                                                                *
!      * File:          residual.f90                                    *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-21-2008                                      *
!      * Last modified: 04-28-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine residualNKPC(wAdj,pAdj,siAdj,sjAdj,skAdj,volAdj,normAdj,&
     sFaceIAdj,sFaceJAdj,sFaceKAdj,&
     radIAdj,radJAdj,radKAdj,&
     dwAdj, iCell, jCell, kCell,  &  
     rotRateAdj, correctForK,nn,level,sps)
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
  use iteration
  implicit none

  !       Subroutine Variables
  integer(kind=intType) :: iCell, jCell, kCell,nn,level,sps
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(inout) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType), dimension(3),intent(in) ::rotRateAdj
  real(kind=realType), dimension(0:0,0:0,0:0,nTimeIntervalsSpectral), intent(inout) :: volAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,nTimeIntervalsSpectral), intent(inout) :: normAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(inout) :: wAdj
  real(kind=realType), dimension(-1:1,-1:1,-1:1,nTimeIntervalsSpectral) :: radIAdj,radJAdj,radKAdj
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral), intent(inout) :: dwAdj
  real(kind=realType),dimension(nw,nTimeIntervalsSpectral):: dwAdj2
  !integer(kind=intType), intent(in) :: discr
  logical, intent(in) :: correctForK

  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),intent(in) :: pAdj
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: fwAdj
  !
  !      Local variables.
  !
  integer(kind=intType) :: discr!sps, nn, discr
  integer(kind=intType) :: i, j, k, l

  logical :: fineGrid
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  !   Come back to this later....
!!$       ! Add the source terms from the level 0 cooling model.
!!$
!!$       call level0CoolingModel

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

  !!massFlowFamilyInv  = zero
  !!massFlowFamilyDiss = zero

  ! Set the value of the discretization, depending on the grid level,
  ! and the logical fineGrid, which indicates whether or not this
  ! is the finest grid level of the current mg cycle.

  discr = spaceDiscrCoarse
  if(currentLevel == 1) discr = spaceDiscr

  fineGrid = .false.
  if(currentLevel == groundLevel) fineGrid = .true.

  call inviscidCentralFluxNKPC(wAdj,  pAdj,  dwAdj,         &
       siAdj, sjAdj, skAdj, volAdj, &
       sFaceIAdj,sFaceJAdj,sFaceKAdj,&
       rotRateAdj,                  &
       iCell, jCell, kCell,nn,level,sps)


  ! Compute the artificial dissipation fluxes.
  ! This depends on the parameter discr.

  select case (discr)

  case (dissScalar) ! Standard scalar dissipation scheme.

     if( fineGrid ) then

        call inviscidDissFluxScalarNKPC(wAdj,  pAdj,  dwadj,&
             radIAdj,radJAdj,radKAdj, &
             iCell, jCell, kCell,nn,level,sps)


     else
        call terminate("residualAdj", &
             "ADjoint does not function on coarse grid level")
        !call inviscidDissFluxScalarCoarse
     endif

     !===========================================================

!!$!             case (dissMatrix) ! Matrix dissipation scheme.
!!$!
!!$!               if( fineGrid ) then
!!$!                 call inviscidDissFluxMatrixAdj()
!!$!               else
!!$!                 call terminate("residualAdj", &
!!$!                        "ADjoint does not function on coarse grid level")
!!$!                 !call inviscidDissFluxMatrixCoarse
!!$!               endif

     !===========================================================

!!$ !            case (dissCusp) ! Cusp dissipation scheme.
!!$!
!!$!               if( fineGrid ) then
!!$!                 call inviscidDissFluxCuspAdj()
!!$!               else
!!$!                 call terminate("residualAdj", &
!!$!                        "ADjoint does not function on coarse grid level")
!!$!                 !call inviscidDissFluxCuspCoarse
!!$!               endif

     !===========================================================

  case (upwind) ! Dissipation via an upwind scheme.


     call inviscidUpwindFluxNKPC(wAdj,  pAdj,  dwAdj, &
          siAdj, sjAdj, skAdj, &
          sFaceIAdj,sFaceJAdj,sFaceKAdj,&
          iCell, jCell, kCell,finegrid,nn,level,sps)

  end select
  !
  !          if( viscous ) call viscousFlux

  ! Add the dissipative and possibly viscous fluxes to the
  ! Euler fluxes. Loop over the owned cells and add fw to dw.
  ! Also multiply by iblank so that no updates occur in holes
  ! or on the overset boundary.

  ! ADD VISCOUS FLUXES HERE
  !            do l=1,nwf
  !               dwAdj(l,sps) =dwAdj(l,sps) + 

  !            enddo


end subroutine residualNKPC

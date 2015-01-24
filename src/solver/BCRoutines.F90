!      ******************************************************************
!      *                                                                *
!      * File:          BCRoutines.F90                                   *
!      * Author:        Gaetan K. W. Kenway                             *
!      * Starting date: 01-23-2015                                      *
!      * Last modified: 01-23-2015                                      *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * This module contains data structures *and* routines used       *
!      * for applying *all* boundary conditions for Navier Stokes part  *
!      * of the code. The reason for using a module to contain the      *
!      * routines is that due to the use of pointers, it eliminates the *
!      * need for using interfaces. All former bc*.f90 routines are now *
!      * included in this module.                                       *
!      *                                                                *
!      ******************************************************************
!
module BCRoutines

  use constants
  implicit none
  save
#ifndef TAPENADE_REVERSE
  real(kind=realType), dimension(:,:,:), pointer :: ww0, ww1, ww2, ww3
  real(kind=realType), dimension(:,:),   pointer :: pp0, pp1, pp2, pp3
  real(kind=realType), dimension(:,:),   pointer :: rlv0, rlv1, rlv2, rlv3
  real(kind=realType), dimension(:,:),   pointer :: rev0, rev1, rev2, rev3
  real(kind=realType), dimension(:,:),   pointer :: gamma0, gamma1, gamma2, gamma3
  real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
  real(kind=realType), dimension(:,:,:), pointer :: ss
  real(kind=realType), dimension(:,:,:), pointer :: xline
#else
  real(kind=realType), dimension(:,:,:), allocatable :: ww0, ww1, ww2, ww3
  real(kind=realType), dimension(:,:)  , allocatable :: pp0, pp1, pp2, pp3
  real(kind=realType), dimension(:,:)  , allocatable :: rlv0, rlv1, rlv2, rlv3
  real(kind=realType), dimension(:,:)  , allocatable :: rev0, rev1, rev2, rev3
  real(kind=realType), dimension(:,:  ), allocatable :: gamma0, gamma1, gamma2, gamma3
#endif
  integer(kind=intType) :: iSize, jSize
contains

  subroutine applyAllBC_block2(secondHalo)

    ! Apply BC's for a single block

    use blockPointers
    use flowVarRefState
    use inputDiscretization
    use inputTimeSpectral
    use iteration
    use bcTypes
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo

    ! Local variables.
    logical :: correctForK
    integer(kind=intType) :: nn
    !
    ! Determine whether or not the total energy must be corrected
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

    ! Apply all the boundary conditions. The order is important!  Only
    ! some of them have been AD'ed

    ! ------------------------------------
    !  Symmetry Boundary Condition 
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nBocos
       if (bcType(nn) == symm) then 
          call setBCPointers2(nn)
          call bcSymm2(nn, secondHalo)
          call resetBCPointers2(nn)
       end if
    end do

#ifndef USE_TAPENADE
    ! ------------------------------------
    !  Symmetry Polar Boundary Condition 
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nBocos
       if (bcType(nn) == symmPolar) then
          call setBCPointers2(nn)
          call bcSymmPolar(nn, secondHalo)
          call resetBCPointers2(nn)
       end if
    end do
#endif

    ! ------------------------------------
    !  Adibatic Wall Boundary Condition 
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nViscBocos
       if (bcType(nn) == NSWallAdiabatic) then 
          call setBCPointers2(nn)
          call bcNSWallAdiabatic2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Isotermal Wall Boundary Condition 
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nViscBocos
       if (bcType(nn) == NSWallIsoThermal) then 
          call setBCPointers2(nn)
          call bcNSWallIsothermal2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Farfield Boundary Condition 
    ! ------------------------------------
    if (precond == Turkel .or. precond == ChoiMerkle) then 
       call terminate("applyAllBC", &
            "Farfield Turkel and Coid/Merkle preconditioners not implemented")
    end if
    !$AD II-LOOP
    do nn=1,nBocos
       if (bcType(nn) == farField) then
          call setBCPointers2(nn)
          call bcFarField2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

#ifndef USE_TAPENADE
    ! ------------------------------------
    !  Subsonic Outflow Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == subSonicOutFlow .or. &
            bcType(nn) == MassBleedOutflow) then 
          call setBCPointers2(nn)
          call bcSubSonicOutFlow2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Subsonic Inflow Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == subSonicInFlow) then 
          call setBCPointers2(nn)
          call bcSubSonicInflow2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Bleed Inflow Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == MassBleedInflow) then
          call setBCPointers2(nn)
          call bcBleedInflow2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Mdot Engine Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == mdot) then 
          call setBCPointers2(nn)
          call bcMDot2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  bcThrust Engine Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == thrust) then 
          call setBCPointers2(nn)
          call bcThrust2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Extrapolation Boundary Condition 
    ! ------------------------------------
    ! Extrapolation boundary conditions; this also includes
    ! the supersonic outflow boundary conditions. The difference
    ! between the two is that the extrap boundary conditions
    ! correspond to singular lines and supersonic outflow
    ! boundaries to physical boundaries. The treatment however
    ! is identical.
    do nn=1,nBocos
       if (bcType(nn) == extrap .or. &
            bcType(nn) == SupersonicOutFlow) then 
          call setBCPointers2(nn)
          call bcExtrap2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do
#endif

    ! ------------------------------------
    !  Euler Wall Boundary Condition 
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1,nBocos
       if (bcType(nn) == EulerWall) then
          call setBCPointers2(nn)
          call bcEulerWall2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

#ifndef USE_TAPENADE
    ! ------------------------------------
    !  Domain Interface Boundary Condition
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == DomainInterfaceAll .or. & 
            bcTYpe(nn) == DomainInterfaceRhoUVW .or. &
            bcType(nn) == DomainInterfaceP .or. &
            bcType(nn) == DomainInterfaceRho .or. &
            bcType(nn) == DomainInterfaceTotal) then
          call setBCPointers2(nn)
          call bcDomainInterface2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do

    ! ------------------------------------
    !  Supersonic inflow condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == SupersonicInflow) then 
          call setBCPointers2(nn)
          call bcSupersonicInflow2(nn, secondHalo, correctForK)
          call resetBCPointers2(nn)
       end if
    end do
#endif
  end subroutine applyAllBC_block2

  ! ===================================================================
  !   Actual implementation of each of the boundary condition routines
  ! ===================================================================
  subroutine bcSymm2(nn, secondHalo)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcSymm applies the symmetry boundary conditions to a block.    *
    ! * It is assumed that the pointers in blockPointers are already   *
    ! * set to the correct block on the correct grid level.            *
    ! *                                                                *
    ! * In case also the second halo must be set the loop over the     *
    ! * boundary subfaces is executed twice. This is the only correct  *
    ! * way in case the block contains only 1 cell between two         *
    ! * symmetry planes, i.e. a 2D problem.                            *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType) :: vn, nnx, nny, nnz

    ! Loop over the generic subface to set the state in the
    ! 1-st level halos 

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1
       
       ! Determine twice the normal velocity component,
       ! which must be substracted from the donor velocity
       ! to obtain the halo velocity.
       
       vn = two*(ww2(i,j,ivx)*BCData(nn)%norm(i,j,1) + &
            ww2(i,j,ivy)*BCData(nn)%norm(i,j,2) + &
            ww2(i,j,ivz)*BCData(nn)%norm(i,j,3))
       
       ! Determine the flow variables in the halo cell.
       
       ww1(i,j,irho) = ww2(i,j,irho)
       ww1(i,j,ivx) = ww2(i,j,ivx) - vn*BCData(nn)%norm(i,j,1)
       ww1(i,j,ivy) = ww2(i,j,ivy) - vn*BCData(nn)%norm(i,j,2)
       ww1(i,j,ivz) = ww2(i,j,ivz) - vn*BCData(nn)%norm(i,j,3)
       ww1(i,j,irhoE) = ww2(i,j,irhoE)
       
       ! Simply copy the turbulent variables.
       
       !$AD II-LOOP
       do l=nt1MG,nt2MG
          ww1(i,j,l) = ww2(i,j,l)
       enddo
       
       ! Set the pressure and gamma and possibly the
       ! laminar and eddy viscosity in the halo.
       
       gamma1(i,j) = gamma2(i,j)
       pp1(i,j)    = pp2(i,j)
       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)
    enddo

    if (secondHalo) then

       ! If we need the second halo, do everything again, but using ww0,
       ! ww3 etc instead of ww2 and ww1. 

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          vn = two*(ww3(i,j,ivx)*BCData(nn)%norm(i,j,1) + &
               ww3(i,j,ivy)*BCData(nn)%norm(i,j,2) + &
               ww3(i,j,ivz)*BCData(nn)%norm(i,j,3))

          ! Determine the flow variables in the halo cell.
          ww0(i,j,irho) = ww3(i,j,irho)
          ww0(i,j,ivx) = ww3(i,j,ivx) - vn*BCData(nn)%norm(i,j,1)
          ww0(i,j,ivy) = ww3(i,j,ivy) - vn*BCData(nn)%norm(i,j,2)
          ww0(i,j,ivz) = ww3(i,j,ivz) - vn*BCData(nn)%norm(i,j,3)

          ww0(i,j,irhoE) = ww3(i,j,irhoE)

          !$AD II-LOOP
          do l=nt1MG,nt2MG
             ww0(i,j,l) = ww3(i,j,l)
          enddo

          ! Set the pressure and gamma and possibly the
          ! laminar and eddy viscosity in the halo.

          gamma0(i,j) = gamma3(i,j)
          pp0(i,j)    = pp3(i,j)
          if( viscous )   rlv0(i,j) = rlv3(i,j)
          if( eddyModel ) rev0(i,j) = rev3(i,j)
       enddo
    end if
  end subroutine bcSymm2

#ifndef USE_TAPENADE
  subroutine bcSymmPolar2(nn, secondHalo)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcSymmPolar applies the polar symmetry boundary conditions     *
    ! * to a singular line of a block. It is assumed that the pointers *
    ! * in blockPointers are already set to the correct block on the   *
    ! * correct grid level.                                            *
    ! * The polar symmetry condition is a special case of a degenerate *
    ! * line, as this line is the axi-symmetric centerline.            *
    ! *                                                                *
    ! ******************************************************************
    !

    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii, mm
    real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

    ! Loop over the generic subface to set the state in the
    ! 1-st level halos 
    if (.not. secondHalo) then 
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Determine the unit vector along the degenerated face.
          ! However it is not known which is the singular
          ! direction and therefore determine the direction along
          ! the diagonal (i,j) -- (i-1,j-1), which is correct for
          ! both singular i and j-direction. Note that due to the
          ! usage of the pointer xline there is an offset of +1
          ! in the indices and therefore (i+1,j+1) - (i,j) must
          ! be used to determine this vector.

          nnx = xline(i+1,j+1,1) - xline(i,j,1)
          nny = xline(i+1,j+1,2) - xline(i,j,2)
          nnz = xline(i+1,j+1,3) - xline(i,j,3)

          ! Determine the unit vector in this direction.

          tmp = one/sqrt(nnx*nnx + nny*nny + nnz*nnz)
          nnx = nnx*tmp
          nny = nny*tmp
          nnz = nnz*tmp

          ! Determine twice the tangential velocity vector of the
          ! internal cell.

          tmp = two*(ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny &
               +      ww2(i,j,ivz)*nnz)
          vtx = tmp*nnx
          vty = tmp*nny
          vtz = tmp*nnz

          ! Determine the flow variables in the halo cell. The
          ! velocity is constructed such that the average of the
          ! internal and the halo cell is along the centerline.
          ! Note that the magnitude of the velocity does not
          ! change and thus the energy is identical.

          ww1(i,j,irho)  = ww2(i,j,irho)
          ww1(i,j,ivx)   = vtx - ww2(i,j,ivx)
          ww1(i,j,ivy)   = vty - ww2(i,j,ivy)
          ww1(i,j,ivz)   = vtz - ww2(i,j,ivz)
          ww1(i,j,irhoE) = ww2(i,j,irhoE)

          ! Simply copy the turbulent variables.

          do l=nt1MG,nt2MG
             ww1(i,j,l) = ww2(i,j,l)
          enddo

          ! Set the pressure and possibly the laminar and
          ! eddy viscosity in the halo.

          pp1(i,j) = pp2(i,j)
          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)
       end do
    else
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Determine the unit vector along the degenerated face.
          ! However it is not known which is the singular
          ! direction and therefore determine the direction along
          ! the diagonal (i,j) -- (i-1,j-1), which is correct for
          ! both singular i and j-direction. Note that due to the
          ! usage of the pointer xline there is an offset of +1
          ! in the indices and therefore (i+1,j+1) - (i,j) must
          ! be used to determine this vector.

          nnx = xline(i+1,j+1,1) - xline(i,j,1)
          nny = xline(i+1,j+1,2) - xline(i,j,2)
          nnz = xline(i+1,j+1,3) - xline(i,j,3)

          ! Determine the unit vector in this direction.

          tmp = one/sqrt(nnx*nnx + nny*nny + nnz*nnz)
          nnx = nnx*tmp
          nny = nny*tmp
          nnz = nnz*tmp

          ! Determine twice the tangential velocity vector of the
          ! internal cell.

          tmp = two*(ww3(i,j,ivx)*nnx + ww3(i,j,ivy)*nny &
               +      ww3(i,j,ivz)*nnz)
          vtx = tmp*nnx
          vty = tmp*nny
          vtz = tmp*nnz

          ! Determine the flow variables in the halo cell. The
          ! velocity is constructed such that the average of the
          ! internal and the halo cell is along the centerline.
          ! Note that the magnitude of the velocity does not
          ! change and thus the energy is identical.

          ww0(i,j,irho)  = ww3(i,j,irho)
          ww0(i,j,ivx)   = vtx - ww3(i,j,ivx)
          ww0(i,j,ivy)   = vty - ww3(i,j,ivy)
          ww0(i,j,ivz)   = vtz - ww3(i,j,ivz)
          ww0(i,j,irhoE) = ww3(i,j,irhoE)

          ! Simply copy the turbulent variables.

          do l=nt1MG,nt2MG
             ww0(i,j,l) = ww3(i,j,l)
          enddo

          ! Set the pressure and possibly the laminar and
          ! eddy viscosity in the halo.

          pp0(i,j) = pp3(i,j)
          if( viscous )   rlv0(i,j) = rlv3(i,j)
          if( eddyModel ) rev0(i,j) = rev3(i,j)
       enddo
    end if
  end subroutine bcSymmPolar2
#endif
  subroutine bcNSWallAdiabatic2(nn, secondHalo, correctForK)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * bcNSWallAdiabatic applies the viscous adiabatic wall           *
    !      * boundary condition the pointers already defined.               *
    !      *                                                                *
    !      ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use iteration
    implicit none

    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn
    integer(kind=intType) :: i, j, ii
    real(kind=realType) :: rhok

    ! Apply the BCWall In case the turbulent transport equations are solved
    ! together with the mean flow equations, aplly the viscous
    ! wall boundary conditions for the turbulent variables.
    ! No need to extrapolate the secondary halo's, because this
    ! is done in extrapolate2ndHalo.

    if( turbCoupled ) call turbBCNSWall(.false.)

    ! Initialize rhok to zero. This will be overwritten if a
    ! correction for k must be applied.

    rhok = zero

    ! Loop over the generic subface to set the state in the
    ! halo cells.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ! Set the value of rhok if a correcton must be applied.
       ! It probably does not matter too much, because k is very
       ! small near the wall.

       if( correctForK ) rhok = ww2(i,j,irho)*ww2(i,j,itu1)

       ! Determine the variables in the halo. As the spacing
       ! is very small a constant pressure boundary condition
       ! (except for the k correction) is okay. Take the slip
       ! velocity into account.

       ww1(i,j,irho) =  ww2(i,j,irho)
       ww1(i,j,ivx)  = -ww2(i,j,ivx) + two*bcData(nn)%uSlip(i,j,1)
       ww1(i,j,ivy)  = -ww2(i,j,ivy) + two*bcData(nn)%uSlip(i,j,2)
       ww1(i,j,ivz)  = -ww2(i,j,ivz) + two*bcData(nn)%uSlip(i,j,3)
       pp1(i,j)      =  pp2(i,j) - four*third*rhok

       ! Set the viscosities. There is no need to test for a
       ! viscous problem of course. The eddy viscosity is
       ! set to the negative value, as it should be zero on
       ! the wall.

       rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = -rev2(i,j)
    enddo

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcNSWallAdiabatic2

  subroutine bcNSWallIsoThermal2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcNSWallAdiabatic applies the viscous isothermal wall          *
    ! * boundary condition to a block. It is assumed that the          *
    ! * BCPointers are already set                                     *
    ! *                                                                *
    ! ******************************************************************

    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, ii

    real(kind=realType) :: rhok, t2, t1

    ! In case the turbulent transport equations are solved
    ! together with the mean flow equations, aplly the viscous
    ! wall boundary conditions for the turbulent variables.
    ! No need to extrapolate the secondary halo's, because this
    ! is done in extrapolate2ndHalo.

    if( turbCoupled ) call turbBCNSWall(.false.)

    ! Initialize rhok to zero. This will be overwritten if a
    ! correction for k must be applied.

    rhok = zero

    ! Loop over the generic subface to set the state in the
    ! halo cells.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ! Set the value of rhok if a correcton must be applied.
       ! It probably does not matter too much, because k is very
       ! small near the wall.

       if( correctForK ) rhok = ww2(i,j,irho)*ww2(i,j,itu1)

       ! Compute the temperature in the internal cell and in the
       ! halo cell such that the average is the wall temperature.

       t2 = pp2(i,j)/(RGas*ww2(i,j,irho))
       t1 = two*bcData(nn)%TNS_Wall(i,j) - t2

       ! Make sure that t1 is within reasonable bounds. These
       ! bounds are such that the clipping is never active in the
       ! converged solution; it is only to avoid instabilities
       ! during the convergence.

       t1 = max(half*bcData(nn)%TNS_Wall(i,j), t1)
       t1 = min(two *bcData(nn)%TNS_Wall(i,j), t1)

       ! Determine the variables in the halo. As the spacing
       ! is very small a constant pressure boundary condition
       ! (except for the k correction) is okay. Take the slip
       ! velocity into account.

       pp1(i,j)      =  pp2(i,j) - four*third*rhok
       ww1(i,j,irho) =  pp1(i,j)/(RGas*t1)
       ww1(i,j,ivx)  = -ww2(i,j,ivx) + two*bcData(nn)%uSlip(i,j,1)
       ww1(i,j,ivy)  = -ww2(i,j,ivy) + two*bcData(nn)%uSlip(i,j,2)
       ww1(i,j,ivz)  = -ww2(i,j,ivz) + two*bcData(nn)%uSlip(i,j,3)

       ! Set the viscosities. There is no need to test for a
       ! viscous problem of course. The eddy viscosity is
       ! set to the negative value, as it should be zero on
       ! the wall.

       rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = -rev2(i,j)
    enddo

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcNSWallIsoThermal2

#ifndef USE_TAPENADE
  subroutine bcSubsonicOutflow2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcSubsonicOutflow applies the subsonic outflow boundary        *
    ! * condition, static pressure prescribed, to a block. It is       *
    ! * assumed that the pointers in blockPointers are already set to  *
    ! * the correct block on the correct grid level.                   *
    ! * Exactly the same boundary condition is also applied for an     *
    ! * outflow mass bleed. Therefore the test is for both a subsonic  *
    ! * outflow and an bleed outflow.                                  *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType), parameter :: twothird = two*third
    real(kind=realType) :: ovg, ovgm1, nnx, nny, nnz
    real(kind=realType) :: pExit, pInt, r, a2, a, ac, ss
    real(kind=realType) :: ue, ve, we, qne, qnh

    ! Loop over the generic subface to set the state in the
    ! halo cells.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1


       ! Store a couple of variables, such as the static
       ! pressure and grid unit outward normal, a bit easier.

       pExit = BCData(nn)%ps(i,j)

       nnx = BCData(nn)%norm(i,j,1)
       nny = BCData(nn)%norm(i,j,2)
       nnz = BCData(nn)%norm(i,j,3)

       ! Abbreviate 1/gamma and 1/(gamma -1) a bit easier.

       ovg   = one/gamma2(i,j)
       ovgm1 = one/(gamma2(i,j)-one)

       ! Store the internal pressure and correct for the
       ! possible presence of a k-equation.

       pInt = pp2(i,j)
       if( correctForK ) &
            pInt = pInt - twothird*ww2(i,j,irho)*ww2(i,j,itu1)

       ! Compute the velocity components, the normal velocity
       ! and the speed of sound for the internal cell.

       r   = one/ww2(i,j,irho)
       a2  = gamma2(i,j)*pInt*r
       a   = sqrt(a2)
       ue  = ww2(i,j,ivx)
       ve  = ww2(i,j,ivy)
       we  = ww2(i,j,ivz)
       qne = ue*nnx + ve*nny + we*nnz

       ! Compute the entropy and the acoustic variable.
       ! These riemann inVariants, as well as the tangential
       ! velocity components, are extrapolated.

       ss = pInt*(r**gamma2(i,j))
       ac = qne + two*a*ovgm1

       ! Compute the state in the halo.

       ww1(i,j,irho) = (pExit/ss)**ovg
       pp1(i,j)      = pExit
       a             = sqrt(gamma2(i,j)*pExit/ww1(i,j,irho))
       qnh           = ac - two*a*ovgm1
       ww1(i,j,ivx)  = ue + (qnh - qne)*nnx
       ww1(i,j,ivy)  = ve + (qnh - qne)*nny
       ww1(i,j,ivz)  = we + (qnh - qne)*nnz

       ! Extrapolate the primitive turbulent variables.

       do l=nt1MG,nt2MG
          ww1(i,j,l) = ww2(i,j,l)
       enddo

       ! Correct the pressure if a k-equation is present.

       if( correctForK )   &
            pp1(i,j) = pp1(i,j) &
            + twothird*ww1(i,j,irho)*ww1(i,j,itu1)

       ! Set the viscosities in the halo to the viscosities
       ! in the donor cell.

       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)

    enddo

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcSubsonicOutflow2

  subroutine bcSubsonicInflow2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcSubsonicInflow applies the subsonic outflow boundary         *
    ! * condition, total pressure, total density and flow direction    *
    ! * prescribed,  to a block. It is assumed that the pointers in    *
    ! * blockPointers are already set to the correct block on the      *
    ! * correct grid level.                                            *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use inputDiscretization
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii

    real(kind=realType) :: gm1, ovgm1
    real(kind=realType) :: ptot, ttot, htot, a2tot, r, alpha, beta
    real(kind=realType) :: aa2, bb, cc, dd, q, q2, a2, m2, scaleFact
    real(kind=realType) :: ssx, ssy, ssz, nnx, nny, nnz
    real(kind=realType) :: rho, velx, vely, velz

    ! Determine the boundary treatment to be used.

    select case (BCData(nn)%subsonicInletTreatment)

    case (totalConditions)

       ! The total conditions have been prescribed.

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Store a couple of variables, such as the total
          ! pressure, total temperature, total enthalpy, flow
          ! direction and grid unit outward normal, a bit easier.

          ptot = BCData(nn)%ptInlet(i,j)
          ttot = BCData(nn)%ttInlet(i,j)
          htot = BCData(nn)%htInlet(i,j)

          ssx  = BCData(nn)%flowXdirInlet(i,j)
          ssy  = BCData(nn)%flowYdirInlet(i,j)
          ssz  = BCData(nn)%flowZdirInlet(i,j)

          nnx = BCData(nn)%norm(i,j,1)
          nny = BCData(nn)%norm(i,j,2)
          nnz = BCData(nn)%norm(i,j,3)

          ! Some abbreviations in which gamma occurs.

          gm1   = gamma2(i,j) - one
          ovgm1 = one/gm1

          ! Determine the acoustic Riemann variable that must be
          ! extrapolated from the domain.

          r    = one/ww2(i,j,irho)
          a2   = gamma2(i,j)*pp2(i,j)*r
          beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
               + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

          ! Correct the value of the Riemann invariant if total
          ! enthalpy scaling must be applied. This scaling may
          ! be needed for stability if large gradients of the
          ! total temperature are prescribed.

          scaleFact = one
          if( hScalingInlet ) &
               scaleFact = sqrt(htot/(r*(ww2(i,j,irhoE) + pp2(i,j))))

          beta = beta*scaleFact

          ! Compute the value of a2 + 0.5*gm1*q2, which is the
          ! total speed of sound for constant cp. However, the
          ! expression below is also valid for variable cp,
          ! although a linearization around the value of the
          ! internal cell is performed.

          q2    = ww2(i,j,ivx)**2 + ww2(i,j,ivy)**2 &
               + ww2(i,j,ivz)**2
          a2tot = gm1*(htot - r*(ww2(i,j,irhoE) + pp2(i,j)) &
               +      half*q2) + a2

          ! Compute the dot product between the normal and the
          ! velocity direction. This value should be negative.

          alpha = nnx*ssx + nny*ssy + nnz*ssz

          ! Compute the coefficients in the quadratic equation
          ! for the magnitude of the velocity.

          aa2 =  half*gm1*alpha*alpha + one
          bb = -gm1*alpha*beta
          cc =  half*gm1*beta*beta - two*ovgm1*a2tot

          ! Solve the equation for the magnitude of the
          ! velocity. As this value must be positive and both aa2
          ! and bb are positive (alpha is negative and beta is
          ! positive up till Mach = 5.0 or so, which is not
          ! really subsonic anymore), it is clear which of the
          ! two possible solutions must be taken. Some clipping
          ! is present, but this is normally not active.

          dd = bb*bb - four*aa2*cc
          dd = sqrt(max(zero,dd))
          q  = (-bb + dd)/(two*aa2)
          q  = max(zero,q)
          q2 = q*q

          ! Compute the speed of sound squared from the total
          ! speed of sound equation (== total enthalpy equation
          ! for constant cp).

          a2 = a2tot - half*gm1*q2

          ! Compute the Mach number squared and cut it between
          ! 0.0 and 1.0. Adapt the velocity and speed of sound
          ! squared accordingly.

          m2 = q2/a2
          m2 = min(one,m2)
          q2 = m2*a2
          q  = sqrt(q2)
          a2 = a2tot - half*gm1*q2

          ! Compute the velocities in the halo cell and use rho,
          ! rhoe and p as temporary buffers to store the total
          ! temperature, total pressure and static temperature.

          ww1(i,j,ivx)  = q*ssx
          ww1(i,j,ivy)  = q*ssy
          ww1(i,j,ivz)  = q*ssz

          ww1(i,j,irho)  = ttot
          pp1(i,j)       = ptot
          ww1(i,j,irhoE) = a2/(gamma2(i,j)*RGas)

          ! Compute the turbulent variables, which are
          ! prescribed.

          do l=nt1MG,nt2MG
             ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the pressure and density for these halo's.

       call pRhoSubsonicInlet(ww1, pp1, correctForK)

       !===========================================================

    case (massFlow)

       ! Density and velocity vector prescribed.

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Store a couple of variables, such as the density,
          ! velocity and grid unit outward normal, a bit easier.

          rho  = BCData(nn)%rho(i,j)
          velx = BCData(nn)%velx(i,j)
          vely = BCData(nn)%vely(i,j)
          velz = BCData(nn)%velz(i,j)

          nnx = BCData(nn)%norm(i,j,1)
          nny = BCData(nn)%norm(i,j,2)
          nnz = BCData(nn)%norm(i,j,3)

          ! Some abbreviations in which gamma occurs.

          gm1   = gamma2(i,j) - one
          ovgm1 = one/gm1

          ! Determine the acoustic Riemann variable that must be
          ! extrapolated from the domain.

          r    = one/ww2(i,j,irho)
          a2   = gamma2(i,j)*pp2(i,j)*r
          beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
               + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

          ! Compute the speed of sound squared in the halo.

          a2 = half*gm1*(beta - velx*nnx - vely*nny - velz*nnz)
          a2 = max(zero,a2)
          a2 = a2*a2

          ! Compute the pressure in the halo, assuming a
          ! constant value of gamma.

          pp1(i,j) = rho*a2/gamma2(i,j)

          ! Simply copy the density and velocities.

          ww1(i,j,irho) = rho
          ww1(i,j,ivx)  = velx
          ww1(i,j,ivy)  = vely
          ww1(i,j,ivz)  = velz

          ! Compute the turbulent variables, which are
          ! prescribed.

          do l=nt1MG,nt2MG
             ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

    end select

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcSubsonicInflow2

  subroutine bcBleedInflow2(nn, secondHalo, correctForK)
    !
    ! *****************************************************************
    ! *                                                                *
    ! * bcBleedInflow applies the boundary conditions to an inflow     *
    ! * bleed region of a block.                                       *
    ! * It is assumed that the pointers in blockPointers are already   *
    ! * set to the correct block on the correct grid level.            *
    ! *                                                                *
    ! ******************************************************************
    !
    use BCTypes
    use bleedFlows
    use blockPointers
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    call terminate("bcBleedInflow", "Not implemented yet")

  end subroutine bcBleedInflow2

  subroutine bcMDot2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcMdot applies the duct inflow boundary condition to a block.  *
    ! * It is assumed that the pointers in blockPointers are already   *
    ! * set to the correct block on the correct grid level.            *
    ! *                                                                *
    ! ******************************************************************

    use BCTypes
    use bleedFlows
    use blockPointers
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    call terminate("bcMdot", "Not implemented yet")

  end subroutine bcMDot2

  subroutine bcThrust2(nn, secondHalo, correctForK)
    ! ******************************************************************
    ! *                                                                *
    ! * bcThrust applies the duct outflow boundary condition to a      *
    ! * block. It is assumed that the pointers in blockPointers are    *
    ! * already set to the correct block on the correct grid level.    *
    ! *                                                                *
    ! ******************************************************************

    use BCTypes
    use bleedFlows
    use blockPointers
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    call terminate("bcThrust", "Not implemented yet")

  end subroutine bcThrust2
#endif

  subroutine bcEulerWall2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcEulerWall applies the inviscid wall boundary condition to    *
    ! * a block. It is assumed that the bcpointers are                 *
    ! * already set to the correct block on the correct grid level.    *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use inputDiscretization
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: j, k, l, ii
    integer(kind=intType) :: jm1, jp1, km1, kp1
    integer(kind=intType) :: wallTreatment

    real(kind=realType) :: sixa, siya, siza, sjxa, sjya, sjza
    real(kind=realType) :: skxa, skya, skza, a1, b1
    real(kind=realType) :: rxj, ryj, rzj, rxk, ryk, rzk
    real(kind=realType) :: dpj, dpk, ri, rj, rk, qj, qk, vn
    real(kind=realType) :: uux, uuy, uuz
    real(kind=realType), dimension(iSize, jSize) :: grad

    ! Make sure that on the coarser grids the constant pressure
    ! boundary condition is used.

    wallTreatment = wallBcTreatment
    if(currentLevel > groundLevel) wallTreatment = constantPressure

    ! **************************************************************
    ! *                                                            *
    ! * Determine the boundary condition treatment and compute the *
    ! * undivided pressure gradient accordingly. This gradient is  *
    ! * temporarily stored in the halo pressure.                   *
    ! *                                                            *
    ! **************************************************************
    !
    BCTreatment: select case (wallTreatment)

    case (constantPressure)

       ! Constant pressure. Set the gradient to zero.
       grad = zero

    case (linExtrapolPressure)

       ! Linear extrapolation. 
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          j = mod(ii, isize) + 1
          k = ii/isize + 1
          grad(j,k) = pp3(j,k) - pp2(j,k)
       end do
#ifndef TAPENADE_REVERSE
    case (normalMomentum)

       ! Pressure gradient is computed using the normal momentum
       ! equation. First set a couple of additional variables for
       ! the normals, depending on the block face. Note that the
       ! construction 1: should not be used in these pointers,
       ! because element 0 is needed. Consequently there will be
       ! an offset of 1 for these normals. This is commented in
       ! the code. For moving faces also the grid velocity of
       ! the 1st cell center from the wall is needed.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          j = mod(ii, isize) + 1
          k = ii/isize + 1

          ! Store the indices k+1, k-1 a bit easier and make
          ! sure that they do not exceed the range of the arrays.

          km1 = k-1; km1 = max(1, km1)
          kp1 = k+1; kp1 = min(jsize, kp1)

          ! Compute the scaling factor for the central difference
          ! in the k-direction.

          b1 = one/max(1_intType,(kp1-km1))

          ! The indices j+1 and j-1. Make sure that they
          ! do not exceed the range of the arrays.

          jm1 = j-1; jm1 = max(1, jm1)
          jp1 = j+1; jp1 = min(isize, jp1)

          ! Compute the scaling factor for the central
          ! difference in the j-direction.

          a1 = one/max(1_intType,(jp1-jm1))

          ! Compute (twice) the average normal in the generic i,
          ! j and k-direction. Note that in j and k-direction
          ! the average in the original indices should be taken
          ! using j-1 and j (and k-1 and k). However due to the
          ! usage of pointers ssj and ssk there is an offset in
          ! the indices of 1 and therefore now the correct
          ! average is obtained with the indices j and j+1
          ! (k and k+1).

          sixa = two*ssi(j,k,1)
          siya = two*ssi(j,k,2)
          siza = two*ssi(j,k,3)

          sjxa = ssj(j,k,1) + ssj(j+1,k,1)
          sjya = ssj(j,k,2) + ssj(j+1,k,2)
          sjza = ssj(j,k,3) + ssj(j+1,k,3)

          skxa = ssk(j,k,1) + ssk(j,k+1,1)
          skya = ssk(j,k,2) + ssk(j,k+1,2)
          skza = ssk(j,k,3) + ssk(j,k+1,3)

          ! Compute the difference of the normal vector and
          ! pressure in j and k-direction. As the indices are
          ! restricted to the 1st halo-layer, the computation
          ! of the internal halo values is not consistent;
          ! however this is not really a problem, because these
          ! values are overwritten in the communication pattern.

          rxj = a1*(BCData(nn)%norm(jp1,k,1) - BCData(nn)%norm(jm1,k,1))
          ryj = a1*(BCData(nn)%norm(jp1,k,2) - BCData(nn)%norm(jm1,k,2))
          rzj = a1*(BCData(nn)%norm(jp1,k,3) - BCData(nn)%norm(jm1,k,3))
          dpj = a1*(pp2(jp1,k)    - pp2(jm1,k))

          rxk = b1*(BCData(nn)%norm(j,kp1,1) - BCData(nn)%norm(j,km1,1))
          ryk = b1*(BCData(nn)%norm(j,kp1,2) - BCData(nn)%norm(j,km1,2))
          rzk = b1*(BCData(nn)%norm(j,kp1,3) - BCData(nn)%norm(j,km1,3))
          dpk = b1*(pp2(j,kp1)    - pp2(j,km1))

          ! Compute the dot product between the unit vector
          ! and the normal vectors in i, j and k-direction.

          ri = BCData(nn)%norm(j,k,1)*sixa + BCData(nn)%norm(j,k,2)*siya &
               + BCData(nn)%norm(j,k,3)*siza
          rj = BCData(nn)%norm(j,k,1)*sjxa + BCData(nn)%norm(j,k,2)*sjya &
               + BCData(nn)%norm(j,k,3)*sjza
          rk = BCData(nn)%norm(j,k,1)*skxa + BCData(nn)%norm(j,k,2)*skya &
               + BCData(nn)%norm(j,k,3)*skza

          ! Store the velocity components in uux, uuy and uuz and
          ! subtract the mesh velocity if the face is moving.

          uux = ww2(j,k,ivx)
          uuy = ww2(j,k,ivy)
          uuz = ww2(j,k,ivz)

          if( addGridVelocities ) then
             uux = uux - ss(j,k,1)
             uuy = uuy - ss(j,k,2)
             uuz = uuz - ss(j,k,3)
          endif

          ! Compute the velocity components in j and
          ! k-direction.

          qj = uux*sjxa + uuy*sjya + uuz*sjza
          qk = uux*skxa + uuy*skya + uuz*skza

          ! Compute the pressure gradient, which is stored
          ! in pp1. I'm not entirely sure whether this
          ! formulation is correct for moving meshes. It could
          ! be that an additional term is needed there.

          grad(j,k) = ((qj*(uux*rxj + uuy*ryj + uuz*rzj)      &
               +   qk*(uux*rxk + uuy*ryk + uuz*rzk))     &
               *  ww2(j,k,irho) - rj*dpj - rk*dpk)/ri
       enddo
#endif
    end select BCTreatment

    ! Determine the state in the halo cell. Again loop over
    ! the cell range for this subface.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       j = mod(ii, isize) + 1
       k = ii/isize + 1

       ! Compute the pressure density and velocity in the
       ! halo cell. Note that rface is the grid velocity
       ! component in the direction of norm, i.e. outward
       ! pointing.

       pp1(j,k) = dim(pp2(j,k), grad(j,k))

       vn = two*(BCData(nn)%rface(j,k) - &
            ww2(j,k,ivx)*BCData(nn)%norm(j,k,1) - &
            ww2(j,k,ivy)*BCData(nn)%norm(j,k,2) - &
            ww2(j,k,ivz)*BCData(nn)%norm(j,k,3))

       ww1(j,k,irho) = ww2(j,k,irho)
       ww1(j,k,ivx)  = ww2(j,k,ivx) + vn*BCData(nn)%norm(j,k,1)
       ww1(j,k,ivy)  = ww2(j,k,ivy) + vn*BCData(nn)%norm(j,k,2)
       ww1(j,k,ivz)  = ww2(j,k,ivz) + vn*BCData(nn)%norm(j,k,3)

       ! Just copy the turbulent variables.

       do l=nt1MG,nt2MG
          ww1(j,k,l) = ww2(j,k,l)
       enddo

       ! The laminar and eddy viscosity, if present.

       if( viscous )    rlv1(j,k) = rlv2(j,k)
       if( eddyModel ) rev1(j,k) = rev2(j,k)

    enddo

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcEulerWall2

#ifndef USE_TAPENADE
  subroutine bcDomainInterface2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcDomainInterface applies the domain-interface boundary        *
    ! * condition, where necessary flow variables are obtained from    *
    ! * the coupler. More options can be added in the future.          *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use couplerParam
    use flowVarRefState
    use inputDiscretization
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType), parameter :: twoThird = two*third

    real(kind=realType) :: gm1, ovg, ovgm1
    real(kind=realType) :: ptot, ttot, htot, a2tot, r, alpha, beta
    real(kind=realType) :: pExit, pInt, a, ac, ss, scaleFact
    real(kind=realType) :: aa2, bb, cc, dd, q, q2, a2, m2
    real(kind=realType) :: ssx, ssy, ssz, nnx, nny, nnz
    real(kind=realType) :: rho, velx, vely, velz
    real(kind=realType) :: ue, ve, we, qne, qnh

    select case (BCType(nn))

    case (DomainInterfaceAll)


       ! Flow variables are already prescribed from the coupler
       ! (see the subroutine setInterfaceData).

       ! Loop over the generic subface and copy the density, 
       ! velocities, pressure and turbulent variables.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ww1(i,j,irho) = BCData(nn)%rho(i,j)
          ww1(i,j,ivx)  = BCData(nn)%velx(i,j)
          ww1(i,j,ivy)  = BCData(nn)%vely(i,j)
          ww1(i,j,ivz)  = BCData(nn)%velz(i,j)
          pp1(i,j)      = BCData(nn)%ps(i,j)

          do l=nt1MG,nt2MG
             ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the energy for these halo's.

       call computeEtot2(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo2(correctForK)

       !============================================================

    case (DomainInterfaceRhoUVW)

       ! Flow variables are already prescribed from the coupler
       ! (see the subroutine setInterfaceData).

       ! Loop over the generic subface and copy the density,
       ! velocities and turbulent variables. The pressure is
       ! computed using Riemann invariants.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          rho  = BCData(nn)%rho(i,j)
          velx = BCData(nn)%velx(i,j)
          vely = BCData(nn)%vely(i,j)
          velz = BCData(nn)%velz(i,j)
          nnx  = BCData(nn)%norm(i,j,1)
          nny  = BCData(nn)%norm(i,j,2)
          nnz  = BCData(nn)%norm(i,j,3)

          ! Some abbreviations in which gamma occurs.

          gm1   = gamma2(i,j) - one
          ovgm1 = one/gm1

          ! Determine the acoustic Riemann variable that must be
          ! extrapolated from the domain.

          r    = one/ww2(i,j,irho)
          a2   = gamma2(i,j)*pp2(i,j)*r
          beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
               + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

          ! Compute the speed of sound squared in the halo.

          a2 = half*gm1*(beta - velx*nnx - vely*nny - velz*nnz)
          a2 = max(zero,a2)
          a2 = a2*a2

          ! Compute the pressure in the halo, assuming a
          ! constant value of gamma.

          pp1(i,j) = rho*a2/gamma2(i,j)

          ! Copy the density, velocities and turbulent variables.

          ww1(i,j,irho) = rho
          ww1(i,j,ivx)  = velx
          ww1(i,j,ivy)  = vely
          ww1(i,j,ivz)  = velz

          do l=nt1MG,nt2MG
             ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the energy for these halo's.

       call computeEtot2(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo2(correctForK)

       !============================================================

    case (DomainInterfaceP)

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Store a couple of variables, such as the static
          ! pressure and grid unit outward normal, a bit easier.

          pExit = BCData(nn)%ps(i,j)

          nnx = BCData(nn)%norm(i,j,1)
          nny = BCData(nn)%norm(i,j,2)
          nnz = BCData(nn)%norm(i,j,3)

          ! Abbreviate 1/gamma and 1/(gamma -1) a bit easier.

          ovg   = one/gamma2(i,j)
          ovgm1 = one/(gamma2(i,j)-one)

          ! Store the internal pressure and correct for the
          ! possible presence of a k-equation.

          pInt = pp2(i,j)
          if( correctForK ) &
               pInt = pInt - twoThird*ww2(i,j,irho)*ww2(i,j,itu1)

          ! Compute the velocity components, the normal velocity
          ! and the speed of sound for the internal cell.

          r   = one/ww2(i,j,irho)
          a2  = gamma2(i,j)*pInt*r
          a   = sqrt(a2)
          ue  = ww2(i,j,ivx)
          ve  = ww2(i,j,ivy)
          we  = ww2(i,j,ivz)
          qne = ue*nnx + ve*nny + we*nnz

          ! Compute the entropy and the acoustic variable.
          ! These Riemann invariants, as well as the tangential
          ! velocity components, are extrapolated.

          ss = pInt*(r**gamma2(i,j))
          ac = qne + two*a*ovgm1

          ! Compute the state in the halo.

          ww1(i,j,irho) = (pExit/ss)**ovg
          pp1(i,j)      = pExit
          a             = sqrt(gamma2(i,j)*pExit/ww1(i,j,irho))
          qnh           = ac - two*a*ovgm1
          ww1(i,j,ivx)  = ue + (qnh - qne)*nnx
          ww1(i,j,ivy)  = ve + (qnh - qne)*nny
          ww1(i,j,ivz)  = we + (qnh - qne)*nnz

          ! Extrapolate the primitive turbulent variables.

          do l=nt1MG,nt2MG
             ww1(i,j,l) = ww2(i,j,l)
          enddo

          ! Correct the pressure if a k-equation is present.

          if( correctForK )     &
               pp1(i,j) = pp1(i,j) &
               + twoThird*ww1(i,j,irho)*ww1(i,j,itu1)

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the energy for these halo's.

       call computeEtot2(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo2(correctForK)

       !============================================================

    case (DomainInterfaceRho)

       call terminate("analyzeString", &
            "DomainInterfaceRho not implemented yet")

       !============================================================

    case (DomainInterfaceTotal)

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Store a couple of variables, such as the total
          ! pressure, total temperature, total enthalpy, flow
          ! direction and grid unit outward normal, a bit easier.

          ptot = BCData(nn)%ptInlet(i,j)
          ttot = BCData(nn)%ttInlet(i,j)
          htot = BCData(nn)%htInlet(i,j)

          ssx  = BCData(nn)%flowXdirInlet(i,j)
          ssy  = BCData(nn)%flowYdirInlet(i,j)
          ssz  = BCData(nn)%flowZdirInlet(i,j)

          nnx = BCData(nn)%norm(i,j,1)
          nny = BCData(nn)%norm(i,j,2)
          nnz = BCData(nn)%norm(i,j,3)

          ! Some abbreviations in which gamma occurs.

          gm1   = gamma2(i,j) - one
          ovgm1 = one/gm1

          ! Determine the acoustic Riemann variable that must be
          ! extrapolated from the domain.

          r    = one/ww2(i,j,irho)
          a2   = gamma2(i,j)*pp2(i,j)*r
          beta = ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny  &
               + ww2(i,j,ivz)*nnz + two*ovgm1*sqrt(a2)

          ! Correct the value of the Riemann invariant if total
          ! enthalpy scaling must be applied. This scaling may
          ! be needed for stability if large gradients of the
          ! total temperature are prescribed.

          scaleFact = one
          if( hScalingInlet ) &
               scaleFact = sqrt(htot/(r*(ww2(i,j,irhoE) + pp2(i,j))))

          beta = beta*scaleFact

          ! Compute the value of a2 + 0.5*gm1*q2, which is the
          ! total speed of sound for constant cp. However, the
          ! expression below is also valid for variable cp,
          ! although a linearization around the value of the
          ! internal cell is performed.

          q2    = ww2(i,j,ivx)**2 + ww2(i,j,ivy)**2 &
               + ww2(i,j,ivz)**2
          a2tot = gm1*(htot - r*(ww2(i,j,irhoE) + pp2(i,j)) &
               +      half*q2) + a2

          ! Compute the dot product between the normal and the
          ! velocity direction. This value should be negative.

          alpha = nnx*ssx + nny*ssy + nnz*ssz

          ! Compute the coefficients in the quadratic equation
          ! for the magnitude of the velocity.

          aa2 =  half*gm1*alpha*alpha + one
          bb = -gm1*alpha*beta
          cc =  half*gm1*beta*beta - two*ovgm1*a2tot

          ! Solve the equation for the magnitude of the
          ! velocity. As this value must be positive and both aa2
          ! and bb are positive (alpha is negative and beta is
          ! positive up till Mach = 5.0 or so, which is not
          ! really subsonic anymore), it is clear which of the
          ! two possible solutions must be taken. Some clipping
          ! is present, but this is normally not active.

          dd = bb*bb - four*aa2*cc
          dd = sqrt(max(zero,dd))
          q  = (-bb + dd)/(two*aa2)
          q  = max(zero,q)
          q2 = q*q

          ! Compute the speed of sound squared from the total
          ! speed of sound equation (== total enthalpy equation
          ! for constant cp).

          a2 = a2tot - half*gm1*q2

          ! Compute the Mach number squared and cut it between
          ! 0.0 and 1.0. Adapt the velocity and speed of sound
          ! squared accordingly.

          m2 = q2/a2
          m2 = min(one,m2)
          q2 = m2*a2
          q  = sqrt(q2)
          a2 = a2tot - half*gm1*q2

          ! Compute the velocities in the halo cell and use rho,
          ! rhoe and p as temporary buffers to store the total
          ! temperature, total pressure and static temperature.

          ww1(i,j,ivx)  = q*ssx
          ww1(i,j,ivy)  = q*ssy
          ww1(i,j,ivz)  = q*ssz

          ww1(i,j,irho)  = ttot
          pp1(i,j)       = ptot
          ww1(i,j,irhoE) = a2/(gamma2(i,j)*RGas)

          ! Compute the turbulent variables, which are
          ! prescribed.

          do l=nt1MG,nt2MG
             ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the pressure and density for these halo's.
       ! This is identical to a subsonic inlet treatment,
       ! so call that routine

       call pRhoSubsonicInlet(ww1, pp1, correctForK)

       ! Compute the energy for these halo's.

       call computeEtot2(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo2(correctForK)

    end select
  end subroutine bcDomainInterface2
#endif

  subroutine bcFarfield2(nn, secondHalo, correctForK)

    !      ******************************************************************
    !      *                                                                *
    !      * bcFarfield applies the farfield boundary condition to a block. *
    !      * It is assumed that the BCPointers are already set              *
    !      *                                                                *
    !      ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK

    ! Local variables.
    integer(kind=intType) :: nn, i, j, k, l, ii

    real(kind=realType) :: nnx, nny, nnz
    real(kind=realType) :: gm1, ovgm1, ac1, ac2
    real(kind=realType) :: r0, u0, v0, w0, qn0, vn0, c0, s0
    real(kind=realType) :: re, ue, ve, we, qne, ce
    real(kind=realType) :: qnf, cf, uf, vf, wf, sf, cc, qq

    ! Some constants needed to compute the riemann inVariants.

    gm1   = gammaInf -one
    ovgm1 = one/gm1

    ! Compute the three velocity components, the speed of sound and
    ! the entropy of the free stream.

    r0  = one/wInf(irho)
    u0  = wInf(ivx)
    v0  = wInf(ivy)
    w0  = wInf(ivz)
    c0  = sqrt(gammaInf*pInfCorr*r0)
    s0  = wInf(irho)**gammaInf/pInfCorr

    ! Loop over the generic subface to set the state in the
    ! halo cells.
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ! Compute the normal velocity of the free stream and
       ! substract the normal velocity of the mesh.

       qn0 = u0*BCData(nn)%norm(i,j,1) + v0*BCData(nn)%norm(i,j,2) + w0*BCData(nn)%norm(i,j,3)
       vn0 = qn0 - BCData(nn)%rface(i,j)

       ! Compute the three velocity components, the normal
       ! velocity and the speed of sound of the current state
       ! in the internal cell.

       re  = one/ww2(i,j,irho)
       ue  = ww2(i,j,ivx)
       ve  = ww2(i,j,ivy)
       we  = ww2(i,j,ivz)
       qne = ue*BCData(nn)%norm(i,j,1) + ve*BCData(nn)%norm(i,j,2) + we*BCData(nn)%norm(i,j,3)
       ce  = sqrt(gamma2(i,j)*pp2(i,j)*re)

       ! Compute the new values of the riemann inVariants in
       ! the halo cell. Either the value in the internal cell
       ! is taken (positive sign of the corresponding
       ! eigenvalue) or the free stream value is taken
       ! (otherwise).

       if(vn0 > -c0) then ! Outflow or subsonic inflow.
          ac1 = qne + two*ovgm1*ce
       else               ! Supersonic inflow.
          ac1 = qn0 + two*ovgm1*c0
       endif

       if(vn0 > c0) then  ! Supersonic outflow.
          ac2 = qne - two*ovgm1*ce
       else                     ! Inflow or subsonic outflow.
          ac2 = qn0 - two*ovgm1*c0
       endif

       qnf = half*  (ac1 + ac2)
       cf  = fourth*(ac1 - ac2)*gm1

       if(vn0 > zero) then ! Outflow.

          uf = ue + (qnf - qne)*BCData(nn)%norm(i,j,1)
          vf = ve + (qnf - qne)*BCData(nn)%norm(i,j,2)
          wf = we + (qnf - qne)*BCData(nn)%norm(i,j,3)

          sf = ww2(i,j,irho)**gamma2(i,j)/pp2(i,j)

          !$AD II-LOOP
          do l=nt1MG,nt2MG
             ww1(i,j,l) = ww2(i,j,l)
          enddo

       else                          
          ! Inflow
          uf = u0 + (qnf - qn0)*BCData(nn)%norm(i,j,1)
          vf = v0 + (qnf - qn0)*BCData(nn)%norm(i,j,2)
          wf = w0 + (qnf - qn0)*BCData(nn)%norm(i,j,3)
          sf = s0

          !$AD II-LOOP
          do l=nt1MG,nt2MG
             ww1(i,j,l) = wInf(l)
          enddo
       endif

       ! Compute the density, velocity and pressure in the
       ! halo cell.

       cc = cf*cf/gamma2(i,j)
       qq = uf*uf + vf*vf + wf*wf
       ww1(i,j,irho) = (sf*cc)**ovgm1
       ww1(i,j,ivx)  = uf
       ww1(i,j,ivy)  = vf
       ww1(i,j,ivz)  = wf
       pp1(i,j)      = ww1(i,j,irho)*cc

       ! Simply set the laminar and eddy viscosity to
       ! the value in the donor cell. Their values do
       ! not matter too much in the far field.

       if( viscous )    rlv1(i,j) = rlv2(i,j)
       if( eddyModel )  rev1(i,j) = rev2(i,j)
    enddo

    ! Compute the energy for these halo's.
    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.
    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcFarfield2
#ifndef USE_TAPENADE
  subroutine bcSupersonicInflow2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * bcSupersonicInflow applies the supersonic inflow boundary      *
    ! * conditions, entire state vector is prescribed, to a block. It  *
    ! * is assumed that the pointers in blockPointers are already set  *
    ! * to the correct block on the correct grid level.                *
    ! *                                                                *
    ! ******************************************************************

    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, kk, mm, ii
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd


    ! Loop over the generic subface to set the state in the
    ! halo cells.
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ww1(i,j,irho) = BCData(nn)%rho(i,j)
       ww1(i,j,ivx)  = BCData(nn)%velx(i,j)
       ww1(i,j,ivy)  = BCData(nn)%vely(i,j)
       ww1(i,j,ivz)  = BCData(nn)%velz(i,j)
       pp1(i,j)      = BCData(nn)%ps(i,j)

       ! The turbulent variables.

       do l=nt1MG,nt2MG
          ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
       enddo

       ! Set the laminar and eddy viscosity in the halo
       ! if needed.

       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)
    end do

    call computeEtot2(ww1, pp1, correctForK)

    if (secondHalo) then 
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ww0(i,j,irho) = BCData(nn)%rho(i,j)
          ww0(i,j,ivx)  = BCData(nn)%velx(i,j)
          ww0(i,j,ivy)  = BCData(nn)%vely(i,j)
          ww0(i,j,ivz)  = BCData(nn)%velz(i,j)
          pp0(i,j)      = BCData(nn)%ps(i,j)

          ! The turbulent variables.

          do l=nt1MG,nt2MG
             ww0(i,j,l) = BCData(nn)%turbInlet(i,j,l)
          enddo

          ! Set the laminar and eddy viscosity in the halo
          ! if needed.

          if( viscous )   rlv0(i,j) = rlv1(i,j)
          if( eddyModel ) rev0(i,j) = rev1(i,j)
       end do

       call computeEtot2(ww0, pp0, correctForK)
    end if

  end subroutine bcSupersonicInflow2

  subroutine bcExtrap2(nn, secondHalo, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * ccExtrap applies the extrapolation boundary condition to a     *
    ! * block. It is assumed that the pointers in blockPointers are    *
    ! * already set to the correct block on the correct grid level.    *
    ! * Extrapolation boundaries are applied to both singular lines or *
    ! * points of a block face and to supersonic outlets. They are     *
    ! * marked differently because of postprocessing reasons, but      *
    ! * their numerical treatment is identical.                        *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use BCTypes
    use constants
    use flowVarRefState
    use inputDiscretization
    use inputPhysics
    use iteration
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local parameter.
    real(kind=realType), parameter :: factor = 0.5

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType) :: fw2, fw3

    ! Set the extrapolation weights, depending on the situation.

    if(BCType(nn) == SupersonicOutflow) then

       ! A physical outflow face. Set the weights depending
       ! on the input parameter.

       select case (outflowTreatment)
       case (constantExtrapol)
          fw2 = one; fw3 = zero
       case (linExtrapol)
          fw2 = two; fw3 = -one
       end select

    else

       ! Singular block boundary. Use linear extrapolation.

       fw2 = two; fw3 = -one

    endif

    ! Loop over the generic subface to set the state in the
    ! 1-st level halos 
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ! Extrapolate the density, velocities and pressure.
       ! Make sure that a certain threshold is kept for the
       ! density and pressure.

       ww1(i,j,irho) = fw2*ww2(i,j,irho) + fw3*ww3(i,j,irho)
       ww1(i,j,irho) = max(factor*ww2(i,j,irho), ww1(i,j,irho))

       ww1(i,j,ivx) = fw2*ww2(i,j,ivx) + fw3*ww3(i,j,ivx)
       ww1(i,j,ivy) = fw2*ww2(i,j,ivy) + fw3*ww3(i,j,ivy)
       ww1(i,j,ivz) = fw2*ww2(i,j,ivz) + fw3*ww3(i,j,ivz)

       pp1(i,j) = fw2*pp2(i,j) + fw3*pp3(i,j)
       pp1(i,j) = max(factor*pp2(i,j), pp1(i,j))

       ! Extrapolate the turbulent variables.

       do l=nt1MG,nt2MG
          ww1(i,j,l) = fw2*ww2(i,j,l) + fw3*ww3(i,j,l)
       enddo

       ! The laminar and eddy viscosity, if present. These
       ! values are simply taken constant. Their values do
       ! not really matter.

       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)

    enddo

    ! Compute the energy for these halo's.

    call computeEtot2(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo2(correctForK)

  end subroutine bcExtrap2

  subroutine pRhoSubsonicInlet(ww, pp, correctForK)
    !
    ! ******************************************************************
    ! *                                                                *
    ! * pRhoSubsonicInlet computes the pressure and density for the    *
    ! * given range of the block to which the pointers in              *
    ! * blockPointers currently point.                                 *
    ! *                                                                *
    ! ******************************************************************
    !
    use blockPointers
    use cpCurveFits
    use flowVarRefState
    use inputPhysics
    implicit none

    ! Local parameter.
    real(kind=realType), parameter :: twoThird = two*third

    ! Subroutine arguments.
    real(kind=realType), dimension(:,:,:) :: ww
    real(kind=realType), dimension(:,:) :: pp
    logical, intent(in)               :: correctForK

    ! Local variables.
    integer(kind=intType) :: i, j, ii, mm, nns, nnt, iii
    real(kind=realType)   :: govgm1, tt, ts, pt, ratio
    real(kind=realType)   :: intTs, intTt, val

    ! Determine the cp model used in the computation.
    select case (cpModel)

    case (cpConstant)

       ! Constant cp and thus constant gamma. Compute the coefficient
       ! gamma/(gamma-1), which occurs in the isentropic expression
       ! for the total pressure.

       govgm1 = gammaConstant/(gammaConstant - one)

       ! Loop over the pointer range
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1

          ! Store the total temperature, total pressure and
          ! static temperature a bit easier.

          tt = ww(i,j,irho)
          pt = pp(i,j)
          ts = ww(i,j,irhoE)

          ! Compute the static pressure from the total pressure
          ! and the temperature ratio. Compute the density using
          ! the gas law.

          ratio        = (ts/tt)**govgm1
          pp(i,j)      = pt*ratio
          ww(i,j,irho) = pp(i,j)/(RGas*ts)

       enddo

       !        ================================================================

    case (cpTempCurveFits)

       ! Cp as function of the temperature is given via curve fits.
       ! The ratio pt/ps is given by exp(a), where a is the integral
       ! from ts to tt of cp/(r*t).

       ! Loop over the pointer range
       !$AD II-LOOP
       do iii=0,isize*jsize-1
          i = mod(iii, isize) + 1
          j = iii/isize + 1

          ! Store the total temperature, total pressure and
          ! static temperature a bit easier. Note that the
          ! temperatures get their dimensional value.

          tt = Tref*ww(i,j,irho)
          pt = pp(i,j)
          ts = Tref*ww(i,j,irhoE)

          ! Determine the integrant of cp/(r*t) for the static
          ! temperature ts and the total temperature tt.

          call cportIntegrant(ts, nns, intTs)
          call cportIntegrant(tt, nnt, intTt)

          ! Compute the value of the integral of cp/(r*t) from
          ! ts to tt. First part is the initialization where it
          ! is assumed that both ts and tt lie in the same
          ! interval of the curve fits.

          val = intTt - intTs

          ! Correct this value if ts and tt belong to different
          ! curve fit intervals.

          do mm=(nns+1),nnt

             ! The contribution from the interval mm-1. Add the
             ! value of the integrant at the upper boundary.

             ii = mm - 1
             if(ii == 0_intType) then
                val = val + (cv0+one)*log(cpTrange(0))
             else
                val = val + cpTempFit(ii)%intCpovrt_2
             endif

             ! The contribution from the interval mm. Substract
             ! the value of integrant at the lower boundary.

             if(mm > cpNparts) then
                val = val - (cvn+one)*log(cpTrange(cpNparts))
             else
                val = val - cpTempFit(mm)%intCpovrt_1
             endif

          enddo

          ! Compute the static pressure from the known
          ! total pressure.

          ratio   = exp(val)
          pp(i,j) = pt/ratio

          ! Compute the density using the gas law.

          ts = ww(i,j,irhoE)
          ww(i,j,irho) = pp(i,j)/(RGas*ts)

       enddo
    end select

    ! Add 2*rho*k/3 to the pressure if a k-equation is present.

    if( correctForK ) then
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1
          pp(i,j) = pp(i,j) &
               + twoThird*ww(i,j,irho)*ww(i,j,itu1)
       enddo
    end if

  contains
    subroutine cportIntegrant(T, nn, int)
      !
      ! ****************************************************************
      ! *                                                              *
      ! * cportIntegrant computes the integrant of the function        *
      ! * cp/(r*t) for the given temperature. It also stores the       *
      ! * correct curve fit interval, which is needed to determine the *
      ! * entire integral in the main subroutine.                      *
      ! *                                                              *
      ! ****************************************************************
      !
      implicit none

      ! Subroutine arguments.
      integer(kind=intType), intent(out) :: nn
      real(kind=realType),   intent(in)  :: t
      real(kind=realType),   intent(out) :: int

      ! Local variables.
      integer(kind=intType) :: mm, ii, start
      real(kind=realType)   :: T2

      ! Determine the situation we are having here for the temperature.

      if(T <= cpTrange(0)) then

         ! Temperature is less than the smallest temperature of the
         ! curve fits. Use extrapolation using constant cp.
         ! Set nn to 0 to indicate this.

         nn  = 0
         int = (cv0+one)*log(T)

      else if(T >= cpTrange(cpNparts)) then

         ! Temperature is larger than the largest temperature of the
         ! curve fits. Use extrapolation using constant cp.
         ! Set nn to cpNparts+1 to indicate this.

         nn  = cpNparts + 1
         int = (cvn+one)*log(T)

      else

         ! Temperature is within the curve fit range. Determine
         ! the correct interval.

         ii    = cpNparts
         start = 1
         interval: do

            ! Next guess for the interval.

            nn = start + ii/2

            ! Determine the situation we are having here.

            if(T > cpTrange(nn)) then

               ! Temperature is larger than the upper boundary of
               ! the current interval. Update the lower boundary.

               start = nn + 1
               ii    = ii - 1

            else if(T >= cpTrange(nn-1)) then

               ! This is the correct range. Exit the do-loop.

               exit

            endif

            ! Modify ii for the next branch to search.

            ii = ii/2

         enddo interval

         ! Nn contains the correct curve fit interval.
         ! Compute the value of the integrant.

         int = zero
         do ii=1,cpTempFit(nn)%nterm

            mm = cpTempFit(nn)%exponents(ii)
            if(mm == 0_intType) then
               int = int + cpTempFit(nn)%constants(ii)*log(T)
            else
               T2  = T**mm
               int = int + cpTempFit(nn)%constants(ii)*T2/mm
            endif

         enddo

      endif
    end subroutine cportIntegrant
  end subroutine pRhoSubsonicInlet
#endif

  subroutine computeEtot2(ww, pp, correctForK)
    ! Simplified total energy computation for boundary conditions.
    ! Only implements the constant cpModel

    use constants
    use flowVarRefState
    use inputPhysics
    implicit none 

    real(kind=realType),  dimension(:,:) :: pp
    real(kind=realType),  dimension(:,:,:) :: ww
    logical :: correctForK
    integer(kind=intType) :: ii, i, j
    real(kind=realType) :: ovgm1, factK

    select case (cpModel) 

    case (cpConstant)

       ! Constant cp and thus constant gamma.
       ! Abbreviate 1/(gamma -1) a bit easier.

       ovgm1 = one/(gammaConstant - one)
       factK = ovgm1*(five*third - gammaConstant)

       ! Loop over the given array and compute the energy, possibly
       ! correcting for K
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + 1
          j = ii/isize + 1
          if( .not. correctForK ) then
             ww(i,j,iRhoE) = ovgm1*pp(i,j) &
                  + half*ww(i,j,irho)*(ww(i,j,ivx)**2 &
                  +                    ww(i,j,ivy)**2 &
                  +                    ww(i,j,ivz)**2)

          else
             ww(i,j, iRhoE) = ovgm1*pp(i,j) &
                  + half*ww(i,j,irho)*(ww(i,j,ivx)**2 &
                  +                    ww(i,j,ivy)**2 &
                  +                    ww(i,j,ivz)**2) & 
                  - factK*ww(i,j,irho)*ww(i,j,itu1)
          end if
       enddo

    case (cpTempCurveFits)

       call terminate("BCRoutines", "CPTempCurveFits not implemented yet.")
    end select
  end subroutine computeEtot2

  subroutine extrapolate2ndHalo2(correctForK)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * extrapolate2ndHalo determines the states of the second layer   *
    !      * halo cells for the given subface of the block. It is assumed   *
    !      * that the appropriate BCPointers are already set
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use constants
    use flowVarRefState
    use iteration
    use inputPhysics
    implicit none
    logical, intent(in) :: correctForK
    real(kind=realType), parameter :: factor = 0.5_realType
    integer(kind=intType) :: i, j, l, ii

    ! Loop over the generic subface to set the state in the
    ! halo cells.
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + 1
       j = ii/isize + 1

       ! Extrapolate the density, momentum and pressure.
       ! Make sure that a certain threshold is kept.

       ww0(i,j,irho) = two*ww1(i,j,irho) - ww2(i,j,irho)
       ww0(i,j,irho) = max(factor*ww1(i,j,irho),ww0(i,j,irho))

       ww0(i,j,ivx) = two*ww1(i,j,ivx) - ww2(i,j,ivx)
       ww0(i,j,ivy) = two*ww1(i,j,ivy) - ww2(i,j,ivy)
       ww0(i,j,ivz) = two*ww1(i,j,ivz) - ww2(i,j,ivz)

       pp0(i,j) = max(factor*pp1(i,j),two*pp1(i,j) - pp2(i,j))

       ! Extrapolate the turbulent variables. Use constant
       ! extrapolation.

       !$AD II-LOOP
       do l=nt1MG,nt2MG
          ww0(i,j,l) = ww1(i,j,l)
       enddo

       ! The laminar and eddy viscosity, if present. These values
       ! are simply taken constant. Their values do not matter.

       if( viscous )   rlv0(i,j) = rlv1(i,j)
       if( eddyModel ) rev0(i,j) = rev1(i,j)
    enddo

    ! Compute the energy for this halo range.
    call computeEtot2(ww0, pp0, correctForK)

  end subroutine extrapolate2ndHalo2

  subroutine setBCPointers2(nn)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setBCPointers sets the pointers needed for the boundary        *
    !      * condition treatment on a general face, such that the boundary  *
    !      * routines are only implemented once instead of 6 times.         *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use flowVarRefState
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables
    integer(kind=intType) :: iStart, iEnd, jStart, jEnd

    ! Determine the sizes of each face and point to just the range we
    ! need on each face. 
    iStart = BCData(nn)%icBeg
    iEnd   = BCData(nn)%icEnd
    jStart = BCData(nn)%jcBeg
    jEnd   = BCData(nn)%jcEnd

    ! Set the size of the subface
    isize = iEnd-iStart + 1
    jsize = jEnd-jStart + 1

    ! Determine the face id on which the subface is located and set
    ! the pointers accordinly.
#ifndef TAPENADE_REVERSE
    select case (BCFaceID(nn))

       !===============================================================
    case (iMin)
       ww3 => w(3, iStart:iEnd, jStart:jend, :)
       ww2 => w(2, iStart:iEnd, jStart:jend, :)
       ww1 => w(1, iStart:iEnd, jStart:jend, :)
       ww0 => w(0, iStart:iEnd, jStart:jend, :)

       pp3 => p(3, iStart:iEnd, jStart:jend)
       pp2 => p(2, iStart:iEnd, jStart:jend)
       pp1 => p(1, iStart:iEnd, jStart:jend)
       pp0 => p(0, iStart:iEnd, jStart:jend)

       rlv3 => rlv(3, iStart:iEnd, jStart:jend)
       rlv2 => rlv(2, iStart:iEnd, jStart:jend)
       rlv1 => rlv(1, iStart:iEnd, jStart:jend)
       rlv0 => rlv(0, iStart:iEnd, jStart:jend)

       rev3 => rev(3, iStart:iEnd, jStart:jend)
       rev2 => rev(2, iStart:iEnd, jStart:jend)
       rev1 => rev(1, iStart:iEnd, jStart:jend)
       rev0 => rev(0, iStart:iEnd, jStart:jend)

       gamma3 => gamma(3, iStart:iEnd, jStart:jend)
       gamma2 => gamma(2, iStart:iEnd, jStart:jend)
       gamma1 => gamma(1, iStart:iEnd, jStart:jend)
       gamma0 => gamma(0, iStart:iEnd, jStart:jend)

       !===============================================================

    case (iMax)

       ww3 => w(nx, iStart:iEnd, jStart:jend, :)
       ww2 => w(il, iStart:iEnd, jStart:jend, :)
       ww1 => w(ie, iStart:iEnd, jStart:jend, :)
       ww0 => w(ib, iStart:iEnd, jStart:jend, :)

       pp3 => p(nx, iStart:iEnd, jStart:jend)
       pp2 => p(il, iStart:iEnd, jStart:jend)
       pp1 => p(ie, iStart:iEnd, jStart:jend)
       pp0 => p(ib, iStart:iEnd, jStart:jend)

       rlv3 => rlv(nx, iStart:iEnd, jStart:jend)
       rlv2 => rlv(il, iStart:iEnd, jStart:jend)
       rlv1 => rlv(ie, iStart:iEnd, jStart:jend)
       rlv0 => rlv(ib, iStart:iEnd, jStart:jend)

       rev3 => rev(nx, iStart:iEnd, jStart:jend)
       rev2 => rev(il, iStart:iEnd, jStart:jend)
       rev1 => rev(ie, iStart:iEnd, jStart:jend)
       rev0 => rev(ib, iStart:iEnd, jStart:jend)

       gamma3 => gamma(nx, iStart:iEnd, jStart:jend)
       gamma2 => gamma(il, iStart:iEnd, jStart:jend)
       gamma1 => gamma(ie, iStart:iEnd, jStart:jend)
       gamma0 => gamma(ib, iStart:iEnd, jStart:jend)

       !===============================================================

    case (jMin)

       ww3 => w(iStart:iEnd, 3, jStart:jend, :)
       ww2 => w(iStart:iEnd, 2, jStart:jend, :)
       ww1 => w(iStart:iEnd, 1, jStart:jend, :)
       ww0 => w(iStart:iEnd, 0, jStart:jend, :)

       pp3 => p(iStart:iEnd, 3, jStart:jend)
       pp2 => p(iStart:iEnd, 2, jStart:jend)
       pp1 => p(iStart:iEnd, 1, jStart:jend)
       pp0 => p(iStart:iEnd, 0, jStart:jend)

       rlv3 => rlv(iStart:iEnd, 3, jStart:jend)
       rlv2 => rlv(iStart:iEnd, 2, jStart:jend)
       rlv1 => rlv(iStart:iEnd, 1, jStart:jend)
       rlv0 => rlv(iStart:iEnd, 0, jStart:jend)

       rev3 => rev(iStart:iEnd, 3, jStart:jend)
       rev2 => rev(iStart:iEnd, 2, jStart:jend)
       rev1 => rev(iStart:iEnd, 1, jStart:jend)
       rev0 => rev(iStart:iEnd, 0, jStart:jend)

       gamma3 => gamma(iStart:iEnd, 3, jStart:jend)
       gamma2 => gamma(iStart:iEnd, 2, jStart:jend)
       gamma1 => gamma(iStart:iEnd, 1, jStart:jend)
       gamma0 => gamma(iStart:iEnd, 0, jStart:jend)

       !===============================================================

    case (jMax)

       ww3 => w(iStart:iEnd, ny, jStart:jend, :)
       ww2 => w(iStart:iEnd, jl, jStart:jend, :)
       ww1 => w(iStart:iEnd, je, jStart:jend, :)
       ww0 => w(iStart:iEnd, jb, jStart:jend, :)

       pp3 => p(iStart:iEnd, ny, jStart:jend)
       pp2 => p(iStart:iEnd, jl, jStart:jend)
       pp1 => p(iStart:iEnd, je, jStart:jend)
       pp0 => p(iStart:iEnd, jb, jStart:jend)

       rlv3 => rlv(iStart:iEnd, ny, jStart:jend)
       rlv2 => rlv(iStart:iEnd, jl, jStart:jend)
       rlv1 => rlv(iStart:iEnd, je, jStart:jend)
       rlv0 => rlv(iStart:iEnd, jb, jStart:jend)

       rev3 => rev(iStart:iEnd, ny, jStart:jend)
       rev2 => rev(iStart:iEnd, jl, jStart:jend)
       rev1 => rev(iStart:iEnd, je, jStart:jend)
       rev0 => rev(iStart:iEnd, jb, jStart:jend)

       gamma3 => gamma(iStart:iEnd, ny, jStart:jend)
       gamma2 => gamma(iStart:iEnd, jl, jStart:jend)
       gamma1 => gamma(iStart:iEnd, je, jStart:jend)
       gamma0 => gamma(iStart:iEnd, jb, jStart:jend)

       !===============================================================

    case (kMin)

       ww3 => w(iStart:iEnd, jStart:jend, 3, :)
       ww2 => w(iStart:iEnd, jStart:jend, 2, :)
       ww1 => w(iStart:iEnd, jStart:jend, 1, :)
       ww0 => w(iStart:iEnd, jStart:jend, 0, :)

       pp3 => p(iStart:iEnd, jStart:jend, 3)
       pp2 => p(iStart:iEnd, jStart:jend, 2)
       pp1 => p(iStart:iEnd, jStart:jend, 1)
       pp0 => p(iStart:iEnd, jStart:jend, 0)

       rlv3 => rlv(iStart:iEnd, jStart:jend, 3)
       rlv2 => rlv(iStart:iEnd, jStart:jend, 2)
       rlv1 => rlv(iStart:iEnd, jStart:jend, 1)
       rlv0 => rlv(iStart:iEnd, jStart:jend, 0)

       rev3 => rev(iStart:iEnd, jStart:jend, 3)
       rev2 => rev(iStart:iEnd, jStart:jend, 2)
       rev1 => rev(iStart:iEnd, jStart:jend, 1)
       rev0 => rev(iStart:iEnd, jStart:jend, 0)

       gamma3 => gamma(iStart:iEnd, jStart:jend, 3)
       gamma2 => gamma(iStart:iEnd, jStart:jend, 2)
       gamma1 => gamma(iStart:iEnd, jStart:jend, 1)
       gamma0 => gamma(iStart:iEnd, jStart:jend, 0)

       !===============================================================

    case (kMax)

       ww3 => w(iStart:iEnd, jStart:jend, nz, :)
       ww2 => w(iStart:iEnd, jStart:jend, kl, :)
       ww1 => w(iStart:iEnd, jStart:jend, ke, :)
       ww0 => w(iStart:iEnd, jStart:jend, kb, :)

       pp3 => p(iStart:iEnd, jStart:jend, nz)
       pp2 => p(iStart:iEnd, jStart:jend, kl)
       pp1 => p(iStart:iEnd, jStart:jend, ke)
       pp0 => p(iStart:iEnd, jStart:jend, kb)

       rlv3 => rlv(iStart:iEnd, jStart:jend, nz)
       rlv2 => rlv(iStart:iEnd, jStart:jend, kl)
       rlv1 => rlv(iStart:iEnd, jStart:jend, ke)
       rlv0 => rlv(iStart:iEnd, jStart:jend, kb)

       rev3 => rev(iStart:iEnd, jStart:jend, nz)
       rev2 => rev(iStart:iEnd, jStart:jend, kl)
       rev1 => rev(iStart:iEnd, jStart:jend, ke)
       rev0 => rev(iStart:iEnd, jStart:jend, kb)

       gamma3 => gamma(iStart:iEnd, jStart:jend, nz)
       gamma2 => gamma(iStart:iEnd, jStart:jend, kl)
       gamma1 => gamma(iStart:iEnd, jStart:jend, ke)
       gamma0 => gamma(iStart:iEnd, jStart:jend, kb)

    end select

    select case (BCFaceID(nn))
    case (iMin)
       ssi => si(1,:,:,:)
       ssj => sj(2,:,:,:)
       ssk => sk(2,:,:,:)
       ss  => s (2,:,:,:)
    case (iMax)
       ssi => si(il,:,:,:)
       ssj => sj(il,:,:,:)
       ssk => sk(il,:,:,:)
       ss  =>  s(il,:,:,:)
    case (jMin)
       ssi => sj(:,1,:,:)
       ssj => si(:,2,:,:)
       ssk => sk(:,2,:,:)
       ss   => s(:,2,:,:)
    case (jMax)
       ssi => sj(:,jl,:,:)
       ssj => si(:,jl,:,:)
       ssk => sk(:,jl,:,:)
       ss  =>  s(:,jl,:,:)
    case (kMin)
       ssi => sk(:,:,1,:)
       ssj => si(:,:,2,:)
       ssk => sj(:,:,2,:)
       ss  =>  s(:,:,2,:)
    case (kMax)
       ssi => sk(:,:,kl,:)
       ssj => si(:,:,kl,:)
       ssk => sj(:,:,kl,:)
       ss  =>  s(:,:,kl,:)
    end select

    select case (BCFaceID(nn))
    case (iMin)
       xline => x(1,:,:,:)
    case (iMax)
       xline => x(il,:,:,:)
    case (jMin)
       xline => x(:,1,:,:)
    case (jMax)
       xline => x(:,jl,:,:)
    case (kMin)
       xline => x(:,:,1,:)
    case (kMax)
       xline => x(:,:,kl,:)
    end select
#else
    select case (BCFaceID(nn))

       !===============================================================
    case (iMin)
       ww3(1:isize,1:jsize,:) = w(3, iStart:iEnd, jStart:jend, :)
       ww2(1:isize,1:jsize,:) = w(2, iStart:iEnd, jStart:jend, :)
       ww1(1:isize,1:jsize,:) = w(1, iStart:iEnd, jStart:jend, :)
       ww0(1:isize,1:jsize,:) = w(0, iStart:iEnd, jStart:jend, :)

       pp3(1:isize,1:jsize) = p(3, iStart:iEnd, jStart:jend)
       pp2(1:isize,1:jsize) = p(2, iStart:iEnd, jStart:jend)
       pp1(1:isize,1:jsize) = p(1, iStart:iEnd, jStart:jend)
       pp0(1:isize,1:jsize) = p(0, iStart:iEnd, jStart:jend)

       rlv3(1:isize,1:jsize) = rlv(3, iStart:iEnd, jStart:jend)
       rlv2(1:isize,1:jsize) = rlv(2, iStart:iEnd, jStart:jend)
       rlv1(1:isize,1:jsize) = rlv(1, iStart:iEnd, jStart:jend)
       rlv0(1:isize,1:jsize) = rlv(0, iStart:iEnd, jStart:jend)

       rev3(1:isize,1:jsize) = rev(3, iStart:iEnd, jStart:jend)
       rev2(1:isize,1:jsize) = rev(2, iStart:iEnd, jStart:jend)
       rev1(1:isize,1:jsize) = rev(1, iStart:iEnd, jStart:jend)
       rev0(1:isize,1:jsize) = rev(0, iStart:iEnd, jStart:jend)

       gamma3(1:isize,1:jsize) = gamma(3, iStart:iEnd, jStart:jend)
       gamma2(1:isize,1:jsize) = gamma(2, iStart:iEnd, jStart:jend)
       gamma1(1:isize,1:jsize) = gamma(1, iStart:iEnd, jStart:jend)
       gamma0(1:isize,1:jsize) = gamma(0, iStart:iEnd, jStart:jend)

       !===============================================================

    case (iMax)

       ww3(1:isize,1:jsize,:) = w(nx, iStart:iEnd, jStart:jend, :)
       ww2(1:isize,1:jsize,:) = w(il, iStart:iEnd, jStart:jend, :)
       ww1(1:isize,1:jsize,:) = w(ie, iStart:iEnd, jStart:jend, :)
       ww0(1:isize,1:jsize,:) = w(ib, iStart:iEnd, jStart:jend, :)

       pp3(1:isize,1:jsize) = p(nx, iStart:iEnd, jStart:jend)
       pp2(1:isize,1:jsize) = p(il, iStart:iEnd, jStart:jend)
       pp1(1:isize,1:jsize) = p(ie, iStart:iEnd, jStart:jend)
       pp0(1:isize,1:jsize) = p(ib, iStart:iEnd, jStart:jend)

       rlv3(1:isize,1:jsize) = rlv(nx, iStart:iEnd, jStart:jend)
       rlv2(1:isize,1:jsize) = rlv(il, iStart:iEnd, jStart:jend)
       rlv1(1:isize,1:jsize) = rlv(ie, iStart:iEnd, jStart:jend)
       rlv0(1:isize,1:jsize) = rlv(ib, iStart:iEnd, jStart:jend)

       rev3(1:isize,1:jsize) = rev(nx, iStart:iEnd, jStart:jend)
       rev2(1:isize,1:jsize) = rev(il, iStart:iEnd, jStart:jend)
       rev1(1:isize,1:jsize) = rev(ie, iStart:iEnd, jStart:jend)
       rev0(1:isize,1:jsize) = rev(ib, iStart:iEnd, jStart:jend)

       gamma3(1:isize,1:jsize) = gamma(nx, iStart:iEnd, jStart:jend)
       gamma2(1:isize,1:jsize) = gamma(il, iStart:iEnd, jStart:jend)
       gamma1(1:isize,1:jsize) = gamma(ie, iStart:iEnd, jStart:jend)
       gamma0(1:isize,1:jsize) = gamma(ib, iStart:iEnd, jStart:jend)

       !===============================================================

    case (jMin)

       ww3(1:isize,1:jsize,:) = w(iStart:iEnd, 3, jStart:jend, :)
       ww2(1:isize,1:jsize,:) = w(iStart:iEnd, 2, jStart:jend, :)
       ww1(1:isize,1:jsize,:) = w(iStart:iEnd, 1, jStart:jend, :)
       ww0(1:isize,1:jsize,:) = w(iStart:iEnd, 0, jStart:jend, :)

       pp3(1:isize,1:jsize) = p(iStart:iEnd, 3, jStart:jend)
       pp2(1:isize,1:jsize) = p(iStart:iEnd, 2, jStart:jend)
       pp1(1:isize,1:jsize) = p(iStart:iEnd, 1, jStart:jend)
       pp0(1:isize,1:jsize) = p(iStart:iEnd, 0, jStart:jend)

       rlv3(1:isize,1:jsize) = rlv(iStart:iEnd, 3, jStart:jend)
       rlv2(1:isize,1:jsize) = rlv(iStart:iEnd, 2, jStart:jend)
       rlv1(1:isize,1:jsize) = rlv(iStart:iEnd, 1, jStart:jend)
       rlv0(1:isize,1:jsize) = rlv(iStart:iEnd, 0, jStart:jend)

       rev3(1:isize,1:jsize) = rev(iStart:iEnd, 3, jStart:jend)
       rev2(1:isize,1:jsize) = rev(iStart:iEnd, 2, jStart:jend)
       rev1(1:isize,1:jsize) = rev(iStart:iEnd, 1, jStart:jend)
       rev0(1:isize,1:jsize) = rev(iStart:iEnd, 0, jStart:jend)

       gamma3(1:isize,1:jsize) = gamma(iStart:iEnd, 3, jStart:jend)
       gamma2(1:isize,1:jsize) = gamma(iStart:iEnd, 2, jStart:jend)
       gamma1(1:isize,1:jsize) = gamma(iStart:iEnd, 1, jStart:jend)
       gamma0(1:isize,1:jsize) = gamma(iStart:iEnd, 0, jStart:jend)

       !===============================================================

    case (jMax)

       ww3(1:isize,1:jsize,:) = w(iStart:iEnd, ny, jStart:jend, :)
       ww2(1:isize,1:jsize,:) = w(iStart:iEnd, jl, jStart:jend, :)
       ww1(1:isize,1:jsize,:) = w(iStart:iEnd, je, jStart:jend, :)
       ww0(1:isize,1:jsize,:) = w(iStart:iEnd, jb, jStart:jend, :)

       pp3(1:isize,1:jsize) = p(iStart:iEnd, ny, jStart:jend)
       pp2(1:isize,1:jsize) = p(iStart:iEnd, jl, jStart:jend)
       pp1(1:isize,1:jsize) = p(iStart:iEnd, je, jStart:jend)
       pp0(1:isize,1:jsize) = p(iStart:iEnd, jb, jStart:jend)

       rlv3(1:isize,1:jsize) = rlv(iStart:iEnd, ny, jStart:jend)
       rlv2(1:isize,1:jsize) = rlv(iStart:iEnd, jl, jStart:jend)
       rlv1(1:isize,1:jsize) = rlv(iStart:iEnd, je, jStart:jend)
       rlv0(1:isize,1:jsize) = rlv(iStart:iEnd, jb, jStart:jend)

       rev3(1:isize,1:jsize) = rev(iStart:iEnd, ny, jStart:jend)
       rev2(1:isize,1:jsize) = rev(iStart:iEnd, jl, jStart:jend)
       rev1(1:isize,1:jsize) = rev(iStart:iEnd, je, jStart:jend)
       rev0(1:isize,1:jsize) = rev(iStart:iEnd, jb, jStart:jend)

       gamma3(1:isize,1:jsize) = gamma(iStart:iEnd, ny, jStart:jend)
       gamma2(1:isize,1:jsize) = gamma(iStart:iEnd, jl, jStart:jend)
       gamma1(1:isize,1:jsize) = gamma(iStart:iEnd, je, jStart:jend)
       gamma0(1:isize,1:jsize) = gamma(iStart:iEnd, jb, jStart:jend)

       !===============================================================

    case (kMin)

       ww3(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, 3, :)
       ww2(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, 2, :)
       ww1(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, 1, :)
       ww0(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, 0, :)

       pp3(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, 3)
       pp2(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, 2)
       pp1(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, 1)
       pp0(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, 0)

       rlv3(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, 3)
       rlv2(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, 2)
       rlv1(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, 1)
       rlv0(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, 0)

       rev3(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, 3)
       rev2(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, 2)
       rev1(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, 1)
       rev0(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, 0)

       gamma3(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, 3)
       gamma2(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, 2)
       gamma1(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, 1)
       gamma0(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, 0)

       !===============================================================

    case (kMax)

       ww3(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, nz, :)
       ww2(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, kl, :)
       ww1(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, ke, :)
       ww0(1:isize,1:jsize,:) = w(iStart:iEnd, jStart:jend, kb, :)

       pp3(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, nz)
       pp2(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, kl)
       pp1(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, ke)
       pp0(1:isize,1:jsize) = p(iStart:iEnd, jStart:jend, kb)

       rlv3(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, nz)
       rlv2(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, kl)
       rlv1(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, ke)
       rlv0(1:isize,1:jsize) = rlv(iStart:iEnd, jStart:jend, kb)

       rev3(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, nz)
       rev2(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, kl)
       rev1(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, ke)
       rev0(1:isize,1:jsize) = rev(iStart:iEnd, jStart:jend, kb)

       gamma3(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, nz)
       gamma2(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, kl)
       gamma1(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, ke)
       gamma0(1:isize,1:jsize) = gamma(iStart:iEnd, jStart:jend, kb)
    end select

    ! Note that the ss{i,j,}, ss and xline pointers are NOT included
    ! here since they are not AD'ed. 
#endif

  end subroutine setBCPointers2


  subroutine resetBCPointers2(nn)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * resetBCPointers nullifyies the boundary pointers. For reverse  *
    !      * mode AD it copies the values back in to the respective arrays  *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use flowVarRefState
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables
    integer(kind=intType) :: iStart, iEnd, jStart, jEnd

    ! Determine the face id on which the subface is located and set
    ! the pointers accordinly.

    iStart = BCData(nn)%icBeg
    iEnd   = BCData(nn)%icEnd
    jStart = BCData(nn)%jcBeg
    jEnd   = BCData(nn)%jcEnd

    ! Set the size of the subface
    isize = iEnd-iStart + 1
    jsize = jEnd-jStart + 1

#ifndef TAPENADE_REVERSE
    ! This is just a no-opt
#else
    select case (BCFaceID(nn))

       !===============================================================
    case (iMin)
       w(3, iStart:iEnd, jStart:jend, :) = ww3(1:isize,1:jsize,:)
       w(2, iStart:iEnd, jStart:jend, :) = ww2(1:isize,1:jsize,:)
       w(1, iStart:iEnd, jStart:jend, :) = ww1(1:isize,1:jsize,:)
       w(0, iStart:iEnd, jStart:jend, :) = ww0(1:isize,1:jsize,:)

       p(3, iStart:iEnd, jStart:jend) = pp3(1:isize,1:jsize)
       p(2, iStart:iEnd, jStart:jend) = pp2(1:isize,1:jsize)
       p(1, iStart:iEnd, jStart:jend) = pp1(1:isize,1:jsize)
       p(0, iStart:iEnd, jStart:jend) = pp0(1:isize,1:jsize)

       rlv(3, iStart:iEnd, jStart:jend) = rlv3(1:isize,1:jsize)
       rlv(2, iStart:iEnd, jStart:jend) = rlv2(1:isize,1:jsize)
       rlv(1, iStart:iEnd, jStart:jend) = rlv1(1:isize,1:jsize)
       rlv(0, iStart:iEnd, jStart:jend) = rlv0(1:isize,1:jsize)

       rev(3, iStart:iEnd, jStart:jend) = rev3(1:isize,1:jsize)
       rev(2, iStart:iEnd, jStart:jend) = rev2(1:isize,1:jsize)
       rev(1, iStart:iEnd, jStart:jend) = rev1(1:isize,1:jsize)
       rev(0, iStart:iEnd, jStart:jend) = rev0(1:isize,1:jsize)

       gamma(3, iStart:iEnd, jStart:jend) = gamma3(1:isize,1:jsize)
       gamma(2, iStart:iEnd, jStart:jend) = gamma2(1:isize,1:jsize)
       gamma(1, iStart:iEnd, jStart:jend) = gamma1(1:isize,1:jsize) 
       gamma(0, iStart:iEnd, jStart:jend) = gamma0(1:isize,1:jsize)

       !===============================================================

    case (iMax)
       w(nx, iStart:iEnd, jStart:jend, :) = ww3(1:isize,1:jsize,:)
       w(il, iStart:iEnd, jStart:jend, :) = ww2(1:isize,1:jsize,:)
       w(ie, iStart:iEnd, jStart:jend, :) = ww1(1:isize,1:jsize,:)
       w(ib, iStart:iEnd, jStart:jend, :) = ww0(1:isize,1:jsize,:)

       p(nx, iStart:iEnd, jStart:jend) = pp3(1:isize,1:jsize)
       p(il, iStart:iEnd, jStart:jend) = pp2(1:isize,1:jsize)
       p(ie, iStart:iEnd, jStart:jend) = pp1(1:isize,1:jsize)
       p(ib, iStart:iEnd, jStart:jend) = pp0(1:isize,1:jsize)

       rlv(nx, iStart:iEnd, jStart:jend) = rlv3(1:isize,1:jsize)
       rlv(il, iStart:iEnd, jStart:jend) = rlv2(1:isize,1:jsize)
       rlv(ie, iStart:iEnd, jStart:jend) = rlv1(1:isize,1:jsize)
       rlv(ib, iStart:iEnd, jStart:jend) = rlv0(1:isize,1:jsize)

       rev(nx, iStart:iEnd, jStart:jend) = rev3(1:isize,1:jsize)
       rev(il, iStart:iEnd, jStart:jend) = rev2(1:isize,1:jsize)
       rev(ie, iStart:iEnd, jStart:jend) = rev1(1:isize,1:jsize)
       rev(ib, iStart:iEnd, jStart:jend) = rev0(1:isize,1:jsize)

       gamma(nx, iStart:iEnd, jStart:jend) = gamma3(1:isize,1:jsize)
       gamma(il, iStart:iEnd, jStart:jend) = gamma2(1:isize,1:jsize)
       gamma(ie, iStart:iEnd, jStart:jend) = gamma1(1:isize,1:jsize) 
       gamma(ib, iStart:iEnd, jStart:jend) = gamma0(1:isize,1:jsize)

       !===============================================================

    case (jMin)

       w(iStart:iEnd, 3, jStart:jend, :) = ww3(1:isize,1:jsize,:)
       w(iStart:iEnd, 2, jStart:jend, :) = ww2(1:isize,1:jsize,:)
       w(iStart:iEnd, 1, jStart:jend, :) = ww1(1:isize,1:jsize,:)
       w(iStart:iEnd, 0, jStart:jend, :) = ww0(1:isize,1:jsize,:)

       p(iStart:iEnd, 3, jStart:jend) = pp3(1:isize,1:jsize)
       p(iStart:iEnd, 2, jStart:jend) = pp2(1:isize,1:jsize)
       p(iStart:iEnd, 1, jStart:jend) = pp1(1:isize,1:jsize)
       p(iStart:iEnd, 0, jStart:jend) = pp0(1:isize,1:jsize)

       rlv(iStart:iEnd, 3, jStart:jend) = rlv3(1:isize,1:jsize)
       rlv(iStart:iEnd, 2, jStart:jend) = rlv2(1:isize,1:jsize)
       rlv(iStart:iEnd, 1, jStart:jend) = rlv1(1:isize,1:jsize)
       rlv(iStart:iEnd, 0, jStart:jend) = rlv0(1:isize,1:jsize)

       rev(iStart:iEnd, 3, jStart:jend) = rev3(1:isize,1:jsize)
       rev(iStart:iEnd, 2, jStart:jend) = rev2(1:isize,1:jsize)
       rev(iStart:iEnd, 1, jStart:jend) = rev1(1:isize,1:jsize)
       rev(iStart:iEnd, 0, jStart:jend) = rev0(1:isize,1:jsize)

       gamma(iStart:iEnd, 3, jStart:jend) = gamma3(1:isize,1:jsize)
       gamma(iStart:iEnd, 2, jStart:jend) = gamma2(1:isize,1:jsize)
       gamma(iStart:iEnd, 1, jStart:jend) = gamma1(1:isize,1:jsize)
       gamma(iStart:iEnd, 0, jStart:jend) = gamma0(1:isize,1:jsize)

       !===============================================================

    case (jMax)

       w(iStart:iEnd, ny, jStart:jend, :) = ww3(1:isize,1:jsize,:)
       w(iStart:iEnd, jl, jStart:jend, :) = ww2(1:isize,1:jsize,:)
       w(iStart:iEnd, je, jStart:jend, :) = ww1(1:isize,1:jsize,:)
       w(iStart:iEnd, jb, jStart:jend, :) = ww0(1:isize,1:jsize,:)

       p(iStart:iEnd, ny, jStart:jend) = pp3(1:isize,1:jsize)
       p(iStart:iEnd, jl, jStart:jend) = pp2(1:isize,1:jsize)
       p(iStart:iEnd, je, jStart:jend) = pp1(1:isize,1:jsize)
       p(iStart:iEnd, jb, jStart:jend) = pp0(1:isize,1:jsize)

       rlv(iStart:iEnd, ny, jStart:jend) = rlv3(1:isize,1:jsize)
       rlv(iStart:iEnd, jl, jStart:jend) = rlv2(1:isize,1:jsize)
       rlv(iStart:iEnd, je, jStart:jend) = rlv1(1:isize,1:jsize)
       rlv(iStart:iEnd, jb, jStart:jend) = rlv0(1:isize,1:jsize)

       rev(iStart:iEnd, ny, jStart:jend) = rev3(1:isize,1:jsize)
       rev(iStart:iEnd, jl, jStart:jend) = rev2(1:isize,1:jsize)
       rev(iStart:iEnd, je, jStart:jend) = rev1(1:isize,1:jsize)
       rev(iStart:iEnd, jb, jStart:jend) = rev0(1:isize,1:jsize)

       gamma(iStart:iEnd, ny, jStart:jend) = gamma3(1:isize,1:jsize)
       gamma(iStart:iEnd, jl, jStart:jend) = gamma2(1:isize,1:jsize)
       gamma(iStart:iEnd, je, jStart:jend) = gamma1(1:isize,1:jsize)
       gamma(iStart:iEnd, jb, jStart:jend) = gamma0(1:isize,1:jsize)

       !===============================================================

    case (kMin)

       w(iStart:iEnd, jStart:jend, 3, :) = ww3(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, 2, :) = ww2(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, 1, :) = ww1(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, 0, :) = ww0(1:isize,1:jsize,:)

       p(iStart:iEnd, jStart:jend, 3) = pp3(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, 2) = pp2(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, 1) = pp1(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, 0) = pp0(1:isize,1:jsize)

       rlv(iStart:iEnd, jStart:jend, 3) = rlv3(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, 2) = rlv2(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, 1) = rlv1(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, 0) = rlv0(1:isize,1:jsize)

       rev(iStart:iEnd, jStart:jend, 3) = rev3(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, 2) = rev2(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, 1) = rev1(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, 0) = rev0(1:isize,1:jsize)

       gamma(iStart:iEnd, jStart:jend, 3) = gamma3(1:isize,1:jsize)
       gamma(iStart:iEnd, jStart:jend, 2) = gamma2(1:isize,1:jsize)
       gamma(iStart:iEnd, jStart:jend, 1) = gamma1(1:isize,1:jsize)
       gamma(iStart:iEnd, jStart:jend, 0) = gamma0(1:isize,1:jsize)
       
       !===============================================================

    case (kMax)

       w(iStart:iEnd, jStart:jend, nz, :) = ww3(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, kl, :) = ww2(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, ke, :) = ww1(1:isize,1:jsize,:)
       w(iStart:iEnd, jStart:jend, kb, :) = ww0(1:isize,1:jsize,:)

       p(iStart:iEnd, jStart:jend, nz) = pp3(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, kl) = pp2(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, ke) = pp1(1:isize,1:jsize)
       p(iStart:iEnd, jStart:jend, kb) = pp0(1:isize,1:jsize)

       rlv(iStart:iEnd, jStart:jend, nz) = rlv3(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, kl) = rlv2(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, ke) = rlv1(1:isize,1:jsize)
       rlv(iStart:iEnd, jStart:jend, kb) = rlv0(1:isize,1:jsize)

       rev(iStart:iEnd, jStart:jend, nz) = rev3(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, kl) = rev2(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, ke) = rev1(1:isize,1:jsize)
       rev(iStart:iEnd, jStart:jend, kb) = rev0(1:isize,1:jsize)

    end select

    ! Note that the ss{i,j,}, ss and xline pointers are NOT included
    ! here since they are not AD'ed. 
#endif

  end subroutine resetBCPointers2
end module BCRoutines


! This module contains routines used for applying *all* boundary
! conditions for Navier Stokes part of the code. Boundary conditions
! pointers from the BCPointers modules are used. The utilty routines
! setBCPointers (and resetBCPointers) are employed.

module BCRoutines

  implicit none
  save

contains

  subroutine applyAllBC_block(secondHalo)

    ! Apply BC's for a single block
    use constants
    use blockPointers , only : nBocos, BCType, nViscBocos, w, dw, x, vol, il, jl, kl, &
         sectionID, wOld, volOld, BCData, &
       si, sj, sk, sfacei, sfacej, sfacek, rlv, gamma, p, rev, &
       bmtj1, bmtj2, scratch, bmtk2, bmtk1, &
       fw, aa, d2wall, bmti1, bmti2, s
    use utils, only : setBCPointers, resetBCPointers, getCorrectForK
    use BCPointers
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo

    ! Local variables.
    logical :: correctForK
    integer(kind=intType) :: nn
    !
    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    correctForK = getCorrectForK()

    ! Apply all the boundary conditions. The order is important!  Only
    ! some of them have been AD'ed

    ! ------------------------------------
    !  Symmetry Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nBocos
       if (BCType(nn) == symm) then
          call setBCPointers(nn, .False.)
          call bcSymm1stHalo(nn)
          call resetBCPointers(nn, .False.)
       end if
    end do

    if (secondHalo) then
       !$AD II-LOOP
       do nn=1, nBocos
          if (BCType(nn) == symm) then
             call setBCPointers(nn, .False.)
             call bcSymm2ndHalo(nn)
             call resetBCPointers(nn, .False.)
          end if
       end do
    end if

    ! ------------------------------------
    !  Symmetry Polar Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nBocos
       if (BCType(nn) == symmPolar) then
          call setBCPointers(nn, .True.)
          call bcSymmPolar1stHalo(nn)
          call resetBCPointers(nn, .True.)
       end if
    end do

    if (secondHalo) then
       !$AD II-LOOP
       do nn=1, nBocos
          if (BCType(nn) == symmPolar) then
             call setBCPointers(nn, .True.)
             call bcSymmPolar2ndHalo(nn)
             call resetBCPointers(nn, .True.)
          end if
       end do
    end if

    ! ------------------------------------
    !  Adibatic Wall Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nViscBocos
       if (BCType(nn) == NSWallAdiabatic) then
          call setBCPointers(nn, .False.)
          call bcNSWallAdiabatic(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

    ! ------------------------------------
    !  Isotermal Wall Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1, nViscBocos
       if (BCType(nn) == NSWallIsoThermal) then
          call setBCPointers(nn, .False.)
          call bcNSWallIsothermal(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

    ! ------------------------------------
    !  Farfield Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1,nBocos
       if (BCType(nn) == farField) then
          call setBCPointers(nn, .False.)
          call bcFarField(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

    ! ------------------------------------
    !  Subsonic Outflow Boundary Condition
    ! ------------------------------------
    do nn=1,nBocos
       if (BCType(nn) == subSonicOutFlow .or. &
            BCType(nn) == MassBleedOutflow) then
          call setBCPointers(nn, .False.)
          call bcSubSonicOutFlow(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

    ! ------------------------------------
    !  Subsonic Inflow Boundary Condition
    ! ------------------------------------
    do nn=1,nBocos
       if (BCType(nn) == subSonicInFlow) then
          call setBCPointers(nn, .False.)
          call bcSubSonicInflow(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
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
       if (BCType(nn) == extrap .or. &
            BCType(nn) == SupersonicOutFlow) then
          call setBCPointers(nn, .False.)
          call bcExtrap(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

    ! ------------------------------------
    !  Euler Wall Boundary Condition
    ! ------------------------------------
    !$AD II-LOOP
    do nn=1,nBocos
       if (BCType(nn) == EulerWall) then
          call setBCPointers(nn, .True.)
          call bcEulerWall(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .True.)
       end if
    end do

#ifndef USE_TAPENADE
    ! ------------------------------------
    !  Domain Interface Boundary Condition
    ! ------------------------------------
    do nn=1,nBocos
       if (BCType(nn) == DomainInterfaceAll .or. &
            bcTYpe(nn) == DomainInterfaceRhoUVW .or. &
            BCType(nn) == DomainInterfaceP .or. &
            BCType(nn) == DomainInterfaceRho .or. &
            BCType(nn) == DomainInterfaceTotal) then
          call setBCPointers(nn, .False.)
          call bcDomainInterface(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do
#endif
    ! ------------------------------------
    !  Supersonic inflow condition
    ! ------------------------------------
    do nn=1,nBocos
       if (BCType(nn) == SupersonicInflow) then
          call setBCPointers(nn, .False.)
          call bcSupersonicInflow(nn, secondHalo, correctForK)
          call resetBCPointers(nn, .False.)
       end if
    end do

  end subroutine applyAllBC_block

  ! ===================================================================
  !   Actual implementation of each of the boundary condition routines
  ! ===================================================================
  subroutine bcSymm1stHalo(nn)

    !  bcSymm1stHalo applies the symmetry boundary conditions to a
    !  block.  * It is assumed that the pointers in blockPointers are
    !  already set to the correct block on the correct grid level.
    !
    !  In case also the second halo must be set, a second loop is
    !  execulted calling bcSymm2ndhalo. This is the only correct way
    !  in case the block contains only 1 cell between two symmetry
    !  planes, i.e. a 2D problem.

    use constants
    use blockPointers, only : BCdata
    use flowVarRefState, only : viscous, eddyModel
    use BCPointers, only : gamma1, gamma2, ww1, ww2, pp1, pp2, rlv1, rlv2, &
         iStart, jStart, iSize, jSize, rev1, rev2
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType) :: vn, nnx, nny, nnz

    ! Loop over the generic subface to set the state in the
    ! 1-st level halos

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

       ! Set the pressure and gamma and possibly the
       ! laminar and eddy viscosity in the halo.

       gamma1(i,j) = gamma2(i,j)
       pp1(i,j)    = pp2(i,j)
       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)
    enddo
  end subroutine bcSymm1stHalo

  subroutine bcSymm2ndHalo(nn)

    !  bcSymm2ndHalo applies the symmetry boundary conditions to a
    !  block for the 2nd halo. This routine is separate as it makes
    !  AD slightly easier.
    use constants
    use blockPointers, only : BCdata
    use flowVarRefState, only : viscous, eddyModel
    use BCPointers, only : gamma0, gamma3, ww0, ww3, pp0, pp3, rlv0, rlv3, &
         rev0, rev3, iStart, jStart, iSize, jSize
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii
    real(kind=realType) :: vn, nnx, nny, nnz

    ! If we need the second halo, do everything again, but using ww0,
    ! ww3 etc instead of ww2 and ww1.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

       vn = two*(ww3(i,j,ivx)*BCData(nn)%norm(i,j,1) + &
            ww3(i,j,ivy)*BCData(nn)%norm(i,j,2) + &
            ww3(i,j,ivz)*BCData(nn)%norm(i,j,3))

       ! Determine the flow variables in the halo cell.
       ww0(i,j,irho) = ww3(i,j,irho)
       ww0(i,j,ivx) = ww3(i,j,ivx) - vn*BCData(nn)%norm(i,j,1)
       ww0(i,j,ivy) = ww3(i,j,ivy) - vn*BCData(nn)%norm(i,j,2)
       ww0(i,j,ivz) = ww3(i,j,ivz) - vn*BCData(nn)%norm(i,j,3)

       ww0(i,j,irhoE) = ww3(i,j,irhoE)

       ! Set the pressure and gamma and possibly the
       ! laminar and eddy viscosity in the halo.

       gamma0(i,j) = gamma3(i,j)
       pp0(i,j)    = pp3(i,j)
       if( viscous )   rlv0(i,j) = rlv3(i,j)
       if( eddyModel ) rev0(i,j) = rev3(i,j)
    enddo

  end subroutine bcSymm2ndHalo

  subroutine bcSymmPolar1stHalo(nn)

    ! bcSymmPolar applies the polar symmetry boundary conditions to a
    ! singular line of a block. It is assumed that the pointers in
    ! blockPointers are already set to the correct block on the
    ! correct grid level.  The polar symmetry condition is a special
    ! case of a degenerate line, as this line is the axi-symmetric
    ! centerline. This routine does just the 1st level halo.

    use constants
    use BCPointers, only: ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2, &
         xx, iStart, jStart, iSize, jSize
    use flowVarRefState, only : viscous, eddyModel
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii, mm
    real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

    ! Loop over the generic subface to set the state in the
    ! 1-st level halos
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

       ! Determine the unit vector along the degenerated face.
       ! However it is not known which is the singular
       ! direction and therefore determine the direction along
       ! the diagonal (i,j) -- (i-1,j-1), which is correct for
       ! both singular i and j-direction. Note that due to the
       ! usage of the pointer xx there is an offset of +1
       ! in the indices and therefore (i+1,j+1) - (i,j) must
       ! be used to determine this vector.

       nnx = xx(i+1,j+1,1) - xx(i,j,1)
       nny = xx(i+1,j+1,2) - xx(i,j,2)
       nnz = xx(i+1,j+1,3) - xx(i,j,3)

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

       ! Set the pressure and possibly the laminar and
       ! eddy viscosity in the halo.

       pp1(i,j) = pp2(i,j)
       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)
    end do
  end subroutine bcSymmPolar1stHalo

  subroutine bcSymmPolar2ndHalo(nn)

    ! bcSymmPolar applies the polar symmetry boundary conditions to a
    ! singular line of a block. It is assumed that the pointers in
    ! blockPointers are already set to the correct block on the
    ! correct grid level.  The polar symmetry condition is a special
    ! case of a degenerate line, as this line is the axi-symmetric
    ! centerline. This routine does just the 2nd level halo.

    use constants
    use BCPointers, only: ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3, &
         xx, iStart, jStart, iSize, jSize
    use flowVarRefState, only : viscous, eddyModel
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, l, ii, mm
    real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

       ! Determine the unit vector along the degenerated face.
       ! However it is not known which is the singular
       ! direction and therefore determine the direction along
       ! the diagonal (i,j) -- (i-1,j-1), which is correct for
       ! both singular i and j-direction. Note that due to the
       ! usage of the pointer xx there is an offset of +1
       ! in the indices and therefore (i+1,j+1) - (i,j) must
       ! be used to determine this vector.

       nnx = xx(i+1,j+1,1) - xx(i,j,1)
       nny = xx(i+1,j+1,2) - xx(i,j,2)
       nnz = xx(i+1,j+1,3) - xx(i,j,3)

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

       ! Set the pressure and possibly the laminar and
       ! eddy viscosity in the halo.

       pp0(i,j) = pp3(i,j)
       if( viscous )   rlv0(i,j) = rlv3(i,j)
       if( eddyModel ) rev0(i,j) = rev3(i,j)
    enddo

  end subroutine bcSymmPolar2ndHalo

  subroutine bcNSWallAdiabatic(nn, secondHalo, correctForK)

    ! bcNSWallAdiabatic applies the viscous adiabatic wall boundary
    ! condition the pointers already defined.

    use constants
    use blockPointers, only : BCData
    use inputDiscretization , only : viscWallBCTreatment
    use BCPointers, only : ww1, ww2, rlv1, rlv2, pp0, pp1, pp2, pp3, rev1, rev2, &
         iStart, jStart, iSize, jSize
    use flowVarRefState, only : viscous, eddyModel
    use iteration, only : currentLevel, groundLevel
    implicit none

    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn
    integer(kind=intType) :: i, j, ii
    real(kind=realType) :: rhok
    integer(kind=intType) :: wallTreatment

    ! Initialize rhok to zero. This will be overwritten if a
    ! correction for k must be applied.

    rhok = zero

    ! Loop over the generic subface to set the state in the
    ! halo cells.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

       ! Set the viscosities. There is no need to test for a
       ! viscous problem of course. The eddy viscosity is
       ! set to the negative value, as it should be zero on
       ! the wall.

       rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = -rev2(i,j)
  
       ! Make sure that on the coarser grids the constant pressure
       ! boundary condition is used.
       
       wallTreatment = viscWallBcTreatment
       if(currentLevel > groundLevel) wallTreatment = constantPressure
       
       BCTreatment: select case (wallTreatment)
          
       case (constantPressure)
          
          ! Constant pressure. Set the gradient to zero.
          pp1(i, j) = pp2(i, j) - four*third*rhok

       case default
          
          pp1(i,j) = 2*pp2(i,j) - pp3(i,j)
          ! Adjust value if pressure is negative
          if (pp1(i,j) .le. zero) pp1(i,j) = pp2(i,j)
          
       end select BCTreatment
    end do

    ! Compute the energy for these halo's.

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcNSWallAdiabatic

  subroutine bcNSWallIsoThermal(nn, secondHalo, correctForK)

    ! bcNSWallAdiabatic applies the viscous isothermal wall boundary
    ! condition to a block. It is assumed that the BCPointers are
    ! already set

    use constants
    use blockPointers, only : BCData
    use inputDiscretization , only : viscWallBCTreatment
    use BCPointers, only : ww1, ww2, rlv1, rlv2, pp1, pp2, pp3, rev1, rev2, &
         iStart, jStart, iSize, jSize
    use flowVarRefState, only : viscous, eddyModel, RGas
    use iteration, only : currentLevel, groundLevel
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo, correctForK
    integer(kind=intType), intent(in) :: nn

    ! Local variables.
    integer(kind=intType) :: i, j, ii
    integer(kind=intType) :: wallTreatment
    real(kind=realType) :: rhok, t2, t1

    ! Initialize rhok to zero. This will be overwritten if a
    ! correction for k must be applied.

    rhok = zero

    ! Loop over the generic subface to set the state in the
    ! halo cells.

    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

       ! PRESSURE EXTRAPOLATION

       ! Make sure that on the coarser grids the constant pressure
       ! boundary condition is used.

       wallTreatment = viscWallBCTreatment
       if(currentLevel > groundLevel) wallTreatment = constantPressure

       BCTreatment: select case (wallTreatment)

       case (constantPressure)

          ! Constant pressure. Set the gradient to zero.
          pp1(i,j) = pp2(i,j) - four*third*rhok

       case default

          ! Linear extrapolation.
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

          pp1(i,j) = 2*pp2(i,j) - pp3(i,j)
          ! Adjust value if pressure is negative
          if (pp1(i,j) .le. zero) pp1(i,j) = pp2(i,j)

       end select BCTreatment

       ! Determine the variables in the halo. As the spacing
       ! is very small a constant pressure boundary condition
       ! (except for the k correction) is okay. Take the slip
       ! velocity into account.

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

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcNSWallIsoThermal

  subroutine bcSubsonicOutflow(nn, secondHalo, correctForK)

    !  bcSubsonicOutflow applies the subsonic outflow boundary
    !  condition, static pressure prescribed, to a block. It is
    !  assumed that the pointers in blockPointers are already set to
    !  the correct block on the correct grid level.  Exactly the same
    !  boundary condition is also applied for an outflow mass
    !  bleed. Therefore the test is for both a subsonic outflow and an
    !  bleed outflow.

    use constants
    use blockPointers, only : BCData
    use BCPointers, only : gamma2, rev2, rlv2, pp2, ww2, &
         rlv1, rev1, pp1, ww1, iSize, jSize, iStart, jStart
    use flowVarRefState, only : eddyModel, viscous
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
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcSubsonicOutflow

  subroutine bcSubsonicInflow(nn, secondHalo, correctForK)

    !  bcSubsonicInflow applies the subsonic outflow boundary
    !  condition, total pressure, total density and flow direction
    !  prescribed, to a block. It is assumed that the pointers in
    !  blockPointers are already set to the correct block on the
    !  correct grid level.

    use constants
    use blockPointers, only : BCData
    use flowVarRefState, only : viscous, eddyModel, RGas
    use inputDiscretization, only : hScalingInlet
    use BCPointers, only : gamma2, ww2, pp2, rlv1, rlv2, rev1, rev2, &
         ww1, pp1, iSize, jSize, iStart, jStart
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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

    end select

    ! Compute the energy for these halo's.

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcSubsonicInflow

  subroutine bcEulerWall(nn, secondHalo, correctForK)

    !  bcEulerWall applies the inviscid wall boundary condition to a
    !  block. It is assumed that the bcpointers are already set to the
    !  correct block on the correct grid level.

    use constants
    use blockPointers, only : BCData, addGridVelocities
    use flowVarRefState, only : viscous, eddyModel, RGas
    use inputDiscretization, only : eulerWallBCTreatment
    use iteration, only : currentLevel, groundLevel
    use utils, only : myDim
    use BCPointers, only : ww1, pp1, rlv1, rev1, ww2, pp2, rlv2, rev2, &
         pp3, ss, ssi, ssj, ssk, iStart, iSize, jStart, jSize, iEnd, jEnd
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
    real(kind=realType), dimension(iStart:iEnd, jStart:jEnd) :: grad

    ! Make sure that on the coarser grids the constant pressure
    ! boundary condition is used.

    wallTreatment = eulerWallBcTreatment
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
          j = mod(ii, isize) + iStart
          k = ii/isize + jStart
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
          j = mod(ii, isize) + iStart
          k = ii/isize + jStart

          ! Store the indices k+1, k-1 a bit easier and make
          ! sure that they do not exceed the range of the arrays.

          km1 = k-1; km1 = max(jStart, km1)
          kp1 = k+1; kp1 = min(jend  , kp1)

          ! Compute the scaling factor for the central difference
          ! in the k-direction.

          b1 = one/max(1_intType,(kp1-km1))

          ! The indices j+1 and j-1. Make sure that they
          ! do not exceed the range of the arrays.

          jm1 = j-1; jm1 = max(iStart, jm1)
          jp1 = j+1; jp1 = min(iEnd  , jp1)

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
       j = mod(ii, isize) + iStart
       k = ii/isize + jStart

       ! Compute the pressure density and velocity in the
       ! halo cell. Note that rface is the grid velocity
       ! component in the direction of norm, i.e. outward
       ! pointing.

       pp1(j,k) = mydim(pp2(j,k), grad(j,k))

       vn = two*(BCData(nn)%rface(j,k) - &
            ww2(j,k,ivx)*BCData(nn)%norm(j,k,1) - &
            ww2(j,k,ivy)*BCData(nn)%norm(j,k,2) - &
            ww2(j,k,ivz)*BCData(nn)%norm(j,k,3))

       ww1(j,k,irho) = ww2(j,k,irho)
       ww1(j,k,ivx)  = ww2(j,k,ivx) + vn*BCData(nn)%norm(j,k,1)
       ww1(j,k,ivy)  = ww2(j,k,ivy) + vn*BCData(nn)%norm(j,k,2)
       ww1(j,k,ivz)  = ww2(j,k,ivz) + vn*BCData(nn)%norm(j,k,3)

       ! The laminar and eddy viscosity, if present.

       if( viscous )    rlv1(j,k) = rlv2(j,k)
       if( eddyModel ) rev1(j,k) = rev2(j,k)

    enddo

    ! Compute the energy for these halo's.

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcEulerWall

#ifndef USE_TAPENADE
  subroutine bcDomainInterface(nn, secondHalo, correctForK)

    !  bcDomainInterface applies the domain-interface boundary
    !  condition, where necessary flow variables are obtained from the
    !  coupler. More options can be added in the future.

    use constants
    use blockPointers, only : BCData, BCType
    use flowVarRefState, only : viscous, eddyModel, PinfCorr
    use inputDiscretization, only : HscalingInlet
    use BCPointers, only : ww1, pp1, rlv1, rev1, ww2, gamma2, rlv2, rev2, &
         pp2, ww2, iSize, jSize, iStart, jStart
    use flowVarRefState, only : viscous, eddyModel, RGas
    use utils, only : terminate
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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

          ww1(i,j,irho) = BCData(nn)%rho(i,j)
          ww1(i,j,ivx)  = BCData(nn)%velx(i,j)
          ww1(i,j,ivy)  = BCData(nn)%vely(i,j)
          ww1(i,j,ivz)  = BCData(nn)%velz(i,j)
          pp1(i,j)      = BCData(nn)%ps(i,j)

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the energy for these halo's.

       call computeEtot(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo(correctForK)

       !============================================================

    case (DomainInterfaceRhoUVW)

       ! Flow variables are already prescribed from the coupler
       ! (see the subroutine setInterfaceData).

       ! Loop over the generic subface and copy the density,
       ! velocities and turbulent variables. The pressure is
       ! computed using Riemann invariants.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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

          ! Set the viscosities in the halo to the viscosities
          ! in the donor cell.

          if( viscous )   rlv1(i,j) = rlv2(i,j)
          if( eddyModel ) rev1(i,j) = rev2(i,j)

       enddo

       ! Compute the energy for these halo's.

       call computeEtot(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo(correctForK)

       !============================================================

    case (DomainInterfaceP)

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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

       call computeEtot(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo(correctForK)

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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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

       call computeEtot(ww1, pp1, correctForK)

       ! Extrapolate the state vectors in case a second halo
       ! is needed.

       if( secondHalo ) call extrapolate2ndHalo(correctForK)

    end select
  end subroutine bcDomainInterface
#endif

  subroutine bcFarfield(nn, secondHalo, correctForK)

    ! bcFarfield applies the farfield boundary condition to a block.
    ! It is assumed that the BCPointers are already set *

    use constants
    use blockPointers, only : BCData
    use flowVarRefState, only : eddyModel, viscous, gammaInf, wInf, pInfCorr
    use BCPointers, only : ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, &
         rev0, rev1, rev2, gamma2, iStart, jStart, iSize, jSize
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
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

       else
          ! Inflow
          uf = u0 + (qnf - qn0)*BCData(nn)%norm(i,j,1)
          vf = v0 + (qnf - qn0)*BCData(nn)%norm(i,j,2)
          wf = w0 + (qnf - qn0)*BCData(nn)%norm(i,j,3)
          sf = s0

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
    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.
    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcFarfield

  subroutine bcSupersonicInflow(nn, secondHalo, correctForK)

    ! bcSupersonicInflow applies the supersonic inflow boundary
    ! conditions, entire state vector is prescribed, to a block. It is
    ! assumed that the pointers in blockPointers are already set to
    ! the correct block on the correct grid level.

    use constants
    use blockPointers, only : BCData
    use flowVarRefState, only : eddyModel, viscous
    use BCPointers, only : ww0, ww1, pp0, pp1, rlv0, rlv1, rlv2, &
         rev0, rev1, rev2, iStart, jStart, iSize, jSize
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
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

       ww1(i,j,irho) = BCData(nn)%rho(i,j)
       ww1(i,j,ivx)  = BCData(nn)%velx(i,j)
       ww1(i,j,ivy)  = BCData(nn)%vely(i,j)
       ww1(i,j,ivz)  = BCData(nn)%velz(i,j)
       pp1(i,j)      = BCData(nn)%ps(i,j)

       ! Set the laminar and eddy viscosity in the halo
       ! if needed.

       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)
    end do

    call computeEtot(ww1, pp1, correctForK)

    if (secondHalo) then
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

          ww0(i,j,irho) = BCData(nn)%rho(i,j)
          ww0(i,j,ivx)  = BCData(nn)%velx(i,j)
          ww0(i,j,ivy)  = BCData(nn)%vely(i,j)
          ww0(i,j,ivz)  = BCData(nn)%velz(i,j)
          pp0(i,j)      = BCData(nn)%ps(i,j)

          ! Set the laminar and eddy viscosity in the halo
          ! if needed.

          if( viscous )   rlv0(i,j) = rlv1(i,j)
          if( eddyModel ) rev0(i,j) = rev1(i,j)
       end do

       call computeEtot(ww0, pp0, correctForK)
    end if

  end subroutine bcSupersonicInflow

  subroutine bcExtrap(nn, secondHalo, correctForK)
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
    use constants
    use blockPointers, only : BCType
    use flowVarRefState, only : viscous, eddyModel
    use inputDiscretization, onlY : outflowTreatment
    !use inputPhysics
    use BCPointers, only : ww1, ww2, ww3, pp1, pp2, pp3, &
         rlv1, rlv2, rev1, rev2, iStart, jStart, iSize, jSize
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
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

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

       ! The laminar and eddy viscosity, if present. These
       ! values are simply taken constant. Their values do
       ! not really matter.

       if( viscous )   rlv1(i,j) = rlv2(i,j)
       if( eddyModel ) rev1(i,j) = rev2(i,j)

    enddo

    ! Compute the energy for these halo's.

    call computeEtot(ww1, pp1, correctForK)

    ! Extrapolate the state vectors in case a second halo
    ! is needed.

    if( secondHalo ) call extrapolate2ndHalo(correctForK)

  end subroutine bcExtrap

  subroutine pRhoSubsonicInlet(ww, pp, correctForK)

    !  pRhoSubsonicInlet computes the pressure and density for the
    !  given range of the block to which the pointers in blockPointers
    !  currently point.

    use constants
    use cpCurveFits
    use flowVarRefState, only : RGas, Tref
    use inputPhysics, only : cpModel, gammaConstant
    use BCPointers, only : iSize, jSize, iStart, jStart
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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart

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

#ifndef USE_TAPENADE
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
#endif
    end select

    ! Add 2*rho*k/3 to the pressure if a k-equation is present.

    if( correctForK ) then
       !$AD II-LOOP
       do ii=0,isize*jsize-1
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart
          pp(i,j) = pp(i,j) &
               + twoThird*ww(i,j,irho)*ww(i,j,itu1)
       enddo
    end if
#ifndef USE_TAPENADAE
  contains
    subroutine cportIntegrant(T, nn, int)

      ! cportIntegrant computes the integrant of the function cp/(r*t)
      ! for the given temperature. It also stores the correct curve
      ! fit interval, which is needed to determine the entire integral
      ! in the main subroutine.

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
#endif
  end subroutine pRhoSubsonicInlet

  subroutine computeEtot(ww, pp, correctForK)

    ! Simplified total energy computation for boundary conditions.
    ! Only implements the constant cpModel

    use constants
    use inputPhysics, only : gammaConstant, cpModel
    use utils, only :terminate
    use BCPointers, only : iSize, jSize, iStart, jStart
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
          i = mod(ii, isize) + iStart
          j = ii/isize + jStart
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
  end subroutine computeEtot

  subroutine extrapolate2ndHalo(correctForK)

    ! extrapolate2ndHalo determines the states of the second layer
    ! halo cells for the given subface of the block. It is assumed
    ! that the appropriate BCPointers are already set

    use constants
    use BCPointers, only : ww0, ww1, ww2, pp0, pp1, pp2, &
         rlv0, rlv1, rlv2, rev0, rev1, rev2, iSize, jSize, iStart, jStart
    use flowVarRefState, only : viscous, eddyModel
    implicit none

    ! Input variables
    logical, intent(in) :: correctForK

    ! Working variables
    real(kind=realType), parameter :: factor = 0.5_realType
    integer(kind=intType) :: i, j, l, ii

    ! Loop over the generic subface to set the state in the
    ! halo cells.
    !$AD II-LOOP
    do ii=0,isize*jsize-1
       i = mod(ii, isize) + iStart
       j = ii/isize + jStart

       ! Extrapolate the density, momentum and pressure.
       ! Make sure that a certain threshold is kept.

       ww0(i,j,irho) = two*ww1(i,j,irho) - ww2(i,j,irho)
       ww0(i,j,irho) = max(factor*ww1(i,j,irho),ww0(i,j,irho))

       ww0(i,j,ivx) = two*ww1(i,j,ivx) - ww2(i,j,ivx)
       ww0(i,j,ivy) = two*ww1(i,j,ivy) - ww2(i,j,ivy)
       ww0(i,j,ivz) = two*ww1(i,j,ivz) - ww2(i,j,ivz)

       pp0(i,j) = max(factor*pp1(i,j),two*pp1(i,j) - pp2(i,j))

       ! The laminar and eddy viscosity, if present. These values
       ! are simply taken constant. Their values do not matter.

       if( viscous )   rlv0(i,j) = rlv1(i,j)
       if( eddyModel ) rev0(i,j) = rev1(i,j)
    enddo

    ! Compute the energy for this halo range.
    call computeEtot(ww0, pp0, correctForK)

  end subroutine extrapolate2ndHalo

end module BCRoutines

! @File    :   cudaBCRoutines.f90
! @Desc    :   This module contains the CUDA Fortran version of the boundary condition routines in BCRoutines.F90

module cudaBCRoutines

    ! --- CUDA FORTRAN module ---
    use cudafor

    implicit none
    save

contains

    subroutine copy_data
        ! Does the same thing as subroutine copycudaBlock in cudaBlock.f90 but does NOT store in cudaBlockType
        ! Eventually, this should get thrown out 
        use blockPointers, only: blckptrBCData => BCData

        implicit none
        
        h_BCData = blckptrBCData
        BCData = blckptrBCData

    end subroutine copy_data

! #ifndef USE_TAPENADE
    attributes(global) subroutine applyAllBC(secondHalo)
        !
        !       applyAllBC applies all boundary conditions for the all
        !       blocks on the grid level currentLevel.
        !
        use constants
        use cudaBlock, only: nDom
        use iteration, only: currentLevel
        use utils, only: setPointers
        use ALEUtils, only: interpLevelALEBC_Block, recoverLevelALEBC_block
        implicit none
        !
        !      Subroutine arguments.
        !
        logical, intent(in) :: secondHalo
        !
        !      Local Variables
        integer(kind=intType) :: sps, nn
        
        ! manually skip time spectral for now
        integer(kind=intType) :: nTimeIntervalsSpectral = 1

        sps = (blockIdx%y-1)*blockDim%y + threadIdx%y
        nn = (blockIdx%z-1)*blockDim%z + threadIdx%z

        ! Loop over the number of spectral solutions and number of blocks.
        
        if (sps <= spse .AND. nn <= nne) then

                ! Set the pointers for this block.

                call setPointers(nn, currentLevel, sps)

                call interpLevelALEBC_block
                call applyAllBC_block(secondHalo)
                call recoverLevelALEBC_block

            end do domains
        end do spectralLoop

    end subroutine applyAllBC
! #endif

    ! ===================================================================
    !   Actual implementation of each of the boundary condition routines
    ! ===================================================================
    attributes(global) subroutine bcSymm1stHalo(nn)

        !  bcSymm1stHalo applies the symmetry boundary conditions to a
        !  block.  * It is assumed that the pointers in blockPointers are
        !  already set to the correct block on the correct grid level.
        !
        !  In case also the second halo must be set, a second loop is
        !  execulted calling bcSymm2ndhalo. This is the only correct way
        !  in case the block contains only 1 cell between two symmetry
        !  planes, i.e. a 2D problem.

        use constants
        use blockPointers, only: BCdata
        use cudaflowVarRefState, only: viscous, eddyModel
        use cudaBCPointers, only: gamma1, gamma2, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              iStart, jStart, iSize, jSize, rev1, rev2
        implicit none

        ! Subroutine arguments.
        integer(kind=intType), intent(in) :: nn

        ! Local variables.
        integer(kind=intType) :: i, j, l, ii
        real(kind=realType) :: vn, nnx, nny, nnz

        ! Loop over the generic subface to set the state in the
        ! 1-st level halos

        ! --- Initialize GPU thread indices Launch params ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then

            ! Determine twice the normal velocity component,
            ! which must be substracted from the donor velocity
            ! to obtain the halo velocity.

            vn = two * (ww2(i, j, ivx) * BCData(nn)%norm(i, j, 1) + &
                        ww2(i, j, ivy) * BCData(nn)%norm(i, j, 2) + &
                        ww2(i, j, ivz) * BCData(nn)%norm(i, j, 3))

            ! Determine the flow variables in the halo cell.

            ww1(i, j, irho) = ww2(i, j, irho)
            ww1(i, j, ivx) = ww2(i, j, ivx) - vn * BCData(nn)%norm(i, j, 1)
            ww1(i, j, ivy) = ww2(i, j, ivy) - vn * BCData(nn)%norm(i, j, 2)
            ww1(i, j, ivz) = ww2(i, j, ivz) - vn * BCData(nn)%norm(i, j, 3)
            ww1(i, j, irhoE) = ww2(i, j, irhoE)

            ! Set the pressure and gamma and possibly the
            ! laminar and eddy viscosity in the halo.

            gamma1(i, j) = gamma2(i, j)
            pp1(i, j) = pp2(i, j)
            if (viscous) rlv1(i, j) = rlv2(i, j)
            if (eddyModel) rev1(i, j) = rev2(i, j)
        end if
    end subroutine bcSymm1stHalo

    ! subroutine bcSymm2ndHalo(nn)

    !     !  bcSymm2ndHalo applies the symmetry boundary conditions to a
    !     !  block for the 2nd halo. This routine is separate as it makes
    !     !  AD slightly easier.
    !     use constants
    !     use blockPointers, only: BCdata
    !     use flowVarRefState, only: viscous, eddyModel
    !     use BCPointers, only: gamma0, gamma3, ww0, ww3, pp0, pp3, rlv0, rlv3, &
    !                           rev0, rev3, iStart, jStart, iSize, jSize
    !     implicit none

    !     ! Subroutine arguments.
    !     integer(kind=intType), intent(in) :: nn

    !     ! Local variables.
    !     integer(kind=intType) :: i, j, l, ii
    !     real(kind=realType) :: vn, nnx, nny, nnz

    !     ! If we need the second halo, do everything again, but using ww0,
    !     ! ww3 etc instead of ww2 and ww1.

    !     !$AD II-LOOP
    !     do ii = 0, isize * jsize - 1
    !         i = mod(ii, isize) + iStart
    !         j = ii / isize + jStart

    !         vn = two * (ww3(i, j, ivx) * BCData(nn)%norm(i, j, 1) + &
    !                     ww3(i, j, ivy) * BCData(nn)%norm(i, j, 2) + &
    !                     ww3(i, j, ivz) * BCData(nn)%norm(i, j, 3))

    !         ! Determine the flow variables in the halo cell.
    !         ww0(i, j, irho) = ww3(i, j, irho)
    !         ww0(i, j, ivx) = ww3(i, j, ivx) - vn * BCData(nn)%norm(i, j, 1)
    !         ww0(i, j, ivy) = ww3(i, j, ivy) - vn * BCData(nn)%norm(i, j, 2)
    !         ww0(i, j, ivz) = ww3(i, j, ivz) - vn * BCData(nn)%norm(i, j, 3)

    !         ww0(i, j, irhoE) = ww3(i, j, irhoE)

    !         ! Set the pressure and gamma and possibly the
    !         ! laminar and eddy viscosity in the halo.

    !         gamma0(i, j) = gamma3(i, j)
    !         pp0(i, j) = pp3(i, j)
    !         if (viscous) rlv0(i, j) = rlv3(i, j)
    !         if (eddyModel) rev0(i, j) = rev3(i, j)
    !     end do

    ! end subroutine bcSymm2ndHalo

    ! subroutine bcSymmPolar1stHalo(nn)

    !     ! bcSymmPolar applies the polar symmetry boundary conditions to a
    !     ! singular line of a block. It is assumed that the pointers in
    !     ! blockPointers are already set to the correct block on the
    !     ! correct grid level.  The polar symmetry condition is a special
    !     ! case of a degenerate line, as this line is the axi-symmetric
    !     ! centerline. This routine does just the 1st level halo.

    !     use constants
    !     use BCPointers, only: ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2, &
    !                           xx, iStart, jStart, iSize, jSize
    !     use flowVarRefState, only: viscous, eddyModel
    !     implicit none

    !     ! Subroutine arguments.
    !     integer(kind=intType), intent(in) :: nn

    !     ! Local variables.
    !     integer(kind=intType) :: i, j, l, ii, mm
    !     real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

    !     ! Loop over the generic subface to set the state in the
    !     ! 1-st level halos
    !     !$AD II-LOOP
    !     do ii = 0, isize * jsize - 1
    !         i = mod(ii, isize) + iStart
    !         j = ii / isize + jStart

    !         ! Determine the unit vector along the degenerated face.
    !         ! However it is not known which is the singular
    !         ! direction and therefore determine the direction along
    !         ! the diagonal (i,j) -- (i-1,j-1), which is correct for
    !         ! both singular i and j-direction. Note that due to the
    !         ! usage of the pointer xx there is an offset of +1
    !         ! in the indices and therefore (i+1,j+1) - (i,j) must
    !         ! be used to determine this vector.

    !         nnx = xx(i + 1, j + 1, 1) - xx(i, j, 1)
    !         nny = xx(i + 1, j + 1, 2) - xx(i, j, 2)
    !         nnz = xx(i + 1, j + 1, 3) - xx(i, j, 3)

    !         ! Determine the unit vector in this direction.

    !         tmp = one / sqrt(nnx * nnx + nny * nny + nnz * nnz)
    !         nnx = nnx * tmp
    !         nny = nny * tmp
    !         nnz = nnz * tmp

    !         ! Determine twice the tangential velocity vector of the
    !         ! internal cell.

    !         tmp = two * (ww2(i, j, ivx) * nnx + ww2(i, j, ivy) * nny &
    !                      + ww2(i, j, ivz) * nnz)
    !         vtx = tmp * nnx
    !         vty = tmp * nny
    !         vtz = tmp * nnz

    !         ! Determine the flow variables in the halo cell. The
    !         ! velocity is constructed such that the average of the
    !         ! internal and the halo cell is along the centerline.
    !         ! Note that the magnitude of the velocity does not
    !         ! change and thus the energy is identical.

    !         ww1(i, j, irho) = ww2(i, j, irho)
    !         ww1(i, j, ivx) = vtx - ww2(i, j, ivx)
    !         ww1(i, j, ivy) = vty - ww2(i, j, ivy)
    !         ww1(i, j, ivz) = vtz - ww2(i, j, ivz)
    !         ww1(i, j, irhoE) = ww2(i, j, irhoE)

    !         ! Set the pressure and possibly the laminar and
    !         ! eddy viscosity in the halo.

    !         pp1(i, j) = pp2(i, j)
    !         if (viscous) rlv1(i, j) = rlv2(i, j)
    !         if (eddyModel) rev1(i, j) = rev2(i, j)
    !     end do
    ! end subroutine bcSymmPolar1stHalo

    ! subroutine bcSymmPolar2ndHalo(nn)

    !     ! bcSymmPolar applies the polar symmetry boundary conditions to a
    !     ! singular line of a block. It is assumed that the pointers in
    !     ! blockPointers are already set to the correct block on the
    !     ! correct grid level.  The polar symmetry condition is a special
    !     ! case of a degenerate line, as this line is the axi-symmetric
    !     ! centerline. This routine does just the 2nd level halo.

    !     use constants
    !     use BCPointers, only: ww0, ww3, pp0, pp3, rlv0, rlv3, rev0, rev3, &
    !                           xx, iStart, jStart, iSize, jSize
    !     use flowVarRefState, only: viscous, eddyModel
    !     implicit none

    !     ! Subroutine arguments.
    !     integer(kind=intType), intent(in) :: nn

    !     ! Local variables.
    !     integer(kind=intType) :: i, j, l, ii, mm
    !     real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

    !     !$AD II-LOOP
    !     do ii = 0, isize * jsize - 1
    !         i = mod(ii, isize) + iStart
    !         j = ii / isize + jStart

    !         ! Determine the unit vector along the degenerated face.
    !         ! However it is not known which is the singular
    !         ! direction and therefore determine the direction along
    !         ! the diagonal (i,j) -- (i-1,j-1), which is correct for
    !         ! both singular i and j-direction. Note that due to the
    !         ! usage of the pointer xx there is an offset of +1
    !         ! in the indices and therefore (i+1,j+1) - (i,j) must
    !         ! be used to determine this vector.

    !         nnx = xx(i + 1, j + 1, 1) - xx(i, j, 1)
    !         nny = xx(i + 1, j + 1, 2) - xx(i, j, 2)
    !         nnz = xx(i + 1, j + 1, 3) - xx(i, j, 3)

    !         ! Determine the unit vector in this direction.

    !         tmp = one / sqrt(nnx * nnx + nny * nny + nnz * nnz)
    !         nnx = nnx * tmp
    !         nny = nny * tmp
    !         nnz = nnz * tmp

    !         ! Determine twice the tangential velocity vector of the
    !         ! internal cell.

    !         tmp = two * (ww3(i, j, ivx) * nnx + ww3(i, j, ivy) * nny &
    !                      + ww3(i, j, ivz) * nnz)
    !         vtx = tmp * nnx
    !         vty = tmp * nny
    !         vtz = tmp * nnz

    !         ! Determine the flow variables in the halo cell. The
    !         ! velocity is constructed such that the average of the
    !         ! internal and the halo cell is along the centerline.
    !         ! Note that the magnitude of the velocity does not
    !         ! change and thus the energy is identical.

    !         ww0(i, j, irho) = ww3(i, j, irho)
    !         ww0(i, j, ivx) = vtx - ww3(i, j, ivx)
    !         ww0(i, j, ivy) = vty - ww3(i, j, ivy)
    !         ww0(i, j, ivz) = vtz - ww3(i, j, ivz)
    !         ww0(i, j, irhoE) = ww3(i, j, irhoE)

    !         ! Set the pressure and possibly the laminar and
    !         ! eddy viscosity in the halo.

    !         pp0(i, j) = pp3(i, j)
    !         if (viscous) rlv0(i, j) = rlv3(i, j)
    !         if (eddyModel) rev0(i, j) = rev3(i, j)
    !     end do

    ! end subroutine bcSymmPolar2ndHalo

    attributes(global) subroutine bcNSWallAdiabatic(nn, secondHalo, correctForK)

        ! bcNSWallAdiabatic applies the viscous adiabatic wall boundary
        ! condition the pointers already defined.

        use constants
        ! use blockPointers, only: BCData
        use cudaBlock, only: d_cudaDoms
        use inputDiscretization, only: viscWallBCTreatment
        use cudaBCPointers, only: ww0, ww1, ww2, rlv0, rlv1, rlv2, pp0, pp1, pp2, pp3, rev0, &
                              rev1, rev2, iStart, jStart, iSize, jSize
        use cudaflowVarRefState, only: viscous, eddyModel
        use iteration, only: currentLevel, groundLevel
        implicit none

        logical, intent(in) :: secondHalo, correctForK
        integer(kind=intType), intent(in) :: nn
        integer(kind=intType) :: sps = 1 ! hard-coding sps = 1 because of how the subroutine signature works
        real(kind=realType) :: rhok
        integer(kind=intType) :: wallTreatment
        integer(kind=intType) :: i, j ! GPU thread indices

        ! Initialize rhok to zero. This will be overwritten if a
        ! correction for k must be applied.

        rhok = zero

        ! Loop over the generic subface to set the state in the
        ! halo cells.

        ! --- Initialize GPU thread indices (assuming 2d launch params) ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then

            ! Set the value of rhok if a correcton must be applied.
            ! It probably does not matter too much, because k is very
            ! small near the wall.

            if (correctForK) rhok = ww2(i, j, irho) * ww2(i, j, itu1)

            ! Determine the variables in the halo. As the spacing
            ! is very small a constant pressure boundary condition
            ! (except for the k correction) is okay. Take the slip
            ! velocity into account.

            ww1(i, j, irho) = ww2(i, j, irho)
            ww1(i, j, ivx) = -ww2(i, j, ivx) + two * d_cudaDoms(nn,sps)%BCData(nn)%uSlip(i, j, 1)
            ww1(i, j, ivy) = -ww2(i, j, ivy) + two * d_cudaDoms(nn,sps)%BCData(nn)%uSlip(i, j, 2)
            ww1(i, j, ivz) = -ww2(i, j, ivz) + two * d_cudaDoms(nn,sps)%BCData(nn)%uSlip(i, j, 3)

            ! Set the viscosities. There is no need to test for a
            ! viscous problem of course. The eddy viscosity is
            ! set to the negative value, as it should be zero on
            ! the wall.

            rlv1(i, j) = rlv2(i, j)
            if (eddyModel) rev1(i, j) = -rev2(i, j)

            ! Make sure that on the coarser grids the constant pressure
            ! boundary condition is used.

            wallTreatment = viscWallBcTreatment
            if (currentLevel > groundLevel) wallTreatment = constantPressure

            BCTreatment:select case(wallTreatment)

            case (constantPressure)

            ! Constant pressure. Set the gradient to zero.
            pp1(i, j) = pp2(i, j) - four * third * rhok

            case default

            pp1(i, j) = 2 * pp2(i, j) - pp3(i, j)
            ! Adjust value if pressure is negative
            if (pp1(i, j) .le. zero) pp1(i, j) = pp2(i, j)

            end select BCTreatment
        end if

        ! Compute the energy for these halo's.

        call computeEtot(ww1, pp1, correctForK)

        ! Extrapolate the state vectors in case a second halo
        ! is needed.

        if (secondHalo) call extrapolate2ndHalo(correctForK)

    end subroutine bcNSWallAdiabatic

    ! subroutine bcNSWallIsoThermal(nn, secondHalo, correctForK)

    !     ! bcNSWallAdiabatic applies the viscous isothermal wall boundary
    !     ! condition to a block. It is assumed that the BCPointers are
    !     ! already set

    !     use constants
    !     use blockPointers, only: BCData
    !     use inputDiscretization, only: viscWallBCTreatment
    !     use BCPointers, only: ww0, ww1, ww2, rlv0, rlv1, rlv2, pp0, pp1, pp2, pp3, &
    !                           rev0, rev1, rev2, iStart, jStart, iSize, jSize
    !     use flowVarRefState, only: viscous, eddyModel, RGas
    !     use iteration, only: currentLevel, groundLevel
    !     implicit none

    !     ! Subroutine arguments.
    !     logical, intent(in) :: secondHalo, correctForK
    !     integer(kind=intType), intent(in) :: nn

    !     ! Local variables.
    !     integer(kind=intType) :: i, j, ii
    !     integer(kind=intType) :: wallTreatment
    !     real(kind=realType) :: rhok, t2, t1

    !     ! Initialize rhok to zero. This will be overwritten if a
    !     ! correction for k must be applied.

    !     rhok = zero

    !     ! Loop over the generic subface to set the state in the
    !     ! halo cells.

    !     !$AD II-LOOP
    !     do ii = 0, isize * jsize - 1
    !         i = mod(ii, isize) + iStart
    !         j = ii / isize + jStart

    !         ! Set the value of rhok if a correcton must be applied.
    !         ! It probably does not matter too much, because k is very
    !         ! small near the wall.

    !         if (correctForK) rhok = ww2(i, j, irho) * ww2(i, j, itu1)

    !         ! Compute the temperature in the internal cell and in the
    !         ! halo cell such that the average is the wall temperature.

    !         t2 = pp2(i, j) / (RGas * ww2(i, j, irho))
    !         t1 = two * bcData(nn)%TNS_Wall(i, j) - t2

    !         ! Make sure that t1 is within reasonable bounds. These
    !         ! bounds are such that the clipping is never active in the
    !         ! converged solution; it is only to avoid instabilities
    !         ! during the convergence.

    !         t1 = max(half * bcData(nn)%TNS_Wall(i, j), t1)
    !         t1 = min(two * bcData(nn)%TNS_Wall(i, j), t1)

    !         ! PRESSURE EXTRAPOLATION

    !         ! Make sure that on the coarser grids the constant pressure
    !         ! boundary condition is used.

    !         wallTreatment = viscWallBCTreatment
    !         if (currentLevel > groundLevel) wallTreatment = constantPressure

    !         BCTreatment:select case(wallTreatment)

    !         case (constantPressure)

    !         ! Constant pressure. Set the gradient to zero.
    !         pp1(i, j) = pp2(i, j) - four * third * rhok

    !         case default

    !         ! Linear extrapolation.
    !         i = mod(ii, isize) + iStart
    !         j = ii / isize + jStart

    !         pp1(i, j) = 2 * pp2(i, j) - pp3(i, j)
    !         ! Adjust value if pressure is negative
    !         if (pp1(i, j) .le. zero) pp1(i, j) = pp2(i, j)

    !         end select BCTreatment

    !         ! Determine the variables in the halo. As the spacing
    !         ! is very small a constant pressure boundary condition
    !         ! (except for the k correction) is okay. Take the slip
    !         ! velocity into account.

    !         ww1(i, j, irho) = pp1(i, j) / (RGas * t1)
    !         ww1(i, j, ivx) = -ww2(i, j, ivx) + two * bcData(nn)%uSlip(i, j, 1)
    !         ww1(i, j, ivy) = -ww2(i, j, ivy) + two * bcData(nn)%uSlip(i, j, 2)
    !         ww1(i, j, ivz) = -ww2(i, j, ivz) + two * bcData(nn)%uSlip(i, j, 3)

    !         ! Set the viscosities. There is no need to test for a
    !         ! viscous problem of course. The eddy viscosity is
    !         ! set to the negative value, as it should be zero on
    !         ! the wall.

    !         rlv1(i, j) = rlv2(i, j)
    !         if (eddyModel) rev1(i, j) = -rev2(i, j)
    !     end do

    !     ! Compute the energy for these halo's.

    !     call computeEtot(ww1, pp1, correctForK)

    !     ! Extrapolate the state vectors in case a second halo
    !     ! is needed.

    !     if (secondHalo) call extrapolate2ndHalo(correctForK)

    ! end subroutine bcNSWallIsoThermal

    attributes(global) subroutine bcSubsonicOutflow(nn, secondHalo, correctForK)
        ! read the subroutine description in BCRoutines.f90 

        use constants
        use cudaBlock, only: BCData
        use cudaBCPointers, only: ww0, ww1, ww2, pp0, pp1, pp2, &
                              rlv0, rlv1, rlv2, rev0, rev1, rev2, gamma2, &
                              iSize, jSize, iStart, jStart
        use cudaflowVarRefState, only: eddyModel, viscous
        implicit none

        ! Subroutine arguments.
        logical, intent(in) :: secondHalo, correctForK
        integer(kind=intType), intent(in) :: nn

        ! Local variables.
        integer(kind=intType) :: i, j, l, ii
        real(kind=realType), parameter :: twothird = two * third
        real(kind=realType) :: ovg, ovgm1, nnx, nny, nnz
        real(kind=realType) :: pExit, pInt, r, a2, a, ac, ss
        real(kind=realType) :: ue, ve, we, qne, qnh

        ! Loop over the generic subface to set the state in the halo cells.
        
        ! --- Initialize GPU thread indices (assuming 2d launch params) ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then

            ! Store a couple of variables, such as the static
            ! pressure and grid unit outward normal, a bit easier.

            pExit = BCData(nn)%ps(i, j)

            nnx = BCData(nn)%norm(i, j, 1)
            nny = BCData(nn)%norm(i, j, 2)
            nnz = BCData(nn)%norm(i, j, 3)

            ! Abbreviate 1/gamma and 1/(gamma -1) a bit easier.

            ovg = one / gamma2(i, j)
            ovgm1 = one / (gamma2(i, j) - one)

            ! Store the internal pressure and correct for the
            ! possible presence of a k-equation.

            pInt = pp2(i, j)
            if (correctForK) &
                pInt = pInt - twothird * ww2(i, j, irho) * ww2(i, j, itu1)

            ! Compute the velocity components, the normal velocity
            ! and the speed of sound for the internal cell.

            r = one / ww2(i, j, irho)
            a2 = gamma2(i, j) * pInt * r
            a = sqrt(a2)
            ue = ww2(i, j, ivx)
            ve = ww2(i, j, ivy)
            we = ww2(i, j, ivz)
            qne = ue * nnx + ve * nny + we * nnz

            ! Compute the entropy and the acoustic variable.
            ! These riemann inVariants, as well as the tangential
            ! velocity components, are extrapolated.

            ss = pInt * (r**gamma2(i, j))
            ac = qne + two * a * ovgm1

            ! Compute the state in the halo.

            ww1(i, j, irho) = (pExit / ss)**ovg
            pp1(i, j) = pExit
            a = sqrt(gamma2(i, j) * pExit / ww1(i, j, irho))
            qnh = ac - two * a * ovgm1
            ww1(i, j, ivx) = ue + (qnh - qne) * nnx
            ww1(i, j, ivy) = ve + (qnh - qne) * nny
            ww1(i, j, ivz) = we + (qnh - qne) * nnz

            ! Correct the pressure if a k-equation is present.

            if (correctForK) &
                pp1(i, j) = pp1(i, j) &
                            + twothird * ww1(i, j, irho) * ww1(i, j, itu1)

            ! Set the viscosities in the halo to the viscosities
            ! in the donor cell.

            if (viscous) rlv1(i, j) = rlv2(i, j)
            if (eddyModel) rev1(i, j) = rev2(i, j)

        end if

        ! Compute the energy for these halo's.

        call computeEtot(ww1, pp1, correctForK)

        ! Extrapolate the state vectors in case a second halo is needed.

        if (secondHalo) call extrapolate2ndHalo(correctForK)

    end subroutine bcSubsonicOutflow

    attributes(global) subroutine computeEtot(ww, pp, correctForK)

        ! Simplified total energy computation for boundary conditions.
        ! Only implements the constant cpModel

        use constants, only: gammaConstant, third, five, half, irho, one
        use cudaBCPointers, only: iSize, jSize, iStart, jStart
        implicit none

        real(kind=realType), dimension(:, :) :: pp
        real(kind=realType), dimension(:, :, :) :: ww
        logical :: correctForK
        integer(kind=intType) :: ii, i, j
        real(kind=realType) :: ovgm1, factK

        ! Constant cp and thus constant gamma.
        ! Abbreviate 1/(gamma -1) a bit easier.

        ovgm1 = one / (gammaConstant - one)
        factK = ovgm1 * (five * third - gammaConstant)

        ! Loop over the given array and compute the energy, possibly
        ! correcting for K

        ! --- Initialize GPU thread indices (assuming 2d launch params) ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then
            i = mod(ii, isize) + iStart
            j = ii / isize + jStart
            if (.not. correctForK) then
                ww(i, j, iRhoE) = ovgm1 * pp(i, j) &
                                    + half * ww(i, j, irho) * (ww(i, j, ivx)**2 &
                                                                + ww(i, j, ivy)**2 &
                                                                + ww(i, j, ivz)**2)

            else
                ww(i, j, iRhoE) = ovgm1 * pp(i, j) &
                                    + half * ww(i, j, irho) * (ww(i, j, ivx)**2 &
                                                                + ww(i, j, ivy)**2 &
                                                                + ww(i, j, ivz)**2) &
                                    - factK * ww(i, j, irho) * ww(i, j, itu1)
            end if
        end if

        ! Including only constant cp for now
        ! Skipping cpTempCurveFits

    end subroutine computeEtot
        
    attributes(global) subroutine extrapolate2ndHalo(corectForK)

        ! extrapolate2ndHalo determines the states of the second layer
        ! halo cells for the given subface of the block. It is assumed
        ! that the appropriate BCPointers are already set

        use constants, only: irho, two
        use cudaBCPointers, only: ww0, ww1, ww2, pp0, pp1, pp2, &
                              rlv0, rlv1, rlv2, rev0, rev1, rev2, iSize, jSize, iStart, jStart
        use cudaFlowVarRefState, only: viscous, eddyModel
        implicit none

        ! Input variables
        logical, intent(in) :: correctForK

        ! Working variables
        real(kind=realType), parameter :: factor = 0.5_realType
        integer(kind=intType) :: i, j, l, ii

        ! Loop over the generic subface to set the state in the halo cells.

        ! --- Initialize GPU thread indices (assuming 2d launch params) ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then

            ! Extrapolate the density, momentum and pressure.
            ! Make sure that a certain threshold is kept.

            ww0(i, j, irho) = two * ww1(i, j, irho) - ww2(i, j, irho)
            ww0(i, j, irho) = max(factor * ww1(i, j, irho), ww0(i, j, irho))

            ww0(i, j, ivx) = two * ww1(i, j, ivx) - ww2(i, j, ivx)
            ww0(i, j, ivy) = two * ww1(i, j, ivy) - ww2(i, j, ivy)
            ww0(i, j, ivz) = two * ww1(i, j, ivz) - ww2(i, j, ivz)

            pp0(i, j) = max(factor * pp1(i, j), two * pp1(i, j) - pp2(i, j))

            ! The laminar and eddy viscosity, if present. These values
            ! are simply taken constant. Their values do not matter.

            if (viscous) rlv0(i, j) = rlv1(i, j)
            if (eddyModel) rev0(i, j) = rev1(i, j)
        end if

        ! Compute the energy for this halo range.
        call computeEtot(ww0, pp0, correctForK)

    end subroutine extrapolate2ndHalo

end module cudaBCRoutines
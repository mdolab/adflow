! @File    :   cudaBCRoutines.f90
! @Desc    :   This module contains the CUDA Fortran version of the boundary condition routines in BCRoutines.F90

module cudaBCRoutines
    ! --- ADflow modules to use ---
    use precision, only: intType, realType
    ! use cudaBlock, only: cudaBCType => d_BCData, h_BCData
    ! use constants
    ! use cudaBCPointers
    ! --- CUDA FORTRAN module ---
    use cudafor
    ! ALL VARIABLES PASSED AROUND IN THIS MODULE
    use cudaBlock, only: BCData => d_BCData, h_BCData

    implicit none
    



    save

contains

    subroutine copydata
        ! This copies all the necessary data from CPU to GPU
        ! Eventually, this should get thrown out 
        use blockPointers, only: nDom, &
            bBCData => BCData,&
            nBocos, bBCType => BCType
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use utils, only: setPointers
        use constants, only: farField, symm
        
        implicit none

        integer(kind=intType) :: nn, sps, &
                mm, iBeg, iEnd, jBeg, jEnd, &
                inodeBeg, inodeEnd, jnodeBeg, jnodeEnd
        
        domainsLoop: do nn = 1, nDom
            spectralLoop: do sps = 1, nTimeIntervalsSpectral
                
                call setPointers(nn, 1, sps)
                
                allocate(BCData(nBocos))
                allocate(h_BCData(nBocos))

                bocoLoop: do mm = 1, nBocos
                    iBeg = bBCData(mm)%icbeg; iEnd = bBCData(mm)%icend
                    jBeg = bBCData(mm)%jcbeg; jEnd = bBCData(mm)%jcend

                    inodeBeg = bBCData(mm)%inbeg; inodeEnd = bBCData(mm)%inend
                    jnodeBeg = bBCData(mm)%jnbeg; jnodeEnd = bBCData(mm)%jnend

                    select case (bBCType(mm))

                    case (farField)
                        ! Only allocate vars that you need for this type of BC
                        allocate(h_BCData(mm)%rface(iBeg:iEnd, jBeg:jEnd))
                        allocate(h_BCData(mm)%norm(iBeg:iEnd, jBeg:jEnd, 3))

                        h_BCData(mm)%rface = bBCData(mm)%rface
                        h_BCData(mm)%norm = bBCData(mm)%norm
                    
                    ! case (symm, symmPolar)

                    !     h_BCData(mm)%rface = bBCData(mm)%rface
                    !     d_BCData(mm)%rface = bBCData(mm)%rface
                    !     h_BCData(mm)%norm = bBCData(mm)%norm

                    end select

                    ! --- Copy data from host to device via the host type ---
                    BCData(mm) = h_BCData(mm)

                end do bocoLoop

            end do spectralLoop
        end do domainsLoop


    end subroutine copydata

! #ifndef USE_TAPENADE
    subroutine applyAllBC(secondHalo)
        !
        !       applyAllBC applies all boundary conditions for the all
        !       blocks on the grid level currentLevel.
        !
        use constants
        use blockPointers, only: nDom
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
        
        ! manually skip time spectral and multigrid for now
        integer(kind=intType) :: nTimeIntervalsSpectral = 1
        integer(kind=intType) :: currentLevel = 1

        ! Loop over the number of spectral solutions and number of blocks.
        
        spectralLoop: do sps = 1, nTimeIntervalsSpectral
            domains: do nn = 1, nDom
                ! Set the pointers for this block.
                
                call setPointers(nn, currentLevel, sps)

                call interpLevelALEBC_block
                call applyAllBC_block(secondHalo)
                call recoverLevelALEBC_block

            end do domains
        end do spectralLoop

    end subroutine applyAllBC
! #endif

    subroutine applyAllBC_block(secondHalo)

        ! Apply BCs for a single block
        use constants
        use blockPointers, only: nDom, nBocos, BCType, nViscBocos, w, dw, x, vol, il, jl, kl, &
                                 sectionID, wOld, volOld, BCData, &
                                 si, sj, sk, sfacei, sfacej, sfacek, rlv, gamma, p, rev, &
                                 bmtj1, bmtj2, scratch, bmtk2, bmtk1, &
                                 fw, aa, d2wall, bmti1, bmti2, s, &
                                 bib=>ib, bjb=>jb, bkb=>kb, bke=>ke, bje=>je, bie=>ie
        use utils, only: setBCPointers, getCorrectForK
        use BCPointers
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        ! subroutines for copying data to device
        use cudafor, only: dim3
        use cudaBlock, only: copyCudaBCData
        use cudaBCPointers, only: copyCudaBC
        use cudaBCDataMod, only: copyCudaBCDataMod
        implicit none

        ! --- GPU VARS ---
        type(dim3) :: grid_size, block_size
        integer(kind = intType) :: ibmax, jbmax, kbmax

        ! Subroutine arguments.
        logical, intent(in) :: secondHalo

        ! Local variables.
        logical :: correctForK 
        integer(kind=intType) :: nn 
        integer(kind=intType) :: istat, dom, sps
        
        ! --- GPU SETUP ---
        ! Copy block data to device for this block
        ! call copyCudaBlock ! this doesn't work right now
        call copydata ! this is the manual way to do it

        ibmax = 0
        jbmax = 0
        kbmax = 0
        ! zeroing
        do dom = 1,nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(dom,1,sps)
                ibmax = max(ibmax,bib)
                jbmax = max(jbmax,bjb)
                kbmax = max(kbmax,bkb)
            end do
        end do 
        ! SYNCHRONIZE
        istat = cudaDeviceSynchronize()
        
        ! metrics
        print *,ibmax,jbmax,kbmax
        block_size = dim3(8, 8, 1)
        grid_size  = dim3(ceiling(real(ibmax+1) / block_size%x), ceiling(real(jbmax+1) / block_size%y), 1)
        
        ! call CPU_TIME(start)
        
        ! Determine whether or not the total energy must be corrected
        ! for the presence of the turbulent kinetic energy.
        correctForK = getCorrectForK()

        ! Apply all the boundary conditions. The order is important!  Only
        ! some of them have been AD'ed

        ! ! ------------------------------------
        ! !  Symmetry Boundary Condition
        ! ! ------------------------------------
        ! !$AD II-LOOP
        ! do nn = 1, nBocos
        !     if (BCType(nn) == symm) then
        !         ! First set pointers on CPU
        !         call setBCPointers(nn, .False.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         ! Run routine on GPU
        !         call bcSymm1stHalo(nn)
        !     end if
        ! end do

        ! if (secondHalo) then
        !     !$AD II-LOOP
        !     do nn = 1, nBocos
        !         if (BCType(nn) == symm) then
        !             ! First set pointers on CPU
        !             call setBCPointers(nn, .False.)
        !             ! Then copy all data to device
        !             call copyCudaBC
        !             call copyCudaBCDataMod
        !             ! Run routine on GPU
        !             call bcSymm2ndHalo(nn)
        !         end if
        !     end do
        ! end if

        ! ! ------------------------------------
        ! !  Symmetry Polar Boundary Condition
        ! ! ------------------------------------
        ! !$AD II-LOOP
        ! do nn = 1, nBocos
        !     if (BCType(nn) == symmPolar) then
        !         call setBCPointers(nn, .True.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         call bcSymmPolar1stHalo(nn)
        !     end if
        ! end do

        ! if (secondHalo) then
        !     !$AD II-LOOP
        !     do nn = 1, nBocos
        !         if (BCType(nn) == symmPolar) then
        !             call setBCPointers(nn, .True.)
        !             ! Then copy all data to device
        !             call copyCudaBC
        !             call copyCudaBCDataMod
        !             call bcSymmPolar2ndHalo(nn)
        !         end if
        !     end do
        ! end if

        ! ! ------------------------------------
        ! !  Adibatic Wall Boundary Condition
        ! ! ------------------------------------
        ! !$AD II-LOOP
        ! do nn = 1, nViscBocos
        !     if (BCType(nn) == NSWallAdiabatic) then
        !         call setBCPointers(nn, .False.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         call bcNSWallAdiabatic(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ! ------------------------------------
        ! !  Isotermal Wall Boundary Condition
        ! ! ------------------------------------
        ! !$AD II-LOOP
        ! do nn = 1, nViscBocos
        !     if (BCType(nn) == NSWallIsoThermal) then
        !         call setBCPointers(nn, .False.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         call bcNSWallIsothermal(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ------------------------------------
        !  Farfield Boundary Condition
        ! ------------------------------------
        !$AD II-LOOP
        do nn = 1, nBocos
            if (BCType(nn) == farField) then
                call setBCPointers(nn, .False.)
                ! Then copy all data to device
                call copyCudaBCData
                call copyCudaBC
                call copyCudaBCDataMod
                ! Run routine on GPU
                call bcFarField<<<grid_size, block_size>>>(nn, secondHalo, correctForK)
            end if
        end do

        ! ! ------------------------------------
        ! !  Subsonic Outflow Boundary Condition
        ! ! ------------------------------------
        ! do nn = 1, nBocos
        !     if (BCType(nn) == subSonicOutFlow .or. &
        !         BCType(nn) == MassBleedOutflow) then
        !         call setBCPointers(nn, .False.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         call bcSubSonicOutFlow(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ! ------------------------------------
        ! !  Subsonic Inflow Boundary Condition
        ! ! ------------------------------------
        ! do nn = 1, nBocos
        !     if (BCType(nn) == subSonicInFlow) then
        !         call setBCPointers(nn, .False.)
        !         ! Then copy all data to device
        !         call copyCudaBC
        !         call copyCudaBCDataMod
        !         call bcSubSonicInflow(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ! ------------------------------------
        ! !  Extrapolation Boundary Condition
        ! ! ------------------------------------
        ! ! Extrapolation boundary conditions; this also includes
        ! ! the supersonic outflow boundary conditions. The difference
        ! ! between the two is that the extrap boundary conditions
        ! ! correspond to singular lines and supersonic outflow
        ! ! boundaries to physical boundaries. The treatment however
        ! ! is identical.
        ! do nn = 1, nBocos
        !     if (BCType(nn) == extrap .or. &
        !         BCType(nn) == SupersonicOutFlow) then
        !         call setBCPointers(nn, .False.)
        !         call bcExtrap(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ! ------------------------------------
        ! !  Euler Wall Boundary Condition
        ! ! ------------------------------------
        ! !$AD II-LOOP
        ! do nn = 1, nBocos
        !     if (BCType(nn) == EulerWall) then
        !         call setBCPointers(nn, .True.)
        !         call bcEulerWall(nn, secondHalo, correctForK)
        !     end if
        ! end do

        ! ! ------------------------------------
        ! !  Supersonic inflow condition
        ! ! ------------------------------------
        ! do nn = 1, nBocos
        !     if (BCType(nn) == SupersonicInflow) then
        !         call setBCPointers(nn, .False.)
        !         call bcSupersonicInflow(nn, secondHalo, correctForK)
        !     end if
        ! end do

    end subroutine applyAllBC_block

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
        ! use cudaBlock, only: BCData=>d_BCdata
        use cudaflowVarRefState, only: viscous, eddyModel
        use cudaBCPointers, only: gamma1, gamma2, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              iStart, jStart, iSize, jSize, rev1, rev2
        implicit none

        ! Subroutine arguments.
        integer(kind=intType), intent(in), value :: nn

        ! Local variables.
        integer(kind=intType) :: i, j, l
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

    attributes(global) subroutine bcSymm2ndHalo(nn)

        !  bcSymm2ndHalo applies the symmetry boundary conditions to a
        !  block for the 2nd halo. This routine is separate as it makes
        !  AD slightly easier.
        use constants
        ! use cudaBlock, only: BCData => d_BCdata
        use cudaflowVarRefState, only: viscous, eddyModel
        use cudaBCPointers, only: gamma0, gamma3, ww0, ww3, pp0, pp3, rlv0, rlv3, &
                              rev0, rev3, iStart, jStart, iSize, jSize
        implicit none

        ! Subroutine arguments.
        integer(kind=intType), intent(in) :: nn

        ! Local variables.
        integer(kind=intType) :: i, j, l
        real(kind=realType) :: vn, nnx, nny, nnz

        ! If we need the second halo, do everything again, but using ww0,
        ! ww3 etc instead of ww2 and ww1.

         ! --- Initialize GPU thread indices Launch params ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart

        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then

            vn = two * (ww3(i, j, ivx) * BCData(nn)%norm(i, j, 1) + &
                        ww3(i, j, ivy) * BCData(nn)%norm(i, j, 2) + &
                        ww3(i, j, ivz) * BCData(nn)%norm(i, j, 3))

            ! Determine the flow variables in the halo cell.
            ww0(i, j, irho) = ww3(i, j, irho)
            ww0(i, j, ivx) = ww3(i, j, ivx) - vn * BCData(nn)%norm(i, j, 1)
            ww0(i, j, ivy) = ww3(i, j, ivy) - vn * BCData(nn)%norm(i, j, 2)
            ww0(i, j, ivz) = ww3(i, j, ivz) - vn * BCData(nn)%norm(i, j, 3)

            ww0(i, j, irhoE) = ww3(i, j, irhoE)

            ! Set the pressure and gamma and possibly the
            ! laminar and eddy viscosity in the halo.

            gamma0(i, j) = gamma3(i, j)
            pp0(i, j) = pp3(i, j)
            if (viscous) rlv0(i, j) = rlv3(i, j)
            if (eddyModel) rev0(i, j) = rev3(i, j)
        end if

    end subroutine bcSymm2ndHalo

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

    attributes(global) subroutine bcFarfield(nn, secondHalo, correctForK)

        ! bcFarfield applies the farfield boundary condition to a block.
        ! It is assumed that the BCPointers are already set *

        use constants
        ! use blockPointers, only: BCData
        use cudaFlowVarRefState, only: eddyModel, viscous, gammaInf, wInf, pInfCorr
        use cudaBCPointers, only: ww0, ww1, ww2, pp0, pp1, pp2, rlv0, rlv1, rlv2, &
                              rev0, rev1, rev2, gamma2, iStart, jStart, iSize, jSize
        implicit none

        ! Subroutine arguments.
        logical, intent(in), value :: secondHalo, correctForK
        integer, intent(in), value :: nn
        ! Local variables.
        integer(kind=intType) ::  i, j, k, l 

        real(kind=realType) :: nnx, nny, nnz
        real(kind=realType) :: gm1, ovgm1, ac1, ac2
        real(kind=realType) :: r0, u0, v0, w0, qn0, vn0, c0, s0
        real(kind=realType) :: re, ue, ve, we, qne, ce
        real(kind=realType) :: qnf, cf, uf, vf, wf, sf, cc, qq

        ! Some constants needed to compute the riemann inVariants.

        gm1 = gammaInf - one
        ovgm1 = one / gm1

        ! Compute the three velocity components, the speed of sound and
        ! the entropy of the free stream.

        r0 = one / wInf(irho)
        u0 = wInf(ivx)
        v0 = wInf(ivy)
        w0 = wInf(ivz)
        c0 = sqrt(gammaInf * pInfCorr * r0)
        s0 = wInf(irho)**gammaInf / pInfCorr

        ! Loop over the generic subface to set the state in the
        ! halo cells.
        !$AD II-LOOP
        ! --- Initialize GPU thread indices Launch params ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + iStart - 1 ! starts at 0 + iStart
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + jStart - 1 ! starts at 0 + jStart
        if ((i - iStart) <= isize .and. (j - jStart) <= jsize) then
        ! do ii = 0, isize * jsize - 1
        !     i = mod(ii, isize) + iStart
        !     j = ii / isize + jStart

            ! Compute the normal velocity of the free stream and
            ! substract the normal velocity of the mesh.

            qn0 = u0 * BCData(nn)%norm(i, j, 1) + v0 * BCData(nn)%norm(i, j, 2) + w0 * BCData(nn)%norm(i, j, 3)
            vn0 = qn0 - BCData(nn)%rface(i, j)

            ! Compute the three velocity components, the normal
            ! velocity and the speed of sound of the current state
            ! in the internal cell.

            re = one / ww2(i, j, irho)
            ue = ww2(i, j, ivx)
            ve = ww2(i, j, ivy)
            we = ww2(i, j, ivz)
            qne = ue * BCData(nn)%norm(i, j, 1) + ve * BCData(nn)%norm(i, j, 2) + we * BCData(nn)%norm(i, j, 3)
            ce = sqrt(gamma2(i, j) * pp2(i, j) * re)

            ! Compute the new values of the riemann inVariants in
            ! the halo cell. Either the value in the internal cell
            ! is taken (positive sign of the corresponding
            ! eigenvalue) or the free stream value is taken
            ! (otherwise).

            if (vn0 > -c0) then ! Outflow or subsonic inflow.
                ac1 = qne + two * ovgm1 * ce
            else               ! Supersonic inflow.
                ac1 = qn0 + two * ovgm1 * c0
            end if

            if (vn0 > c0) then  ! Supersonic outflow.
                ac2 = qne - two * ovgm1 * ce
            else                     ! Inflow or subsonic outflow.
                ac2 = qn0 - two * ovgm1 * c0
            end if

            qnf = half * (ac1 + ac2)
            cf = fourth * (ac1 - ac2) * gm1

            if (vn0 > zero) then ! Outflow.

                uf = ue + (qnf - qne) * BCData(nn)%norm(i, j, 1)
                vf = ve + (qnf - qne) * BCData(nn)%norm(i, j, 2)
                wf = we + (qnf - qne) * BCData(nn)%norm(i, j, 3)

                sf = ww2(i, j, irho)**gamma2(i, j) / pp2(i, j)

            else
                ! Inflow
                uf = u0 + (qnf - qn0) * BCData(nn)%norm(i, j, 1)
                vf = v0 + (qnf - qn0) * BCData(nn)%norm(i, j, 2)
                wf = w0 + (qnf - qn0) * BCData(nn)%norm(i, j, 3)
                sf = s0

            end if

            ! Compute the density, velocity and pressure in the
            ! halo cell.

            cc = cf * cf / gamma2(i, j)
            qq = uf * uf + vf * vf + wf * wf
            ww1(i, j, irho) = (sf * cc)**ovgm1
            ww1(i, j, ivx) = uf
            ww1(i, j, ivy) = vf
            ww1(i, j, ivz) = wf
            pp1(i, j) = ww1(i, j, irho) * cc

            ! Simply set the laminar and eddy viscosity to
            ! the value in the donor cell. Their values do
            ! not matter too much in the far field.

            if (viscous) rlv1(i, j) = rlv2(i, j)
            if (eddyModel) rev1(i, j) = rev2(i, j)
        end if

        ! Compute the energy for these halos.
        call computeEtot(ww1, pp1, correctForK)

        ! Extrapolate the state vectors in case a second halo
        ! is needed.
        if (secondHalo) call extrapolate2ndHalo(correctForK)

    end subroutine bcFarfield

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

    ! ========================================
    !  DEVICE SUBROUTINES
    ! ========================================
    attributes(device) subroutine computeEtot(ww, pp, correctForK)

        ! Simplified total energy computation for boundary conditions.
        ! Only implements the constant cpModel

        
        use constants, only: gammaConstant, third, five, half, irho, one
        use cudaBCPointers, only: iSize, jSize, iStart, jStart
        implicit none

        real(kind=realType) :: pp
        real(kind=realType) :: ww
        logical :: correctForK
        integer(kind=intType) :: i, j
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
        
    attributes(device) subroutine extrapolate2ndHalo(corectForK)

        ! extrapolate2ndHalo determines the states of the second layer
        ! halo cells for the given subface of the block. It is assumed
        ! that the appropriate BCPointers are already set

        use constants, only: irho, two
        use cudaBCPointers, only: ww0, ww1, ww2, pp0, pp1, pp2, &
                              rlv0, rlv1, rlv2, rev0, rev1, rev2, iSize, jSize, iStart, jStart
        use cudaFlowVarRefState, only: viscous, eddyModel
        implicit none

        ! Input variables
        logical :: correctForK

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
module smoothers

contains
    subroutine RungeKuttaSmoother
        !
        !       RungeKuttaSmoother performs one multi-stage runge kutta
        !       explicit time step for the current multigrid level. On
        !       entrance it is assumed that the residual and time step are
        !       already computed. On exit the solution in the halo's contain
        !       the latest values. However, the residual corresponding to
        !       these values is not computed.
        !
        use constants
        use blockPointers, only: w, p, wn, pn, il, jl, kl, nDom
        use flowVarRefState, only: nwf
        use inputIteration, only: nRKStages
        use inputTimeSpectral, only: ntimeIntervalsSpectral
        use iteration, only: currentLevel, rkStage
        use utils, only: setPointers
        use residuals, only: initRes, residual, residualAveraging, sourceTerms
        use BCRoutines, only: applyAllBC
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: sps, nn, i, j, k, l

        ! Store the variables of the zeroth runge kutta stage.

        spectralLoop: do sps = 1, nTimeIntervalsSpectral
            domains: do nn = 1, nDom

                ! Set the pointers to this block.

                call setPointers(nn, currentLevel, sps)

                ! The variables stored in w.

                do l = 1, nwf
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                wn(i, j, k, l) = w(i, j, k, l)
                            end do
                        end do
                    end do
                end do

                ! And the pressure.

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            pn(i, j, k) = p(i, j, k)
                        end do
                    end do
                end do

            end do domains
        end do spectralLoop

        ! Loop over all but the last Runge Kutta stages. Note that the
        ! counter variable rkStage is defined in the module iteration.

        do rkStage = 1, (nRKStages - 1)

            ! Execute a Runge Kutta stage and exchange the externals.

            call executeRkStage

            ! Compute the residuals for the next stage.

            call initres(1_intType, nwf)
            call sourceTerms()
            call residual

        end do

        ! Execute the last RK stage. Set rkStage to nRKStages, for
        ! clarity; after the previous loop rkStage == nRKStages.

        rkStage = nRKStages

        call executeRkStage

    end subroutine RungeKuttaSmoother

    !      ==================================================================

    subroutine executeRkStage
        !
        !       executeRkStage executes one runge kutta stage. The stage
        !       number, rkStage, is defined in the local module iteration.
        !
        use blockPointers
        use constants
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use inputTimeSpectral
        use inputUnsteady
        use iteration
        use inputDiscretization
        use haloExchange, only: whalo1, whalo2
        use utils, only: setPointers
        use flowUtils, only: computeEtotBlock, computeLamViscosity
        use turbutils, only: computeeddyviscosity
        use residuals, only: residualAveraging
        use BCRoutines, only: applyAllBC

        implicit none
        !
        !      Local parameter.
        !
        real(kind=realType), parameter :: fiveThird = five * third
        !
        !      Local variables.
        !
        integer(kind=intType) :: sps, nn, i, j, k, l

        real(kind=realType) :: tmp, unsteadyImpl, mult
        real(kind=realType) :: dt, currentCfl, gm1, gm53
        real(kind=realType) :: v2, ovr, dp, factK, ru, rv, rw

        logical :: secondHalo, smoothResidual, correctForK

        ! Set the value of secondHalo and the current cfl number,
        ! depending on the situation. On the finest grid in the mg cycle
        ! the second halo is computed, otherwise not.

        if (currentLevel <= groundLevel) then
            secondHalo = .true.
        else
            secondHalo = .false.
        end if

        currentCfl = cflCoarse
        if (currentLevel == 1) then
            currentCfl = cfl
        end if

        ! Determine whether or not residual averaging must be applied.

        if (resAveraging == noResAveraging) then
            smoothResidual = .false.
        else if (resAveraging == alwaysResAveraging) then
            smoothResidual = .true.
        else if (mod(rkStage, 2_intType) == 1) then
            smoothResidual = .true.
        else
            smoothResidual = .false.
        end if

        ! Determine whether or not the total energy must be corrected
        ! for the presence of the turbulent kinetic energy.

        if (kPresent) then
            if ((currentLevel <= groundLevel)) then
                correctForK = .true.
            else
                correctForK = .false.
            end if
        else
            correctForK = .false.
        end if
        !
        !       Compute the updates of the conservative variables.
        !
        ! Loop over the local number of blocks.

        domainsUpdate: do nn = 1, nDom

            ! Determine the equation mode solved.

            select case (equationMode)

            case (steady, timeSpectral)

                ! Steady equations, including time spectral. Everything is
                ! solved explicitly. Store the cfl number times the RK
                ! coefficient in tmp.

                tmp = currentCfl * etaRk(rkStage)

                ! Loop over the number of spectral solutions. Note that
                ! for the steady mode this value is 1.

                spectralSteady: do sps = 1, nTimeIntervalsSpectral

                    ! Set the pointers to this block.

                    call setPointers(nn, currentLevel, sps)

                    ! Loop over the owned cells of this block.

                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! Determine the local time step (multiplied by the
                                ! current rk coefficient).
                                if (lowspeedpreconditioner) then
                                    dt = 0.8 * tmp * dtl(i, j, k)
                                else
                                    dt = tmp * dtl(i, j, k)
                                end if

                                ! Compute the updates of the flow field variables.

                                dw(i, j, k, irho) = dw(i, j, k, irho) * dt
                                dw(i, j, k, imx) = dw(i, j, k, imx) * dt
                                dw(i, j, k, imy) = dw(i, j, k, imy) * dt
                                dw(i, j, k, imz) = dw(i, j, k, imz) * dt
                                dw(i, j, k, irhoE) = dw(i, j, k, irhoE) * dt

                            end do
                        end do
                    end do

                end do spectralSteady

                !=============================================================

            case (unsteady)

                ! Unsteady equations are solved via the dual time
                ! stepping technique. The leading term of the time
                ! integrator must be treated implicitly for stability
                ! reasons. This leads to a different multiplication
                ! factor of the residual compared to the steady case.
                ! Compute this additional term and store the cfl number
                ! times the rk coefficient in tmp.

                unsteadyImpl = coefTime(0) * timeRef / deltaT
                tmp = currentCfl * etaRk(rkStage)

                ! Loop over the number of spectral modes, although this is
                ! always 1 for the unsteady mode. The loop is executed for
                ! consistency reasons.

                spectralUnsteady: do sps = 1, nTimeIntervalsSpectral

                    ! Set the pointers to this block.

                    call setPointers(nn, currentLevel, sps)

                    ! Determine the updates of the flow field variables.
                    ! Owned cells only. The result is stored in dw.

                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! Determine the local time step (multiplied by the
                                ! current rk coefficient) and the multiplication
                                ! factor for the residuals.

                                dt = tmp * dtl(i, j, k)
                                mult = dt / (dt * unsteadyImpl * vol(i, j, k) + one)

                                ! Compute the updates of the flow field variables.

                                dw(i, j, k, irho) = dw(i, j, k, irho) * mult
                                dw(i, j, k, imx) = dw(i, j, k, imx) * mult
                                dw(i, j, k, imy) = dw(i, j, k, imy) * mult
                                dw(i, j, k, imz) = dw(i, j, k, imz) * mult
                                dw(i, j, k, irhoE) = dw(i, j, k, irhoE) * mult

                            end do
                        end do
                    end do

                end do spectralUnsteady

            end select

        end do domainsUpdate
        !
        !       Compute the new state vector.
        !
        ! Loop over the number of spectral solutions and local blocks.

        spectralLoop: do sps = 1, nTimeIntervalsSpectral
            domainsState: do nn = 1, nDom

                ! Set the pointers to this block.

                call setPointers(nn, currentLevel, sps)

                ! Possibility to smooth the updates.

                if (smoothResidual) then
                    call residualAveraging
                end if
                ! Flow variables.

                factK = zero
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! Store gamma -1 and gamma - 5/3 a bit easier.

                            gm1 = gamma(i, j, k) - one
                            gm53 = gamma(i, j, k) - fiveThird

                            ! Compute the pressure update from the conservative
                            ! updates. The expression used below is valid even if
                            ! cp is not constant. For the calorically perfect case,
                            ! cp is constant, it can be simplified to the usual
                            ! expression, but not for variable cp.

                            ovr = one / w(i, j, k, irho)
                            v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
                            if (correctForK) factK = gm53 * w(i, j, k, itu1)

                            dp = (ovr * p(i, j, k) + factK &
                                  - gm1 * (ovr * w(i, j, k, irhoE) - v2)) * dw(i, j, k, irho) &
                                 + gm1 * (dw(i, j, k, irhoE) - w(i, j, k, ivx) * dw(i, j, k, imx) &
                                          - w(i, j, k, ivy) * dw(i, j, k, imy) &
                                          - w(i, j, k, ivz) * dw(i, j, k, imz))

                            ! Compute the density. Correct for negative values.

                            w(i, j, k, irho) = wn(i, j, k, irho) - dw(i, j, k, irho)
                            w(i, j, k, irho) = max(w(i, j, k, irho), 1.e-4_realType * rhoInf)

                            ! Compute the velocities.

                            ru = wn(i, j, k, irho) * wn(i, j, k, ivx) - dw(i, j, k, imx)
                            rv = wn(i, j, k, irho) * wn(i, j, k, ivy) - dw(i, j, k, imy)
                            rw = wn(i, j, k, irho) * wn(i, j, k, ivz) - dw(i, j, k, imz)

                            ovr = one / w(i, j, k, irho)
                            w(i, j, k, ivx) = ovr * ru
                            w(i, j, k, ivy) = ovr * rv
                            w(i, j, k, ivz) = ovr * rw

                            ! Compute the pressure. Correct for negative values.

                            p(i, j, k) = pn(i, j, k) - dp
                            p(i, j, k) = max(p(i, j, k), 1.e-4_realType * pInfCorr)

                        end do
                    end do
                end do

                ! Compute the total energy and possibly the laminar and eddy
                ! viscosity in the owned cells.

                call computeEtotBlock(2_intType, il, 2_intType, jl, &
                                      2_intType, kl, correctForK)
                call computeLamViscosity(.False.)
                call computeEddyViscosity(.False.)

            end do domainsState
        end do spectralLoop

        ! Exchange the pressure if the pressure must be exchanged early.
        ! Only the first halo's are needed, thus whalo1 is called.
        ! Only on the fine grid.

        if (exchangePressureEarly .and. currentLevel <= groundLevel) &
            call whalo1(currentLevel, 1_intType, 0_intType, .true., &
                        .false., .false.)

        ! Apply all boundary conditions to all blocks on this level.

        call applyAllBC(secondHalo)

        ! Exchange the solution. Either whalo1 or whalo2
        ! must be called.

        if (secondHalo) then
            call whalo2(currentLevel, 1_intType, nwf, .true., &
                        .true., .true.)
        else
            call whalo1(currentLevel, 1_intType, nwf, .true., &
                        .true., .true.)
        end if

    end subroutine executeRkStage
    subroutine DADISmoother
        !
        !       RungeKuttaSmoother performs one multi-stage runge kutta
        !       explicit time step for the current multigrid level. On
        !       entrance it is assumed that the residual and time step are
        !       already computed. On exit the solution in the halo's contain
        !       the latest values. However, the residual corresponding to
        !       these values is not computed.
        !
        use blockPointers
        use flowVarRefState
        use inputIteration
        use inputTimeSpectral
        use iteration
        use residuals, only: initRes, residual, sourceTerms
        implicit none

        if (groundLevel == 1) then
            do Subit = 1, nSubiterations - 1

                ! Execute a DADI step and exchange the externals.

                call executeDADIStep

                ! Compute the residuals for the next stage.
                call initres(1_intType, nwf)
                call sourceTerms()
                call residual

            end do

            ! Set Subit to nSubiterations, for clarity;
            Subit = nSubiterations
        end if

        ! Execute the last subiteration.
        call executeDADIStep

    end subroutine DADISmoother

    !      ==================================================================

    subroutine executeDADIStep
        !
        !       executeDADIStep executes one DADI step.
        !
        use blockPointers
        use constants
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use inputTimeSpectral
        use inputUnsteady
        use iteration
        use utils, only: getCorrectForK, setPointers
        use haloExchange, only: whalo1, whalo2
        use flowUtils, only: computeETotBlock, computeLamViscosity
        use turbutils, only: computeeddyviscosity
        use residuals, only: residualAveraging, computeDwDADI
        use BCRoutines, only: applyAllBC
        implicit none
        !
        !      Local parameter.
        !
        real(kind=realType), parameter :: fiveThird = five * third
        !
        !      Local variables.
        !
        integer(kind=intType) :: sps, nn, i, j, k, l

        real(kind=realType) :: unsteadyImpl, mult
        real(kind=realType) :: dt, currentCfl, gm1, gm53
        real(kind=realType) :: v2, ovr, dp, factK, ru, rv, rw

        logical :: secondHalo, smoothResidual, correctForK

        ! Set the value of secondHalo and the current cfl number,
        ! depending on the situation. On the finest grid in the mg cycle
        ! the second halo is computed, otherwise not.

        if (currentLevel <= groundLevel) then
            secondHalo = .true.
        else
            secondHalo = .false.
        end if

        currentCfl = cflCoarse
        if (currentLevel == 1) then
            currentCfl = cfl
        end if

        ! Determine whether or not residual averaging must be applied.

        if (resAveraging == noResAveraging) then
            smoothResidual = .false.
        else if (resAveraging == alwaysResAveraging) then
            smoothResidual = .true.
        else if (mod(rkStage, 2_intType) == 1) then
            smoothResidual = .true.
        else
            smoothResidual = .false.
        end if

        ! Determine whether or not the total energy must be corrected
        ! for the presence of the turbulent kinetic energy.
        correctForK = getCorrectForK()

        !
        !       Compute the updates of the conservative variables.
        !
        ! Loop over the local number of blocks.

        domainsUpdate: do nn = 1, nDom

            ! Determine the equation mode solved.

            select case (equationMode)

            case (steady, timeSpectral)

                ! Loop over the number of spectral solutions. Note that
                ! for the steady mode this value is 1.

                spectralSteady: do sps = 1, nTimeIntervalsSpectral

                    ! Set the pointers to this block.

                    call setPointers(nn, currentLevel, sps)

                    ! Loop over the owned cells of this block.

                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! Determine the local time step

                                dt = -currentCfl * dtl(i, j, k) * vol(i, j, k)

                                ! Compute the updates of the flow field variables.

                                dw(i, j, k, irho) = dw(i, j, k, irho) * dt
                                dw(i, j, k, imx) = dw(i, j, k, imx) * dt
                                dw(i, j, k, imy) = dw(i, j, k, imy) * dt
                                dw(i, j, k, imz) = dw(i, j, k, imz) * dt
                                dw(i, j, k, irhoE) = dw(i, j, k, irhoE) * dt

                            end do
                        end do
                    end do

                    call computedwDADI

                end do spectralSteady

                !=============================================================

            case (unsteady)

                ! Unsteady equations are solved via the dual time
                ! stepping technique. The leading term of the time
                ! integrator must be treated implicitly for stability
                ! reasons. This leads to a different multiplication
                ! factor of the residual compared to the steady case.

                unsteadyImpl = coefTime(0) * timeRef / deltaT

                ! Loop over the number of spectral modes, although this is
                ! always 1 for the unsteady mode. The loop is executed for
                ! consistency reasons.

                spectralUnsteady: do sps = 1, nTimeIntervalsSpectral

                    ! Set the pointers to this block.

                    call setPointers(nn, currentLevel, sps)

                    ! Determine the updates of the flow field variables.
                    ! Owned cells only. The result is stored in dw.

                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! Determine the local time step

                                dt = currentCfl * dtl(i, j, k)
                                mult = dt / (dt * unsteadyImpl * vol(i, j, k) + one)
                                mult = -mult * vol(i, j, k)

                                dw(i, j, k, irho) = dw(i, j, k, irho) * mult
                                dw(i, j, k, imx) = dw(i, j, k, imx) * mult
                                dw(i, j, k, imy) = dw(i, j, k, imy) * mult
                                dw(i, j, k, imz) = dw(i, j, k, imz) * mult
                                dw(i, j, k, irhoE) = dw(i, j, k, irhoE) * mult

                            end do
                        end do
                    end do

                    call computedwDADI

                end do spectralUnsteady

            end select

        end do domainsUpdate
        !
        !       Compute the new state vector.
        !
        ! Loop over the number of spectral solutions and local blocks.

        spectralLoop: do sps = 1, nTimeIntervalsSpectral
            domainsState: do nn = 1, nDom

                ! Set the pointers to this block.

                call setPointers(nn, currentLevel, sps)

                ! Possibility to smooth the updates.

                if (smoothResidual) call residualAveraging

                ! Flow variables.

                factK = zero
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! Store gamma -1 and gamma - 5/3 a bit easier.

                            gm1 = gamma(i, j, k) - one
                            gm53 = gamma(i, j, k) - fiveThird

                            ! Compute the pressure update from the conservative
                            ! updates. The expression used below is valid even if
                            ! cp is not constant. For the calorically perfect case,
                            ! cp is constant, it can be simplified to the usual
                            ! expression, but not for variable cp.

                            ovr = one / w(i, j, k, irho)
                            v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
                            if (correctForK) factK = gm53 * w(i, j, k, itu1)

                            dp = (ovr * p(i, j, k) + factK &
                                  - gm1 * (ovr * w(i, j, k, irhoE) - v2)) * dw(i, j, k, irho) &
                                 + gm1 * (dw(i, j, k, irhoE) - w(i, j, k, ivx) * dw(i, j, k, imx) &
                                          - w(i, j, k, ivy) * dw(i, j, k, imy) &
                                          - w(i, j, k, ivz) * dw(i, j, k, imz))

                            ! Compute the velocities.

                            ru = w(i, j, k, irho) * w(i, j, k, ivx) - dw(i, j, k, imx)
                            rv = w(i, j, k, irho) * w(i, j, k, ivy) - dw(i, j, k, imy)
                            rw = w(i, j, k, irho) * w(i, j, k, ivz) - dw(i, j, k, imz)

                            ! Compute the density. Correct for negative values.

                            w(i, j, k, irho) = w(i, j, k, irho) - dw(i, j, k, irho)
                            w(i, j, k, irho) = max(w(i, j, k, irho), 1.e-4_realType * rhoInf)

                            ovr = one / w(i, j, k, irho)
                            w(i, j, k, ivx) = ovr * ru
                            w(i, j, k, ivy) = ovr * rv
                            w(i, j, k, ivz) = ovr * rw

                            ! Compute the pressure. Correct for negative values.

                            p(i, j, k) = p(i, j, k) - dp
                            p(i, j, k) = max(p(i, j, k), 1.e-4_realType * pInfCorr)

                        end do
                    end do
                end do

                ! Compute the total energy and possibly the laminar and eddy
                ! viscosity in the owned cells.

                call computeEtotBlock(2_intType, il, 2_intType, jl, &
                                      2_intType, kl, correctForK)
                call computeLamViscosity(.False.)
                call computeEddyViscosity(.False.)

            end do domainsState
        end do spectralLoop

        ! Exchange the pressure if the pressure must be exchanged early.
        ! Only the first halo's are needed, thus whalo1 is called.
        ! Only on the fine grid.

        if (exchangePressureEarly .and. currentLevel <= groundLevel) &
            call whalo1(currentLevel, 1_intType, 0_intType, .true., &
                        .false., .false.)

        ! Apply all boundary conditions to all blocks on this level.

        call applyAllBC(secondHalo)

        ! Exchange the solution. Either whalo1 or whalo2
        ! must be called.

        if (secondHalo) then
            call whalo2(currentLevel, 1_intType, nwf, .true., &
                        .true., .true.)
        else
            call whalo1(currentLevel, 1_intType, nwf, .true., &
                        .true., .true.)
        end if

    end subroutine executeDADIStep

end module smoothers

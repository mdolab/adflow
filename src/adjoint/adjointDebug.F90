! This module is used for debugging and testing only

module adjointDebug

contains

#ifndef USE_COMPLEX

    subroutine computeMatrixFreeProductFwdFD(xvdot, extradot, wdot, bcDataValuesdot, &
                                             useSpatial, useState, famLists, &
                                             bcDataNames, bcDataValues, bcDataFamLists, bcVarsEmpty, &
                                             dwdot, funcsDot, fDot, &
                                             costSize, fSize, nTime, h)

        ! This routine is used to debug master_d. It uses the forward seeds to set perturbations
        ! and then computes the value of the derivatives using forward finite diffenece

        use constants
        use adjointvars
        use blockPointers, only: nDom
        use communication, only: adflow_comm_world
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputPhysics, only: pointRef, alpha, beta, equations, machCoef, &
                                mach, machGrid, rgasdim
        use iteration, only: currentLevel, groundLevel
        use flowVarRefState, only: pInfDim, rhoInfDim, TinfDim
        use blockPointers, only: nDom, il, jl, kl, wd, x, w, dw, dwd, nBocos, nViscBocos

        use adjointUtils, only: allocDerivativeValues, zeroADSeeds
        use masterRoutines, only: master
        use utils, only: isWallType, setPointers, setPointers_d, EChk
        use flowVarRefState, only: nw, nwf
        use wallDistanceData, only: xSurf, xSurfVec
        use wallDistance, only: updateXSurf
        implicit none
        !
        ! Input Variables
        !
        ! derivative seeds
        real(kind=realType), dimension(:), intent(in) :: xvDot
        real(kind=realType), dimension(:), intent(in) :: extraDot
        real(kind=realType), dimension(:), intent(in) :: wDot

        ! size data
        logical, intent(in) :: useSpatial, useState
        integer(kind=intType), dimension(:, :) :: famLists
        integer(kind=intType) :: costSize, fSize, nTime

        character, dimension(:, :), intent(in) :: bcDataNames
        real(kind=realType), dimension(:), intent(inout) :: bcDataValues, bcDataValuesDot
        integer(kind=intType), dimension(:, :) :: bcDataFamLists
        logical, intent(in) :: BCVarsEmpty

        ! Finite difference parameters
        real(kind=realType), intent(in) :: h ! step size for Finite Difference

        !
        ! Ouput Variables
        !
        ! Output derivative seeds
        real(kind=realType), dimension(size(wDot)), intent(out) :: dwDot
        real(kind=realType), dimension(costSize, size(famLists, 1)), intent(out) :: funcsDot
        real(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot

        !
        ! Working Variables
        !
        integer(kind=intType) :: nn, sps, level
        integer(kind=intType) :: ierr, mm, i, j, k, l, ii, jj, iRegion

        real(kind=realType), dimension(costSize, size(famLists, 1)) :: funcs

        ! Input Arguments for master:
        real(kind=realType), dimension(costSize, size(famLists, 1)) :: funcValues
        real(kind=realType), dimension(:, :, :), allocatable :: forces

        fSize = size(fDot, 2)
        allocate (forces(3, fSize, nTimeIntervalsSpectral))

        ! Need to trick the residual evalution to use coupled (mean flow and
        ! turbulent) together.
        level = 1
        currentLevel = level
        groundLevel = level

        ! Allocate the memory we need for derivatives if not done so
        ! already. Note this isn't deallocated until the adflow is
        ! destroyed.

        ! a bit over kill for our needs, but so values are stored in the seeds
        if (.not. derivVarsAllocated) then
            call allocDerivativeValues(level)
        end if

        ! Zero all AD seesd.
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call zeroADSeeds(nn, level, sps)
            end do
        end do

        ! ----------------------------- Run Master ---------------------------------
        ! Run the super-dee-duper master rotuine
        if (bcVarsEmpty) then
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces)

        else
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces, &
                        bcDataNames, bcDataValues, bcDataFamLists)
        end if

        ! Copy base val (f) into variables for the final vals (f(x+dx) - f)/dx
        ! we add the negative sign here instead of doing it later

        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers_d(nn, 1, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                ii = ii + 1
                                dwd(i, j, k, l) = -dw(i, j, k, l)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        fDot = -forces
        funcsDot = -funcValues

        !  --------------------- apply the perturbations ----------------------------
        ! Set the extra seeds now do the extra ones. Note that we are assuming the
        ! machNumber used for the coefficients follows the Mach number,
        ! not the grid mach number.

        alpha = alpha + h * extraDot(iAlpha)
        beta = beta + h * extraDot(iBeta)
        mach = mach + h * extraDot(iMach)
        machCoef = machCoef + h * extraDot(iMach)
        machGrid = machGrid + h * extraDot(iMachGrid)
        PinfDim = PinfDim + h * extraDot(iPressure)
        rhoinfDim = rhoinfDim + h * extraDot(iDensity)
        tinfdim = tinfdim + h * extraDot(iTemperature)
        pointref(1) = pointref(1) + h * extraDot(iPointRefX)
        pointref(2) = pointref(2) + h * extraDot(iPointRefY)
        pointref(3) = pointref(3) + h * extraDot(iPointRefZ)
        rgasdim = rgasdim + h * zero

        ! Set the provided w and x seeds:
        ii = 0
        jj = 0
        domainLoop1: do nn = 1, nDom
            spectalLoop1: do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do l = 1, 3
                                ii = ii + 1
                                x(i, j, k, l) = x(i, j, k, l) + xvDot(ii) * h
                            end do
                        end do
                    end do
                end do
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                jj = jj + 1
                                w(i, j, k, l) = w(i, j, k, l) + wDot(jj) * h
                            end do
                        end do
                    end do
                end do
            end do spectalLoop1
        end do domainLoop1

        bcDataValues = bcDataValues + bcDataValuesdot * h

        ! The xvolume that is set by used by update Xsurf is only allocated with
        ! rans
        if (equations == RANSEquations) then
            call updateXSurf(level)
        end if

        ! ----------------------------- Run Master ---------------------------------
        ! Run the super-dee-duper master rotuine
        if (bcVarsEmpty) then
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces)

        else
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces, &
                        bcDataNames, bcDataValues, bcDataFamLists)

        end if

        ! ------------------------- set the output vectors ----------------------
        ! Set the extra seeds now do the extra ones. Note that we are assuming the
        ! machNumber used for the coefficients follows the Mach number,
        ! not the grid mach number.

        alpha = alpha - h * extraDot(iAlpha)
        beta = beta - h * extraDot(iBeta)
        mach = mach - h * extraDot(iMach)
        machCoef = machCoef - h * extraDot(iMach)
        machGrid = machGrid - h * extraDot(iMachGrid)
        PinfDim = PinfDim - h * extraDot(iPressure)
        rhoinfDim = rhoinfDim - h * extraDot(iDensity)
        tinfdim = tinfdim - h * extraDot(iTemperature)
        pointref(1) = pointref(1) - h * extraDot(iPointRefX)
        pointref(2) = pointref(2) - h * extraDot(iPointRefY)
        pointref(3) = pointref(3) - h * extraDot(iPointRefZ)
        rgasdim = rgasdim - h * zero

        ! Copy out the residual derivative into the provided dwDot and remove the
        ! perturbation
        ii = 0
        jj = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers_d(nn, 1, sps)
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do l = 1, 3
                                ii = ii + 1
                                x(i, j, k, l) = x(i, j, k, l) - xvDot(ii) * h
                            end do
                        end do
                    end do
                end do
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                jj = jj + 1
                                w(i, j, k, l) = w(i, j, k, l) - wDot(jj) * h
                                dwd(i, j, k, l) = (dwd(i, j, k, l) + dw(i, j, k, l)) / h
                                dwDot(jj) = dwd(i, j, k, l) ! copy values to output
                            end do
                        end do
                    end do
                end do

            end do
        end do

        if (equations == RANSEquations) then
            call updateXSurf(level)
        end if
        bcDataValues = bcDataValues - bcDataValuesdot * h

        fDot = (fDot + forces) / h
        funcsDot = (funcsDot + funcValues) / h

    end subroutine computeMatrixFreeProductFwdFD

    subroutine printADSeeds(nn, level, sps)
!DIR$ NOOPTIMIZE
        ! this routine is used for debugging master_d, and master_b.
        ! it prints all the AD seeds used
        ! it is useful to save the output and compare it with a diff tool

        use constants
        use block, only: flowDomsd, flowDoms
        use blockPointers
        use inputTimeSpectral
        use flowVarRefState
        use inputPhysics
        use BCPointers_b
        use communication
        use oversetData, only: oversetPresent
        use cgnsGrid, only: cgnsDoms, cgnsDomsd, cgnsNDom
        use actuatorRegionData, only: nActuatorRegions, actuatorRegionsd
        implicit none

        ! Input parameters
        integer(kind=intType) :: nn, level, sps

        ! Working parameters
        integer(kind=intType) :: mm, i, iDom
        integer(kind=intType) :: iBoco, iData, iDirichlet
        write (*, *) 'ADSeeds for block', nn, ' at level ', level, ' at sps ', sps
        write (*, *) 'd2wall ', minval(flowDomsd(nn, level, sps)%d2wall), &
            maxval(flowDomsd(nn, level, sps)%d2wall), &
            norm2(flowDomsd(nn, level, sps)%d2wall)
        write (*, *) 'x ', minval(flowDomsd(nn, level, sps)%x), &
            maxval(flowDomsd(nn, level, sps)%x), &
            norm2(flowDomsd(nn, level, sps)%x)
        write (*, *) 'si ', minval(flowDomsd(nn, level, sps)%si), &
            maxval(flowDomsd(nn, level, sps)%si), &
            norm2(flowDomsd(nn, level, sps)%si)
        write (*, *) 'sj ', minval(flowDomsd(nn, level, sps)%sj), &
            maxval(flowDomsd(nn, level, sps)%sj), &
            norm2(flowDomsd(nn, level, sps)%sj)
        write (*, *) 'sk ', minval(flowDomsd(nn, level, sps)%sk), &
            maxval(flowDomsd(nn, level, sps)%sk), &
            norm2(flowDomsd(nn, level, sps)%sk)
        write (*, *) 'vol ', minval(flowDomsd(nn, level, sps)%vol), &
            maxval(flowDomsd(nn, level, sps)%vol), &
            norm2(flowDomsd(nn, level, sps)%vol)

        write (*, *) 's ', minval(flowDomsd(nn, level, sps)%s), &
            maxval(flowDomsd(nn, level, sps)%s), &
            norm2(flowDomsd(nn, level, sps)%s)
        write (*, *) 'sFaceI ', minval(flowDomsd(nn, level, sps)%sFaceI), &
            maxval(flowDomsd(nn, level, sps)%sFaceI), &
            norm2(flowDomsd(nn, level, sps)%sFaceI)
        write (*, *) 'sFaceJ ', minval(flowDomsd(nn, level, sps)%sFaceJ), &
            maxval(flowDomsd(nn, level, sps)%sFaceJ), &
            norm2(flowDomsd(nn, level, sps)%sFaceJ)
        write (*, *) 'sFaceK ', minval(flowDomsd(nn, level, sps)%sFaceK), &
            maxval(flowDomsd(nn, level, sps)%sFaceK), &
            norm2(flowDomsd(nn, level, sps)%sFaceK)

        write (*, *) 'w ', minval(flowDomsd(nn, level, sps)%w), &
            maxval(flowDomsd(nn, level, sps)%w), &
            norm2(flowDomsd(nn, level, sps)%w)
        write (*, *) 'dw ', minval(flowDomsd(nn, level, sps)%dw), &
            maxval(flowDomsd(nn, level, sps)%dw), &
            norm2(flowDomsd(nn, level, sps)%dw)
        write (*, *) 'fw ', minval(flowDomsd(nn, level, sps)%fw), &
            maxval(flowDomsd(nn, level, sps)%fw), &
            norm2(flowDomsd(nn, level, sps)%fw)
        write (*, *) 'scratch ', minval(flowDomsd(nn, level, sps)%scratch), &
            maxval(flowDomsd(nn, level, sps)%scratch), &
            norm2(flowDomsd(nn, level, sps)%scratch)

        write (*, *) 'p ', minval(flowDomsd(nn, level, sps)%p), &
            maxval(flowDomsd(nn, level, sps)%p), &
            norm2(flowDomsd(nn, level, sps)%p)
        write (*, *) 'gamma ', minval(flowDomsd(nn, level, sps)%gamma), &
            maxval(flowDomsd(nn, level, sps)%gamma), &
            norm2(flowDomsd(nn, level, sps)%gamma)
        write (*, *) 'aa ', minval(flowDomsd(nn, level, sps)%aa), &
            maxval(flowDomsd(nn, level, sps)%aa), &
            norm2(flowDomsd(nn, level, sps)%aa)

        write (*, *) 'rlv ', minval(flowDomsd(nn, level, sps)%rlv), &
            maxval(flowDomsd(nn, level, sps)%rlv), &
            norm2(flowDomsd(nn, level, sps)%rlv)
        write (*, *) 'rev ', minval(flowDomsd(nn, level, sps)%rev), &
            maxval(flowDomsd(nn, level, sps)%rev), &
            norm2(flowDomsd(nn, level, sps)%rev)

        write (*, *) 'radI ', minval(flowDomsd(nn, level, sps)%radI), &
            maxval(flowDomsd(nn, level, sps)%radI), &
            norm2(flowDomsd(nn, level, sps)%radI)
        write (*, *) 'radJ ', minval(flowDomsd(nn, level, sps)%radJ), &
            maxval(flowDomsd(nn, level, sps)%radJ), &
            norm2(flowDomsd(nn, level, sps)%radJ)
        write (*, *) 'radK ', minval(flowDomsd(nn, level, sps)%radK), &
            maxval(flowDomsd(nn, level, sps)%radK), &
            norm2(flowDomsd(nn, level, sps)%radK)

        write (*, *) 'ux ', minval(flowDomsd(nn, level, sps)%ux), &
            maxval(flowDomsd(nn, level, sps)%ux), &
            norm2(flowDomsd(nn, level, sps)%ux)
        write (*, *) 'uy ', minval(flowDomsd(nn, level, sps)%uy), &
            maxval(flowDomsd(nn, level, sps)%uy), &
            norm2(flowDomsd(nn, level, sps)%uy)
        write (*, *) 'uz ', minval(flowDomsd(nn, level, sps)%uz), &
            maxval(flowDomsd(nn, level, sps)%uz), &
            norm2(flowDomsd(nn, level, sps)%uz)
        write (*, *) 'vx ', minval(flowDomsd(nn, level, sps)%vx), &
            maxval(flowDomsd(nn, level, sps)%vx), &
            norm2(flowDomsd(nn, level, sps)%vx)
        write (*, *) 'vy ', minval(flowDomsd(nn, level, sps)%vy), &
            maxval(flowDomsd(nn, level, sps)%vy), &
            norm2(flowDomsd(nn, level, sps)%vy)
        write (*, *) 'vz ', minval(flowDomsd(nn, level, sps)%vz), &
            maxval(flowDomsd(nn, level, sps)%vz), &
            norm2(flowDomsd(nn, level, sps)%vz)
        write (*, *) 'wx ', minval(flowDomsd(nn, level, sps)%wx), &
            maxval(flowDomsd(nn, level, sps)%wx), &
            norm2(flowDomsd(nn, level, sps)%wx)
        write (*, *) 'wy ', minval(flowDomsd(nn, level, sps)%wy), &
            maxval(flowDomsd(nn, level, sps)%wy), &
            norm2(flowDomsd(nn, level, sps)%wy)
        write (*, *) 'wz ', minval(flowDomsd(nn, level, sps)%wz), &
            maxval(flowDomsd(nn, level, sps)%wz), &
            norm2(flowDomsd(nn, level, sps)%wz)
        write (*, *) 'qx ', minval(flowDomsd(nn, level, sps)%qx), &
            maxval(flowDomsd(nn, level, sps)%qx), &
            norm2(flowDomsd(nn, level, sps)%qx)
        write (*, *) 'qy ', minval(flowDomsd(nn, level, sps)%qy), &
            maxval(flowDomsd(nn, level, sps)%qy), &
            norm2(flowDomsd(nn, level, sps)%qy)
        write (*, *) 'qz ', minval(flowDomsd(nn, level, sps)%qz), &
            maxval(flowDomsd(nn, level, sps)%qz), &
            norm2(flowDomsd(nn, level, sps)%qz)

        write (*, *) 'bmti1 ', minval(flowDomsd(nn, level, sps)%bmti1), &
            maxval(flowDomsd(nn, level, sps)%bmti1), &
            norm2(flowDomsd(nn, level, sps)%bmti1)
        write (*, *) 'bmti2 ', minval(flowDomsd(nn, level, sps)%bmti2), &
            maxval(flowDomsd(nn, level, sps)%bmti2), &
            norm2(flowDomsd(nn, level, sps)%bmti2)
        write (*, *) 'bmtj1 ', minval(flowDomsd(nn, level, sps)%bmtj1), &
            maxval(flowDomsd(nn, level, sps)%bmtj1), &
            norm2(flowDomsd(nn, level, sps)%bmtj1)
        write (*, *) 'bmtj2 ', minval(flowDomsd(nn, level, sps)%bmtj2), &
            maxval(flowDomsd(nn, level, sps)%bmtj2), &
            norm2(flowDomsd(nn, level, sps)%bmtj2)
        write (*, *) 'bmtk1 ', minval(flowDomsd(nn, level, sps)%bmtk1), &
            maxval(flowDomsd(nn, level, sps)%bmtk1), &
            norm2(flowDomsd(nn, level, sps)%bmtk1)
        write (*, *) 'bmtk2 ', minval(flowDomsd(nn, level, sps)%bmtk2), &
            maxval(flowDomsd(nn, level, sps)%bmtk2), &
            norm2(flowDomsd(nn, level, sps)%bmtk2)
        write (*, *) 'bvti1 ', minval(flowDomsd(nn, level, sps)%bvti1), &
            maxval(flowDomsd(nn, level, sps)%bvti1), &
            norm2(flowDomsd(nn, level, sps)%bvti1)
        write (*, *) 'bvti2 ', minval(flowDomsd(nn, level, sps)%bvti2), &
            maxval(flowDomsd(nn, level, sps)%bvti2), &
            norm2(flowDomsd(nn, level, sps)%bvti2)
        write (*, *) 'bvtj1 ', minval(flowDomsd(nn, level, sps)%bvtj1), &
            maxval(flowDomsd(nn, level, sps)%bvtj1), &
            norm2(flowDomsd(nn, level, sps)%bvtj1)
        write (*, *) 'bvtj2 ', minval(flowDomsd(nn, level, sps)%bvtj2), &
            maxval(flowDomsd(nn, level, sps)%bvtj2), &
            norm2(flowDomsd(nn, level, sps)%bvtj2)
        write (*, *) 'bvtk1 ', minval(flowDomsd(nn, level, sps)%bvtk1), &
            maxval(flowDomsd(nn, level, sps)%bvtk1), &
            norm2(flowDomsd(nn, level, sps)%bvtk1)
        write (*, *) 'bvtk2 ', minval(flowDomsd(nn, level, sps)%bvtk2), &
            maxval(flowDomsd(nn, level, sps)%bvtk2), &
            norm2(flowDomsd(nn, level, sps)%bvtk2)

        bocoLoop: do mm = 1, flowDoms(nn, level, sps)%nBocos
            write (*, *) 'mm', mm, 'BCData(mm)%norm', minval(flowDomsd(nn, level, sps)%BCData(mm)%norm), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%norm), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%norm)
            write (*, *) 'mm', mm, 'bcData(mm)%rface ', minval(flowDomsd(nn, level, sps)%bcData(mm)%rface), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%rface), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%rface)
            write (*, *) 'mm', mm, 'bcData(mm)%Fv ', minval(flowDomsd(nn, level, sps)%bcData(mm)%Fv), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%Fv), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%Fv)
            write (*, *) 'mm', mm, 'bcData(mm)%Fp ', minval(flowDomsd(nn, level, sps)%bcData(mm)%Fp), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%Fp), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%Fp)
            write (*, *) 'mm', mm, 'bcData(mm)%Tv ', minval(flowDomsd(nn, level, sps)%bcData(mm)%Tv), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%Tv), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%Tv)
            write (*, *) 'mm', mm, 'bcData(mm)%Tp ', minval(flowDomsd(nn, level, sps)%bcData(mm)%Tp), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%Tp), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%Tp)
            write (*, *) 'mm', mm, 'bcData(mm)%area ', minval(flowDomsd(nn, level, sps)%bcData(mm)%area), &
                maxval(flowDomsd(nn, level, sps)%bcData(mm)%area), &
                norm2(flowDomsd(nn, level, sps)%bcData(mm)%area)
            write (*, *) 'mm', mm, 'BCData(mm)%uSlip ', minval(flowDomsd(nn, level, sps)%BCData(mm)%uSlip), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%uSlip), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%uSlip)
            write (*, *) 'mm', mm, 'BCData(mm)%TNS_Wall ', minval(flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall)
            write (*, *) 'mm', mm, 'BCData(mm)%ptInlet ', minval(flowDomsd(nn, level, sps)%BCData(mm)%ptInlet), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%ptInlet), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%ptInlet)
            write (*, *) 'mm', mm, 'BCData(mm)%htInlet ', minval(flowDomsd(nn, level, sps)%BCData(mm)%htInlet), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%htInlet), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%htInlet)
            write (*, *) 'mm', mm, 'BCData(mm)%ttInlet ', minval(flowDomsd(nn, level, sps)%BCData(mm)%ttInlet), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%ttInlet), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%ttInlet)
            write (*, *) 'mm', mm, 'BCData(mm)%turbInlet ', minval(flowDomsd(nn, level, sps)%BCData(mm)%turbInlet), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%turbInlet), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%turbInlet)
            write (*, *) 'mm', mm, 'BCData(mm)%ps ', minval(flowDomsd(nn, level, sps)%BCData(mm)%ps), &
                maxval(flowDomsd(nn, level, sps)%BCData(mm)%ps), &
                norm2(flowDomsd(nn, level, sps)%BCData(mm)%ps)

        end do bocoLoop

        viscbocoLoop: do mm = 1, flowDoms(nn, level, sps)%nViscBocos
            write (*, *) 'mm', mm, 'viscSubface(mm)%tau ', minval(flowDomsd(nn, level, sps)%viscSubface(mm)%tau), &
                maxval(flowDomsd(nn, level, sps)%viscSubface(mm)%tau), &
                norm2(flowDomsd(nn, level, sps)%viscSubface(mm)%tau)
            write (*, *) 'mm', mm, 'viscSubface(mm)%q ', minval(flowDomsd(nn, level, sps)%viscSubface(mm)%q), &
                maxval(flowDomsd(nn, level, sps)%viscSubface(mm)%q), &
                norm2(flowDomsd(nn, level, sps)%viscSubface(mm)%q)
        end do viscbocoLoop

        ! For overset, the weights may be active in the comm structure. We
        ! need to zero them before we can accumulate.
        if (oversetPresent) then
            ! Pointers to the overset comms to make it easier to read
            sends: do i = 1, commPatternOverset(level, sps)%nProcSend
                write (*, *) 'commPatternOverset(level, sps)%sendList(i)%interpd ', &
                    minval(commPatternOverset(level, sps)%sendList(i)%interpd), &
                    maxval(commPatternOverset(level, sps)%sendList(i)%interpd), &
                    norm2(commPatternOverset(level, sps)%sendList(i)%interpd)
            end do sends
           write (*, *) 'internalOverset(level, sps)%donorInterpd ', minval(internalOverset(level, sps)%donorInterpd), &
                maxval(internalOverset(level, sps)%donorInterpd), &
                norm2(internalOverset(level, sps)%donorInterpd)
        end if

        write (*, *) 'alphad ', alphad
        write (*, *) 'betad ', betad
        write (*, *) 'machd ', machd
        write (*, *) 'machGridd ', machGridd
        write (*, *) 'machCoefd ', machCoefd
        write (*, *) 'pinfdimd ', pinfdimd
        write (*, *) 'tinfdimd ', tinfdimd
        write (*, *) 'rhoinfdimd ', rhoinfdimd
        write (*, *) 'rgasdimd ', rgasdimd
        write (*, *) 'pointrefd ', pointrefd
        write (*, *) 'prefd ', prefd
        write (*, *) 'rhoRefd ', rhoRefd
        write (*, *) 'Trefd ', Trefd
        write (*, *) 'murefd ', murefd
        write (*, *) 'urefd ', urefd
        write (*, *) 'hrefd ', hrefd
        write (*, *) 'timerefd ', timerefd
        write (*, *) 'pinfd ', pinfd
        write (*, *) 'pinfCorrd ', pinfCorrd
        write (*, *) 'rhoinfd ', rhoinfd
        write (*, *) 'uinfd ', uinfd
        write (*, *) 'rgasd ', rgasd
        write (*, *) 'muinfd ', muinfd
        write (*, *) 'gammainfd ', gammainfd
        write (*, *) 'winfd ', winfd
        write (*, *) 'veldirfreestreamd ', veldirfreestreamd
        write (*, *) 'liftdirectiond ', liftdirectiond
        write (*, *) 'dragdirectiond ', dragdirectiond

        ! Zero all the reverse seeds in the dirichlet input arrays
        do iDom = 1, cgnsNDom
            do iBoco = 1, cgnsDoms(iDom)%nBocos
                if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)) then
                    do iData = 1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet)
                        if (associated(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)) then
                            do iDirichlet = 1, size(cgnsDoms(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays)
                                write (*, *) iDom, iBoco, iData, iDirichlet, 'dataArr(:) ' &
                                 , cgnsDomsd(iDom)%bocoInfo(iBoco)%dataSet(iData)%dirichletArrays(iDirichlet)%dataArr(:)
                            end do
                        end if
                    end do
                end if
            end do
        end do

        ! And the reverse seeds in the actuator zones
        do i = 1, nActuatorRegions
            write (*, *) 'actuatorRegionsd(i)%Force ', actuatorRegionsd(i)%force
            write (*, *) 'actuatorRegionsd(i)%Torque ', actuatorRegionsd(i)%torque
        end do

    end subroutine printADSeeds

#else

    subroutine computeMatrixFreeProductFwdCS(xvdot, extradot, wdot, bcDataValuesdot, &
                                             useSpatial, useState, famLists, &
                                             bcDataNames, bcDataValues, bcDataFamLists, bcVarsEmpty, &
                                             dwdot, funcsDot, fDot, &
                                             costSize, fSize, nTime, h_mag)

        ! This routine is used to debug master_d. It uses the forward seeds to set perturbations
        ! and then computes the value of the derivatives using forward complex step

        use constants
        use adjointvars
        use blockPointers, only: nDom
        use communication, only: adflow_comm_world
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputPhysics, only: pointRef, alpha, beta, equations, machCoef, &
                                mach, machGrid, rgasdim
        use iteration, only: currentLevel, groundLevel
        use flowVarRefState, only: pInfDim, rhoInfDim, TinfDim
        use blockPointers, only: nDom, il, jl, kl, wd, x, w, dw, dwd, nBocos, nViscBocos

        use adjointUtils, only: allocDerivativeValues, zeroADSeeds
        use masterRoutines, only: master
        use utils, only: isWallType, setPointers, setPointers_d, EChk
        use flowVarRefState, only: nw, nwf
        use wallDistanceData, only: xSurf, xSurfVec
        use wallDistance, only: updateXSurf

        implicit none

        ! Input Variables
        complex(kind=realType), dimension(:), intent(in) :: xvdot
        complex(kind=realType), dimension(:), intent(in) :: extradot
        complex(kind=realType), dimension(:), intent(in) :: wdot

        logical, intent(in) :: useSpatial, useState
        integer(kind=intType), dimension(:, :) :: famLists
        integer(kind=intType) :: costSize, fSize, nTime

        character, dimension(:, :), intent(in) :: bcDataNames
        real(kind=realType), dimension(:), intent(inout) :: bcDataValues, bcDataValuesDot
        integer(kind=intType), dimension(:, :) :: bcDataFamLists
        logical, intent(in) :: BCVarsEmpty

        ! step parameters
        real(kind=alwaysRealType), intent(in) :: h_mag ! step size for step

        ! Ouput Variables
        complex(kind=realType), dimension(size(wdot)), intent(out) :: dwDot
        complex(kind=realType), dimension(costSize, size(famLists, 1)), intent(out) :: funcsDot
        complex(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot

        ! Working Variables
        integer(kind=intType) :: nn, sps, level
        integer(kind=intType) :: ierr, mm, i, j, k, l, ii, jj, iRegion

        complex(kind=realType), dimension(costSize, size(famLists, 1)) :: funcs

        ! Input Arguments for master:
        complex(kind=realType), dimension(costSize, size(famLists, 1)) :: funcValues

        ! Working Variables
        complex(kind=realType), dimension(:, :, :), allocatable :: forces
        complex(kind=realType) :: h ! step size for Finite Difference

        ! note that h_mag does not have to be complex
        ! it is just the magnitude of the complex perturbation
        h = cmplx(0, h_mag)

        fSize = size(fDot, 2)
        allocate (forces(3, fSize, nTimeIntervalsSpectral))

        ! Need to trick the residual evalution to use coupled (mean flow and
        ! turbulent) together.
        level = 1
        currentLevel = level
        groundLevel = level

        ! Allocate the memory we need for derivatives if not done so
        ! already. Note this isn't deallocated until the adflow is
        ! destroyed.
        if (.not. derivVarsAllocated) then
            call allocDerivativeValues(level)
        end if

        ! Zero all AD seesd.
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call zeroADSeeds(nn, level, sps)
            end do
        end do

        ! Set the extra seeds now do the extra ones. Note that we are assuming the
        ! machNumber used for the coefficients follows the Mach number,
        ! not the grid mach number.

        alpha = alpha + h * extraDot(iAlpha)
        beta = beta + h * extraDot(iBeta)
        mach = mach + h * extraDot(iMach)
        machCoef = machCoef + h * extraDot(iMach)
        machGrid = machGrid + h * extraDot(iMachGrid)
        PinfDim = PinfDim + h * extraDot(iPressure)
        rhoinfDim = rhoinfDim + h * extraDot(iDensity)
        tinfdim = tinfdim + h * extraDot(iTemperature)
        pointref(1) = pointref(1) + h * extraDot(iPointRefX)
        pointref(2) = pointref(2) + h * extraDot(iPointRefY)
        pointref(3) = pointref(3) + h * extraDot(iPointRefZ)
        rgasdim = rgasdim + h * zero

        !  --------------------- apply the perturbations ----------------------------
        ! Set the provided w and x seeds:
        ii = 0
        jj = 0
        domainLoop1: do nn = 1, nDom
            spectalLoop1: do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do l = 1, 3
                                ii = ii + 1
                                x(i, j, k, l) = x(i, j, k, l) + xvdot(ii) * h
                            end do
                        end do
                    end do
                end do
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                jj = jj + 1
                                w(i, j, k, l) = w(i, j, k, l) + wDot(jj) * h
                            end do
                        end do
                    end do
                end do
            end do spectalLoop1
        end do domainLoop1

        bcDataValues = bcDataValues + bcDataValuesdot * h

        if (equations == RANSEquations) then
            call updateXSurf(level)
        end if

        ! ----------------------------- Run Master ---------------------------------
        ! Run the super-dee-duper master rotuine
        if (bcVarsEmpty) then
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces)
        else
            call master(useSpatial, &
                        famLists, funcValues, &
                        forces, &
                        bcDataNames, bcDataValues, bcDataFamLists)
        end if

        ! Copy out the residual derivative into the provided dwDot and remove the
        ! perturbation
        ii = 0
        jj = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers_d(nn, 1, sps)
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do l = 1, 3
                                ii = ii + 1
                                x(i, j, k, l) = x(i, j, k, l) - xvdot(ii) * h
                            end do
                        end do
                    end do
                end do
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                jj = jj + 1
                                w(i, j, k, l) = w(i, j, k, l) - wDot(jj) * h
                                dwd(i, j, k, l) = aimag(dw(i, j, k, l)) / aimag(h)
                                dwdot(jj) = dwd(i, j, k, l) ! copy values to output
                            end do
                        end do
                    end do
                end do

            end do
        end do

        alpha = alpha - h * extraDot(iAlpha)
        beta = beta - h * extraDot(iBeta)
        mach = mach - h * extraDot(iMach)
        machCoef = machCoef - h * extraDot(iMach)
        machGrid = machGrid - h * extraDot(iMachGrid)
        PinfDim = PinfDim - h * extraDot(iPressure)
        rhoinfDim = rhoinfDim - h * extraDot(iDensity)
        tinfdim = tinfdim - h * extraDot(iTemperature)
        pointref(1) = pointref(1) - h * extraDot(iPointRefX)
        pointref(2) = pointref(2) - h * extraDot(iPointRefY)
        pointref(3) = pointref(3) - h * extraDot(iPointRefZ)
        rgasdim = rgasdim - h * zero

        bcDataValues = bcDataValues - bcDataValuesdot * h
        if (equations == RANSEquations) then
            call updateXSurf(level)
        end if

        fDot = aimag(forces) / aimag(h)
        funcsDot = aimag(funcValues) / aimag(h)

    end subroutine computeMatrixFreeProductFwdCS

#endif

end module adjointDebug

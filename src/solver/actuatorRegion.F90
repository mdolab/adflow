module actuatorRegion

    use constants
    use communication, only: commType, internalCommType
    use actuatorRegionData
    implicit none

contains
    subroutine computeCellSpatialMetrics(i, j, k, centerPoint, thrustVector, distance2plane, distance2axis, tangent)
        use constants
        use blockPointers, only: x
        implicit none

        ! Inputs
        real(kind=realType), dimension(3), intent(in) :: centerPoint, thrustVector
        integer(kind=intType), intent(in) :: i, j, k

        ! Outputs
        real(kind=realType), intent(out) :: distance2axis, distance2plane
        real(kind=realType), dimension(3), intent(out) :: tangent

        ! Working
        real(kind=realType)  :: thrustVectorNorm, dotProduct, tangentNorm
        real(kind=realType), dimension(3) ::  xCen, distance2center, rawTangent

        ! Compute the cell center
        xCen = eighth * (x(i - 1, j - 1, k - 1, :) + x(i, j - 1, k - 1, :) &
                         + x(i - 1, j, k - 1, :) + x(i, j, k - 1, :) &
                         + x(i - 1, j - 1, k, :) + x(i, j - 1, k, :) &
                         + x(i - 1, j, k, :) + x(i, j, k, :))

        thrustVectorNorm = sqrt(thrustVector(1)**2 + thrustVector(2)**2 + thrustVector(3)**2)
        distance2center = xCen - centerPoint

        ! compute distance to plane
        dotProduct = thrustVector(1)*distance2center(1) + thrustVector(2)*distance2center(2) + thrustVector(3)*distance2center(3)
        distance2plane = dotProduct / thrustVectorNorm

        ! compute distance to axis
        rawTangent(1) = distance2center(2)*thrustVector(3) - distance2center(3)*thrustVector(2)
        rawTangent(2) = distance2center(3)*thrustVector(1) - distance2center(1)*thrustVector(3)
        rawTangent(3) = distance2center(1)*thrustVector(2) - distance2center(2)*thrustVector(1)

        tangentNorm = sqrt(rawTangent(1)**2 + rawTangent(2)**2 + rawTangent(3)**2)

        distance2axis = tangentNorm / thrustVectorNorm

        ! compute tangential vector
        tangent = rawTangent / tangentNorm


    end subroutine computeCellSpatialMetrics


    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE
    subroutine computeInitialSpatialMetrics(centerPoint, thrustVector, n, distance2plane, distance2axis, tangent)

        use constants
        use blockPointers, only: x, il, jl, kl, nDom, iBlank
        use utils, only: setPointers
        implicit none

        ! Input variables
        real(kind=realType), dimension(3), intent(in) :: centerPoint, thrustVector
        integer(kind=intType), intent(in) :: n
        real(kind=realType), dimension(n), intent(out) :: distance2plane, distance2axis
        real(kind=realType), dimension(3, n), intent(out) :: tangent


        ! Working variables
        integer(kind=intType) :: ii, i, j, k, nn, iDim, level, sps
        real(kind=realType), dimension(3) ::  xCen

        ! Loop through all the cells and compute the distance of each cell center from the thrust axis (centerPoint + thrustVector)
        ! also compute the distance from the plane formed by the centerPoint and the normal trhustVector
        distance2axis = large
        distance2plane = large

        sps = 1
        level = 1

        ii = 0
        do nn = 1, nDom
            call setPointers(nn, level, sps)
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ii = ii + 1

                        ! cycle if it is not a compute cell
                        if (iblank(i, j, k) /= 1) then
                            cycle
                        end if

                        call computeCellSpatialMetrics(i, j, k, &
                            centerPoint, thrustVector, &
                            distance2plane(ii), distance2axis(ii), tangent(:, ii))
                    end do
                end do
            end do
        end do

    end subroutine computeInitialSpatialMetrics


    subroutine addActuatorRegion(flag, n, famName, famID, relaxStart, relaxEnd, iRegion, nLocalCells)
        ! Add a user-supplied integration surface.

        use communication, only: myID, adflow_comm_world
        use constants
        use blockPointers, only: x, il, jl, kl, nDom, iBlank, vol
        use adjointVars, only: nCellsLocal
        use utils, only: setPointers, EChk
        implicit none

        ! Input variables
        integer(kind=intType), intent(in) :: famID
        character(len=*), intent(in) :: famName
        real(kind=realType), intent(in) :: relaxStart, relaxEnd

        integer(kind=intType), intent(in) :: n
        integer(kind=intType), dimension(n), intent(in) :: flag

        ! Output Variables
        integer(kind=intType), intent(out) :: iRegion
        integer(kind=intType), intent(out) :: nLocalCells

        ! Working variables
        integer(kind=intType) :: i, j, k, nn, sps, level, ii, iii, ierr
        type(actuatorRegionType), pointer :: region
        integer(kind=intType), dimension(:, :), pointer :: tmp

        nActuatorRegions = nActuatorRegions + 1
        if (nActuatorRegions > nActuatorRegionsMax) then
            print *, "Error: Exceeded the maximum number of actuatorDiskRegions.  "&
                 &"Increase nActuatorDiskRegionsMax"
            stop
        end if

        ! Save the input information
        region => actuatorRegions(nActuatorRegions)
        iRegion = nActuatorRegions
        region%famName = famName
        region%famID = famID
        region%relaxStart = relaxStart
        region%relaxEnd = relaxEnd


        allocate (region%blkPtr(0:nDom))
        region%blkPtr(0) = 0
        allocate (region%cellIDs(3, nCellsLocal(1)))

        ! Now search for all the coordinate. Note that We have explictly
        ! set sps to 1 becuase it is only implemented for single grid.
        sps = 1
        level = 1

        ii = 0
        do nn = 1, nDom
            call setPointers(nn, level, sps)
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ii = ii + 1

                        ! cycle if it is not a compute cell
                        if (iblank(i, j, k) /= 1) then
                            cycle
                        end if

                        ! cycle if cell was not tagged
                        if (flag(ii) /= 1) then
                            cycle
                        end if

                        region%nCellIDs = region%nCellIDs + 1
                        region%cellIDs(:, region%nCellIDs) = (/i, j, k/)
                    end do
                end do
            end do
            ! Since we're doing all the blocks in order, simply store the
            ! current counter into blkPtr which gives up the range of cells
            ! we have found on this block
            region%blkPtr(nn) = region%nCellIDs

        end do


        ! Resize the cellIDs to the correct size now that we know the
        ! correct exact number.
        tmp => region%cellIDs
        allocate (region%cellIDs(3, region%nCellIDs))
        region%cellIDs = tmp(:, 1:region%nCellIDs)
        deallocate (tmp)
        nLocalCells = region%nCellIDs


        ! Allocate variables for the source terms
        allocate (region%force(3, region%nCellIDs))
        allocate (region%heat(region%nCellIDs))

    end subroutine addActuatorRegion

    subroutine computeSpatialMetrics(iRegion, nLocalCells, centerPoint, thrustVector, &
            distance2plane, distance2axis, tangent, volume)
        use constants
        use blockPointers, only: nDom, vol
        use utils, only: setPointers
        implicit none

        ! Inputs
        real(kind=realType), dimension(3), intent(in) :: centerPoint, thrustVector
        integer(kind=intType), intent(in) :: nLocalCells
        integer(kind=intType), intent(in) :: iRegion

        ! Outputs
        real(kind=realType), intent(out), dimension(nLocalCells) :: distance2axis, distance2plane, volume
        real(kind=realType), dimension(3, nLocalCells), intent(out) :: tangent

        ! Working
        type(actuatorRegionType), pointer :: region
        integer(kind=intType) :: sps, level, nn, i, j, k, iii

        region => actuatorRegions(iRegion)


        sps = 1
        level = 1
        do nn = 1, nDom
            call setPointers(nn, level, sps)
            do iii = actuatorRegions(iRegion)%blkPtr(nn - 1) + 1, actuatorRegions(iRegion)%blkPtr(nn)
                i = actuatorRegions(iRegion)%cellIDs(1, iii)
                j = actuatorRegions(iRegion)%cellIDs(2, iii)
                k = actuatorRegions(iRegion)%cellIDs(3, iii)


                call computeCellSpatialMetrics(i, j, k, &
                    centerPoint, thrustVector, &
                    distance2plane(iii), distance2axis(iii), tangent(:, iii))

                volume(iii) = vol(i, j, k)
            end do
        end do

    end subroutine computeSpatialMetrics

    subroutine populateBCValues(iRegion, nLocalCells, force, heat)
        use constants
        use blockPointers, only: nDom
        use utils, only: setPointers
        implicit none

        ! Inputs
        integer(kind=intType), intent(in) :: iRegion
        integer(kind=intType), intent(in) :: nLocalCells
        real(kind=realType), dimension(3, nLocalCells), intent(in) :: force
        real(kind=realType), dimension(nLocalCells), intent(in) :: heat


        ! Working
        type(actuatorRegionType), pointer :: region
        integer(kind=intType) :: sps, level, nn, i, j, k, iii


        region => actuatorRegions(iRegion)


        sps = 1
        level = 1
        do nn = 1, nDom
            call setPointers(nn, level, sps)
            do iii = actuatorRegions(iRegion)%blkPtr(nn - 1) + 1, actuatorRegions(iRegion)%blkPtr(nn)
                actuatorRegions(iRegion)%force(:, iii) = force(:, iii)
                actuatorRegions(iRegion)%heat(iii) = heat(iii)
            end do
        end do

    end subroutine populateBCValues


    subroutine writeActuatorRegions(fileName)

        ! This a (mostly) debug routine that is used to verify to the user
        ! the that the cells that the user thinks should be specified as
        ! being inside the actuator region actually are. We will dump a
        ! hex unstructured ascii tecplot file with all the zones we
        ! found. We won't be super concerned about efficiency here.

        use constants
        use utils, only: EChk, pointReduce, setPointers
        use communication, only: myID, adflow_comm_world, nProc
        use blockPointers, only: x, nDom
        use commonFormats, only: sci12
        implicit none

        ! Input
        character(len=*) :: fileName

        ! Working
        integer(kind=intType) :: iRegion, nn, i, j, k, ii, jj, kk, iii, kkk, iDim
        integer(kind=intType) :: level, sps, iProc, ierr, totalCount, offset, nUnique
        integer(kind=intType), dimension(:), allocatable :: sizesProc, cumSizesProc
        real(kind=realType), dimension(:), allocatable :: pts, allPts
        real(kind=realType), dimension(:, :), allocatable :: tmp, uniquePts
        real(kind=realType), parameter :: tol = 1e-8
        integer(kind=intType), dimension(:), allocatable :: conn, allConn, link
        character(80) :: zoneName
        type(actuatorRegionType), pointer :: region

        ! Before we start the main region loop the root procesoor has to
        ! open up the tecplot file and write the header

        if (myid == 0) then
            open (unit=101, file=trim(fileName), form='formatted')
            write (101, *) 'TITLE = "Actuator Regions"'
            write (101, *) 'Variables = "CoordinateX", "CoordinateY", "CoordinateZ"'
        end if

        ! Region Loop
        regionLoop: do iRegion = 1, nActuatorRegions

            ! Only for the finest grid level.
            level = 1
            sps = 1

            ! Do an allgather with the number of actuator cells on each
            ! processor so that everyone knows the sizes and can compute the offsets,
            region => actuatorRegions(iRegion)

            allocate (sizesProc(nProc), cumSizesProc(0:nProc))

            call mpi_allgather(region%nCellIDs, 1, adflow_integer, sizesProc, 1, &
                               adflow_integer, adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            cumSizesProc(0) = 0
            do iProc = 1, nProc
                cumSizesProc(iProc) = cumSizesProc(iProc - 1) + sizesProc(iProc)
            end do

            ! Fill up our own nodes/conn with the nodes we have here.
            allocate (conn(8 * region%nCellIDs), pts(24 * region%nCellIDs))

            kkk = 0
            do nn = 1, nDom
                call setPointers(nn, level, sps)
                ! Loop over the ranges for this block
                do iii = region%blkPtr(nn - 1) + 1, region%blkPtr(nn)

                    ! Carful with the conn values! They need to be in counter clock wise ordering!
                    offset = (iii - 1) * 8 + cumSizesProc(myID) * 8
                    conn((iii - 1) * 8 + 1) = 1 + offset
                    conn((iii - 1) * 8 + 2) = 2 + offset
                    conn((iii - 1) * 8 + 3) = 4 + offset
                    conn((iii - 1) * 8 + 4) = 3 + offset
                    conn((iii - 1) * 8 + 5) = 5 + offset
                    conn((iii - 1) * 8 + 6) = 6 + offset
                    conn((iii - 1) * 8 + 7) = 8 + offset
                    conn((iii - 1) * 8 + 8) = 7 + offset

                    ! Add in the 24 values for the nodal coordinates in coordinate
                    ! ordering. Do all the coordinates interlaced
                    do kk = -1, 0
                        do jj = -1, 0
                            do ii = -1, 0
                                do iDim = 1, 3
                                    i = region%cellIDs(1, iii)
                                    j = region%cellIDs(2, iii)
                                    k = region%cellIDs(3, iii)
                                    kkk = kkk + 1
                                    pts(kkk) = x(i + ii, j + jj, k + kk, iDim)
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            ! Now that we've filled up our array, we can allocate the total
            ! space we need on the root proc and it
            if (myid == 0) then
                totalCount = sum(sizesProc)
                allocate (allConn(8 * totalCount), allPts(24 * totalCount))
            end if

            ! Perform the two gatherV's
            call mpi_gatherV(pts, region%nCellIDs * 24, adflow_real, &
                             allPts, 24 * sizesProc, 24 * cumSizesProc, adflow_real, &
                             0, adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            call mpi_gatherV(conn, region%nCellIDs * 8, adflow_integer, &
                             allConn, 8 * sizesProc, 8 * cumSizesProc, adflow_integer, &
                             0, adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            ! We can deallocate all the per-proc memory now
            deallocate (sizesProc, cumSizesProc, pts, conn)

            if (myid == 0) then

                ! Now the poor root processor dumps everything out to a
                ! file. To help cut down on the already bloated file size,
                ! we'll point reduce it which will help tecplot display them
                ! better as well.

                allocate (tmp(3, totalCount * 8))
                do i = 1, totalCount * 8
                    do iDim = 1, 3
                        tmp(iDim, i) = allPts((i - 1) * 3 + iDim)
                    end do
                end do
                deallocate (allPts)
                allocate (uniquePts(3, totalCount * 8), link(totalCount * 8))

                ! Get unique set of nodes.
                call pointReduce(tmp, totalCount * 8, tol, uniquePts, link, nUnique)

                write (zoneName, "(a,a,a)") 'Zone T="', trim(region%famName), ' Region"'
                write (101, *) trim(zoneName)
                write (101, *) "Nodes = ", nUnique, " Elements= ", totalCount, " ZONETYPE=FEBRICK"
                write (101, *) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"

                ! Write all the coordinates...this is horrendously slow...
                do iDim = 1, 3
                    do i = 1, nUnique
                        write (101, sci12) uniquePts(iDim, i)
                    end do
                end do

                ! Write out the connectivity
                do i = 1, totalCount
                    do j = 1, 8
                        write (101, "(I8)", advance='no') link(allConn((i - 1) * 8 + j))
                    end do
                    write (101, "(1x)")
                end do

                ! Ditch the memory only allocated on this proc
                deallocate (allConn, link, tmp, uniquePts)
            end if
        end do regionLoop

        ! Close the output file on the root proc
        if (myid == 0) then
            close (101)
        end if
    end subroutine writeActuatorRegions

    subroutine integrateActuatorRegions(localValues, famList, sps)
        !--------------------------------------------------------------
        ! Manual Differentiation Warning: Modifying this routine requires
        ! modifying the hand-written forward and reverse routines.
        ! --------------------------------------------------------------

        ! Perform volume integrals over the actuator region.
        use constants
        use blockPointers, only: vol, dw, w, nDom
        use flowVarRefState, only: Pref, uRef
        use utils, only: setPointers
        use sorting, only: famInList
        use residuals, only: sourceTerms_block
        use actuatorRegionData
        use communication
        ! Input/output Variables
        real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
        integer(kind=intType), dimension(:), intent(in) :: famList
        integer(kind=intType), intent(in) :: sps

        ! Working
        integer(kind=intType) :: nn, iRegion
        real(kind=realType) :: PLocal, PLocald
        ! Zero the accumulation variable. We comptue flow power across
        ! 'all' actuaor zones. The famInclude is used to section out which
        ! one we want.
        PLocal = zero

        domainLoop: do nn = 1, nDom
            call setPointers(nn, 1, sps)

            ! Loop over each region
            regionLoop: do iRegion = 1, nActuatorRegions

                ! Check if this needs to be included:
                famInclude: if (famInList(actuatorRegions(iRegion)%famID, famList)) then

                    ! If so, call the regular sourceTerms_block routine
                    call sourceTerms_block(nn, .False., iRegion, pLocal)

                end if famInclude
            end do regionLoop
        end do domainLoop

        ! Add in the contribution from this processor.
        localValues(iPower) = localValues(iPower) + PLocal

    end subroutine integrateActuatorRegions
#ifndef USE_COMPLEX
    subroutine integrateActuatorRegions_d(localValues, localValuesd, famList, sps)
        !--------------------------------------------------------------
        ! Manual Differentiation Warning: Modifying this routine requires
        ! modifying the hand-written forward and reverse routines.
        ! --------------------------------------------------------------

        ! Perform volume integrals over the actuator region.
        use constants
        use blockPointers, only: vol, dw, w, nDom
        use flowVarRefState, only: Pref, uRef
        use utils, only: setPointers_d
        use sorting, only: famInList
        use actuatorRegionData
        use residuals_d, only: sourceTerms_block_d

        ! Input/output Variables
        real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
        integer(kind=intType), dimension(:), intent(in) :: famList
        integer(kind=intType), intent(in) :: sps

        ! Working
        integer(kind=intType) :: nn, iRegion
        real(kind=realType) :: PLocal, PLocald

        ! Zero the accumulation variable. We comptue flow power across
        ! 'all' actuaor zones. The famInclude is used to section out which
        ! one we want.
        PLocal = zero
        PLocald = zero

        domainLoop: do nn = 1, nDom
            call setPointers_d(nn, 1, sps)

            ! Loop over each region
            regionLoop: do iRegion = 1, nActuatorRegions

                ! Check if this needs to be included:
                famInclude: if (famInList(actuatorRegions(iRegion)%famID, famList)) then

                    ! If so, call the regular sourceTerms_block routine
                    call sourceTerms_block_d(nn, .False., iRegion, pLocal, pLocald)

                end if famInclude
            end do regionLoop
        end do domainLoop

        ! Add in the contribution from this processor.
        localValues(iPower) = localValues(iPower) + PLocal
        localValuesd(iPower) = localValuesd(iPower) + PLocald

    end subroutine integrateActuatorRegions_d

    subroutine integrateActuatorRegions_b(localValues, localValuesd, famList, sps)
        !--------------------------------------------------------------
        ! Manual Differentiation Warning: Modifying this routine requires
        ! modifying the hand-written forward and reverse routines.
        ! --------------------------------------------------------------

        ! Perform volume integrals over the actuator region.
        use constants
        use blockPointers, only: vol, dw, w, nDom
        use flowVarRefState, only: Pref, uRef
        use utils, only: setPointers_b
        use sorting, only: famInList
        use actuatorRegionData
        use residuals_b, only: sourceTerms_block_b

        ! Input/output Variables
        real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
        integer(kind=intType), dimension(:), intent(in) :: famList
        integer(kind=intType), intent(in) :: sps

        ! Working
        integer(kind=intType) :: nn, iRegion
        real(kind=realType) :: PLocal, PLocald

        ! Zero the accumulation variable. We comptue flow power across
        ! 'all' actuaor zones. The famInclude is used to section out which
        ! one we want.
        PLocal = zero
        PLocald = zero

        ! Pull out the seed
        PLocald = localValuesd(iPower)

        domainLoop: do nn = 1, nDom
            call setPointers_b(nn, 1, sps)

            ! Loop over each region
            regionLoop: do iRegion = 1, nActuatorRegions

                ! Check if this needs to be included:
                famInclude: if (famInList(actuatorRegions(iRegion)%famID, famList)) then

                    ! If so, call the regular sourceTerms_block routine
                    call sourceTerms_block_b(nn, .False., iRegion, pLocal, pLocald)

                end if famInclude
            end do regionLoop
        end do domainLoop
    end subroutine integrateActuatorRegions_b

#endif
#endif
end module actuatorRegion

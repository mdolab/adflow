module actuatorRegion

    use constants
    use communication, only: commType, internalCommType
    use actuatorRegionData
    implicit none

contains

    subroutine computeActuatorRegionVolume(nn, iRegion)
        use blockPointers, only: nDom, vol
        implicit none

        ! Inputs
        integer(kind=intType), intent(in) :: nn, iRegion

        ! Working
        integer(kind=intType) :: iii
        integer(kind=intType) :: i, j, k

        ! Loop over the region for this block
        do iii = actuatorRegions(iRegion)%blkPtr(nn - 1) + 1, actuatorRegions(iRegion)%blkPtr(nn)
            i = actuatorRegions(iRegion)%cellIDs(1, iii)
            j = actuatorRegions(iRegion)%cellIDs(2, iii)
            k = actuatorRegions(iRegion)%cellIDs(3, iii)

            ! Sum the volume of each cell within the region on this proc
            actuatorRegions(iRegion)%volLocal = actuatorRegions(iRegion)%volLocal + vol(i, j, k)
        end do

    end subroutine computeActuatorRegionVolume

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE
    subroutine addActuatorRegion(pts, conn, axis1, axis2, famName, famID, &
                                 thrust, torque, heat, relaxStart, relaxEnd, nPts, nConn)
        ! Add a user-supplied integration surface.

        use communication, only: myID, adflow_comm_world
        use constants
        use adtBuild, only: buildSerialQuad, destroySerialQuad
        use adtLocalSearch, only: minDistanceTreeSearchSinglePoint
        use ADTUtils, only: stack
        use ADTData
        use blockPointers, only: x, il, jl, kl, nDom, iBlank, vol
        use adjointVars, only: nCellsLocal
        use utils, only: setPointers, EChk
        implicit none

        ! Input variables
        integer(kind=intType), intent(in) :: nPts, nConn, famID
        real(kind=realType), dimension(3, nPts), intent(in), target :: pts
        integer(kind=intType), dimension(4, nConn), intent(in), target :: conn
        real(kind=realType), intent(in), dimension(3) :: axis1, axis2
        character(len=*) :: famName
        real(kind=realType) :: thrust, torque, heat, relaxStart, relaxEnd

        ! Working variables
        integer(kind=intType) :: i, j, k, nn, iDim, cellID, intInfo(3), sps, level, iii, ierr
        real(kind=realType) :: dStar, frac, volLocal
        type(actuatorRegionType), pointer :: region
        real(kind=realType), dimension(3) :: minX, maxX, v1, v2, v3, xCen, axisVec
        type(adtType) :: ADT
        real(kind=realType) :: axisVecNorm
        real(kind=realType), dimension(:, :), allocatable :: norm
        integer(kind=intType), dimension(:), allocatable :: normCount
        integer(kind=intType), dimension(:, :), pointer :: tmp

        ! ADT Type required data
        integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew
        type(adtBBoxTargetType), dimension(:), pointer :: BB
        real(kind=realType) :: coor(4), uvw(5)
        real(kind=realType) :: dummy(3, 2)

        nActuatorRegions = nActuatorRegions + 1
        if (nActuatorRegions > nActuatorRegionsMax) then
            print *, "Error: Exceeded the maximum number of actuatorDiskRegions.  "&
                 &"Increase nActuatorDiskRegionsMax"
            stop
        end if

        ! Save the input information
        region => actuatorRegions(nActuatorRegions)
        region%famName = famName
        region%famID = famID
        region%torque = torque
        region%heat = heat
        region%relaxStart = relaxStart
        region%relaxEnd = relaxEnd
        ! We use the axis to define the direction of F. Since we are
        ! dealing with rotating machinary, it is pretty good approximation
        ! to assume that the thrust is going to be in the direction of the
        ! axis.
        axisVec = axis2 - axis1
        axisVecNorm = sqrt((axisVec(1)**2 + axisvec(2)**2 + axisVec(3)**2))
        if (axisVecNorm < 1e-12) then
            print *, "Error: Axis cannot be determined by the supplied points. They are too close"
            stop
        end if

        axisVec = axisVec / axisVecNorm

        region%force = axisVec * thrust
        region%axisVec = axisVec

        allocate (region%blkPtr(0:nDom))
        region%blkPtr(0) = 0

        ! Next thing we need to do is to figure out if any of our cells
        ! are inside the actuator disk region. If so we will save them in
        ! the actuatorRegionType data structure

        ! Since this is effectively a wall-distance calc it gets super
        ! costly for the points far away. Luckly, we can do a fairly
        ! simple shortcut: Just compute the bounding box of the region and
        ! use that as the "already found" distance in the cloest point
        ! search. This will eliminate all the points further away
        ! immediately and this should be sufficiently fast.

        ! So...compute that bounding box:
        do iDim = 1, 3
            minX(iDim) = minval(pts(iDim, :))
            maxX(iDim) = maxval(pts(iDim, :))
        end do

        ! Get the max distance. This should be quite conservative.
        dStar = (maxX(1) - minx(1))**2 + (maxX(2) - minX(2))**2 + (maxX(3) - minX(3))**2

        ! Now build the tree.
        call buildSerialQuad(size(conn, 2), size(pts, 2), pts, conn, ADT)

        ! Compute the (averaged) unique nodal vectors:
        allocate (norm(3, size(pts, 2)), normCount(size(pts, 2)))

        norm = zero
        normCount = 0

        do i = 1, size(conn, 2)

            ! Compute cross product normal and normalize
            v1 = pts(:, conn(3, i)) - pts(:, conn(1, i))
            v2 = pts(:, conn(4, i)) - pts(:, conn(2, i))

            v3(1) = (v1(2) * v2(3) - v1(3) * v2(2))
            v3(2) = (v1(3) * v2(1) - v1(1) * v2(3))
            v3(3) = (v1(1) * v2(2) - v1(2) * v2(1))
            v3 = v3 / sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)

            ! Add to each of the four pts and increment the number added
            do j = 1, 4
                norm(:, conn(j, i)) = norm(:, conn(j, i)) + v3
                normCount(conn(j, i)) = normCount(conn(j, i)) + 1
            end do
        end do

        ! Now just divide by the norm count
        do i = 1, size(pts, 2)
            norm(:, i) = norm(:, i) / normCount(i)
        end do

        ! Norm count is no longer needed
        deallocate (normCount)

        ! Allocate the extra data the tree search requires.
        allocate (stack(100), BB(20), frontLeaves(25), frontLeavesNew(25))

        ! Allocate sufficient space for the maximum possible number of cellIDs
        allocate (region%cellIDs(3, nCellsLocal(1)))

        ! Now search for all the coordinate. Note that We have explictly
        ! set sps to 1 becuase it is only implemented for single grid.
        sps = 1
        level = 1

        do nn = 1, nDom
            call setPointers(nn, level, sps)
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ! Only check real cells
                        if (iblank(i, j, k) == 1) then
                            ! Compute the cell center
                            xCen = eighth * (x(i - 1, j - 1, k - 1, :) + x(i, j - 1, k - 1, :) &
                                             + x(i - 1, j, k - 1, :) + x(i, j, k - 1, :) &
                                             + x(i - 1, j - 1, k, :) + x(i, j - 1, k, :) &
                                             + x(i - 1, j, k, :) + x(i, j, k, :))

                            ! The current point to search for and continually
                            ! reset the "closest point already found" variable.
                            coor(1:3) = xCen
                            coor(4) = dStar
                            intInfo(3) = 0
                            call minDistancetreeSearchSinglePoint(ADT, coor, intInfo, &
                                                                  uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
                            cellID = intInfo(3)
                            if (cellID > 0) then
                                ! Now check if this was successful or now:
                                if (checkInside()) then
                                    ! Whoohoo! We are inside the region. Add this cell
                                    ! to the list.
                                    region%nCellIDs = region%nCellIDs + 1
                                    region%cellIDs(:, region%nCellIDs) = (/i, j, k/)
                                end if
                            end if
                        end if
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

        ! Now go back and generate the total volume of the the cells we've flagged
        volLocal = zero

        do nn = 1, nDom
            call setPointers(nn, level, sps)

            ! Loop over the region for this block
            do iii = region%blkPtr(nn - 1) + 1, region%blkPtr(nn)
                i = region%cellIDs(1, iii)
                j = region%cellIDs(2, iii)
                k = region%cellIDs(3, iii)
                volLocal = volLocal + vol(i, j, k)
            end do
        end do

        call mpi_allreduce(volLocal, region%volume, 1, adflow_real, &
                           MPI_SUM, adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Final memory cleanup
        deallocate (stack, norm, frontLeaves, frontLeavesNew, BB)
        call destroySerialQuad(ADT)
    contains

        function checkInside()

            implicit none
            logical :: checkInside
            integer(kind=intType) :: jj
            real(kind=realType) :: shp(4), xp(3), normal(3), v1(3), dp

            ! bi-linear shape functions (CCW ordering)
            shp(1) = (one - uvw(1)) * (one - uvw(2))
            shp(2) = (uvw(1)) * (one - uvw(2))
            shp(3) = (uvw(1)) * (uvw(2))
            shp(4) = (one - uvw(1)) * (uvw(2))

            xp = zero
            normal = zero
            do jj = 1, 4
                xp = xp + shp(jj) * pts(:, conn(jj, cellID))
                normal = normal + shp(jj) * norm(:, conn(jj, cellID))
            end do

            ! Compute the dot product of normal with cell center
            ! (stored in coor) with the point on the surface.
            v1 = coor(1:3) - xp
            dp = normal(1) * v1(1) + normal(2) * v1(2) + normal(3) * v1(3)

            if (dp < zero) then
                checkInside = .True.
            else
                checkInside = .False.
            end if
        end function checkInside
    end subroutine addActuatorRegion

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

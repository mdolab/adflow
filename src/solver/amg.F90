module amg

    use constants
    use utils, only: EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Structure used for storing the interpolation indices
    type arr1int4
        integer, dimension(:), allocatable :: arr
    end type arr1int4

    ! Structure used for storing the interpolation indices
    type arr3int4
        integer, dimension(:, :, :), allocatable :: arr
    end type arr3int4

    ! Structure used for storing the interpolation indices
    type arr4int4
        integer, dimension(:, :, :, :), allocatable :: arr
    end type arr4int4

    ! The number of levels
    integer(kind=intType) amgLevels

    ! The number of outer iterations
    integer(kind=intType) amgOuterIts

    ! The number of smoothing iterations
    integer(kind=intType) amgNSmooth

    ! ASM overlap for the solver/smoother
    integer(kind=intType) :: amgASMOverlap

    ! ILU fill for the solver/smoother
    integer(kind=intType) :: amgFillLevel

    ! Ordering
    character(len=maxStringLen) :: amgMatrixOrdering

    ! Arrays for matrices/vectors/solvers on each level.
    Mat, dimension(:), allocatable :: A
    KSP, dimension(:), allocatable :: kspLevels
    Vec, dimension(:), allocatable :: res, rhs, sol, sol2
    PC shellPC

    Mat, pointer :: fineMat

    ! The interpolation arrays
    type(arr1int4), dimension(:), allocatable :: interps

    ! The coarse level indices
    type(arr3int4), dimension(:, :), allocatable, target :: coarseIndices
    type(arr4int4), dimension(:, :), allocatable, target :: coarseOversetIndices

    logical :: amgSetup = .False.
    integer :: bs
contains

    subroutine setupAMG(inputMat, nCell, blockSize, levels, nSmooth)

        use blockPointers
        use communication
        use utils, only: setPointers
        use haloExchange, only: whalo1to1intgeneric
        use adjointVars, only: nCellsLocal
        integer(kind=intType), intent(in) :: nCell, blockSize, levels, nSmooth
        Mat, target :: inputMat

        VecScatter :: scat
        Vec :: recvVec, indexVec
        IS :: IS1

        integer(kind=intTYpe) :: nnx, nny, nnz, l, nn, i, j, k, ii, jj, kk, lvl, n, m
        integer(kind=intType) :: idim, count, cursize, ierr, coarseIndex
        integer(kind=intType) :: level, nVar, sps, nextLvl, ncoarse
        integer(kind=intType), dimension(:, :, :), allocatable :: sizes
        integer(kind=intType), dimension(:), allocatable :: offsets, procStarts, nnzOn, nnzOff
        real(kind=realType), dimension(:), pointer :: indPtr
        integer(kind=intType), dimension(:), allocatable :: indicesToGet

        if (amgSetup) then
            return
        end if

        ! Set some module variables
        amgLevels = levels
        amgNSmooth = nSmooth

        ! Block size here refers to the number of states in each cell
        bs = blockSize

        ! Set the pointer to the fine level the AMG will be working on
        fineMat => inputMat

        ! Allocate the list of the mats/vects/ksp
        allocate ( &
            A(2:amgLevels), &
            kspLevels(1:amgLevels), &
            res(1:amgLevels), &
            rhs(1:amgLevels), &
            sol(1:amgLevels), &
            sol2(1:amgLevels))

        ! First allocate the coarse level indices.
        allocate (coarseIndices(nDom, 1:amgLevels - 1))
        allocate (coarseOversetIndices(nDom, 1:amgLevels - 1))
        allocate (sizes(3, nDom, 1:amgLevels))
        allocate (offsets(0:nProc), procStarts(1:amgLevels))
        allocate (interps(1:amgLevels - 1))
        do lvl = 1, amgLevels

            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                if (lvl == 1) then

                    ! Set the sizes for the finest level.
                    ! This is the number of owned cells.
                    sizes(1, nn, 1) = nx
                    sizes(2, nn, 1) = ny
                    sizes(3, nn, 1) = nz
                end if
            end do

            if (lvl > 1) then
                do nn = 1, nDom
                    do iDim = 1, 3
                        curSize = sizes(iDim, nn, lvl - 1)
                        if (curSize == 1) then
                            sizes(iDim, nn, lvl) = curSize
                        else if (mod(curSize, 2) == 0) then
                            ! Evenly divides
                            sizes(iDim, nn, lvl) = curSize / 2
                        else if (mod(curSize, 2) == 1) then
                            ! Odd, so have an extra
                            sizes(iDim, nn, lvl) = (curSize - 1) / 2 + 1
                        end if
                    end do
                end do
            end if

            ! Get the number of cells on this proc:
            i = 0
            do nn = 1, nDom
                i = i + sizes(1, nn, lvl) * sizes(2, nn, lvl) * sizes(3, nn, lvl)
            end do

            call MPI_Allgather(i, 1, mpi_integer, offsets(1:nProc), 1, mpi_integer, &
                               adflow_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Prefix sum
            offsets(0) = 0
            do i = 1, nProc
                offsets(i) = offsets(i - 1) + offsets(i)
            end do

            ! Get my starting index.
            procStarts(lvl) = offsets(myid)

        end do

        ! Now set up the interp and the coarse-level indices. Note that this
        ! loop is the number of levels MINUS 1 since we are generating the
        ! interpolations between levels and we already have the 1st level
        ! of the  coarseIndices

        call VecCreateMPI(adflow_comm_world, nCellsLocal(1_intType), &
                          PETSC_DETERMINE, indexVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        do lvl = 1, amgLevels - 1

            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                allocate (coarseIndices(nn, lvl)%arr(0:ib, 0:jb, 0:kb))
                allocate (coarseOversetIndices(nn, lvl)%arr(8, 0:ib, 0:jb, 0:kb))
                coarseIndices(nn, lvl)%arr = -1
                coarseOversetIndices(nn, lvl)%arr = -1
            end do

            ! Allocate the linear algebra interpolation array for this level
            ! (first count the number of nodes to be restricted on level lvl)
            n = 0
            do nn = 1, nDom
                n = n + sizes(1, nn, lvl) * sizes(2, nn, lvl) * sizes(3, nn, lvl)
            end do
            allocate (interps(lvl)%arr(n))

            ! Interps uses the local ordering so we always start at 0.
            n = 0
            count = 0
            ! Loop over the blocks
            do nn = 1, nDom

                ! Sizes for next level
                nnx = sizes(1, nn, lvl + 1)
                nny = sizes(2, nn, lvl + 1)
                nnz = sizes(3, nn, lvl + 1)

                ! Loop over the sizes of this level
                do k = 1, sizes(3, nn, lvl)
                    do j = 1, sizes(2, nn, lvl)
                        do i = 1, sizes(1, nn, lvl)

                            ! These are the indices on the next level
                            ii = (i - 1) / 2 + 1
                            jj = (j - 1) / 2 + 1
                            kk = (k - 1) / 2 + 1

                            coarseIndex = (kk - 1) * nnx * nny + (jj - 1) * nnx + ii + count

                            ! Linear algebra info
                            n = n + 1
                            interps(lvl)%arr(n) = coarseIndex

                        end do
                    end do
                end do
                count = count + nnx * nny * nnz
            end do

            ! We are not done yet; We need to fill in the block-based
            ! coarse indices and then do a halo-exchange on it so procs
            ! know where to put their off-proc entries on the coarser levels.

            ! Loop over the blocks.
            n = 0
            ii = 0

            call VecGetArrayF90(indexVec, indPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                ii = ii + (kb + 1) * (jb + 1) * (ib + 1) ! Count maximum double halos
                ! Loop over the sizes of this level
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            n = n + 1

                            ! Coarse Index: This is the first coarse index for
                            ! the current finest level element.
                            coarseIndex = interps(1)%arr(n)

                            ! For levels higher than 2, we need to trace
                            ! through the subsequent levels to find what the
                            ! coarse index is for level lvl.
                            do nextLvl = 2, lvl
                                coarseIndex = interps(nextLvl)%arr(coarseIndex)
                            end do

                            ! Now we set the coarse level index into the
                            ! coarseIndices array. Note that we make the index
                            ! global here by adding procStarts. We also
                            ! subtract 1 to make it zero-based for petsc.
                            coarseIndex = coarseIndex + procStarts(lvl + 1) - 1
                            coarseIndices(nn, lvl)%arr(i, j, k) = coarseIndex

                            ! Now put that index into the array
                            indPtr(n) = transfer(int(coarseIndex, 8), one)

                        end do
                    end do
                end do
            end do
            call VecRestoreArrayF90(indexVec, indPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            allocate (indicesToGet(8 * ii))
            ii = 0
            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                do k = 0, kb
                    do j = 0, jb
                        do i = 0, ib
                            ! If this cell is is an interpolated cell, record the indices we need to get
                            if (iblank(i, j, k) == -1) then
                                do m = 1, 8
                                    if (flowDoms(nn, 1, 1)%gInd(m, i, j, k) >= 0) then
                                        ii = ii + 1
                                        indicesToGet(ii) = flowDoms(nn, 1, 1)%gInd(m, i, j, k)
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do

            ! Now create the scatter to retrieve the "indicesToGet" from
            ! the indexVec. Petsc is always annoying for this.

            call ISCreateGeneral(adflow_comm_world, ii, indicesToGet(1:ii), PETSC_COPY_VALUES, &
                                 IS1, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            deallocate (indicesToGet)

            ! Create array to dump the result
            call VecCreateMPI(adflow_comm_world, ii, PETSC_DETERMINE, recvVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Create the scatter
            call VecScatterCreate(indexVec, IS1, recvVec, PETSC_NULL_IS, scat, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Do the actual scatter
            call VecScatterBegin(scat, indexVec, recvVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call VecScatterEnd(scat, indexVec, recvVec, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecGetArrayF90(recvVec, indPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Loop back over at set the coarse indices
            ii = 0
            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                ! Loop over the sizes of this level
                do k = 0, kb
                    do j = 0, jb
                        do i = 0, ib
                            if (iblank(i, j, k) == -1) then
                                do m = 1, 8
                                    if (flowDoms(nn, 1, 1)%gInd(m, i, j, k) >= 0) then
                                        ii = ii + 1
                                        coarseOversetIndices(nn, lvl)%arr(m, i, j, k) = &
                                            int(transfer(indPtr(ii), 1_8), intType)
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
            end do

            call VecRestoreArrayF90(recvVec, indPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecScatterDestroy(scat, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(recvVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call ISDestroy(IS1, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Now we need to exchange the coarse level indices. Set
            ! pointers to the coarseIndices arrays and then call the
            ! generic integer halo exchange.
            level = 1
            sps = 1
            nVar = 1

            do nn = 1, nDom
                flowDoms(nn, level, sps)%intCommVars(1)%var => coarseIndices(nn, lvl)%arr(:, :, :)
            end do

            call whalo1to1IntGeneric(nVar, level, sps, commPatternCell_2nd, internalCell_2nd)

        end do

        call VecDestroy(indexVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Next we need to set up the matrices and vectors

        do lvl = 1, amgLevels

            call KSPCreate(adflow_comm_world, kspLevels(lvl), ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Create the coarse levels
            if (lvl >= 2) then
                ncoarse = maxval(interps(lvl - 1)%arr)

                allocate (nnzOn(1:ncoarse), nnzOff(1:ncoarse))
                nnzOn = 14
                nnzOff = 7
                call MatCreateBAIJ(adflow_comm_world, bs, ncoarse * bs, ncoarse * bs, &
                                   PETSC_DETERMINE, PETSC_DETERMINE, 0, nnzOn, 0, nnzOff, &
                                   A(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)
                deallocate (nnzOn, nnzOff)

                call MatSetOption(A(lvl), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call MatSetOption(A(lvl), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call MatCreateVecs(A(lvl), res(lvl), sol(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

            else
                call vecCreateMPI(adflow_comm_world, nCell * bs, PETSC_DETERMINE, res(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDuplicate(res(lvl), sol(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

            end if
            call VecDuplicate(res(lvl), rhs(lvl), ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(res(lvl), sol2(lvl), ierr)
            call EChk(ierr, __FILE__, __LINE__)

        end do
        amgSetup = .True.

    end subroutine setupAMG

    subroutine destroyAMG

        integer(kind=intType) :: lvl, ierr, i, j
        if (amgSetup) then
            do lvl = 1, amgLevels
                ! Destroy all of our vectors/matrices
                if (lvl > 1) then
                    call MatDestroy(A(lvl), ierr)
                    call EChk(ierr, __FILE__, __LINE__)
                end if

                call VecDestroy(res(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(sol(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(sol2(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(rhs(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call KSPDestroy(kspLevels(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)
            end do

            deallocate (A, res, rhs, sol, sol2, kspLevels)
            deallocate (coarseIndices, coarseOversetIndices, interps)
            amgSetup = .False.

        end if
    end subroutine destroyAMG

    subroutine applyShellPC(pc, x, y, ierr)
        use communication
        ! Input/Output
        PC pc
        Vec x, y
        integer(kind=intType) ierr

        ! Working
        integer(kind=intType) :: i

        if (amgLevels > 1) then

            ! Copy the input vector x into rhs because x is read-only
            call VecCopy(x, rhs(1), ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Zero the output vector y
            call VecSet(y, zero, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (amgOuterIts == 1) then
                call MGPreCon(rhs(1), y, 1) ! y is the new approximate solution
            else

                do i = 1, amgOuterIts

                    call MGPreCon(rhs(1), sol(1), 1) ! The solution update is stored in sol(1)

                    ! Update the solution
                    call VecAYPX(y, one, sol(1), ierr) ! y = y + sol(1)
                    call EChk(ierr, __FILE__, __LINE__)

                    if (i < amgOuterIts) then
                        ! Compute new residual
                        call matMult(fineMat, y, rhs(1), ierr)
                        call EChk(ierr, __FILE__, __LINE__)

                        call VecAYPX(rhs(1), -one, x, ierr)
                        call EChk(ierr, __FILE__, __LINE__)
                    end if

                end do
            end if
        else
            ! Solve the fine level
            ! This is equivalent to not using multigrid
            call KSPSolve(kspLevels(1), x, y, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

    end subroutine ApplyShellPC

    subroutine setupShellPC(pc, ierr)

        use communication, only: adflow_comm_world
        use inputAdjoint

        ! Input/Output
        PC pc
        integer(kind=intTYpe) :: ierr

        ! Working
        PC globalPC, subpc
        KSP subksp
        integer(kind=intType) :: lvl
        integer(kind=intType) :: nlocal, first

        do lvl = 1, amgLevels
            if (lvl == 1) then
                call KSPSetOperators(kspLevels(lvl), fineMat, fineMat, ierr)
                call EChk(ierr, __FILE__, __LINE__)
            else
                call KSPSetOperators(kspLevels(lvl), A(lvl), A(lvl), ierr)
                call EChk(ierr, __FILE__, __LINE__)
            end if

            call kspsetnormtype(kspLevels(lvl), KSP_NORM_NONE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPSetType(kspLevels(lvl), 'richardson', ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPSetTolerances(kspLevels(lvl), PETSC_DEFAULT_REAL, &
                                  PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
                                  amgNSmooth, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPgetPC(kspLevels(lvl), globalPC, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call PCSetType(globalPC, 'asm', ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call PCASMSetOverlap(globalPC, amgASMOverlap, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set up the main ksp context before extracting the subdomains
            call KSPSetUp(kspLevels(lvl), ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Extract the ksp objects for each subdomain
            call PCASMGetSubKSP(globalPC, nlocal, first, subksp, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPSetType(subksp, 'preonly', ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Extract the preconditioner for subksp object
            call KSPGetPC(subksp, subpc, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! The subpc type is always ILU
            call PCSetType(subpc, 'ilu', ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set up the matrix ordering for the subpc object
            call PCFactorSetMatOrderingtype(subpc, amgMatrixOrdering, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the ILU parameters
            call PCFactorSetLevels(subpc, amgFillLevel, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end do

    end subroutine setupShellPC

    subroutine destroyShellPC(pc, ierr)

        ! Input/Ouput
        PC pc
        integer(kind=intType) :: ierr

        ! Working
        integer(kind=intType) :: lvl

    end subroutine destroyShellPC

    subroutine restrictVec(x, y, interp)

        ! Input/Output
        Vec x, y
        integer(kind=intType), dimension(:), intent(in) :: interp

        ! Working
        real(kind=realType), pointer :: yPtr(:), xPtr(:)
        real(kind=realType), pointer :: xPtrBlk(:, :), yPtrBlk(:, :)
        integer(kind=intType) :: ierr, n, i, j

        ! Restrict x -> y
        call VecGetArrayF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayF90(y, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Number of block nodes
        n = size(interp)

        ! Convenience block-based pointers
        xPtrBlk(1:bs, 1:size(xPtr) / bs) => xPtr
        yPtrBlk(1:bs, 1:size(yPtr) / bs) => yPtr

        ! Zero the output array
        yPtr = zero

        ! Loop over the interpolation array, summing into the coarse array
        do i = 1, n
            j = interp(i)
            yPtrBlk(:, j) = yPtrBlk(:, j) + xPtrBlk(:, i)
        end do

        call VecRestoreArrayF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(y, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine restrictVec

    subroutine prolongVec(x, y, interp)

        ! Input/Output
        Vec x, y
        integer(kind=intType), dimension(:), intent(in) :: interp

        ! Working
        real(kind=realType), pointer :: yPtr(:), xPtr(:)
        real(kind=realType), pointer :: xPtrBlk(:, :), yPtrBlk(:, :)
        integer(kind=intType) :: ierr, n, i, j

        ! Prolong vector x -> y

        call VecGetArrayF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayF90(y, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Number of block nodes
        n = size(interp)

        ! Convenience pointers
        yPtrBlk(1:bs, 1:size(yPtr) / bs) => yPtr
        xPtrBlk(1:bs, 1:size(xPtr) / bs) => xPtr

        ! Loop over the interpolation array, injecting into the fine array
        do i = 1, n
            j = interp(i)
            yPtrBlk(:, i) = xPtrBlk(:, j)
        end do

        call VecRestoreArrayF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(y, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine prolongVec

    recursive subroutine MGPreCon(r, y, k)

        ! r is the residual and y is the approximate solution
        ! k is the level

        ! Input/Output
        Vec y, r
        integer(kind=intType), intent(in) :: k

        ! Working
        integer(kind=intType) :: ierr, i

        ! Step 1: Restrict the residual
        call restrictVec(r, rhs(k + 1), interps(k)%arr)

        ! Step 2: Solve the coarse level directly or recursively
        if (k == amglevels - 1) then
            ! The next level down is the bottom
            ! Break the recursion by solving
            call kspSolve(KSPLevels(k + 1), rhs(k + 1), sol(k + 1), ierr)
            call EChk(ierr, __FILE__, __LINE__)
        else
            ! Call the next level down recursively
            call MGPreCon(rhs(k + 1), sol(k + 1), k + 1)
        end if

        ! Step 3: Prolongate the solution
        call prolongVec(sol(k + 1), y, interps(k)%arr)

        ! Step 4: Compute the new residual
        if (k == 1) then
            call matMult(fineMat, y, res(k), ierr)  ! res = A(i) * z1(k)
        else
            call matMult(A(k), y, res(k), ierr)  ! res = A(i) * z1(k)
        end if
        call EChk(ierr, __FILE__, __LINE__)

        call VecAYPX(res(k), -one, r, ierr)    ! res = -res + r
        call EChk(ierr, __FILE__, __LINE__)

        ! Step 5: Relax using the smoother
        call kspSolve(KSPLevels(k), res(k), sol2(k), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call vecAXPY(y, one, sol2(k), ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine MGPreCon

end module amg

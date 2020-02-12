module cartMesh

  use oversetData
  use communication
  use utils
  use haloExchange
  use oversetPackingRoutines
  use su_cgns
  implicit none

 contains

  subroutine createCartMesh(level, sps)

    use constants
    use blockPointers
    use surfaceFamilies, only : BCFamGroups
    use su_cgns
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Params
    integer(kind=intType), intent(in) :: level, sps

    ! Working params
    integer(kind=intType) :: i, j, k, l, ii, jj, kk, iBeg, iEnd, jBeg
    integer(kind=intType) :: jEnd, kBeg, kEnd, nn, mm, iDim, symOnFace(6)
    integer(kind=intType) :: count, countLocal , symOnFaceLocal(6),procDims(3)
    integer(kind=intType) :: cellDims(3), nNodes, nCells, ierr, globalInd
    integer(kind=intType) :: nSubI, nSubJ, nSeed, iSeed, nChanged, loopIter
    integer(kind=intType) :: iSize, jSize, kSize, nChangedLocal, stackPointer
    integer(kind=intType), dimension(:), allocatable :: indices
    integer(kind=intType), dimension(:,:), allocatable :: lSizes
    integer(kind=intType), dimension(:, :), allocatable :: stack, floodSeeds

    real(kind=realType), dimension(:), pointer :: cartPointer
    real(kind=realType), dimension(:, :, :), pointer :: xx
    real(kind=realType), dimension(:, :, :), pointer :: arrVals, changed
    real(kind=realType), dimension(:), allocatable :: values
    real(kind=realType), dimension(3) :: xMinLocal, xMaxLocal, xMin, xMax, sss
    real(kind=realType), dimension(3) ::  pt1, pt2, pt3, pt4, newPt, v1, v2
    real(kind=realType) :: areaLocal, area, areaAvg, err1, err2, h
    real(kind=realType) :: coorAvg, scaleSize, length, u, v

    DM cartArray
    AO cartOrdering
    Vec cartVecGlobal, cartVecLocal, blankVec, blankVecLocal, changedVecLocal, changedVecGlobal
    IS IS1, IS2
    VecScatter blankScatter, blankScatterLocal


    xMinLocal = large
    xMaxLocal = -large
    areaLocal = zero
    countLocal = 0
    ! First we have to determine the bounding box of our surfaces.

    do nn=1,nDom
       call setPointers(nn, level, sps)
       do mm=1,nBocos
          if(  BCType(mm) == NSWallAdiabatic .or. &
               BCType(mm) == NSWallIsothermal .or. &
               BCType(mm) == EulerWall) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             select case (BCFaceID(mm))
             case (iMin)
                xx   => x(1,:,:,:)
             case (iMax)
                xx   => x(il,:,:,:)
             case (jMin)
                xx   => x(:,1,:,:)
             case (jMax)
                xx   => x(:,jl,:,:)
             case (kMin)
                xx   => x(:,:,1,:)
             case (kMax)
                xx   => x(:,:,kl,:)
             end select

             do j=jBeg+1, jEnd+1
                do i=iBeg+1,iEnd+1
                   do iDim=1,3
                      xMinLocal(iDim) = min(xMinLocal(iDim), xx(i, j, iDim))
                      xMaxLocal(iDim) = max(xMaxLocal(iDim), xx(i, j, iDim))
                   end do
                end do
             end do

             ! Determine the total area of the patch:
             do j=jBeg+1, jEnd
                do i=iBeg+1, iEnd
                   v1 = xx(i, j, :) - xx(i-1, j-1, :)
                   v2 = xx(i-1, j, :) - xx(i, j-1, :)
                   ! Cross Product
                   sss(1) = (v1(2)*v2(3) - v1(3)*v2(2))
                   sss(2) = (v1(3)*v2(1) - v1(1)*v2(3))
                   sss(3) = (v1(1)*v2(2) - v1(2)*v2(1))

                   areaLocal = areaLocal + half*(sss(1)**2 + sss(2)**2 + sss(3)**2)
                   countLocal = countLocal + 1
                end do
             end do
          end if
       end do
    end do

    call MPI_Allreduce(xMinLocal, xMin, 3, adflow_real, MPI_MIN, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MPI_Allreduce(xMaxLocal, xMax, 3, adflow_real, MPI_MAX, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MPI_Allreduce(areaLocal, area, 1, adflow_real, MPI_SUM, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MPI_Allreduce(countLocal, count, 1, adflow_integer, MPI_SUM, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Compute the average area:
    areaAvg = area/count

    ! Compute the master size 'h':
    h = sqrt(areaAvg)
    h = h * 6
    ! Before we setup the cartesian mesh, we will expand the bounds
    ! slightly to ensure that there isn't a wall cell right at the
    ! boundary. This is actually vastly tricker than it sounds becuase
    ! we have to correctly account for symmetry planes. That is we
    ! CANNOT
    ! outside-in flood would just "go around the back" of the symmetry
    ! plane and flood out the inside. Not good. So what we do is now
    ! that we know the exact bounds of our surface, we loop over the
    ! symmetry planes and flag values in an array of lenght 6
    ! corresponding to (iLow, iHigh, jLow, jHigh, kLow, kHigh). We set
    ! the value to 1 if the symmetry plane I'm looking at corresponds
    ! to that value. Then we all reduce so that everyone knows which
    ! of the 6 faces cannot be extended.

    scaleSize = mynorm2(xMax-xMin)

    symOnFaceLocal = 0
    do nn=1,nDom
       call setPointers(nn, level, sps)
       do mm=1,nBocos
          if(  BCType(mm) ==symm) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             select case (BCFaceID(mm))
             case (iMin)
                xx   => x(1,:,:,:)
             case (iMax)
                xx   => x(il,:,:,:)
             case (jMin)
                xx   => x(:,1,:,:)
             case (jMax)
                xx   => x(:,jl,:,:)
             case (kMin)
                xx   => x(:,:,1,:)
             case (kMax)
                xx   => x(:,:,kl,:)
             end select

             ! First we have to determine the principle coordinate
             ! direction of the face. Do this with symNorm which is
             ! already computed.

             ! Location, ie coordiante direction of dominate direction
             iDim = maxloc(abs(bcData(mm)%symNorm), 1)

             ! Now determine the average value in "iDim" dimension
             coorAvg = sum(xx(iBeg+1:iEnd+1, jBeg+1:jEnd+1, iDim))/((iEnd-iBeg+1)*(jEnd-jBeg+1))

             ! Check if it is sufficently close to the bounding box:
             err1 = abs(coorAvg - xMin(iDim))/scaleSize
             err2 = abs(coorAvg - xMax(iDim))/scaleSize

             if (err1 < 1e-8) then
                symOnFaceLocal(2*iDim-1) = 1
             end if

             if (err2< 1e-8) then
                symOnFaceLocal(2*iDim  ) = 1
             end if
          end if
       end do
    end do

    ! Now we all reduce with max. Values that are not zero have
    ! symmetry planes on them.

    call MPI_Allreduce(symOnFaceLocal, symOnFace, 6, adflow_integer, MPI_MAX, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Now we will expand xmin/xmax by 2h (just to be sure) but only
    ! on the planes without symmetry conditions

    do iDim=1,3
       if (symOnFace(2*iDim-1) == 0) then
          xMin(iDim) = xMin(iDim) - 2*h
       end if

       if (symOnFace(2*iDim  ) == 0) then
          xMax(iDim) = xMax(iDim) + 2*h
       end if
    end do

    ! Ok. Now we have the final size of our cartesian block. The
    ! next step is to determine how it is to be partitioned.
    call tripleFactor(nProc, procDims)

    ! And the dimensions in each direction. Remember the
    ! distributed array is cell-centered based.
    cellDims = floor((xMax-xMin)/h) + 1

    ! Now reorder the procDims so that the largest dimension
    ! corresponds to largest number of procs. Since we know procDims
    ! are sorted, we can do the 6 cases:

    if (cellDims(1) >= cellDims(2) .and. cellDims(2) >= cellDims(3)) then
       procDims = (/procDims(3), procDims(2), procDims(1)/)
    else if (cellDims(1) >= cellDims(3) .and. cellDims(3) >= cellDims(2)) then
       procDims = (/procDims(3), procDims(1), procDims(2)/)

    else if (cellDims(2) >= cellDims(1) .and. cellDims(1) >= cellDims(3)) then
       procDims = (/procDims(2), procDims(3), procDims(1)/)
    else if (cellDims(2) >= cellDims(3) .and. cellDims(3) >= cellDims(1)) then
       procDims = (/procDims(1), procDims(3), procDims(2)/)

    else if (cellDims(3) >= cellDims(1) .and. cellDims(1) >= cellDims(2)) then
       procDims = (/procDims(2), procDims(1), procDims(3)/)
    else if (cellDims(3) >= cellDims(2) .and. cellDims(2) >= cellDims(1)) then
       procDims = (/procDims(1), procDims(2), procDims(3)/)
    end if

    ! Since we can't pass in null array objects to petsc (this is
    ! still broken) we have to do up all the arguments ourself. In
    ! particular, the lx, ly and lz one which are quite annonying.

    allocate(lSizes(maxval(procDims), 3))
    do iDim=1, 3
       do j=1, procDims(iDim)

          ii = cellDims(iDim)/procDims(iDim)
          if (mod(cellDims(iDim), procDims(iDim)) > j-1) then
             ii = ii + 1
          end if
          lSizes(j, iDim) = ii
       end do
    end do
    if (myid == 0) then
       print *,'xmin:', xmin
       print *,'xmax:', xmax
       print *,'dims:', cellDims
       print *, 'H:', h
       print *, 'procDims:', procDims
       print *,'I sizes:', lSizes(1:procDims(1), 1)
       print *,'J sizes:', lSizes(1:procDims(2), 2)
       print *,'K sizes:', lSizes(1:procDims(3), 3)
    end if

    call DMDAcreate3d(adflow_comm_world, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, &
         DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, cellDims(1), cellDims(2), &
         cellDims(3),  procDims(1), procDims(2), procDims(3), 1, 1, &
         lSizes(1:procDims(1), 1), lSizes(1:procDims(2), 2), lSizes(1:procDims(3), 3), &
         cartArray, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    deallocate(lSizes)

    ! Now loop back over the surfaces. For each surface, determine
    ! it's global index of the point that will be flagged as a cut
    ! cell. Due to the stupid PETSC ordering crap, we have to store
    ! all the indices as we go.

    i = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)
       call getWallSize(BCFamGroups(iBCGroupWalls)%famList, nNodes, nCells, .False.)
       i = i + nNodes
    end do
    i = i * 10
    allocate(indices(i))
    count = 0
    do nn=1,nDom
       call setPointers(nn, level, sps)
       do mm=1,nBocos
          if(  BCType(mm) == NSWallAdiabatic .or. &
               BCType(mm) == NSWallIsothermal .or. &
               BCType(mm) == EulerWall) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             select case (BCFaceID(mm))
             case (iMin)
                xx   => x(1,:,:,:)
             case (iMax)
                xx   => x(il,:,:,:)
             case (jMin)
                xx   => x(:,1,:,:)
             case (jMax)
                xx   => x(:,jl,:,:)
             case (kMin)
                xx   => x(:,:,1,:)
             case (kMax)
                xx   => x(:,:,kl,:)
             end select


             ! Loop over Raw Nodes
             do j=jBeg, jEnd
                do i=iBeg, iEnd

                   ! Note that ii, jj, kk are zero based since we are
                   ! working with PETSc indices.
                   ii = int((xx(i+1, j+1, 1) - xMin(1))/h)
                   jj = int((xx(i+1, j+1, 2) - xMin(2))/h)
                   kk = int((xx(i+1, j+1, 3) - xMin(3))/h)

                   globalInd = cellDims(1)*cellDims(2)*kk + cellDims(1)*jj + ii
                   count = count + 1
                   indices(count) = globalInd
                end do
             end do

             ! Loop over iEdges
             do j=jBeg, jEnd
                do i=iBeg+1, iEnd

                   ! Check the Length of this edge
                   pt1 = xx(i  , j+1, :)
                   pt2 = xx(i+1, j+1, :)

                   length = mynorm2(pt1-pt2)
                   nSubI = int(length/h)

                   do k=1, nSubI

                      newPt = pt1 + dble(k)/(nSubI+1)*(pt2-pt1)
                      ii = int((newPt(1) - xMin(1))/h)
                      jj = int((newPt(2) - xMin(2))/h)
                      kk = int((newPt(3) - xMin(3))/h)

                      globalInd = cellDims(1)*cellDims(2)*kk + cellDims(1)*jj + ii
                      count = count + 1
                      indices(count) = globalInd
                   end do
                end do
             end do

             ! Loop over jEdges
             do j=jBeg+1, jEnd
                do i=iBeg, iEnd

                   ! Check the Length of this edge
                   pt1 = xx(i+1, j  , :)
                   pt2 = xx(i+1, j+1, :)

                   length = mynorm2(pt1-pt2)
                   nSubJ = int(length/h)
                   do k=1, nSubJ

                      newPt = pt1 + dble(k)/(nSubJ+1)*(pt2-pt1)
                      ii = int((newPt(1) - xMin(1))/h)
                      jj = int((newPt(2) - xMin(2))/h)
                      kk = int((newPt(3) - xMin(3))/h)

                      globalInd = cellDims(1)*cellDims(2)*kk + cellDims(1)*jj + ii
                      count = count + 1
                      indices(count) = globalInd
                   end do
                end do
             end do

             ! Loop over faces
             do j=jBeg+1, jEnd
                do i=iBeg+1, iEnd

                   ! Extract the 4 pts. CCW ordering
                   pt1 = xx(i  , j  , :)
                   pt2 = xx(i+1, j  , :)

                   pt3 = xx(i+1, j+1, :)
                   pt4 = xx(i  , j+1, :)

                   ! Sub pts in I
                   length = mynorm2(pt2-pt1)
                   nSubI = int(length/h)

                   length = mynorm2(pt3-pt4)
                   nSubI = max(nSubI, int(length/h))

                   ! Sub Pts in J
                   length = mynorm2(pt4-pt1)
                   nSubJ = int(length/h)

                   length = mynorm2(pt3-pt2)
                   nSubJ = max(nSubJ, int(length/h))

                   do l=1, nSubJ
                      do k=1, nSubI
                         u = dble(k)/(nSubI+1)
                         v = dble(l)/(nSubJ+1)

                         newPt = (one-u)*(one-v)*pt1 + u*(one-v)*pt2 + &
                              u*v*pt3 + (one-u)*v*pt4

                         ii = int((newPt(1) - xMin(1))/h)
                         jj = int((newPt(2) - xMin(2))/h)
                         kk = int((newPt(3) - xMin(3))/h)

                         globalInd = cellDims(1)*cellDims(2)*kk + cellDims(1)*jj + ii
                         count = count + 1
                         indices(count) = globalInd
                      end do
                   end do
                end do
             end do
          end if
       end do
    end do

    call DMCreateGlobalVector(cartArray, cartVecGlobal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Initialize all values to 1 ("Compute")
    call vecSet(cartVecGlobal, one, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call DMDAGetAO(cartArray, cartOrdering, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Convert the indices from application ordering to petsc
    call AOApplicationToPetsc(cartOrdering, size(indices), indices, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Now we set all the values, with a simple single vecSet call
    allocate(values(count))
    values = -three ! -3 is flood seed. Use the same notation here.

    call vecSetValues(cartVecGlobal, count, indices, values, INSERT_VALUES, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    deallocate(values, indices)

    ! Don't forget to assemble
    call vecAssemblyBegin(cartVecGlobal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call vecAssemblyEnd(cartVecGlobal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! These are "get" vecs. We have to restore them later.
    call DMGetLocalVector(cartArray, cartVecLocal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call DMGetLocalVector(cartArray, changedVecLocal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call VecSet(changedVecLocal, zero, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call DMGetGlobalVector(cartArray, changedVecGlobal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Now the next step is to perform the flooding from the outside in.
    call DMGlobalToLocalBegin(cartArray, cartVecGlobal, INSERT_VALUES, cartVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call DMGlobalToLocalEnd(cartArray, cartVecGlobal, INSERT_VALUES, cartVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Determine the bounds of the arrays we will be getting back.
    call DMDAVecGetArrayF90(cartArray, cartVecLocal, arrVals, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    iBeg = lbound(arrVals, 1)
    jBeg = lbound(arrVals, 2)
    kBeg = lbound(arrVals, 3)

    iEnd = ubound(arrVals, 1)
    jEnd = ubound(arrVals, 2)
    kEnd = ubound(arrVals, 3)

    iSize = iEnd - iBeg + 1
    jSize = jEnd - jBeg + 1
    kSize = kEnd - kBeg + 1

    call DMDAVecRestoreArrayF90(cartArray, cartVecLocal, arrVals, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    loopIter = 1
    allocate(stack(3, 6*iSize*jSize*kSize + 1))
    allocate(floodSeeds(3, 2*iSize*jSize + 2*iSize*kSize + 2*jSize*kSize))

    parallelSyncLoop: do

       call DMDAVecGetArrayF90(cartArray, cartVecLocal, arrVals, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call DMDAVecGetArrayF90(cartArray, changedVecLocal, changed, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       ! Keep track of the total number of fringes we've modified
       nChangedLocal = 0

       ! Allocate space for our queue (stack). It needs to be 6*nx*ny*nz + 1:
       ! 6 for each of the 6 coordinate directions plus our extra
       ! seed. It should never come close to this unless the entire
       ! block will be blanked.

       nSeed = 0
       if (loopIter == 1) then

          if (myid == 0) then

             ! Set the single seed on the bottom corner of the root proc

             call addSeed(0,0,0)

          end if

       else

          ! On the second and subsequent passes, check each 1st
          ! non-corner halos in the 6 faces to see if we received
          ! "changed" info from neighbour proc. This will allow us to
          ! continue the flooding on this processor/block. Note that
          ! even in a single processor case, the halo exchange in
          ! necessary to communicate between two local blocks

          ! iMin/iMax
          do k=kBeg+1, kEnd-1
             do j=jBeg+1, jEnd-1
                if (int(changed(iBeg  , j, k)) == 1) then
                   call addSeed(iBeg+1, j, k)
                end if
                if (int(changed(iEnd  , j, k)) == 1) then
                   call addSeed(iEnd-1, j, k)
                end if
             end do
          end do

          ! jMin/jMax
          do k=kBeg+1, kEnd-1
             do i=iBeg+1, iEnd-1
                if (int(changed(i, jBeg,   k)) == 1) then
                   call addSeed(i, jBeg+1, k)
                end if
                if (int(changed(i, jEnd,   k)) == 1) then
                   call addSeed(i, jEnd-1, k)
                end if
             end do
          end do

          ! kMin:
          do j=jBeg+1, jEnd-1
             do i=iBeg+1, iEnd-1
                if (int(changed(i, j, kBeg )) == 1) then
                   call addSeed(i, j, kBeg+1)
                end if
                if (int(changed(i, j, kEnd )) == 1) then
                   call addSeed(i, j, kEnd-1)
                end if
             end do
          end do
       end if

       ! Loop over our seeds we currently have
       do iSeed = 1, nSeed

          ! Put the particular seed in the first slot of the stack
          stack(:, 1) = floodSeeds(:, iSeed)

          ! Reset the stack pointer length back to just the one seed
          ! we have
          stackPointer = 1

          ! Start the flooding (stacked based, not recursive)
          do while (stackPointer > 0 )

             ! 'Pop' the current point off the stack
             i = stack(1, stackPointer)
             j = stack(2, stackPointer)
             k = stack(3, stackPointer)
             stackPointer = stackPointer - 1

             if (int(arrVals(i, j, k)) == 1) then

                ! Flag the cell (using changed) as being changed
                changed(i, j, k) = one

                ! Keep track of the total number we've changed.  For
                ! reporting purposes...only count the ones that are
                ! on actual compute cells:
                if (onBlock(i, j, k)) then
                   nChangedLocal = nChangedLocal + 1
                end if

                ! Set the value as flooded. Use 0 as per usual.
                arrVals(i, j, k) = zero

                ! Now add the six nearest neighbours to the stack
                ! provided they are in the owned cell range:

                if (i-1 >= iBeg) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i-1, j  , k  /)
                end if

                if (i+1 <= iEnd) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i+1, j  , k  /)
                end if

                if (j-1 >= jBeg) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i  , j-1, k  /)
                end if

                if (j+1 <= jEnd) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i  , j+1, k  /)
                end if

                if (k-1 >= kBeg) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i  , j  , k-1/)
                end if

                if (k+1 <= kEnd) then
                   stackPointer = stackPointer + 1
                   stack(:, stackPointer) = (/i  , j , k+1 /)
                end if
             end if
          end do
       end do

       call DMDAVecRestoreArrayF90(cartArray, cartVecLocal, arrVals, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call DMDAVecRestoreArrayF90(cartArray, changedVecLocal, changed, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       ! Exchange "changed"
       call DMLocalToGlobalBegin(cartArray, changedVecLocal, INSERT_VALUES, changedVecGlobal, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call DMLocalToGlobalEnd(cartArray, changedVecLocal, INSERT_VALUES, changedVecGlobal, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call DMGlobalToLocalBegin(cartArray, changedVecGlobal, INSERT_VALUES, changedVecLocal, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call DMGlobalToLocalEnd(cartArray, changedVecGlobal, INSERT_VALUES, changedVecLocal, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       ! Determine if cells got changd. If so do another loop.
       call mpi_allreduce(nChangedLocal, nChanged, 1, adflow_integer, MPI_SUM, &
            adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
       if (myid == 0) then
          print *, 'Cart Flood Iteration:', loopIter, 'Blanked ', nChanged, 'Interior Cells.'
       end if

       if (nChanged == 0) then
          exit parallelSyncLoop
       end if

       loopIter = loopIter + 1
    end do parallelSyncLoop

    deallocate(stack, floodSeeds)

    ! Now that we have flooded everything, any cells left over must
    ! be *inside. Do one last pass through and flip those.

    call DMDAVecGetArrayF90(cartArray, cartVecLocal, arrVals, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    do k=kBeg+1, kEnd-1
       do j=jBeg+1, jEnd-1
          do i=iBeg+1, iEnd-1
             if (int(arrVals(i, j, k)) == 1) then
                arrVals(i, j, k) = -three
             end if
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(cartArray, cartVecLocal, arrVals, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Restore the vectors obtained with "get"
    call DMRestoreLocalVector(cartArray, cartVecLocal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call DMRestoreLocalVector(cartArray, changedVecLocal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call DMRestoreGlobalVector(cartArray, changedVecGlobal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Now update the final global variables in cartArray.
    call DMLocalToGlobalBegin(cartArray, cartVecLocal, INSERT_VALUES, cartVecGlobal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call DMLocalToGlobalEnd(cartArray, cartVecLocal, INSERT_VALUES, cartVecGlobal, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Now that we are done with the flooding and any local
    ! communication, we can create a global vector that is in the
    ! ordering that we actually want.

    i = cellDims(1)*cellDims(2)*cellDims(3)
    call VecCreateMPI(adflow_comm_world, PETSC_DECIDE, i, blankVec, ierr)
    call EChk(ierr, __FILE__, __LINE__)


    ! Create the index for the real global Vector
    call VecGetOwnershipRange(blankVec, i, j, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call ISCreateStride(adflow_comm_world, j-i, i, 1, IS2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call ISDuplicate(IS2, IS1, ierr)
    call AOApplicationToPetscIS(cartOrdering, IS1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call vecScatterCreate(cartVecGlobal, IS1, blankVec, IS2, blankScatter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterBegin(blankScatter, cartVecGlobal, blankVec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(blankScatter, cartVecGlobal, blankVec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! We are now done with everything petsc related except for
    ! blankVec. That's all we  need now.
    call VecScatterDestroy(blankScatter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call ISDestroy(IS1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call ISDestroy(IS2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecDestroy(cartVecGlobal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call DMDestroy(cartArray, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Write our mesh
    call writeCartMesh(blankVec, cellDims, xMin, h)


    ! Now what we have to do is to distribute parts of the cart mesh
    ! to the processors that need it. For now, just do a create to all.

    call vecScatterCreateToAll(blankVec, blankScatterLocal, blankVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterBegin(blankScatterLocal, blankVec, blankVecLOCAL, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(blankScatterLocal, blankVec, blankVecLOCAL, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterDestroy(blankScatterLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(blanKVecLocal, cartPointer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    do nn=1,nDom
       call setPointers(nn, level, sps)
       do k=2, kl
          do j=2, jl
             do i=2, il
                ii = ii + 1
                do iDim=1, 3
                   pt1(iDim) = eighth*(&
                        x(i-1, j-1, k-1, iDim) + &
                        x(i  , j-1, k-1, iDim) + &
                        x(i-1, j  , k-1, iDim) + &
                        x(i  , j  , k-1, iDim) + &
                        x(i-1, j-1, k  , iDim) + &
                        x(i  , j-1, k  , iDim) + &
                        x(i-1, j  , k  , iDim) + &
                        x(i  , j  , k  , iDim))
                end do

                ! Now Simply check if the cell this point is in has a value of -3
                ii = int((pt1(1) - xMin(1))/h)
                jj = int((pt1(2) - xMin(2))/h)
                kk = int((pt1(3) - xMin(3))/h)

                ! Clip the bounds to the actual ranges:
                ii = min(max(0, ii), cellDims(1)-1)
                jj = min(max(0, jj), cellDims(2)-1)
                kk = min(max(0, kk), cellDims(3)-1)

                ! GlobalInd is 0 based
                globalInd = cellDims(1)*cellDims(2)*kk + cellDims(1)*jj + ii

                ! Only blank if we are more than sqrt(3)*h from our own wall:
                if (xSeed(i, j, k, 1) < large) then
                   ! We have a wall: See how far we are away form it
                   length = mynorm2(pt1 - xSeed(i, j, k, :))
                   if (length > 1.73205080*h) then
                      if (int(cartPointer(globalInd+1) ) == -3) then
                         iblank(i, j, k) = -3
                      end if
                   end if
                else
                   ! No wall...no wall check.
                   if (int(cartPointer(globalInd+1) ) == -3) then
                      iblank(i, j, k) = -3
                   end if
                end if
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(blankVecLocal, cartPointer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Clean up the remaining PETScmemory
    call VecDestroy(blankVec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecDestroy(blankVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Update the iblank info.
    domainLoop:do nn=1, nDom
       flowDoms(nn, level, sps)%intCommVars(1)%var => &
            flowDoms(nn, level, sps)%iblank(:, :, :)
    end do domainLoop

    ! Run the generic integer exchange
    call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)


  contains
    ! Simple routine to make code easier to read above
    subroutine addSeed(i, j, k)
      use constants
      implicit none
      integer(kind=intType), intent(in) :: i, j, k
      nSeed = nSeed + 1
      floodSeeds(:, nSeed) = (/i, j, k/)
    end subroutine addSeed

    function onBlock(i, j, k)

      use constants
      implicit none

      integer(kind=intType), intent(in) :: i, j, k
      logical :: onBlock

      if (i >= 0 .and. i <= iEnd-1 .and. j >= 0 .and. j<= jEnd-1 .and. k >= 0 .and. k <= kEnd-1) then
         onBlock = .True.
      else
         onBlock = .False.
      end if
    end function onBlock
  end subroutine createCartMesh

  subroutine writeCartMesh(blankVec, cellDims, xMin, h)

#include <petsc/finclude/petsc.h>
   use petsc
   implicit none

    ! Input
    integer(kind=intType), intent(in), dimension(3) :: cellDims
    real(kind=realType), intent(in), dimension(3) :: xMin
    real(kind=realType), intent(in) :: h
    Vec blankVec
    ! Working
    integer(kind=intType) :: ierr

    ! CGNS
    character*32 :: coorNames(3)
    integer(kind=intType) :: base, zoneID, coordID, cg, zone,  iField, iSol, i, j, k, iDim
    real(kind=realType), dimension(:, :, :, :), allocatable :: xTmp
    real(kind=realType), dimension(:), pointer :: cartPointer

    Vec blankVecLocal
    IS IS1, IS2
    VecScatter blankScatterLocal

    coorNames(1) = "CoordinateX"
    coorNames(2) = "CoordinateY"
    coorNames(3) = "CoordinateZ"

    call VecScatterCreateToZero(blankVec, blankScatterLocal, blankVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterBegin(blankScatterLocal, blankVec, blankVecLocal, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterEnd(blankScatterLocal, blankVec, blankVecLocal, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (myid == 0) then
       ! Open the CGNS File
       call cg_open_f("cartblock.cgns", mode_write, cg, ierr)
       base = 1
       call cg_base_write_f(cg, "Base#1", 3, 3, base, ierr)

       call cg_zone_write_f(cg, base, "cartblock", int((/cellDims(1)+1, cellDims(2)+1, cellDims(3)+1, &
            cellDims(1), cellDims(2), cellDims(3), 0, 0, 0/), cgsize_t), Structured, zoneID, ierr)

       allocate(xtmp(cellDims(1)+1, cellDims(2)+1, cellDims(3)+1, 3))

       do k=1, cellDims(3)+1
          do j=1, cellDims(2)+1
             do i=1, cellDims(1)+1
                xTmp(i, j, k, 1) = xMin(1) + (i-1)*h
             end do
          end do
       end do

       do k=1, cellDims(3)+1
          do j=1, cellDims(2)+1
             do i=1, cellDims(1)+1
                xTmp(i, j, k, 2) = xMin(2) + (j-1)*h
             end do
          end do
       end do

       do k=1, cellDims(3)+1
          do j=1, cellDims(2)+1
             do i=1, cellDims(1)+1
                xTmp(i, j, k, 3) = xMin(3) + (k-1)*h
             end do
          end do
       end do

       do idim=1, 3
          call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
               xtmp(:, :, :, idim), coordID, ierr)
       end do
       deallocate(xtmp)

       call cg_sol_write_f(cg, base, zoneID, "flowSolution", CellCenter, iSol, ierr)
       call VecGetArrayF90(blanKVecLocal, cartPointer, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "iBlank", &
            cartPointer, iField, ierr)

       call VecRestoreArrayF90(blankVecLocal, cartPointer, ierr)
       call EChk(ierr,__FILE__,__LINE__)
       call cg_close_f(cg, ierr)
    end if

    call VecDestroy(blankVecLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecScatterDestroy(blankScatterLocal, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine writeCartMesh

  subroutine tripleFactor(N, s)
    use constants
    use sorting, only : qsortIntegers
    implicit none
    ! Input/Output
    integer(kind=intType), intent(in) :: N
    integer(kind=intType), intent(out) :: s(3)

    ! Working
    integer(kind=intType) :: a, b, a1, b1, a2, b2, s1(3), s2(3)

    ! Determine a set of triple factors for integer N

    call largeFactor(N, a, b)
    call largeFactor(b, a1, a2)
    call largeFactor(a, b1, b2)

    ! Our options are a, a1 and a2 OR b, b1, and b2
    s1 = (/a, a1, a2/)
    s2 = (/b, b1, b2/)

    ! Sort them
    call qsortIntegers(s1, 3)
    call qsortIntegers(s2, 3)

    ! And take the set that has the largest, smallest value.
    if (s1(1) > s2(1)) then
       s = s1
    else
       s = s2
    end if
  end subroutine tripleFactor

  subroutine largeFactor(N, f1, f2)

    implicit none
    integer(kind=intType) :: N, f1, f2, i, j, s
    ! Return the two factors that are closest to the sqrt(N)

    s = int(sqrt(dble(N)))
    do j=s, 1, -1
       if (mod(N, j) == 0) then
          f1 = N/j
          f2 = j
          exit
       end if
    end do
  end subroutine largeFactor


end module cartMesh

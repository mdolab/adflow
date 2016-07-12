subroutine writeIsoSurface(isoName , sps, nIsoSurfVar, isoSurfSolNames)

  ! Implements a marching cubes algrorithm which can be used to
  ! extract iso surfaces or slcies from a solution and store them in a
  ! CGNS surface file.

  use communication
  use blockPointers
  use flowVarRefState
  use inputPhysics
  use su_cgns
  use BCTypes
  use inputIO
  use outputMod
  use cgnsNames

  implicit none
  
  ! Input param
  character(len=*), intent(in)                   :: isoName
  integer(kind=intType) :: sps
  integer(kind=intType) :: nIsoSurfVar
  character(len=*), dimension(*), intent(in) :: isoSurfSolNames
  ! Working param
  integer(kind=intType) :: i, j, k, nn, kk, nMax, icon, iCoor, indexCube, num1, num2
  real(kind=realType), dimension(:, :, :), pointer :: fn
  real(kind=realType), dimension(:, :), pointer :: coords, uniqueCoords, weights
  real(kind=realType), dimension(:), allocatable :: buffer

  integer(kind=intType), dimension(:, :), pointer :: indices
  integer(kind=intType), dimension(:, :), allocatable :: connBuffer, allConn
  integer(kind=intType), dimension(:), allocatable :: link, nPtsProc, nConnProc

  integer(kind=intType) :: ccwOrdering(3, 8), n1(3), n2(3)
  integer(kind=intType) :: matCon1(256, 16), matCon2(12,2), nUnique, ivar
  integer(kind=intType) :: ierr, iProc, tag, cumNodes, cumConn, iCorner, curBlock, idim
  real(kind=realType) :: f(8)

  logical :: logic1
  integer, dimension(mpi_status_size) :: status
  integer(kind=intType) :: cgnsInd, cgnsBase, cgnsZOne, coordID, secID, solID, fieldID
  real(kind=realType) :: tol=1e-8 ! Node tol for isosurf pointReduce
  interface
     subroutine reallocateReal2(realArray,           &
          newSize1, newSize2, &
          oldSize1, oldSize2, &
          alwaysFreeMem)
       use precision
       implicit none

       real(kind=realType), dimension(:,:), pointer :: realArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
            oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
     end subroutine reallocateReal2

     subroutine reallocateInteger2(intArray, newSize1, newSize2, &
          oldSize1, oldSize2,           &
          alwaysFreeMem)
       use precision
       implicit none
       integer(kind=intType), dimension(:,:), pointer :: intArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                            oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
     end subroutine reallocateInteger2
  end interface

  ! Fill up the connecivity matrices
  call getMatCons(matCon1, matCon2, ccwOrdering)

  ! Generate a guess for the size of the iso surface -- sum of the
  ! number of nodes on the faces
  nMax = 0
  do nn=1, nDom
     call setPointers(nn, 1, sps)
     nMax = nMax + 2*il*jl + 2*il*kl + 2*jl*kl
  end do

  ! Allocate size nMax
  allocate(weights(2, nMax), indices(7, nMax))

  iCoor = 0
  do nn=1, nDom
     call setPointers(nn, 1, sps)
     fn => flowDoms(nn, 1, sps)%fn
     ! Now loop over the Cells:
     do k=2, kl
        do j=2, jl
           do i=2, il

              ! Extract value at the corners
              do iCorner=1,8
                 f(iCorner) = fn(i+ccwOrdering(1, iCorner), &
                      j + ccwOrdering(2, iCorner), &
                      k + ccwOrdering(3, iCorner))
              end do

              ! Based on the values at each corner, determine which
              ! type surface we have
              indexcube = 1
              if (f(1) .lt. zero) indexcube = indexcube + 1
              if (f(2) .lt. zero) indexcube = indexcube + 2
              if (f(3) .lt. zero) indexcube = indexcube + 4
              if (f(4) .lt. zero) indexcube = indexcube + 8
              if (f(5) .lt. zero) indexcube = indexcube + 16
              if (f(6) .lt. zero) indexcube = indexcube + 32
              if (f(7) .lt. zero) indexcube = indexcube + 64
              if (f(8) .lt. zero) indexcube = indexcube + 128

              logic1 = .true.

              kk = 1
              do while (logic1)
                 icon = matcon1(indexcube, kk)

                 if (icon == 0) then
                    logic1=.false.
                 else

                    iCoor = iCoor + 1
                    if (iCoor > nMax) then
                       ! Need to realloc the coord array. Make it double the size
                       call reallocateReal2(weights, 2, 2*nMax, 2, nMax, .true.)
                       call reallocateInteger2(indices, 7, 2*nMax, 7, nMax, .true.)
                       nMax = nMax * 2
                    end if

                    num1 = matcon2(icon,1)
                    num2 = matcon2(icon,2)

                    ! Weight factors
                    weights(2, iCoor) = (zero - f(num1))/(f(num2) - f(num1))
                    weights(1, iCoor) = one - weights(2, icoor)

                    ! Indices of nodes
                    n1 = (/i, j, k/) + ccwOrdering(: ,num1)
                    n2 = (/i, j, k/) + ccwOrdering(:, num2) 
                    indices(:, iCoor) = (/nn, n1(1), n1(2), n1(3), n2(1), n2(2), n2(3)/)              

                    kk = kk + 1
                 end if
              end do

           end do ! I loop
        end do ! J loop
     end do ! K loop
  end do ! Domain loop

  ! We have not actually stored the coordintes; only the positions and
  ! the weights. To compute the coordinates we pass back through and assemble
  allocate(Coords(3, iCoor))

  ! Set pointer to first block
  call setPointers(1, 1, sps)
  curBlock = 1
  do i=1,iCoor

     ! If we've switched blocks, reset points. This stil only calls
     ! setPointer nDom times since there are at most that many
     ! switches
     if (indices(1, i) /= curBlock) then
        call setPointers(indices(1, i), 1, sps)
        curBlock = indices(1, i)
     end if
     
     ! Computing coordinates is easy; we just juse the weights and the
     ! indices on x
     do idim=1,3
        coords(idim, i) = &
             weights(1, i) * &
             X(indices(2,i), indices(3,i), indices(4,i), idim) + &
             weights(2, i) * &
             X(indices(5,i), indices(6,i), indices(7,i), idim)
     end do
  end do

  ! Now we know the maximum number of coordinates so we can allocate
  ! the unique set and the link array
  allocate(uniqueCoords(3, icoor))
  allocate(link(icoor))

  ! Compute the reduced set of coordinates. The sole purpose of this
  ! is to reduce the filesize. This will typicaly reduce the number of
  ! coordinates by about a factor of 4. 

  call pointReduce(coords, iCoor, tol, uniqueCoords, link, nUnique)
  
  ! Now that we have produced the desired isosurface on each
  ! processor. Communicate the number of number of coordinates and the
  ! number of triangles each proc is going to send to the root:

  allocate(nPtsProc(nProc), nConnProc(nProc))
  nPtsProc(:) = 0_intType
  nConnProc(:) = 0_intType

  call  MPI_Allgather(nUnique, 1, mpi_integer4, nPtsProc, 1, mpi_integer4, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call  MPI_Allgather(iCoor/3, 1, mpi_integer4, nConnProc, 1, mpi_integer4, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  if (sum(nPtsProc) > 0) then 
     
     if (myid == 0) then ! Root proc does the writing
     
        ! Write a new zone:
        cgnsInd = fileIDs(sps) 
        cgnsBase = cgnsIsoSurfBases(sps) 

        ! Write the unstructured zone
        call cg_zone_write_f(cgnsInd, cgnsBase, isoName, (/sum(nPtsProc), sum(nConnProc), 0/), &
             Unstructured, cgnsZone, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        if(ierr /= CG_OK)                    &
             call returnFail("writeIsoSurface", &
             "Something wrong when calling cg_zone_write_f")
     end if

  else
     if (myid == 0) then
        ! We don't actually have an isosurface. We will create a zone
        ! that the same structure, but contains only a single triangle
        ! with all the coordinates at zero. This way the zone still
        ! exists and yields a uniform structure which can make
        ! processing easier

        ! Write a new zone:
        cgnsInd = fileIDs(sps) 
        cgnsBase = cgnsIsoSurfBases(sps) 

        call writeEmptyZone

     end if
     ! Don't forget to deallocate the stuff allocated so far:
     deallocate(nPtsProc, nConnProc, link, uniqueCoords, coords, weights, indices)
     return
  end if

  ! We need to keep track of the cumulative number of nodes since each
  ! proc has done its own ordering
  cumNodes = 0

  ! Communicate and write the coordinates
  do iproc=0, nProc-1

     dataOnProc: if (myid == iproc) then
        allocate(buffer(3*nPtsProc(iProc+1)))

        ! We will swap the order of the coordinates to packed format
        ! since this is what we need for CGNS
    
        do i=1,nPtsProc(iProc+1)
           buffer(i)                      = uniqueCoords(1, i)
           buffer(1*nPtsProc(iProc+1)+ i) = uniqueCoords(2, i)
           buffer(2*nPtsProc(iProc+1)+ i) = uniqueCoords(3, i)
        end do

     end if dataOnProc

     if (iproc .ne. 0) then
        tag = 13
        if (myid == 0) then
           ! allocate space for the recv
           allocate(buffer(3*nPtsProc(iProc+1)))

           call mpi_recv(buffer, nPtsProc(iProc+1)*3, sumb_real, iProc, tag, &
                sumb_comm_world, status, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if

        if (myid == iProc) then
           call mpi_send(buffer, nPtsProc(iProc+1)*3, sumb_real, 0, tag, &
                sumb_comm_world, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if
     end if

     if (myid == 0) then
        ! Now do partial writes on the root proc with points we've
        ! received from iProc
        if (nPtsProc(iProc+1) > 0) then
           
           call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                'CoordinateX', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                buffer(1:nPtsProc(iProc+1)), coordID, ierr)
           
           call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                'CoordinateY', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                buffer(nPtsProc(iProc+1)+1:2*nPtsProc(iProc+1)), coordID, ierr)
           
           call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                'CoordinateZ', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                buffer(2*nPtsProc(iProc+1)+1:3*nPtsProc(iProc+1)), coordID, ierr)

           if(ierr /= CG_OK)                    &
                call returnFail("writeIsoSurface", &
                "Something wrong when calling cg_coord_write_f")
           
           ! Increment by the number of nodes on this proc
           cumNodes = cumNodes + nPtsProc(iProc+1)
        end if
     end if

     ! Buffer was only allocated on root and current iProc
     if (myid == iProc .or. myid == 0) then
        deallocate(buffer)
     end if
  end do

  ! We need to keep track of the cumulative number of nodes since each
  ! proc has done its own ordering
  cumNodes = 0
  cumConn  = 0

  ! The partial write functionality is different between versions 2.5
  ! and 3.1, so we will just gather all the connectivities and do a
  ! final write at the end
  if (myid == 0) then
     allocate(allConn(3, sum(nConnProc)))
  endif
    
  ! Communicate and write the connectivity
  do iProc=0, nProc-1
     connOnProc: if (myid == iProc) then
        allocate(connBuffer(3,nConnProc(iProc+1)))
        do i=1,nConnProc(iProc+1)
           connBuffer(1, i) = link(3*i-2)
           connBuffer(2, i) = link(3*i-1)
           connBuffer(3, i) = link(3*i  )
        end do
     end if connOnProc

     ! Communication is only necessary if we are not dealing with root
     ! proc:
     if (iproc .ne. 0) then
        tag = 13
        if (myid == 0) then
           ! allocate space for the recv
           allocate(connBuffer(3,nConnProc(iProc+1)))
           call mpi_recv(connBuffer, nConnProc(iProc+1)*3, sumb_integer, iProc, tag, &
                sumb_comm_world, status, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if

        if (myid == iProc) then
           call mpi_send(connBuffer, nConnProc(iProc+1)*3, sumb_integer, 0, tag, &
                sumb_comm_world, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if
     endif

     if (myid == 0) then

        ! Copy into the allCon array and increment the received
        ! local connectivity by cummNodes
        
        do i=1,nConnProc(iProc+1)
           allConn(:, cumConn + i) = connBuffer(:, i) + cumNodes
        end do
   
        cumNodes = cumNodes + nPtsProc(iProc+1)
        cumConn  = cumConn  + nConnProc(iProc+1)    
     end if
        
     ! Buffer was only allocated on roto and iProc
     if (myid == iProc .or. myid == 0) then
        deallocate(connBuffer)
     end if
  end do

  ! Finally do the (full) connectivity write
  if (myid == 0) then
     ! Now write on root proc:
  
     ! Write just the connectively we have in buffer
     call cg_section_write_f(cgnsInd, cgnsBase, cgnsZone, "ELEM", TRI_3, &
          1, sum(nConnProc), 0, allConn, secID, ierr)
     if(ierr /= CG_OK)                    &
          call returnFail("writeIsoSurface", &
          "Something wrong when calling cg_section_partial_write_f")

     ! Also free allConn
     deallocate(allConn)
  end if
  
  ! Finally we have to write solution data for the iso surface
  ! iself. The main reason is that the same code is used for "slices"
  ! as well and in that case, you want to have other data interpolated
  ! on the "isoSurafce" (slice)

  ! Write the solution node:
  if (myid == 0) then
     call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, "isoSurfSolution", Vertex, solID, ierr)
     if(ierr /= CG_OK)                    &
          call returnFail("writeIsoSurface", &
          "Something wrong when calling cg_sol_write_f")
  end if

  ! Make the buffer large enough
  allocate(buffer(maxval(nPtsProc)))

  ! Loop over variables to write:
  do iVar=1,nIsoSurfVar

     ! We will reuse the same code as was used for computing the value
     ! onwhich we did the interpolation. However, set 'zero' for the
     ! isovalue such that we get the true value back
     
     call computeIsoVariable(isoSurfSolNames(iVar), sps, zero) 
    
     ! Set points to first block:
     call setPointers(1, 1, sps)
     curBlock = 1
     fn => flowDoms(1, 1, sps)%fn
     
     do i=1,iCoor

        ! If we've switched blocks, reset points. This stil only calls
        ! setPointer nDom times since there are at most that many
        ! switches
        if (indices(1, i) /= curBlock) then
           call setPointers(indices(1, i), 1, sps)
           curBlock = indices(1, i)
           fn => flowDoms(curBlock, 1, sps)%fn
        end if
        
        ! Computing interpolated value is easy using weights:
        ! indices on x
        buffer(link(i)) = weights(1, i) * fn(indices(2,i), indices(3,i), indices(4,i)) + &
             weights(2, i) * fn(indices(5,i), indices(6,i), indices(7,i))
     end do
     
     cumNodes = 0
     ! Communicate and write the solutions
     do iproc=0, nProc-1
        
        if (iproc .ne. 0) then
           tag = 13
           if (myid == 0) then
              call mpi_recv(buffer, nPtsProc(iProc+1), sumb_real, iProc, tag, &
                   sumb_comm_world, status, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end if
           
           if (myid == iProc) then
              call mpi_send(buffer, nPtsProc(iProc+1), sumb_real, 0, tag, &
                   sumb_comm_world, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end if
        end if
        
        if (myid == 0) then
           ! Now do partial writes on the root proc with points we've
           ! received from iProc
           if (nPtsProc(iProc+1) > 0) Then
              call cg_field_partial_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
                   isoSurfSolNames(iVar), cumNodes+1, cumNodes + nPtsProc(iProc+1), &
                   buffer, fieldID, ierr)
              
              if(ierr /= CG_OK)                    &
                   call returnFail("writeIsoSurface", &
                   "Something wrong when calling cg_field_partial_write_f")
              
              ! Increment by the number of nodes on this proc
              cumNodes = cumNodes + nPtsProc(iProc+1)
           end if
        end if
     end do
  end do
  ! Everyone deallocs buffer
  deallocate(buffer)

  ! Clear up temporary allocatable data.
  deallocate(nPtsProc, nConnProc)
  deallocate(coords, uniqueCoords, link)
  deallocate(weights, indices)

contains

  subroutine writeEmptyZone
    
    call cg_zone_write_f(cgnsInd, cgnsBase, isoName, (/3, 1, 0/), &
         Unstructured, cgnsZone, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f

    call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
         'CoordinateX', (/zero, zero, zero/), coordID, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f

    call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
         'CoordinateY', (/zero, zero, zero /), coordID, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f

    call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
         'CoordinateZ', (/zero, zero, zero/), coordID, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f

    call cg_section_write_f(cgnsInd, cgnsBase, cgnsZone, "ELEM", TRI_3, &
         1, 1, 0, (/1, 2, 3/), secID, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, "isoSurfSolution", Vertex, solID, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    do iVar = 1, nIsoSurfVar
       call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, isoSurfSolNames(iVar), &
            (/zero, zero, zero/), fieldID, ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
    end do
  end subroutine writeEmptyZone

end subroutine writeIsoSurface

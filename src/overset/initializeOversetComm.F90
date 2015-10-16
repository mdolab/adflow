subroutine initializeOversetComm

  ! This routine setups the data structures needed for overset
  ! communication. Specificially it creates the required petsc scatter
  ! context.
  use communication
  use overset
  use blockPointers
  use adjointVars
  use inputOverset
  implicit none
  ! Working Variables
  integer(kind=intType) :: ii, jj, kk, i, j, k, nn, ierr, tag1, tag2, dest
  integer(kind=intType) :: iDom, nFirstHalos, counter, donorID, iFringe
  integer(kind=intType) :: nFringeProc, nDonorsPerCell
  integer(kind=intType), dimension(:), allocatable :: offsetBlock
  integer(kind=intType), dimension(:), allocatable :: fringesProc
  integer(kind=intType), dimension(:), allocatable :: nFringePerBlock
  integer(kind=intType), dimension(:), allocatable :: intBuffer
  real(kind=realType), dimension(:), allocatable :: realBuffer
  integer(kind=intType) :: nFringeOnBlock, iStart, iEnd
  integer status(MPI_STATUS_SIZE) 
  ! Number of donors depends on the type of interpolation
  if (oversetInterpolation == linear) then 
     nDonorsPerCell = 8
  else
     nDonorsPerCell = 27
  end if

  ! First we determine the total number of fringes on each
  ! block. Currently on the root processor with the oBlocks knows this
  ! information. 

  ! Let everyone know the number of fringes from each block
  allocate(nFringePerBlock(nDomTotal))

  if (myid == 0) then 
     do iDom=1, nDomTotal
        oBlocks(iDom)%nFringe = 0
        do k=2, oBlocks(iDom)%kl
           do j=2, oBlocks(iDom)%jl
              do i=2, oBlocks(iDom)%il
                 if (oBlocks(iDom)%iBlank(i,j,k) == -1) then 
                    oBlocks(iDom)%nFringe = oBlocks(iDom)%nFringe + 1

                    ! Also check that we actuall have info for
                    ! it...totherwise this is an orphan that we can't
                    ! deal with yer. 

                    if (oBlocks(iDom)%donors(i, j, k)%donorProcID == -1) then 
                       print *,'We have an orphan. Cannot do this yet.'
                       print *,'Bad block, cell:', iDom, i, j, k 
                       stop
                    end if
                 end if
              end do
           end do
        end do
        nFringePerBlock(iDom) = oBlocks(iDom)%nFringe
     end do
  end if

  ! Broadcast nFringePerBlock to everyone
  call mpi_bcast(nFringePerBlock, nDomTotal, sumb_integer, 0, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Now go back through the oBlocks, pack up donorBlockID, ind and
  ! frac and fire it to the correct processor.
  if (myid == 0) then 
     do iDom=1, nDomTotal

        if (oBlocks(iDom)%procID /= 0) then 

           ! Allocate buffer space, pack and send. 
           allocate(realBuffer(nFringePerBlock(iDom)*3), &
                intBuffer(nFringePerBlock(iDom)*8))

           ii = 0
           do k=2, oBlocks(iDom)%kl
              do j=2, oBlocks(iDom)%jl
                 do i=2, oBlocks(iDom)%il
                    if (oBlocks(iDom)%iBlank(i,j,k) == -1) then 

                       ! Frac for the donor
                       realBuffer(ii*3+1) = oBlocks(iDom)%donors(i,j,k)%frac(1)
                       realBuffer(ii*3+2) = oBlocks(iDom)%donors(i,j,k)%frac(2)
                       realBuffer(ii*3+3) = oBlocks(iDom)%donors(i,j,k)%frac(3)

                       ! Block/index for the donor
                       intBuffer(ii*8+1) =  oBlocks(iDom)%donors(i,j,k)%donorBlockID
                       intBuffer(ii*8+2) =  oBlocks(iDom)%donors(i,j,k)%ind(1)
                       intBuffer(ii*8+3) =  oBlocks(iDom)%donors(i,j,k)%ind(2)
                       intBuffer(ii*8+4) =  oBlocks(iDom)%donors(i,j,k)%ind(3)

                       ! Block/index for the receiver
                       intBuffer(ii*8+5) =  oBlocks(iDom)%localBlockID
                       intBuffer(ii*8+6) =  i
                       intBuffer(ii*8+7) =  j
                       intBuffer(ii*8+8) =  k

                       ii = ii + 1
                    end if
                 end do
              end do
           end do

           tag1 = 2*(iDom-1) + 1
           tag2 = tag1 + 1
           dest = oBlocks(iDom)%procID

           call mpi_send(intBuffer, size(intBuffer), sumb_integer, dest, tag1, &
                sumb_comm_world, ierr)

           call mpi_send(realBuffer, size(realBuffer), sumb_real, dest, tag2, &
                sumb_comm_world, ierr)

           deallocate(realBuffer, intBuffer)

        else

           ! Its a local block so we can allocate and copy
           allocate(flowDoms(iDom, 1, 1)%fringeFrac(3, nFringePerBlock(iDom)), &
                flowDoms(iDom, 1, 1)%fringeIndices(8, nFringePerBlock(iDom)))

           flowDoms(iDom, 1, 1)%nFringe = nFringePerBlock(iDom)

           ! Set pointers to make it easier to read
           fringeFrac => flowDoms(iDom, 1, 1)%fringeFrac
           fringeIndices => flowDoms(iDom, 1, 1)%fringeIndices

           ii = 0
           do k=2, oBlocks(iDom)%kl
              do j=2, oBlocks(iDom)%jl
                 do i=2, oBlocks(iDom)%il
                    if (oBlocks(iDom)%iBlank(i,j,k) == -1) then 
                       ii = ii + 1
                       fringeFrac(:, ii) = oBlocks(iDom)%donors(i,j,k)%frac
                       fringeIndices(1  , ii) = oBlocks(iDom)%donors(i,j,k)%donorBlockID
                       fringeIndices(2:4, ii) = oBlocks(iDom)%donors(i,j,k)%ind

                       fringeIndices(5  , ii) = oBlocks(iDom)%localBlockID
                       fringeIndices(6:8, ii) = (/i, j, k/)

                    end if
                 end do
              end do
           end do
        end if
     end do

  else ! Not root...need to receive.

     do nn=1, nDom
        call setPointers(nn, 1, 1)

        iDom = cumDomProc(myid) + nn
        allocate(realBuffer(nFringePerBlock(iDom)*3), &
             intBuffer(nFringePerBlock(iDom)*8))

        tag1 = 2*(iDom-1) + 1
        tag2 = tag1 + 1

        call MPI_Recv(intBuffer, size(intBuffer), sumb_integer, 0, tag1, &
             sumb_comm_world, status, ierr)

        call MPI_Recv(realBuffer, size(realBuffer), sumb_real, 0, tag2, &
             sumb_comm_world, status, ierr)

        ! Allocate space in flowDoms for this
        allocate(flowDoms(nn, 1, 1)%fringeFrac(3, nFringePerBlock(iDom)), &
             flowDoms(nn, 1, 1)%fringeIndices(8, nFringePerBlock(iDom)))
        flowDoms(nn, 1, 1)%nFringe = nFringePerBlock(iDom)

        ! Set pointers to make it easier to read
        fringeFrac => flowDoms(nn, 1, 1)%fringeFrac
        fringeIndices => flowDoms(nn, 1, 1)%fringeIndices

        do ii=1, nFringePerBlock(iDom)
           fringeFrac(:, ii) = realBuffer(3*ii-2:3*ii)
           fringeIndices(:, ii) = intBuffer(8*ii-7:8*ii)
        end do

        ! Free the buffers
        deallocate(intBuffer, realBuffer)
     end do
  end if

  nFringeProc = 0
  do nn=1, nDom
     ! Count up the total number of fringes (times 8) I need. This is
     ! quite inefficient but we'll fix later
     nFringeProc = nFringeProc + nDonorsPerCell*flowDoms(nn, 1, 1)%nFringe
  end do
 
  ! One last thing...we need the single-level halo based onset per proc 
  allocate(offsetBlock(nDomTotal))
  offsetBlock(1) = 0
  do nn=2, nDomTotal
     offsetBlock(nn) = offsetBlock(nn-1) + &
          (dims(1, nn-1) + 1) * (dims(2, nn-1) + 1) * (dims(3, nn-1) + 1)
  end do


  ! Allocate spaces for the fringes indices
  allocate(fringesProc(nFringeProc))

  nFirstHalos = 0
  counter = 0

  ! Copy in the indices each block needs. 
  do nn=1, nDom
     call setPointers(nn, 1, 1)
     do iFringe=1, nFringe

        ! What we will do now is back out global index of the donor
        ! cells
        donorID = fringeIndices(1, iFringe)
        i = fringeIndices(2, iFringe)
        j = fringeIndices(3, iFringe)
        k = fringeIndices(4, iFringe)
        
        ! Note that the -1 here makes the ii,jj,kk effectively start
        ! at 0. This is necessary since the fringesProc must be
        ! zero-based for PETSc. 
        do kk=k-1, k
           do jj=j-1, j
              do ii=i-1, i
                 counter = counter + 1
                 fringesProc(counter) = offsetBlock(donorID) + &
                      kk*( (dims(1, donorID)+1)*(dims(2, donorID) + 1)) + &
                      jj*( dims(1, donorID)+1 ) + ii
              end do
           end do
        end do
     end do
     
     ! Also count up the number of (local) cells plus first halos here
     ! as well.
     nFirstHalos = nFirstHalos + ie*je*ke
  end do

  ! Create the two arrays needed for the overset comm. The first one
  ! is just a scalar version of w. That is just the number of local cells. 
  call VecCreateMPI(sumb_comm_world, nFirstHalos, PETSC_DETERMINE, oversetDonors, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! The second is equal to the nDonorsPerCell*the number of fringes we have. 
  call VecCreateMPI(sumb_comm_world, nFringeProc, PETSC_DETERMINE, oversetFringes, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! The only other thing we need to do is create the actual scatter context. 
  call ISCreateGeneral(sumb_comm_world, nFringeProc, fringesProc, PETSC_COPY_VALUES, IS1, ierr)
  call ECHK(ierr, __FILE__, __LINE__)  

  ! We are putting all values into oversetFringes in order, hence the
  ! strided index set for IS2 
  call VecGetOwnershipRange(oversetFringes, iStart, iEnd, ierr)
  call ISCreateStride(sumb_comm_world, nFringeProc, iStart, 1, IS2, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  call VecScatterCreate(oversetDonors, IS1, oversetFringes, IS2, oversetScatter, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Clean up the temporary IS1
  deallocate(fringesProc, offsetBlock, nFringePerBLock)

end subroutine initializeOversetComm

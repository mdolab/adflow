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
  integer(kind=intType) :: nFringeProc, nDonorsPerCell, level, sps
  integer(kind=intType), dimension(:), allocatable :: offsetBlock
  integer(kind=intType), dimension(:), allocatable :: fringesProc
  integer(kind=intType) :: nFringeOnBlock, iStart, iEnd
  integer status(MPI_STATUS_SIZE) 

  ! Explictly set level and sps to 1. This will be removed in the future.
  level = 1
  sps = 1 

  ! Number of donors depends on the type of interpolation
  if (oversetInterpolation == linear) then 
     nDonorsPerCell = 8
  else
     nDonorsPerCell = 27
  end if

  ! My own proc flowDoms already knows my fringe cells donors info
  ! Now just update the 'nFringe, 'fringeFrac' and 'fringeIndices' info.

  ! First find out nFringe
  loop_nDom: do nn=1, nDom
     iDom = nn + cumDomProc(myid)
     call setPointers(nn, level, sps)

     nFringe = 0 ! = flowDoms(nn, level, 1)%nFringe = 0 
     do k=2, kl
        do j=2, jl
           do i=2, il
              
              if (iBlank(i, j, k) == -1) then
                 nFringe = nFringe + 1

                 ! Also check that we actuall have info for
                 ! it...totherwise this is an orphan that we can't
                 ! deal with yer. 
                 if ( donors(i, j, k)%donorBlockId == -1 ) then

                    print *,'We have an orphan. Cannot do this yet.'
                    print *,'Bad global block, local cell:', iDom, i, j, k 
                    stop
                 end if
              end if 

           end do
        end do
     end do
     flowDoms(nn, level, sps)%nFringe = nFringe

     ! Allocate 'fringeFrac' and 'fringeIndices'
     ! fringeFrac(3, nFringe)      : donor cell fractions
     ! fringeIndices( 1   , nFringe) : donor cell globalBlockId
     ! fringeIndices( 2:4 , nFringe) : donor cell indices
     ! fringeIndices( 5:12, nFringe) : donor cell global cells i
     ! fringeIndices(13:15, nFringe) : current or receiver/fringe cell indices 
     
     allocate(flowDoms(nn, level, sps)%fringeFrac(3, nFringe), &
              flowDoms(nn, level, sps)%fringeIndices(8, nFringe))
              !flowDoms(nn, level, sps)%fringeIndices(15, nFringe))

     fringeFrac => flowDoms(nn, level, sps)%fringeFrac
     fringeIndices => flowDoms(nn, level, sps)%fringeIndices
  
     ! Now, actually save fringeFrac and fringeIndices from flowDoms()%donors
     ii = 0
     do k=2, kl
        do j=2, jl
           do i=2, il
              
              if (iBlank(i, j, k) == -1) then
                 ii = ii + 1

                 fringeFrac(:, ii)      = donors(i, j, k)%frac

                 fringeIndices( 1   , ii) = donors(i, j, k)%donorBlockID !globalBlockId of donor
                 fringeIndices( 2:4 , ii) = donors(i, j, k)%ind !donor indices

                 !a fringeIndices( 5:12, ii) = donors(i, j, k)%gInd !donor global cell indices

                 !a fringeIndices(13:15, ii) = (/i, j, k/) !receiver/fringe indices

                 fringeIndices(  5, ii)   = nn !localBlockId
                 fringeIndices(6:8, ii)   = (/i, j, k/) !receiver/fringe indices

              end if
     
           end do
        end do
     end do

  end do loop_nDom

  nFringeProc = 0
  do nn=1, nDom
     ! Count up the total number of fringes (times 8) I need. 
     nFringeProc = nFringeProc + nDonorsPerCell*flowDoms(nn, level, sps)%nFringe
  end do

  ! One last thing...we need the single-level halo based onset per proc 
  allocate(offsetBlock(nDomTotal))
  offsetBlock(1) = 0
  do nn=2, nDomTotal
     offsetBlock(nn) = offsetBlock(nn-1) + &
          (dims(1, nn-1) + 1) * (dims(2, nn-1) + 1) * (dims(3, nn-1) + 1)
  end do

  call MPi_barrier(sumb_comm_world, ierr)

  ! Allocate spaces for the fringes indices
  allocate(fringesProc(nFringeProc))

  nFirstHalos = 0
  counter = 0
  
  ! Copy in the indices each block needs. 
  do nn=1, nDom
     call setPointers(nn, level, sps)
     !fringeIndices => flowDoms(nn, level, sps)%fringeIndices

     do iFringe=1, nFringe  !flowDoms(nn, level, sps)%nFringe

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
  deallocate(fringesProc, offsetBlock)

end subroutine initializeOversetComm

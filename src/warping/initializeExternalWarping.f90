!
! ***********************************
! *  File: initializeExternalWarping.f90
! *  Author: Gaetan Kenway
! *  Started: 05-29-2010
! *  Modified: 05-29-2010
! ***********************************

subroutine initializeExternalWarping(ndofcgns)

  ! The purpose of this routine is to create the PETSc index sets and
  ! the vector scatter operators used to transfer the grid between the
  ! original cgns ordering and the (possibly) block split ordering.

  use precision 
  use blockPointers
  use communication
  use monitor
  use cgnsGrid
  use warpingPetsc
  use mddata
  use mddatalocal
  use bctypes
  implicit none

  integer(kind=intType) :: ndofcgns

  ! Local Variables

  integer(kind=intType) :: nn,ndof
  integer(kind=intType) ,allocatable,dimension(:) :: indices
  integer(kind=intType) ,allocatable,dimension(:) :: indices2
  integer(kind=intType) :: il_cgns,jl_cgns,kl_cgns
  integer(kind=intType) :: i,j,k,ii,mm
  integer(kind=intType) :: indx,indy,indz
  integer(kind=intTYpe) :: iBeg,iEnd,jBeg,jEnd,kBeg,kEnd
  integer(kind=intType) :: ist,ien,jst,jen,kst,ken
  integer(kind=intType) :: il_cg,jl_cg,kl_cg
  integer(kind=intType) :: ierr
  integer(kind=intType) ,allocatable,dimension(:) :: dof_offset
  integer(kind=intType) ,allocatable,dimension(:) :: temp

  logical :: storesubface
  ! -------------------------------------------------------------
  !                   Grid Scatter Context 
  ! -------------------------------------------------------------

  !Number of sumb dof per proc and cumulative version
  allocate(temp(nProc))
  ndof = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     ndof = ndof + il*jl*kl*3
  end do
  allocate(indices(ndof))

  call mpi_allgather(ndof,1,sumb_integer,temp,1,sumb_integer, &
       sumb_comm_world,ierr)
  allocate(cumdofproc(0:nproc))
  cumdofproc(0) = 0
  do i=1,nProc
     cumdofproc(i) = cumdofproc(i-1) + temp(i)
  end do
  deallocate(temp)

  ! Determine the local block offset (sumb)
  allocate(cumdofblock(ndom+1))
  cumdofblock(1) = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     cumdofblock(nn+1) = cumdofblock(nn) + il*jl*kl*3
  end do


  ! Determine the global offset list for the cgns blocks (cumulative)

  allocate(dof_offset(cgnsNDom))
  dof_offset(1) = 0
  do nn=2,cgnsNDom
     dof_offset(nn) = dof_offset(nn-1) + cgnsDoms(nn-1)%il*cgnsDoms(nn-1)%jl*cgnsDoms(nn-1)%kl*3
  end do

  ii = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
 
     do k=1,kl
        do j=1,jl
           do i=1,il
              il_cg = cgnsDoms(nbkGlobal)%il
              jl_cg = cgnsDoms(nbkGlobal)%jl
              kl_cg = cgnsDoms(nbkGlobal)%kl

              ii = ii + 1

              indx = iBegOr + i - 1
              indy = jBegOr + j - 1
              indz = kBegOr + k - 1

              indices(ii*3-2) = dof_offset(nbkGlobal)  + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3

              indices(ii*3-1) = dof_offset(nbkGlobal)  + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3 + 1

              indices(ii*3 ) = dof_offset(nbkGlobal)   + &
                   (indz-1)*jl_cg*il_cg*3 + &
                   (indy-1)*il_cg*3       + &
                   (indx-1)*3 + 2

           end do ! k loop
        end do ! j loop
     end do ! i loop
  end do ! domain loop

  ! Create the two index sets
  call ISCreateGeneral(sumb_comm_world,ndof,indices,IScgns,ierr)
  call ISCreateStride (sumb_comm_world,ndof,cumdofproc(myid),1_intType,ISsumb,ierr)

  ! Indices no longer required
  deallocate(indices)

  ! Setup the cgnsGrid and sumbGrid Petsc Vectors
  ! cgnsGrid Vector
  call VecCreate(sumb_comm_world,cgnsGridVec,ierr)
  call VecSetType(cgnsGridVec,"mpi",ierr)
  call VecSetSizes(cgnsGridVec,ndofcgns,PETSC_DETERMINE,ierr)

  ! sumbGrid Vector
  call VecCreate(sumb_comm_world,sumbGridVec,ierr)
  call VecSetType(sumbGridVec,"mpi",ierr)
  call VecSetSizes(sumbGridVec,ndof,PETSC_DETERMINE,ierr)

  ! Create the scatter context; cgnsTOsumbGrid is the actual scatter context
  call VecScatterCreate(cgnsGridVec,IScgns,sumbGridVec,ISsumb,cgnsTOsumbGrid,ierr)

  ! The index sets are no longer required
  call ISDestroy(ISsumb,ierr)
  call ISDestroy(IScgns,ierr)

  ! -------------------------------------------------------------
  !                Surface Force Scatter Context
  ! -------------------------------------------------------------

  ! We need to transfer  mdSumNSurfNodesLocal*3 values from this processor

  allocate(indices(mdSumNSurfNodesLocal*3))   ! -> These are CGNS indices
  allocate(indices2(mdSumNSurfNodesLocal*3))  ! -> These are SUMB indices

  ! Loop in the same way we do in the force calculation to get the
  ! original CGNS dof cooresponding to the point we're calculating the
  ! force for

  ii = 0

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     il_cg = cgnsDoms(nbkGlobal)%il
     jl_cg = cgnsDoms(nbkGlobal)%jl
     kl_cg = cgnsDoms(nbkGlobal)%kl

     bocos: do mm=1,nBocos
        ! Check if the data of this subface must be stored.
        storeSubface = .false.
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
      


           do k= knBeg(mm),knEnd(mm)
              do j= jnBeg(mm),jnEnd(mm)
                 do i= inBeg(mm),inEnd(mm)
                    ii = ii + 1

                    indx = iBegOr + i - 1
                    indy = jBegOr + j - 1
                    indz = kBegOr + k - 1

                    indices(ii*3-2) = dof_offset(nbkGlobal)  + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3

                    indices(ii*3-1) = dof_offset(nbkGlobal)  + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3 + 1

                    indices(ii*3 ) = dof_offset(nbkGlobal)   + &
                         (indz-1)*jl_cg*il_cg*3 + &
                         (indy-1)*il_cg*3       + &
                         (indx-1)*3 + 2

                    indices2(ii*3-2) = cumdofproc(myID) + cumdofblock(nn) + &
                         (k-1)*jl*il*3 + &
                         (j-1)*il*3    + &
                         (i-1)*3

                    indices2(ii*3-1) = cumdofproc(myID) + cumdofblock(nn) + &
                         (k-1)*jl*il*3 + &
                         (j-1)*il*3    + &
                         (i-1)*3 + 1

                    indices2(ii*3  ) = cumdofproc(myID) + cumdofblock(nn) + &
                         (k-1)*jl*il*3 + &
                         (j-1)*il*3    + &
                         (i-1)*3 + 2
                 end do
              end do
           end do
        end if
     end do bocos
  end do domains

  ! Create the two index sets
  call ISCreateGeneral(sumb_comm_world,mdSumNSurfNodesLocal*3,indices ,IScgns,ierr)
  call ISCreateGeneral(sumb_comm_world,mdSumNSurfNodesLocal*3,indices2,ISsumb,ierr)

  ! Indices no longer required
  deallocate(indices,indices2)

  ! Create the scatter context; cgnsTOsumbGrid is the actual scatter context
  call VecScatterCreate(sumbGridVec,ISsumb,cgnsGridVec,IScgns,sumbTOcgnsForce,ierr)

  ! The index sets are no longer required
  call ISDestroy(ISsumb,ierr)
  call ISDestroy(IScgns,ierr)

  deallocate(dof_offset)

end subroutine initializeExternalWarping




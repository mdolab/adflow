!
!      ******************************************************************
!      *                                                                *
!      * determineClusterArea determine the average cell surface area   *
!      * for all blocks in a particular cluster. This is used for       *
!      * determine blanking preference for overlapping surface cells.   *
!      *                                                                *
!      ******************************************************************

subroutine determineClusterAreas

  use BCTypes
  use blockPointers
  use communication
  use overset
  implicit none

  ! Working
  integer(kind=intType) :: i, j, mm, nn, clusterID, ierr, nPts, nCells
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii
  integer(kind=intType) :: lower_left, lower_right, upper_left, upper_right
  real(kind=realType), dimension(:), allocatable :: localAreas
  integer(kind=intType), dimension(:), allocatable :: localCount, globalCount
  real(kind=realType), dimension(:, :), allocatable :: pts
  real(kind=realType) :: fact , v1(3), v2(3), sss(3), da

  logical :: isWallType

  if (allocated(clusterAreas)) then 
     ! We only ever do this once!
     return
  end if

  allocate(clusterAreas(nClusters), localAreas(nClusters), &
       localCount(nClusters), globalCount(nClusters))
  
  localAreas = zero
  localCount = 0

  call getForceSize(nPts, nCells)
  allocate(pts(3, nPts))
  call getForcePoints(pts, nPts, 1_intType)

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)

     clusterID = clusters(cumDomProc(myid) + nn)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        if( isWalltype(BCType(mm))) then            
              
           ! Store the cell range of the subfaces a bit easier.
           ! As only owned faces must be considered the nodal range
           ! in BCData must be used to obtain this data.
           
           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
           
           ! Compute the dual area at each node. Just store in first dof
           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 
                 ! Compute Normal
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 v1(:) = pts(:, upper_right) -   pts(:, lower_left)
                 v2(:) = pts(:, upper_left) -  pts(:, lower_right)

                 ! Cross Product
                 call cross_prod(v1, v2, sss)
                 da = fourth*(sss(1)**2 + sss(2)**2 + sss(3)**2)
                 localAreas(clusterID) = localAreas(clusterID) + da
                 localCount(clusterID) = localCount(clusterID) + 1
              end do
           end do
           
           ! Note how iBeg,iBeg is defined above... it is one MORE
           ! then the starting node (used for looping over faces, not
           ! nodes)
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)
           
        end if
     end do bocos
  end do domains

  ! All reduce sum for the localAreas to get clusterAreas and
  ! localCount to get globalCount

  call mpi_allreduce(localAreas, clusterAreas, nClusters, sumb_real, &
       MPI_SUM, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  call mpi_allreduce(localCount, globalCount, nClusters, sumb_integer, &
       MPI_SUM, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Final get the average global area
  do i=1, nClusters
     clusterAreas(i) = clusterAreas(i)/max(globalCount(i), 1)
  end do
    
  deallocate(localAreas, localCount, globalCount)
end subroutine determineClusterAreas

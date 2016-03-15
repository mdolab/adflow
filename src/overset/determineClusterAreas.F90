!
!      ******************************************************************
!      *                                                                *
!      * determineClusterArea determine the average cell surface area   *
!      * for all blocks in a particular cluster. This is used for       *
!      * determine blanking preference for overlapping surface cells.   *
!      *                                                                *
!      ******************************************************************

subroutine determineClusterAreas(clusters, N, totalClusters)

  use BCTypes
  use blockPointers
  use communication
  use overset
  implicit none

  ! Input/output variables
  integer(kind=intType), dimension(N), intent(in) :: clusters
  integer(kind=intType), intent(in) :: N, totalClusters

  ! Working
  integer(kind=intType) :: i, j, mm, nn, clusterID, ierr
  real(kind=realType), dimension(:), allocatable :: localAreas
  integer(kind=intType), dimension(:), allocatable :: localCount, globalCount
  real(kind=realType), dimension(:, :, :), pointer :: ssi
  if (allocated(clusterAreas)) then 
     ! We only ever do this once!
     return
  end if

  allocate(clusterAreas(totalClusters), localAreas(totalClusters), &
       localCount(totalClusters), globalCount(totalClusters))
  
  localAreas = zero
  localCount = 0

  do nn=1, nDom
     call setPointers(nn, 1_intType, 1_intType)
     clusterID = clusters(cumDomProc(myid) + nn)
     do mm=1, nBocos

        if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsoThermal) then 
           
           select case (BCFaceID(mm))
           case (iMin)
              ssi => si(1,:,:,:)
           case (iMax)
              ssi => si(il,:,:,:)
           case (jMin)
              ssi => sj(:,1,:,:)
           case (jMax)
              ssi => sj(:,jl,:,:)
           case (kMin)
              ssi => sk(:,:,1,:)
           case (kMax)
              ssi => sk(:,:,kl,:)
           end select
           
           ! Loop over the owned cells
           do j = BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
              do i = BCData(mm)%inBeg+1, BCData(mm)%inEnd
                 ! No need to check for orientation since we just need
                 ! the area.
                 
                 localAreas(clusterID) = localAreas(clusterID) + &
                      sqrt(ssi(i, j, 1)**2 + ssi(i, j, 2)**2 + ssi(i, j, 3)**2)
                 localCount(clusterID) = localCount(clusterID) + 1
              end do
           end do
        end if
     end do
  end do

  ! All reduce sum for the localAreas to get clusterAreas and
  ! localCount to get globalCount

  call mpi_allreduce(localAreas, clusterAreas, totalClusters, sumb_real, &
       MPI_SUM, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  call mpi_allreduce(localCount, globalCount, totalClusters, sumb_integer, &
       MPI_SUM, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Final get the average global area
  do i=1, totalClusters
     clusterAreas(i) = clusterAreas(i)/max(globalCount(i), 1)
  end do
    
  deallocate(localAreas, localCount, globalCount)
end subroutine determineClusterAreas

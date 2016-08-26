subroutine determineClusterMarchDist
  use constants
  use blockPointers, only : nDom, x, xSeed, il, jl, kl
  use communication, only : sumb_comm_world, myid
  use overset, only : clusterMarchDist, cumDomProc, nClusters, clusters
  implicit none

  ! Working
  integer(kind=intType) :: i, j, k, nn, clusterID, ierr, nPts, nCells
  real(kind=realType), dimension(:), allocatable :: localDist
  real(kind=realType) :: xp(3)
  if (allocated(clusterMarchDist)) then 
     ! We only ever do this once!
     return
  end if

  allocate(clusterMarchDist(nClusters), localDist(nClusters))
  
  localDist = zero

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)

     clusterID = clusters(cumDomProc(myid) + nn)

     do k=2, kl
        do j=2, jl
           do i=2, il
              xp = eighth*(&
                   x(i-1, j-1, k-1, :) + &
                   x(i  , j-1, k-1, :) + &
                   x(i-1, j  , k-1, :) + &
                   x(i  , j  , k-1, :) + &
                   x(i-1, j-1, k  , :) + &
                   x(i  , j-1, k  , :) + &
                   x(i-1, j  , k  , :) + &
                   x(i  , j  , k  , :))
              if (xSeed(i, j, k, 1) < large) then 
                 localDist(clusterID) = max(localDist(clusterID), &
                      norm2(xp - xSeed(i, j, k, :)))
              end if
           end do
        end do
     end do
  end do domains
  ! All reduce sum for the localDist to get clusterMarchDist

  call mpi_allreduce(localDist, clusterMarchDist, nClusters, sumb_real, &
       MPI_SUM, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
    
  deallocate(localDist)

  clusterMarchDist = maxval(clusterMarchDist)

end subroutine determineClusterMarchDist

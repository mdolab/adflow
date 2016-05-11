subroutine computeCellWallPoint(level, sps)

  !  This routine is vastly more complex that it really should
  !  be. Essentially what we want to do is to determine the "wall
  !  point" for every cell. Essentially for every wall surface, we wan
  !  to record the coordintes of the wall surface cell center along
  !  cell centers eminating from the surface. The reason why this gets
  !  complex is that block can get cut in the off-wall direction which
  !  breaks the propagation. If this propatation isn't continued, the
  !  overset hole cut will be dependent on the block distribution and
  !  hense the numbe rof processors. 

  use blockPointers
  use BCTypes
  use communication
  use kdtree2_module
  use overset
  implicit none 

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps

  ! Working paramters
  integer(kind=intType) :: i, j, k, nn, cluster, ii, ind
  type(oversetWall), pointer :: wall
  logical, dimension(:), allocatable :: treeBuilt
  type(kdtree2_result), dimension(1) :: results  
  real(kind=realType) :: xp(3)

  ! We already have clusterWalls. Build the KD tree from the nodes
  ! only for the clusters we have.
  allocate(treeBuilt(nClusters))
  treeBuilt = .False.

  do nn=1, nDom
     call setPointers(nn, level, sps)
     cluster = clusters(cumDomProc(myid) + nn)
     wall => clusterWalls(cluster)

     ! If tree for this cluster is not built
     if (treeBuilt(cluster) .eqv. .False. .and. wall%nNodes > 0) then 
        
        ! Only build tree for real surface nodes. Copy these and the
        !  indices out. This is an overestimate of the size.
        allocate(wall%xPrimalCen(3, 1:wall%nNodes), wall%indPrimal(1:wall%nNodes))
  
        j = 0
        do i=1, wall%nNodes
           if (wall%ind(i) >= 0) then 
              j = j + 1
              wall%xPrimalCen(:, j) = wall%x(:, i)
              wall%indPrimal(j) = wall%ind(i)
           end if
        end do

        wall%tree => kdtree2_create(wall%xPrimalCen(:, 1:j), sort=.False.)
     end if

     if (.not. associated(flowDoms(nn, level, sps)%xSeed)) then 
        allocate(flowDoms(nn, level, sps)%XSeed(0:ib, 0:jb, 0:kb, 3))
        allocate(flowDoms(nn, level, sps)%wallInd(2:il, 2:jl, 2:kl))
        ! Manaully set the pointer for xSeed so we don't call
        ! setPointers again
        xSeed => flowDoms(nn, level, sps)%xSeed
        wallInd => flowDoms(nn, level, sps)%wallInd
     end if

     ! Initialize to large to indicate that nothing has been changed. 
     xSeed = large
     wallInd = -1

     if (wall%nNodes > 0) then 
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

                 call kdtree2_n_nearest(wall%tree, xp, 1, results)
                 
                 ! Need to store the value in xseed and wall ind
                 xseed(i, j, k, :) = wall%xPrimalCen(:, results(1)%idx)
                 wallInd(i, j, k) = wall%indPrimal(results(1)%idx)
                 
              end do
           end do
        end do
     end if
  end do

  ! Loop back over the blocks destroying the kd_trees as necessary
  do nn=1, nDom
     call setPointers(nn, level, sps)
     cluster = clusters(cumDomProc(myid) + nn)
     wall => clusterWalls(cluster)

     if (treeBuilt(cluster)) then 
        call kdtree2_destroy(wall%tree)
        deallocate(wall%xPrimalCen, wall%indPrimal)
     end if
  end do

  ! Exchange the xSeeds. Need these for building oBlock.
  do nn=1, nDom
     flowDoms(nn, level, sps)%realCommVars(1)%var => &
          flowDoms(nn, level, sps)%xSeed(:, :, :, 1)
     flowDoms(nn, level, sps)%realCommVars(2)%var => & 
          flowDoms(nn, level, sps)%xSeed(:, :, :, 2) 
     flowDoms(nn, level, sps)%realCommVars(3)%var => &    
          flowDoms(nn, level, sps)%xSeed(:, :, :, 3) 
  end do
  
  ! Run the generic halo exchange.
  call wHalo1to1RealGeneric(3, level, sps, commPatternCell_2nd, internalCell_2nd)

end subroutine computeCellWallPoint


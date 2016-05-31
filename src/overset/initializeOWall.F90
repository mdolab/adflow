subroutine initializeOWall(oWall, dualMesh, cluster)

  ! This routine builds the ADT tree for any wall surfaces for the
  ! block currently being pointed to by block Pointers.
  use overset
  use blockPointers
  use adtAPI
  use BCTypes
  use kdtree2_module
  implicit none 

  ! Input Params
  type(oversetWall), intent(inout) :: oWall
  logical, intent(in) :: dualMesh 
  integer(kind=intType), intent(in) :: cluster

  ! Working paramters
  integer(kind=intType) :: i, j, k, n, ii, jj, jjj, mm, ni, nj, nodeCount
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, nNodes, maxCells, nCells, iNode
  logical :: isWallType

  ! Set all the sizes for this block.
  oWall%il = il
  oWall%jl = jl
  oWall%kl = kl

  call getWallSize(nNodes, maxCells, dualMesh)

  oWall%nNodes = nNodes
  oWall%maxCells = maxCells
  oWall%cluster = cluster
  ! Allocate space for the x array and connectivity array. cellPtr is
  ! larger than necessary.
  allocate(oWall%x(3, nNodes), oWall%conn(4, maxCells), &
       oWall%cellPtr(maxCells), oWall%iBlank(maxCells), &
       oWall%delta(nNodes), oWall%nte(4, nNodes))
 
  ii = 0 ! Cumulative node counter
  jj = 0 ! Cumulative cell counter (with iblanks)
  jjj = 0 !Cumulative cell counter (without iblanks)
  nodeCount = 0

  do mm=1, nBocos
     if (isWallType(BCType(mm))) then 

        ! THIS IS SUPER IMPORTANT: It is absolutely critical that the
        ! wall be built *FROM THE DUAL MESH!!* IT WILL NOT WORK IF YOU
        ! USE THE PRIMAL MESH! The -1 for the node ranges below gives
        ! the extra '1' node for the mesh formed from the dual cells. 

        dualCheck: if (dualMesh) then 
           jBeg = BCData(mm)%jnBeg-1 ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg-1 ; iEnd = BCData(mm)%inEnd
           ! Now fill up the point array
           do j=jBeg, jEnd
              do i=iBeg, iEnd 
                 ii = ii +1
                 select case(BCFaceID(mm))
                 case(imin)
                    oWall%x(:,ii) = fourth*(x(1, i, j, :) + x(1, i+1, j, :) + &
                         x(1, i, j+1, :) + x(1, i+1, j+1, :))
                 case(imax)
                    oWall%x(:,ii) = fourth*(x(il, i, j, :) + x(il, i+1, j, :) + &
                         x(il, i, j+1, :) + x(il, i+1, j+1, :))
                 case(jmin) 
                    oWall%x(:,ii) = fourth*(x(i, 1, j, :) + x(i+1, 1, j, :) + &
                         x(i, 1, j+1, :) + x(i+1, 1, j+1, :))
                 case(jmax) 
                    oWall%x(:,ii) = fourth*(x(i, jl, j, :) + x(i+1, jl, j, :) + &
                         x(i, jl, j+1, :) + x(i+1, jl, j+1, :))
                 case(kmin) 
                    oWall%x(:,ii) = fourth*(x(i, j, 1, :) + x(i+1, j, 1, :) + &
                      x(i, j+1, 1, :) + x(i+1, j+1, 1, :))
                 case(kmax) 
                    oWall%x(:,ii) = fourth*(x(i, j, kl, :) + x(i+1, j, kl, :) + &
                         x(i, j+1, kl, :) + x(i+1, j+1, kl, :))
                 end select
              end do
           end do

           ! Fill up the conn array. Note that don't take the
           ! surface`normal direction (in or out) or the cell
           ! handed-ness into account...it is not necessary since we
           ! are just getting distance to the wall, which is
           ! independent of the orientation.
           
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0, nj-2
              do i=0, ni-2
                 jj = jj + 1
                 oWall%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                 oWall%conn(2, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                 oWall%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                 oWall%conn(4, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
              end do
           end do
           nodeCount = nodeCount + ni*nj

           ! We don't care about iBlank, cellPtr or delta for the dual
           ! mesh
           oWall%iBlank = 1
           oWall%cellPtr = 0
           oWall%delta = zero
        else ! Using the primal mesh
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           ! Now fill up the point array. Owned node loop.
           do j=jBeg, jEnd
              do i=iBeg, iEnd 
                 ii = ii +1
                 select case(BCFaceID(mm))
                 case(imin)
                    oWall%x(:,ii) = x(1, i, j, :)
                 case(imax)
                    oWall%x(:,ii) = x(il, i, j, :) 
                 case(jmin) 
                    oWall%x(:,ii) = x(i, 1, j, :) 
                 case(jmax) 
                    oWall%x(:,ii) = x(i, jl, j, :) 
                 case(kmin) 
                    oWall%x(:,ii) = x(i, j, 1, :) 
                 case(kmax) 
                    oWall%x(:,ii) = x(i, j, kl, :) 
                 end select
                 oWall%delta(ii) = BCData(mm)%deltaNode(i, j)
              end do
           end do

           ! Fill up the conn array being careful to *only* adding
           ! cells that are not already blanked. 
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           do j=0, nj-2
              do i=0, ni-2
                 jjj = jjj + 1
                 oWall%iBlank(jjj) = BCData(mm)%iblank(iBeg+i+1,jBeg+j+1)
                 if (oWall%iBlank(jjj) == 1) then 
                    jj = jj + 1
                    oWall%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                    oWall%conn(2, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                    oWall%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                    oWall%conn(4, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
                    oWall%cellPtr(jj) = jjj
                 end if
              end do
           end do
           nodeCount = nodeCount + ni*nj
        end if dualCheck
     end if
  end do

  ! Set the actual number of cells
  oWall%nCells = jj

  ! Build the tree itself.
  call buildSerialQuad(oWall%nCells, nNodes, oWall%x, oWall%conn, owall%ADT)

  ! Build the KDTree
  if (oWall%nNodes > 0) then 
     oWall%tree => kdtree2_create(oWall%x)
  end if
  
  ! Build the inverse of the connectivity, the nodeToElem array. 
  oWall%nte = 0
  do i=1, oWall%nCells
     do j=1, 4
        n = oWall%conn(j, i)
        inner:do k=1, 4
           if (oWall%nte(k, n) == 0) then 
              oWall%nte(k, n) = i
              exit inner
           end if
        end do inner
     end do
  end do

  ! Flag this wall as being allocated
  oWall%allocated = .True.

end subroutine initializeOWall

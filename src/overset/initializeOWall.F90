subroutine initializeOWall(oWall)

  ! This routine builds the ADT tree for any wall surfaces for the
  ! block currently being pointed to by block Pointers.
  use overset
  use blockPointers
  use adtAPI
  use BCTypes

  implicit none 

  ! Input Params
  type(oversetWall), intent(inout) :: oWall

  ! Working paramters
  integer(kind=intType) :: i, j, ii, jj, mm, ni, nj, nodeCount
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, nNodes, nCells

  ! Set all the sizes for this block.
  oWall%il = il
  oWall%jl = jl
  oWall%kl = kl

  call getWallSize(nNodes, nCells)

  oWall%nNodes = nNodes
  oWall%nCells = nCells

  ! Allocate space for the x array and connectivity array
  allocate(oWall%x(3, nNodes), oWall%conn(4, nCells))

  ii = 0
  jj = 0
  nodeCount = 0
  do mm=1, nBocos
     if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic .or. &
          BCType(mm) == NSWallIsoThermal) then 

        ! THIS IS SUPER IMPORTANT: It is absolutely critical that the
        ! wall be built *FROM THE DUAL MESH!!* IT WILL NOT WORK IF YOU
        ! USE THE PRIMAL MESH! The -1 for the node ranges below gives
        ! the extra '1' node for the mesh fromed from the dual cells. 

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
        ! surface`normal direction (in or out) or the cell handed-ness
        ! into account...it is not necessary since we are just getting
        ! distance to the wall, which is independent of the orientation. 

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
     end if
  end do
  
  ! Build the tree itself.
  call buildSerialQuad(nCells, nNodes, oWall%x, oWall%conn, owall%ADT)

  ! Flag this wall as being allocated
  oWall%allocated = .True.

end subroutine initializeOWall

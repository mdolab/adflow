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

        ! NODE Loop
        jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
        iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

        ! Now fill up the point array
        do j=jBeg, jEnd
           do i=iBeg, iEnd ! This is a node loop
              ii = ii +1
              select case(BCFaceID(mm))
              case(imin)
                 oWall%x(:,ii) = x(1,i,j,:)
              case(imax)
                 oWall%x(:,ii) = x(il,i,j,:)
              case(jmin) 
                 oWall%x(:,ii) = x(i,1,j,:)
              case(jmax) 
                 oWall%x(:,ii) = x(i,jl,j,:)
              case(kmin) 
                 oWall%x(:,ii) = x(i,j,1,:)
              case(kmax) 
                 oWall%x(:,ii) = x(i,j,kl,:)
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

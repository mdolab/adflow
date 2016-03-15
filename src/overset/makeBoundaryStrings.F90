subroutine makeBoundaryStrings(level, sps)

  use blockPointers
  use bctypes
  use communication
  implicit none
  
  ! Input Params
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(kind=intType) :: i, j, k, nn, mm, curElem, jj
  integer(kind=intType) :: iStart, iEnd, iSize, ierr, iProc, n(2), m(2)
  integer(kind=intType) :: below, above, left, right, nEdge, nNode, nUnique
  logical :: duplicateElement
  real(kind=realType), dimension(:, :, :), pointer :: xx
  integer(kind=intType), dimension(:, :), allocatable :: edgeConn, globalConn, nodeToElem
  integer(kind=intType), dimension(:), allocatable :: link
  real(kind=realType), dimension(:, :), allocatable :: edgeNodes, uniqueNodes, globalNodes
  integer(kind=intType), dimension(:), allocatable :: nEdgeProc, nNodeProc
  integer(kind=intType), dimension(:), allocatable :: tmp
  integer(kind=intType), dimension(:, :), pointer :: gcp
  logical, dimension(:), allocatable :: duplicated
  character(80) :: fileName, zoneName
  integer status(MPI_STATUS_SIZE) 

  interface
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none

       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduce

  end interface

  ! Loop over the wall faces counting up the edges that stradle a
  ! compute cell and a blanked (or interpolated) cell. 

  nEdge = 0
  domainLoop: do nn=1, nDom
     call setPointers(nn, level, sps)
          
     ! Push the surface iblank back to the volume:
     bocoLoop: do mm=1, nBocos
        wallType: if (BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           select case (BCFaceID(mm))
           case (iMin)
              gcp => globalCell(2, :, :)
           case (iMax)
              gcp => globalCell(il, :, :)
           case (jMin)
              gcp => globalCell(:, 2, :)
           case (jMax)
              gcp => globalCell(:, jl, :)
           case (kMin)
              gcp => globalCell(:, :, 2)
           case (kMax)
              gcp => globalCell(:, :, kl)
           end select

           ! Check the i-edges
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd       ! <------- Node loop
              do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd  ! <------- Face Loop

                 if (gcp(i+1, j+1) > 0 .and. gcp(i+1, j+2) > 0) then 
                    
                    below = max(BCData(mm)%iBlank(i, j  ), 0)
                    above = max(BCData(mm)%iBlank(i, j+1), 0)
                    
                    if ((below == 0 .and. above == 1) .or. &
                         (below == 1 .and. above == 0)) then 
                       nEdge = nEdge + 1
                    end if
                 end if
              end do
           end do
           
           ! Check the j-edges
           do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd   ! <------- Face loop
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd  ! <------- Node Loop
                 if (gcp(i+1, j+1) > 0 .and. gcp(i+2, j+1)> 0) then 

                    left = max(BCData(mm)%iBlank(i, j), 0)
                    right = max(BCData(mm)%iBlank(i+1,  j), 0)
                    
                    if ((left == 0 .and. right == 1) .or. &
                         (left == 1 .and. right == 0)) then 
                       nEdge = nEdge + 1
                    end if
                 end if
              end do
           end do
        end if wallType
     end do bocoLoop
  end do domainLoop
  
  ! Allocate the spae we need
  allocate(edgeConn(2, nEdge), edgeNodes(3, nEdge*2))

  ! And do the same thing again, but this time add the actual nodes
  ! and conn. 
  nEdge = 0
  domainLoop2: do nn=1, nDom
     call setPointers(nn, level, sps)
          
     ! Push the surface iblank back to the volume:
     bocoLoop2: do mm=1, nBocos
        wallType2: if (BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           select case (BCFaceID(mm))
           case (iMin)
              xx => x(1, :, :, :)
              gcp => globalCell(2, :, :)
           case (iMax)
              xx => x(il, :, :, :)
              gcp => globalCell(il, :, :)
           case (jMin)
              xx => x(:, 1, :, :)
              gcp => globalCell(:, 2, :)
           case (jMax)
              xx => x(:, jl, :, :)
              gcp => globalCell(:, jl, :)
           case (kMin)
              xx => x(:, :, 1, :)
              gcp => globalCell(:, :, 2)
           case (kMax)
              xx => x(:, :, kl, :)
              gcp => globalCell(:, :, kl)
           end select

           ! Check the i-edges
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd       ! <------- Node loop
              do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd  ! <------- Face Loop
                 if (gcp(i+1, j+1) > 0 .and. gcp(i+1, j+2) > 0) then 
                    below = max(BCData(mm)%iBlank(i, j), 0)
                    above = max(BCData(mm)%iBlank(i, j+1), 0)

                    if (below == 0 .and. above == 1) then 
                       nEdge = nEdge + 1

                       ! Make sure the real cell is on the LEFT of the
                       ! directed edge an don't forget pointer offset
                       edgeNodes(:, nEdge*2-1) = xx(i  , j+1, :)
                       edgeNodes(:, nEdge*2 ) = xx(i+1, j+1, :)
                       edgeConn(:, nEdge) = (/2*nEdge-1, 2*nEdge/)

                    else if (below == 1 .and. above == 0) then 
                       nEdge = nEdge + 1
                       edgeNodes(:, nEdge*2-1) = xx(i+1, j+1, :)
                       edgeNodes(:, nEdge*2 ) = xx(i  , j+1, :)
                       edgeConn(:, nEdge) = (/2*nEdge-1, 2*nEdge/)
                    end if
                 end if
              end do
           end do
           
           ! Check the j-edges
           do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd   ! <------- Face loop
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd  ! <------- Node Loop
                 if (gcp(i+1, j+1) > 0 .and. gcp(i+2, j+1)> 0)then 
                    left = max(BCData(mm)%iBlank(i, j), 0)
                    right = max(BCData(mm)%iBlank(i+1,  j), 0)
                    
                    if (left == 0 .and. right == 1) then 
                       nEdge = nEdge + 1
                       ! Again, make sure the real cell is on the LEFT
                       ! of the directed edge
                       edgeNodes(:, nEdge*2-1) = xx(i+1, j+1, :)
                       edgeNodes(:, nEdge*2  ) = xx(i+1, j  , :) 
                       edgeConn(:, nEdge) = (/2*nEdge-1, 2*nEdge/)

                    else if (left == 1 .and. right == 0) then 
                       nEdge = nEdge + 1
                       edgeNodes(:, nEdge*2-1) = xx(i+1, j  , :)
                       edgeNodes(:, nEdge*2  ) = xx(i+1, j+1, :) 
                       edgeConn(:, nEdge) = (/2*nEdge-1, 2*nEdge/)
                    end if
                 end if
              end do
           end do
        end if wallType2
     end do bocoLoop2
  end do domainLoop2

  ! Currently we have approximately twice as may nodes as we should
  ! have since most nodes are duplicated on each edge. To help reduce
  ! the amount of work the poor root processor has to do, we will
  ! unique-ify the nodes on each processor before we send then so it
  ! has fewer nodes to deal with. 

  allocate(uniqueNodes(3, nEdge*2), link(nEdge*2))
  call pointReduce(edgeNodes, nEdge*2, 1e-8, uniqueNodes, link, nUnique)

  ! Update the connectivity to use the new set of nodes:
  do i=1,nEdge
     edgeConn(1, i) = link(edgeConn(1, i))
     edgeConn(2, i) = link(edgeConn(2, i))
  end do
  
  ! Now let the root processor know how many nodes/elements my
  ! processor will be providing:
  allocate(nEdgeProc(0:nProc), nNodeProc(0:nProc))
  nEdgeProc(0) = 0
  nNodeProc(0) = 0

  call MPI_Gather(nEdge, 1, sumb_integer, nEdgeProc(1:nProc), 1, sumb_integer, 0, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  call MPI_Gather(nUnique, 1, sumb_integer, nNodeProc(1:nProc), 1, sumb_integer, 0, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  if (myid == 0) then 

     ! Before we can receive stuff, we need to determine the node
     ! off-sets such that the conn from the strings on each processor
     ! don't overlap. 
     
     do i=2, nProc
        ! The 0 and 1st entry of the nEdgeProc and nNodeProc arrays are already correct:
        nNodeProc(i) = nNodeProc(i) + nNodeProc(i-1)
        nEdgeProc(i) = nEdgeProc(i) + nEdgeProc(i-1)
     end do
        
     allocate(globalNodes(3, nNodeProc(nProc)), globalConn(2, nEdgeProc(nProc)))

     ! Put proc 0's node/elements in the global list:
     do i=1, nUnique
        globalNodes(:, i) = uniqueNodes(:, i)
     end do
     
     do i=1,nEdge
        globalConn(:, i) = edgeConn(:, i)
     end do
     
     nNode = nNodeProc(nProc)
     nEdge = nEdgeProc(nProc)

     ! Now receive from each of the other procs. 
     do iProc=1, nProc-1
        ! Check if this proc actually has anything to send:
        if ((nEdgeProc(iProc+1) - nEdgeProc(iProc))>0) then 
           iStart = nNodeProc(iProc) + 1
           iEnd =   nNodeProc(iProc+1)
           iSize = iEnd - iStart + 1

           call MPI_Recv(globalNodes(:, iStart:iEnd), iSize*3, sumb_real, iProc, iProc, &
                sumb_comm_world, status, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           iStart = nEdgeProc(iProc) + 1
           iEnd =   nEdgeProc(iProc+1)
           iSize = iEnd - iStart + 1
           call MPI_Recv(globalConn(:, iStart:iEnd), iSize*2, sumb_integer, iProc, iProc, &
                sumb_comm_world, status, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           ! Increment the conn we just received by the node offset:
           do i=iStart, iEnd
              globalConn(:, i) = globalConn(:, i) + nNodeProc(iProc-1)
           end do

        end if
     end do
  else
     ! Not root proc so send my stuff if we have anything:
     if (nEdge > 0) then 
        call MPI_Send(uniqueNodes, 3*nUnique, sumb_real, 0, myid, &
             sumb_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
        
        call MPI_Send(edgeConn, 2*nEdge, sumb_integer, 0, myid, &
             sumb_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end if
  end if
  
  if (myid == 0) then 
     ! This is where we start the real serial bottleneck where only
     ! the root procssor does stuff. First thing we will do is
     ! uniquify the nodes onces more to ensure that any duplicates
     ! that occured across boundaries are elminated.
     deallocate(uniqueNodes, link)
     allocate(uniqueNodes(3, nNode), link(nNode))

     call pointReduce(globalNodes, nNode, 1e-8, uniqueNodes, link, nUnique)

     ! Update the global connectivity to use the new set of nodes:
     do i=1, nEdge
        globalConn(1, i) = link(globalConn(1, i))
        globalConn(2, i) = link(globalConn(2, i))
     end do

     ! Next thing we have to do is produce the inverse of the
     ! connectivity...the nodeToElem array. Normally, each node should
     ! point to 1 element (at a boundary) or two elements for a normal
     ! part of a chain. There are however, some pathalogical cases
     ! that will result in more: 

     ! 1. The same edge may be added twice if it is on a boundary
     ! between two block and each add the same edge. 

     ! 2. Chains that are "just" touching such as middle node of the
     ! following: (double lines are chains)
     !
     !     +----------++------------+
     !     +          ||            |
     !     +  ib=1    /\   ib=0     |
     !     +          ||            |
     !     +===>>=====++=====<<=====+
     !     +          ||            |
     !     +  ib=0    \/   ib=1     |
     !     +          ||            |
     !     +----------++------------+
     ! 
     ! The center node has 4 edges eminating from it. It is part of
     ! the chain that enters from the left and cointinues to the top
     ! as well as the chain that enters from the right and continues
     ! to the bottom. 







     ! This algorithm should take linear time, but has lots of
     ! indirect addressing (welcome to unstructured grids). It looks
     ! like it might be quadratic, but there should never be more than
     ! 4 edges attached to each node. 
     
     allocate(nodeToElem(4, nUnique))
     nodeToElem = 0 ! Only non-zero entries will be valid

     allocate(tmp(nUnique), duplicated(nEdge))
     tmp = 0
     duplicated = .False.

     do i=1, nEdge
        ! Node numbers we're working with:
        n(1) = globalConn(1, i)
        n(2) = globalConn(2, i)

        ! For each node check which elements (if any) are already
        ! connected. We need to check them again the node numbers n1 and n2

        duplicateElement = .False.
        do jj=1,2
           
           do j=1, tmp(n(jj)) ! Loop over the element numbers already here:
              curElem = nodeToElem(n(jj), j)
              m(1) = globalConn(1, curElem)
              m(2) = globalConn(1, curElem)

              if (m(1) == n(1) .and. m(2) == n(2)) then 
                 duplicateElement = .True.
              else if(m(1) == n(2) .and. m(2) == n(1)) then 
                 ! Element exists, but it is the wrong order...don't
                 ! know what to do with this, probably an error or
                 ! maybe a corner case I haven't thought of.
                 call terminate("makeBoundaryString", "Inconsistnet duplicate edge.")
              end if
           end do
        end do
        
        if (.not. duplicateElement) then 
           do jj=1,2
              tmp(n(jj)) = tmp(n(jj)) + 1
              nodeToElem(tmp(n(jj)), n(jj)) = i
           end do
        else
           ! Well, we've figured out that this element is actually a
           ! duplicate so we'll make a note of that
           duplicated(i) = .True.
        end if
     end do
     deallocate(tmp)
  end if

  ! Write the global gap strings
  if (myid == 0) then 

     !write (fileName,"(a,I1,a)") "gapStrings_", myid, ".dat"
  
    open(unit=101,file='global_gapstrings.dat',form='formatted')
    write(101,*) 'TITLE = "Gap Strings Data" '
    write(101,*) 'Variables = "X", "Y", "Z"'

    write (101, *) "Zone"
    write (101,*) "Nodes = ", nNode, " Elements= ", nEdge, " ZONETYPE=FELINESEG"
    write (101,*) "DATAPACKING=POINT"
13  format (E14.6)
       
    do i=1, nNode
       ! Write the coordinates
       do j=1,3
          write(101,13, advance='no') globalNodes(j, i)
       end do
       write(101,"(1x)")
    end do
       
15     format(I5, I5)
    do i=1, nEdge
       write(101, 15) globalConn(1, i), globalConn(2, i)
    end do
 end if

end subroutine makeBoundaryStrings


subroutine writeWalls

  use communication
  use overset
  use constants
  use blockPointers
  use BCTypes
  implicit none

  character(80) :: fileName, zoneName
  integer(kind=intType) :: i, j, nn, iDom, iBeg, iEnd, jBeg, jEnd, mm
  real(kind=realType), dimension(:, :, :), pointer :: xx
 
  write (fileName,"(a,I1,a)") "wall_", myid, ".dat"

  open(unit=101,file=trim(fileName),form='formatted')
  write(101,*) 'TITLE = "mywalls"'
  write(101,*) 'Variables = "X", "Y", "Z", "CellIBlank"'

  do nn=1,nDom
     iDom = nn + cumDomProc(myid)
     call setPointers(nn, 1, 1)
     if (nBocos > 0) then 
        do mm=1, nBocos
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           if (BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              select case (BCFaceID(mm))
              case (iMin)
                 xx => x(1,:,:,:)
              case (iMax)
                 xx => x(il,:,:,:)
              case (jMin)
                 xx => x(:,1,:,:)
              case (jMax)
                 xx => x(:,jl,:,:)
              case (kMin)
                 xx => x(:,:,1,:)
              case (kMax)
                 xx => x(:,:,kl,:)
              end select

              write(zoneName, "(a,I1,a,I1)") "Zone", iDom, "_Proc_", myid
110           format('ZONE T=',a, " I=", i5, " J=", i5)
              write(101, 110), trim(zoneName), iEnd-iBeg+1, jEnd-jBeg+1
              write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"

13            format (E14.6)

              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    write(101, *) xx(i+1, j+1, 1)
                 end do
              end do

              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    write(101, *) xx(i+1, j+1, 2)
                 end do
              end do

              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    write(101, *) xx(i+1, j+1, 3)
                 end do
              end do

              do j=jBeg+1, jEnd
                 do i=iBeg+1, iEnd
                    write(101, *) BCData(mm)%iBlank(i, j)
                 end do
              end do
           end if
        end do
     else
        ! Write dummy zone
        write(zoneName, "(a,I1,a,I1)") "Zone", 0, "_Proc_", myid
        write(101, 110), trim(zoneName), 1, 1
        write (101,*) "DATAPACKING=POINT"
        write(101, *) zero, zero, zero, one, one
     end if
  end do
  close(101)
end subroutine writeWalls


! subroutine writeOWalls(oWalls, n)

!   use communication
!   use overset
!   use constants
!   use blockPointers
!   implicit none

!   type(oversetWall), intent(inout), dimension(n) :: oWalls
!   integer(kind=intType), intent(in) :: n
!   character(80) :: fileName, zoneName
!   integer(kind=intType) :: i, j, nn, iDom
  
  
!   write (fileName,"(a,I1,a)") "oWall_", myid, ".dat"

!   open(unit=101,file=trim(fileName),form='formatted')
!   write(101,*) 'TITLE = "OWalls"'
!   write(101,*) 'Variables = "X", "Y", "Z", "IBlank"'

!   do nn=1,nDom
!      iDom = nn + cumDomProc(myid)
!      if (oWalls(iDom)%nNOdes > 0) then 
!         write(zoneName, "(a,I1,a,I1)") "Zone", iDom, "_Proc_", myid

! 110  format('ZONE T=',a, " Nodes=", i5, " Elements=", i5, "  ZONETYPE=FEQUADRILATERAL")
        
!         write(101, 110), trim(zoneName), oWalls(iDom)%nNodes, oWalls(iDom)%nCells
!         write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"
! 13      format (E14.6)
!         do j=1,3
!            do i=1,oWalls(iDom)%nNodes
!               write(101,13) oWalls(iDom)%x(j, i)
!            end do
!         end do
        
!         do i=1,oWalls(iDom)%nCells
!            write(101,*) min(oWalls(iDom)%iBlank(i), oWalls(iDom)%iBlankNew(i))
!            !write(101,*) oWalls(iDom)%iBlank(i)
!         end do
        
! 15      format(I7, I7, I7, I7)
        
!         do i=1,oWalls(iDom)%nCells
!            write(101,15) oWalls(iDom)%conn(1, i), oWalls(iDom)%conn(2, i), &
!                 oWalls(iDom)%conn(3, i), oWalls(iDom)%conn(4, i)
!         end do
!      else
!         ! Do dummy stuff   
!         write(101, *) "Zone T=dummy Nodes=4, Elements=1, ZONETYPE=FEQUADRILATERAL"
!         write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"
!         do j=1,3
!            do i=1,4
!               write(101,13) zero
!            end do
!         end do
!         write(101,*) 1
!         write(101,*) '1 2 3 4'
!      end if
        

!   end do
!   close(101)
! end subroutine writeOWalls

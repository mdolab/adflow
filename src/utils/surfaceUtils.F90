module surfaceUtils

contains

  subroutine getSurfaceSize(sizeNode, sizeCell, famList, n, includeZipper)
    ! Compute the number of points that will be returned from getForces
    ! or getForcePoints
    use constants
    use blockPointers, only : BCData, nDom, nBocos
    use utils, only : setPointers
    use sorting, only : famInList
    use surfaceFamilies, only  : BCFamGroups
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent

    implicit none

    integer(kind=intType),intent(out) :: sizeNode, sizeCell
    logical, intent(in) :: includeZipper
    integer(kind=intType) :: nn, mm, i, j,iimax, shp(1)
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iBCGroup
    integer(kind=intType), intent(in) :: famList(n), n
    type(zipperMesh), pointer :: zipper
    logical :: BCGroupNeeded
    sizeNode = 0_intType
    sizeCell = 0_intType

    domains: do nn=1,nDom
       call setPointers(nn,1_intType,1_intType)
       bocos: do mm=1,nBocos
          ! Check if this surface should be included or not:
          famInclude: if (famInList(BCData(mm)%famID, famList)) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
             sizeNode = sizeNode + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)

             ! If we don't care about blanking, it's easy:
             blanking:if (.not. includeZipper) then
                sizeCell = sizeCell + (iEnd - iBeg)*(jEnd - jBeg)
             else
                ! Otherwise we have to consider the iBlank
                do j=jBeg+1, jEnd
                   do i=iBeg+1, iEnd
                      if (BCData(mm)%iBlank(i,j) == 1) then
                         sizeCell = sizeCell + 1
                      end if
                   end do
                end do
             end if blanking
          end if famInclude
       end do bocos
    end do domains

    ! We know must consider additional nodes that are required by the
    ! zipper mesh triangles on the root proc.

    ! No overset or we don't want to include the zipper, return immediately
    if (.not. oversetPresent .or. .not. includeZipper) then
       return
    end if

    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    do iBCGroup=1, nfamExchange
       BCGroupNeeded = .False.
       BCGroupFamLoop: do j=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(j), famList)) then
             BCGroupNeeded = .True.
             exit BCGroupFamLoop
          end if
       end do BCGroupFamLoop

       if (.not. BCGroupNeeded) then
          cycle
       end if

       ! Pointer for easier reading.
       zipper => zipperMeshes(iBCGroup)

       ! If we don't have a zipper for this BCGroup, just keep going.
       if (.not. zipper%allocated) then
          cycle
       end if

       ! Include the total extra number of nodes. Not necessairly all
       ! nodes are needed, but they will be returned anyway.
       sizeNode = sizeNode + size(zipper%indices)

       ! Include the extra number of cells. Not necessairly all cells
       ! are needed, but here we have to check indvidually.

       do i=1,size(zipper%fam)
          if (famInList(zipper%fam(i), famList)) then
             sizeCell = sizeCell + 1
          end if
       end do
    end do

  end subroutine getSurfaceSize

  subroutine getSurfaceConnectivity(conn, cgnsBlockID, ncell, famList, nFamList, includeZipper)
    ! Return the connectivity list for the each of the patches
    ! cgnsBlockID is the domain number of the CGNS file that each face patch corresponds to.
    ! the zipper mesh patches will have cgnsBlockID = -1
    use constants
    use blockPointers, only : nDom, nBocos, BCData, BCFaceID, rightHanded, nbkGlobal
    use utils, only : setPointers
    use sorting, only : famInList
    use surfaceFamilies, only : BCFamGroups
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent

    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: ncell
    integer(kind=intType), intent(inout) :: conn(4*ncell), cgnsBlockID(ncell)
    integer(kind=intType), intent(in) :: nFamList, famList(nFamList)
    logical, intent(in) :: includeZipper

    ! Working
    integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iBCGroup
    logical regularOrdering, BCGroupNeeded
    type(zipperMesh), pointer :: zipper
    cellCount = 0
    nodeCount = 0

    ! Initialize all cgns IDs to -1. This will take care of the zipper meshes IDs, which should
    ! be -1 since they do not belong to the CGNS file.
    cgnsBlockID = -1

    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)
       bocos: do mm=1,nBocos
          famInclude: if (famInList(BCData(mm)%famID, famList)) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             ni = iEnd - iBeg + 1
             nj = jEnd - jBeg + 1

             ! We want to ensure that all the normals of the faces are
             ! consistent. To ensure this, we enforce that all normals
             ! are "into" the domain. Therefore we must treat difference
             ! faces of a block differently. For example for an iLow
             ! face, when looping over j-k in the regular way, results
             ! in in a domain inward pointing normal for iLow but
             ! outward pointing normal for iHigh. The same is true for
             ! kMin and kMax. However, it is reverse for the J-faces:
             ! This is becuase the way the pointers are extracted i then
             ! k is the reverse of what "should" be for consistency. The
             ! other two, the pointers are cyclic consistent: i,j->k,
             ! j,k (wrap) ->i, but for the j-direction is is i,k->j when
             ! to be consistent with the others it should be
             ! k,i->j. Hope that made sense.

             select case(BCFaceID(mm))
             case(iMin, jMax, kMin)
                regularOrdering = .True.
             case default
                regularOrdering = .False.
             end select

             ! Now this can be reversed *again* if we have a block that
             ! is left handed.
             if (.not. rightHanded) then
                regularOrdering = .not. (regularOrdering)
             end if

             if (regularOrdering) then
                ! Do regular ordering.

                ! Loop over generic face size...Note we are doing zero
                ! based ordering!

                ! This cartoon of a generic cell might help:
                !
                ! i, j+1 +-----+ i+1, j+1
                !   n4   |     | n3
                !        +-----+
                !       i,j    i+1, j
                !       n1     n2
                !

                do j=0,nj-2
                   do i=0,ni-2
                      if (.not. includeZipper .or. BCData(mm)%iBlank(i+iBeg+1, j+jBeg+1) ==1 ) then
                         conn(4*cellCount+1) = nodeCount + (j  )*ni + i       + 1! n1
                         conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1   + 1! n2
                         conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1   + 1! n3
                         conn(4*cellCount+4) = nodeCount + (j+1)*ni + i       + 1! n4

                         ! Assign the corresponding block IDs.
                         ! Remember that nbkGlobal is the CGNS ID of the current domain,
                         ! and this is set when we call setPointers.
                         cgnsBlockID(cellCount+1) = nbkGlobal

                         cellCount = cellCount + 1
                      end if
                   end do
                end do
             else
                ! Do reverse ordering:
                do j=0,nj-2
                   do i=0,ni-2
                      if (.not. includeZipper .or. BCData(mm)%iBlank(i+iBeg+1, j+JBeg+1) ==1 ) then
                         conn(4*cellCount+1) = nodeCount + (j  )*ni + i        + 1! n1
                         conn(4*cellCount+2) = nodeCount + (j+1)*ni + i        + 1! n4
                         conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1    + 1! n3
                         conn(4*cellCount+4) = nodeCount + (j  )*ni + i + 1    + 1! n2

                         ! Assign the corresponding block IDs.
                         ! Remember that nbkGlobal is the CGNS ID of the current domain,
                         ! and this is set when we call setPointers.
                         cgnsBlockID(cellCount+1) = nbkGlobal

                         cellCount = cellCount + 1
                      end if
                   end do
                end do
             end if
             nodeCount = nodeCount + ni*nj
          end if famInclude
       end do bocos
    end do domains

    ! We know must consider additional connectivity required by the
    ! zipper mesh triangles on the root proc

    ! No overset or don't want zipper return immediately
    if (.not. oversetPresent .or. .not. includeZipper) then
       return
    end if

    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    BCGroupLoop: do iBCGroup=1, nFamExchange

       BCGroupNeeded = .False.
       BCGroupFamLoop: do i=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList)) then
             BCGroupNeeded = .True.
             exit BCGroupFamLoop
          end if
       end do BCGroupFamLoop

       if (.not. BCGroupNeeded) then
          cycle
       end if

       ! Pointer for easier reading.
       zipper => zipperMeshes(iBCGroup)

       ! If the zipper isn't done yet, don't do anything
       if (.not. zipper%allocated) then
          cycle
       end if

       ! Include the extra number of cells. Not necessairly all cells
       ! are needed, but ehre we have to check indvidually.

       do i=1,size(zipper%fam)
          if (famInList(zipper%fam(i), famList)) then
             ! This triangle should be included. Note that we use
             ! degenerate quads for the triangles.
             conn(4*cellCount+1) = nodeCount + zipper%conn(1, i)
             conn(4*cellCount+2) = nodeCount + zipper%conn(2, i)
             conn(4*cellCount+3) = nodeCount + zipper%conn(3, i)
             conn(4*cellCount+4) = nodeCount + zipper%conn(3, i)
             cellCount = cellCount + 1
          end if
       end do
    end do BCGroupLoop

  end subroutine getSurfaceConnectivity

  subroutine getSurfaceFamily(elemFam, ncell, famList, nFamList, includeZipper)

    use constants
    use blockPointers, only : nDom, nBocos, BCData
    use utils, only : setPointers
    use sorting, only : famInList
    use surfaceFamilies, only : BCFamGroups
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: ncell
    integer(kind=intType), intent(inout) :: elemFam(nCell)
    integer(kind=intType), intent(in) :: famList(nFamList), nFamList
    logical, intent(in) :: includeZipper

    ! Working
    integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iBCGroup
    logical BCGroupNeeded
    type(zipperMesh), pointer :: zipper

    cellCount = 0
    nodeCount = 0

    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)
       bocos: do mm=1,nBocos
          famInclude: if (famInList(BCData(mm)%famID, famList)) then

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             ni = iEnd - iBeg + 1
             nj = jEnd - jBeg + 1
             do j=0,nj-2
                do i=0,ni-2
                   if (.not. includeZipper .or. BCData(mm)%iBlank(i+iBeg+1, j+JBeg+1)==1 ) then
                      cellCount = cellCount + 1
                      elemFam(cellCount) = BCdata(mm)%famID
                   end if
                end do
             end do
             nodeCount = nodeCount + ni*nj
          end if famInclude
       end do bocos
    end do domains

    ! We know must consider additional elements quired by the zipper
    ! mesh triangles on the root proc

    ! No overset or don't want zipper
    if (.not. oversetPresent .or. .not. includeZipper) then
       return
    end if

    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    BCGroupLoop: do iBCGroup=1, nFamExchange

       BCGroupNeeded = .False.
       BCGroupFamLoop: do i=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList)) then
             BCGroupNeeded = .True.
             exit BCGroupFamLoop
          end if
       end do BCGroupFamLoop

       if (.not. BCGroupNeeded) then
          cycle
       end if

       ! Pointer for easier reading.
       zipper => zipperMeshes(iBCGroup)

       ! If the zipper isn't done yet, don't do anything
       if (.not. zipper%allocated) then
          cycle
       end if

       ! Include the extra number of cells. Not necessairly all cells
       ! are needed, but ehre we have to check indvidually.

       do i=1, size(zipper%fam)
          if (famInList(zipper%fam(i), famList)) then
             ! This triangle should be included. Note that we use
             ! degenerate quads for the triangles.
             cellCount = cellCount + 1
             elemFam(cellCount) = zipper%fam(i)
          end if
       end do
    end do BCGroupLoop
  end subroutine getSurfaceFamily

  subroutine getSurfacePoints(points, npts, sps_in, famList, nFamList, includeZipper)
    use constants
    use communication, only : myid
    use blockPointers, only : nDom, BCData, nBocos, x, BCFaceID, il, jl, kl
    use BCPointers, only : xx
    use surfaceFamilies, only : BCFamGroups, familyExchange, BCFamExchange
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent
    use sorting, only : famInList
    use utils, only : setPointers, EChk, setBCPointers
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    !
    !      Local variables.
    !
    integer(kind=intType), intent(in) :: npts,sps_in
    real(kind=realType), intent(inout) :: points(3,npts)
    integer(kind=intType), intent(in) :: nFamList, famList(nFamList)
    logical, intent(in) :: includeZipper

    integer(kind=intType) :: mm, nn, i, j, ii,sps, iDim, jj, ierr
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iBCGroup
    type(zipperMesh), pointer :: zipper
    type(familyexchange), pointer :: exch
    logical :: BCGroupNeeded
    real(kind=realType), dimension(:), pointer :: localPtr
    sps = sps_in

    ii = 0
    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, sps)

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,nBocos

          famInclude: if (famInList(BCData(mm)%famID, famList)) then

             ! NODE Based
             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             do j=jBeg, jEnd ! This is a node loop
                do i=iBeg, iEnd ! This is a node loop
                   ii = ii +1
                   select case(BCFaceID(mm))

                   case(imin)
                      points(:,ii) = x(1,i,j,:)
                   case(imax)
                      points(:,ii) = x(il,i,j,:)
                   case(jmin)
                      points(:,ii) = x(i,1,j,:)
                   case(jmax)
                      points(:,ii) = x(i,jl,j,:)
                   case(kmin)
                      points(:,ii) = x(i,j,1,:)
                   case(kmax)
                      points(:,ii) = x(i,j,kl,:)
                   end select
                end do
             end do
          end if famInclude
       end do bocos
    end do domains

    ! No overset or not zipper, return
    if (.not. oversetPresent .or. .not. includeZipper) then
       return
    end if

    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    do iBCGroup=1, nFamExchange

       zipper => zipperMeshes(iBCGroup)

       if (.not. zipper%allocated) then
          cycle
       end if

       exch => BCFamExchange(iBCGroup, sps)
       BCGroupNeeded = .False.
       BCGroupFamLoop: do i=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList)) then
             BCGroupNeeded = .True.
             exit BCGroupFamLoop
          end if
       end do BCGroupFamLoop

       if (.not. BCGroupNeeded) then
          cycle
       end if

       ! Now we know we *actually* need something from this BCGroup.

       ! Loop over each dimension individually since we have a scalar
       ! scatter.
       dimLoop: do iDim=1,3

          call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          localPtr = zero

          ! jj is the running counter through the pointer array.
          jj = 0
          do nn=1, nDom
             call setPointers(nn, 1_intType, sps)
             do mm=1, nBocos
                famInclude2: if (famInList(BCData(mm)%famID, exch%famList)) then
                   iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                   jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                   call setBCPointers(mm, .True.)
                   do j=jBeg, jEnd
                      do i=iBeg, iEnd
                         jj = jj+ 1
                         localPtr(jj) = xx(i+1, j+1, iDim)
                      end do
                   end do
                end if famInclude2
             end do
          end do

          ! Restore the pointer
          call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Now scatter this to the zipper
          call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! The values we need are precisely what is in zipper%localVal
          call vecGetArrayF90(zipper%localVal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Just copy the received data into the points array. Only root proc.
          if (myid == 0) then
             points(iDim, ii+1:ii+size(localPtr)) = localPtr
          end if
          ! The values we need are precisely what is in zipper%localVal
          call vecRestoreArrayF90(zipper%localVal, localPtr, ierr)
          call EChk(ierr,__FILE__,__LINE__)

       end do dimLoop

       ! Increcment the running ii counter.
       ii = ii + size(localPtr)

    end do
  end subroutine getSurfacePoints

  subroutine mapVector(vec1, n1, famList1, nf1, vec2, n2, famList2, nf2, includeZipper,nDim)

    ! Map one vector, vec1 of size (3,n1) defined on family list 'famList1' onto
    ! vector, vec2, of size (3, n2) defined on family list 'famList2'

    ! This operation is actually pretty fast since it just requires a
    ! single copy of surface-based data.
    use constants
    use blockPointers, onlY :nDom, flowDoms
    use sorting, only : famInList
    use surfaceFamilies, only : BCFamGroups
    use oversetData, only : zipperMeshes, zipperMesh, oversetPresent

    implicit none

    ! Input/Output
    integer(kind=intType) :: n1, n2, nf1, nf2,nDim
    integer(kind=intType), intent(in) :: famList1(nf1), famList2(nf2)
    real(kind=realType), intent(in) :: vec1(nDim, n1)
    real(kind=realType), intent(inout) :: vec2(nDim, n2)
    logical, intent(in) :: includeZipper

    ! Working
    integer(kind=intType) :: i, k, ii, jj, nn, mm, iSize
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, famID, iBCGroup
    logical :: fam1Included, fam2Included
    type(zipperMesh), pointer :: zipper
    logical :: BCGroupNeeed
    ii = 0
    jj = 0
    domains: do nn=1,nDom
       ! Don't set pointers for speed

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,flowDoms(nn, 1, 1)%nBocos
          famId = flowDoms(nn, 1, 1)%BCdata(mm)%famID

          fam1Included = famInList(famID, famList1)
          fam2Included = famInList(famid, famList2)

          jBeg = flowDoms(nn, 1, 1)%bcData(mm)%jnBeg
          jEnd = flowDoms(nn, 1, 1)%bcData(mm)%jnEnd

          iBeg = flowDoms(nn, 1, 1)%bcData(mm)%inBeg
          iEnd = flowDoms(nn, 1, 1)%bcData(mm)%inEnd
          iSize = (iEnd-iBeg+1)*(jEnd-jBeg+1)

          if (fam1Included .and. fam2Included) then
             ! The two lists overlap so copy:
             do k=1, iSize
                vec2(:, k+jj) = vec1(:, k+ii)
             end do
          end if

          ! Finally increment the counters if the face had been inclded
          if (fam1Included) then
             ii = ii + iSize
          end if

          if (fam2Included) then
             jj =jj + iSize
          end if

       end do bocos
    end do domains

    ! As with the rest of the code we have to account for the zipper
    ! mesh on the root proc.

    ! We know must consider additional nodes that are required by the
    ! zipper mesh triangles on the root proc.

    ! No overset or don't want to include zipper, return immediately
    if (.not. oversetPresent .or. .not. includeZipper) then
       return
    end if

    ! If there are zipper meshes, we must include the nodes that the
    ! zipper triangles will use.
    BCGroupLoop: do iBCGroup=1, nFamExchange

       zipper => zipperMeshes(iBCGroup)

       fam1Included = .False.
       fam2Included = .False.
       BCGroupFamLoop: do i=1, size(BCFamGroups(iBCGroup)%famList)
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList1)) then
             fam1Included = .True.
          end if
          if (famInList(BCFamGroups(iBCGroup)%famList(i), famList2)) then
             fam2Included = .True.
          end if
       end do BCGroupFamLoop

       ! This is the total number of nodes that this BCGroup has. It
       ! is not further broken down by family group.
       iSize = size(zipper%indices)

       if (fam1Included .and. fam2Included) then
          ! The two lists overlap so copy:
          do k=1, iSize
             if (k+ii <= n1) then
                vec2(:, k+jj) = vec1(:, k+ii)
             end if
          end do
       end if

       ! Finally increment the counters if this BCGroup had been included.
       if (fam1Included) then
          ii = ii + iSize
       end if

       if (fam2Included) then
          jj = jj + iSize
       end if
    end do BCGroupLoop

  end subroutine mapVector

  subroutine getWallList(wallList, nWallList, nFamTotal)

    ! Python wrapped utility function to return the list of families
    ! that are walls to Python since we need that information in
    ! Python for a few default values.

    use constants
    use surfaceFamilies, only :BCFamGroups
    implicit none

    integer(kind=intType), intent(in) :: nFamtotal
    integer(kind=intType), dimension(nFamTotal), intent(out) :: wallList
    integer(kind=intType), intent(out) :: nWallList

    nWallList = size(BCFamGroups(iBCGroupWalls)%famList)
    wallList(1:nWallList) = BCfamGroups(iBCGroupWalls)%famList

  end subroutine getWallList

end module surfaceUtils

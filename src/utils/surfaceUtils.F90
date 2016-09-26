module surfaceUtils

contains

  subroutine getSurfaceSize(size, sizeCell, famList, n)
    ! Compute the number of points that will be returned from getForces
    ! or getForcePoints
    use constants
    use blockPointers, only : BCData, nDom, nBocos
    use utils, only : setPointers
    use sorting, only : bsearchIntegers
    implicit none

    integer(kind=intType),intent(out) :: size, sizeCell
    integer(kind=intType) :: nn,mm
    integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd
    integer(kind=intType), intent(in) :: famList(n), n

    size = 0_intType
    sizeCell = 0_intType

    domains: do nn=1,nDom
       call setPointers(nn,1_intType,1_intType)
       bocos: do mm=1,nBocos
          ! Check if this surface should be included or not:
          famInclude: if (bsearchIntegers(BCdata(mm)%famID, famList) > 0) then 

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
             size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
             sizeCell = sizeCell + (iEnd - iBeg)*(jEnd - jBeg)
          end if famInclude
       end do bocos
    end do domains
  end subroutine getSurfaceSize

  subroutine getSurfaceConnectivity(conn, ncell)
    ! Return the connectivity list for the each of the patches
    use constants
    use blockPointers, only : nDom, nBocos, BCData, BCFaceID, rightHanded
    use surfaceFamilies, only : famGroups
    use utils, only : setPointers
    use sorting, only : bsearchIntegers
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: ncell
    integer(kind=intType), intent(inout) :: conn(4*ncell)

    ! Working
    integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
    logical regularOrdering
    cellCount = 0
    nodeCount = 0

    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)
       bocos: do mm=1,nBocos
          famInclude: if (bsearchIntegers(BCdata(mm)%famID, famGroups) > 0) then 

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
                      conn(4*cellCount+1) = nodeCount + (j  )*ni + i       + 1! n1
                      conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1   + 1! n2
                      conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1   + 1! n3
                      conn(4*cellCount+4) = nodeCount + (j+1)*ni + i       + 1! n4
                      cellCount = cellCount + 1
                   end do
                end do
             else
                ! Do reverse ordering:
                do j=0,nj-2
                   do i=0,ni-2
                      conn(4*cellCount+1) = nodeCount + (j  )*ni + i        + 1! n1
                      conn(4*cellCount+2) = nodeCount + (j+1)*ni + i        + 1! n4
                      conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1    + 1! n3
                      conn(4*cellCount+4) = nodeCount + (j  )*ni + i + 1    + 1! n2
                      cellCount = cellCount + 1
                   end do
                end do
             end if
             nodeCount = nodeCount + ni*nj
          end if famInclude
       end do bocos
    end do domains
  end subroutine getSurfaceConnectivity


  subroutine getSurfaceFamily(elemFam, ncell)

    use constants
    use blockPointers, only : nDom, nBocos, BCData
    use surfaceFamilies, only : famGroups
    use utils, only : setPointers
    use sorting, only : bsearchIntegers
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: ncell
    integer(kind=intType), intent(inout) :: elemFam(nCell)

    ! Working
    integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
    logical regularOrdering
    cellCount = 0
    nodeCount = 0

    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)
       bocos: do mm=1,nBocos
          famInclude: if (bsearchIntegers(BCdata(mm)%famID, famGroups) > 0) then 

             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             ni = iEnd - iBeg + 1
             nj = jEnd - jBeg + 1
             do j=0,nj-2
                do i=0,ni-2
                   cellCount = cellCount + 1
                   elemFam(cellCount) = BCdata(mm)%famID
                end do
             end do
             nodeCount = nodeCount + ni*nj
          end if famInclude
       end do bocos
    end do domains
  end subroutine getSurfaceFamily

  subroutine getSurfacePoints(points, npts, sps_in)
    use constants
    use blockPointers, only : nDom, BCData, nBocos, x, BCFaceID, il, jl, kl
    use surfaceFamilies, only : famGroups
    use utils, only : setPointers
    use sorting, only : bsearchIntegers
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType), intent(in) :: npts,sps_in
    real(kind=realType), intent(inout) :: points(3,npts)

    integer(kind=intType) :: mm, nn, i, j, ii,sps
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

    sps = sps_in

    ii = 0 
    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, sps)

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,nBocos

          famInclude: if (bsearchIntegers(BCdata(mm)%famID, famGroups) > 0) then 

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

  end subroutine getSurfacePoints

  subroutine mapVector(vec1, n1, famList1, nf1, vec2, n2, famList2, nf2)

    ! Map one vector, vec1 of size (3,n1) defined on family list 'famList1' onto
    ! vector, vec2, of size (3, n2) defined on family list 'famList2'

    ! This operation is actually pretty fast since it just requires a
    ! single copy of surface-based data. 
    use constants
    use blockPointers
    use sorting, only : bsearchIntegers
    implicit none

    ! Input/Output
    integer(kind=intType) :: n1, n2, nf1, nf2
    integer(kind=intType), intent(in) :: famList1(nf1), famList2(nf2)
    real(kind=realType), intent(in) :: vec1(3, n1)
    real(kind=realType), intent(inout) :: vec2(3, n2)

    ! Working
    integer(kind=intType) :: k, ii, jj, nn, mm, iSize, iBeg, iEnd, jBeg, jEnd, famID
    logical :: fam1Included, fam2Included

    ii = 0 
    jj = 0
    domains: do nn=1,nDom
       ! Don't set pointers for speed

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,flowDoms(nn, 1, 1)%nBocos
          famId = flowDoms(nn, 1, 1)%BCdata(mm)%famID

          fam1Included = bsearchIntegers(famID, famList1) > 0
          fam2Included = bsearchIntegers(famID, famList2) > 0

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
  end subroutine mapVector
  ! The only two groups we need to deal with in fortran directly is the
  ! fullFamilyList and all walls. We write special routines for them. 
  subroutine setFullFamilyList()

    use surfaceFamilies
    implicit none

    integer(kind=intType) :: i
    if (allocated(famGroups)) then 
       deallocate(famGroups)
    end if
    allocate(famGroups(totalFamilies))
    do i=1, totalFamilies
       famGroups(i) = i
    end do

  end subroutine setFullFamilyList

  subroutine setWallFamilyList()

    use surfaceFamilies
    implicit none
    integer(kind=intType) :: i

    if (allocated(famGroups)) then 
       deallocate(famGroups)
    end if
    allocate(famGroups(totalWallFamilies))
    do i=1, totalWallFamilies
       famGroups(i) = wallFamilies(i)
    end do

  end subroutine setWallFamilyList

  subroutine setFamilyInfo(famList, n)

    use surfaceFamilies
    implicit none
    integer(kind=intType), intent(in) :: famList(n), n

    if (allocated(famGroups)) then 
       deallocate(famGroups)
    end if
    allocate(famGroups(n))
    famGroups = famList

  end subroutine setFamilyInfo


end module surfaceUtils

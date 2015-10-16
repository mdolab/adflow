subroutine readInterpolation(baseName)

  use constants
  use communication
  use blockPointers
  use BCTypes

  ! Input variables
  character*(*), intent(in) :: baseName

  ! ! Local Variables
  ! integer(Kind=intType) :: nn, i, j, k, fileID
  ! character *128 integer_string, inter_file
  ! integer(kind=intType) :: nFringe, nDonor, iieptr, iisptr, level, sps
  ! integer(kind=intType), allocatable, dimension(:, :, :) :: iBlankTmp
  ! integer(kind=intTYpe), allocatable, dimension(:) :: tmpInt
  ! real(kind=realType), allocatable, dimension(:) :: tmpReal

  ! domains: do nn=1, nDom
     
  !    ! Get the name of the interpolation file
  !    write(integer_string,*) nn-1
  !    inter_file = 'inter.'//adjustl(integer_string)

  !    fileID = 2
  !    open(unit=fileID, file=inter_file, form='unformatted', status='unknown')
  !    read(fileID) nfringe, ndonor, iieptr, iisptr

  !    ! Allocate the required space in flowDoms for the interpolation
  !    ! information for this block. 
  !    allocate(flowDoms(nn, 1, 1)%iDonor(3, ndonor), flowDoms(nn, 1, 1)%frac(3, ndonor))
  !    allocate(flowDoms(nn, 1, 1)%iMesh(3, nfringe), flowDoms(nn, 1, 1)%iBC(nfringe))

  !    ! Now set the pointers for the memory allocated above to make the
  !    ! code easier to read
  !    call setPointers(nn, 1, 1)
         
  !    ! Since the donor information is currently stored in field as
  !    ! opposed to block format we cannot reliable read
  !    ! directly. Allocate a temp integer array and use that. 

  !    allocate(tmpInt(3*nDonor))
  !    read(fileID) (tmpInt(i), i=1, 3*nDonor)
     
  !    ! Distribute
  !    do i=1, nDonor
  !       idonor(1, i) = tmpInt(3*i-2)
  !       idonor(2, i) = tmpInt(3*i-1)
  !       idonor(3, i) = tmpInt(3*i  )
  !    end do

  !    ! Release memory
  !    deallocate(tmpInt)

  !    ! Same logic for the fraction..tmpReal and redistribute
  !    allocate(tmpReal(3*nDonor))
  !    read(fileID) (tmpReal(i), i=1, 3*nDonor)
     
  !    ! Distribute
  !    do i=1, nDonor
  !       frac(1, i) = tmpReal(3*i-2)
  !       frac(2, i) = tmpReal(3*i-1)
  !       frac(3, i) = tmpReal(3*i  )
  !    end do

  !    ! Release memory
  !    deallocate(tmpReal)

  !    ! Again same logic for iMesh
  !    allocate(tmpInt(3*nFringe))
  !    read(fileID) (tmpInt(i), i=1, 3*nFringe)

  !    ! Distribute
  !    do i=1, nFringe
  !       imesh(1, i) = tmpInt(3*i-2)
  !       imesh(2, i) = tmpInt(3*i-1)
  !       imesh(3, i) = tmpInt(3*i  )
  !    end do

  !    ! Release tmp memory
  !    deallocate(tmpInt)

  !    ! ibc is scalar so we can read directly
  !    read(fileID) (ibc(i), i=1, nFringe)

  !    ! We also need to get the iblank array from the interpolation
  !    ! file. We will only read the blanks on the owned cells, not the
  !    ! halos...we'll communicate using 1to1 connectivity to iblank on
  !    ! the halos. 
         
  !    ! Since this is scalar we can do this one directly as well
  !    read(fileID) (((iblank(i, j, k), i=2, il), j=2, jl), k=2, kl)

  !    ! Set the blanking value to -1 for all fringe cells on this block
  !    do id=1, nfringe
  !       i = imesh(1, id)
  !       j = imesh(2, id)
  !       k = imesh(3, id)
  !       iblank(i, j, k) = -1 
  !    enddo

  !    ! We also ned to refine the iblanking based on the boundary
  !    ! conditions. For now just set the last two layers next to a
  !    ! overSetOuterBoundary condition. 
     
  !    do mm=1,nBocos
  !       if (bcType(mm) == OversetOuterBoundary) then 
  !          select case (BCFaceID(mm))
  !          case (iMin)
  !             iblank(2, :, :) = -1
  !             iblank(3, :, :) = -1
  !          case (iMax)
  !             iblank(il, :, :) = -1
  !             iblank(nx, :, :) = -1
  !          case (jMin)
  !             iblank(:, 2, :) = -1
  !             iblank(:, 3, :) = -1
  !          case (jMax)
  !             iblank(:, jl, :) = -1
  !             iblank(:, ny, :) = -1
  !          case (kMin)
  !             iblank(:, :, 2) = -1
  !             iblank(:, :, 3) = -1
  !          case (kMax)
  !             iblank(:, :, kl) = -1
  !             iblank(:, :, nz) = -1
  !          end select
  !       end if
  !    end do
  ! end do domains

  ! ! Now we need to do a halo exchange for hte iblank arrays. 
  ! level = 1
  ! sps = 1 
  ! call exchangeIblanks(level, sps, commPatternCell_2nd, internalCell_2nd)

end subroutine readInterpolation

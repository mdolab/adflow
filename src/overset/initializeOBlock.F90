
subroutine initializeOBlock(oBlock, nn)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use inputOverset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  use communication
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType) :: nn, kk

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  real(kind=realType) :: factor, frac, exponent, wallEdge, avgEdge, dist
  logical :: wallsPresent, isWallType
  real(kind=realType), dimension(:, :, :, :), allocatable, target :: tmpX
  type(Xplane), dimension(:), allocatable :: planes

  ! Set all the sizes for this block.
  oBlock%il = il
  oBlock%jl = jl
  oBlock%kl = kl

  oBlock%proc = myID
  oBlock%block = nn

  call wallsOnBlock(wallsPresent)

  allocate( &
       oBlock%qualDonor(1, ie*je*ke), &
       oBlock%globalCell(0:ib, 0:jb, 0:kb), &
       oBlock%nearWall(1:ie, 1:je, 1:ke), &
       oBlock%invalidDonor(1:ie, 1:je, 1:ke))

  oBlock%nearWall = 0
  oBlock%invalidDonor = 0

  call flagForcedReceivers(oBlock%invalidDonor)
  
  ! Add to the invalid donor list:
  kk = 0
  do k=1, ke
     do j=1, je
        do i=1, ie
           if (iblank(i,j,k) == -2 .or. iblank(i,j,k)==-3 .or. iblank(i,j,k) == 0) then 
              oBlock%invalidDonor(i,j,k) = 1
              kk = kk + 1
           end if
        end do
     end do
  end do

  ! Copy Volume to qualDonor and do minVol while we're at it
  oBlock%minVol = Large
  mm = 0
  exponent = third
  do k=1,ke
     do j=1,je
        do i=1,ie
           mm = mm + 1
           if (wallsPresent) then 

              wallEdge = fourth*(&
                   norm2(x(i-1, j-1, k-1, :) - x(i-1, j-1, k, :)) + &
                   norm2(x(i  , j-1, k-1, :) - x(i  , j-1, k, :)) + &
                   norm2(x(i-1, j  , k-1, :) - x(i-1, j  , k, :)) + &
                   norm2(x(i  , j  , k-1, :) - x(i  , j  , k, :)))
              
              avgEdge = vol(i, j, k)**exponent
              !oBlock%qualDonor(1, mm) = half*(avgEdge + wallEdge)
              !oBlock%qualDonor(1, mm) = min(avgEdge, wallEdge)
              oBlock%qualDonor(1, mm) = avgEdge
           else
              factor = 4.0
              oBlock%qualDonor(1, mm) = (vol(i, j, k)*factor)**exponent
           end if

           oBlock%minVol = min(oBlock%minVol,  oBlock%qualDonor(1, mm))
        end do
     end do
  end do
  
  !Copy over global cell
  oBlock%globalCell = globalCell
  
  ! Now setup the data for the ADT
  nHexa = il * jl * kl
  nADT = ie * je * ke

  allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa))
  allocate(tmpX(3, ie, je, ke))
  ! Fill up the xADT using cell centers (dual mesh)
  mm = 0
  do k=1, ke
     do j=1, je
        do i=1, ie
           mm = mm + 1
           oBlock%xADT(:, mm) = eighth*(&
                x(i-1, j-1, k-1, :) + &
                x(i  , j-1, k-1, :) + &
                x(i-1, j  , k-1, :) + &
                x(i  , j  , k-1, :) + &
                x(i-1, j-1, k  , :) + &
                x(i  , j-1, k  , :) + &
                x(i-1, j  , k  , :) + &
                x(i  , j  , k  , :))
           tmpX(:, i, j, k) = oBlock%xADT(:, mm)
        end do
     end do
  end do

  ! Flag all nodes that are within nearWallDist as being nearWall
  do mm=1, nBocos
     if (isWallType(BCType(mm))) then 
        select case (BCFaceID(mm))
        case (iMin)
           allocate(planes(ie))
           do i=1, ie
              planes(i)%xx => tmpX(:, i, :, :)
              planes(i)%nearWall => oBlock%nearWall(i, :, :)
           end do
        case (iMax)
           allocate(planes(ie))
           do i=1, ie
              planes(i)%xx => tmpX(:, ie-i+1, :, :)
              planes(i)%nearWall => oBlock%nearWall(ie-i+1, :, :)
           end do
        case (jMin)
           allocate(planes(je))
           do j=1, je
              planes(j)%xx => tmpX(:, :, j, :)
              planes(j)%nearWall => oBlock%nearWall(:, j, :)
           end do
        case (jMax)
           allocate(planes(je))
           do j=1, je
              planes(j)%xx => tmpX(:, :, je-j+1, :)
              planes(j)%nearWall => oBlock%nearWall(:, je-j+1, :)
           end do
        case (kMin)
           allocate(planes(ke))
           do k=1, ke
              planes(k)%xx => tmpX(:, :, :, k)
              planes(k)%nearWall => oBlock%nearWall(:, :, k)
           end do
        case (kMax)
           allocate(planes(ke))
           do k=1, ke
              planes(k)%xx => tmpX(:, :, :, ke-k+1)
              planes(k)%nearWall => oBlock%nearWall(:, :, ke-k+1)
           end do
        end select
        
        ! Loop over the size of the generalized plane
        do j=1, size(planes(1)%xx, 3)
           do i=1, size(planes(1)%xx, 2)
              
              ! Loop over the 'k' ie offwall direction
              do k=1, size(planes) 
                 
                 dist = norm2(planes(k)%xx(:, i, j) - planes(1)%xx(:, i, j))
                 if (dist < nearWallDist) then 
                    planes(k)%nearWall(i, j) = 1
                 end if
              end do
           end do
        end do
     deallocate(planes)
  end if
end do ! BocoLoop
deallocate(tmpX)



  mm = 0
  ! These are the 'elements' of the dual mesh.
  planeOffset = ie * je
  do k=2, ke
     do j=2, je
        do i=2, ie
           mm = mm + 1
           oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
           oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
           oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
           oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

           oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
           oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
           oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
           oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
        end do
     end do
  end do

  ! Call the custom build routine -- Serial only, only Hexa volumes,
  ! we supply our own ADT Type

  call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)

  ! Flag this block as being allocated
  oBlock%allocated = .True.

end subroutine initializeOBlock


subroutine initializeOBlock(oBlock, nn, level, sps)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use inputOverset
  use blockPointers
  use adtBuild, only : buildSerialHex
  use cgnsGrid
  use communication
  use stencils
  use utils, only : mynorm2
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType) :: nn, level, sps, kk

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  real(kind=realType) :: factor, frac,  dist, xp(3)
  integer(kind=intType) :: i_stencil, ii, jj, iii
  logical :: wallsPresent, isWallType
  logical, allocatable, dimension(:, :, :)  :: nearWallTmp

  ! Set all the sizes for this block.
  oBlock%il = il
  oBlock%jl = jl
  oBlock%kl = kl

  oBlock%proc = myID
  oBlock%block = nn
  oBlock%cluster = clusters(cumDomProc(myid) + nn)
  call wallsOnBlock(wallsPresent)

  ! Do the reset of the allocs
  allocate( &
       oBlock%qualDonor(1, ie*je*ke), &
       oBlock%globalCell(0:ib, 0:jb, 0:kb), &
       oBlock%invalidDonor(1:ie, 1:je, 1:ke))

  oBlock%invalidDonor = 0
  call flagForcedReceivers(oBlock%invalidDonor)
  
  ! Add to the invalid donor list if it got flooded with iblank of -2 or -3:
  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! This is a hard interior cell. Flag EVERY cell it it's
           ! stencil as a invalid donor. 
           if (iblank(i, j, k) ==-3 .or. iBlank(i, j,k) == -2) then 
          
              stencilLoop: do i_stencil=1, N_visc_drdw
                 ii = visc_drdw_stencil(i_stencil, 1) + i
                 jj = visc_drdw_stencil(i_stencil, 2) + j
                 kk = visc_drdw_stencil(i_stencil, 3) + k
                 
                 ! Make sure we're at least at 1-level halos
                 if (ii >= 1 .and. ii <= ie .and. jj >= 1 .and. jj<= je .and. &
                      kk >= 1 .and. kk <= ke) then 
                    oBlock%invalidDonor(ii, jj, kk) = 1
                
                 end if
              end do stencilLoop
           end if
        end do
     end do
  end do

  ! Copy Volume to qualDonor and do minVol while we're at it
  oBlock%minVol = Large
  mm = 0

  do k=1,ke
     do j=1,je
        do i=1,ie
           mm = mm + 1
           if (wallsPresent) then 
              
              ii = i
              jj = j
              kk = k
              ! If the cell is a boundary halo, use the real cell
              if (globalCell(i, j, k) < 0) then 

                 ii = min(max(2, i), il)
                 jj = min(max(2, j), jl)
                 kk = min(max(2, k), kl)
              end if
              
              xp = eighth*(&
                   x(ii-1, jj-1, kk-1, :) + &
                   x(ii  , jj-1, kk-1, :) + &
                   x(ii-1, jj  , kk-1, :) + &
                   x(ii  , jj  , kk-1, :) + &
                   x(ii-1, jj-1, kk  , :) + &
                   x(ii  , jj-1, kk  , :) + &
                   x(ii-1, jj  , kk  , :) + &
                   x(ii  , jj  , kk  , :))

              ! dist = mynorm2(xp - xSeed(i, j, k, :))
              ! frac = dist/clusterMarchDist(oBlock%cluster)
              frac = one
              oBlock%qualDonor(1, mm) = frac*vol(i, j, k)**third 
              
           else
              oBlock%qualDonor(1, mm) = (backGroundVolScale*vol(i, j, k))**third
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
  ! Fill up the xADT using cell centers (dual mesh)
  mm = 0

  ! Allocate the nearWall
  allocate(oBlock%nearWall(1:il, 1:jl, 1:kl))
  oBlock%nearWall = 0

  allocate(nearWallTmp(1:ie, 1:je, 1:ke))
  nearWallTmp = .False.

  do k=1, ke
     do j=1, je
        do i=1, ie
           mm = mm + 1
           xp = eighth*(&
                x(i-1, j-1, k-1, :) + &
                x(i  , j-1, k-1, :) + &
                x(i-1, j  , k-1, :) + &
                x(i  , j  , k-1, :) + &
                x(i-1, j-1, k  , :) + &
                x(i  , j-1, k  , :) + &
                x(i-1, j  , k  , :) + &
                x(i  , j  , k  , :))
           oBlock%xADT(:, mm) = xp
           
           ! Determine if this point is near wall. Note that the
           ! boundary halos sill have xSeed as "large" so these won't
           ! be flagged as nearWall. We will account for this below. 
           dist = mynorm2(xp - xSeed(i, j, k, :))
           if (dist < nearWallDist) then 
              nearWallTmp(i, j, k) = .True. 
           end if
        end do
     end do
  end do

  ! Now finally set the nearwall for the dual mesh cells. It is
  ! considered a near wall if all "nodes" of the dual mesh cell are
  ! also near wall. Have to be carful not to count boundary halos
  ! since they do not have nearWallTmp Values.
  
  do k=1, kl
     do j=1, jl
        do i=1, il
           if (&
                (nearWallTmp(i  , j  , k  ) .or. globalCell(i  , j,   k  ) < 0) .and. &
                (nearWallTmp(i+1, j  , k  ) .or. globalCell(i+1, j,   k  ) < 0) .and. &
                (nearWallTmp(i  , j+1, k  ) .or. globalCell(i  , j+1, k  ) < 0) .and. &
                (nearWallTmp(i+1, j+1, k  ) .or. globalCell(i+1, j+1, k  ) < 0) .and. &
                (nearWallTmp(i  , j  , k+1) .or. globalCell(i  , j,   k+1) < 0) .and. &
                (nearWallTmp(i+1, j  , k+1) .or. globalCell(i+1, j,   k+1) < 0) .and. &
                (nearWallTmp(i  , j+1, k+1) .or. globalCell(i  , j+1, k+1) < 0) .and. &
                (nearWallTmp(i+1, j+1, k+1) .or. globalCell(i+1, j+1, k+1) < 0)) then 
              oBlock%nearWall(i, j, k) = 1
           end if
        end do
     end do
  end do
  
  deallocate(nearWallTmp)
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

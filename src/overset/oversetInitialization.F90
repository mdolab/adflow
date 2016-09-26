module oversetInitialization

contains
  subroutine initializeFringes(nn, level, sps)

    ! This subroutine initializes the fringe information for the given
    ! block, level and spectral instance. It is assumed that
    ! blockPointers are already set. 
    use constants
    use communication
    use blockPointers
    use overset
    use stencils
    use inputOverset
    use utils, only : isWallType
    use oversetUtilities, only : flagForcedReceivers, emptyFringe, setIsCompute, &
         wallsOnBlock
    implicit none

    ! Input Params
    integer(kind=intType), intent(in) :: nn, level, sps

    ! Working Params
    integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj
    integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
    logical :: wallsPresent
    integer(kind=intType) :: i_stencil, clusterID
    integer(kind=intType), dimension(:, :, :), allocatable :: tmp
    real(kind=realType) :: frac, dist, xp(3)
    ! Allocate space for the double halo fringes. 
    if (.not. associated(flowDoms(nn, level, sps)%fringes)) then 
       allocate(flowDoms(nn, level, sps)%fringes(0:ib, 0:jb, 0:kb))
    end if

    clusterID = clusters(cumDomProc(myid) + nn)

    ! Check if we have walls:
    call wallsOnBLock(wallsPresent)

    ! Manually set the pointer to fringes we (possibly) just allocated
    ! instead of calling setPointers again
    fringes => flowDoms(nn, level, sps)%fringes

    ! Loop over all cells, including halos, setting all isCompute to
    ! false. This is important. Halos cells on boundaries are *not*
    ! considered compute cells. 

    do k=0, kb
       do j=0, jb
          do i=0, ib
             call emptyFringe(fringes(i, j, k))
             fringes(i, j, k)%myI = i
             fringes(i, j, k)%myJ = j
             fringes(i, j, k)%myK = k
             fringes(i, j, k)%myBlock = nn
             call setIsCompute(fringes(i, j, k)%status, .False.)
          end do
       end do
    end do

    ii = 0
    do k=2, kl
       do j=2, jl
          do i=2, il
             call setIsCompute(fringes(i, j, k)%status, .True.)
             ii = ii + 1
             if (wallsPresent) then 

                xp = eighth*(&
                     x(i-1, j-1, k-1, :) + &
                     x(i  , j-1, k-1, :) + &
                     x(i-1, j  , k-1, :) + &
                     x(i  , j  , k-1, :) + &
                     x(i-1, j-1, k  , :) + &
                     x(i  , j-1, k  , :) + &
                     x(i-1, j  , k  , :) + &
                     x(i  , j  , k  , :))

                ! dist = norm2(xp - xSeed(i, j, k, :))
                ! frac = dist/clusterMarchDist(clusterID)
                frac = one

                fringes(i, j, k)%quality = frac*vol(i, j, k)**third
             else
                fringes(i, j, k)%quality = (backgroundVolScale*vol(i, j, k))**third
             end if
          end do
       end do
    end do

    ! Now loop over this block's boundary condiitons and we need to set
    ! a litle info: We need to flag the two levels next to an overset
    ! outer boundary as being a "forced receiver'. To implement this, we
    ! set the volume of these cells to "large".  We use the generic
    ! flagForcedReceiver routine for this since the same information is
    ! used elsewhere.

    allocate(tmp(1:ie, 1:je, 1:ke))
    call flagForcedReceivers(tmp)
    do k=2, kl
       do j=2, jl
          do i=2, il
             if (tmp(i,j,k) == 1) then
                fringes(i, j, k)%quality = large
             end if
          end do
       end do
    end do
    deallocate(tmp)

    ! Flag all the interior hole cells with a *negative* quality since
    ! this means it won't try to get a stencil:

    do k=2, kl
       do j=2, jl
          do i=2, il
             if (iBlank(i, j, k) == -2 .or. iblank(i, j, k)==-3) then 
                fringes(i, j, k)%quality = -large
             end if
          end do
       end do
    end do

    ! Flag the cells *surrounding* the hold cells with large
    ! quality. That forces them to get a donor. Note that we *don't*
    ! overwrite the exisitng -2 and -3 cells. 

    do k=0, kb
       do j=0, jb
          do i=0, ib
             if (iBlank(i,j,k) == -2 .or. iblank(i,j,k)==-3) then 

                stencilLoop: do i_stencil=1, N_visc_drdw
                   ii = visc_drdw_stencil(i_stencil, 1) + i
                   jj = visc_drdw_stencil(i_stencil, 2) + j
                   kk = visc_drdw_stencil(i_stencil, 3) + k

                   ! Make sure we're on-block
                   if (ii >=2 .and. ii <= il .and. jj >= 2 .and. jj<= jl .and. &
                        kk >=2 .and. kk <= kl) then 
                      if (iblank(ii, jj, kk) /= -3 .and. iblank(ii, jj, kk) /=-2) then 
                         fringes(ii, jj, kk)%quality = large
                      end if
                   end if
                end do stencilLoop
             end if
          end do
       end do
    end do

    ! Flag the actual halo cells behind an overset outer boundary as holes.
    do mm=1,nBocos

       select case (BCFaceID(mm))
       case (iMin)
          iStart=0; iEnd=1;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (iMax)
          iStart=ie; iEnd=ib;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=0; jEnd=1;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=je; jEnd=jb;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (kMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=0; kEnd=1;
       case (kMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=ke; kEnd=kb;
       end select

       if (BCType(mm) == OversetOuterBound) then
          do k=kStart, kEnd
             do j=jStart, jEnd
                do i=iStart, iEnd
                   iblank(i,j,k) = 0
                end do
             end do
          end do
       end if
    end do

    ! Set the original quality
    do k=2, kl
       do j=2, jl
          do i=2, il
             fringes(i, j, k)%origQuality = fringes(i, j, k)%quality
          end do
       end do
    end do
  end subroutine initializeFringes


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
    use utils, only : mynorm2, isWallType
    use oversetUtilities, only : flagForcedReceivers, wallsOnBlock
    implicit none 

    ! Input Params
    type(oversetBlock), intent(inout) :: oBlock
    integer(kind=intType) :: nn, level, sps, kk

    ! Working paramters
    integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
    integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
    real(kind=realType) :: factor, frac,  dist, xp(3)
    integer(kind=intType) :: i_stencil, ii, jj, iii
    logical :: wallsPresent
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

  subroutine initializeOFringes(oFringe, nn)

    ! This subroutine initializes the fringe information for the given
    ! block, level and spectral instance. It is assumed that
    ! blockPointers are already set. 
    use communication
    use blockPointers
    use overset
    use stencils
    use inputOverset
    use utils
    use oversetUtilities, only : flagForcedReceivers, wallsOnBlock
    implicit none

    ! Input Params
    type(oversetFringe), intent(inout) :: oFringe
    integer(kind=intType), intent(in) :: nn

    ! Working Params
    integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj
    integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
    logical :: wallsPresent
    integer(kind=intType) :: i_stencil
    integer(kind=intType), dimension(:, :, :), allocatable :: tmp
    real(kind=realType) :: dist, frac, xp(3)
    ! Check if we have walls:
    call wallsOnBLock(wallsPresent)

    ! Set the sizes for the oFringe and allocate the required space. 
    oFringe%il = il
    oFringe%jl = jl
    oFringe%kl = kl
    oFringe%cluster = clusters(cumDomProc(myid) + nn)

    mm = nx*ny*nz
    allocate(oFringe%x(3, mm))
    allocate(oFringe%quality(mm))
    allocate(oFringe%origQuality(mm))
    allocate(oFringe%myBlock(mm))
    allocate(oFringe%myIndex(mm))
    allocate(oFringe%donorProc(mm))
    allocate(oFringe%donorBlock(mm))
    allocate(oFringe%dI(mm))
    allocate(oFringe%dJ(mm))
    allocate(oFringe%dK(mm))
    allocate(oFringe%donorFrac(3, mm))
    allocate(oFringe%gInd(8, mm))
    allocate(oFringe%xSeed(3, mm))
    allocate(oFringe%wallInd(mm))
    allocate(oFringe%isWall(mm))

    ! These default set the entire array
    oFringe%myBlock = nn
    oFringe%donorProc = -1
    oFringe%dI = -1
    oFringe%dJ = -1
    oFringe%dK = -1
    oFringe%donorFrac = -one
    oFringe%gInd = -1
    oFringe%isWall = 0
    oFringe%xSeed = large
    oFringe%wallInd = 0

    ! Now loop over the actual compute cells, setting the cell center
    ! value 'x', the volume and flag these cells as compute
    ii = 0
    do k=2, kl
       do j=2, jl
          do i=2, il
             ii = ii + 1
             do iDim=1, 3
                oFringe%x(iDim, ii) = eighth*(&
                     x(i-1, j-1, k-1, iDim) + &
                     x(i  , j-1, k-1, iDim) + &
                     x(i-1, j  , k-1, iDim) + &
                     x(i  , j  , k-1, iDim) + &
                     x(i-1, j-1, k  , iDim) + &
                     x(i  , j-1, k  , iDim) + &
                     x(i-1, j  , k  , iDim) + &
                     x(i  , j  , k  , iDim))
             end do

             if (wallsPresent) then 

                ! dist = norm2(oFringe%x(:, ii) - xSeed(i, j, k, :))
                ! frac = dist/clusterMarchDist(oFringe%cluster)
                frac = one
                oFringe%quality(ii) = frac*vol(i, j, k)**third
             else
                oFringe%quality(ii) = (backgroundVolScale*vol(i, j, k))**third
             end if

             oFringe%myIndex(ii) = ii

             oFringe%xSeed(:, ii) = xSeed(i, j, k, :)
             oFringe%wallInd(ii) = wallInd(i, j, k)

          end do
       end do
    end do

    ! Now loop over this block's boundary condiitons and we need to set
    ! a litle info: We need to flag the two levels next to an overset
    ! outer boundary as being a "forced receiver'. To implement this, we
    ! set the volume of these cells to "large".  We use the generic
    ! flagForcedReceiver routine for this since the same information is
    ! used elsewhere.

    allocate(tmp(1:ie, 1:je, 1:ke))
    call flagForcedReceivers(tmp)
    do k=2, kl
       do j=2, jl
          do i=2, il
             if (tmp(i,j,k) == 1) then
                ii = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
                oFringe%quality(ii) = large
             end if

          end do
       end do
    end do
    deallocate(tmp)

    ! Flag all the interior hole cells with a *negative* quality since
    ! this means it won't try to get a stencil:

    do k=2, kl
       do j=2, jl
          do i=2, il
             if (iblank(i,j,k)==-3 .or. iblank(i,j,k)==-2) then 
                iii = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
                oFringe%quality(iii) = -large
             end if
          end do
       end do
    end do

    ! Flag the cells *surrounding* the hold cells with large
    ! quality. That forces them to get a donor. 

    do k=0, kb
       do j=0, jb
          do i=0, ib
             if (iblank(i,j,k)==-3 .or. iblank(i,j,k)==-2) then 

                stencilLoop: do i_stencil=1, N_visc_drdw
                   ii = visc_drdw_stencil(i_stencil, 1) + i
                   jj = visc_drdw_stencil(i_stencil, 2) + j
                   kk = visc_drdw_stencil(i_stencil, 3) + k

                   ! Make sure we're on-block
                   if (ii >=2 .and. ii <= il .and. jj >= 2 .and. jj<= jl .and. &
                        kk >=2 .and. kk <= kl) then 
                      if (iblank(ii, jj, kk) /= -3 .and. iblank(ii,jj,kk) /= -2) then 
                         iii = (kk-2)*nx*ny + (jj-2)*nx + (ii-2) + 1
                         oFringe%quality(iii) = large
                      end if
                   end if
                end do stencilLoop
             end if
          end do
       end do
    end do

    ! Set the original quality. 
    do ii=1, mm
       oFringe%origQuality(ii) = oFringe%quality(ii)
    end do

    ! We also need to flag a single layer of cells next a wall
    ! boundary condition as being "isWall". Knowing the fringes
    ! next to walls will be necessary for determine the overap
    ! wall distance correction as well as the flooding procedure.

    do mm=1, nBocos
       select case(BCFaceID(mm))
       case (iMin)
          iStart=2; iEnd=2;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (iMax)
          iStart=il; iEnd=il;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=2; jEnd=2;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=jl; jEnd=jl;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (kMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=2; kEnd=2;
       case (kMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=kl; kEnd=kl;
       end select

       if (isWallType(BCType(mm))) then 
          do k=kStart, kEnd
             do j=jStart, jEnd
                do i=iStart, iEnd
                   ! Recompute the index
                   ii = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
                   oFringe%isWall(ii) = bcFaceID(mm)
                end do
             end do
          end do
       end if
    end do ! BocoLoop

    ! Flag this set of fringes as being allocated
    oFringe%allocated = .True.


  end subroutine initializeOFringes

  subroutine initializeOWall(oWall, dualMesh, cluster)

    ! This routine builds the ADT tree for any wall surfaces for the
    ! block currently being pointed to by block Pointers.
    use overset
    use blockPointers
    use adtBuild, only : buildSerialQuad
    use kdtree2_module
    use utils
    use oversetPackingRoutines, only : getWallSize
    implicit none 

    ! Input Params
    type(oversetWall), intent(inout) :: oWall
    logical, intent(in) :: dualMesh 
    integer(kind=intType), intent(in) :: cluster

    ! Working paramters
    integer(kind=intType) :: i, j, k, n, ii, jj, jjj, mm, ni, nj, nodeCount
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, nNodes, maxCells, nCells, iNode

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

end module oversetInitialization

subroutine initializeFringes(nn, level, sps)

  ! This subroutine initializes the fringe information for the given
  ! block, level and spectral instance. It is assumed that
  ! blockPointers are already set. 
  use communication
  use blockPointers
  use overset
  use stencils
  use inputOverset
  implicit none
  
  ! Input Params
  integer(kind=intType), intent(in) :: nn, level, sps

  ! Working Params
  integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj
  integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  logical :: wallsPresent, isWallType
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

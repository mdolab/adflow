subroutine initializeFringes(nn, level, sps)

  ! This subroutine initializes the fringe information for the given
  ! block, level and spectral instance. It is assumed that
  ! blockPointers are already set. 
  use communication
  use blockPointers
  use overset
  use BCTypes
  use stencils
  implicit none
  
  ! Input Params
  integer(kind=intType), intent(in) :: nn, level, sps

  ! Working Params
  integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj
  integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  real(kind=realType) :: factor, frac, exponent, avgEdge, wallEdge
  logical :: wallsPresent, isWallType
  integer(kind=intType) :: i_stencil
  integer(kind=intType), dimension(:, :, :), allocatable :: tmp
  ! Allocate space for the double halo fringes. 
  if (.not. associated(flowDoms(nn, level, sps)%fringes)) then 
     allocate(flowDoms(nn, level, sps)%fringes(0:ib, 0:jb, 0:kb))
  end if

  ! Check if we have walls:
  call wallsOnBLock(wallsPresent)

  ! Manually set the pointer to fringes we just allocated instead of
  ! calling setPointers again
  fringes => flowDoms(nn, level, sps)%fringes

  ! Loop over all cells, including halos, setting all isCompute to
  ! false

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
  exponent = third
  do k=2, kl
     do j=2, jl
        do i=2, il
           call setIsCompute(fringes(i, j, k)%status, .True.)
           ii = ii + 1
           if (wallsPresent) then 

              wallEdge = fourth*(&
                   norm2(x(i-1, j-1, k-1, :) - x(i-1, j-1, k, :)) + &
                   norm2(x(i  , j-1, k-1, :) - x(i  , j-1, k, :)) + &
                   norm2(x(i-1, j  , k-1, :) - x(i-1, j  , k, :)) + &
                   norm2(x(i  , j  , k-1, :) - x(i  , j  , k, :)))

              avgEdge = vol(i, j, k)**exponent
              
              !fringes(i, j, k)%quality = half*(avgEdge + wallEdge)
              !fringes(i, j, k)%quality = min(avgEdge, wallEdge)
              fringes(i, j, k)%quality = avgEdge
           else
              factor = 4.0
              fringes(i, j, k)%quality = (factor*vol(i, j, k))**exponent
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

  ! jjj = 0
  ! do k=0, kb
  !    do j=0, jb
  !       do i=0, ib
  !          if (iBlank(i,j,k) == -2 .or. iblank(i,j,k)==-3 .or. iblank(i,j,k)==0) then 
  !             ! Any cell in this cell's stencil 
              
  !             stencilLoop: do i_stencil=1, N_visc_drdw
  !                ii = visc_drdw_stencil(i_stencil, 1) + i
  !                jj = visc_drdw_stencil(i_stencil, 2) + j
  !                kk = visc_drdw_stencil(i_stencil, 3) + k

  !                ! Make sure we're on-block
  !                if (ii >=2 .and. ii <= il .and. jj >= 2 .and. jj<= jl .and. &
  !                     kk >=2 .and. kk <= kl) then 
  !                   if (iblank(ii, jj, kk) == 1) then 
  !                      iii = (kk-2)*nx*ny + (jj-2)*nx + (ii-2) + 1
  !                      fringes(ii, jj, kk)%quality = large
  !                      jjj = jjj + 1
  !                   end if
  !                end if
  !             end do stencilLoop
  !          end if
  !       end do
  !    end do
  ! end do
end subroutine initializeFringes

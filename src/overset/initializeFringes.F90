subroutine initializeFringes(oFringe, nn, level, sps)

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
  type(oversetFringe), intent(inout) :: oFringe
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


  ! Set the sizes for the oFringe and allocate the required space. 
  oFringe%il = il
  oFringe%jl = jl
  oFringe%kl = kl

  mm = nx*ny*nz
  allocate(& 
       oFringe%x(3, mm), &
       oFringe%quality(mm), &
       oFringe%myBlock(mm), &
       oFringe%myIndex(mm), &
       oFringe%donorProc(mm), &
       oFringe%donorBlock(mm), &
       oFringe%dI(mm), &
       oFringe%dJ(mm), &
       oFringe%dK(mm), &
       oFringe%donorFrac(3, mm), &
       oFringe%gInd(8, mm), &    
       oFringe%isWall(mm))

  ! These default set the entire array
  oFringe%myBlock = nn
  oFringe%donorProc = -1
  oFringe%dI = -1
  oFringe%dJ = -1
  oFringe%dK = -1
  oFringe%donorFrac = -one
  oFringe%gInd = -1
  oFringe%isWall = 0

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
           fringes(i, j, k)%isCompute = .False.
        end do
     end do
  end do

  ! Now loop over the actual compute cells, setting the cell center
  ! value 'x', the volume and flag these cells as compute
  ii = 0
  exponent = third
  do k=2, kl
     do j=2, jl
        do i=2, il
           fringes(i,j,k)%isCompute = .True. 
           ii = ii + 1
           do iDim=1, 3
              fringes(i, j, k)%x(iDim) = eighth*(&
                   x(i-1, j-1, k-1, iDim) + &
                   x(i  , j-1, k-1, iDim) + &
                   x(i-1, j  , k-1, iDim) + &
                   x(i  , j  , k-1, iDim) + &
                   x(i-1, j-1, k  , iDim) + &
                   x(i  , j-1, k  , iDim) + &
                   x(i-1, j  , k  , iDim) + &
                   x(i  , j  , k  , iDim))
              oFringe%x(iDim, ii) = fringes(i, j, k)%x(iDim)
           end do

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
              factor = 5.0
              fringes(i, j, k)%quality = (factor*vol(i, j, k))**exponent
           end if

           oFringe%quality(ii) = fringes(i, j, k)%quality
           oFringe%myIndex(ii) = ii
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
              ii = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
              oFringe%quality(ii) = large
           end if

        end do
     end do
  end do
  deallocate(tmp)


  jjj = 0
  do k=0, kb
     do j=0, jb
        do i=0, ib
           if (iBlank(i,j,k) == -2 .or. iblank(i,j,k)==-3 .or. iblank(i,j,k)==0) then 
              ! Any cell in this cell's stencil 
              
              stencilLoop: do i_stencil=1, N_visc_drdw
                 ii = visc_drdw_stencil(i_stencil, 1) + i
                 jj = visc_drdw_stencil(i_stencil, 2) + j
                 kk = visc_drdw_stencil(i_stencil, 3) + k

                 ! Make sure we're on-block
                 if (ii >=2 .and. ii <= il .and. jj >= 2 .and. jj<= jl .and. &
                      kk >=2 .and. kk <= kl) then 
                    if (iblank(ii, jj, kk) == 1) then 
                       iii = (kk-2)*nx*ny + (jj-2)*nx + (ii-2) + 1
                       oFringe%quality(iii) = large
                       fringes(ii, jj, kk)%quality = large
                       jjj = jjj + 1
                    end if
                 end if
              end do stencilLoop
           end if
        end do
     end do
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
                 fringes(i, j, k)%isWall  = .True. 
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

end subroutine initializeFringes

subroutine initializeOFringes(oFringe, nn)

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
  type(oversetFringe), intent(inout) :: oFringe
  integer(kind=intType), intent(in) :: nn

  ! Working Params
  integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj
  integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  logical :: wallsPresent, isWallType
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

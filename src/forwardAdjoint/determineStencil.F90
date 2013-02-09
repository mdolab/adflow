! This function attemps to "experimentally" determine the stencil for
! a given configuration. Basically, we're just going to peturb a cell
! in the "middle" of a block with FD, re-run the residual and see
! which cells in the output are change. This is (by defination) the
! stencil.

subroutine determineStencil(lumped)
  use inputDiscretization 
  use blockPointers
  use flowvarrefstate
  use inputPhysics
  use iteration
  implicit none
  logical :: lumped
  integer(kind=intType) :: ncellx,ncelly,ncellz
  integer(kind=intType) :: nnodex,nnodey,nnodez

  integer(kind=intType) :: i,j,k,l
  integer(kind=intType) :: icell,jcell,kcell,inode,jnode,knode
  real(kind=realType) :: delta_x
  real(kind=realType) :: dw_0(-3:3,-3:3,-3:3,nw)
  real(kind=realType) :: fw_0(-3:3,-3:3,-3:3,nw)
  integer(kind=intType) :: stencil(-3:3,-3:3,-3:3)
  logical :: different
  delta_x = 1.0
 
  if (lumped) then
     lumpedDiss=.True.
  end if

  stencil = 0
  call setPointers(1,1,1)
  ! Make sure the first block is big enough
  ncellx = il-2+1
  ncelly = jl-2+1
  ncellz = kl-2+1
  lumpedDiss = .True.

  if (ncellx < 5 .or. ncelly < 5 .or. ncellz < 5) then
     print *,'Block 1 is too small'
     stop
  end if
  
  ! Get the center cell
  icell = int((il + 2)/2)
  jcell = int((il + 2)/2)
  kcell = int((il + 2)/2)

  groundLevel = 1
  currentLevel = 1
  ! Get current dw
  !call block_res(1,1)
  
  ! Copy out dw
  do k=-3,3
     do j=-3,3
        do i=-3,3
           dw_0(i,j,k,:) = dw(icell+i,jcell+j,kcell+k,:)
           fw_0(i,j,k,:) = fw(icell+i,jcell+j,kcell+k,:)
        end do
     end do
  end do

  ! Set a peturbation in all states
  do l=1,nw
     w(icell,jcell,kcell,l) = w(icell,jcell,kcell,l) + delta_x
  end do

  ! Re-run dw
  !call block_res(1,1)

  do i=-3,3
     do j=-3,3
        do k=-3,3
           different = .False.
           do l=1,nw

              if (abs(dw(icell+i,jcell+j,kcell+k,l) - dw_0(i,j,k,l)) > 1e-12) then
                 different = .true. 
              end if
           end do
           if (different) then
              stencil(i,j,k) = 1
           end if
        end do
     end do
  end do

  print *,'State Stencil Consists of the following cells:',sum(stencil)
  do i=-3,3
     do j=-3,3
        do k=-3,3
           if (stencil(i,j,k) == 1) then
              print *,i,j,k
           end if
        end do
     end do
  end do

  ! Now also do it for the spatial stencil

  stencil = 0
  delta_x = .001
  ! Make sure the first block is big enough
  nnodex = il
  nnodey = jl
  nnodez = kl
  lumpedDiss = .False.

  if (nnodex < 5 .or. nnodey < 5 .or. nnodez < 5) then
     print *,'Block 1 is too small'
     stop
  end if
  
  ! Get the center cell
  inode = int((il + 1)/2)
  jnode = int((il + 1)/2)
  knode = int((il + 1)/2)

  groundLevel = 1
  currentLevel = 1

  ! Get current dw -> We can use the normal version here
  !call block_res(1,1)
  
  ! Copy out dw
  do k=-3,3
     do j=-3,3
        do i=-3,3
           dw_0(i,j,k,:) = dw(inode+i,jnode+j,knode+k,:)
        end do
     end do
  end do

  ! Set a peturbation in all DOF at coordinate inode,jnode,knode
  do l=1,3
     x(inode,jnode,knode,l) = x(inode,jnode,knode,l) + delta_x
  end do

  ! Re-run dw
  !call block_res_spatial(1,1)

  do i=-3,3
     do j=-3,3
        do k=-3,3
           different = .False.
           do l=1,nw

              if (abs(dw(inode+i,jnode+j,knode+k,l) - dw_0(i,j,k,l)) > 1e-12) then
                 different = .true. 
              end if
           end do
           if (different) then
              stencil(i,j,k) = 1
           end if
        end do
     end do
  end do

  print *,'Spatial Stencil Consists of the following cells:',sum(stencil)
  do i=-3,3
     do j=-3,3
        do k=-3,3
           if (stencil(i,j,k) == 1) then
              print *,i,j,k
           end if
        end do
     end do
  end do

  if (lumped) then
     lumpedDiss=.False.
  end if

end subroutine determineStencil

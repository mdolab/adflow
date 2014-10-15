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

  integer(kind=intType) :: i,j,k,l, nn
  integer(kind=intType) :: icell,jcell,kcell,inode,jnode,knode
  real(kind=realType) :: delta_x
  real(kind=realType) :: dw_0(-3:3,-3:3,-3:3,nw)
  real(kind=realType) :: fw_0(-3:3,-3:3,-3:3,nw)
  integer(kind=intType) :: stencil(-3:3,-3:3,-3:3)
  real(kind=realType) :: deriv(-3:3,-3:3,-3:3,nw,2)
  logical :: different
  real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor, Cavitation
  integer(kind=intType) :: liftIndex
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

  delta_x = one
 
  if (lumped) then
     lumpedDiss=.True.
  end if

  stencil = 0
  nn = 1
  call setPointers(nn,1,1)
  ! Make sure the first block is big enough
  ncellx = il-2+1
  ncelly = jl-2+1
  ncellz = kl-2+1
  lumpedDiss = .True.
  print *,'Block size:',il, jl, kl
  if (ncellx < 7 .or. ncelly < 7 .or. ncellz < 7) then
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
  call block_res(nn, 1, .False., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
  
! Copy out dw
  do k=3,3
     do j=3,3
        do i=3,3
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
  call block_res(nn, 1, .False., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
 
  do i=3,3
     do j=3,3
        do k=-3,3
           different = .False.
           do l=1,nw

              if (abs(dw(icell+i,jcell+j,kcell+k,l) - dw_0(i,j,k,l)) > 1e-13)then
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
  delta_x = 1.0
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
  jnode = int((jl + 1)/2)
  knode = int((kl + 1)/2)

  groundLevel = 1
  currentLevel = 1

  ! Get current dw -> We can use the normal version here
  call block_res(nn, 1, .True., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
  call block_res(nn, 1, .False., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
  dw_0=zero
  ! Copy out dw
  do k=-3,3
     do j=-3,3
        do i=-3,3
           dw_0(i,j,k,:) = dw(inode+i,jnode+j,knode+k,:)
        end do
     end do
  end do

  ! Set a peturbation in all DOF at coordinate inode,jnode,knode
  x(inode,jnode,knode,1) = x(inode,jnode,knode,1) + delta_x
  x(inode,jnode,knode,2) = x(inode,jnode,knode,2) + delta_x
  x(inode,jnode,knode,3) = x(inode,jnode,knode,3) + delta_x
  

  ! Re-run dw
  call block_res(nn, 1, .True., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
  do k=-3,3
     do j=-3,3
        do i=-3,3

           if (inode+i > il .or. inode + i < 2 .or. &
                jnode+j > jl .or. jnode + j < 2 .or. &
                knode+k > kl .or. knode +k < 2) then
              print *,'Block too small'
              stop
           end if
    
           different = .False.
           do l=1,nw

              if (abs(dw(inode+i,jnode+j,knode+k,l) - dw_0(i,j,k,l)) > 1e-12) then
                 different = .true. 
                 print *,'diff becuase of:'
                 print *, dw(inode+i,jnode+j,knode+k,l)
                 print *, dw_0(i,j,k,l)
              end if
           end do
           if (different) then
              stencil(i,j,k) = 1
              deriv(i,j,k,1,1) = sum((abs(dw(inode+i, jnode+j,knode+k,:) - dw_0(i,j,k,:)))/delta_x)
           end if
        end do
     end do
  end do

  print *,'Spatial Stencil Consists of the following cells:',sum(stencil)
  do i=-3,3
     do j=-3,3
        do k=-3,3
           if (stencil(i,j,k) == 1) then
              print *,i,j,k,deriv(i,j,k,1,1)
           end if
        end do
     end do
  end do

  if (lumped) then
     lumpedDiss=.False.
  end if

end subroutine determineStencil

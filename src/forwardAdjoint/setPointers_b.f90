!
! Set the pointers for the derivative values AND the normal pointers
subroutine setPointers_b(nn, level, sps)

  use blockPointers_b
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: nn,level,sps

  ! Set normal pointers
  call setPointers(nn, level, sps)
  
  viscSubfaceb => flowDomsb(nn,1,sps)%viscSubface

  xb    => flowDomsb(nn,1,sps)%x

  sib     => flowDomsb(nn,1,sps)%si
  sjb     => flowDomsb(nn,1,sps)%sj
  skb     => flowDomsb(nn,1,sps)%sk

  volb    => flowDomsb(nn,1,sps)%vol

  rotMatrixIb => flowDomsb(nn,1,sps)%rotMatrixI
  rotMatrixJb => flowDomsb(nn,1,sps)%rotMatrixJ
  rotMatrixKb => flowDomsb(nn,1,sps)%rotMatrixK

  sFaceIb => flowDomsb(nn,1,sps)%sFaceI
  sFaceJb => flowDomsb(nn,1,sps)%sFaceJ
  sFaceKb => flowDomsb(nn,1,sps)%sFaceK

  ! Flow variables. Note that wOld, gamma and the laminar viscosity
  ! point to the entries on the finest mesh. The reason is that
  ! they are computed from the other variables. For the eddy
  ! viscosity this is not the case because in a segregated solver
  ! its values are obtained from the fine grid level.

  wb     => flowDomsb(nn,1,sps)%w
  wOldb  => flowDomsb(nn,1,sps)%wOld
  pb     => flowDomsb(nn,1,sps)%p

  gammab => flowDomsb(nn,1,sps)%gamma
  rlvb   => flowDomsb(nn,1,sps)%rlv
  revb   => flowDomsb(nn,1,sps)%rev
  sb     => flowDomsb(nn,1,sps)%s

  ! Residual and multigrid variables. The residual point to the
  ! finest grid entry, the multigrid variables to their own level.

  dwb => flowDomsb(nn,1,sps)%dw
  fwb => flowDomsb(nn,1,sps)%fw

  w1b => flowDomsb(nn,1,sps)%w1
  wrb => flowDomsb(nn,1,sps)%wr

  ! Time-stepping variables and spectral radIi.
  ! They asps point to the fine mesh entry.

  wnb  => flowDomsb(nn,1,sps)%wn
  dtlb => flowDomsb(nn,1,sps)%dtl

  radIb => flowDomsb(nn,1,sps)%radI
  radJb => flowDomsb(nn,1,sps)%radJ
  radKb => flowDomsb(nn,1,sps)%radK

  d2Wallb => flowDomsb(nn,1,sps)%d2Wall

  ! Arrays used for the implicit treatment of the turbulent wasps
  ! boundary conditions. As these variables are only aspocated for
  ! the 1st spectral solution of the fine mesh, the pointers point
  ! to those arrays.

  bmti1b => flowDomsb(nn,1,1)%bmti1
  bmti2b => flowDomsb(nn,1,1)%bmti2
  bmtj1b => flowDomsb(nn,1,1)%bmtj1
  bmtj2b => flowDomsb(nn,1,1)%bmtj2
  bmtk1b => flowDomsb(nn,1,1)%bmtk1
  bmtk2b => flowDomsb(nn,1,1)%bmtk2

  bvti1b => flowDomsb(nn,1,1)%bvti1
  bvti2b => flowDomsb(nn,1,1)%bvti2
  bvtj1b => flowDomsb(nn,1,1)%bvtj1
  bvtj2b => flowDomsb(nn,1,1)%bvtj2
  bvtk1b => flowDomsb(nn,1,1)%bvtk1
  bvtk2b => flowDomsb(nn,1,1)%bvtk2

  !BCData Array
  BCDatab => flowDomsb(nn,1,sps)%BCdata
  
end subroutine setPointers_b



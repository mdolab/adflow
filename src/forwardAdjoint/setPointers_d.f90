!
! Set the pointers for the derivative values AND the normal pointers
subroutine setPointers_d(nn, level, sps)

  use blockPointers
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: nn,level,sps

  ! Set normal pointers
  call setPointers(nn, level, sps)
  
  ! Set shockSensor pointer IN FLOWDOMS!
  shockSensor => flowDoms(nn,1,sps)%shockSensor

  viscSubfaced => flowDomsd(nn,1,sps)%viscSubface

  xd    => flowDomsd(nn,1,sps)%x

  sid     => flowDomsd(nn,1,sps)%si
  sjd     => flowDomsd(nn,1,sps)%sj
  skd     => flowDomsd(nn,1,sps)%sk

  vold    => flowDomsd(nn,1,sps)%vol

  rotMatrixId => flowDomsd(nn,1,sps)%rotMatrixI
  rotMatrixJd => flowDomsd(nn,1,sps)%rotMatrixJ
  rotMatrixKd => flowDomsd(nn,1,sps)%rotMatrixK

  sFaceId => flowDomsd(nn,1,sps)%sFaceI
  sFaceJd => flowDomsd(nn,1,sps)%sFaceJ
  sFaceKd => flowDomsd(nn,1,sps)%sFaceK

  ! Flow variables. Note that wOld, gamma and the laminar viscosity
  ! point to the entries on the finest mesh. The reason is that
  ! they are computed from the other variables. For the eddy
  ! viscosity this is not the case because in a segregated solver
  ! its values are obtained from the fine grid level.

  wd     => flowDomsd(nn,1,sps)%w
  wOldd  => flowDomsd(nn,1,sps)%wOld
  pd     => flowDomsd(nn,1,sps)%p

  gammad => flowDomsd(nn,1,sps)%gamma
  rlvd   => flowDomsd(nn,1,sps)%rlv
  revd   => flowDomsd(nn,1,sps)%rev
  sd     => flowDomsd(nn,1,sps)%s

  uxd => flowDomsd(nn,1,sps)%ux
  uyd => flowDomsd(nn,1,sps)%uy
  uzd => flowDomsd(nn,1,sps)%uz
  
  vxd => flowDomsd(nn,1,sps)%vx
  vyd => flowDomsd(nn,1,sps)%vy
  vzd => flowDomsd(nn,1,sps)%vz
  
  wxd => flowDomsd(nn,1,sps)%wx
  wyd => flowDomsd(nn,1,sps)%wy
  wzd => flowDomsd(nn,1,sps)%wz
  
  qxd => flowDomsd(nn,1,sps)%qx
  qyd => flowDomsd(nn,1,sps)%qy
  qzd => flowDomsd(nn,1,sps)%qz

  ! Residual and multigrid variables. The residual point to the
  ! finest grid entry, the multigrid variables to their own level.

  dwd => flowDomsd(nn,1,sps)%dw
  fwd => flowDomsd(nn,1,sps)%fw

  w1d => flowDomsd(nn,1,sps)%w1
  wrd => flowDomsd(nn,1,sps)%wr

  ! Time-stepping variables and spectral radIi.
  ! They asps point to the fine mesh entry.

  wnd  => flowDomsd(nn,1,sps)%wn
  dtld => flowDomsd(nn,1,sps)%dtl

  radId => flowDomsd(nn,1,sps)%radI
  radJd => flowDomsd(nn,1,sps)%radJ
  radKd => flowDomsd(nn,1,sps)%radK

  d2Walld => flowDomsd(nn,1,sps)%d2Wall

  ! Arrays used for the implicit treatment of the turbulent wasps
  ! boundary conditions. As these variables are only aspocated for
  ! the 1st spectral solution of the fine mesh, the pointers point
  ! to those arrays.

  bmti1d => flowDomsd(nn,1,1)%bmti1
  bmti2d => flowDomsd(nn,1,1)%bmti2
  bmtj1d => flowDomsd(nn,1,1)%bmtj1
  bmtj2d => flowDomsd(nn,1,1)%bmtj2
  bmtk1d => flowDomsd(nn,1,1)%bmtk1
  bmtk2d => flowDomsd(nn,1,1)%bmtk2

  bvti1d => flowDomsd(nn,1,1)%bvti1
  bvti2d => flowDomsd(nn,1,1)%bvti2
  bvtj1d => flowDomsd(nn,1,1)%bvtj1
  bvtj2d => flowDomsd(nn,1,1)%bvtj2
  bvtk1d => flowDomsd(nn,1,1)%bvtk1
  bvtk2d => flowDomsd(nn,1,1)%bvtk2

  !BCData Array
  BCDatad => flowDomsd(nn,1,sps)%BCdata
  
end subroutine setPointers_d



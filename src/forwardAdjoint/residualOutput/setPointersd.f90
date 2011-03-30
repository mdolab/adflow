!
!      ******************************************************************
!      *                                                                *
!      * File:          setPointers.f90                                 *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine setPointersd(sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * setPointers makes the variables in blockPointers point to      *
  !      * block nn for grid level mm and spectral solution ll.           *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers_d
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: sps

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Store the info of the current block, such that inside the
  ! module blockPointers it is known to which block the data
  ! belongs.

  viscSubfaced => flowDomsd(sps)%viscSubface


  xd    => flowDomsd(sps)%x

  sid     => flowDomsd(sps)%si
  sjd     => flowDomsd(sps)%sj
  skd     => flowDomsd(sps)%sk

  rotMatrixId => flowDomsd(sps)%rotMatrixI
  rotMatrixJd => flowDomsd(sps)%rotMatrixJ
  rotMatrixKd => flowDomsd(sps)%rotMatrixK

  sFaceId => flowDomsd(sps)%sFaceI
  sFaceJd => flowDomsd(sps)%sFaceJ
  sFaceKd => flowDomsd(sps)%sFaceK

  ! Flow variables. Note that wOld, gamma and the laminar viscosity
  ! point to the entries on the finest mesh. The reason is that
  ! they are computed from the other variables. For the eddy
  ! viscosity this is not the case because in a segregated solver
  ! its values are obtained from the fine grid level.

  wd     => flowDomsd(sps)%w
  wOldd  => flowDomsd(sps)%wOld
  pd     => flowDomsd(sps)%p

  gammad => flowDomsd(sps)%gamma
  rlvd   => flowDomsd(sps)%rlv
  revd   => flowDomsd(sps)%rev
  sd     => flowDomsd(sps)%s

  ! Residual and multigrid variables. The residual point to the
  ! finest grid entry, the multigrid variables to their own level.

  dwd => flowDomsd(sps)%dw
  fwd => flowDomsd(sps)%fw

  w1d => flowDomsd(sps)%w1
  wrd => flowDomsd(sps)%wr

  ! Time-stepping variables and spectral radIi.
  ! They asps point to the fine mesh entry.

  wnd  => flowDomsd(sps)%wn
  dtld => flowDomsd(sps)%dtl

  radId => flowDomsd(sps)%radI
  radJd => flowDomsd(sps)%radJ
  radKd => flowDomsd(sps)%radK

  d2Walld => flowDomsd(sps)%d2Wall

  ! Arrays used for the implicit treatment of the turbulent wasps
  ! boundary conditions. As these variables are only aspocated for
  ! the 1st spectral solution of the fine mesh, the pointers point
  ! to those arrays.

  bmti1d => flowDomsd(1)%bmti1
  bmti2d => flowDomsd(1)%bmti2
  bmtj1d => flowDomsd(1)%bmtj1
  bmtj2d => flowDomsd(1)%bmtj2
  bmtk1d => flowDomsd(1)%bmtk1
  bmtk2d => flowDomsd(1)%bmtk2

  bvti1d => flowDomsd(1)%bvti1
  bvti2d => flowDomsd(1)%bvti2
  bvtj1d => flowDomsd(1)%bvtj1
  bvtj2d => flowDomsd(1)%bvtj2
  bvtk1d => flowDomsd(1)%bvtk1
  bvtk2d => flowDomsd(1)%bvtk2

end subroutine setPointersd

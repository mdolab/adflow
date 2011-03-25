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
  use blockPointers
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

  viscSubfaced => flowDoms(sps)%viscSubface


  xd    => flowDomsd(sps)%x

  sid     => flowDoms(sps)%si
  sjd     => flowDoms(sps)%sj
  skd     => flowDoms(sps)%sk

  rotMatrixId => flowDoms(sps)%rotMatrixI
  rotMatrixJd => flowDoms(sps)%rotMatrixJ
  rotMatrixKd => flowDoms(sps)%rotMatrixK

  sFaceId => flowDoms(sps)%sFaceI
  sFaceJd => flowDoms(sps)%sFaceJ
  sFaceKd => flowDoms(sps)%sFaceK

  ! Flow variables. Note that wOld, gamma and the laminar viscosity
  ! point to the entries on the finest mesh. The reason is that
  ! they are computed from the other variables. For the eddy
  ! viscosity this is not the case because in a segregated solver
  ! its values are obtained from the fine grid level.

  wd     => flowDoms(sps)%w
  wOldd  => flowDoms(sps)%wOld
  pd     => flowDoms(sps)%p

  gammad => flowDoms(sps)%gamma
  rlvd   => flowDoms(sps)%rlv
  revd   => flowDoms(sps)%rev
  sd     => flowDoms(sps)%s

  ! Residual and multigrid variables. The residual point to the
  ! finest grid entry, the multigrid variables to their own level.

  dwd => flowDoms(sps)%dw
  fwd => flowDoms(sps)%fw

  w1d => flowDoms(sps)%w1
  wrd => flowDoms(sps)%wr

  ! Time-stepping variables and spectral radIi.
  ! They asps point to the fine mesh entry.

  wnd  => flowDoms(sps)%wn
  dtld => flowDoms(sps)%dtl

  radId => flowDoms(sps)%radI
  radJd => flowDoms(sps)%radJ
  radKd => flowDoms(sps)%radK

  d2Walld => flowDoms(sps)%d2Wall

  ! Arrays used for the implicit treatment of the turbulent wasps
  ! boundary conditions. As these variables are only aspocated for
  ! the 1st spectral solution of the fine mesh, the pointers point
  ! to those arrays.

  bmti1d => flowDoms(1)%bmti1
  bmti2d => flowDoms(1)%bmti2
  bmtj1d => flowDoms(1)%bmtj1
  bmtj2d => flowDoms(1)%bmtj2
  bmtk1d => flowDoms(1)%bmtk1
  bmtk2d => flowDoms(1)%bmtk2

  bvti1d => flowDoms(1)%bvti1
  bvti2d => flowDoms(1)%bvti2
  bvtj1d => flowDoms(1)%bvtj1
  bvtj2d => flowDoms(1)%bvtj2
  bvtk1d => flowDoms(1)%bvtk1
  bvtk2d => flowDoms(1)%bvtk2

end subroutine setPointersd

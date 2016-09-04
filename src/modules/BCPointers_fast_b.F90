
module BCPointers_fast_b

! Thiss module contains data structures used to apply BCs. 

  use constants, only : intType, realType
  implicit none
  save

  real(kind=realType), dimension(:,:,:), pointer :: ww0, ww1, ww2, ww3
  real(kind=realType), dimension(:,:)  , pointer :: pp0, pp1, pp2, pp3
  real(kind=realType), dimension(:,:)  , pointer :: rlv0, rlv1, rlv2, rlv3
  real(kind=realType), dimension(:,:)  , pointer :: rev0, rev1, rev2, rev3
  real(kind=realType), dimension(:,:  ), pointer :: gamma0, gamma1, gamma2, gamma3
  real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
  real(kind=realType), dimension(:,:,:), pointer :: ss, xx
  real(kind=realType), dimension(:,:), pointer :: dd2wall, sface
  integer(kind=intType), dimension(:,:), pointer :: gcp

  integer(kind=intType) :: iStart, iEnd, iSize
  integer(kind=intType) :: jStart, jEnd, jSize

  real(kind=realType), dimension(:,:,:), pointer :: ww0d, ww1d, ww2d, ww3d
  real(kind=realType), dimension(:,:),   pointer :: pp0d, pp1d, pp2d, pp3d
  real(kind=realType), dimension(:,:),   pointer :: rlv0d, rlv1d, rlv2d, rlv3d
  real(kind=realType), dimension(:,:),   pointer :: rev0d, rev1d, rev2d, rev3d
  real(kind=realType), dimension(:,:),   pointer :: gamma0d, gamma1d, gamma2d, gamma3d
  real(kind=realType), dimension(:,:,:), pointer :: ssid, ssjd, sskd
  real(kind=realType), dimension(:,:,:), pointer :: ssd, xxd
  real(kind=realType), dimension(:,:),   pointer :: dd2walld, sfaced

end module BCPointers_fast_b


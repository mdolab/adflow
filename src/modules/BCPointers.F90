
module BCPointers

! Thiss module contains data structures used to apply BCs.

  use constants, only : intType, realType
  implicit none
  save
#if !defined USE_TAPENADE || defined TAPENADE_POINTERS || defined TAPENADE_FORWARD
  real(kind=realType), dimension(:,:,:), pointer :: ww0, ww1, ww2, ww3
  real(kind=realType), dimension(:,:),   pointer :: pp0, pp1, pp2, pp3
  real(kind=realType), dimension(:,:),   pointer :: rlv0, rlv1, rlv2, rlv3
  real(kind=realType), dimension(:,:),   pointer :: rev0, rev1, rev2, rev3
  real(kind=realType), dimension(:,:),   pointer :: gamma0, gamma1, gamma2, gamma3
  real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
  real(kind=realType), dimension(:,:,:), pointer :: ss, xx
  real(kind=realType), dimension(:,:),   pointer :: dd2wall, sFace
  integer(kind=intType), dimension(:,:), pointer :: gcp
#else
  real(kind=realType), dimension(:,:,:), allocatable :: ww0, ww1, ww2, ww3
  real(kind=realType), dimension(:,:)  , allocatable :: pp0, pp1, pp2, pp3
  real(kind=realType), dimension(:,:)  , allocatable :: rlv0, rlv1, rlv2, rlv3
  real(kind=realType), dimension(:,:)  , allocatable :: rev0, rev1, rev2, rev3
  real(kind=realType), dimension(:,:  ), allocatable :: gamma0, gamma1, gamma2, gamma3
  real(kind=realType), dimension(:,:,:), allocatable :: ssi, xx
  integer(kind=intType), dimension(:,:), allocatable :: gcp
#endif
  integer(kind=intType) :: iStart, iEnd, iSize
  integer(kind=intType) :: jStart, jEnd, jSize

end module BCPointers

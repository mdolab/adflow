module turbMod
  !
  !       This local module contains variables used when the turbulence  
  !       equations are solved.                                          
  !
  use precision
  implicit none
  save

  ! secondOrd:  whether or not a second order discretization for
  !             the advective terms must be used.
  ! sig1, sig2: Sigma coefficients in the diffusion terms of the
  !             different turbulence models.

  logical :: secondOrd
  real(kind=realType) :: sig1, sig2

  ! dvt:     Pointer, which points to an unused part of dw. It is
  !          used for temporary storage of residual.
  ! vort:    Pointer, which points to an unused part of dw. It is
  !          used for temporary storage of the magnitude of
  !          vorticity squared.
  ! prod:    Pointer, which points to an unused part of dw. It is
  !          used for temporary storage of the unscaled production
  !          term.
  ! f1:      F1 blending function in the SST model.
  ! kwCD:   Cross diffusion term in the k-omega type models.
  ! ktCD:   Cross diffusion term in the k-tau model
  ! sct:     Time scale in the v2-f model.
  ! scl2:    Length scale in the v2-f model.
  ! strain2: Square of the strain.

  real(kind=realType), dimension(:,:,:,:), pointer :: dvt
  real(kind=realType), dimension(:,:,:),   pointer :: vort
  real(kind=realType), dimension(:,:,:),   pointer :: prod
  real(kind=realType), dimension(:,:,:),   pointer :: f1
  real(kind=realType), dimension(:,:,:),   pointer :: kwCD
  real(kind=realType), dimension(:,:,:),   pointer :: ktCD
  real(kind=realType), dimension(:,:,:),   pointer :: sct
  real(kind=realType), dimension(:,:,:),   pointer :: scl2
  real(kind=realType), dimension(:,:,:),   pointer :: strain2

#ifndef USE_TAPENADE
  real(kind=realType), dimension(:,:,:,:), pointer :: dvtd
  real(kind=realType), dimension(:,:,:),   pointer :: vortd
  real(kind=realType), dimension(:,:,:),   pointer :: prodd
  real(kind=realType), dimension(:,:,:),   pointer :: f1d
  real(kind=realType), dimension(:,:,:),   pointer :: kwCDd
  real(kind=realType), dimension(:,:,:),   pointer :: ktCDd
  real(kind=realType), dimension(:,:,:),   pointer :: sctd
  real(kind=realType), dimension(:,:,:),   pointer :: scl2d
  real(kind=realType), dimension(:,:,:),   pointer :: strain2d
#endif
end module turbMod

#ifndef USE_TAPENADE
module turbmod_d
  use turbmod
end module turbmod_d
module turbmod_b
  use turbmod
end module turbmod_b

#endif

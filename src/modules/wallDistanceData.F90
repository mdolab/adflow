module wallDistanceData

  ! This module stores some additional data required for the fast wall
  ! distance computations.

  use constants
#ifndef USE_TAPENADE
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#endif

  ! xVolume: flatten 1-D vector of all volume nodes for all
  ! spectral instances. The dimension is the number of levels

  ! xSurf: flatten 1-D vector of the surface nodes for the faces
  ! that individual processors require for doing its own wall
  ! distance calculation

  ! wallScatter: The vecScatter objects that select the nodes
  ! from xVolume and desposit them into xSurf. This is the
  ! forward operation. The reverse operation is used for the
  ! sensitivities.

  ! wallDistanceDataAllocated : Logical array keeping track of
  ! whether or not the petsc data is allocated.

  ! indicesForSPS: A simple derived type for keeping track of
  ! indices while doing wall distance computation.


  real(kind=realType), dimension(:), pointer :: xSurf

  logical, dimension(:), allocatable  :: wallDistanceDataAllocated
  logical, dimension(:), allocatable ::  updateWallAssociation

#ifndef USE_TAPENADE
  real(kind=realType), dimension(:), pointer :: xSurfd
  Vec , dimension(:), allocatable :: xVolumeVec
  Vec , dimension(:, :), allocatable :: xSurfVec
  Vec , dimension(:), allocatable :: xSurfVecd
  VecScatter, dimension(:, :), allocatable :: wallScatter
  IS IS1, IS2
  real(kind=realType), dimension(:), pointer :: xVolume
#endif

end module wallDistanceData

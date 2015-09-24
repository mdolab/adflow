module overset

  use precision
  implicit none

#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Helper data type for overset communication
  type variablePointer
     real(kind=realType), dimension(:, :, :), pointer :: arr
  end type variablePointer

  ! Helper dataType for communicated overset grid points. This data
  ! structure mirrros the blockType structure in block.F90, but only
  ! contains minimum amount of information required for computing
  ! overset connectivities. 

  type oversetBlock
     integer(kind=intType) :: globalID
     integer(kind=intType) :: ie, je, ke
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: nx, ny, nz
     integer(kind=intType), dimension(:, :), pointer :: hexaConn
     real(kind=realType), dimension(:, :), pointer :: xDual, xSearch
     real(kind=realType), dimension(:, :,:, :), pointer :: xDual2, xSearch2
     real(kind=realType), dimension(:, :), pointer :: qualDonor, qualRecv
     integer(kind=intType), dimension(:,:,:), pointer :: iblank
     integer(kind=intType), dimension(:, :, :), pointer :: globalCell
     logical, dimension(:, :, :), pointer :: forceRecv
     character(len=15) :: adtName
     integer(kind=intTYpe) :: nFringe
     integer(kind=intTYpe) :: nDonor
     integer(kind=intType), dimension(:, :), pointer :: fringeIndices
     integer(kind=intType), dimension(:, :), pointer :: donorIndices
     integer(kind=intType), dimension(:, :, :), pointer :: donorIndices2
     real(kind=realType), dimension(:, :), pointer :: donorFrac
     integer(kind=intType), dimension(:, :), pointer :: ind
     
  end type oversetBlocK

  type(oversetBlock), dimension(:), allocatable :: oBlocks

  real(kind=realType), dimension(:, :), allocatable :: nodesVisc, normVisc
  type(variablePointer), dimension(:, :), allocatable :: variables

  ! Two vectors for overset communication
  Vec oversetDonors
  Vec oversetFringes
  
  ! The vecscatter context for the communication
  VecScatter oversetScatter
  
  ! Temporary index sets
  IS IS1, IS2

end module overset

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

  type oversetDonor

     ! Make everything in here static such that we can potentially use
     ! MPI to send these types around direclty

     integer(kind=intType) :: donorProcID
     integer(kind=intType) :: donorBlockID
     real(kind=realType), dimension(3) :: frac
     integer(kind=intType), dimension(3) :: ind

  end type oversetDonor

  type oversetBlock
     integer(kind=intType) :: cluster
     integer(kind=intType) :: ib, jb, kb
     integer(kind=intType) :: ie, je, ke
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: nx, ny, nz
     integer(kind=intType), dimension(:, :), pointer :: hexaConn
     real(kind=realType), dimension(:, :, :, :), pointer :: x ! Orignal X coordinates
     real(kind=realType), dimension(:, :), pointer :: xDual, xPrimal
     real(kind=realType), dimension(:, :, :, :), pointer :: xSearch
     real(kind=realType), dimension(:, :, :), pointer :: qualDonor, qualRecv
     integer(kind=intType), dimension(:,:,:), pointer :: iblank
     integer(kind=intType), dimension(:,:,:), pointer :: cellStatus
     integer(kind=intType), dimension(:, :, :), pointer :: globalCell
     logical, dimension(:, :, :), pointer :: recvStatus
     logical, dimension(:, :, :), pointer :: donorStatus
     logical, dimension(:, :, :), pointer :: forceRecv
     logical, dimension(:, :, :), pointer :: invalidDonor
     real(kind=realType) :: minVolume
     integer :: globalBlockID, localBlockID
     character(len=15) :: adtName
     integer(kind=intTYpe) :: nFringe
     integer(kind=intTYpe) :: nDonor
     integer(kind=intType), dimension(:, :), pointer :: fringeIndices
     integer(kind=intType), dimension(:, :), pointer :: donorIndices
     real(kind=realType), dimension(:, :), pointer :: donorFrac
     integer(kind=intType) :: procID
     type(oversetDonor), dimension(:, :, :), pointer :: donors
  end type oversetBlocK

  type(oversetBlock), dimension(:), allocatable :: oBlocks
  type(variablePointer), dimension(:, :), allocatable :: variables

  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  integer(kind=intType) :: nDomTotal
  ! Two vectors for overset communication
  Vec oversetDonors
  Vec oversetFringes
  
  ! The vecscatter context for the communication
  VecScatter oversetScatter
  
  ! Temporary index sets
  IS IS1, IS2

  ! Variales for tracking cots
  integer(kind=intType) :: totalSearches
  real(kind=realType) :: searchCosts

end module overset

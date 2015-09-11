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
  
  type(variablePointer), dimension(:, :), allocatable :: variables

  ! Two vectors for overset communication
  Vec oversetDonors
  Vec oversetFringes
  
  ! The vecscatter context for the communication
  VecScatter oversetScatter
  
  ! Temporary index sets
  IS IS1, IS2
end module overset

module surfaceFamilies

  use constants
#ifndef USE_TAPENADE
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none


  type familyExchange
     ! Vectors for global traction calc

     ! Parallel vector of un-uniqufied values concatenated for all
     ! included surfaces
     Vec nodeValLocal

     ! Parallel vector of uniqueifed values.
     Vec nodeValGlobal

     ! Sum global. Same size as nodeValGlobal
     Vec sumGlobal

     ! Scatter from nodeValLocal to NodeValGlobal
     VecScatter scatter

     ! Flag for allocated petsc variables
     logical :: allocated =.False.
     integer(kind=intType), dimension(:) , allocatable :: famList

     integer(Kind=intType) :: sps
     integer(kind=intType) :: nNodes
  end type familyExchange

  type BCGroupType
     integer(kind=intType), pointer, dimension(:) :: famList
  end type BCGroupType

  ! Generic PETSc scatters
  IS IS1, IS2

  ! The list of exchanges based on boundary condition type
  type(familyExchange), dimension(:, :), allocatable, target :: BCFamExchange

  ! List of familis grouped by BC. See constants.F90 for the indices
  ! to use for this array.
  type(BCGroupType), dimension(nFamExchange) :: BCFamGroups

  ! The full list of the family names
  character(len=maxCGNSNameLen), dimension(:), allocatable :: famNames

  ! List of all families. This is just 1,2,3,4...nFam. It is just used
  ! in fortran when a specific family is not required.
  integer(kind=intType), dimension(:), allocatable :: fullFamList
#endif

  ! Special BC array's that are sometime required for reducitons.
  real(kind=realType), dimension(:, :), allocatable, target :: zeroCellVal
  real(kind=realType), dimension(:, :), allocatable, target :: oneCellVal
  real(kind=realType), dimension(:, :), allocatable, target :: zeroNodeVal


#ifndef USE_TAPENADE
  contains

    subroutine getnfam(nfam)

      implicit none
      integer(kind=inttype), intent(out) :: nfam
      if (allocated(famnames)) then
         nfam = size(famnames)
      else
         nfam = 0
      end if

    end subroutine getnfam

    subroutine getfam(i, fam)
      implicit none
      character(len=maxCGNSNameLen), intent(out) :: fam
      integer(kind=intType), intent(in) :: i

      if (allocated(famnames)) then
         fam = famnames(i)
      end if

    end subroutine getfam

    subroutine destroyFamilyExchange(exch)

      use constants
      type(familyExchange) :: exch
      integer(kind=intType) :: ierr
      if (exch%allocated) then

           call vecDestroy(exch%nodeValLocal, ierr)
           call vecDestroy(exch%nodeValGlobal, ierr)
           call vecDestroy(exch%sumGlobal, ierr)
           call vecScatterDestroy(exch%scatter, ierr)

        end if

        exch%allocated = .False.

      end subroutine destroyFamilyExchange
#endif
end module surfaceFamilies

module surfaceFamilies

  use constants
  implicit none
  save
#ifndef USE_TAPENADE
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif
  type familyExchange 
     ! Vectors for global traction calc
     Vec nodeValLocal
     Vec nodevalGlobal
     Vec localWeight
     Vec sumGlobal
     VecScatter scatter
     logical :: allocated =.False.

     integer(kind=intType), dimension(:) , allocatable :: famGroups
     integer(kind=intType) :: nFam
     integer(Kind=intType) :: sps
     ! Data required for exchange
     integer(kind=intType) :: nNodes
     real(kind=realType), dimension(:, :), allocatable :: nodalValues
     integer(kind=intType), allocatable, dimension(:, :) :: conn
     integer(kind=intType), allocatable, dimension(:) :: elemFam
     real(kind=realType),   allocatable, dimension(:) :: fc

  end type familyExchange

  IS IS1, IS2


  ! The list of family exchanges for each individual family
  type(familyExchange) , dimension(:, :), allocatable, target :: familyExchanges

  ! All the walls get their own exchange that is reduced over all
  ! walls. This is used for the lift distribution/slice data to make
  ! sure all data is consistently reduced. 
  type(familyExchange), dimension(:), allocatable, target :: wallExchange

#endif
  integer(kind=intType), dimension(:), allocatable :: famGroups
  character(len=maxCGNSNameLen), dimension(:), allocatable :: famNames
  integer(kind=intType), dimension(:), allocatable :: wallFamilies, famIsWall
  integer(kind=intType) :: totalFamilies, totalWallFamilies

#ifndef USE_TAPENADE
  contains

    subroutine destroyFamilyExchange(exch) 

      use constants
      type(familyExchange) :: exch
      integer(kind=intType) :: ierr
      if (exch%allocated) then

           call vecDestroy(exch%nodeValLocal, ierr)
           call vecDestroy(exch%localWeight, ierr)
           call vecDestroy(exch%nodeValGlobal, ierr)
           call vecDestroy(exch%sumGlobal, ierr)
           call vecScatterDestroy(exch%scatter, ierr)
        end if

        if (allocated(exch%nodalvalues)) &
             deallocate(exch%nodalValues)
        
        if (allocated(exch%conn)) &
             deallocate(exch%conn)
        
        if (allocated(exch%elemFam)) & 
             deallocate(exch%elemFam)
        
        if (allocated(exch%fc)) & 
             deallocate(exch%fc)
        
        if (allocated(famGroups)) &
             deallocate(famGroups)
        
        exch%allocated = .False.
        
      end subroutine destroyFamilyExchange
#endif
end module surfaceFamilies

module overset

  use precision
  use adtData
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

  ! Store the coordinates from a block that will need to be searched.

  ! x , size(3, arrSize) : Seach point. 
  ! fInd, size(3, arrSize) : i,j,k block indices of the point
  ! dInd, size(4, arrSize) : blockID, i,j,k indices of the point's donor
  ! frac, size(3, arrSize) : Fraction weights in the cell
  ! gInd, size(8, arrSize) : Global indices of the 8 donor cells
  ! qualRecv, size(1, arraySize) : qualRecv of point i,j,k, usually
  !                                donorQual of the donor to point i,j,k
  type oversetSearchCoords
     real(kind=realType), dimension(:, :), pointer :: x, frac
     integer(kind=intType), dimension(:, :), pointer :: fInd, dInd, gInd
     integer(kind=intType) :: n, arrSize
     real(kind=realType), dimension(:, :), pointer :: qualRecv

     ! Buffer space
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer
  end type oversetSearchCoords

  type specialSearchCoords
     ! procId, size(3, arrSize) : size(fInd), saves procId of
     !                            incoming oversetSearchCoords blocks
     integer(kind=intType), dimension(:, :), pointer :: procId
  end type specialSearchCoords

  ! A simple generic sparse matrix storage container for storing the
  ! (sparse) overlap structure of an overset mesh
  type CSRMatrix
     integer(kind=intType) :: nRow, nCol, nnz
     integer(kind=intType), dimension(:), pointer :: colInd, rowPtr
     real(kind=realType), dimension(:), pointer :: data
     integer(Kind=intType), dimension(:), pointer :: assignedProc
  end type CSRMatrix

  ! This derived type contains sufficient information to perfom ADT
  ! donor searches. The information is com
  type oversetBlock

     ! Sizes for the block
     integer(kind=intType) :: ib, jb, kb
     integer(kind=intType) :: ie, je, ke
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: nx, ny, nz
     
     ! This is the cell volume of the donor
     real(kind=realType), dimension(:, :), pointer :: qualDonor

     ! Connectivity for the ADT
     integer(kind=intType), dimension(:, :), pointer :: hexaConn

     ! Coordinates for the ADT
     real(kind=realType), dimension(:, :), pointer :: xADT

     ! Flag eliminating a cell as a potential donor. Using an integer
     ! since it is easier to transfer
     integer, dimension(:, :, :), pointer :: invalidDonor

     ! Copy of global cell
     integer(kind=intType), dimension(:, :, :), pointer :: globalCell
     
     ! Minimum volume for this block
     real(kind=realType) :: minVolume

     ! The ADT for this block
     type(adtType) :: ADT

     ! Buffer space
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! Additional variables for donor info
     ! ------------------------------------
     integer(kind=intType), dimension(:,:,:), pointer :: iblank

     ! Status of a cell based on oversetbc info
     integer(kind=intType), dimension(:,:,:), pointer :: cellStatus

     ! Valid status of donor cells
     logical, dimension(:, :, :), pointer :: recvStatus
     logical, dimension(:, :, :), pointer :: donorStatus
     logical, dimension(:, :, :), pointer :: forceRecv
     ! ------------------------------------

  end type oversetBlocK

  type(oversetBlock), dimension(:), allocatable :: oBlocks
  type(variablePointer), dimension(:, :), allocatable :: variables

  type(oversetSearchCoords), dimension(:), allocatable :: searchCoords
  type(specialSearchCoords), dimension(:), allocatable :: specialArray

  ! temporary search coordinates to save. size= (nDom,0:nProc-1),
  ! nDom=total compute blocks in my proc
  type(oversetSearchCoords), dimension(:,:), allocatable :: tmpsearchCoords

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

end module overset

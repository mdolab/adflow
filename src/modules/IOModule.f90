!
!      ******************************************************************
!      *                                                                *
!      * File:          IOModule.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-02-2005                                      *
!      * Last modified: 10-31-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module IOModule
!
!      ******************************************************************
!      *                                                                *
!      * Constants and variables used in the IO routines.               *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the parameters.                                  *
!      *                                                                *
!      ******************************************************************
!
       ! The maximum amount of data a processor can read in one call
       ! to su_file_read. This value is 2Gbyte.
 
       integer(kind=mpi_offset_kind), parameter :: &
                                  maxSizeIO = 2147483648_mpi_offset_kind

       ! Definition of the possibilities for storing data.

       integer(kind=intType), parameter :: nodeData          = 1
       integer(kind=intType), parameter :: cellDataNoHalo    = 2
       integer(kind=intType), parameter :: cellDataPlusHalo  = 3

!
!      ******************************************************************
!      *                                                                *
!      * Definition of the derived datatype IOType, which is used to    *
!      * to make the IO as general as needed.                           *
!      *                                                                *
!      ******************************************************************
!
       type IOType

         ! pointerOffset: offset due to the usage of a pointer to a
         !                subarray. The Fortran standard is such that
         !                the starting indices of the pointer array
         !                is 1, no matter what the original starting
         !                index is. This can lead to a shift in the
         !                indices.
         ! w:             The variable(s) to be read/written.

         integer(kind=intType) :: pointerOffset

         real(kind=realType), dimension(:,:,:,:), pointer :: w

       end type IOType

!
!      ******************************************************************
!      *                                                                *
!      * Definition of the variables used for both CGNS                 *
!      *                                                                *
!      ******************************************************************
!
       ! IOVar(nDom,nIOFiles): Array of the derived datatype IOType to
       !                       facilitate a general IO implementation.

       type(IOType), dimension(:,:), allocatable :: IOVar

       end module IOModule

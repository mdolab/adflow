!
!      ******************************************************************
!      *                                                                *
!      * File:          determineIndirectHalos.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-31-2003                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineIndirectHalos(nHalo, iihalo, entityHalo,  &
                                         transform, entityIndex,     &
                                         start, nLevel, offset,      &
                                         gridLevel)
!
!      ******************************************************************
!      *                                                                *
!      * determineIndirectHalos determines the indirect halo's via a    *
!      * recursive algorithm.                                           *
!      *                                                                *
!      * Step 1.                                                        *
!      * =======                                                        *
!      * Determine for every indirect halo the closest face halo.       *
!      * If several options exist choose the one that does not          *
!      * correspond to a boundary halo, if possible. If several         *
!      * non-boundary halo's exist, just pick one. If all the closest   *
!      * face halo's are boundary halo's then there is no corresponding *
!      * halo and the state is determined by the boundary conditions    *
!      * and/or extrapolation.                                          *
!      * Store the direction from the face halo to the indirect halo.   *
!      *                                                                *
!      * Step 2.                                                        *
!      * =======                                                        *
!      * Determine the level of indirectness of every indirect halo.    *
!      * This is the sum of the absolute values of the elements of the  *
!      * direction vector. For 1st level halo's the maximum level of    *
!      * of indirectness is 2; for 2nd level halo's it is 5. These      *
!      * numbers are for 3 space dimensions.                            *
!      *                                                                *
!      * Step 3.                                                        *
!      * =======                                                        *
!      * Loop over the number of indirect levels.                       *
!      *   For every halo of the current level of indirectness do:      *
!      *     - apply the transformation matrix of its corresponding     *
!      *       face halo to the direction vector.                       *
!      *     - start in the donor cell of the face halo and travel in   *
!      *       the direction of the transformed direction vector.       *
!      *     - you either end up in an internal cell/node or in a halo. *
!      *       Case internal: you're done. Internal cell/node is the    *
!      *                      donor.                                    *
!      *       Case halo: this is guarenteed to be a halo of at least   *
!      *                  one level of indirectness less than the       *
!      *                  current level. Thus the donor is known and    *
!      *                  you're done too. It is possible that it is a  *
!      *                  boundary halo, but this is allowed.           *
!      * End loop over the number of indirect levels.                   *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       use indirectHalo
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nHalo, start, nLevel, offset
       integer(kind=intType), intent(in) :: gridLevel
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType), dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i
!
!      Interfaces
!
       interface
       subroutine closestDirectHalos(entityHalo, entityIndex, &
                                     start, nLevel, offset, gridLevel)
       use block
       use bcHalo
       use haloList
       use indirectHalo
       implicit none
       integer(kind=intType), intent(in) :: start, nLevel, offset
       integer(kind=intType), intent(in) :: gridLevel
       type(haloListType), dimension(:), intent(in) :: entityHalo
       type(indexListType), dimension(:), intent(in) :: entityIndex

       end subroutine closestDirectHalos

       !=================================================================

       subroutine indirectHalosPerLevel(level, iihalo, entityHalo, &
                                        transform, entityIndex)
       use haloList
       use indirectHalo
       use communication
       use BCTypes
       implicit none
       integer(kind=intType), intent(in)    :: level
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType),  dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
       end subroutine indirectHalosPerLevel
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine for every indirect halo the closest direct halo.

       call closestDirectHalos(entityHalo, entityIndex, start, &
                               nLevel, offset, gridLevel)

       ! Sort the indirect halo's.

       call qsortIndHaloType(indHalo, nIndHalo)

       ! Abbreviate the number of indirect levels a bit easier, allocate
       ! the memory for nHaloPerLev and nHaloPerProc, and determine
       ! the values of nHaloPerLev.

       nLevOfInd = indHalo(nIndHalo)%levOfInd
       allocate(nHaloPerLev(0:nlevOfInd), nHaloPerProc(0:nProc), &
                stat=ierr)
       if(ierr /= 0)                                          &
         call terminate("determineIndirectHalos",             &
                        "Allocation error for nHaloPerLev and &
                        &nHaloPerProc")

       nHaloPerLev = 0
       do i=1,nIndHalo
         nHaloPerLev(indHalo(i)%levOfInd) = &
              nHaloPerLev(indHalo(i)%levOfInd) + 1
       enddo

       ! Put nHaloPerLev in cumulative storage format.

       do i=1,nlevOfInd
         nHaloPerLev(i) = nHaloPerLev(i) + nHaloPerLev(i-1)
       enddo

       ! Loop over the number of levels of indirectness.

       do i=1,nLevOfInd

         ! Determine the halo's for this level of indirectness.

         call indirectHalosPerLevel(i, iihalo, entityHalo, transform, &
                                    entityIndex)

         ! Synchronize the processors to avoid possible problems.

         call mpi_barrier(SUmb_comm_world, ierr)
       enddo

       ! Release the memory of the module indirectHalo.

       deallocate(indHalo, nHaloPerLev, nHaloPerProc, stat=ierr)
       if(ierr /= 0) &
         call terminate("determineIndirectHalos",        &
                        "Deallocation error for indHalo, &
                        &nHaloPerLev and nHaloPerProc")

       ! Check in debug mode if iihalo equals nHalo, as it should be.

       if( debug ) then
         if(iihalo /= nHalo)                        &
           call terminate("determineIndirectHalos", &
                          "iihalo differs from nHalo")
       endif

       end subroutine determineIndirectHalos

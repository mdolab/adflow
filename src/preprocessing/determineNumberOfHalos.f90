!
!      ******************************************************************
!      *                                                                *
!      * File:          determineCommPattern.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-20-2003                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineNumberOfHalos(level)
!
!      ******************************************************************
!      *                                                                *
!      * determineNumberOfHalos determines the amount of 1st and 2nd    *
!      * level cell halo's as well as the number of 1st level node      *
!      * halo's stored on this processor.                               *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use haloList
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer(kind=intType) :: i, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the amount of halo cells and nodes to 0.

       nCellHalo1st = 0
       nCellHalo2nd = 0
       nNodeHalo1st = 0

       ! Loop over the number of blocks stored on this processor.

       do i=1,nDom

         ! Set the pointers for this block. The halo construction is
         ! the same for all time spectral solutions, so only the 1st
         ! needs to be considered.

         call setPointers(i,level,1_intType)

         ! Determine the number of 1st level halo cells for this block
         ! and add it to nCellHalo1st. Note the ie == nx + 2, etc.

         nn = nx*ny*nz

         nCellHalo1st = nCellHalo1st - nn + ie*je*ke

         ! Idem for the second level halo's. However there is no variable
         ! which stores nx + 4, so it is computed.

         nCellHalo2nd = nCellHalo2nd - nn + (nx+4)*(ny+4)*(nz+4)

         ! Idem for the 1st level node halo's. Use is made of the fact
         ! that ib == il + 2, etc.

         nNodeHalo1st = nNodeHalo1st + ib*jb*kb - il*jl*kl
       enddo

       end subroutine determineNumberOfHalos

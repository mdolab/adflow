!
!     ******************************************************************
!     *                                                                *
!     * File:          haloIndexCommunication.f90                      *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 07-20-2006                                      *
!     * Last modified: 07-26-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine haloIndexCommunication(level, nHalo)
!
!     ******************************************************************
!     *                                                                *
!     * haloIndexCommunication communicates the global node and cell   *
!     * indices of halo cells between the blocks.                      *
!     * nHalo layers are exchanged.                                    *
!     * These variables are the same for all spectral modes, therefore *
!     * only the 1st mode needs to be communicated.                    * 
!     *                                                                *
!     * This routine was based on                                      *
!     * /utils/whalo.f90 and /preprocessing/xhalo.f90, except that     *
!     * overset grids are not supported.                               *
!     *                                                                *
!     ******************************************************************
!
      use communication
      use flowVarRefState
      use inputPhysics
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level, nHalo
!
!     Local variables.
!

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution                                                *
!     *                                                                *
!     ******************************************************************
!
 
!
!     ******************************************************************
!     *                                                                *
!     *                  The cell halo communication.                  *
!     *                                                                *
!     ******************************************************************
!     
      call indexHalo1(level)
      
      if(nHalo ==2) call indexHalo2(level)

!     
!     ******************************************************************
!     *                                                                *
!     *                  The node halo communication.                  *
!     *                                                                *
!     ******************************************************************


    end subroutine haloIndexCommunication

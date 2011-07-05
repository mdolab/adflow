!
!      ******************************************************************
!      *                                                                *
!      * File:          bl.f90                                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-27-2004                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bl(resOnly)
!
!      ******************************************************************
!      *                                                                *
!      * bl is a sort of dummy routine to be consistent with the other  *
!      * turbulence models. It perfroms the loop over the numer of      *
!      * local blocks and computes the eddy viscosity.                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine argument.
!
       logical, intent(in) :: resOnly
!
!      Local variables.
!
       integer(kind=intType) :: nn, sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if only the residual needs to be computed.
       ! Baldwin-lomax is an algebraic model and therefore it does not
       ! have a residual.

       if( resOnly ) return

       ! Loop over the number of spectral modes and local blocks and
       ! compute the eddy viscosity.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           ! Set the pointers for this block and compute the eddy
           ! viscosity.

           call setPointers(nn, currentLevel, sps)
           call blEddyViscosity

         enddo
       enddo

       end subroutine bl

!      ==================================================================

       subroutine blEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * blEddyViscosity computes the eddy-viscosity according to the   *
!      * Baldwin-Lomax turbulence model for the block given in          *
!      * blockPointers.                                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call terminate("blEddyViscosity", "Not implemented yet")

       end subroutine blEddyViscosity

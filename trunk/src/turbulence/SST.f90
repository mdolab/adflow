!
!      ******************************************************************
!      *                                                                *
!      * File:          SST.f90                                         *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 07-31-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine SST(resOnly)
!
!      ******************************************************************
!      *                                                                *
!      * SST solves the transport equations for menter's SST-variant of *
!      * the k-omega turbulence model in a segregated manner using a    *
!      * diagonal dominant ADI-scheme.                                  *
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
       ! Compute the time derivative for the time spectral mode.

       call unsteadyTurbSpectral(itu1,itu2)

       ! Loop over the number of spectral modes and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, currentLevel, sps)

           ! Set the arrays for the boundary condition treatment.

           call bcTurbTreatment

           ! Solve the transport equations for k and omega.

           call SSTSolve(resOnly)

           ! The eddy viscosity and the boundary conditions are only
           ! applied if an actual update has been computed in SSTSolve.

           if(.not. resOnly ) then

             ! Compute the corresponding eddy viscosity.

             call SSTEddyViscosity

             ! Set the halo values for the turbulent variables.
             ! We are on the finest mesh, so the second layer of halo
             ! cells must be computed as well.

             call applyAllTurbBCThisBlock(.true.)

             ! Write the loglaw for a flat plate boundary layer.

             ! call writeLoglaw

           endif

         enddo domains
       enddo spectralLoop

       end subroutine SST

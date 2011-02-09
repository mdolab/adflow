!
!      ******************************************************************
!      *                                                                *
!      * File:          vf.f90                                          *
!      * Author:        Georgi Kalitzin                                 *
!      * Starting date: 04-14-2004                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine vf(resOnly)
!
!      ******************************************************************
!      *                                                                *
!      * vf solves the transport equations for the v2-f model           *
!      * in a coupled manner using a diagonal dominant ADI-scheme.      *
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
       ! Note that the f-equation, itu4, does not have a time derivative
       ! and therefore should not be included.

       call unsteadyTurbSpectral(itu1,itu3)

       ! Loop over the number of spectral modes and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, currentLevel, sps)

           ! Set the arrays for the boundary condition treatment.

           call bcTurbTreatment

           ! Compute time and length scale

           call vfScale

           ! Solve the transport equations for k and epsilon.

           call keSolve(resOnly)

           ! Solve the transport equation for v2 and the elliptic
           ! equation for f.

           call vfSolve(resOnly)

           ! The eddy viscosity and the boundary conditions are only
           ! applied if actual updates have been computed in keSolve
           ! and vfSolve.

           if(.not. resOnly ) then

             ! Compute the corresponding eddy viscosity.

             call vfEddyViscosity

             ! Set the halo values for the turbulent variables.
             ! We are on the finest mesh, so the second layer of halo
             ! cells must be computed as well.

             call applyAllTurbBCThisBlock(.true.)

             ! Write the loglaw for a flat plate boundary layer.

             ! call writeLoglaw

           endif

         enddo domains
       enddo spectralLoop

       end subroutine vf

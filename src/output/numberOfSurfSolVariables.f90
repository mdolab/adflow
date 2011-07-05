!
!      ******************************************************************
!      *                                                                *
!      * File:          numberOfSurfSolVariables.f90                    *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 05-15-2003                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine numberOfSurfSolVariables(nSolVar)
!
!      ******************************************************************
!      *                                                                *
!      * numberOfSurfSolVariables determines the number of surface      *
!      * variables to be written to the surface solution file.          *
!      *                                                                *
!      ******************************************************************
!
       use precision
       use extraOutput
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: nSolVar
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the number of solution variables to zero.

       nSolVar = 0

       ! Determine the number of surface variables to be written.

       if( surfWriteRho )      nSolVar = nSolVar +1
       if( surfWriteP   )      nSolVar = nSolVar +1
       if( surfWriteTemp )     nSolVar = nSolVar +1
       if( surfWriteVx )       nSolVar = nSolVar +1
       if( surfWriteVy )       nSolVar = nSolVar +1
       if( surfWriteVz )       nSolVar = nSolVar +1
       if( surfWriteCp )       nSolVar = nSolVar +1
       if( surfWritePtotloss ) nSolVar = nSolVar +1
       if( surfWriteMach )     nSolVar = nSolVar +1
       if( surfWriteCf )       nSolVar = nSolVar +1
       if( surfWriteCh )       nSolVar = nSolVar +1
       if( surfWriteYplus )    nSolVar = nSolVar +1
       if( surfWriteCfx )      nSolVar = nSolVar +1
       if( surfWriteCfy )      nSolVar = nSolVar +1
       if( surfWriteCfz )      nSolVar = nSolVar +1
       if( surfWriteBlank )    nSolVar = nSolVar +1

       end subroutine numberOfSurfSolVariables

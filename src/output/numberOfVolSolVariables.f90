!
!      ******************************************************************
!      *                                                                *
!      * File:          numberOfVolVariables.f90                        *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 04-14-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine numberOfVolSolVariables(nVolSolvar, nVolDiscrVar)
!
!      ******************************************************************
!      *                                                                *
!      * numberOfVolSolVariables determines the number of volume        *
!      * variables to be written to the solution file. A distinction is *
!      * made between solution variables and discrete variables. The    *
!      * former discribes the actual solution, the latter is additional *
!      * info such as equation residuals.                               *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       use extraOutput
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: nVolSolvar, nVolDiscrVar
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the number of solution variables to the number of
       ! independent variables and the number of discrete variables to 0.

       nVolSolvar   = nw
       nVolDiscrVar = 0

       ! Check whether or not some additional solution variables must
       ! be written.

       if( volWriteMx )           nVolSolvar = nVolSolvar + 1
       if( volWriteMy )           nVolSolvar = nVolSolvar + 1
       if( volWriteMz )           nVolSolvar = nVolSolvar + 1
       if( volWriteRhoe )         nVolSolvar = nVolSolvar + 1
       if( volWriteTemp )         nVolSolvar = nVolSolvar + 1
       if( volWriteCp )           nVolSolvar = nVolSolvar + 1
       if( volWriteMach )         nVolSolvar = nVolSolvar + 1
       if( volWriteMachTurb )     nVolSolvar = nVolSolvar + 1
       if( volWriteEddyVis )      nVolSolvar = nVolSolvar + 1
       if( volWriteRatioEddyVis ) nVolSolvar = nVolSolvar + 1
       if( volWriteDist )         nVolSolvar = nVolSolvar + 1
       if( volWriteVort )         nVolSolvar = nVolSolvar + 1
       if( volWriteVortx )        nVolSolvar = nVolSolvar + 1
       if( volWriteVorty )        nVolSolvar = nVolSolvar + 1
       if( volWriteVortz )        nVolSolvar = nVolSolvar + 1
       if( volWritePtotloss )     nVolSolvar = nVolSolvar + 1

       ! Check the discrete variables.

       if( volWriteResRho )  nVolDiscrVar  = nVolDiscrVar + 1
       if( volWriteResMom )  nVolDiscrVar  = nVolDiscrVar + 3
       if( volWriteResRhoe ) nVolDiscrVar  = nVolDiscrVar + 1
       if( volWriteResTurb ) nVolDiscrVar  = nVolDiscrVar   &
                                           + (nw - nwf)

       if( volWriteBlank ) nVolDiscrVar  = nVolDiscrVar + 1

       end subroutine numberOfVolSolVariables

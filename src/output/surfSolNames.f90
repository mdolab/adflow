!
!      ******************************************************************
!      *                                                                *
!      * File:          surfSolNames.f90                                *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 05-15-2003                                      *
!      * Last modified: 07-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine surfSolNames(solNames)
!
!      ******************************************************************
!      *                                                                *
!      * surfSolNames sets the names for the surface variables to be    *
!      * written to the surface solution file. Sids convention names    *
!      * are used as much as possible.                                  *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use extraOutput
       implicit none
!
!      Subroutine argument.
!
       character(len=*), dimension(*), intent(out) :: solNames
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize nn to 0.

       nn = 0

       ! Check which surfaces variables must be written and set
       ! solNames accordingly.

       if( surfWriteRho ) then
         nn = nn + 1
         solNames(nn) = cgnsDensity
       endif

       if( surfWriteP   ) then
         nn = nn + 1
         solNames(nn) = cgnsPressure
       endif

       if( surfWriteTemp ) then
         nn = nn + 1
         solNames(nn) = cgnsTemp
       endif

       if( surfWriteVx ) then
         nn = nn + 1
         solNames(nn) = cgnsVelx
       endif

       if( surfWriteVy ) then
         nn = nn + 1
         solNames(nn) = cgnsVely
       endif

       if( surfWriteVz ) then
         nn = nn + 1
         solNames(nn) = cgnsVelz
       endif

       if( surfWriteCp ) then
         nn = nn + 1
         solNames(nn) = cgnsCp
       endif

       if( surfWritePtotloss ) then
         nn = nn + 1
         solNames(nn) = cgnsPtotloss
       endif

       if( surfWriteMach ) then
         nn = nn + 1
         solNames(nn) = cgnsMach
       endif

       if( surfWriteCf ) then
         nn = nn + 1
         solNames(nn) = cgnsSkinFmag
       endif

       if( surfWriteCh ) then
         nn = nn + 1
         solNames(nn) = cgnsStanton
       endif

       if( surfWriteYplus ) then
         nn = nn + 1
         solNames(nn) = cgnsYplus
       endif

       if( surfWriteCfx )  then
         nn = nn + 1
         solNames(nn) = cgnsSkinFx
       endif

       if( surfWriteCfy ) then
         nn = nn + 1
         solNames(nn) = cgnsSkinFy
       endif

       if( surfWriteCfz ) then
         nn = nn + 1
         solNames(nn) = cgnsSkinFz
       endif

       if( surfWriteBlank ) then
         nn = nn + 1
         solNames(nn) = cgnsBlank
       endif

       end subroutine surfSolNames

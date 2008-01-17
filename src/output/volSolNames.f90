!
!      ******************************************************************
!      *                                                                *
!      * File:          volSolNames.f90                                 *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 04-14-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine volSolNames(solNames)
!
!      ******************************************************************
!      *                                                                *
!      * volSolNames sets the names for the volume variables to be      *
!      * written to the volume solution file. Sids convention names are *
!      * used as much as possible.                                      *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use inputPhysics
       use flowVarRefState
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
       ! First store the names of the independent flow variables.

       solNames(1) = cgnsDensity
       solNames(2) = cgnsVelx
       solNames(3) = cgnsVely
       solNames(4) = cgnsVelz
       solNames(5) = cgnsPressure

       ! The turbulent variables if the RANS equations are solved.
       ! Note that these are the primitive variables and not the
       ! conservative ones. The reason is that the sids conventions only
       ! defines these names and not the conservative ones.

       if(equations == RANSEquations) then

         select case(turbModel)

           case(spalartAllmaras, spalartAllmarasEdwards)
             solNames(itu1) = cgnsTurbSaNu

           case(komegaWilcox, komegaModified, menterSST)
             solNames(itu1) = cgnsTurbK
             solNames(itu2) = cgnsTurbOmega

           case(ktau)
             solNames(itu1) = cgnsTurbK
             solNames(itu2) = cgnsTurbTau

           case(v2f)
             solNames(itu1) = cgnsTurbK
             solNames(itu2) = cgnsTurbEpsilon
             solNames(itu3) = cgnsTurbV2
             solNames(itu4) = cgnsTurbF

         end select

       endif

       ! Initialize nn to the number of independent variables.

       nn = nw

       ! Check the additional variables to be written.

       if( volWriteMx )  then
         nn = nn + 1
         solNames(nn) = cgnsMomx
       endif

       if( volWriteMy )  then
         nn = nn + 1
         solNames(nn) = cgnsMomy
       endif

       if( volWriteMz )  then
         nn = nn + 1
         solNames(nn) = cgnsMomz
       endif

       if( volWriteRhoe ) then
         nn = nn + 1
         solNames(nn) = cgnsEnergy
       endif

       if( volWriteTemp ) then
         nn = nn + 1
         solNames(nn) = cgnsTemp
       endif

       if( volWriteCp ) then
         nn = nn + 1
         solNames(nn) = cgnsCp
       endif

       if( volWriteMach ) then
         nn = nn + 1
         solNames(nn) = cgnsMach
       endif

       if( volWriteMachTurb ) then
         nn = nn + 1
         solNames(nn) = cgnsMachTurb
       endif

       if( volWriteEddyVis ) then
         nn = nn + 1
         solNames(nn) = cgnsEddy
       endif

       if( volWriteRatioEddyVis ) then
         nn = nn + 1
         solNames(nn) = cgnsEddyRatio
       endif

       if( volWriteDist ) then
         nn = nn + 1
         solNames(nn) = cgNSWallDist
       endif

       if( volWriteVort ) then
         nn = nn + 1
         solNames(nn) = cgnsVortMagn
       endif

       if( volWriteVortx ) then
         nn = nn + 1
         solNames(nn) = cgnsVortx
       endif

       if( volWriteVorty ) then
         nn = nn + 1
         solNames(nn) = cgnsVorty
       endif

       if( volWriteVortz ) then
         nn = nn + 1
         solNames(nn) = cgnsVortz
       endif

       if( volWritePtotloss ) then
         nn = nn + 1
         solNames(nn) = cgnsPtotloss
       endif

       if( volWriteResRho ) then
         nn = nn + 1
         solNames(nn) = cgnsResRho
       endif

       if( volWriteResMom ) then
         nn = nn + 1
         solNames(nn) = cgnsResMomx

         nn = nn + 1
         solNames(nn) = cgnsResMomy

         nn = nn + 1
         solNames(nn) = cgnsResMomz
       endif

       if( volWriteResRhoe) then
         nn = nn + 1
         solNames(nn) = cgnsResRhoe
       endif

       if( volWriteResTurb ) then

         select case(turbModel)

           case(spalartAllmaras, spalartAllmarasEdwards)
             nn = nn + 1
             solNames(nn) = cgnsResNu

           case(komegaWilcox, komegaModified, menterSST)
             nn = nn + 1
             solNames(nn) = cgnsResK

             nn = nn + 1
             solNames(nn) = cgnsResOmega

           case(ktau)
             nn = nn + 1
             solNames(nn) = cgnsResK

             nn = nn + 1
             solNames(nn) = cgnsResTau

           case(v2f)
             nn = nn + 1
             solNames(nn) = cgnsResK

             nn = nn + 1
             solNames(nn) = cgnsResEpsilon

             nn = nn + 1
             solNames(nn) = cgnsResV2

             nn = nn + 1
             solNames(nn) = cgnsResF

         end select

       endif

       if( volWriteBlank) then
         nn = nn + 1
         solNames(nn) = cgnsBlank
       endif

       end subroutine volSolNames

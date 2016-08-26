!
!       File:          volSolNames.f90                                 
!       Author:        Gaetan Kenway                                   
!       Starting date: 07-21-2013                                      
!       Last modified: 07-21-2013                                      
!
       subroutine isoSurfNames(solNames)
!
!       isoNames sets the names for the volume variables to be         
!       written to the isosurfaces. Sids convention names are          
!       used as much as possible.                                      
!
       use constants
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
!       Begin execution                                                
!

       ! Check the additional variables to be written -- there are no
       ! default variables already written
       nn = 0
       if (isoWriteRho) then
          nn = nn + 1
          solNames(nn) = cgnsDensity
       end if

       if (isoWriteVx) then
          nn = nn + 1
          solNames(nn) = cgnsVelx
       end if

       if (isoWriteVy) then
          nn = nn + 1
          solNames(nn) = cgnsVely
       end if

       if (isoWriteVz) then
          nn = nn + 1
          solNames(nn) = cgnsVelz
       end if

       if (isoWriteP) then
          nn = nn + 1
          solNames(nn) = cgnsPressure
       end if

       if( isoWriteTurb ) then

         select case(turbModel)

           case(spalartAllmaras, spalartAllmarasEdwards)
              nn = nn + 1
             solNames(nn) = cgnsTurbSaNu

           case(komegaWilcox, komegaModified, menterSST)
              nn = nn + 1
              solNames(nn) = cgnsTurbK
              nn = nn + 1
              solNames(nn) = cgnsTurbOmega
              
           case(ktau)
              nn = nn + 1
              solNames(nn) = cgnsTurbK
              nn = nn + 1
              solNames(nn) = cgnsTurbTau

           case(v2f)
              nn = nn + 1
             solNames(nn) = cgnsTurbK
             nn = nn + 1
             solNames(nn) = cgnsTurbEpsilon
             nn = nn + 1
             solNames(nn) = cgnsTurbV2
             nn = nn + 1
             solNames(nn) = cgnsTurbF

         end select

       endif
          
       if( isoWriteMx )  then
         nn = nn + 1
         solNames(nn) = cgnsMomx
       endif

       if( isoWriteMy )  then
         nn = nn + 1
         solNames(nn) = cgnsMomy
       endif

       if( isoWriteMz )  then
         nn = nn + 1
         solNames(nn) = cgnsMomz
       endif

       if( isoWriteRVx )  then
         nn = nn + 1
         solNames(nn) = cgnsRelVelx
       endif

       if( isoWriteRVy )  then
         nn = nn + 1
         solNames(nn) = cgnsRelVely
       endif

       if( isoWriteRVz )  then
         nn = nn + 1
         solNames(nn) = cgnsRelVelz
       endif

       if( isoWriteRhoe ) then
         nn = nn + 1
         solNames(nn) = cgnsEnergy
       endif

       if( isoWriteTemp ) then
         nn = nn + 1
         solNames(nn) = cgnsTemp
       endif

       if( isoWriteCp ) then
         nn = nn + 1
         solNames(nn) = cgnsCp
       endif

       if( isoWriteMach ) then
         nn = nn + 1
         solNames(nn) = cgnsMach
       endif

       if( isoWriteRMach ) then
         nn = nn + 1
         solNames(nn) = cgnsRelMach
       endif

       if( isoWriteMachTurb ) then
         nn = nn + 1
         solNames(nn) = cgnsMachTurb
       endif

       if( isoWriteEddyVis ) then
         nn = nn + 1
         solNames(nn) = cgnsEddy
       endif

       if( isoWriteRatioEddyVis ) then
         nn = nn + 1
         solNames(nn) = cgnsEddyRatio
       endif

       if( isoWriteDist ) then
         nn = nn + 1
         solNames(nn) = cgNSWallDist
       endif

       if( isoWriteVort ) then
         nn = nn + 1
         solNames(nn) = cgnsVortMagn
       endif

       if( isoWriteVortx ) then
         nn = nn + 1
         solNames(nn) = cgnsVortx
       endif

       if( isoWriteVorty ) then
         nn = nn + 1
         solNames(nn) = cgnsVorty
       endif

       if( isoWriteVortz ) then
         nn = nn + 1
         solNames(nn) = cgnsVortz
       endif

       if( isoWritePtotloss ) then
         nn = nn + 1
         solNames(nn) = cgnsPtotloss
       endif

       if( isoWriteResRho ) then
         nn = nn + 1
         solNames(nn) = cgnsResRho
       endif

       if( isoWriteResMom ) then
         nn = nn + 1
         solNames(nn) = cgnsResMomx

         nn = nn + 1
         solNames(nn) = cgnsResMomy

         nn = nn + 1
         solNames(nn) = cgnsResMomz
       endif

       if( isoWriteResRhoe) then
         nn = nn + 1
         solNames(nn) = cgnsResRhoe
       endif

       if( isoWriteResTurb ) then

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

       if (isoWriteShock) then
          nn = nn + 1
          solNames(nn) = cgnsShock
       end if

       if (isoWriteFilteredShock) then
          nn = nn + 1
          solNames(nn) = cgnsFilteredShock
       end if

       if( isoWriteBlank) then
         nn = nn + 1
         solNames(nn) = cgnsBlank
       endif

     end subroutine isoSurfNames
       

!
!      ******************************************************************
!      *                                                                *
!      * File:          writeInputParam.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-07-2004                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeInputParam
!
!      ******************************************************************
!      *                                                                *
!      * writeInputParam writes the input parameters used for the       *
!      * current computation to stdout. This can be used as a reference *
!      * for later purposes by the user.                                *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use flowVarRefState
       use iteration
       use allInputParam
       implicit none
!
!      Local variables.
!
       character(len=20) :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is not processor 0. The parameters
       ! need to be written only once and this is done by processor 0.

       IF(myID /= 0) return
!
!      ******************************************************************
!      *                                                                *
!      * Header such that everything looks nice.                        *
!      *                                                                *
!      ******************************************************************
!
       print "(a)", "#"
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#"
       print "(a)", "# Summary of the relevant input parameters &
                     &used for this computation."
       print "(a)", "#"
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#"
!
!      ******************************************************************
!      *                                                                *
!      * IO parameters.                                                 *
!      *                                                                *
!      ******************************************************************
!
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#   IO parameters"
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#"

       print "(2a)",    "#              Parameter file: ", &
                        trim(paramFile)

       select case (fileFormatRead)
         case (cgnsFormat)
           print "(a)", "#            File format read: CGNS"
           print "(2a)","#                   Grid file: ", &
                        trim(gridFile)

         case (plot3DFormat)
           print "(a)", "#            File format read: PLOT3D"
           print "(2a)","#    PLOT3D Connectivity file: ", &
                         trim(plot3DConnFile)
           print "(2a)","#                   Grid file: ", &
                        trim(gridFile)
       end select

       if( restart ) then
         print "(a)",   "#                     Restart: yes"
         print "(2a)",  "#                Restart file: ", &
                        trim(restartFile)
         if( checkRestartSol) then
           print "(a)", "# Check nonDimensionalization: yes"
         else
           print "(a)", "# Check nonDimensionalization: no"
         endif
       else
         print "(a)",   "#          Restart: no"
       endif
       print "(a)",     "#"

       select case (fileFormatWrite)
         case (cgnsFormat)
           print "(a)", "#           File format write: CGNS"

         case (plot3DFormat)
           print "(a)", "#           File format write: PLOT3D"
       end select

       if(changing_Grid .or. gridMotionSpecified) &
         print "(2a)",  "#               New grid file: ", &
                        trim(newGridFile)
       print "(2a)",    "#               Solution file: ", &
                        trim(solFile)
       print "(2a)",    "#       Surface solution file: ", &
                        trim(surfaceSolFile)
       if( storeRindLayer ) then
         print "(a)",   "#Rind layer in solution files: yes"
       else
         print "(a)",   "#Rind layer in solution files: no"
       endif
       if( writeCoorMeter ) then
         print "(a)",   "#  Write coordinates in meter: yes"
       else
         print "(a)",   "#  Write coordinates in meter: no"
       endif

       if( autoParameterUpdate ) then
         print "(a)",   "#  Automatic parameter update: yes"
       else
         print "(a)",   "#  Automatic parameter update: no"
       endif

       if(cpModel == cpTempCurveFits) &
         print "(2a)",  "#           Cp curve fit file: ", trim(cpFile)

       print "(a)",     "#"
!
!      ******************************************************************
!      *                                                                *
!      * Physics parameters.                                            *
!      *                                                                *
!      ******************************************************************
!
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#   Physics parameters"
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#"

       write(*,"(a)",advance="no") "#                         Equations:"
       select case (equations)
         case (EulerEquations)
           print "(a)", " Euler"
         case (NSEquations)
           print "(a)", " Laminar NS"
         case (RANSEquations)
           print "(a)", " RANS"
       end select

       write(*,"(a)",advance="no") "#                              Mode:"
       select case (equationMode)
         case (steady)
           print "(a)", " Steady"
         case (unsteady)
           print "(a)", " Unsteady"
         case (timeSpectral)
           print "(a)", " Time spectral"
       end select

       write(*,"(a)",advance="no") "#                         Flow type:"
       select case (flowType)
         case (internalFlow)
           print "(a)", " Internal flow"
         case (externalFlow)
           print "(a)", " External flow"
       end select

       write(*,"(a)",advance="no") "#                          Cp model:"
       select case (cpModel)
         case (cpConstant)
           print "(a)", " Constant"
         case (cpTempCurveFits)
           print "(a)", " Temperature curve fits"
       end select

       if(equations == RANSEquations) then
         print "(a)", "#"
         write(*,"(a)",advance="no") "#                  Turbulence &
                                        &model:"
         select case (turbModel)
           case (baldwinLomax)
             print "(a)", " Baldwin Lomax"
           case (spalartAllmaras)
             print "(a)", " Spalart Allmaras"
           case (spalartAllmarasEdwards)
             print "(a)", " Spalart Allmaras Edwards"
           case (komegaWilcox)
             print "(a)", " KOmega Wilcox"
           case (komegaModified)
             print "(a)", " KOmega Modified"
           case (ktau)
             print "(a)", " KTau"
           case (menterSST)
             print "(a)", " Menter SST"
           case (v2f)
             print "(a)", " v2f"
             print "(a,i1)", "#            v2f version (n1 or n6): ", &
                             rvfN
             if( rvfB ) then
               print "(a)", "#              v2f with upper bound: yes"
             else
               print "(a)", "#              v2f with upper bound: no"
             endif
         end select

         select case (turbModel)
           case (spalartAllmaras, spalartAllmarasEdwards, &
                 komegaWilcox,    komegaModified,          &
                 ktau,             v2f)

           write(*,"(a)",advance="no") "#        Turbulence &
                                          &production term:"
           select case (turbProd)
             case (strain)
               print "(a)", " Strain"
             case (vorticity)
               print "(a)", " Vorticity"
             case (katoLaunder)
               print "(a)", " Kato-Launder"
           end select
         end select

         if( wallFunctions ) then
           print "(a)", "#                Use wall functions: yes"
           if(wallOffset > zero) then
             write(string,100) wallOffset
             call removeZerosEformat(string)
             print "(2a)", "# Offset from wall in wall functions: ", &
                           trim(string)
           endif
         else
           print "(a)", "#                Use wall functions: no"
         endif

         write(string,100) eddyVisInfRatio
         call removeZerosEformat(string)
         print "(2a)", "#  Free stream eddy viscosity ratio: ", &
                       trim(string)

         select case (turbModel)
           case (komegaWilcox, komegaModified, ktau, v2f)
             write(string,100) turbIntensityInf
             call removeZerosEformat(string)
             print "(2a)", "#   Free stream turbulent &
                           &intensity: ", trim(string)
         end select

       endif

       print "(a)", "#"
       if(cpModel == cpConstant) then
         write(string,100) gammaConstant
         call removeZerosEformat(string)
         print "(2a)", "#       Constant specific heat ratio: ", &
                       trim(string)
       endif
       write(string,100) RGasDim
       call removeZerosEformat(string)
       print "(2a)", "#            Gas constant (J/(kg k)): ", &
                     trim(string)

       select case (equations)
         case (NSEquations, RANSEquations)
           write(string,100) prandtl
           call removeZerosEformat(string)
           print "(2a)", "#                     Prandtl number: ", &
                         trim(string)
       end select

       if(equations == RANSEquations) then
         write(string,100) prandtlTurb
         call removeZerosEformat(string)
         print "(2a)", "#           Turbulent Prandtl number: ", &
                       trim(string)

         if( kPresent ) then
           write(string,100) pklim
           call removeZerosEformat(string)
           print "(2a)", "#              Max ratio k-prod/dest: ", &
                         trim(string)
         endif
       endif

       if(flowType == externalFlow) then

         write(string,100) Mach
         call removeZerosEformat(string)
         print "(2a)", "#                               Mach: ", &
                       trim(string)

         write(*,"(a)",advance="no") "#     Free stream velocity &
                                     &direction:"
         write(string,100) velDirFreestream(1)
         call removeZerosEformat(string)
         write(*,"(2a)",advance="no") " ", trim(string)
         write(string,100) velDirFreestream(2)
         call removeZerosEformat(string)
         write(*,"(2a)",advance="no") " ", trim(string)
         write(string,100) velDirFreestream(3)
         call removeZerosEformat(string)
         print "(2a)", " ", trim(string)

         write(*,"(a)",advance="no") "#                     Lift &
                                     &direction:"
         write(string,100) liftDirection(1)
         call removeZerosEformat(string)
         write(*,"(2a)",advance="no") " ", trim(string)
         write(string,100) liftDirection(2)
         call removeZerosEformat(string)
         write(*,"(2a)",advance="no") " ", trim(string)
         write(string,100) liftDirection(3)
         call removeZerosEformat(string)
         print "(2a)", " ", trim(string)

         select case (equations)
           case (NSEquations, RANSEquations)

             write(string,100) Reynolds
             call removeZerosEformat(string)
             print "(2a)", "#                           Reynolds: ", &
                           trim(string)

             write(string,100) ReynoldsLength
             call removeZerosEformat(string)
             print "(2a)", "#         Reynolds length (in meter): ", &
                           trim(string)

             write(string,100) tempFreestream
             call removeZerosEformat(string)
             print "(2a)", "#     Free stream temperature (in K): ", &
                           trim(string)

         end select

       endif

       write(string,100) MachCoef
       call removeZerosEformat(string)
       print "(2a)", "#              Mach for coefficients: ", &
                     trim(string)

       write(string,100) surfaceRef
       call removeZerosEformat(string)
       print "(2a)", "#                  Reference surface: ", &
                     trim(string)

       write(string,100) lengthRef
       call removeZerosEformat(string)
       print "(2a)", "#                   Reference length: ", &
                     trim(string)

       write(string,100) pointRef(1)
       call removeZerosEformat(string)
       print "(2a)", "#           Moment reference point x: ", &
                     trim(string)
       write(string,100) pointRef(2)
       call removeZerosEformat(string)
       print "(2a)", "#           Moment reference point y: ", &
                     trim(string)
       write(string,100) pointRef(3)
       call removeZerosEformat(string)
       print "(2a)", "#           Moment reference point z: ", &
                     trim(string)

       print "(a)",  "#"
!
!      ******************************************************************
!      *                                                                *
!      * Reference state.                                               *
!      *                                                                *
!      ******************************************************************
!
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#   Reference State"
       print "(a)", "#---------------------------------------------&
                    &---------------------------------"
       print "(a)", "#"





       ! Format statements.

 100   format(e12.5)

       end subroutine writeInputParam

!      ==================================================================

       subroutine removeZerosEformat(string)
!
!      ******************************************************************
!      *                                                                *
!      * removeZerosEformat removes the non-meaning zero's from the     *
!      * given string, which contains a real number in scientific       *
!      * format, i.e. written with the e edit descriptor.               *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(inout) :: string
!
!      Local variables.
!
       integer :: i, pos, pos1, lenString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Get rid of the leading spaces, convert string to lower case
       ! and determine its trimmed length.

       string = adjustl(string)
       call convertToLowerCase(string)
       lenString = len_trim(string)

       ! Determine the position of "e" in the string and check if it
       ! was found.

       pos = index(string,"e")
       if(pos > 0) then

         ! Determine the 1st nonzero digit to the left of pos.
         ! Store this position in pos1.

         pos1 = pos - 1
         do
           if(string(pos1:pos1) /= "0") exit
           if(pos1 <= 3)                exit
           pos1 = pos1 - 1
         enddo

         ! Add 1 to pos1, such that it contains the position where to
         ! store string(pos1:). Remove the zero's by overwriting them.

         pos1 = pos1 + 1
         do i=pos,lenString
           string(pos1:pos1) = string(i:i)
           pos1 = pos1 + 1
         enddo

         ! Overwrite the rest with blanks.

         do i=pos1,lenString
           string(i:i) = " "
         enddo

       endif

       end subroutine removeZerosEformat

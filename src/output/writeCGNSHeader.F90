!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSHeader.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-11-2003                                      *
!      * Last modified: 07-03-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSHeader(cgnsInd, base)
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSHeader writes a descriptive header to the given base  *
!      * of the given CGNS file. Only processor 0 performs this task.   *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use flowVarRefState
       use su_cgns
       use inputPhysics
       use inputTimeSpectral
       use monitor
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in) :: cgnsInd, base
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS

       real(kind=cgnsRealType) :: val

       character(len=2048) :: message
       character(len=7)    :: integerString
       character(len=12)   :: realString
!
!      Function definition.
!
       integer :: setCGNSRealType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("writeCGNSHeader", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else

       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Go to the correct position in the CGNS file.

       call cg_goto_f(cgnsInd, base, ierr, "end")
       if(ierr /= all_ok)                    &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling cg_goto_f")

       ! Create a data class type node to indicate that nonDimensional
       ! solution data is written for which the reference state
       ! is known.

       call cg_dataclass_write_f(NormalizedByDimensional,ierr)
       if(ierr /= all_ok)                  &
         call terminate("writeCGNSHeader", &
                        "Something wrong when calling &
                        &cg_dataclass_write_f")

       ! Write the info about the solver used.

       call cg_descriptor_write_f("SolverInfo", &
                                  "SUmb multiblock code", ierr)
       if(ierr /= all_ok)                  &
         call terminate("writeCGNSHeader", &
                        "Something wrong when calling &
                        &cg_descriptor_write_f")

       ! Write the info about the scheme used; message is used as
       ! storage for the string containing the scheme description.

       call describeScheme(message)
       call cg_descriptor_write_f("DiscretizationScheme", message, ierr)
       if(ierr /= all_ok)                  &
         call terminate("writeCGNSHeader", &
                        "Something wrong when calling &
                        &cg_descriptor_write_f")

       ! Write the similation type to the CGNS file.

       select case (equationMode)

         case (steady)

           ! Steady mode. Just write this info.

           call cg_simulation_type_write_f(cgnsInd, base, &
                                           nonTimeaccurate, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_simulation_type_write_f")

         !===============================================================

         case (unsteady)

           ! Unsteady mode. First write the simulation type.
 
           call cg_simulation_type_write_f(cgnsInd, base, &
                                           timeaccurate, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_simulation_type_write_f")

           ! Write some additional stuff, like time step and
           ! physical time. First store it in the big string message.

           write(integerString,"(i7)") timeStepUnsteady + &
                                       nTimeStepsRestart
           write(realString,"(e12.5)") timeUnsteady + &
                                       timeUnsteadyRestart

           integerString = adjustl(integerString)
           realString    = adjustl(realString)

           write(message,101) trim(integerString), trim(realString)
 101       format("Unsteady time step ",a,", physical time ",a, &
                  " seconds")

           ! And write the info.

           call cg_descriptor_write_f("UnsteadyInfo", message, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_descriptor_write_f")

         !===============================================================

         case (timeSpectral)

           ! Time spectral mode. This is not a predefined mode in CGNS
           ! and therefore use userDefined.

           call cg_simulation_type_write_f(cgnsInd, base, &
                                           UserDefined, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_simulation_type_write_f")

           ! Write some info to the string message.

           write(integerString,"(i7)") nTimeIntervalsSpectral
           integerString = adjustl(integerString)

           write(message,102) trim(integerString)
 102       format("Time spectral mode for periodic problems; ",a,  &
                  " spectral solutions have been used to model the &
                  &problem.")

           ! And write the info.

           call cg_descriptor_write_f("PeriodicInfo", message, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_descriptor_write_f")
       end select

       ! Go back to the given base in the cgns file.

       call cg_goto_f(cgnsInd, base, ierr, "end")
       if(ierr /= all_ok)                    &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling cg_goto_f")

       ! Create a flow equation set.

       call cg_equationset_write_f(cgnsPhysDim, ierr)
       if(ierr /= all_ok)                    &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_equationset_write_f")

       ! Write the rest of the physical model under the flow
       ! equation set just created.

       call cg_goto_f(cgnsInd, base, ierr, &
                      "FlowEquationSet_t", 1, "end")
       if(ierr /= all_ok)                  &
         call terminate("writeCGNSHeader", &
                        "Something wrong when calling cg_goto_f")

       ! Write the governing equations solved.

       select case (equations)
         case (EulerEquations)
           call cg_governing_write_f(Euler, ierr)

         case (NSEquations)
           call cg_governing_write_f(nsLaminar, ierr)

         case (RANSEquations)
           call cg_governing_write_f(nsTurbulent, ierr)
       end select

       if(ierr /= all_ok)                  &
         call terminate("writeCGNSHeader", &
                        "Something wrong when calling &
                        &cg_governing_write_f")

       ! Write the information about the gas model used.
       ! Determine the cp model used in the computation.

       select case (cpModel)

         case (cpConstant)

           ! Constant cp and thus constant gamma.

           call cg_model_write_f("GasModel_t", Ideal, ierr)
           if(ierr /= all_ok)                &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_model_write_f")

           ! Write the actual value of gamma; this must be done under
           ! gas model type, which explains the goto statement.

           call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", &
                          1, "GasModel_t", 1, "end")
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling cg_goto_f")

           val = gammaConstant
           call cg_array_write_f(cgnsHeatRatio, realTypeCGNS, &
                                   1, 1, val, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_array_write_f")

           ! And create a data class under SpecificHeatRatio to tell that
           ! this is a nonDimensional parameter.

           call cg_goto_f(cgnsInd, base, ierr,       &
                          "FlowEquationSet_t", 1,    &
                          "GasModel_t", 1, "DataArray_t", 1,"end")
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling cg_goto_f")

           call cg_dataclass_write_f(NonDimensionalParameter,ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_dataclass_write_f")

         !===============================================================

         case (cpTempCurveFits)

           ! Cp as function of the temperature is given via curve fits.

           call cg_model_write_f("GasModel_t", ThermallyPerfect, ierr)
           if(ierr /= all_ok)                  &
             call terminate("writeCGNSHeader", &
                            "Something wrong when calling &
                            &cg_model_write_f")

       end select

       ! The rest of physical model description is only
       ! for viscous flows.

       viscousTest: if( viscous ) then

         ! Write the info of the viscosity model. Under the flow
         ! equation set.

         call cg_goto_f(cgnsInd, base, ierr, &
                        "FlowEquationSet_t", 1, "end")
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling cg_goto_f")

         call cg_model_write_f("ViscosityModel_t", sutherlandlaw, ierr)
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_model_write_f")

         ! Write the info about the thermal conductivity, i.e.
         ! Constant Prandtl number. Write the used value as well.

         call cg_model_write_f("ThermalConductivityModel_t", &
                               constantPrandtl, ierr)
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_model_write_f")

         call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", 1,&
                        "ThermalConductivityModel_t", 1, "end")
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling cg_goto_f")

         val = prandtl
         call cg_array_write_f(cgnsPrandtl, realTypeCGNS, 1, 1, &
                               val, ierr)
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_array_write_f")

         ! And create a data class under Prandtl number to tell that
         ! this is a nonDimensional parameter.

         call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", &
                        1, "ThermalConductivityModel_t", 1,       &
                        "DataArray_t", 1,"end")
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(NonDimensionalParameter,ierr)
         if(ierr /= all_ok)                  &
           call terminate("writeCGNSHeader", &
                          "Something wrong when calling &
                          &cg_dataclass_write_f")

         ! The rest of the physical model description is only for the
         ! RANS equations.

         turbulentTest: if(equations == RANSEquations) then

           select case(turbModel)

             case (baldwinLomax)
               call writeCGNSBaldwinLomaxInfo(cgnsInd, base)

             case (spalartAllmaras)
               call writeCGNSSaInfo(cgnsInd, base)

             case (spalartAllmarasEdwards)
               call writeCGNSSaeInfo(cgnsInd, base)

             case (komegaWilcox)
               call writeCGNSKomegaWilcoxInfo(cgnsInd, base)

             case (komegaModified)
               call writeCGNSKomegaModifiedInfo(cgnsInd, base)

             case (ktau)
               call writeCGNSKtauInfo(cgnsInd, base)

             case (menterSST)
               call writeCGNSMenterSSTInfo(cgnsInd, base)

             case (v2f)
               call writeCGNSV2fInfo(cgnsInd, base)
           end select

         endif turbulentTest

       endif viscousTest

       ! Write the reference state.

       call writeCGNSReferenceState(cgnsInd, base)

#endif

       end subroutine writeCGNSHeader


!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSReferenceState.F90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-22-2003                                      *
!      * Last modified: 06-29-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSReferenceState(cgnsInd, cgnsBase)
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSReferenceState writes the reference state to the      *
!      * cgns file. Enough info is specified such that a restart can be *
!      * performed by a different solver, which uses a different        *
!      * nonDimensionalization.                                         *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use su_cgns
       use inputPhysics
       use flowVarRefState
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsInd, cgnsBase
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS, ii

       integer(kind=intType) :: i

       real(kind=cgnsRealType) :: val
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

       call terminate("writeReferenceState", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Go to the base.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "end")
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling cg_goto_f")

       ! Create the reference state node with a nice description.

       call cg_state_write_f("Reference state variables for &
                             &nonDimensional data. Variables are &
                             &nonDimensionalized using the reference &
                             &density, pressure and temperature.", ierr)
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling cg_state_write_f")

       ! The actual data should be written below the reference state
       ! node. So go there first.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                      "ReferenceState_t", 1, "end")
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling cg_goto_f")

       ! Write the Mach number and indicate that it is a nonDimensional
       ! parameter

       val = Mach
       call cg_array_write_f(cgnsMach, realTypeCGNS, 1, 1, val, ierr)
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling cg_array_write_f")

       ii = 1
       call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
                      "DataArray_t", ii, "end")
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling cg_goto_f")

       call cg_dataclass_write_f(NonDimensionalParameter, ierr)
       if(ierr /= all_ok)                      &
         call terminate("writeReferenceState", &
                        "Something wrong when calling &
                        &cg_dataclass_write_f")

       ! Write the 3 flow angles. The units are degrees.

       velocityDir: do i=1,3

         ! Go to the reference state node.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                        "ReferenceState_t", 1, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         ! Store component i of the direction in val.

         val = velDirFreestream(i)

         select case(i)
           case (1_intType)
             call cg_array_write_f(cgnsVelVecX, realTypeCGNS, &
                                   1, 1, val, ierr)

           case (2_intType)
             call cg_array_write_f(cgnsVelVecY, realTypeCGNS, &
                                   1, 1, val, ierr)

           case (3_intType)
             call cg_array_write_f(cgnsVelVecZ, realTypeCGNS, &
                                   1, 1, val, ierr)
         end select

         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_array_write_f")

         ! Write the info that the unit vector is nondimensional.

         ii = ii + 1
         call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
                        "DataArray_t", ii, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(NonDimensionalParameter, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_dataclass_write_f")

       enddo velocityDir

       ! Write Reynolds number, etc. In case of an external viscous flow.

       testViscous: if(viscous .and. flowType == externalFlow) then

         ! Go to the reference state node and write the Reynolds number.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                        "ReferenceState_t", 1, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         val = Reynolds
         call cg_array_write_f(cgnsReyn, realTypeCGNS, 1, 1, &
                               val, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_array_write_f")

         ii = ii + 1
         call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
                        "DataArray_t", ii, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(NonDimensionalParameter, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_dataclass_write_f")

         ! Go to the reference state and write the length on which the
         ! Reynolds number is based.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                        "ReferenceState_t", 1, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         val = ReynoldsLength
         call cg_array_write_f(cgnsReynLen, realTypeCGNS, &
                               1, 1, val, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_array_write_f")

         ii = ii + 1
         call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
                        "DataArray_t", ii, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(NormalizedByDimensional, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_dataclass_write_f")

       endif testViscous

       ! Write some reference values of the density, pressure, temperature,
       ! velocity and length.

       refLoop: do i=1,5

         ! Go to the reference state node.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                        "ReferenceState_t", 1, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         ! Write a value, depending on i.

         select case(i)
           case (1_intType)
             val = rhoref
             call cg_array_write_f(cgnsDensity, realTypeCGNS, &
                                   1, 1, val, ierr)
           case (2_intType)
             val = pref
             call cg_array_write_f(cgnsPressure, realTypeCGNS, &
                                   1, 1, val, ierr)

           case (3_intType)
             val = Tref
             call cg_array_write_f(cgnsTemp, realTypeCGNS, &
                                   1, 1, val, ierr)

           case (4_intType)
             val = sqrt(pref/rhoref)
             call cg_array_write_f(cgnsVelocity, realTypeCGNS, &
                                   1, 1, val, ierr)

           case (5_intType)
             val = one
             call cg_array_write_f(cgnsLength, realTypeCGNS, &
                                   1, 1, val, ierr)
         end select

         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_array_write_f")

         ! Write the info that the this reference value is dimensional
         ! and based on si units.

         ii = ii + 1
         call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
                        "DataArray_t", ii, "end")
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(Dimensional, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_dataclass_write_f")

         call cg_units_write_f(Kilogram, Meter, Second, Kelvin, &
                               Null, ierr)
         if(ierr /= all_ok)                      &
           call terminate("writeReferenceState", &
                          "Something wrong when calling &
                          &cg_units_write_f")

       enddo refLoop

#endif

       end subroutine writeCGNSReferenceState

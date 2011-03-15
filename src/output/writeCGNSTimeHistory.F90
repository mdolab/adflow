!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCgnsTimeHistory.F90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-20-2004                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSTimeHistory
!
!      ******************************************************************
!      *                                                                *
!      * WriteCGNSTimeHistory writes for unsteady computations          *
!      * the time history of the monitoring variables to the            *
!      * cgns file.                                                     *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("writeCGNSTimeHistory", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use cgnsNames
       use inputIO
       use monitor
       use su_cgns
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: i, nn, mm, cgnsInd, base, ierr, realTypeCGNS

       real(kind=4), dimension(:), allocatable :: buf4
       real(kind=8), dimension(:), allocatable :: buf8

       character(len=maxCGNSNameLen) :: cgnsName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the file index and base a bit easier. Note that this info
       ! is only written to the 1st file, because this is an unsteady
       ! computation.

       cgnsInd = fileIDs(1)
       base    = cgnsBases(1)

       ! Store the total number of time steps in nn.

       nn = timeStepUnsteady + nTimeStepsRestart

       ! Depending on the input option, set the CGNS type and allocate
       ! the memory for either buf4 or buf8.

       ! Set the cgns real type depending on the input option.

       select case (precisionSol)
         case (precisionSingle)
           realTypeCGNS = RealSingle
           allocate(buf4(nn), stat=ierr)

         !===============================================================

         case (precisionDouble)
           realTypeCGNS = RealDouble
           allocate(buf8(nn), stat=ierr)
       end select

       if(ierr /= 0)                            &
         call terminate("writeCGNSTimeHistory", &
                        "Memory allocation failure for either buf4 &
                        &or buf8")

       ! Go to the correct position in the cgns file.

       call cg_goto_f(cgnsInd, base, ierr, "end")
       if(ierr /= all_ok)                       &
         call terminate("writeCGNSTimeHistory", &
                        "Something wrong when calling cg_goto_f")

       ! Create the name of the base iterative data node.

       cgnsName = "TimeHistory"

       ! Create the base iterative node and check if everything
       ! went okay.

       call cg_biter_write_f(cgnsInd, base, cgnsName, nn, ierr)
       if(ierr /= all_ok)                       &
         call terminate("writeCGNSTimeHistory", &
                        "Something wrong when calling cg_biter_write_f")

       ! The time history must be written under the node just created.
       ! Go there.

       call cg_goto_f(cgnsInd, base, ierr, &
                      "BaseIterativeData_t" , 1, "end")
       if(ierr /= all_ok)                          &
         call terminate("writeCGNSTimeHistory", &
                        "Something wrong when calling cg_goto_f")

       ! Write the time values.

       cgnsName = cgnsTimeValue
       call cg_array_write_f(cgnsName, realTypeCGNS, 1, nn, &
                             timeArray, ierr)
       if(ierr /= all_ok)                       &
         call terminate("writeCGNSTimeHistory", &
                        "Something wrong when calling cg_array_write_f")

       ! Loop over the number of monitoring variables and write
       ! their time history.

       monLoop: do i=1,nMon

         ! Copy the time history to either buf4 or buf8 and write it
         ! to file.

         select case (precisionSol)
           case (precisionSingle)
             do mm=1,nn
               buf4(mm) = timeDataArray(mm,i)
             enddo

             call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
                                   buf4, ierr)

           !=============================================================

           case (precisionDouble)
             do mm=1,nn
               buf8(mm) = timeDataArray(mm,i)
             enddo

             call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
                                   buf8, ierr)
         end select

         if(ierr /= all_ok)                       &
           call terminate("writeCGNSTimeHistory", &
                          "Something wrong when calling &
                          &cg_array_write_f")
       enddo monLoop

       ! Release the memory of buf4 or buf8.

       select case (precisionSol)
         case (precisionSingle)
           deallocate(buf4, stat=ierr)

         case (precisionDouble)
           deallocate(buf8, stat=ierr)
       end select

       if(ierr /= 0)                            &
         call terminate("writeCGNSTimeHistory", &
                        "Deallocation failure for either buf4 or buf8") 

#endif

       end subroutine writeCGNSTimeHistory

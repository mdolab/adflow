!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSConvInfo.F90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-24-2003                                      *
!      * Last modified: 10-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSConvInfo
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSConvInfo writes the convergence info to the           *
!      * cgns file(s).                                                  *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("writeCGNSConvInfo", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use inputIO
       use inputPhysics
       use monitor
       use su_cgns
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: nConv, base, cgnsInd, conv
       integer :: i, nn, mm, ierr, realTypeCGNS

       real(kind=4), dimension(:), allocatable :: buf4
       real(kind=8), dimension(:), allocatable :: buf8
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if the convergence history (of the inner
       ! iterations) does not need to be stored. This logical can
       ! only be .false. for an unsteady computation.

       if(.not. storeConvInnerIter) return

       ! Store the number of iterations to be written in nn.
       ! This is nIterCur + 1, because the array starts at 0.

       nn = nIterCur + 1

       ! Depending on the input option, set the CGNS type and allocate
       ! the memory for either buf4 or buf8.

       ! Set the cgns real type depending on the input option.

       select case (precisionSol)
         case (precisionSingle)
           realTypeCGNS = RealSingle
           allocate(buf4(0:nIterCur), stat=ierr)

         !===============================================================

         case (precisionDouble)
           realTypeCGNS = RealDouble
           allocate(buf8(0:nIterCur), stat=ierr)
       end select

       if(ierr /= 0)                         &
         call terminate("writeCGNSConvInfo", &
                        "Memory allocation failure for either buf4 &
                        &or buf8")

       ! Determine the number of convergence histories to be written.
       ! This depends on the equation mode.

       select case (equationMode)

         case (steady,unsteady)
           nConv = 1

         case (timeSpectral)
           nConv = nVolSolToWrite    ! == number of spectral solutions.

       end select

       ! Loop over the number of convergence histories to be written.

       convLoop: do conv=1,nConv

         ! Abbreviate the corresponding file and base a bit easier.

         cgnsInd = fileIDs(conv)
         base    = cgnsBases(conv)

         ! Go to the correct position in the CGNS file.

         call cg_goto_f(cgnsInd, base, ierr, "end")
         if(ierr /= all_ok)                    &
           call terminate("writeCGNSConvInfo", &
                          "Something wrong when calling cg_goto_f")

         ! Create the convergence history node. Add a small description.

         call cg_convergence_write_f(nn,"L2 norms are computed by taking &
                                        &the square root of the quotient &
                                        &the sum of the square of the &
                                        &residuals and the total number &
                                        &of cells in the grid.", ierr)
         if(ierr /= all_ok)                    &
           call terminate("writeCGNSConvInfo", &
                          "Something wrong when calling &
                          &cg_convergence_write_f")

         ! The convergence history must be written under the node just
         ! created. Go there.

         call cg_goto_f(cgnsInd, base, ierr, &
                        "ConvergenceHistory_t", 1, "end")
         if(ierr /= all_ok)                    &
           call terminate("writeCGNSConvInfo", &
                          "Something wrong when calling cg_goto_f")

         ! Loop over the number of monitoring variables.

         monLoop: do i=1,nMon

           ! Copy the convergence info to either buf4 or buf8 and write
           ! it to file.

           select case (precisionSol)
             case (precisionSingle)
               do mm=0,nIterCur
                 buf4(mm) = convArray(mm,conv,i)
               enddo

               call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
                                     buf4, ierr)

             !===========================================================

             case (precisionDouble)
               do mm=0,nIterCur
                 buf8(mm) = convArray(mm,conv,i)
               enddo

               call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
                                     buf8, ierr)
           end select

           if(ierr /= all_ok)                    &
             call terminate("writeCGNSConvInfo", &
                            "Something wrong when calling &
                            &cg_array_write_f")
         enddo monLoop
       enddo convLoop

       ! Release the memory of buf4 or buf8.

       select case (precisionSol)
         case (precisionSingle)
           deallocate(buf4, stat=ierr)

         case (precisionDouble)
           deallocate(buf8, stat=ierr)
       end select

       if(ierr /= 0)                         &
         call terminate("writeCGNSConvInfo", &
                        "Deallocation failure for either buf4 or buf8")

#endif

       end subroutine writeCGNSConvInfo

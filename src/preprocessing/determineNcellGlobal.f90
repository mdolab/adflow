!
!      ******************************************************************
!      *                                                                *
!      * File:          determineNcellGlobal.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-28-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineNcellGlobal(level)
!      ******************************************************************
!      *                                                                *
!      * determineNcellGlobal determines the global number of cells     *
!      * the given grid level. This info is needed to compute the L2    *
!      * norm of the residuals in the flow solver.                      *
!      * Only the 1st spectral solution needs to be considered, because *
!      * this info is identical for all of them.                        *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn, nCellLocal

       character(len=12) :: int1String, int2String
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the local number of cells by looping over the blocks.

       nCellLocal = 0
       do nn=1,nDom
         nCellLocal = nCellLocal + flowDoms(nn,level,1)%nx &
                    *              flowDoms(nn,level,1)%ny &
                    *              flowDoms(nn,level,1)%nz
       enddo

       ! And determine the global sum.

       call mpi_allreduce(nCellLocal, nCellGlobal(level), 1, &
                          sumb_integer, mpi_sum, SUmb_comm_world, ierr)

       ! Write the total number of cells to stdout; only done by
       ! processor 0 to avoid a messy output.

       if(myID == 0) then

         write(int1String,"(i12)") level
         write(int2String,"(i12)") nCellGlobal(level)
         int1String = adjustl(int1String)
         int2String = adjustl(int2String)

         print "(a)", "#"
         print 101, trim(int1String), trim(int2String)
         print "(a)", "#"
 101     format("# Grid level: ", a,", Total number of cells: ", a)

       endif

       end subroutine determineNcellGlobal

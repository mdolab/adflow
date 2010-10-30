!
!     ******************************************************************
!     *                                                                *
!     * File:          getADjoint.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-18-2008                                      *
!     * Last modified: 10-04-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getADjoint(nnodes,functionGradLocal)
      use ADjointPETSc
      use ADjointVars
      use blockpointers !globalnode

      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType),intent(in):: nnodes
      real(kind=realType),dimension(nNodes),intent(out) :: functionGradLocal

!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

!      integer(kind=intType) :: idx!tmp for fd

      character(len=2*maxStringLen) :: errorMessage

      integer(kind=intType) :: idxmg, iLow, iHigh, n 		

#ifndef USE_NO_PETSC		
    
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

!     ******************************************************************
!     *                                                                *
!     * Transfer solution from PETSc context.                          *
!     *                                                                *
!     ******************************************************************
!

      ! Query about the ownership range.
      ! iHigh is one more than the last element stored locally.

      call VecGetOwnershipRange(psi, iLow, iHigh, PETScIerr)
	
      !print *,'irange',ilow,ihigh,petscrank
!      print *,'irange',iLow,iHigh

      if( PETScIerr/=0 ) &
        call terminate("getADjoint", &
                       "Error in VecGetOwnershipRange psi")

	
      ! VecGetValues - Gets values from certain locations of a vector.
      !           Currently can only get values on the same processor.

      n = 0

      do idxmg=iLow, iHigh-1

        n = n + 1
      
        !print *,'getting value', idxmg,n
        call VecGetValues(psi, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)
	

        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("getADjoint", errorMessage)
        endif

     enddo

     !print *,'flushing'
      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

    end subroutine getADjoint

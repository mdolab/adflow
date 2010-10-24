!
!     ******************************************************************
!     *                                                                *
!     * File:          setADjoint.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 10-04-2010                                      *
!     * Last modified: 10-04-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setADjoint(ncells,functionGradLocal)
      use ADjointPETSc
      use ADjointVars

      use blockpointers !globalnode

      implicit none
!
!     Subroutine arguments.
! 
      integer(kind=intType),intent(in):: ncells
      real(kind=realType),dimension(ncells),intent(in) :: functionGradLocal
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

      ! Send some feedback to screen.

    !   if( PETScRank==0 ) &
!         write(*,10) "Setting ADjoint Vector..."	
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
	
      if( PETScIerr/=0 ) &
        call terminate("setADjoint", &
                       "Error in VecGetOwnershipRange psi")
      ! VecGetValues - Gets values from certain locations of a vector.
      !           Currently can only get values on the same processor.

      n = 0

      do idxmg=iLow, iHigh-1

        n = n + 1
       

        call VecSetValue(psi, idxmg, &
                          functionGradLocal(n),INSERT_VALUES, PETScIerr)
	
        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("getADjoint", errorMessage)
        endif

     enddo


     !assemble vector

     call VecAssemblyBegin(psi, PETScIerr)

     if( PETScIerr/=0 ) &
          call terminate("setADjoint", "Error in VecAssemblyBegin") 

     call VecAssemblyEnd(psi,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setADjoint", "Error in VecAssemblyEnd")
      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif
    end subroutine setADjoint

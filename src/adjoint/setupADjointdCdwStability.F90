!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointdCdwStability.F90                   *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-27-2009                                      *
!     * Last modified: 11-27-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupADjointdCdwStability(level,costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * This routine computes the dCdw portion of the right hand side  *
!     * of the discrete ADjoint problem when using the Time Spectral   *
!     * stability derivative formulation.                              *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use BCTypes          ! i/j/k/max/min
      use blockPointers    ! il,jl,kl,nViscBocos,nInvBocos
                           ! w,x,BCFaceID,groupNum,BCData,d2Wall,muLam
                           ! nBocos,viscousSubface
                           ! globalNode
      use inputTimeSpectral   ! nTimeIntervalsSpectral
      use cgnsGrid         ! cgnsFamilies
      use flowVarRefState  ! nw,LRef, pInf, gammaInf
      use inputPhysics     ! pointRef, equations, RANSEquations, 
                           ! MachCoef,surfaceRef
      use monitor          ! timeUnsteadyRestart
      use section          ! nSections,section%
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level, costFunction!, sps
!
!     Local variables.
!
      real(kind=realType), dimension(3) :: cFp, cFv
      real(kind=realType), dimension(3) :: cMp, cMv

      real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
      real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

      integer(kind=intType) :: nn, mm, i, j,liftIndex,nnn
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
      integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End

      integer(kind=intType) :: iCell, jCell, kCell

      real(kind=realType), dimension(nTimeIntervalsSpectral) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,&
                             &CmxAdj,CmyAdj,CmzAdj 

      real(kind=realType), dimension(nTimeIntervalsSpectral) :: ClAdjB,CdAdjB,CfxAdjB,CfyAdjB,CfzAdjB,&
                             &CmxAdjB,CmyAdjB,CmzAdjB  

      real(kind=realType) :: yplusMax, test

      real(kind=realType), dimension(3) :: refPoint

      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,xAdjB
      real(kind=realType), dimension(:,:,:,:), allocatable :: wAdj,wAdjB
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj

      REAL(KIND=REALTYPE) :: machadj, machcoefadj, uinfadj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb,machgridadj, machgridadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb
      REAL(KIND=REALTYPE) :: rotcenteradj(3), rotrateadj(3), rotrateadjb(3)

      logical :: contributeToForce, viscousSubface,secondHalo,righthanded

      real(kind=realType), dimension(nSections) :: t

      ! dJ/dw row block
      
      real(kind=realType), dimension(nw) :: dJdWlocal

      ! idxmgb - global block row index

      integer(kind=intType) :: idxmgb

      character(len=2*maxStringLen) :: errorMessage
	
      integer :: ierr,sps,n

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj
      
      !
      !     ******************************************************************
      !     *                                                                *
      !     * Begin execution.                                               *
      !     *                                                                *
      !     ******************************************************************
      !

      !Specify the existence of a second halo.
      secondhalo = .true.
#ifndef USE_NO_PETSC

     ! Get the initial time.

      call cpu_time(time(1))
      ! Determine the reference point for the moment computation in
      ! meters.

      refPoint(1) = LRef*pointRef(1)
      refPoint(2) = LRef*pointRef(2)
      refPoint(3) = LRef*pointRef(3)

      ! Initialize the force and moment coefficients to 0 as well as
      ! yplusMax.

      ClAdj = zero
      CDAdj = zero
      CmxAdj = zero 
      CmyAdj = zero
      CmzAdj = zero
      CfxAdj = zero
      CfyAdj = zero
      CfzAdj = zero

      yplusMax = zero

      spectralLoopAdj: do sps=1,nTimeIntervalsSpectral

         ! Initialize the force and moment coefficients to 0.

         cFp(1) = zero; cFp(2) = zero; cFp(3) = zero
         cFv(1) = zero; cFv(2) = zero; cFv(3) = zero
         cMp(1) = zero; cMp(2) = zero; cMp(3) = zero
         cMv(1) = zero; cMv(2) = zero; cMv(3) = zero

         domainLoopAD: do nn=1,nDom

            ! Set some pointers to make the code more readable.
            !print *,'setting pointers'
            call setPointersAdj(nn,level,sps)
            !print *,'allocating memory'
            allocate(xAdj(0:ie,0:je,0:ke,3), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for xAdj.")

            allocate(xAdjB(0:ie,0:je,0:ke,3), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for xAdjB.")

            allocate(wAdj(0:ib,0:jb,0:kb,nw), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for wAdj.")

            allocate(wAdjB(0:ib,0:jb,0:kb,nw), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for wAdjB.")

            allocate(pAdj(0:ib,0:jb,0:kb), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for pAdj.")

            !print *,'finished allocating',nn,level,sps
            righthanded = flowDoms(nn,level,sps)%righthanded
            
            ! Copy the coordinates into xAdj and
            ! Compute the face normals on the subfaces
            call copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
           rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
           timerefAdj,pInfCorrAdj,nn,level,sps,liftIndex)

        
            bocoLoop: do mm=1,nBocos
	    
            ! Initialize the seed for reverse mode. Select case based on
            ! desired cost function.
		    
	    !print *,'selecting case',costfunction
            select case (costFunction)
		               
            !case (costFuncLiftCoef)
            case (costFunccl0,costFuncclalpha)
		!print *,'computing lift derivatives'
               ClAdjB = 0
               ClAdjB(sps) = 1
               CDAdjB = 0
               CfxAdjB = 0
               CfyAdjB = 0
               CfzAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0     
		!print *,'Case CL',clAdjb
	 
             
            case (costFunccd0,costFunccdalpha)
	       !print *,'computing drag derivatives'

               ClAdjB = 0
               CDAdjB = 0
               CDAdjB(sps) = 1
               CfxAdjB = 0
               CfyAdjB = 0
               CfzAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0 

!!$            case (costFuncForceXCoef)
!!$
!!$               ClAdjB = 0
!!$               CDAdjB = 0
!!$	       CfxAdjB(sps) = 1
!!$               CfyAdjB = 0
!!$               CfzAdjB = 0
!!$               CmxAdjB = 0
!!$               CmyAdjB = 0
!!$               CmzAdjB = 0
               
!!$            case (costFuncForceYCoef)
!!$
!!$               ClAdjB = 0
!!$               CDAdjB = 0
!!$   	       CfxAdjB = 0
!!$               CfyAdjB(sps) = 1
!!$               CfzAdjB = 0
!!$               CmxAdjB = 0
!!$               CmyAdjB = 0
!!$               CmzAdjB = 0      
               
!!$            case (costFuncForceZCoef)
!!$
!!$               ClAdjB = 0
!!$               CDAdjB = 0
!!$	       CfxAdjB = 0
!!$               CfyAdjB = 0
!!$               CfzAdjB(sps) = 1
!!$               CmxAdjB = 0
!!$               CmyAdjB = 0
!!$               CmzAdjB = 0                
!!$               
!!$            case (costFuncMomXCoef)
!!$
!!$               ClAdjB = 0
!!$               CDAdjB = 0
!!$               CfxAdjB = 0
!!$               CfyAdjB = 0
!!$               CfzAdjB = 0
!!$               CmxAdjB(sps) = 1
!!$               CmyAdjB = 0
!!$               CmzAdjB = 0
!!$               
!!$            case (costFuncMomYCoef)
!!$
!!$               ClAdjB = 0
!!$               CDAdjB = 0
!!$               CfxAdjB = 0
!!$               CfyAdjB = 0
!!$               CfzAdjB = 0
!!$               CmxAdjB = 0
!!$               CmyAdjB(sps) = 1
!!$               CmzAdjB = 0      
               
            case (costFunccm0,costFunccmzalpha)

               ClAdjB = 0
               CDAdjB = 0
               CfxAdjB = 0
               CfyAdjB = 0
               CfzAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB =0
               CmzAdjB(sps) = 1      
                              
            end select
            
            wAdjB(:,:,:,:) = zero ! > return dCf/dw
            xAdjB(:,:,:,:) = zero ! > return dCf/dx
       

           
               
               ! Determine the range of cell indices of the owned cells
               ! Notice these are not the node indices
               iiBeg = BCData(mm)%icBeg
               iiEnd = BCData(mm)%icEnd
               jjBeg = BCData(mm)%jcBeg
               jjEnd = BCData(mm)%jcEnd
               
               i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
               j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd

               t = timeUnsteadyRestart
               
               if(equationMode == timeSpectral) then
                  do nnn=1,nSections
                     !t(nnn) = t(nnn) + (sps2-1)*sections(nnn)%timePeriod &
                     !     /         real(nTimeIntervalsSpectral,realType)
                     t(nnn) = t(nnn) + (sps-1)*sections(nnn)%timePeriod &
                          /         (nTimeIntervalsSpectral*1.0)!to make denomenator a real number...
                  enddo
               endif
               !print *,'cladjb',cladjb,cdadjb
               !print *,'before',sum(wadjb)
	       call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfxadjb&
&  , cfyadj, cfyadjb, cfzadj, cfzadjb, cmxadj, cmxadjb, cmyadj, cmyadjb&
&  , cmzadj, cmzadjb, yplusmax, refpoint, cladj, cladjb, cdadj, cdadjb, &
&  nn, level, sps, cfpadj, cmpadj, righthanded, secondhalo, alphaadj, &
&  alphaadjb, betaadj, betaadjb, machadj, machadjb, machcoefadj, &
&  machcoefadjb, machgridadj, machgridadjb, prefadj, rhorefadj, &
&  pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, murefadj, timerefadj, &
&  pinfcorradj, rotcenteradj, rotrateadj, rotrateadjb, liftindex, t)



		!print *,'CL',cladj,cladjb,p(1,1,1),wadj(1,1,1,1)
           
        
            !*******************************************************
            !                                                      *
            ! Transfer the block Jacobians to the PETSc vector.    *
            !                                                      *
            !*******************************************************
            
            ! When using the block vector format, one can insert
            ! elements more efficiently using the block variant.
            
            ! VecSetValuesBlocked - Inserts or adds blocks of values
            !                    into certain locations of a vector.
            ! Synopsis
            ! 
            ! #include "petscvec.h" 
            ! call VecSetValuesBlocked(Vec x,PetscInt ni,          &
            !            const PetscInt ix[],const PetscScalar y[],&
            !            InsertMode iora,PetscErrorCode ierr) 
            !
            ! Not Collective
            !
            ! Input Parameters
            !   x    - vector to insert in
            !   ni   - number of blocks to add
            !   ix   - indices where to add in block count,
            !          rather than element count
            !   y    - array of values
            !   iora - either INSERT_VALUES or ADD_VALUES,
            !          where ADD_VALUES adds values to any existing
            !          entries, and INSERT_VALUES replaces existing
            !          entries with new values
            ! Notes
            ! VecSetValuesBlocked() sets x[bs*ix[i]+j] = y[bs*i+j],
            !   for j=0,...,bs, for i=0,...,ni-1. where bs was set
            !   with VecSetBlockSize().
                            !
            ! Calls to VecSetValuesBlocked() with the INSERT_VALUES
            !   and ADD_VALUES options cannot be mixed without
            !   intervening calls to the assembly routines.
            !
            ! These values may be cached, so VecAssemblyBegin() and
            !   VecAssemblyEnd() MUST be called after all calls to
            !   VecSetValuesBlocked() have been completed.
            !
            ! VecSetValuesBlocked() uses 0-based indices in Fortran
            !   as well as in C. 
            !
            ! see .../petsc/docs/manualpages/Vec/VecSetValuesBlocked.html
            
            ! Global vector block row mgb function of node indices.
            !
            ! VecSetValuesBlocked() uses 0-based row numbers but the
            ! global node numbering already accounts for that since
            ! it starts at node 0.

            do kcell = 0,kb
               do jcell = 0,jb
                  do icell = 0,ib  
                     do n = 1,nw
                        !idxmgb = globalCell(icell,jcell,kcell)
                        idxmgb = globalCell(iCell,jCell,kCell)*nw+n
                        !print *,'index',idxmgb,sps,nn,icell,jcell,kcell,n,globalCell(iCell,jCell,kCell)
                        test = sum(wAdjB(icell,jcell,kcell,:))
                        !if(abs(test)>1e-12)then
                        !if(abs(wAdjB(icell,jcell,kcell,n))>1e-12)then
                        !   print *,'wadjb',wAdjB(icell,jcell,kcell,n),idxmgb,sps,nn,globalcell(icell,jcell,kcell),icell,jcell,kcell
                        !endif
                        !if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal*nTimeIntervalsSpectral) then
                        if (idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<=nCellsGlobal*nTimeIntervalsSpectral*nw) then
                           !print *,'test',test,idxmgb
                           !print *,'index',idxmgb,sps,nn,icell,jcell,kcell,n,globalCell(iCell,jCell,kCell)
                           !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
                           if (wAdjb(icell,jcell,kcell,n).ne.0.0)then
                              !print *,'mat set',wadjb(icell,jcell,kcell,n),idxmgb
                              call MatSetValues(dCdw, 1, sps-1, 1, idxmgb-1,   &
                                   wAdjb(icell,jcell,kcell,n), ADD_VALUES, PETScIerr)
                              
                              if( PETScIerr/=0 ) then
                                 write(errorMessage,99) &
                                      "Error in MatSetValues for global cell", &
                                      idxmgb
                                 call terminate("setupADjointdCdwStability", &
                                      errorMessage)
                              endif
                           endif
                        endif
                     enddo
                  end do
               enddo
            enddo
	    
            enddo bocoLoop
            
            !===============================================================
            
            !print *,' deallocating'
            ! Deallocate the xAdj.
            deallocate(pAdj, stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("verifydCfdx", &
                 "Deallocation failure for xAdj.")
            
            ! Deallocate the xAdj.
            deallocate(wAdj, stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("verifydCfdx", &
                 "Deallocation failure for xAdj.") 
            
            deallocate(wAdjB, stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("verifydCfdx", &
                 "Deallocation failure for xAdj.") 
            ! Deallocate the xAdj.
            deallocate(xAdj, stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("verifydCfdx", &
                 "Deallocation failure for xAdj.") 
            
            deallocate(xAdjB, stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("verifydCfdx", &
                 "Deallocation failure for xAdj.") 
            !print *,'finishhed deallocating'
            
         enddo domainLoopAD
         
      enddo spectralLoopAdj


!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc matrix assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! MatAssemblyBegin - Begins assembling the matrix. This routine
      !  should be called after completing all calls to MatSetValues().
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatAssemblyBegin(Mat mat, &
      !                                            MatAssemblyType type)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat  - the matrix
      !   type - type of assembly, either MAT_FLUSH_ASSEMBLY or
      !          MAT_FINAL_ASSEMBLY
      ! Notes
      ! MatSetValues() generally caches the values. The matrix is ready
      !  to use only after MatAssemblyBegin() and MatAssemblyEnd() have
      !  been called. Use MAT_FLUSH_ASSEMBLY when switching between
      !  ADD_VALUES and INSERT_VALUES in MatSetValues(); use
      !  MAT_FINAL_ASSEMBLY for the final assembly before using the
      !  matrix.
      !
      ! see .../petsc/docs/manualpages/Mat/MatAssemblyBegin.html

      call MatAssemblyBegin(dCdw,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointdCdwStability", &
                       "Error in MatAssemblyBegin dCdw")

      ! MatAssemblyEnd - Completes assembling the matrix. This routine
      !                  should be called after MatAssemblyBegin().
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatAssemblyEnd(Mat mat,&
      !                                            MatAssemblyType type)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat  - the matrix
      !   type - type of assembly, either MAT_FLUSH_ASSEMBLY or
      !          MAT_FINAL_ASSEMBLY
      !
      ! see .../petsc/docs/manualpages/Mat/MatAssemblyEnd.html

      call MatAssemblyEnd  (dCdw,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointdCdwStability", &
                       "Error in MatAssemblyEnd dCdw")

      ! Let PETSc know that the dRda matrix retains the same nonzero 
      ! pattern, in case the matrix is assembled again, as for a new
      ! point in the design space.

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted,row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetOption(Mat mat,MatOption op,PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   option - the option, one of those listed below (and possibly
      !     others), e.g., MAT_ROWS_SORTED, MAT_NEW_NONZERO_LOCATION_ERR
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetOption.html
      ! or PETSc users manual, pp.52
#ifndef USE_PETSC_3
      call MatSetOption(dCdw,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointdcdwStability", &
                       "Error in MatSetOption dCdw")
#endif
      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dC/dw matrix time (s) =", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled matrix.                                *
!     *                                                                *
!     ******************************************************************
!
      ! MatView - Visualizes a matrix object.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatView(Mat mat, &
      !                                              PetscViewer viewer)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   viewer - visualization context
      !
      ! Notes
      ! The available visualization contexts include
      !  PETSC_VIEWER_STDOUT_SELF  - standard output (default)
      !  PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where
      !                         only the first processor opens the file.
      !                         All other processors send their data to
      !                         the first processor to print.
      !  PETSC_VIEWER_DRAW_WORLD- graphical display of nonzero structure
      !
      ! see .../petsc/docs/manualpages/Mat/MatView.html
      ! or PETSc users manual, pp.57,148

      if( debug ) then
        !call MatView(dCdw,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
        call MatView(dCdw,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupADjointdCdwStability", "Error in MatView")
        !pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)
      ! Output format.

20    format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

    end subroutine setupADjointdCdwStability

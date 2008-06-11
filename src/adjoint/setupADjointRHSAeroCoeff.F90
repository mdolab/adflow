!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointRHSAeroCoeff.F90                    *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 10-04-2006                                      *
!     * Last modified: 06-09-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupADjointRHSAeroCoeff(level,costFunction)
!subroutine setupADjointRHSAeroCoeff(level,sps,costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * This routine computes the right hand side of the discrete      *
!     * ADjoint problem when using aerodynamic coefficients (CD,CL,CM) *
!     * as cost functions:                                             *
!     *                         dPartial(J)                            *
!     *   J = CD, CL or CM  =>  -----------                            *
!     *                         dPartial(W)                            *
!     *                                                                *
!     *                                                                *
!     * This routine is based on /solver/convergenceInfo.f90 and       *
!     * /solver/forcesAndMoments.f90 (as of 05-16-2006).               *
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

      integer(kind=intType) :: nn, mm, i, j
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
      integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End

      integer(kind=intType) :: iCell, jCell, kCell

      real(kind=realType) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,&
                             &CmxAdj,CmyAdj,CmzAdj 

      real(kind=realType) :: ClAdjB,CdAdjB,CfxAdjB,CfyAdjB,CfzAdjB,&
                             &CmxAdjB,CmyAdjB,CmzAdjB  

      real(kind=realType) :: yplusMax, test

      real(kind=realType), dimension(3) :: refPoint

      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,xAdjB
      real(kind=realType), dimension(:,:,:,:), allocatable :: wAdj,wAdjB
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj

      REAL(KIND=REALTYPE) :: machadj, machcoefadj, uinfadj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb

      logical :: contributeToForce, viscousSubface,secondHalo,righthanded

      ! dJ/dw row block
      
      real(kind=realType), dimension(nw) :: dJdWlocal

      ! idxmgb - global block row index

      integer(kind=intType) :: idxmgb

      character(len=2*maxStringLen) :: errorMessage
	
      integer :: ierr,sps

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
              MachAdj,machCoefAdj,prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
              rhoinfAdj, pinfAdj,murefAdj, timerefAdj,pInfCorrAdj,nn,level,sps)
!copyADjointForcesStencil(wAdj,xAdj,nn,level,sps)
        
            bocoLoop: do mm=1,nBocos
	    
            ! Initialize the seed for reverse mode. Select case based on
            ! desired cost function.
		    
	    !print *,'selecting case',costfunction
            select case (costFunction)
		               
            case (costFuncLiftCoef)
		!print *,'computing lift derivatives'
               ClAdjB = 1
               CDAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0     
		!print *,'Case CL',clAdjb
	 
             
            case (costFuncDragCoef)
	       !print *,'computing drag derivatives'

               ClAdjB = 0
               CDAdjB = 1
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0 

            case (costFuncForceXCoef)

               ClAdjB = 0
               CDAdjB = 0
	       !CfxAdjB = 1
               !CfyAdjB = 0
               !CfzAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0
               
            case (costFuncForceYCoef)

               ClAdjB = 0
               CDAdjB = 0
   	       !CfxAdjB = 0
               !CfyAdjB = 1
               !CfzAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0      
               
            case (costFuncForceZCoef)

               ClAdjB = 0
               CDAdjB = 0
	       !CfxAdjB = 0
               !CfyAdjB = 0
               !CfzAdjB = 1
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 0                
               
            case (costFuncMomXCoef)

               ClAdjB = 0
               CDAdjB = 0
               CmxAdjB = 1
               CmyAdjB = 0
               CmzAdjB = 0
               
            case (costFuncMomYCoef)

               ClAdjB = 0
               CDAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 1
               CmzAdjB = 0      
               
            case (costFuncMomZCoef)

               ClAdjB = 0
               CDAdjB = 0
               CmxAdjB = 0
               CmyAdjB = 0
               CmzAdjB = 1      
                              
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
               !print *,'cladjb',cladjb
               
	       call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, prefadj, &
&  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, murefadj, &
&  timerefadj, pinfcorradj)

!		print *,'output',xadj, xadjb, wadj, wadjb, padj, iibeg, &
!&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
!&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
!&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
!&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
!&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, prefadj, &
!&  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, murefadj, &
!&  timerefadj, pinfcorradj

!COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
!&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
!&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
!&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
!&  cmpadj, righthanded, alphaadj, alphaadjb, betaadj, betaadjb, machadj&
!&  , machadjb, machcoefadj, machcoefadjb, prefadj, rhorefadj, pinfdimadj&
!&  , rhoinfdimadj, rhoinfadj, pinfadj, murefadj, timerefadj, pinfcorradj)
!               call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
!&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
!&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
!&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
!&  cmpadj, righthanded, alphaadj, alphaadjb, betaadj, betaadjb, machadj&
!&  , machadjb, machcoefadj, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj&
!&  , rhoinfadj, pinfadj, murefadj, timerefadj, pinfcorradj)

!		print *,'setup',xadj, xadjb, wadj, wadjb, padj, iibeg, &
!&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
!&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
!&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
!&  cmpadj, righthanded, alphaadj, alphaadjb, betaadj, betaadjb, machadj&
!&  , machadjb, machcoefadj, machcoefadjb, prefadj, rhorefadj, pinfdimadj&
!&  , rhoinfdimadj, rhoinfadj, pinfadj, murefadj, timerefadj, pinfcorradj

!COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
!                    &  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
!                    &  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
!                    &  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
!                    &  cmpadj, righthanded)

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
                     idxmgb = globalCell(icell,jcell,kcell)
                     
                     test = sum(wAdjB(icell,jcell,kcell,:))
                     !print *,'test',wAdjB(icell,jcell,kcell,:)!test
                     if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
	                !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
                       dJdWlocal(:) = wAdjB(icell,jcell,kcell,:)
                     
                       call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
                                                ADD_VALUES, PETScIerr)
       !                 call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
       !                      ADD_VALUES, PETScIerr)
            
                        if( PETScIerr/=0 ) then
                           write(errorMessage,99) &
                                "Error in VecSetValuesBlocked for global node", &
                                idxmgb
                           call terminate("setupADjointRHSAeroCoeff", &
                                errorMessage)
                        endif
		     endif
                  enddo
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

      ! Output format.

   99 format(a,1x,i6)

#endif

      end subroutine setupADjointRHSAeroCoeff

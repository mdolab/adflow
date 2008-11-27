!
!     ******************************************************************
!     *                                                                *
!     * File:          setupGradientRHSExtra.F90                       *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 08-23-2005                                      *
!     * Last modified: 03-01-2007                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupGradientRHSExtra(level,costFunction)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the cost function partial sensitivity with respect to  *
  !     * the extra design variables:                                    *
  !     *                                                                *
  !     *                                dpartial J                      *
  !     *   dJda(nDesignExtra) = ----------------------------            *
  !     *                        dpartial alpha(nDesignExtra)            *
  !     *                                                                *
  !     *   where nDesignExtra includes:                                 *
  !     *   > angle of attack                                            *
  !     *   > side slip angle                                            *
  !     *   > Mach Number                                                *
  !     *                                                                *
  !     * This routine has to be consistent with the design variable     *
  !     * ordering done in "setDesignVaribles".                          *
  !     *                                                                *
  !     * Supported cost functions J:                                    *
  !     *   - lift coefficient Cl  if "costFunction" = costFuncLiftCoef  *
  !     *   - drag coefficient Cd  if "costFunction" = costFuncDragCoef  *
  !     *   - x-force coef.   Cfx if "costFunction" = costFuncForceXCoef *
  !     *   - y-force coef.   Cfy if "costFunction" = costFuncForceYCoef *
  !     *   - z-force coef.   Cfz if "costFunction" = costFuncForceZCoef *
  !     *   - x-moment coef.   Cmx if "costFunction" = costFuncMomXCoef  *
  !     *   - y-moment coef.   Cmy if "costFunction" = costFuncMomYCoef  *
  !     *   - z-moment coef.   Cmz if "costFunction" = costFuncMomZCoef  *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use ADjointVars
  use blockPointers  ! il,jl,kl, globalNode
  use communication  ! myID, nProc
  use inputPhysics   ! liftDirection, dragDirection
  use flowVarRefState!nw
  use inputTimeSpectral !nTimeIntervalsSpectral
  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, costFunction
  !
  !     Local variables.
  !
  integer(kind=intType) :: sps,ierr
  integer(kind=intType) :: n

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj

  ! dJ/da local vector at node (iNode,jNode,kNode)

  real(kind=realType) :: dJdaLocal

  ! idxmg - global row index

  integer(kind=intType) :: idxmg,liftIndex

  ! auxiliar variables to compute dJ/dalpha and dJ/dbeta
  real(kind=realType), dimension(3) :: cFp, cFv
  real(kind=realType), dimension(3) :: cMp, cMv

  real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
  real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

  integer(kind=intType) :: nn, mm, i, j
  integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
  integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End

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
  REAL(KIND=REALTYPE) :: rotcenteradj(3), rotrateadj(3), rotrateadjb(3)

  logical :: secondHalo,exchangeTurb,correctfork,finegrid,righthanded
  integer(kind=intType):: discr

  character(len=2*maxStringLen) :: errorMessage
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

  ! Reset the RHS vector dJ/da by assigning the value zero to all
  ! its components.

  ! VecSet - Sets all components of vector to a single scalar value.

  call VecSet(dJda,PETScZero,PETScIerr)

  if( PETScIerr/=0 ) &
       call terminate("setupGradientRHSExtra", "Error in VecSet")


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

     !zero the force components
     cFpAdj(1) = zero; cFpAdj(2) = zero; cFpAdj(3) = zero
     cFvAdj(1) = zero; cFvAdj(2) = zero; cFvAdj(3) = zero
     cMpAdj(1) = zero; cMpAdj(2) = zero; cMpAdj(3) = zero
     cMvAdj(1) = zero; cMvAdj(2) = zero; cMvAdj(3) = zero
 
     domainLoopAD: do nn=1,nDom

        ! Set some pointers to make the code more readable.
 
        call setPointersAdj(nn,level,sps)
  
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

        righthanded = flowDoms(nn,level,sps)%righthanded

        ! Copy the coordinates into xAdj 
       
        call copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
           rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj, timerefAdj,&
           pInfCorrAdj,nn,level,sps,liftIndex)

 
        wAdjB(:,:,:,:) = zero ! > return dCf/dw
        xAdjB(:,:,:,:) = zero ! > return dCf/dx
        !alphaadjb = zero
        !betaadjb = zero
        !machadjb = zero
	!machcoefadjb = zero	
	!rotrateadjb(:)=zero

	bocoLoop: do mm=1,nBocos
	   ! Initialize the seed for reverse mode. Select case based on
           ! desired cost function.
        
           ClAdjB = 0
           CDAdjB = 0
           !CfxAdjB = 0
           !CfyAdjB = 0
           !CfzAdjB = 0
           CmxAdjB = 0
           CmyAdjB = 0
           CmzAdjB = 0  
           
           select case (costFunction)
           case (costFuncLiftCoef)
              
              ClAdjB = 1   
              
 
           case (costFuncDragCoef)
              
              CDAdjB = 1

           case (costFuncForceXCoef)

              !CfxAdjB = 1
 

           case (costFuncForceYCoef)

              !CfyAdjB = 1

           case (costFuncForceZCoef)

              !CfzAdjB = 1

           case (costFuncMomXCoef)

              CmxAdjB = 1

           case (costFuncMomYCoef)

              CmyAdjB = 1

           case (costFuncMomZCoef)
              CmzAdjB = 1      

           end select
	
	   
           ! Determine the range of cell indices of the owned cells
           ! Notice these are not the node indices
           
	   iiBeg = BCData(mm)%icBeg
           iiEnd = BCData(mm)%icEnd
           jjBeg = BCData(mm)%jcBeg
           jjEnd = BCData(mm)%jcEnd

           i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
           j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd
           

	   call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfyadj, &
&  cfzadj, cmxadj, cmxadjb, cmyadj, cmyadjb, cmzadj, cmzadjb, yplusmax, &
&  refpoint, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, prefadj, &
&  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, murefadj, &
&  timerefadj, pinfcorradj, rotcenteradj, rotrateadj, rotrateadjb, &
&  liftindex)


	   enddo bocoLoop

           !***************************************************
           !        set angle of of attack (AOA) derivative
           !***********************************************

           dJdaLocal = alphaadjb
           
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignAOA - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
	   call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif


           !
           !     ******************************************************************
           !     *                                                                *
           !     * Side slip angle > beta.                                        *
           !     *                                                                *
           !     ******************************************************************
           !

           dJdaLocal = betaadjb
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignSSA - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
	   call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif

           !
           !     ******************************************************************
           !     *                                                                *
           !     * Mach Number derivative.                                        *
           !     *                                                                *
           !     ******************************************************************
           !

           dJdaLocal = machadjb+machcoefadjb
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignMach - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
	   call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif

	   ! 
           !     ******************************************************************
           !     *                                                                *
           !     * Rot X derivative.                                              *
           !     *                                                                *
           !     ******************************************************************
           !

           dJdaLocal = rotrateadjb(1)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignRotX - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
	   call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif

           !
           !     ******************************************************************
           !     *                                                                *
           !     * Rot Y derivative.                                              *
           !     *                                                                *
           !     ******************************************************************
           !

           dJdaLocal = rotrateadjb(2)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignRotY - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
           call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif
 
           !
           !     ******************************************************************
           !     *                                                                *
           !     * Rot Z derivative.                                              *
           !     *                                                                *
           !     ******************************************************************
           !

           dJdaLocal = rotrateadjb(3)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignMach - 1

           ! Transfer data to PETSc vector

           !call VecSetValue(dJda, idxmg, dJdaLocal, &
           !     ADD_VALUES, PETScIerr)
           call VecSetValue(dJda, idxmg, dJdaLocal, &
                INSERT_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in VecSetValue for global node", idxmg
              call terminate("setupGradientRHSExtra", errorMessage)
           endif

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
        

     enddo domainLoopAD

  enddo spectralLoopAdj


!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc vector assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! VecAssemblyBegin - Begins assembling the vector. This routine
      ! should be called after completing all calls to VecSetValues().

      call VecAssemblyBegin(dJda,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientRHSExtra", &
                       "Error in VecAssemblyBegin")

      ! VecAssemblyEnd - Completes assembling the vector. This routine
      ! should be called after VecAssemblyBegin().

      call VecAssemblyEnd  (dJda,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientRHSExtra", &
                       "Error in VecAssemblyEnd")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dJ/da vector time (s) =", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled vector.                                *
!     *                                                                *
!     ******************************************************************
!
      ! VecView - Views a vector object.

      if( debug ) then
        call VecView(dJda,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupGradientRHS", "Error in VecView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

      end subroutine setupGradientRHSExtra

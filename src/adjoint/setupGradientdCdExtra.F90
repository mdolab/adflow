!
!     ******************************************************************
!     *                                                                *
!     * File:          setupGradientdCdExtra.F90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-30-2009                                      *
!     * Last modified: 11-30-2009                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupGradientdCdExtra(level,costFunction)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the force coefficient partial sensitivity with respect *
  !     * to the extra design variables:                                 *
  !     *                                                                *
  !     *                                  dpartial C(sps)               *
  !     *   dCda(sps,nDesignExtra) = ----------------------------        *
  !     *                            dpartial alpha(nDesignExtra)        *
  !     *                                                                *
  !     *   where nDesignExtra includes:                                 *
  !     *   > angle of attack                                            *
  !     *   > side slip angle                                            *
  !     *   > Mach Number                                                *
  !     *   > Rotation Rate(p,q,r)                                       *
  !     *                                                                *
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
  use monitor        !timeunsteadyrestart
  use section        !nsection,sections%
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

  integer(kind=intType) :: idxmg,liftIndex,nnn

  ! auxiliar variables to compute dJ/dalpha and dJ/dbeta
  real(kind=realType), dimension(3) :: cFp, cFv
  real(kind=realType), dimension(3) :: cMp, cMv

  real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
  real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

  integer(kind=intType) :: nn, mm, i, j
  integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
  integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End

  real(kind=realType),dimension(nTimeIntervalsSpectral) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,&
       &CmxAdj,CmyAdj,CmzAdj 
  
  real(kind=realType),dimension(nTimeIntervalsSpectral) :: ClAdjB,CdAdjB,CfxAdjB,CfyAdjB,CfzAdjB,&
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
  REAL(KIND=REALTYPE) :: rotcenteradj(3), rotcenteradjb(3), rotrateadj(3&
       &  ), rotrateadjb(3)
  REAL(KIND=REALTYPE) :: pointrefadj(3), pointrefadjb(3), rotpointadj(3)&
       &  , rotpointadjb(3)

  logical :: secondHalo,exchangeTurb,correctfork,finegrid,righthanded
  integer(kind=intType):: discr

  real(kind=realType), dimension(nSections) :: t     

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
!print *,'setupgradientdcdextra'
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
     !print *,'spectral loop',sps
     ! Initialize the force and moment coefficients to 0.

     !zero the force components
     cFpAdj(1) = zero; cFpAdj(2) = zero; cFpAdj(3) = zero
     cFvAdj(1) = zero; cFvAdj(2) = zero; cFvAdj(3) = zero
     cMpAdj(1) = zero; cMpAdj(2) = zero; cMpAdj(3) = zero
     cMvAdj(1) = zero; cMvAdj(2) = zero; cMvAdj(3) = zero
 
     domainLoopAD: do nn=1,nDom
      !  print *,'domain loop',nn
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
           MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
           rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
           timerefAdj,pInfCorrAdj,pointRefAdj,rotPointAdj,nn,level,sps,&
           liftIndex)

 
        wAdjB(:,:,:,:) = zero ! > return dCf/dw
        xAdjB(:,:,:,:) = zero ! > return dCf/dx
        !alphaadjb = zero
        !betaadjb = zero
        !machadjb = zero
	!machcoefadjb = zero	
	!rotrateadjb(:)=zero

	bocoLoop: do mm=1,nBocos
           !print *,'bc loop',mm
	   ! Initialize the seed for reverse mode. Select case based on
           ! desired cost function.
        
           ClAdjB = 0
           CDAdjB = 0
           CfxAdjB = 0
           CfyAdjB = 0
           CfzAdjB = 0
           CmxAdjB = 0
           CmyAdjB = 0
           CmzAdjB = 0  
           
           select case (costFunction)
           case (costFuncLiftCoef)
              
              ClAdjB(sps) = 1   
              
 
           case (costFuncDragCoef)
              
              CDAdjB(sps) = 1

           case (costFuncForceXCoef)

              CfxAdjB(sps) = 1
 

           case (costFuncForceYCoef)

              CfyAdjB(sps) = 1

           case (costFuncForceZCoef)

              CfzAdjB(sps) = 1

           case (costFuncMomXCoef)

              CmxAdjB(sps) = 1

           case (costFuncMomYCoef)

              CmyAdjB(sps) = 1

           case (costFuncMomZCoef)
              CmzAdjB(sps) = 1      

           case (costFunccl0,costFuncclalpha)
              ClAdjB(sps) = 1

           case (costFunccd0,costFunccdalpha)
              CdAdjB(sps) = 1

           case (costFunccm0,costFunccmzalpha)
              CmzAdjB(sps) = 1     
              
           end select
	
	   
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

	   call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfxadjb&
&  , cfyadj, cfyadjb, cfzadj, cfzadjb, cmxadj, cmxadjb, cmyadj, cmyadjb&
&  , cmzadj, cmzadjb, yplusmax, pointrefadj, pointrefadjb, rotpointadj, &
&  rotpointadjb, cladj, cladjb, cdadj, cdadjb, nn, level, sps, cfpadj, &
&  cmpadj, righthanded, secondhalo, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machcoefadjb, machgridadj, &
&  machgridadjb, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj&
&  , pinfadj, murefadj, timerefadj, pinfcorradj, rotcenteradj, &
&  rotcenteradjb, rotrateadj, rotrateadjb, liftindex, t)


	   enddo bocoLoop
           !print *,'end bc loop'
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
     !print *,'end domain loop'
 
	
           !***************************************************
           !        set angle of of attack (AOA) derivative
           !***********************************************

           dJdaLocal = alphaadjb
           
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignAOA - 1

           ! Transfer data to PETSc matrix
           !print *,'alpha',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValue for alpha", idxmg
              call terminate("setupGradientdcdExtra", errorMessage)
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
           !print *,'sideslip',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for beta", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
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
           !print *,'mach',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for mach", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

!
!     ******************************************************************
!     *                                                                *
!     * Mach Number Grid derivative.                                   *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = machcoefadjb+machgridadjb
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignMachGrid - 1

           ! Transfer data to PETSc vector
           !print *,'machgrid',idxmg,sps-1,dJdaLocal
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)
     
           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for MachGrid", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
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
           !print *,'rotx',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for roterate(1)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
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
           !print *,'roty',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for roterate(2)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
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

           idxmg = nDesignRotZ - 1

           ! Transfer data to PETSc vector
           !print *,'rotz',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)
           
           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for roterate(3)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

! 
!     ******************************************************************
!     *                                                                *
!     * Rotcen X derivative.                                           *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = rotcenteradjb(1)+rotpointadjb(1)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignRotCenX - 1

           ! Transfer data to PETSc vector
           !print *,'rotx',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for rotcen(1)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

!
!     ******************************************************************
!     *                                                                *
!     * RotCen Y derivative.                                           *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = rotcenteradjb(2)+rotpointadjb(2)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignRotCenY - 1

           ! Transfer data to PETSc vector
           !print *,'roty',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for rotcen(2)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif
 
!
!     ******************************************************************
!     *                                                                *
!     * RotCen Z derivative.                                           *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = rotcenteradjb(3)+rotpointadjb(3)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignRotCenZ - 1

           ! Transfer data to PETSc vector
           !print *,'rotz',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)
           
           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for rotcenter(3)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

! 
!     ******************************************************************
!     *                                                                *
!     * PointRef X derivative.                                         *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = pointrefadjb(1)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignPointRefX - 1

           ! Transfer data to PETSc vector
           !print *,'rotx',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for pointref(1)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

!
!     ******************************************************************
!     *                                                                *
!     * PointRef Y derivative.                                         *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = pointrefadjb(2)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignPointRefY - 1

           ! Transfer data to PETSc vector
           !print *,'roty',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)

           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for pointref(2)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif
 
!
!     ******************************************************************
!     *                                                                *
!     * PointRef Z derivative.                                         *
!     *                                                                *
!     ******************************************************************
!

           dJdaLocal = pointrefadjb(3)
	   
           ! Set the corresponding single entry of the PETSc vector dJda.

           ! Global vector row idxmg function of design variable index.

           idxmg = nDesignPointRefZ - 1

           ! Transfer data to PETSc vector
           !print *,'rotz',idxmg,sps-1
           call MatSetValues(dCda,1,sps-1,1,idxmg, dJdaLocal, &
                ADD_VALUES, PETScIerr)
           
           if( PETScIerr/=0 ) then
              write(errorMessage,99) &
                   "Error in MatSetValues for pointref(3)", idxmg
              call terminate("setupGradientdCdExtra", errorMessage)
           endif

        enddo spectralLoopAdj
        !print *,'finished loops'
!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc vector assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! VecAssemblyBegin - Begins assembling the vector. This routine
      ! should be called after completing all calls to VecSetValues().
        !print *,'assemblybegin',petscierr
      call MatAssemblyBegin(dCda,MAT_FINAL_ASSEMBLY,PETScIerr)
      !print *,'assemblybegin1',petscierr
      if( PETScIerr/=0 ) &
        call terminate("setupGradientdCdExtra", &
                       "Error in MatAssemblyBegin")
!      print *,'assemblybegin2',petscierr
      ! VecAssemblyEnd - Completes assembling the vector. This routine
      ! should be called after VecAssemblyBegin().

      call MatAssemblyEnd(dCda,MAT_FINAL_ASSEMBLY,PETScIerr)
 !     print *,'assemblyend',petscierr
      if( PETScIerr/=0 ) &
        call terminate("setupGradientdCdExtra", &
                       "Error in MatAssemblyEnd")

  !    print *,'end assembly'
      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)

   !   print *,'time',petscrank,timeadj
      if( PETScRank==0 ) &
        write(*,20) "Assembling dC/da vector time (s) =", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled vector.                                *
!     *                                                                *
!     ******************************************************************
!
      ! VecView - Views a vector object.

      if( debug ) then
        !call MatView(dCda,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
	call MatView(dCda,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupGradientdCdExtra", "Error in MatView")
        !pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

    end subroutine setupGradientdCdExtra

!
!     ******************************************************************
!     *                                                                *
!     * File:          computObjPartial.f90                            *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 11-02-2010                                      *
!     * Last modified: 11-02-2010                                      *
!     *                                                                *
!     ******************************************************************

subroutine computeObjPartials(costFunction,pts,npts,nTS,usedJdw,usedJdx)
#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * This routine computes the partial derivative of objective,     *
  !     * costFunction wrt the surface points pts, and the states w for  *
  !     * time spectral level sps.                                       *
  !     *                                                                *
  !     *         dPartial(J)        dPartial(J)                         *
  !     *         -----------        -----------                         *
  !     *         dPartial(W)        dPartial(X)                         *
  !     *                                                                *
  !     * dJ/dW is set in PETSc as this the RHS vector for the adjoint.  *
  !     * dI/dx is returned as this must be run through the (external)   *
  !     * warping first before it can be used with the actual shape      *
  !     * design variables                                               *
  !     ******************************************************************

  use ADjointVars ! include costFunctions
  use ADjointPETSc
  use BCTypes          ! i/j/k/max/min
  use blockPointers    ! il,jl,kl,nViscBocos,nInvBocos
  use flowVarRefState  ! nw,LRef, pInf, gammaInf
  use inputPhysics     ! pointRef, equations, RANSEquations, 
  use inputTimeSpectral 
  use communication    !myID
  use costFunctions
  use section          !sections
  use monitor          !TimeUnsteady
  implicit none
  !
  ! Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: costFunction, npts,nTS
  real(kind=realType), intent(in) :: pts(3,npts,nTS)
  real(kind=realType) :: ptsb(3,npts,nTS)
  logical, intent(in) :: usedjdx,usedjdw
  ! Variables for computeforceandmomentadj_b
  real(kind=realtype) :: force(3), cforce(3)
  real(kind=realtype) :: forceb(3), cforceb(3)
  real(kind=realtype) :: lift, drag, cl, cd
  real(kind=realtype) :: liftb, dragb, clb, cdb
  real(kind=realtype) :: moment(3), cmoment(3)
  real(kind=realtype) :: momentb(3), cmomentb(3)

  real(kind=realType) :: alphaadj,alphaadjb,betaadj,betaadjb
  integer(kind=intType) :: liftindex
  real(kind=realType) :: machcoefadj,machcoefadjb
  real(kind=realType) :: pointrefadj(3),pointrefadjb(3)
  logical :: righthandedadj

  real(kind=realType), dimension(:,:,:,:),allocatable :: wblock,wblockb

  ! Variables for TimeSptectral Derivatives
!!$  real(kind=realType), dimension(nTimeIntervalsSpectral)::       &
!!$       ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,   &
!!$       CmxAdj,CmyAdj,CmzAdj
!!$  real(kind=realType), dimension(nTimeIntervalsSpectral)::       &
!!$       ClAdjb,CdAdjb,CfxAdjb,CfyAdjb,CfzAdjb,   &
!!$       CmxAdjb,CmyAdjb,CmzAdjb
  integer(kind=inttype)::level
  real(kind=realType),dimension(nTimeIntervalsSpectral,8)::BaseCoef
  real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
  real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,dcdbetadot,dcdMach,dcdMachdot
  real(kind=realType),dimension(8)::Coef0,Coef0dot
  real(kind=realType),dimension(nTimeIntervalsSpectral,8)::basecoefb
  real(kind=realType),dimension(8)::dcdpb,dcdpdotb,dcdqb,dcdqdotb,dcdrb,dcdrdotb
  real(kind=realType),dimension(8)::dcdalphab,dcdalphadotb,dcdbetab,dcdbetadotb,dcdMachb,dcdMachdotb
  real(kind=realType),dimension(8)::Coef0b,Coef0dotb
  real(kind=realType), dimension(nCostFunction)::globalCFVals
  !bending derivatives
  real(kind=realType), dimension(nCostFunction)::globalCFValsb
  real(kind=realType)::bendingMoment,bendingMomentb

  real(kind=realType) :: lengthRefAdj,lengthRefAdjb
  real(kind=realType) :: surfaceRefAdj,surfaceRefAdjb
!!$  real(kind=realType) :: cl0,cd0,cmz0,dcldalpha,dcddalpha,dcmzdalpha
!!$  real(kind=realType) :: cl0b,cd0b,cmz0b,dcldalphab,dcddalphab,dcmzdalphab
!!$  real(kind=realType) :: dcmzdqb,dcmzdq
!!$  real(kind=realType) :: dcmzdalphadot,dcmzdalphadotb
  ! Working Variables
  integer(kind=intTYpe) :: sps,ii,iInc,ierr,n
  integer(kind=intTYpe) :: i,j,icell,jcell,kcell,nn,mm,idxmgb,faceID,ibeg,iend,jbeg,jend
  real(kind=realType) :: dIdctemp,val
  real(kind=realType) :: dJdc(nTimeIntervalsSpectral)
  integer(kind=intType) :: row_start,row_end

  !rotation matrix variables
  real(kind=realType),dimension(3):: RpXCorrection,RpYCorrection,RpZCorrection
  real(kind=realType)::rotpointxcorrection,rotpointycorrection,rotpointzcorrection
  !real(kind=realType), dimension(3)   :: RpCorrection, rotPointCorrection
  real(kind=realType), dimension(3)   :: rotationPoint,r
  real(kind=realType), dimension(3,3) :: rotationMatrix  
  real(kind=realType) :: t(nSections),dt(nSections)
  real(kind=realType) :: tOld,tNew

  ! Copy over values we need for the computeforcenadmoment call:
  MachCoefAdj = MachCoef
  pointRefAdj = pointRef

  call getDirAngle(velDirFreestream,LiftDirection,liftIndex,alphaAdj,betaAdj)
  pointRefAdj(1) = pointRef(1)
  pointRefAdj(2) = pointRef(2)
  pointRefAdj(3) = pointRef(3)

  dIda = 0.0
  select case(costFunction)
  case(costFuncLift,costFuncDrag, &
       costFuncLiftCoef,costFuncDragCoef, &
       costFuncForceX,costFuncForceY,costFuncForceZ, &
       costFuncForceXCoef,costFuncForceYCoef,costFuncForceZCoef, &
       costFuncMomX,costFuncMomY,costFuncMomZ,&
       costFuncMomXCoef,costFuncMomYCoef,costFuncMomZCoef,&
       costFuncBendingCoef)

     ! For non-timeSpectral type functions, just time average the
     ! objectives

     dJdc(:) = 1.0/nTimeIntervalsSpectral
     !dJdc(:) = 0.0!1.0/nTimeIntervalsSpectral
     !dJdc(1) = 1.0!/nTimeIntervalsSpectral
     
  case(costFuncCl0,costFuncCd0,costFuncCm0, &
       costFuncClAlpha,costFuncCdAlpha,costFuncCmzAlpha,&
       costFuncClAlphaDot,costFuncCdAlphaDot,costFuncCmzAlphaDot,&
       costFuncClq,costFuncCdq,costFuncCmzq,&
       costFuncClqDot,costFuncCdqDot,costFuncCmzqDot)
     
     ! We have stability derivative cost functions, so there is a more
     ! complex dependance of J on the values computed at each time instance

     ! Get the (sumed) solution values by running getSolution

     do sps =1,nTimeIntervalsSpectral
     
        level = 1
        call computeAeroCoef(globalCFVals,sps)
        
        BaseCoef(sps,1) = globalCFVals(costFuncLiftCoef)
        BaseCoef(sps,2) = globalCFVals(costFuncDragCoef)
        BaseCoef(sps,3) = globalCFVals(costFuncForceXCoef)
        BaseCoef(sps,4) = globalCFVals(costFuncForceYCoef)
        BaseCoef(sps,5) = globalCFVals(costFuncForceZCoef)
        BaseCoef(sps,6) = globalCFVals(costFuncMomXCoef)
        BaseCoef(sps,7) = globalCFVals(costFuncMomYCoef)
        BaseCoef(sps,8) = globalCFVals(costFuncMomZCoef)
        
     end do

     lengthRefAdj = lengthRef

     ! Set Reverse Mode Seeds
     coef0b= 0.0
     dcdalphab= 0.0
     dcdalphadotb= 0.0
     dcdqb= 0.0
     dcdqdotb= 0.0
     lengthrefadjb = 0.0

     select case(costFunction)
     case(costfunccl0)
        coef0b(1)=1.0
     case(costfuncclalpha)
        dcdalphab(1) = 1.0
     case(costfunccd0)
        coef0b(2)=1.0
     case(costfunccdalpha)
        dcdalphab(2) = 1.0
     case(costfunccm0)
        coef0b(8)=1.0
     case(costfunccmzalpha)
        dcdalphab(8) = 1.0
     case(costfunccmzalphadot)
        dcdalphadotb(8) =1.0
     case(costfuncclq)
        dcdqb(1) = 1.0
     case(costfunccmzq)
        dcdqb(8) = 1.0 
     end select

     call COMPUTETSSTABILITYDERIVADJ_B(basecoef, basecoefb, coef0, &
&  coef0b, dcdalpha, dcdalphab, dcdalphadot, dcdalphadotb, dcdq, dcdqb, &
&  dcdqdot, dcdqdotb, lengthrefadj, lengthrefadjb)


     do sps = 1,nTimeIntervalsSpectral
        select case(costFunction)
        case(costfunccl0,costfuncclalpha, costFuncClAlphaDot,costFuncClq,costFuncClqDot)
           dIdctemp = basecoefb(sps,1)!Cladjb(sps)
        case(costfunccd0,costfunccdalpha,costFuncCdAlphaDot,costFuncCdq,costFuncCdqDot)
           dIdctemp = basecoefb(sps,2)!Cdadjb(sps)
        case(costfunccm0,costfunccmzalpha,costfunccmzalphadot,costfunccmzq,costFuncCmzqDot)
           dIdctemp = basecoefb(sps,8)!cmzAdjb(sps)
        end select

        dJdc(sps) = dIdctemp
       
     end do
     if (nDesignLengthRef >=0) then
        !Because the above calculation is based on globally reduced coef., 
        !we only need to store the derivative on the root process,
        !otherwise we end up with nProc times the derivative
        if (myID==0) then

           !print *,'lengthref',lengthRefAdjb,dIda(nDesignLengthRef+1),myid
           dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + lengthRefAdjb!*dJdc(sps)
        end if
     end if
  end select

  ! Now we have dJdc on each processor...when we go through the
  ! reverse mode AD we can take the dot-products on the fly SUM the
  ! entries into dJdw
  if (usedJdw) then
     call VecZeroEntries(dJdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if
  
  if (usedJdx) then
     call VecZeroEntries(dJdx,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecGetOwnershipRange(dJdx,row_start,row_end,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  !for correction to rotPoint derivatives
  do nn=1,nSections
     dt(nn) = sections(nn)%timePeriod &
          / real(nTimeIntervalsSpectral,realType)
  enddo
  
  timeUnsteady = zero


  spectralLoopAdj: do sps=1,nTimeIntervalsSpectral

     do nn=1,nSections
        t(nn) = (sps-1)*dt(nn)
     enddo
     
     ! Compute the displacements due to the rigid motion of the mesh.
     
     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)
     
     call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
     !r = (/-1,-1,-1/)
     !RpCorrection = matmul(rotationMatrix,r)
     r = (/-1,0,0/)
     RpXCorrection = matmul(rotationMatrix,r)
     r = (/0,-1,0/)
     RpYCorrection = matmul(rotationMatrix,r)
     r = (/0,0,-1/)
     RpZCorrection = matmul(rotationMatrix,r)
     ii = 0 

     !zero out the pointrefb value in the module
     pointRefb(:) = 0.0
     lengthRefb = 0.0
     if (costfunction==costFuncBendingCoef)then
        level = 1
        call computeAeroCoef(globalCFVals,sps)
        bendingmomentb = 1.0
        call COMPUTEROOTBENDINGMOMENT_B(globalCFVals, globalCFValsb, bendingmoment, &
             &  bendingmomentb)

        if (nDesignPointRefX >=0) then
           if (myID==0)then
              dIda(nDesignPointRefX + 1) = dIda(nDesignPointRefX + 1) + pointrefb(1)*dJdc(sps)
           end if
        end if
        
        if (nDesignPointRefY >=0) then
           if (myID==0)then
              dIda(nDesignPointRefY + 1) = dIda(nDesignPointRefY + 1) +pointrefb(2)*dJdc(sps)
           end if
        end if
        
        if (nDesignPointRefZ >=0) then
           if (myID==0)then
              dIda(nDesignPointRefZ + 1) = dIda(nDesignPointRefZ + 1) +pointrefb(3)*dJdc(sps)
           end if
        end if
        if (nDesignLengthRef >=0) then
           if (myID==0)then
              dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + lengthRefb*dJdc(sps)
           end if
        end if
     endif

     domainLoopAD: do nn=1,nDom
        call setPointersadj(nn,1_intType,sps)
        bocos: do mm=1,nBocos
           rotpointxcorrection = 0.0
           rotpointycorrection = 0.0
           rotpointzcorrection = 0.0
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
                BCType(mm) == NSWallIsothermal) then

              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

              lengthRefAdj = lengthRef
              SurfaceRefAdj = SurfaceRef
              ! Zero all the backward-mode seeds
              forceb = 0.0
              cforceb = 0.0
              liftb = 0.0
              dragb = 0.0
              clb = 0.0
              cdb = 0.0
              momentb = 0.0
              cmomentb = 0.0
              alphaadjb = 0.0
              betaadjb = 0.0
              machcoefadjb = 0.0
              pointrefadjb = 0.0
              ptsb = 0.0
              lengthRefAdjb= 0.0
              SurfaceRefAdjb = 0.0

              ! Seed the correct value based on the cost function
              select case (costFunction)

              case (costFuncLift)
                 liftb = 1.0
              case (costFuncDrag)
                 dragb = 1.0
              case (costFuncLiftCoef,costFuncCL0,costFuncCLalpha,costFuncClAlphaDot,costFuncClq,costFuncClqDot)
                 clb = 1.0
              case (costFuncDragCoef,costFuncCd0,costFuncCdAlpha,costFuncCdAlphaDot,costFuncCdq,costFuncCdqDot)
                 cdb = 1.0
              case (costFuncForceX)
                 forceb(1) = 1.0
              case (costFuncForceY)
                 forceb(2) = 1.0
              case (costFuncForceZ)
                 forceb(3) = 1.0
              case (costFuncForceXCoef)
                 cforceb(1) = 1.0
              case (costFuncForceYCoef)
                 cforceb(2) = 1.0
              case (costFuncForceZCoef)
                 cforceb(3) = 1.0
              case (costFuncMomX)
                 momentb(1) = 1.0
              case (costFuncMomY)
                 momentb(2) = 1.0
              case (costFuncMomZ)
                 momentb(3) = 1.0
              case (costFuncMomXCoef)
                 cmomentb(1) = 1.0
              case (costFuncMomYCoef)
                 cmomentb(2) = 1.0
              case (costFuncMomZCoef,costFuncCm0,costFuncCMzAlpha,costFuncCmzalphadot,&
                   costFuncCmzq,costFuncCmzqDot)
                
                 cmomentb(3) = 1.0
              case(costFuncBendingCoef)
                 cforceb(1) = globalCFValsb(costFuncForceXCoef)
                 cforceb(2) = globalCFValsb(costFuncForceYCoef)
                 cforceb(3) = globalCFValsb(costFuncForceZCoef)
                 cmomentb(1) = globalCFValsb(costFuncMomXCoef)
                 cmomentb(2) = globalCFValsb(costFuncMomYCoef)
                 cmomentb(3) = globalCFValsb(costFuncMomZCoef)
              end select

              allocate(wblock(0:ib,0:jb,0:kb,nw),&
                   wblockb(0:ib,0:jb,0:kb,nw),stat=ierr)
    
              wblock(:,:,:,:) = w(:,:,:,:)
              wblockb(:,:,:,:) = 0.0
              righthandedadj = righthanded

              faceID = bcfaceid(mm)
 
              call COMPUTEFORCEANDMOMENTADJ_B(force, forceb, cforce, cforceb, &
                   &  lift, liftb, drag, dragb, cl, clb, cd, cdb, moment, momentb, cmoment&
                   &  , cmomentb, alphaadj, alphaadjb, betaadj, betaadjb, liftindex, &
                   &  machcoefadj, machcoefadjb, pointrefadj, pointrefadjb, lengthrefadj, &
                   &  lengthrefadjb, surfacerefadj, surfacerefadjb, pts(:,:,sps), ptsb(:,:,sps), npts, wblock&
                   &  , wblockb, righthandedadj, faceid, ibeg, iend, jbeg, jend, ii, &
                   &  sps)
              
              ! Set the w-values derivatives in dJdw
              if (usedJdw) then
                 do kcell = 2,kl
                    do jcell = 2,jl
                       do icell = 2,il
                          idxmgb = globalCell(icell,jcell,kcell)
                          call VecSetValuesBlocked(dJdw,1,idxmgb,&
                               wblockb(icell,jcell,kcell,:)*dJdc(sps),&
                               ADD_VALUES,PETScIerr)
                          call EChk(PETScIerr,__FILE__,__LINE__)
                       enddo
                    enddo
                 enddo
              end if
              
              if (usedJdx) then
                 ! Set the pt derivative values in dIdpt
                 do j=jBeg,jEnd
                    do i=iBeg,iEnd
                       ! This takes care of the ii increments -- 
                       ! DO NOT NEED INCREMENT ON LINE BELOW
                       ii = ii + 1
                       call VecSetValues(dJdx,3,&
                            
                            (/row_start+3*ii-3,row_start+3*ii-2,row_start+3*ii-1/)+(sps-1)*npts*3,&
                            ptsb(:,ii,sps)*dJdc(sps),ADD_VALUES,PETScIerr)
                       
                       rotpointxcorrection = rotpointxcorrection+DOT_PRODUCT((ptsb(:,ii,sps)*dJdc(sps)),((/1,0,0/)+RpXCorrection))
                       rotpointycorrection = rotpointycorrection+DOT_PRODUCT((ptsb(:,ii,sps)*dJdc(sps)),((/0,1,0/)+RpYCorrection))
                       rotpointzcorrection = rotpointzcorrection+DOT_PRODUCT((ptsb(:,ii,sps)*dJdc(sps)),((/0,0,1/)+RpZCorrection))
                       
                       call EChk(PETScIerr,__file__,__line__)
                    end do
                 end do
              end if
              !ii = ii + (iEnd-iBeg+1)*(jEnd-jBeg+1)

              ! We also have the derivative of the Objective wrt the
              ! "AeroDVs" intrinsic aero design variables, alpha, beta etc

              if (nDesignAoA >=0) then
                 dIda(nDesignAoA+1) = dIda(nDesignAoA+1) + alphaAdjb*dJdc(sps)

              end if

              if (nDesignSSA >= 0) then
                 dIda(nDesignSSA+1) = dIda(nDesignSSA+1) + betaAdjb*dJdc(sps)
              end if

              if (nDesignMach >= 0) then
                 dIda(nDesignMach+1) = dIda(nDesignMach+1) + machCoefAdjb*dJdc(sps)
              end if

              if (nDesignMachGrid >= 0) then
                 dIda(nDesignMachGrid+1) = dIda(nDesignMachGrid+1) + machCoefAdjb*dJdc(sps)
              end if
              
              if (nDesignPointRefX >=0) then
                 
                 dIda(nDesignPointRefX + 1) = dIda(nDesignPointRefX + 1) + pointrefAdjb(1)*dJdc(sps)
                
              end if

              if (nDesignPointRefY >=0) then
                 dIda(nDesignPointRefY + 1) = dIda(nDesignPointRefY + 1) + pointrefAdjb(2)*dJdc(sps)
              end if

              if (nDesignPointRefZ >=0) then
                 dIda(nDesignPointRefZ + 1) = dIda(nDesignPointRefZ + 1) + pointrefAdjb(3)*dJdc(sps)
              end if

              if (nDesignRotCenX >= 0) then
                 dIda(nDesignRotCenX+1) = dIda(nDesignRotCenX+1)+rotpointxcorrection
              endif
              
              if (nDesignRotCenY >= 0) then
                  dIda(nDesignRotCenY+1) = dIda(nDesignRotCenY+1)+rotpointycorrection
              end if
              if (nDesignRotCenZ >= 0) then
                  dIda(nDesignRotCenZ+1) = dIda(nDesignRotCenZ+1)+rotpointzcorrection
              endif

              if (nDesignLengthRef >=0) then

                 dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + lengthRefAdjb*dJdc(sps)
                 
              end if
              if (nDesignSurfaceRef >=0) then
                 dIda(nDesignSurfaceRef+1) = dIda(nDesignSurfaceRef+1) + SurfaceRefAdjb*dJdc(sps)
              end if
              if (nDesignDissError >=0) then
                 dIda(nDesignDissError+1) = 0
              end if
                  
           end if
           deallocate(wblock,wblockb,stat=ierr)

        end do bocos
     end do domainLoopAD
  end do spectralLoopAdj
  
  ! Assemble the petsc vectors
  if (usedJdw) then
     call VecAssemblyBegin(dJdw,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     call VecAssemblyEnd(dJdw,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  end if

  if (usedJdx) then
     call VecAssemblyBegin(dJdx,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     call VecAssemblyEnd(dJdx,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  end if
#endif
end subroutine computeObjPartials


! Add two functions to return dIdw and dIdx. dIda is available
! directly in python in dIda in the adjointVars module. 

subroutine getdIdw(ndof,output)
#ifndef USE_NO_PETSC
  use ADjointPETSc
  use ADjointVars
  use precision 
  implicit none

  integer(kind=intType),intent(in) :: ndof
  real(kind=realType),intent(out)  :: output(ndof)

  integer(kind=intType) :: ilow,ihigh,i

  call VecGetOwnershipRange(dJdw,ilow,ihigh,PETScIerr)

  do i=1,(ihigh-ilow)
     call VecGetValues(dJdw,1,ilow+i-1,output(i),PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  end do
#endif
end subroutine getdIdw

subroutine getdIdx(ndof,output)
#ifndef USE_NO_PETSC
  use ADjointPETSc
  use ADjointVars
  use inputTimeSpectral
  use section
  use monitor

  implicit none

  integer(kind=intType),intent(in) :: ndof
  real(kind=realType),intent(out)  :: output(ndof)
  integer(kind=intType),dimension(3) :: idx
  real(kind=realType),dimension(3) ::temp
  !rotation matrix variables
  real(kind=realType), dimension(3)   :: rotationPoint,r
  real(kind=realType), dimension(3,3) :: rotationMatrix  
  real(kind=realType) :: t(nSections),dt(nSections)
  real(kind=realType) :: tOld,tNew
  integer(kind=intType) :: ilow,ihigh,i,sps,nn
  output(:) = 0.0

  call VecGetOwnershipRange(dJdx,ilow,ihigh,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Comput rotation from each time instance back to a single surface
  
  do nn=1,nSections
     dt(nn) = sections(nn)%timePeriod &
          / real(nTimeIntervalsSpectral,realType)
  enddo
  
  timeUnsteady = zero
  
  do sps = 1,nTimeIntervalsSpectral
     do nn=1,nSections
        t(nn) = (sps-1)*dt(nn)
     enddo
     
     ! Compute the displacements due to the rigid motion of the mesh.
     
     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)
     
     call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

     ! Take rotation Matrix Transpose
     rotationMatrix = transpose(rotationMatrix)

     do i=1,ndof/3

        idx = (/ilow+ndof*(sps-1)+(3*(i-1)),&
                ilow+ndof*(sps-1)+(3*(i-1))+1,&
                ilow+ndof*(sps-1)+(3*(i-1))+2/)
        call VecGetValues(dJdx,3,idx,temp,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)

        output((i-1)*3+1:(i-1)*3+3) = output((i-1)*3+1:(i-1)*3+3)+ &
             matmul(rotationMatrix,temp)

     end do
  end do
#endif
end subroutine getdIdx

subroutine zeroObjPartials
#ifndef USE_NO_PETSC
  use precision 
  use ADjointVars ! include costFunctions
  use ADjointPETSc

  integer(kind=intType) :: ierr

  call VecZeroEntries(dJdw,ierr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecZeroEntries(dJdx,ierr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif
end subroutine zeroObjPartials

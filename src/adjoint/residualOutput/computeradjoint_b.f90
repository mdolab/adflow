   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
   !  
   !  Differentiation of computeradjoint in reverse (adjoint) mode:
   !   gradient, with respect to input variables: pointrefadj rotrateadj
   !                machadj alphaadj rotpointadj xadj xblockcorneradj
   !                dwadj wadj betaadj machgridadj rotcenteradj
   !   of linear combination of output variables: dwadj
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          computeRAdj.f90                                 *
   !      * Author:        C.A.(Sandy) Mader                               *
   !      * Starting date: 02-01-2008                                      *
   !      * Last modified: 04-23-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
   &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
   &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
   &  icell, jcell, kcell, nn, level, sps, correctfork, secondhalo, prefadj&
   &  , rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj&
   &  , rotrateadjb, rotcenteradj, rotcenteradjb, pointrefadj, pointrefadjb&
   &  , rotpointadj, rotpointadjb, murefadj, timerefadj, pinfcorradj, &
   &  liftindex)
   USE blockpointers
   USE flowvarrefstate
   USE inputphysics
   USE inputtimespectral
   USE monitor
   USE section
   IMPLICIT NONE
   REAL(KIND=REALTYPE) :: alphaadj, alphaadjb, betaadj, betaadjb
   LOGICAL :: correctfork, secondhalo
   REAL(KIND=REALTYPE) :: dwadj(nw, ntimeintervalsspectral), dwadjb(nw, &
   &  ntimeintervalsspectral)
   INTEGER(KIND=INTTYPE), INTENT(IN) :: icell
   INTEGER(KIND=INTTYPE), INTENT(IN) :: jcell
   INTEGER(KIND=INTTYPE), INTENT(IN) :: kcell
   INTEGER(KIND=INTTYPE), INTENT(IN) :: level
   INTEGER(KIND=INTTYPE) :: liftindex
   REAL(KIND=REALTYPE) :: machadj, machadjb, machcoefadj, machgridadj, &
   &  machgridadjb, pinfcorradj
   REAL(KIND=REALTYPE) :: murefadj, timerefadj
   INTEGER(KIND=INTTYPE), INTENT(IN) :: nn
   REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
   REAL(KIND=REALTYPE) :: pinfadj, rhoinfadj
   REAL(KIND=REALTYPE) :: prefadj, rhorefadj
   REAL(KIND=REALTYPE), DIMENSION(3), INTENT(IN) :: rotcenteradj
   REAL(KIND=REALTYPE) :: rotcenteradjb(3), rotrateadjb(3)
   REAL(KIND=REALTYPE) :: pointrefadj(3), pointrefadjb(3), rotpointadj(3)&
   &  , rotpointadjb(3)
   REAL(KIND=REALTYPE), DIMENSION(3), INTENT(IN) :: rotrateadj
   INTEGER(KIND=INTTYPE), INTENT(IN) :: sps
   REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral), INTENT(IN) :: wadj
   REAL(KIND=REALTYPE) :: wadjb(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE), DIMENSION(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: xadj
   REAL(KIND=REALTYPE) :: xadjb(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: xblockcorneradj(2, 2, 2, 3, &
   &  ntimeintervalsspectral), xblockcorneradjb(2, 2, 2, 3, &
   &  ntimeintervalsspectral)
   INTEGER :: branch
   REAL(KIND=REALTYPE) :: dragdirectionadj(3)
   INTEGER(KIND=INTTYPE) :: i, ii, j, jj, k, kk, nnn, sps2
   INTEGER(KIND=INTTYPE) :: iend, istart, jend, jstart, kend, kstart
   REAL(KIND=REALTYPE) :: liftdirectionadj(3)
   REAL(KIND=REALTYPE) :: normadj(nbocos, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral), normadjb(nbocos, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: padj(-2:2, -2:2, -2:2, ntimeintervalsspectral)&
   &  , padjb(-2:2, -2:2, -2:2, ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: radiadj(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral), radiadjb(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral), radjadj(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral), radjadjb(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral), radkadj(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral), radkadjb(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: rfaceadj(nbocos, -2:2, -2:2, &
   &  ntimeintervalsspectral), rfaceadjb(nbocos, -2:2, -2:2, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: sadj(-2:2, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral), sadjb(-2:2, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: sfaceiadj(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral), sfaceiadjb(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral), sfacejadj(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral), sfacejadjb(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral), sfacekadj(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral), sfacekadjb(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: siadj(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), siadjb(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), sjadj(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), sjadjb(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), skadj(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), skadjb(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: t(nsections)
   REAL(KIND=REALTYPE) :: tempb(nw)
   REAL(KIND=REALTYPE) :: pinfcorradjb, uinfadj, uinfadjb
   LOGICAL :: useoldcoor=.false.
   REAL(KIND=REALTYPE) :: veldirfreestreamadj(3), veldirfreestreamadjb(3)
   REAL(KIND=REALTYPE) :: voladj(ntimeintervalsspectral), voladjb(&
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE) :: winfadj(nw), winfadjb(nw)
   !      Set Use Modules
   !nTimeIntervalsSpectral
   !nsection
   !      Set Passed in Variables
   !      Set Local Variables
   !variables for test loops
   !  real(kind=realType), dimension(-2:2,-2:2,-2:2,3) :: siAdj, sjAdj, skAdj
   !!$!File Parameters remove for AD
   !!$      integer :: unitxAD = 15,ierror
   !!$      integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart2,jstart2,kstart2,iend2,jend2,kend2,n
   !!$      character(len = 16)::outfile
   !!$      
   !!$      outfile = "xAD.txt"
   !!$      
   !!$      open (UNIT=unitxAD,File=outfile,status='old',position='append',action='write',iostat=ierror)
   !!$      if(ierror /= 0)                        &
   !!$           call terminate("verifyResiduals", &
   !!$           "Something wrong when &
   !!$           &calling open")
   ! *************************************************************************
   !      Begin Execution
   ! *************************************************************************
   !print *,'in computeRadj',wadj(:,:,:,irho)!
   !      call the initialization routines to calculate the effect of Mach and alpha
   CALL ADJUSTINFLOWANGLEADJ(alphaadj, betaadj, veldirfreestreamadj, &
   &                      liftdirectionadj, dragdirectionadj, liftindex)
   CALL PUSHREAL8ARRAY(veldirfreestreamadj, 3)
   CALL CHECKINPUTPARAMADJ(veldirfreestreamadj, liftdirectionadj, &
   &                    dragdirectionadj, machadj, machcoefadj)
   CALL PUSHREAL8(gammainf)
   CALL PUSHREAL8(rhorefadj)
   CALL PUSHREAL8(prefadj)
   CALL REFERENCESTATEADJ(machadj, machcoefadj, uinfadj, prefadj, &
   &                   rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, &
   &                   pinfadj, murefadj, timerefadj)
   CALL PUSHREAL8ARRAY(winfadj, nw)
   !call referenceStateAdj(velDirFreestreamAdj,liftDirectionAdj,&
   !     dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj,prefAdj,&
   !     rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
   !     murefAdj, timerefAdj)
   !(velDirFreestreamAdj,liftDirectionAdj,&
   !     dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj)
   CALL SETFLOWINFINITYSTATEADJ(veldirfreestreamadj, liftdirectionadj, &
   &                         dragdirectionadj, machadj, machcoefadj, &
   &                         uinfadj, winfadj, prefadj, rhorefadj, &
   &                         pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, &
   &                         murefadj, timerefadj, pinfcorradj)
   ! sadj = 0.0
   DO sps2=1,ntimeintervalsspectral
   CALL PUSHREAL8ARRAY(xadj, 6**3*3*ntimeintervalsspectral)
   !      Call the metric routines to generate the areas, volumes and surface normals for the stencil.
   CALL XHALOADJ(xadj, xblockcorneradj, icell, jcell, kcell, nn, level&
   &            , sps, sps2)
   CALL PUSHREAL8ARRAY(normadj, nbocos*5**2*3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(voladj, ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(skadj, 6**3*3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(sjadj, 6**3*3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(siadj, 6**3*3*ntimeintervalsspectral)
   CALL METRICADJ(xadj, siadj, sjadj, skadj, voladj, normadj, icell, &
   &             jcell, kcell, nn, level, sps, sps2)
   CALL PUSHREAL8ARRAY(t, nsections)
   !call the gridVelocities function to get the cell center ,face center and boundary mesh velocities.
   ! Compute the time, which corresponds to this spectral solution.
   ! For steady and unsteady mode this is simply the restart time;
   ! for the spectral mode the periodic time must be taken into
   ! account, which can be different for every section.
   t = timeunsteadyrestart
   IF (equationmode .EQ. timespectral) THEN
   DO nnn=1,nsections
   !t(nnn) = t(nnn) + (sps2-1)*sections(nnn)%timePeriod &
   !     /         real(nTimeIntervalsSpectral,realType)
   !to make denomenator a real number...
   t(nnn) = t(nnn) + (sps2-1)*sections(nnn)%timeperiod/(&
   &          ntimeintervalsspectral*1.0)
   END DO
   CALL PUSHINTEGER4(1)
   ELSE
   CALL PUSHINTEGER4(0)
   END IF
   CALL PUSHREAL8ARRAY(sfacekadj, 5**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(sfacejadj, 5**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(sfaceiadj, 5**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(sadj, 5**3*3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(rotrateadj, 3)
   !first two arguments needed for time spectral.just set to initial values for the current steady case...
   CALL GRIDVELOCITIESFINELEVELADJ(useoldcoor, t, sps, xadj, siadj, &
   &                              sjadj, skadj, rotcenteradj, rotrateadj, &
   &                              sadj, sfaceiadj, sfacejadj, sfacekadj, &
   &                              machgridadj, veldirfreestreamadj, &
   &                              liftdirectionadj, alphaadj, betaadj, &
   &                              liftindex, icell, jcell, kcell, &
   &                              pointrefadj, rotpointadj, nn, level, sps2&
   &                             )
   CALL PUSHREAL8ARRAY(rfaceadj, nbocos*5**2*ntimeintervalsspectral)
   CALL NORMALVELOCITIESALLLEVELSADJ(sps, icell, jcell, kcell, &
   &                                sfaceiadj, sfacejadj, sfacekadj, siadj&
   &                                , sjadj, skadj, rfaceadj, nn, level, &
   &                                sps2)
   CALL PUSHREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   !needed for uSlip in Viscous Calculations
   !call slipVelocitiesFineLevel(.false., t, mm)
   !      Mimic the Residual calculation in the main code
   !Compute the Pressure in the stencil based on the current 
   !States
   !print *,'Calling computepressure',il,jl,kl,nn,secondhalo!,wadj(:,:,:,irho)!
   ! replace with Compute Pressure Adjoint!
   CALL COMPUTEPRESSUREADJ(wadj, padj, nn, level, sps, sps2)
   CALL PUSHBOOLEAN(secondhalo)
   CALL PUSHREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   !!$       !print out pAdj
   !!$       istart2 = -2
   !!$       jstart2 = -2
   !!$       kstart2 = -2
   !!$       iend2 = 2
   !!$       jend2 = 2
   !!$       kend2 = 2 
   !!$       if(icell==2) istart2=-1
   !!$       if(jcell==2) jstart2=-1
   !!$       if(kcell==2) kstart2=-1
   !!$       if(icell==il) iend2=1
   !!$       if(jcell==jl) jend2=1
   !!$       if(kcell==kl) kend2=1
   !!$       do iiii = istart2,iend2
   !!$          do jjjj = jstart2,jend2
   !!$             do kkkk = kstart2,kend2
   !!$                !do n = 1,3!nw
   !!$                   do n = 1,1!nw
   !!$                   !do n = 1,nw 
   !!$                   !do sps2 = 1,nTimeIntervalsSpectral
   !!$                      i = icell+iiii
   !!$                      j = jcell+jjjj
   !!$                      k = kcell+kkkk
   !!$                                 !write(unitxAD,11) i,j,k,n,nnn,sps,sps2,wAdj(iiii,jjjj,kkkk,n,sps2) 
   !!$                      !write(unitxAD,11) i,j,k,n,nn,sps,sps2,pAdj(iiii,jjjj,kkkk,sps2)
   !!$                      write(unitxAD,11) i,j,k,n,nn,sps,sps2,sFaceIAdj(iiii,jjjj,kkkk,sps2)
   !!$                      !write(unitxAD,11) i,j,k,n,nn,sps,sps2,sAdj(iiii,jjjj,kkkk,n,sps2)
   !!$
   !!$11                    format(1x,'wadj',7I8,f20.14)
   !!$                   !enddo
   !!$                enddo
   !!$             enddo
   !!$          enddo
   !!$       enddo
   ! Apply all boundary conditions to stencil.
   ! In case of a full mg mode, and a segegated turbulent solver,
   ! first call the turbulent boundary conditions, such that the
   ! turbulent kinetic energy is properly initialized in the halo's.
   !###! Ignore Viscous for now
   !###!       if(turbSegregated .and. (.not. corrections)) &
   !###!         call applyAllTurbBCAdj(secondHalo)
   ! Apply all boundary conditions of the mean flow.
   !print *,'applying bcs',nn,secondhalo,sps2
   !******************************************
   CALL APPLYALLBCADJ(winfadj, pinfcorradj, wadj, padj, sadj, siadj, &
   &                 sjadj, skadj, voladj, normadj, rfaceadj, icell, jcell&
   &                 , kcell, secondhalo, nn, level, sps, sps2)
   CALL PUSHREAL8ARRAY(radkadj, 3**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(radjadj, 3**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(radiadj, 3**3*ntimeintervalsspectral)
   !print *,'bcsiAdj',siAdj(:,0,0,1,:)
   !#!#Shouldn't need this section for derivatives...
   !#!$       ! In case this routine is called in full mg mode call the mean
   !#!$       ! flow boundary conditions again such that the normal momentum
   !#!$       ! boundary condition is treated correctly.
   !#!$
   !#!$       if(.not. corrections) call applyAllBCAdj(wAdj, pAdj, &
   !#!$                              siAdj, sjAdj, skAdj, volAdj, normAdj, &
   !#!$                              iCell, jCell, kCell,secondHalo)
   !Leave out State exchanges for now. If there are discrepancies 
   !Later, this may be a source...
   !#!$       ! Exchange the solution. Either whalo1 or whalo2
   !#!$       ! must be called.
   !#!$
   !#!$       if( secondHalo ) then
   !#!$         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
   !#!$                     .true., .true.)
   !#!$       else
   !#!$         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
   !#!$                     .true., .true.)
   !#!$       endif
   !Again this should not be required, so leave out for now...
   ! For full multigrid mode the bleeds must be determined, the
   ! boundary conditions must be applied one more time and the
   ! solution must be exchanged again.
   !#!$       if(.not. corrections) then
   !#!$         call BCDataMassBleedOutflowAdj(.true., .true.)
   !#!$         call applyAllBCAdj(secondHalo)
   !#!$
   !#!$       !Leave out State exchanges for now. If there are discrepancies 
   !#!$       !Later, this may be a source...
   !#!$!         if( secondHalo ) then
   !#!$!           call whalo2(currentLevel, 1_intType, nVarInt, .true., &
   !#!$!                       .true., .true.)
   !#!$!         else
   !#!!$           call whalo1(currentLevel, 1_intType, nVarInt, .true., &
   !#!!$                       .true., .true.)
   !#!!$         endif
   !#!$       endif
   !#!$
   !#!$
   !#!$
   !#!$       ! Reset the values of rkStage and currentLevel, such that
   !#!$       ! they correspond to a new iteration.
   !#!$
   !#!$       rkStage = 0
   !#!$       currentLevel = groundLevel
   !#!$
   !#!$       ! Compute the latest values of the skin friction velocity.
   !#!$       ! The currently stored values are of the previous iteration.
   !#!$
   !#!$       call computeUtauAdj
   !#!$
   !#!$       ! Apply an iteration to the turbulent transport equations in
   !#!$       ! case these must be solved segregatedly.
   !#!$
   !#!$       if( turbSegregated ) call turbSolveSegregatedAdj
   !#!$
   ! Compute the time step.
   !call timeStepAdj(.false.)
   !print *,'nntimestep',nn
   CALL TIMESTEPADJ(.true., wadj, padj, siadj, sjadj, skadj, sfaceiadj&
   &               , sfacejadj, sfacekadj, voladj, radiadj, radjadj, &
   &               radkadj, icell, jcell, kcell, pinfcorradj, rhoinfadj, nn&
   &               , level, sps, sps2)
   END DO
   !print *,'tstepsiAdj',siAdj(:,0,0,1,:)
   !#!$
   !#!$       ! Compute the residual of the new solution on the ground level.
   !#!$
   !#!$       if( turbCoupled ) then
   !#!$         call initresAdj(nt1MG, nMGVar)
   !#!$         call turbResidualAdj
   !#!$       endif
   !#!$
   !!$   dwadj(:,sps) = wadj(0,0,0,:,sps)*voladj(sps)
   !!$     !do sps2 = 1,nTimeIntervalsSpectral
   !!$     !  print *,'calculating initres',nn
   !!$       !call initresAdj(1_intType, nwf,sps,dwAdj)
   !    dwAdj(:,sps) = 0.0
   !    dwAdj(1:3,sps) = sAdj(0,0,0,:,sps)!xAdj(0,0,0,:,sps)
   !dwAdj(4,sps) = volAdj(sps)
   CALL INITRESADJ(1, nwf, wadj, voladj, dwadj, nn, level, sps)
   CALL PUSHREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   !print *,'dwadj',dwadj,icell,jcell,kcell
   !  print *,'calculating residuals',nn
   CALL RESIDUALADJ(wadj, padj, siadj, sjadj, skadj, voladj, normadj, &
   &             sfaceiadj, sfacejadj, sfacekadj, radiadj, radjadj, radkadj&
   &             , dwadj, icell, jcell, kcell, rotrateadj, correctfork, nn&
   &             , level, sps)
   !end do
   ! print *,'nn end',nn
   !stop
   !close (UNIT=unitxAD)
   DO sps2=1,ntimeintervalsspectral
   CALL PUSHREAL8ARRAY(dwadj(:, sps2), nw)
   dwadj(:, sps2) = dwadj(:, sps2)/voladj(sps2)
   END DO
   voladjb(1:ntimeintervalsspectral) = 0.0
   DO sps2=ntimeintervalsspectral,1,-1
   CALL POPREAL8ARRAY(dwadj(:, sps2), nw)
   tempb = dwadjb(:, sps2)/voladj(sps2)
   voladjb(sps2) = voladjb(sps2) + SUM(-(dwadj(:, sps2)*tempb/voladj(&
   &      sps2)))
   dwadjb(:, sps2) = tempb
   END DO
   CALL POPREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   CALL RESIDUALADJ_B(wadj, wadjb, padj, padjb, siadj, siadjb, sjadj, &
   &               sjadjb, skadj, skadjb, voladj, voladjb, normadj, &
   &               sfaceiadj, sfaceiadjb, sfacejadj, sfacejadjb, sfacekadj&
   &               , sfacekadjb, radiadj, radiadjb, radjadj, radjadjb, &
   &               radkadj, radkadjb, dwadj, dwadjb, icell, jcell, kcell, &
   &               rotrateadj, rotrateadjb, correctfork, nn, level, sps)
   CALL INITRESADJ_B(1, nwf, wadj, wadjb, voladj, voladjb, dwadj, dwadjb&
   &              , nn, level, sps)
   pointrefadjb(1:3) = 0.0
   alphaadjb = 0.0
   pinfcorradjb = 0.0
   rotpointadjb(1:3) = 0.0
   xadjb(-3:2, -3:2, -3:2, 1:3, 1:ntimeintervalsspectral) = 0.0
   xblockcorneradjb(1:2, 1:2, 1:2, 1:3, 1:ntimeintervalsspectral) = 0.0
   betaadjb = 0.0
   machgridadjb = 0.0
   rotcenteradjb(1:3) = 0.0
   winfadjb(1:nw) = 0.0
   rfaceadjb(1:nbocos, -2:2, -2:2, 1:ntimeintervalsspectral) = 0.0
   sadjb(-2:2, -2:2, -2:2, 1:3, 1:ntimeintervalsspectral) = 0.0
   veldirfreestreamadjb(1:3) = 0.0
   normadjb(1:nbocos, -2:2, -2:2, 1:3, 1:ntimeintervalsspectral) = 0.0
   DO sps2=ntimeintervalsspectral,1,-1
   CALL POPREAL8ARRAY(radiadj, 3**3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(radjadj, 3**3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(radkadj, 3**3*ntimeintervalsspectral)
   CALL TIMESTEPADJ_B(.true., wadj, wadjb, padj, padjb, siadj, siadjb, &
   &                 sjadj, sjadjb, skadj, skadjb, sfaceiadj, sfaceiadjb, &
   &                 sfacejadj, sfacejadjb, sfacekadj, sfacekadjb, voladj, &
   &                 radiadj, radiadjb, radjadj, radjadjb, radkadj, &
   &                 radkadjb, icell, jcell, kcell, pinfcorradj, &
   &                 pinfcorradjb, rhoinfadj, nn, level, sps, sps2)
   CALL POPREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   CALL POPBOOLEAN(secondhalo)
   CALL APPLYALLBCADJ_B(winfadj, winfadjb, pinfcorradj, pinfcorradjb, &
   &                   wadj, wadjb, padj, padjb, sadj, sadjb, siadj, siadjb&
   &                   , sjadj, sjadjb, skadj, skadjb, voladj, normadj, &
   &                   normadjb, rfaceadj, rfaceadjb, icell, jcell, kcell, &
   &                   secondhalo, nn, level, sps, sps2)
   CALL POPREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   CALL COMPUTEPRESSUREADJ_B(wadj, wadjb, padj, padjb, nn, level, sps, &
   &                        sps2)
   CALL POPREAL8ARRAY(rfaceadj, nbocos*5**2*ntimeintervalsspectral)
   CALL NORMALVELOCITIESALLLEVELSADJ_B(sps, icell, jcell, kcell, &
   &                                  sfaceiadj, sfaceiadjb, sfacejadj, &
   &                                  sfacejadjb, sfacekadj, sfacekadjb, &
   &                                  siadj, siadjb, sjadj, sjadjb, skadj, &
   &                                  skadjb, rfaceadj, rfaceadjb, nn, &
   &                                  level, sps2)
   CALL POPREAL8ARRAY(rotrateadj, 3)
   CALL POPREAL8ARRAY(sadj, 5**3*3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(sfaceiadj, 5**3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(sfacejadj, 5**3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(sfacekadj, 5**3*ntimeintervalsspectral)
   CALL GRIDVELOCITIESFINELEVELADJ_B(useoldcoor, t, sps, xadj, xadjb, &
   &                                siadj, siadjb, sjadj, sjadjb, skadj, &
   &                                skadjb, rotcenteradj, rotcenteradjb, &
   &                                rotrateadj, rotrateadjb, sadj, sadjb, &
   &                                sfaceiadj, sfaceiadjb, sfacejadj, &
   &                                sfacejadjb, sfacekadj, sfacekadjb, &
   &                                machgridadj, machgridadjb, &
   &                                veldirfreestreamadj, &
   &                                veldirfreestreamadjb, liftdirectionadj&
   &                                , alphaadj, alphaadjb, betaadj, &
   &                                betaadjb, liftindex, icell, jcell, &
   &                                kcell, pointrefadj, pointrefadjb, &
   &                                rotpointadj, rotpointadjb, nn, level, &
   &                                sps2)
   CALL POPINTEGER4(branch)
   IF (.NOT.branch .LT. 1) nnn = 0
   CALL POPREAL8ARRAY(t, nsections)
   CALL POPREAL8ARRAY(siadj, 6**3*3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(sjadj, 6**3*3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(skadj, 6**3*3*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(voladj, ntimeintervalsspectral)
   CALL POPREAL8ARRAY(normadj, nbocos*5**2*3*ntimeintervalsspectral)
   CALL METRICADJ_B(xadj, xadjb, siadj, siadjb, sjadj, sjadjb, skadj, &
   &               skadjb, voladj, voladjb, normadj, normadjb, icell, jcell&
   &               , kcell, nn, level, sps, sps2)
   CALL POPREAL8ARRAY(xadj, 6**3*3*ntimeintervalsspectral)
   CALL XHALOADJ_B(xadj, xadjb, xblockcorneradj, xblockcorneradjb, &
   &              icell, jcell, kcell, nn, level, sps, sps2)
   END DO
   CALL POPREAL8ARRAY(winfadj, nw)
   CALL SETFLOWINFINITYSTATEADJ_B(veldirfreestreamadj, &
   &                           veldirfreestreamadjb, liftdirectionadj, &
   &                           dragdirectionadj, machadj, machcoefadj, &
   &                           uinfadj, uinfadjb, winfadj, winfadjb, &
   &                           prefadj, rhorefadj, pinfdimadj, rhoinfdimadj&
   &                           , rhoinfadj, pinfadj, murefadj, timerefadj, &
   &                           pinfcorradj, pinfcorradjb)
   CALL POPREAL8(prefadj)
   CALL POPREAL8(rhorefadj)
   CALL POPREAL8(gammainf)
   CALL REFERENCESTATEADJ_B(machadj, machadjb, machcoefadj, uinfadj, &
   &                     uinfadjb, prefadj, rhorefadj, pinfdimadj, &
   &                     rhoinfdimadj, rhoinfadj, pinfadj, murefadj, &
   &                     timerefadj)
   CALL POPREAL8ARRAY(veldirfreestreamadj, 3)
   CALL CHECKINPUTPARAMADJ_B(veldirfreestreamadj, veldirfreestreamadjb, &
   &                      liftdirectionadj, dragdirectionadj, machadj, &
   &                      machcoefadj)
   CALL ADJUSTINFLOWANGLEADJ_B(alphaadj, alphaadjb, betaadj, betaadjb, &
   &                        veldirfreestreamadj, veldirfreestreamadjb, &
   &                        liftdirectionadj, dragdirectionadj, liftindex)
   dwadjb(1:nw, 1:ntimeintervalsspectral) = 0.0
   END SUBROUTINE COMPUTERADJOINT_B

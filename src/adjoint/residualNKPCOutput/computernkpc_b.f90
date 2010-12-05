   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of computernkpc in reverse (adjoint) mode:
   !   gradient     of useful results: dwadj
   !   with respect to varying inputs: dwadj wadj
   !   RW status of diff variables: dwadj:in-out wadj:out
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          computeRAdj.f90                                 *
   !      * Author:        C.A.(Sandy) Mader                               *
   !      * Starting date: 02-01-2008                                      *
   !      * Last modified: 04-23-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   SUBROUTINE COMPUTERNKPC_B(wadj, wadjb, dwadj, dwadjb, siadj, sjadj, &
   &  skadj, sadj, voladj, sfaceiadj, sfacejadj, sfacekadj, rotrateadj, &
   &  icell, jcell, kcell, nn, level, sps)
   USE BLOCKPOINTERS
   USE SECTION
   USE INPUTTIMESPECTRAL
   USE INPUTPHYSICS
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   ! Passed in Variables
   ! Input Variables --- Note no intent(in) --- this can cause problems
   ! with reverse mode AD
   INTEGER(kind=inttype), INTENT(IN) :: icell, jcell, kcell, nn, level, &
   &  sps
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral) :: wadj
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral) :: wadjb
   REAL(kind=realtype), DIMENSION(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral) :: siadj, sjadj, skadj
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral) :: sadj
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral) :: voladj
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral) :: sfaceiadj, sfacejadj, sfacekadj
   REAL(kind=realtype), DIMENSION(3) :: rotrateadj
   ! Ouptut Variables
   REAL(kind=realtype), DIMENSION(nw, ntimeintervalsspectral) :: dwadj
   REAL(kind=realtype), DIMENSION(nw, ntimeintervalsspectral) :: dwadjb
   !      Set Local Variables
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral) :: padj
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, -2:2, &
   &  ntimeintervalsspectral) :: padjb
   REAL(kind=realtype), DIMENSION(nbocos, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral) :: normadj
   REAL(kind=realtype), DIMENSION(nbocos, -2:2, -2:2, &
   &  ntimeintervalsspectral) :: rfaceadj
   REAL(kind=realtype), DIMENSION(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral) :: radiadj, radjadj, radkadj
   REAL(kind=realtype), DIMENSION(-1:1, -1:1, -1:1, &
   &  ntimeintervalsspectral) :: radiadjb, radjadjb, radkadjb
   REAL(kind=realtype), DIMENSION(nsections) :: t
   ! Integer variables
   INTEGER(kind=inttype) :: sps2
   ! Logical 
   LOGICAL :: secondhalo, correctfork
   DO sps2=1,ntimeintervalsspectral
   CALL PUSHREAL8ARRAY(padj, realtype*5**3*ntimeintervalsspectral/8)
   !Compute the Pressure in the stencil based on the current 
   !States
   CALL COMPUTEPRESSURENKPC(wadj, padj, nn, level, sps, sps2)
   CALL PUSHREAL8ARRAY(radkadj, realtype*3**3*ntimeintervalsspectral/8)
   CALL PUSHREAL8ARRAY(radjadj, realtype*3**3*ntimeintervalsspectral/8)
   CALL PUSHREAL8ARRAY(radiadj, realtype*3**3*ntimeintervalsspectral/8)
   ! Apply all boundary conditions to stencil.
   ! In case of a full mg mode, and a segegated turbulent solver,
   ! first call the turbulent boundary conditions, such that the
   ! turbulent kinetic energy is properly initialized in the halo's.
   !      call applyAllBCNKPC(wInf,pInfCorr,wAdj, pAdj,sAdj, &
   !           siAdj, sjAdj, skAdj, volAdj, normAdj, &
   !           rFaceAdj,iCell, jCell, kCell,secondHalo,nn,level,sps,sps2)
   CALL TIMESTEPNKPC(.true., wadj, padj, siadj, sjadj, skadj, &
   &                   sfaceiadj, sfacejadj, sfacekadj, voladj, radiadj, &
   &                   radjadj, radkadj, icell, jcell, kcell, pinfcorr, &
   &                   rhoinf, nn, level, sps, sps2)
   END DO
   DO sps2=ntimeintervalsspectral,1,-1
   dwadjb(:, sps2) = dwadjb(:, sps2)/voladj(sps2)
   END DO
   CALL RESIDUALNKPC_B(wadj, wadjb, padj, padjb, siadj, sjadj, skadj, &
   &                voladj, normadj, sfaceiadj, sfacejadj, sfacekadj, &
   &                radiadj, radiadjb, radjadj, radjadjb, radkadj, radkadjb&
   &                , dwadj, dwadjb, icell, jcell, kcell, rotrateadj, &
   &                correctfork, nn, level, sps)
   CALL INITRESNKPC_B(1, nwf, wadj, wadjb, voladj, dwadj, dwadjb, nn, &
   &               level, sps)
   DO sps2=ntimeintervalsspectral,1,-1
   CALL POPREAL8ARRAY(radiadj, realtype*3**3*ntimeintervalsspectral/8)
   CALL POPREAL8ARRAY(radjadj, realtype*3**3*ntimeintervalsspectral/8)
   CALL POPREAL8ARRAY(radkadj, realtype*3**3*ntimeintervalsspectral/8)
   CALL TIMESTEPNKPC_B(.true., wadj, wadjb, padj, padjb, siadj, sjadj, &
   &                  skadj, sfaceiadj, sfacejadj, sfacekadj, voladj, &
   &                  radiadj, radiadjb, radjadj, radjadjb, radkadj, &
   &                  radkadjb, icell, jcell, kcell, pinfcorr, rhoinf, nn, &
   &                  level, sps, sps2)
   CALL POPREAL8ARRAY(padj, realtype*5**3*ntimeintervalsspectral/8)
   CALL COMPUTEPRESSURENKPC_B(wadj, wadjb, padj, padjb, nn, level, sps&
   &                         , sps2)
   END DO
   END SUBROUTINE COMPUTERNKPC_B

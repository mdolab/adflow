   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
   !  
   !  Differentiation of applyallbcnkpc in reverse (adjoint) mode:
   !   gradient, with respect to input variables: padj wadj
   !   of linear combination of output variables: padj wadj
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          applyAllBCAdj.f90                               *
   !      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
   !      *                C.A.(Sandy) Mader                               *
   !      * Starting date: 04-16-2008                                      *
   !      * Last modified: 04-17-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE APPLYALLBCNKPC_B(winfadj, pinfcorradj, wadj, wadjb, padj, &
   &  padjb, sadj, siadj, sjadj, skadj, voladj, normadj, rfaceadj, icell, &
   &  jcell, kcell, secondhalo, nn, level, sps, sps2)
   USE bctypes
   USE blockpointers
   USE flowvarrefstate
   USE inputdiscretization
   USE inputtimespectral
   IMPLICIT NONE
   INTEGER(KIND=INTTYPE) :: icell, jcell, kcell
   INTEGER(KIND=INTTYPE) :: level, nn, sps, sps2
   REAL(KIND=REALTYPE), DIMENSION(nbocos, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: normadj
   REAL(KIND=REALTYPE) :: padj(-2:2, -2:2, -2:2, ntimeintervalsspectral)&
   &  , padjb(-2:2, -2:2, -2:2, ntimeintervalsspectral)
   REAL(KIND=REALTYPE), INTENT(IN) :: pinfcorradj
   REAL(KIND=REALTYPE), DIMENSION(nbocos, -2:2, -2:2, &
   &  ntimeintervalsspectral), INTENT(IN) :: rfaceadj
   REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: sadj
   LOGICAL :: secondhalo
   REAL(KIND=REALTYPE), DIMENSION(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: siadj
   REAL(KIND=REALTYPE), DIMENSION(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: sjadj
   REAL(KIND=REALTYPE), DIMENSION(-3:2, -3:2, -3:2, 3, &
   &  ntimeintervalsspectral), INTENT(IN) :: skadj
   REAL(KIND=REALTYPE), DIMENSION(ntimeintervalsspectral), INTENT(IN) :: &
   &  voladj
   REAL(KIND=REALTYPE) :: wadj(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral), wadjb(-2:2, -2:2, -2:2, nw, &
   &  ntimeintervalsspectral)
   REAL(KIND=REALTYPE), DIMENSION(nw), INTENT(IN) :: winfadj
   INTEGER :: branch
   LOGICAL :: correctfork
   INTEGER(KIND=INTTYPE) :: i, ii, j, jj, k, kk, l
   INTEGER(KIND=INTTYPE) :: iend, istart, jend, jstart, kend, kstart
   CALL PUSHBOOLEAN(secondhalo)
   CALL PUSHREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * applyAllBCAdj applies the possible boundary conditions for the *
   !      * halo cells adjacent to the cell for which the residual needs   *
   !      * to be computed.                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   !, only : ie, ib, je, jb, ke, kb, nBocos, &
   !         BCFaceID, BCType, BCData,p,w
   !precond,choimerkle, etc...
   !nIntervalTimespectral
   !
   !      Subroutine arguments.
   !
   !
   !      Local variables.
   !
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   CALL BCSYMMNKPC(wadj, padj, normadj, icell, jcell, kcell, secondhalo, &
   &            nn, level, sps, sps2)
   SELECT CASE  (precond) 
   CASE (noprecond) 
   CALL PUSHREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   CALL PUSHREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   CALL PUSHBOOLEAN(secondhalo)
   CALL BCFARFIELDNKPC(secondhalo, winfadj, pinfcorradj, wadj, padj, &
   &                  siadj, sjadj, skadj, normadj, rfaceadj, icell, jcell&
   &                  , kcell, nn, level, sps, sps2)
   CALL PUSHINTEGER4(3)
   CASE (turkel) 
   CALL PUSHINTEGER4(1)
   CASE (choimerkle) 
   CALL PUSHINTEGER4(0)
   CASE DEFAULT
   CALL PUSHINTEGER4(2)
   END SELECT
   CALL BCEULERWALLNKPC_B(secondhalo, wadj, wadjb, padj, padjb, sadj, &
   &                   siadj, sjadj, skadj, normadj, rfaceadj, icell, jcell&
   &                   , kcell, nn, level, sps, sps2)
   CALL POPINTEGER4(branch)
   IF (.NOT.branch .LT. 2) THEN
   IF (.NOT.branch .LT. 3) THEN
   CALL POPBOOLEAN(secondhalo)
   CALL POPREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   CALL POPREAL8ARRAY(padj, 5**3*ntimeintervalsspectral)
   CALL BCFARFIELDNKPC_B(secondhalo, winfadj, pinfcorradj, wadj, &
   &                      wadjb, padj, padjb, siadj, sjadj, skadj, normadj&
   &                      , rfaceadj, icell, jcell, kcell, nn, level, sps, &
   &                      sps2)
   END IF
   END IF
   CALL POPREAL8ARRAY(wadj, 5**3*nw*ntimeintervalsspectral)
   CALL POPBOOLEAN(secondhalo)
   CALL BCSYMMNKPC_B(wadj, wadjb, padj, padjb, normadj, icell, jcell, &
   &              kcell, secondhalo, nn, level, sps, sps2)
   END SUBROUTINE APPLYALLBCNKPC_B

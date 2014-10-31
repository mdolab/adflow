   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of sa_block in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *rev *dw *w *rlv (global)timeref
   !   with respect to varying inputs: *rev *dw *w *rlv (global)timeref
   !   Plus diff mem management of: rev:in dw:in w:in rlv:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          sa.f90                                          *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 06-11-2003                                      *
   !      * Last modified: 04-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE SA_BLOCK_B(resonly)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * sa solves the transport equation for the Spalart-Allmaras      *
   !      * turbulence model in a segregated manner using a diagonal       *
   !      * dominant ADI-scheme.                                           *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_B
   USE INPUTTIMESPECTRAL
   USE ITERATION
   IMPLICIT NONE
   !
   !      Subroutine argument.
   !
   LOGICAL, INTENT(IN) :: resonly
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, sps
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   ! Set the arrays for the boundary condition treatment.
   ! Solve the transport equation for nuTilde.
   CALL PUSHREAL8ARRAY(w, SIZE(w, 1)*SIZE(w, 2)*SIZE(w, 3)*SIZE(w, 4))
   CALL PUSHREAL8ARRAY(dw, SIZE(dw, 1)*SIZE(dw, 2)*SIZE(dw, 3)*SIZE(dw, 4&
   &               ))
   CALL SASOLVE(resonly)
   ! The eddy viscosity and the boundary conditions are only
   ! applied if an actual update has been computed in saSolve.
   IF (.NOT.resonly) CALL SAEDDYVISCOSITY_B()
   CALL POPREAL8ARRAY(dw, SIZE(dw, 1)*SIZE(dw, 2)*SIZE(dw, 3)*SIZE(dw, 4)&
   &             )
   CALL POPREAL8ARRAY(w, SIZE(w, 1)*SIZE(w, 2)*SIZE(w, 3)*SIZE(w, 4))
   CALL SASOLVE_B(resonly)
   END SUBROUTINE SA_BLOCK_B

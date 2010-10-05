!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBCForcesAdj.f90                         *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      *                C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine applyAllBCForceCouplingAdj(wInfAdj,pInfCorrAdj,wAdj, pAdj, &
                              siAdj, sjAdj, skAdj, normAdj, &
                              iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End,&
                              secondHalo,mm)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllBCAdj applies the possible boundary conditions for the *
!      * halo cells adjacent to the cell for which the residual needs   *
!      * to be computed.                                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers!, only : ie, ib, je, jb, ke, kb, nBocos, &
                         !         BCFaceID, BCType, BCData,p,w
       use flowVarRefState
       use inputDiscretization !precond,choimerkle, etc...
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo
       integer(kind=intType), intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
       integer(kind=intType), intent(in) :: i2Beg,i2End,j2Beg,j2End
       integer(kind=intType), intent(in) :: mm

       real(kind=realType), dimension(0:ib,0:jb,0:kb,nw), intent(in) :: wAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb), intent(in) :: pAdj
       real(kind=realType), intent(in) ::pInfCorrAdj
       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd,3) :: siAdj 
       ! notice the range of y dim is set 1:2 which corresponds to 1/jl
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd,3) :: sjAdj
       ! notice the range of z dim is set 1:2 which corresponds to 1/kl
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2,3) :: skAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: normAdj
       real(kind=realType), dimension(nw),intent(in)::wInfAdj

!
!      Local variables.
!
       integer(kind=intType) :: nn, sps
       integer(kind=intType)::i,j,k,ii,jj,kk,l
       integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       
       call bcSymmForceCouplingAdj(secondHalo, wAdj,pAdj,normAdj,mm,&
            iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End)

       select case (precond)
          
       case (noPrecond)
          call bcFarfieldForceCouplingAdj(secondHalo,wInfAdj,pInfCorrAdj, wAdj,pAdj,&
               siAdj, sjAdj, skAdj, normAdj,mm,&
               iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End)

       case (Turkel)
          call terminate("applyAllBC", "Farfield boundary conditions for Turkel preconditioner not implemented")
          
       case (ChoiMerkle)
          call terminate("applyAllBC", "Farfield boundary conditions for Choi and Merkle preconditioner not implemented")

       end select

       ! Inviscid wall boundary conditions.
       call bcEulerWallForceCouplingAdj(secondHalo, wAdj,pAdj,      &
            siAdj, sjAdj, skAdj, normAdj,mm,&
            iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End)

     end subroutine applyAllBCForceCouplingAdj

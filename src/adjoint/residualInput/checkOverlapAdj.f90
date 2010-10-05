!
!      ******************************************************************
!      *                                                                *
!      * File:          checkOverlapAdj.f90                             *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-23-2008                                      *
!      * Last modified: 04-23-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkOverlapAdj(nn,icell,jcell,kcell,isbeg,jsbeg,&
            ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
            computeBC)

         
       !modules
       use BCTypes
       use blockPointers, only : ie, ib, je, jb, ke, kb, nBocos, &
                                  BCFaceID, BCType, BCData
       use flowVarRefState
       implicit none


       !subroutine arguments

       integer(kind=intType), intent(in) :: nn!, offset
       integer(kind=intType), intent(in) ::iCell, jCell,kCell
       integer(kind=intType) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
       integer(kind=intType) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
       logical, intent(out) :: computeBC
         

       !local Variables
       integer(kind=intType) :: iOffset, jOffset, kOffset
       logical :: iOverlap, jOverlap, kOverlap

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the range of the stencil for the given cell.

       iSBeg = iCell - 2; iSEnd = iCell + 2
       jSBeg = jCell - 2; jSEnd = jCell + 2
       kSBeg = kCell - 2; kSEnd = kCell + 2

       iOffset = 0!zero
       jOffset = 0!zero
       kOffset = 0!zero

       ! Determine the range of halo cells which this boundary subface
       ! will change.
       
       select case (BCFaceID(nn))
       case (iMin)
          iBBeg = 0;                iBEnd = 1
          jBBeg = BCData(nn)%icBeg; jBEnd = BCData(nn)%icEnd
          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
          
          !=============================================================
          
       case (iMax)
          iBBeg = ie;               iBEnd = ib
          jBBeg = BCData(nn)%icBeg; jBEnd = BCData(nn)%icEnd
          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
          
          !=============================================================
          
       case (jMin)
          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
          jBBeg = 0;                jBEnd = 1
          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
          
          !============================================================= 
      case (jMax)
          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
          jBBeg = je;               jBEnd = jb
          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
          
          !=============================================================
          
       case (kMin)
          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
          jBBeg = BCData(nn)%jcBeg; jBEnd = BCData(nn)%jcEnd
          kBBeg = 0;                kBEnd = 1
          
          !=============================================================
          
       case (kMax)
          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
          jBBeg = BCData(nn)%jcBeg; jBEnd = BCData(nn)%jcEnd
          kBBeg = ke;               kBEnd = kb
          
       end select


       iOverlap = .false.
       if(iSBeg <= iBEnd .and. iSEnd >= iBBeg) iOverlap = .true.

       jOverlap = .false.
       if(jSBeg <= jBEnd .and. jSEnd >= jBBeg) jOverlap = .true.
       
       kOverlap = .false.
       if(kSBeg <= kBEnd .and. kSEnd >= kBBeg) kOverlap = .true.

       
       computeBC = .false.
       checkOverlap: if(iOverlap .and. jOverlap .and. kOverlap) then
          computeBC = .true.
     
       endif checkOverlap
       
       
     end subroutine checkOverlapAdj

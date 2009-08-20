!
!      ******************************************************************
!      *                                                                *
!      * File:          checkXOverlapAdj.f90                            *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 08-10-2009                                      *
!      * Last modified: 08-10-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkXOverlapAdj(icell,jcell,kcell,ioffset,joffset,&
            koffset,iMinOverlap, jMinOverlap, kMinOverlap ,&
            iMaxOverlap, jMaxOverlap, kMaxOverlap)

         
       !modules
       use BCTypes
       use blockPointers, only : ie, ib, je, jb, ke, kb, nBocos, &
                                  BCFaceID, BCType, BCData
       use flowVarRefState
       implicit none


       !subroutine arguments
       integer(kind=intType), intent(in) ::iCell, jCell,kCell
       integer(kind=intType) :: iOffset, jOffset, kOffset
       logical :: iMinOverlap, jMinOverlap, kMinOverlap 
       logical :: iMaxOverlap, jMaxOverlap, kMaxOverlap 

       !local Variables
       integer(kind=intType) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
       integer(kind=intType) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the range of the stencil for the given cell.

       iSBeg = iCell - 3; iSEnd = iCell + 2
       jSBeg = jCell - 3; jSEnd = jCell + 2
       kSBeg = kCell - 3; kSEnd = kCell + 2

       iOffset = zero
       jOffset = zero
       kOffset = zero

       iBBeg = 0; iBEnd = ie
       jBBeg = 0; jBEnd = je
       kBBeg = 0; kBEnd = ke


!!$       ! Determine the range of halo cells which this boundary subface
!!$       ! will change.
!!$       
!!$       select case (BCFaceID(nn))
!!$       case (iMin)
!!$          iBBeg = 0;                iBEnd = 0
!!$          jBBeg = BCData(nn)%inBeg; jBEnd = BCData(nn)%inEnd
!!$          kBBeg = BCData(nn)%jnBeg; kBEnd = BCData(nn)%jnEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (iMax)
!!$          iBBeg = ie;               iBEnd = ie
!!$          jBBeg = BCData(nn)%inBeg; jBEnd = BCData(nn)%inEnd
!!$          kBBeg = BCData(nn)%jnBeg; kBEnd = BCData(nn)%jnEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (jMin)
!!$          iBBeg = BCData(nn)%inBeg; iBEnd = BCData(nn)%inEnd
!!$          jBBeg = 0;                jBEnd = 0
!!$          kBBeg = BCData(nn)%jnBeg; kBEnd = BCData(nn)%jnEnd
!!$          
!!$          !============================================================= 
!!$      case (jMax)
!!$          iBBeg = BCData(nn)%inBeg; iBEnd = BCData(nn)%inEnd
!!$          jBBeg = je;               jBEnd = je
!!$          kBBeg = BCData(nn)%jnBeg; kBEnd = BCData(nn)%jnEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (kMin)
!!$          iBBeg = BCData(nn)%inBeg; iBEnd = BCData(nn)%inEnd
!!$          jBBeg = BCData(nn)%jnBeg; jBEnd = BCData(nn)%jnEnd
!!$          kBBeg = 0;                kBEnd = 0
!!$          
!!$          !=============================================================
!!$          
!!$       case (kMax)
!!$          iBBeg = BCData(nn)%inBeg; iBEnd = BCData(nn)%inEnd
!!$          jBBeg = BCData(nn)%jnBeg; jBEnd = BCData(nn)%jnEnd
!!$          kBBeg = ke;               kBEnd = ke
!!$          
!!$       end select


       iMinOverlap = .false.
       if(iSBeg <= iBBeg .and. iSEnd >= iBBeg) then
          iMinOverlap = .true.
          ioffset= iBBeg-icell
       end if

       iMaxOverlap = .false.
       if(iSBeg <= iBEnd .and. iSEnd >= iBEnd)then
          iMaxOverlap = .true.
          ioffset= iBEnd-icell
       endif

       jMinOverlap = .false.
       if(jSBeg <= jBBeg .and. jSEnd >= jBBeg)then
          jMinOverlap = .true.
          joffset = jBBeg-jcell
       endif
       
       jMaxOverlap = .false.
       if(jSBeg <= jBEnd .and. jSEnd >= jBEnd) then
          jMaxOverlap = .true.
          joffset = jBEnd-jcell
       endif
       
       kMinOverlap = .false.
       if(kSBeg <= kBBeg .and. kSEnd >= kBBeg) then
          kMinOverlap = .true.
          koffset = kBBeg-kcell
       endif

       kMaxOverlap = .false.
       if(kSBeg <= kBEnd .and. kSEnd >= kBEnd) then
          kMaxOverlap = .true.
          koffset = kBEnd-kcell
       endif
       

       
       
  end subroutine checkXOverlapAdj

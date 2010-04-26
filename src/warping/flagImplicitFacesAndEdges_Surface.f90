!
! ***********************************
! *  File: flagImplicitEdgesAndFaces_Surface.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 26-04-2010
! *  Modified: 26-04-2008
! ***********************************

subroutine flagImplicitEdgesAndFacesSurface

!*****************************************************
!
!  This routine identifies what edges and faces are explicitly and 
!  implicitly perturbed for a give block.
!
!*****************************************************

use blockPointers
use mdDataLocal
implicit none

!Subroutine Arguments


!Local Variables
integer(kind=intType)::mm,ll,ii,jj,kk,i,j,k,m,n
integer(kind=intType)::edge,face,counter
integer(kind=intType),dimension(3)::ijk_num
logical::cornerPoint,on_edge,on_face,is_corner,edgepoint,facepoint

integer(kind=intType),dimension(3,8)::relatedFaces,relatedEdges
integer(kind=intType),dimension(2,12)::edgeRelatedFaces
integer(kind=intType),dimension(12,6)::searchPattern

integer(kind=intType)::level=1,sps=1
integer(kind=intType)::which_corner,on_which_edge
    
!
! Begin Execution
!  
do nn = 1,nDom
   if(.not. associated(flowdoms(i,level,sps)%ifaceptb))then
      allocate(flowdoms(i,level,sps)%ifaceptb(6),status=ierr)
   endif
   if(.not. associated(flowdoms(i,level,sps)%iedgeptb))then
      allocate(flowdoms(i,level,sps)%iedgeptb(12),status=ierr)
   endif
end do

do nn = 1,nDom
   call setpointers(nn,level,sps)
   IFACEPTB(:)=0
   IEDGEPTB(:)=0
enddo
call blockRelations(relatedFaces,relatedEdges,edgeRelatedFaces,&
     searchPattern)

do i = 1,size(mdSurfGlobalIndLocal(5,:))
   !update if they match the current global surface node
   call setPointers(mdSurfGlobalIndLocal(4,i),level,sps)
   ii = mdSurfGlobalIndLocal(1,i)
   jj = mdSurfGlobalIndLocal(2,i)
   kk = mdSurfGlobalIndLocal(3,i)
   IJK_NUM = (/ii,jj,kk/)
   !check the location of this point
   cornerPoint = IS_CORNER(IJK_NUM)
   edgePoint = ON_EDGE(IJK_NUM)
   facePoint = ON_FACE(IJK_NUM)
   
   if (cornerPoint) then
      !print *,'corner perturbed',i
      mm = which_corner(IJK_NUM)
      do n =1,3
         edge = relatedEdges(n,mm)
         !print *,'edge',edge
         if (.not. iedgeptb(edge)==2)then
            IEDGEPTB(edge)= 1
         endif
         face = relatedFaces(n,mm)
         !print *,'face',face
         if (.not. IFACEPTB(face)==2) then
            IFACEPTB(face)=1
         endif
      end do
   endif
   if (edgepoint )then
      !print *,'edgeperturbed',i,j,k,local,local0!,n,mm
      mm = on_which_edge(IJK_NUM)
      IEDGEPTB(i) = 2
      do m=1,2
         face = edgeRelatedFaces(m,mm)
         if (.not. IFACEPTB(face)==2)then
            IFACEPTB(face)=1
         endif
      end do
      exit!break
   else
      IEDGEPTB(mm) = IEDGEPTB(mm) 
   endif
      
   if ((.not. cornerPoint) .and. (.not. edgePoint) .and.&
        facePoint)then
      if (IJK_NUM(1) == 1) then
         IFACEPTB(1) = 2
         exit!break
      elseif (IJK_NUM(1) == il) then
         IFACEPTB(2) = 2
         exit!break
      elseif (IJK_NUM(2) == 1) then
         IFACEPTB(3) = 2
         exit!break
      elseif (IJK_NUM(2) == jl) then
         IFACEPTB(4) = 2
         exit!break
      elseif (IJK_NUM(3) == 1) then
         IFACEPTB(5) = 2
         exit!break
      elseif (IJK_NUM(3) == kl) then
         IFACEPTB(6) = 2
         exit!break
      END IF
   elseif (.not. facePoint)then
      print *,'WARNING - ILLEGAL PERTURBATION!  PERTURBATION PASSED TO '
      print *,'WARP IS IN THE INTERIOR OF THE BLOCK, RATHER THAN ON A'
      print *,'BLOCK FACE AS IT SHOULD BE - WARP WILL LIKELY FAIL'
      print *,'ILLEGAL BLOCK PERTURBATION OF A INTERIOR POINT IN A BLOCK'
      stop
   endif
end do

end subroutine flagImplicitEdgesAndFacesSurface



 

!
! ***********************************
! *  File: flagImplicitEdgesAndFaces.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-11-2008
! *  Modified: 12-11-2008
! ***********************************

subroutine flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)

!*****************************************************
!
!  This routine identifies what edges and faces are explicitly and 
!  implicitly perturbed for a give block. Converted from python...
!
!*****************************************************

use blockPointers
implicit none

!Subroutine Arguments

integer(kind=intType),dimension(6)::IFACEPTB
integer(kind=intType),dimension(12)::IEDGEPTB

!Local Variables
integer(kind=intType)::mm,ll,ii,jj,kk,i,j,k,m,n
integer(kind=intType)::edge,face,counter
integer(kind=intType),dimension(3)::ijk_num
logical::cornerPoint,on_edge,on_face,is_corner,edgepoint,facepoint
real(kind=realType):: local,local0,tolerance

integer(kind=intType),dimension(3,8)::relatedFaces,relatedEdges
integer(kind=intType),dimension(2,12)::edgeRelatedFaces
integer(kind=intType),dimension(12,6)::searchPattern

integer(kind=intType),dimension(nSubface)::incrementI,&
     incrementJ,incrementK

!Create storage for corner perturbation info
logical,dimension(8)::perturbedCorner = (/.False.,.False.,.False.,&
     .False.,.False.,.False.,.False.,.False./)
    
!
! Begin Execution
!

! Determine whether the coordinates are increasing or
! decreasing in each direction for each subface

do i =1,nSubface
   !check for +ve vs -ve increment
   if (inend(i) >=inbeg(i)) then
      incrementI(i) = 1
   else
       incrementI(i) = -1
   endif
                
   if ( jnend(i) >= jnbeg(i)) then
      incrementJ(i) = 1
   else
      incrementJ(i) = -1
   endif
      
   if ( knend(i) >= knbeg(i)) then
      incrementK(i) = 1
   else
      incrementK(i) = -1
   endif
end do
  

!Now for a single block only!!

!loop through blocks and faces checking faces
!do nn=1,nDom
!   call setPointer(nn,level,sps)
      !set the lists for this block to zero
IFACEPTB(:)=0
IEDGEPTB(:)=0
do mm= 1,nSubface
!!$   do ii=inbeg(mm),inend(mm)+incrementI(mm),incrementI(mm)
!!$      do jj=jnbeg(mm),jnend(mm)+incrementJ(mm),incrementJ(mm)
!!$         do kk=knbeg(mm),knend(mm)+incrementK(mm),incrementK(mm)
   do ii=inbeg(mm),inend(mm),incrementI(mm)
      do jj=jnbeg(mm),jnend(mm),incrementJ(mm)
         do kk=knbeg(mm),knend(mm),incrementK(mm)
            !check whether this point has moved
            !print *,'indicies',ii,jj,kk,mm
            do ll=1,3
               local = x(ii,jj,kk,ll)
               local0 = xInit(ii,jj,kk,ll)
               tolerance = 1.0e-12
               if (abs(local -local0)/max(abs(local0),abs(local),tolerance)>1e-12.and. abs(local -local0)>tolerance )then
                  !print *,'ifcheck movement',abs(local -local0)/max(abs(local0),abs(local),tolerance)>1e-12,local,local0,abs(local -local0),max(abs(local0),abs(local),tolerance),abs(local -local0)/max(abs(local0),abs(local),tolerance)
                  !point has moved
                  !set the index for this point
                  IJK_NUM = (/ii,jj,kk/)
                     !check the location of this point
                  cornerPoint = IS_CORNER(IJK_NUM)
                  edgePoint = ON_EDGE(IJK_NUM)
                  facePoint = ON_FACE(IJK_NUM)
                  !print *,'logicals',cornerPoint,edgePoint,facePoint,(.not. cornerPoint) .and. (.not. edgePoint) .and. facePoint,ii,jj,kk,mm,local -local0,local ,local0
                  ! Flag the explicitly perturbed faces
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
               endif
            end do
         end do
      end do
   end do
end do
!end do

!print *,'ifaceptb',ifaceptb,iedgeptb
!!identify perturbed corners on each block
!do nn=1,nDom
!   call setPointer(nn,level,sps)
counter = 1
do k=1,kl,kl-1!ijk(3),ijk(3)!jump corner to corner
   do j=1,jl,jl-1!ijk(2),ijk(2)
      do i=1,il,il-1!ijk(1),ijk(1)
         do n =1,3
            !print *,'i,j,k',i,j,k,n
            local = x(i,j,k,n)
            local0 = xInit(i,j,k,n)
            tolerance = 1.0e-12
            !#             if (abs(local.real -local0.real)/max(abs(local0.real),abs(local.real),tolerance)>1e-6 or abs(local.imag -local0.imag)/max(abs(local0.imag),abs(local.imag),tolerance)>1e-6):
            if (abs(local -local0)/max(abs(local0),abs(local),tolerance)>1e-12.and. abs(local -local0)>tolerance)then
               !print *,'corner perturbed',local,local0,abs(local -local0)/max(abs(local0),abs(local),tolerance),i,j,k,n
               perturbedCorner(counter) = .True.
               exit!break
            else
               perturbedCorner(counter) = .False.
            endif
         end do
         counter= counter+1
      end do
   end do
end do
!end do
!print *,'perturbedCorner',perturbedCorner
!get the Block relations for the next part of the algorithm
call blockRelations(relatedFaces,relatedEdges,edgeRelatedFaces,&
     searchPattern)

!identify the edges of the block

! Flag the edges and faces  connected with each
! perturbed corner as implicitly perturbed
!do nn=1,nDom
!   call setPointer(nn,level,sps)
do i = 1,8!(len(perturbedCorner))
   if (perturbedCorner(i)) then
      !print *,'corner perturbed',i
      do n =1,3
         edge = relatedEdges(n,i)
         !print *,'edge',edge
         IEDGEPTB(edge)= 1
         face = relatedFaces(n,i)
         !print *,'face',face
         if (.not. IFACEPTB(face)==2) then
            IFACEPTB(face)=1
         endif
      end do
   endif
end do
!end do
!print *,'ifaceptb',ifaceptb,iedgeptb        
! Flag the explicitly perturbed edges
!do nn=1,nDom
!   call setPointer(nn,level,sps)
do mm = 1,12!len(searchPattern)!Loop over edges
   !print *,'mm',mm,'i',searchPattern(mm,1),searchPattern(mm,2),'j',searchPattern(mm,3),searchPattern(mm,4),'k',searchPattern(mm,5),searchPattern(mm,6),shape(x),ie,je,ke
   do i=searchPattern(mm,1),searchPattern(mm,2)
      do j=searchPattern(mm,3),searchPattern(mm,4)
         do k=searchPattern(mm,5),searchPattern(mm,6)
            !print *,'i,j,k',i,j,k
            do n=1,3
               local = x(i,j,k,n)
               local0 = xInit(i,j,k,n)
               tolerance = 1.0e-12
!#if (abs(local.real -local0.real)/max(abs(local0.real),abs(local.real),tolerance)>1e-6 or abs(local.imag -local0.imag)/max(abs(local0.imag),abs(local.imag),tolerance)>1e-6):
               if (abs(local -local0)/max(abs(local0),abs(local),tolerance)>1e-12.and. abs(local -local0)>tolerance )then
                  !print *,'edgeperturbed',i,j,k,local,local0!,n,mm
                  IEDGEPTB(mm) = 2
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
            end do
         end do
      end do
   end do
end do
!end do

end subroutine flagImplicitEdgesAndFaces



 

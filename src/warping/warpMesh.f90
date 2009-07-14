!
! ***********************************
! *  File: warpMesh.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 07-13-2009
! *  Modified: 07-13-2009
! ***********************************

subroutine warpMesh

  ! Warp the current mesh to take into acccout updates from
  ! updateGlobalFaces

  use blockpointers
  implicit none
  !Subroutine Arguments
 
  
  ! Local Arguments
  
  integer(kind=intType)::level=1
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB  

  !begin execution

  ! Now warp the blocks based on new face coordinates
  
  do nn=1,nDom
     !print*,'warpingblock',nn
     call setPointers(nn,level,sps)     
     
     !print *,'flag implicites'
     !determine the explicitly and implicitly perturbed faces and edges
     call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)

     ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
     
     ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
     IMAX = IL
     JMAX = JL
     KMAX = KL
     !print *,'allocate xyz0'
     ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
     xyz0 = 0
     xyznew = 0
     
     XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,1)
     XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,2)
     XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,3)
     XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
     XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
     XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)
     !print *,'call warp local'
     call warp_local(xyznew,xyz0,ifaceptb,iedgeptb,imax,jmax,kmax)
     
     !print *,'update coordinates'
     ! ASSIGN THESE NEW XYZ VALUES TO THE MESH ITSELF
     DO I=1,IMAX
        DO J=1,JMAX
           DO K=1,KMAX
              X(I,J,K,1) = XYZNEW(1,I,J,K)
              X(I,J,K,2) = XYZNEW(2,I,J,K)
              X(I,J,K,3) = XYZNEW(3,I,J,K)
           END DO
        END DO
     END DO
      deALLOCATE(XYZ0,XYZNEW)
  end do
 
end subroutine warpMesh

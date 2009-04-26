!
! ***********************************
! *  File: warp_local.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-07-2008
! *  Modified: 12-07-2008
! ***********************************

subroutine warp_local(xyznew,xyz0,ifaceptb,iedgeptb,imax,jmax,kmax)

  use blockpointers
  implicit none

  ! Subroutine Arguments
  real(kind=realType),dimension(3,0:IMAX+1,0:JMAX+1,0:KMAX+1)::xyznew,xyz0
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB
  integer(kind=intType)::imax,jmax,kmax
  ! Local Variables

  real(kind=realType),allocatable,dimension(:,:,:,:,:)::dFaceI,dFaceJ,dFaceK
  real(kind=realType),allocatable,dimension(:,:,:,:)::s0


  !**************************
  ! Begin execution
  !**************************

  ! SAVE DFACE VALUES TO BE PASSED TO WARPBLK
  !          ALLOCATE(DFACEI(3,0:JMAX+1,0:KMAX+1,2,4),DFACEJ(3,0:IMAX+1,0:KMAX+1,2,4),&
  !               DFACEK(3,0:IMAX+1,0:JMAX+1,2,4))
  print *,'allocating dface'
  ALLOCATE(DFACEI(3,1:JMAX,1:KMAX,2,4),DFACEJ(3,1:IMAX,1:KMAX,2,4),&
       DFACEK(3,1:IMAX,1:JMAX,2,4))
  
  DFACEI = 0.0
  DFACEJ = 0.0
  DFACEK = 0.0
  DFACEI(1,1:JMAX,1:KMAX,1,1) = Xyznew(1,1,1:JMAX,1:KMAX)- Xyz0(1,1,1:JMAX,1:KMAX)
  DFACEI(2,1:JMAX,1:KMAX,1,1) = Xyznew(2,1,1:JMAX,1:KMAX)- Xyz0(2,1,1:JMAX,1:KMAX)
  DFACEI(3,1:JMAX,1:KMAX,1,1) = Xyznew(3,1,1:JMAX,1:KMAX)- Xyz0(3,1,1:JMAX,1:KMAX)
  DFACEI(1,1:JMAX,1:KMAX,2,1) = Xyznew(1,IMAX,1:JMAX,1:KMAX)-Xyz0(1,IMAX,1:JMAX,1:KMAX)
  DFACEI(2,1:JMAX,1:KMAX,2,1) =Xyznew(2,IMAX,1:JMAX,1:KMAX)-Xyz0(2,IMAX,1:JMAX,1:KMAX)
  DFACEI(3,1:JMAX,1:KMAX,2,1) =Xyznew(3,IMAX,1:JMAX,1:KMAX)-Xyz0(3,IMAX,1:JMAX,1:KMAX)
  DFACEJ(1,1:IMAX,1:KMAX,1,1) =Xyznew(1,1:IMAX,1,1:KMAX)-Xyz0(1,1:IMAX,1,1:KMAX)
  DFACEJ(2,1:IMAX,1:KMAX,1,1) =Xyznew(2,1:IMAX,1,1:KMAX)-Xyz0(2,1:IMAX,1,1:KMAX)
  DFACEJ(3,1:IMAX,1:KMAX,1,1) =Xyznew(3,1:IMAX,1,1:KMAX)-Xyz0(3,1:IMAX,1,1:KMAX)
  DFACEJ(1,1:IMAX,1:KMAX,2,1) =Xyznew(1,1:IMAX,JMAX,1:KMAX)-Xyz0(1,1:IMAX,JMAX,1:KMAX)
  DFACEJ(2,1:IMAX,1:KMAX,2,1) =Xyznew(2,1:IMAX,JMAX,1:KMAX)-Xyz0(2,1:IMAX,JMAX,1:KMAX)
  DFACEJ(3,1:IMAX,1:KMAX,2,1) =Xyznew(3,1:IMAX,JMAX,1:KMAX)-Xyz0(3,1:IMAX,JMAX,1:KMAX)
  DFACEK(1,1:IMAX,1:JMAX,1,1) =Xyznew(1,1:IMAX,1:JMAX,1)-Xyz0(1,1:IMAX,1:JMAX,1)
  DFACEK(2,1:IMAX,1:JMAX,1,1) =Xyznew(2,1:IMAX,1:JMAX,1)-Xyz0(2,1:IMAX,1:JMAX,1)
  DFACEK(3,1:IMAX,1:JMAX,1,1) =Xyznew(3,1:IMAX,1:JMAX,1)-Xyz0(3,1:IMAX,1:JMAX,1)
  DFACEK(1,1:IMAX,1:JMAX,2,1) =Xyznew(1,1:IMAX,1:JMAX,KMAX)-Xyz0(1,1:IMAX,1:JMAX,KMAX)
  DFACEK(2,1:IMAX,1:JMAX,2,1) =Xyznew(2,1:IMAX,1:JMAX,KMAX)-Xyz0(2,1:IMAX,1:JMAX,KMAX)
  DFACEK(3,1:IMAX,1:JMAX,2,1) =Xyznew(3,1:IMAX,1:JMAX,KMAX)-Xyz0(3,1:IMAX,1:JMAX,KMAX)
  print *,'allocating s0'
  ALLOCATE(S0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
  s0 = 0.0
  ! CALL WARPBLK FOR THE CURRENT BLOCK IF ANY OF THE FACES OR EDGES IN THAT 
  ! BLOCK ARE IMPLICITLY OR EXPLICITLY PERTURBED
  print *,'if maxval',ifaceptb,'edge',iedgeptb
  IF (MAXVAL(IFACEPTB(1:6)) >= 1 .OR. MAXVAL(IEDGEPTB(1:12)) >= 1) THEN
     print *,'calling warpblk'
     CALL WARPBLK(IFACEPTB(1:6),IEDGEPTB(1:12),-3,0,1,IMAX,JMAX,&
          KMAX,XYZ0,S0,DFACEI,DFACEJ,DFACEK,XYZNEW)
     
  END IF
  print *,'deallocating'
  DEALLOCATE(S0,DFACEI,DFACEJ,DFACEK)
END SUBROUTINE WARP_LOCAL


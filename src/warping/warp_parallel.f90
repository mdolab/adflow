!
!      ******************************************************************
!      *                                                                *
!      * FILE:          warp_parallel.f90                               *
!      * AUTHOR:        C.A.(Sandy) Mader                               *
!      * STARTING DATE: 11-01-2007                                      *
!      * LAST MODIFIED: 11-01-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine warp_parallel(ifaceptb, iedgeptb, ncall, ig, igo, il, jl, kl,xyz0, s0, dfacei, dfacej, dfacek, xyz)
!
!      ******************************************************************
!      *                                                                *
!      * This is and intrface rountine between the main python driver   *
!      * Code and the fortran single block warping code WARPBLK         *
!      *                                                                *
!      ******************************************************************
  use complexify
!
!      SUBROUTINE ARGUMENTS.
!
  integer*4,intent(in) ::ifaceptb(6)
  integer*4,intent(in) ::iedgeptb(12)
  integer*4,intent(in) ::ncall
  integer*4,intent(in) ::ig
  integer*4,intent(in) ::igo                     
  integer*4,intent(in) ::il, jl, kl   
!         real(kind=real_type), dimension(3,il,jl,kl), intent(in) :: x0   
!  real, dimension(3,0:il+1,0:jl+1,0:kl+1) :: x0            
  complex*16, dimension(3,0:il+1,0:jl+1,0:kl+1) :: xyz0
  complex*16, dimension(3,0:il+1,0:jl+1,0:kl+1) :: s0
!         
  complex*16, dimension(3,0:jl+1,0:kl+1,2,4) :: dfacei 
  complex*16, dimension(3,0:il+1,0:kl+1,2,4) :: dfacej         
  complex*16, dimension(3,0:il+1,0:jl+1,2,4) :: dfacek         
  complex*16, dimension(3,0:il+1,0:jl+1,0:kl+1):: xyz 
         
!      LOCAL VARIABLES.

!!$  do i=0,il+1
!!$     do j = 0,jl+1
!!$        do k=0,kl+1
!!$           do n=1,3
!!$              if((real(xyz(n,I,J,K))-real(xyz0(n,I,J,K)))>1e-15)then
!!$                 print *,'warp parallel coords',real(xyz(n,I,J,K)),&
!!$                      real(xyz0(n,I,J,K)),i,j,k,n
!!$              endif
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo
  !print *,'xyz',xyz
  !print *,'in warp_parallel',IFACEPTB,'edge', IEDGEPTB,'ncall', NCALL,&
  !     IG, IGO, IL, JL, KL
  ! BLOCK ARE IMPLICITLY OR EXPLICITLY PERTURBED
  IF (MAXVAL(IFACEPTB(1:6)) >= 1 .OR. MAXVAL(IEDGEPTB(1:12)) >= 1) THEN

     CALL WARPBLK(IFACEPTB(1:6),IEDGEPTB(1:12),ncall,ig,igo,IL,JL,KL,XYZ0,S0,&
                 DFACEI,DFACEJ,DFACEK,XYZ)
  endif

  
  
end subroutine warp_parallel




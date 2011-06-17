
!
!      ******************************************************************
!      *                                                                *
!      * File:          storeReferenceMesh.f90                          *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 12-11-2008                                      *
!      * Last modified: 12-11-2008                                      *
!      *                                                                *
!      ******************************************************************
!

subroutine storeReferenceMesh

! use blockPointers
! use inputTimeSpectral
! implicit none

! !
! ! Local Variables
! !

! integer(kind=intType)::i,j,k,nn,sps

! !
! !Begin Execution
! !
! !perhaps call xhalo(level) here?
! do sps = 1, nTimeIntervalsSpectral
!    do nn=1,nDom
!       !print *,'prepointers',flowdoms(nn,1,sps)%x(1,1,1,1),nn,sps
!       call setPointers(nn,1,sps)
!       do i = 0,ie
!          do j=0,je
!             do k=0,ke
!                xInit(i,j,k,:) = x(i,j,k,:)
!                !print *,'storex',xInit(i,j,k,1), x(i,j,k,1),i,j,k,flowdoms(nn,1,sps)%x(i,j,k,1)
!             enddo
!          enddo
!       enddo
!       !print *,'refxinit',sum(xinit),nn,sps,sum(x)
!    enddo
! end do

end subroutine storeReferenceMesh

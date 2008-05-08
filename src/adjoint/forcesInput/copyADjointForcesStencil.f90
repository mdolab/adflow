!
!     ******************************************************************
!     *                                                                *
!     * File:          copyADjointForcesStencil.f90                    *
!     * Author:        C.A.(Sandy) Mader                               *
!     *                Seongim Choi                                    *
!     * Starting date: 01-15-2007                                      *
!     * Last modified: 05-06-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine copyADjointForcesStencil(wAdj,xAdj,nn,level,sps)

!
!     ******************************************************************
!     *                                                                *
!     * Transfer the variable  x to the auxiliary stencil              *
!     * xAdj is used by the Tapenade differentiated routines.          *
!     *                                                                *
!     * And compute boundary face normals (siAdj, sjAdj, skAdj)        *
!     *                                                                *
!     * It is assumed that the pointers in blockPointers have already  *
!     * been set.                                                      *
!     *                                                                *
!     ******************************************************************
!
      use blockPointers   ! w,il,jl,kl,ie,je,ke
      use communication   ! myID for debug
      use flowvarrefstate ! nw
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: nn,level,sps    
      real(kind=realType), dimension(0:ie,0:je,0:ke,3), intent(out) :: xAdj
      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw), intent(out) :: wAdj

!
!     Local variables.
!
      integer(kind=intType) :: ii, jj, kk, i, j, k, l, m, n

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

       ! Initialize the auxiliar array wAdj.

      do l=1,nw
        do kk=0,kb
          do jj=0,jb
            do ii=0,ib
              wAdj(ii,jj,kk,l) = zero
            enddo
          enddo
        enddo
      enddo

      ! Set the values of the array xAdj from the original matrix x
      do l=1,nw
         do kk=0,kb
            do jj=0,jb
               do ii=0,ib
                  wAdj(ii,jj,kk,l) = w(ii,jj,kk,l)
               enddo
            enddo
         enddo
      enddo

      ! Initialize the auxiliar array xAdj.

      do l=1,3
        do kk=0,ke
          do jj=0,je
            do ii=0,ie
              xAdj(ii,jj,kk,l) = zero
            enddo
          enddo
        enddo
      enddo

      ! Set the values of the array xAdj from the original matrix x
      do l=1,3
         do kk=0,ke
            do jj=0,je
               do ii=0,ie
                  xAdj(ii,jj,kk,l) = x(ii,jj,kk,l)
               enddo
            enddo
         enddo
      enddo
      
      end subroutine copyADjointForcesStencil

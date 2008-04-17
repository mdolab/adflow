!
!     ******************************************************************
!     *                                                                *
!     * File:          copyADjointStencil.f90                          *
!     * Author:        Andre C. Marta                                  *
!     *                Seongim Choi                                    *
!     * Starting date: 08-03-2006                                      *
!     * Last modified: 11-18-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine copyADjointStencil(wAdj, xAdj, iCell, jCell, kCell)
!
!     ******************************************************************
!     *                                                                *
!     * Transfer the state variable w to the auxiliary stencil wAdj    *
!     * used by the Tapenade diferentiated routines. It takes into     *
!     * account whether or not the stencil is centered close to a      *
!     * physical block face (not an internal boundary created by block *
!     * splitting) since those do not have halo nodes.                 *
!     *                                                                *
!     * It is assumed that the pointers in blockPointers have already  *
!     * been set.                                                      *
!     *                                                                *
!     ******************************************************************
!
      use blockPointers   ! w, il, jl, kl
!      use indices         ! nw
      use flowVarRefState
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: iCell, jCell, kCell
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
                                                     intent(out) :: wAdj
      real(kind=realType), dimension(-2:3,-2:3,-2:3,3), &
                                                     intent(out) :: xAdj

!
!     Local variables.
!
      integer(kind=intType) :: ii, jj, kk, i1, j1, k1, i2, j2, k2, l
      integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

      ! Initialize the auxiliary array wAdj 
      do l=1,nw
        do kk=-2,2
          do jj=-2,2
            do ii=-2,2
              wAdj(ii,jj,kk,l) = 0.0
            enddo
          enddo
        enddo
      enddo

      ! Initialize the auxiliary array xAdj 
      do l=1,3
        do kk=-2,3
          do jj=-2,3
            do ii=-2,3
              xAdj(ii,jj,kk,l) = 0.0
            enddo
          enddo
        enddo
      enddo

      ! Copy the wAdj from w
      do l=1,nw
        do kk=-2,2
          do jj=-2,2
            do ii=-2,2
              wAdj(ii,jj,kk,l) = w(iCell+ii, jCell+jj, kCell+kk,l)
            enddo
          enddo
        enddo
      enddo


      ! Copy xAdj from x

      iStart=-2; iEnd=3
      jStart=-2; jEnd=3
      kStart=-2; kEnd=3

!!$      ! Special care needs to be done for subfaces. 
!!$      ! There're no points for -3 and 2 indices
!!$
!!$      if(iCell==2) iStart=-2; if(iCell==il) iEnd=1
!!$      if(jCell==2) jStart=-2; if(jCell==jl) jEnd=1
!!$      if(kCell==2) kStart=-2; if(kCell==kl) kEnd=1

      do l=1,3
        do kk=kStart,kEnd
          do jj=jStart,jEnd
            do ii=iStart,iEnd
              xAdj(ii,jj,kk,l) = x(iCell+ii, jCell+jj, kCell+kk,l)
            enddo
          enddo
        enddo
      enddo

    end subroutine copyADjointStencil

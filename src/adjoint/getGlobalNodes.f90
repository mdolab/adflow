!
!     ******************************************************************
!     *                                                                *
!     * File:          getGlobalNodes.f90                              *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 06-19-2007                                      *
!     * Last modified: 06-19-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine getGlobalNodes(blocknum,idim,jdim,kdim,nodenumbers)
!
!     ******************************************************************
!     *                                                                *
!     * Retrieve the global node numbering that is used to assemble    *
!     * the adjoint system of equations. Use this routine to ensure    *
!     * consitancy between python and fortran                          *
!     *                                                                *
!     * The nodes are numbered according to the following sequence:    *
!     *                                                                *
!     * loop processor = 1, nProc                                      *
!     *   loop domain = 1, nDom                                        *
!     *     loop k = 1, kl                                             *
!     *       loop j = 1, jl                                           *
!     *         loop i = 1, il                                         *
!     *                                                                *
!     * Only the onwned nodes are numbered, meaning i/j/k span from 1  *
!     * to il/jl/kl. The halo nodes receive the numbering from the     *
!     * neighboring block that owns them.                              *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars ! nNodesGlobal, nNodesLocal, nOffsetLocal
      use blockpointers
      use communication
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: blocknum
      integer(kind=intType), intent(in) :: idim,jdim,kdim
      integer(kind=intType), intent(out):: nodenumbers(idim,jdim,kdim)
!
!     Local variables.
!
      integer(kind=intType) :: nn, i, j, k
      integer :: ierr

      integer(kind=intType) ::level = 1, sps = 1

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution                                                *
!     *                                                                *
!     ******************************************************************
!
      !print *,'in getGlobalNodes'
      call setPointersAdj(blocknum,level,sps)

      !print *,'shapes',blocknum,shape(nodenumbers),shape(globalnode),'i',ie,ib,il,'j',je,jb,jl,'k',ke,kb,kl
      nodenumbers(:,:,:) = globalNode(:,:,:)
      
!!$      do i = 1,idim
!!$         do j = 1,jdim
!!$            do k = 1,kdim
!!$               print *,'node numbers', nodenumbers(i,j,k), globalNode(i,j,k),i,j,k
!!$            enddo
!!$         enddo
!!$      enddo

      end subroutine getGlobalNodes

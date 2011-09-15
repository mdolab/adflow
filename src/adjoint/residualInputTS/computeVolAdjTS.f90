!
!      ******************************************************************
!      *                                                                *
!      * File:          computeVolTS.f90                                *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 09-15-2011                                      *
!      * Last modified: 09-15-2011                                      *
!      *                                                                *
!      ******************************************************************
!

       subroutine computeVolTS(xAdj,volAdj,iCell,jCell,kCell,&
            nn,level,sps)
!
!      ******************************************************************
!      *                                                                *
!      * computes the volume necessary for the given grid level for all *
!      * spectral solutions. First the volumes are                      *
!      * computed assuming that the block is right handed. Then the     *
!      * number of positive and negative volumes are determined. If all *
!      * volumes are positive the block is indeed right handed; if all  *
!      * volumes are negative the block is left handed and both the     *
!      * volumes and the normals must be negated (for the normals this  *
!      * is done by the introduction of fact, which is either -0.5 or   *
!      * 0.5); if there are both positive and negative volumes the mesh *
!      * is not valid.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use inputTimeSpectral !nTimeIntervalsSpectral
       use section
       use constants

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(in) :: xAdj
       real(kind=realType),dimension(nTimeIntervalsSpectral) :: volAdj

       integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,level,sps

!
!      Local parameter.
!
       real(kind=realType), parameter :: thresVolume = 1.e-2_realType
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, n, m, l,iSt,iEn,jSt,jEn,jj,kk
       integer(kind=intType) :: iSBeg,iSEnd,jSBeg,jSEnd,kSBeg,kSEnd
       integer(kind=intType) :: iBBeg,iBEnd,jBBeg,jBEnd,kBBeg,kBEnd
       integer(kind=intType) :: iRBeg,iREnd,jRBeg,jREnd,kRBeg,kREnd
       integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
       integer(kind=intType) :: mm, nTime
       integer(kind=intType) :: nVolBad,   nVolBadGlobal
       integer(kind=intType) :: nVolNeg,   nVolPos

       real(kind=realType) :: fact, mult
       real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6

       real(kind=realType), dimension(3) :: v1, v2

       real(kind=realType), dimension(-2:2,-2:2,3) :: ss

       character(len=10) :: integerString

       logical :: checkK, checkJ, checkI, checkAll
       logical :: badVolume

       logical :: volumeIsNeg, iOverlap, jOverlap, kOverlap, secondHalo
       
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************

!
!      **************************************************************
!      *                                                            *
!      * Volume computation
!      *                                                            *
!      **************************************************************
!

       ! Initialize the number of positive and negative volumes for
       ! this block to 0.Needed to catch right/lefthanded blocks
       
       nVolNeg = 0
       nVolPos = 0


       ! Compute the volumes. The hexahedron is split into 6 pyramids
       ! whose volumes are computed. The volume is positive for a
       ! right handed block.
       ! Initialize the volumes to zero. The reasons is that the second
       ! level halo's must be initialized to zero and for convenience
       ! all the volumes are set to zero.
       
       k=0  ; j=0  ; i=0
       n=k-1; m=j-1; l=i-1;

       checkAll = .true. ! always check the volume of changed cell

       ! Compute the coordinates of the center of gravity.

       xp = eighth*(xAdj(i,j,k,1,sps) + xAdj(i,m,k,1,sps) &
            +         xAdj(i,m,n,1,sps) + xAdj(i,j,n,1,sps) &
            +         xAdj(l,j,k,1,sps) + xAdj(l,m,k,1,sps) &
            +         xAdj(l,m,n,1,sps) + xAdj(l,j,n,1,sps))
       yp = eighth*(xAdj(i,j,k,2,sps) + xAdj(i,m,k,2,sps) &
            +         xAdj(i,m,n,2,sps) + xAdj(i,j,n,2,sps) &
            +         xAdj(l,j,k,2,sps) + xAdj(l,m,k,2,sps) &
            +         xAdj(l,m,n,2,sps) + xAdj(l,j,n,2,sps))
       zp = eighth*(xAdj(i,j,k,3,sps) + xAdj(i,m,k,3,sps) &
            +         xAdj(i,m,n,3,sps) + xAdj(i,j,n,3,sps) &
            +         xAdj(l,j,k,3,sps) + xAdj(l,m,k,3,sps) &
            +         xAdj(l,m,n,3,sps) + xAdj(l,j,n,3,sps))

       ! Compute the volumes of the 6 sub pyramids. The
       ! arguments of volpym must be such that for a (regular)
       ! right handed hexahedron all volumes are positive.
       
       call volpym3(xAdj(i,j,k,1,sps), xAdj(i,j,k,2,sps), xAdj(i,j,k,3,sps), &
            xAdj(i,j,n,1,sps), xAdj(i,j,n,2,sps), xAdj(i,j,n,3,sps), &
            xAdj(i,m,n,1,sps), xAdj(i,m,n,2,sps), xAdj(i,m,n,3,sps), &
            xAdj(i,m,k,1,sps), xAdj(i,m,k,2,sps), xAdj(i,m,k,3,sps),xp,yp,zp,vp1)
       
       call volpym3(xAdj(l,j,k,1,sps), xAdj(l,j,k,2,sps), xAdj(l,j,k,3,sps), &
            xAdj(l,m,k,1,sps), xAdj(l,m,k,2,sps), xAdj(l,m,k,3,sps), &
            xAdj(l,m,n,1,sps), xAdj(l,m,n,2,sps), xAdj(l,m,n,3,sps), &
            xAdj(l,j,n,1,sps), xAdj(l,j,n,2,sps), xAdj(l,j,n,3,sps),xp,yp,zp,vp2)
       
       call volpym3(xAdj(i,j,k,1,sps), xAdj(i,j,k,2,sps), xAdj(i,j,k,3,sps), &
            xAdj(l,j,k,1,sps), xAdj(l,j,k,2,sps), xAdj(l,j,k,3,sps), &
            xAdj(l,j,n,1,sps), xAdj(l,j,n,2,sps), xAdj(l,j,n,3,sps), &
            xAdj(i,j,n,1,sps), xAdj(i,j,n,2,sps), xAdj(i,j,n,3,sps),xp,yp,zp,vp3)
       
       call volpym3(xAdj(i,m,k,1,sps), xAdj(i,m,k,2,sps), xAdj(i,m,k,3,sps), &
            xAdj(i,m,n,1,sps), xAdj(i,m,n,2,sps), xAdj(i,m,n,3,sps), &
            xAdj(l,m,n,1,sps), xAdj(l,m,n,2,sps), xAdj(l,m,n,3,sps), &
            xAdj(l,m,k,1,sps), xAdj(l,m,k,2,sps), xAdj(l,m,k,3,sps),xp,yp,zp,vp4)
       
       call volpym3(xAdj(i,j,k,1,sps), xAdj(i,j,k,2,sps), xAdj(i,j,k,3,sps), &
            xAdj(i,m,k,1,sps), xAdj(i,m,k,2,sps), xAdj(i,m,k,3,sps), &
            xAdj(l,m,k,1,sps), xAdj(l,m,k,2,sps), xAdj(l,m,k,3,sps), &
            xAdj(l,j,k,1,sps), xAdj(l,j,k,2,sps), xAdj(l,j,k,3,sps),xp,yp,zp,vp5)
       
       call volpym3(xAdj(i,j,n,1,sps), xAdj(i,j,n,2,sps), xAdj(i,j,n,3,sps), &
            xAdj(l,j,n,1,sps), xAdj(l,j,n,2,sps), xAdj(l,j,n,3,sps), &
            xAdj(l,m,n,1,sps), xAdj(l,m,n,2,sps), xAdj(l,m,n,3,sps), &
            xAdj(i,m,n,1,sps), xAdj(i,m,n,2,sps), xAdj(i,m,n,3,sps),xp,yp,zp,vp6)

       ! Set the volume to 1/6 of the sum of the volumes of the
       ! pyramid. Remember that volpym computes 6 times the
       ! volume.
       

       volAdj(sps) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)
       
       ! Check the volume and update the number of positive
       ! and negative volumes if needed.
       
       if( checkAll ) then
          
          ! Update either the number of negative or positive
          ! volumes. Negative volumes should only occur for left
          ! handed blocks. This is checked later.
          ! Set the logical volumeIsNeg accordingly.
               
          if(volAdj(sps) < zero) then
             nVolNeg            = nVolNeg + 1
             volumeIsNeg = .true.
          else
             nVolPos            = nVolPos + 1
             volumeIsNeg = .false.
          endif
          
          ! terminate if negative volume is located
          !if(volumeIsNeg) &
          !     write(*,*)"VOLUME NEGATIVE",voladj(sps),nbkglobal,iCell,jCell,kCell
!            call terminate("negative volume located")

          ! Set the threshold for the volume quality.
          
          fact = thresVolume*abs(voladj(sps))
          
          ! Check the quality of the volume.
          
          badVolume = .false.

          if(vp1*volAdj(sps) < zero .and. &
               abs(vp1)       > fact) badVolume = .true.
          if(vp2*volAdj(sps) < zero .and. &
               abs(vp2)       > fact) badVolume = .true.
          if(vp3*volAdj(sps) < zero .and. &
               abs(vp3)       > fact) badVolume = .true.
          if(vp4*volAdj(sps) < zero .and. &
               abs(vp4)       > fact) badVolume = .true.
          if(vp5*volAdj(sps) < zero .and. &
               abs(vp5)       > fact) badVolume = .true.
          if(vp6*volAdj(sps) < zero .and. &
               abs(vp6)       > fact) badVolume = .true.
          
          ! Update nVolBad if this is a bad volume.
          
          if( badVolume .and. myID==0) then
             write(*,'(a)')"bad quality volumes found"
             write(*,'(a)')"Computation will continue, but be aware of this"
          endif
       
       end if
          ! Set the volume to the absolute value.
          
       
       volAdj(sps) = abs(volAdj(sps))
          
!!$           ! Some additional safety stuff for halo volumes.
!!$
!!$           do k=2,kl
!!$             do j=2,jl
!!$               if(vol(1, j,k) <= eps) vol(1, j,k) = vol(2, j,k)
!!$               if(vol(ie,j,k) <= eps) vol(ie,j,k) = vol(il,j,k)
!!$             enddo
!!$           enddo
!!$
!!$           do k=2,kl
!!$             do i=1,ie
!!$               if(vol(i,1, k) <= eps) vol(i,1, k) = vol(i,2, k)
!!$               if(vol(i,je,k) <= eps) vol(i,je,k) = vol(i,jl,k)
!!$             enddo
!!$           enddo
!!$
!!$           do j=1,je
!!$             do i=1,ie
!!$               if(vol(i,j,1)  <= eps) vol(i,j,1)  = vol(i,j,2)
!!$               if(vol(i,j,ke) <= eps) vol(i,j,ke) = vol(i,j,kl)
!!$             enddo
!!$           enddo

       ! Determine the orientation of the block. For the fine level
       ! this is based on the number of positive and negative
       ! volumes; on the coarse levels the corresponding fine level
       ! value is taken. If both positive and negative volumes are
       ! present it is assumed that the block was intended to be
       ! right handed. The code will terminate later on anyway.
       
       if(level == 1) then
          if(nVolPos == 0) then       ! Left handed block.
             rightHanded = .false.
          else                        ! Right handed (or bad) block.
             rightHanded = .true.
          endif
       else
          print *,'ADjoint not setup on lower levels'
          stop
       endif

  
     end subroutine computeVolTS


!      ==================================================================
!      ================================================================

  subroutine volpym3(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,xp,yp,zp,vp)
!
!        ****************************************************************
!        *                                                              *
!        * volpym computes 6 times the volume of a pyramid. Node p,     *
!        * whose coordinates are set in the subroutine metric itself,   *
!        * is the top node and a-b-c-d is the quadrilateral surface.    *
!        * It is assumed that the cross product vCa * vDb points in     *
!        * the direction of the top node. Here vCa is the diagonal      *
!        * running from node c to node a and vDb the diagonal from      *
!        * node d to node b.                                            *
!        *                                                              *
!        ****************************************************************
!
       use precision
       use constants
       implicit none

!
!        Subroutine arguments.
!
       real(kind=realType), intent(in) :: xa, ya, za, xb, yb, zb
       real(kind=realType), intent(in) :: xc, yc, zc, xd, yd, zd
       real(kind=realType), intent(in) :: xp, yp, zp

       real(kind=realType), intent(out) :: vp
       
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
       vp = (xp - fourth*(xa + xb  + xc + xd))              &
          * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
             (yp - fourth*(ya + yb  + yc + yd))              &
          * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
             (zp - fourth*(za + zb  + zc + zd))              &
          * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

     end subroutine volpym3



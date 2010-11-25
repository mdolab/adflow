!
!     ******************************************************************
!     *                                                                *
!     * File:          copyADjointStencil.f90                          *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     *                Seongim Choi                                    *
!     * Starting date: 08-03-2006                                      *
!     * Last modified: 11-30-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
           betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
           nn,level,sps,pointRefAdj,rotPointAdj,&
           prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
           rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
           murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
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
      use flowVarRefState  !timeref,nw
      use inputPhysics
      use inputTimeSpectral !nTimeIntervalsSpectral
      use cgnsgrid    !cgnsdoms 
      use inputMotion     ! rotPoint
      implicit none

!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: iCell, jCell, kCell
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
                                                     intent(out) :: wAdj
!      real(kind=realType), dimension(-2:3,-2:3,-2:3,3), &
!                                                     intent(out) :: xAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), &
                                                     intent(out) :: xAdj

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,MachCoefAdj,machGridAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
      integer(kind=intType)::liftIndex

      real(kind=realType), dimension(3),intent(out) ::rotRateAdj,rotCenterAdj
      real(kind=realType), dimension(3),intent(out) ::pointRefAdj,rotPointAdj
      real(kind=realType), dimension(2,2,2,3,nTimeIntervalsSpectral),intent(out) ::xBlockCornerAdj

      integer(kind=intType) :: nn,level,sps

!
!     Local variables.
!
      integer(kind=intType) :: sps2
      integer(kind=intType) :: ii, jj, kk, i1, j1, k1, i2, j2, k2, l,j
      integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
      integer(kind=intType) :: ind(3,56),cellstodo,iiCell
      
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

!      call five_point_cell_stencil(ind,cellstodo)

      ! Zero w and x
      wAdj = 0.0
      xAdj = 0.0

!       do sps2 = 1,nTimeIntervalsSpectral
!          do l=1,nw
!             do iiCell = 1,cellsToDo
!                ii = ind(1,iiCell)
!                jj = ind(2,iiCell)
!                kk = ind(3,iiCell)
!                wAdj(ii,jj,kk,l,sps2) = flowDoms(nn,level,sps2)%w(iCell+ii,jCell+jj,kCell+kk,l)
!             end do
!          end do
!       end do


      ! Copy the wAdj from w
      do sps2 = 1,nTimeIntervalsSpectral
         call setPointers(nn,level,sps2)
         do l=1,nw
            do kk=-2,2
               do jj=-2,2
                  do ii=-2,2
                     !print *,'wadj',wAdj(ii,jj,kk,l), iCell+ii, jCell+jj, kCell+kk
                     wAdj(ii,jj,kk,l,sps2) = w(iCell+ii, jCell+jj, kCell+kk,l)
                  enddo
               enddo
            enddo
         enddo
         call setPointers(nn,level,sps)
      end do

      ! Copy xAdj from x
      

!      call five_pt_node_stencil_all(icell,jcell,kcell,ind,CellsToDo)
      
 !      do sps2 = 1,nTimeIntervalsSpectral
!           do l=1,3
!             do iiCell = 1,cellsToDo
!                ii = ind(1,iiCell)
!                jj = ind(2,iiCell)
!                kk = ind(3,iiCell)
               
!                xAdj(ii,jj,kk,l,sps2) = flowDoms(nn,level,sps2)%x(iCell+ii, jCell+jj, kCell+kk,l)
!             end do
!          end do
!       end do


      iStart=-3; iEnd=2
      jStart=-3; jEnd=2
      kStart=-3; kEnd=2

      ! Special care needs to be done for subfaces. 
      ! There're no points for -3 and 2 indices

      if(iCell==2) iStart=-2; if(iCell==il) iEnd=1
      if(jCell==2) jStart=-2; if(jCell==jl) jEnd=1
      if(kCell==2) kStart=-2; if(kCell==kl) kEnd=1
     do sps2 = 1,nTimeIntervalsSpectral
        !call setPointers(nn,level,sps2)
         do l=1,3
            do kk=kStart,kEnd
               do jj=jStart,jEnd
                  do ii=iStart,iEnd
                     xAdj(ii,jj,kk,l,sps2) = flowDoms(nn,level,sps2)%x(iCell+ii, jCell+jj, kCell+kk,l)
                  enddo
               enddo
            enddo
         enddo
         !call setPointers(nn,level,sps)
      end do
    
      do sps2 = 1,nTimeIntervalsSpectral
         !call setPointers(nn,level,sps2)     
         xBlockCornerAdj(1,1,1,:,sps2) = flowDoms(nn,level,sps2)%x(1 , 1, 1,:)
         xBlockCornerAdj(2,1,1,:,sps2) = flowDoms(nn,level,sps2)%x(il, 1, 1,:)
         xBlockCornerAdj(1,2,1,:,sps2) = flowDoms(nn,level,sps2)%x(1 ,jl, 1,:)
         xBlockCornerAdj(2,2,1,:,sps2) = flowDoms(nn,level,sps2)%x(il,jl, 1,:)
         xBlockCornerAdj(1,2,2,:,sps2) = flowDoms(nn,level,sps2)%x(1 ,jl,kl,:)
         xBlockCornerAdj(2,2,2,:,sps2) = flowDoms(nn,level,sps2)%x(il,jl,kl,:)
         xBlockCornerAdj(1,1,2,:,sps2) = flowDoms(nn,level,sps2)%x(1 , 1,kl,:)
         xBlockCornerAdj(2,1,2,:,sps2) = flowDoms(nn,level,sps2)%x(il, 1,kl,:)
         !call setPointers(nn,level,sps)
      end do
      
      MachAdj = Mach
      MachCoefAdj = MachCoef
      MachGridAdj = MachGrid

      call getDirAngle(velDirFreestream,liftDirection,liftIndex,alphaAdj,betaAdj)


!      call getDirAngle(velDirFreestream,velDirFreestream,liftIndex,alphaAdj,betaAdj)
      !call getDirAngle(velDirFreestream(1), velDirFreestream(2),&
      !     velDirFreestream(3), alphaAdj, betaAdj)

      prefAdj = pRef
      rhorefAdj = rhoref
      pinfdimAdj = pinfdim
      rhoinfdimAdj = rhoinfdim
      rhoinfAdj = rhoinf
      pinfAdj = pInf
      murefAdj = muref
      timerefAdj = timeref
      pInfCorrAdj = pInfCorr
      
      ! Store the rotation center and determine the 
      ! nonDimensional rotation rate of this block. As the
      ! reference length is 1 timeRef == 1/uRef and at the end
      ! the nonDimensional velocity is computed.

      j = nbkGlobal
      
      rotCenterAdj = cgnsDoms(j)%rotCenter
      rotRateAdj   = timeRef*cgnsDoms(j)%rotRate
!      rotRateAdj   = cgnsDoms(j)%rotRate
      pointRefAdj = pointRef
      rotPointAdj = rotPoint
    end subroutine copyADjointStencil




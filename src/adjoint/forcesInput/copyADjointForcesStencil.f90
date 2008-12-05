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
      subroutine copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
           rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
           timerefAdj,pInfCorrAdj,nn,level,sps,liftIndex)

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
      use inputPhysics    ! Mach,veldirfreestream
      use cgnsgrid        ! cgnsdoms
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: nn,level,sps    
      real(kind=realType), dimension(0:ie,0:je,0:ke,3), intent(out) :: xAdj
      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw), intent(out) :: wAdj

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,MachCoefAdj,MachGridAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
      integer(kind=intType)::liftIndex

      real(kind=realType), dimension(3),intent(out) ::rotRateAdj,rotCenterAdj

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

      MachAdj = Mach
      MachCoefAdj = MachCoef
      MachGridAdj = MachGrid
      call getDirAngle(velDirFreestream,LiftDirection,liftIndex,alphaAdj,betaAdj)
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
      
      end subroutine copyADjointForcesStencil

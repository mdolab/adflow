!
!      ******************************************************************
!      *                                                                *
!      * File:          normalVelocitiesAdj.f90                         *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 10-25-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine normalVelocitiesAllLevelsAdj(sps,iCell, jCell, kCell,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj,siAdj, sjAdj, skAdj,rFaceAdj,nn,level,sps2)
!
!      ******************************************************************
!      *                                                                *
!      * normalVelocitiesAllLevels computes the normal grid             *
!      * velocities of some boundary faces of the moving blocks for     *
!      * spectral mode sps. All grid levels from ground level to the    *
!      * coarsest level are considered.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use iteration
       use inputTimeSpectral !nIntervalTimespectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn,level,sps,sps2
       integer(kind=intType), intent(in) :: iCell, jCell, kCell
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(in) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(nBocos,-2:2,-2:2,nTimeIntervalsSpectral), intent(out) :: rFaceAdj
       
!
!      Local variables.
!

       integer(kind=intType) :: nLevels, mm
       integer(kind=intType) :: i, j,l,m

       integer(kind=intType) :: iSt,iEn,jSt,jEn,ii,jj
       integer(kind=intType) :: iSBeg,iSEnd,jSBeg,jSEnd,kSBeg,kSEnd
       integer(kind=intType) :: iBBeg,iBEnd,jBBeg,jBEnd,kBBeg,kBEnd
       integer(kind=intType) :: iRBeg,iREnd,jRBeg,jREnd,kRBeg,kREnd

       real(kind=realType) :: weight, mult

       real(kind=realType), dimension(-2:2,-2:2)::sFaceAdj
       real(kind=realType), dimension(-3:2,-3:2,3):: ss

       logical :: secondHalo,computeBC
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! since rface normally isn't allocate for all BC's set to zero by 
       ! default, then let it get filled in.

       rFaceAdj(:,:,:,sps2) = zero

!!Governed Outside computeRAdj!
!!$       ! Loop over the number of grid levels, starting at groundLevel,
!!$       ! the currently finest mesh.
!!$
!!$       nLevels = ubound(flowDoms,2)
!!$       levelLoop: do level=groundLevel,nLevels
!!$
!!$         ! Loop over the number of local blocks.
!!$
!!$         domains: do nn=1,nDom
!!$
!!$           ! Set the pointers for this block.
!!$
!!$           call setPointersAdj(nn, level, sps)

           ! Check for a moving block. As it is possible that in a
           ! multidisicplinary environment additional grid velocities
           ! are set, the test should be done on addGridVelocities
           ! and not on blockIsMoving.

           testMoving: if( addGridVelocities ) then
!
!            ************************************************************
!            *                                                          *
!            * Determine the normal grid velocities of the boundaries.  *
!            * As these values are based on the unit normal. A division *
!            * by the length of the normal is needed.                   *
!            * Furthermore the boundary unit normals are per definition *
!            * outward pointing, while on the iMin, jMin and kMin       *
!            * boundaries the face normals are inward pointing. This    *
!            * is taken into account by the factor mult.                *
!            *                                                          *
!            ************************************************************
!
             ! Loop over the boundary subfaces.

             bocoLoop: do mm=1,nBocos

                call checkOverlapAdj(mm,icell,jcell,kcell,isbeg,jsbeg,&
                     ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
                     jbend,kbend,computeBC)

                if (computeBC) then
                   !There is overlap....

                   iRBeg = max(iSBeg, iBBeg); iREnd = min(iSEnd, iBEnd)
                   jRBeg = max(jSBeg, jBBeg); jREnd = min(jSEnd, jBEnd)
                   kRBeg = max(kSBeg, kBBeg); kREnd = min(kSEnd, kBEnd)
                   
                   iSt=-2; iEn=2
                   jSt=-2; jEn=2


               ! Check whether rFace is allocated.
               !testAssoc: if( associated(BCData(mm)%rFace) ) then
                 testrface: if(BCType(mm) == FarField .or. BCType(mm) == EulerWall) then

                 ! Determine the block face on which the subface is
                 ! located and set some variables accordingly.

                 select case (BCFaceID(mm))

                   case (iMin)
                      if(jCell==2)  iSt=-1
                      if(jCell==jl) iEn=1
                      if(kCell==2)  jSt=-1
                      if(kCell==kl)  jEn=1
                      l=jcell
                      m=kcell
                      secondHalo = .true.
                      if(iRBeg == iREnd) secondHalo = .false.
                      if(secondHalo)then
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = siAdj(-1,iSt:iEn,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn)= sFaceIAdj(-1,iSt:iEn,jSt:jEn,sps2)
                      else 
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = siAdj(-2,iSt:iEn,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn)= sFaceIAdj(-2,iSt:iEn,jSt:jEn,sps2)
                      end if
                   
                   case (iMax)
                      if(jCell==2)  iSt=-1
                      if(jCell==jl) iEn=1
                      if(kCell==2)  jSt=-1
                      if(kCell==kl)  jEn=1
                      l=jcell
                      m=kcell
                      secondHalo = .true.
                      if(iRBeg == iREnd) secondHalo = .false.
                      if(secondHalo) then
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = siAdj(0,iSt:iEn,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceIAdj(0,iSt:iEn,jSt:jEn,sps2)
                      else
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = siAdj(1,iSt:iEn,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceIAdj(1,iSt:iEn,jSt:jEn,sps2)
                      end if
                   
                   case (jMin)
                      if(iCell==2)  iSt=-1
                      if(iCell==il) iEn=1
                      if(kCell==2)  jSt=-1
                      if(kCell==kl) jEn=1
                      l=icell
                      m=kcell
                      secondHalo = .true.
                      if(jRBeg == jREnd) secondHalo = .false.
                      if(secondHalo) then
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-1,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceJAdj(iSt:iEn,-1,jSt:jEn,sps2)
                      else
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-2,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceJAdj(iSt:iEn,-2,jSt:jEn,sps2)
                      end if

                   case (jMax)
                      if(iCell==2)  iSt=-1
                      if(iCell==il) iEn=1
                      if(kCell==2)  jSt=-1
                      if(kCell==kl) jEn=1
                      l=icell
                      m=kcell
                      secondHalo = .true.
                      if(jRBeg == jREnd) secondHalo = .false.
                      if(secondHalo) then
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,0,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceJAdj(iSt:iEn,0,jSt:jEn,sps2)
                      else
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,1,jSt:jEn,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceJAdj(iSt:iEn,1,jSt:jEn,sps2)
                      end if
                     
                   case (kMin)
                      if(iCell==2)  iSt=-1
                      if(iCell==il) iEn=1
                      if(jCell==2)  jSt=-1
                      if(jCell==jl) jEn=1
                      l=icell
                      m=jcell
                      secondHalo = .true.
                      if(kRBeg == kREnd) secondHalo = .false.
                      if(secondHalo) then
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-1,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceKAdj(iSt:iEn,jSt:jEn,-1,sps2)
                      else
                         mult = -one
                         ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-2,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceKAdj(iSt:iEn,jSt:jEn,-2,sps2)
                      end if
                     
                   case (kMax)
                      if(iCell==2)  iSt=-1
                      if(iCell==il) iEn=1
                      if(jCell==2)  jSt=-1
                      if(jCell==jl) jEn=1   
                      l=icell
                      m=jcell
                      secondHalo = .true.
                      if(kRBeg == kREnd) secondHalo = .false.
                      if(secondHalo) then
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,0,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceKAdj(iSt:iEn,jSt:jEn,0,sps2)
                      else
                         mult = one
                         ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,1,:,sps2)
                         sFaceAdj(iSt:iEn,jSt:jEn) = sFaceKAdj(iSt:iEn,jSt:jEn,1,sps2)
                      end if
                      
                   end select

                 ! Loop over the faces of the subface.
                 do jj=jSt,jEn
                    do ii=iSt,iEn
                 !do j=jcBeg, jcEnd
                 !  do i=icBeg, icEnd

                     ! Compute the inverse of the length of the normal
                     ! vector and possibly correct for inward pointing.
                       if( ss(ii,jj,1)**2>zero .or. ss(ii,jj,2)**2>zero&
                            .or. ss(ii,jj,3)**2>zero)then
                          weight = sqrt(ss(ii,jj,1)**2 + ss(ii,jj,2)**2 &
                               +      ss(ii,jj,3)**2)
                          !if(weight > zero) weight = mult/weight
                          weight = mult/weight
                          ! Compute the normal velocity based on the outward
                          ! pointing unit normal.
                          
                          !BCData(mm)%rFace(i,j) = weight*sFace(i,j)
                          rFaceAdj(mm,ii,jj,sps2) = weight*sFaceAdj(ii,jj)
!!$                     if (abs(BCData(mm)%rFace(l+ii,m+jj)-rFaceAdj(mm,ii,jj))>1e-16)then
!!$                        print *,'indices',mm,ii,jj,l,m
!!$                        print *,'rface',BCData(mm)%rFace(l+ii,m+jj),rFaceAdj(mm,ii,jj),BCData(mm)%rFace(l+ii,m+jj)-rFaceAdj(mm,ii,jj)
!!$                     endif
                       else
                          rFaceAdj(mm,ii,jj,sps2)= zero
                       endif

                   enddo
                 enddo
              endif testrface
 
           endif
           enddo bocoLoop

           else testMoving

             ! Block is not moving. Loop over the boundary faces and set
             ! the normal grid velocity to zero if allocated.

             do mm=1,nBocos
                !if( associated(BCData(mm)%rFace) ) &
                !    BCData(mm)%rFace = zero 
                if(BCType(mm) == FarField .or. BCType(mm) == EulerWall) then
                   rFaceAdj(mm,:,:,sps2) = zero
                endif
             enddo

           endif testMoving
!         enddo domains

!       enddo levelLoop

       end subroutine normalVelocitiesAllLevelsAdj

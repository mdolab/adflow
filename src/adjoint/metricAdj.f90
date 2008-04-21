!
!      ******************************************************************
!      *                                                                *
!      * File:          metricAdj.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi                                    *
!      * Starting date: 12-15-2007                                      *
!      * Last modified: 12-26-2007                                      *
!      *                                                                *
!      ******************************************************************
!

       subroutine metricAdj(xAdj,siAdj,sjAdj,skAdj,volAdj,normAdj, &
                            iCell,jCell,kCell)
!
!      ******************************************************************
!      *                                                                *
!      * metric computes the face normals and the volume for the given  *
!      * grid level for all spectral solutions. First the volumes are   *
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
       use inputTimeSpectral
       use section
       use constants

       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(-2:3,-2:3,-2:3,3), intent(in) :: xAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,3), intent(out) :: siAdj, sjAdj, skAdj
       !real(kind=realType), dimension(0:0,0:0,0:0), intent(out) :: volAdj
       real(kind=realType), intent(out) :: volAdj
       real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(out) :: normAdj
       integer(kind=intType), intent(in) :: iCell, jCell, kCell

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
       integer(kind=intType) :: mm, sps, nTime
       integer(kind=intType) :: nVolBad,   nVolBadGlobal

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

!       print *,'in metric',xAdj,siAdj,sjAdj,skAdj,volAdj,normAdj, &
!            iCell,jCell,kCell

       ! Some initialization for siAdj,sjAdj,skAdj,normAdj 
       ! Volume needs only one stencil so it does not need initialization

       siAdj = zero; sjAdj = zero; skAdj = zero
       normAdj = zero

       ! print *,'in metric',xAdj, &!,siAdj,sjAdj,skAdj,volAdj,normAdj, &
       !     iCell,jCell,kCell
       

!
!
!      **************************************************************
!      *                                                            *
!      * Volume computation
!      *                                                            *
!      **************************************************************
!
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

       xp = eighth*(xAdj(i,j,k,1) + xAdj(i,m,k,1) &
            +         xAdj(i,m,n,1) + xAdj(i,j,n,1) &
            +         xAdj(l,j,k,1) + xAdj(l,m,k,1) &
            +         xAdj(l,m,n,1) + xAdj(l,j,n,1))
       yp = eighth*(xAdj(i,j,k,2) + xAdj(i,m,k,2) &
            +         xAdj(i,m,n,2) + xAdj(i,j,n,2) &
            +         xAdj(l,j,k,2) + xAdj(l,m,k,2) &
            +         xAdj(l,m,n,2) + xAdj(l,j,n,2))
       zp = eighth*(xAdj(i,j,k,3) + xAdj(i,m,k,3) &
            +         xAdj(i,m,n,3) + xAdj(i,j,n,3) &
            +         xAdj(l,j,k,3) + xAdj(l,m,k,3) &
            +         xAdj(l,m,n,3) + xAdj(l,j,n,3))
       !print *,'xp',xp,yp,zp
       ! Compute the volumes of the 6 sub pyramids. The
       ! arguments of volpym must be such that for a (regular)
       ! right handed hexahedron all volumes are positive.
       
       call volpym2(xAdj(i,j,k,1), xAdj(i,j,k,2), xAdj(i,j,k,3), &
            xAdj(i,j,n,1), xAdj(i,j,n,2), xAdj(i,j,n,3), &
            xAdj(i,m,n,1), xAdj(i,m,n,2), xAdj(i,m,n,3), &
            xAdj(i,m,k,1), xAdj(i,m,k,2), xAdj(i,m,k,3),xp,yp,zp,vp1)
       
       call volpym2(xAdj(l,j,k,1), xAdj(l,j,k,2), xAdj(l,j,k,3), &
            xAdj(l,m,k,1), xAdj(l,m,k,2), xAdj(l,m,k,3), &
            xAdj(l,m,n,1), xAdj(l,m,n,2), xAdj(l,m,n,3), &
            xAdj(l,j,n,1), xAdj(l,j,n,2), xAdj(l,j,n,3),xp,yp,zp,vp2)
       
       call volpym2(xAdj(i,j,k,1), xAdj(i,j,k,2), xAdj(i,j,k,3), &
            xAdj(l,j,k,1), xAdj(l,j,k,2), xAdj(l,j,k,3), &
            xAdj(l,j,n,1), xAdj(l,j,n,2), xAdj(l,j,n,3), &
            xAdj(i,j,n,1), xAdj(i,j,n,2), xAdj(i,j,n,3),xp,yp,zp,vp3)
       
       call volpym2(xAdj(i,m,k,1), xAdj(i,m,k,2), xAdj(i,m,k,3), &
            xAdj(i,m,n,1), xAdj(i,m,n,2), xAdj(i,m,n,3), &
            xAdj(l,m,n,1), xAdj(l,m,n,2), xAdj(l,m,n,3), &
            xAdj(l,m,k,1), xAdj(l,m,k,2), xAdj(l,m,k,3),xp,yp,zp,vp4)
       
       call volpym2(xAdj(i,j,k,1), xAdj(i,j,k,2), xAdj(i,j,k,3), &
            xAdj(i,m,k,1), xAdj(i,m,k,2), xAdj(i,m,k,3), &
            xAdj(l,m,k,1), xAdj(l,m,k,2), xAdj(l,m,k,3), &
            xAdj(l,j,k,1), xAdj(l,j,k,2), xAdj(l,j,k,3),xp,yp,zp,vp5)
       
       call volpym2(xAdj(i,j,n,1), xAdj(i,j,n,2), xAdj(i,j,n,3), &
            xAdj(l,j,n,1), xAdj(l,j,n,2), xAdj(l,j,n,3), &
            xAdj(l,m,n,1), xAdj(l,m,n,2), xAdj(l,m,n,3), &
            xAdj(i,m,n,1), xAdj(i,m,n,2), xAdj(i,m,n,3),xp,yp,zp,vp6)

       ! Set the volume to 1/6 of the sum of the volumes of the
       ! pyramid. Remember that volpym computes 6 times the
       ! volume.
       
!       volAdj(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)
       volAdj = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)
       !print *,'VolAdj',volAdj,checkAll
       ! Check the volume and update the number of positive
       ! and negative volumes if needed.
       
       if( checkAll ) then
          
          ! Update either the number of negative or positive
          ! volumes. Negative volumes should only occur for left
          ! handed blocks. This is checked later.
          ! Set the logical volumeIsNeg accordingly.
          
          !if(volAdj(i,j,k) < zero) then
          if(volAdj < zero) then
             volumeIsNeg = .true.
          else
             volumeIsNeg = .false.
          endif
          
          ! terminate if negative volume is located
          if(volumeIsNeg) &
               write(*,*)"VOLUME NEGATIVE"
!            call terminate("negative volume located")

          ! Set the threshold for the volume quality.
          
          !fact = thresVolume*abs(voladj(i,j,k))
          fact = thresVolume*abs(voladj)
          
          ! Check the quality of the volume.
          
          badVolume = .false.
!!$          if(vp1*volAdj(i,j,k) < zero .and. &
!!$               abs(vp1)       > fact) badVolume = .true.
!!$          if(vp2*volAdj(i,j,k) < zero .and. &
!!$               abs(vp2)       > fact) badVolume = .true.
!!$          if(vp3*volAdj(i,j,k) < zero .and. &
!!$               abs(vp3)       > fact) badVolume = .true.
!!$          if(vp4*volAdj(i,j,k) < zero .and. &
!!$               abs(vp4)       > fact) badVolume = .true.
!!$          if(vp5*volAdj(i,j,k) < zero .and. &
!!$               abs(vp5)       > fact) badVolume = .true.
!!$          if(vp6*volAdj(i,j,k) < zero .and. &
!!$               abs(vp6)       > fact) badVolume = .true.
          if(vp1*volAdj < zero .and. &
               abs(vp1)       > fact) badVolume = .true.
          if(vp2*volAdj < zero .and. &
               abs(vp2)       > fact) badVolume = .true.
          if(vp3*volAdj < zero .and. &
               abs(vp3)       > fact) badVolume = .true.
          if(vp4*volAdj < zero .and. &
               abs(vp4)       > fact) badVolume = .true.
          if(vp5*volAdj < zero .and. &
               abs(vp5)       > fact) badVolume = .true.
          if(vp6*volAdj < zero .and. &
               abs(vp6)       > fact) badVolume = .true.
          
          ! Update nVolBad if this is a bad volume.
          
          if( badVolume .and. myID==0) then
             write(*,'(a)')"bad quality volumes found"
             write(*,'(a)')"Computation will continue, but be aware of this"
          endif
       
       end if
          ! Set the volume to the absolute value.
          
       !volAdj(i,j,k) = abs(volAdj(i,j,k))
       volAdj = abs(volAdj)
       !print *,'absvol',volAdj
          
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


!
!          **************************************************************
!          *                                                            *
!          * Computation of the face normals in i-, j- and k-direction. *
!          * Formula's are valid for a right handed block; for a left   *
!          * handed block the correct orientation is obtained via fact. *
!          * The normals point in the direction of increasing index.    *
!          * The absolute value of fact is 0.5, because the cross       *
!          * product of the two diagonals is twice the normal vector.   *
!          *                                                            *
!          * Note that also the normals of the first level halo cells   *
!          * are computed. These are needed for the viscous fluxes.     *
!          *                                                            *
!          **************************************************************
!

!s           if( flowDoms(nn,level,sps)%rightHanded ) then
             fact =  half
!s           else
!s             fact = -half
!s           endif

           ! Projected areas of cell faces in the i direction.

          kStart=-2; kEnd=2
          jStart=-2; jEnd=2
          iStart=-2; iEnd=2

          if(iCell==il) iEnd=1

          if(jCell==2)  jStart=-1
          if(jCell==jl) jEnd=1 

          if(kCell==2) kStart=-1
          if(kCell==kl) kEnd=1

           do k=kStart,kEnd !-2,2
             n = k -1
             do j=jStart,jEnd !-2,2
               m = j -1
               do i=iStart,iEnd !-2,2

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = xAdj(i,j,n,1) - xAdj(i,m,k,1)
                 v1(2) = xAdj(i,j,n,2) - xAdj(i,m,k,2)
                 v1(3) = xAdj(i,j,n,3) - xAdj(i,m,k,3)

                 v2(1) = xAdj(i,j,k,1) - xAdj(i,m,n,1)
                 v2(2) = xAdj(i,j,k,2) - xAdj(i,m,n,2)
                 v2(3) = xAdj(i,j,k,3) - xAdj(i,m,n,3)

                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 siAdj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 siAdj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 siAdj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

               enddo
             enddo
           enddo
           
           ! Projected areas of cell faces in the j direction.
          kStart=-2; kEnd=2
          jStart=-2; jEnd=2
          iStart=-2; iEnd=2

          if(iCell==2)  iStart=-1
          if(iCell==il) iEnd=1

          if(jCell==jl) jEnd=1 

          if(kCell==2) kStart=-1
          if(kCell==kl) kEnd=1

           do k=kStart,kEnd !-2,2
             n = k -1
             do j=jStart,jEnd !-2,2
               do i=iStart,iEnd !-2,2
                 l = i -1

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = xAdj(i,j,n,1) - xAdj(l,j,k,1)
                 v1(2) = xAdj(i,j,n,2) - xAdj(l,j,k,2)
                 v1(3) = xAdj(i,j,n,3) - xAdj(l,j,k,3)

                 v2(1) = xAdj(l,j,n,1) - xAdj(i,j,k,1)
                 v2(2) = xAdj(l,j,n,2) - xAdj(i,j,k,2)
                 v2(3) = xAdj(l,j,n,3) - xAdj(i,j,k,3)

                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 sjAdj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 sjAdj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 sjAdj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
 
               enddo
             enddo
           enddo

           ! Projected areas of cell faces in the k direction.
           ! Projected areas of cell faces in the j direction.
          kStart=-2; kEnd=2
          jStart=-2; jEnd=2
          iStart=-2; iEnd=2

          if(iCell==2)  iStart=-1
          if(iCell==il) iEnd=1

          if(jCell==2)  jStart=-1
          if(jCell==jl) jEnd=1 

          if(kCell==kl) kEnd=1

           do k=kStart,kEnd !-2,2
             do j=jStart,jEnd !-2,2
               m = j -1
               do i=iStart,iEnd !-2,2
                 l = i -1

                 ! Determine the two diagonal vectors of the face.

                 v1(1) = xAdj(i,j,k,1) - xAdj(l,m,k,1)
                 v1(2) = xAdj(i,j,k,2) - xAdj(l,m,k,2)
                 v1(3) = xAdj(i,j,k,3) - xAdj(l,m,k,3)

                 v2(1) = xAdj(l,j,k,1) - xAdj(i,m,k,1)
                 v2(2) = xAdj(l,j,k,2) - xAdj(i,m,k,2)
                 v2(3) = xAdj(l,j,k,3) - xAdj(i,m,k,3)

                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact is
                 ! either -0.5 or 0.5.

                 skAdj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                 skAdj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                 skAdj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

               enddo
             enddo
           enddo
!
!        **************************************************************
!        *                                                            *
!        * If the considering cell is on the subfaces, then compute   *
!        * normAdj(-2:2,-2:2,-2:2,3). Otherwise return!               *
!        *                                                            *
!        * The unit normals on the boundary faces. These always point *
!        *  out of the domain, so a multiplication by -1 is needed for *
!        * the iMin, jMin and kMin boundaries.                        *
!        *                                                            *
!        **************************************************************
!

           ! Determine the range of the stencil for the given cell.
           
           iSBeg = iCell - 2; iSEnd = iCell + 2
           jSBeg = jCell - 2; jSEnd = jCell + 2
           kSBeg = kCell - 2; kSEnd = kCell + 2
           

           ! Loop over the number of physical boundary subfaces of the block.
           
           bocoLoop: do mm=1,nBocos
              
              ! Determine the range of halo cells which this boundary subface
              ! will change.
              
              select case (BCFaceID(mm))
              case (iMin)
                 iBBeg = 0;                iBEnd = 1
                 jBBeg = BCData(mm)%icBeg; jBEnd = BCData(mm)%icEnd
                 kBBeg = BCData(mm)%jcBeg; kBEnd = BCData(mm)%jcEnd
                 
                 !=============================================================
                 
              case (iMax)
                 iBBeg = ie;               iBEnd = ib
                 jBBeg = BCData(mm)%icBeg; jBEnd = BCData(mm)%icEnd
                 kBBeg = BCData(mm)%jcBeg; kBEnd = BCData(mm)%jcEnd
                 
                 !=============================================================
                 
              case (jMin)
                 iBBeg = BCData(mm)%icBeg; iBEnd = BCData(mm)%icEnd
                 jBBeg = 0;                jBEnd = 1
                 kBBeg = BCData(mm)%jcBeg; kBEnd = BCData(mm)%jcEnd
                 
                 !=============================================================
                 
              case (jMax)
                 iBBeg = BCData(mm)%icBeg; iBEnd = BCData(mm)%icEnd
                 jBBeg = je;               jBEnd = jb
                 kBBeg = BCData(mm)%jcBeg; kBEnd = BCData(mm)%jcEnd
                 
                 !=============================================================
                 
              case (kMin)
                 iBBeg = BCData(mm)%icBeg; iBEnd = BCData(mm)%icEnd
                 jBBeg = BCData(mm)%jcBeg; jBEnd = BCData(mm)%jcEnd
                 kBBeg = 0;                kBEnd = 1
                 
                 !=============================================================
                 
              case (kMax)
                 iBBeg = BCData(mm)%icBeg; iBEnd = BCData(mm)%icEnd
                 jBBeg = BCData(mm)%jcBeg; jBEnd = BCData(mm)%jcEnd
                 kBBeg = ke;               kBEnd = kb
                 
              end select
              
              ! Check for an overlap between the stencil range and the
              ! halo cells influenced by this boundary subface.
              
              iOverlap = .false.
              if(iSBeg <= iBEnd .and. iSEnd >= iBBeg) iOverlap = .true.
              
              jOverlap = .false.
              if(jSBeg <= jBEnd .and. jSEnd >= jBBeg) jOverlap = .true.
              
              kOverlap = .false.
              if(kSBeg <= kBEnd .and. kSEnd >= kBBeg) kOverlap = .true.
              
              checkOverlap: if(iOverlap .and. jOverlap .and. kOverlap) then
                 iRBeg = max(iSBeg, iBBeg); iREnd = min(iSEnd, iBEnd)
                 jRBeg = max(jSBeg, jBBeg); jREnd = min(jSEnd, jBEnd)
                 kRBeg = max(kSBeg, kBBeg); kREnd = min(kSEnd, kBEnd)
                 
                 iSt=-2; iEn=2
                 jSt=-2; jEn=2
                
                 select case (BCFaceID(mm))

                 case (iMin)
                    if(jCell==2)  iSt=-1
                    if(jCell==jl) iEn=1
                    if(kCell==2)  jSt=-1
                    if(kCell==kl)  jEn=1
                    secondHalo = .true.
                    if(iRBeg == iREnd) secondHalo = .false.
                    if(secondHalo)then
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = siAdj(-1,iSt:iEn,jSt:jEn,:)
                    else 
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = siAdj(-2,iSt:iEn,jSt:jEn,:)                       
                    end if
                    
                 case (iMax)
                    if(jCell==2)  iSt=-1
                    if(jCell==jl) iEn=1
                    if(kCell==2)  jSt=-1
                    if(kCell==kl)  jEn=1
                    secondHalo = .true.
                    if(iRBeg == iREnd) secondHalo = .false.
                    if(secondHalo) then
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = siAdj(0,iSt:iEn,jSt:jEn,:)
                    else
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = siAdj(1,iSt:iEn,jSt:jEn,:)
                    end if

                    
                 case (jMin)
                    if(iCell==2)  iSt=-1
                    if(iCell==il) iEn=1
                    if(kCell==2)  jSt=-1
                    if(kCell==kl) jEn=1
                    secondHalo = .true.
                    if(jRBeg == jREnd) secondHalo = .false.
                   if(secondHalo) then
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-1,jSt:jEn,:)
                    else
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-2,jSt:jEn,:)
                    end if                    

                 case (jMax)
                    if(iCell==2)  iSt=-1
                    if(iCell==il) iEn=1
                    if(kCell==2)  jSt=-1
                    if(kCell==kl) jEn=1                    
                    secondHalo = .true.
                    if(jRBeg == jREnd) secondHalo = .false.
                   if(secondHalo) then
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,0,jSt:jEn,:)
                    else
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,1,jSt:jEn,:)
                    end if                    
                    
                 case (kMin)
                    if(iCell==2)  iSt=-1
                    if(iCell==il) iEn=1
                    if(jCell==2)  jSt=-1
                    if(jCell==jl) jEn=1
                    secondHalo = .true.
                    if(kRBeg == kREnd) secondHalo = .false.
                    if(secondHalo) then
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-1,:)
                    else
                       mult = -one; ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-2,:)
                    end if

                    
                 case (kMax)
                    if(iCell==2)  iSt=-1
                    if(iCell==il) iEn=1
                    if(jCell==2)  jSt=-1
                    if(jCell==jl) jEn=1                    
                    secondHalo = .true.
                    if(kRBeg == kREnd) secondHalo = .false.
                    if(secondHalo) then
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,0,:)
                    else
                       mult = one;  ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,1,:)
                    end if                    
                 end select
                 
                 do kk=jSt,jEn
                    do jj=iSt,iEn

                       ! Compute the inverse of the length of the normal vector
                       ! and possibly correct for inward pointing.
                       
                       xp = ss(jj,kk,1);  yp = ss(jj,kk,2);  zp = ss(jj,kk,3)
                       fact = sqrt(xp*xp + yp*yp + zp*zp)
                       if(fact > zero) fact = mult/fact
                       
                       ! Compute the unit normal.
                       
                       normAdj(mm,jj,kk,1) = fact*xp
                       normAdj(mm,jj,kk,2) = fact*yp
                       normAdj(mm,jj,kk,3) = fact*zp
                       
                    enddo
                 enddo


              end if checkOverlap
           end do bocoLoop


!
!          **************************************************************
!          *                                                            *
!          * Check in debug mode the sum of the normals of the cells.   *
!          * If everything is correct this should sum up to zero.       *
!          *                                                            *
!          **************************************************************
!
           debugging: if( debug ) then

              ! Loop over the cells including the 1st level halo's.
              
              i=0;   j=0;   k=0
              l=i-1; m=j-1; n=k-1

              ! Store the sum of the outward pointing surrounding
              ! normals in v1. Due to the outward convention the
              ! normals with the lowest index get a negative sign;
              ! normals point in the direction of the higher index.
              
              v1(1) = siAdj(i,j,k,1) + sjAdj(i,j,k,1) + skAdj(i,j,k,1) &
                    - siAdj(l,j,k,1) - sjAdj(i,m,k,1) - skAdj(i,j,n,1)
              v1(2) = siAdj(i,j,k,2) + sjAdj(i,j,k,2) + skAdj(i,j,k,2) &
                    - siAdj(l,j,k,2) - sjAdj(i,m,k,2) - skAdj(i,j,n,2)
              v1(3) = siAdj(i,j,k,3) + sjAdj(i,j,k,3) + skAdj(i,j,k,3) &
                    - siAdj(l,j,k,3) - sjAdj(i,m,k,3) - skAdj(i,j,n,3)
              
              ! Store the inverse of the sum of the areas of the
              ! six faces in fact.

              fact = one/(sqrt(siAdj(i,j,k,1)*siAdj(i,j,k,1)  &
                   +           siAdj(i,j,k,2)*siAdj(i,j,k,2)  &
                   +           siAdj(i,j,k,3)*siAdj(i,j,k,3)) &
                   +      sqrt(siAdj(l,j,k,1)*siAdj(l,j,k,1)  &
                   +           siAdj(l,j,k,2)*siAdj(l,j,k,2)  &
                   +           siAdj(l,j,k,3)*siAdj(l,j,k,3)) &
                   +      sqrt(sjAdj(i,j,k,1)*sjAdj(i,j,k,1)  &
                   +           sjAdj(i,j,k,2)*sjAdj(i,j,k,2)  &
                   +           sjAdj(i,j,k,3)*sjAdj(i,j,k,3)) &
                   +      sqrt(sjAdj(i,m,k,1)*sjAdj(i,m,k,1)  &
                   +           sjAdj(i,m,k,2)*sjAdj(i,m,k,2)  &
                   +           sjAdj(i,m,k,3)*sjAdj(i,m,k,3)) &
                   +      sqrt(skAdj(i,j,k,1)*skAdj(i,j,k,1)  &
                   +           skAdj(i,j,k,2)*skAdj(i,j,k,2)  &
                   +           skAdj(i,j,k,3)*skAdj(i,j,k,3)) &
                   +      sqrt(skAdj(i,j,n,1)*skAdj(i,j,n,1)  &
                   +           skAdj(i,j,n,2)*skAdj(i,j,n,2)  &
                   +           skAdj(i,j,n,3)*skAdj(i,j,n,3)))
              
              ! Multiply v1 by fact to obtain a nonDimensional
              ! quantity and take tha absolute value of it.
              
              v1(1) = abs(v1(1)*fact)
              v1(2) = abs(v1(2)*fact)
              v1(3) = abs(v1(3)*fact)
              
              ! Check if the control volume is closed.
              
              if(v1(1) > thresholdReal .or. &
                   v1(2) > thresholdReal .or. &
                   v1(3) > thresholdReal)     &
                   write(*,*)'NORMALS DO NOT SUM UP TO 0',icell,jcell,kcell
!                   call terminate("metric", &
!                   "Normals do not sum up to 0")

     
  endif debugging
  
end subroutine metricAdj


!      ==================================================================
!      ================================================================

  subroutine volpym2(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,xp,yp,zp,vp)
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

     end subroutine volpym2



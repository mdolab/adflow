
subroutine computeNormNKPC(siAdj,sjAdj,skAdj,normAdj, &
     iCell,jCell,kCell,nn,level,sps,sps2)
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
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,&
       nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,&
       nTimeIntervalsSpectral)  :: normAdj
  integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,level,sps,sps2

  !!      Local variables.
  !
  integer(kind=intType) :: ierr

  integer(kind=intType) :: i, j, k, n, m, l,iSt,iEn,jSt,jEn,jj,kk,mm
  integer(kind=intType) :: iSBeg,iSEnd,jSBeg,jSEnd,kSBeg,kSEnd
  integer(kind=intType) :: iBBeg,iBEnd,jBBeg,jBEnd,kBBeg,kBEnd
  integer(kind=intType) :: iRBeg,iREnd,jRBeg,jREnd,kRBeg,kREnd
  integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd
  logical ::  iOverLap,jOverlap,kOverlap
  real(kind=realType) :: mult
  logical :: secondHalo
  real(kind=realType), dimension(-2:2,-2:2,3) :: ss
  real(kind=realType) :: xp,yp,zp,fact
  secondHalo = .true.
  normAdj(:,:,:,:,sps2) = zero
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
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = siAdj(-1,iSt:iEn,jSt:jEn,:,sps2)
           else 
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = siAdj(-2,iSt:iEn,jSt:jEn,:,sps2)                       
           end if

        case (iMax)
           if(jCell==2)  iSt=-1
           if(jCell==jl) iEn=1
           if(kCell==2)  jSt=-1
           if(kCell==kl)  jEn=1
           secondHalo = .true.
           if(iRBeg == iREnd) secondHalo = .false.
           if(secondHalo) then
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = siAdj(0,iSt:iEn,jSt:jEn,:,sps2)
           else
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = siAdj(1,iSt:iEn,jSt:jEn,:,sps2)
           end if


        case (jMin)
           if(iCell==2)  iSt=-1
           if(iCell==il) iEn=1
           if(kCell==2)  jSt=-1
           if(kCell==kl) jEn=1
           secondHalo = .true.
           if(jRBeg == jREnd) secondHalo = .false.
           if(secondHalo) then
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-1,jSt:jEn,:,sps2)
           else
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,-2,jSt:jEn,:,sps2)
           end if

        case (jMax)
           if(iCell==2)  iSt=-1
           if(iCell==il) iEn=1
           if(kCell==2)  jSt=-1
           if(kCell==kl) jEn=1                    
           secondHalo = .true.
           if(jRBeg == jREnd) secondHalo = .false.
           if(secondHalo) then
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,0,jSt:jEn,:,sps2)
           else
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = sjAdj(iSt:iEn,1,jSt:jEn,:,sps2)
           end if

        case (kMin)
           if(iCell==2)  iSt=-1
           if(iCell==il) iEn=1
           if(jCell==2)  jSt=-1
           if(jCell==jl) jEn=1
           secondHalo = .true.
           if(kRBeg == kREnd) secondHalo = .false.
           if(secondHalo) then
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-1,:,sps2)
           else
              mult = -one; ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,-2,:,sps2)
           end if


        case (kMax)
           if(iCell==2)  iSt=-1
           if(iCell==il) iEn=1
           if(jCell==2)  jSt=-1
           if(jCell==jl) jEn=1                    
           secondHalo = .true.
           if(kRBeg == kREnd) secondHalo = .false.
           if(secondHalo) then
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,0,:,sps2)
           else
              mult = one;  ss(iSt:iEn,jSt:jEn,:) = skAdj(iSt:iEn,jSt:jEn,1,:,sps2)
           end if
        end select

        do kk=jSt,jEn
           do jj=iSt,iEn

              ! Compute the inverse of the length of the normal vector
              ! and possibly correct for inward pointing.
                       
              xp = ss(jj,kk,1);  yp = ss(jj,kk,2);  zp = ss(jj,kk,3)

              !alternate form to allow inclusion of degenrate halos???
              if( xp**2>zero .or. yp**2>zero .or. zp**2>zero)then
                 !if (fact > zero)then
                 !compute length
                 fact = sqrt(xp*xp + yp*yp + zp*zp)
                 !set factor to 1/length
                 fact = mult/fact
                 !compute unit normal...
                 normAdj(mm,jj,kk,1,sps2) = fact*xp
                 normAdj(mm,jj,kk,2,sps2) = fact*yp
                 normAdj(mm,jj,kk,3,sps2) = fact*zp
              else
                 !Length is zero
                 normAdj(mm,jj,kk,:,sps2) = zero
              endif
           enddo
        enddo
     end if checkOverlap
 end do bocoLoop


end subroutine computeNormNKPC


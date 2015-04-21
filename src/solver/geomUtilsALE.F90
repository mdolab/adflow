!
!      ******************************************************************
!      *                                                                *
!      * This set of utilities operates on geometry-related data        *
!      * members for ALE scheme. They are seperated out so that         *
!      * cumbersome loops over spectral and domains can be removed from *
!      * other subroutines.                                             *
!      * Subroutines labeled _block operates only on one block, so      *
!      * setPointer has to be used in advance. In those cases, block    *
!      * and boundary data are treated seperately.                      *
!      *                                                                *
!      * The data members considered here include:                      *
!      * sFace[I,J,K], s[I,J,K], norm, rFace, uSlip,                    *
!      * and their ALE counterparts                                     *
!      *                                                                *
!      * Added by HDN                                                   *
!      *                                                                *
!      ******************************************************************
!
! ===========================================================
subroutine setLevelALE(setType)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * setLevelALE sets specified ALE level(s) with current data.     *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Input variables.
  !
  integer(kind=intType), intent(in) :: setType
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,l, aleBeg,aleEnd, nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  select case (setType)
  case (-1_intType)
     aleBeg = 0_intType
     aleEnd = nALEsteps
  case (-2_intType)
     aleBeg = 1_intType
     aleEnd = nALEMeshes
  case default
     aleBeg = setType
     aleEnd = setType
  end select

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom

        ! Set the pointers for this block on the ground level.

        call setPointers(nn, groundLevel, kk)

        blkALE : do l = aleBeg,aleENd
           fillI2 : do k = 1,ke
              do j = 1,je
                 do i = 0,ie
                    sFaceIALE(l,i,j,k) = sFaceI(i,j,k)
                    sIALE(l,i,j,k,1)   = sI(i,j,k,1)
                    sIALE(l,i,j,k,2)   = sI(i,j,k,2)
                    sIALE(l,i,j,k,3)   = sI(i,j,k,3)
                 enddo
              enddo
           enddo fillI2

           fillJ2 : do k = 1,ke
              do j = 0,je
                 do i = 1,ie
                    sFaceJALE(l,i,j,k) = sFaceJ(i,j,k)
                    sJALE(l,i,j,k,1)   = sJ(i,j,k,1)
                    sJALE(l,i,j,k,2)   = sJ(i,j,k,2)
                    sJALE(l,i,j,k,3)   = sJ(i,j,k,3)
                 enddo
              enddo
           enddo fillJ2

           fillK2 : do k = 0,ke
              do j = 1,je
                 do i = 1,ie
                    sFaceKALE(l,i,j,k) = sFaceK(i,j,k)
                    sKALE(l,i,j,k,1)   = sK(i,j,k,1)
                    sKALE(l,i,j,k,2)   = sK(i,j,k,2)
                    sKALE(l,i,j,k,3)   = sK(i,j,k,3)
                 enddo
              enddo
           enddo fillK2
        enddo blkALE

        normLoop: do mm=1,nBocos
           do l = aleBeg,aleENd
              do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%normALE(l,i,j,1) = BCData(mm)%norm(i,j,1)
                    BCData(mm)%normALE(l,i,j,2) = BCData(mm)%norm(i,j,2)
                    BCData(mm)%normALE(l,i,j,3) = BCData(mm)%norm(i,j,3)
                 enddo
              enddo
           enddo
        enddo normLoop

        rFaceLoop: do mm=1,nBocos
           testAssoc: if( associated(BCData(mm)%rFace) ) then
              do l = aleBeg,aleENd
                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                       BCData(mm)%rFaceALE(l,i,j) = BCData(mm)%rFace(i,j)
                    enddo
                 enddo
              enddo
           endif testAssoc
        enddo rFaceLoop

        uSlipLoop: do mm=1,nViscBocos
           do l = aleBeg,aleENd
              do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%uSlipALE(l,i,j,1) = BCData(mm)%uSlip(i,j,1)
                    BCData(mm)%uSlipALE(l,i,j,2) = BCData(mm)%uSlip(i,j,2)
                    BCData(mm)%uSlipALE(l,i,j,3) = BCData(mm)%uSlip(i,j,3)
                 enddo
              enddo
           enddo
        enddo uSlipLoop

     end do domains
  end do spectralLoop

end subroutine setLevelALE

! ===========================================================
subroutine shiftLevelALE
  !
  !      ******************************************************************
  !      *                                                                *
  !      * shiftLevelALE move current ALE levels to older levels and      *
  !      * update them with current data.                                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,l,lo,nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom

        ! Set the pointers for this block on the ground level.

        call setPointers(nn, groundLevel, kk)

        blkALE : do l = nALEsteps, nALEMeshes+1, -1
           lo = l - nALEMeshes

           fillI2 : do k = 1,ke
              do j = 1,je
                 do i = 0,ie
                    sFaceIALE(l,i,j,k) = sFaceIALE(lo,i,j,k)
                    sIALE(l,i,j,k,1)   = sIALE(lo,i,j,k,1)
                    sIALE(l,i,j,k,2)   = sIALE(lo,i,j,k,2)
                    sIALE(l,i,j,k,3)   = sIALE(lo,i,j,k,3)
                 enddo
              enddo
           enddo fillI2

           fillJ2 : do k = 1,ke
              do j = 0,je
                 do i = 1,ie
                    sFaceJALE(l,i,j,k) = sFaceJALE(lo,i,j,k)
                    sJALE(l,i,j,k,1)   = sJALE(lo,i,j,k,1)
                    sJALE(l,i,j,k,2)   = sJALE(lo,i,j,k,2)
                    sJALE(l,i,j,k,3)   = sJALE(lo,i,j,k,3)
                 enddo
              enddo
           enddo fillJ2

           fillK2 : do k = 0,ke
              do j = 1,je
                 do i = 1,ie
                    sFaceKALE(l,i,j,k) = sFaceKALE(lo,i,j,k)
                    sKALE(l,i,j,k,1)   = sKALE(lo,i,j,k,1)
                    sKALE(l,i,j,k,2)   = sKALE(lo,i,j,k,2)
                    sKALE(l,i,j,k,3)   = sKALE(lo,i,j,k,3)
                 enddo
              enddo
           enddo fillK2
        enddo blkALE

        normLoop: do mm=1,nBocos
           do l = nALEsteps, nALEMeshes+1, -1
              lo = l - nALEMeshes
              do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%normALE(l,i,j,1) = BCData(mm)%normALE(lo,i,j,1)
                    BCData(mm)%normALE(l,i,j,2) = BCData(mm)%normALE(lo,i,j,2)
                    BCData(mm)%normALE(l,i,j,3) = BCData(mm)%normALE(lo,i,j,3)
                 enddo
              enddo
           enddo
        enddo normLoop

        rFaceLoop: do mm=1,nBocos
           testAssoc: if( associated(BCData(mm)%rFace) ) then
              do l = nALEsteps, nALEMeshes+1, -1
                 lo = l - nALEMeshes
                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                       BCData(mm)%rFaceALE(l,i,j) = BCData(mm)%rFaceALE(lo,i,j)
                    enddo
                 enddo
              enddo
           endif testAssoc
        enddo rFaceLoop

        uSlipLoop: do mm=1,nViscBocos
           do l = nALEsteps, nALEMeshes+1, -1
              lo = l - nALEMeshes
              do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%uSlipALE(l,i,j,1) = BCData(mm)%uSlipALE(lo,i,j,1)
                    BCData(mm)%uSlipALE(l,i,j,2) = BCData(mm)%uSlipALE(lo,i,j,2)
                    BCData(mm)%uSlipALE(l,i,j,3) = BCData(mm)%uSlipALE(lo,i,j,3)
                 enddo
              enddo
           enddo
        enddo uSlipLoop

     end do domains
  end do spectralLoop

  ! Set latest levels with current data
  call setLevelALE(-2_intType)

end subroutine shiftLevelALE

! ===========================================================
subroutine interpLevelALE_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * interpLevelALE_block interpolates geometric data over the      *
  !      * latest time step.                                              *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use iteration
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,l,nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  ! --------------------------------
  ! First store then clear current data
  ! --------------------------------
  clearI : do k = 1,ke
     do j = 1,je
        do i = 0,ie
           sFaceIALE(0,i,j,k) = sFaceI(i,j,k)
           sIALE(0,i,j,k,1)   = sI(i,j,k,1)
           sIALE(0,i,j,k,2)   = sI(i,j,k,2)
           sIALE(0,i,j,k,3)   = sI(i,j,k,3)
           sFaceI(i,j,k)      = zero
           sI(i,j,k,1)        = zero
           sI(i,j,k,2)        = zero
           sI(i,j,k,3)        = zero
        enddo
     enddo
  enddo clearI

  clearJ : do k = 1,ke
     do j = 0,je
        do i = 1,ie
           sFaceJALE(0,i,j,k) = sFaceJ(i,j,k)
           sJALE(0,i,j,k,1)   = sJ(i,j,k,1)
           sJALE(0,i,j,k,2)   = sJ(i,j,k,2)
           sJALE(0,i,j,k,3)   = sJ(i,j,k,3)
           sFaceJ(i,j,k)      = zero
           sJ(i,j,k,1)        = zero
           sJ(i,j,k,2)        = zero
           sJ(i,j,k,3)        = zero
        enddo
     enddo
  enddo clearJ

  clearK : do k = 0,ke
     do j = 1,je
        do i = 1,ie
           sFaceKALE(0,i,j,k) = sFaceK(i,j,k)
           sKALE(0,i,j,k,1)   = sK(i,j,k,1)
           sKALE(0,i,j,k,2)   = sK(i,j,k,2)
           sKALE(0,i,j,k,3)   = sK(i,j,k,3)
           sFaceK(i,j,k)      = zero
           sK(i,j,k,1)        = zero
           sK(i,j,k,2)        = zero
           sK(i,j,k,3)        = zero
        enddo
     enddo
  enddo clearK

  ALEloop : do l = 1,nALEsteps
     ! --------------------------------
     ! Then average surface normal and normal velocity from array of old variables
     ! --------------------------------
     updateI : do k = 1,ke
        do j = 1,je
           do i = 0,ie
              sFaceI(i,j,k) = sFaceI(i,j,k) + coefTimeALE(l) * sFaceIALE(l,i,j,k)
              sI(i,j,k,1)   = sI(i,j,k,1)   + coefTimeALE(l) * sIALE(l,i,j,k,1)
              sI(i,j,k,2)   = sI(i,j,k,2)   + coefTimeALE(l) * sIALE(l,i,j,k,2)
              sI(i,j,k,3)   = sI(i,j,k,3)   + coefTimeALE(l) * sIALE(l,i,j,k,3)
           enddo
        enddo
     enddo updateI

     updateJ : do k = 1,ke
        do j = 0,je
           do i = 1,ie
              sFaceJ(i,j,k) = sFaceJ(i,j,k) + coefTimeALE(l) * sFaceJALE(l,i,j,k)
              sJ(i,j,k,1)   = sJ(i,j,k,1)   + coefTimeALE(l) * sJALE(l,i,j,k,1)
              sJ(i,j,k,2)   = sJ(i,j,k,2)   + coefTimeALE(l) * sJALE(l,i,j,k,2)
              sJ(i,j,k,3)   = sJ(i,j,k,3)   + coefTimeALE(l) * sJALE(l,i,j,k,3)
           enddo
        enddo
     enddo updateJ

     updateK : do k = 0,ke
        do j = 1,je
           do i = 1,ie
              sFaceK(i,j,k) = sFaceK(i,j,k) + coefTimeALE(l) * sFaceKALE(l,i,j,k)
              sK(i,j,k,1)   = sK(i,j,k,1)   + coefTimeALE(l) * sKALE(l,i,j,k,1)
              sK(i,j,k,2)   = sK(i,j,k,2)   + coefTimeALE(l) * sKALE(l,i,j,k,2)
              sK(i,j,k,3)   = sK(i,j,k,3)   + coefTimeALE(l) * sKALE(l,i,j,k,3)
           enddo
        enddo
     enddo updateK
  enddo ALEloop


end subroutine interpLevelALE_block

! ===========================================================
subroutine recoverLevelALE_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * recoverLevelALE_block recovers current geometric data from     *
  !      * temporary interpolation                                        *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  recoverI : do k = 1,ke
     do j = 1,je
        do i = 0,ie
           sFaceI(i,j,k) = sFaceIALE(0,i,j,k)
           sI(i,j,k,1)   = sIALE(0,i,j,k,1)
           sI(i,j,k,2)   = sIALE(0,i,j,k,2)
           sI(i,j,k,3)   = sIALE(0,i,j,k,3)
        enddo
     enddo
  enddo recoverI

  recoverJ : do k = 1,ke
     do j = 0,je
        do i = 1,ie
           sFaceJ(i,j,k) = sFaceJALE(0,i,j,k)
           sJ(i,j,k,1)   = sJALE(0,i,j,k,1)
           sJ(i,j,k,2)   = sJALE(0,i,j,k,2)
           sJ(i,j,k,3)   = sJALE(0,i,j,k,3)
        enddo
     enddo
  enddo recoverJ

  recoverK : do k = 0,ke
     do j = 1,je
        do i = 1,ie
           sFaceK(i,j,k) = sFaceKALE(0,i,j,k)
           sK(i,j,k,1)   = sKALE(0,i,j,k,1)
           sK(i,j,k,2)   = sKALE(0,i,j,k,2)
           sK(i,j,k,3)   = sKALE(0,i,j,k,3)
        enddo
     enddo
  enddo recoverK


end subroutine recoverLevelALE_block

! ===========================================================
subroutine interpLevelALEBC_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * interpLevelALEBC_block interpolates geometric data on boundary *
  !      * over the latest time step.                                     *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use iteration
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,l,nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  ! --------------------------------
  ! First store then clear current data
  ! --------------------------------
  clearNM: do mm=1,nBocos
     do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
        do i=BCData(mm)%icBeg, BCData(mm)%icEnd
           BCData(mm)%normALE(0,i,j,1) = BCData(mm)%norm(i,j,1)
           BCData(mm)%normALE(0,i,j,2) = BCData(mm)%norm(i,j,2)
           BCData(mm)%normALE(0,i,j,3) = BCData(mm)%norm(i,j,3)
           BCData(mm)%norm(i,j,1)      = zero
           BCData(mm)%norm(i,j,2)      = zero
           BCData(mm)%norm(i,j,3)      = zero
        enddo
     enddo
  enddo clearNM

  clearRF: do mm=1,nBocos
     testAssoc1: if( associated(BCData(mm)%rFace) ) then
        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd
              BCData(mm)%rFaceALE(0,i,j) = BCData(mm)%rFace(i,j)
              BCData(mm)%rFace(i,j)      = zero
           enddo
        enddo
     endif testAssoc1
  enddo clearRF

  clearUS: do mm=1,nViscBocos
     do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
        do i=BCData(mm)%icBeg, BCData(mm)%icEnd
           BCData(mm)%uSlipALE(0,i,j,1) = BCData(mm)%uSlip(i,j,1)
           BCData(mm)%uSlipALE(0,i,j,2) = BCData(mm)%uSlip(i,j,2)
           BCData(mm)%uSlipALE(0,i,j,3) = BCData(mm)%uSlip(i,j,3)
           BCData(mm)%uSlip(i,j,1)      = zero
           BCData(mm)%uSlip(i,j,2)      = zero
           BCData(mm)%uSlip(i,j,3)      = zero
        enddo
     enddo
  enddo clearUS

  ALEloop : do l = 1,nALEsteps
     ! --------------------------------
     ! Then average surface normal and normal velocity from array of old variables
     ! --------------------------------
     updateNM: do mm=1,nBocos
        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd
              BCData(mm)%norm(i,j,1) = BCData(mm)%norm(i,j,1) &
                   + coefTimeALE(l) * BCData(mm)%normALE(l,i,j,1)
              BCData(mm)%norm(i,j,2) = BCData(mm)%norm(i,j,2) &
                   + coefTimeALE(l) * BCData(mm)%normALE(l,i,j,2)
              BCData(mm)%norm(i,j,3) = BCData(mm)%norm(i,j,3) &
                   + coefTimeALE(l) * BCData(mm)%normALE(l,i,j,3)
           enddo
        enddo
     enddo updateNM

     updateRF: do mm=1,nBocos
        testAssoc2: if( associated(BCData(mm)%rFace) ) then
           do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
              do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                 BCData(mm)%rFace(i,j) = BCData(mm)%rFace(i,j) &
                      + coefTimeALE(l) * BCData(mm)%rFaceALE(0,i,j)
              enddo
           enddo
        endif testAssoc2
     enddo updateRF

     updateUS: do mm=1,nViscBocos
        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd
              BCData(mm)%uSlip(i,j,1) = BCData(mm)%uSlip(i,j,1) &
                   + coefTimeALE(l) * BCData(mm)%uSlipALE(l,i,j,1)
              BCData(mm)%uSlip(i,j,2) = BCData(mm)%uSlip(i,j,2) &
                   + coefTimeALE(l) * BCData(mm)%uSlipALE(l,i,j,2)
              BCData(mm)%uSlip(i,j,3) = BCData(mm)%uSlip(i,j,3) &
                   + coefTimeALE(l) * BCData(mm)%uSlipALE(l,i,j,3)
           enddo
        enddo
     enddo updateUS
  enddo ALEloop

end subroutine interpLevelALEBC_block

! ===========================================================
subroutine recoverLevelALEBC_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * recoverLevelALEBC_block recovers current geometric data on     *
  !      * boundary from temporary interpolation                          *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,mm,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  recoverNM: do mm=1,nBocos
     do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
        do i=BCData(mm)%icBeg, BCData(mm)%icEnd
           BCData(mm)%norm(i,j,1) = BCData(mm)%normALE(0,i,j,1)
           BCData(mm)%norm(i,j,2) = BCData(mm)%normALE(0,i,j,2)
           BCData(mm)%norm(i,j,3) = BCData(mm)%normALE(0,i,j,3)
        enddo
     enddo
  enddo recoverNM

  recoverRF: do mm=1,nBocos
     testAssoc: if( associated(BCData(mm)%rFace) ) then
        do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
           do i=BCData(mm)%icBeg, BCData(mm)%icEnd
              BCData(mm)%rFace(i,j) = BCData(mm)%rFaceALE(0,i,j)
           enddo
        enddo
     endif testAssoc
  enddo recoverRF

  recoverUS: do mm=1,nViscBocos
     do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
        do i=BCData(mm)%icBeg, BCData(mm)%icEnd
           BCData(mm)%uSlip(i,j,1) = BCData(mm)%uSlipALE(0,i,j,1)
           BCData(mm)%uSlip(i,j,2) = BCData(mm)%uSlipALE(0,i,j,2)
           BCData(mm)%uSlip(i,j,3) = BCData(mm)%uSlipALE(0,i,j,3)
        enddo
     enddo
  enddo recoverUS

end subroutine recoverLevelALEBC_block

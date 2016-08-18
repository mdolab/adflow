!
!      ******************************************************************
!      *                                                                *
!      * This set of utilities operates on coordinate-related data      *
!      * members for ALE scheme. They are seperated out so that         *
!      * cumbersome loops over spectral and domains can be removed from *
!      * other subroutines. The subroutines include:                    *
!      * fillCoor:    Initialize xOld and volOld with current x and vol *
!      * storeCoor:   Store x to temporary variable xtmp                *
!      * interpCoor:  Interpolate x over xALE and xOld                  *
!      * recoverCoor: Recover x from xALE                               *
!      *                                                                *
!      * Added by HDN                                                   *
!      *                                                                *
!      ******************************************************************
!
! ===========================================================
subroutine fillCoor
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
             
        ! Set the pointers for this block on the ground level.
        
        call setPointers(nn, groundLevel,kk)
        
        do k=0,ke
           do j=0,je
              do i=0,ie
                 xOld(:,i,j,k,1) = x(i,j,k,1)
                 xOld(:,i,j,k,2) = x(i,j,k,2)
                 xOld(:,i,j,k,3) = x(i,j,k,3)
              enddo
           enddo
        enddo

        do k=2,kl
           do j=2,jl
              do i=2,il
                 volOld(:,i,j,k) = vol(i,j,k)
              enddo
           enddo
        enddo

     end do domains
  end do spectralLoop
  
end subroutine fillCoor

! ===========================================================
subroutine storeCoor
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then 
     return
  end if

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        
        ! Set the pointers for this block on the ground level.
        
        call setPointers(nn, groundLevel,kk)
        
        storex : do k = 0,ke
           do j = 0,je
              do i = 0,ie
                 xALE(i,j,k,1) = x(i,j,k,1)
                 xALE(i,j,k,2) = x(i,j,k,2)
                 xALE(i,j,k,3) = x(i,j,k,3)
              enddo
           enddo
        enddo storex

     end do domains
  end do spectralLoop
  
end subroutine storeCoor

! ===========================================================
subroutine interpCoor(lale)
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Input variables.
  !
  integer(kind=intType), intent(in) :: lale
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,kk

  if (.not. useALE .or. equationMode .ne. unsteady)  then
     return
  end if

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom

        ! Set the pointers for this block on the ground level.
        ! This eq. 11a, found paper by C.Farhat http://dx.doi.org/10.1016/S0021-9991(03)00311-5

        call setPointers(nn, groundLevel,kk)

        interpmesh : do k = 0,ke
           do j = 0,je
              do i = 0,ie
                 x(i,j,k,1) = coefMeshALE(lale,1)*xALE(i,j,k,1) &
                            + coefMeshALE(lale,2)*xOld(1,i,j,k,1)
                 x(i,j,k,2) = coefMeshALE(lale,1)*xALE(i,j,k,2) &
                            + coefMeshALE(lale,2)*xOld(1,i,j,k,2)
                 x(i,j,k,3) = coefMeshALE(lale,1)*xALE(i,j,k,3) &
                            + coefMeshALE(lale,2)*xOld(1,i,j,k,3)
              enddo
           enddo
        enddo interpmesh

     end do domains
  end do spectralLoop

end subroutine interpCoor

! ===========================================================
subroutine recoverCoor
  use blockPointers
  use iteration
  use inputTimeSpectral
  use inputUnsteady
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: i,j,k,nn,kk


  if (.not. useALE .or. equationMode .ne. unsteady)  then
     return
  end if

  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom

        ! Set the pointers for this block on the ground level.

        call setPointers(nn, groundLevel,kk)

        recoverx : do k = 0,ke
           do j = 0,je
              do i = 0,ie
                 x(i,j,k,1) = xALE(i,j,k,1)
                 x(i,j,k,2) = xALE(i,j,k,2)
                 x(i,j,k,3) = xALE(i,j,k,3)
              enddo
           enddo
        enddo recoverx

     end do domains
  end do spectralLoop

end subroutine recoverCoor


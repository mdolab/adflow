!
!      ******************************************************************
!      *                                                                *
!      * wOverset controls the communication between overset halos      *
!      * for the cell-centered variables by interpolating the solution  *
!      * from other blocks consistent with the chimera approach. A tri- *
!      * linear interpolation is used. It is possible to send only a    *
!      * range of variables.                                            *
!      * set, e.g. only the flow variables or only the turbulent        *
!      * variables. This is controlled by the arguments start, end,     *
!      * commPressure and commViscous. The exchange takes place for     *
!      * the given grid level.                                          *
!      *                                                                *
!      ******************************************************************

subroutine wOverset(level, start, end, commPressure, commVarGamma, &
     commLamVis, commEddyVis)
  
  use blockPointers
  use communication
  use inputTimeSpectral
  use overset
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: level, start, end
  logical, intent(in) :: commPressure, commVarGamma
  logical, intent(in) :: commLamVis, commEddyVis

  ! Working variables
  integer :: iVar, nn, l, id, ii, jj, kk, nComm, ierr, sps
  real(kind=realType), dimension(:), pointer :: donorPtr, fringePtr
  real(kind=realType) :: f1, f2, f3, f4, f5, f6, f7, f8
  real(kind=realType) :: di0, di1, dj0, dj1, dk0, dk1

  ! Assume sps=1 for now
  sps = 1
  
  ! Will do each of the required communications individually for now.
  
  ! Determine the number of variables to communicate  
  iVar = 1
  do l=start, end
     iVar = iVar + 1
  end do

  if (commPressure) then 
     iVar = iVar + 1
  end if

  if (commPressure) then 
     iVar = iVar + 1
  end if

  if (commLamVis) then 
     iVar = iVar + 1
  end if

  if (commEddyVis) then 
     iVar = iVar + 1
  end if

  ! Total number of required communiations
  nComm = iVar

  ! Allocate the 2D "variables" derived tpe
  allocate(variables(nComm, nDom))

  ! Now set all the required pointers. We will need these when filling
  ! up the donors as well as when we set the fringes.
  do nn=1, nDom
     iVar = 1
     do l=start, end
        variables(iVar, nn)%arr => w(:, :, :, l)
        iVar = iVar + 1
     end do

     if (commPressure) then 
        variables(iVar, nn)%arr => P(:, :, :)
        iVar = iVar + 1
     end if

     if (commPressure) then 
        variables(iVar, nn)%arr => gamma(:, :, :)
        iVar = iVar + 1
     end if

     if (commLamVis) then 
        variables(iVar, nn)%arr => rlv(:, :, :)
        iVar = iVar + 1
     end if
     
     if (commEddyVis) then 
        variables(iVar, nn)%arr => rev(:, :, :)
        iVar = iVar + 1
     end if
  end do

  ! Master loop over the required number of communications
  masterLoop: do iVar = 1, nComm

     call VecGetArrayF90(oversetDonors, donorPtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Counter for donorPtr array
     ii = 0

     ! Loop over the number of blocks on this processor
     do nn=1, nDom
     
        ! No setPointers here for speed

        ! Loop of the number of donors on this block
        do id=1, flowDoms(nn, level, sps)%nDonor
           
           ii = flowDoms(nn, level, sps)%idonor(1, id)
           jj = flowDoms(nn, level, sps)%idonor(2, id)
           kk = flowDoms(nn, level, sps)%idonor(3, id)
          
           ! Reconstruct the weights from the fractions. We are using
           ! linear interpolation here
           di0 = one - flowDoms(nn, level, sps)%frac(1, id)
           dj0 = one - flowDoms(nn, level, sps)%frac(2, id)
           dk0 = one - flowDoms(nn, level, sps)%frac(3, id)
           
           di1 = flowDoms(nn, level, sps)%frac(1, id)
           dj1 = flowDoms(nn, level, sps)%frac(2, id)
           dk1 = flowDoms(nn, level, sps)%frac(3, id)
           
           ! The corresponding linear weights
           f1   = di0*dj0*dk0
           f2   = di1*dj0*dk0
           f3   = di0*dj1*dk0
           f4   = di1*dj1*dk0
           f5   = di0*dj0*dk1
           f6   = di1*dj0*dk1
           f7   = di0*dj1*dk1
           f8   = di1*dj1*dk1
     
           ! Set this donor point into the arrary from the petsc vector
           ii = ii + 1
           donorPtr(ii) = &
                f1*variables(iVar, nn)%arr(ii  , jj,   kk  ) + &
                f2*variables(iVar, nn)%arr(ii+1, jj,   kk  ) + &
                f3*variables(iVar, nn)%arr(ii  , jj+1, kk  ) + &
                f4*variables(iVar, nn)%arr(ii+1, jj+1, kk  ) + &
                f5*variables(iVar, nn)%arr(ii  , jj,   kk+1) + &
                f6*variables(iVar, nn)%arr(ii+1, jj,   kk+1) + &
                f7*variables(iVar, nn)%arr(ii  , jj+1, kk+1) + &
                f8*variables(iVar, nn)%arr(ii+1, jj+1, kk+1)
        end do
     end do

     ! Return the pointer to petsc before the scatter
     call VecRestoreArrayF90(oversetDonors, donorPtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     
     ! Now perform the actual vec scatter
     call VecScatterBegin(oversetScatter, oversetDonors, oversetFringes, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Now we pull out a pointer from oversetFringe and set the
     ! required data back into the required array
     call VecGetArrayF90(oversetFringes, fringePtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Counter for fringePtr array
     ii = 0

     ! Loop over the number of blocks on this processor
     do nn=1, nDom
     
        ! No setPointers here for speed

        ! Loop of the number of fringe points on this block
        do id=1, flowDoms(nn, level, sps)%nFringe
           
           ii = flowDoms(nn, level, sps)%iMesh(1, id)
           jj = flowDoms(nn, level, sps)%iMesh(2, id)
           kk = flowDoms(nn, level, sps)%iMesh(3, id)
                
           ! Set this donor point into the arrary from the petsc vector
           ii = ii + 1
           variables(iVar, nn)%arr(ii, jj, kk) = fringePtr(ii)
        end do
     end do
     
     call vecRestoreArrayF90(oversetFringes, fringePtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     
  end do masterLoop

  ! Free up the variables array...doesnt' have any data in it, just pointers
  deallocate(variables)

end subroutine wOverset
 

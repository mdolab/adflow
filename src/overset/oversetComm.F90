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
  integer(kind=intType) :: iVar, nn, i, j,k, l, ii, iFringe
  integer(kind=intType) :: iCell, jCell, kCell, nComm, ierr, sps, gid
  real(kind=realType), dimension(:), pointer :: donorPtr, fringePtr
  real(kind=realType) :: f(8)
  real(kind=realType) :: di0, di1, dj0, dj1, dk0, dk1

  ! Assume sps=1 for now
  sps = 1
  
  ! Will do each of the required communications individually for now.
  
  ! Determine the number of variables to communicate  
  iVar = 0
  do l=start, end
     iVar = iVar + 1
  end do

  if (commPressure) then 
     iVar = iVar + 1
  end if

  if (commVarGamma) then 
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

  ! Allocate the 2D "variables" derived type
  allocate(variables(nComm, nDom))

  ! Now set all the required pointers. We will need these when filling
  ! up the donors as well as when we set the fringes.

  do nn=1, nDom
     iVar = 0
     do l=start, end
        iVar = iVar + 1
        variables(iVar, nn)%arr => flowDoms(nn, level, sps)%w(:, :, :, l)
     end do

     if (commPressure) then 
        iVar = iVar + 1
        variables(iVar, nn)%arr => flowDoms(nn, level, sps)%P(:, :, :)
     end if

     if (commVarGamma) then 
        iVar = iVar + 1
        variables(iVar, nn)%arr => flowDoms(nn, level, sps)%gamma(:, :, :)
     end if

     if (commLamVis) then 
        iVar = iVar + 1
        variables(iVar, nn)%arr => flowDoms(nn, level, sps)%rlv(:, :, :)
     end if
     
     if (commEddyVis) then 
        iVar = iVar + 1
        variables(iVar, nn)%arr => flowDoms(nn, level, sps)%rev(:, :, :)
     end if
  end do

  ! Master loop over the required number of communications
  masterLoop: do iVar=1, nComm

     call VecGetArrayF90(oversetDonors, donorPtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Counter for donorPtr array
     ii = 0

     ! Loop over the number of blocks on this processor
     do nn=1, nDom
        ! We just have to put our required variable into donorPtr
        call setPointers(nn, 1, 1)
        do k=2, kl
           do j=2, jl
              do i=2, il 
                 ii = ii + 1
                 ! the plus a's are due to the pointer offset. All the
                 ! interpolation variables are double haloed and they
                 ! start at zero.
                 donorPtr(ii) = variables(iVar, nn)%arr(i+1, j+1, k+1)

              end do
           end do
        end do
     end do

     ! Return the pointer to petsc before the scatter
     call VecRestoreArrayF90(oversetDonors, donorPtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     
     ! Now perform the actual vec scatter
     call VecScatterBegin(oversetScatter, oversetDonors, oversetFringes, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     call VecScatterEnd(oversetScatter, oversetDonors, oversetFringes, &
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
     
        ! Here is where we actually do the overset interpolation. Loop
        ! over the number of fringes on my block
        
        do iFringe=1, oBlocks(nn)%nFringe
           ii = ii + 1
          
           di1 = oBlocks(nn)%donorFrac(1, iFringe)
           dj1 = oBlocks(nn)%donorFrac(2, iFringe)
           dk1 = oBlocks(nn)%donorFrac(3, iFringe)

           di0 = one - di1
           dj0 = one - dj1
           dk0 = one - dk1

           ! The corresponding linear weights
           f(1)   = di0*dj0*dk0
           f(2)   = di1*dj0*dk0
           f(3)   = di0*dj1*dk0
           f(4)   = di1*dj1*dk0
           f(5)   = di0*dj0*dk1
           f(6)   = di1*dj0*dk1
           f(7)   = di0*dj1*dk1
           f(8)   = di1*dj1*dk1

           iCell = oBlocks(nn)%fringeIndices(1, iFringe)
           jCell = oBlocks(nn)%fringeIndices(2, iFringe)
           kCell = oBlocks(nn)%fringeIndices(3, iFringe)
           
           ! The +1 is due to the pointer offset
           variables(iVar, nn)%arr(iCell+1, jCell+1, kCell+1) = &
                f(1)*fringePtr(8*ii-7) + &
                f(2)*fringePtr(8*ii-6) + &
                f(3)*fringePtr(8*ii-5) + &
                f(4)*fringePtr(8*ii-4) + &
                f(5)*fringePtr(8*ii-3) + &
                f(6)*fringePtr(8*ii-2) + &
                f(7)*fringePtr(8*ii-1) + &
                f(8)*fringePtr(8*ii  ) 
        end do
     end do

     call vecRestoreArrayF90(oversetFringes, fringePtr, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     
  end do masterLoop

  ! Free up the variables array...doesnt' have any data in it, just pointers
  deallocate(variables)

end subroutine wOverset
 

subroutine getFreeStreamResidual(rhoRes,totalRRes)
  use communication
  use ADjointVars
  use precision
  use inputTimeSpectral
  use flowVarRefState
  use inputIteration
  use blockpointers
  use iteration

  implicit none

  real(kind=realType), intent(out) :: rhoRes,totalRRes
  real(kind=realType),dimension(:), allocatable :: wtemp
  real(kind=realType) :: temp1,temp2
  integer(kind=intType) :: nDimW,ierr,counter
  integer(kind=intType) :: tempStartLevel,tempCurrentLevel,tempMGStartLevel
  integer(kind=intType) :: nn,sps,i,j,k,l

  ! Get the residual cooresponding to the free-stream on the fine grid-level

  ! We need to copy the current wvector temporirly since it may be a
  ! restart and actually useful
  nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
  allocate(wtemp(nDimW))

  ! Copy w to wTemp
  counter = 0
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    counter = counter + 1
                    wtemp(counter) = w(i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
     end do domains
  end do spectralLoop
  
  tempMGStartLevel = mgStartLevel
  tempCurrentLevel = currentLevel
  !tempStartLevel   = startLevel

  mgStartLevel = 1
  currentlevel = 1
  !startLevel   = 1

  call setUniformFlow
  call getCurrentResidual(rhoRes,totalRRes)
  counter = 0
  spectralLoop2: do sps=1,nTimeIntervalsSpectral
     domains2: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    counter = counter + 1
                    w(i,j,k,l) = wtemp(counter) 
                 enddo
              enddo
           enddo
        enddo
     end do domains2
  end do spectralLoop2

  mgStartLevel = tempMGStartLevel
  currentLevel = tempCurrentLevel
  !startLevel   = tempStartLevel

  deallocate(wtemp)
  call getCurrentResidual(temp1,temp2)

end subroutine getFreeStreamResidual

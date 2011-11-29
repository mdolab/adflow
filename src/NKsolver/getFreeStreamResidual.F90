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
  real(kind=realType),dimension(:), allocatable :: wtemp,ptemp
  real(kind=realType) :: temp1,temp2,v2,gm1
  integer(kind=intType) :: nDimW,nDimP,ierr,counter
  integer(kind=intType) :: tempStartLevel,tempCurrentLevel,tempMGStartLevel
  integer(kind=intType) :: nn,sps,i,j,k,l

  ! Get the residual cooresponding to the free-stream on the fine grid-level

  ! We need to copy the current wvector temporirly since it may be a
  ! restart and actually useful

  ! Copy out ALL states and ALL pressures including the halos. 
  nDimW = 0
  nDimP = 0
  do nn=1,nDom
     call setPointers(nn,1,1)
     nDimp = nDimp + (ib+1)*(jb+1)*(kb+1)
     nDimW = nDimw + (ib+1)*(jb+1)*(kb+1)
  end do

  nDimp = nDimp * nTimeIntervalsSpectral
  nDimw = nDimw * nTimeIntervalsSpectral * nw 

  allocate(wtemp(nDimW),ptemp(nDimP))

  ! Copy w to wTemp
  counter = 0
  spectralLoop1: do sps=1,nTimeIntervalsSpectral
     domains1: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    wtemp(counter) = w(i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
     end do domains1
  end do spectralLoop1

  ! Copy p to pTemp
  counter = 0
  spectralLoop2: do sps=1,nTimeIntervalsSpectral
     domains2: do nn=1,nDom
        call setPointers(nn,1,sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 counter = counter + 1
                 ptemp(counter) = p(i,j,k)
              enddo
           enddo
        enddo
     end do domains2
  end do spectralLoop2

  
  tempMGStartLevel = mgStartLevel
  tempCurrentLevel = currentLevel
  !tempStartLevel   = startLevel

  mgStartLevel = 1
  currentlevel = 1
  !startLevel   = 1
  
  call setUniformFlow
  call getCurrentResidual(rhoRes,totalRRes)

  counter = 0
  spectralLoop3: do sps=1,nTimeIntervalsSpectral
     domains3: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    w(i,j,k,l) = wtemp(counter) 
                 enddo
              enddo
           enddo
        enddo
     end do domains3
  end do spectralLoop3

  counter = 0
  spectralLoop4: do sps=1,nTimeIntervalsSpectral
     domains4: do nn=1,nDom
        call setPointers(nn,1,sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 counter = counter + 1
                 p(i,j,k) = ptemp(counter) 
              enddo
           enddo
        enddo
     end do domains4
  end do spectralLoop4

  mgStartLevel = tempMGStartLevel
  currentLevel = tempCurrentLevel
  !startLevel   = tempStartLevel

  deallocate(wtemp,ptemp)

end subroutine getFreeStreamResidual

subroutine getFreeStreamResidual(rhoRes,totalRRes)

  use communication
  use inputTimeSpectral
  use flowVarRefState
  use inputIteration
  use blockpointers
  use iteration

  implicit none

  real(kind=realType), intent(out) :: rhoRes,totalRRes
  real(kind=realType),dimension(:), allocatable :: wtemp, ptemp, rlvtemp, revtemp
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
  if (viscous) then
     allocate(rlvtemp(nDimP))
  end if

  if (eddyModel) then 
     allocate(revtemp(nDimP))
  end if

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

  if (viscous) then
     counter = 0
     spectralLoop3: do sps=1,nTimeIntervalsSpectral
        domains3: do nn=1,nDom
           call setPointers(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    rlvtemp(counter) = rlv(i,j,k)
                 enddo
              enddo
           enddo
        end do domains3
     end do spectralLoop3
  end if

  if (eddyModel) then
     counter = 0
     spectralLoop4: do sps=1,nTimeIntervalsSpectral
        domains4: do nn=1,nDom
           call setPointers(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    revtemp(counter) = rev(i,j,k)
                 enddo
              enddo
           enddo
        end do domains4
     end do spectralLoop4
  end if
  
  tempMGStartLevel = mgStartLevel
  tempCurrentLevel = currentLevel

  mgStartLevel = 1
  currentlevel = 1

  call setUniformFlow
  call getCurrentResidual(rhoRes,totalRRes)

  counter = 0
  redospectralLoop1: do sps=1,nTimeIntervalsSpectral
     redomains1: do nn=1,nDom
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
     end do redomains1
  end do redospectralLoop1

  counter = 0
  respectralLoop2: do sps=1,nTimeIntervalsSpectral
     redomains2: do nn=1,nDom
        call setPointers(nn,1,sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 counter = counter + 1
                 p(i,j,k) = ptemp(counter) 
              enddo
           enddo
        enddo
     end do redomains2
  end do respectralLoop2

  if (viscous) then
     counter = 0
     redospectralLoop3: do sps=1,nTimeIntervalsSpectral
        redodomains3: do nn=1,nDom
           call setPointers(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    rlv(i,j,k) = rlvtemp(counter)
                 enddo
              enddo
           enddo
        end do redodomains3
     end do redospectralLoop3
  end if

  if (eddyModel) then
     counter = 0
     redospectralLoop4: do sps=1,nTimeIntervalsSpectral
        redodomains4: do nn=1,nDom
           call setPointers(nn,1,sps)
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    counter = counter + 1
                    rev(i,j,k) = revtemp(counter)
                 enddo
              enddo
           enddo
        end do redodomains4
     end do redospectralLoop4
  end if

  mgStartLevel = tempMGStartLevel
  currentLevel = tempCurrentLevel

  deallocate(wtemp,ptemp)

  if (viscous) then
     deallocate(rlvtemp)
  end if

  if (eddyModel) then 
     deallocate(revtemp)
  end if

end subroutine getFreeStreamResidual

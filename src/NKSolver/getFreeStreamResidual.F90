subroutine getFreeStreamResidual(rhoRes, totalRRes)

  use constants
  use blockPointers, only : nDom, ib, jb, kb, w
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowVarRefState, only : nw, winf
  implicit none

  real(kind=realType), intent(out) :: rhoRes, totalRRes
  real(kind=realType), dimension(:), allocatable :: tmp
  integer(kind=intType) :: nDimW, nDimP, counter
  integer(kind=intType) :: nn, sps, i, j, k, l, n

  ! Get the residual cooresponding to the free-stream on the fine
  ! grid-level --- This saves the current values in W, P, rlv and rev
  ! and restores them when finished. 

  call getInfoSize(n)
  allocate(tmp(n))
  call getInfo(tmp, n)

  ! Set the w-variables to the ones of the uniform flow field.
  spectralLoop4b: do sps=1, nTimeIntervalsSpectral
     domains4b: do nn=1, nDom
        call setPointers(nn, 1, sps)
        do l=1,nw
           do k=0, kb
              do j=0, jb
                 do i=0, ib
                    w(i,j,k,l) = winf(l)
                 enddo
              enddo
           enddo
        end do
     end do domains4b
  end do spectralLoop4b
  
  ! Evaluate the residual now
  call computeResidualNK()
  call getCurrentResidual(rhoRes, totalRRes)

  ! Put everything back
  call setInfo(tmp, n)

  deallocate(tmp)

end subroutine getFreeStreamResidual

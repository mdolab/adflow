subroutine setReferenceVolume(level)

  use constants
  use blockPointers, only : nDom, flowDoms, il, jl, kl
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use utils, only : setPointers
  implicit none
  integer :: ierr

  integer(kind=intType), intent(in) :: level
  integer(kind=intType) :: nn, sps
  integer(kind=intType) :: i,j, k
  
  spectral: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        call setPointers(nn, level, sps)
        allocate(flowDoms(nn, level, sps)%volRef(2:il, 2:jl, 2:kl))

        do k=2, kl
           do j=2, jl
              do i=2, il
                 flowDoms(nn, level, sps)%volRef(i, j, k) = &
                      flowDoms(nn, level, sps)%vol(i, j, k)
              end do
           end do
        end do
     end do domains
  end do spectral
end subroutine setReferenceVolume

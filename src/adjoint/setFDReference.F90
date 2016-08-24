subroutine setFDReference(level)

  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use inputPhysics
  use utils, only : EChk, setPointers
  implicit none

  ! Input Parameters
  integer(kind=intType) :: level
  ! Working Parameters
  integer(kind=intType) :: i, j, k, l, nn, sps, liftIndex
  real(kind=realType) :: alpha, beta

  ! Compute the reference values for doing jacobian with FD
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        
        call setPointers(nn, level, sps)
        shockSensor => flowDoms(nn,level,sps)%shockSensor
        call block_res(nn, sps, .False., alpha, beta, liftIndex, .False.)

        ! Set the values
        do l=1, nw
           do k=0, kb 
              do j=0, jb
                 do i=0, ib
                    flowdoms(nn,1,sps)%wtmp(i,j,k,l)  = w(i, j, k, l)
                    flowdoms(nn, 1, sps)%dwtmp(i, j, k, l) = dw(i, j, k, l)
                 end do
              end do
           end do
        end do

        call initRes_block(1, nwf, nn, sps)
        
        ! Note: we have to divide by the volume for dwtmp2 since
        ! normally, dw would have been mulitpiled by 1/Vol in block_res 
        
        do l=1, nw
           do k=0, kb 
              do j=0, jb
                 do i=0, ib
                    flowdoms(nn, 1, sps)%dwtmp2(i, j, k, l) = &
                         dw(i, j, k, l)/vol(i, j, k)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine setFDReference

subroutine resetFDReference(level)

  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use inputPhysics
  use utils, only : setPointers
  implicit none

  ! Input Parameters
  integer(kind=intType) :: level
  
  ! Working Parameters
  integer(kind=intType) :: i, j, k, l, nn, sps, liftIndex
  real(kind=realType) :: alpha, beta, sepSensor, Cavitation

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Reset w and dw                            
        do l=1, nw
           do k=0, kb 
              do j=0, jb
                 do i=0, ib
                    w(i, j, k, l) = flowdoms(nn, 1, sps)%wtmp(i, j, k, l)
                    dw(i, j, k, l) = flowdoms(nn, 1, sps)%dwtmp(i, j, k, l)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine resetFDReference

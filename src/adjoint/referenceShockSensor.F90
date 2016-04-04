subroutine referenceShockSensor

  ! Compute the reference shock sensor for PC computations
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use inputPhysics
  use inputDiscretization
  implicit none

  ! Working variables
  integer(kind=intType) :: nn, level, sps, i, j, k
  
  level = 1

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)

        ! Allocate shockSensor in flowDoms *NOT* flowDomsd....and
        ! compute the value depending on equations/dissipation
        ! type. Note we are just doing all cells including corners
        ! halos..those values are not used anyway. 
        
        allocate(flowDoms(nn, 1, sps)%shockSensor(0:ib, 0:jb, 0:kb))
        shockSensor => flowDoms(nn,1,sps)%shockSensor

        if (equations == EulerEquations .or. spaceDiscr == dissMatrix) then 
           !shockSensor is Pressure
           do k=0, kb
              do j=0, jb
                 do i=0, ib
                    shockSensor(i,j,k) = P(i,j,k)
                 end do
              end do
           end do
        else
           ! Enthalpy is used instead
           do k=0, kb
              do j=2, jl
                 do i=2, il
                    shockSensor(i, j, k) = p(i, j, k)/(w(i, j, k, irho)**gamma(i, j, k))
                 enddo
              enddo
           enddo
           
           do k=2, kl
              do j=2, jl
                 shockSensor(0,  j, k) = p(0,  j, k)/(w(0,  j, k, irho)**gamma(0,  j, k))
                 shockSensor(1,  j, k) = p(1,  j, k)/(w(1,  j, k, irho)**gamma(1,  j, k))
                 shockSensor(ie, j, k) = p(ie, j, k)/(w(ie, j, k, irho)**gamma(ie, j, k))
                 shockSensor(ib, j, k) = p(ib, j, k)/(w(ib, j, k, irho)**gamma(ib, j, k))
              enddo
           enddo
           
           do k=2, kl
              do i=2, il
                 shockSensor(i, 0,  k) = p(i, 0,  k)/(w(i, 0,  k, irho)**gamma(i, 0,  k))
                 shockSensor(i, 1,  k) = p(i, 1,  k)/(w(i, 1,  k, irho)**gamma(i, 1,  k))
                 shockSensor(i, je, k) = p(i, je, k)/(w(i, je, k, irho)**gamma(i, je, k))
                 shockSensor(i, jb, k) = p(i, jb, k)/(w(i, jb, k, irho)**gamma(i, jb, k))
              enddo
           enddo
        end if
     end do
  end do
end subroutine referenceShockSensor

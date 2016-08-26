! Compute the pressure on a block with the pointers already set. This
! routine is used by the forward mode AD code only. 

subroutine computePressureSimple
  use constants
  use blockPointers
  use flowVarRefState
  use inputPhysics
  implicit none

  ! Local Variables
  integer(kind=intType) :: i, j, k, ii
  real(kind=realType) :: gm1, v2

  ! Compute the pressures
  gm1 = gammaConstant - one

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,(ib+1)*(jb+1)*(kb+1)-1
     i = mod(ii, ib+1) 
     j = mod(ii/(ib+1), (jb+1))
     k = ii/((ib+1)*(jb+1))
#else
     do k=0,kb
        do j=0,jb
           do i=0,ib
#endif             
              v2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
              p(i, j, k) = gm1*(w(i, j, k, irhoE) - half*w( i, j, k, irho)*v2)
              p(i, j, k) = max(p(i, j, k), 1.e-4_realType*pInfCorr)

#ifdef TAPENADE_FAST
           end do
#else
        end do
     end do
  end do
#endif
end subroutine computePressureSimple

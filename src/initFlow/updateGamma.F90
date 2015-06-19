subroutine updateGamma
!
!      ******************************************************************
!      *                                                                *
!      * This is a utility routine to update the gamma variable from    *
!      * from gammaConstant if gammaConstant has changed.               *
!      *                                                                *
!      ******************************************************************
!
  use blockPointers
  use inputtimespectral
  use inputPhysics
  implicit none

  integer :: nn, sps

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn, 1_intType, sps)
        gamma = gammaConstant
     end do
  end do
end subroutine updateGamma

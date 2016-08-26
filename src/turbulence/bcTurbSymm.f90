!
!       File:          bcTurbSymm.F90                                  
!       Author:        Georgi Kalitzin, Edwin van der Weide            
!       Starting date: 06-11-2003                                      
!       Last modified: 06-12-2005                                      
!
subroutine bcTurbSymm(nn)
  !
  !       bcTurbSymm applies the implicit treatment of the symmetry      
  !       boundary condition (or inviscid wall) to subface nn. As the    
  !       symmetry boundary condition is independent of the turbulence   
  !       model, this routine is valid for all models. It is assumed     
  !       that the pointers in blockPointers are already set to the      
  !       correct block on the correct grid level.                       
  !
  use constants
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, l
  !
  !       Begin execution                                                
  !
  ! Loop over the faces of the subfaces and set the values of bmt
  ! for an implicit treatment. For a symmetry face this means
  ! that the halo value is set to the internal value.

  do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
     do i=BCData(nn)%icBeg, BCData(nn)%icEnd
        do l=nt1,nt2
           select case (BCFaceID(nn))
           case (iMin)
              bmti1(i,j,l,l) = -one
           case (iMax)
              bmti2(i,j,l,l) = -one
           case (jMin)
              bmtj1(i,j,l,l) = -one
           case (jMax)
              bmtj2(i,j,l,l) = -one
           case (kMin)
              bmtk1(i,j,l,l) = -one
           case (kMax)
              bmtk2(i,j,l,l) = -one
           end select
        enddo
     enddo
  enddo
end subroutine bcTurbSymm


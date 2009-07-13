!File containing the block logical test functions for the warping algorithm

function IS_CORNER(IJK)
  
  use blockPointers
  implicit none
  !Function Arguments
  
  integer(kind=intType),dimension(3),intent(in)::ijk
  
  !Function Value
  logical:: is_corner        
  
  !Local Variables
  
  integer(kind=intType)::imax,jmax,kmax
  !
  !      ******************************************************************
  !      *                                                                *
  !      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
  !      * I,J,K VALUES IS A CORNER OF THIS BLOCK                         *
  !      *                                                                *
  !      ******************************************************************
  !         
  IMAX = il! ijk(1)
  JMAX = jl! ijk(2)
  KMAX = kl! ijk(3)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * BEGIN EXECUTION                                                *
  !      *                                                                *
  !      ******************************************************************
  !      
  IS_CORNER = .False.
  !print *,'ijk',ijk
  ! CHECKS THE EIGHT POSSIBILITIES TO SEE IF THE POINT IS A CORNER
  if (IJK(1) == 1 .and. IJK(2) == 1 .and. IJK(3) == 1) then
     IS_CORNER = .True.
  elseif (IJK(1) == IMAX .and. IJK(2) == 1 .and. IJK(3) == 1) then
     IS_CORNER = .True.
  elseif (IJK(1) == 1 .and. IJK(2) == JMAX .and. IJK(3) == 1) then
     IS_CORNER = .True.
  elseif (IJK(1) == IMAX .and. IJK(2) == JMAX .and. IJK(3) == 1) then
     IS_CORNER = .True.
  elseif (IJK(1) == 1 .and. IJK(2) == 1 .and. IJK(3) == KMAX) then
     IS_CORNER = .True.
  elseif (IJK(1) == IMAX .and. IJK(2) == 1 .and. IJK(3) == KMAX) then
     IS_CORNER = .True.
  elseif (IJK(1) == 1 .and. IJK(2) == JMAX .and. IJK(3) == KMAX ) then
     IS_CORNER = .True.
  elseif( IJK(1) == IMAX .and. IJK(2) == JMAX .and. IJK(3) == KMAX) then
     IS_CORNER = .True.
  endif
          
  return !IS_CORNER
end function IS_CORNER
   
function ON_FACE(IJK)
  use blockPointers
  implicit none
  !Function Arguments
  
  integer(kind=intType),dimension(3),intent(in)::ijk
  
  !Function Value
  logical:: on_face
  
  !Local Variables
  
  integer(kind=intType)::imax,jmax,kmax
  !
  !      ******************************************************************
  !      *                                                                *
  !      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
  !      * I,J,K VALUES IS ON A FACE OF THIS BLOCK                        *
  !      *                                                                *
  !      ******************************************************************
  !         
  !         
  IMAX = il! ijk(1)! iend-1
  JMAX = jl! ijk(2)! jend-1
  KMAX = kl! ijk(3)! kend-1
  !
  !      ******************************************************************
  !      *                                                                *
  !      * BEGIN EXECUTION                                                *
  !      *                                                                *
  !      ******************************************************************
  !      
  ON_FACE = .False.
        
  if ((IJK(1) == 1) .or. (IJK(1) == IMAX) .or. (IJK(2) == 1)&
       .or. (IJK(2) == JMAX) .or.(IJK(3) == 1) .or. (IJK(3) == KMAX))then
     ON_FACE = .True.
  endif
  
  return !ON_FACE
end function ON_FACE


function ON_EDGE(IJK)
  use blockPointers
  implicit none  
  !Function Arguments
  
  integer(kind=intType),dimension(3),intent(in)::ijk
  
  !Function Value
  logical:: on_edge
  
  !Local Variables
  
  integer(kind=intType)::imax,jmax,kmax
  !
  !      ******************************************************************
  !      *                                                                *
  !      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
  !      * BLOCK NUMBER AND I,J,K VALUES IS ON AN EDGE OF THAT BLOCK      *
  !      *                                                                *
  !      ******************************************************************
  !         
  IMAX = il! ijk(1)! iend-1
  JMAX = jl! ijk(2)! jend-1
  KMAX = kl! ijk(3)! kend-1
  !
  !      ******************************************************************
  !      *                                                                *
  !      * BEGIN EXECUTION                                                *
  !      *                                                                *
  !      ******************************************************************
  !      
  ON_EDGE = .False.

 
  if (IJK(1) == 1 .and. IJK(2) == 1)then
     ON_EDGE = .True.
  elseif (IJK(1) == IMAX .and. IJK(2) == 1)then  
     ON_EDGE = .True.
  elseif( IJK(1) == 1 .and. IJK(2) == JMAX )then
     ON_EDGE = .True.
  elseif( IJK(1) == IMAX .and. IJK(2) == JMAX)then
     ON_EDGE = .True.
  elseif( IJK(1) == 1 .and. IJK(3) == 1 )then
     ON_EDGE = .True.
  elseif( IJK(1) == IMAX .and. IJK(3) == 1)then  
     ON_EDGE = .True.
  elseif( IJK(1) == 1 .and. IJK(3) == KMAX )then
     ON_EDGE = .True.
  elseif( IJK(1) == IMAX .and. IJK(3) == KMAX)then 
     ON_EDGE = .True.
  elseif( IJK(2) == 1 .and. IJK(3) == 1)then 
     ON_EDGE = .True.
  elseif( IJK(2) == JMAX .and. IJK(3) == 1)then  
     ON_EDGE = .True.
  elseif( IJK(2) == 1 .and. IJK(3) == KMAX )then
     ON_EDGE = .True.
  elseif( IJK(2) == JMAX .and. IJK(3) == KMAX)then  
     ON_EDGE = .True.   
  END IF
          
  return !ON_EDGE
end function ON_EDGE


subroutine ON_WHICH_FACE(Ijk,which_face)
  use blockPointers
  use BCTypes
  implicit none
  !Function Arguments
  
  integer(kind=intType),dimension(3),intent(in)::ijk
  
  !Function Value
  logical,dimension(6),intent(out):: which_face
  
  !Local Variables
  integer :: i,j,k
  !
  !      ******************************************************************
  !      *                                                                *
  !      * DETERMINES WHICH FACE(S) a given point is on.                  *
  !      *                                                                *
  !      ******************************************************************
  !         
  !         
  i = ijk(1) 
  j = ijk(2)
  k = ijk(3)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * BEGIN EXECUTION                                                *
  !      *                                                                *
  !      ******************************************************************
  !      
  WHICH_FACE(:) = .False.
        
  if (I == 1)then
     which_face(imin) = .True.
  endif
  
  if (I == IL)then
     which_face(imax) = .True.
  endif
  
  if (J == 1)then
     which_face(jmin) = .True.
  endif
  
  if (J == JL)then
     which_face(jmax) = .True.
  endif
  
  if (k == 1)then
     which_face(kmin) = .True.
  endif
  
  if (K == KL)then
     which_face(kmax) = .True.
  endif
  
  
end subroutine ON_WHICH_FACE


function on_wall_surface(ijk)  
  use blockPointers
  use BCTypes
  implicit none
  !Function Arguments
  
  integer(kind=intType),dimension(3),intent(in)::ijk
  
  !Function Value
  integer:: on_wall_surface
  
  !Local Variables
  logical,dimension(6):: blockFaces
  integer::n
  on_wall_surface = 0 !false

  call on_which_face(ijk,blockFaces)
  do n=1,nSubface
     if (BCFaceID(n)== imin) then !imin
        if(blockFaces(imin))then
           if(BCType(n) == Eulerwall .or. &
                BCType(n) == NSWallAdiabatic .or. &
                BCType(n) == NSWallIsothermal) then
              on_wall_surface =1! .true.
           endif
        endif
     elseif (BCFaceID(n)== imax) then!imax
     if(blockFaces(imax))then
        if(BCType(n) == Eulerwall   .or. &
             BCType(n) == NSWallAdiabatic .or. &
             BCType(n) == NSWallIsothermal) then
           on_wall_surface =1! .true.
        endif
     endif
  elseif (BCFaceID(n)== jmin)then!:#jmin
     if(blockFaces(jmin))then
        if(BCType(n) == Eulerwall   .or. &
             BCType(n) == NSWallAdiabatic .or. &
             BCType(n) == NSWallIsothermal) then
           on_wall_surface =1! .true.
        endif
     endif
  elseif (BCFaceID(n)== jmax)then!:#jmax
     if(blockFaces(jmax))then
        if(BCType(n) == Eulerwall   .or. &
             BCType(n) == NSWallAdiabatic .or. &
             BCType(n) == NSWallIsothermal) then
           on_wall_surface = 1!.true.
        endif
     endif
  elseif(BCFaceID(n)== kmin)then!:#kmin
     if(blockFaces(kmin))then
        if(BCType(n) == Eulerwall   .or. &
             BCType(n) == NSWallAdiabatic .or. &
             BCType(n) == NSWallIsothermal) then
           on_wall_surface =1! .true.
        endif
     endif
  elseif (BCFaceID(n)== kmax)then!:#kmax
     if(blockFaces(kmax))then
        if(BCType(n) == Eulerwall   .or. &
             BCType(n) == NSWallAdiabatic .or. &
             BCType(n) == NSWallIsothermal) then
           on_wall_surface =1! .true.
        endif
     endif
  else
     print *,'Error:Not a valid face type',BCFaceID(n)
  endif
end do
return
end function on_wall_surface

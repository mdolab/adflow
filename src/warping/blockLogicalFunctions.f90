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

!
! ***********************************
! *  File: blockRelations.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-11-2008
! *  Modified: 12-11-2008
! ***********************************

subroutine blockRelations(relatedFaces, relatedEdges, edgeRelatedFaces, &
     searchPattern)
  use blockPointers
  implicit none

  ! Input/Output
  integer(kind=intType), dimension(3, 8), intent(out)::relatedFaces, relatedEdges
  integer(kind=intType), dimension(2, 12), intent(out)::edgeRelatedFaces
  integer(kind=intType), dimension(12, 6), intent(out)::searchPattern

  ! Working
  integer(kind=intType), dimension(3)::ijk

  !set the face - corner relations for WARPBLK
  relatedFaces(:, 1)=(/1, 3, 5/)
  relatedFaces(:, 2)=(/2, 3, 5/)
  relatedFaces(:, 3)=(/1, 4, 5/)
  relatedFaces(:, 4)=(/2, 4, 5/)
  relatedFaces(:, 5)=(/1, 3, 6/)
  relatedFaces(:, 6)=(/2, 3, 6/)
  relatedFaces(:, 7)=(/1, 4, 6/)
  relatedFaces(:, 8)=(/2, 4, 6/)


  !set the edge - corner relations for WARPBLK
  relatedEdges(:, 1) = (/1, 5, 9/)
  relatedEdges(:, 2) = (/1, 6, 10/)
  relatedEdges(:, 3) = (/2, 5, 11/)
  relatedEdges(:, 4) = (/2, 6, 12/)
  relatedEdges(:, 5) = (/3, 7, 9/)
  relatedEdges(:, 6) = (/3, 8, 10/)
  relatedEdges(:, 7) = (/4, 7, 11/)
  relatedEdges(:, 8) = (/4, 8, 12/)

  !set the face -edge relations for WARPBLK
  edgeRelatedFaces(:, 1) = (/3, 5/)
  edgeRelatedFaces(:, 2) = (/4, 5/)
  edgeRelatedFaces(:, 3) = (/3, 6/)
  edgeRelatedFaces(:, 4) = (/4, 6/)
  edgeRelatedFaces(:, 5) = (/1, 5/)
  edgeRelatedFaces(:, 6) = (/2, 5/)
  edgeRelatedFaces(:, 7) = (/1, 6/)
  edgeRelatedFaces(:, 8) = (/2, 6/)
  edgeRelatedFaces(:, 9) = (/1, 3/)
  edgeRelatedFaces(:, 10) = (/2, 3/)
  edgeRelatedFaces(:, 11) = (/1, 4/)
  edgeRelatedFaces(:, 12) = (/2, 4/)

  !Setup block dimensions
  ijk(1) = il
  ijk(2) = jl
  ijk(3) = kl

  !Edge search pattern
  !set to match the edge pattern in flagimplicits
  searchPattern(1, :)=(/1+1,  ijk(1)-1, 1, 1, 1, 1/)
  searchPattern(2, :)=(/1+1,  ijk(1)-1,  ijk(2),  ijk(2), 1, 1/)
  searchPattern(3, :)=(/1+1,  ijk(1)-1, 1, 1,  ijk(3),  ijk(3)/)
  searchPattern(4, :)=(/1+1,  ijk(1)-1,  ijk(2),  ijk(2),  ijk(3),  ijk(3)/)
  searchPattern(5, :)=(/1, 1, 1+1,  ijk(2)-1, 1, 1/)
  searchPattern(6, :)=(/ ijk(1),  ijk(1), 1+1,  ijk(2)-1, 1, 1/)
  searchPattern(7, :)=(/1, 1, 1+1,  ijk(2)-1,  ijk(3),  ijk(3)/)
  searchPattern(8, :)=(/ ijk(1),  ijk(1), 1+1,  ijk(2)-1,  ijk(3),  ijk(3)/)
  searchPattern(9, :)=(/1, 1, 1, 1, 1+1,  ijk(3)-1/)
  searchPattern(10, :)=(/ ijk(1),  ijk(1), 1, 1, 1+1,  ijk(3)-1/)
  searchPattern(11, :)=(/1, 1,  ijk(2),  ijk(2), 1+1,  ijk(3)-1/)
  searchPattern(12, :)=(/ ijk(1),  ijk(1),  ijk(2),  ijk(2), 1+1,  ijk(3)-1/)

end subroutine blockRelations

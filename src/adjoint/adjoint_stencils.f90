subroutine five_pt_node_stencil_i(icell,jcell,kcell,ind,CellsToDo)

  ! Return The incdices you need for a 5-pt node stencil in i-direction 

  use blockPointers

  integer(kind=intType), intent(in) :: icell,jcell,kcell
  integer(kind=intType), intent(out) :: ind(3,22),CellsToDo
  integer(kind=intType) :: counter

  ! We want i: -3 -> 2
  !         j: -2 -> 2
  !         k: -2 -> 2

  ind(:,1) = (/-2, 0, 0/)
  ind(:,2) = (/-1, 0, 0/)
  ind(:,3) = (/ 0, 0, 0/)
  ind(:,4) = (/ 1, 0, 0/)

  ind(:,5) = (/ -1,-1, 0/) ! j - 1's
  ind(:,6) = (/ 0 ,-1, 0/)

  ind(:,7) = (/ -1, 1, 0/) ! j + 1's
  ind(:,8) = (/ 0 , 1, 0/)

  ind(:, 9) = (/ -1, 0,-1/) ! k - 1's
  ind(:,10) = (/ 0 , 0,-1/)

  ind(:,11) = (/ -1, 0, 1/) ! k + 1's
  ind(:,12) = (/ 0 , 0, 1/)

  counter = 12

  if (iCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-3, 0, 0/)
  end if

  if (iCell .ne. il) then
     counter = counter + 1
     ind(:,counter) = (/ 2, 0, 0/)
  end if

  if (jCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-1,-2, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-2, 0/)
  end if

  if (jCell .ne. jl) then
     counter = counter + 1
     ind(:,counter) = (/-1, 2, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 2, 0/)
  end if

  if (kCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-1, 0,-2/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0,-2/)
  end if

  if (kCell .ne. kl) then
     counter = counter + 1
     ind(:,counter) = (/-1, 0, 2/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0, 2/)
  end if
  CellsToDo = counter

end subroutine five_pt_node_stencil_i

subroutine five_pt_node_stencil_j(icell,jcell,kcell,ind,CellsToDo)

  ! Return The incdices you need for a 5-pt node stencil in i-direction 

  use blockPointers

  integer(kind=intType), intent(in) :: icell,jcell,kcell
  integer(kind=intType), intent(out) :: ind(3,22),CellsToDo
  integer(kind=intType) :: counter

  ! We want i: -2 -> 2
  !         j: -3 -> 2
  !         k: -2 -> 2

  ind(:,1) = (/ 0,-2, 0/)
  ind(:,2) = (/ 0,-1, 0/)
  ind(:,3) = (/ 0, 0, 0/)
  ind(:,4) = (/ 0, 1, 0/)

  ind(:,5) = (/ -1,-1, 0/) ! i - 1's
  ind(:,6) = (/ -1, 0, 0/)

  ind(:,7) = (/  1,-1, 0/) ! i + 1's
  ind(:,8) = (/  1, 0, 0/)

  ind(:, 9) = (/ 0 ,-1,-1/) ! k - 1's
  ind(:,10) = (/ 0 , 0,-1/)

  ind(:,11) = (/ 0 ,-1, 1/) ! k + 1's
  ind(:,12) = (/ 0 , 0, 1/)

  counter = 12

  if (jCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/ 0, -3, 0/)
  end if

  if (jCell .ne. jl) then
     counter = counter + 1
     ind(:,counter) = (/ 0, 2, 0/)
  end if

  if (iCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-2,-1, 0/)
     counter = counter + 1
     ind(:,counter) = (/-2, 0, 0/)
  end if

  if (iCell .ne. jl) then
     counter = counter + 1
     ind(:,counter) = (/ 2,-1, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 2, 0, 0/)
  end if

  if (kCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/ 0,-1,-2/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0,-2/)
  end if

  if (kCell .ne. kl) then
     counter = counter + 1
     ind(:,counter) = (/ 0,-1, 2/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0, 2/)
  end if
  CellsToDo = counter

end subroutine five_pt_node_stencil_j

subroutine five_pt_node_stencil_k(icell,jcell,kcell,ind,CellsToDo)

  ! Return The incdices you need for a 5-pt node stencil in i-direction 

  use blockPointers

  integer(kind=intType), intent(in) :: icell,jcell,kcell
  integer(kind=intType), intent(out) :: ind(3,22),CellsToDo
  integer(kind=intType) :: counter

  ! We want i: -2 -> 2
  !         j: -2 -> 2
  !         k: -3 -> 2

  ind(:,1) = (/ 0, 0,-2/)
  ind(:,2) = (/ 0, 0,-1/)
  ind(:,3) = (/ 0, 0, 0/)
  ind(:,4) = (/ 0, 0, 1/)

  ind(:,5) = (/ -1, 0,-1/) ! i - 1's
  ind(:,6) = (/ -1, 0, 0/)

  ind(:,7) = (/  1, 0,-1/) ! i + 1's
  ind(:,8) = (/  1, 0, 0/)

  ind(:, 9) = (/ 0 ,-1,-1/) ! j - 1's
  ind(:,10) = (/ 0 ,-1, 0/)

  ind(:,11) = (/ 0 , 1,-1/) ! j + 1's
  ind(:,12) = (/ 0 , 1, 0/)

  counter = 12

  if (kCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/ 0, 0,-3/)
  end if

  if (kCell .ne. kl) then
     counter = counter + 1
     ind(:,counter) = (/ 0, 0, 2/)
  end if

  if (iCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-2, 0,-1/)
     counter = counter + 1
     ind(:,counter) = (/-2, 0, 0/)
  end if

  if (iCell .ne. jl) then
     counter = counter + 1
     ind(:,counter) = (/ 2, 0,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 2, 0, 0/)
  end if

  if (jCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/ 0,-2,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-2, 0/)
  end if

  if (jCell .ne. jl) then
     counter = counter + 1
     ind(:,counter) = (/ 0, 2,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 2, 0/)
  end if
  CellsToDo = counter

end subroutine five_pt_node_stencil_k

subroutine five_pt_node_stencil_all(icell,jcell,kcell,ind,CellsToDo)

  ! Return The incdices you need for a 5-pt node stencil in
  ! "all-directions"

  use blockPointers

  integer(kind=intType), intent(in) :: icell,jcell,kcell
  integer(kind=intType), intent(out) :: ind(3,56),CellsToDo
  integer(kind=intType) :: counter

  ! We want i: -3 -> 2
  !         j: -3 -> 2
  !         k: -3 -> 2
  ! If possible

  ! Lets try to do everything along the i-direction to help with
  ! caching

  ! Coordinates for the three cells in i-direction

  ind(:,1) = (/ -2,-1,-1/)
  ind(:,2) = (/ -1,-1,-1/)
  ind(:,3) = (/ 0 ,-1,-1/)
  ind(:,4) = (/ 1 ,-1,-1/)

  ind(:,5) = (/ -2, 0,-1/)
  ind(:,6) = (/ -1, 0,-1/)
  ind(:,7) = (/ 0 , 0,-1/)
  ind(:,8) = (/ 1 , 0,-1/)

  ind(:,9) = (/ -2,-1, 0/)
  ind(:,10) = (/ -1,-1, 0/)
  ind(:,11) = (/ 0 ,-1, 0/)
  ind(:,12) = (/ 1 ,-1, 0/)

  ind(:,13) = (/ -2, 0, 0/)
  ind(:,14) = (/ -1, 0, 0/)
  ind(:,15) = (/ 0 , 0, 0/)
  ind(:,16) = (/ 1 , 0, 0/)

  ! Coordinates for j = -1 cell

  ind(:,17) = (/ -1,-2,-1/)
  ind(:,18) = (/  0,-2,-1/)
  ind(:,19) = (/ -1,-2, 0/)
  ind(:,20) = (/  0,-2, 0/)

  ! Coordinates for j = 1 cell

  ind(:,21) = (/ -1, 1,-1/)
  ind(:,22) = (/  0, 1,-1/)
  ind(:,23) = (/ -1, 1, 0/)
  ind(:,24) = (/  0, 1, 0/)

  ! Coordinates for k = -1 cell

  ind(:,25) = (/ -1,-1,-2/)
  ind(:,26) = (/  0,-1,-2/)
  ind(:,27) = (/ -1, 0,-2/)
  ind(:,28) = (/  0, 0,-2/)

  ! Coordinates for k = 1 cell

  ind(:,29) = (/ -1,-1, 1/)
  ind(:,30) = (/  0,-1, 1/)
  ind(:,31) = (/ -1, 0, 1/)
  ind(:,32) = (/  0, 0, 1/)

  counter = 32

  if (iCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-3,-1,-1/)
     counter = counter + 1
     ind(:,counter) = (/-3, 0,-1/)
     counter = counter + 1
     ind(:,counter) = (/-3,-1, 0/)
     counter = counter + 1
     ind(:,counter) = (/-3, 0, 0/)
  end if

  if (iCell .ne.il) then
     counter = counter + 1
     ind(:,counter) = (/ 2,-1,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 2, 0,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 2,-1, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 2, 0, 0/)
  end if

  if (jCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-1,-3,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-3,-1/)
     counter = counter + 1
     ind(:,counter) = (/-1,-3, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-3, 0/)
  end if

  if (jCell .ne.jl) then
     counter = counter + 1
     ind(:,counter) = (/-1, 2,-1/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 2,-1/)
     counter = counter + 1
     ind(:,counter) = (/-1, 2, 0/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 2, 0/)
  end if

  if (kCell .ne. 2) then
     counter = counter + 1
     ind(:,counter) = (/-1,-1,-3/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-1,-3/)
     counter = counter + 1
     ind(:,counter) = (/-1, 0,-3/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0,-3/)
  end if

  if (kCell .ne.kl) then
     counter = counter + 1
     ind(:,counter) = (/-1,-1, 2/)
     counter = counter + 1
     ind(:,counter) = (/ 0,-1, 2/)
     counter = counter + 1
     ind(:,counter) = (/-1, 0, 2/)
     counter = counter + 1
     ind(:,counter) = (/ 0, 0, 2/)
  end if

  CellsToDo = counter

end subroutine five_pt_node_stencil_all

subroutine three_point_cell_stencil(ind,cellstodo)

  use blockPointers

  integer(kind=intType), intent(out) :: ind(3,7),CellsToDo
  integer(kind=intType) :: counter

  ind(:,1) = (/0,0,0/)
  ind(:,2) = (/-1,0,0/)
  ind(:,3) = (/ 1,0,0/)
  ind(:,4) = (/0,-1,0/)
  ind(:,5) = (/0, 1,0/)
  ind(:,6) = (/0,0,-1/)
  ind(:,7) = (/0,0, 1/)

  CellsToDo = 7
end subroutine three_point_cell_stencil


subroutine three_point_cell_stencil_small(icell,jcell,kcell,ind,cellstodo)

  use blockPointers
  integer(kind=intType), intent(in) :: icell,jcell,kcell
  integer(kind=intType), intent(out) :: ind(3,7),CellsToDo
  integer(kind=intType) :: counter

  ind(:,1) = (/0,0,0/)
  counter = 1
  
  ! i - 1
  if (icell > 2) then
     counter = counter + 1
     ind(:,counter) = (/-1,0,0/)
  end if

  ! i + 1
  if (icell < il) then
     counter = counter + 1
     ind(:,counter) = (/1,0,0/)
  end if

  ! j - 1
  if (jcell > 2) then
     counter = counter + 1
     ind(:,counter) = (/0,-1,0/)
  end if

  ! j + 1
  if (jcell < jl) then
     counter = counter + 1
     ind(:,counter) = (/0,1,0/)
  end if

  ! k - 1
  if (kcell > 2) then
     counter = counter + 1
     ind(:,counter) = (/0,0,-1/)
  end if

  ! k + 1
  if (kcell < kl) then
     counter = counter + 1
     ind(:,counter) = (/0,0,1/)
  end if

  CellsToDo = counter
end subroutine three_point_cell_stencil_small

subroutine five_point_cell_stencil(ind,cellstodo)

  use blockPointers

  integer(kind=intType), intent(out) :: ind(3,13),CellsToDo

  ind(:,1) = (/ 0,0,0/)

  ind(:,2) = (/-2,0,0/)
  ind(:,3) = (/-1,0,0/)
  ind(:,4) = (/ 1,0,0/)
  ind(:,5) = (/ 2,0,0/)

  ind(:,6) = (/0,-2,0/)
  ind(:,7) = (/0,-1,0/)
  ind(:,8) = (/0, 1,0/)
  ind(:,9) = (/0, 2,0/)

  ind(:,10) = (/0,0,-2/)
  ind(:,11) = (/0,0,-1/)
  ind(:,12) = (/0,0, 1/)
  ind(:,13) = (/0,0, 2/)

  CellsToDo = 13
end subroutine five_point_cell_stencil

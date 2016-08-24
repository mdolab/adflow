!
!      ******************************************************************
!      *                                                                *
!      * File:          turbBCNSWall.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-30-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine turbBCNSWall(secondHalo)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * turbBCNSWall applies the viscous wall boundary conditions      *
  !      * of the turbulent transport equations to a block. It is assumed *
  !      * that the pointers in blockPointers are already set to the      *
  !      * correct block on the correct grid level.                       *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine argument.
  !
  logical, intent(in) :: secondHalo
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, i, j, l, m

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the viscous subfaces of this block.

  bocos: do nn=1,nViscBocos

     ! Set the corresponding arrays.

     call BCTurbWall(nn)

     ! Loop over the faces and set the state in
     ! the turbulent halo cells.

     select case (BCFaceID(nn))
     case (iMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(1,i,j,l) = bvti1(i,j,l)
                 do m=nt1,nt2
                    w(1,i,j,l) = w(1,i,j,l) - bmti1(i,j,l,m)*w(2,i,j,m)
                 enddo
                 if (secondHalo) w(0,i,j,l) = w(1,i,j,l)
              end do

              if (eddyModel) then 
                 rev(1,i,j) = -rev(2,i,j)
                 if (secondHalo) then 
                    rev(0,i,j) = rev(1,i,j)
                 end if
              end if
           end do
        end do
     case (iMax) 
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(ie,i,j,l) = bvti2(i,j,l)
                 do m=nt1,nt2
                    w(ie,i,j,l) = w(ie,i,j,l) - bmti2(i,j,l,m)*w(il,i,j,m)
                 enddo
                 if (secondHalo) w(ib,i,j,l) = w(ie,i,j,l)
              end do

              if (eddyModel) then 
                 rev(ie,i,j) = -rev(il,i,j)
                 if (secondHalo) then 
                    rev(ib,i,j) = rev(ie,i,j)
                 end if
              end if
           end do
        end do
     case (jMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(i,1,j,l) = bvtj1(i,j,l)
                 do m=nt1,nt2
                    w(i,1,j,l) = w(i,1,j,l) - bmtj1(i,j,l,m)*w(i,2,j,m)
                 enddo
                 if (secondHalo) w(i,0,j,l) = w(i,1,j,l)
              end do

              if (eddyModel) then 
                 rev(i,1,j) = -rev(i,2,j)
                 if (secondHalo) then 
                    rev(i,0,j) = rev(i,1,j)
                 end if
              end if
           end do
        end do
     case (jMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(i,je,j,l) = bvtj2(i,j,l)
                 do m=nt1,nt2
                    w(i,je,j,l) = w(i,je,j,l) - bmtj2(i,j,l,m)*w(i,jl,j,m)
                 enddo
                 if (secondHalo) w(i,jb,j,l) = w(i,je,j,l)
              end do

              if (eddyModel) then 
                 rev(i,je,j) = -rev(i,jl,j)
                 if (secondHalo) then 
                    rev(i,jb,j) = rev(i,je,j)
                 end if
              end if
           end do
        end do
     case (kMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(i,j,1,l) = bvtk1(i,j,l)
                 do m=nt1,nt2
                    w(i,j,1,l) = w(i,j,1,l) - bmtk1(i,j,l,m)*w(i,j,2,m)
                 enddo
                 if (secondHalo) w(i,j,0,l) = w(i,j,1,l)
              end do

              if (eddyModel) then 
                 rev(i,j,1) = -rev(i,j,2)
                 if (secondHalo) then 
                    rev(i,j,0) = rev(i,j,1)
                 end if
              end if
           end do
        end do
     case (kMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd

              do l=nt1,nt2
                 w(i,j,ke,l) = bvtk2(i,j,l)
                 do m=nt1,nt2
                    w(i,j,ke,l) = w(i,j,ke,l) - bmtk2(i,j,l,m)*w(i,j,kl,m)
                 enddo
                 if (secondHalo) w(i,j,kb,l) = w(i,j,ke,l)
              end do

              if (eddyModel) then 
                 rev(i,j,ke) = -rev(i,j,kl)
                 if (secondHalo) then 
                    rev(i,j,kb) = rev(i,j,ke)
                 end if
              end if
           end do
        end do
     end select
  enddo bocos
end subroutine turbBCNSWall

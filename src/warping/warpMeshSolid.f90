!
! ***********************************
! *  File: warpMeshSolid.f90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************

subroutine warpMeshSolid

  ! An alternative implentation using a combination of explicit/implict aproach
  
  use blockpointers
  use communication, only: myID,sumb_comm_world
  implicit none
  integer(kind=intType)::level=1
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k
  real(kind=realType)::  points(8,3),stiffness(24,24)

  print *,'My id is:',myID,' I have ',nDom,'blocks'
  do nn=1,nDom
       call setPointers(nn,level,sps)


       IMAX = IL
       JMAX = JL
       KMAX = KL    

       points(1,:) = X(1    ,1  ,1   ,:)
       points(2,:) = X(IMAX,1   ,1   ,:)
       points(3,:) = X(1   ,JMAX,1   ,:)
       points(4,:) = X(IMAX,JMAX,1   ,:)
       points(5,:) = X(1    ,1  ,KMAX,:)
       points(6,:) = X(IMAX,1   ,KMAX,:)
       points(7,:) = X(1   ,JMAX,KMAX,:)
       points(8,:) = X(IMAX,JMAX,KMAX,:)

       call calcstiffness(points,stiffness)
    end do
      

end subroutine warpMeshSolid


subroutine calcstiffness(points,Kstif)

  implicit none

  ! Input/Output
  double precision  points(8,3)
  double precision  Kstif(24,24)

  ! Working 
  double precision c,G,E,volume,nu
  double precision Cm(6,6),Jac(3,3),Jinv(3,3),Bi(6,3),Bj(6,3)
  double precision r,s,t
  integer          i,j,k,l,m,n,ii,jj,kk,jm1,im1
  double precision Nr(8),Ns(8),Nt(8),Nx(8),Ny(8),Nz(8),intM(3,3)
  double precision nnr(2),nns(2),nnt(2),ndr(2),nds(2),ndt(2)
  double precision g_points(2),det,invdet

  ! First calculate the volume needed for the stiffness

  call volume_hexa(points,volume)

  E = 1.0/volume
  nu = 0.3
  Cm(:,:) = 0.0
  c = E/((1+nu)*(1-2*nu))
  G = E/(2*(1+nu))

  Cm(1,1) = (1-nu)*c
  Cm(1,2) = nu*c
  Cm(1,3) = nu*c

  Cm(2,1) = nu*c
  Cm(2,2) = (1-nu)*c
  Cm(2,3) = nu*c

  Cm(3,1) = nu*c
  Cm(3,2) = nu*c
  Cm(3,3) =(1-nu)*c

  Cm(4,4) = G
  Cm(5,5) = G
  Cm(6,6) = G

  ! Now here is where we do the integration numerically

  Nr(:) = 0.0
  Ns(:) = 0.0
  Nt(:) = 0.0
  Nx(:) = 0.0
  Ny(:) = 0.0
  Nz(:) = 0.0

  g_points(1) = -0.5773502691
  g_points(2) = 0.5773502691

  do k=1,2 ! Gauss in r
     do l=1,2 ! Gauss in s
        do m=1,2 ! Gauss in t

           r = g_points(k)
           s = g_points(l)
           t = g_points(m)

           nnr(1) = 0.5*( 1.0 - r )
           nnr(2) = 0.5*( 1.0 + r )
           ndr(1) = - 0.5
           ndr(2) =   0.5

           nns(1) = 0.5*( 1.0 - s )
           nns(2) = 0.5*( 1.0 + s )
           nds(1) = - 0.5
           nds(2) =   0.5

           nnt(1) = 0.5*( 1.0 - t )
           nnt(2) = 0.5*( 1.0 + t )
           ndt(1) = - 0.5
           ndt(2) =   0.5

           Nr(1) = ndr(1)*ns(1)*nt(1)
           Nr(2) = ndr(2)*ns(1)*nt(1)
           Nr(3) = ndr(2)*ns(2)*nt(1)
           Nr(4) = ndr(1)*ns(2)*nt(1)
           Nr(5) = ndr(1)*ns(1)*nt(2)
           Nr(6) = ndr(2)*ns(1)*nt(2)
           Nr(7) = ndr(2)*ns(2)*nt(2)
           Nr(8) = ndr(1)*ns(2)*nt(2)


           Ns(1) = nr(1)*nds(1)*nt(1)
           Ns(2) = nr(2)*nds(1)*nt(1)
           Ns(3) = nr(2)*nds(2)*nt(1)
           Ns(4) = nr(1)*nds(2)*nt(1)
           Ns(5) = nr(1)*nds(1)*nt(2)
           Ns(6) = nr(2)*nds(1)*nt(2)
           Ns(7) = nr(2)*nds(2)*nt(2)
           Ns(8) = nr(1)*nds(2)*nt(2)

           Nt(1) = nr(1)*ns(1)*ndt(1)
           Nt(2) = nr(2)*ns(1)*ndt(1)
           Nt(3) = nr(2)*ns(2)*ndt(1)
           Nt(4) = nr(1)*ns(2)*ndt(1)
           Nt(5) = nr(1)*ns(1)*ndt(2)
           Nt(6) = nr(2)*ns(1)*ndt(2)
           Nt(7) = nr(2)*ns(2)*ndt(2)
           Nt(8) = nr(1)*ns(2)*ndt(2)


           do n=1,8
              Jac(1,1) = Jac(1,1) + Nr(n)*points(n,1)
              Jac(1,2) = Jac(1,2) + Nr(n)*points(n,2)
              Jac(1,3) = Jac(1,3) + Nr(n)*points(n,3)

              Jac(2,1) = Jac(2,1) + Ns(n)*points(n,1)
              Jac(2,2) = Jac(2,2) + Ns(n)*points(n,2)
              Jac(2,3) = Jac(2,3) + Ns(n)*points(n,3)

              Jac(3,1) = Jac(3,1) + Nt(n)*points(n,1)
              Jac(3,2) = Jac(3,2) + Nt(n)*points(n,2)
              Jac(3,3) = Jac(3,3) + Nt(n)*points(n,3)
           end do

           ! Express the inverse of the jacobian explicitly 
           det = Jac(1,1)*Jac(2,2)*Jac(3,3)  &
           +   Jac(1,2)*Jac(2,3)*Jac(3,1)&
           +   Jac(1,3)*Jac(2,1)*Jac(3,2)&
           -   Jac(1,1)*Jac(2,3)*Jac(3,2)&
           -   Jac(1,2)*Jac(2,1)*Jac(3,3)&
           -   Jac(1,3)*Jac(2,2)*Jac(3,1)
           invdet = 1/det

           !             [ |a22 a23|   |a12 a13|  |a12 a13|]     */
           !             [ |a32 a33|  -|a32 a33|  |a22 a23|]     */
           !             [                                 ]     */
           !             [ |a21 a23|   |a11 a13|  |a11 a13|]     */
           !    A^(-1) = [-|a31 a33|   |a31 a33| -|a21 a23|] / d */
           !             [                                 ]     */
           !             [ |a21 a22|   |a11 a12|  |a11 a12|]     */
           !             [ |a31 a32|  -|a31 a32|  |a21 a22|]     */

           Jinv(1,1) =  invdet*(Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2))
           Jinv(1,2) = -invdet*(Jac(1,2)*Jac(3,3)-Jac(1,3)*Jac(3,2))
           Jinv(1,3) =  invdet*(Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2))

           Jinv(2,1) = -invdet*(Jac(2,1)*Jac(3,3)-Jac(2,3)*Jac(3,1))
           Jinv(2,2) =  invdet*(Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1))
           Jinv(2,3) = -invdet*(Jac(1,1)*Jac(2,3)-Jac(1,3)*Jac(2,1))

           Jinv(3,1) =  invdet*(Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1))
           Jinv(3,2) = -invdet*(Jac(1,1)*Jac(3,2)-Jac(1,2)*Jac(3,1))
           Jinv(3,3) =  invdet*(Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1))


           ! this calc is [die Ni/die x, die Ni/die y]T = J^-1 * [die Ni/die s , die Ni/die t]T
           do n=1,8
              Nx(n) = Jinv(1,1)*Nr(n) + Jinv(1,2)*Ns(n) + Jinv(1,3)*Nt(n)
              Ny(n) = Jinv(2,1)*Nr(n) + Jinv(2,2)*Ns(n) + Jinv(2,3)*Nt(n)
              Nz(n) = Jinv(3,1)*Nr(n) + Jinv(3,2)*Ns(n) + Jinv(3,3)*Nt(n)
           end do

         !   do i=1,8
!               do j=i,8
                 
!                  Bi(:,:) = 0.0
!                  ! Bi
!                  Bi(1,1) = Nx(i)
!                  Bi(2,2) = Ny(i)
!                  Bi(3,3) = Nz(i)

!                  Bi(4,1) = Ny(i)
!                  Bi(4,2) = Nx(i)

!                  Bi(5,2) = Nz(i)
!                  Bi(5,3) = Ny(i)

!                  Bi(6,1) = Nz(i)
!                  Bi(6,3) = Nx(i)

!                  ! Bj 
!                  Bj(1,1) = Nx(j)
!                  Bj(2,2) = Ny(j)
!                  Bj(3,3) = Nz(j)

!                  Bj(4,1) = Ny(j)
!                  Bj(4,2) = Nx(j)

!                  Bj(5,2) = Nz(j)
!                  Bj(5,3) = Ny(j)

!                  Bj(6,1) = Nz(j)
!                  Bj(6,3) = Nx(j)

!                  intM =matmul(transpose(Bi),MATMUL(Cm,Bj))
!                  im1 = i-1
!                  jm1 = j-1

!                  Kstif(3*im1 + 1, 3*jm1 + 1) = Kstif(3*im1 + 1, 3*jm1 + 1) + intM(1,1)*det
!                  Kstif(3*im1 + 1, 3*jm1 + 2) = Kstif(3*im1 + 1, 3*jm1 + 2) + intM(1,2)*det
!                  Kstif(3*im1 + 1, 3*jm1 + 3) = Kstif(3*im1 + 1, 3*jm1 + 3) + intM(1,3)*det

!                  Kstif(3*im1 + 2, 3*jm1 + 1) = Kstif(3*im1 + 2, 3*jm1 + 1) + intM(2,1)*det
!                  Kstif(3*im1 + 2, 3*jm1 + 2) = Kstif(3*im1 + 2, 3*jm1 + 2) + intM(2,2)*det
!                  Kstif(3*im1 + 2, 3*jm1 + 3) = Kstif(3*im1 + 2, 3*jm1 + 3) + intM(2,3)*det

!                  Kstif(3*im1 + 3, 3*jm1 + 1) = Kstif(3*im1 + 3, 3*jm1 + 1) + intM(3,1)*det
!                  Kstif(3*im1 + 3, 3*jm1 + 2) = Kstif(3*im1 + 3, 3*jm1 + 2) + intM(3,2)*det
!                  Kstif(3*im1 + 3, 3*jm1 + 3) = Kstif(3*im1 + 3, 3*jm1 + 3) + intM(3,3)*det

!               end do
!            end do
        end do
     end do
  end do

  ! fill up the remiander
!   do i=24,2
!      do j=i-1,1
!         Kstif(i,j) = Kstif(j,i)
!      end do
!   end do


end subroutine calcstiffness

subroutine volume_hexa(points,volume)
  implicit none
  double precision points(8,3),volume
  double precision center(3),sum,volpymrid
  integer i,idim
  ! Compute the center of the points
  ! Average x
  do idim =1,3
     sum = 0.0
     do i=1,8
        sum = sum + points(i,idim)
     end do
     center(idim) = sum / 8
  end do
    
  ! Compute the volumes of the 6 sub pyramids. The
  ! arguments of volpym must be such that for a (regular)
  ! right handed hexahedron all volumes are positive.
  
  volume = (volpymrid(center,points(1,:),points(2,:),points(4,:),points(2,:)) + &
       volpymrid(center,points(7,:),points(8,:),points(6,:),points(5,:)) + &
       volpymrid(center,points(1,:),points(3,:),points(7,:),points(5,:)) + & 
       volpymrid(center,points(4,:),points(2,:),points(6,:),points(8,:)) + &
       volpymrid(center,points(2,:),points(1,:),points(5,:),points(4,:)) + & 
       volpymrid(center,points(3,:),points(4,:),points(8,:),points(7,:)))/6

end subroutine volume_hexa

function volpymrid(p,a,b,c,d)
  implicit none
  double precision p(3),a(3),b(3),c(3),d(3),volpymrid
  
  ! 6*Volume of a pyrimid -> Counter clockwise ordering
  volpymrid = (p(1) - 0.25*(a(1) + b(1)  + c(1) + d(1))) *&
       ((a(2) - c(2))*(b(3) - d(3)) - (a(3) - c(3))*(b(2) - d(2)))   + &
       (p(2) - .25*(a(2) + b(2)  + c(2) + d(2)))*&
       ((a(3) - c(3))*(b(1) - d(1)) - (a(1) - c(1))*(b(3) - d(3)))   + &
       (p(3) - .25*(a(3) + b(3)  + c(3) + d(3)))* &
       ((a(1) - c(1))*(b(2) - d(2)) - (a(2) - c(2))*(b(1) - d(1)))
  
end function volpymrid

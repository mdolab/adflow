!
! ***********************************
! *  File: warpMeshSolid.f90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************

subroutine warpMeshSolid(g_index,gptr,l_index,lptr,l_sizes,ngi,ngptr,nli,nblock)

  ! An alternative implentation using a combination of explicit/implict aproach

  use blockpointers
  use communication, only: myID,sumb_comm_world
  implicit none
  integer(kind=intType) :: g_index(ngi),gptr(ngptr),l_index(nli),lptr(nblock+1)
  integer(kind=intType) :: l_sizes(nblock,3),indx,indy,indz,indxp1,indyp1,indzp1
  integer(kind=intType) :: ngi,nli,ngptr,nblock,nelemx,nelemy,nelemz
  integer(kind=intType) :: elemPtr(nBlock+1)

  integer(kind=intType)::level=1,norm_ind(1),maxloc,ierr
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,ii,jj,j,k,npts,dof(3),iface
  ! Temporary Variables
  real(kind=realType)::  points(8,3)
  integer(kind=intType), dimension(:,:,:,:),allocatable::bc_block
  ! Sending Variables 
  real(kind=realType) :: normal(3),delta(24)
  integer(kind=intType) :: bcs(8,3)
  real(kind=realType),  dimension(:,:,:), allocatable :: allK
  real(kind=realType),  dimension(:,:)  , allocatable :: allBCVal
  integer(kind=intType),dimension(:,:,:), allocatable :: allBC

  integer(kind=intType) ::  elemCounter ,elemID

  ! Compute total number of points
  print *,'My id is:',myID,' I have ',nDom,'blocks'

  ! FIRST compute the number of elements on each blcok (each processor does this)
  elemPtr(1) = 1
  do i=1,nBlock
     elemPtr(i+1) = elemPtr(i)+(l_sizes(i,1)-1)*(l_sizes(i,2)-1)*(l_sizes(i,3)-1)
  end do

  allocate(allK(elemPtr(nBlock+1)-1,24,24))
  allocate(allBC(elemPtr(nBlock+1)-1,24,3))
  allocate(allBCVal(elemPtr(nBlock+1)-1,24))

  elemCounter = 0
  do nn=1,nDom
     call setPointers(nn,level,sps)
     ! nbkGlobal is the original CGNS block ID we're on...if we
     ! haven't been block split this should be exactly the same
     nelemx = l_sizes(nbkGlobal,1)-1
     nelemy = l_sizes(nbkGlobal,2)-1
     nelemz = l_sizes(nbkGlobal,3)-1
     
     allocate(bc_block(l_sizes(nbkGlobal,1),l_sizes(nbkGlobal,2),l_sizes(nbkGlobal,3),3))
      
     ! Now determine the BCs for all nodes on the domain
     do i=1,l_sizes(nbkGlobal,1)
        do j =1,l_sizes(nbkGlobal,2)
           do k=1,l_sizes(nbkGlobal,3)

              ! Now we must figure out the indices
              indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
              indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
              indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

              indxp1 = int(floor(dble(i)/(nelemx)*(nx))) + 1
              indyp1 = int(floor(dble(j)/(nelemy)*(ny))) + 1
              indzp1 = int(floor(dble(k)/(nelemz)*(nz))) + 1

              bc_block(i,j,k,:) = 0

              do iface=1,nSubFace
                 if ( (indx == inBeg(iface) .and. indx == inEnd(iface)) .or. &
                      (indy == jnBeg(iface) .and. indy == jnEnd(iface)) .or. & 
                      (indz == knBeg(iface) .and. indz == knEnd(iface))) then

                    if (BCType(iface) == -1) then ! Symmetry
                       if (indx == inBeg(iface) .or. indx == inEnd(iface)) then
                          normal = BCData(iface)%norm(indy,indz,:)
                       else if (indy == jnBeg(iface) .or. indy == jnEnd(iface)) then
                          normal = BCData(iface)%norm(indx,indz,:)
                       else if (indz == knBeg(iface) .or. indz == knEnd(iface)) then
                          normal = BCData(iface)%norm(indx,indy,:)
                       end if

                       ! Now check the normal --> *SHOULD* be cartesian
                       do ii=1,3
                          normal(ii) = abs(normal(ii))
                       end do
                       norm_ind = MAXLOC(normal) !1 for x,2 for y, 3 for z
                       bc_block(i,j,k,norm_ind(1)) = 1

                    else if (BCType(iface) == -5) then
                       if (bc_block(i,j,k,1) == 0) then
                          bc_block(i,j,k,1) = -1
                       end if
                       if (bc_block(i,j,k,2) == 0) then
                          bc_block(i,j,k,2) = -1
                       end if
                       if (bc_block(i,j,k,3)== 0)then
                          bc_block(i,j,k,3)= -1
                       end if

                    else if (BCType(iface) == -6) then
                       bc_block(i,j,k,:) = 1
                    end if
                 end if
              end do ! subface loop
           end do ! k loop
        end do ! j loop
     end do ! i loop

     do i=1,nelemx   
        do j =1,nelemy
           do k=1,nelemz
              elemcounter =elemcounter + 1
              ! Now we must figure out the indices
              indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
              indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
              indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

              indxp1 = int(floor(dble(i)/(nelemx)*(nx))) + 1
              indyp1 = int(floor(dble(j)/(nelemy)*(ny))) + 1
              indzp1 = int(floor(dble(k)/(nelemz)*(nz))) + 1

              points(1,:) = X(indx  ,indy  ,indz  ,:)
              points(2,:) = X(indxp1,indy  ,indz  ,:)
              points(3,:) = X(indx  ,indyp1,indz  ,:)
              points(4,:) = X(indxp1,indyp1,indz  ,:)
              points(5,:) = X(indx  ,indy  ,indzp1,:)
              points(6,:) = X(indxp1,indy  ,indzp1,:)
              points(7,:) = X(indx  ,indyp1,indzp1,:)
              points(8,:) = X(indxp1,indyp1,indzp1,:)

              ! Compute Deltas for constrained DOF
              delta(:) = 0.0 ! Zeros the deltas by default

              ! Find the bc's that have been moved by updateFacesGlobal

              do jj=1,3
                 if (bc_block(i  ,j  ,k  ,jj) == -1) then
                    delta(0*3 + jj) = X(indx  ,indy  ,indz  ,jj)-XInit(indx  ,indy  ,indz  ,jj)
                 end if
                 if (bc_block(i+1,j  ,k  ,jj) == -1) then
                    delta(1*3 + jj) = X(indxp1,indy  ,indz  ,jj)-XInit(indxp1,indy  ,indz  ,jj)
                 end if
                 if (bc_block(i  ,j+1,k  ,jj) == -1) then
                    delta(2*3 + jj) = X(indx  ,indyp1,indz  ,jj)-XInit(indx  ,indyp1,indz  ,jj)
                 end if
                 if (bc_block(i+1,j+1,k  ,jj) == -1) then
                    delta(3*3 + jj) = X(indxp1,indyp1,indz  ,jj)-XInit(indxp1,indyp1,indz  ,jj)
                 end if
                 if (bc_block(i  ,j  ,k+1,jj) == -1) then
                    delta(4*3 + jj) = X(indx  ,indy  ,indzp1,jj)-XInit(indx  ,indy  ,indzp1,jj)
                 end if
                 if (bc_block(i+1,j  ,k+1,jj) == -1) then
                    delta(5*3 + jj) = X(indxp1,indy  ,indzp1,jj)-XInit(indxp1,indy  ,indzp1,jj)
                 end if
                 if (bc_block(i  ,j+1,k+1,jj) == -1) then
                    delta(6*3 + jj) = X(indx  ,indyp1,indzp1,jj)-XInit(indx  ,indyp1,indzp1,jj)
                 end if
                 if (bc_block(i+1,j+1,k+1,jj) == -1) then
                    delta(7*3 + jj) = X(indxp1,indyp1,indzp1,jj)-XInit(indxp1,indyp1,indzp1,jj)
                 end if
              end do
          
              !Calculate the element ID for the Bcast
              elemID = elemPtr(nBkGlobal) + (i-1)*(nelemy*nelemz) + (j-1)*(nelemz) + (k-1)

              ! Compute Stiffness
              call calcStiffness(points,allK(elemID,:,:)
              allBCval(elemID,:) = delta

              call mpi_bcast(allK(elemID,:,:),24*24,sumb_real,myID,SUmb_comm_world,ierr)
              call mpi_bcast(allBCVal(elemID,:),24,sumb_real,myID,SUmb_comm_world,ierr)
              
              if (ierr .ne. 0) then
                 print *,'Error in mpi_bast. Number is:',ierr
              end if
              
           end do ! k loop
        end do ! j loop
     end do ! iloop
      deallocate(bc_block)
   end do !nDoms loop
   call mpi_barrier(sumb_comm_world, ierr) ! Make sure everything gets caught up

   ! Now we have all the data we need on every processor so we can go
   ! right to town to assemble and solve the structural system.


   deallocate(allK,allBC,allBCVal)
end subroutine warpMeshSolid

subroutine calcstiffness(points,Kstif)
  use precision 
  implicit none

  ! Input/Output
  real(kind=realType) :: points(8,3)
  real(kind=realType) :: Kstif(24,24)

  ! Working 
  real(kind=realType) :: c,G,E,volume,nu,onemnuc,nuc
  real(kind=realType) :: Jac(3,3),Jinv(3,3)
  integer(kind=intType):: i,j,k,l,m,n,ii,jj,kk,jm1_3,im1_3
  real(kind=realType) :: r,s,t
  real(kind=realType) :: Nr(8),Ns(8),Nt(8),Nx(8),Ny(8),Nz(8)
  real(kind=realType) :: nnr(2),nns(2),nnt(2),ndr(2),nds(2),ndt(2)
  real(kind=realType) :: g_points(2),det,invdet

  ! First calculate the volume needed for the stiffness

  call volume_hexa(points,volume)

  E = 1.0/volume
  nu = 0.3
  c = E/((1+nu)*(1-2*nu))
  G = E/(2*(1+nu))

  onemnuc = (1-nu)*c
  nuc = nuc

  ! Now here is where we do the integration numerically

  g_points(1) = -0.5773502691
  g_points(2) = 0.5773502691

  ! 2x2x2 Gaussian Integration
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

           Nr(1) = ndr(1)*nns(1)*nnt(1)
           Nr(2) = ndr(2)*nns(1)*nnt(1)
           Nr(3) = ndr(1)*nns(2)*nnt(1)
           Nr(4) = ndr(2)*nns(2)*nnt(1)
           Nr(5) = ndr(1)*nns(1)*nnt(2)
           Nr(6) = ndr(2)*nns(1)*nnt(2)
           Nr(7) = ndr(1)*nns(2)*nnt(2)
           Nr(8) = ndr(2)*nns(2)*nnt(2)

           Ns(1) = nnr(1)*nds(1)*nnt(1)
           Ns(2) = nnr(2)*nds(1)*nnt(1)
           Ns(3) = nnr(1)*nds(2)*nnt(1)
           Ns(4) = nnr(2)*nds(2)*nnt(1)
           Ns(5) = nnr(1)*nds(1)*nnt(2)
           Ns(6) = nnr(2)*nds(1)*nnt(2)
           Ns(7) = nnr(1)*nds(2)*nnt(2)
           Ns(8) = nnr(2)*nds(2)*nnt(2)

           Nt(1) = nnr(1)*nns(1)*ndt(1)
           Nt(2) = nnr(2)*nns(1)*ndt(1)
           Nt(3) = nnr(1)*nns(2)*ndt(1)
           Nt(4) = nnr(2)*nns(2)*ndt(1)
           Nt(5) = nnr(1)*nns(1)*ndt(2)
           Nt(6) = nnr(2)*nns(1)*ndt(2)
           Nt(7) = nnr(1)*nns(2)*ndt(2)
           Nt(8) = nnr(2)*nns(2)*ndt(2)

           Jac(:,:) = 0.0
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

           ! This does transpose(Bi)*C*Bj
           do i=1,8
              do j=i,8

                 im1_3 = (i-1)*3
                 jm1_3 = (j-1)*3

                 Kstif(im1_3 + 1, jm1_3 + 1) = Kstif(im1_3 + 1, jm1_3 + 1) + &
                      (Nx(i)*onemnuc*Nx(j)+Ny(i)*G*Ny(j)+Nz(i)*G*Nz(j))*det
                 Kstif(im1_3 + 1, jm1_3 + 2) = Kstif(im1_3 + 1, jm1_3 + 2) + &
                      (Nx(i)*nuc*Ny(j)+Ny(i)*G*Nx(j))*det
                 Kstif(im1_3 + 1, jm1_3 + 3) = Kstif(im1_3 + 1, jm1_3 + 3) + &
                      (Nx(i)*nuc*Nz(j)+Nz(i)*G*Nx(j))*det

                 Kstif(im1_3 + 2, jm1_3 + 1) = Kstif(im1_3 + 2, jm1_3 + 1) + &
                      (Ny(i)*nuc*Nx(j)+Nx(i)*G*Ny(j))*det
                 Kstif(im1_3 + 2, jm1_3 + 2) = Kstif(im1_3 + 2, jm1_3 + 2) + &
                      (Ny(i)*onemnuc*Ny(j)+Nx(i)*G*Nx(j)+Nz(i)*G*Nz(j))*det
                 Kstif(im1_3 + 2, jm1_3 + 3) = Kstif(im1_3 + 2, jm1_3 + 3) + &
                      (Ny(i)*nuc*Nz(j)+Nz(i)*G*Ny(j))*det

                 Kstif(im1_3 + 3, jm1_3 + 1) = Kstif(im1_3 + 3, jm1_3 + 1) + &
                      (Nz(i)*nuc*Nx(j)+Nx(i)*G*Nz(j))*det
                 Kstif(im1_3 + 3, jm1_3 + 2) = Kstif(im1_3 + 3, jm1_3 + 2) + &
                      (Nz(i)*nuc*Ny(j)+Ny(i)*G*Nz(j))*det
                 Kstif(im1_3 + 3, jm1_3 + 3) = Kstif(im1_3 + 3, jm1_3 + 3) + &
                      (Nz(i)*onemnuc*Nz(j)+Ny(i)*G*Ny(j)+Nx(i)*G*Nx(j))*det

              end do
           end do
        end do
     end do
  end do
end subroutine calcstiffness

subroutine volume_hexa(points,volume)
  use precision 

  implicit none
  real(kind=realType) :: points(8,3),volume
  real(kind=realType) :: center(3),sum,volpymrid
  integer(kind=intType) ::  i,idim
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
  use precision 
  implicit none
  real(kind=realType) :: p(3),a(3),b(3),c(3),d(3),volpymrid

  ! 6*Volume of a pyrimid -> Counter clockwise ordering
  volpymrid = (p(1) - 0.25*(a(1) + b(1)  + c(1) + d(1))) *&
       ((a(2) - c(2))*(b(3) - d(3)) - (a(3) - c(3))*(b(2) - d(2)))   + &
       (p(2) - .25*(a(2) + b(2)  + c(2) + d(2)))*&
       ((a(3) - c(3))*(b(1) - d(1)) - (a(1) - c(1))*(b(3) - d(3)))   + &
       (p(3) - .25*(a(3) + b(3)  + c(3) + d(3)))* &
       ((a(1) - c(1))*(b(2) - d(2)) - (a(2) - c(2))*(b(1) - d(1)))

end function volpymrid

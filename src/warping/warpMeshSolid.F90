!
! ***********************************
! *  File: warpMeshSolid.F90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************

subroutine warpMeshSolid(g_index,gptr,l_index,lptr,l_sizes,ngi,ngptr,nli,nblock, &
     nFree,nSurface,nBoundary)
#ifndef USE_NO_PETSC

  use blockpointers
  use communication, only: myID,sumb_comm_world,sumb_comm_self
  use solidwarpmodule

  implicit none
  ! Input Data
  integer(kind=intType) :: g_index(ngi),gptr(ngptr)
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: ngi,nli,ngptr,nblock
  integer(kind=intType) :: nFree,nSurface,nBoundary

  ! Working Data
  integer(kind=intType) :: nelemx,nelemy,nelemz,nelem
  integer(kind=intType) :: indx,indy,indz,indxp1,indyp1,indzp1
  integer(kind=intType) :: elemPtr(nBlock+1),n

  integer(kind=intType)::level=1,ierr
  integer(kind=intType)::nn,sps=1,i,ii,j,jj,k,npts,idim,jdim,irow,jcol
  integer(kind=intType) :: indices(8,3),rows(3),cols(3)
  ! Temporary Variables
  real(kind=realType)::  points(8,3)

  ! Sending Variables 
  real(kind=realType),  dimension(:,:,:), allocatable :: allK
  real(kind=realType),  dimension(:,:)  , allocatable :: allBCVal
  real(kind=realType),  dimension(:)    , allocatable :: solution
  integer(kind=intType) ::  elemCounter ,elemID

  ! Petsc sizes
  integer(kind=intType) :: nonz
  integer(kind=intType),dimension(:), allocatable :: nnz
  ! Compute total number of points

  external MyKSPMonitor


  print *,'My id is:',myID,' I have ',nDom,'blocks'
  call initPETScWrap()

  ! FIRST compute the number of elements on each blcok (each processor does this)
  elemPtr(1) = 1
  do i=1,nBlock
     elemPtr(i+1) = elemPtr(i)+(l_sizes(i,1)-1)*(l_sizes(i,2)-1)*(l_sizes(i,3)-1)
  end do

  allocate(allK(elemPtr(nBlock+1)-1,24,24))
  allocate(allBCVal(elemPtr(nBlock+1)-1,24))

  elemCounter = 0
  do nn=1,nDom
     call setPointers(nn,level,sps)
     ! nbkGlobal is the original CGNS block ID we're on...if we
     ! haven't been block split this should be exactly the same
     
     nelemx = l_sizes(nbkGlobal,1)-1
     nelemy = l_sizes(nbkGlobal,2)-1
     nelemz = l_sizes(nbkGlobal,3)-1
     
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

              !Calculate the element ID for the Bcast
              elemID = elemPtr(nBkGlobal) + (i-1)*(nelemy*nelemz) + (j-1)*(nelemz) + (k-1)
              ! Compute Deltas for constrained DOF
              do jj=1,3
                 allBCVal(elemID,0*3 + jj) = X(indx  ,indy  ,indz  ,jj)-XInit(indx  ,indy  ,indz  ,jj)
                 allBCVal(elemID,1*3 + jj) = X(indxp1,indy  ,indz  ,jj)-XInit(indxp1,indy  ,indz  ,jj)
                 allBCVal(elemID,2*3 + jj) = X(indx  ,indyp1,indz  ,jj)-XInit(indx  ,indyp1,indz  ,jj)
                 allBCVal(elemID,3*3 + jj) = X(indxp1,indyp1,indz  ,jj)-XInit(indxp1,indyp1,indz  ,jj)
                 allBCVal(elemID,4*3 + jj) = X(indx  ,indy  ,indzp1,jj)-XInit(indx  ,indy  ,indzp1,jj)
                 allBCVal(elemID,5*3 + jj) = X(indxp1,indy  ,indzp1,jj)-XInit(indxp1,indy  ,indzp1,jj)
                 allBCVal(elemID,6*3 + jj) = X(indx  ,indyp1,indzp1,jj)-XInit(indx  ,indyp1,indzp1,jj)
                 allBCVal(elemID,7*3 + jj) = X(indxp1,indyp1,indzp1,jj)-XInit(indxp1,indyp1,indzp1,jj)
              end do
          
              ! Compute Stiffness
              call calcStiffness(points,allK(elemID,:,:))

              call mpi_bcast(allK(elemID,:,:),24*24,sumb_real,myID,SUmb_comm_world,ierr)
              call mpi_bcast(allBCVal(elemID,:),24,sumb_real,myID,SUmb_comm_world,ierr)
              
              if (ierr .ne. 0) then
                 print *,'Error in mpi_bast. Number is:',ierr
              end if
              
           end do ! k loop
        end do ! j loop
     end do ! iloop
   end do !nDoms loop
   call mpi_barrier(sumb_comm_world, ierr) ! Make sure everything gets caught up

   ! Now we have all the data we need on every processor so we can go
   ! right to town to assemble and solve the structural system.

   ! Initialize the PETSc Matrices

   ! --------------- Kuu ----------------
   nrow = nfree*3
   allocate(nnz(nfree*3))
   nnz(:) = 27 ! Petsc is screwed up...we can't just pass in the singe,
               ! nz value we MUST pass in the full nnz array
   nonz = 27

   !call MatCreateSeqSBAIJ(PETSC_COMM_SELF,3,nrow,nrow,nonz,nnz,kuu,ierr)
   call MatCreate(PETSC_COMM_SELF,kuu,ierr)
   call MatSetSizes(kuu,nrow,nrow,nrow,nrow,ierr)
   call MatSetType(kuu,'seqsbaij',ierr)
   call MatSeqSBAIJSetPreallocation(kuu,3,nonz,nnz,ierr)
   call MatSetOption(kuu,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)

   ! --------------- Kus ----------------
   ncol = nSurface*3
   !call MatCreateSeqBAIJ(PETSC_COMM_SELF,3,nrow,ncol,nonz,nnz,kus,ierr)
   call MatCreate(PETSC_COMM_SELF,kus,ierr)
   call MatSetSizes(kus,nrow,ncol,nrow,ncol,ierr)
   call MatSetType(kus,'seqbaij',ierr)
   call MatSeqBAIJSetPreallocation(kus,3,nonz,nnz,ierr)
   call MatSetOption(kus,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)

   deallocate(nnz)

   ! --------------- uu,us,fu -----------------
   call VecCreate(PETSC_COMM_SELF,uu,ierr)
   call VecSetSizes(uu,nrow,nrow,ierr)
   call VecSetType(uu,'seq',ierr)

   call VecCreate(PETSC_COMM_SELF,us,ierr)
   call VecSetSizes(us,ncol,ncol,ierr)
   call VecSetType(us,'seq',ierr)
   call VecSetBlockSize(us,3,ierr)

   call VecCreate(PETSC_COMM_SELF,fu,ierr)
   call VecSetSizes(fu,nrow,nrow,ierr)
   call VecSetType(fu,'seq',ierr)

   print *,'Assembling'
   print *,'nFree,nSurface,nBoundary:',nfree,nsurface,nboundary
   
   do nn=1,nblock
      nelemx = l_sizes(nn,1)-1
      nelemy = l_sizes(nn,2)-1
      nelemz = l_sizes(nn,3)-1
      
      do i=1,nelemx   
         do j =1,nelemy
            do k=1,nelemz
               elemID = elemPtr(nn)+(i-1)*(nelemy*nelemz)+(j-1)*(nelemz)+(k-1)
               call hexa_index(i-1,j-1,k-1,indices) ! Make this zero based
              
               do ii=1,8
                  irow = l_index(lptr(nn) +&
                       indices(ii,1)*l_sizes(nn,2)*l_sizes(nn,3) + &
                       indices(ii,2)*l_sizes(nn,3) +&
                       indices(ii,3)+1)
                  
                  if (irow .ge. nFree .and. irow .lt. (nFree+Nsurface)) then
                     call VecSetValuesBlocked(us,1,irow-nFree, &
                     allBCVal(elemID,(ii-1)*3+1:ii*3),INSERT_VALUES,ierr)
                  end if

                  do jj=ii,8
                     jcol = l_index(lptr(nn) +&
                          indices(jj,1)*l_sizes(nn,2)*l_sizes(nn,3) + &
                          indices(jj,2)*l_sizes(nn,3) +&
                          indices(jj,3)+1)

                     if (irow .lt. nFree .and. jcol .lt. nFree) then
                        call MatSetValuesBlocked(Kuu,1,irow,1,jcol, &
                             allK(elemID,(ii-1)*3+1:ii*3,(jj-1)*3+1:jj*3), &
                             ADD_VALUES,ierr)
                     end if

                     if (irow .lt. nFree .and. jcol .ge. nFree .and. &
                          jcol .lt. (nFree+nSurface)) then
                        call MatSetValuesBlocked(Kus,1,irow,1,jcol-nFree, &
                             allK(elemID,(ii-1)*3+1:ii*3,(jj-1)*3+1:jj*3), &
                             ADD_VALUES,ierr)
                     end if
                  end do ! j node loop
               end do ! inode loop
            end do !z elem loop
         end do ! yelem loop
      end do ! xelem loop
   end do !ndomain loop

   call MatAssemblyBegin(Kuu,MAT_FINAL_ASSEMBLY,ierr)
   call MatAssemblyEnd(Kuu,MAT_FINAL_ASSEMBLY,ierr)
   call MatAssemblyBegin(Kus,MAT_FINAL_ASSEMBLY,ierr)
   call MatAssemblyEnd(Kus,MAT_FINAL_ASSEMBLY,ierr)
   call VecAssemblyBegin(us,ierr)
   call VecAssemblyEnd(us,ierr)
   print *,'Done Assmebly'
   call mpi_barrier(sumb_comm_world, ierr) 
   if (ierr .ne. 0) then 
      print *,'Error creating Stiffness Matrix' 
   end if

   ! Multiply kus by us to get -Fu
   PETScNegOne = -1.0
   call MatMult(Kus,us,fu,ierr)
   call VecScale(fu,PETScNegOne,ierr)

   ! Ksp 
    call KSPCreate(PETSC_COMM_SELF,ksp, ierr)
    call KSPSetOperators(ksp,kuu,kuu,DIFFERENT_NONZERO_PATTERN,ierr)
    call KSPSetFromOptions(ksp, ierr)
    call KSPSetType(ksp, "preonly", ierr)
    call KSPSetTolerances(ksp, 1e-12,1e-16,1e3,250,ierr) !reltol,abstol,
                                                         !divergence tol,
                                                         !max iter
    call KSPGetPC(ksp, pc,ierr)
    call PCSetType( pc, "cholesky",ierr)
    !call PCSetType( pc, "icc",ierr)
    !print *,'pc factor ordering'
    !call PCFactorSetMatOrderingtype(pc,"rcm",ierr)
        call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
                         PETSC_NULL_FUNCTION, ierr)
    
    call KSPSetUp(ksp,ierr)
    call KSPSolve(ksp,fu,uu,ierr)

    allocate(solution(nfree*3))
    do i=1,nfree*3
       call VecGetValues(uu,1,i-1,ierr,solution(i))
    end do

    call writeFEAP('solid_test.dat',g_index,gptr,l_index,lptr,l_sizes,ngi,&
         ngptr,nli,nblock,elemPtr(nBlock+1)-1,solution,nfree)

    !call VecView(uu,PETSC_VIEWER_DRAW_WORLD,ierr)
    !call VecView(uu,PETSC_VIEWER_STDOUT_SELF,ierr)
      !pause

    deallocate(allK,allBCVal,solution)
    call MatDestroy(Kuu,ierr)
    call MatDestroy(Kus,ierr)
    call VecDestroy(uu,ierr)
    call VecDestroy(us,ierr)
    call VecDestroy(fu,ierr)
   
#endif
end subroutine warpMeshSolid

subroutine writeFEAP(file_name,g_index,gptr,l_index,lptr,l_sizes,ngi,ngptr,&
  nli,nblock,nelem,solution,nfree)
  use blockpointers

  implicit none
  character*32 file_name
  integer(kind=intType) :: g_index(ngi),gptr(ngptr)
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: ngi,nli,ngptr,nblock,nelem,nfree
  real(kind=realType)   :: solution(nfree*3)
  integer(kind=intType) :: ii,i,j,k,blockID,nelemx,nelemy,nelemz,indx,indy,indz
  integer(kind=intType) :: indices(8,3),nn

  OPEN (7, FILE = file_name)
  WRITE(7,*) 'TITLE = "SOLID WARP Finite Element Data"'
  WRITE(7,*) 'VARIABLES = "CoordinateX", "CoordinateY",', &
       '"CoordinateZ"'
  WRITE(7,*) 'Zone N=',ngptr-1,'E=', nelem
  WRITE(7,*)  'DATAPACKING=POINT, ZONETYPE = FEBRICK'
    
  ! ONLY WORKS ON ONE PROCESSOR
  print *,'About to do x'
  do ii=1,ngptr-1 ! Global number of points
     blockID = g_index(gptr(ii)+1)+1 ! gptr is zero based
     i       = g_index(gptr(ii)+2)+1
     j       = g_index(gptr(ii)+3)+1
     k       = g_index(gptr(ii)+4)+1
     call setPointers(blockID,1,1)
     nelemx = l_sizes(nbkGlobal,1)-1
     nelemy = l_sizes(nbkGlobal,2)-1
     nelemz = l_sizes(nbkGlobal,3)-1

     indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
     indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
     indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

     !write(7,*) Xinit(indx,indy,indz,1),Xinit(indx,indy,indz,2),Xinit(indx,indy,indz,3)

     if (ii .le. nfree) then
        write(7,*) Xinit(indx,indy,indz,1)+solution(ii*3-2),&
             Xinit(indx,indy,indz,2)+solution(ii*3-1),&
             Xinit(indx,indy,indz,3)+solution(ii*3  )
     else
        write(7,*) X(indx,indy,indz,1),&
             X(indx,indy,indz,2),&
             X(indx,indy,indz,3)
     end if
  end do

  print *,'About to do connectivity'
  ! Connectivity (CCW Ordered)
  do nn=1,nblock
 
     do i=1,l_sizes(nn,1)-1
        do j =1,l_sizes(nn,2)-1
           do k=1,l_sizes(nn,3)-1

              call hexa_index_ccw(i-1,j-1,k-1,indices)
              do ii=1,8
                 write(7,'(I5)') l_index(lptr(nn) + &
                      indices(ii,1)*l_sizes(nn,2)*l_sizes(nn,3) + &
                      indices(ii,2)*l_sizes(nn,3) + &
                      indices(ii,3)+1)+1
              end do
           end do
        end do
     end do
  end do
  close(7)
end subroutine writeFEAP


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
  E = 1.0
  nu = 0.3
  c = E/((1+nu)*(1-2*nu))
  G = E/(2*(1+nu))

  onemnuc = (1-nu)*c
  nuc = nu*c

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
   ! Fill up the remainder
   do i=1,24
      do j = 1,i
         Kstif(i,j) = Kstif(j,i)
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

subroutine hexa_index(i,j,k,indices)
  use precision
  integer(kind=intType) :: i,j,k,indices(8,3)
  
  indices(1,:) = (/i,j,k/)
  indices(2,:) = (/i+1,j,k/)
  indices(3,:) = (/i,j+1,k/)
  indices(4,:) = (/i+1,j+1,k/)
  indices(5,:) = (/i,j,k+1/)
  indices(6,:) = (/i+1,j,k+1/)
  indices(7,:) = (/i,j+1,k+1/)
  indices(8,:) = (/i+1,j+1,k+1/)
end subroutine hexa_index

subroutine hexa_index_ccw(i,j,k,indices)
  use precision
  integer(kind=intType) :: i,j,k,indices(8,3)
  
  indices(1,:) = (/i,j,k/)
  indices(2,:) = (/i+1,j,k/)
  indices(3,:) = (/i+1,j+1,k/)
  indices(4,:) = (/i,j+1,k/)
  indices(5,:) = (/i,j,k+1/)
  indices(6,:) = (/i+1,j,k+1/)
  indices(7,:) = (/i+1,j+1,k+1/)
  indices(8,:) = (/i,j+1,k+1/)
end subroutine hexa_index_ccw

!   allocate(bc_block(l_sizes(nbkGlobal,1),l_sizes(nbkGlobal,2),l_sizes(nbkGlobal,3),3))
      
!      ! Now determine the BCs for all nodes on the domain
!      do i=1,l_sizes(nbkGlobal,1)
!         do j =1,l_sizes(nbkGlobal,2)
!            do k=1,l_sizes(nbkGlobal,3)

!               ! Now we must figure out the indices
!               indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
!               indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
!               indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

!               indxp1 = int(floor(dble(i)/(nelemx)*(nx))) + 1
!               indyp1 = int(floor(dble(j)/(nelemy)*(ny))) + 1
!               indzp1 = int(floor(dble(k)/(nelemz)*(nz))) + 1

!               bc_block(i,j,k,:) = 0

!               do iface=1,nSubFace
!                  if ( (indx == inBeg(iface) .and. indx == inEnd(iface)) .or. &
!                       (indy == jnBeg(iface) .and. indy == jnEnd(iface)) .or. & 
!                       (indz == knBeg(iface) .and. indz == knEnd(iface))) then

!                     if (BCType(iface) == -1) then ! Symmetry
!                        if (indx == inBeg(iface) .or. indx == inEnd(iface)) then
!                           normal = BCData(iface)%norm(indy,indz,:)
!                        else if (indy == jnBeg(iface) .or. indy == jnEnd(iface)) then
!                           normal = BCData(iface)%norm(indx,indz,:)
!                        else if (indz == knBeg(iface) .or. indz == knEnd(iface)) then
!                           normal = BCData(iface)%norm(indx,indy,:)
!                        end if

!                        ! Now check the normal --> *SHOULD* be cartesian
!                        do ii=1,3
!                           normal(ii) = abs(normal(ii))
!                        end do
!                        norm_ind = MAXLOC(normal) !1 for x,2 for y, 3 for z
!                        bc_block(i,j,k,norm_ind(1)) = 1

!                     else if (BCType(iface) == -5) then
!                        if (bc_block(i,j,k,1) == 0) then
!                           bc_block(i,j,k,1) = -1
!                        end if
!                        if (bc_block(i,j,k,2) == 0) then
!                           bc_block(i,j,k,2) = -1
!                        end if
!                        if (bc_block(i,j,k,3)== 0)then
!                           bc_block(i,j,k,3)= -1
!                        end if

!                     else if (BCType(iface) == -6) then
!                        bc_block(i,j,k,:) = 1
!                     end if
!                  end if
!               end do ! subface loop
!            end do ! k loop
!         end do ! j loop
!      end do ! i loop



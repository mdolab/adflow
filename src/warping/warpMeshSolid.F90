!
! ***********************************
! *  File: warpMeshSolid.F90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************

subroutine warpMeshSolid(g_index,gptr,l_index,lptr,l_sizes,ngi,ngptr,nli,nblock, &
     nFree,nSurface,nBoundary,free_dof)
#ifndef USE_NO_PETSC

  use blockpointers
  use communication, only: myID,sumb_comm_world,sumb_comm_self
  use solidwarpmodule

  implicit none
  ! Input Data
  integer(kind=intType) :: g_index(ngi),gptr(ngptr)
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: ngi,nli,ngptr,nblock
  integer(kind=intType) :: nFree,nSurface,nBoundary,counter,free_dof(3*nfree)

  ! Working Data
  integer(kind=intType) :: nelemx,nelemy,nelemz,nelem
  integer(kind=intType) :: indx,indy,indz,indxp1,indyp1,indzp1
  integer(kind=intType) :: elemPtr(nBlock+1),n

  integer(kind=intType) :: level=1,ierr
  integer(kind=intType) :: nn,sps=1,i,ii,iii,j,jj,jjj,k,kk,kkk,iset,jset
  integer(kind=intType) :: npts,idim,jdim,irow,jcol
  integer(kind=intType) :: indices(8,3),rows(3),cols(3),lenx,leny,lenz,indices2(8)
  ! Temporary Variables
  real(kind=realType)::  points(8,3),value(1,1),value1(1,1),value2(1,1)
  real(kind=realType):: deltas(8,3),volume

  ! Sending Variables 
  real(kind=realType),  dimension(:,:,:), allocatable :: allK
  real(kind=realType),  dimension(:,:)  , allocatable :: allBCVal
  real(kind=realType),  dimension(:)    , allocatable :: solution
  integer(kind=intType) :: elemID,blockID

  ! Petsc sizes
  integer(kind=intType) :: nonz
  integer(kind=intType),dimension(:), allocatable :: nnz
  ! Compute total number of points

  ! Temporary Mesh Arrays
  real(kind=realType), dimension(:,:,:,:), allocatable :: Xstart
  real(kind=realType), dimension(:,:,:,:), allocatable :: Xupdated
  real(kind=realType), dimension(:,:,:,:), allocatable :: SS
  real(kind=realType) ::  shp(8),pt_delta(3),nns(2),nnt(2),nnr(2)

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

  open(7,FILE='FEAP_input.inp')
  write(7,*)"FEAP * * Solid Element Element Example"
  write(7,*) nfree+nsurface+nboundary,elemPtr(nBlock+1)-1,1,3,3,8
  write(7,*) " "
 

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
              ! Now we must figure out the indices
              indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1 ! note nx is number of elems, or mesh size-1
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

              call volume_hexa(points,volume)
            !   write(7,*) "MATErial ",elemID
!               write(7,*) "SOLID"
!               write(7,*) "ELAStic ISOtripoic ",1/volume,0.3
!               write(7,*) " "


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

  write(7,*) "MATErial ",1
  write(7,*) "SOLID"
  write(7,*) "ELAStic ISOtropic ",1,0.3
  write(7,*) " "

  call mpi_barrier(sumb_comm_world, ierr) ! Make sure everything gets caught up

  write(7,*) "COORdinate ALL"
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

     write(7,*) ii,0,Xinit(indx,indy,indz,1),Xinit(indx,indy,indz,2),Xinit(indx,indy,indz,3)
  
  end do
  write(7,*) " "
  write(7,*) "ELEMents"
        ! Connectivity (CCW Ordered)
  counter = 1
  do nn=1,nblock
     do i=1,l_sizes(nn,1)-1
        do j =1,l_sizes(nn,2)-1
           do k=1,l_sizes(nn,3)-1
              call hexa_index_ccw(i-1,j-1,k-1,indices)
              do ii=1,8
                 indices2(ii) = l_index(lptr(nn) + &
                      indices(ii,1)*l_sizes(nn,2)*l_sizes(nn,3) + &
                      indices(ii,2)*l_sizes(nn,3) + &
                      indices(ii,3)+1)+1
              end do
              write(7,*) counter,1,1,indices2(1),indices2(2),indices2(3),indices2(4),&
                   indices2(5),indices2(6),indices2(7),indices2(8)
              counter = counter + 1
           end do
        end do
     end do
  end do

  write(7,*) " "
  write(7,*) "BOUNdary restraints"
  
  ! Now do the boundary restraints
  ! Set the prescribed displacements and the fixed ones
  do i=nfree+nsurface+1,nfree+nsurface+nboundary
     write(7,*) i,0,1,1,1
  end do
  write(7,*)
  write(7,*) "DISPlacement restraints"
  
  ! Now do the boundary restraints
  ! Set the prescribed displacements and the fixed ones
  
  do ii=nfree+1,nfree+nsurface ! Global number of points
     blockID = g_index(gptr(ii)+1)+1 ! gptr is zero based
     i       = g_index(gptr(ii)+2)+1
     j       = g_index(gptr(ii)+3)+1
     k       = g_index(gptr(ii)+4)+1

     call setPointers(blockID,1,1)
     pt_delta = X(i,j,k,:)-Xinit(i,j,k,:)
     write(7,*) ii,0,pt_delta(1),pt_delta(2),pt_delta(3)
  end do

  write(7,*) " "
  write(7,*) "END"
  write(7,*) "BATCh"
  write(7,*) "TANGent"
  write(7,*) "FORM"
  write(7,*) "SOLV"
  write(7,*) "DISPlacement all"
  write(7,*) "END"
  write(7,*) "STOP"
  close(7)

  ! Now we have all the data we need on every processor so we can go
  ! right to town to assemble and solve the structural system.

  ! Initialize the PETSc Matrices

  ! --------------- Kuu ----------------
  nrow = nfree*3
  allocate(nnz(nfree*3))
  nnz(:) = min(27,nfree) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(27,nfree)

  !call MatCreateSeqSBAIJ(PETSC_COMM_SELF,3,nrow,nrow,nonz,nnz,kuu,ierr)
  call MatCreate(PETSC_COMM_SELF,kuu,ierr)
  call MatSetSizes(kuu,nrow,nrow,nrow,nrow,ierr)
  call MatSetType(kuu,'seqaij',ierr)
  call MatSeqBAIJSetPreallocation(kuu,3,nonz,nnz,ierr)
  call MatSetOption(kuu,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  deallocate(nnz)
  ! --------------- Kus ----------------
  ncol = nSurface*3
  !call MatCreateSeqBAIJ(PETSC_COMM_SELF,3,nrow,ncol,nonz,nnz,kus,ierr)

  allocate(nnz(nfree*3))
  nnz(:) = min(27,nsurface) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(27,nsurface)
  call MatCreate(PETSC_COMM_SELF,Kus,ierr)
  call MatSetSizes(Kus,nrow,ncol,nrow,ncol,ierr)
  call MatSetType(Kus,'seqbaij',ierr)
  call MatSeqBAIJSetPreallocation(Kus,3,nonz,nnz,ierr)
  call MatSetOption(Kus,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
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
  print *,'fre dof'
  do i=1,3*nfree
     print *,i,free_dof(i)
  end do
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
                         allBCVal(elemID,3*ii-2:ii*3),INSERT_VALUES,ierr)
                 end if

                 do jj=1,8
                    jcol = l_index(lptr(nn) +&
                         indices(jj,1)*l_sizes(nn,2)*l_sizes(nn,3) + &
                         indices(jj,2)*l_sizes(nn,3) +&
                         indices(jj,3)+1)

                    if (irow .lt. nFree .and. jcol .lt. nFree .and. jcol .ge. jcol) then
                       do iii=1,3
                          do jjj=1,3
                             if (free_dof(3*irow+iii) == -1) then
                                iset = -1
                                print *,'iset -1'
                             else
                                iset = 3*irow+iii-1
                             end if

                             if (free_dof(3*jcol+jjj) == -1) then
                                jset = -1
                             else
                                jset = 3*jcol+jjj-1
                             end if

                             call MatSetValues(Kuu,1,iset,1,jset, &
                                  allK(elemID,3*(ii-1)+iii,(jj-1)*3+jjj), &
                                  ADD_VALUES,ierr)
                          end do
                       end do
                    end if

                    if (irow .lt. nFree .and. jcol .ge. nFree .and. &
                         jcol .lt. (nFree+nSurface)) then
                       call MatSetValuesBlocked(Kus,1,irow,1,jcol-nFree, &
                            allK(elemID,ii*3-2:ii*3,jj*3-2:jj*3), &
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

  ! Multiply kus by us to get -Fu
  PETScNegOne = -1.0
  call MatMult(Kus,us,fu,ierr)
  call VecScale(fu,PETScNegOne,ierr)

  !    open(7,file="rhs.txt")
  !    do i=1,nfree*3
  !       call VecGetValues(fu,1,i-1,value1,ierr)
  !       write(7,*) value1
  !    end do
  !    close(7)

  !     open(7,file='matrix.txt')
  !     do i=1,nfree*3
  !        do j=1,nfree*3
  !           call MatGetValues(kuu,1,i-1,1,j-1,value1,ierr)
  !           write(7,*) value1
  !        end do
  !     end do
  !     close(7)
  !     call VecNorm(fu,NORM_2,value1,ierr)

  ! Do the cholesky directly, no ksp
  call MatGetOrdering(kuu,"rcm",row_perm,col_perm,ierr)
  call MatLUFactor(kuu,row_perm,col_perm,info,ierr)
  call MatSolve(kuu,fu,uu,ierr)

  open(7,file="uu.txt")
  do i=1,nfree*3
     call VecGetValues(uu,1,i-1,value1,ierr)
     write(7,*) value1
  end do
  close(7)

  print *,'Done Solution'

  OPEN (7, FILE = 'zones.dat')
  WRITE(7,*) 'TITLE = "SOLID WARP Finite Element Data"'
  WRITE(7,*) 'VARIABLES = "X", "Y","Z"'


  do nn=1,nDom ! This is now done for each processor
     call setPointers(nn,1,1)

     nelemx = l_sizes(nbkGlobal,1)-1
     nelemy = l_sizes(nbkGlobal,2)-1
     nelemz = l_sizes(nbkGlobal,3)-1
     !print *,'nn:',nn


     do i=1,nelemx   
        do j =1,nelemy
           do k=1,nelemz

              elemID = elemPtr(nbkglobal)+(i-1)*(nelemy*nelemz)+(j-1)*(nelemz)+(k-1)
              ! Now we must figure out the indices
              indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
              indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
              indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

              indxp1 = int(floor(dble(i)/(nelemx)*(nx))) + 1
              indyp1 = int(floor(dble(j)/(nelemy)*(ny))) + 1
              indzp1 = int(floor(dble(k)/(nelemz)*(nz))) + 1

              lenx = indxp1-indx+1
              leny = indyp1-indy+1
              lenz = indzp1-indz+1

              write(7,*) 'Zone I=',lenx,' J=',leny,' K=',lenz
              write(7,*) 'DATAPACKING=POINT\n'

              allocate(Xstart(lenx,leny,lenz,3))
              allocate(SS(lenx,leny,lenz,3))

              XStart = Xinit(indx:indxp1,indy:indyp1,indz:indzp1,:)

              call para3d(Xstart,lenx,leny,lenz,3,SS)

              ! Now get the deltas from the FE solution
              call hexa_index(i-1,j-1,k-1,indices)
              do ii=1,8
                 irow = l_index(lptr(nbkglobal) +&
                      indices(ii,1)*l_sizes(nbkGlobal,2)*l_sizes(nbkGlobal,3)+&
                      indices(ii,2)*l_sizes(nbkGlobal,3) +&
                      indices(ii,3)+1) !! 0 based

                 if (irow .lt. nFree) then
                    call VecGetValues(uu,1,3*irow,deltas(ii,1),ierr)
                    call VecGetValues(uu,1,3*irow+1,deltas(ii,2),ierr)
                    call VecGetValues(uu,1,3*irow+2,deltas(ii,3),ierr)

                 elseif (irow .lt. nFree+nSurface) then
                    deltas(ii,:) = allBCval(elemID,3*ii-2:ii*3)
                 else
                    deltas(ii,:) = 0.0
                 end if
              end do
              ! Now actuall update the Mesh Variables 'X'
              do ii=1,lenx
                 do jj=1,leny
                    do kk=1,lenz

                       nnr(1) =  1.0 - SS(ii,jj,kk,1)
                       nnr(2) =  SS(ii,jj,kk,1)

                       nns(1) = 1.0 - SS(ii,jj,kk,2)
                       nns(2) = SS(ii,jj,kk,2)

                       nnt(1) = 1.0 - SS(ii,jj,kk,3)
                       nnt(2) = SS(ii,jj,kk,3)

                       shp(1) = nnr(1)*nns(1)*nnt(1)
                       shp(2) = nnr(2)*nns(1)*nnt(1)
                       shp(3) = nnr(1)*nns(2)*nnt(1)
                       shp(4) = nnr(2)*nns(2)*nnt(1)
                       shp(5) = nnr(1)*nns(1)*nnt(2)
                       shp(6) = nnr(2)*nns(1)*nnt(2)
                       shp(7) = nnr(1)*nns(2)*nnt(2)
                       shp(8) = nnr(2)*nns(2)*nnt(2)

                       pt_delta = matmul(shp,deltas)

                       X(indx+ii-1,indy+jj-1,indz+kk-1,:) =  &
                            Xinit(indx+ii-1,indy+jj-1,indz+kk-1,:) + pt_delta

                    end do ! kk loop
                 end do ! jj loop
              end do !ii loop

              do kk=1,lenz
                 do jj=1,leny
                    do ii=1,lenx
                       write(7,*) X(indx+ii-1,indy+jj-1,indz+kk-1,1),&
                            X(indx+ii-1,indy+jj-1,indz+kk-1,2),&
                            X(indx+ii-1,indy+jj-1,indz+kk-1,3)
                    end do
                 end do
              end do

              deallocate(Xstart,SS)

           end do ! k loop
        end do ! j loop
     end do ! i loop
  end do

  close(7)


  allocate(solution(nfree*3))
  do i=1,nfree*3
     call VecGetValues(uu,1,i-1,solution(i),ierr)
     print *,free_dof(i),solution(i)
  end do

  call writeFEAP(g_index,gptr,l_index,lptr,l_sizes,ngi,&
       ngptr,nli,nblock,elemPtr(nBlock+1)-1,solution,nfree)

  ! Deallocation
  deallocate(allK,allBCVal,solution)
  call MatDestroy(Kuu,ierr)
  call MatDestroy(Kus,ierr)
  call VecDestroy(uu,ierr)
  call VecDestroy(us,ierr)
  call VecDestroy(fu,ierr)
  call mpi_barrier(sumb_comm_world, ierr) 
#endif
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
  integer(kind=intType):: i,j,k,l,m,n,ii,jj,kk
  real(kind=realType) :: r,s,t
  real(kind=realType) :: Nr(8),Ns(8),Nt(8),Nx(8),Ny(8),Nz(8)
  real(kind=realType) :: nnr(2),nns(2),nnt(2),ndr(2),nds(2),ndt(2)
  real(kind=realType) :: g_points(2),det,invdet


  real(kind=realType) :: Bi(6,3),Bj(6,3),Cm(6,6),intM(3,3)

  ! First calculate the volume needed for the stiffness

  call volume_hexa(points,volume)

  E = 1.0/(volume)!*volume)
  !E = 1
  nu = .2
  c = E/((1+nu)*(1-2*nu))
  G = c/(2-nu)

  onemnuc = (1.0-nu)*c
  nuc = nu*c
  Cm(:,:) = 0.0
  Cm(1,1) = onemnuc
  Cm(1,2) = nuc
  Cm(1,3) = nuc

  Cm(2,1) = nuc
  Cm(2,2) = onemnuc
  Cm(2,3) = nuc

  Cm(3,1) = nuc
  Cm(3,2) = nuc
  Cm(3,3) = onemnuc

  Cm(4,4) = G
  Cm(5,5) = G
  Cm(6,6) = G
  Kstif(:,:) = 0.0

  !print *,'E,nu,c,g,onemnuc,nuc:',E,nu,c,G,onemnuc,nuc
  ! Now here is where we do the integration numerically

  g_points(1) = -0.5773502691
  g_points(2) =  0.5773502691

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
           !             [ |a21\ a23|   |a11 a13|  |a11 a13|]     */
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

           do i=1,8
              do j=i,8

                 Bi(:,:) = 0.0
                 Bj(:,:) = 0.0

                 Bi(1,1) = Nx(i)
                 Bi(2,2) = Ny(i)
                 Bi(3,3) = Nz(i)

                 Bi(4,1) = Ny(i)
                 Bi(4,2) = Nx(i)

                 Bi(5,2) = Nz(i)
                 Bi(5,3) = Ny(i)

                 Bi(6,1) = Nz(i)
                 Bi(6,3) = Nx(i)

                 !Bj 
                 Bj(1,1) = Nx(j)
                 Bj(2,2) = Ny(j)
                 Bj(3,3) = Nz(j)

                 Bj(4,1) = Ny(j)
                 Bj(4,2) = Nx(j)

                 Bj(5,2) = Nz(j)
                 Bj(5,3) = Ny(j)

                 Bj(6,1) = Nz(j)
                 Bj(6,3) = Nx(j)


                 intM = matmul(transpose(Bi), matmul(Cm,Bj))

                 Kstif(3*i-2,3*j-2) = Kstif(3*i-2,3*j-2) +   intM(1,1)*det
                 Kstif(3*i-2,3*j-1) = Kstif(3*i-2,3*j-1) +   intM(1,2)*det
                 Kstif(3*i-2,3*j  ) = Kstif(3*i-2,3*j  ) +   intM(1,3)*det

                 Kstif(3*i-1,3*j-2) = Kstif(3*i-1,3*j-2) +   intM(2,1)*det
                 Kstif(3*i-1,3*j-1) = Kstif(3*i-1,3*j-1) +   intM(2,2)*det
                 Kstif(3*i-1,3*j  ) = Kstif(3*i-1,3*j  ) +   intM(2,3)*det

                 Kstif(3*i  ,3*j-2) = Kstif(3*i  ,3*j-2) +   intM(3,1)*det
                 Kstif(3*i  ,3*j-1) = Kstif(3*i  ,3*j-1) +   intM(3,2)*det
                 Kstif(3*i  ,3*j  ) = Kstif(3*i  ,3*j  ) +   intM(3,3)*det
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

subroutine writeFEAP(g_index,gptr,l_index,lptr,l_sizes,ngi,ngptr,&
     nli,nblock,nelem,solution,nfree)
  use blockpointers

  implicit none
  integer(kind=intType) :: g_index(ngi),gptr(ngptr)
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: ngi,nli,ngptr,nblock,nelem,nfree
  real(kind=realType)   :: solution(nfree*3)
  integer(kind=intType) :: ii,i,j,k,blockID,nelemx,nelemy,nelemz,indx,indy,indz
  integer(kind=intType) :: indices(8,3),nn

  OPEN (7, FILE ='solid_test.dat')
  WRITE(7,*) 'TITLE = "SOLID WARP Finite Element Data"'
  WRITE(7,*) 'VARIABLES = "X", "Y","Z"'
  WRITE(7,*) 'Zone N=',ngptr-1,'E=', nelem
  WRITE(7,*)  'DATAPACKING=POINT, ZONETYPE = FEBRICK'

  ! ONLY WORKS ON ONE PROCESSOR
  print *,'About to do x'
  print *,'nfree:',nfree
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
        write(7,*) Xinit(indx,indy,indz,1)+solution(ii*3-2)/1.,&
             Xinit(indx,indy,indz,2)+solution(ii*3-1)/1.,&
             Xinit(indx,indy,indz,3)+solution(ii*3  )/1.
     else
        write(7,*) X(indx,indy,indz,1),&
             X(indx,indy,indz,2),&
             X(indx,indy,indz,3)

     end if
  end do

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

  ! Zone with DOF it THINKS are free
  WRITE(7,*) 'Zone T=free I=',nfree
  write(7,*) 'DATAPACKING=POINT'

  do ii=1,nfree
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

     write(7,*) Xinit(indx,indy,indz,1),Xinit(indx,indy,indz,2),Xinit(indx,indy,indz,3)
  end do

  close(7)
end subroutine writeFEAP

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



!   open(7,file='kus.txt')
!    do j=1,nsurface*3
!       do i=1,nfree*3
!          call MatGetValues(Kus,1,i-1,1,j-1,value1,ierr)
!          write(7,*) value1
!       end do 
!    end do
!    close(7)


subroutine para3d(X,n,m,l,ndim,S)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract para3d calculates the parametric locations for a 3d block
  !
  !     Description of Arguments
  !     Input
  !     X       - Real,size(n,m,l,ndim): Coordiantes
  !     Output
  !     S       - Real,size(n,m,l,3): The u,v,w parametric positions

  implicit none

  ! Input
  integer          , intent(in)   :: n,m,l,ndim
  double precision , intent(in)   :: X(n,m,l,ndim)


  ! Output
  double precision , intent(out)  :: S(n,m,l,ndim)

  ! Working 
  integer                         :: i,j,k

  double precision DELI,DELJ,DELK

  DELI(I,J,K) = SQRT ((X(I,J,K,1) - X(I-1,J,K,1)) ** 2 + &
       (X(I,J,K,2) - X(I-1,J,K,2)) ** 2 + &
       (X(I,J,K,3) - X(I-1,J,K,3)) ** 2)

  DELJ(I,J,K) = SQRT ((X(I,J,K,1) - X(I,J-1,K,1)) ** 2 + &
       (X(I,J,K,2) - X(I,J-1,K,2)) ** 2 + &
       (X(I,J,K,3) - X(I,J-1,K,3)) ** 2)

  DELK(I,J,K) = SQRT ((X(I,J,K,1) - X(I,J,K-1,1)) ** 2 + &
       (X(I,J,K,2) - X(I,J,K-1,2)) ** 2 + &
       (X(I,J,K,3) - X(I,J,K-1,3)) ** 2)


  !     Zero the three low-end faces (or edges if one plane is specified).

  DO K = 1, l
     DO J = 1, m
        S(1,J,K,1) = 0.0
     END DO

     DO I = 1, n
        S(I,1,K,2) = 0.0
     END DO
  END DO

  DO J = 1, m
     DO I = 1, n
        S(I,J,1,3) = 0.0
     END DO
  END DO

  !     Set up the low-end edge lines because they are missed by the
  !     following loops over most of the low-end faces:

  DO I = 2, n
     S(I,1,1,1) = S(I-1,1,1,1) + DELI(I,1,1)
  END DO

  DO J = 2, m
     S(1,J,1,2) = S(1,J-1,1,2) + DELJ(1,J,1)
  END DO

  DO K = 2, l
     S(1,1,K,3) = S(1,1,K-1,3) + DELK(1,1,K)
  END DO

  !     Set up the rest of the low-end face lines because they are
  !     missed by the the main loop over most of the volume.

  DO K = 2, l
     DO J = 2, m
        S(1,J,K,2) = S(1,J-1,K,2) + DELJ(1,J,K)
        S(1,J,K,3) = S(1,J,K-1,3) + DELK(1,J,K)
     END DO
     DO I = 2, n
        S(I,1,K,1) = S(I-1,1,K,1) + DELI(I,1,K)
        S(I,1,K,3) = S(I,1,K-1,3) + DELK(I,1,K)
     END DO

  END DO

  DO J = 2, m
     DO I = 2, n
        S(I,J,1,1) = S(I-1,J,1,1) + DELI(I,J,1)
        S(I,J,1,2) = S(I,J-1,1,2) + DELJ(I,J,1)
     END DO
  END DO

  !     Traverse the block just once for all lines except those within
  !     the low-end faces.

  DO K = 2, l
     DO J = 2, m
        DO I = 2, n
           S(I,J,K,1) = S(I-1,J,K,1) + DELI(I,J,K)
           S(I,J,K,2) = S(I,J-1,K,2) + DELJ(I,J,K)
           S(I,J,K,3) = S(I,J,K-1,3) + DELK(I,J,K)
        END DO
     END DO
  END DO

  !     Normalizing requires another pass through the volume.
  !     Handle lines of zero length first by inserting uniform
  !     distributions.  Then the standard normalization can be
  !     applied safely everywhere.

  DO K = 1, l

     !        Zero-length lines in the I direction?

     DO J = 1, m
        IF (S(n,J,K,1) == 0.0) THEN
           DO I = 2, n
              S(I,J,K,1) = I - 1
           END DO
        END IF
     END DO

     !        Zero-length lines in the J direction?

     DO I = 1, n
        IF (S(I,m,K,2) == 0.0) THEN
           DO J = 2, m
              S(I,J,K,2) = J - 1
           END DO
        END IF
     END DO
  END DO

  !     Zero-length lines in the K direction?

  DO J = 1, m
     DO I = 1, n
        IF (S(I,J,l,3) == 0.0) THEN
           DO K = 2, l
              S(I,J,K,3) = K - 1
           END DO
        END IF
     END DO
  END DO

  !     Normalize:

  DO K = 1, l
     DO J = 1, m
        DO I = 1, n
           S(I,J,K,1) = S(I,J,K,1) / S(n,J,K,1)
           S(I,J,K,2) = S(I,J,K,2) / S(I,m,K,2)
           S(I,J,K,3) = S(I,J,K,3) / S(I,J,l,3)
        END DO
     END DO
  END DO

  !     Finally, precise 1s for the three high-end faces:

  DO K = 1, l
     DO J = 1, m
        S(n,J,K,1) = 1.0
     END DO

     DO I = 1, n
        S(I,m,K,2) = 1.0
     END DO
  END DO

  DO J = 1, m
     DO I = 1, n
        S(I,J,l,3) = 1.0
     END DO
  END DO



end subroutine para3d


!    call MatCholeskyFactor(kuu,col_perm,10,ierr)
!    ! Ksp 
!     call KSPCreate(PETSC_COMM_SELF,ksp, ierr)
!     call KSPSetOperators(ksp,kuu,kuu,DIFFERENT_NONZERO_PATTERN,ierr)
!     call KSPSetFromOptions(ksp, ierr)
!     call KSPSetType(ksp, "preonly", ierr) !preonly
!     call KSPSetTolerances(ksp, 1e-12,1e-16,1e3,250,ierr) !reltol,abstol,
!                                                          !divergence tol,
!                                                          !max iter
!     call KSPGetPC(ksp, pc,ierr)
!     call PCSetType( pc, "cholesky",ierr)
!     !call PCSetType( pc, "icc",ierr)
!     call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
!          PETSC_NULL_FUNCTION, ierr)

!     call KSPSolve(ksp,fu,uu,ierr)

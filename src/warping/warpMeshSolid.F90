
! ***********************************
! *  File: warpMeshSolid.F90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************
subroutine warpMeshSolid(nuu,nus,l_index,lptr,l_sizes,nli,nblock)
#ifndef USE_NO_PETSC

  use blockpointers
  use communication, only: myID,sumb_comm_world,sumb_comm_self,nProc
  use solidwarpmodule
  use mdData    

  implicit none
  ! Input Data
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: nli,nblock,nuu,nus

  ! Working Data
  integer(kind=intType) :: nelemx,nelemy,nelemz
  integer(kind=intType) :: indx,indy,indz,indxp1,indyp1,indzp1
  integer(kind=intType) :: lenx,leny,lenz
  integer(kind=intType) :: iproc,index
  integer(kind=intType) :: level=1,sps=1,ierr,imax,jmax,kmax

  integer(kind=intType) :: nn,i,ii,iii,j,jj,jjj,k,kk,kkk
  integer(kind=intType) :: irow,jcol,blockID,counter,indices(8,3)

  !integer(kind=intType),dimension(6)::IFACEPTB
  !integer(kind=intType),dimension(12)::IEDGEPTB  

  ! Temporary Variables
  real(kind=realType)   ::  points(8,3),deltas(8,3)
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0

  ! Temporary Mesh Arrays
  real(kind=realType), dimension(:,:,:,:), allocatable :: Xstart,Xupdated,SS
  real(kind=realType) ::  shp(8),pt_delta(3),nns(2),nnt(2),nnr(2)

  external MyKSPMonitor

  call initPETScWrap()
  call initializeWarpMeshSolid(nuu,nus,l_index,lptr,l_sizes,nli,nblock)

  ! --------------- Compute stiffness matrices on each processor -----------
  counter = 0
  do nn=1,nDom
     call setPointers(nn,level,sps)
     ! nbkGlobal is the original CGNS block ID we're on

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
              points(1,:) = Xinit(indx  ,indy  ,indz  ,:)
              points(2,:) = Xinit(indxp1,indy  ,indz  ,:)
              points(3,:) = Xinit(indx  ,indyp1,indz  ,:)
              points(4,:) = Xinit(indxp1,indyp1,indz  ,:)
              points(5,:) = Xinit(indx  ,indy  ,indzp1,:)
              points(6,:) = Xinit(indxp1,indy  ,indzp1,:)
              points(7,:) = Xinit(indx  ,indyp1,indzp1,:)
              points(8,:) = Xinit(indxp1,indyp1,indzp1,:)

              ! Compute Deltas for constrained DOF
              do jj=1,3
                 localBCVal(counter*24      + jj) = X(indx  ,indy  ,indz  ,jj)-XInit(indx  ,indy  ,indz  ,jj)
                 localBCVal(counter*24 + 3  + jj) = X(indxp1,indy  ,indz  ,jj)-XInit(indxp1,indy  ,indz  ,jj)
                 localBCVal(counter*24 + 6  + jj) = X(indx  ,indyp1,indz  ,jj)-XInit(indx  ,indyp1,indz  ,jj)
                 localBCVal(counter*24 + 9  + jj) = X(indxp1,indyp1,indz  ,jj)-XInit(indxp1,indyp1,indz  ,jj)
                 localBCVal(counter*24 + 12 + jj) = X(indx  ,indy  ,indzp1,jj)-XInit(indx  ,indy  ,indzp1,jj)
                 localBCVal(counter*24 + 15 + jj) = X(indxp1,indy  ,indzp1,jj)-XInit(indxp1,indy  ,indzp1,jj)
                 localBCVal(counter*24 + 18 + jj) = X(indx  ,indyp1,indzp1,jj)-XInit(indx  ,indyp1,indzp1,jj)
                 localBCVal(counter*24 + 21 + jj) = X(indxp1,indyp1,indzp1,jj)-XInit(indxp1,indyp1,indzp1,jj)
              end do

              ! Compute Stiffness
              call calcStiffness(points,localK(counter*24*24+1))
              counter = counter + 1
           end do ! k loop
        end do ! j loop
     end do ! iloop
  end do !nDoms loop
  print *,'done matrices'
  ! ------------------ Communicate all stifness matrices to each other processor ----------

  call mpi_allgatherv(localK,allNelem(myID+1)*24*24,sumb_real,allK,Krecvcount,Kdispls,sumb_real,SUmb_comm_world,ierr)
  call mpi_allgatherv(localBCVal,allNelem(myID+1)*24,sumb_real,allBCVal,BCrecvcount,BCdispls,sumb_real,SUmb_comm_world,ierr)

  ! ------------------ Assemble Global Stiffness Matrix on each Processor ----------------

  counter = 0
  do iproc =1,nProc
     do nn=1,allnDom(iproc)
        blockID = allBlockIDs(cumNDom(iproc)+nn)
        nelemx = l_sizes(blockID,1)-1
        nelemy = l_sizes(blockID,2)-1
        nelemz = l_sizes(blockID,3)-1
        do i=1,nelemx   
           do j =1,nelemy
              do k=1,nelemz
                 call hexa_index(i-1,j-1,k-1,indices) ! Make this zero based
                 do ii=1,8
                    do iii=1,3
                       irow = l_index(lptr(blockID) +&
                            indices(ii,1)*l_sizes(blockID,2)*l_sizes(blockID,3)*3 + &
                            indices(ii,2)*l_sizes(blockID,3)*3 +&
                            indices(ii,3)*3+ iii )

                       if (irow .ge. nuu) then
                          call VecSetValues(us,1,irow-nuu,allBCVal(counter*24+3*(ii-1)+iii),&
                               INSERT_VALUES,ierr)
                       end if

                       do jj=1,8
                          do jjj =1,3
                             jcol = l_index(lptr(blockID) +&
                                  indices(jj,1)*l_sizes(blockID,2)*l_sizes(blockID,3)*3 + &
                                  indices(jj,2)*l_sizes(blockID,3)*3 +&
                                  indices(jj,3)*3 + jjj )
                             if (irow .lt. nuu) then
                                if (jcol .lt. nuu) then
                                   call MatSetValues(Kuu,1,irow,1,jcol, &
                                        allK(counter*24*24+(3*(ii-1)+iii-1)*24+3*(jj-1)+jjj), &
                                        ADD_VALUES,ierr)
                                else
                                   call MatSetValues(Kus,1,irow,1,jcol-nuu, &
                                        allK(counter*24*24+(3*(ii-1)+iii-1)*24+3*(jj-1)+jjj), &
                                        ADD_VALUES,ierr)
                                end if
                             end if
                          end do ! jjj dof loop
                       end do ! j node loop
                    end do ! iii dof loop
                 end do ! inode loop
                 counter = counter + 1
              end do !z elem loop
           end do ! yelem loop
        end do ! xelem loop
     end do !ndomain loop
  end do !iproc loop

  ! ------------------- Do petsc assmebly and solve -------------
  call MatAssemblyBegin(Kuu,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Kuu,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(Kus,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Kus,MAT_FINAL_ASSEMBLY,ierr)
  call VecAssemblyBegin(us,ierr)
  call VecAssemblyEnd(us,ierr)

  ! Multiply kus by us to get -Fu
  call MatMult(Kus,us,fu,ierr)
  call VecScale(fu,PETScNegOne,ierr)
  
  call MatGetOrdering(kuu,"nd",row_perm,col_perm,ierr)
  call MatLUFactor(kuu,row_perm,col_perm,info,ierr)  
  call MatSolve(kuu,fu,uu,ierr)

  ! ------------------ Set the Volume Coordinates from Solution -----------------

  counter = cumNElem(myID+1)
  do nn=1,nDom ! This is now done for each processor
     call setPointers(nn,level,sps)

     nelemx = l_sizes(nbkGlobal,1)-1
     nelemy = l_sizes(nbkGlobal,2)-1
     nelemz = l_sizes(nbkGlobal,3)-1
     
     do i=1,nelemx   
        do j =1,nelemy
           do k=1,nelemz

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

              allocate(Xstart(lenx,leny,lenz,3))
              allocate(SS(lenx,leny,lenz,3))

              XStart = Xinit(indx:indxp1,indy:indyp1,indz:indzp1,:)

              call para3d(Xstart,lenx,leny,lenz,3,SS)

              ! Now get the deltas from the FE solution
              call hexa_index(i-1,j-1,k-1,indices)
              do ii=1,8
                 do iii=1,3
                    irow = l_index(lptr(nbkGlobal) +&
                         indices(ii,1)*l_sizes(nbkGlobal,2)*l_sizes(nbkGlobal,3)*3 + &
                         indices(ii,2)*l_sizes(nbkGlobal,3)*3 +&
                         indices(ii,3)*3+ iii )
                    
                    if (irow .lt. nuu) then
                       call VecGetValues(uu,1,irow,deltas(ii,iii),ierr)
                    else
                       deltas(ii,iii) = allBCval(counter*24+3*(ii-1)+iii)
                    end if
                 end do
              end do
              ! Now actually update the Mesh Variables 'X'
              do ii=1,lenx
                 do jj=1,leny
                    do kk=1,lenz

                       nnr(1) = 1.0 - SS(ii,jj,kk,1)
                       nnr(2) = SS(ii,jj,kk,1)

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

                       XSW(indx+ii-1,indy+jj-1,indz+kk-1,:) =  &
                            Xinit(indx+ii-1,indy+jj-1,indz+kk-1,:) + pt_delta
                       X  (indx+ii-1,indy+jj-1,indz+kk-1,:) =  &
                            Xinit(indx+ii-1,indy+jj-1,indz+kk-1,:) + pt_delta

                    end do ! kk loop
                 end do ! jj loop
              end do !ii loop
              counter = counter + 1 
              deallocate(Xstart,SS)
           end do ! k loop
        end do ! j loop
     end do ! i loop
  end do

  ! Re-update the global faces

  call updateFacesGlobal(mdNSurfNodesCompact,mdGlobalSurfxx,.false.)
  
  do nn=1,nDom
     call setPointers(nn,level,sps)     
     !call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)

     ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
     ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
     IMAX = IL
     JMAX = JL
     KMAX = KL

     ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
     xyz0 = 0
     xyznew = 0
     
     XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XSW(1:IMAX,1:JMAX,1:KMAX,1)
     XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XSW(1:IMAX,1:JMAX,1:KMAX,2)
     XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XSW(1:IMAX,1:JMAX,1:KMAX,3)
     XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
     XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
     XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)

     call warp_local(xyznew,xyz0,imax,jmax,kmax)
     
     ! ASSIGN THESE NEW XYZ VALUES TO THE MESH ITSELF
     DO I=1,IMAX
        DO J=1,JMAX
           DO K=1,KMAX
              X(I,J,K,1) = XYZNEW(1,I,J,K)
              X(I,J,K,2) = XYZNEW(2,I,J,K)
              X(I,J,K,3) = XYZNEW(3,I,J,K)
           END DO
        END DO
     END DO
     deALLOCATE(XYZ0,XYZNEW)
   end do


call destroyWarpMeshSolid()
#endif
end subroutine warpMeshSolid


subroutine initializeWarpMeshSolid(nuu,nus,l_index,lptr,l_sizes,nli,nblock)
  ! This subroutine, sets up the preprocessing information required
  ! for warpMeshSolid
  use blockpointers
  use communication, only: myID,sumb_comm_world,sumb_comm_self,nProc
  use solidwarpmodule
  implicit none

  ! Input Data
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: nli,nblock,nuu,nus


  ! Working Data
  integer(kind=intType) :: i,nn

 ! Petsc non-zero sizes
  integer(kind=intType) :: nonz
  integer(kind=intType),dimension(:), allocatable :: nnz
  integer(kind=intType) :: nelem_local,blockID_local(nDom)
  integer(kind=intType) :: level=1,sps=1,ierr

  ! Allocate some size arrays
  allocate(allNDom(nProc),cumNDom(nProc+1),allNElem(nProc),cumNElem(nProc+1))
  allocate(allBlockIDs(nBlock))
  allocate(Kdispls(nProc),BCdispls(nProc),Krecvcount(nProc),BCrecvcount(nProc))


  ! Distribute the number of domains on each processor
  call mpi_allgather(nDom,1,sumb_integer,allnDom,1,sumb_integer,PETSC_COMM_WORLD,ierr)
  cumNDom(1) = 0
  do i=1,nProc
     cumNDom(i+1) = cumNDom(i) + allNDom(i)
  end do
  ! Compute the number of elements on this processor
  nelem_local = 0
  do nn=1,nDom
     call setPointers(nn,level,sps)
     nelem_local = nelem_local + (l_sizes(nbkGlobal,1)-1)*(l_sizes(nbkGlobal,2)-1)*(l_sizes(nBkGlobal,3)-1)
     blockID_local(nn) = nbkGlobal
  end do
  call mpi_allgather(nelem_local,1,sumb_integer,allNElem,1,sumb_integer,PETSC_COMM_WORLD,ierr)
  cumNElem(1) = 0
  do i=1,nProc
     cumNElem(i+1) = cumNElem(i) + allNElem(i)
  end do
  
  nElem = cumNElem(nProc+1)
  ! Get the global BlockID list
  call mpi_allgatherv(blockID_local,nDom,sumb_integer,allBlockIDs,allNDom,cumNDom(1:nProc),sumb_integer,SUmb_comm_world,ierr)

  ! Now we can allocate Allk,allBc,Klocal ect

  allocate(allK(24*24*nElem),allBCVal(24*nElem))
  allocate(localK(allNelem(myID+1)*24*24),localBCVal(allNelem(myID+1)*24))

  ! Set up K and BC sending data
  Kdispls = cumNElem(1:nProc)*24*24
  BCdispls = cumNElem(1:nProc)*24
  Krecvcount = allNElem(:)*24*24
  BCrecvcount = allNelem(:)*24
  ! --------- Setup PETSc Arrays --------------
  PETScNegOne = -1.0
  PETScOne = 1
  ! Kuu
  allocate(nnz(nuu))
  nnz(:) = min(81,nuu) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(27,nuu)

  call MatCreate(PETSC_COMM_SELF,kuu,ierr)
  call MatSetSizes(kuu,nuu,nuu,nuu,nuu,ierr)
  call MatSetType(kuu,'seqaij',ierr)
  call MatSeqAIJSetPreallocation(kuu,nonz,nnz,ierr)
  call MatSetOption(kuu,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  deallocate(nnz)

  ! Kus 
  allocate(nnz(nuu))
  nnz(:) = min(81,nus) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(81,nus)
  call MatCreate(PETSC_COMM_SELF,Kus,ierr)
  call MatSetSizes(Kus,nuu,nus,nuu,nus,ierr)
  call MatSetType(Kus,'seqaij',ierr)
  call MatSeqAIJSetPreallocation(Kus,nonz,nnz,ierr)
  call MatSetOption(Kus,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  deallocate(nnz)

  ! --------------- uu,us,fu -----------------
  call VecCreate(PETSC_COMM_SELF,uu,ierr)
  call VecSetSizes(uu,nuu,nuu,ierr)
  call VecSetType(uu,'seq',ierr)

  call VecCreate(PETSC_COMM_SELF,us,ierr)
  call VecSetSizes(us,nus,nus,ierr)
  call VecSetType(us,'seq',ierr)

  call VecCreate(PETSC_COMM_SELF,fu,ierr)
  call VecSetSizes(fu,nuu,nuu,ierr)
  call VecSetType(fu,'seq',ierr)
  print *,'Done initialization'
end subroutine initializeWarpMeshSolid

subroutine destroyWarpMeshSolid ! Deallocation
  use solidwarpmodule
  use communication, only: sumb_comm_world
  implicit none
  integer(kind=intType) :: ierr
  deallocate(allK,allBCVal,allBlockIDs,localK,localBCVal)
  deallocate(allNDom,cumNDom,allNElem,cumNElem)
  deallocate(Kdispls,BCdispls,Krecvcount,BCrecvcount)

  call MatDestroy(Kuu,ierr)
  call MatDestroy(Kus,ierr)
  call VecDestroy(uu,ierr)
  call VecDestroy(us,ierr)
  call VecDestroy(fu,ierr)
  call mpi_barrier(sumb_comm_world, ierr)
end subroutine destroyWarpMeshSolid

subroutine calcstiffness(points,Kstif)
  use precision 
  implicit none

  ! Input/Output
  real(kind=realType) :: points(8,3)
  real(kind=realType) :: Kstif(24*24)

  ! Working 
  real(kind=realType) :: c,G,E,volume,nu,onemnuc,nuc
  real(kind=realType) :: Jac(3,3),Jinv(3,3)
  integer(kind=intType):: i,j,k,l,m,n,ii,jj,kk
  real(kind=realType) :: r,s,t
  real(kind=realType) :: Nr(8),Ns(8),Nt(8),Nx(8),Ny(8),Nz(8)
  real(kind=realType) :: nnr(2),nns(2),nnt(2),ndr(2),nds(2),ndt(2)
  real(kind=realType) :: g_points(2),det,invdet

  ! First calculate the volume needed for the stiffness

  call volume_hexa(points,volume)

  E = 1.0/(volume)!*volume)!*volume)
  nu = -.2
  c = E/((1+nu)*(1-2*nu))
  G = c/(2-nu)

  onemnuc = (1.0-nu)*c
  nuc = nu*c

  Kstif(:) = 0.0

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

                 !                  intM = matmul(transpose(Bi), matmul(Cm,Bj)) ==->

                 ![ Nxi*onemnuc*Nxj+Nyi*G*Nyj+Nzi*G*Nzj,               Nxi*nuc*Nyj+Nyi*G*Nxj,               Nxi*nuc*Nzj+Nzi*G*Nxj]
                 ![               Nyi*nuc*Nxj+Nxi*G*Nyj, Nyi*onemnuc*Nyj+Nxi*G*Nxj+Nzi*G*Nzj,               Nyi*nuc*Nzj+Nzi*G*Nyj]
                 ![               Nzi*nuc*Nxj+Nxi*G*Nzj,               Nzi*nuc*Nyj+Nyi*G*Nzj, Nzi*onemnuc*Nzj+Nyi*G*Nyj+Nxi*G*Nxj]

                 ! USE ROW ORIENTED 1D array for Kstif
                 Kstif((3*i-3)*24  + 3*j-2) =Kstif((3*i-3)*24  + 3*j-2) + &
                      (Nx(i)*onemnuc*Nx(j)+Ny(i)*G*Ny(j)+Nz(i)*G*Nz(j))*det
                 Kstif((3*i-3)*24  + 3*j-1) = Kstif((3*i-3)*24  + 3*j-1) + &
                      (Nx(i)*nuc*Ny(j)+Ny(i)*G*Nx(j))*det
                 Kstif((3*i-3)*24  + 3*j  ) =Kstif((3*i-3)*24  + 3*j)  + &
                      (Nx(i)*nuc*Nz(j)+Nz(i)*G*Nx(j))*det

                 Kstif((3*i-2)*24 + 3*j-2) = Kstif((3*i-2)*24 + 3*j-2) + & 
                      (Ny(i)*nuc*Nx(j)+Nx(i)*G*Ny(j))*det
                 Kstif((3*i-2)*24 + 3*j-1) = Kstif((3*i-2)*24 + 3*j-1) + &
                      (Ny(i)*onemnuc*Ny(j)+Nx(i)*G*Nx(j)+Nz(i)*G*Nz(j))*det
                 Kstif((3*i-2)*24 + 3*j  ) = Kstif((3*i-2)*24 + 3*j  ) + & 
                      (Ny(i)*nuc*Nz(j)+Nz(i)*G*Ny(j))*det

                 Kstif((3*i-1)*24  +3*j-2) = Kstif((3*i-1)*24  +3*j-2) + &
                      (Nz(i)*nuc*Nx(j)+Nx(i)*G*Nz(j))*det
                 Kstif((3*i-1)*24  +3*j-1) = Kstif((3*i-1)*24  +3*j-1) + &
                      (Nz(i)*nuc*Ny(j)+Ny(i)*G*Nz(j))*det
                 Kstif((3*i-1)*24  +3*j  ) = Kstif((3*i-1)*24  +3*j  ) + &
                      (Nz(i)*onemnuc*Nz(j)+Ny(i)*G*Ny(j)+Nx(i)*G*Nx(j))*det
              end do
           end do
        end do
     end do
  end do
  ! Fill up the remainder
  do i=1,24
     do j = 1,i
        Kstif((24*(i-1))+j) = Kstif((j-1)*24 + i)
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

! !    ! Ksp 
!      call KSPCreate(PETSC_COMM_SELF,ksp, ierr)
!      call KSPSetOperators(ksp,kuu,kuu,DIFFERENT_NONZERO_PATTERN,ierr)
!      call KSPSetFromOptions(ksp, ierr)
!      call KSPSetType(ksp, "gmres", ierr) !preonly
!      call KSPSetPreconditionerSide(ksp, PC_LEFT, ierr)
!      reltol = 1e-12
!      abstol = 1.0e-16
!      divtol = 1.0e6
!      iterations = 10000
!      call KSPSetTolerances(ksp,reltol,abstol,divtol,iterations,ierr)
!      call KSPGetPC(ksp, pc,ierr)
!      call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
!           PETSC_NULL_FUNCTION, ierr)
!      call PCSetType(pc,"ilu",ierr)
!      call PCFactorSetLevels(pc, 0, ierr)
!      call PCFactorSetMatOrderingType(pc,'rcm',ierr)
!      call KspView(ksp,PETSC_VIEWER_STDOUT_SELF,ierr)
!      call KSPSolve(ksp,fu,uu,ierr)

!   call MatSetType(kuu,'seqsbaij',ierr)
!   call MatSeqSBAIJSetPreallocation(kuu,1,nonz,nnz,ierr)
!   call MatSetOption(kuu,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE,ierr)

!   print *,'cholesky factor'
!   call ISCreateStride(PETSC_COMM_SELF,nuu,0,1,row_perm,ierr)
!   call MatCholeskyFactor(kuu,row_perm,info,ierr)

!  print *,'solve'

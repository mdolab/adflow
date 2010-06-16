
! ***********************************
! *  File: warpMeshSolid.F90
! *  Author: Gaetan Kenway
! *  Started: 03-23-2010
! *  Modified: 02-23-2010
! ***********************************
subroutine warpMeshSolid(nuu,nus,l_index,lptr,l_sizes,nli,nblock,ifaceptb,iedgeptb)
#ifndef USE_NO_PETSC

  use blockpointers
  use communication, only: myID,sumb_comm_world,sumb_comm_self,nProc
  use solidwarpmodule
  use mdData    

  implicit none
  ! Input Data
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: ifaceptb(nblock,6),iedgeptb(nblock,12)
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
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,paraS
  real(kind=realType)   :: reltol,abstol,divtol
  integer(kind=intType) :: max_iter
  !real(kind=realType),dimension(:), allocatable   :: uu_local,us_local,values


  ! Temporary Mesh Arrays
  real(kind=realType), dimension(:,:,:,:), allocatable :: SS
  real(kind=realType) ::  shp(8),pt_delta(3),nns(2),nnt(2),nnr(2),value

  logical :: print_it
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

              ! This is NOT indx+1, its the next mesh index for the next node
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

  if (.not. use_parallel) then
     call mpi_allgatherv(localK,allNelem(myID+1)*24*24,sumb_real,allK,Krecvcount,Kdispls,sumb_real,SUmb_comm_world,ierr)
     call mpi_allgatherv(localBCVal,allNelem(myID+1)*24,sumb_real,allBCVal,BCrecvcount,BCdispls,sumb_real,SUmb_comm_world,ierr)
     call checkError(ierr,'warpMeshSolid','MPI allgatherv error')
  end if
  ! ------------------ Assemble Global Stiffness Matrix on each Processor ----------------

  if (.not. use_parallel) then 
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
  else ! We are using parallel assembly
     counter = 0
     do nn=1,nDom
        call setPointers(nn,1,1)
        blockID = nbkGlobal
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
                          call VecSetValues(us,1,irow-nuu,localBCVal(counter*24+3*(ii-1)+iii),&
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
                                        localK(counter*24*24+(3*(ii-1)+iii-1)*24+3*(jj-1)+jjj), &
                                        ADD_VALUES,ierr)
                                else
                                   call MatSetValues(Kus,1,irow,1,jcol-nuu, &
                                        localK(counter*24*24+(3*(ii-1)+iii-1)*24+3*(jj-1)+jjj), &
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
  end if

  call checkError(ierr,'warpMeshSolid','Error in Matrix assmebly')

  ! ------------------- Do petsc assmebly and solve -------------
  call MatAssemblyBegin(Kuu,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Kuu,MAT_FINAL_ASSEMBLY,ierr)
  call checkError(ierr,'warpMeshSolid','Error in kuu assembly')
  call MatAssemblyBegin(Kus,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Kus,MAT_FINAL_ASSEMBLY,ierr)
  call checkError(ierr,'warpMeshSolid','Error in kus assembly')
  call VecAssemblyBegin(us,ierr)
  call VecAssemblyEnd(us,ierr)
  call checkError(ierr,'warpMeshSolid','Error in us assembly')
  ! Multiply kus by us to get -Fu
  call MatMult(Kus,us,fu,ierr)
  call checkError(ierr,'warpMeshSolid','Kus*us error')
  call VecScale(fu,PETScNegOne,ierr)

  ! Ksp 
  if (use_parallel) then
!      call KSPCreate(PETSC_COMM_WORLD,ksp, ierr)
!      call KSPSetOperators(ksp,kuu,kuu,DIFFERENT_NONZERO_PATTERN,ierr)
!      call KSPSetFromOptions(ksp, ierr)
!      call KSPSetType(ksp, "bicg", ierr) !preonly
!      reltol = 1e-10
!      abstol = 1e-16
!      divtol = 1e3
!      max_iter = 10000
!      call KSPSetTolerances(ksp,reltol,abstol,divtol,max_iter,ierr)
!      call KSPGetPC(ksp, pc,ierr)
!      call PCSetType( pc, "asm",ierr)
!      call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
!           PETSC_NULL_FUNCTION, ierr)
!      call KSPSolve(ksp,fu,uu,ierr)

!      Solve direct with spooles
     call KSPCreate(PETSC_COMM_WORLD,ksp, ierr)
     call KSPSetOperators(ksp,kuu,kuu,DIFFERENT_NONZERO_PATTERN,ierr)
     call KSPSetFromOptions(ksp, ierr)
     call KSPSetType(ksp, "preonly", ierr) !preonly
     reltol = 1e-10
     abstol = 1e-16
     divtol = 1e3
     max_iter = 10000
     call KSPSetTolerances(ksp,reltol,abstol,divtol,max_iter,ierr)
     call KSPGetPC(ksp, pc,ierr)
     call PCSetType( pc, "lu",ierr)
     call PCFactorSetMatSolverPackage(pc,"spooles",ierr)
     call PCFactorSetMatOrderingType(pc, "nd",ierr)
     call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
          PETSC_NULL_FUNCTION, ierr)
     call KSPSolve(ksp,fu,uu,ierr)
  else 
     call MatGetOrdering(kuu,"nd",row_perm,col_perm,ierr)
     call checkError(ierr,'warpMeshSolid','kuu matordering error')
     call MatLUFactor(kuu,row_perm,col_perm,info,ierr)  
     call checkError(ierr,'warpMeshSolid','MatLUfactor error')
     call MatSolve(kuu,fu,uu,ierr)
     call checkError(ierr,'warpMeshSolid','MatSolve error')
  end if

  ! ------------------ Set the Volume Coordinates from Solution -----------------

  if (.not. use_parallel) then
     counter = cumNelem(myID+1)
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

                 allocate(SS(lenx,leny,lenz,3))

                 call para3d(Xinit(indx:indxp1,indy:indyp1,indz:indzp1,:),lenx,leny,lenz,3,SS)

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
                 call checkError(ierr,'warpMeshSolid','VecGetValues error')
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
                 deallocate(SS)
              end do ! k loop
           end do ! j loop
        end do ! i loop
     end do
  else
    
     ! Scatter the vector to all the procs
     call VecScatterCreateToAll(uu,uu_scatter,uu_local,ierr)
     call VecScatterBegin(uu_scatter,uu,uu_local,ADD_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (uu_scatter,uu,uu_local,ADD_VALUES,SCATTER_FORWARD,ierr)

     call VecScatterCreateToAll(us,us_scatter,us_local,ierr)
     call VecScatterBegin(us_scatter,us,us_local,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (us_scatter,us,us_local,INSERT_VALUES,SCATTER_FORWARD,ierr)

  call updateFacesGlobal(mdNSurfNodesCompact,mdGlobalSurfxx,.false.)
  
  do nn=1,nDom
     call setPointers(nn,level,sps)     
     !call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)

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

                 allocate(SS(lenx,leny,lenz,3))

                 call para3d(Xinit(indx:indxp1,indy:indyp1,indz:indzp1,:),lenx,leny,lenz,3,SS)

                 ! Now get the deltas from the FE solution
                 call hexa_index(i-1,j-1,k-1,indices)
                 do ii=1,8
                    do iii=1,3
                       irow = l_index(lptr(nbkGlobal) +&
                            indices(ii,1)*l_sizes(nbkGlobal,2)*l_sizes(nbkGlobal,3)*3 + &
                            indices(ii,2)*l_sizes(nbkGlobal,3)*3 +&
                            indices(ii,3)*3+ iii )

                       if (irow .lt. nuu) then
                          call VecGetValues(uu_local,1,irow,deltas(ii,iii),ierr)
                       else
                          call VecGetValues(us_local,1,irow-nuu,deltas(ii,iii),ierr)
                       end if
                    end do
                 end do
                 !call checkError(ierr,'warpMeshSolid','VecGetValues error')
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
                 deallocate(SS)
              end do ! k loop
           end do ! j loop
        end do ! i loop
     end do
     call VecScatterDestroy(uu_scatter,ierr)
     call VecScatterDestroy(us_scatter,ierr)
     call VecDestroy(us_local,ierr)
     call VecDestroy(uu_local,ierr)
     
  end if
  ! Re-update the global faces
  call updateFacesGlobal(mdNSurfNodesCompact,mdGlobalSurfxx,.false.)

  do nn=1,nDom
     call setPointers(nn,level,sps)   

     ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
     ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK

     ALLOCATE(XYZ0(3,0:IL+1,0:JL+1,0:KL+1),XYZNEW(3,0:IL+1,0:JL+1,0:KL+1),&
          paraS(3,0:IL+1,0:JL+1,0:KL+1))

     XYZ0(1,1:IL,1:JL,1:KL) = XSW(1:IL,1:JL,1:KL,1)
     XYZ0(2,1:IL,1:JL,1:KL) = XSW(1:IL,1:JL,1:KL,2)
     XYZ0(3,1:IL,1:JL,1:KL) = XSW(1:IL,1:JL,1:KL,3)

     XYZNEW(1,1:IL,1:JL,1:KL) = X(1:IL,1:JL,1:KL,1)
     XYZNEW(2,1:IL,1:JL,1:KL) = X(1:IL,1:JL,1:KL,2)
     XYZNEW(3,1:IL,1:JL,1:KL) = X(1:IL,1:JL,1:KL,3)

     ! Only warp if we have at least one implicitly perturbed face or edge
     IF (MAXVAL(IFACEPTB(nbkGlobal,:)) >= 1 .OR. MAXVAL(IEDGEPTB(nbkGlobal,:)) >= 1) THEN
        CALL WARPBLK_solid(IFACEPTB(nbkGlobal,:),IEDGEPTB(nbkGlobal,:),-3,0,1,IL,JL,KL,&
             XYZ0,paraS,XYZnew)
     END IF

     ! ASSIGN THESE NEW XYZ VALUES TO THE MESH ITSELF
     DO I=1,IL
        DO J=1,JL
           DO K=1,KL
              X(I,J,K,1) = XYZNEW(1,I,J,K)
              X(I,J,K,2) = XYZNEW(2,I,J,K)
              X(I,J,K,3) = XYZNEW(3,I,J,K)
           END DO
        END DO
     END DO
     deALLOCATE(XYZ0,XYZNEW,paraS)
  end do

  call mpi_barrier(sumb_comm_world, ierr)
  !call calculateSolidWarpDeriv(nuu,nus)
  call destroyWarpMeshSolid()
  call mpi_barrier(sumb_comm_world, ierr)
#endif
end subroutine warpMeshSolid

subroutine calculateSolidWarpDeriv(nuu,nus)

  use precision 
  use communication, only: myID,sumb_comm_world
  use blockpointers
  use solidwarpmodule
  !use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use mdData          !mdNSurfNodes,mdNSurfNodesCompact

  implicit none

  ! After we have the back solve solution to the FE stiffness problem,
  ! we can compute the (dense) back-solve solution for each DOF
  integer(kind=intType) :: nuu,nus,ierr
  integer(kind=intType) :: n,nn,m,mm,i,j,k,idim
  integer(kind=intType) :: level=1,sps=1
  real(kind=realType)    :: values(nuu),val
  integer(kind=intType) :: indices(nuu)
  integer(kind=intType) :: nDimX

  ! This technically performs MORE backs solves then strictly
  ! necessary, since it includes the constrained DOF on the
  ! boundary. However, since this is EMBARASSINGLY parrallel, it makes
  ! little difference

  call MatCreate(sumb_comm_world,dXvFEdXsFE,ierr)
  call MatSetSizes(dXvFEdXsFE,PETSC_DECIDE,PETSC_DECIDE,nuu,nus,ierr) ! mat,m,n,M,n,ierr
  call MatSetType(dXvFEdXsFE,"mpidense",ierr)
  call MatSetOption(dXvFEdXsFE,MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
  call MatGetOwnershipRangeColumn(dXvFEdXsFE,m,n,ierr)
  call checkError(ierr,'warpMeshSolid','Error creating dXvFEdXsFE')

  do i=1,nuu
     indices(i) = i-1
  end do
  print *,'before back solve'

  do i=m,n-1
     call VecZeroEntries(us,ierr)
     call VecSetValues(us,1,i,one,INSERT_VALUES,ierr)
     call VecAssemblyBegin(us,ierr)
     call VecAssemblyEnd(us,ierr)
     call MatMult(Kus,us,fu,ierr)
     call VecScale(fu,PETScNegOne,ierr)
     call MatSolve(kuu,fu,uu,ierr)
     call VecGetValues(uu,nuu,indices,values,ierr)
     call MatSetValues(dXvFEdXsFE,nuu,indices,1,i,values,INSERT_VALUES,ierr)
  end do

  call MatAssemblyBegin(dXvFEdXsFE,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(dXvFEdXsFE,MAT_FINAL_ASSEMBLY,ierr)

  do nn = 1,nDom
     call setPointers(nn,level,sps)
     call setPointersAdj(nn,level,sps)

     do mm = 1,mdNGlobalSurfNodesLocal(myID+1)
        ! Is this coordinate on this block?
        if( mdSurfGlobalIndLocal(4,mm)==nn)then
           ! We need to determine IF this node is on the structural FE grid
           i = mdSurfGlobalIndLocal(1,mm)
           j = mdSurfGlobalIndLocal(2,mm)
           k = mdSurfGlobalIndLocal(3,mm)
           do idim=1,3
              if (FE_id(i,j,k,idim) .ne. -1) then
                 !                call VecZeroEntries(us,ierr)
                 !                  call VecSetValues(us,1,FE_ID(i,j,k,idim),one,INSERT_VALUES,ierr)
                 !                  call VecAssemblyBegin(us,ierr)
                 !                  call VecAssemblyEnd(us,ierr)
                 !                  call MatMult(Kus,us,fu,ierr)
                 !                  call VecScale(fu,PETScNegOne,ierr)
                 !                  call MatSolve(kuu,fu,uu,ierr)
                 !                  call VecGetValues(uu,nuu,indices,values,ierr)

                 !                 print *,'nn,mm,i,j,k:',nn,mm,i,j,k,FE_ID(i,j,k,:)
              end if
           end do
        end if
     end do
  end do

end subroutine calculateSolidWarpDeriv


subroutine initializeWarpMeshSolid(nuu,nus,l_index,lptr,l_sizes,nli,nblock)
  ! This subroutine, sets up the preprocessing information required
  ! for warpMeshSolid
  use blockpointers
  !use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication, only: myID,sumb_comm_world,sumb_comm_self,nProc
  use solidwarpmodule
  use mdData          !mdNSurfNodes,mdNSurfNodesCompact

  implicit none

  ! Input Data
  integer(kind=intType) :: l_index(nli),lptr(nblock+1),l_sizes(nblock,3)
  integer(kind=intType) :: nli,nblock,nuu,nus

  ! Working Data

  character             :: mat_type*6,vec_type*3
  integer               :: comm_to_use
  integer(kind=intType) :: n,nn,i,ii,iii,j,jj,jjj,k,kk,kkk
  integer(kind=intType) :: istart,jstart,kstart
  integer(kind=intType) :: indx,indy,indz,indxp1,indyp1,indzp1
  integer(kind=intType) :: lenx,leny,lenz
  integer(kind=intType) :: nelemx,nelemy,nelemz
  integer(kind=intType) :: col_indices(8,3),indices(8,3)
  integer(kind=intType) :: nDimX,idxvol
  real(kind=realType) ::  shp(8),pt_delta(3),nns(2),nnt(2),nnr(2)
  real(kind=realType), dimension(:,:,:,:), allocatable :: SS
  integer(kind=intType) :: counter
  ! Petsc non-zero sizes
  integer(kind=intType) :: nonz
  integer(kind=intType),dimension(:), allocatable :: nnz
  integer(kind=intType) :: nelem_local,blockID_local(nDom)
  integer(kind=intType) :: level=1,sps=1,ierr

  use_parallel = .True.
  if (.not. use_parallel) then

     ! Allocate some size arrays
     allocate(allNDom(nProc),cumNDom(nProc+1),allNElem(nProc),cumNElem(nProc+1))
     allocate(allBlockIDs(nBlock))
     allocate(Kdispls(nProc),BCdispls(nProc),Krecvcount(nProc),BCrecvcount(nProc))

     ! Distribute the number of domains on each processor
     call mpi_allgather(nDom,1,sumb_integer,allnDom,1,sumb_integer,PETSC_COMM_WORLD,ierr)
     call checkError(ierr,'warpMeshSolid','Error cathering allnDom')

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
     call checkError(ierr,'warpMeshSolid','Error cathering allNElem')

     cumNElem(1) = 0
     do i=1,nProc
        cumNElem(i+1) = cumNElem(i) + allNElem(i)
     end do

     nElem = cumNElem(nProc+1)
     ! Get the global BlockID list
     call mpi_allgatherv(blockID_local,nDom,sumb_integer,allBlockIDs,allNDom,cumNDom(1:nProc),sumb_integer,SUmb_comm_world,ierr)
     call checkError(ierr,'warpMeshSolid','Error cathering allBlockIDs')

     ! Now we can allocate Allk,allBc,Klocal ect

     allocate(allK(24*24*nElem),allBCVal(24*nElem))
     allocate(localK(allNelem(myID+1)*24*24),localBCVal(allNelem(myID+1)*24))
     ! Set up K and BC sending data
     Kdispls = cumNElem(1:nProc)*24*24
     BCdispls = cumNElem(1:nProc)*24
     Krecvcount = allNElem(:)*24*24
     BCrecvcount = allNelem(:)*24

  else
     nelem_local = 0
     do nn=1,nDom
        call setPointers(nn,level,sps)
        nelem_local = nelem_local + (l_sizes(nbkGlobal,1)-1)*(l_sizes(nbkGlobal,2)-1)*(l_sizes(nBkGlobal,3)-1)
     end do
     allocate(localK(nelem_local*24*24),localBCVal(nelem_local*24)) ! This is all we need for parallel case
  end if

  ! --------- Setup PETSc Arrays --------------
  PETScNegOne = -1.0
  PETScOne = 1

  if (use_parallel) then
     mat_type = 'mpiaij'
     vec_type = 'mpi'
     comm_to_use = PETSC_COMM_WORLD
  else
     mat_type = 'seqaij'
     vec_type = 'seq'
     comm_to_use = PETSC_COMM_SELF
  end if


  ! Kuu
  allocate(nnz(nuu))
  nnz(:) = min(81,nuu) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(81,nuu)

  

  call MatCreate(comm_to_use,kuu,ierr)
  if (use_parallel) then
     call MatSetSizes(kuu,PETSC_DECIDE,PETSC_DECIDE,nuu,nuu,ierr)
  else
     call MatSetSizes(kuu,nuu,nuu,nuu,nuu,ierr)
  end if
  call MatSetType(kuu,mat_type,ierr)
  if (use_parallel) then
     call MatMPIAIJSetPreallocation(kuu,nonz,nnz,nonz,nnz,ierr)
  else
     call MatSeqAIJSetPreallocation(kuu,nonz,nnz,ierr)
  end if

  call MatSetOption(kuu,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  call checkError(ierr,'warpMeshSolid','Error creating kuu')
  deallocate(nnz)

  ! Kus 
  allocate(nnz(nuu))
  nnz(:) = min(81,nus) ! Petsc is screwed up...we can't just pass in the singe,
  ! nz value we MUST pass in the full nnz array
  nonz = min(81,nus)
  call MatCreate(comm_to_use,Kus,ierr)
  if (use_parallel) then
     call MatSetSizes(Kus,PETSC_DECIDE,PETSC_DECIDE,nuu,nus,ierr)
  else
     call MatSetSizes(Kus,nuu,nus,nuu,nus,ierr)
  end if
  call MatSetType(Kus,mat_type,ierr)
  if (use_parallel) then
     call MatMPIAIJSetPreallocation(kus,nonz,nnz,nonz,nnz,ierr)
  else
     call MatSeqAIJSetPreallocation(kus,nonz,nnz,ierr)
  end if
  call MatSetOption(Kus,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  call checkError(ierr,'warpMeshSolid','Error creating kus')
  deallocate(nnz)

  ! --------------- uu,us,fu -----------------
  call VecCreate(comm_to_use,uu,ierr)
  call VecSetType(uu,vec_type,ierr)

  call VecCreate(comm_to_use,us,ierr)
  call VecSetType(us,vec_type,ierr)

  call VecCreate(comm_to_use,fu,ierr)
  call VecSetType(fu,vec_type,ierr)

  if (use_parallel) then
     call VecSetSizes(uu,PETSC_DECIDE,nuu,ierr)
     call VecSetSizes(us,PETSC_DECIDE,nus,ierr)
     call VecSetSizes(fu,PETSC_DECIDE,nuu,ierr)
  else
     call VecSetSizes(uu,nuu,nuu,ierr)
     call VecSetSizes(us,nus,nus,ierr)
     call VecSetSizes(fu,nuu,nuu,ierr)
  end if
  call checkError(ierr,'warpMeshSolid','Error creating uu,us,or fu')

!   if (use_parallel) then
!      ! Create the data we need to be able to scatter uu and us to all
!      ! processors EVEN for the parallel version
!      allocate(uu_sizes(nProc),cum_uu(nProc+1))
!      call VecGetOwnershipRanges(uu,cum_uu,ierr)
!      do i=1,nProc
!         uu_sizes(i) = cum_uu(i+1)-cum_uu(i)
!      end do
!      print *,'myid,cum_uu:',myid,cum_uu
!      allocate(us_sizes(nProc),cum_us(nProc+1))
!      call VecGetOwnershipRanges(us,cum_us,ierr)
!      do i=1,nProc
 !         us_sizes(i) = cum_us(i+1)-cum_us(i)
!      end do
!      print *,'myid,cum_us:',myid,cum_us
!   end if

  !  nDimX = 3 * nNodesLocal*1
  !   allocate(nnz(nDimX))
  !   nnz(:) = 8
  !   nonz = 1 ! Ignored anyway

  !   call MatCreate(sumb_comm_world,dXvdXvFE, ierr)
  !   call MatSetSizes(dXvdXvFE,nDimX,PETSC_DECIDE,PETSC_DETERMINE,nuu+nus,ierr)
  !   call MatSetType(dXvdXvFE,"mpiaij",ierr)
  !   call MatMPIAIJSetPreallocation(dXvdXvFE,nonz,nnz,nonz,nnz,ierr)
  !   call MatSetOption(dXvdXvFE,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
  !   call checkError(ierr,'warpMeshSolid','Error creating dXvdXvFE')

  !   deallocate(nnz)
  !   counter= 1
  !   call MatGetSize(dXvdXvFE,i,j,ierr)
  !   print *,'damn size:',i,j
  !   do nn=1,nDom
  !      call setPointers(nn,level,sps)
  !      call setPointersAdj(nn,level,sps)
  !      FE_ID(:,:,:,:) = -1 ! Set them to -1 by default...not on FE_grid
  !      nelemx = l_sizes(nbkGlobal,1)-1
  !      nelemy = l_sizes(nbkGlobal,2)-1
  !      nelemz = l_sizes(nbkGlobal,3)-1
  !      do i=1,nelemx   
  !         do j =1,nelemy
  !            do k=1,nelemz

  !               ! Now we must figure out the indices
  !               indx = int(floor(dble(i-1.0)/(nelemx)*(nx))) + 1
  !               indy = int(floor(dble(j-1.0)/(nelemy)*(ny))) + 1
  !               indz = int(floor(dble(k-1.0)/(nelemz)*(nz))) + 1

  !               indxp1 = int(floor(dble(i-0)/(nelemx)*(nx))) + 1
  !               indyp1 = int(floor(dble(j-0)/(nelemy)*(ny))) + 1
  !               indzp1 = int(floor(dble(k-0)/(nelemz)*(nz))) + 1

  !               lenx = indxp1-indx+1
  !               leny = indyp1-indy+1
  !               lenz = indzp1-indz+1

  !               allocate(SS(lenx,leny,lenz,3))

  !               call para3d(Xinit(indx:indxp1,indy:indyp1,indz:indzp1,:),lenx,leny,lenz,3,SS)
  !               call hexa_index(i-1,j-1,k-1,indices)

  !               do ii=1,8
  !                  do iii =1,3
  !                     col_indices(ii,iii) = l_index(lptr(nbkGlobal) +&
  !                          indices(ii,1)*l_sizes(nbkGlobal,2)*l_sizes(nbkGlobal,3)*3 + &
  !                          indices(ii,2)*l_sizes(nbkGlobal,3)*3 +&
  !                          indices(ii,3)*3+ iii )
  !                  end do
  !               end do

  !               do iii=1,3
  !                  FE_ID(indx  ,indy  ,indz  ,iii) = col_indices(1,iii)
  !                  FE_ID(indxp1,indy  ,indz  ,iii) = col_indices(2,iii)
  !                  FE_ID(indx  ,indyp1,indz  ,iii) = col_indices(3,iii)
  !                  FE_ID(indxp1,indyp1,indz  ,iii) = col_indices(4,iii)
  !                  FE_ID(indx  ,indy  ,indzp1,iii) = col_indices(5,iii)
  !                  FE_ID(indxp1,indy  ,indzp1,iii) = col_indices(6,iii)
  !                  FE_ID(indx  ,indyp1,indzp1,iii) = col_indices(7,iii)
  !                  FE_ID(indxp1,indyp1,indzp1,iii) = col_indices(8,iii)
  !               end do


  !               if (i == 1) then
  !                  istart = 1
  !               else
  !                  istart = 2
  !               end if

  !               if (j == 1) then
  !                  jstart = 1
  !               else
  !                  jstart = 2
  !               end if
  !               if (k == 1) then
  !                  kstart = 1
  !               else
  !                  kstart = 2
  !               end if

  !               do ii=istart,lenx
  !                  do jj=jstart,leny
  !                     do kk=kstart,lenz

  !                        nnr(1) = 1.0 - SS(ii,jj,kk,1)
  !                        nnr(2) = SS(ii,jj,kk,1)

  !                        nns(1) = 1.0 - SS(ii,jj,kk,2)
  !                        nns(2) = SS(ii,jj,kk,2)

  !                        nnt(1) = 1.0 - SS(ii,jj,kk,3)
  !                        nnt(2) = SS(ii,jj,kk,3)

  !                        shp(1) = nnr(1)*nns(1)*nnt(1)
  !                        shp(2) = nnr(2)*nns(1)*nnt(1)
  !                        shp(3) = nnr(1)*nns(2)*nnt(1)
  !                        shp(4) = nnr(2)*nns(2)*nnt(1)
  !                        shp(5) = nnr(1)*nns(1)*nnt(2)
  !                        shp(6) = nnr(2)*nns(1)*nnt(2)
  !                        shp(7) = nnr(1)*nns(2)*nnt(2)
  !                        shp(8) = nnr(2)*nns(2)*nnt(2)

  !                        do n=1,3
  !                           counter = counter + 1
  !                           idxvol = globalNode(indx+ii-1,indy+jj-1,indz+kk-1)*3+n
  !                           do iii=1,8
  !                              call MatSetValue(dXvdXvFE,idxvol-1,col_indices(iii,n),&
  !                                   shp(iii),INSERT_VALUES,ierr)
  !                           end do
  !                        end do
  !                     end do ! kk loop
  !                  end do ! jj loop
  !               end do !ii loop
  !               deallocate(SS)
  !            end do ! k loop
  !         end do ! j loop
  !      end do ! i loop
  !   end do ! nndomain loop

  !   call checkError(ierr,'warpMeshSolid','Error setting dXvdXvFE values')
  !   call MatAssemblyBegin(dXvdXvFE,MAT_FINAL_ASSEMBLY,ierr)
  !   call MatAssemblyEnd(dXvdXvFE,MAT_FINAL_ASSEMBLY,ierr)
  !   call checkError(ierr,'warpMeshSolid','Error assembling dXvdXvFE')

  print *,'Done initialization'
end subroutine initializeWarpMeshSolid

subroutine destroyWarpMeshSolid ! Deallocation
  use solidwarpmodule
  use communication, only: sumb_comm_world
  implicit none
  integer(kind=intType) :: ierr
  if (.not. use_parallel) then
     deallocate(allK,allBCVal,allBlockIDs,localK,localBCVal)
     deallocate(allNDom,cumNDom,allNElem,cumNElem)
     deallocate(Kdispls,BCdispls,Krecvcount,BCrecvcount)
  else
     deallocate(localK,localBCVal)!,uu_sizes,us_sizes,cum_uu,cum_us)
  end if

  call MatDestroy(Kuu,ierr)
  call checkError(ierr,'warpMeshSolid','Error destroying kuu')

  call MatDestroy(Kus,ierr)
  call checkError(ierr,'warpMeshSolid','Error destroying kus')

  !call MatDestroy(dXvFEdXsFE,ierr)
  !call checkError(ierr,'warpMeshSolid','Error destroying dXvFEdXsFE')

  !call MatDestroy(dXvdXvFE,ierr)
  !call checkError(ierr,'warpMeshSolid','Error destroying dXvdXsFE')

  !call MatDestroy(dXvdXsFE,ierr)
  !call checkError(ierr,'warpMeshSolid','Error destroying dXvdXsFE')

  call VecDestroy(uu,ierr)
  call checkError(ierr,'warpMeshSolid','Error destroying uu')

  call VecDestroy(us,ierr)
  call checkError(ierr,'warpMeshSolid','Error destroying us')

  call VecDestroy(fu,ierr)
  call checkError(ierr,'warpMeshSolid','Error destroying fu')

  call mpi_barrier(sumb_comm_world, ierr)
end subroutine destroyWarpMeshSolid


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



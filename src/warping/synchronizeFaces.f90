
!==========================================================
!
!       File: synchronizeFaces.f90
!       Author: C.A.(Sandy) Mader
!       Start Date: 08-19-2008
!       Modified Date: 12-11-2008
!
!==========================================================

!     ******************************************************************
!     *   SyncrhonizeFaces syncs the face only perturbartions of the   *
!     *   mesh blocks so that the single block warping algorithm can   *
!     *   be used.                                                     *
!     *                                                                *
!     *   08/19/08  C.A.Mader   Converted to *.f90 from python         *
!     *   12/11/08      "       Integrated into the SUmb flow solver   *
!     *                         Now uses the SUmb block.f90 module to  *
!     *                         store data                             *
!     *                                                                *
!     *   C.A.(Sandy) Mader, UTIAS,Toronto,Canada                      *
!     ******************************************************************



subroutine synchronizeBlockFaces(level,sps)
  !""" Synchronize the block faces to ensure that the perturbations are
  !      present everywhere
  !      """
  use blockPointers
  use communication, only: myID,sumb_comm_world
  use BCTypes
  implicit none

  !Subroutine Arguments
  integer(kind=intType),intent(in)::level,sps

  !Local Variables


  logical :: notSynchronized
  integer(kind=inttype) :: count,i,j,k,l,m,n,counter,ierr,ii,jj,kk,nn
  integer(kind=inttype), dimension(3)::step
  integer(kind=inttype) :: counterI,counterJ,counterK,length
  integer(kind=inttype) :: inmin,jnmin,knmin,inmax,jnmax,knmax
  integer(kind=inttype) :: sendtag,recvtag
  real(KIND=realTYPE) :: neighbour,neighbour0,local,local0,tmp,a
  real(KIND=REALTYPE)    :: eps2,realbuf,imagbuf
  integer(kind=inttype) :: destblock,index,neighbourindex
  logical ::test
  complex(KIND=REALTYPE), dimension(3)::coords
  integer(kind=inttype),dimension(3)::indices
  integer(kind=intType),dimension(:),allocatable::incrementI,&
     incrementJ,incrementK,  incrementdI,&
     incrementdJ,incrementdK  
  !complex, allocatable, dimension(:) :: sendBuffer
  !complex, allocatable, dimension(:) :: recvBuffer

  !************************
  ! Begin Execution
  !************************

  !Set the inital state to not syncronized
  
  notSynchronized = .True.
  count = 0
  !print *,'starting loop'
  do while (notSynchronized)
     ! Set state to synchronized and the check for truth at the end of the 
     ! loop.
     notSynchronized = .False.

   !  print *,'looping over blocks'
     !loop over blocks and subfaces
     do i=1,nDom
    !    print *,'setting pointers'
        call setPointers(i,level,sps)

     !   print *,'checkking for allocation of nNodesSubface'
        ! Check to see if the memory is allocated to store the total number
        ! of nodes on each subface. If not, allocate.
        !if(.not. allocated(nNodesSubface)) then
        if(.not. associated(nNodesSubface)) then
           allocate(nNodesSubface(nsubface), stat=ierr)
        endif

        !allocate the memory for the block face communicators
        if (.not. associated(warp_comm)) &
                allocate(warp_comm(nsubface))

        !allocate the memory for the block increments
        allocate(incrementI(nsubface),incrementJ(nsubface),incrementK(nsubface))
        allocate(incrementdI(nsubface),incrementdJ(nsubface),incrementdK(nsubface))
        !print *,'i',i
        call getIncrement(nSubface)
        !print *,'iinter',i
        call getIncrementD(nSubface)
        !print *,'iafter',i,incrementI,'j',incrementJ,'k',incrementK

        !print *,'looping over subfaces'
        do j=1,nSubface
           !print *,'bcfaceID',bcfaceid(j)
           ! determine the number of nodes on this subface
           ! Determine which face is active for this subface
           if (BCFaceID(j)== imin) then !imin
              inmin =  inbeg(j); inmax = inbeg(j)
              jnmin =  jnbeg(j); jnmax = jnend(j)
              knmin =  knbeg(j); knmax = knend(j)
              nNodesSubface(j) =(abs(jnmax-jnmin)+1)*(abs(knmax-knmin)+1)
           elseif (BCFaceID(j)== imax) then!imax
              inmin = inend(j); inmax = inend(j)
              jnmin = jnbeg(j); jnmax = jnend(j)
              knmin = knbeg(j); knmax = knend(j)
              nNodesSubface(j) = (abs(jnmax-jnmin)+1)*(abs(knmax-knmin)+1)
           elseif (BCFaceID(j)== jmin)then!:#jmin
              inmin = inbeg(j); inmax =  inend(j)
              jnmin = jnbeg(j); jnmax =  jnbeg(j)
              knmin = knbeg(j); knmax =  knend(j)
              nNodesSubface(j) = (abs(inmax-inmin)+1)*(abs(knmax-knmin)+1)
           elseif (BCFaceID(j)== jmax)then!:#jmax
              inmin = inbeg(j); inmax = inend(j)
              jnmin = jnend(j); jnmax = jnend(j)
              knmin = knbeg(j); knmax = knend(j)
              nNodesSubface(j) = (abs(inmax-inmin)+1)*(abs(knmax-knmin)+1)
           elseif(BCFaceID(j)== kmin)then!:#kmin
              inmin = inbeg(j); inmax = inend(j)
              jnmin = jnbeg(j); jnmax = jnend(j)
              knmin = knbeg(j); knmax = knbeg(j)
              nNodesSubface(j) = (abs(jnmax-jnmin)+1)*(abs(inmax-inmin)+1)
           elseif (BCFaceID(j)== kmax)then!:#kmax
              inmin =  inbeg(j); inmax =  inend(j)
              jnmin =  jnbeg(j); jnmax =  jnend(j)
              knmin =  knend(j); knmax =  knend(j)
              nNodesSubface(j) = (abs(jnmax-jnmin)+1)*(abs(inmax-inmin)+1)
           else
              print *,'Error:Not a valid face type',BCFaceID(j)
           endif

           ! Set the length for the communication buffer 3 coords + 3 indices + 1 blocknum
           length = 7*nNodesSubface(j)
           
           !print *,'checking allocation of buffers'
           !allocate the communicators for this subface
           !Generate the receive buffer
           ! Dimension 1-3 are coords., 4-6 are indices, 7 is remote block number
           if (.not. allocated(warp_comm(j)%recvBuffer)) &
                allocate(warp_comm(j)%recvBuffer(length))
           if (.not. allocated(warp_comm(j)%sendBuffer)) &
                allocate(warp_comm(j)%sendBuffer(length))
           !print *,'buffers allocated',shape(meshblocks(i)%comm(j)%sendBuffer),shape(meshblocks(i)%comm(j)%recvBuffer)
           
           !print *,'synchronizing faces'
           
           !         #Check for one to one matching internal faces
           !print *,'ifs',meshblocks(i)%BCType(j),meshblocks(i)%neigh_proc(j),meshblocks(i)%neigh_block(j)
           !stop
           if(BCType(j) == B2BMatch)then
           !if(meshblocks(i)%BCType(j) == -16)then
              !print *,'boundary match',i,j
              !      #Check that face requires interprocessor communication
              if(neighproc(j) /= myid)then
                 !#need communication with blocks on another processor
                           
                 !#post non blocking recieves for face communicator
                
                 !Generate a unique tag for this subface
                 recvtag = 10000*neighproc(j)+100*neighblock(j)+i
                 !recvtag = 10000*meshblocks(i)%neigh_proc(j)+100*meshblocks(i)%neigh_block(j)+i+1
                 !print *,'recvtag',recvtag,myid
                 
                 !call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                 !     myID, SUmb_comm_world, recvRequests(i), ierr)
                 
                 !post the non blocking recv
                 !print *,'posting receive'
                 call mpi_irecv(warp_comm(j)%recvBuffer(1), length, sumb_real,&
                      neighproc(j), recvtag, sumb_comm_world, &
                      warp_comm(j)%recvreq, ierr)
                 !recvreq = MPI.WORLD.Irecv( test, meshblocks(i).neigh_proc(j), recvtag)
                                                  
                 !#set a placement counter for the send buffer
                 counter=1
                 
                 !print *,'i',imin,imax+meshblocks(i)%incrementI(j),meshblocks(i)%incrementI(j)
                 !print *,'j',jmin,jmax+meshblocks(i)%incrementJ(j),meshblocks(i)%incrementJ(j)
                 !print *,'k',kmin,kmax+meshblocks(i)%incrementK(j),meshblocks(i)%incrementK(j)
                 
                 !print *,'looping over local indices'
                 !#Loop over the local indices
                 do l =inmin,inmax,incrementI(j)
                    do m =jnmin,jnmax,incrementJ(j)
                       do n=knmin,knmax,incrementK(j)
!!$ do l =imin,imax+meshblocks(i)%incrementI(j),meshblocks(i)%incrementI(j)
!!$                    do m =jmin,jmax+meshblocks(i)%incrementJ(j),meshblocks(i)%incrementJ(j)
!!$                       do n=kmin,kmax+meshblocks(i)%incrementK(j),meshblocks(i)%incrementK(j)
                          !#Set the counter step based on the coordinate transformation for the face
                          step(1) = (abs(l-inmin))
                          step(2) = (abs(m-jnmin))
                          step(3) = (abs(n-knmin))
                          !print *,'step',step(:),l,m,n
                          !print *,'di',meshblocks(i)%dibeg(j),meshblocks(i)%djbeg(j),meshblocks(i)%dkbeg(j)
                          !print *,'increment',meshblocks(i)%incrementdI(j),meshblocks(i)%incrementdJ(j),meshblocks(i)%incrementdK(j)
                          !print *,'abs(meshblocks(i)%l1(j))-1',abs(meshblocks(i)%l1(j)),abs(meshblocks(i)%l2(j)),abs(meshblocks(i)%l3(j))
                          counterI =dinbeg(j)+step(abs(l1(j)))*incrementdI(j)
                          !print *,'counterI'
                          counterJ =djnbeg(j)+step(abs(l2(j)))*incrementdJ(j)
                          !print *,'counterJ'
                          counterK =dknbeg(j)+step(abs(l3(j)))*incrementdK(j)
                          !print *,'counters',counterI,counterJ,counterK
                          !#set the value in the send buffer
                          !print *,'meshblocks',meshblocks(i)%x(1,l,m,n)
                          warp_comm(j)%sendBuffer(counter)=x(l,m,n,1)
                          !print *,'meshblockx'
                          warp_comm(j)%sendBuffer(counter+1)=x(l,m,n,2)
                          warp_comm(j)%sendBuffer(counter+2)=x(l,m,n,3)
                          !print *,'setting counters'
                          warp_comm(j)%sendBuffer(counter+3)=float(counterI)
                          warp_comm(j)%sendBuffer(counter+4)=float(counterJ)
                          warp_comm(j)%sendBuffer(counter+5)=float(counterK)
                          warp_comm(j)%sendBuffer(counter+6)=float(neighblock(j))!-1
                          !print *,'counter',counter,counter+6,shape(meshblocks(i)%comm(j)%sendBuffer)
                          !print *,'sendbuffer',meshblocks(i)%comm(j)%sendBuffer(counter:counter+7)
                                                               
                          counter=counter+7
                       enddo
                    enddo
                 enddo
                 
                 !#Determine the number of elements in the send buffer
                 !length = len(sbuf(:,0))*len(sbuf(0,:))
                            
                 !#Generate the send buffer
                 !test = MPI.Buffer(sbuf, length, datatype)
                 
                 !#Post the buffer send
                 sendtag = 10000*myid+100*(i)+neighblock(j)
                 !sendtag = 10000*myid+100*(i+1)+meshblocks(i)%neigh_block(j)
                 !print *,'sendtag',sendtag,myid        
                 !sendreq = MPI.WORLD.Isend( test, meshblocks(i).neigh_proc(j),sendtag)
                 call mpi_isend(warp_comm(j)%sendBuffer(1), length,sumb_real,&
                      neighproc(j),sendtag, sumb_comm_world, &
                      warp_comm(j)%sendreq, ierr)
                 
                 !#Save the request for later checking
                 !reqList.append(sendreq)
                        
              else              
                 !print *,'internal communication'
                 !internal communication
                 !communicate the blocks that are on the local processor

                 !print *,'looping over local indices'
                 !#Loop over the local indices
                 !print *,'increments',incrementI(j),incrementJ(j),incrementK(j),'loopzize',inmin,inmax,jnmin,jnmax,knmin,knmax
                 do ii =inmin,inmax,incrementI(j)
                    do jj =jnmin,jnmax,incrementJ(j)
                       do kk=knmin,knmax,incrementK(j)
                 !for ii in xrange(imin,imax+meshblocks(i).incrementI(j),meshblocks(i).incrementI(j)):
                 !               for jj in xrange(jmin,jmax+meshblocks(i).incrementJ(j),meshblocks(i).incrementJ(j)):
                 !                   for kk in xrange(kmin,kmax+meshblocks(i).incrementK(j),meshblocks(i).incrementK(j)):
                          !print *,'setting step',ii,jj,kk
                          
                          !#set the counter step based on the face coordinate transformation
                          step(1) = (abs(ii-inmin))
                          step(2) = (abs(jj-jnmin))
                          step(3) = (abs(kk-knmin))
                          !print *,'step is:',step
                          !step = (abs(ii-inmin),abs(jj-jnmin),abs(kk-knmin))
                          !#compute the neighbouring face counters
                          counterI =dinbeg(j)+step(abs(l1(j)))*incrementdI(j)
                          counterJ =djnbeg(j)+step(abs(l2(j)))*incrementdJ(j)
                          counterK =dknbeg(j)+step(abs(l3(j)))*incrementdK(j)
                          !print *,'indices',ii,jj,kk,counterI,counterJ,counterK
                          do nn=1,3
                             neighbourindex = neighblock(j)
                             !reset pointers and get neighbour data
                             call setpointers(neighbourindex,level,sps)
                             !neighbour = meshblocks(meshblocks(i).neigh_block(j)).x(nn,counterI,counterJ,counterK)
                             neighbour = x(counterI,counterJ,counterK,nn)
                             !neighbour0 = meshblocks(meshblocks(i).neigh_block(j)).xInit(nn,counterI,counterJ,counterK)
                             neighbour0 = xInit(counterI,counterJ,counterK,nn)
                             !Reset pointers to the local block.
                             call setpointers(i,level,sps)
                             local = x(ii,jj,kk,nn)
                             local0 = xInit(ii,jj,kk,nn)
                             eps2 = 1.0e-12
                             !print *,'coordcheck',neighbour0,local0
                             !print *,'neighbours',neighbour,neighbour0,abs(neighbour -neighbour0),max(abs(neighbour0),abs(neighbour),eps2)
                             !print *,'locals',local,local0
                             !check the various options and act accordingly
                             if( abs(neighbour -neighbour0)/max(abs(neighbour0),abs(neighbour),eps2)<1e-12) then
                                !#Then do nothing
                                !#print 'No Change in donor Block'
                                a=1
                             elseif(abs(neighbour -local)/max(abs(local),abs(neighbour),eps2)<1e-12)then
                                !#do nothing
                                a=1
                             elseif( (abs(neighbour-local)/max(abs(local),abs(neighbour),eps2)>1e-12) .and. (abs(local -local0)/max(abs(local0),abs(local),eps2)<1e-12))then
                                call setPointers(neighbourindex,level,sps)
                                tmp = x(counterI,counterJ,counterK,nn)
                                call setPointers(i,level,sps)
                                x(ii,jj,kk,nn)=tmp
                                !meshblocks(i)%x(nn,ii,jj,kk) = meshblocks(meshblocks(i).neigh_block(j)).x(nn,counterI,counterJ,counterK)
                                
                                notSynchronized = .True.
                                
                             elseif( (abs(neighbour -local)/max(abs(local),abs(neighbour),eps2)>1e-12 ) .and. (abs(neighbour -local0)/max(abs(local0),abs(neighbour),eps2)<1e-12 ))then
                                
                                !#do nothing
                                a=1

                             elseif ((abs(neighbour -local)/max(abs(local),abs(neighbour),eps2)>1e-12 ) .and.   (abs(local -local0)/max(abs(local0),abs(local),eps2)>1e-12) .and. (abs(neighbour -local0)/max(abs(local0),abs(neighbour),eps2)>1e-12 )) then

                                !print *,'neighbour',neighbour, local,local0
                                !#average new and current
                                print *,'not matching int!!!', myid!MPI.rank
                                call setPointers(neighbourindex,level,sps)     
                                !tmp = (meshblocks(meshblocks(i).neigh_block(j)).x(nn,counterI,counterJ,counterK) +meshblocks(i).x(nn,ii,jj,kk))/2.0
                                tmp = x(counterI,counterJ,counterK,nn) 
                                call setPointers(i,level,sps)
                                tmp = (tmp +x(ii,jj,kk,nn))/2.0
                                                
                                !                #meshblocks(meshblocks(i).neigh_block(j)-1).x(nn,counterI,counterJ,counterK)=tmp
                                x(ii,jj,kk,nn)=tmp
                                
                                notSynchronized = .True.
                             else !#x_current = x_new, no change
                                x(ii,jj,kk,nn)=x(ii,jj,kk,nn)
                             endif
                          end do
                       end do
                    end do
                 end do
                 
              endif
           else
              !#no communicaton required
              !print *,'no communication required'
              !a=1
           endif
              

        end do
        !allocate the memory for the block increments
        deallocate(incrementI,incrementJ,incrementK)
        deallocate(incrementdI,incrementdJ,incrementdK)
     end do
     
     call mpi_barrier(sumb_comm_world, ierr)
       
     !print *,'waiting for communication to finish...', myid!MPI.rank
     !stop
     !loop over blocks and subfaces
     do i=1,nDom
        call setPointers(i,level,sps)
     !do i=1,nmeshblocks!nBlocksLocal
        do j=1,nSubface 
           !only check for receives on procesors that posted them
           !if(BCType(j) == -16)then
           if(BCType(j) == B2BMatch)then
              !      #Check that face requires interprocessor communication
              if(neighproc(j) /= myid)then
                 !#need communication with blocks on another processor
                 !print *,'i',i,j,meshblocks(i)%comm(j)%recvreq,myid
                 call mpi_wait(warp_comm(j)%recvreq,MPI_STATUS_IGNORE,ierr)
              endif
           endif
        enddo
     enddo
     call mpi_barrier(sumb_comm_world, ierr)
     !stop

!            #Wait for all of the faces to be comunicated
!            MPI.Request.Waitall(reqList)
            

     !Now loop over the recieved buffers and update the local blocks
     !do for buf in recvList:
     !loop over blocks and subfaces
     do l=1,nDom
        call setPointers(l,level,sps)
        !do l=1,nmeshblocks!nBlocksLocal
        do m=1,nSubface
           !         #Check for one to one matching internal faces
           !if(BCType(m) == -16)then
           if(BCType(m) == B2BMatch)then
              !      #Check that face requires interprocessor communication
              if(neighproc(m) /= myid)then
                 !#need communication with blocks on another processor   
                 length = nNodesSubface(m)!shape(meshblocks(l)%comm(m)%recvbuffer(:))
                 do i=1,int(length)
                    index = 1+(i-1)*7
                    destblock = int(warp_comm(m)%recvbuffer(index+6))
                    !print *,'destblock',destblock
                    coords(:) =  warp_comm(m)%recvbuffer(index:index+2)
                    indices(:) = int(warp_comm(m)%recvbuffer(index+3:index+5))
                    !print *,'indices',indices,l,m,i
                    !Set pointers for the given destblock
                    call setPointers(destblock,level,sps)
                    
                    do nn=1,3
                      
                       local = x(indices(1),indices(2),indices(3),nn)
                       local0 = xInit(indices(1),indices(2),indices(3),nn)
                       eps2 = 1.0e-12
                       realbuf = real(coords(nn))
                       !if (abs(local-local0)/max(abs(local0),abs(local),eps2)>1e-12) print*,'local',local,local0,realbuf
                       if ( abs(realbuf -local)/max(abs(local),abs(realbuf),eps2)<1e-12 )then
                          
                          !do nothing
                          !print *,'no perturbation'
                          a=1
                        
                       elseif ((abs(realbuf -local)/max(abs(local),abs(realbuf),eps2)>1e-12 ) .and. (abs(local -local0)/max(abs(local0),abs(local),eps2)<1e-12 )) then
                          !meshblocks(destblock).x(nn,indices(1),indices(2),indices(3))=coords(nn)
                          x(indices(1),indices(2),indices(3),nn)=coords(nn)
                          !print *,'synching'
                          !set synchronization flag
                          notSynchronized = .True.

                       elseif ((abs(realbuf -local)/max(abs(local),abs(realbuf),eps2)>1e-12) .and. (abs(realbuf -local0)/max(abs(local0),abs(realbuf),eps2)<1e-12 ))then
                          !print *,'already moved'
                          !print *,'realbuf',realbuf,local,local0
                          !Do nothing
                          a=1

                       elseif ((abs(realbuf -local)/max(abs(local),abs(realbuf),eps2)>1e-12) .and.  (abs(local-local0)/max(abs(local0),abs(local),eps2)>1e-12) .and. (abs(realbuf -local0)/max(abs(local0),abs(realbuf),eps2)>1e-12 )) then

                          !#Donor and recipient do not match
                          print *,'not matching ext!!!', myid!MPI.rank
                        
                          !#average new and current
                          !       #if MPI.COMM_WORLD.rank ==0:
                          !       #    print "doing averaging"
                          tmp = (x(indices(1),indices(2),indices(3),nn)+coords(nn))/2.0
                                  
                          x(indices(1),indices(2),indices(3),nn)=tmp
                            
                          notSynchronized = .True.
                       else !#x_current = x_new, no change
                          if (myid==0)then !myid
                             print *,"doing else"!,meshblocks(int(buf(9,i))).x(nn,int(buf(6,i)),int(buf(7,i)),int(buf(8,i))),buf(nn,i),buf(nn+3,i)
                            !meshblocks(int(buf(9,i))).x(nn,int(buf(6,i)),int(buf(7,i)),int(buf(8,i)))=meshblocks(int(buf(9,i))).x(nn,int(buf(6,i)),int(buf(7,i)),int(buf(8,i)))
                          endif
                       endif
                    end do
                    !return pointers to normal
                    call setPointers(l,level,sps)
                 end do
              end if
           endif
        end do
     end do
     
     call mpi_barrier(sumb_comm_world, ierr)
     !call mpi_allreduce(notSynchronized,test,1,MPI_INTEGER,MPI_SUM,warp_comm_world,ierr)
     call mpi_allreduce(notSynchronized,test,1,MPI_LOGICAL,MPI_LOR,sumb_comm_world,ierr)
     !test = MPI.WORLD.Allreduce(notSynchronized,MPI.SUM)
     if (myid ==0) then !myid
        print *,'Synchronization test',test
     endif
     !notSynchronized = .True.
!!$     if( test>0) then
!!$        notSynchronized = .True.
!!$     endif
     if( test) then
        notSynchronized = .True.
     else
        notSynchronized = .false.
     endif
     if (count >6)then
     !if (count >1)then
        print *,'count',count
        exit!break
     endif
     if (myid ==0) then !myid
        print *,'count',count
        !count+=1
     endif
     count=count+1
     !endif
  end do !do while
       
  !stop

  contains
    
    subroutine getIncrement(nSubface)

      integer(kind=intType),intent(in):: nSubface
      integer(kind=intType)::i
      
      !begin execution
      
      ! Determine whether the coordinates are increasing or
      ! decreasing in each direction for each subface
      
      do i =1,nSubface
         !check for +ve vs -ve increment
         if (inend(i) >=inbeg(i)) then
            incrementI(i) = 1
         else
            incrementI(i) = -1
         endif
         
         if ( jnend(i) >= jnbeg(i)) then
            incrementJ(i) = 1
         else
            incrementJ(i) = -1
         endif
         
         if ( knend(i) >= knbeg(i)) then
            incrementK(i) = 1
         else
            incrementK(i) = -1
         endif
      end do
    end subroutine getIncrement

    subroutine getIncrementD(nSubface)

      integer(kind=intType),intent(in):: nSubface
      integer(kind=intType)::i
      
      !begin execution
      
      ! Determine whether the coordinates are increasing or
      ! decreasing in each direction for each subface
      
      do i =1,nSubface
         !check for +ve vs -ve increment
         if (dinend(i) >=dinbeg(i)) then
            incrementdI(i) = 1
         else
            incrementdI(i) = -1
         endif
         
         if ( djnend(i) >= djnbeg(i)) then
            incrementdJ(i) = 1
         else
            incrementdJ(i) = -1
         endif
         
         if ( dknend(i) >= dknbeg(i)) then
            incrementdK(i) = 1
         else
            incrementdK(i) = -1
         endif
      end do
    end subroutine getIncrementD

  end subroutine synchronizeBlockFaces

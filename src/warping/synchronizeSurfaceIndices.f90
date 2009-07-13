
!==========================================================
!
!       File: synchronizeSurfaceIndices.f90
!       Author: C.A.(Sandy) Mader
!       Start Date: 05-20-2009
!       Modified Date: 05-20-2009
!
!==========================================================

!     ******************************************************************
!     *   SynchronizeSurfaceIndices syncs global surface indices for   *
!     *   the parallel derivative computation                          *
!     *                                                                *
!     *   05/20/09  C.A.Mader   Initial Implementation                 *
!     *                                                                *
!     *   C.A.(Sandy) Mader, UTIAS,Toronto,Canada                      *
!     ******************************************************************



subroutine synchronizeSurfaceIndices(level,sps)
  !""" Synchronize surface indices for the parallel derivatives
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
  real(KIND=REALTYPE), dimension(3)::coords
  integer(kind=inttype),dimension(3)::indices
  integer(kind=intType),dimension(:),allocatable::incrementI,&
     incrementJ,incrementK,  incrementdI,&
     incrementdJ,incrementdK  
 

  !************************
  ! Begin Execution
  !************************

  !Set the inital state to not syncronized
  
  notSynchronized = .True.
  count = 0

  !allocate warp_comm storage
  do i=1,nDom
     call setPointers(i,level,sps)

     !allocate the memory for the block face communicators
     
     allocate(flowdoms(i,1,1)%warp_comm(nsubface), stat=ierr)
       
  end do

  do while (notSynchronized)
     ! Set state to synchronized and the check for truth at the end of the 
     ! loop.
     notSynchronized = .False.

     !loop over blocks and subfaces
     do i=1,nDom
    
        call setPointers(i,level,sps)

        ! Check to see if the memory is allocated to store the total number
        ! of nodes on each subface. If not, allocate.
        if(.not. associated(nNodesSubface)) then
           allocate(nNodesSubface(nsubface), stat=ierr)
        endif


        !allocate the memory for the block increments
        allocate(incrementI(nsubface),incrementJ(nsubface),incrementK(nsubface))
        allocate(incrementdI(nsubface),incrementdJ(nsubface),incrementdK(nsubface))
        !get the increment size and direction for each face
        call getIncrement(nSubface)
        ! and its donor
        call getIncrementD(nSubface)
        
        !Loop over the subfaces
        do j=1,nSubface
           
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
           
           !allocate the communicators for this subface
           !Generate the receive buffer
           ! Dimension 1-3 are coords., 4-6 are indices, 7 is remote block number
           if (.not. allocated(warp_comm(j)%recvBuffer)) &
                allocate(warp_comm(j)%recvBuffer(length))
           if (.not. allocated(warp_comm(j)%sendBuffer)) &
                allocate(warp_comm(j)%sendBuffer(length))
                      
           !Check for one to one matching internal faces
           if(BCType(j) == B2BMatch)then
              !Check that face requires interprocessor communication
              if(neighproc(j) /= myid)then
                 !need communication with blocks on another processor
                          
                 !post non blocking recieves for face communicator
                
                 !Generate a unique tag for this subface
                 recvtag = 10000*neighproc(j)+100*neighblock(j)+i
                                  
                 !post the non blocking recv
                 call mpi_irecv(warp_comm(j)%recvBuffer(1), length, sumb_real,&
                      neighproc(j), recvtag, sumb_comm_world, &
                      warp_comm(j)%recvreq, ierr)
                                                                   
                 !set a placement counter for the send buffer
                 counter=1
                                                   
                 !Loop over the local indices
                 do l =inmin,inmax,incrementI(j)
                    do m =jnmin,jnmax,incrementJ(j)
                       do n=knmin,knmax,incrementK(j)
                          !Set the counter step based on the coordinate transformation for the face
                          step(1) = (abs(l-inmin))
                          step(2) = (abs(m-jnmin))
                          step(3) = (abs(n-knmin))
                          
                          !Determine the indices in the recieving block
                          counterI =dinbeg(j)+step(abs(l1(j)))*incrementdI(j)
                          counterJ =djnbeg(j)+step(abs(l2(j)))*incrementdJ(j)
                          counterK =dknbeg(j)+step(abs(l3(j)))*incrementdK(j)
                          
                          !#set the value in the send buffer
                          warp_comm(j)%sendBuffer(counter)=x(l,m,n,1)
                          warp_comm(j)%sendBuffer(counter+1)=float(flowDoms(i,1,1)%cgnsBlockID)!x(l,m,n,2)
                          warp_comm(j)%sendBuffer(counter+2)=x(l,m,n,3)
                          warp_comm(j)%sendBuffer(counter+3)=float(counterI)
                          warp_comm(j)%sendBuffer(counter+4)=float(counterJ)
                          warp_comm(j)%sendBuffer(counter+5)=float(counterK)
                          warp_comm(j)%sendBuffer(counter+6)=float(neighblock(j))!-1
                          
                          counter=counter+7
                       enddo
                    enddo
                 enddo
                 
                 !Post the buffer send
                 sendtag = 10000*myid+100*(i)+neighblock(j)
                
                 call mpi_isend(warp_comm(j)%sendBuffer(1), length,sumb_real,&
                      neighproc(j),sendtag, sumb_comm_world, &
                      warp_comm(j)%sendreq, ierr)
                 
                                      
              else              
                 
                 !internal communication
                 !communicate the blocks that are on the local processor
                 
                 !#Loop over the local indices
                 do ii =inmin,inmax,incrementI(j)
                    do jj =jnmin,jnmax,incrementJ(j)
                       do kk=knmin,knmax,incrementK(j)
               
                          !#set the counter step based on the face coordinate transformation
                          step(1) = (abs(ii-inmin))
                          step(2) = (abs(jj-jnmin))
                          step(3) = (abs(kk-knmin))
                         
                          !#compute the neighbouring face counters
                          counterI =dinbeg(j)+step(abs(l1(j)))*incrementdI(j)
                          counterJ =djnbeg(j)+step(abs(l2(j)))*incrementdJ(j)
                          counterK =dknbeg(j)+step(abs(l3(j)))*incrementdK(j)
                          
                          !since this is only index communication only need to do 1st coord

                          do nn=1,1
                             neighbourindex = neighblock(j)
                             !reset pointers and get neighbour data
                             call setpointers(neighbourindex,level,sps)
                            
                             neighbour = x(counterI,counterJ,counterK,nn)
                            
                             !Reset pointers to the local block.
                             call setpointers(i,level,sps)
                             local = x(ii,jj,kk,nn)
                           
                             !check the various options and act accordingly
                             if( int(local) ==-5 )then
                                !#Then local face is not on surface
                                
                                !check if neighbour is...
                                if (int(neighbour)/=-5)then
                                   !neighbout is on surface, update my Global index.
                                   x(ii,jj,kk,nn) = neighbour
                                   notSynchronized = .True.
                                else
                                   !neither point is on surface, cycle
                                   a=1
                                endif
                             else
                                ! I am a surface point
                                if(int(neighbour)==-5)then
                                   !neighbour is not a surface,
                                   !keep current index
                                   local = local
                                  
                                else
                                   !neighbour is on a surface,
                                   !set both indices to lower value 
                                   if (local==int(neighbour))then
                                      !Both have correct index
                                      !do nothing
                                      a=1
                                   elseif (local<int(neighbour))then
                                      !I am the correct index
                                      ! other index will be corrected
                                      local = local
                                      notSynchronized = .True.
                                  
                                   else
                                      !neighbour block is lower
                                      x(ii,jj,kk,nn) = neighbour
                                  
                                      notSynchronized = .True.
                                   endif
                                end if
                             endif
                             
                          end do
                       end do
                    end do
                 end do
                 
              endif
           else
              !#no communicaton required
              !print *,'no communication required'
              a=1
           endif
           
           
        end do
        !allocate the memory for the block increments
        deallocate(incrementI,incrementJ,incrementK)
        deallocate(incrementdI,incrementdJ,incrementdK)
     end do
     
     call mpi_barrier(sumb_comm_world, ierr)
       
 
     !loop over blocks and subfaces
     do i=1,nDom
        call setPointers(i,level,sps)
 
        do j=1,nSubface 
           !only check for receives on procesors that posted them
           if(BCType(j) == B2BMatch)then
              !      #Check that face requires interprocessor communication
              if(neighproc(j) /= myid)then
                 !#need communication with blocks on another processor
                 call mpi_wait(warp_comm(j)%recvreq,MPI_STATUS_IGNORE,ierr)
              endif
           endif
        enddo
     enddo
     
     call mpi_barrier(sumb_comm_world, ierr)

     
!            #Wait for all of the faces to be comunicated
!            MPI.Request.Waitall(reqList)
            

     !Now loop over the recieved buffers and update the local blocks
    
     !loop over blocks and subfaces
     do l=1,nDom
        call setPointers(l,level,sps)
      
        do m=1,nSubface
           ! Check for one to one matching internal faces
           
           if(BCType(m) == B2BMatch)then
              !#Check that face requires interprocessor communication
              if(neighproc(m) /= myid)then

                 !#need communication with blocks on another processor   
                 length = nNodesSubface(m)
                 
                 do i=1,int(length)
                    index = 1+(i-1)*7
                    !Set the destination block index
                    destblock = int(warp_comm(m)%recvbuffer(index+6))
                    
                    !get the surface indices (coord(1)) and location 
                    !(indices(:)) from the receive buffer
                    coords(:) =  warp_comm(m)%recvbuffer(index:index+2)        
                    indices(:) = int(warp_comm(m)%recvbuffer(index+3:index+5))
                    
                    !Set pointers for the given destblock
                    call setPointers(destblock,level,sps)
                    !need only the first coordinate because this 
                    !is only index communication
                    
                    do nn=1,1
                      
                       local = x(indices(1),indices(2),indices(3),nn)
                       realbuf = coords(nn)!real(coords(nn))
                       
                       !check the various options and act accordingly
                       if( int(local) ==-5 )then
                          !#Then local face is not on surface
                          !check if neighbour is...
                          if (int(realbuf)/=-5)then
                             !neighbour is on surface, update my Global index.
                             x(indices(1),indices(2),indices(3),nn) = realbuf!neighbour
                             notSynchronized = .True.
                             
                          else
                             !neither point is on surface, cycle
                             a=1
                          endif
                       else
                          ! I am a surface point
                          if(int(realbuf)==-5)then
                             !neighbour is not a surface,
                             !keep current index
                             local = local
                             notSynchronized = .True.
                             
                          else
                             !neighbour is on a surface,
                             !set both indices to value of
                             !lower numbered block
                             if (int(local)==int(realbuf))then
                               !Both have correct index
                                !do nothing
                                a=1
                             elseif (int(local)<int(realbuf))then
                                !I am the correct index
                                ! other index will be corrected
                                local = local
                                notSynchronized = .True.

                             else
                                !neighbour block is lower
                                x(indices(1),indices(2),indices(3),nn) = int(realbuf)!neighbour
                                notSynchronized = .True.
                     
                             endif
                          end if
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
     !Collect logicals for synchronization check.
     call mpi_allreduce(notSynchronized,test,1,MPI_LOGICAL,MPI_LOR,sumb_comm_world,ierr)
     
     if (myid ==0) then !myid
        print *,'Synchronization test',test
     endif

     !if synchrinized then set logical to break loop...
     if( test) then
        notSynchronized = .True.
     else
        notSynchronized = .false.
     endif
     !loop limiter...
     if (count >6)then
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
  
  !deallocate warp_comm storage
  do i=1,nDom
     
     call setPointers(i,level,sps)
     !allocate the memory for the block face communicators
     
     deallocate(flowdoms(i,1,1)%warp_comm, stat=ierr)
     
  end do

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

  end subroutine synchronizeSurfaceIndices

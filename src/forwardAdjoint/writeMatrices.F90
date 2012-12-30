! Debugging routines for dumping dRdw and dRdx matrices to files for
! comparison. These routines are to be called from python. Basically,
! the idea is to run routine, assembling dRdw or dRdx, write the file,
! assemble another and then dump a new one. 

subroutine writedRdx(file_name)

  ! This subroutine write the current state of dRdx to a file

  use ADjointPETSc
  use blockPointers      
  use inputDiscretization 
  USE inputTimeSpectral 
  use iteration         
  use flowVarRefState    
  use stencils
  implicit none
  character*(*) :: file_name

  integer(kind=intType) :: icell,jcell,kcell,inode,jnode,knode
  integer(kind=intType) :: nn,sps,i,j,ii,jj,kk,ierr,irow,icol
  real(kind=realType) :: vals(nw,3)

  integer(kind=intType) :: n_stencil,i_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: rowInds(nw), colInds(3)
  logical :: boundary 

  ! Open the file we want to write. 
  open (UNIT=13,File=trim(file_name),status='replace',action='write')

  if( viscous ) then
     !stencil => visc_drdx_stencil
     !n_stencil = N_visc_drdx
  else 
     stencil => euler_drdx_stencil
     n_stencil = N_euler_drdx
  end if

  do nn=1,1
     do sps =1,nTimeIntervalsSpectral
     call setPointers(nn,1,sps)

     write(13,10),nn,sps
     write(13,15),il,jl,kl
10   format(1x,'Block:',I4, 'Spectral Instance:',I4)
15   format(2x,'Block Bounds: ',I4,' ',I4, ' ',I4)

       do knode=1,3
          do jnode=1,3
             do inode=1,3
                boundary = .False.
20              format(3x,'Node:',I3,' ',I3,' ',I3)
                
                if (inode == 1 .or. inode == il)  boundary = .True. 
                if (jnode == 1 .or. jnode == jl)  boundary = .True. 
                if (knode == 1 .or. knode == kl)  boundary = .True. 

                write(13,20),inode,jnode,knode

                if (boundary) then
                   write(13,*), 'Boundary Node'
                end if

                icol = globalNode(inode,jnode,knode)

                do i_stencil=1,n_stencil
                   ii = stencil(i_stencil,1)
                   jj = stencil(i_stencil,2)
                   kk = stencil(i_stencil,3)

                   if ( inode+ii >= 2 .and. icell+ii <= il .and. &
                        jnode+jj >= 2 .and. jnode+jj <= jl .and. &
                        knode+kk >= 2 .and. knode+kk <= kl) then 

                   write(13,30),ii,jj,kk
30                 format(5x,'Sten:',I3,' ',I3,' ',I3)

                      irow = globalCell(inode+ii,jnode+jj,knode+kk)
                      ! Extract all nw*3 values from the MatGetValues

                      ! Row Indices
                      do i=1,nw
                         rowInds(i) = iRow*nw+i-1
                      end do

                      ! Colume Indices
                      do j=1,3
                         colInds(j) = iCol*3+j-1
                      end do

                      ! NOTE: dRdx is actually stored as the
                      ! TRANSPOSE!. We flip row and col into the call

                      call MatGetValues(dRdx,3,colInds,nw,rowInds,vals,ierr)
                      call EChk(ierr,__FILE__,__LINE__)

40                    format(7x,'x ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ')
50                    format(7x,'y ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ')
60                    format(7x,'z ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ',g15.8,' ')
                        
                      ! Zero out really small values
                      do i=1,5
                         do j=1,3
                            if (abs(vals(i,j))< 1e-10)then
                               vals(i,j) = zero
                            end if
                         end do
                      end do
                      ! Write out all 15  values. 
                      write(13,40),vals(:,1)
                      write(13,50),vals(:,2)
                      write(13,60),vals(:,3)

                      
                   end if
                end do
             end do
          end do
       end do
    end do
 end do
 ! Finally close the file
 close(13)
end subroutine writedRdx
  

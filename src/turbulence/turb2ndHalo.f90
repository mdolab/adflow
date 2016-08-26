!
!       File:          turb2ndHalo.f90                                 
!       Author:        Edwin van der Weide                             
!       Starting date: 06-16-2003                                      
!       Last modified: 06-12-2005                                      
!
       subroutine turb2ndHalo(nn)
!
!       turb2ndHalo sets the turbulent variables in the second halo    
!       cell for the given subface. Simple constant extrapolation is   
!       used to avoid problems.                                        
!
       use constants
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l
!
!       Begin execution                                                
!
       ! Determine the face on which this subface is located and set
       ! some pointers accordingly.

       ! Loop over the turbulent variables and set the second halo
       ! value. If this is an eddy model, also set the eddy viscosity.

       select case (BCFaceID(nn))
         case (iMin)
            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(0,i,j,l) = w(1,i,j,l)
                  enddo
                  if( eddyModel ) rev(0,i,j) = rev(1,i,j)
               enddo
            enddo
            
         !===============================================================

         case (iMax)
            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(ib,i,j,l) = w(ie,i,j,l)
                  enddo
                  if( eddyModel ) rev(ib,i,j) = rev(ie,i,j)
               enddo
            enddo
            
         !===============================================================

         case (jMin)
            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(i,0,j,l) = w(i,1,j,l)
                  enddo
                  if( eddyModel ) rev(i,0,j) = rev(i,1,j)
               enddo
            enddo

         !===============================================================

         case (jMax)
            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(i,jb,j,l) = w(i,je,j,l)
                  enddo
                  if( eddyModel ) rev(i,jb,j) = rev(i,je,j)
               enddo
            enddo

         !===============================================================

         case (kMin)
            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(i,j,0,l) = w(i,j,1,l)
                  enddo
                  if( eddyModel ) rev(i,j,0) = rev(i,j,1)
               enddo
            enddo

         !===============================================================

         case (kMax)

            do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                  do l=nt1,nt2
                     w(i,j,kb,l) = w(i,j,ke,l)
                  enddo
                  if( eddyModel ) rev(i,j,kb) = rev(i,j,ke)
               enddo
            enddo

       end select

       end subroutine turb2ndHalo

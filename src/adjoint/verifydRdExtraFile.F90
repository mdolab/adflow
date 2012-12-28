!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdExtraFile.f90                          *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 13-07-2011                                      *
!     * Last modified: 13-07-2011                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdExtraFile(level)
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine prints the values of dRda for comparison      *
!     *  to CS results                                                 *
!     *                                                                *
!     ******************************************************************

#ifndef USE_NO_PETSC
      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate !nw
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr,nTimeIntervalsSpectral
      use inputIO

      use ADjointPETSc, only: drda,insert_values,petscierr,&
           mat_final_assembly,petsc_viewer_draw_world
      use ADjointVars
      use precision

      implicit none

!
!     Subroutine arguments
      integer(kind=intType), intent(in) :: level

!
!     Local Variables
!
      integer(kind=intType) ::icolor,sps,nn
      integer(kind=intType) ::i,j,k,n,idxres
      real(kind=realType):: value1

      !File Parameters
      integer :: unitdRda = 8,ierr
      character(len = 128)::outfile,testfile

      write(testfile,100) myid
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("ADdRdafile",a,".out")
      
      unitdRda = 8+myID
            
      open (UNIT=unitdRda,File=outfile,status='replace',action='write',iostat=ierr)
      if(ierr /= 0)                        &
           call terminate("verifydRdExtraFile", &
           "Something wrong when &
           &calling open")

      ! VerifydRdExtraFile portion of script - printing to file ADdRdaFile0.out

      do icolor = 1,nDesignExtra
         do sps = 1,nTimeIntervalsSpectral
            do nn = 1,ndom
               call setPointersAdj(nn,1,sps)
               DO I=2,Il
                  DO J=2,Jl
                     DO K=2,Kl
                        do n = 1,nw
                           idxres = globalCell(i,j,k)*nw+n                                   
                          
                           call MatGetValues(drda,1,idxres-1,1,icolor-1,value1,ierr)
                           !call MatGetValues(drdw,1,idxres-1,1,idxstate-1,value1,ierr)
                           call EChk(ierr,__FILE__,__LINE__)
                           if(abs(value1)>1e-10)then!-7)then!
                              !print *,'value',value1,icolor,idxres
                              write(unitdrda,13) icolor,idxres,n,k,j,i,nn,sps,real(value1)
13                            format(1x,'drda',8I8,f18.10)
                           endif
                        enddo
                     END DO
                  END DO
               END DO
            end do
         end do
      end do
      
      close(unitdRda)
#endif

    end subroutine verifydRdExtraFile


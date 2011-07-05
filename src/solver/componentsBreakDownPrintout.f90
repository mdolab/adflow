!
!      ******************************************************************
!      *                                                                *
!      * File:          componentBreakDownPrintout.f90                  *
!      * Author:        Eran Arad                                       *
!      * Starting date: 09-12-2006                                      *
!      * Last modified:                                                 *
!      *                                                                *
!      ******************************************************************
!
       subroutine componentsBreakDownPrintout(writeToFile)
!
!      ******************************************************************
!      *                                                                *
!      * Called by convergenceInfo computes and writes the components   *
!      * contribution to forces/moments to standard output.             *  
!      * input variable: writeToFile = 1 .... CBD data written on file  *
!      *                                 cbdFileName                    *
!      *                               0 .... CBD data written on output*
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use cgnsNames
       use inputTimeSpectral
       use monitor
       use iteration
       use cgnsGrid 
       use inputIO
 
       implicit none

       integer :: writeToFile
!
!      Local variables.
!
       integer :: ierr, lenParamFile, iPoint,nameLength, CBDUnit=37, ios=0, &
            outUnit 

       integer(kind=intType) :: sps, nn, i, j, dataTransSize

       character(len=maxStringLen) :: cbdFileName

!
!  cFpA .... value in domain
!  cFpL .... value in processor
!  cFpG .... Global value
!
       real(kind=realType), allocatable, dimension(:,:) :: cfpA, cfvA, cmpA, cmvA 
       real(kind=realType), allocatable, dimension(:,:) :: cfpL, cfvL, cmpL, cmvL 
       real(kind=realType), allocatable, dimension(:,:) :: cfpG, cfvG, cmpG, cmvG 


!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!
!  allocate memory for the forces and moments coefficients up
!  to the number of wall surfaces in the CGNS file, to enable 
!  components contribution break down. 
!  The 0 location serves for the total value (sum of contribution of all surfaces)
!
       allocate(cfpA(3,0:cgnsNWallSurfaces))
       allocate(cfvA(3,0:cgnsNWallSurfaces))
       allocate(cmpA(3,0:cgnsNWallSurfaces))
       allocate(cmvA(3,0:cgnsNWallSurfaces))
       allocate(cfpL(3,0:cgnsNWallSurfaces))
       allocate(cfvL(3,0:cgnsNWallSurfaces))
       allocate(cmpL(3,0:cgnsNWallSurfaces))
       allocate(cmvL(3,0:cgnsNWallSurfaces))
       allocate(cfpG(3,0:cgnsNWallSurfaces))
       allocate(cfvG(3,0:cgnsNWallSurfaces))
       allocate(cmpG(3,0:cgnsNWallSurfaces))
       allocate(cmvG(3,0:cgnsNWallSurfaces))


       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

!
!  Initialize the force and moment coefficients to 0
!
          direction : do j=1,3
             surfacesCounter : do i=0,cgnsNWallSurfaces
                cfpL(j,i)= zero ; cfvL(j,i) = zero
                cmpL(j,i)= zero ; cmvL(j,i) = zero
             end do surfacesCounter
          end do direction


         ! Loop over the blocks.

          domains: do nn=1,nDom

           ! Set the pointers for this block.

             call setPointers(nn, groundLevel, sps)

           ! Compute the forces and moments for this block.

             call forcesAndMomentsCBD(nn,cfpA, cfvA, cmpA, cmvA) 

             direction2 : do j=1,3
                surfacesCounter2 : do i=0,cgnsNWallSurfaces  
                   cFpL(j,i) = cFpL(j,i) + cFpA(j,i)
                   cFvL(j,i) = cFvL(j,i) + cFvA(j,i)
                   cMpL(j,i) = cMpL(j,i) + cMpA(j,i)
                   cMvL(j,i) = cMvL(j,i) + cMvA(j,i)
                end do surfacesCounter2
             end do direction2

          enddo domains

         GroundLevelCh : if(groundLevel == 1)then

            dataTransSize = 3*(cgnsNWallSurfaces+1)
            call mpi_reduce(cFpL,cFpG,dataTransSize,sumb_real, &
                 mpi_sum, 0, SUmb_comm_world, ierr)
            call mpi_reduce(cFvL,cFvG,dataTransSize,sumb_real, &
                 mpi_sum, 0, SUmb_comm_world, ierr)
            call mpi_reduce(cMpL,cMpG,dataTransSize,sumb_real, &
                 mpi_sum, 0, SUmb_comm_world, ierr)
            call mpi_reduce(cMvL,cMvG,dataTransSize,sumb_real, &
                 mpi_sum, 0, SUmb_comm_world, ierr)

         end if GroundLevelCh
         

         ! Write the convergence info; only processor 0 does this.

         testRootProc: if(myID == 0) then
!
! determine output media: file or screen
!
            if(writeToFile == 0)then
               outUnit = 6 ! print on screen
            else
!
! define a name and create a CBD out file
!
               outUnit = CBDUnit
               lenParamFile = len_trim(paramFile)
               iPoint=index(paramFile,'.')
               nameLength=lenParamFile
               if (iPoint > 0 )then
                  nameLength=iPoint-1
               end if

               cbdFilename=paramFile(1:nameLength)//'_CBDOUT.dat'
               open(unit=CBDUnit, file=trim(cbdFilename), status="unknown", &
                    action="write", position='rewind',iostat=ios)
               if(ios /= 0)call terminate("componentsBreakDownPrintout",&
                    "failed to open CBDOut file") 
            end if

            nIterCur = nIterCur + nIterOld

!
!----------------- Components Break Down printout
!
              write(outUnit,'("#     ======== COMPONENTS BREAK DOWN: @Iteration ",i5," ============")')&
                   nIterCur
              write(outUnit,'("    convergenceQuality : ",i2)')convergenceQuality
              write(outUnit,'("#====================================================================")') 
              write(outUnit,'("# ")')

              surfacesCounter3 : do i=1,cgnsNWallSurfaces
                 write(outUnit,'("#   Wall surface ",i3," : ",a)')i,wallBCNames(i)
                 write(outUnit,'("# ",50(1h-))')
                 write(outUnit,'("#   cFx          cFy            cFz  ",&
                   "       cMx          cMy          cMz   ")')
              write(outUnit,'("#  --------     --------       ------ ",&
                   "     --------     --------     -------- ")')

                 write(outUnit,'("  ",6(1pg12.5,1x))')&
                    (cFpG(j,i)+cFvG(j,i),j=1,3),(cMpG(j,i)+cMvG(j,i),j=1,3  )
                 write(outUnit,'("# ")')
              end do surfacesCounter3


              write(outUnit,'("# ")')
              write(outUnit,'("# ")')
            

              write(outUnit,'("#   Forces/moments sum (total)")')
              write(outUnit,'("# ",50(1h-))')
              i=0
              write(outUnit,'("#    cFx          cFy            cFz  ",&
                   "       cMx          cMy          cMz   ")')
              write(outUnit,'("#  --------     --------       ------ ",&
                   "     --------     --------     -------- ")')
              write(outUnit,'("  ",6(1pg12.5,1x))')&
                    (cFpG(j,i)+cFvG(j,i),j=1,3),(cMpG(j,i)+cMvG(j,i),j=1,3  )

              if(writeToFile == 1)then
                 close(unit=CBDUnit,status='keep', iostat=ios)
                 if( ios == 0)then
                    write(*,*)"CBD file ",trim(cbdFilename),' created'
                 else
                    call terminate("componentsBreakDownPrintout",&
                         "failed to close CBDOut file") 
                 end if
              end if

         endif testRootProc

       enddo spectralLoop


       deallocate(cfpA,cfvA,cmpA,cmvA,cFpL,cFvL,cMpL,cMvL,&
            cFpG,cFvG,cMpG,cMvG) ! eran-CBD


       end subroutine componentsBreakDownPrintout


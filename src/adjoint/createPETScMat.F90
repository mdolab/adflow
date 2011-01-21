
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScMat.F90                              *
!     * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine createPETScMat
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matries: dRdWT,dRdwTPre,dRdx,dRda,dFdw,dFdx         *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState ! 
  use inputADjoint    !ApproxPC

  implicit none
  !
  !     Local variables.
  !
  integer       :: nDimW, nDimX,nDimS
  integer :: i
  integer, dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer, dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  ! PETSc macros are lost and have to be redefined.
  ! They were extracted from: <PETSc root dir>/include/petscmat.h
#define MATMPIDENSE        "mpidense"
#define MATMPIAIJ          "mpiaij"
  ! Define matrix dRdW local size, taking into account the total
  ! number of Cells owned by the processor and the number of 
  ! equations.

  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral

  ! Define matrix dRdx local size (number of columns) for the
  ! spatial derivatives.

  nDimX = 3 * nNodesLocal*nTimeIntervalsSpectral

  ! Define matrix dFdx local size (number of Rows) for the
  ! Coupling derivatives.

  call getForceSize(nDimS)
  nDimS = nDimS * 3 ! Multiply by 3 for each dof on each point

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdW that define the adjoint linear system of    *
  !     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDimW,nDimW]*
  !     * but is very sparse because of the computational stencil R=R(W).*
  !     *                                                                *
  !     ******************************************************************
  !

  allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
       nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )
  
  call drdwPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal*nTimeIntervalsSpectral)
  if( nw <= 7 ) then

     PETScBlockMatrix = .true.

     ! Create a Block AIJ Matrix with block size nw, (number of states)
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  else

     PETScBlockMatrix = .false.
     
     allocate(nnzDiagonal2(nDimw),nnzOffDiag2(nDimw))
     ! The drdw prealloc function is done per block, which we can use
     ! to compute the correct preallocation for the non-block matrix:

     do i=1,nCellsLocal*nTimeIntervalsSpectral
        nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
        nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
     end do
     call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          8, nnzDiagonal2,         &
          8, nnzOffDiag2,            &
          dRdWT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     deallocate(nnzDiagonal2,nnzOffDiag2)
  endif

  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set the matrix dRdW options.

  ! Warning: The array values is logically two-dimensional, 
  ! containing the values that are to be inserted. By default the
  ! values are given in row major order, which is the opposite of
  ! the Fortran convention, meaning that the value to be put in row
  ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
  ! the insertion of values in column major order, one can call the
  ! command MatSetOption(Mat A,MAT COLUMN ORIENTED);

#ifdef USE_PETSC_3
  call MatSetOption(dRdWt, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#else
  call MatSetOption(dRdWt, MAT_COLUMN_ORIENTED, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif
  !****************
  !create dRdWPre
  !***************

  if (ApproxPC) then

     !     ******************************************************************
     !     *                                                                *
     !     * Create matrix dRdWPre for the ADjoint approximate preconitioner*
     !     * This matrix is sparse with a narrower bandwidth than drdw      *
     !     *                                                                *
     !     ******************************************************************
     !

     allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
          nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )
     
     call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal*nTimeIntervalsSpectral)

     if( nw <= 7 ) then

        PETScBlockMatrix = .true.

        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWPreT, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     else

        PETScBlockMatrix = .false.
        allocate(nnzDiagonal2(nDimw),nnzOffDiag2(nDimw))
        ! The drdwPC prealloc function is done per block, which we can use
        ! to compute the correct preallocation for the non-block matrix:
        
        do i=1,nCellsLocal*nTimeIntervalsSpectral
           nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
           nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
        end do

        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal2,         &
             0,nnzOffDiag2,            &
             dRdWPret, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        deallocate(nnzDiagonal2,nnzOffDiag2)
     endif

     deallocate( nnzDiagonal, nnzOffDiag )

     ! Set the matrix dRdWPre options.

#ifdef USE_PETSC_3
     call MatSetOption(dRdWPret, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
#else
     call MatSetOption(dRdWPret, MAT_COLUMN_ORIENTED, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
#endif
  end if ! Approx PC

  ! dRda

  call MatCreate(SUMB_PETSC_COMM_WORLD, dRda, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetSizes(dRda, nDimW, PETSC_DECIDE, &
       PETSC_DETERMINE, nDesignExtra, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetType(dRda,MATMPIDENSE,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
  call MatSetOption(dRda, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#else
  call MatSetOption(dRda, MAT_COLUMN_ORIENTED, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdx that is used to compute the total cost /    *
  !     * constraint function sensitivity with respect to the spatial    *
  !     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
  !     *                                                                *
  !     * Matrix dRdx has size [nDimW,nDimX] and is generally            *
  !     * sparse for the coordinate design variables.                    *
  !     *                                                                *
  !     * The local dimensions are specified so that the spatial         *
  !     * coordinates x (a) are placed in the local processor. This has  *
  !     * to be consistent with the vectors dIdx and dJdx.               *
  !     *                                                                *
  !     ******************************************************************

  allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )
  ! Create the matrix dRdx.
  call drdxPreAllocation(nnzDiagonal,nnzOffDiag,nDimW)
  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
       nDimW, nDimX,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       8, nnzDiagonal,     &
       8, nnzOffDiag,            &
       dRdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )

  call MatSetFromOptions(dRdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Set column major order for the matrix dRdx.
#ifdef USE_PETSC_3
  call MatSetOption(dRdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetOption(dRdx,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,PETScIErr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#else
  call MatSetOption(dRdx, MAT_COLUMN_ORIENTED, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif

  ! Create dFdx and dFdw

  ! Each nodal force is contribed by (nominally) 4 quadrilateral cells
  ! surrounding it. The pressure on each of these cells depend on the
  ! average of the pressure of the cell above and the halo below. The
  ! 1st halo is computed by linear pressure extrapolation from the
  ! first two cells ABOVE the surface. The results in each node being
  ! affected by 8 cells. All of these cells on on-processors so there
  ! should be zero offdiag entries. For dFdx, the coordinates of each
  ! of the 4 quadrilateral affects the force, so this results in 9
  ! points spatial points affecting the force. Each coordiante has 3
  ! dimension which results in 3x3x3=27 nonzeros per row

  ! don't know where the non-zeros will end up

  ! Create the matrix dFdw

  allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )
  nnzDiagonal = 8*nw
  nnzOffDiag  = 8*nw! Make the off diagonal the same, since we
  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
       nDimS, nDimW,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       0, nnzDiagonal,         &
       0, nnzOffDiag,            &
       dFdw, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Create the matrix dFdx
  nnzDiagonal = 27
  nnzOffDiag = 27
  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
       nDimS, nDimS,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       0, nnzDiagonal,         &
       0, nnzOffDiag,            &
       dFdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dFdw.
#ifdef USE_PETSC_3
  call MatSetOption(dFdw, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetOption(dFdx, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#else
  call MatSetOption(dFdw, MAT_ROW_ORIENTED, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetOption(dFdx, MAT_ROW_ORIENTED, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif

  call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)  

#endif

end subroutine createPETScMat

subroutine drdwPreAllocation(onProc,offProc,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use BCTypes
  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: onProc(wSize),offProc(Wsize)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,ii,sps

  ii = 0
  onProc(:) = 1+(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  offProc(:) = 0_intType 
  
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 ii = ii + 1
                 if (i-2 < 2) then
                    call checkCell(iMin,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (i-1 < 2) then 
                    call checkCell(iMin,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (i+1 > il) then 
                    call checkCell(iMax,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if
                 if (i+2 > il) then
                    call checkCell(iMax,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j-2 < 2) then 
                    call checkCell(jMin,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j-1 < 2) then 
                    call checkCell(jMin,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j+1 > jl) then 
                    call checkCell(jMax,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j+2 > jl) then 
                    call checkCell(jMax,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k-2 < 2) then 
                    call checkCell(kMin,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k-1 < 2) then 
                    call checkCell(kMin,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k+1 > kl) then 
                    call checkCell(kMax,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k+2 > kl) then 
                    call checkCell(kMax,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps Loop
  end do ! Domain Loop


end subroutine drdwPreAllocation

subroutine drdwPCPreAllocation(onProc,offProc,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use BCTypes

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: onProc(wSize),offProc(wSize)


  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,sps,ii
  ii = 0
  onProc(:) = 1+(nTimeIntervalsSpectral) ! ALWAYS have the center cell ON-PROCESSOR
  offProc(:) = 0_intType 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 ii = ii + 1
                 if (i-1 < 2) then 
                    call checkCell(iMin,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (i+1 > il) then 
                    call checkCell(iMax,j,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j-1 < 2) then 
                    call checkCell(jMin,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (j+1 > jl) then 
                    call checkCell(jMax,i,k,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k-1 < 2) then 
                    call checkCell(kMin,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if

                 if (k+1 > kl) then 
                    call checkCell(kMax,i,j,onProc(ii),offProc(ii),1)
                 else
                    onProc(ii) = onProc(ii) + 1
                 end if
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop
end subroutine drdwPCPreAllocation

subroutine drdxPreAllocation(onProc,offProc,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use BCTypes

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: onProc(wSize),offProc(wSize)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,l,sps,ii

  onProc(:) = 8*3+8*3*(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  offProc(:) = 0_intType 

  ii = 0 

  ! This is for the "Regular" drdx calculation. i.e. xadjb
 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 do l=1,nw
                    ii = ii + 1
                    if (i-1 < 2) then 
                       call checkCell(iMin,j,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (i+1 > il) then 
                       call checkCell(iMax,j,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (j-1 < 2) then 
                       call checkCell(jMin,i,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (j+1 > jl) then 
                       call checkCell(jMax,i,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (k-1 < 2) then 
                       call checkCell(kMin,i,j,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (k+1 > kl) then 
                       call checkCell(kMax,i,j,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if
                 end do ! l loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

  ! However, drdx is more complex since we ALSO have
  ! xblockcorners. These however, only show up for the cells that are
  ! along a symmetry plane. Lets try to estimate those. 

  ! THIS MAY NOT WORK FOR SPS CASE!!! 
  ii = 0
 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    ii = ii + 1
                    call checkCellSym(i,j,k,onProc(ii),24)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine drdxPreAllocation

subroutine checkCell(iface,i,j,onProc,offProc,addVal)

  use blockPointers
  use BCTypes
  use communication
  implicit None

  ! Subroutine Arguments

  integer(kind=intType), intent(in) :: iface,i,j
  integer(kind=intType), intent(inout) :: onProc,offProc
  integer(kind=intType), intent(in) :: addVal

  !local Variables
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd,mm,ll

  ! It is assumed blockPointers are already set for this block

  ! Basically what we want to do is take the cell defined by index i
  ! and j on face defined by iface and determine:
  ! 1. What block subface it is on.
  ! 2. Determine if this is a boundary condition or a block-match
  ! 3. If its a block-match is the connecting block on- or off-processor


  ! If its a BC condition, nothing is added onProc and offProc
  ! If its a block-match on-proc, addVal is added to onProc
  ! If its a block-match off-proc, addVal is added to offProc

  n1to1Loop: do mm=1,nsubFace

     ! Store the correct index for this subface, i.e. add the
     ! offset from the boundary subfaces.

     if (BCFaceID(mm) == iface) then ! Check Face
        if (iface == iMin .or. iface == iMax) then
           iBeg = jcBeg(mm) ; iEnd = jcEnd(mm)
           jBeg = kcBeg(mm) ; jEnd = kcEnd(mm)
        else if(iface == jMin .or. iface == jMax) then
           iBeg = icBeg(mm) ; iEnd = icEnd(mm)
           jBeg = kcBeg(mm) ; jEnd = kcEnd(mm)
        else
           iBeg = icBeg(mm) ; iEnd = icEnd(mm)
           jBeg = jcBeg(mm) ; jEnd = jcEnd(mm)
        end if

        ! Check to make sure cell is on this (possible) sub-face
        if (i>=iBeg .and. i<=iEnd .and. j>=jBeg .and. j<= jEnd) then

           if (neighproc(mm) == myid) then
              onProc = onProc + addVal
           else
              if (neighproc(mm) >=0) then
                 offProc = offProc + addVal
              end if
           end if
        end if

     end if
  end do n1to1Loop
end subroutine checkCell

subroutine checkCellSym(i,j,k,onProc,addVal)

  use blockPointers
  use BCTypes
  use communication
  implicit None

  ! Subroutine Arguments

  integer(kind=intType), intent(in) :: i,j,k
  integer(kind=intType), intent(inout) :: onProc
  integer(kind=intType), intent(in) :: addVal

  !local Variables
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd,mm,ll

  n1to1Loop: do mm=1,nsubFace
     if (BCType(mm) == Symm) then 
        ! Is cell i,j,k "on" the symmetry plane
        select case (BCFaceID(mm))
        case (iMin)
           if (i == 2) then
              onProc = onProc + addVal
           end if
        case (iMax)
           if (i == il) then
              onProc = onProc + addVal
           end if
        case (jMin)
           if (j == 2) then
              onProc = onProc + addVal
           end if
        case (jMax)
           if (j == jl) then
              onProc = onProc + addVal
           end if
        case (kMin)
           if (k == 2) then
              onProc = onProc + addVal
           end if
        case (kMax)
           if (k == kl) then
              onProc = onProc + addVal
           end if
        end select
     end if
  end do n1to1Loop
end subroutine checkCellSym

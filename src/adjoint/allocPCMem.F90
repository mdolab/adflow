subroutine allocPCMem(level)

  ! This routine allocates memory for the fortran-based PC. It is
  ! currently not used anywhere, but it become useful in the future.
  
  use constants
  use blockPointers, only : nDom, nx, ny, nz, il, jl, kl, ie, je, ke, flowDoms
  use inputTimeSpectral, onlY : nTimeIntervalsSpectral
  use flowVarRefState, only : nw
  implicit none

  integer(kind=intType), intent(in) :: level
  integer(kind=intType) :: nn, sps
  
  
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)

        if (.not. associated(flowDoms(nn, level, sps)%PCMat)) then 

           allocate(flowDoms(nn, level, sps)%PCMat(nw, 7*nw, 2:il, 2:jl, 2:kl))
           
           ! I direciton data
           allocate(flowDoms(nn, level, sps)%i_D_Fact(nw, nx*nw, 2:jl, 2:kl))
           allocate(flowDoms(nn, level, sps)%i_L_Fact(nw, (nx-1)*nw, 2:jl, 2:kl))
           allocate(flowDoms(nn, level, sps)%i_U_Fact(nw, (nx-1)*nw, 2:jl, 2:kl))
           allocate(flowDoms(nn, level, sps)%i_U2_Fact(nw, (nx-2)*nw, 2:jl, 2:kl))
           
           ! J direction data
           allocate(flowDoms(nn, level, sps)%j_D_Fact(nw, ny*nw, 2:il, 2:kl))
           allocate(flowDoms(nn, level, sps)%j_L_Fact(nw, (ny-1)*nw, 2:il, 2:kl))
           allocate(flowDoms(nn, level, sps)%j_U_Fact(nw, (ny-1)*nw, 2:il, 2:kl))
           allocate(flowDoms(nn, level, sps)%j_U2_Fact(nw, (ny-2)*nw, 2:il, 2:kl))

           ! K direciton data
           allocate(flowDoms(nn, level, sps)%k_D_Fact(nw, nz*nw, 2:il, 2:jl))
           allocate(flowDoms(nn, level, sps)%k_L_Fact(nw, (nz-1)*nw, 2:il, 2:jl))
           allocate(flowDoms(nn, level, sps)%k_U_Fact(nw, (nz-1)*nw, 2:il, 2:jl))
           allocate(flowDoms(nn, level, sps)%k_U2_Fact(nw, (nz-2)*nw, 2:il, 2:jl))
           
           ! iPIV arrays
           allocate(flowDoms(nn, level, sps)%i_ipiv(nw, nx, 2:jl, 2:kl))
           allocate(flowDoms(nn, level, sps)%j_ipiv(nw, ny, 2:il, 2:kl))
           allocate(flowDoms(nn, level, sps)%k_ipiv(nw, nz, 2:il, 2:jl))
           
           ! Vectors
           allocate(flowDoms(nn, level, sps)%pcVec1(nw, 1:ie, 1:je, 1:ke))
           allocate(flowDoms(nn, level, sps)%pcVec2(nw, 1:ie, 1:je, 1:ke))

        end if
     end do
  end do
end subroutine allocPCMem

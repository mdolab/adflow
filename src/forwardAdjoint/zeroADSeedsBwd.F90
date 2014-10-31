! This function zeros all 'd' values of on the specifified block/level

subroutine zeroADSeedsBwd(nn, level, sps)

  use blockPointers_b ! This modules includes blockPointers
  use inputTimeSpectral
  use flowVarRefState
  implicit none

  ! Input parameters
  integer(kind=intType) :: nn, level, sps

  ! Working parameters
  integer(kind=intType) :: mm

  flowDomsb(nn, level, sps)%d2wall = zero	
  flowDomsb(nn, level, sps)%x = zero
  flowDomsb(nn, level, sps)%si = zero
  flowDomsb(nn, level, sps)%sj = zero
  flowDomsb(nn, level, sps)%sk = zero
  flowDomsb(nn, level, sps)%vol = zero
  
  flowDomsb(nn, level, sps)%s = zero
  flowDomsb(nn, level, sps)%sFaceI = zero
  flowDomsb(nn, level, sps)%sFaceJ = zero
  flowDomsb(nn, level, sps)%sFaceK = zero
  
  flowDomsb(nn, level, sps)%w = zero
  flowDomsb(nn, level, sps)%dw = zero
  flowDomsb(nn, level, sps)%fw = zero
  
  flowDomsb(nn, level, sps)%p = zero
  flowDomsb(nn, level, sps)%gamma = zero
  
  flowDomsb(nn, level, sps)%rlv = zero
  flowDomsb(nn, level, sps)%rev = zero
  
  flowDomsb(nn, level, sps)%dtl  = zero
  flowDomsb(nn, level, sps)%radI = zero
  flowDomsb(nn, level, sps)%radJ = zero
  flowDomsb(nn, level, sps)%radK = zero
  
  bocoLoop: do mm=1,nBocos
     flowDomsb(nn, level, sps)%BCData(mm)%norm = zero
     flowDomsb(nn, level, sps)%bcData(mm)%rface = zero
     flowDomsb(nn, level, sps)%bcData(mm)%Fp = zero
     flowDomsb(nn, level, sps)%bcData(mm)%Fv = zero
     flowDomsb(nn, level, sps)%bcData(mm)%M = zero
     flowDomsb(nn, level, sps)%bcData(mm)%oArea = zero
     flowDomsb(nn, level, sps)%BCData(mm)%uSlip = zero
     flowDomsb(nn, level, sps)%BCData(mm)%TNS_Wall = zero
     flowDomsb(nn, level, sps)%BCData(mm)%Cavitation = zero
     flowDomsb(nn, level, sps)%BCData(mm)%SepSensor = zero
  end do bocoLoop
  
   viscbocoLoop: do mm=1,nviscBocos
      flowDomsb(nn, level, sps)%viscSubface(mm)%tau = zero 
      flowDomsb(nn, level, sps)%viscSubface(mm)%q = zero
   end  do viscbocoLoop

  end subroutine zeroADSeedsBwd

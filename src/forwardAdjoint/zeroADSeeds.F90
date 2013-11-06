! This function zeros all 'd' values of on the specifified block/level

subroutine zeroADSeeds(nn, level, sps)

  use blockPointers_d ! This modules includes blockPointers
  use inputTimeSpectral
  use flowVarRefState
  implicit none

  ! Input parameters
  integer(kind=intType) :: nn, level, sps

  ! Working parameters
  integer(kind=intType) :: mm

  flowDomsd(nn, level, sps)%x = zero
  flowDomsd(nn, level, sps)%si = zero
  flowDomsd(nn, level, sps)%sj = zero
  flowDomsd(nn, level, sps)%sk = zero
  flowDomsd(nn, level, sps)%vol = zero
  
  flowDomsd(nn, level, sps)%s = zero
  flowDomsd(nn, level, sps)%sFaceI = zero
  flowDomsd(nn, level, sps)%sFaceJ = zero
  flowDomsd(nn, level, sps)%sFaceK = zero
  
  flowDomsd(nn, level, sps)%w = zero
  flowDomsd(nn, level, sps)%dw = zero
  flowDomsd(nn, level, sps)%fw = zero
  
  flowDomsd(nn, level, sps)%p = zero
  flowDomsd(nn, level, sps)%gamma = zero
  
  flowDomsd(nn, level, sps)%rlv = zero
  flowDomsd(nn, level, sps)%rev = zero
  
  flowDomsd(nn, level, sps)%dtl  = zero
  flowDomsd(nn, level, sps)%radI = zero
  flowDomsd(nn, level, sps)%radJ = zero
  flowDomsd(nn, level, sps)%radK = zero
  
  bocoLoop: do mm=1,nBocos
     flowDomsd(nn, level, sps)%BCData(mm)%norm = zero
     flowDomsd(nn, level, sps)%bcData(mm)%rface = zero
     flowDomsd(nn, level, sps)%bcData(mm)%Fp = zero
     flowDomsd(nn, level, sps)%bcData(mm)%Fv = zero
     flowDomsd(nn, level, sps)%bcData(mm)%M = zero
     flowDomsd(nn, level, sps)%bcData(mm)%oArea = zero
     flowDomsd(nn, level, sps)%BCData(mm)%uSlip = zero
     flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall = zero
  end do bocoLoop
  
  if (viscous) then
     viscbocoLoop: do mm=1,nviscBocos
        flowDomsd(nn, level, sps)%viscSubface(mm)%tau = zero 
        flowDomsd(nn, level, sps)%viscSubface(mm)%q = zero
     end do viscbocoLoop
  end if

end subroutine zeroADSeeds

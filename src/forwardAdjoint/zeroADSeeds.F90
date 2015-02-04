! This function zeros all 'd' values of on the specifified block/level

subroutine zeroADSeeds(nn, level, sps)

  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  use bcroutines_b
  implicit none

  ! Input parameters
  integer(kind=intType) :: nn, level, sps

  ! Working parameters
  integer(kind=intType) :: mm

  flowDomsd(nn, level, sps)%d2wall = zero
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
  flowDomsd(nn, level, sps)%scratch = zero

  flowDomsd(nn, level, sps)%p = zero
  flowDomsd(nn, level, sps)%gamma = zero
  flowDomsd(nn, level, sps)%aa = zero
  
  flowDomsd(nn, level, sps)%rlv = zero
  flowDomsd(nn, level, sps)%rev = zero
  
  flowDomsd(nn, level, sps)%radI = zero
  flowDomsd(nn, level, sps)%radJ = zero
  flowDomsd(nn, level, sps)%radK = zero
  
  flowDomsd(nn, level, sps)%ux = zero
  flowDomsd(nn, level, sps)%uy = zero
  flowDomsd(nn, level, sps)%uz = zero
  flowDomsd(nn, level, sps)%vx = zero
  flowDomsd(nn, level, sps)%vy = zero
  flowDomsd(nn, level, sps)%vz = zero
  flowDomsd(nn, level, sps)%wx = zero
  flowDomsd(nn, level, sps)%wy = zero
  flowDomsd(nn, level, sps)%wz = zero
  flowDomsd(nn, level, sps)%qx = zero
  flowDomsd(nn, level, sps)%qy = zero
  flowDomsd(nn, level, sps)%qz = zero

  bocoLoop: do mm=1,nBocos
     flowDomsd(nn, level, sps)%BCData(mm)%norm = zero
     flowDomsd(nn, level, sps)%bcData(mm)%rface = zero
     flowDomsd(nn, level, sps)%bcData(mm)%Fp = zero
     flowDomsd(nn, level, sps)%bcData(mm)%Fv = zero
     flowDomsd(nn, level, sps)%bcData(mm)%M = zero
     flowDomsd(nn, level, sps)%bcData(mm)%oArea = zero
     flowDomsd(nn, level, sps)%BCData(mm)%uSlip = zero
     flowDomsd(nn, level, sps)%BCData(mm)%TNS_Wall = zero
     flowDomsd(nn, level, sps)%BCData(mm)%Cavitation = zero
     flowDomsd(nn, level, sps)%BCData(mm)%SepSensor = zero
  end do bocoLoop

  if (sps == 1) then
     flowDomsd(nn,1,sps)%bmti1 = zero
     flowDomsd(nn,1,sps)%bmti2 = zero
     flowDomsd(nn,1,sps)%bmtj1 = zero
     flowDomsd(nn,1,sps)%bmtj2 = zero
     flowDomsd(nn,1,sps)%bmtk1 = zero
     flowDomsd(nn,1,sps)%bmtk2 = zero
     flowDomsd(nn,1,sps)%bvti1 = zero
     flowDomsd(nn,1,sps)%bvti2 = zero
     flowDomsd(nn,1,sps)%bvtj1 = zero
     flowDomsd(nn,1,sps)%bvtj2 = zero
     flowDomsd(nn,1,sps)%bvtk1 = zero
     flowDomsd(nn,1,sps)%bvtk2 = zero
  end if

  viscbocoLoop: do mm=1,nviscBocos
     flowDomsd(nn, level, sps)%viscSubface(mm)%tau = zero 
     flowDomsd(nn, level, sps)%viscSubface(mm)%q = zero
  end do viscbocoLoop

  ! Now zero these
  ww0d = zero
  ww1d = zero
  ww2d = zero
  ww3d = zero

  pp0d = zero
  pp1d = zero
  pp2d = zero
  pp3d = zero

  rlv0d = zero
  rlv1d = zero
  rlv2d = zero
  rlv3d = zero

  rev0d = zero
  rev1d = zero
  rev2d = zero
  rev3d = zero
  ssid = zero
  xxd = zero

end subroutine zeroADSeeds

MODULE TAPENADE_DP_TEST

  IMPLICIT NONE
  REAL*8 fwdarray(1000), revarray(1000)
  CHARACTER*20 dppoints(0:1000)
  INTEGER :: indexinarray = 1
  INTEGER :: maxindexinarray = 0
  REAL*8 :: dotp=0.0
  INTEGER :: unitcallcount = 0
  INTEGER :: dptriggercount = 1
  REAL*8 :: dpdelta

CONTAINS

  subroutine DPFWDUPDATECOUNT(unit)
    CHARACTER unit*(*)
    unitcallcount = unitcallcount + 1
  end subroutine DPFWDUPDATECOUNT

  subroutine DPREVUPDATECOUNT(unit)
    CHARACTER unit*(*)
    unitcallcount = unitcallcount - 1
  end subroutine DPREVUPDATECOUNT

  LOGICAL FUNCTION DPFWDTESTCOUNT(unit, pos)
    CHARACTER unit*(*)
    INTEGER pos
    DPFWDTESTCOUNT = (unitcallcount.EQ.dptriggercount)
  end FUNCTION DPFWDTESTCOUNT

  LOGICAL FUNCTION DPREVTESTCOUNT(unit, pos)
    CHARACTER unit*(*)
    INTEGER pos
    DPREVTESTCOUNT = (unitcallcount.EQ.dptriggercount)
  end FUNCTION DPREVTESTCOUNT

  subroutine DPFWDREAL8(xd)
    REAL*8 xd
    call DPPUSHREAL8(xd)
    dotp = dotp + xd*xd
  end subroutine DPFWDREAL8

  subroutine DPFWDREAL8ARRAY(xd,sz)
    INTEGER sz,i
    REAL*8 xd(sz)
    call DPPUSHREAL8ARRAY(xd,sz)
    do i=1,sz
       dotp = dotp + xd(i)*xd(i)
    enddo
  end subroutine DPFWDREAL8ARRAY

  subroutine DPFWDREAL4(xd)
    REAL xd
    call DPPUSHREAL4(xd)
    dotp = dotp + xd*xd
  end subroutine DPFWDREAL4

  subroutine DPFWDREAL4ARRAY(xd,sz)
    INTEGER sz,i
    REAL xd(sz)
    call DPPUSHREAL4ARRAY(xd,sz)
    do i=1,sz
       dotp = dotp + xd(i)*xd(i)
    enddo
  end subroutine DPFWDREAL4ARRAY

  subroutine DPCUTFWDREAL8(xd)
    REAL*8 xd
    xd = 0.0
    call DPPUSHREAL8(xd)
  end subroutine DPCUTFWDREAL8

  subroutine DPCUTFWDREAL8ARRAY(xd,sz)
    INTEGER sz,i
    REAL*8 xd(sz)
    do i=1,sz
       xd(i) = 0.0
    enddo
    call DPPUSHREAL8ARRAY(xd,sz)
  end subroutine DPCUTFWDREAL8ARRAY

  subroutine DPREVREAL8(xb)
    REAL*8 xb, xd
    call DPPOPREAL8(xd)
    dotp = dotp + xd*xb
    xb = xd
  end subroutine DPREVREAL8

  subroutine DPREVREAL8ARRAY(xb,sz)
    INTEGER sz,i
    REAL*8 xb(sz), xd(100000)
    call DPPOPREAL8ARRAY(xd,sz)
    do i=1,sz
       dotp = dotp + xd(i)*xb(i)
       xb(i) = xd(i)
    enddo
  end subroutine DPREVREAL8ARRAY

  subroutine DPREVREAL4(xb)
    REAL xb, xd
    call DPPOPREAL4(xd)
    dotp = dotp + xd*xb
    xb = xd
  end subroutine DPREVREAL4

  subroutine DPREVREAL4ARRAY(xb,sz)
    INTEGER sz,i
    REAL xb(sz), xd(100000)
    call DPPOPREAL4ARRAY(xd,sz)
    do i=1,sz
       dotp = dotp + xd(i)*xb(i)
       xb(i) = xd(i)
    enddo
  end subroutine DPREVREAL4ARRAY

  subroutine DPFWDINITDISPLAY(point)
    CHARACTER point*(*)
    dotp = 0.0
    dppoints(0)=point
    print *, point,' call#',unitcallcount
  end subroutine DPFWDINITDISPLAY

  subroutine DPREVINITDISPLAY(point)
    CHARACTER point*(*)
    dotp = 0.0
    print *, point,' call#',unitcallcount
  end subroutine DPREVINITDISPLAY

  subroutine DPFWDDISPLAY(point)
    CHARACTER point*(*)
    print *, 'DotProduct = ', dotp
    dppoints(indexinarray)=point
    fwdarray(indexinarray)=dotp
    maxindexinarray = indexinarray
    indexinarray = indexinarray+1
    dotp = 0.0
    print *, point,' call#',unitcallcount
  end subroutine DPFWDDISPLAY

  subroutine DPREVDISPLAY(point)
    CHARACTER point*(*)
    print *, 'DotProduct = ', dotp
    indexinarray = indexinarray-1
    revarray(indexinarray)=dotp
    dotp = 0.0
    print *, point,' call#',unitcallcount
  end subroutine DPREVDISPLAY

  subroutine DPPRINTSUMMARY()
    INTEGER i
    print *, 'Dot Product test summary:'
    print *,dppoints(0)
    do i=1,maxindexinarray,1
       if (.not.(fwdarray(i).eq.revarray(i))) then
          print *,fwdarray(i)
          print *,revarray(i)
       else
          print *,fwdarray(i),'OK!'
       endif
       print *,dppoints(i)
    enddo
  end subroutine DPPRINTSUMMARY

  subroutine DPPRINTPARTIALSUM(point)
    CHARACTER point*(*)
    print *, 'Partial sum at ',point,': ',dotp
  end subroutine DPPRINTPARTIALSUM

  subroutine DPINITDELTA()
    dpdelta = 0.0_8
  end subroutine DPINITDELTA

  subroutine DPPRINTDELTA(point)
    CHARACTER point*(*)
    print *, 'Delta at ',point,': ',dotp-dpdelta
    dpdelta = dotp
  end subroutine DPPRINTDELTA

END MODULE TAPENADE_DP_TEST

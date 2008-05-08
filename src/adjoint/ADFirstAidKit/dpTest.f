C$Id: dpTest.f 367 2006-11-29 00:47:02Z acmarta $

      BLOCK DATA callcounter
      INTEGER unitcallcount,triggercount
      COMMON /callcount/unitcallcount,triggercount
      DATA unitcallcount/0/
      DATA triggercount/10/
      END

      subroutine DPFWDUPDATECOUNT(unit)
      CHARACTER unit*(*)
      INTEGER unitcallcount,triggercount
      COMMON /callcount/unitcallcount,triggercount
      unitcallcount = unitcallcount + 1
      end

      subroutine DPREVUPDATECOUNT(unit)
      CHARACTER unit*(*)
      INTEGER unitcallcount,triggercount
      COMMON /callcount/unitcallcount,triggercount
      unitcallcount = unitcallcount - 1
      end

      LOGICAL FUNCTION DPFWDTESTCOUNT(unit)
      CHARACTER unit*(*)
      INTEGER unitcallcount,triggercount
      COMMON /callcount/unitcallcount,triggercount
      DPFWDTESTCOUNT = (unitcallcount.EQ.triggercount)
      end

      LOGICAL FUNCTION DPREVTESTCOUNT(unit)
      CHARACTER unit*(*)
      INTEGER unitcallcount,triggercount
      COMMON /callcount/unitcallcount,triggercount
      DPREVTESTCOUNT = (unitcallcount.EQ.triggercount)
      end

      subroutine DPFWDREAL8(xd)
      REAL*8 xd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPUSHREAL8(xd)
      dotp = dotp + xd*xd
      end

      subroutine DPFWDREAL8ARRAY(xd,sz)
      INTEGER sz,i
      REAL*8 xd(sz)
      REAL*8 dotp
      common /cdotp/dotp
      call DPPUSHREAL8ARRAY(xd,sz)
      do i=1,sz
        dotp = dotp + xd(i)*xd(i)
      enddo
      end

      subroutine DPFWDREAL4(xd)
      REAL xd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPUSHREAL4(xd)
      dotp = dotp + xd*xd
      end

      subroutine DPFWDREAL4ARRAY(xd,sz)
      INTEGER sz,i
      REAL xd(sz)
      REAL*8 dotp
      common /cdotp/dotp
      call DPPUSHREAL4ARRAY(xd,sz)
      do i=1,sz
        dotp = dotp + xd(i)*xd(i)
      enddo
      end

      subroutine DPCUTFWDREAL8(xd)
      REAL*8 xd
      xd = 0.0
      call DPPUSHREAL8(xd)
      end

      subroutine DPCUTFWDREAL8ARRAY(xd,sz)
      INTEGER sz,i
      REAL*8 xd(sz)
      do i=1,sz
         xd(i) = 0.0
      enddo
      call DPPUSHREAL8ARRAY(xd,sz)
      end

      subroutine DPREVREAL8(xb)
      REAL*8 xb, xd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL8(xd)
      dotp = dotp + xd*xb
      xb = xd
      end

      subroutine DPREVREAL8ARRAY(xb,sz)
      INTEGER sz,i
      REAL*8 xb(sz), xd(100000)
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL8ARRAY(xd,sz)
      do i=1,sz
        dotp = dotp + xd(i)*xb(i)
        xb(i) = xd(i)
      enddo
      end

      subroutine DPREVREAL4(xb)
      REAL xb, xd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL4(xd)
      dotp = dotp + xd*xb
      xb = xd
      end

      subroutine DPREVREAL4ARRAY(xb,sz)
      INTEGER sz,i
      REAL xb(sz), xd(100000)
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL4ARRAY(xd,sz)
      do i=1,sz
        dotp = dotp + xd(i)*xb(i)
        xb(i) = xd(i)
      enddo
      end

      subroutine DPFWDINITDISPLAY(point)
      CHARACTER point*(*) 
      REAL*8 dotp
      common /cdotp/dotp
      dotp = 0.0
      print *, point
      end

      subroutine DPREVINITDISPLAY(point)
      CHARACTER point*(*) 
      REAL*8 dotp
      common /cdotp/dotp
      dotp = 0.0
      print *, point
      end

      subroutine DPFWDDISPLAY(point)
      CHARACTER point*(*) 
      REAL*8 dotp
      common /cdotp/dotp
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*10 dppoints(1000)
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints
      DATA indexinarray/1/
      DATA maxindexinarray/0/
      print *, 'DotProduct = ', dotp
      dppoints(indexinarray)=point
      fwdarray(indexinarray)=dotp
      maxindexinarray = indexinarray
      indexinarray = indexinarray+1
      dotp = 0.0
      print *, point
      end

      subroutine DPREVDISPLAY(point)
      CHARACTER point*(*) 
      REAL*8 dotp
      common /cdotp/dotp
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*10 dppoints(1000)
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints
      print *, 'DotProduct = ', dotp
      indexinarray = indexinarray-1
      revarray(indexinarray)=dotp
      dotp = 0.0
      print *, point
      end

      subroutine DPPRINTSUMMARY()
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*10 dppoints(1000)
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints
      INTEGER i
      print *, 'Dot Product test summary:'
      do i=1,maxindexinarray,1
         if (.not.(fwdarray(i).eq.revarray(i))) then
            print *,fwdarray(i)
            print *,revarray(i)
         endif
         print *, ' ====> ',dppoints(i)
      enddo
      end

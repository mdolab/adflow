C$Id: dpTest.f 2595 2008-10-15 09:04:22Z llh $

      BLOCK DATA DPWORKZONE
      REAL *8 arrayd(1000000)
      common /dpwz/arrayd
      END

      BLOCK DATA DPCALLCOUNTER
      INTEGER dpunitcallcount,dptriggercount
      COMMON /DPCALLCOUNT/dpunitcallcount,dptriggercount
      DATA dpunitcallcount/0/
      DATA dptriggercount/3/
      END

      BLOCK DATA DPSUMM
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
      DATA indexinarray/1/
      DATA maxindexinarray/0/
      END

      subroutine DPFWDUPDATECOUNT(unit)
      CHARACTER unit*(*)
      INTEGER dpunitcallcount,dptriggercount
      COMMON /DPCALLCOUNT/dpunitcallcount,dptriggercount
      dpunitcallcount = dpunitcallcount + 1
      end

      subroutine DPREVUPDATECOUNT(unit)
      CHARACTER unit*(*)
      INTEGER dpunitcallcount,dptriggercount
      COMMON /DPCALLCOUNT/dpunitcallcount,dptriggercount
      dpunitcallcount = dpunitcallcount - 1
      end

      LOGICAL FUNCTION DPFWDTESTCOUNT(unit, place)
      CHARACTER unit*(*)
      INTEGER place
      INTEGER dpunitcallcount,dptriggercount
      COMMON /DPCALLCOUNT/dpunitcallcount,dptriggercount
      DPFWDTESTCOUNT = (dpunitcallcount.EQ.dptriggercount)
      end

      LOGICAL FUNCTION DPREVTESTCOUNT(unit, place)
      CHARACTER unit*(*)
      INTEGER place
      INTEGER dpunitcallcount,dptriggercount
      COMMON /DPCALLCOUNT/dpunitcallcount,dptriggercount
      DPREVTESTCOUNT = (dpunitcallcount.EQ.dptriggercount)
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
      dotp = dotp - xd*xb
      xb = xd
      end

      subroutine DPREVREAL8ARRAY(xb,sz)
      INTEGER sz,i
      REAL*8 xb(sz)
      REAL *8 arrayd(1000000)
      common /dpwz/arrayd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL8ARRAY(arrayd,sz)
      do i=sz,1,-1
        dotp = dotp - arrayd(i)*xb(i)
        xb(i) = arrayd(i)
      enddo
      end

      subroutine DPREVREAL4(xb)
      REAL xb, xd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL4(xd)
      dotp = dotp - xd*xb
      xb = xd
      end

      subroutine DPREVREAL4ARRAY(xb,sz)
      INTEGER sz,i
      REAL xb(sz)
      REAL*4 arrayd(2000000)
      common /dpwz/arrayd
      REAL*8 dotp
      common /cdotp/dotp
      call DPPOPREAL4ARRAY(arrayd,sz)
      do i=sz,1,-1
        dotp = dotp - arrayd(i)*xb(i)
        xb(i) = arrayd(i)
      enddo
      end

      subroutine DPFWDINITDISPLAY(point)
      CHARACTER point*(*)
      REAL*8 dotp
      common /cdotp/dotp
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
      initdppoint=point
      dotp = 0.0
      print *, point
      end

      subroutine DPREVINITDISPLAY(point)
      CHARACTER point*(*)
      REAL*8 dotp
      common /cdotp/dotp
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
c      dotp = 0.0
      dotp = fwdarray(indexinarray-1)
      print *, point
      end

      subroutine DPSONDE(point)
      CHARACTER point*(*)
      REAL*8 dotp
      common /cdotp/dotp
      write (6,1002) 'sonde DotProduct ', point, ': ', dotp
 1002 format(a,a,a,e26.20)
      end

      subroutine DPFWDDISPLAY(point)
      CHARACTER point*(*)
      REAL*8 dotp
      common /cdotp/dotp
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
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
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
      print *, 'DotProduct error = ', dotp
      indexinarray = indexinarray-1
      revarray(indexinarray)=dotp
c      dotp = 0.0
      if (indexinarray.gt.1) then
         dotp = fwdarray(indexinarray-1)
      else
         dotp = 0.0
      endif
      print *, point
      end

      subroutine DPPRINTSUMMARY()
      REAL*8 fwdarray(1000), revarray(1000)
      CHARACTER*20 dppoints(1000),initdppoint
      INTEGER indexinarray,maxindexinarray
      COMMON /dpsummary/ fwdarray,revarray,
     +     indexinarray,maxindexinarray,dppoints,initdppoint
      INTEGER i
      write (6,*) ''
      write (6,*) '-------------------------'
      write (6,*) 'Dot Product test summary:'
      write (6,*) ' ====> ',initdppoint

      do i=1,maxindexinarray,1
         if (SERIOUSERROR(fwdarray(i),revarray(i)).gt.-35) then
            write (6,1001) '  INVALID Dot product: ',
     +            fwdarray(i),'  err:',revarray(i)
         else
            write (6,1001) '    valid fragment:    ',
     +            fwdarray(i),'  err:',revarray(i)
         endif
         write (6,*) ' ====> ',dppoints(i)
      enddo
 1001 format(a,e26.20,a,e26.20)
      end

c Measures if this "error" on this "value"
c  may be significant wrt machine precision.
c If returns < -35, then error may be only due to machine precision.
      REAL*8 function SERIOUSERROR(value, delta)
      REAL*8 value, delta, logvalue, logdelta
      if (delta.eq.0.0) then
         seriouserror = 0.0
      else
         logdelta = LOG(ABS(delta))
         if (value.eq.0.0) then
            logvalue = 0.0
         else
            logvalue = LOG(ABS(value))
         endif
         seriouserror = logdelta-logvalue
      endif
c      print *,'seriousness:',seriouserror
      return
      end


      subroutine DPINITREAL8(xd)
      REAL*8 xd
      xd = 0.0
      end

      subroutine DPINITREAL8ARRAY(xd,sz)
      INTEGER sz,i
      REAL*8 xd(sz)
      do i=1,sz
        xd(i) = 0.0
      enddo
      end

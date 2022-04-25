C
C $Id: ddTest.f $
C

c Utility subroutines for validation of tangent AD derivatives
c compared to Divided Differences.

      subroutine DDFWDREAL4(varname, var, vard)
      IMPLICIT NONE
      real*4 var, vard, ddvar, dd, diff, varwr
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname
      integer ddphase, ddfile
      LOGICAL areNaNs
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         WRITE(ddfile, *) var
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         call DDPICKTWO4(var, varwr, ddfile, ddvar, areNaNs)
         if (.not.areNaNs) then
            dd = (ddvar-varwr)/ddeps
            if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
               diff = (abs(vard-dd)*100.0)/
     +              max(abs(vard),abs(dd))
               if (diff.gt.1.0) then
                  diffstr = 'DIFFERENCE!!'
!                  print *,'DIFF ',var,varwr,ddvar,ddeps                    !!Debug level 2
               else
                  diffstr = '            '
               endif
               write (*, 300) varname, vard, dd, diffstr
            endif
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
      end

      subroutine DDFWDREAL4ARRAY
     +     (varname, var, vard, length)
      IMPLICIT NONE
      integer length
      real*4 var(length)
      real*4 vard(length)
      real*4 ddvar, dd, diff, varwr
      real*4 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader
      integer ddphase, ddfile
      LOGICAL areNaNs
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         do i1=1,length
            WRITE(ddfile, *) var(i1)
         enddo
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         notprintedheader=.true.
         ibuf = 1
         do i1=1,length
            call DDPICKTWO4(var(i1), varwr, ddfile, ddvar, areNaNs)
            if (.not.areNaNs) then
               dd = (ddvar-varwr)/ddeps
               if ((abs(vard(i1)).gt.epszero)
     +              .or.(abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  diff = (abs(vard(i1)-dd)*100.0)/
     +                 max(abs(vard(i1)),abs(dd))
                  if (diff.gt.1.0) then
                     diffbuf(ibuf) = '  DIFFERENCE!!'
!                     print *,'DIFF ',var(i1),varwr,ddvar,ddeps               !!Debug level 2
                     ibuf = ibuf+1
                  else
c                     diffbuf(ibuf) = '              '
c                     ibuf = ibuf+1
                  endif
               endif
            endif
            if(ibuf.gt.10.or.(i1.eq.length.and.ibuf.gt.1)) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)
     +              (indexbuf1(j),valbuf(j), j=1,ibuf-1)
               write (*, 250) (ddbuf(j), j=1,ibuf-1)
               write (*, 270) (diffbuf(j), j=1,ibuf-1)
               ibuf = 1
            endif
         end do
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end

      subroutine DDPICKTWO4(var, varwr, ddfile, ddvar, areNaNs)
      IMPLICIT NONE
      REAL*4 var,varwr,ddvar
      LOGICAL areNaNs
      INTEGER ddfile,stat1,stat2

      OPEN(38, FILE='ddwrfile')
      WRITE(38, *) var
      REWIND(38)
      READ(38, *,IOSTAT=stat1) varwr
      CLOSE(38)
      READ(ddfile, *,IOSTAT=stat2) ddvar
      areNaNs = stat1.eq.225.and.stat2.eq.225
      end



      subroutine DDFWDREAL8(varname, var, vard)
      IMPLICIT NONE
      real*8 var, vard, ddvar, dd, diff, varwr
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname
      integer ddphase, ddfile
      LOGICAL areNaNs
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         WRITE(ddfile, *) var
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         call DDPICKTWO8(var, varwr, ddfile, ddvar, areNaNs)
         if (.not.areNaNs) then
            dd = (ddvar-varwr)/ddeps
            if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
               diff = (abs(vard-dd)*100.0)/
     +              max(abs(vard),abs(dd))
!               print *,'DIFF ',var,varwr,ddvar,ddeps !!Debug level 2
               if (diff.gt.1.0) then
                  diffstr = 'DIFFERENCE!!'
               else
                  diffstr = '            '
               endif
               write (*, 300) varname, vard, dd, diffstr
            endif
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
      end

      subroutine DDFWDREAL8ARRAY
     +     (varname, var, vard, length)
      IMPLICIT NONE
      integer length
      real*8 var(length)
      real*8 vard(length)
      real*8 ddvar, dd, diff, varwr
      real*8 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader
      integer ddphase, ddfile
      LOGICAL areNaNs
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         do i1=1,length
            WRITE(ddfile, *) var(i1)
         enddo
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         notprintedheader=.true.
         ibuf = 1
         do i1=1,length
            call DDPICKTWO8(var(i1), varwr, ddfile, ddvar, areNaNs)
            if (.not.areNaNs) then
               dd = (ddvar-varwr)/ddeps
               if ((abs(vard(i1)).gt.epszero)
     +              .or.(abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  diff = (abs(vard(i1)-dd)*100.0)/
     +                 max(abs(vard(i1)),abs(dd))
                  if (diff.gt.1.0) then
                     diffbuf(ibuf) = '  DIFFERENCE!!'
                     ibuf = ibuf+1
                  else
c                     diffbuf(ibuf) = '              '
c                     ibuf = ibuf+1
                  endif
               endif
            endif
            if(ibuf.gt.10.or.(i1.eq.length.and.ibuf.gt.1)) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)
     +              (indexbuf1(j),valbuf(j), j=1,ibuf-1)
               write (*, 250) (ddbuf(j), j=1,ibuf-1)
               write (*, 270) (diffbuf(j), j=1,ibuf-1)
               ibuf = 1
            endif
         end do
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end

      subroutine DDFWDREAL8ARRAYREPLACE
     +     (varname, var, vard, length)
      IMPLICIT NONE
      integer length
      real*8 var(length)
      real*8 vard(length)
      real*8 ddvar, dd, diff, varwr
      real*8 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader
      integer ddphase, ddfile
      LOGICAL areNaNs
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         do i1=1,length
            WRITE(ddfile, *) var(i1)
         enddo
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         notprintedheader=.true.
         ibuf = 1
         do i1=1,length
            call DDPICKTWO8(var(i1), varwr, ddfile, ddvar, areNaNs)
            if (.not.areNaNs) then
               dd = (ddvar-varwr)/ddeps
               if ((abs(vard(i1)).gt.epszero)
     +              .or.(abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  diff = (abs(vard(i1)-dd)*100.0)/
     +                 max(abs(vard(i1)),abs(dd))
                  if (diff.gt.1.0) then
                     vard(i1) = dd
                  endif
               endif
            endif
            if(ibuf.gt.10.or.(i1.eq.length.and.ibuf.gt.1)) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)
     +              (indexbuf1(j),valbuf(j), j=1,ibuf-1)
               write (*, 250) (ddbuf(j), j=1,ibuf-1)
               write (*, 270) (diffbuf(j), j=1,ibuf-1)
               ibuf = 1
            endif
         end do
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end

      subroutine DDPICKTWO8(var, varwr, ddfile, ddvar, areNaNs)
      IMPLICIT NONE
      REAL*8 var,varwr,ddvar
      LOGICAL areNaNs
      INTEGER ddfile,stat1,stat2

      OPEN(38, FILE='ddwrfile')
      WRITE(38, *) var
      REWIND(38)
      READ(38, *,IOSTAT=stat1) varwr
      CLOSE(38)
      READ(ddfile, *,IOSTAT=stat2) ddvar
      areNaNs = stat1.eq.225.and.stat2.eq.225
      end

      subroutine DDCHECKVARNAME(varname)
      IMPLICIT NONE
      character varname*(*)
      character*50 ddvarname
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero
      integer linesskip

      linesskip = 0
 100  CONTINUE
      if (linesskip.GT.200) THEN
         write(*,*)
     +        'ERROR: Too many lines skipped. Bad DD program control ?'
         write(*,*) 'Was looking for variable:',varname
         STOP
      ENDIF
      READ(ddfile, '(a)') ddvarname
      if (ddvarname.ne.varname) then
c         write(*,*) 'ERROR: mismatch in DD program control !!!',
c     +        ' read ', ddvarname, ' expecting ', varname
         linesskip = linesskip+1
         GOTO 100
      endif
      end

      BLOCK DATA DDFWDCALLCOUNTER
      INTEGER sz, ddnum
      PARAMETER (sz=99)
      CHARACTER*33 ddunames(sz), ddcontainernames(sz)
      INTEGER dducounts(sz) ,ddufroms(sz) ,ddutos(sz)
      LOGICAL dduactives(sz)
      COMMON /DDCALLCOUNT/ddnum,ddunames,ddcontainernames,
     +     dducounts,ddufroms,ddutos,dduactives

      DATA ddnum/0/
      END

      subroutine DDFWDADDCHECKEDUNIT(unitname, containername, from, to)
      INTEGER sz, ddnum
      PARAMETER (sz=99)
      CHARACTER*33 ddunames(sz), ddcontainernames(sz)
      INTEGER dducounts(sz) ,ddufroms(sz) ,ddutos(sz)
      LOGICAL dduactives(sz)
      COMMON /DDCALLCOUNT/ddnum,ddunames,ddcontainernames,
     +     dducounts,ddufroms,ddutos,dduactives

      CHARACTER unitname*(*)
      CHARACTER containername*(*)
      INTEGER from, to

      ddnum = ddnum+1
      ddunames(ddnum) = unitname
      ddcontainernames(ddnum) = containername
      ddufroms(ddnum) = from
      ddutos(ddnum) = to
      dducounts(ddnum) = 0
      dduactives(ddnum) = .FALSE.
      end

      subroutine DDFWDUPDATECOUNT(unitname)
      INTEGER sz, ddnum
      PARAMETER (sz=99)
      CHARACTER*33 ddunames(sz), ddcontainernames(sz)
      INTEGER dducounts(sz) ,ddufroms(sz) ,ddutos(sz)
      LOGICAL dduactives(sz)
      COMMON /DDCALLCOUNT/ddnum,ddunames,ddcontainernames,
     +     dducounts,ddufroms,ddutos,dduactives

      CHARACTER unitname*(*)
      CHARACTER*33 containername
      INTEGER i, index
      LOGICAL DDFWDTESTCOUNT
      index = -1
      do i=1,ddnum
         if (ddunames(i)==unitname) then
            index = i
            GOTO 100
         endif
      enddo
 100  if (index.ne.-1) then
         containername = ddcontainernames(index)
         if (containername.eq.""
     +        .OR.DDFWDTESTCOUNT(containername,0)) then
            dducounts(index) = dducounts(index)+1
         endif
      endif
      end

C   place: -1 => entry ; 0 => somewhere inside ; 1 => exit
      LOGICAL function DDFWDTESTCOUNT(unitname, place)
      INTEGER sz, ddnum
      PARAMETER (sz=99)
      CHARACTER*33 ddunames(sz), ddcontainernames(sz)
      INTEGER dducounts(sz) ,ddufroms(sz) ,ddutos(sz)
      LOGICAL dduactives(sz)
      COMMON /DDCALLCOUNT/ddnum,ddunames,ddcontainernames,
     +     dducounts,ddufroms,ddutos,dduactives

      CHARACTER unitname*(*)
      CHARACTER*33 containername
      INTEGER i, index, place, index2
      LOGICAL containerisactive
      DDFWDTESTCOUNT = .FALSE.
      index = -1
      do i=1,ddnum
         if (ddunames(i)==unitname) then
            index = i
            GOTO 100
         endif
      enddo
 100  if (index.ne.-1) then
         containername = ddcontainernames(index)
         if (containername.eq."") then
            containerisactive = .TRUE.
         else
            index2 = -1
            do i=1,ddnum
               if (ddunames(i)==containername) then
                  index2 = i
                  GOTO 200
               endif
            enddo
 200        containerisactive = (index2.ne.-1).AND.
     +           dduactives(index2)
         endif
         if (containerisactive) then
            if (place==-1) then
               if (dducounts(index).ge.ddufroms(index)
     +              .AND.(ddutos(index).le.0
     +                    .OR.dducounts(index).le.ddutos(index))) then
                  dduactives(index) = .TRUE.
                  DDFWDTESTCOUNT = .TRUE.
                  write(6,777) '>>> AT ENTRY INTO ',unitname,
     +                 ' #',dducounts(index)
               else
                  dduactives(index) = .FALSE.
                  DDFWDTESTCOUNT = .FALSE.
               endif
            else if (place==1) then
               DDFWDTESTCOUNT = dduactives(index)
               if (DDFWDTESTCOUNT) then
                  write(6,777) ' <<< AT EXIT FROM ',unitname,
     +                 ' #',dducounts(index)
               endif
               dduactives(index) = .FALSE.
            else
               DDFWDTESTCOUNT = dduactives(index)
            endif
         else
            dduactives(index) = .FALSE.
            DDFWDTESTCOUNT = .FALSE.
         endif
      endif
 777  format(a,a15,a,i4)
      end

PROGRAM example2
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-09  Time: 23:45:36

USE Kernel_Regression
IMPLICIT NONE

INTEGER, PARAMETER  :: n = 250, m = 300, ihom = 0, nue = 0
INTEGER, PARAMETER  :: irnd = 0, m1 = 400
INTEGER    :: ismo = 0, i, j, kord = nue+2
REAL (dp)  :: t(n), x(n), tt(m), b
REAL (dp)  :: yg(m), yl(m), ban(m)

REAL (dp)  :: tl = 1.0_dp, tu = 0.0_dp, sig = -1.0_dp,  &
              s(0:n) = (/ (0.0_dp,i=0,n-1), -1.0_dp /)

OPEN (81, FILE='raw.dat')
DO  i = 1, n
  READ (81,*) t(i), x(i)
END DO
CLOSE (81)
DO  j = 1, m
  tt(j) = DBLE(j-1) * (t(n)-t(1)) / DBLE(m-1) + t(1)
END DO

CALL glkern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,b,yg)
CALL lokern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,ban,yl)

OPEN (82,FILE='ex2.out')
WRITE(82, *) '      tt        Smoothed Y    Smoothed Y      Local'
WRITE(82, *) '                 (glkern)      (lokern)     Bandwidth'
WRITE (82,'(4G14.6)') (tt(j),yg(j),yl(j),ban(j),j = 1,m)
CLOSE (82)

STOP
END PROGRAM example2

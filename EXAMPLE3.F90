PROGRAM example3
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-09  Time: 23:45:41

USE Kernel_Regression
IMPLICIT NONE

INTEGER, PARAMETER  :: n = 250, m = 300, ihom = 1, nue = 0
INTEGER, PARAMETER  :: irnd = 0, m1 = 400
INTEGER    :: ismo = 1, i, j, kord = nue+2
REAL (dp)  :: t(n), x(n), tt(m)
REAL (dp)  :: y(m), ban(m)

REAL (dp)  :: tl = 1.0_dp, tu = 0.0_dp, sig = -1.0_dp,  &
              s(0:n) = (/ (0.0_dp,i=0,n-1), -1.0_dp /)

OPEN (81, FILE='raw.dat')
DO  i = 1, n
  READ (81,*) t(i), x(i)
END DO
CLOSE (81)

OPEN (82, FILE='band.dat')
DO  j = 1, m
  READ (82,*) tt(j), ban(j)
END DO
CLOSE (82)

CALL lokern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,ban,y)

OPEN (82, FILE='ex3.out')
WRITE (82, '(2G14.6)') (tt(j),y(j),j = 1,m)
CLOSE (82)
STOP
END PROGRAM example3

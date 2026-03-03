PROGRAM example1
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-09  Time: 23:45:31

USE Kernel_Regression
IMPLICIT NONE

INTEGER    :: nue, kord, i, ismo = 0, j
INTEGER, PARAMETER  :: n = 75, m = 300, ihom = 0, irnd = 1, m1 = 400
REAL (dp)  :: t(n), tt(m), tl = 1.0_dp, tu = 0.0_dp, sig = -1.0_dp
REAL (dp)  :: b0, b1, b2, y0(m), y1(m), y2(m)
REAL (dp)  :: s(0:n) = (/ (0.0_dp, i=0,n-1), -1.0_dp /)

REAL (dp)  :: x(n) =  (/  &
    1.9666, 1.9000, 1.6449, 1.4275, 2.4000, 1.8487, 1.3383,  &
    1.9514, 1.5722, 1.1978, 1.3109, 1.6940, 1.0816, 1.1070,  &
    1.2650, 1.1143, 0.5717, 0.9492, 0.3116, 0.7942, 0.5670, 0.3555,  &
    0.3243, 0.3242, 0.7647, 0.4167, 0.7583, 0.9462, 1.4337, 1.8820,  &
    2.3054, 2.1040, 3.8500, 3.7052, 4.3047, 4.4505, 4.3425,  &
    4.3516, 3.7578, 4.0112, 3.3286, 2.9024, 2.4331, 1.5544,  &
    1.1867, 0.6870, -0.3630, 0.0436, -0.8197, -0.2883, -1.0559, -0.9795,  &
    -1.7802, -2.2206, -1.3922, -1.6511, -1.7694, -1.1309, -2.1912,  &
    -1.5785, -2.6189, -1.8125, -2.6155, -1.4585, -2.0951, -2.1428,  &
    -2.4827, -2.4171, -2.6610, -2.6509, -2.6475, -3.1040, -2.7631,  &
    -2.9486, -3.1323 /)

OPEN (81, FILE='ex1.dat')

DO  i = 1, n
  t(i) = (i-0.5) / DBLE(n)
  WRITE (81, '(2F9.5)') t(i), x(i)
END DO
CLOSE (81)
DO  j = 1, m
  tt(j) = DBLE(j-1) * (t(n)-t(1)) / DBLE(m-1) + t(1)
END DO

nue = 0
kord = nue + 2
CALL glkern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,b0,y0)
nue = 1
kord = nue + 2
s(0) = 1.0
s(n) = 0.0
CALL glkern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,b1,y1)
nue = 2
kord = nue + 2
s(0) = 1.0
s(n) = 0.0
CALL glkern(t,x,n,tt,m,ihom,nue,kord,irnd,ismo,m1,tl,tu,s,sig,b2,y2)
WRITE(*, *) ' Plugin bandwidths used in glkern:'
WRITE(*, *) ' regression function:', b0
WRITE(*, *) ' first derivative:   ', b1
WRITE(*, *) ' second derivative:  ', b2

OPEN (82,FILE='ex1.out')
WRITE (82,'(4G14.6)') (tt(j),y0(j),y1(j),y2(j),j = 1,m)
CLOSE (82)
STOP
END PROGRAM example1

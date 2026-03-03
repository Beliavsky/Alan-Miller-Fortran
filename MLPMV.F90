MODULE lpmv_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
 
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp

CONTAINS
 

SUBROUTINE lpmv(v, m, x, pmv)

!    =======================================================
!    Purpose: Compute the associated Legendre function
!             Pmv(x) with an integer order and an arbitrary
!             nonnegative degree v
!    Input :  x   --- Argument of Pm(x)  ( -1 ó x ó 1 )
!             m   --- Order of Pmv(x)
!             v   --- Degree of Pmv(x)
!    Output:  PMV --- Pmv(x)
!    Routine called:  PSI for computing Psi function
!    =======================================================

REAL (dp), INTENT(IN)   :: v
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: pmv

REAL (dp)  :: c0, eps, pa, pss, psv, pv0, qr, r, r0, r1, r2, rg,  &
              s, s0, s1, s2, v0, vs, xq
INTEGER    :: j, k, nv

eps = 1.0D-14
nv = v
v0 = v - nv
IF (x == -1.0D0 .AND. v /= nv) THEN
  IF (m == 0) pmv = -1.0D+300
  IF (m /= 0) pmv = 1.0D+300
  RETURN
END IF
c0 = 1.0D0
IF (m /= 0) THEN
  rg = v * (v+m)
  DO  j = 1, m - 1
    rg = rg * (v*v-j*j)
  END DO
  xq = SQRT(1.0D0-x*x)
  r0 = 1.0D0
  DO  j = 1, m
    r0 = .5D0 * r0 * xq / j
  END DO
  c0 = r0 * rg
END IF
IF (v0 == 0.0D0) THEN
  pmv = 1.0D0
  r = 1.0D0
  DO  k = 1, nv - m
    r = 0.5D0 * r * (-nv+m+k-1.0D0) * (nv+m+k) / (k*(k+m)) * (1.0D0+x)
    pmv = pmv + r
  END DO
  pmv = (-1) ** nv * c0 * pmv
ELSE
  IF (x >= -0.35D0) THEN
    pmv = 1.0D0
    r = 1.0D0
    DO  k = 1, 100
      r = 0.5D0 * r * (-v+m+k-1.0D0) * (v+m+k) / (k*(m+k)) * (1.0D0-x)
      pmv = pmv + r
      IF (k > 12 .AND. ABS(r/pmv) < eps) EXIT
    END DO
    pmv = (-1) ** m * c0 * pmv
  ELSE
    vs = SIN(v*pi) / pi
    pv0 = 0.0D0
    IF (m /= 0) THEN
      qr = SQRT((1.0D0-x)/(1.0D0+x))
      r2 = 1.0D0
      DO  j = 1, m
        r2 = r2 * qr * j
      END DO
      s0 = 1.0D0
      r1 = 1.0D0
      DO  k = 1, m - 1
        r1 = 0.5D0 * r1 * (-v+k-1) * (v+k) / (k*(k-m)) * (1.0D0+x)
        s0 = s0 + r1
      END DO
      pv0 = -vs * r2 / m * s0
    END IF
    CALL psi(v, psv)
    pa = 2.0D0 * (psv+el) + pi / TAN(pi*v) + 1.0D0 / v
    s1 = 0.0D0
    DO  j = 1, m
      s1 = s1 + (j*j+v*v) / (j*(j*j-v*v))
    END DO
    pmv = pa + s1 - 1.0D0 / (m-v) + LOG(0.5D0*(1.0D0+x))
    r = 1.0D0
    DO  k = 1, 100
      r = 0.5D0 * r * (-v+m+k-1.0D0) * (v+m+k) / (k*(k+m)) * (1.0D0+x)
      s = 0.0D0
      DO  j = 1, m
        s = s + ((k+j)**2+v*v) / ((k+j)*((k+j)**2 - v*v))
      END DO
      s2 = 0.0D0
      DO  j = 1, k
        s2 = s2 + 1.0D0 / (j*(j*j - v*v))
      END DO
      pss = pa + s + 2.0D0*v*v*s2 - 1.0D0/(m+k-v) + LOG(0.5D0*(1.0D0+x))
      r2 = pss * r
      pmv = pmv + r2
      IF (ABS(r2/pmv) < eps) EXIT
    END DO
    pmv = pv0 + pmv * vs * c0
  END IF
END IF
RETURN
END SUBROUTINE lpmv


SUBROUTINE psi(x, ps)

!    ======================================
!    Purpose: Compute psi function
!    Input :  x  --- Argument of psi(x)
!    Output:  PS --- psi(x)
!    ======================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ps

REAL (dp), PARAMETER  ::  a1 = -.8333333333333D-01, a2 = .83333333333333333D-02, &
      a3 = -.39682539682539683D-02, a4 = .41666666666666667D-02,  &
      a5 = -.75757575757575758D-02, a6 = .21092796092796093D-01,  &
      a7 = -.83333333333333333D-01, a8 = .4432598039215686_dp
REAL (dp)  :: s, x2, xa
INTEGER    :: k, n

xa = ABS(x)
s = 0.0D0
IF (x == INT(x) .AND. x <= 0.0) THEN
  ps = 1.0D+300
  RETURN
ELSE IF (xa == INT(xa)) THEN
  n = xa
  DO  k = 1, n - 1
    s = s + 1.0_dp / k
  END DO
  ps = -el + s
ELSE IF (xa+0.5 == INT(xa+0.5)) THEN
  n = xa - .5
  DO  k = 1, n
    s = s + 1.0_dp / (2*k-1)
  END DO
  ps = -el + 2.0D0 * s - 1.386294361119891_dp
ELSE
  IF (xa < 10.0) THEN
    n = 10 - INT(xa)
    DO  k = 0, n - 1
      s = s + 1.0_dp / (xa+k)
    END DO
    xa = xa + n
  END IF
  x2 = 1.0D0 / (xa*xa)
  ps = LOG(xa) - .5D0 / xa + x2 * (((((((a8*x2 + a7)*x2 + a6)*x2 + a5)*  &
                 x2 + a4)*x2 + a3)*x2 + a2)*x2 + a1)
  ps = ps - s
END IF
IF (x < 0.0) ps = ps - pi * COS(pi*x) / SIN(pi*x) - 1.0D0 / x
RETURN
END SUBROUTINE psi
 
END MODULE lpmv_func
 
 
 
PROGRAM mlpmv
USE lpmv_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    =========================================================
!    Purpose: This program computes the associated Legendre
!             function Pmv(x) with an integer order and an
!             arbitrary nonnegative degree using subroutine
!             LPMV
!    Input :  x   --- Argument of Pm(x)  ( -1 ó x ó 1 )
!             m   --- Order of Pmv(x)
!             v   --- Degree of Pmv(x)
!    Output:  PMV --- Pmv(x)
!    Example:    m = 4,  x = 0.5
!                 v          Pmv(x)
!              -----------------------
!                1.5       .46218726
!                1.6       .48103143
!                1.7       .45031429
!                1.8       .36216902
!                1.9       .21206446
!                2.0       .00000000
!                2.5     -1.51996235
!    =========================================================

REAL (dp)  :: pmv, v, x
INTEGER    :: m

WRITE (*,*) 'Please enter m,v,x = ? '
READ (*,*) m, v, x
WRITE (*,5100) m, x
WRITE (*,*)
WRITE (*,*) '     v        Pmv(x)'
WRITE (*,*) '  -----------------------'
CALL lpmv(v, m, x, pmv)
WRITE (*,5000) v, pmv
STOP

5000 FORMAT (t4, f5.1, g16.8)
5100 FORMAT (t4, ' m =', i2, ',    x =', f6.2)
END PROGRAM mlpmv

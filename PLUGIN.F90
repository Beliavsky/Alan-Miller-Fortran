SUBROUTINE plugin(x,n,z,m,f,h)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-11  Time: 12:41:50
 
!-----------------------------------------------------------------------
!       Version: 1995
!       fortran77

!       Purpose:

!       Simple  Subroutine for kernel density estimation
!       with iterative plug-in bandwidth selection

!       This version only uses the gauss kernel and estimates only
!       the density itself and not its derivatives.

!  INPUT    X(N)   DOUBLE PRECISION   sorted data
!  INPUT    N      INTEGER            length of  X
!  INPUT    Z(M)   DOUBLE PRECISION   output grid (sorted)
!  INPUT    M      INTEGER            length of Z
!  OUTPUT   F(M)   DOUBLE PRECISION   estimated density
!  OUTPUT   H      DOUBLE PRECISION   estimated iterative plugin bandwidth

!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: z(:)
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(OUT)  :: f(:)
REAL (dp), INTENT(OUT)  :: h

! Local variables
INTEGER    :: i, it, iter, j, jbegin, jend
REAL (dp)  :: a, co1, const, d2, d3, e2, e3, h2, h3, pi, r2, rhat2, rhat3, &
              rt2pi, s, s2, s3, t, xiqr, xn

xn = n
xiqr = x(nint(.75*xn)) - x(nint(.25*xn)+1)
pi = 3.14159265358979324_dp
rt2pi = SQRT(2._dp*pi)
iter = 5
!-------  Estimate inflation constant C
h2 = (.920*xiqr) / (xn**(1._dp/7._dp))
h3 = (.912*xiqr) / (xn**(1._dp/9._dp))
s2 = 0._dp
s3 = 0._dp
loop20:  DO  i = 1, n - 1
  DO  j = i + 1, n
    d2 = ((x(i)-x(j))/h2) ** 2
    d3 = ((x(i)-x(j))/h3) ** 2
    IF (d2 > 50 .AND. d3 > 60) CYCLE loop20
    e2 = EXP(-0.5_dp*d2)
    e3 = EXP(-0.5_dp*d3)
    s2 = s2 + (d2**2 - 6*d2 + 3._dp) * e2
    s3 = s3 + (d3**3 - 15*d3**2 + 45*d3 - 15) * e3
  END DO
END DO loop20
rhat2 = (2.0_dp*s2) / ((xn**2)*(h2**5)*rt2pi)
rhat2 = rhat2 + 3._dp / (rt2pi*xn*(h2**5))
rhat3 = (-2.0_dp*s3) / ((xn**2)*(h3**7)*rt2pi)
rhat3 = rhat3 + 15.0_dp / (rt2pi*xn*(h3**7))
co1 = 1.357_dp * (rhat2/rhat3) ** (1.0_dp/7.0_dp)
!-
!-------  Compute constant of asymptotic formula
!-
const = 1._dp / (2._dp*SQRT(pi))
a = 1.132795764 / rhat3 ** (1._dp/7._dp) * xn ** (-1.0_dp/2.0_dp)
!-
!------  Loop over iterations
!-
DO  it = 1, iter
!-
!-------  Estimate functional
  
  s = 0._dp
  loop40:  DO  i = 1, n - 1
    DO  j = i + 1, n
      d2 = ((x(i) - x(j))/a) ** 2
      IF (d2 > 50) CYCLE loop40
      e2 = EXP(-0.5_dp*d2)
      s = s + (d2**2 - 6._dp*d2 + 3._dp) * e2
    END DO
  END DO loop40
  r2 = (2.0_dp*s) / ((xn**2)*(a**5)*rt2pi)
  r2 = r2 + 3._dp / (rt2pi*xn*(a**5))
!-
!-------  Estimate bandwidth by asymptotic formula
!-
  h = (const/(r2*xn)) ** (0.2_dp)
  a = co1 * h ** (5.0_dp/7.0_dp)
END DO
!-
!------- Estimate density with plugin bandwidth
!-
jbegin = 1
jend = 1
DO  i = 1, m
  s = 0._dp
  DO  j = jbegin, jend
    t = (z(i) - x(j)) / h
    IF (t > 5.0 .AND. jbegin < n) THEN
      jbegin = jbegin + 1
      CYCLE
    END IF
    s = s + EXP(-t*t/2.)
  END DO
  DO  jend = j, n
    t = (z(i) - x(jend)) / h
    IF (t < -5.0_dp) EXIT
    s = s + EXP(-t*t/2.)
  END DO
  f(i) = s / (xn*h*rt2pi)
  jend = jend - 1
!-
END DO

RETURN
END SUBROUTINE plugin

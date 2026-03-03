PROGRAM fit_poly
! A simple program to fit a polynomial in one variable.
! Data must be store in the form of pairs, either (x,y) or (y,x)
! Polynomial to be fitted:

! Y = a(0) + a(1).X + a(2).X^2 + ... + a(m).X^m

USE lsq
IMPLICIT NONE

CHARACTER (LEN=50)  :: fname
CHARACTER (LEN= 1)  :: ans
REAL (dp)           :: x(1000), y(1000), xrow(0:20), wt = 1.0_dp, beta(0:20), &
                       var, covmat(231), sterr(0:20), totalSS
INTEGER             :: i, ier, iostatus, j, m, n
LOGICAL             :: fit_const = .TRUE., lindep(0:20), xfirst


WRITE(*, '(a)', ADVANCE='NO') ' Enter name of data set: '
READ(*, *) fname
OPEN(UNIT=8, FILE=fname, STATUS='OLD')
OPEN(9, FILE='fit_poly.out')

WRITE(*, '(a)', ADVANCE='NO') ' Is X value before Y in data file? (Y/N): '
READ(*, *) ans
xfirst = (ans == 'Y') .OR. (ans == 'y')

n = 1
DO
  IF (xfirst) THEN
    READ(8, *, IOSTAT=iostatus) x(n), y(n)
  ELSE
    READ(8, *, IOSTAT=iostatus) y(n), x(n)
  END IF
  IF (iostatus > 0) CYCLE
  IF (iostatus < 0) EXIT
  n = n + 1
END DO
n = n - 1
WRITE(9, '(i4, a, a)') n, ' cases read from file: ', fname

WRITE(*, '(a)', ADVANCE='NO') ' Enter highest power to fit: '
READ(*, *) m

! Least-squares calculations

CALL startup(m, fit_const)
DO i = 1, n
  xrow(0) = 1.0_dp
  DO j = 1, m
    xrow(j) = x(i) * xrow(j-1)
  END DO
  CALL includ(wt, xrow, y(i))
END DO

CALL sing(lindep, ier)
IF (ier /= 0) THEN
  DO i = 0, m
    IF (lindep(i)) WRITE(*, '(a, i3)') ' Singularity detected for power: ', i
    IF (lindep(i)) WRITE(9, '(a, i3)') ' Singularity detected for power: ', i
  END DO
END IF

! Calculate progressive residual sums of squares
CALL ss()
var = rss(m+1) / (n - m - 1)

! Calculate least-squares regn. coeffs.
CALL regcf(beta, m+1, ier)

! Calculate covariance matrix, and hence std. errors of coeffs.
CALL cov(m+1, var, covmat, 231, sterr, ier)

WRITE(*, *) 'Least-squares coefficients & std. errors'
WRITE(9, *) 'Least-squares coefficients & std. errors'
WRITE(*, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
WRITE(9, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
DO i = 0, m
  WRITE(*, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
        i, beta(i), sterr(i), rss(i+1)
  WRITE(9, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
        i, beta(i), sterr(i), rss(i+1)
END DO

WRITE(*, *)
WRITE(9, *)
WRITE(*, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)
WRITE(9, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)
totalSS = rss(1)
WRITE(*, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS
WRITE(9, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS

STOP
END PROGRAM fit_poly

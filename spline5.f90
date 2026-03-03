PROGRAM Quintic_Spline

! Fit a quintic spline with user control of knot positions.
! If the knots are at tk1, tk2, ..., then the fitted spline is
! b0 + b1.t + b2.t^2 + b3.t^3 + b4.t^4 + b5.t^5    for t <= tk1
! b0 + ... + b5.t^5 + b6.(t-tk1)^5                 for tk1 < t <= tk2
! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5  for tk2 < t <= tk3
! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5 + b8.(t-tk3)^5
!                                                  for tk3 < t <= tk4, etc.

! In this version, the knots are evenly spaced.
! Also calculates first & 2nd derivatives of the spline.

! Uses the author's least-squares package in file lsq.f90
! Latest revision - 2 November 2003
! Alan Miller (amiller @ bigpond.net.au)

USE lsq
IMPLICIT NONE

CHARACTER (LEN=50)      :: infile, outfile
INTEGER                 :: i, ier, iostatus, j, n, next_knot, nk, pos
REAL (dp)               :: t, t1, t2, y, d1, d2, dist, fitted
REAL (dp), PARAMETER    :: one = 1.0_dp
REAL (dp), ALLOCATABLE  :: knot(:), xrow(:), b(:)

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter file name: '
  READ(*, *) infile
  OPEN(UNIT=8, FILE=infile, STATUS='OLD', IOSTAT=iostatus)
  IF (iostatus == 0) EXIT
  WRITE(*, '(a, a, a)') ' *** File: ', TRIM(infile), ' not found -- TRY AGAIN!'
END DO

! Name of output file = name of input file with extension .OUT
pos = INDEX(infile, '.')
IF (pos > 0) THEN
  outfile = infile(:pos) // 'out'
ELSE
  outfile = TRIM(infile) // '.out'
END IF
OPEN(UNIT=9, FILE=outfile)

! Display first & last cases in data file

n = 0
DO
  READ(UNIT=8, *, IOSTAT=iostatus) t, y
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  IF (n == 0) THEN
    WRITE(*, '(a, 2g13.4)') ' First t, y = ', t, y
    t1 = t
  END IF
  n = n + 1
END DO
t2 = t
WRITE(*, '(a, 2g13.4)') ' Last  t, y = ', t, y
WRITE(*, '(a, i6)') ' Number of cases read = ', n
REWIND (UNIT=8)

DO
  WRITE(*, '(a)', ADVANCE='NO') ' How many knots do you want? '
  READ(*, *) nk
  IF (nk > n/5) THEN
    WRITE(*, *) ' ** Too many knots requested - TRY AGAIN'
    CYCLE
  ELSE IF (nk < 0) THEN
    WRITE(*, *) ' ** Must be >= 0 - TRY AGAIN'
    CYCLE
  ELSE
    EXIT
  END IF
END DO

! Calculate knot positions, evenly spaced.
ALLOCATE ( knot(nk) )
dist = (t2 - t1) / (nk + 1)
DO i = 1, nk
  knot(i) = t1 + dist * i
END DO
WRITE(UNIT=9, '(a, i4)') 'Number of knots = ', nk

! Read in the data & update the least-squares calculations.
ALLOCATE ( xrow(0:5+nk), b(0:5+nk) )
next_knot = 1

! Initialize the least-squares calculations
CALL startup(6+nk, .FALSE.)

DO
  READ(UNIT=8, *, IOSTAT=iostatus) t, y
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  xrow(0) = one
  xrow(1) = (t - t1)
  xrow(2) = (t - t1) * xrow(1)
  xrow(3) = (t - t1) * xrow(2)
  xrow(4) = (t - t1) * xrow(3)
  xrow(5) = (t - t1) * xrow(4)
  IF (t > knot(next_knot)) THEN
    next_knot = MIN(nk, next_knot + 1)
  END IF
  DO j = 1, next_knot-1
    xrow(5+j) = (t - knot(j))**5
  END DO
  xrow(5+next_knot:5+nk) = 0.0_dp
  CALL includ(one, xrow, y)
END DO

CALL regcf(b, 6+nk, ier)

WRITE(*, *) ' Coefficient   Value'
WRITE(*, '(a, g13.5)') ' Constant   ', b(0)
WRITE(*, '(a, g13.5)') ' Linear     ', b(1)
WRITE(*, '(a, g13.5)') ' Quadratic  ', b(2)
WRITE(*, '(a, g13.5)') ' Cubic      ', b(3)
WRITE(*, '(a, g13.5)') ' Quartic    ', b(4)
WRITE(*, '(a, g13.5)') ' Quintic    ', b(5)
WRITE(*, *) ' Knot position   Quintic Coefficient'
DO j = 1, nk
  WRITE(*, '(g13.5, t17, g13.5)') knot(j), b(5+j)
END DO

! Calculate fitted values and derivatives
next_knot = 1
REWIND (UNIT=8)
WRITE(UNIT=9, '(a)') '   Time  1st deriv.   2nd deriv.   Fitted value   Actual'
DO i = 1, n
  READ(UNIT=8, *) t, y
  xrow(0) = one
  xrow(1) = (t - t1)
  xrow(2) = (t - t1) * xrow(1)
  xrow(3) = (t - t1) * xrow(2)
  xrow(4) = (t - t1) * xrow(3)
  xrow(5) = (t - t1) * xrow(4)
  IF (t > knot(next_knot)) THEN
    next_knot = next_knot + 1
    IF (next_knot <= nk) THEN
      WRITE(UNIT=9, '(a, g13.5)') 'New knot at t = ', knot(next_knot-1)
    ELSE
      next_knot = nk + 1
    END IF
  END IF
  DO j = 1, next_knot-1
    xrow(5+j) = (t - knot(j))**5
  END DO
  fitted = DOT_PRODUCT( b(0:5+next_knot-1), xrow(0:5+next_knot-1) )
  d2 = ((20*b(5)*(t-t1) + 12*b(4))*(t-t1) + 6*b(3))*(t-t1) + 2*b(2)
  d1 = (((5*b(5)*(t-t1) + 4*b(4))*(t-t1) + 3*b(3))*(t-t1) + 2*b(2))*(t-t1) + b(1)
  DO j = 1, next_knot-1
    d1 = d1 + 5*b(j+5)*(t - knot(j))**4
    d2 = d2 + 20*b(j+5)*(t - knot(j))**3
  END DO
  WRITE(UNIT=9, '(f8.3, 4g13.4)') t, d1, d2, fitted, y
END DO

STOP
END PROGRAM Quintic_Spline

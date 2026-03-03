PROGRAM Weighted_Quintic

! Fit a quintic polynomial using exponential weighting.
! Full weight = 1.0 is given to the latest case, weight = lambda
! is given to the previous case, lambda^2 to the previous case, etc.

! Also calculates first & 2nd derivatives of the polynomial.

! Uses the author's least-squares package in file lsq.f90
! Latest revision - 11 November 2001
! Alan Miller (amiller @ bigpond.net.au)

USE lsq
IMPLICIT NONE

CHARACTER (LEN=50)    :: infile, outfile
CHARACTER (LEN=10)    :: text
INTEGER               :: ier, iostatus, pos
REAL (dp)             :: lambda, xrow(0:5), t, t0, y, b(0:5), d1, d2, fitted
REAL (dp), PARAMETER  :: one = 1.0_dp

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

! Enter the exponential weight
DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter exponential weight (default = 0.95): '
  READ(*, '(a)') text
  IF (LEN_TRIM(text) == 0) THEN
    lambda = 0.95_dp
    EXIT
  ELSE
    READ(text, '(f9.3)', IOSTAT=iostatus) lambda
    IF (iostatus /= 0) CYCLE
    IF (lambda > 0.0_dp .AND. lambda < one) EXIT
    WRITE(*, *) ' Exponential weight must be between 0 and 1'
  END IF
END DO
WRITE(UNIT=9, '(a, f9.4)') 'Exponential weight = ', lambda

! Initialize the least-squares calculations
CALL startup(6, .FALSE.)

! Read in the data, one line at a time, and update the quintic
WRITE(UNIT=9, '(a)') '   Time  1st deriv.   2nd deriv.   Fitted value   Actual'
DO
  READ(UNIT=8, *, IOSTAT=iostatus) t, y
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  IF (nobs == 0) t0 = t
  xrow(0) = one
  xrow(1) = t - t0
  xrow(2) = (t - t0) * xrow(1)
  xrow(3) = (t - t0) * xrow(2)
  xrow(4) = (t - t0) * xrow(3)
  xrow(5) = (t - t0) * xrow(4)
  CALL includ(one, xrow, y)
  IF (nobs > 5) THEN
    CALL regcf(b, 6, ier)
    d2 = ((20*b(5)*(t-t0) + 12*b(4))*(t-t0) + 6*b(3))*(t-t0) + 2*b(2)
    d1 = (((5*b(5)*(t-t0) + 4*b(4))*(t-t0) + 3*b(3))*(t-t0) + 2*b(2))*(t-t0) + b(1)
    fitted = ((((b(5)*(t-t0) + b(4))*(t-t0) + b(3))*(t-t0) + b(2))*(t-t0) + b(1))*(t-t0) + b(0)
    WRITE(UNIT=9, '(f8.3, 4g13.4)') t, d1, d2, fitted, y
  END IF

! Now downweight the calculations up to this point
  d = lambda * d
END DO

STOP
END PROGRAM Weighted_Quintic

MODULE My_Data
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), SAVE  :: x(100), y(100)
END MODULE My_Data



PROGRAM Test_AS290
USE Nonlin_Conf_Regions
USE My_Data
IMPLICIT NONE

CHARACTER (LEN=20)  :: fname
CHARACTER (LEN=80)  :: text
INTEGER             :: dim1, ifault, fstatus, i, length, line1, ncols, nconf, &
                       ngrid1, ngrid2, nobs, np, np1, np2, row
REAL (dp)           :: conf(6), f(6), grid(31,31), p1lo, p1hi, p2lo, p2hi,  &
                       param(4)
LOGICAL             :: lindep(4)

INTERFACE
  SUBROUTINE Logistic4(param, np, ncols, iobs, grad, resid, wt)
    USE My_Data
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: param(:)
    INTEGER, INTENT(IN)     :: np, ncols, iobs
    REAL (dp), INTENT(OUT)  :: grad(:), resid, wt
  END SUBROUTINE Logistic4
END INTERFACE

DO
  WRITE(*, *) 'Enter name of data file: '
  READ(*, *) fname
  OPEN(10, file=fname, IOSTAT=fstatus)
  IF (fstatus == 0) EXIT
  WRITE(*, *) ' !! File not found -- try again !!'
END DO

WRITE(*, *) ' This is the start of the file:'
DO row = 1, 4
  READ(10, '(a)') text
  length = LEN_TRIM(text)
  WRITE(*, '(" ", a75)') text(1:length)
END DO
REWIND(10)
WRITE(*, *)
WRITE(*, *) ' What line does data start on? '
READ(*, *) line1

IF (line1 > 1) THEN
  DO row = 1, line1-1
    READ(10, *)
  END DO
END IF

nobs = 1
DO
  READ(10, *, iostat=fstatus) x(nobs), y(nobs)
  IF (fstatus < 0) EXIT
  nobs = nobs + 1
END DO
nobs = nobs - 1
WRITE(*, *) ' No. of cases read = ', nobs

! Set up the 2_D grid of variance ratios.
! The values below are for the case in the Applied Statistics paper.
! Change for your needs.

param(1) = 0.88_dp
param(2) = 21.3_dp
np = 4
ncols = 4
np1 = 3
np2 = 4
p1lo = 0.4_dp
p1hi = 1.5_dp
p2lo = 5.8_dp
p2hi = 7.2_dp
conf(1:6) = (/ 10._dp, 50._dp, 90._dp, 95._dp, 98._dp, 99._dp /)
nconf = 6
dim1 = 31
ngrid1 = 31
ngrid2 = 31
CALL halprn(param, np, ncols, nobs, np1, np2, p1lo, p1hi, p2lo, p2hi,  &
                  conf, nconf, f, grid, dim1, ngrid1, ngrid2, Logistic4,  &
                  lindep, ifault)
WRITE(*, *) ' Percentage points of the F-distribution'
WRITE(*, *) '   %    F-value'
DO i = 1, nconf
  WRITE(*, '(f6.1, f8.2)') conf(i), f(i)
END DO

! Write the grid to an output file for later use with a contouring program.
OPEN(11, file='confgrid.dat')
DO row = 1, ngrid1
  WRITE(11, '(10g13.5)') grid(row,1:ngrid2)
END DO

STOP
END PROGRAM Test_AS290



SUBROUTINE Logistic4(param, np, ncols, iobs, grad, resid, wt)

! Calculate residual and first derivatives for the 4-parameter logistic model:
!      Y = a + (b-a) / [1 + exp{-c(X-d)}]
! for one case IOBS.

USE My_Data
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: param(:)
INTEGER, INTENT(IN)     :: np, ncols, iobs
REAL (dp), INTENT(OUT)  :: grad(:), resid, wt

! Local variables
REAL (dp)  :: denom, diff, expntl
REAL (dp), PARAMETER  :: one = 1.0_dp

wt = one
diff = x(iobs) - param(4)
expntl = EXP(-param(3)*diff)
denom = one + expntl
resid = y(iobs) - param(1) - (param(2) - param(1))/denom
grad(1) = -one + one/denom
grad(2) = -one/denom
grad(3) = -(param(2) - param(1))*diff*expntl/denom**2
grad(4) = (param(2) - param(1))*param(3)*expntl/denom**2

RETURN
END SUBROUTINE Logistic4

MODULE data_storage
! The data is stored in this module.
! This is equivalent to a COMMON area in old Fortran.

IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), SAVE  :: x(100), y(100)

END MODULE data_storage



PROGRAM logistic4
! Fit the 4-parameter logistic curve:

!     Y = A + C / [1 + exp{-B(X - D)}]

! by unweighted least squares.

! Latest revision - 25 October 2001

USE Levenberg_Marquardt
USE data_storage
IMPLICIT NONE

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE

REAL (dp), ALLOCATABLE  :: fvec(:), fjac(:,:)
INTEGER, PARAMETER      :: n = 4                 ! The number of parameters
INTEGER                 :: info, iostatus, ipvt(n), m
REAL (dp)               :: p(n), tol = 1.0E-04_dp

! Read in the data, counting the number of cases.

OPEN(UNIT=8, FILE='lgstic4.dat', STATUS='OLD')
m = 1
DO
  READ(8, *, IOSTAT=iostatus) x(m), y(m)
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE    ! Skips cases with errors - e.g. a heading
  m = m + 1
END DO
m = m - 1
ALLOCATE( fvec(m), fjac(m,n) )

! Set starting values for parameters.
! A = p(1), B = p(2), C = p(3), D = p(4)

p(1) = MINVAL( y(1:m) )
p(3) = MAXVAL( y(1:m) ) - p(1)
p(2) = 2.0_dp / ( MAXVAL( x(1:m) ) - MINVAL( x(1:m) ) )
p(4) = 0.5_dp * ( MAXVAL( x(1:m) ) + MINVAL( x(1:m) ) )

CALL lmder1(fcn, m, n, p, fvec, fjac, tol, info, ipvt)

SELECT CASE (info)
  CASE (:-1)
    WRITE(*, *) 'Users FCN returned INFO = ', -info
  CASE (0)
    WRITE(*, *) 'Improper values for input parameters'
  CASE (1:3)
    WRITE(*, *) 'Convergence'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (4)
    WRITE(*, *) 'Residuals orthogonal to the Jacobian'
    WRITE(*, *) 'There may be an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (5)
    WRITE(*, *) 'Too many calls to FCN'
    WRITE(*, *) 'Either slow convergence, or an error in FCN'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE (6:7)
    WRITE(*, *) 'TOL was set too small'
    WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
  CASE DEFAULT
    WRITE(*, *) 'INFO =', info, ' ???'
END SELECT

STOP
END PROGRAM logistic4



SUBROUTINE fcn(m, n, p, fvec, fjac, iflag)
! Calculate either residuals or the Jacobian matrix.
! A = p(1), B = p(2), C = p(3), D = p(4)
! m = no. of cases, n = no. of parameters (4)

USE data_storage
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: p(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag

! Local variables

REAL (dp), PARAMETER  :: one = 1.0_dp
INTEGER               :: i
REAL (dp)             :: expntl, temp

IF (iflag == 1) THEN
  fvec = y(1:m) - p(1) - p(3)/(one + EXP(-p(2)*(x(1:m) - p(4))))
ELSE IF (iflag == 2) THEN
  fjac(1:m,1) = -one
  DO i = 1, m
    expntl = EXP(-p(2)*(x(i) - p(4)))
    temp = one / (one + expntl)
    fjac(i,2) = - p(3) * temp**2 * (x(i) - p(4)) * expntl
    fjac(i,3) = - temp
    fjac(i,4) =   p(3) * temp**2 * expntl * p(2)
  END DO
END IF

RETURN
END SUBROUTINE fcn

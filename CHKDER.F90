SUBROUTINE chkder(m, n, x, fvec, fjac, xp, fvecp, mode, ERR)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-16  Time: 10:36:21

! N.B. Argument LDFJAC has been removed.

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: fvec(:)
REAL (dp), INTENT(IN)   :: fjac(:,:)
REAL (dp), INTENT(OUT)  :: xp(:)
REAL (dp), INTENT(IN)   :: fvecp(:)
INTEGER, INTENT(IN)     :: mode
REAL (dp), INTENT(OUT)  :: ERR(:)


!     **********

!     subroutine chkder

!     this subroutine checks the gradients of m nonlinear functions
!     in n variables, evaluated at a point x, for consistency with
!     the functions themselves. the user must call chkder twice,
!     first with mode = 1 and then with mode = 2.

!     mode = 1. on input, x must contain the point of evaluation.
!               on output, xp is set to a neighboring point.

!     mode = 2. on input, fvec must contain the functions and the rows of fjac
!                         must contain the gradients of the respective
!                         functions each evaluated at x, and fvecp must contain
!                         the functions evaluated at xp.
!               on output, err contains measures of correctness of the
!                          respective gradients.

!     the subroutine does not perform reliably if cancellation or rounding
!     errors cause a severe loss of significance in the evaluation of a
!     function.  Therefore, none of the components of x should be unusually
!     small (in particular, zero) or any other value which may cause loss of
!     significance.

!     the subroutine statement is

!       subroutine chkder(m, n, x, fvec, fjac, xp, fvecp, mode, err)

!     where

!       m is a positive integer input variable set to the number of functions
!         (i.e. the number of cases in most applications).

!       n is a positive integer input variable set to the number of variables.

!       x is an input array of length n.

!       fvec is an array of length m.  On input when mode = 2,
!         fvec must contain the functions evaluated at x.

!       fjac is an m by n array. on input when mode = 2,
!         the rows of fjac must contain the gradients of
!         the respective functions evaluated at x.

!       ldfjac is a positive integer input parameter not less than m
!         which specifies the leading dimension of the array fjac.

!       xp is an array of length n.  On output when mode = 1,
!         xp is set to a neighboring point of x.

!       fvecp is an array of length m.  On input when mode = 2,
!         fvecp must contain the functions evaluated at xp.

!       mode is an integer input variable set to 1 on the first call and 2 on
!         the second.  Other values of mode are equivalent to mode = 1.

!       err is an array of length m. on output when mode = 2, err contains
!         measures of correctness of the respective gradients.  If there is
!         no severe loss of significance, then if err(i) is 1.0 the i-th
!         gradient is correct, while if err(i) is 0.0 the i-th gradient is
!         incorrect.  For values of err between 0.0 and 1.0, the categorization
!         is less certain.  In general, a value of err(i) greater than 0.5
!         indicates that the i-th gradient is probably correct, while a value
!         of err(i) less than 0.5 indicates that the i-th gradient is probably
!         incorrect.

!     subprograms called

!       minpack supplied ... dpmpar

!       fortran supplied ... ABS,LOG10,SQRT

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, j
REAL (dp) :: eps, epsf, epslog, epsmch, temp
REAL (dp), PARAMETER :: factor = 100._dp, one = 1.0_dp, zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(one)

eps = SQRT(epsmch)

IF (mode /= 2) THEN
  
!        mode = 1.
  
  DO  j = 1, n
    temp = eps * ABS(x(j))
    IF (temp == zero) temp = eps
    xp(j) = x(j) + temp
  END DO
ELSE
  
!        mode = 2.
  
  epsf = factor * epsmch
  epslog = LOG10(eps)
  ERR(1:m) = zero
  DO  j = 1, n
    temp = ABS(x(j))
    IF (temp == zero) temp = one
    DO  i = 1, m
      ERR(i) = ERR(i) + temp * fjac(i,j)
    END DO
  END DO
  DO  i = 1, m
    temp = one
    IF (fvec(i) /= zero.AND.fvecp(i) /= zero .AND.  &
        ABS(fvecp(i)-fvec(i)) >= epsf*ABS(fvec(i))) temp = eps *  &
        ABS((fvecp(i)-fvec(i))/eps-ERR(i)) / (ABS(fvec(i)) + ABS(fvecp(i)))
    ERR(i) = one
    IF (temp > epsmch .AND. temp < eps) ERR(i) = (LOG10(temp) - epslog) / epslog
    IF (temp >= eps) ERR(i) = zero
  END DO
END IF

RETURN

!     last card of subroutine chkder.

END SUBROUTINE chkder

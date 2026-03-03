PROGRAM test_cube_root
USE quadruple_precision
IMPLICIT NONE

TYPE (quad)           :: x, lhs, rhs, diff
REAL (r10)            :: half = 0.5_r10, small = EPSILON(half), r
INTEGER, ALLOCATABLE  :: seed(:)
INTEGER               :: k, i

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

DO i = 1, 20
  CALL RANDOM_NUMBER(x%hi)
  CALL RANDOM_NUMBER(r)
  x%hi = (x%hi - half) / r
  x%lo = x%hi * small
  lhs = cube_root(x)
  IF (x%hi >= 0._r10) THEN
    rhs = EXP( LOG(x) / 3._r10)
  ELSE
    rhs = - EXP( LOG(- x) / 3._r10)
  END IF
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
END DO

STOP

CONTAINS


FUNCTION cube_root(y) RESULT(x)
!     Calculates the cube root of y in quadruple-precision.

! Programmer : Alan Miller, CSIRO Mathematical Information Sciences
! e-mail: alan @ mel.dms.csiro.au   WWW: http://www.mel.dms.csiro.au/~alan
! Fax: (+61) 3-9545-8080

! Latest revision - 6 December 1996

IMPLICIT NONE
TYPE (quad), INTENT(IN) :: y
TYPE (quad)             :: x

!     Local variables
TYPE (quad) :: f, xsq
REAL (r10)  :: one_third = 0.3333333333333333_r10

!     First approximation is the double-precision cube root.
!     Only one iteration is necessary.

x%hi = SIGN(ABS(y%hi) ** one_third, y%hi)

! Put f(x) = x**3 - y
! Then Newton's iterative method uses:
!
! new x = x - f(x) / f'(x)
!
! Halley's method uses:
!
! new x = x - f / (f' - 0.5 * f * f'' / f')
!
! In this case
!
! new x = x - x * f / (3.x**3 - f)

x%lo = 0._r10
xsq  = exactmul2(x%hi, x%hi)
f    = x * xsq - y
x    = x - x * f / (3._r10 * x * xsq - f)

RETURN
END FUNCTION cube_root

END PROGRAM test_cube_root

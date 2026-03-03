MODULE Toms715_Utilities

! Code used by many of the special functions in TOMS algorithm 715

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC :: machar, ren

CONTAINS


!      INCLUDES CHANGES GIVEN IN REMARK BY PRICE, TOMS 22 (2)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01

SUBROUTINE machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp,  &
                  maxexp, eps, epsneg, xmin, xmax)
!----------------------------------------------------------------------
!  This Fortran 77 subroutine is intended to determine the parameters of
!   the floating-point arithmetic system specified below.  The determination
!   of the first three uses an extension of an algorithm due to
!   M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some, but not all,
!   of the improvements suggested by M. Gentleman and S. Marovich, CACM 17
!   (1974), pp. 276-277.  An earlier version of this program was published
!   in the book Software Manual for the Elementary Functions by W. J. Cody
!   and W. Waite, Prentice-Hall, Englewood Cliffs, NJ, 1980.

!  Parameter values reported are as follows:

!     IBETA   - the radix for the floating-point representation
!     IT      - the number of base IBETA digits in the floating-point
!               significand
!     IRND    - 0 if floating-point addition chops
!               1 if floating-point addition rounds, but not in the IEEE style
!               2 if floating-point addition rounds in the IEEE style
!               3 if floating-point addition chops, and there is partial
!                 underflow
!               4 if floating-point addition rounds, but not in the
!                 IEEE style, and there is partial underflow
!               5 if floating-point addition rounds in the IEEE style,
!                 and there is partial underflow
!     NGRD    - the number of guard digits for multiplication with
!               truncating arithmetic.  It is
!               0 if floating-point arithmetic rounds, or if it truncates and
!                 only  IT  base  IBETA digits participate in the post-
!                 normalization shift of the floating-point significand in
!                 multiplication;
!               1 if floating-point arithmetic truncates and more than  IT
!                 base  IBETA  digits participate in the post-normalization
!                 shift of the floating-point significand in multiplication.
!     MACHEP  - the largest negative integer such that
!               1.0 + FLOAT(IBETA)**MACHEP .NE. 1.0, except that
!               MACHEP is bounded below by  -(IT+3)
!     NEGEPS  - the largest negative integer such that
!               1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
!               NEGEPS is bounded below by  -(IT+3)
!     IEXP    - the number of bits (decimal places if IBETA = 10)
!               reserved for the representation of the exponent
!               (including the bias or sign) of a floating-point number
!     MINEXP  - the largest in magnitude negative integer such that
!               FLOAT(IBETA)**MINEXP is positive and normalized
!     MAXEXP  - the smallest positive power of  BETA  that overflows
!     EPS     - FLOAT(IBETA)**MACHEP.
!     EPSNEG  - FLOAT(IBETA)**NEGEPS.
!     XMIN    - the smallest non-vanishing normalized floating-point
!               power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
!     XMAX    - the largest finite floating-point number.  In particular
!               XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
!               Note - on some machines  XMAX  will be only the second, or
!               perhaps third, largest number, being too small by 1 or 2 units
!               in the last digit of the significand.

!  Latest modification: May 30, 1989

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(OUT)    :: ibeta
INTEGER, INTENT(OUT)    :: it
INTEGER, INTENT(OUT)    :: irnd
INTEGER, INTENT(OUT)    :: ngrd
INTEGER, INTENT(OUT)    :: machep
INTEGER, INTENT(OUT)    :: negep
INTEGER, INTENT(OUT)    :: iexp
INTEGER, INTENT(OUT)    :: minexp
INTEGER, INTENT(OUT)    :: maxexp
REAL (dp), INTENT(OUT)  :: eps
REAL (dp), INTENT(OUT)  :: epsneg
REAL (dp), INTENT(OUT)  :: xmin
REAL (dp), INTENT(OUT)  :: xmax

INTEGER    :: i, itemp, iz, j, k,  mx, nxres
REAL (dp)  :: a, b, beta, betain, betah, one, t, temp, tempa, temp1, two,  &
              y, z, zero
!----------------------------------------------------------------------
one = 1.0_dp
two = one + one
zero = one - one
!----------------------------------------------------------------------
!  Determine IBETA, BETA ala Malcolm.
!----------------------------------------------------------------------
a = one
10 a = a + a
temp = a + one
temp1 = temp - a
IF (temp1-one == zero) GO TO 10
b = one
20 b = b + b
temp = a + b
itemp = INT(temp-a)
IF (itemp == 0) GO TO 20
ibeta = itemp
beta = ibeta
!----------------------------------------------------------------------
!  Determine IT, IRND.
!----------------------------------------------------------------------
it = 0
b = one
30 it = it + 1
b = b * beta
temp = b + one
temp1 = temp - b
IF (temp1-one == zero) GO TO 30
irnd = 0
betah = beta / two
temp = a + betah
IF (temp-a /= zero) irnd = 1
tempa = a + beta
temp = tempa + betah
IF (irnd == 0 .AND. temp-tempa /= zero) irnd = 2
!----------------------------------------------------------------------
!  Determine NEGEP, EPSNEG.
!----------------------------------------------------------------------
negep = it + 3
betain = one / beta
a = one
DO  i = 1, negep
  a = a * betain
END DO
b = a
50 temp = one - a
IF (temp-one == zero) THEN
  a = a * beta
  negep = negep - 1
  GO TO 50
END IF
negep = -negep
epsneg = a
!----------------------------------------------------------------------
!  Determine MACHEP, EPS.
!----------------------------------------------------------------------
machep = -it - 3
a = b
60 temp = one + a
IF (temp-one == zero) THEN
  a = a * beta
  machep = machep + 1
  GO TO 60
END IF
eps = a
!----------------------------------------------------------------------
!  Determine NGRD.
!----------------------------------------------------------------------
ngrd = 0
temp = one + eps
IF (irnd == 0 .AND. temp*one-one /= zero) ngrd = 1
!----------------------------------------------------------------------
!  Determine IEXP, MINEXP, XMIN.

!  Loop to determine largest I and K = 2**I such that
!         (1/BETA) ** (2**(I))
!  does not underflow.
!  Exit from loop is signaled by an underflow.
!----------------------------------------------------------------------
i = 0
k = 1
z = betain
t = one + eps
nxres = 0
70 y = z
z = y * y
!----------------------------------------------------------------------
!  Check for underflow here.
!----------------------------------------------------------------------
a = z * one
temp = z * t
IF (.NOT.(a+a == zero .OR. ABS(z) >= y)) THEN
  temp1 = temp * betain
  IF (temp1*beta /= z) THEN
    i = i + 1
    k = k + k
    GO TO 70
  END IF
END IF
IF (ibeta /= 10) THEN
  iexp = i + 1
  mx = k + k
ELSE
!----------------------------------------------------------------------
!  This segment is for decimal machines only.
!----------------------------------------------------------------------
  iexp = 2
  iz = ibeta
  80 IF (k >= iz) THEN
    iz = iz * ibeta
    iexp = iexp + 1
    GO TO 80
  END IF
  mx = iz + iz - 1
END IF
!----------------------------------------------------------------------
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!----------------------------------------------------------------------
90 xmin = y
y = y * betain
!----------------------------------------------------------------------
!  Check for underflow here.
!----------------------------------------------------------------------
a = y * one
temp = y * t
IF (.NOT.((a+a) == zero .OR. ABS(y) >= xmin)) THEN
  k = k + 1
  temp1 = temp * betain
  IF (temp1*beta /= y .OR. temp == y) THEN
    GO TO 90
  ELSE
    nxres = 3
    xmin = y
  END IF
END IF
minexp = -k
!----------------------------------------------------------------------
!  Determine MAXEXP, XMAX.
!----------------------------------------------------------------------
IF (.NOT.(mx > k+k-3 .OR. ibeta == 10)) THEN
  mx = mx + mx
  iexp = iexp + 1
END IF
maxexp = mx + minexp
!----------------------------------------------------------------------
!  Adjust IRND to reflect partial underflow.
!----------------------------------------------------------------------
irnd = irnd + nxres
!----------------------------------------------------------------------
!  Adjust for IEEE-style machines.
!----------------------------------------------------------------------
IF (irnd >= 2) maxexp = maxexp - 2
!----------------------------------------------------------------------
!  Adjust for machines with implicit leading bit in binary significand,
!  and machines with radix point at extreme right of significand.
!----------------------------------------------------------------------
i = maxexp + minexp
IF ((ibeta == 2) .AND. (i == 0)) maxexp = maxexp - 1
IF (i > 20) maxexp = maxexp - 1
IF (a /= y) maxexp = maxexp - 2
xmax = one - epsneg
IF (xmax*one /= xmax) xmax = one - beta * epsneg
xmax = xmax / (beta*beta*beta*xmin)
i = maxexp + minexp + 3
IF (i > 0) THEN
  DO  j = 1, i
    IF (ibeta == 2) xmax = xmax + xmax
    IF (ibeta /= 2) xmax = xmax * beta
  END DO
END IF
RETURN
!---------- Last line of MACHAR ----------
END SUBROUTINE machar



FUNCTION REN() RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01

! N.B. The argument K has been removed; it was not used.

!---------------------------------------------------------------------
!  Random number generator - based on Algorithm 266 by Pike and
!   Hill (modified by Hansson), Communications of the ACM,
!   Vol. 8, No. 10, October 1965.

!  This subprogram is intended for use on computers with fixed point
!   wordlength of at least 29 bits.  It is best if the floating-point
!   significand has at most 29 bits.

!  Latest modification: March 13, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!---------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: fn_val

! Local variables

! DATA iy /100001/
INTEGER, SAVE  :: iy = 100001
REAL (dp), PARAMETER  :: one = 1.0_dp, c1 = 2796203.0_dp, c2 = 1.0E-6_dp,  &
                         c3 = 1.0E-12_dp
!---------------------------------------------------------------------
iy = iy * 125
iy = iy - (iy/2796203) * 2796203
fn_val = (DBLE(iy) / c1) * (one + c2+c3)
RETURN
!---------- Last card of REN ----------
END FUNCTION REN

END MODULE Toms715_Utilities

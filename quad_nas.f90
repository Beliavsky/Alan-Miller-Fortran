MODULE quadruple_precision
!
! This version is for NAS FortranPlus, and uses its 10-byte reals to give
! about 38 decimal digit accuracy in quadruple precision.
!
! N.B. This is quadruple precision implemented in SOFTWARE, hence it is SLOW.
! This package has NOT been tested with other Fortran 90 compilers.
! The operations +-*/ & sqrt will probably work correctly with other compilers
! for PCs.   Functions and routines which initialize quadruple-precision
! constants, particularly the trigonometric and exponential functions will
! only work if the compiler initializes double-precision constants from
! decimal code in exactly the same way as the Lahey compilers.

! Programmer : Alan Miller, CSIRO Mathematical & Information Sciences
! e-mail: amiller @ bigpond.net.au   URL:  www.ozemail.com.au/~milleraj

! Latest revision - 18 February 2002

! 4 Aug 97.  Added new algorithm for exponential.  Takes half the time of the
!            previous Taylor series, but errors 2.5 times larger.   The old
!            exponential is included here under the name exp_taylor for anyone
!            who needs the extra accuracy.  To use it instead of the new
!            algorithm, change the module procedure name in INTERFACE exp from
!            longexp to exp_taylor.
! 5 Aug 97.  Found way to reduce cancellation errors in yesterday's algorithm
!            for the exponential.   Removed exp_taylor.
! 18 Aug 97. Added table of quadruple-precision constants.
! 8 Sept 97. Added str_quad which reads a character string and converts it to
!            quadruple-precision form.
! 15 Oct 97. Added quad_str which takes a quadruple-precision value and outputs
!            a character string containing its decimal value to 30 significant
!            digits.
!            Added overlays to the ** operator for quadruple-precision values.
! 15 Jan 98. Added ATAN2.
! 19 Jan 98. Added <, <=, ==, >= and >.
! 27 Dec 98. Added quad/real, quad/integer, integer/quad, real/quad, dp/quad,
!            int+quad, real+quad, dp+quad, int-quad, real-quad, dp-quad,
!            SUM, DOT_PRODUCT & MATMUL.
! 29 Dec 98. Correction to routine string_quad for strings of < 5 characters.
! 10 May 99. Added EPSILON for quad type.
! 5 Oct 99.  log10e corrected.
! 15 Oct 99. Corrected function quad_pow_int.
! 18 Oct 99. Rewrote quad_pow_int to use the binary power method.
! 11 Nov 99. Added overlaying of assignments, e.g. quad = int, etc.
! 21 Jan 00. Added interface for EPSILON.
! 17 Feb 02. Corrected ACOS & ASIN for arguments close to +1 or -1.
! 18 Feb 02. Further improvement to ACOS.

IMPLICIT NONE

INTEGER, PARAMETER    :: r10 = 3        ! 10-byte real KIND
INTEGER, PARAMETER    :: nbits = DIGITS(1.0_r10)
REAL (r10), PARAMETER :: constant = 2._r10**(nbits - nbits/2) + 1._r10
                        !
                        ! Special for 10-byte reals
REAL (r10), PARAMETER :: zero = 0.0_r10, one = 1.0_r10

PRIVATE :: zero, one, constant, nbits

TYPE    :: quad
  REAL (r10) :: hi, lo
END TYPE quad

INTERFACE OPERATOR (*)
  MODULE PROCEDURE longmul
  MODULE PROCEDURE mult_quad_int
  MODULE PROCEDURE mult_int_quad
  MODULE PROCEDURE mult_quad_r10
  MODULE PROCEDURE mult_r10_quad
  MODULE PROCEDURE mult_quad_real
  MODULE PROCEDURE mult_real_quad
END INTERFACE

INTERFACE OPERATOR (/)
  MODULE PROCEDURE longdiv
  MODULE PROCEDURE div_quad_r10
  MODULE PROCEDURE div_quad_real
  MODULE PROCEDURE div_quad_int
  MODULE PROCEDURE div_int_quad
  MODULE PROCEDURE div_real_quad
  MODULE PROCEDURE div_r10_quad
END INTERFACE

INTERFACE OPERATOR (+)
  MODULE PROCEDURE longadd
  MODULE PROCEDURE quad_add_int
  MODULE PROCEDURE quad_add_real
  MODULE PROCEDURE quad_add_r10
  MODULE PROCEDURE int_add_quad
  MODULE PROCEDURE real_add_quad
  MODULE PROCEDURE r10_add_quad
END INTERFACE

INTERFACE OPERATOR (-)
  MODULE PROCEDURE longsub
  MODULE PROCEDURE quad_sub_int
  MODULE PROCEDURE quad_sub_real
  MODULE PROCEDURE quad_sub_r10
  MODULE PROCEDURE int_sub_quad
  MODULE PROCEDURE real_sub_quad
  MODULE PROCEDURE r10_sub_quad
  MODULE PROCEDURE negate_quad
END INTERFACE

INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE quad_eq_int
  MODULE PROCEDURE quad_eq_real
  MODULE PROCEDURE quad_eq_dp
  MODULE PROCEDURE int_eq_quad
  MODULE PROCEDURE real_eq_quad
  MODULE PROCEDURE dp_eq_quad
END INTERFACE

INTERFACE OPERATOR (**)
  MODULE PROCEDURE quad_pow_int
  MODULE PROCEDURE quad_pow_real
  MODULE PROCEDURE quad_pow_r10
  MODULE PROCEDURE quad_pow_quad
END INTERFACE

INTERFACE OPERATOR (<)
  MODULE PROCEDURE quad_lt
END INTERFACE

INTERFACE OPERATOR (<=)
  MODULE PROCEDURE quad_le
END INTERFACE

INTERFACE OPERATOR (==)
  MODULE PROCEDURE quad_eq
END INTERFACE

INTERFACE OPERATOR (>=)
  MODULE PROCEDURE quad_ge
END INTERFACE

INTERFACE OPERATOR (>)
  MODULE PROCEDURE quad_gt
END INTERFACE

INTERFACE SCALE
  MODULE PROCEDURE qscale
END INTERFACE

INTERFACE ABS
  MODULE PROCEDURE qabs
END INTERFACE

INTERFACE SQRT
  MODULE PROCEDURE longsqrt
END INTERFACE

INTERFACE LOG
  MODULE PROCEDURE longlog
END INTERFACE

INTERFACE EXP
  MODULE PROCEDURE longexp
END INTERFACE

INTERFACE SIN
  MODULE PROCEDURE longsin
END INTERFACE

INTERFACE COS
  MODULE PROCEDURE longcos
END INTERFACE

INTERFACE TAN
  MODULE PROCEDURE longtan
END INTERFACE

INTERFACE ASIN
  MODULE PROCEDURE longasin
END INTERFACE

INTERFACE ACOS
  MODULE PROCEDURE longacos
END INTERFACE

INTERFACE ATAN
  MODULE PROCEDURE longatan
END INTERFACE

INTERFACE ATAN2
  MODULE PROCEDURE qatan2
END INTERFACE

INTERFACE SUM
  MODULE PROCEDURE quad_sum
END INTERFACE

INTERFACE DOT_PRODUCT
  MODULE PROCEDURE quad_dot_product
END INTERFACE

INTERFACE MATMUL
  MODULE PROCEDURE q_matmul12
  MODULE PROCEDURE q_matmul21
  MODULE PROCEDURE q_matmul22
END INTERFACE

INTERFACE EPSILON
  MODULE PROCEDURE q_epsilon
END INTERFACE


TYPE (quad), PARAMETER ::  &
          pi = quad( 3.141592653589793239_r10, -0.4838464451208850938E-18_r10 ), &
       piby2 = quad( 1.570796326794896619_r10,  0.3001778636823096701E-18_r10 ), &
       piby3 = quad( 1.047197551196597746_r10,  0.1278384309558394845E-18_r10 ), &
       piby4 = quad( 0.7853981633974483096_r10, 0.4166871459260439164E-19_r10 ), &
       piby6 = quad( 0.5235987755982988731_r10, -0.4450100177063070117E-19_r10 ), &
       twopi = quad( 6.283185307179586477_r10, -0.1003311522533666405E-18_r10 ), &
       ln_pi = quad( 1.144729885849400174_r10,  0.1828709585926350065E-18_r10 ), &
      sqrtpi = quad( 1.772453850905516027_r10,  0.3124892482128265356E-18_r10 ), &
    fact_pt5 = quad( 0.8862269254527580137_r10, -0.6059581039068761899E-19_r10 ), &
     sqrt2pi = quad( 2.506628274631000502_r10,  0.3780482144527760417E-18_r10 ), &
   lnsqrt2pi = quad( 0.9189385332046727418_r10, -0.2271391431563230654E-19_r10 ), &
    two_onpi = quad( 0.6366197723675813431_r10, -0.8144891411629536075E-21_r10 ), &
 two_on_rtpi = quad( 1.128379167095512574_r10, -0.6021123612828221731E-19_r10 ), &
     deg2rad = quad( 0.1745329251994329577E-01_r10, -0.1596304451988263418E-20_r10 ), &
     rad2deg = quad( 57.29577951308232088_r10, -0.3325910540161179126E-17_r10 ), &
         ln2 = quad( 0.6931471805599453094_r10,  0.4275175589747648889E-19_r10 ), &
        ln10 = quad( 2.302585092994045684_r10, -0.1834740662695737210E-19_r10 ), &
       log2e = quad( 1.442695040888963407_r10,  0.4059537671977504638E-18_r10 ), &
      log10e = quad( 0.4342944819032518276_r10,  0.5985961437619392648E-19_r10 ), &
     log2_10 = quad( 3.321928094887362348_r10, -0.1548615619171673199E-18_r10 ), &
     log10_2 = quad( 0.3010299956639811952_r10,  0.1519752067714118941E-19_r10 ), &
       euler = quad( 0.5772156649015328606_r10, -0.9795267621599925471E-20_r10 ), &
           e = quad( 2.718281828459045235_r10,  0.3658002323529239324E-18_r10 ), &
       sqrt2 = quad( 1.414213562373095049_r10, -0.1789397833192357452E-18_r10 ), &
       sqrt3 = quad( 1.732050807568877294_r10, -0.4799598230828858348E-18_r10 ), &
      sqrt10 = quad( 3.162277660168379332_r10,  0.3141265028276078512E-19_r10 )

CONTAINS


FUNCTION exactmul2(a, c) RESULT(ac)
!  Procedure exactmul2, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL (r10), INTENT(IN)  :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (r10)  :: a1, a2, c1, c2, t

t = constant * a
a1 = (a - t) + t             ! Lahey's optimization removes the brackets
                             ! and sets a1 = a which defeats the whole point.
a2 = a - a1
t = constant * c
c1 = (c - t) + t
c2 = c - c1
ac%hi = a * c
ac%lo = (((a1 * c1 - ac%hi) + a1 * c2) + c1 * a2) + c2 * a2

RETURN
END FUNCTION exactmul2



FUNCTION longmul(a, c) RESULT(ac)
!  Procedure longmul, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (r10  ) :: zz
TYPE (quad)  :: z

z = exactmul2(a%hi, c%hi)
zz = ((a%hi + a%lo) * c%lo + a%lo * c%hi) + z%lo
ac%hi = z%hi + zz
ac%lo = (z%hi - ac%hi) + zz

RETURN
END FUNCTION longmul



FUNCTION mult_quad_int(a, i) RESULT(c)
!     Multiply quadruple-precision number (a) by an integer (i).

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: c

IF (i == 0) THEN
  c = quad(zero, zero)
ELSE IF (i == 1) THEN
  c = a
ELSE IF (i == -1) THEN
  c = -a
ELSE
  c = exactmul2(a%hi, REAL(i, KIND=r10)) + exactmul2(a%lo, REAL(i, KIND=r10))
END IF

RETURN
END FUNCTION mult_quad_int



FUNCTION mult_int_quad(i, a) RESULT(c)
!     Multiply quadruple-precision number (a) by an integer (i).

INTEGER, INTENT(IN)     :: i
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: c

IF (i == 0) THEN
  c = quad(zero, zero)
ELSE IF (i == 1) THEN
  c = a
ELSE IF (i == -1) THEN
  c = -a
ELSE
  c = exactmul2(a%hi, REAL(i, KIND=r10)) + exactmul2(a%lo, REAL(i, KIND=r10))
END IF

RETURN
END FUNCTION mult_int_quad



FUNCTION mult_quad_r10(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (r10)              :: zz

z = exactmul2(a%hi, b)
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_quad_r10



FUNCTION mult_quad_real(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a real number (b).

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (r10)              :: zz

z = exactmul2(a%hi, REAL(b, KIND=r10))
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_quad_real



FUNCTION mult_r10_quad(b, a) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (r10)              :: zz

z = exactmul2(a%hi, b)
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_r10_quad



FUNCTION mult_real_quad(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (r10)              :: zz

z = exactmul2(b%hi, REAL(a, KIND=r10))
zz = b%lo * a + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_real_quad



FUNCTION longdiv(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (r10  ) :: z, zz
TYPE (quad)  :: q

z = a%hi / c%hi
q = exactmul2(c%hi, z)
zz = ((((a%hi - q%hi) - q%lo) + a%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longdiv



FUNCTION div_quad_r10(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by a double-precision number (b)

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: b
TYPE (quad)             :: c

!     Local variables
REAL (r10  ) :: z, zz
TYPE (quad) :: q

z = a%hi / b
q = exactmul2(b, z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_r10


FUNCTION div_quad_real(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by a real number (b)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: b
TYPE (quad)             :: c

!     Local variables
REAL (r10  ) :: z, zz
TYPE (quad)  :: q

z = a%hi / b
q = exactmul2(REAL(b, KIND=r10), z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_real


FUNCTION div_quad_int(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by an integer (b)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: b
TYPE (quad)             :: c

!     Local variables
REAL (r10 ) :: z, zz
TYPE (quad) :: q

z = a%hi / b
q = exactmul2(REAL(b, KIND=r10), z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_int



FUNCTION div_int_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

INTEGER, INTENT(IN)     :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, zz
TYPE (quad) :: q

z = REAL(a, KIND=r10) / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_int_quad



FUNCTION div_real_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, zz
TYPE (quad) :: q

z = REAL(a, KIND=r10) / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_real_quad



FUNCTION div_r10_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL (r10), INTENT(IN)  :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, zz
TYPE (quad) :: q

z = a / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_r10_quad



FUNCTION longadd(a, c) RESULT(ac)
!  Procedure longadd, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c%hi
q = a%hi - z
zz = (((q + c%hi) + (a%hi - (q + z))) + a%lo) + c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longadd



FUNCTION quad_add_int(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_int



FUNCTION quad_add_real(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_real



FUNCTION quad_add_r10(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_r10



FUNCTION int_add_quad(c, a) RESULT(ac)

INTEGER, INTENT(IN)     :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION int_add_quad



FUNCTION real_add_quad(c, a) RESULT(ac)

REAL, INTENT(IN)        :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION real_add_quad



FUNCTION r10_add_quad(c, a) RESULT(ac)

REAL (r10), INTENT(IN)  :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION r10_add_quad



FUNCTION longsub(a, c) RESULT(ac)
!  Adapted from longadd by changing signs of c.
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi - c%hi
q = a%hi - z
zz = (((q - c%hi) + (a%hi - (q + z))) + a%lo) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longsub



FUNCTION quad_sub_int(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_int



FUNCTION quad_sub_real(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_real



FUNCTION quad_sub_r10(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_r10



FUNCTION int_sub_quad(a, c) RESULT(ac)

INTEGER, INTENT(IN)     :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = REAL(a, KIND=r10) - c%hi
q = REAL(a, KIND=r10) - z
zz = ((q - c%hi) + (DBLE(a) - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION int_sub_quad



FUNCTION real_sub_quad(a, c) RESULT(ac)

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a - c%hi
q = a - z
zz = ((q - c%hi) + (a - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION real_sub_quad



FUNCTION r10_sub_quad(a, c) RESULT(ac)

REAL (r10), INTENT(IN)  :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (r10 ) :: z, q, zz

z = a - c%hi
q = a - z
zz = ((q - c%hi) + (a - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION r10_sub_quad



FUNCTION negate_quad(a) RESULT(b)
!     Change the sign of a quadruple-precision number.
!     In many cases, a & b will occupy the same locations.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

b%hi = -a%hi
b%lo = -a%lo

RETURN
END FUNCTION negate_quad



SUBROUTINE quad_eq_int(a, i)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
INTEGER, INTENT(IN)      :: i

a%hi = i
a%lo = 0

RETURN
END SUBROUTINE quad_eq_int



SUBROUTINE quad_eq_real(a, r)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
REAL, INTENT(IN)         :: r

a%hi = r
a%lo = zero

RETURN
END SUBROUTINE quad_eq_real



SUBROUTINE quad_eq_dp(a, d)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
REAL (r10), INTENT(IN)   :: d

a%hi = d
a%lo = zero

RETURN
END SUBROUTINE quad_eq_dp



SUBROUTINE int_eq_quad(i, a)
!     Assignment

INTEGER, INTENT(OUT)     :: i
TYPE (quad), INTENT(IN)  :: a

i = a%hi

RETURN
END SUBROUTINE int_eq_quad



SUBROUTINE real_eq_quad(r, a)
!     Assignment

REAL, INTENT(OUT)        :: r
TYPE (quad), INTENT(IN)  :: a

r = a%hi

RETURN
END SUBROUTINE real_eq_quad



SUBROUTINE dp_eq_quad(d, a)
!     Assignment

REAL (r10), INTENT(OUT)  :: d
TYPE (quad), INTENT(IN)  :: a

d = a%hi

RETURN
END SUBROUTINE dp_eq_quad



FUNCTION quad_lt(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi < zero)

RETURN
END FUNCTION quad_lt



FUNCTION quad_le(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = .NOT. (diff%hi > zero)

RETURN
END FUNCTION quad_le



FUNCTION quad_eq(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi == zero)

RETURN
END FUNCTION quad_eq



FUNCTION quad_ge(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = .NOT. (diff%hi < zero)

RETURN
END FUNCTION quad_ge



FUNCTION quad_gt(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi > zero)

RETURN
END FUNCTION quad_gt



FUNCTION quad_pow_int(a, i) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: b

! Local variables

INTEGER     :: ia, j, emax
TYPE (quad) :: power
LOGICAL     :: first

SELECT CASE (i)
  CASE (0)
    b = quad(one, zero)
  CASE (-1023:-1, 1:1023)
    ia = ABS(i)
    first = .TRUE.
    power = a
    emax = MAXEXPONENT(one)
    DO j = 0, 9
      IF (BTEST(ia, j)) THEN
        IF (first) THEN
          b = power
          first = .FALSE.
        ELSE
          IF (EXPONENT(power%hi) + EXPONENT(b%hi) > emax) THEN
            WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
            RETURN
          END IF
          b = b * power
        END IF
        ia = IBCLR(ia, j)
        IF (ia == 0) EXIT
      END IF
      IF (EXPONENT(power%hi) > emax/2) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      END IF
      power = power * power
    END DO
    IF (i < 0) b = quad(one, zero) / b
  CASE DEFAULT
    IF (a%hi < zero) THEN
      IF (i * LOG(-a%hi) > LOG( HUGE(one) )) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      ELSE IF (i * LOG(-a%hi) < LOG( TINY(one) )) THEN
        b = quad(zero, zero)
      ELSE
        ia = ABS(i)
        IF (BTEST(ia,0)) THEN
          b = - EXP( LOG(-a)*i )         ! i is odd (test least sign. bit)
        ELSE
          b = EXP( LOG(-a)*i )           ! i is even
        END IF
      END IF
    ELSE IF (a%hi > zero) THEN
      IF (i * LOG(a%hi) > LOG( HUGE(one) )) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      ELSE IF (i * LOG(a%hi) < LOG( TINY(one) )) THEN
        b = quad(zero, zero)
      ELSE
        b = EXP( LOG(a)*i )
      END IF
    ELSE
      b = quad(zero, zero)
    END IF
END SELECT

RETURN
END FUNCTION quad_pow_int



FUNCTION quad_pow_real(a, r) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: r
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*r )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_real



FUNCTION quad_pow_r10(a, d) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
REAL (r10), INTENT(IN)  :: d
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*d )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_r10



FUNCTION quad_pow_quad(a, q) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a, q
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*q )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_quad



FUNCTION qscale(a, i) RESULT(b)
!     Multiply a by 2^i

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: b

b%hi = SCALE(a%hi, i)
b%lo = SCALE(a%lo, i)

RETURN
END FUNCTION qscale



FUNCTION qabs(a) RESULT(b)
!     Absolute value of a quadruple-precision number

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  b%hi = - a%hi
  b%lo = - a%lo
ELSE
  b%hi = a%hi
  b%lo = a%lo
END IF

RETURN
END FUNCTION qabs



FUNCTION longsqrt(a) RESULT(b)
! This is modified from procedure sqrt2 of:
!    Dekker, T.J. (1971). 'A floating-point technique for extending the
!    available precision', Numer. Math., 18, 224-242.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

!     Local variables
REAL (r10)  :: t, res
TYPE (quad) :: tt

! Check that ahi >= 0.

IF (a%hi < 0._r10) THEN
  WRITE(*, *) ' *** Negative argument for longsqrt ***'
  RETURN
ELSE IF (a%hi == 0.d0) THEN
  b%hi = 0._r10
  b%lo = 0._r10
  RETURN
END IF

! First approximation is  t = sqrt(a).

t = SQRT(a%hi)
tt = exactmul2(t, t)
res = t + (((a%hi - tt%hi) - tt%lo) + a%lo) * 0.5_r10 / t
b%lo = (t - res) + (((a%hi - tt%hi) - tt%lo) + a%lo) * 0.5_r10 / t
b%hi = res

RETURN
END FUNCTION longsqrt



FUNCTION longlog(x) RESULT(y)
!  Quadruple-precision logarithm to base e
!  Halley's algorithm using double-precision logarithm as starting value.
!  Solves:         y
!          f(y) = e  - x = 0

TYPE (quad), INTENT(IN) :: x
TYPE (quad)             :: y

!     Local variables
TYPE (quad) :: expy, f

y%hi = LOG(x%hi)
y%lo = 0._r10
expy = EXP(y)
f = expy - x
f = SCALE(f, 1)
y = y - f / (expy + x)

RETURN
END FUNCTION longlog



FUNCTION longexp(x) RESULT(y)
!  Calculate a quadruple-precision exponential
!  Method:
!   x    x.log2(e)    nint[x.log2(e)] + frac[x.log2(e)]
!  e  = 2          = 2
!
!                     iy    fy
!                  = 2   . 2
!  Then
!   fy    y.ln(2)
!  2   = e
!
!  Now y.ln(2) will be less than 0.3466 in absolute value.
!  This is halved and a Pade approximation is used to approximate e^x over
!  the region (-0.1733, +0.1733).   This approximation is then squared.

!  WARNING: No overflow checks!

TYPE (quad), INTENT(IN) :: x
TYPE (quad)             :: y

! Local variables
TYPE (quad)  :: temp, ysq, sum1, sum2
INTEGER      :: iy

y = x / ln2
iy = NINT(REAL(y%hi))
y = (y - REAL(iy, KIND=r10)) * ln2
y = SCALE(y, -1)

! The Pade series is:
!     p0 + p1.y + p2.y^2 + p3.y^3 + ... + p9.y^9
!     ------------------------------------------
!     p0 - p1.y + p2.y^2 - p3.y^3 + ... - p9.y^9
!
! sum1 is the sum of the odd powers, sum2 is the sum of the even powers

ysq = y * y
sum1 = y * ((((ysq + 3960.)*ysq + 2162160._r10)*ysq + 302702400._r10)*ysq +   &
               8821612800._r10)
sum2 = (((90.*ysq + 110880.)*ysq + 30270240._r10)*ysq + 2075673600._r10)*ysq +  &
       17643225600._r10

!                     sum2 + sum1         2.sum1
! Now approximation = ----------- = 1 + ----------- = 1 + 2.temp
!                     sum2 - sum1       sum2 - sum1
!
! Then (1 + 2.temp)^2 = 4.temp.(1 + temp) + 1

temp = sum1 / (sum2 - sum1)
y = temp * (temp + 1._r10)
y = SCALE(y, 2)
y = y + 1._r10
y = SCALE(y, iy)

RETURN
END FUNCTION longexp



SUBROUTINE longmodr(a, b, n, rem)

! Extended arithmetic calculation of the 'rounded' modulus:
!  a = n.b + rem
! where all quantities are in quadruple-precision, except the integer
! number of multiples, n.   The absolute value of the remainder (rem)
! is not greater than b/2.
! The result is exact.   rem may occupy the same location as either input.

! Programmer: Alan Miller

! Latest revision - 11 September 1986
! Fortran version - 4 December 1996

TYPE (quad), INTENT(IN)  :: a, b
INTEGER, INTENT(OUT)     :: n
TYPE (quad), INTENT(OUT) :: rem

! Local variables

TYPE (quad) :: temp

! Check that b%hi .ne. 0

IF (b%hi == 0._r10) THEN
  WRITE(*, *) ' *** Error in longmodr - 2nd argument zero ***'
  RETURN
END IF

! Calculate n.

temp = a / b
n = NINT(REAL(temp%hi))

! Calculate remainder.

rem = a - b*n

RETURN
END SUBROUTINE longmodr



! Extended accuracy arithmetic sine, cosine & tangent (about 38 decimals).
! Calculates  b = sin, cos or tan (a), where all quantities are in
! quadruple-precision, using table look-up and a Taylor series expansion.
! The result may occupy the same locations as the input value.
! Much of the code is common to all three functions, and this is in a
! subroutine longcst.


FUNCTION longsin(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .true.
cosine = .false.
tangent = .false.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longsin



FUNCTION longcos(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .false.
cosine = .true.
tangent = .false.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longcos



FUNCTION longtan(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .false.
cosine = .false.
tangent = .true.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longtan



SUBROUTINE longcst(a, b, sine, cosine, tangent)

TYPE (quad), INTENT(IN)  :: a
TYPE (quad), INTENT(OUT) :: b
LOGICAL, INTENT(IN)      :: sine, cosine, tangent

! Local variables

LOGICAL     :: pos
TYPE (quad) :: d, term, temp, angle, piby40, sum1, sum2, sin
INTEGER     :: npi, ipt, i
REAL (r10)  :: tol19 = 1.E-19_r10, tol38 = 1.E-38_r10

! sin(i.pi/40), i = 0(1)20
REAL (r10), PARAMETER   :: table(2, 0:20) = RESHAPE( (/    &
   0.0000000000000000000_r10 ,  0.000000000000000000          , &
   0.07845909572784494503_r10,  0.7243166770598514154E-21_r10 , &
   0.1564344650402308690_r10 ,  0.9822383634853578774E-20_r10 , &
   0.2334453638559054118_r10 , -0.2727149579911293712E-19_r10 , &
   0.3090169943749474241_r10 , -0.1311597650643346617E-20_r10 , &
   0.3826834323650897717_r10 ,  0.3230750765360341853E-19_r10 , &
   0.4539904997395467916_r10 , -0.4542956404774908511E-19_r10 , &
   0.5224985647159488650_r10 , -0.1417786591071146915E-19_r10 , &
   0.5877852522924731292_r10 , -0.2007497904583745515E-19_r10 , &
   0.6494480483301836557_r10 ,  0.2557134761526789203E-19_r10 , &
   0.7071067811865475244_r10 ,  0.1895032558893257080E-19_r10 , &
   0.7604059656000309382_r10 ,  0.5224249890853254941E-20_r10 , &
   0.8090169943749474241_r10 , -0.1311597650643346617E-20_r10 , &
   0.8526401643540922215_r10 ,  0.4989584747625498625E-20_r10 , &
   0.8910065241883678624_r10 , -0.1994650995938756338E-19_r10 , &
   0.9238795325112867561_r10 ,  0.2676178144722950291E-19_r10 , &
   0.9510565162951535721_r10 ,  0.5868888759500729479E-20_r10 , &
   0.9723699203976766018_r10 ,  0.6578667829614901267E-19_r10 , &
   0.9876883405951377262_r10 ,  0.2105846888637170833E-19_r10 , &
   0.9969173337331279762_r10 , -0.8374704372080533265E-20_r10 , &
   1.000000000000000000_r10  ,  0.000000000000000000_r10 /),  (/ 2, 21 /) )

! pi/40

piby40%hi = 0.7853981633974483096E-01_r10
piby40%lo = -0.1254139403167083006E-20_r10

! Reduce angle to range (-pi/2, +pi/2) by subtracting an integer multiple of pi.

CALL longmodr(a, pi, npi, angle)

! Find nearest multiple of pi/40 to angle.

CALL longmodr(angle, piby40, ipt, d)

! Sum 1 = 1 - d**2/2! + d**4/4! - d**6/6! + ...
! Sum 2 = d - d**3/3! + d**5/5! - d**7/7! + ...

sum1 = quad(zero, zero)
sum2 = quad(zero, zero)
pos = .false.
term = d
i = 2
20 IF (ABS(term%hi) > tol19) THEN
  term = term * d                                ! Use quad. precision
  IF (i == 2 .OR. i == 4 .OR. i == 8) THEN
    term%hi = term%hi / i
    term%lo = term%lo / i
  ELSE
    term = term / i
  END IF
  IF (pos) THEN
    sum1 = sum1 + term
  ELSE
    sum1 = sum1 - term
  END IF
ELSE
  term%hi = term%hi * d%hi / i                   ! Double prec. adequate
  IF (pos) THEN
    sum1%lo = sum1%lo + term%hi
  ELSE
    sum1%lo = sum1%lo - term%hi
  END IF
END IF

! Repeat for sum2

i = i + 1
IF (ABS(term%hi) > tol19) THEN
  term = term * d / i                            ! Use quad. precision
  IF (pos) THEN
    sum2 = sum2 + term
  ELSE
    sum2 = sum2 - term
  END IF
ELSE
  term%hi = term%hi * d%hi / i                   ! Double prec. adequate
  IF (pos) THEN
    sum2%lo = sum2%lo + term%hi
  ELSE
    sum2%lo = sum2%lo - term%hi
  END IF
END IF

i = i + 1
pos = .NOT. pos
IF (ABS(term%hi) > tol38) GO TO 20

sum1 = sum1 + 1._r10                             ! Now add the 1st terms
sum2 = sum2 + d                                  ! for max. accuracy

! Construct sine, cosine or tangent.
! Sine first.    sin(angle + d) = sin(angle).cos(d) + cos(angle).sin(d)

IF (sine .OR. tangent) THEN
  IF (ipt >= 0) THEN
    temp%hi = table(1, ipt)
    temp%lo = table(2, ipt)
  ELSE
    temp%hi = - table(1, -ipt)
    temp%lo = - table(2, -ipt)
  END IF
  b = sum1 * temp
  IF (ipt >= 0) THEN
    temp%hi = table(1, 20-ipt)
    temp%lo = table(2, 20-ipt)
  ELSE
    temp%hi = table(1, 20+ipt)
    temp%lo = table(2, 20+ipt)
  END IF
  b = b + sum2 * temp
  IF (npi /= 2*(npi/2)) THEN
    b = -b
  END IF
  IF (tangent) THEN
    sin = b
  END IF
END IF

! Cosine or tangent.

IF (cosine .OR. tangent) THEN
  IF (ipt >= 0) THEN
    temp%hi = table(1, ipt)
    temp%lo = table(2, ipt)
  ELSE
    temp%hi = - table(1, -ipt)
    temp%lo = - table(2, -ipt)
  END IF
  b = sum2 * temp
  IF (ipt >= 0) THEN
    temp%hi = table(1, 20-ipt)
    temp%lo = table(2, 20-ipt)
  ELSE
    temp%hi = table(1, 20+ipt)
    temp%lo = table(2, 20+ipt)
  END IF
  b = sum1 * temp - b
  IF (npi /= 2*(npi/2)) THEN
    b = -b
  END IF
END IF

! Tangent.

IF (tangent) THEN

! Check that bhi .ne. 0

  IF (b%hi == zero) THEN
    WRITE(*, *) ' *** Infinite tangent - routine longcst ***'
    b%hi = HUGE(one)
    b%lo = zero
    RETURN
  END IF
  b = sin / b
END IF

RETURN
END SUBROUTINE longcst



FUNCTION longasin(a) RESULT(b)

! Quadratic-precision arc sine (about 31 decimals).
! One Newton-Raphson iteration to solve:  f(b) = sin(b) - a = 0,
! except when a close to -1 or +1.
! The result (b) may occupy the same location as the input values (a).
! Use ACOS when |a| is close to 1.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y, c

! Check that -1 <= a%hi <= +1.

IF (a%hi < -one .OR. a%hi > one) THEN
  WRITE(*, *) ' *** Argument outside range for longasin ***'
  RETURN
END IF

IF (ABS(a%hi) < 0.866) THEN
  ! First approximation is  y = asin(a).
  ! Quadruple-precision result is  y - [sin(y) - a]/cos(y).

  y%hi = ASIN(a%hi)
  y%lo = zero
  b = y + (a - SIN(y)) / COS(y%hi)
ELSE
  ! Calculate acos(c) where c = sqrt(1 - a^2)
  c = SQRT(one - a*a)
  y%hi = ACOS(c%hi)
  y%lo = zero
  b = y + (COS(y) - c) / SIN(y%hi)
  IF (a%hi < zero) b = -b
END IF


RETURN
END FUNCTION longasin



FUNCTION longacos(a) RESULT(b)

! Quadratic-precision arc cosine (about 31 decimals).
! Newton-Raphson iteration to solve: f(b) = cos(b) - a = 0.
! The result (b) may occupy the same location as the input values (a).
! When |a| is near 1, use formula from p.175 of
! `Software Manual for the Elementary Functions' by W.J. Cody, Jr. &
! W. Waite, Prentice-Hall, 1980.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y, c

! Check that -1 <= a%hi <= +1.

IF (a%hi < -one .OR. a%hi > one) THEN
  WRITE(*, *) ' *** Argument outside range for longacos ***'
  RETURN
END IF

IF (ABS(a%hi) < 0.866) THEN
  ! First approximation is  y = acos(a).
  ! Quadruple-precision result is  y + [cos(y) - a]/sin(y).

  y%hi = ACOS(a%hi)
  y%lo = zero
  b = y + (COS(y) - a) / SIN(y%hi)
ELSE
  ! Calculate 2.asin(c) where c = sqrt([1 - |a|]/2)
  c = SQRT((one - ABS(a))/2)
  y%hi = ASIN(c%hi)
  y%lo = zero
  b = (y - (SIN(y) - c) / COS(y%hi))*2
  IF (a%hi < zero) b = pi - b
END IF

RETURN
END FUNCTION longacos



FUNCTION longatan(a) RESULT(b)

! Quadratic-precision arc tangent (about 31 decimals).
! Newton-Raphson iteration to solve: f(b) = tan(b) - a = 0.
! The result (b) may occupy the same location as the input values (a).

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y

! First approximation is  y = atan(a).
! Quadruple-precision result is  y - [tan(y) - a] * cos(y)**2.

y%hi = ATAN(a%hi)
y%lo = 0._r10
b = y - (TAN(y) - a) * (COS(y%hi))**2

RETURN
END FUNCTION longatan



FUNCTION qatan2(y, x) RESULT(b)

! Quadratic-precision arc tangent (about 31 decimals).
! As for arc tangent (y/x) except that the result is in the range
!       -pi < ATAN2 <= pi.
! The signs of x and y determine the quadrant.

TYPE (quad), INTENT(IN) :: y, x
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: z

! First approximation is  z = atan2(y, x).
! Quadruple-precision result is  z - [tan(z) - (y/x)] * cos(z)**2.

z%hi = ATAN2(y%hi, x%hi)
z%lo = 0._r10
IF (x%hi == zero) THEN
  b = z
ELSE
  b = z - (TAN(z) - y/x) * (COS(z%hi))**2
END IF

RETURN
END FUNCTION qatan2



FUNCTION quad_sum(a) RESULT(s)

! Quadruple-precision SUM

TYPE (quad), INTENT(IN), DIMENSION(:) :: a
TYPE (quad)                           :: s

! Local variables

INTEGER :: i

s = quad(zero, zero)

DO i = 1, SIZE(a)
  s = s + a(i)
END DO

RETURN
END FUNCTION quad_sum



FUNCTION quad_dot_product(a, b) RESULT(ab)

! Quadruple-precision DOT_PRODUCT

TYPE (quad), INTENT(IN), DIMENSION(:) :: a, b
TYPE (quad)                           :: ab

! Local variables

INTEGER :: i, n

ab = quad(zero, zero)
n = SIZE(a)
IF (n /= SIZE(b)) THEN
  WRITE(*, *) ' ** Error invoking DOT_PRODUCT - different argument sizes **'
  WRITE(*, '(a, i10, a, i10)') ' Size of 1st argument = ', n,   &
                               '   Size of 2nd argument = ', SIZE(b)
  RETURN
END IF

DO i = 1, n
  ab = ab + a(i)*b(i)
END DO

RETURN
END FUNCTION quad_dot_product



FUNCTION q_matmul12(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 1, rank of B = 2

TYPE (quad), INTENT(IN), DIMENSION(:)   :: a
TYPE (quad), INTENT(IN), DIMENSION(:,:) :: b
TYPE (quad), DIMENSION( SIZE(b,2) )     :: ab

! Local variables

INTEGER :: j, na, nb1, nb2

! Check dimensions

na = SIZE(a)
nb1 = SIZE(b, 1)
nb2 = SIZE(b, 2)
IF (na /= nb1) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO j = 1, nb2
  ab(j) = quad_dot_product( a, b(:,j) )
END DO

RETURN
END FUNCTION q_matmul12



FUNCTION q_matmul21(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 2, rank of B = 1

TYPE (quad), INTENT(IN), DIMENSION(:,:) :: a
TYPE (quad), INTENT(IN), DIMENSION(:)   :: b
TYPE (quad), DIMENSION( SIZE(a,1) )     :: ab

! Local variables

INTEGER :: i, na1, na2, nb

! Check dimensions

na1 = SIZE(a, 1)
na2 = SIZE(a, 2)
nb = SIZE(b)

IF (na2 /= nb) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO i = 1, na1
  ab(i) = quad_dot_product( a(i,:), b )
END DO

RETURN
END FUNCTION q_matmul21



FUNCTION q_matmul22(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 2, rank of B = 2

TYPE (quad), INTENT(IN), DIMENSION(:,:)     :: a
TYPE (quad), INTENT(IN), DIMENSION(:,:)     :: b
TYPE (quad), DIMENSION(SIZE(a,1),SIZE(b,2)) :: ab

! Local variables

INTEGER :: i, j, na1, na2, nb1, nb2

! Check dimensions

na1 = SIZE(a, 1)
na2 = SIZE(a, 2)
nb1 = SIZE(b, 1)
nb2 = SIZE(b, 2)

IF (na2 /= nb1) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO i = 1, na1
  DO j = 1, nb2
    ab(i,j) = quad_dot_product( a(i,:), b(:,j) )
  END DO
END DO

RETURN
END FUNCTION q_matmul22



SUBROUTINE string_quad(string, value, ier)
! Convert a character string to a quadruple-precision quantity.
! Error indicator ier = 0 if value OK
!                     = 1 if string has > 50 characters and no decimal point
!                            in about the first 45
!                     = 2 if a non-numeric character is read in the mantissa
!                            part

CHARACTER (LEN=*), INTENT(IN) :: string
TYPE (quad), INTENT(OUT)      :: value
INTEGER, INTENT(OUT)          :: ier

! Local variables
CHARACTER (LEN=50)  :: str
INTEGER             :: length, pos, decpt, status, power10, i, i1, i2
REAL (r10)          :: sgn, temp
CHARACTER (LEN= 5)  :: expnt

str = ADJUSTL(string)
length = LEN_TRIM( str )
IF (length > 50) THEN
                             ! Truncate string to 50 characters preserving
                             ! the exponent, if present.
  pos = SCAN(str(length-4:length), 'DdEe')
  IF (pos > 5) pos = 0
  IF (pos > 0) THEN
    str(45+pos:50) = str(length-5+pos:length)
  END IF
                             ! Check that the string contains a '.'
  IF (INDEX(str, '.') == 0) THEN
    ier = 1
    value = quad(zero, zero)
    RETURN
  END IF
END IF

ier = 0

! Determine sign of mantissa

sgn = +1._r10
IF (str(1:1) == '-') THEN
  sgn = -1._r10
  str = str(2:)
ELSE IF (str(1:1) == '+') THEN
  str = str(2:)
END IF

! Separate the exponent, if there is one

length = LEN_TRIM(str)
IF (length > 4) THEN
  pos = SCAN(str(length-4:length), 'DdEe')
  IF (pos > 5) pos = 0
  IF (pos > 0) THEN
    expnt = str(length-5+pos:)
    str(length-5+pos:) = '     '
  ELSE
    expnt = ' '
  END IF
ELSE
  expnt = ' '
END IF

! decpt = position of decimal point

decpt = INDEX(str, '.')
IF (decpt > 0) str = str(:decpt-1) // str(decpt+1:)

! Read str in blocks of up to 9 digits at a time, as a large integer.

i1 = 1
length = LEN_TRIM(str)
value = quad(zero, zero)
DO
  i2 = MIN(i1+8, length)
  READ(str(i1:i2), '(f9.0)', IOSTAT=status) temp
  IF (status /= 0) THEN
    ier = 2
    RETURN
  END IF
  IF (i1 > 1) THEN
    value = value * (10._r10 ** (i2+1-i1)) + temp
  ELSE
    value = quad(temp, zero)
  END IF
  i1 = i2 + 1
  IF (i1 > length) EXIT
END DO

IF (sgn < zero) value = -value

! Multiply by appropriate power of 10 for position of the decimal point,
! and the exponent.

IF (expnt == ' ') THEN
  power10 = 0
ELSE
  i = LEN_TRIM(expnt)
  READ(expnt(2:i), '(i4)') power10
END IF
IF (decpt > 0) power10 = power10 - (length + 1 - decpt)

! As 1.D+18 is represented exactly in 10-byte arithmetic,
! multiply by multiples of 1.D+18 or 1.D-18.

IF (power10 == 0) RETURN
IF (power10 < 0) THEN
  DO i = 1, -power10/18
    value = value / 1.E+18_r10
  END DO
  i = MOD(-power10, 18)
  IF (i /= 0) value = value / (10._r10 ** i)
ELSE
  DO i = 1, power10/18
    value = value * 1.E+18_r10
  END DO
  i = MOD(power10, 18)
  IF (i /= 0) value = value * (10._r10 ** i)
END IF

RETURN
END SUBROUTINE string_quad


SUBROUTINE quad_string(value, string, ier)
! Convert a quadruple-precision quantity to a decimal character string.
! Error indicator ier = 0 if conversion OK
!                     = 1 if the length of the string < 42 characters.

TYPE (quad), INTENT(IN)        :: value
CHARACTER (LEN=*), INTENT(OUT) :: string
INTEGER, INTENT(OUT)           :: ier

! Local variables
CHARACTER (LEN= 1)   :: sgn
CHARACTER (LEN=20)   :: str1, str2
TYPE (quad)          :: val
INTEGER              :: dec_expnt, i
REAL (r10)           :: tmp

IF (LEN(string) < 42) THEN
  ier = 1
  string = '***'
  RETURN
END IF
ier = 0

! Check if value = zero.
IF (value%hi == zero) THEN
  string = ' 0.00'
  RETURN
END IF

IF (value%hi < zero) THEN
  sgn = '-'
  val = - value
ELSE
  sgn = ' '
  val = value
END IF

! Use LOG10 to set the exponent.
dec_expnt = FLOOR( LOG10(val%hi) )

! Get the first 18 decimal digits
IF (dec_expnt /= 17) THEN
  val = val * EXP( LOG(quad(10._r10, zero)) * REAL(17 - dec_expnt) )
END IF
WRITE(str1, '(f19.0)') val%hi

! Calculate the remainder
READ(str1, '(f19.0)') tmp
val = val - tmp

! If val is -ve, subtract 1 from the last digit of str1, and add 1 to val.
IF (val%hi < -0.5D-19) THEN
  tmp = tmp - one
  WRITE(str1, '(f19.0)') tmp
  val = val + one
END IF
val = val * 1.E18_r10

! write the second 18 digits
WRITE(str2, '(f19.0)') val%hi

! If str2 consists of asterisks, add 1 in the last digit to str1.
! Set str2 to zeroes.
IF (str2(2:2) == '*') THEN
  tmp = tmp + one
  WRITE(str1, '(f20.0)') tmp
  IF (str1(1:1) /= ' ') THEN
    dec_expnt = dec_expnt + 1
  ELSE
    str1 = str1(2:20)
  END IF
  str2 = '000000000000000000.'
END IF

! Replace leading blanks with zeroes
DO i = 1, 18
  IF (str2(i:i) /= ' ') EXIT
  str2(i:i) = '0'
END DO

! Combine str1 & str2, removing decimal points & adding exponent.
i = INDEX(str1, '.')
str1(i:i) = ' '
str2(19:19) = ' '
string = '.' // TRIM(ADJUSTL(str1)) // TRIM(ADJUSTL(str2)) // 'E'
WRITE(str1, '(i4.2)') dec_expnt+1
string = TRIM(string) // ADJUSTL(str1)

! Restore the sign.
IF (sgn == '-') THEN
  string = '-' // ADJUSTL(string)
ELSE
  string = ADJUSTL(string)
END IF

RETURN
END SUBROUTINE quad_string


FUNCTION q_epsilon() RESULT(eps)
! Returns the machine accuracy.
! eps is the smallest value such that x.(1 + eps) > x for all x.
! This value is machine dependent.   It has been checked only for the
! the following PC compilers:
! Lahey's ELF90
! Lahey/Fujitsi LF95
! Compaq (Digital) DVF5

TYPE (quad) :: eps

eps = quad(0.702e-37_r10, 0.0_r10)

RETURN
END FUNCTION q_epsilon

END MODULE quadruple_precision

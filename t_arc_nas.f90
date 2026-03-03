PROGRAM test_inv_trig
!  1. a + b = asin [ sin(a).cos(b) + cos(a).sin(b) ]
!  2. a + b = acos [ cos(a).cos(b) - sin(a).sin(b) ]
!  3. a + b = atan [ tan(a) + tan(b) ] / [ 1 - tan(a).tan(b) ]

USE quadruple_precision
IMPLICIT NONE

TYPE (quad)          :: a, b, lhs, rhs, diff, sina, sinb, cosa, cosb,  &
                        tana, tanb, qone
REAL (r10)           :: half = 0.5_r10, small = EPSILON(half)
INTEGER, ALLOCATABLE :: seed(:)
INTEGER              :: k, i

qone = quad(1.0_r10, 0.0_r10)

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

DO i = 1, 10
  CALL RANDOM_NUMBER(a%hi)
  a%hi = (a%hi - half) / a%hi
  a%lo = a%hi * small
  CALL RANDOM_NUMBER(b%hi)
  b%hi = (b%hi - half) / b%hi
  b%lo = b%hi * small
                             ! asin [ sin(a).cos(b) + cos(a).sin(b) ]
                             ! Range (-pi/2, +pi/2)
  diff = (a + b)
  CALL longmodr(diff, pi, k, lhs)
  IF (k /= 2*(k/2)) lhs = -lhs
  sina = SIN(a)
  sinb = SIN(b)
  cosa = COS(a)
  cosb = COS(b)
  rhs = longasin( sina*cosb + cosa*sinb )
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! acos [ cos(a).cos(b) - sin(a).sin(b) ]
                             ! Range (0, pi)
  diff = (a + b)
  CALL longmodr(diff, pi, k, lhs)
  IF (k /= 2*(k/2)) THEN
    IF (lhs%hi > 0._r10) THEN
      lhs = pi - lhs
    ELSE
      lhs = pi + lhs
    END IF
  END IF
  IF (lhs%hi < 0._r10) lhs = -lhs
  rhs = longacos( cosa*cosb - sina*sinb )
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
                             ! atan [ tan(a) + tan(b) ] / [ 1 - tan(a).tan(b) ]
                             ! Range (-pi/2, +pi/2)
  diff = (a + b)
  CALL longmodr(diff, pi, k, lhs)
  tana = TAN(a)
  tanb = TAN(b)
  rhs = longatan( (tana + tanb) / (qone - tana*tanb) )
  diff = lhs - rhs
  WRITE(*, '(" lhs =", g13.5, "  Diff. =", g12.4)') lhs%hi, diff%hi
END DO

STOP
END PROGRAM test_inv_trig

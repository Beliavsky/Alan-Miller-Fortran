MODULE Lin_Feedback_Shift_Reg
! L'Ecuyer's 1999 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! This version is for 64-bit integers and assumes that KIND=8 identifies them.
! It has only been tested using Lahey's LF95 compiler
! http://users.bigpond.net.au/amiller/
! http://www.ozemail.com.au/~milleraj
! Latest revision - 12 January 2001

IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(14, 60)
! These are unsigned integers in the C version
INTEGER (KIND=8), SAVE  :: s1 = 153587801, s2 = -759022222, s3 = 1288503317, &
			   s4 = -1718083407, s5 = -123456789

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3, i4, i5)
IMPLICIT NONE

INTEGER (KIND=8), INTENT(IN) :: i1, i2, i3, i4, i5

s1 = i1
s2 = i2
s3 = i3
s4 = i4
s5 = i5
IF (IAND(s1,      -2) == 0) s1 = i1 - 8388607
IF (IAND(s2,    -512) == 0) s2 = i2 - 8388607
IF (IAND(s3,   -4096) == 0) s3 = i3 - 8388607
IF (IAND(s4, -131072) == 0) s4 = i4 - 8388607
IF (IAND(s5,-8388608) == 0) s5 = i5 - 8388607

RETURN
END SUBROUTINE init_seeds



FUNCTION lfsr258() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
! generators', Math. of Comput., 68, 261-269.

! The cycle length is claimed to be about 2^(258) or about 4.6 x 10^77.
! Actually - (2^63 - 1).(2^55 - 1).(2^52 - 1).(2^47 - 1).(2^41 - 1)

IMPLICIT NONE
REAL (dp) :: random_numb

INTEGER (KIND=8)   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,1), s1), -53)
s1 = IEOR( ISHFT( IAND(s1,-2), 10), b)
b  = ISHFT( IEOR( ISHFT(s2,24), s2), -50)
s2 = IEOR( ISHFT( IAND(s2,-512), 5), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -23)
s3 = IEOR( ISHFT( IAND(s3,-4096), 29), b)
b  = ISHFT( IEOR( ISHFT(s4,5), s4), -24)
s4 = IEOR( ISHFT( IAND(s4,-131072), 23), b)
b  = ISHFT( IEOR( ISHFT(s5,3), s5), -33)
s5 = IEOR( ISHFT( IAND(s5,-8388608), 8), b)

! The constant below is the reciprocal of (2^64 - 1)
random_numb = IEOR( IEOR( IEOR( IEOR(s1,s2), s3), s4), s5)  &
              * 5.4210108624275221E-20_dp + 0.5_dp

RETURN
END FUNCTION lfsr258

END MODULE Lin_Feedback_Shift_Reg



PROGRAM t_lfsr258
USE Lin_Feedback_Shift_Reg
IMPLICIT NONE
INTEGER (KIND=8)               :: i1, i2, i3, i4, i5
INTEGER (KIND=4)               :: i, j, k, l, m, freq(0:15,0:15,0:15,0:15)
INTEGER (KIND=8), ALLOCATABLE  :: seed(:)
REAL (dp)                      :: x, chi_sq, expctd, deg_freedom, stdev,  &
                                  lower, upper

WRITE(*, *)'Enter 5 integers as seeds: '
READ(*, *) i1, i2, i3, i4, i5
CALL init_seeds(i1, i2, i3, i4, i5)

x = lfsr258()
i = 16*x
x = lfsr258()
j = 16*x
x = lfsr258()
k = 16*x
freq = 0
DO m= 1, 10000000
  x = lfsr258()
  l = 16*x
  freq(i,j,k,l) = freq(i,j,k,l) + 1
  i = j
  j = k
  k = l
END DO

expctd = REAL( SUM(freq) ) / (16. * 16. * 16. * 16.)
chi_sq = 0.0_dp
DO i = 0, 15
  DO j = 0, 15
    DO k = 0, 15
      chi_sq = chi_sq + SUM( (REAL(freq(i,j,k,:)) - expctd)**2 ) / expctd
    END DO
  END DO
END DO

deg_freedom = (16. * 16. * 16. * 16.) - 1.
WRITE(*, '(a, f10.1, a, f6.0, a)') ' Chi-squared = ', chi_sq, ' with ',  &
                                  deg_freedom, ' deg. of freedom'

! Now repeat the exercise with the compiler's random number generator.

CALL RANDOM_SEED(size=k)
ALLOCATE( seed(k) )
IF (i1 /= 0) THEN
  seed(1) = i1
ELSE
  seed(1) = 1234567
END IF
IF (k >= 2) seed(2) = i2
IF (k >= 3) seed(3) = i3
DO i = 4, k
  seed(i) = seed(i-3)
END DO
CALL RANDOM_SEED(put=seed)

CALL RANDOM_NUMBER(x)
i = 16*x
CALL RANDOM_NUMBER(x)
j = 16*x
CALL RANDOM_NUMBER(x)
k = 16*x
freq = 0
DO m= 1, 10000000
  CALL RANDOM_NUMBER(x)
  l = 16*x
  freq(i,j,k,l) = freq(i,j,k,l) + 1
  i = j
  j = k
  k = l
END DO

expctd = REAL( SUM(freq) ) / (16. * 16. * 16. * 16.)
chi_sq = 0.0_dp
DO i = 0, 15
  DO j = 0, 15
    DO k = 0, 15
      chi_sq = chi_sq + SUM( (REAL(freq(i,j,k,:)) - expctd)**2 ) / expctd
    END DO
  END DO
END DO

WRITE(*, '(a, f10.1, a, f6.0, a)') ' Chi-squared = ', chi_sq, ' with ',  &
                                  deg_freedom, ' deg. of freedom'

! Calculate rough limits for chi-squared based upon a normal approximation.
stdev = SQRT(2. * deg_freedom)
upper = deg_freedom + 2.*stdev
lower = deg_freedom - 2.*stdev
WRITE(*, '(a, f8.1, a, f8.1)') ' Approx. 95% limits for chi-squared: ',  &
            lower, ' - ', upper

STOP
END PROGRAM t_lfsr258

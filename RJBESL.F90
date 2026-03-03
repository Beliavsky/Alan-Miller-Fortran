SUBROUTINE rjbesl(x, alpha, nb, b, ncalc)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!---------------------------------------------------------------------
! This routine calculates Bessel functions J sub(N+ALPHA) (X)
!   for non-negative argument X, and non-negative order N+ALPHA.


!  Explanation of variables in the calling sequence.

!   X     - working precision non-negative real argument for which
!           J's are to be calculated.
!   ALPHA - working precision fractional part of order for which
!           J's or exponentially scaled J'r (J*exp(X)) are
!           to be calculated.  0 <= ALPHA < 1.0.
!   NB  - integer number of functions to be calculated, NB > 0.
!           The first function calculated is of order ALPHA, and the
!           last is of order (NB - 1 + ALPHA).
!   B  - working precision output vector of length NB.  If RJBESL
!           terminates normally (NCALC=NB), the vector B contains the
!           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
!           corresponding exponentially scaled functions.
!   NCALC - integer output variable indicating possible errors.
!           Before using the vector B, the user should check that
!           NCALC=NB, i.e., all orders have been calculated to
!           the desired accuracy.  See Error Returns below.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   it     = Number of bits in the mantissa of a working precision variable
!   NSIG   = Decimal significance desired.  Should be set to
!            INT(LOG10(2)*it+1).  Setting NSIG lower will result
!            in decreased accuracy while setting NSIG higher will
!            increase CPU time without increasing accuracy.  The
!            truncation error is limited to a relative error of
!            T=.5*10**(-NSIG).

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   ENTEN  = 10.0 ** K, where K is the largest integer such that
!            ENTEN is machine-representable in working precision
!   ENSIG  = 10.0 ** NSIG
!   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
!            K .GE. NSIG/4
!   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
!   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
!            then at least N iterations of the backward recursion will be
!            executed.  The value of 10.0 ** 4 is used on every machine.

!     Approximate values for some important machines are:

!                            it    NSIG    ENTEN       ENSIG

!   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
!   Cyber 180/855
!     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
!   IEEE (IBM/XT,
!     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
!   IEEE (IBM/XT,
!     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
!   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
!   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
!   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
!   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16


!                           RTNSIG      ENMTEN      XLARGE

!   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
!   Cyber 180/855
!     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
!   IEEE (IBM/XT,
!     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
!   IEEE (IBM/XT,
!     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
!   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
!   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
!   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
!   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4

!*******************************************************************

!  Error returns

!    In case of an error,  NCALC .NE. NB, and not all J's are
!    calculated to the desired accuracy.

!    NCALC .LT. 0:  An argument is out of range. For example,
!       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
!       In this case, B(1) is set to zero, the remainder of the
!       B-vector is not calculated, and NCALC is set to
!       MIN(NB,0)-1 so that NCALC .NE. NB.

!    NB .GT. NCALC .GT. 0: Not all requested function values could
!       be calculated accurately.  This usually occurs because NB is
!       much larger than ABS(X).  In this case, B(N) is calculated
!       to the desired accuracy for N .LE. NCALC, but precision
!       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
!       for N .GT. NCALC (because it is too small to be represented),
!       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!       significant figures of B(N) can be trusted.

!  Intrinsic and other functions required are:

!     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,

!     REAL, SIN, SQRT


!  Acknowledgement

!   This program is based on a program written by David J. Sookne
!   (2) that computes values of the Bessel functions J or I of real
!   argument and integer order.  Modifications include the restriction
!   of the computation to the J Bessel function of non-negative real
!   argument, the extension of the computation to arbitrary positive
!   order, and the elimination of most underflow.

!  References: "A Note on Backward Recurrence Algorithms," Olver, F. W. J.,
!               and Sookne, D. J., Math. Comp. 26, 1972, pp 941-947.

!              "Bessel Functions of Real Argument and Integer Order,"
!               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp 125-132.

!  Latest modification: March 15, 1992

!  Author: W. J. Cody
!          Applied Mathematics Division
!          Argonne National Laboratory
!          Argonne, IL  60439, USA

!---------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: alpha
INTEGER, INTENT(IN)     :: nb
REAL (dp), INTENT(OUT)  :: b(:)
INTEGER, INTENT(OUT)    :: ncalc

INTEGER    :: i, j, k, l, m, magx, n, nbmx, nend, nstart
REAL (dp)  :: alpem, alp2em, capp, capq, em, en, gnu, halfx, p,  &
              plast, pold, psave, psavel, s, sum, t, t1, tempa, tempb,   &
              tempc, test, tover, xc, xin, xk, xm, vcos, vsin, z
!---------------------------------------------------------------------
!  Mathematical constants

!   PI2    - 2 / PI
!   TWOPI1 - first few significant digits of 2 * PI
!   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
!            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: pi2 = 0.636619772367581343075535_dp,  &
           twopi1 = 6.28125_dp, twopi2 = 1.935307179586476925286767D-3,  &
           zero = 0.0_dp, eighth = 0.125_dp, half = 0.5_dp, one = 1.0_dp,  &
           two = 2.0_dp, three = 3.0_dp, four = 4.0_dp, twofiv = 25.0_dp,  &
           one30 = 130.0_dp, three5 = 35.0_dp
!---------------------------------------------------------------------
!  Machine-dependent parameters
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: enten = 1.0D308, ensig = 1.0D16, rtnsig = 1.0D-4,  &
           enmten = 8.90D-308, xlarge = 1.0D4
!---------------------------------------------------------------------
!     Factorial(N)
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: FACT(25) = (/ 1.0_dp, 1.0_dp, 2.0_dp, 6.0_dp,   &
           24.0_dp, 1.2D2, 7.2D2, 5.04D3, 4.032D4, 3.6288D5, 3.6288D6,   &
           3.99168D7, 4.790016D8, 6.2270208D9, 8.71782912D10, 1.307674368D12,  &
           2.0922789888D13, 3.55687428096D14, 6.402373705728D15,   &
           1.21645100408832D17, 2.43290200817664D18, 5.109094217170944D19,  &
           1.12400072777760768D21, 2.585201673888497664D22,   &
           6.2044840173323943936D23 /)

INTERFACE
  FUNCTION GAMMA(X) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION GAMMA
END INTERFACE

!---------------------------------------------------------------------
! Check for out of range arguments.
!---------------------------------------------------------------------
magx = INT(x)
IF (nb > 0 .AND. x >= zero .AND. x <= xlarge .AND. alpha >= zero .AND.  &
    alpha < one) THEN
!---------------------------------------------------------------------
! Initialize result array to zero.
!---------------------------------------------------------------------
  ncalc = nb
  b(1:nb) = zero
!---------------------------------------------------------------------
! Branch to use 2-term ascending series for small X and asymptotic
! form for large X when NB is not too large.
!---------------------------------------------------------------------
  IF (x < rtnsig) THEN
!---------------------------------------------------------------------
! Two-term ascending series for small X.
!---------------------------------------------------------------------
    tempa = one
    alpem = one + alpha
    halfx = zero
    IF (x > enmten) halfx = half * x
    IF (alpha /= zero) tempa = halfx ** alpha / (alpha* gamma(alpha))
    tempb = zero
    IF ((x+one) > one) tempb = -halfx * halfx
    b(1) = tempa + tempa * tempb / alpem
    IF ((x /= zero) .AND. (b(1) == zero)) ncalc = 0
    IF (nb /= 1) THEN
      IF (x <= zero) THEN
        DO  n = 2, nb
          b(n) = zero
        END DO
      ELSE
!---------------------------------------------------------------------
! Calculate higher order functions.
!---------------------------------------------------------------------
        tempc = halfx
        tover = (enmten+enmten) / x
        IF (tempb /= zero) tover = enmten / tempb
        DO  n = 2, nb
          tempa = tempa / alpem
          alpem = alpem + one
          tempa = tempa * tempc
          IF (tempa <= tover*alpem) tempa = zero
          b(n) = tempa + tempa * tempb / alpem
          IF ((b(n) == zero) .AND. (ncalc > n)) ncalc = n - 1
        END DO
      END IF
    END IF
  ELSE IF ((x > twofiv) .AND. (nb <= magx+1)) THEN
!---------------------------------------------------------------------
! Asymptotic series for X .GT. 21.0.
!---------------------------------------------------------------------
    xc = SQRT(pi2/x)
    xin = (eighth/x) ** 2
    m = 11
    IF (x >= three5) m = 8
    IF (x >= one30) m = 4
    xm = four * m
!---------------------------------------------------------------------
! Argument reduction for SIN and COS routines.
!---------------------------------------------------------------------
    t = INT(x/(twopi1+twopi2)+half)
    z = ((x-t*twopi1)-t*twopi2) - (alpha+half) / pi2
    vsin = SIN(z)
    vcos = COS(z)
    gnu = alpha + alpha
    DO  i = 1, 2
      s = ((xm-one)-gnu) * ((xm-one)+gnu) * xin * half
      t = (gnu-(xm-three)) * (gnu+(xm-three))
      capp = s * t / fact(2*m+1)
      t1 = (gnu-(xm+one)) * (gnu+(xm+one))
      capq = s * t1 / fact(2*m+2)
      xk = xm
      k = m + m
      t1 = t
      DO  j = 2, m
        xk = xk - four
        s = ((xk-one)-gnu) * ((xk-one)+gnu)
        t = (gnu-(xk-three)) * (gnu+(xk-three))
        capp = (capp+one/fact(k-1)) * s * t * xin
        capq = (capq+one/fact(k)) * s * t1 * xin
        k = k - 2
        t1 = t
      END DO
      capp = capp + one
      capq = (capq+one) * (gnu*gnu-one) * (eighth/x)
      b(i) = xc * (capp*vcos-capq*vsin)
      IF (nb == 1) GO TO 180
      t = vsin
      vsin = -vcos
      vcos = t
      gnu = gnu + two
    END DO
!---------------------------------------------------------------------
! If  NB > 2, compute J(X,ORDER+I)  I = 2, NB-1
!---------------------------------------------------------------------
    IF (nb > 2) THEN
      gnu = alpha + alpha + two
      DO  j = 3, nb
        b(j) = gnu * b(j-1) / x - b(j-2)
        gnu = gnu + two
      END DO
    END IF
!---------------------------------------------------------------------
! Use recurrence to generate results.
! First initialize the calculation of P*S.
!---------------------------------------------------------------------
  ELSE
    nbmx = nb - magx
    n = magx + 1
    en = n+n + (alpha+alpha)
    plast = one
    p = en / x
!---------------------------------------------------------------------
! Calculate general significance test.
!---------------------------------------------------------------------
    test = ensig + ensig
    IF (nbmx >= 3) THEN
!---------------------------------------------------------------------
! Calculate P*S until N = NB-1.  Check for possible overflow.
!---------------------------------------------------------------------
      tover = enten / ensig
      nstart = magx + 2
      nend = nb - 1
      en = nstart + nstart - two + (alpha+alpha)
      DO  k = nstart, nend
        n = k
        en = en + two
        pold = plast
        plast = p
        p = en * plast / x - pold
        IF (p > tover) THEN
!---------------------------------------------------------------------
! To avoid overflow, divide P*S by TOVER.  Calculate P*S until ABS(P) > 1.
!---------------------------------------------------------------------
          tover = enten
          p = p / tover
          plast = plast / tover
          psave = p
          psavel = plast
          nstart = n + 1
          70 n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x - pold
          IF (p <= one) GO TO 70
          tempb = en / x
!---------------------------------------------------------------------
! Calculate backward test and find NCALC, the highest N such that
! the test is passed.
!---------------------------------------------------------------------
          test = pold * plast * (half-half/(tempb*tempb))
          test = test / ensig
          p = plast * tover
          n = n - 1
          en = en - two
          nend = MIN(nb,n)
          DO  l = nstart, nend
            pold = psavel
            psavel = psave
            psave = en * psavel / x - pold
            IF (psave*psavel > test) THEN
              ncalc = l - 1
              GO TO 110
            END IF
          END DO
          ncalc = nend
          GO TO 110
        END IF
      END DO
      n = nend
      en = n+n + (alpha+alpha)
!---------------------------------------------------------------------
! Calculate special significance test for NBMX > 2.
!---------------------------------------------------------------------
      test = MAX(test, SQRT(plast*ensig)*SQRT(p+p))
    END IF
!---------------------------------------------------------------------
! Calculate P*S until significance test passes.
!---------------------------------------------------------------------
    100 n = n + 1
    en = en + two
    pold = plast
    plast = p
    p = en * plast / x - pold
    IF (p < test) GO TO 100
!---------------------------------------------------------------------
! Initialize the backward recursion and the normalization sum.
!---------------------------------------------------------------------
    110 n = n + 1
    en = en + two
    tempb = zero
    tempa = one / p
    m = 2 * n - 4 * (n/2)
    sum = zero
    em = n/2
    alpem = (em-one) + alpha
    alp2em = (em+em) + alpha
    IF (m /= 0) sum = tempa * alpem * alp2em / em
    nend = n - nb
    IF (nend > 0) THEN
!---------------------------------------------------------------------
! Recur backward via difference equation, calculating (but not
! storing) B(N), until N = NB.
!---------------------------------------------------------------------
      DO  l = 1, nend
        n = n - 1
        en = en - two
        tempc = tempb
        tempb = tempa
        tempa = (en*tempb) / x - tempc
        m = 2 - m
        IF (m /= 0) THEN
          em = em - one
          alp2em = (em+em) + alpha
          IF (n == 1) GO TO 130
          alpem = (em-one) + alpha
          IF (alpem == zero) alpem = one
          sum = (sum+tempa*alp2em) * alpem / em
        END IF
      END DO
    END IF
!---------------------------------------------------------------------
! Store B(NB).
!---------------------------------------------------------------------
    130 b(n) = tempa
    IF (nend >= 0) THEN
      IF (nb <= 1) THEN
        alp2em = alpha
        IF ((alpha+one) == one) alp2em = one
        sum = sum + b(1) * alp2em
        GO TO 160
      ELSE
!---------------------------------------------------------------------
! Calculate and store B(NB-1).
!---------------------------------------------------------------------
        n = n - 1
        en = en - two
        b(n) = (en*tempa) / x - tempb
        IF (n == 1) GO TO 150
        m = 2 - m
        IF (m /= 0) THEN
          em = em - one
          alp2em = (em+em) + alpha
          alpem = (em-one) + alpha
          IF (alpem == zero) alpem = one
          sum = (sum+b(n)*alp2em) * alpem / em
        END IF
      END IF
    END IF
    nend = n - 2
    IF (nend /= 0) THEN
!---------------------------------------------------------------------
! Calculate via difference equation and store B(N), until N = 2.
!---------------------------------------------------------------------
      DO  l = 1, nend
        n = n - 1
        en = en - two
        b(n) = (en*b(n+1)) / x - b(n+2)
        m = 2 - m
        IF (m /= 0) THEN
          em = em - one
          alp2em = (em+em) + alpha
          alpem = (em-one) + alpha
          IF (alpem == zero) alpem = one
          sum = (sum+b(n)*alp2em) * alpem / em
        END IF
      END DO
    END IF
!---------------------------------------------------------------------
! Calculate B(1).
!---------------------------------------------------------------------
    b(1) = two * (alpha+one) * b(2) / x - b(3)
    150 em = em - one
    alp2em = (em+em) + alpha
    IF (alp2em == zero) alp2em = one
    sum = sum + b(1) * alp2em
!---------------------------------------------------------------------
! Normalize.  Divide all B(N) by sum.
!---------------------------------------------------------------------
    160 IF ((alpha+one) /= one) sum = sum * gamma(alpha) * (x*half) ** (-alpha)
    tempa = enmten
    IF (sum > one) tempa = tempa * sum
    DO  n = 1, nb
      IF (ABS(b(n)) < tempa) b(n) = zero
      b(n) = b(n) / sum
    END DO
  END IF
!---------------------------------------------------------------------
! Error return -- X, NB, or ALPHA is out of range.
!---------------------------------------------------------------------
ELSE
  b(1) = zero
  ncalc = MIN(nb,0) - 1
END IF
!---------------------------------------------------------------------
! Exit
!---------------------------------------------------------------------
180 RETURN
! ---------- Last line of RJBESL ----------
END SUBROUTINE rjbesl

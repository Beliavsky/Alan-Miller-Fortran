SUBROUTINE ribesl(x, alpha, nb, ize, b, ncalc)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!-------------------------------------------------------------------

!  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
!  for non-negative argument X, and non-negative order N+ALPHA,
!  with or without exponential scaling.

! Explanation of variables in the calling sequence

! X     - Working precision non-negative real argument for which
!         I's or exponentially scaled I's (I*EXP(-X))
!         are to be calculated.  If I's are to be calculated,
!         X must be less than EXPARG (see below).
! ALPHA - Working precision fractional part of order for which I's or
!         exponentially scaled I's (I*EXP(-X)) are to be calculated.
!         0 <= ALPHA < 1.0.
! NB    - Integer number of functions to be calculated, NB > 0.
!         The first function calculated is of order ALPHA, and the
!         last is of order (NB - 1 + ALPHA).
! IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
!         and 2 if exponentially scaled I's are to be calculated.
! B     - Working precision output vector of length NB.  If the routine
!         terminates normally (NCALC=NB), the vector B contains the
!         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
!         corresponding exponentially scaled functions.
! NCALC - Integer output variable indicating possible errors.
!         Before using the vector B, the user should check that
!         NCALC=NB, i.e., all orders have been calculated to
!         the desired accuracy.  See error returns below.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   it     = Number of bits in the mantissa of a working precision variable

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   NSIG   = Decimal significance desired.  Should be set to
!            INT(LOG10(2)*it+1).  Setting NSIG lower will result
!            in decreased accuracy while setting NSIG higher will
!            increase CPU time without increasing accuracy.  The
!            truncation error is limited to a relative error of
!            T=.5*10**(-NSIG).
!   ENTEN  = 10.0 ** K, where K is the largest integer such that
!            ENTEN is machine-representable in working precision
!   ENSIG  = 10.0 ** NSIG
!   RTNSIG = 10.0 ** (-K) for the smallest integer K such that K >= NSIG/4
!   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
!   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
!            in mind that if ABS(X)=N, then at least N iterations
!            of the backward recursion will be executed.  The value
!            of 10.0 ** 4 is used on every machine.
!   EXPARG = Largest working precision argument that the library EXP routine
!            EXP routine can handle and upper limit on the magnitude of X
!            when IZE=1; approximately LOG(beta**maxexp)


!     Approximate values for some important machines are:

!                        beta       minexp      maxexp       it

!  CRAY-1        (S.P.)    2        -8193        8191        48
!  Cyber 180/855
!    under NOS   (S.P.)    2         -975        1070        48
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2         -126         128        24
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2        -1022        1024        53
!  IBM 3033      (D.P.)   16          -65          63        14
!  VAX           (S.P.)    2         -128         127        24
!  VAX D-Format  (D.P.)    2         -128         127        56
!  VAX G-Format  (D.P.)    2        -1024        1023        53


!                        NSIG       ENTEN       ENSIG      RTNSIG

! CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
! Cyber 180/855
!   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
! IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
! VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
! VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4


!                         ENMTEN      XLARGE   EXPARG

! CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
! Cyber 180/855
!   under NOS   (S.P.)   1.25E-293    1.0E+4     741
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
! IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
! VAX           (S.P.)   1.17E-38     1.0E+4      88
! VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
! VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709

!*******************************************************************

! Error returns

!  In case of an error,  NCALC .NE. NB, and not all I's are
!  calculated to the desired accuracy.

!  NCALC < 0:  An argument is out of range. For example, NB <= 0,
!     IZE is not 1 or 2, or IZE=1 and ABS(X) >= EXPARG.
!     In this case, the B-vector is not calculated, and NCALC is
!     set to MIN(NB,0)-1 so that NCALC .NE. NB.

!  NB > NCALC > 0: Not all requested function values could be calculated
!     accurately.  This usually occurs because NB is much larger than
!     ABS(X).  In this case, B(N) is calculated to the desired accuracy for
!     N <= NCALC, but precision is lost for NCALC < N <= NB.
!     If B(N) does not vanish for N > NCALC (because it is too small to be
!     represented), and B(N)/B(NCALC) = 10**(-K), then only the first
!     NSIG-K significant figures of B(N) can be trusted.

! Intrinsic functions required are:

!     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT


! Acknowledgement

!  This program is based on a program written by David J. Sookne (2) that
!  computes values of the Bessel functions J or I of real argument and
!  integer order.  Modifications include the restriction of the computation
!  to the I Bessel function of non-negative real argument, the extension of
!  the computation to arbitrary positive order, the inclusion of optional
!  exponential scaling, and the elimination of most underflow.
!  An earlier version was published in (3).

! References: "A Note on Backward Recurrence Algorithms," Olver, F. W. J.,
!              and Sookne, D. J., Math. Comp. 26, 1972, pp 941-947.

!             "Bessel Functions of Real Argument and Integer Order,"
!              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp 125-132.

!             "ALGORITHM 597, Sequence of Modified Bessel Functions of the
!              First Kind," Cody, W. J., Trans. Math. Soft., 1983, pp. 242-245.

!  Latest modification: March 14, 1992

!  Modified by: W. J. Cody and L. Stoltz
!               Applied Mathematics Division
!               Argonne National Laboratory
!               Argonne, IL  60439, USA

!-------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: alpha
INTEGER, INTENT(IN)     :: nb
INTEGER, INTENT(IN)     :: ize
REAL (dp), INTENT(OUT)  :: b(:)
INTEGER, INTENT(OUT)    :: ncalc

INTEGER    :: k, l, magx, n, nbmx, nend, nstart
REAL (dp)  :: em, empal, emp2al, en, halfx, p, plast, pold, &
              psave, psavel, sum, tempa, tempb, tempc, test, tover

!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: one = 1.0_dp, two = 2.0_dp, zero = 0.0_dp,  &
                         half = 0.5_dp, const = 1.585_dp
!-------------------------------------------------------------------
!  Machine-dependent parameters
!-------------------------------------------------------------------
INTEGER, PARAMETER    :: nsig = 16
REAL (dp), PARAMETER  :: xlarge = 1.0D4, exparg = 709.0_dp, enten = 1.0D308, &
                         ensig = 1.0D16, rtnsig = 1.0D-4, ENMTEN = 8.9D-308

INTERFACE
  FUNCTION GAMMA(X) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION GAMMA
END INTERFACE

!-------------------------------------------------------------------
! Check for X, NB, OR IZE out of range.
!-------------------------------------------------------------------
IF (nb > 0 .AND. x >= zero .AND. alpha >= zero .AND. alpha < one .AND.  &
    (ize == 1 .AND. x <= exparg) .OR. (ize == 2 .AND. x <= xlarge)) THEN
!-------------------------------------------------------------------
! Use 2-term ascending series for small X
!-------------------------------------------------------------------
  ncalc = nb
  magx = INT(x)
  IF (x >= rtnsig) THEN
!-------------------------------------------------------------------
! Initialize the forward sweep, the P-sequence of Olver
!-------------------------------------------------------------------
    nbmx = nb - magx
    n = magx + 1
    en = (n+n) + (alpha+alpha)
    plast = one
    p = en / x
!-------------------------------------------------------------------
! Calculate general significance test
!-------------------------------------------------------------------
    test = ensig + ensig
    IF (2*magx > 5*nsig) THEN
      test = SQRT(test*p)
    ELSE
      test = test / const ** magx
    END IF
    IF (nbmx >= 3) THEN
!-------------------------------------------------------------------
! Calculate P-sequence until N = NB-1.  Check for possible overflow.
!-------------------------------------------------------------------
      tover = enten / ensig
      nstart = magx + 2
      nend = nb - 1
      DO  k = nstart, nend
        n = k
        en = en + two
        pold = plast
        plast = p
        p = en * plast / x + pold
        IF (p > tover) THEN
!-------------------------------------------------------------------
! To avoid overflow, divide P-sequence by TOVER.  Calculate
! P-sequence until ABS(P) > 1.
!-------------------------------------------------------------------
          tover = enten
          p = p / tover
          plast = plast / tover
          psave = p
          psavel = plast
          nstart = n + 1
          10 n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x + pold
          IF (p <= one) GO TO 10
          tempb = en / x
!-------------------------------------------------------------------
! Calculate backward test, and find NCALC, the highest N
! such that the test is passed.
!-------------------------------------------------------------------
          test = pold * plast / ensig
          test = test * (half - half/(tempb*tempb))
          p = plast * tover
          n = n - 1
          en = en - two
          nend = MIN(nb, n)
          DO  l = nstart, nend
            ncalc = l
            pold = psavel
            psavel = psave
            psave = en * psavel / x + pold
            IF (psave*psavel > test) GO TO 30
          END DO
          ncalc = nend + 1
          30 ncalc = ncalc - 1
          GO TO 60
        END IF
      END DO
      n = nend
      en = (n+n) + (alpha+alpha)
!-------------------------------------------------------------------
! Calculate special significance test for NBMX > 2.
!-------------------------------------------------------------------
      test = MAX(test,SQRT(plast*ensig)*SQRT(p+p))
    END IF
!-------------------------------------------------------------------
! Calculate P-sequence until significance test passed.
!-------------------------------------------------------------------
    50 n = n + 1
    en = en + two
    pold = plast
    plast = p
    p = en * plast / x + pold
    IF (p < test) GO TO 50
!-------------------------------------------------------------------
! Initialize the backward recursion and the normalization sum.
!-------------------------------------------------------------------
    60 n = n + 1
    en = en + two
    tempb = zero
    tempa = one / p
    em = n - one
    empal = em + alpha
    emp2al = (em-one) + (alpha+alpha)
    sum = tempa * empal * emp2al / em
    nend = n - nb
    IF (nend < 0) THEN
!-------------------------------------------------------------------
! N < NB, so store B(N) and set higher orders to zero.
!-------------------------------------------------------------------
      b(n) = tempa
      nend = -nend
      DO  l = 1, nend
        b(n+l) = zero
      END DO
    ELSE
      IF (nend > 0) THEN
!-------------------------------------------------------------------
! Recur backward via difference equation, calculating (but
! not storing) B(N), until N = NB.
!-------------------------------------------------------------------
        DO  l = 1, nend
          n = n - 1
          en = en - two
          tempc = tempb
          tempb = tempa
          tempa = (en*tempb) / x + tempc
          em = em - one
          emp2al = emp2al - one
          IF (n == 1) GO TO 90
          IF (n == 2) emp2al = one
          empal = empal - one
          sum = (sum + tempa*empal) * emp2al / em
        END DO
      END IF
!-------------------------------------------------------------------
! Store B(NB)
!-------------------------------------------------------------------
      90 b(n) = tempa
      IF (nb <= 1) THEN
        sum = (sum+sum) + tempa
        GO TO 120
      END IF
!-------------------------------------------------------------------
! Calculate and Store B(NB-1)
!-------------------------------------------------------------------
      n = n - 1
      en = en - two
      b(n) = (en*tempa) / x + tempb
      IF (n == 1) GO TO 110
      em = em - one
      emp2al = emp2al - one
      IF (n == 2) emp2al = one
      empal = empal - one
      sum = (sum+b(n)*empal) * emp2al / em
    END IF
    nend = n - 2
    IF (nend > 0) THEN
!-------------------------------------------------------------------
! Calculate via difference equation and store B(N), until N = 2.
!-------------------------------------------------------------------
      DO  l = 1, nend
        n = n - 1
        en = en - two
        b(n) = (en*b(n+1)) / x + b(n+2)
        em = em - one
        emp2al = emp2al - one
        IF (n == 2) emp2al = one
        empal = empal - one
        sum = (sum+b(n)*empal) * emp2al / em
      END DO
    END IF
!-------------------------------------------------------------------
! Calculate B(1)
!-------------------------------------------------------------------
    b(1) = two * empal * b(2) / x + b(3)
    110 sum = (sum+sum) + b(1)
!-------------------------------------------------------------------
! Normalize.  Divide all B(N) by sum.
!-------------------------------------------------------------------
    120 IF (alpha /= zero) sum = sum * gamma(one+alpha) * (x*half) ** (-alpha)
    IF (ize == 1) sum = sum * EXP(-x)
    tempa = enmten
    IF (sum > one) tempa = tempa * sum
    DO  n = 1, nb
      IF (b(n) < tempa) b(n) = zero
      b(n) = b(n) / sum
    END DO
    RETURN
!-------------------------------------------------------------------
! Two-term ascending series for small X.
!-------------------------------------------------------------------
  ELSE
    tempa = one
    empal = one + alpha
    halfx = zero
    IF (x > enmten) halfx = half * x
    IF (alpha /= zero) tempa = halfx ** alpha / gamma(empal)
    IF (ize == 2) tempa = tempa * EXP(-x)
    tempb = zero
    IF ((x+one) > one) tempb = halfx * halfx
    b(1) = tempa + tempa * tempb / empal
    IF (x /= zero .AND. b(1) == zero) ncalc = 0
    IF (nb > 1) THEN
      IF (x == zero) THEN
        DO  n = 2, nb
          b(n) = zero
        END DO
      ELSE
!-------------------------------------------------------------------
! Calculate higher-order functions.
!-------------------------------------------------------------------
        tempc = halfx
        tover = (enmten+enmten) / x
        IF (tempb /= zero) tover = enmten / tempb
        DO  n = 2, nb
          tempa = tempa / empal
          empal = empal + one
          tempa = tempa * tempc
          IF (tempa <= tover*empal) tempa = zero
          b(n) = tempa + tempa * tempb / empal
          IF (b(n) == zero .AND. ncalc > n) ncalc = n - 1
        END DO
      END IF
    END IF
  END IF
ELSE
  ncalc = MIN(nb,0) - 1
END IF
RETURN
!---------- Last line of RIBESL ----------
END SUBROUTINE ribesl

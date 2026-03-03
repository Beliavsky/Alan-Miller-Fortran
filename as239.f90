FUNCTION gammad(x, p) RESULT(fn_val)

!      ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!      Computation of the Incomplete Gamma Integral

!      Auxiliary functions required: ALOGAM = logarithm of the gamma
!      function, and ALNORM = algorithm AS66

! ELF90-compatible version by Alan Miller
! Latest revision - 27 October 2000

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN) :: x, p
REAL (dp)             :: fn_val

! Local variables
REAL (dp)             :: pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an
REAL (dp), PARAMETER  :: zero = 0.d0, one = 1.d0, two = 2.d0, &
                         oflo = 1.d+37, three = 3.d0, nine = 9.d0, &
                         tol = 1.d-14, xbig = 1.d+8, plimit = 1000.d0, &
                         elimit = -88.d0
! EXTERNAL alogam, alnorm

fn_val = zero

!      Check that we have valid values for X and P

IF (p <= zero .OR. x < zero) THEN
  WRITE(*, *) 'AS239: Either p <= 0 or x < 0'
  RETURN
END IF
IF (x == zero) RETURN

!      Use a normal approximation if P > PLIMIT

IF (p > plimit) THEN
  pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one /(nine * p) - one)
  fn_val = alnorm(pn1, .false.)
  RETURN
END IF

!      If X is extremely large compared to P then set fn_val = 1

IF (x > xbig) THEN
  fn_val = one
  RETURN
END IF

IF (x <= one .OR. x < p) THEN
  
!      Use Pearson's series expansion.
!      (Note that P is not large enough to force overflow in ALOGAM).
!      No need to test IFAULT on exit since P > 0.
  
  arg = p * LOG(x) - x - alogam(p + one, ifault)
  c = one
  fn_val = one
  a = p
  40   a = a + one
  c = c * x / a
  fn_val = fn_val + c
  IF (c > tol) GO TO 40
  arg = arg + LOG(fn_val)
  fn_val = zero
  IF (arg >= elimit) fn_val = EXP(arg)
  
ELSE
  
!      Use a continued fraction expansion
  
  arg = p * LOG(x) - x - alogam(p, ifault)
  a = one - p
  b = a + x + one
  c = zero
  pn1 = one
  pn2 = x
  pn3 = x + one
  pn4 = x * b
  fn_val = pn3 / pn4
  60   a = a + one
  b = b + two
  c = c + one
  an = a * c
  pn5 = b * pn3 - an * pn1
  pn6 = b * pn4 - an * pn2
  IF (ABS(pn6) > zero) THEN
    rn = pn5 / pn6
    IF (ABS(fn_val - rn) <= MIN(tol, tol * rn)) GO TO 80
    fn_val = rn
  END IF
  
  pn1 = pn3
  pn2 = pn4
  pn3 = pn5
  pn4 = pn6
  IF (ABS(pn5) >= oflo) THEN
    
!      Re-scale terms in continued fraction if terms are large
    
    pn1 = pn1 / oflo
    pn2 = pn2 / oflo
    pn3 = pn3 / oflo
    pn4 = pn4 / oflo
  END IF
  GO TO 60
  80   arg = arg + LOG(fn_val)
  fn_val = one
  IF (arg >= elimit) fn_val = one - EXP(arg)
END IF

RETURN
END FUNCTION gammad

PROGRAM test395

!     Test translated versions of CACM algorithms 395 & 396

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp)          :: ndf, prob, t

INTERFACE
  FUNCTION student(t, ndf, normal) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: t, ndf
    REAL (dp)             :: fn_val
    INTERFACE
      FUNCTION normal(x) RESULT(fx)
        IMPLICIT NONE
        INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
        REAL (dp), INTENT(IN) :: x
        REAL (dp)             :: fx
      END FUNCTION normal
    END INTERFACE
  END FUNCTION student

  FUNCTION t_quantile(p, ndf, normdev) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: p, ndf
    REAL (dp)             :: fn_val
    INTERFACE
      FUNCTION normdev(p) RESULT(x)
        IMPLICIT NONE
        INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
        REAL (dp), INTENT(IN) :: p
        REAL (dp)             :: x
      END FUNCTION normdev
    END INTERFACE
  END FUNCTION t_quantile

  FUNCTION normal(x) RESULT(fx)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fx
  END FUNCTION normal

  FUNCTION normdev(p) RESULT(x)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: p
    REAL (dp)             :: x
  END FUNCTION normdev
END INTERFACE

10 WRITE(*, *) 'Enter no. of deg. of freedom & t-value: '
READ(*, *) ndf, t
prob = student(t, ndf, normal)
WRITE(*, 900) prob
900 FORMAT(' Prob = ', g15.7)
t = t_quantile(prob, ndf, normdev)
WRITE(*, 910) t
910 FORMAT(' t-value = ', f10.6/)
GO TO 10

END PROGRAM test395



FUNCTION normal(z) RESULT(prob)
IMPLICIT NONE

INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN) :: z
REAL (dp)             :: prob

CALL nprob(z, p=prob)

RETURN

CONTAINS


SUBROUTINE nprob(z, p, q, pdf)
! For a given number (z) of standard deviations from the mean, the
! probabilities to the left (p) and right (q) of z are calculated,
! and the probability density (pdf).

! Programmer - Alan J. Miller
! Latest revision - 7 August 1997
! This version is compatible with the ELF90 subset of Fortran 90.

! REFERENCE: ADAMS,A.G. AREAS UNDER THE NORMAL CURVE,
! ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969.

IMPLICIT NONE

INTEGER, PARAMETER :: doubleprec = SELECTED_REAL_KIND(15, 60)
REAL(doubleprec), INTENT(IN)            :: z
REAL(doubleprec), INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL(doubleprec) :: a0 = 0.5D0, a1 = 0.398942280444D0, a2 = 0.399903438504D0, &
                    a3 = 5.75885480458D0, a4 = 29.8213557808D0,  &
                    a5 = 2.62433121679D0, a6 = 48.6959930692D0,  &
                    a7 = 5.92885724438D0, b0 = 0.398942280385D0, &
                    b1 = 3.8052D-8, b2 = 1.00000615302D0, b3 = 3.98064794D-4, &
                    b4 = 1.98615381364D0, b5 = 0.151679116635D0, &
                    b6 = 5.29330324926D0, b7 = 4.8385912808D0,   &
                    b8 = 15.1508972451D0, b9 = 0.742380924027D0, &
                    b10 = 30.789933034D0, b11 = 3.99019417011D0, &
                    zero = 0.D0, one = 1.D0, zabs, y, pp, qq, ppdf


zabs = ABS(z)
IF (zabs < 12.7D0) THEN
  y = a0*z*z
  ppdf = EXP(-y)*b0
  IF (PRESENT(pdf)) pdf = ppdf

!     Z BETWEEN -1.28 AND +1.28

  IF (zabs < 1.28) THEN
    qq = a0 - zabs*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
    IF(z < zero) THEN
      pp = qq
      qq = one - pp
    ELSE
      pp = one - qq
    END IF
    IF (PRESENT(p)) p = pp
    IF (PRESENT(q)) q = qq
    RETURN
  END IF

!     ZABS BETWEEN 1.28 AND 12.7

  qq = ppdf/(zabs-b1+b2/(zabs+b3+b4/(zabs-b5+b6/(zabs+b7-b8/ &
            (zabs+b9+b10/(zabs+b11))))))
  IF(z < zero) THEN
    pp = qq
    qq = one - pp
  ELSE
    pp = one - qq
  END IF
  IF (PRESENT(p)) p = pp
  IF (PRESENT(q)) q = qq
  RETURN

!     Z FAR OUT IN TAIL

ELSE
  ppdf = zero
  IF(z < zero) THEN
    pp = zero
    qq = one
  ELSE
    pp = one
    qq = zero
  END IF
  IF (PRESENT(p)) p = pp
  IF (PRESENT(q)) q = qq
  IF (PRESENT(pdf)) pdf = ppdf
  RETURN
END IF

RETURN
END SUBROUTINE nprob

END FUNCTION normal



FUNCTION normdev(p) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN) :: p
REAL (dp)             :: fn_val

! Local variable
INTEGER :: ifault

CALL ppnd16(p, fn_val, ifault)
IF (ifault /= 0) WRITE(*, *) 'Error in ppnd16: ifault =', ifault

RETURN

CONTAINS


SUBROUTINE ppnd16 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE

INTEGER, PARAMETER            :: doub_prec = SELECTED_REAL_KIND(15, 60)
REAL (doub_prec), INTENT(IN)  :: p
INTEGER, INTENT(OUT)          :: ifault
REAL (doub_prec), INTENT(OUT) :: normal_dev

! Local variables

REAL (doub_prec) :: zero = 0.d0, one = 1.d0, half = 0.5d0,  &
                    split1 = 0.425d0, split2 = 5.d0, const1 = 0.180625d0, &
                    const2 = 1.6d0, q, r

! Coefficients for P close to 0.5

REAL (doub_prec) :: a0 = 3.3871328727963666080D0, &
                    a1 = 1.3314166789178437745D+2, &
                    a2 = 1.9715909503065514427D+3, &
                    a3 = 1.3731693765509461125D+4, &
                    a4 = 4.5921953931549871457D+4, &
                    a5 = 6.7265770927008700853D+4, &
                    a6 = 3.3430575583588128105D+4, &
                    a7 = 2.5090809287301226727D+3, &
                    b1 = 4.2313330701600911252D+1, &
                    b2 = 6.8718700749205790830D+2, &
                    b3 = 5.3941960214247511077D+3, &
                    b4 = 2.1213794301586595867D+4, &
                    b5 = 3.9307895800092710610D+4, &
                    b6 = 2.8729085735721942674D+4, &
                    b7 = 5.2264952788528545610D+3
! HASH SUM AB           55.8831928806149014439

! Coefficients for P not close to 0, 0.5 or 1.

REAL (doub_prec) :: c0 = 1.42343711074968357734D0, &
                    c1 = 4.63033784615654529590D0, &
                    c2 = 5.76949722146069140550D0, &
                    c3 = 3.64784832476320460504D0, &
                    c4 = 1.27045825245236838258D0, &
                    c5 = 2.41780725177450611770D-1, &
                    c6 = 2.27238449892691845833D-2, &
                    c7 = 7.74545014278341407640D-4, &
                    d1 = 2.05319162663775882187D0, &
                    d2 = 1.67638483018380384940D0, &
                    d3 = 6.89767334985100004550D-1, &
                    d4 = 1.48103976427480074590D-1, &
                    d5 = 1.51986665636164571966D-2, &
                    d6 = 5.47593808499534494600D-4, &
                    d7 = 1.05075007164441684324D-9
! HASH SUM CD           49.33206503301610289036

! Coefficients for P near 0 or 1.

REAL (doub_prec) :: e0 = 6.65790464350110377720D0, &
                    e1 = 5.46378491116411436990D0, &
                    e2 = 1.78482653991729133580D0, &
                    e3 = 2.96560571828504891230D-1, &
                    e4 = 2.65321895265761230930D-2, &
                    e5 = 1.24266094738807843860D-3, &
                    e6 = 2.71155556874348757815D-5, &
                    e7 = 2.01033439929228813265D-7, &
                    f1 = 5.99832206555887937690D-1, &
                    f2 = 1.36929880922735805310D-1, &
                    f3 = 1.48753612908506148525D-2, &
                    f4 = 7.86869131145613259100D-4, &
                    f5 = 1.84631831751005468180D-5, &
                    f6 = 1.42151175831644588870D-7, &
                    f7 = 2.04426310338993978564D-15
! HASH SUM EF           47.52583317549289671629

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
  ELSE
    r = r - split2
    normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF

RETURN
END SUBROUTINE ppnd16

END FUNCTION normdev

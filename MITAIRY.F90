MODULE itairy_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
 
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE itairy(x, apt, bpt, ant, bnt)

!    ======================================================
!    Purpose: Compute the integrals of Airy fnctions with
!             respect to t from 0 and x ( x ò 0 )
!    Input  : x   --- Upper limit of the integral
!    Output : APT --- Integration of Ai(t) from 0 and x
!             BPT --- Integration of Bi(t) from 0 and x
!             ANT --- Integration of Ai(-t) from 0 and x
!             BNT --- Integration of Bi(-t) from 0 and x
!    ======================================================

REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(OUT)     :: apt
REAL (dp), INTENT(OUT)     :: bpt
REAL (dp), INTENT(OUT)     :: ant
REAL (dp), INTENT(OUT)     :: bnt

REAL (dp), PARAMETER  :: a(16) = (/   &
        .569444444444444_dp, .891300154320988_dp,  &
        .226624344493027D+01, .798950124766861D+01,  &
        .360688546785343D+02, .198670292131169D+03,  &
        .129223456582211D+04, .969483869669600D+04,  &
        .824184704952483D+05, .783031092490225D+06,  &
        .822210493622814D+07, .945557399360556D+08,  &
        .118195595640730D+10, .159564653040121D+11,  &
        .231369166433050D+12, .358622522796969D+13 /)

REAL (dp)  :: eps, pi, c1, c2, sr3, fx, gx, r, q0, q1, q2, xe, xp6, xr1,  &
              xr2, su1, su2, su3, su4, su5, su6
INTEGER    :: k, l

eps = 1.0D-15
pi = 3.141592653589793_dp
c1 = .355028053887817_dp
c2 = .258819403792807_dp
sr3 = 1.732050807568877_dp
IF (x == 0.0_dp) THEN
  apt = 0.0_dp
  bpt = 0.0_dp
  ant = 0.0_dp
  bnt = 0.0_dp
ELSE
  IF (ABS(x) <= 9.25_dp) THEN
    DO  l = 0, 1
      x = (-1) ** l * x
      fx = x
      r = x
      DO  k = 1, 40
        r = r * (3.0*k - 2.0_dp) / (3.0*k + 1.0_dp) * x / (3.0*k) * x /  &
            (3.0*k - 1.0_dp) * x
        fx = fx + r
        IF (ABS(r) < ABS(fx)*eps) EXIT
      END DO

      gx = .5_dp * x * x
      r = gx
      DO  k = 1, 40
        r = r * (3.0*k-1.0_dp) / (3.0*k+2.0_dp) * x / (3.0*k) * x /  &
            (3.0*k + 1.0_dp) * x
        gx = gx + r
        IF (ABS(r) < ABS(gx)*eps) EXIT
      END DO

      ant = c1 * fx - c2 * gx
      bnt = sr3 * (c1*fx + c2*gx)
      IF (l == 0) THEN
        apt = ant
        bpt = bnt
      ELSE
        ant = -ant
        bnt = -bnt
        x = -x
      END IF
    END DO
  ELSE
    q2 = 1.414213562373095_dp
    q0 = .3333333333333333_dp
    q1 = .6666666666666667_dp
    xe = x * SQRT(x) / 1.5_dp
    xp6 = 1.0_dp / SQRT(6.0_dp*pi*xe)
    su1 = 1.0_dp
    r = 1.0_dp
    xr1 = 1.0_dp / xe
    DO  k = 1, 16
      r = -r * xr1
      su1 = su1 + a(k) * r
    END DO
    su2 = 1.0_dp
    r = 1.0_dp
    DO  k = 1, 16
      r = r * xr1
      su2 = su2 + a(k) * r
    END DO
    apt = q0 - EXP(-xe) * xp6 * su1
    bpt = 2.0_dp * EXP(xe) * xp6 * su2
    su3 = 1.0_dp
    r = 1.0_dp
    xr2 = 1.0_dp / (xe*xe)
    DO  k = 1, 8
      r = -r * xr2
      su3 = su3 + a(2*k) * r
    END DO
    su4 = a(1) * xr1
    r = xr1
    DO  k = 1, 7
      r = -r * xr2
      su4 = su4 + a(2*k+1) * r
    END DO
    su5 = su3 + su4
    su6 = su3 - su4
    ant = q1 - q2 * xp6 * (su5*COS(xe)-su6*SIN(xe))
    bnt = q2 * xp6 * (su5*SIN(xe)+su6*COS(xe))
  END IF
END IF
RETURN
END SUBROUTINE itairy
 
END MODULE itairy_func
 
 
 
PROGRAM mitairy
USE itairy_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!    ===========================================================
!    Purpose: This program computes the integrals of Airy
!             functions using subroutine ITAIRY
!    Input  : x   --- Upper limit of the integral
!    Output : APT --- Integration of Ai(t) from 0 and x
!             BPT --- Integration of Bi(t) from 0 and x
!             ANT --- Integration of Ai(-t) from 0 and x
!             BNT --- Integration of Bi(-t) from 0 and x
!    Example:

!      x      Ai(t)dt       Bi(t)dt       Ai(-t)dt     Bi(-t)dt
!     ----------------------------------------------------------
!      5    .33328759   .32147832D+03    .71788220    .15873094
!     10    .33333333   .14780980D+09    .76569840    .01504043
!     15    .33333333   .49673090D+16    .68358063    .07202621
!     20    .33333333   .47447423D+25    .71173925   -.03906173
!     25    .33333333   .78920820D+35    .70489539    .03293190
!    ===========================================================

REAL (dp)  :: x, apt, bpt, ant, bnt

DO
  WRITE (*,*) 'Please enter x: '
  READ (*,*) x
  WRITE (*,5100)
  WRITE (*,5200)
  CALL itairy(x, apt, bpt, ant, bnt)
  WRITE (*,5000) x, apt, bpt, ant, bnt
END DO
STOP

5000 FORMAT (' ', f5.1, f14.8, '  ', g15.8, 2F14.8)
5100 FORMAT (t4, 'x        Ai(t)dt       Bi(t)dt         Ai(-t)dt      Bi(-t)dt')
5200 FORMAT (t3, '----------------------------------------------------------------')
END PROGRAM mitairy

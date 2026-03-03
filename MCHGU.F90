MODULE chgu_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 27 December 2001
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variables h0 & r0 were used without values assigned to them

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE chgu(a, b, x, hu, md)

!       =======================================================
!       Purpose: Compute the confluent hypergeometric function
!                U(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument  ( x > 0 )
!       Output:  HU --- U(a,b,x)
!                MD --- Method code
!       Routines called:
!            (1) CHGUS for small x ( MD=1 )
!            (2) CHGUL for large x ( MD=2 )
!            (3) CHGUBI for integer b ( MD=3 )
!            (4) CHGUIT for numerical integration ( MD=4 )
!       =======================================================

REAL (dp), INTENT(IN OUT)  :: a
REAL (dp), INTENT(IN OUT)  :: b
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(IN OUT)  :: hu
INTEGER, INTENT(OUT)       :: md

LOGICAL    :: il1, il2, il3, bl1, bl2, bl3, bn
REAL (dp)  :: a00, aa, b00, hu1
INTEGER    :: id, id1

aa = a - b + 1.0_dp
il1 = a == INT(a) .AND. a <= 0.0
il2 = aa == INT(aa) .AND. aa <= 0.0
il3 = ABS(a*(a-b+1.0)) / x <= 2.0
bl1 = x <= 5.0 .OR. (x <= 10.0 .AND. a <= 2.0)
bl2 = (x > 5.0 .AND. x <= 12.5) .AND. (a >= 1.0 .AND. b >= a+4.0)
bl3 = x > 12.5 .AND. a >= 5.0 .AND. b >= a + 5.0
bn = b == INT(b) .AND. b /= 0.0
id1 = -100
IF (b /= INT(b)) THEN
  CALL chgus(a, b, x, hu, id1)
  md = 1
  IF (id1 >= 6) RETURN
  hu1 = hu
END IF
IF (il1 .OR. il2 .OR. il3) THEN
  CALL chgul(a, b, x, hu, id)
  md = 2
  IF (id >= 6) RETURN
  IF (id1 > id) THEN
    md = 1
    id = id1
    hu = hu1
  END IF
END IF
IF (a >= 0.0) THEN
  IF (bn .AND. (bl1 .OR. bl2 .OR. bl3)) THEN
    CALL chgubi(a, b, x, hu, id)
    md = 3
  ELSE
    CALL chguit(a, b, x, hu, id)
    md = 4
  END IF
ELSE
  IF (b <= a) THEN
    a00 = a
    b00 = b
    a = a - b + 1.0_dp
    b = 2.0_dp - b
    CALL chguit(a, b, x, hu, id)
    hu = x ** (1.0_dp-b00) * hu
    a = a00
    b = b00
    md = 4
  ELSE IF (bn .AND. (.NOT.il1)) THEN
    CALL chgubi(a, b, x, hu, id)
    md = 3
  END IF
END IF
IF (id < 6) WRITE (*,*) 'No accurate result obtained'
RETURN
END SUBROUTINE chgu



SUBROUTINE chgus(a, b, x, hu, id)

!       ======================================================
!       Purpose: Compute confluent hypergeometric function
!                U(a,b,x) for small argument x
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,...)
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       Routine called: GAMMA for computing gamma function
!       ======================================================

REAL (dp), INTENT(IN)      :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: hu
INTEGER, INTENT(OUT)       :: id

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: d1, d2, ga, gab, gb, gb2, h0, hmax, hmin, hu0, hua,  &
              r1, r2, xg1, xg2
INTEGER    :: j

id = -100
CALL gamma(a, ga)
CALL gamma(b, gb)
xg1 = 1.0_dp + a - b
CALL gamma(xg1, gab)
xg2 = 2.0_dp - b
CALL gamma(xg2, gb2)
hu0 = pi / SIN(pi*b)
r1 = hu0 / (gab*gb)
r2 = hu0 * x ** (1.0_dp-b) / (ga*gb2)
hu = r1 - r2
h0 = hu
hmax = 0.0_dp
hmin = 1.0D+300
DO  j = 1, 150
  r1 = r1 * (a+j-1.0_dp) / (j*(b+j-1.0_dp)) * x
  r2 = r2 * (a-b+j) / (j*(1.0_dp-b+j)) * x
  hu = hu + r1 - r2
  hua = ABS(hu)
  IF (hua > hmax) hmax = hua
  IF (hua < hmin) hmin = hua
  IF (ABS(hu-h0) < ABS(hu)*1.0D-15) EXIT
  h0 = hu
END DO

d1 = LOG10(hmax)
IF (hmin /= 0.0) d2 = LOG10(hmin)
id = 15 - ABS(d1-d2)
RETURN
END SUBROUTINE chgus



SUBROUTINE chgul(a, b, x, hu, id)

!       =======================================================
!       Purpose: Compute the confluent hypergeometric function
!                U(a,b,x) for large argument x
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       =======================================================

REAL (dp), INTENT(IN)      :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: hu
INTEGER, INTENT(OUT)       :: id

LOGICAL    :: il1, il2
REAL (dp)  :: aa, r, r0, ra
INTEGER    :: k, nm

id = -100
aa = a - b + 1.0_dp
il1 = a == INT(a) .AND. a <= 0.0
il2 = aa == INT(aa) .AND. aa <= 0.0
IF (il1) nm = ABS(a)
IF (il2) nm = ABS(aa)
IF (il1 .OR. il2) THEN
  hu = 1.0_dp
  r = 1.0_dp
  DO  k = 1, nm
    r = -r * (a+k-1.0_dp) * (a-b+k) / (k*x)
    hu = hu + r
  END DO
  hu = x ** (-a) * hu
  id = 10
ELSE
  hu = 1.0_dp
  r = 1.0_dp
  r0 = r
  DO  k = 1, 25
    r = -r * (a+k-1.0_dp) * (a-b+k) / (k*x)
    ra = ABS(r)
    IF (k > 5 .AND. ra >= r0 .OR. ra < 1.0D-15) EXIT
    r0 = ra
    hu = hu + r
  END DO
  id = ABS(LOG10(ra))
  hu = x ** (-a) * hu
END IF
RETURN
END SUBROUTINE chgul



SUBROUTINE chgubi(a, b, x, hu, id)

!       ======================================================
!       Purpose: Compute confluent hypergeometric function
!                U(a,b,x) with integer b ( b = ס1,ס2,... )
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument
!       Output:  HU --- U(a,b,x)
!                ID --- Estimated number of significant digits
!       Routines called:
!            (1) GAMMA for computing gamma function ג(x)
!            (2) PSI for computing psi function
!       ======================================================

REAL (dp), INTENT(IN)      :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: hu
INTEGER, INTENT(OUT)       :: id

REAL (dp), PARAMETER  :: el = 0.5772156649015329_dp
REAL (dp)  :: a0, a1, a2, da1, da2, db1, db2, ga, ga1, h0, hm1, hm2, hm3,  &
              hmax, hmin, hu1, hu2, hw, ps, r, rn, rn1, s0, s1, s2,  &
              sa, sb, ua, ub
INTEGER    :: id1, id2, j, k, m, n

id = -100
n = ABS(b-1)
rn1 = 1.0_dp
rn = 1.0_dp
DO  j = 1, n
  rn = rn * j
  IF (j == n-1) rn1 = rn
END DO
CALL psi(a, ps)
CALL gamma(a, ga)
IF (b > 0.0) THEN
  a0 = a
  a1 = a - n
  a2 = a1
  CALL gamma(a1, ga1)
  ua = (-1) ** (n-1) / (rn*ga1)
  ub = rn1 / ga * x ** (-n)
ELSE
  a0 = a + n
  a1 = a0
  a2 = a
  CALL gamma(a1, ga1)
  ua = (-1) ** (n-1) / (rn*ga) * x ** n
  ub = rn1 / ga1
END IF
hm1 = 1.0_dp
r = 1.0_dp
hmax = 0.0_dp
hmin = 1.0D+300
DO  k = 1, 150
  r = r * (a0+k-1.0_dp) * x / ((n+k)*k)
  hm1 = hm1 + r
  hu1 = ABS(hm1)
  IF (hu1 > hmax) hmax = hu1
  IF (hu1 < hmin) hmin = hu1
  IF (ABS(hm1-h0) < ABS(hm1)*1.0D-15) EXIT
  h0 = hm1
END DO

da1 = LOG10(hmax)
IF (hmin /= 0.0) da2 = LOG10(hmin)
id = 15 - ABS(da1-da2)
hm1 = hm1 * LOG(x)
s0 = 0.0_dp
DO  m = 1, n
  IF (b >= 0.0) s0 = s0 - 1.0_dp / m
  IF (b < 0.0) s0 = s0 + (1.0_dp-a) / (m*(a+m-1.0_dp))
END DO
hm2 = ps + 2.0_dp * el + s0
r = 1.0_dp
hmax = 0.0_dp
hmin = 1.0D+300
DO  k = 1, 150
  s1 = 0.0_dp
  s2 = 0.0_dp
  IF (b > 0.0) THEN
    DO  m = 1, k
      s1 = s1 - (m+2.0_dp*a-2.0_dp) / (m*(m+a-1.0_dp))
    END DO
    DO  m = 1, n
      s2 = s2 + 1.0_dp / (k+m)
    END DO
  ELSE
    DO  m = 1, k + n
      s1 = s1 + (1.0_dp-a) / (m*(m+a-1.0_dp))
    END DO
    DO  m = 1, k
      s2 = s2 + 1.0_dp / m
    END DO
  END IF
  hw = 2.0_dp * el + ps + s1 - s2
  r = r * (a0+k-1.0_dp) * x / ((n+k)*k)
  hm2 = hm2 + r * hw
  hu2 = ABS(hm2)
  IF (hu2 > hmax) hmax = hu2
  IF (hu2 < hmin) hmin = hu2
  IF (ABS((hm2-h0)/hm2) < 1.0D-15) EXIT
  h0 = hm2
END DO

db1 = LOG10(hmax)
IF (hmin /= 0.0) db2 = LOG10(hmin)
id1 = 15 - ABS(db1-db2)
IF (id1 < id) id = id1
hm3 = 1.0_dp
IF (n == 0) hm3 = 0.0_dp
r = 1.0_dp
DO  k = 1, n - 1
  r = r * (a2+k-1.0_dp) / ((k-n)*k) * x
  hm3 = hm3 + r
END DO
sa = ua * (hm1+hm2)
sb = ub * hm3
hu = sa + sb
IF (sa /= 0.0) id1 = INT(LOG10(ABS(sa)))
IF (hu /= 0.0) id2 = INT(LOG10(ABS(hu)))
IF (sa*sb < 0.0) id = id - ABS(id1-id2)
RETURN
END SUBROUTINE chgubi



SUBROUTINE chguit(a, b, x, hu, id)

!       ======================================================
!       Purpose: Compute hypergeometric function U(a,b,x) by
!                using Gaussian-Legendre integration (n=60)
!       Input  : a  --- Parameter ( a > 0 )
!                b  --- Parameter
!                x  --- Argument ( x > 0 )
!       Output:  HU --- U(a,b,z)
!                ID --- Estimated number of significant digits
!       Routine called: GAMMA for computing ג(x)
!       ======================================================

REAL (dp), INTENT(IN)      :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: hu
INTEGER, INTENT(OUT)       :: id

REAL (dp), PARAMETER  :: t(30) = (/ 0.259597723012478D-01,  &
    .778093339495366D-01, .129449135396945_dp, .180739964873425_dp,  &
    .231543551376029_dp, .281722937423262_dp, .331142848268448_dp,  &
    .379670056576798_dp, .427173741583078_dp, .473525841761707_dp,  &
    .518601400058570_dp, .562278900753945_dp, .604440597048510_dp,  &
    .644972828489477_dp, .683766327381356_dp, .720716513355730_dp,  &
    .755723775306586_dp, .788693739932264_dp, .819537526162146_dp,  &
    .848171984785930_dp, .874519922646898_dp, .898510310810046_dp,  &
    .920078476177628_dp, .939166276116423_dp, .955722255839996_dp,  &
    .969701788765053_dp, .981067201752598_dp, .989787895222222_dp,  &
    .995840525118838_dp, .999210123227436_dp /)
REAL (dp), PARAMETER  :: w(30) = (/ 0.519078776312206D-01,  &
    .517679431749102D-01, .514884515009810D-01, .510701560698557D-01,  &
    .505141845325094D-01, .498220356905502D-01, .489955754557568D-01,  &
    .480370318199712D-01, .469489888489122D-01, .457343797161145D-01,  &
    .443964787957872D-01, .429388928359356D-01, .413655512355848D-01,  &
    .396806954523808D-01, .378888675692434D-01, .359948980510845D-01,  &
    .340038927249464D-01, .319212190192963D-01, .297524915007890D-01,  &
    .275035567499248D-01, .251804776215213D-01, .227895169439978D-01,  &
    .203371207294572D-01, .178299010142074D-01, .152746185967848D-01,  &
    .126781664768159D-01, .100475571822880D-01, .738993116334531D-02,  &
    .471272992695363D-02, .202681196887362D-02 /)

REAL (dp)  :: a1, b1, c, d, f1, f2, g, ga, hu0, hu1, hu2,   &
              s, t1, t2, t3, t4
INTEGER    :: j, k, m

id = 7
a1 = a - 1.0_dp
b1 = b - a - 1.0_dp
c = 12.0 / x
DO  m = 10, 100, 5
  hu1 = 0.0_dp
  g = 0.5_dp * c / m
  d = g
  DO  j = 1, m
    s = 0.0_dp
    DO  k = 1, 30
      t1 = d + g * t(k)
      t2 = d - g * t(k)
      f1 = EXP(-x*t1) * t1 ** a1 * (1.0_dp+t1) ** b1
      f2 = EXP(-x*t2) * t2 ** a1 * (1.0_dp+t2) ** b1
      s = s + w(k) * (f1+f2)
    END DO
    hu1 = hu1 + s * g
    d = d + 2.0_dp * g
  END DO
  IF (ABS(1.0_dp-hu0/hu1) < 1.0D-7) EXIT
  hu0 = hu1
END DO

CALL gamma(a, ga)
hu1 = hu1 / ga
DO  m = 2, 10, 2
  hu2 = 0.0_dp
  g = 0.5_dp / m
  d = g
  DO  j = 1, m
    s = 0.0_dp
    DO  k = 1, 30
      t1 = d + g * t(k)
      t2 = d - g * t(k)
      t3 = c / (1.0_dp-t1)
      t4 = c / (1.0_dp-t2)
      f1 = t3 * t3 / c * EXP(-x*t3) * t3 ** a1 * (1.0_dp+t3) ** b1
      f2 = t4 * t4 / c * EXP(-x*t4) * t4 ** a1 * (1.0_dp+t4) ** b1
      s = s + w(k) * (f1+f2)
    END DO
    hu2 = hu2 + s * g
    d = d + 2.0_dp * g
  END DO
  IF (ABS(1.0_dp - hu0/hu2) < 1.0D-7) EXIT
  hu0 = hu2
END DO

CALL gamma(a, ga)
hu2 = hu2 / ga
hu = hu1 + hu2
RETURN
END SUBROUTINE chguit



SUBROUTINE gamma(x, ga)

!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GA --- ג(x)
!       ==================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ga

REAL (dp), PARAMETER  :: g(26) = (/1.0_dp, 0.5772156649015329_dp, -0.6558780715202538_dp,  &
      -0.420026350340952D-1, 0.1665386113822915_dp, -  &
      .421977345555443D-1, -.96219715278770D-2,  &
      .72189432466630D-2, -.11651675918591D-2, -.2152416741149D-3,  &
      .1280502823882D-3, -.201348547807D-4, -.12504934821D-5,  &
      .11330272320D-5, -.2056338417D-6, .61160950D-8,  &
      .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: gr, r, z
INTEGER    :: k, m, m1

IF (x == INT(x)) THEN
  IF (x > 0.0_dp) THEN
    ga = 1.0_dp
    m1 = x - 1
    DO  k = 2, m1
      ga = ga * k
    END DO
  ELSE
    ga = 1.0D+300
  END IF
ELSE
  IF (ABS(x) > 1.0_dp) THEN
    z = ABS(x)
    m = INT(z)
    r = 1.0_dp
    DO  k = 1, m
      r = r * (z-k)
    END DO
    z = z - m
  ELSE
    z = x
  END IF

  gr = g(26)
  DO  k = 25, 1, -1
    gr = gr * z + g(k)
  END DO
  ga = 1.0_dp / (gr*z)
  IF (ABS(x) > 1.0_dp) THEN
    ga = ga * r
    IF (x < 0.0_dp) ga = -pi / (x*ga*SIN(pi*x))
  END IF
END IF
RETURN
END SUBROUTINE gamma



SUBROUTINE psi(x, ps)

!       ======================================
!       Purpose: Compute Psi function
!       Input :  x  --- Argument of psi(x)
!       Output:  PS --- psi(x)
!       ======================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ps

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp)  :: a1, a2, a3, a4, a5, a6, a7, a8, s, x2, xa
INTEGER    :: k, n

xa = ABS(x)
s = 0.0_dp
IF (x == INT(x) .AND. x <= 0.0) THEN
  ps = 1.0D+300
  RETURN
ELSE IF (xa == INT(xa)) THEN
  n = xa
  DO  k = 1, n - 1
    s = s + 1.0_dp / k
  END DO
  ps = -el + s
ELSE IF (xa+.5 == INT(xa+.5)) THEN
  n = xa - .5
  DO  k = 1, n
    s = s + 1.0 / (2.0_dp*k-1.0_dp)
  END DO
  ps = -el + 2.0_dp * s - 1.386294361119891_dp
ELSE
  IF (xa < 10.0) THEN
    n = 10 - INT(xa)
    DO  k = 0, n - 1
      s = s + 1.0_dp / (xa+k)
    END DO
    xa = xa + n
  END IF
  x2 = 1.0_dp / (xa*xa)
  a1 = -.8333333333333D-01
  a2 = .83333333333333333D-02
  a3 = -.39682539682539683D-02
  a4 = .41666666666666667D-02
  a5 = -.75757575757575758D-02
  a6 = .21092796092796093D-01
  a7 = -.83333333333333333D-01
  a8 = .4432598039215686_dp
  ps = LOG(xa) - .5_dp / xa + x2 * (((((((a8*x2 + a7)*x2 + a6)*x2 + a5)*  &
       x2 + a4)*x2 + a3)*x2 + a2)*x2 + a1)
  ps = ps - s
END IF
IF (x < 0.0) ps = ps - pi * COS(pi*x) / SIN(pi*x) - 1.0_dp / x
RETURN
END SUBROUTINE psi
 
END MODULE chgu_func
 
 
 
PROGRAM mchgu
USE chgu_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!       =======================================================
!       Purpose: This program computes the confluent
!                hypergeometric function U(a,b,x) using
!                subroutine CHGU
!       Input  : a  --- Parameter
!                b  --- Parameter
!                x  --- Argument  ( x ע 0 )
!       Output:  HU --- U(a,b,x)
!                MD --- Method code
!       Example:
!                a       b       x        U(a,b,x)
!             --------------------------------------
!              -2.5     2.5     5.0     -9.02812446
!              -1.5     2.5     5.0      2.15780560
!               -.5     2.5     5.0      1.76649370
!                .0     2.5     5.0      1.00000000
!                .5     2.5     5.0       .49193496
!               1.5     2.5     5.0       .08944272
!               2.5     2.5     5.0       .01239387

!                a       b       x        U(a,b,x)
!             --------------------------------------
!              -2.5     5.0    10.0     -2.31982196
!              -1.5     5.0    10.0      8.65747115
!               -.5     5.0    10.0      2.37997143
!                .0     5.0    10.0      1.00000000
!                .5     5.0    10.0       .38329536
!               1.5     5.0    10.0       .04582817
!               2.5     5.0    10.0       .00444535
!       =======================================================

REAL (dp)  :: a, b, hu, x
INTEGER    :: md

DO
  WRITE (*,*) 'Please enter a, b and x: '
  READ (*,*) a, b, x
  WRITE (*,*) '   a       b       x        U(a,b,x)'
  WRITE (*,*) '--------------------------------------'
  CALL chgu(a, b, x, hu, md)
  WRITE (*,5000) a, b, x, hu
END DO
STOP

5000 FORMAT (' ', f5.1, '   ', f5.1, '   ', f5.1, g15.8)
END PROGRAM mchgu

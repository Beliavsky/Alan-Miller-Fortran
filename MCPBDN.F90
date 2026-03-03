MODULE cpbdn_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
 
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

CONTAINS
 

SUBROUTINE cpbdn(n, z, cpb, cpd)

!     ==================================================
!     Purpose: Compute the parabolic cylinder functions
!               Dn(z) and Dn'(z) for a complex argument
!     Input:   z --- Complex argument of Dn(z)
!              n --- Order of Dn(z)  ( n=0,ס1,ס2,תתת )
!     Output:  CPB(|n|) --- Dn(z)
!              CPD(|n|) --- Dn'(z)
!     Routines called:
!          (1) CPDSA for computing Dn(z) for a small |z|
!          (2) CPDLA for computing Dn(z) for a large |z|
!     ==================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cpb(0:)
COMPLEX (dp), INTENT(OUT)  :: cpd(0:)

COMPLEX (dp)  :: c0, ca0, cf, cf0, cf1, cfa, cfb, cs0, z1
REAL (dp)     :: a0, x
INTEGER       :: k, m, n0, n1, nm1

x = REAL(z, KIND=dp)
a0 = ABS(z)
c0 = (0.0_dp,0.0_dp)
ca0 = EXP(-0.25_dp*z*z)
IF (n >= 0) THEN
  cf0 = ca0
  cf1 = z * ca0
  cpb(0) = cf0
  cpb(1) = cf1
  DO  k = 2, n
    cf = z * cf1 - (k-1) * cf0
    cpb(k) = cf
    cf0 = cf1
    cf1 = cf
  END DO
ELSE
  n0 = -n
  IF (x <= 0.0 .OR. ABS(z) == 0.0) THEN
    cf0 = ca0
    cpb(0) = cf0
    z1 = -z
    IF (a0 <= 7.0) THEN
      CALL cpdsa(-1, z1, cf1)
    ELSE
      CALL cpdla(-1, z1, cf1)
    END IF
    cf1 = SQRT(2.0_dp*pi) / ca0 - cf1
    cpb(1) = cf1
    DO  k = 2, n0
      cf = (-z*cf1 + cf0) / (k-1)
      cpb(k) = cf
      cf0 = cf1
      cf1 = cf
    END DO
  ELSE
    IF (a0 <= 3.0) THEN
      CALL cpdsa(-n0, z, cfa)
      cpb(n0) = cfa
      n1 = n0 + 1
      CALL cpdsa(-n1, z, cfb)
      cpb(n1) = cfb
      nm1 = n0 - 1
      DO  k = nm1, 0, -1
        cf = z * cfa + (k+1) * cfb
        cpb(k) = cf
        cfb = cfa
        cfa = cf
      END DO
    ELSE
      m = 100 + ABS(n)
      cfa = c0
      cfb = (1.0D-30,0.0_dp)
      DO  k = m, 0, -1
        cf = z * cfb + (k+1) * cfa
        IF (k <= n0) cpb(k) = cf
        cfa = cfb
        cfb = cf
      END DO
      cs0 = ca0 / cf
      DO  k = 0, n0
        cpb(k) = cs0 * cpb(k)
      END DO
    END IF
  END IF
END IF
cpd(0) = -0.5_dp * z * cpb(0)
IF (n >= 0) THEN
  DO  k = 1, n
    cpd(k) = -0.5_dp * z * cpb(k) + k * cpb(k-1)
  END DO
ELSE
  DO  k = 1, n0
    cpd(k) = 0.5_dp * z * cpb(k) - cpb(k-1)
  END DO
END IF
RETURN
END SUBROUTINE cpbdn



SUBROUTINE cpdsa(n, z, cdn)

!     ===========================================================
!     Purpose: Compute complex parabolic cylinder function Dn(z)
!              for small argument
!     Input:   z   --- complex argument of D(z)
!              n   --- Order of D(z) (n = 0,-1,-2,תתת)
!     Output:  CDN --- Dn(z)
!     Routine called: GAIH for computing ג(x), x=n/2 (n=1,2,...)
!     ===========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cdn

COMPLEX (dp)  :: ca0, cb0, cdw, cr
REAL (dp)     :: eps, g0, g1, ga0, gm, pd, sq2, va0, vm, vt, xn
INTEGER       :: m

eps = 1.0D-15
sq2 = SQRT(2.0_dp)
ca0 = EXP(-.25_dp*z*z)
va0 = 0.5_dp * (1.0_dp-n)
IF (n == 0.0) THEN
  cdn = ca0
ELSE
  IF (ABS(z) == 0.0) THEN
    IF (va0 <= 0.0 .AND. va0 == INT(va0)) THEN
      cdn = 0.0_dp
    ELSE
      CALL gaih(va0, ga0)
      pd = SQRT(pi) / (2.0_dp**(-.5_dp*n)*ga0)
      cdn = CMPLX(pd, 0.0_dp, KIND=dp)
    END IF
  ELSE
    xn = -n
    CALL gaih(xn, g1)
    cb0 = 2.0_dp ** (-0.5_dp*n-1.0_dp) * ca0 / g1
    vt = -.5_dp * n
    CALL gaih(vt, g0)
    cdn = CMPLX(g0, 0.0_dp, KIND=dp)
    cr = (1.0_dp,0.0_dp)
    DO  m = 1, 250
      vm = .5_dp * (m-n)
      CALL gaih(vm, gm)
      cr = -cr * sq2 * z / m
      cdw = gm * cr
      cdn = cdn + cdw
      IF (ABS(cdw) < ABS(cdn)*eps) EXIT
    END DO
    cdn = cb0 * cdn
  END IF
END IF
RETURN
END SUBROUTINE cpdsa



SUBROUTINE cpdla(n, z, cdn)

!     ===========================================================
!     Purpose: Compute complex parabolic cylinder function Dn(z)
!              for large argument
!     Input:   z   --- Complex argument of Dn(z)
!              n   --- Order of Dn(z) (n = 0,ס1,ס2,תתת)
!     Output:  CDN --- Dn(z)
!     ===========================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cdn

COMPLEX (dp)  :: cb0, cr
INTEGER       :: k

cb0 = z ** n * EXP(-.25_dp*z*z)
cr = (1.0_dp,0.0_dp)
cdn = (1.0_dp,0.0_dp)
DO  k = 1, 16
  cr = -0.5_dp * cr * (2*k-n-1) * (2*k-n-2) / (k*z*z)
  cdn = cdn + cr
  IF (ABS(cr) < ABS(cdn)*1.0D-12) EXIT
END DO
cdn = cb0 * cdn
RETURN
END SUBROUTINE cpdla



SUBROUTINE gaih(x, ga)

!     =====================================================
!     Purpose: Compute gamma function ג(x)
!     Input :  x  --- Argument of ג(x), x = n/2, n=1,2,תתת
!     Output:  GA --- ג(x)
!     =====================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ga

INTEGER    :: k, m, m1

IF (x == INT(x) .AND. x > 0.0) THEN
  ga = 1.0_dp
  m1 = INT(x-1.0)
  DO  k = 2, m1
    ga = ga * k
  END DO
ELSE IF (x + .5_dp == INT(x + .5_dp) .AND. x > 0.0) THEN
  m = INT(x)
  ga = SQRT(pi)
  DO  k = 1, m
    ga = 0.5_dp * ga * (2*k-1)
  END DO
END IF
RETURN
END SUBROUTINE gaih
 
END MODULE cpbdn_func
 
 
 
PROGRAM mcpbdn
USE cpbdn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!     =============================================================
!     Purpose: This program computes parabolic cylinder functions
!              Dn(z) for an integer order and a complex argument
!              using subroutine CPBDN
!     Input :  x --- Real part of z
!              y --- Imaginary part of z
!              n --- Order of Dn(z)
!     Output:  CPB(|n|) --- Dn(z)
!              CPD(|n|) --- Dn'(z)
!     Example:
!              z = 5.0 + 5.0 i

!     n     Re[Dn(z)]      Im[Dn(z)]      Re[Dn'(z)]     Im[Dn'(z)]
!   -----------------------------------------------------------------
!     0   .99779828D+00  .66321897D-01 -.23286910D+01 -.26603004D+01
!     1   .46573819D+01  .53206009D+01  .26558457D+01 -.24878635D+02
!     2  -.43138931D+01  .49823592D+02  .14465848D+03 -.10313305D+03
!     3  -.28000219D+03  .21690729D+03  .12293320D+04  .30720802D+03
!     4  -.24716057D+04 -.46494526D+03  .38966424D+04  .82090067D+04
!    -1   .10813809D+00 -.90921592D-01 -.50014908D+00 -.23280660D-01
!    -2   .24998820D-02 -.19760577D-01 -.52486940D-01  .47769856D-01
!    -3  -.15821033D-02 -.23090595D-02 -.68249161D-03  .10032670D-01
!    -4  -.37829961D-03 -.10158757D-03  .89032322D-03  .11093416D-02
!       =============================================================

COMPLEX (dp)  :: cpb(0:100), cpd(0:100), z
REAL (dp)     :: x, y
INTEGER       :: i, n, n0

WRITE(*,*) 'Please enter n, x and y '
READ(*,*) n, x, y
WRITE(*,5100) n, x, y
z = CMPLX(x, y, KIND=dp)
n0 = ABS(n)
CALL cpbdn(n, z, cpb, cpd)
WRITE (*,*)
IF (n >= 0) THEN
  WRITE(*,*)'  n     Re[Dn(z)]       Im[Dn(z)]       Re[Dn''(Z)]      IM[DN''(Z)]'
ELSE
  WRITE(*,*)' -n     Re[Dn(z)]       Im[Dn(z)]       Re[Dn''(Z)]      IM[DN''(Z)]'
END IF
WRITE(*,*) '--------------------------------------------------------------------'
DO  i = 0, n0
  WRITE(*,5000) i, cpb(i), cpd(i)
END DO
STOP

5000 FORMAT(' ', i3, 4g16.8)
5100 FORMAT(' N =', i3, ',   z =x+iy :', f6.2, '+', f6.2, ' i')
END PROGRAM mcpbdn

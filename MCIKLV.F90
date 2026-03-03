MODULE ciklv_func
 
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
 


SUBROUTINE ciklv(v, z, cbiv, cdiv, cbkv, cdkv)

!       =====================================================
!       Purpose: Compute modified Bessel functions Iv(z) and
!                Kv(z) and their derivatives with a complex
!                argument and a large order
!       Input:   v --- Order of Iv(z) and Kv(z)
!                z --- Complex argument
!       Output:  CBIV --- Iv(z)
!                CDIV --- Iv'(z)
!                CBKV --- Kv(z)
!                CDKV --- Kv'(z)
!       Routine called:
!                CJK to compute the expansion coefficients
!       ====================================================

REAL (dp), INTENT(IN)      :: v
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cbiv
COMPLEX (dp), INTENT(OUT)  :: cdiv
COMPLEX (dp), INTENT(OUT)  :: cbkv
COMPLEX (dp), INTENT(OUT)  :: cdkv

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)     :: a(91), v0, vr
COMPLEX (dp)  :: ceta, cf(12), cfi, cfk, csi, csk, ct, ct2, cws
INTEGER       :: i, k, km, l, l0, lf

km = 12
CALL cjk(km, a)
DO  l = 1, 0, -1
  v0 = v - l
  cws = SQRT(1.0D0 + (z/v0)*(z/v0))
  ceta = cws + LOG(z/v0/(1.0D0 + cws))
  ct = 1.0D0 / cws
  ct2 = ct * ct
  DO  k = 1, km
    l0 = k * (k+1) / 2 + 1
    lf = l0 + k
    cf(k) = a(lf)
    DO  i = lf - 1, l0, -1
      cf(k) = cf(k) * ct2 + a(i)
    END DO
    cf(k) = cf(k) * ct ** k
  END DO
  vr = 1.0D0 / v0
  csi = (1.0D0,0.0D0)
  DO  k = 1, km
    csi = csi + cf(k) * vr ** k
  END DO
  cbiv = SQRT(ct/(2.0D0*pi*v0)) * EXP(v0*ceta) * csi
  IF (l == 1) cfi = cbiv
  csk = (1.0D0,0.0D0)
  DO  k = 1, km
    csk = csk + (-1) ** k * cf(k) * vr ** k
  END DO
  cbkv = SQRT(pi*ct/(2.0D0*v0)) * EXP(-v0*ceta) * csk
  IF (l == 1) cfk = cbkv
END DO
cdiv = cfi - v / z * cbiv
cdkv = -cfk - v / z * cbkv
RETURN
END SUBROUTINE ciklv



SUBROUTINE cjk(km, a)

!       ========================================================
!       Purpose: Compute the expansion coefficients for the
!                asymptotic expansion of Bessel functions
!                with large orders
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L
!                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
!       ========================================================

INTEGER, INTENT(IN)        :: km
REAL (dp), INTENT(OUT)     :: a(:)

REAL (dp)  :: f, f0, g, g0
INTEGER    :: j, k, l1, l2, l3, l4

a(1) = 1.0D0
f0 = 1.0D0
g0 = 1.0D0
DO  k = 0, km - 1
  l1 = (k+1) * (k+2) / 2 + 1
  l2 = (k+1) * (k+2) / 2 + k + 2
  f = (0.5D0*k + 0.125D0/(k+1)) * f0
  g = -(1.5D0*k + 0.625D0/(3*(k+1))) * g0
  a(l1) = f
  a(l2) = g
  f0 = f
  g0 = g
END DO
DO  k = 1, km - 1
  DO  j = 1, k
    l3 = k * (k+1) / 2 + j + 1
    l4 = (k+1) * (k+2) / 2 + j + 1
    a(l4) = (j + 0.5D0*k + 0.125D0/(2*j+k+1)) * a(l3) - (j + 0.5D0*k  &
            -1.0 + 0.625D0/(2*j+k+1)) * a(l3-1)
  END DO
END DO
RETURN
END SUBROUTINE cjk

END MODULE ciklv_func
 
 
 
PROGRAM mciklv
USE ciklv_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:35

!       =========================================================
!       Purpose: This program computes modified Bessel functions
!                Iv(z) and Kv(z) and their derivatives for a
!                large order and a complex argument using
!                subroutine CIKLV
!       Input:   v --- Order of Iv(z) and Kv(z)
!                z --- Complex argument
!       Output:  CBIV --- Iv(z)
!                CDIV --- Iv'(z)
!                CBKV --- Kv(z)
!                CDKV --- Kv'(z)
!       Examples:
!                v =100.00,    z =   4.00 + i   2.00

!       Iv(z) = -.7373606617-123 + .6461109082-123 i
!       Iv'(z)= -.8307094243-122 + .2030132500-121 i
!       Kv(z) = -.3836166007+121 - .3356017795+121 i
!       Kv'(z)=  .1103271276+123 + .2886519240+122 i

!                v =100.50,    z =   4.00 + i   2.00
!       Iv(z) = -.1289940051-123 + .6845756182-124 i
!       Iv'(z)= -.1907996261-122 + .2672465997-122 i
!       Kv(z) = -.3008779281+122 - .1593719779+122 i
!       Kv'(z)=  .7653781978+123 + .1857772148+122 i
!       =========================================================

COMPLEX (dp)  :: z, cbiv, cdiv, cbkv, cdkv
REAL (dp)     :: v, x, y

WRITE (*,*) 'Please enter v,x,y ( z = x+iy )'
READ (*,*) v, x, y
WRITE (*,5000) v, x, y
z = CMPLX(x, y, KIND=dp)
CALL ciklv(v, z, cbiv, cdiv, cbkv, cdkv)
WRITE (*,*)
WRITE (*,5100) cbiv
WRITE (*,5200) cdiv
WRITE (*,*)
WRITE (*,5300) cbkv
WRITE (*,5400) cdkv
STOP

5000 FORMAT (t9, 'v =', f6.2, ',    z =', f7.2, ' + i', f7.2)
5100 FORMAT (t9, 'Iv(z) =', g17.10, ' + i ', g17.10)
5200 FORMAT (t9, 'Iv''(Z)=', g17.10, ' + I ', g17.10)
5300 FORMAT (t9, 'Kv(z) =', g17.10, ' + i ', g17.10)
5400 FORMAT (t9, 'Kv''(Z)=', g17.10, ' + I ', g17.10)
END PROGRAM mciklv

MODULE cjylv_func
 
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


SUBROUTINE cjylv(v, z, cbjv, cdjv, cbyv, cdyv)

!       ===================================================
!       Purpose: Compute Bessel functions Jv(z) and Yv(z)
!                and their derivatives with a complex
!                argument and a large order
!       Input:   v --- Order of Jv(z) and Yv(z)
!                z --- Complex argument
!       Output:  CBJV --- Jv(z)
!                CDJV --- Jv'(z)
!                CBYV --- Yv(z)
!                CDYV --- Yv'(z)
!       Routine called:
!                CJK to compute the expansion coefficients
!       ===================================================

REAL (dp), INTENT(IN)      :: v
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cbjv
COMPLEX (dp), INTENT(OUT)  :: cdjv
COMPLEX (dp), INTENT(OUT)  :: cbyv
COMPLEX (dp), INTENT(OUT)  :: cdyv

COMPLEX (dp)  :: ceta, cf(12), cfj, cfy, csj, csy, ct, ct2, cws
REAL (dp)     :: a(91), pi, v0, vr
INTEGER       :: i, k, km, l, l0, lf

km = 12
CALL cjk(km, a)
pi = 3.141592653589793_dp
DO  l = 1, 0, -1
  v0 = v - l
  cws = SQRT(1.0_dp - (z/v0)*(z/v0))
  ceta = cws + LOG(z/v0/(1.0_dp+cws))
  ct = 1.0_dp / cws
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
  vr = 1.0_dp / v0
  csj = (1.0_dp,0.0_dp)
  DO  k = 1, km
    csj = csj + cf(k) * vr ** k
  END DO
  cbjv = SQRT(ct/(2.0_dp*pi*v0)) * EXP(v0*ceta) * csj
  IF (l == 1) cfj = cbjv
  csy = (1.0_dp,0.0_dp)
  DO  k = 1, km
    csy = csy + (-1) ** k * cf(k) * vr ** k
  END DO
  cbyv = -SQRT(2.0_dp*ct/(pi*v0)) * EXP(-v0*ceta) * csy
  IF (l == 1) cfy = cbyv
END DO
cdjv = -v / z * cbjv + cfj
cdyv = -v / z * cbyv + cfy
RETURN
END SUBROUTINE cjylv



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

a(1) = 1.0_dp
f0 = 1.0_dp
g0 = 1.0_dp
DO  k = 0, km - 1
  l1 = (k+1) * (k+2) / 2 + 1
  l2 = (k+1) * (k+2) / 2 + k + 2
  f = (0.5_dp*k + 0.125_dp/(k+1)) * f0
  g = -(1.5_dp*k + 0.625_dp/(3*(k+1))) * g0
  a(l1) = f
  a(l2) = g
  f0 = f
  g0 = g
END DO
DO  k = 1, km - 1
  DO  j = 1, k
    l3 = k * (k+1) / 2 + j + 1
    l4 = (k+1) * (k+2) / 2 + j + 1
    a(l4) = (j + 0.5_dp*k + 0.125_dp/(2*j+k+1)) * a(l3) - (j+0.5_dp*k  &
            -1.0 + 0.625_dp/(2*j+k+1)) * a(l3-1)
  END DO
END DO
RETURN
END SUBROUTINE cjk

END MODULE cjylv_func
 
 
 
PROGRAM mcjylv
USE cjylv_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ========================================================
!       Purpose: This program computes Bessel functions Jv(z)
!                and Yv(z) and their derivatives with a large
!                order and complex argument using subroutine
!                CJYLV
!       Input:   v --- Order of Jv(z) and Yv(z)
!                z --- Complex argument
!       Output:  CBJV --- Jv(z)
!                CDJV --- Jv'(z)
!                CBYV --- Yv(z)
!                CDYV --- Yv'(z)
!       Examples:
!                v = 100.00,    z = 4.00 + 2.00 i

!                Jv(z) = -.6444792518-123 + .6619157435-123 i
!                Jv'(z)= -.6251103777-122 + .1967638668-121 i
!                Yv(z) =  .2403065353+121 + .2472039414+121 i
!                Yv'(z)= -.7275814786+122 - .2533588851+122 i

!                v =100.5,     z = 4.00 + 2.00 i

!                Jv(z) = -.1161315754-123 + .7390127781-124 i
!                Jv'(z)= -.1588519437-122 + .2652227059-122 i
!                Yv(z) =  .1941381412+122 + .1237578195+122 i
!                Yv'(z)= -.5143285247+123 - .5320026773+122 i
!       ========================================================

COMPLEX (dp)  :: z, cbjv, cdjv, cbyv, cdyv
REAL (dp)     :: v, x, y

WRITE (*,*) 'Please enter v,x and y ( z = x+iy )'
READ (*,*) v, x, y
WRITE (*,5000) v, x, y
z = CMPLX(x, y, KIND=dp)
CALL cjylv(v, z, cbjv, cdjv, cbyv, cdyv)
WRITE (*,*)
WRITE (*,5100) cbjv
WRITE (*,5200) cdjv
WRITE (*,*)
WRITE (*,5300) cbyv
WRITE (*,5400) cdyv
STOP

5000 FORMAT (t9, 'v = ', f6.2, ',    z =', f7.2, ' + i ', f7.2)
5100 FORMAT (t9, 'Jv(z) =', g17.10, ' + i', g17.10)
5200 FORMAT (t9, 'Jv''(Z)=', g17.10, ' + I', g17.10)
5300 FORMAT (t9, 'Yv(z) =', g17.10, ' + i', g17.10)
5400 FORMAT (t9, 'Yv''(Z)=', g17.10, ' + I', g17.10)
END PROGRAM mcjylv

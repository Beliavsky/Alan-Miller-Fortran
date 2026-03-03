MODULE comelp_func
 
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


SUBROUTINE comelp(hk, ck, ce)

!       ==================================================
!       Purpose: Compute complete elliptic integrals K(k)
!                and E(k)
!       Input  : K  --- Modulus k ( 0 ó k ó 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       ==================================================

REAL (dp), INTENT(IN)      :: hk
REAL (dp), INTENT(OUT)     :: ck
REAL (dp), INTENT(OUT)     :: ce

REAL (dp)  :: ae, ak, be, bk, pk

pk = 1.0_dp - hk * hk
IF (hk == 1.0) THEN
  ck = 1.0D+300
  ce = 1.0_dp
ELSE
  ak = (((.01451196212_dp*pk + .03742563713_dp)*pk + .03590092383_dp)*pk +   &
      .09666344259_dp) * pk + 1.38629436112_dp
  bk = (((.00441787012_dp*pk + .03328355346_dp)*pk + .06880248576_dp)*pk +   &
      .12498593597_dp) * pk + .5_dp
  ck = ak - bk * LOG(pk)
  ae = (((.01736506451_dp*pk + .04757383546_dp)*pk + .0626060122_dp)*pk +   &
      .44325141463_dp) * pk + 1.0_dp
  be = (((.00526449639_dp*pk + .04069697526_dp)*pk + .09200180037_dp)*pk +   &
      .2499836831_dp) * pk
  ce = ae - be * LOG(pk)
END IF
RETURN
END SUBROUTINE comelp
 
END MODULE comelp_func
 
 
 
PROGRAM mcomelp
USE comelp_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!       ===================================================
!       Purpose: This program computes complete elliptic
!                integrals K(k) and E(k) using subroutine
!                COMELP
!       Input  : K  --- Modulus k ( 0 ó k ó 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       Example:
!                  k         K(k)          E(K)
!                ---------------------------------
!                 .00      1.570796      1.570796
!                 .25      1.596242      1.545957
!                 .50      1.685750      1.467462
!                 .75      1.910990      1.318472
!                1.00       ì            1.000000
!       ===================================================

REAL (dp)  :: hk, ck, ce

WRITE (*,*) 'Please enter the modulus k '
READ (*,*) hk
WRITE (*,*) '    k         K(k)          E(K)'
WRITE (*,*) '  ---------------------------------'
CALL comelp(hk, ck, ce)
IF (hk /= 1.0) WRITE (*,5000) hk, ck, ce
IF (hk == 1.0) WRITE (*,5100) hk, ce
STOP

5000 FORMAT (t3, f5.2, 2F14.6)
5100 FORMAT (t3, f5.2, t15, 'ì', t22, f14.6)
END PROGRAM mcomelp

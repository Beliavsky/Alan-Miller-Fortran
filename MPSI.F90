MODULE psi_func
 
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
 

SUBROUTINE psi(x, ps)

!    ======================================
!    Purpose: Compute the psi function
!    Input :  x  --- Argument of psi(x)
!    Output:  PS --- psi(x)
!    ======================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ps

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = .5772156649015329_dp
REAL (dp), PARAMETER  :: a1 = -.8333333333333D-01, a2 = .83333333333333333D-02, &
     a3 = -.39682539682539683D-02, a4 = .41666666666666667D-02,  &
     a5 = -.75757575757575758D-02, a6 = .21092796092796093D-01,  &
     a7 = -.83333333333333333D-01, a8 = .4432598039215686_dp
REAL (dp)  :: s, x2, xa
INTEGER    :: k, n

xa = ABS(x)
s = 0.0D0
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
    s = s + 1.0_dp / (2*k-1)
  END DO
  ps = -el + 2.0_dp * s - 1.386294361119891_dp
ELSE
  IF (xa < 10.0) THEN
    n = 10 - xa
    DO  k = 0, n - 1
      s = s + 1.0_dp / (xa+k)
    END DO
    xa = xa + n
  END IF
  x2 = 1.0D0 / (xa*xa)
  ps = LOG(xa) - .5D0 / xa + x2 * (((((((a8*x2 + a7)*x2 + a6)*x2 + a5)*  &
       x2 + a4)*x2 + a3)*x2 + a2)*x2 + a1)
  ps = ps - s
END IF
IF (x < 0.0) ps = ps - pi * COS(pi*x) / SIN(pi*x) - 1.0D0 / x
RETURN
END SUBROUTINE psi
 
END MODULE psi_func
 
 
 
PROGRAM mpsi
USE psi_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:45

!    ==================================================
!    Purpose: This program computes the psi function
!             using subroutine PSI
!    Input :  x  --- Argument of psi(x)
!    Output:  PS --- psi(x)
!    Examples:
!                x          Psi(x)
!             ------------------------
!               .25      -4.227453533
!               .50      -1.963510026
!               .75      -1.085860880
!              1.00       -.577215665
!              1.25       -.227453533
!              1.50        .036489974
!              1.75        .247472454
!              2.00        .422784335
!    ==================================================

REAL (dp) :: x, ps

WRITE (*,*) 'Please enter x '
READ (*,*) x
WRITE (*,*) '    x          Psi(x)'
WRITE (*,*) ' ------------------------'
CALL psi(x, ps)
WRITE (*,5000) x, ps
STOP

5000 FORMAT (' ', f6.2, f18.9)
END PROGRAM mpsi

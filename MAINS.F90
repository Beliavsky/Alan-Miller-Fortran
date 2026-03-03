PROGRAM No_Bounds

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-14  Time: 16:13:44
 
! EASY TO USE, NO BOUNDS
!***********************************************************************
! MAIN PROGRAM TO MINIMIZE A FUNCTION (REPRESENTED BY THE ROUTINE SFUN)
! OF N VARIABLES X

USE Truncated_Newton
IMPLICIT NONE

REAL (dp)  :: x(50), f, g(50)
INTEGER    :: i, ierror, n
! EXTERNAL    sfun

INTERFACE
  SUBROUTINE sfun(n, x, f, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f, g(:)
  END SUBROUTINE sfun
END INTERFACE

! DEFINE SUBROUTINE PARAMETERS
! N  - NUMBER OF VARIABLES
! X  - INITIAL ESTIMATE OF THE SOLUTION
! F  - ROUGH ESTIMATE OF FUNCTION VALUE AT SOLUTION

n  = 10
DO  i = 1,n
  x(i) = i / REAL(n+1, KIND=dp)
END DO
f  = 1.d0
CALL tn (ierror, n, x, f, g, sfun)
STOP
END PROGRAM No_Bounds


SUBROUTINE sfun (n, x, f, g)
IMPLICIT NONE
INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f
REAL (dp), INTENT(OUT)  :: g(:)

REAL (dp)  :: t
INTEGER    :: i

! ROUTINE TO EVALUATE FUNCTION (F) AND GRADIENT (G) OF THE OBJECTIVE
! FUNCTION AT THE POINT X

f = 0.d0
DO  i = 1,n
  t    = x(i) - i
  f    = f + t*t
  g(i) = 2 * t
END DO
RETURN
END SUBROUTINE sfun

PROGRAM Bounded

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-14  Time: 16:13:49
 
! EASY TO USE, BOUNDS ON VARIABLES
!***********************************************************************
! MAIN PROGRAM TO MINIMIZE A FUNCTION (REPRESENTED BY THE ROUTINE SFUN)
! SUBJECT TO BOUNDS ON THE VARIABLES X

USE Truncated_Newton
IMPLICIT NONE

REAL (dp)  :: x(50), f, g(50), low(50), up(50)
INTEGER    :: i, ierror, ipivot(50), n
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
! N   - NUMBER OF VARIABLES
! X   - INITIAL ESTIMATE OF THE SOLUTION
! LOW - LOWER BOUNDS
! UP  - UPPER BOUNDS
! F   - ROUGH ESTIMATE OF FUNCTION VALUE AT SOLUTION

n  = 10
DO  i = 1,n
  low(i) = 0.d0
  up(i)  = 6.d0
  x(i)   = i / REAL(n+1, KIND=dp)
END DO
f  = 1.d0
CALL tnbc (ierror, n, x, f, g, sfun, low, up, ipivot)
STOP
END PROGRAM Bounded


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

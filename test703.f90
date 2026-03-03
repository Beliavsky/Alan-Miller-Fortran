PROGRAM runmeb

USE common703
USE toms703
IMPLICIT NONE

INTERFACE
  SUBROUTINE fcn(t,y,ydot)
    USE common703, ONLY: dp
    IMPLICIT NONE
    REAL (dp), INTENT(IN)   :: t
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: ydot(:)
  END SUBROUTINE fcn
END INTERFACE

REAL (dp) :: y(20)

!    THE INCLUSION OF THE COMMON BLOCK  /MEBDF2/  ALLOWS THE USER TO
!    ACCESS VARIOUS COUNTERS THAT MAY BE OF INTEREST TO HIM.
!    THESE VARIABLES ARE EXPLAINED IN THE COMMENTS IN SUBROUTINE MEBDF.

INTEGER, PARAMETER :: lout = 6
REAL (dp)          :: err1, err2, err3, err4, hstart, h0, tol, x, xend, xout
INTEGER            :: index, n, nn

n = 4
DO nn=3,10
  INDEX = 1
  xend = 20.0_dp
  x = 0.0_dp
  y(1) = 1.0_dp
  y(2) = 1.0_dp
  y(3) = 1.0_dp
  y(4) = 1.0_dp
  IF (nn == 3 .OR. nn == 4) hstart = 1.0D-4
  IF (nn == 5 .OR. nn == 6) hstart = 1.0D-5
  IF (nn == 7 .OR. nn == 8) hstart = 1.0D-6
  IF (nn == 9 .OR. nn == 10) hstart = 1.0D-7

!   THESE ARE THE INITIAL STEPS CHOSEN BY THE DETEST PACKAGE FOR THIS PROBLEM

  h0 = hstart
  tol = 10.0_dp**(-nn)
  xout = 20.0_dp

  10 CALL mebdf(n, x, h0, y, xout, xend, tol, 21, INDEX, lout, fcn)
  IF ((INDEX /= 0) .AND. (INDEX /= 3)) THEN
    IF(INDEX == 1) THEN
      INDEX = 0
      GO TO 10
    END IF
    WRITE( 6 ,15) INDEX
    15 FORMAT('   ***INTEGRATION HAS FAILED***   WITH INDEX = ', i3)
    STOP
  END IF

!       HAVE COMPLETED ONE STEP

!       HAVE WE FINISHED YET?

  IF( INDEX == 0) THEN

!          THEN WE HAVE EFFECTIVELY HIT TOUT

    x = xout
  ELSE
    INDEX = 3
    GO TO 10
  END IF
  WRITE(6, 690) tol
  690 FORMAT(' REQUESTED TOLERANCE     ', G10.3)
  WRITE(6,700) nfe, nje

  700 FORMAT(' ', i5, '  FUNCTION EVALUATIONS          ',  &
             i5, '  JACOBIAN EVALUATIONS')
  err1 = y(1) - 4.539992969929191D-05
  err2 = y(2) - 2.061153036149920D-09
  err3 = y(3) - 2.823006338263857D-56
  err4 = y(4) - 5.235792540515384D-52

!       THESE VERY SMALL CONSTANTS SHOULD BE SET TO ZERO IF THERE
!       ARE LIKELY TO BE DIFFICULTIES DUE TO UNDERFLOW.

  WRITE(6,710)
  710 FORMAT(' ERRORS AT ENDPOINT')
  WRITE(6,720) err1, err2, err3, err4
  720 FORMAT(' ERRS', 5G14.5)
END DO
STOP

END PROGRAM runmeb


SUBROUTINE fcn(t,y,ydot)

USE common703, ONLY: dp
IMPLICIT NONE
REAL (dp), INTENT(IN)     :: t
REAL (dp), INTENT(IN)     :: y(:)
REAL (dp), INTENT(OUT)    :: ydot(:)

ydot(1) = -0.5_dp*y(1)
ydot(2) = -y(2)
ydot(3) = -100.0_dp*y(3)
ydot(4) = -90.0_dp*y(4)

! N.B. This application does not use argument t, but for ELF90:
IF (t < -9.9E99_dp) STOP

RETURN
END SUBROUTINE fcn


SUBROUTINE pderv(t, y, pw)

USE common703, ONLY: dp
IMPLICIT NONE
REAL (dp), INTENT(IN)     :: t
REAL (dp), INTENT(IN)     :: y(:,:)
REAL (dp), INTENT(OUT)    :: pw(:,:)     ! pw(4,4)

pw(1:4,1:4) = 0.0_dp
pw(1,1) = -0.5_dp
pw(2,2) = -1.0_dp
pw(3,3) = -100.0_dp
pw(4,4) = -90.0_dp

! N.B. This application does not use arguments t & y, but for ELF90:
IF (t < -9.9E99_dp) STOP
IF (y(1,1) < -9.9E99_dp) STOP

RETURN
END SUBROUTINE pderv

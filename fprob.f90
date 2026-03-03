SUBROUTINE fprob(num, denom, ratio, p, q, ier)
! Calculates the probability p that a variable with the F-distribution
! with num d.f. for the numerator, and denom d.f. for the denominator
! is less than the value ratio.   q gives the probability in the other
! tail, so that p + q = 1.

! Most of the code is translated from the Fortran 66 code for the incomplete
! beta function from the Naval Surface Warfare Center library and was written
! by Alfred Morris.

! This Fortran 90 version by Alan Miller
! e-mail: amiller@bigpond.net.au
! WWW:    http://users.bigpond.net.au/amiller

! Latest revision - 24 June 1997

USE constants_NSWC
USE inc_beta
IMPLICIT NONE
INTEGER, INTENT(IN)    :: num, denom
REAL (dp), INTENT(IN)  :: ratio
INTEGER, INTENT(OUT)   :: ier
REAL (dp), INTENT(OUT) :: p, q

!     Local variables
REAL (dp) :: zero = 0.D0, one = 1.D0, two = 2.D0, x

ier = 0
IF (num <= 0 .OR. denom <= 0) THEN
  ier = 1
  RETURN
ELSE IF (ratio < zero) THEN
  ier = 2
  RETURN
END IF

x = DBLE(denom) / (DBLE(denom) + DBLE(num)*ratio)
CALL bratio(DBLE(denom)/two, DBLE(num)/two, x, one-x, q, p, ier)

RETURN
END SUBROUTINE fprob



PROGRAM test_fprob
USE constants_NSWC
IMPLICIT NONE
INTEGER   :: num, denom, ier
REAL (dp) :: ratio, p, q

INTERFACE
  SUBROUTINE fprob(num, denom, ratio, p, q, ier)
    USE constants_NSWC
    USE inc_beta
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: num, denom
    REAL (dp), INTENT(IN)  :: ratio
    INTEGER, INTENT(OUT)   :: ier
    REAL (dp), INTENT(OUT) :: p, q
  END SUBROUTINE fprob
END INTERFACE

DO
  WRITE(*, *) 'Enter numerator then denominator deg. of freedom: '
  READ(*, *) num, denom
  WRITE(*, *) 'Enter variance ratio: '
  READ(*, *) ratio
  CALL fprob(num, denom, ratio, p, q, ier)
  IF (ier /= 0) THEN
    WRITE(*, *) 'Error from fprob: ', ier
  ELSE
    WRITE(*, '(a, 2g13.5)') ' Probabilities in each tail: ', p, q
  END IF
  WRITE(*, *)
END DO

STOP
END PROGRAM test_fprob

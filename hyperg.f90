MODULE hypergeometric

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER, SAVE      :: ier
PUBLIC             :: hyperg, lnfact, lngamma, ier

CONTAINS

FUNCTION hyperg(n) RESULT(prob)

!     Calculate hypergeometric probability for a 2 x 2 table:
!         n(1,1)  n(1,2)
!         n(2,1)  n(2,2)

!     Programmer: Alan Miller
!     http://users.bigpond.net.au
!     e-mail: amiller @ bigpond.net.au
!     Latest revision - 22 May 1997

INTEGER, INTENT(IN)  :: n(:,:)
REAL (dp)     :: prob

!     Local variables

INTEGER          :: ncell(1:4), ntot(1:4), nall
LOGICAL          :: ended
REAL (dp) :: zero = 0.d0, one = 1.d0

!     Check that all elements are positive or zero.

ier = 0
IF (n(1,1) < 0 .OR. n(1,2) < 0 .OR. n(2,1) < 0 .OR. n(2,2) < 0) THEN
  ier = 1
  prob = zero
  RETURN
END IF

!     Copy the cell frequencies into ncell and order them.
!     Calculate the row & column totals and order them.

ncell(1) = n(1,1)
ncell(2) = n(1,2)
ncell(3) = n(2,1)
ncell(4) = n(2,2)
CALL order4(ncell)

ntot(1) = n(1,1) + n(1,2)
ntot(2) = n(2,1) + n(2,2)
ntot(3) = n(1,1) + n(2,1)
ntot(4) = n(1,2) + n(2,2)
CALL order4(ntot)

nall = ncell(1) + ncell(2) + ncell(3) + ncell(4)
IF (nall == 0) THEN
  ier = 2
  prob = zero
  RETURN
END IF

!     Calculate hyperg as:

!       ntot(1)!   ntot(2)!    ntot(3)!    ntot(4)!        1
!       -------- . --------- . --------- . --------- . ---------
!        nall!     ncell(1)!   ncell(2)!   ncell(3)!   ncell(4)!

!     The order of multiplications and divisions minimizes the chance
!     of overflows or underflows without using logarithms.

IF (ncell(2) < 10) THEN                ! Use divisions & multiplications
  prob = one
  DO
    ended = .true.
    IF (nall > ntot(1)) THEN
      prob = prob / nall
      nall = nall - 1
      ended = .false.
    END IF

    IF (ntot(2) > ncell(1)) THEN
      prob = prob * ntot(2)
      ntot(2) = ntot(2) - 1
      ended = .false.
    END IF

    IF (ntot(3) > ncell(2)) THEN
      prob = prob * ntot(3)
      ntot(3) = ntot(3) - 1
      ended = .false.
    END IF

    IF (ntot(4) > ncell(3)) THEN
      prob = prob * ntot(4)
      ntot(4) = ntot(4) - 1
      ended = .false.
    END IF

    IF (ncell(4) > 1) THEN
      prob = prob / ncell(4)
      ncell(4) = ncell(4) - 1
      ended = .false.
    END IF
    IF (ended) EXIT
  END DO
ELSE                                   ! Use logs of factorials for speed
  prob = EXP( lnfact(ntot(1)) + lnfact(ntot(2)) + lnfact(ntot(3)) +   &
              lnfact(ntot(4)) - lnfact(nall) - lnfact(ncell(1)) -     &
              lnfact(ncell(2)) - lnfact(ncell(3)) - lnfact(ncell(4)) )
END IF

RETURN
END FUNCTION hyperg



SUBROUTINE order4(n)

!     Order 4 integers; largest first.

INTEGER, INTENT(IN OUT) :: n(:)

!     Local variable
INTEGER :: nn

IF (n(1) < n(2)) THEN
  nn = n(1)
  n(1) = n(2)
  n(2) = nn
END IF

IF (n(3) < n(4)) THEN
  nn = n(3)
  n(3) = n(4)
  n(4) = nn
END IF

IF (n(2) < n(3)) THEN
  nn = n(2)
  n(2) = n(3)
  n(3) = nn
  IF (n(3) < n(4)) THEN
    nn = n(3)
    n(3) = n(4)
    n(4) = nn
  END IF
  IF (n(1) < n(2)) THEN
    nn = n(1)
    n(1) = n(2)
    n(2) = nn
    IF (n(2) < n(3)) THEN
      nn = n(2)
      n(2) = n(3)
      n(3) = nn
    END IF
  END IF
END IF

RETURN
END SUBROUTINE order4



FUNCTION lnfact(n) RESULT(fn_val)
! Calculate the logarithm of n!

INTEGER, INTENT(IN) :: n
REAL (dp)    :: fn_val

! Local variable
REAL (dp) :: zero = 0.0D0

SELECT CASE (n)
  CASE (0:1)
    fn_val = zero
  CASE DEFAULT
    fn_val = lngamma(DBLE(n+1))
END SELECT

RETURN
END FUNCTION lnfact



FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma function',
!               J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
		       -1259.139216722289D0, 771.3234287757674D0, &
		       -176.6150291498386D0, 12.50734324009056D0, &
		       -0.1385710331296526D0, 0.9934937113930748D-05, &
			0.1659470187408462D-06 /), zero = 0.D0,   &
			one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
			half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER   :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma

END MODULE hypergeometric



PROGRAM t_hyperg

!     Test function hyperg for calculating hypergeometric probabilities

USE hypergeometric
IMPLICIT NONE
INTEGER  :: n(2,2)

DO
  WRITE(*, *)'Enter 4 integers: '
  READ(*, *) n
  IF (ier == 0) THEN
    WRITE(*, 900) hyperg(n)
    900 FORMAT(' Prob. = ', g14.6/)
  ELSE
    WRITE(*, *) 'Error in input values'
  END IF
END DO

STOP
END PROGRAM t_hyperg

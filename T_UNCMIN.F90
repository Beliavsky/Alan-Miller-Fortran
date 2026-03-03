SUBROUTINE fcn(n, x, f)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp), INTENT(OUT) :: f

REAL (dp), PARAMETER :: one = 1.0_dp

IF (n /= 2) STOP
f = 100._dp*(x(2) - x(1)**2)**2 + (one - x(1))**2
RETURN
END SUBROUTINE fcn


SUBROUTINE d1fcn(n, x, g)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp), INTENT(OUT) :: g(:)

REAL (dp), PARAMETER :: one = 1.0_dp

IF (n /= 2) STOP
g(1) = -400._dp*(x(2) - x(1)**2)*x(1) - 2._dp*(one - x(1))
g(2) =  200._dp*(x(2) - x(1)**2)
RETURN
END SUBROUTINE d1fcn


SUBROUTINE d2fcn(nr, n, x, h)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
INTEGER, INTENT(IN)    :: nr
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp), INTENT(OUT) :: h(:,:)

IF (nr /= 2) STOP
IF (n /= 2) STOP
h(1,1) = -400._dp*(x(2) - 3._dp*x(1)**2) + 2._dp
h(1,2) = -400._dp*x(1)
h(2,1) = -400._dp*x(1)
h(2,2) =  200._dp
RETURN
END SUBROUTINE d2fcn


PROGRAM test_uncmin
! Test the unconstrained minimization package UNCMIN

USE unconstrained_min
IMPLICIT NONE
REAL (dp) :: x(2), typsiz(2), fscale, dlt, gradtl, stepmx, steptl, xpls(2), &
             fpls, gpls(2)
INTEGER   :: nr, n, m, method, iexp, msg, ndigit, itnlim, iagflg, iahflg,  &
             itrmcd

REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

INTERFACE
  SUBROUTINE fcn(n, x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: x(:)
    REAL (dp), INTENT(OUT) :: f
  END SUBROUTINE fcn

  SUBROUTINE d1fcn(n, x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: x(:)
    REAL (dp), INTENT(OUT) :: g(:)
  END SUBROUTINE d1fcn

  SUBROUTINE d2fcn(nr, n, x, h)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: nr
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: x(:)
    REAL (dp), INTENT(OUT) :: h(:,:)
  END SUBROUTINE d2fcn
END INTERFACE


DO m = 1, 3
  method = m
  nr = 2
  n = 2
  x(1) = zero
  x(2) = zero
  typsiz = one
  fscale = 10._dp
  iexp = 0
  msg = 0
  ndigit = 12
  itnlim = 500
  iagflg = 1
  iahflg = 1
  dlt = one
  gradtl = 1.E-4_dp
  stepmx = one
  steptl = 1.E-4_dp
  WRITE(*, *) 'METHOD:', method
  CALL optdrv(nr, n, x, fcn, d1fcn, d2fcn, typsiz, fscale,  &
              method, iexp, msg, ndigit, itnlim, iagflg, iahflg, &
              dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd)
  WRITE(*, *) 'ITRMCD =', itrmcd
  WRITE(*, *) 'MSG    =', msg
  WRITE(*, '(a, 2f12.8)') ' Final X: ', xpls
  WRITE(*, '(a, g13.5)') ' Function value: ', fpls
  WRITE(*, '(a, 2g13.5)') ' Gradients: ', gpls
  WRITE(*, *)
END DO


STOP

END PROGRAM test_uncmin

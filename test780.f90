MODULE common780
IMPLICIT NONE

REAL, SAVE, PUBLIC :: x

END MODULE common780



PROGRAM DRIVER
!
!     The return variable, x, is in COMMON to make certain that no
!     overly-enthuiastic optimizer eliminates the function calls.
!
USE common780
IMPLICIT NONE

INTEGER :: i, mtime0, mtime1, n
REAL    :: dt_lg, dta_lg, dt_eau, dta_eau, dt_eas, dta_eas

INTERFACE
  FUNCTION REXPU() RESULT(fn_val)
    IMPLICIT NONE
    REAL :: fn_val
  END FUNCTION REXPU

  FUNCTION REXPS() RESULT(fn_val)
    IMPLICIT NONE
    REAL :: fn_val
  END FUNCTION REXPS
END INTERFACE

!
n = 1000000
call msec(mtime0)
do i=1,n
  x = rexp_lg()
END DO
call msec(mtime1)
dt_lg = 0.001 * REAL(mtime1 - mtime0)
dta_lg = 1.E6 * dt_lg / REAL(n)
WRITE(*, *) 'Old logarithmic method:'
WRITE(*, *) '  LG Run time =', dt_lg, ' seconds'
WRITE(*, *) '  LG Average time =', dta_lg
!
call msec(mtime0)
do i=1,n
  x = rexpu()
END DO
call msec(mtime1)
dt_eau = 0.001 * REAL(mtime1 - mtime0)
dta_eau = 1.E6 * dt_eau / REAL(n)
WRITE(*, *) 'Unstructured:'
WRITE(*, *) '  EA Run time =', dt_eau, ' seconds'
WRITE(*, *) '  EA Average time =', dta_eau
!
call msec(mtime0)
do i=1,n
  x = rexps()
END DO
call msec(mtime1)
dt_eas = 0.001 * REAL(mtime1 - mtime0)
dta_eas = 1.E6 * dt_eas / REAL(n)
WRITE(*, *) 'Structured:'
WRITE(*, *) '  EA Run time =', dt_eas, ' seconds'
WRITE(*, *) '  EA Average time =', dta_eas
!
stop

CONTAINS


SUBROUTINE MSEC(M)

INTEGER, INTENT(OUT) :: m

! Local variable

INTEGER :: iv(8)

call date_and_time(values=iv)                    ! F90 intrinsic
m = iv(8) + 1000*(iv(7) + 60*(iv(6) + 60*iv(5)))
return
end SUBROUTINE MSEC



FUNCTION REXP_LG() RESULT(fn_val)                ! Old method for comparison

REAL :: fn_val

! Local variable

REAL :: u

call random_number(u)
fn_val = -alog(u)

return
end FUNCTION REXP_LG

end PROGRAM DRIVER


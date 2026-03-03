PROGRAM hh_test
! Test of half-Hermite function
! The exact value for the integral is 1.0 for all values of a

IMPLICIT DOUBLE PRECISION (a-h,o-z)
EXTERNAL func
COMMON a

a = 0.5D0
DO
  DO npt = 2,20
    z = hh(npt,func,ier)
    IF(ier == 0) WRITE(*, 910) a, npt, z
    910 FORMAT(' A =', f5.1, '  NPT =',i3,'  INTEGRAL =',f20.15)
  END DO
  a = a + 0.5D0
  IF(a > 5.01D0) EXIT
END DO
STOP
END PROGRAM hh_test



DOUBLE PRECISION FUNCTION func(x)
IMPLICIT DOUBLE PRECISION (a-h,o-z)
COMMON a

func = 2.d0*(a+x)*EXP(-2.0*a*x)
RETURN
END FUNCTION func

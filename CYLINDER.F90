PROGRAM cylinder
! Fit a cylinder:

! A vertical cylinder centred at zero, x^2 + y^2 = r^2 (any z)
! is first rotated in the (x,z) plane, giving
! [x.cos(theta) + z.sin(theta)]^2 + y^2 = r^2
! then rotated in the (y,z) plane giving
! [x.cos(theta) + z.sin(theta).cos(alpha) - y.sin(theta).sin(alpha)]^2 +
!                 [y.cos(alpha) + z.sin(alpha)]^2 = r^2
! The cylinder is then moved in the (x,y) plane by replacing (x,y)
! with (x-x0, y-y0).

! Latest revision - 2 August 2000
! Alan Miller
! amiller @ bigpond.net.au

USE Errors_in_Variables
IMPLICIT NONE

INTEGER, PARAMETER  :: ndat = 50, npar = 5, nvar = 3
INTEGER    :: i, ifault(2)
REAL (dp)  :: bs(1,nvar), DATA(ndat,nvar), est(ndat,nvar), g1(npar,npar), &
              g3(npar,npar), theta(npar), v(nvar,nvar), sine1, cos1,  &
              sine2, cos2, e1, e2, e3, scale, x, y, z, resid(ndat,nvar)

v = 0.0_dp
DO i = 1, nvar
  v(i,i) = 1.0_dp
END DO

! Create artificial data
! True parameter values: 1.0, 2.0, 3.0, 0.6, 1.2

sine1 = SIN(0.6)
cos1  = COS(0.6)
sine2 = SIN(1.2)
cos2  = COS(1.2)

DO i = 1, ndat
! Generate e1, e2 such that e1^2 + e2^2 = 9
  CALL RANDOM_NUMBER(e1)
  CALL RANDOM_NUMBER(e2)
  e1 = e1 - 0.5_dp
  e2 = e2 - 0.5_dp
  scale = SQRT(9.0 / (e1**2 + e2**2))
  e1 = e1*scale
  e2 = e2*scale
! From e1, e2, calculate (x,y,z) which will be stored in DATA(i,1:3)
! e1 = x.cos1 + (z.cos2 - y.sine2).sine1
! e2 = y.cos2 + z.sine2
  CALL RANDOM_NUMBER(z)
  z = 5.0_dp * z
  y = (e2 - z*sine2) / cos2
  x = (e1 - (z*cos2 - y*sine2)*sine1) / cos1
  DATA(i,1) = x + 1.0
  DATA(i,2) = y + 2.0
  DATA(i,3) = z
END DO

! Now perturb the data slightly
DO i = 1, ndat
  CALL RANDOM_NUMBER(e1)
  CALL RANDOM_NUMBER(e2)
  CALL RANDOM_NUMBER(e3)
  DATA(i,1) = DATA(i,1) + 0.25_dp*(e1 - 0.5_dp)
  DATA(i,2) = DATA(i,2) + 0.25_dp*(e2 - 0.5_dp)
  DATA(i,3) = DATA(i,3) + 0.25_dp*(e3 - 0.5_dp)
END DO

! Set starting values

theta = (/ 1.5_dp, 2.5_dp, 3.5_dp, 1.0_dp, 1.0_dp /)

CALL evms (bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)

WRITE(*, '(a, 2i5)') ' Error status = ', ifault
WRITE(*, '(a / 5f9.4)') ' Parameter estimates: ', theta
WRITE(*, '(/a)') ' Information matrix:'
WRITE(*, '(5f9.2)') g1
WRITE(*, '(/a)') ' Inverse of Information matrix:'
WRITE(*, '(5f10.4)') g3
WRITE(*, *)
WRITE(*, *) '   X      Y      Z           Estimates              Residuals'
DO i = 1, ndat
  WRITE(*, '(3f7.2, "  ", 3f7.2, "  ", 3f7.3)')  &
           DATA(i,1:3), est(i,1:3), resid(i,1:3)
END DO

STOP
END PROGRAM cylinder



SUBROUTINE bf (b, f, theta, xi)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(OUT)  :: b(:,:)
REAL (dp), INTENT(OUT)  :: f(:)
REAL (dp), INTENT(IN)   :: theta(:)
REAL (dp), INTENT(IN)   :: xi(:)

! The parameters in theta are:
! theta(1) = x0
! theta(2) = y0
! theta(3) = r
! theta(4) = angle 1
! theta(5) = angle 2

! The (x,y,z) co-ordinates are xi(1), xi(2) & xi(3)

REAL (dp)  :: e1, e2, sine1, cos1, sine2, cos2, x, y, z
REAL (dp), PARAMETER  :: two = 2.0_dp

! Calculate F
sine1 = SIN(theta(4))
cos1  = COS(theta(4))
sine2 = SIN(theta(5))
cos2  = COS(theta(5))
x = xi(1) - theta(1)
y = xi(2) - theta(2)
z = xi(3)
e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
e2 = y*cos2 + z*sine2
f(1) = e1**2 + e2**2 - theta(3)**2

! Calculate derivatives w.r.t. x, y, z
b(1,1) = two*e1*cos1
b(1,2) = two*e2*cos2 - two*e1*sine2*sine1
b(1,3) = two*e1*cos2*sine1 + two*e2*sine2

RETURN
END SUBROUTINE bf



SUBROUTINE zed (b, f, theta, xi, grad)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: b(:,:)
REAL (dp), INTENT(IN)   :: f(:)
REAL (dp), INTENT(IN)   :: theta(:)
REAL (dp), INTENT(IN)   :: xi(:)
REAL (dp), INTENT(OUT)  :: grad(:,:)

REAL (dp)  :: e1, e2, sine1, cos1, sine2, cos2, x, y, z
REAL (dp), PARAMETER  :: two = 2.0_dp

sine1 = SIN(theta(4))
cos1  = COS(theta(4))
sine2 = SIN(theta(5))
cos2  = COS(theta(5))
x = xi(1) - theta(1)
y = xi(2) - theta(2)
z = xi(3)
e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
e2 = y*cos2 + z*sine2

! Calculate derivatives w.r.t. parameters
grad(1,1) = -two*e1*cos1
grad(1,2) = -two*e2*cos2 + two*e1*sine2*sine1
grad(1,3) = -two*theta(3)
grad(1,4) = two*e1*(-x*sine1 + (z*cos2 - y*sine2)*cos1)
grad(1,5) = two*e2*(-y*sine2 + z*cos2) + two*e1*(-z*sine2 - y*cos2)*sine1

RETURN
END SUBROUTINE zed


MODULE common_nf
IMPLICIT NONE

!     .. Common blocks ..
! COMMON /nf/nfun

INTEGER, SAVE  :: nfun

END MODULE common_nf



! ***********************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-12  Time: 17:50:27

!                                              **
! FUNCTION FZ CALLED BY THE INVERSION SOFTWARE **
!                                              **
! ***********************************************

FUNCTION fz(z) RESULT(fn_val)
!     .. Scalar Arguments ..

USE common_nf
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

COMPLEX (dp), INTENT(IN)  :: z
COMPLEX (dp)              :: fn_val

!     ..
!     .. Local Scalars ..
COMPLEX (dp)  :: b, ci, cvar, arg, z2, z4, a1, a2, bb
REAL (dp)     :: r, c
!     ..
!     .. External Functions ..
! INTRINSIC EXP,CDLOG,SQRT
!     ..

! GO TO(1,2,3,4,5,6,7,8,10,11,12,13,14,15, 16,17,18,19,20,  &
!       21,22,23,24,25,26,27,28,29,30,31,32, 33,34,35,36,37,38) nfun

! ************************************************
! FUNCTIONS ORDER AND NUMBERING REFLECT FUNCTIONS
! ORDER AND NUMBERING AS IN THE PAPER:

!   D'AMORE L., LACCETTI G., MURLI A.,  -

!  "ALGORITHM XXX: A FORTRAN SOFTWARE
!   PACKAGE FOR THE NUMERICAL INVERSION OF THE
!   LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"

!   ACM TRANS. MATH. SOFTWARE, VOL. ##,
!   NO. #, MONTH YEAR, PP. ##-##.
! *******************************************************

SELECT CASE (nfun)
  CASE (1)
    fn_val = (1.d0,0.d0)/z

  CASE (2)
    fn_val = (2.d0,0.d0)*(SQRT(z + (1.d0,0.d0)) - SQRT(z))

  CASE (3)
    fn_val = (1.d0,0.d0)/SQRT(z)

  CASE (4)
    fn_val = (z*z - (1.d0,0.d0))/ ((z*z + (1.d0,0.d0))**2)

  CASE (5)
    fn_val = (1.d0,0.d0)/ (z + (1.d0,0.d0))**2

  CASE (6)
   fn_val = (1.d0,0.d0)/ (z**2)

  CASE (7)
    fn_val = (1.d0,0.d0)/ (z**2 + (1.d0,0.d0))

  CASE (8)
    fn_val = (1.d0,0.d0)/ (z + (0.5D0,0.d0))

! ************************************************
!       Such inverse function is the Bessel function J_0.
!       Actually we compute it by using the Nag library.  That's why
!       it does not appear here.

! 9  fn_val = 1.d0/SQRT(z*z + (1.d0,0.d0))
!

  CASE (10)
    fn_val = EXP((-1.d0,0.d0)/z)/SQRT(z)

  CASE (11)
    fn_val = EXP((-4.d0,0.d0)*SQRT(z))

  CASE (12)
    ci = (0.d0,1.d0)
    cvar = 1.d0 / (0.d0, 2.d0)
    fn_val = cvar*LOG((z+ci)/ (z-ci))

  CASE (13)
    fn_val = (1.d0,0.d0) / ((z +.2_dp)**2 + (1.d0,0.d0))

  CASE (14)
    fn_val =  1.d0/z**3

  CASE (15)
    fn_val = EXP(-2*z)/z

  CASE (16)
    fn_val = (1.d0,0.d0)/ (z* ((1.d0,0.d0) + EXP(-z)))

  CASE (17)
    fn_val = (1.d0,0.d0)/ (z*z + z + (1.d0,0.d0))

  CASE (18)
    fn_val = (3.d0,0.d0)/ (z**2 - (9.d0,0.d0))

  CASE (19)
    fn_val = (120.d0,0.d0)/z**6

  CASE (20)
    fn_val = z / (z**2 + (1.d0,0.d0))**2

  CASE (21)
    fn_val = (1.d0,0.d0) / (z + (1.d0,0.d0)) - (1.d0,0.d0)/ (z + (1000.d0,0.d0))

  CASE (22)
    fn_val = z / (z*z + (1.d0,0.d0))

  CASE (23)
    fn_val = (1.d0,0.d0) / ((z - (0.25D0,0.d0))**2)

  CASE (24)
    fn_val = (1.d0,0.d0) / (z*SQRT(z))

  CASE (25)
    fn_val = (1.d0,0.d0) / SQRT(z + (1.d0,0.d0))

  CASE (26)
    fn_val = (z + (2.d0,0.d0)) / (z*SQRT(z))

  CASE (27)
    fn_val = (1.d0,0.d0)/ ((z*z + (1.d0,0.d0))**2)

  CASE (28)
    fn_val = (1.d0,0.d0)/ (z* (z + (1.d0,0.d0))**2)

  CASE (29)
    fn_val = (1.d0,0.d0)/ (z**3 - (8.d0,0.d0))

  CASE (30)
    fn_val = LOG((z*z + (1.d0,0.d0))/ (z*z + (4.d0,0.d0)))

  CASE (31)
    fn_val = LOG((z + (1.d0,0.d0))/z)

  CASE (32)
    fn_val = LOG(z)/z

  CASE (33)
    fn_val = ((1.d0,0.d0) - EXP(-z)) / (z*z)

  CASE (34)
    fn_val = (1.d0,0.d0) / (z* ((1.d0,0.d0) + EXP(z)))

  CASE (35)
    b = (1.d0,0.d0)/ (2.*z) - (EXP(-2.*z)/ (1. - EXP(-2.*z)))
    fn_val = (1./ (z*z + z))*b

  CASE (36)
    r = 0.5D0
    c = 0.4D0
    b = -r*SQRT((z*(1 + z)) / (1 + c*z))
    fn_val = 1./z*EXP(b)

  CASE (37)
!         FZ = EXP(-2.*PSI)/Z
!   COSH(PSI)= SQRT(1 + Z**2 + (Z**4)/16)

    z2 = z*z
    z4 = z2*z2
    arg = 1.d0 + z2 + z4/16
    a1 = SQRT(arg)
    a2 = (z/4.d0)*SQRT(16.d0 + z2)
    fn_val = 1.d0/(z*(a1 + a2)**2)

  CASE (38)
    b = z - SQRT(z*z - 1.d0)
    bb = SQRT(z)*SQRT(z*z - 1.d0)*SQRT(z - 0.5D0*(SQRT(z*z - 1.d0)))
    fn_val = b/bb

END SELECT

RETURN
END FUNCTION fz


! *********************************************************************

!    FZ's COMPANION FUNCTION FEX TO COMPUTE THE EXACT VALUE OF THE   **
!    INVERSE TRANSFORM                                               **

! *********************************************************************

FUNCTION fex(x) RESULT(fn_val)

USE common_nf
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

!     .. Scalar Arguments ..

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), SAVE  :: comx
INTEGER, SAVE    :: kount

!     ..
!     .. Local Scalars ..
REAL (dp) :: a, b, pi, pi2, sum, t
REAL (dp) :: eulero
INTEGER   :: k, n
LOGICAL   :: trov
!     .. Intrinsic Functions ..
! INTRINSIC DATAN, DCOS, EXP, DSIN, DSINH, SQRT, dfloat
!     ..

! REAL (dp) :: fst3,fst4,fst52,fst51
! EXTERNAL     fst3,fst4,fst52,fst51

pi = 4.d0*ATAN(1.d0)
comx = x
eulero = .5772156649015329_dp

! GO TO(1,2,3,4,5,6,7,8,10,11,12,13,14, 15,16,17,18,19,20,21,22,  &
!    23,24,25,26,27,28,29,30,31,32,33,34,  35,36,37,38) NFUN

! ********************************************************************
! FUNCTIONS ORDER AND NUMBERING REFLECT FUNCTIONS
! ORDER AND NUMBERING AS IN THE PAPER:

!   D'AMORE L., LACCETTI G., MURLI A.,  -

!  "ALGORITHM XXX: A FORTRAN SOFTWARE
!   PACKAGE FOR THE NUMERICAL INVERSION OF THE
!   LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"

!   ACM TRANS. MATH. SOFTWARE, VOL. ##,
!   NO. #, MONTH YEAR, PP. ##-##.
! ********************************************************************

SELECT CASE (nfun)
  CASE (1)
    fn_val = 1.d0

  CASE (2)
    fn_val = (1.d0 - EXP(-x))/ (x*SQRT(pi*x))

  CASE (3)
    fn_val = 1.d0/SQRT(pi*x)

  CASE (4)
    fn_val = x*COS(x)

  CASE (5)
    fn_val = x*EXP(-x)

  CASE (6)
    fn_val = x

  CASE (7)
    fn_val = SIN(x)

  CASE (8)
    fn_val = EXP(-.5D0*x)

! ************************************************************
!       Such inverse function is the Bessel function J_0.
!       Actually we compute it by using the Nag library.  That's why
!       it does not appear here .

! 9   fn_val = S17AEF(X,IFAIL)
!
! ************************************************************

  CASE (10)
    fn_val = COS(2.d0*SQRT(x))/SQRT(pi*x)


  CASE (11)
    fn_val = 2.d0*EXP(-4.d0/x) / (x*SQRT(pi*x))

  CASE (12)
    fn_val = SIN(x)/x

  CASE (13)
    fn_val = EXP(-.2D0*x)*SIN(x)

  CASE (14)
    fn_val = 0.5D0 * x**2

  CASE (15)
   IF (x > 2.d0) THEN
     fn_val = 1.d0

   ELSE IF (x < 2.d0) THEN
     fn_val = 0.d0

   ELSE
     fn_val = 0.5D0
   END IF

  CASE (16)
    trov = .false.
    k = 0
    310 IF (x > 2*k .AND. x < 2*k+1) THEN
      fn_val = 1.d0
      trov = .true.

    ELSE IF (x > 2*k+1 .AND. x < 2*k+2) THEN
      fn_val = 0.d0
      trov = .true.
    END IF

    IF (x == DBLE(k)) THEN
      fn_val = 0.5D0
      trov = .true.
    END IF

    k = k + 1
    IF (.NOT.trov .AND. k <= 49) GO TO 310

  CASE (17)
    fn_val = 2.d0/SQRT(3.d0)*EXP(-x/2.d0)*SIN(x*SQRT(3.d0)/2.d0)

  CASE (18)
    fn_val = SINH(3.d0*x)

  CASE (19)
    fn_val = x**5

  CASE (20)
    fn_val = x/2.d0*SIN(x)

  CASE (21)
    fn_val = EXP(-x) - EXP(-1000.d0*x)

  CASE (22)
    fn_val = COS(x)

  CASE (23)
    fn_val = x*EXP(x/4.)

  CASE (24)
    fn_val = 2.d0*SQRT(x/pi)

  CASE (25)
    fn_val = EXP(-x)/SQRT(pi*x)

  CASE (26)
    fn_val = (1.d0 + 4.d0*x)/SQRT(pi*x)

  CASE (27)
    fn_val = (SIN(x) - x*COS(x))/2.d0

  CASE (28)
    fn_val = 1. - EXP(-x)* (1. + x)

  CASE (29)
    fn_val = EXP(-x)/12.d0* (EXP(3.d0*x) - COS(SQRT(3.d0)*x) -  &
             SQRT(3.d0)*SIN(SQRT(3.d0)*x))

  CASE (30)
    fn_val = 2.d0* (COS(2.d0*x) - COS(x))/x

  CASE (31)
    fn_val = (1. - EXP(-x))/x


  CASE (32)
    fn_val = -eulero - LOG(x)

  CASE (33)
    IF (x >= 0.d0 .AND. x <= 1.d0) THEN
      fn_val = x
    ELSE
      fn_val = 1.d0
    END IF

  CASE (34)
    trov = .false.
    k = 0

    210 IF (x > 2*k .AND. x < 2*k+1) THEN
      fn_val = 0.d0
      trov = .true.

    ELSE IF (x > 2*k+1 .AND. x < 2*k+2) THEN
      fn_val = 1.d0
      trov = .true.
    END IF

    IF (x == DBLE(k)) THEN
      fn_val = 0.5D0
      trov = .true.
    END IF

    k = k + 1
    IF (.NOT.trov .AND. k <= 49) GO TO 210

  CASE (35)
    t = x
    pi2 = pi*pi
    a = (0.5D0- (EXP(2.d0)/ (EXP(2.d0)-1.))*EXP(-t)) + 0.5D0
    n = 0.d0
    sum = 0.d0

    370 n = n + 1.d0
    b = SIN((n*pi*t) - ATAN(n*pi))
    b = b/ (n* (SQRT((n*n*pi2) + 1.d0)))
    sum = sum + b
    IF (n < 2.d+4) GO TO 370
    sum = sum/pi
    sum = sum*EXP(-t)
    fn_val = a - sum

!
!36    EPSABS = 0.0D0
!       EPSREL = 1.0D-04
!       A = 0.0D0
!       INF = 1
!       KOUNT = 0
!       IFAIL = -1
!       CALL D01AMF(FST3,A,INF,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
!     +            IFAIL)
!       IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
!       fn_val=(RESULT/PI)+0.5d0
!       goto 380
!99996 FORMAT (1X,A,I4)


!37    EPSABS = 0.0D0
!      EPSREL = 1.0D-04
!      U1=2.*SQRT(2.D0-SQRT(3.D0))
!      U2=2.*SQRT(2.D0+SQRT(3.D0))
!      A = .0D0
!      B =U1
!      KOUNT = 0
!      IFAIL = -1
!      CALL D01AKF(FST4,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
!     *            IFAIL)
!      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
!      DUFFY6=RESULT
!      A=U2
!      B=4.D0
!      IFAIL=-1
!       KOUNT=0
!      CALL D01AKF(FST4,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
!     *            IFAIL)
!      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
!      fn_val=(-DUFFY6+RESULT)/PI +1.D0
!
!

!38    EPSABS = 0.0D0
!      EPSREL = 1.0D-04
!      A = 0.0D0
!      B = (1.-0.5d0)/0.5d0
!      KOUNT = 0
!      IFAIL = -1
!      CALL D01AJF(FST51,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
!    * IFAIL)
!      IF (IFAIL.NE.0) WRITE (NOUT,99996) 'IFAIL = ', IFAIL
!      PARTE1=RESULT
!      A=0.d0
!      B=SQRT((1.-0.5d0)/1.5d0)
!      KOUNT=0
!      IFAIL=-1
!      CALL D01AJF(FST52,A,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,LIW,
!    * IFAIL)
!      fn_val=(RESULT+PARTE1)*(2.D0/PI)

END SELECT

RETURN

CONTAINS



FUNCTION fst51(u) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: u
REAL (dp)              :: fn_val

REAL (dp) :: n, c
REAL (dp) :: r

! INTRINSIC   SIN, SQRT
! COMMON      /telnum/comx,kount

kount = kount + 1
n = 0.5D0
c = (1. - n)/n
r = SQRT(u*u + n*n*(c*c - u*u))
fn_val = COSH(comx*u)* ((u*SQRT((r + u)/2.)) +  &
         SQRT(c**2 - u**2)*SQRT((r-u)/2.)) / (r*SQRT(c**2 - u**2)*SQRT(u))

RETURN
END FUNCTION fst51



FUNCTION fst52(u) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: u
REAL (dp)              :: fn_val

REAL (dp) :: n, c

! INTRINSIC   SIN, SQRT
! COMMON      /telnum/comx,kount

kount = kount + 1
n = 0.5D0
c = (1. - n)/n
fn_val = COS(comx*u)*((u - SQRT(c**2 + u**2))/(SQRT(u)*  &
         SQRT(c**2 + u**2)*SQRT(n*SQRT(c**2 + u**2) - u)))

RETURN
END FUNCTION fst52



FUNCTION fst3(u) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: u
REAL (dp)              :: fn_val

REAL (dp) :: m, theta, arg1, arg2, r, c

! INTRINSIC  SQRT
! COMMON     /telnum/comx,kount

kount = kount + 1
r = 0.5D0
c = 0.4D0
m = (1.d0 + u*u)/(1.d0 + c*c*u*u)
m = m**0.25
theta = ATAN(u) - ATAN(u*c)
theta = theta/2
arg1 = -r*m*SQRT(u/2)*(COS(theta) - SIN(theta))
arg2 = comx*u - r*m*SQRT(u/2)*(COS(theta) + SIN(theta))
fn_val = EXP(arg1)*SIN(arg2)/u

RETURN
END FUNCTION fst3




FUNCTION fst4(u) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: u
REAL (dp)              :: fn_val

REAL (dp) :: u1, u2, k

! INTRINSIC  SIN, SQRT,ACOS
! COMMON     /telnum/comx,kount

kount = kount + 1
u1 = 4.*(2.d0 - SQRT(3.d0))
u2 = 4.*(2.d0 + SQRT(3.d0))
IF((u1 - u**2 )*( u2 - u**2) >= 0.d0) THEN
  k = ACOS(.25D0*SQRT((u1 - u**2)*(u2 - u**2)))
  fn_val = (SIN(u*comx + 2.*k) - SIN(u*comx - 2*k))/u
END IF

RETURN
END FUNCTION fst4

END FUNCTION fex


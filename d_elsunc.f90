MODULE common_ex25
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /DATA/ u(99)
REAL (dp), SAVE  :: u(99)

END MODULE common_ex25



PROGRAM ex25
! Example driver for ELSUNC
! EX25 HOCK,SCHITTKOWSKI
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-04-29  Time: 15:46:48
 
!     THIS IS A MAIN PROGRAM FOR SOLVING LEAST SQUARES PROBLEMS
!     WHERE THE UNKNOWNS ARE RESTRICTED BY LOWER AND UPPER BOUNDS.
!     THE ROUTINE ELSUNC IS USED AS SOLVER.

USE bounded_nonlin_LS
USE common_ex25, ONLY: u
IMPLICIT NONE

! COMMON /DATA/ u(99)
! EXTERNAL bndfun

INTERFACE
  SUBROUTINE bndfun(x, n, f, m, ctrl, c, mdc)
    IMPLICIT NONE
    INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)    :: x(:)
    INTEGER, INTENT(IN)      :: n
    REAL (dp), INTENT(OUT)   :: f(:)
    INTEGER, INTENT(IN)      :: m
    INTEGER, INTENT(IN OUT)  :: ctrl
    REAL (dp), INTENT(OUT)   :: c(:,:)
    INTEGER, INTENT(IN)      :: mdc
  END SUBROUTINE bndfun
END INTERFACE

INTEGER   :: EXIT,i,mdc,mdw,n,m,nout,bnd
INTEGER   :: p(17)
REAL (dp) :: x(3),f(99),c(99,4),w(327),bl(3),bu(3),two3rd

mdc=99
mdw=327
m=99
n=3
nout=10
OPEN(UNIT=nout, FILE='ex25.out')

!     INITIATE THE PROBLEM VECTOR U

two3rd=2.0_dp/3.0_dp
DO  i=1,m
  u(i)=(-50.0*LOG(REAL(i)*0.01_dp))**two3rd + 25.0_dp
END DO

!     GIVE LOWER AND UPPER BOUNDS

bnd=2
bl(1)=0.1
bu(1)=100.0
bl(2)=0.0
bu(2)=25.6
bl(3)=0.0
bu(3)=5.0

!     GIVE STARTING POINT

x(1)=100.0
x(2)=12.5
x(3)=3.0

!     SET THE FIRST ELEMENTS OF THE ARRAYS P AND W TO SOMETHING
!     NEGATIVE TO INDICATE THAT DEFAULT VALUES SHOULD BE USED

DO  i=1,5
  p(i)=-1
  w(i)=-1.0D0
END DO
WRITE(*,*) 'How many iterations?'
READ(*,*) p(3)

!     DO SOME INITIAL PRINTING

WRITE(nout,1000)
1000 FORMAT(/t11, 'EX25 hock schittkowski' /)
WRITE(nout,1010) x(1:n)
1010 FORMAT(' ', 'STARTING point =', 5g12.3)
WRITE(nout,1020) bl(1:n)
1020 FORMAT(' ', 'LOWER bounds =', 5g12.3)
WRITE(nout,1030) bu(1:n)
1030 FORMAT(' ', 'UPPER bounds =', 5g12.3)

!     INVOKE ELSUNC TO MINIMIZE

CALL elsunc(x,n,mdc,m,bndfun,bnd,bl,bu,p,w,EXIT,f,c)

!     WRITE FINAL RESULTS

CALL PRINT(nout,EXIT,x,n,f,m,c,mdc,p,w)
WRITE(nout,1070) itotal
1070 FORMAT(t36, 'NUMBER of restarts', i5)
STOP

CONTAINS


!PRINT

SUBROUTINE PRINT(iut,EXIT,x,n,f,m,c,mdc,p,w)

INTEGER, INTENT(IN)    :: iut
INTEGER, INTENT(IN)    :: EXIT
REAL (dp), INTENT(IN)  :: x(:)
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: f(:)
INTEGER, INTENT(IN)    :: m
REAL (dp), INTENT(IN)  :: c(:,:)
INTEGER, INTENT(IN)    :: mdc
INTEGER, INTENT(IN)    :: p(:)
REAL (dp), INTENT(IN)  :: w(:)

! Internal variables

REAL (dp)  :: fsum, fl2nrm
INTEGER    :: np11

!     PRINT FINAL RESULTS ON UNIT IUT

WRITE(iut,1000) EXIT,p(6)

1000 FORMAT(/' EXIT =', i8, '    no. of iterations =', i5)
WRITE(iut,1010) p(11),n,m,w(6)
1010 FORMAT(' PRANK =', i4, '   n =', i4, '   m =', i4, '  speed =', g12.3)
WRITE(iut,1020) p(7),p(8),p(9),p(10)
1020 FORMAT(' FUNCEV =', i5, '  jacev =', i5, '  secev =', i5, '  linev =', i5)
fsum=2.0D0*w(5)
fl2nrm=SQRT(fsum)
WRITE(iut,1030) w(5),fsum,fl2nrm
1030 FORMAT(' PHI(x) =', g15.6, '  fsum =', g15.6, '  fl2nrm =', g15.6)
WRITE(iut,1040) x(1:n)
1040 FORMAT(' AT the point'/ (' ', 6g13.5))
WRITE(iut,1090)
1090 FORMAT(' ACTIVE set indicators - NOT AVAILABLE')
IF(EXIT < 0) GO TO 20
WRITE(iut,1050)
1050 FORMAT(' ESSENTIAL part of covariance matrix')
DO  i=1,n
  WRITE(iut,1060) c(i,1:n)
END DO
1060 FORMAT(' ', 5g12.4)

!     PRINT DEFAULT VALUES WHICH ARE RETURNED FROM ELSUNC

20 WRITE(iut,1070)
1070 FORMAT(/t11, 'DEFAULT values')
WRITE(iut,1080) p(1),p(2),p(3),w(1),w(2),w(3),w(4)
1080 FORMAT(t11, 'IPRINT =', i5/ t11, '  nout =', i5/  &
            t11, ' maxit =', i5/ t11, '   tol =', g15.6/  &
            t11, 'EPSREL =', g15.6/  &
            t11, 'EPSABS =', g15.6/ t11, '  epsx =', g15.6)

RETURN
END SUBROUTINE PRINT

END PROGRAM ex25


!BNDFUN

SUBROUTINE bndfun(x,n,f,m,ctrl,c,mdc)
USE common_ex25
IMPLICIT NONE

REAL (dp), INTENT(IN)    :: x(:)
INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(OUT)   :: f(:)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: ctrl
REAL (dp), INTENT(OUT)   :: c(:,:)
INTEGER, INTENT(IN)      :: mdc

!   THIS IS AN EXAMPLE OF A USER WRITTEN SUBROUTINE THAT SHOULD BE USED
!   IN THE BOUNDARY CONSTRAINED NONLINEAR LEAST SQUARES SOLVER ELSUNC.
!   (WRITTEN BY PER LINDSTROM, UNIVERSITY OF UMEA, SWEDEN)

!   THE ACTUAL TEST EXAMPLE IS
!   EX25 BY HOCK AND SCHITTKOWSKI IN LECTURE NOTES IN ECONOMICS AND
!   MATH. SYSTEMS, V. 187, "Test examples for nonlinear programming codes"

! COMMON /DATA/ u(99)

REAL (dp), PARAMETER  :: argbnd = 88.0_dp
INTEGER               :: i, lctrl
REAL (dp)             :: e1, term

IF(ABS(ctrl) /= 1) GO TO 50

!     EVALUATE THE RESIDUALS

DO  i=1,m
  IF((u(i)-x(2)) <= 0.0) GO TO 20
  term=-(u(i)-x(2))**x(3)/x(1)
  IF(term < -argbnd) e1=0.0
  IF(term > argbnd) GO TO 20
  IF(term >= -argbnd) e1=EXP(term)
  f(i)=-REAL(i)*0.01 + e1
END DO
RETURN

!     IF CTRL=1 ON ENTRY AND UNCOMPUTABILITY SEEMS LIKELY THE USER
!     SHOULD SIGNAL BY SETTING CTRL=-1 ON RETURN.
!     IF CTRL=-1 ON ENTRY AND UNCOMPUTABILITY SEEMS LIKELY THE USER
!     SHOULD SIGNAL BY SETTING CTRL TO A VALUE <-10 WHICH WILL
!     TERMINATE THE ITERATION.

20 lctrl=ctrl
IF(lctrl == 1) ctrl=-1
IF(-1 == lctrl) ctrl=-20
RETURN

!     THE JACOBIAN IS REQUESTED.  IT IS COMPUTED
!     NUMERICALLY BY SETTING CTRL=0

50 DO  i=1,m
  term=(u(i)-x(2))**x(3)/x(1)
  IF(-term < -argbnd) e1=0.0
  IF(-term >= -argbnd) e1=EXP(-term)
  c(i,1)=e1*term/x(1)
  c(i,2)=(u(i)-x(2))**(x(3)-1.0)*x(3)/x(1)*e1
  c(i,3)=-e1*term*LOG(u(i) - x(2))
END DO

RETURN
END SUBROUTINE bndfun

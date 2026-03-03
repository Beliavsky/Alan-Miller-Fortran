!   DRIVER PROGRAM FOR TESTING SOFTWARE INVLTF  IMPLEMENTING
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-12  Time: 17:34:52

!   FOURIER BASED METHOD FOR THE NUMERICAL INVERSION

!                  OF LAPLACE TRANSFORMS


!                     JANUARY , 1998


!   AUTHORS: D'AMORE LUISA, GIULIANO LACCETTI, ALMERICO MURLI



!   REFERENCES
!   ==========

!   D'AMORE L., LACCETTI G., MURLI A.,  -    "ALGORITHM XXX: A FORTRAN
!      SOFTWARE PACKAGE FOR THE NUMERICAL INVERSION OF THE
!      LAPLACE TRANSFORM BASE ON FOURIER SERIES' METHOD"


PROGRAM main

!   ***********************************************************
!   DRIVER  PROGRAM TO TEST  THE ROUTINE INVLTF
!   FOR THE INVERSION OF A LAPLACE TRANSFORM FUNCTION.
!   THIS VERSION USES BOTH REAL AND COMPLEX DOUBLE PRECISION OPERATIONS.
!   THE MAIN PROGRAM ALLOWS THE INVERSION OF A SET OF LAPLACE TRANSFORM
!   FUNCTIONS WHICH CAN BE ADDRESSED  BY A NATURAL NUMBER  BETWEEN 1 AND 34.
!   FOR THE COMPLETE LIST SEE THE TABLE IN THE COMPANION PAPER.

!   THIS IS A SELF-CONTAINED DRIVER FOR THE INVLTF ROUTINE
!   COMPRISING A MAIN PROGRAM AND SUBPROGRAMS QDACC, BACKCF,
!   AND THE LAPLACE TRANSFORM PAIRS FZ AND FEX

!   THIS DRIVER ALLOWS TO OBTAIN UP TO 50 VALUES OF THE
!   INVERSE LAPLACE FUNCTION

!   ************************************************************

USE invert_Laplace
USE common_nf
IMPLICIT NONE

!     .. Parameters ..
INTEGER, PARAMETER  :: nmax = 550, fmax = 50, lastfn = 34
!     ..
!     .. Local Scalars ..
REAL (dp)           :: exf, sigma0, tol, ssbar
INTEGER             :: i, nt
CHARACTER (LEN=77)  :: aa, bb
LOGICAL             :: STOP
!     ..
!     .. Local Arrays ..
REAL (dp)    :: difabs(fmax), difrel(fmax), fzinv(fmax),  &
                tarray(fmax), error(3,fmax), ERR(3)
INTEGER      :: ifail(fmax), ifzeval(fmax)
!     ..
!     .. External Functions ..
! COMPLEX (dp) fz
! REAL (dp) :: fex
! EXTERNAL fz, fex

INTERFACE
  FUNCTION fz(z) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER        :: dp = SELECTED_REAL_KIND(12, 60)
    COMPLEX (dp), INTENT(IN)  :: z
    COMPLEX (dp)              :: fn_val
  END FUNCTION fz

  FUNCTION fex(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION fex
END INTERFACE

!     ..
!     .. External Subroutines ..
! EXTERNAL invltf
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC ABS
!     ..
!     .. Common blocks ..
! COMMON /nf/nfun
!     ..
!     *****************************************************************
!                   SET UP THE OUTPUT
!     *****************************************************************


! Loop through all the tests

DO  nfun = 1,lastfn
  
  
! Skip the tests that require the use of NAG library routines
  
! Commented out lines need to be reintroduced in the routines
! FZ and FEX after labels 9, 36, 37, 38
  
  IF (nfun == 9 .OR. nfun == 36 .OR. nfun == 37 .OR. nfun == 38) CYCLE
  aa = '*****************************************************************************'
  bb = '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
  WRITE (*,FMT=9000) aa, aa, aa, aa
  WRITE (*,FMT=9090)
!     READ (*,FMT=*) NFUN
  WRITE (*,FMT=9050) nfun
  WRITE (*,FMT=9110)
!     READ (*,FMT=*) SIGMA0
  sigma0 = 0.0
  WRITE (*,FMT=9120) sigma0
  IF (nfun == 18) sigma0 = 3.0
  IF (nfun == 23) sigma0 = 0.25D0
  IF (nfun == 29) sigma0 = 2.d0
  WRITE (*,FMT=9080) bb, bb, bb, bb

!     *****************************************************************
!          NOW, SET INPUT PARAMETERS FOR INVLTF
!     *****************************************************************
  
!      NT IS THE NUMBER OF T-VALUES USED HERE FOR THE TEST  FUNCTIONS.
!      THE MAXIMUM ALLOWED IS THE DIMENSION OF THE ARRAY TVAL
!     *****************************************************************
  
!     IN THIS TEST PROGRAM THE INVERSE FUNCTION IS REQUESTED IN THE
!     FOLLOWING VALUES OF T:
!     T=1,20, STEP=0.5 AND T=20,65 STEP=5
!     REMARK:
!     ======
!     FOR T SMALL, THAT IS FOR  T=1,20, STEP=0.5, ALL THE RESULTS ARE
!     QUITE ACCURATE, WHILE FOR T LARGE, THAT IS FOR T=20,65 STEP=5,
!     IN SOME CASES THE RESULTS COULDN'T BE ACCURATE. IN SUCH CASES
!     THE AUTHORS SUGGEST TO USE AN ASYMPTOTIC INVERSION METHOD.
!     *****************************************************************
  
  nt = 44
  DO  i = 1,39
    tarray(i) = 1.d0 + (i-1)*0.5D0
  END DO
  tarray(40) = 30.
  DO  i = 41,44
    tarray(i) = 30.d0 + 5.d0*(i-40)
  END DO
  tol = .1D-5
  ssbar = 0.d0
  WRITE (*,FMT=9010) aa, aa
  WRITE (*,FMT=9030) aa, aa
  WRITE (*,FMT=9020) tol
  
! ************************************************************
!                    CALL OF THE ROUTINE INVLTF
! ************************************************************
  
  DO  i=1,nt
    CALL invltf(tol, tarray(i), fz, sigma0, ssbar, nmax, fzinv(i), ERR,  &
                ifzeval(i), ifail(i))
    error(1,i) = ERR(1)
    error(2,i) = ERR(2)
    error(3,i) = ERR(3)
  END DO
  
  
  WRITE (*,FMT=9040) aa, aa
  
  DO  i = 1,nt
    STOP = .false.
    IF (ifail(i) /= 0 ) THEN
      WRITE (*,FMT=9100) i, ifail(i)
      STOP =.true.
    END IF
    
    
    IF (.NOT. STOP) THEN
      
      
      exf = fex(tarray(i))
      difabs(i) = ABS(exf - fzinv(i))
      IF (exf /= 0.d0) THEN
        difrel(i) = difabs(i)/ABS(exf)
        
      ELSE
        difrel(i) = difabs(i)
      END IF
      
      exf = fex(tarray(i))
      IF (ifail(i) /= -2) THEN
        WRITE (*,FMT=9070) tarray(i), fzinv(i), exf, error(1,i), difrel(i), &
                           error(3,i), error(2,i), difabs(i), ifzeval(i),  &
                           ifail(i)
        
      ELSE
        WRITE (*,FMT=9060) tarray(i), ifail(i)
      END IF
    END IF
    
  END DO
END DO

STOP

9000 FORMAT (' ', a45, a45//  &
             t16, '           SUBROUTINE INVLTF             '/  &
             t16, 'NUMERICAL INVERSION OF A LAPLACE TRANSFORM:'/  &
             t16, ' THIS VERSION USES BOTH REAL AND COMPLEX'/  &
             t16, ' REAL (dp) OPERATIONS'// ' ', a45, a45/)
9010 FORMAT (/a45, a45// t17, '        OUTPUT'/)
9020 FORMAT (/' TOLL --> ', e15.7)
9030 FORMAT (' T      : POINT AT WHICH THE INVERSE TRANSFORM IS',  &
             ' COMPUTED;'//  &
             ' FEX    : EXACT VALUE OF THE INVERSE TRANSFORM;'//   &
             ' FCAL   : COMPUTED VALUE OF THE INVERSE TRANSFORM;'//  &
             ' ESTREL : ESTIMATED RELATIVE ERROR;'//   &
             ' RELERR : ACTUAL RELATIVE ERROR;'//   &
             ' ESTABS : ESTIMATED ABSOLUTE ERROR ;'//   &
             ' ABSERR : ACTUAL ABSOLUTE ERROR;'//   &
             ' N      : # OF FUNCTION EVALUATIONS;'//   &
             ' IFAIL  : = 0 NO INPUT ERRORS; SUCCESSFUL RUN',  &
             '    (ACCURACY REACHED AND IFZEVAL<NMAX),'/   &
             t11, '= 1 TOL >= 1,'/   &
             t11, '= 2 VALT LESS THAN ZERO,'/  &
             t11, '= -1 ACCURACY NOT REACHED AND IFZEVAL > NMAX;,'/  &
             t11, '= -2 THE CHOICE FOR SSBAR  MAY BE NOT OPTIMAL;'/  &
             t11, '        IN SUCH A CASE THE USER MAY SLIGHTLY'/  &
             t11, '        CHANGE THE DEFAULT VALUE;'//a45, a45//)
9040 FORMAT (/ a86/  &
             '   T         FCAL         FEX         ESTREL',  &
             '   RELERR  TRUNERR  ESTABS   ABSERR   N  IFAIL'/ a100//)
9050 FORMAT (/ ' TEST FUNCTION ----->  ', i2/)
9060 FORMAT (f5.1, t96, i2)
9070 FORMAT (f5.1, ' ', e14.8, ' ', e14.8, ' ', e8.3, ' ', e8.3, ' ', e8.3, &
             ' ', e8.3, ' ', e8.3,' ', i3, '   ', i2)
9080 FORMAT (//' ',a45,a45//   &
             '     THE T-VALUES AT WHICH THE INVERSE IS REQUIRED ARE T=1,20', &
             ' STEP=0.5 AND T=20,100 STEP=10.'// ' ', a45, a45)
9090 FORMAT (' TEST  FUNCTION : ')
9100 FORMAT (/ ' ERROR DETECTED , I = ', i3, ' IFAIL=', i3)
9110 FORMAT (/ ' ABSCISSA OF CONVERGENCE  ---> '/)
9120 FORMAT (/ ' ABSCISSA OF CONVERGENCE  ---> ', f5.1/)

END PROGRAM main

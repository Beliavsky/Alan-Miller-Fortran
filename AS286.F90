MODULE Errors_in_Variables

! Reference:
! Park M. Reilly, S.E. Keeler, and H.V. Reilly (1993),
! Algorithm AS 286, Parameter Estimation in the Error-in-Variables Model
! Appl. Statist, vol. 42, 693-709.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-11  Time: 12:08:30


IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! Interfaces to the user routines BF & ZED

INTERFACE
  SUBROUTINE bf(b, f, theta, xi)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(OUT)  :: b(:,:)
    REAL (dp), INTENT(OUT)  :: f(:)
    REAL (dp), INTENT(IN)   :: theta(:), xi(:)
  END SUBROUTINE bf

  SUBROUTINE zed(b, f, theta, xi, z)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: b(:,:)
    REAL (dp), INTENT(IN)   :: f(:)
    REAL (dp), INTENT(IN)   :: theta(:), xi(:)
    REAL (dp), INTENT(OUT)  :: z(:,:)
  END SUBROUTINE zed
END INTERFACE

CONTAINS


SUBROUTINE evm (b, DATA, est, g1, g3, ifault, ndat, neq, npar, nvar, phi, q1, &
                resid, theta, v)

! N.B. Arguments BV, F, G2, H, Q2, S, T, XI & Z
!      have been removed.   They are all workspaces.

!     ALGORITHM AS 286  APPL. STATIST. (1993) VOL.42, NO.4

!  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL

REAL (dp), INTENT(IN OUT)  :: b(:,:)      ! b(neq,nvar)
REAL (dp), INTENT(IN OUT)  :: DATA(:,:)   ! DATA(ndat,nvar)
REAL (dp), INTENT(IN OUT)  :: est(:,:)    ! est(ndat,nvar)
REAL (dp), INTENT(OUT)     :: g1(:,:)     ! g1(npar,npar)
REAL (dp), INTENT(OUT)     :: g3(:,:)     ! g3(npar,npar)
INTEGER, INTENT(OUT)       :: ifault(2)
INTEGER, INTENT(IN)        :: ndat
INTEGER, INTENT(IN)        :: neq
INTEGER, INTENT(IN)        :: npar
REAL (dp), INTENT(OUT)     :: phi
REAL (dp), INTENT(OUT)     :: q1(:)       ! q1(npar)
INTEGER, INTENT(IN)        :: nvar
REAL (dp), INTENT(IN OUT)  :: resid(:,:)  ! resid(ndat,nvar)
REAL (dp), INTENT(IN OUT)  :: theta(:)
REAL (dp), INTENT(IN)      :: v(:,:)      ! v(nvar,nvar)

! Automatic arrays providing workspaces

REAL (dp) :: g2(npar,npar), q2(npar,1)

! Local variables

REAL (dp) :: wdiff, w1, w2
INTEGER   :: ipar, ipar1, ipar2, iter

REAL (dp), PARAMETER :: crit1=1.e-6, zero = 0.0_dp, one = 1.0_dp
INTEGER, PARAMETER   :: niter1=50

ifault(1)=0
ifault(2)=0
DO  iter=1,niter1
  
!  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
  
  CALL inner (b, DATA, est, g1, ifault, ndat, neq, npar,  &
              nvar, phi, q1, resid, theta, v)
  
!  PUT Q1 INTO Q2 AND G1 INTO G2
  
  DO  ipar1=1,npar
    q2(ipar1,1)=q1(ipar1)
    DO  ipar2=ipar1,npar
      g2(ipar1,ipar2)=g1(ipar1,ipar2)
    END DO
  END DO
  
!  CALCULATE CORRECTION FOR THETA
  
  CALL chol (g2, q2(:,1), npar)
  CALL bksub (g2, q2, 1, npar, 1)
  
! CHECK FOR CONVERGENCE
  
  wdiff=zero
  DO  ipar=1,npar
    w1=q2(ipar,1)
    w2=theta(ipar)
    wdiff=MAX(wdiff,ABS(w1/w2))
    theta(ipar)=w2 - w1
  END DO
  IF (wdiff <= crit1) GO TO 50
END DO
ifault(1)=1

!  FILL LOWER TRIANGLE OF G1 AND PREPARE G2
!  FOR THE INVERSION OF G1

50   DO  ipar1=1,npar
  DO  ipar2=ipar1,npar
    w1=g1(ipar1,ipar2)
    g1(ipar2,ipar1)=w1
    g2(ipar1,ipar2)=w1
    g3(ipar1,ipar2)=zero
  END DO
END DO

!  INVERT G1 INTO G3

CALL chol (g2, q2(:,1), npar)
DO  ipar1=1,npar
  g3(ipar1,ipar1)=one/g2(ipar1,ipar1)
END DO
CALL bksub (g2, g3, npar, npar, 3)

RETURN
END SUBROUTINE evm



SUBROUTINE inner (b, DATA, est, g1, ifault, ndat, neq,  &
                  npar, nvar, phi, q1, resid, theta, v)

! N.B. Arguments BV, F, H, S, T, XI & Z have been removed.

REAL (dp), INTENT(OUT)     :: b(:,:)     ! b(neq,nvar)
REAL (dp), INTENT(IN)      :: DATA(:,:)  ! DATA(ndat,nvar)
REAL (dp), INTENT(OUT)     :: est(:,:)   ! est(ndat,nvar)
REAL (dp), INTENT(OUT)     :: g1(:,:)    ! g1(npar,npar)
INTEGER, INTENT(OUT)       :: ifault(2)
INTEGER, INTENT(IN)        :: ndat
INTEGER, INTENT(IN)        :: neq
INTEGER, INTENT(IN)        :: npar
INTEGER, INTENT(IN)        :: nvar
REAL (dp), INTENT(OUT)     :: phi
REAL (dp), INTENT(OUT)     :: q1(:)
REAL (dp), INTENT(OUT)     :: resid(:,:) ! resid(ndat,nvar)
REAL (dp), INTENT(IN OUT)  :: theta(:)
REAL (dp), INTENT(IN)      :: v(:,:)     ! v(nvar,nvar)

! Local variables

REAL (dp) :: bv(neq,nvar), f(neq), h(neq), s(neq,neq), t(neq,1),  &
             xi(nvar), z(neq,npar)
REAL (dp) :: wdiff, w1
INTEGER   :: idat, ieq1, ieq2, ipar1, ipar2, iter, ivar1

REAL (dp), PARAMETER :: crit2=1.e-6, half = 0.5_dp, zero = 0.0_dp
INTEGER, PARAMETER   :: niter2=20

!  INITIALIZE PHI, Q1, AND G1

phi=zero
DO  ipar1=1,npar
  q1(ipar1)=zero
  DO  ipar2=ipar1,npar
    g1(ipar1,ipar2)=zero
  END DO
END DO

!  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT

DO  idat=1,ndat
  xi(1:nvar)=DATA(idat,1:nvar)
  DO  iter=1,niter2
    
!  THIS LOOP PERFORMS A SINGLE INNER ITERATION
    
    wdiff=zero
    CALL bf (b, f, theta, xi)
    
!  CALCULATE BV
    
    DO  ieq1=1,neq
      DO  ivar1=1,nvar
        w1=DOT_PRODUCT( b(ieq1,1:nvar), v(ivar1,1:nvar) )
        bv(ieq1,ivar1)=w1
      END DO
    END DO
    
!  CALCULATE BVB' AND STORE IT IN S
    
    DO  ieq1=1,neq
      DO  ieq2=ieq1,neq
        w1=DOT_PRODUCT( bv(ieq1,1:nvar), b(ieq2,1:nvar) )
        s(ieq1,ieq2)=w1
      END DO
    END DO
    
!  PUT F+B(X-XI) INTO H
    
    DO  ieq1=1,neq
      w1=f(ieq1)
      DO  ivar1=1,nvar
        w1=w1 + b(ieq1,ivar1)*(DATA(idat,ivar1) - xi(ivar1))
      END DO
      h(ieq1)=w1
    END DO
    
!  REPLACE S BY ITS TRIANGULAR FACTORIZATION
!  AND SOLVE ST = H FOR T
    
    CALL chol (s, h, neq)
    DO  ieq1=1,neq
      t(ieq1,1)=h(ieq1)
      DO  ieq2=ieq1,neq
        s(ieq2,ieq1)=s(ieq1,ieq2)
      END DO
    END DO
    CALL bksub (s, t, 1, neq, 1)
    
!  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
    
    DO  ivar1=1,nvar
      w1=DATA(idat,ivar1)
      DO  ieq1=1,neq
        w1=w1 - bv(ieq1,ivar1)*t(ieq1,1)
      END DO
      wdiff=MAX(wdiff, ABS((w1-xi(ivar1))/w1))
      xi(ivar1)=w1
    END DO
    
!  CHECK FOR CONVERGENCE
    
    IF (wdiff <= crit2) GO TO 170
  END DO
  ifault(2)=1
  
!  ACCUMULATE PHI
  
  170  DO  ieq1=1,neq
    phi=phi + h(ieq1)**2
  END DO
  
!  CALCULATE FITTED VALUES AND RESIDUALS
  
  DO  ivar1=1,nvar
    w1=xi(ivar1)
    est(idat,ivar1)=w1
    resid(idat,ivar1)=DATA(idat,ivar1) - w1
  END DO
  
!  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
  
  CALL zed (b, f, theta, xi, z)
  CALL bksub (s, z, npar, neq, 2)
  DO  ipar1=1,npar
    w1=q1(ipar1)
    DO  ieq1=1,neq
      w1=w1 + z(ieq1,ipar1)*h(ieq1)
    END DO
    q1(ipar1)=w1
    DO  ipar2=ipar1,npar
      w1=g1(ipar1,ipar2)
      DO  ieq1=1,neq
        w1=w1 + z(ieq1,ipar1)*z(ieq1,ipar2)
      END DO
      g1(ipar1,ipar2)=w1
    END DO
  END DO
END DO
phi=phi*half

RETURN
END SUBROUTINE inner



SUBROUTINE chol (a, b, n)

!  PERFORM CHOLESKY DECOMPOSITION OF A AND PRELIMINARY TREATMENT OF B

REAL (dp), INTENT(IN OUT)  :: a(:,:)
REAL (dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(IN)        :: n

REAL (dp)  :: w1, w2
INTEGER    :: ic1, ir1, ir2

DO  ir1=1,n
  w1=a(ir1,ir1)
  DO  ir2=1,ir1-1
    w1=w1 - a(ir2,ir1)**2
  END DO
  w1=SQRT(w1)
  a(ir1,ir1)=w1
  DO  ic1=ir1+1,n
    w2=a(ir1,ic1)
    DO  ir2=1,ir1-1
      w2=w2 - a(ir2,ir1)*a(ir2,ic1)
    END DO
    a(ir1,ic1)=w2/w1
  END DO
  w2=b(ir1)
  DO  ir2=1,ir1-1
    w2=w2 - a(ir2,ir1)*b(ir2)
  END DO
  b(ir1)=w2/w1
END DO

RETURN
END SUBROUTINE chol



SUBROUTINE bksub (a, b, m, n, job)

!  PERFORM BACK-SOLUTION ACCORDING AS JOB = 1, 2, OR 3

REAL (dp), INTENT(IN)      :: a(:,:)   ! a(n,n)
REAL (dp), INTENT(IN OUT)  :: b(:,:)   ! b(n,m)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: job

REAL (dp)  :: w1
INTEGER    :: ic1, ir1, ir2, k1, k2, k3

DO  ic1=m,1,-1
  IF (job == 2) THEN
    k1=1
    k2=n
    k3=1
  ELSE
    k2=1
    k3=-1
    IF (job == 1) THEN
      k1=n
    ELSE
      k1=ic1
    END IF
  END IF
  DO  ir1=k1,k2,k3
    w1=b(ir1,ic1)
    IF (job == 2) THEN
      k2=ir1-1
    ELSE
      k1=n
      k2=ir1+1
    END IF
    DO  ir2=k1,k2,k3
      w1=w1 - a(ir1,ir2)*b(ir2,ic1)
    END DO
    b(ir1,ic1)=w1/a(ir1,ir1)
    IF (job == 3) b(ic1,ir1)=b(ir1,ic1)
  END DO
END DO

RETURN
END SUBROUTINE bksub



SUBROUTINE evms (bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, &
                 theta, v)

! N.B. Arguments BVS, G2, PHI, Q1, Q2, XI & ZS have been removed.

!  PARAMETER ESTIMATION FOR THE ERROR-IN-VARIABLES MODEL WITH M = 1
!  91/3/21 PMR

REAL (dp), INTENT(IN OUT)  :: bs(:,:)
REAL (dp), INTENT(IN)      :: DATA(:,:)  ! DATA(ndat,nvar)
REAL (dp), INTENT(IN OUT)  :: est(:,:)   ! est(ndat,nvar)
REAL (dp), INTENT(OUT)     :: g1(:,:)
REAL (dp), INTENT(OUT)     :: g3(:,:)
INTEGER, INTENT(OUT)       :: ifault(2)
INTEGER, INTENT(IN)        :: ndat
INTEGER, INTENT(IN)        :: npar
INTEGER, INTENT(IN)        :: nvar
REAL (dp), INTENT(OUT)     :: resid(:,:) ! resid(ndat,nvar)
REAL (dp), INTENT(OUT)     :: theta(:)
REAL (dp), INTENT(IN)      :: v(:,:)     ! v(nvar,nvar)

! Automatic arrays providing workspaces

REAL (dp)  :: g2(npar,npar), q1(npar), q2(npar,1)

! Local variables

REAL (dp)  :: wdiff, w1, w2
INTEGER    :: ipar, ipar1, ipar2, iter

REAL (dp), PARAMETER :: crit1=1.e-6, one = 1.0_dp, zero = 0.0_dp
INTEGER, PARAMETER   :: niter1=50

ifault(1)=0
ifault(2)=0
DO  iter=1,niter1
  
!  EACH PASS THROUGH THIS LOOP PERFORMS A SINGLE OUTER ITERATION
  
  CALL inners (bs, DATA, est, g1, ifault, ndat, npar, nvar,  &
               q1, resid, theta, v)
  
!  PUT Q1 INTO Q2 AND G1 INTO G2
  
  DO  ipar1=1,npar
    q2(ipar1,1)=q1(ipar1)
    DO  ipar2=ipar1,npar
      g2(ipar1,ipar2)=g1(ipar1,ipar2)
    END DO
  END DO
  
!  CALCULATE CORRECTION FOR THETA
  
  CALL chol (g2, q2(:,1), npar)
  CALL bksub (g2, q2, 1, npar, 1)
  
! CHECK FOR CONVERGENCE
  
  wdiff=zero
  DO  ipar=1,npar
    w1=q2(ipar,1)
    w2=theta(ipar)
    wdiff=MAX(wdiff,ABS(w1/w2))
    theta(ipar)=w2 - w1
  END DO
  IF (wdiff <= crit1) GO TO 50
END DO
ifault(1)=1

!  FILL LOWER TRIANGLE OF G1 AND PREPARE G2 FOR THE INVERSION OF G1

50 DO  ipar1=1,npar
  DO  ipar2=ipar1,npar
    w1=g1(ipar1,ipar2)
    g1(ipar2,ipar1)=w1
    g2(ipar1,ipar2)=w1
    g3(ipar1,ipar2)=zero
  END DO
END DO

!  INVERT G1 INTO G3

CALL chol (g2, q2(:,1), npar)
DO  ipar1=1,npar
  g3(ipar1,ipar1)=one/g2(ipar1,ipar1)
END DO
CALL bksub (g2, g3, npar, npar, 3)

RETURN
END SUBROUTINE evms



SUBROUTINE inners (bs, DATA, est, g1, ifault, ndat, npar,  &
                   nvar, q1, resid, theta, v)

! N.B. Arguments BVS, PHI, XI & ZS have been removed.

!  INNER ITERATION FOR M = 1

REAL (dp), INTENT(OUT)     :: bs(:,:)
REAL (dp), INTENT(IN)      :: DATA(:,:)  ! DATA(ndat,nvar)
REAL (dp), INTENT(OUT)     :: est(:,:)   ! est(ndat,nvar)
REAL (dp), INTENT(OUT)     :: g1(:,:)    ! g1(npar,npar)
INTEGER, INTENT(OUT)       :: ifault(2)
INTEGER, INTENT(IN)        :: ndat
INTEGER, INTENT(IN)        :: npar
INTEGER, INTENT(IN)        :: nvar
REAL (dp), INTENT(OUT)     :: q1(:)
REAL (dp), INTENT(OUT)     :: resid(:,:) ! resid(ndat,nvar)
REAL (dp), INTENT(IN OUT)  :: theta(:)
REAL (dp), INTENT(IN)      :: v(:,:)

! Local variables

REAL (dp)  :: bvs(nvar), phi, xi(nvar), zs(1,npar)
REAL (dp)  :: bvb, fs(1), hs, wdiff, w1
INTEGER    :: idat, ipar1, ipar2, iter, ivar1

REAL (dp), PARAMETER :: crit2=1.e-5, half = 0.5_dp, zero = 0.0_dp
INTEGER, PARAMETER   :: niter2=20

!  INITIALIZE PHI, Q1, AND G1

phi=zero
DO  ipar1=1,npar
  q1(ipar1)=zero
  DO  ipar2=ipar1,npar
    g1(ipar1,ipar2)=zero
  END DO
END DO

!  PERFORM INNER ITERATION TO CONVERGENCE FOR EACH IDAT

DO  idat=1,ndat
  
!  INITIALIZE XI
  
  xi(1:nvar)=DATA(idat,1:nvar)
  DO  iter=1,niter2
    
!  THIS LOOP PERFORMS A SINGLE INNER ITERATION
    
    wdiff=zero
    CALL bf (bs, fs, theta, xi)
    
!  CALCULATE BVS
    
    DO  ivar1=1,nvar
      w1=DOT_PRODUCT( bs(1,1:nvar), v(ivar1,1:nvar) )
      bvs(ivar1)=w1
    END DO
    
!  CALCULATE BVB
    
    bvb=DOT_PRODUCT( bvs(1:nvar), bs(1,1:nvar) )
    
!  PUT F+B(X-XI) INTO HS
    
    hs=fs(1)
    DO  ivar1=1,nvar
      hs=hs + bs(1,ivar1)*(DATA(idat,ivar1) - xi(ivar1))
    END DO
    
!  CALCULATE NEW XI AND PUT RELATIVE CHANGE IN XI INTO WDIFF
    
    DO  ivar1=1,nvar
      w1=DATA(idat,ivar1) - bvs(ivar1)*hs/bvb
      wdiff=MAX(wdiff, ABS((w1-xi(ivar1))/w1))
      xi(ivar1)=w1
    END DO
    
!  CHECK FOR CONVERGENCE
    
    IF (wdiff <= crit2) GO TO 100
  END DO
  ifault(2)=1
  
!  ACCUMULATE PHI
  
  100 phi=phi + hs*hs/bvb
  
!  CALCULATE FITTED VALUES AND RESIDUALS
  
  DO  ivar1=1,nvar
    w1=xi(ivar1)
    est(idat,ivar1)=w1
    resid(idat,ivar1)=DATA(idat,ivar1) - w1
  END DO
  
!  CALCULATE CONTRIBUTIONS TO Q1 AND G1 AND COMPLETE THE LOOPS
  
  CALL zed (bs, fs, theta, xi, zs)
  DO  ipar1=1,npar
    w1=zs(1,ipar1)
    q1(ipar1)=q1(ipar1) + w1*hs/bvb
    DO  ipar2=ipar1,npar
      g1(ipar1,ipar2)=g1(ipar1,ipar2) + w1*zs(1,ipar2)/bvb
    END DO
  END DO
END DO
phi=phi*half

RETURN
END SUBROUTINE inners

END MODULE Errors_in_Variables

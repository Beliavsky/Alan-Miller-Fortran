MODULE common_cnmn1
IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /cnmn1/ delfun, dabfun, fdch, fdchm, ct, ctmin, ctl,  &
!     ctlmin, alphax, abobj1, theta, obj, ndv, ncon, nside, iprint,  &
!     nfdg, nscal, linobj, itmax, itrm, icndir, igoto, nac, info, infog, iter

REAL (dp), SAVE  :: delfun, dabfun, fdch, fdchm, ct, ctmin, ctl,  &
                    ctlmin, alphax, abobj1, theta, obj
INTEGER, SAVE    :: ndv, ncon, nside, iprint, nfdg, nscal, linobj,  &
                    itmax, itrm, icndir, igoto, nac, info, infog, iter
END MODULE common_cnmn1



MODULE common_consav
IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /consav/ dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9, dm10,  &
!     dm11, dm12, dct, dctl, phi, abobj, cta, ctam, ctbm, obj1,  &
!     slope, dx, dx1, fi, xi, dftdf1, alp, fff, a1, a2, a3, a4, f1,  &
!     f2, f3, f4, cv1, cv2, cv3, cv4, app, alpca, alpfes, alpln,  &
!     alpmin, alpnc, alpsav, alpsid, alptot, rspace, idm1, idm2,  &
!     idm3, jdir, iobj, kobj, kcount, ncal(2), nfeas, mscal, ncobj,  &
!     nvc, kount, icount, igood1, igood2, igood3, igood4, ibest,  &
!     iii, nlnc, jgoto, ispace(2)

REAL (dp), SAVE  :: dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9, dm10, dm11, &
                    dm12, dct, dctl, phi, abobj, cta, ctam, ctbm, obj1,  &
                    slope, dx, dx1, fi, xi, dftdf1, alp, fff, a1, a2, a3,  &
                    a4, f1, f2, f3, f4, cv1, cv2, cv3, cv4, app, alpca,   &
                    alpfes, alpln, alpmin, alpnc, alpsav, alpsid, alptot, &
                    rspace
INTEGER, SAVE    :: idm1, idm2, idm3, jdir, iobj, kobj, kcount, ncal(2),  &
                    nfeas, mscal, ncobj, nvc, kount, icount, igood1, igood2, &
                    igood3, igood4, ibest, iii, nlnc, jgoto, ispace(2)
END MODULE common_consav



MODULE Constrained_minimization
IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

!----- CONMIN double precision version
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-14  Time: 15:24:16


CONTAINS


SUBROUTINE conmin(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)
!
!  ROUTINE TO SOLVE CONSTRAINED OR UNCONSTRAINED FUNCTION MINIMIZATION.
!  BY G. N. VANDERPLAATS                          APRIL, 1972.
!  * * * * * * * * * * *   JUNE, 1979 VERSION   * * * * * * * * * * *
!  NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!  REFERENCE;  CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION
!      MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS,
!      NASA TM X-62,282, AUGUST, 1973.
!  STORAGE REQUIREMENTS:
!      PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER)
!      ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON, WHERE N3 = NDV+2.
!  RE-SCALE VARIABLES IF REQUIRED.

USE common_cnmn1
USE common_consav

INTEGER, INTENT(IN)        :: n1, n2, n3, n4, n5
REAL (dp), INTENT(IN OUT)  :: x(n1), vlb(n1), vub(n1), g(n2), scal(n1),  &
                              df(n1), a(n1,n3), s(n1), g1(n2), g2(n2),  &
                              b(n3,n3), c(n4)
INTEGER, INTENT(IN OUT)    :: isc(n2), ic(n3), ms1(n5)

! Local variables
REAL (dp)  :: alp1, alp11, alp12, c1, ct1, ctc, ff1, gi, objb, objd,  &
              scj, si, sib, x1, x12, xid, xx
INTEGER    :: i, ii, j, k, m1, m2, m3, mcn1, nci, ndv1, ndv2, nfeasct, nic, nnac

IF (nscal /= 0 .AND. igoto /= 0) THEN
  DO  i = 1, ndv
    x(i) = c(i)
  END DO
END IF
!     CONSTANTS.
ndv1 = ndv + 1
ndv2 = ndv + 2
IF (igoto /= 0) THEN
!     ------------------------------------------------------------------
!                     CHECK FOR UNBOUNDED SOLUTION
!     ------------------------------------------------------------------
!     STOP IF OBJ IS LESS THAN -1.0D+40
  IF (obj <= -1.0D+40) THEN
    WRITE (6,5100)
    GO TO 520
  END IF
  SELECT CASE ( igoto )
    CASE (    1)
      GO TO 60
    CASE (    2)
      GO TO 210
    CASE (    3)
      GO TO 200
    CASE (    4)
      GO TO 410
    CASE (    5)
      GO TO 430
  END SELECT
END IF
!     ------------------------------------------------------------------
!                      SAVE INPUT CONTROL PARAMETERS
!     ------------------------------------------------------------------
IF (iprint > 0) WRITE (6,7500)
IF (.NOT.(linobj == 0 .OR. (ncon > 0 .OR. nside > 0))) THEN
!     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
!     SOLUTION IS UNBOUNDED.
  WRITE (6,5000) linobj, ncon, nside
  RETURN
END IF
idm1 = itrm
idm2 = itmax
idm3 = icndir
dm1 = delfun
dm2 = dabfun
dm3 = ct
dm4 = ctmin
dm5 = ctl
dm6 = ctlmin
dm7 = theta
dm8 = phi
dm9 = fdch
dm10 = fdchm
dm11 = abobj1
dm12 = alphax
!     ------------------------------------------------------------------
!                                DEFAULTS
!     ------------------------------------------------------------------
IF (itrm <= 0) itrm = 3
IF (itmax <= 0) itmax = 20
ndv1 = ndv + 1
IF (icndir == 0) icndir = ndv1
IF (delfun <= 0.) delfun = .0001
ct = -ABS(ct)
IF (ct >= 0.) ct = -.1
ctmin = ABS(ctmin)
IF (ctmin <= 0.) ctmin = .004
ctl = -ABS(ctl)
IF (ctl >= 0.) ctl = -0.01
ctlmin = ABS(ctlmin)
IF (ctlmin <= 0.) ctlmin = .001
IF (theta <= 0.) theta = 1.
IF (abobj1 <= 0.) abobj1 = .1
IF (alphax <= 0.) alphax = .1
IF (fdch <= 0.) fdch = .01
IF (fdchm <= 0.) fdchm = .01
!     ------------------------------------------------------------------
!                     INITIALIZE INTERNAL PARAMETERS
!     ------------------------------------------------------------------
infog = 0
iter = 0
jdir = 0
iobj = 0
kobj = 0
ndv2 = ndv + 2
kcount = 0
ncal(1) = 0
ncal(2) = 0
nac = 0
nfeas = 0
mscal = nscal
ct1 = itrm
ct1 = 1. / ct1
dct = (ctmin/ABS(ct)) ** ct1
dctl = (ctlmin/ABS(ctl)) ** ct1
phi = 5.
abobj = abobj1
ncobj = 0
ctam = ABS(ctmin)
ctbm = ABS(ctlmin)
!     CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC.
nlnc = 0
IF (ncon /= 0) THEN
  DO  i = 1, ncon
    IF (isc(i) > 0) nlnc = nlnc + 1
  END DO
END IF
!     ------------------------------------------------------------------
!          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
!     ------------------------------------------------------------------
IF (nside /= 0) THEN
  DO  i = 1, ndv
    IF (vlb(i) > vub(i)) THEN
      xx = .5 * (vlb(i)+vub(i))
      x(i) = xx
      vlb(i) = xx
      vub(i) = xx
      WRITE (6,6500) i
    END IF
    xx = x(i) - vlb(i)
    IF (xx < 0.) THEN
!     LOWER BOUND VIOLATED.
      WRITE (6,6600) x(i), vlb(i), i
      x(i) = vlb(i)
    ELSE
      xx = vub(i) - x(i)
      IF (xx < 0.) THEN
        WRITE (6,6700) x(i), vub(i), i
        x(i) = vub(i)
      END IF
    END IF
  END DO
END IF
!     ------------------------------------------------------------------
!                        INITIALIZE SCALING VECTOR, SCAL
!     ------------------------------------------------------------------
IF (nscal /= 0) THEN
  IF (nscal >= 0) THEN
    scal(1:ndv) = 1.
  ELSE
    DO  i = 1, ndv
      si = ABS(scal(i))
      IF (si < 1.0D-20) si = 1.0D-5
      scal(i) = si
      si = 1. / si
      x(i) = x(i) * si
      IF (nside /= 0) THEN
        vlb(i) = vlb(i) * si
        vub(i) = vub(i) * si
      END IF
    END DO
  END IF
END IF
!     ------------------------------------------------------------------
!     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
!     ------------------------------------------------------------------
info = 1
ncal(1) = 1
igoto = 1
GO TO 580

60 obj1 = obj
IF (dabfun <= 0.) dabfun = .001 * ABS(obj)
IF (dabfun < 1.0D-10) dabfun = 1.0D-10
IF (iprint > 0) THEN
!     ------------------------------------------------------------------
!                    PRINT INITIAL DESIGN INFORMATION
!     ------------------------------------------------------------------
  IF (iprint > 1) THEN
    IF (nside == 0 .AND. ncon == 0) WRITE (6,8200)
    IF (nside /= 0 .OR. ncon > 0) WRITE (6,7600)
    WRITE (6,7700) iprint, ndv, itmax, ncon, nside, icndir, nscal,  &
                   nfdg, linobj, itrm, n1, n2, n3, n4, n5
    WRITE (6,7900) ct, ctmin, ctl, ctlmin, theta, phi, delfun, dabfun
    WRITE (6,7800) fdch, fdchm, alphax, abobj1
    IF (nside /= 0) THEN
      WRITE (6,8000)
      DO  i = 1, ndv, 6
        m1 = MIN(ndv,i+5)
        WRITE (6,5400) i, vlb(i:m1)
      END DO
      WRITE (6,8100)
      DO  i = 1, ndv, 6
        m1 = MIN(ndv,i+5)
        WRITE (6,5400) i, vub(i:m1)
      END DO
    END IF
    IF (nscal < 0) THEN
      WRITE (6,8300)
      WRITE (6,9900) scal(1:ndv)
    END IF
    IF (ncon /= 0) THEN
      IF (nlnc /= 0 .AND. nlnc /= ncon) THEN
        WRITE (6,5500)
        DO  i = 1, ncon, 15
          m1 = MIN(ncon,i+14)
          WRITE (6,5600) i, isc(i:m1)
        END DO
      ELSE
        IF (nlnc == ncon) WRITE (6,5700)
        IF (nlnc == 0) WRITE (6,5800)
      END IF
    END IF
  END IF
  WRITE (6,9700) obj
  WRITE (6,9800)
  DO  i = 1, ndv
    x1 = 1.
    IF (nscal /= 0) x1 = scal(i)
    g1(i) = x(i) * x1
  END DO
  DO  i = 1, ndv, 6
    m1 = MIN(ndv,i+5)
    WRITE (6,5400) i, g1(i:m1)
  END DO
  IF (ncon /= 0) THEN
    WRITE (6,10000)
    DO  i = 1, ncon, 6
      m1 = MIN(ncon,i+5)
      WRITE (6,5400) i, g(i:m1)
    END DO
  END IF
END IF
IF (iprint > 1) WRITE (6,8900)
!     ------------------------------------------------------------------
!     ********************  BEGIN MINIMIZATION  ************************
!     ------------------------------------------------------------------
130 iter = iter + 1
IF (abobj1 < .0001) abobj1 = .0001
IF (abobj1 > .2) abobj1 = .2
IF (alphax > 1.) alphax = 1.
IF (alphax < .001) alphax = .001
!
!  THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81
!
!     NFEAS=NFEAS+1
!     IF (NFEAS.GT.10) GO TO 810
IF (iprint > 2) WRITE (6,8400) iter
IF (iprint > 3 .AND. ncon > 0) WRITE (6,8500) ct, ctl, phi
cta = ABS(ct)
IF (ncobj /= 0) THEN
!     ------------------------------------------------------------------
!     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO LONGER ACTIVE.
!     ------------------------------------------------------------------
  nnac = nac
  DO  i = 1, nnac
    IF (ic(i) > ncon) nac = nac - 1
  END DO
  IF (nac <= 0) GO TO 250
  nnac = nac
  DO  i = 1, nnac
    150 nic = ic(i)
    ct1 = ct
    IF (isc(nic) > 0) ct1 = ctl
    IF (g(nic) <= ct1) THEN
      nac = nac - 1
      IF (i > nac) GO TO 250
      DO  k = i, nac
        ii = k + 1
        DO  j = 1, ndv2
          a(j,k) = a(j,ii)
        END DO
        ic(k) = ic(ii)
      END DO
      GO TO 150
    END IF
  END DO
  GO TO 250
END IF
IF (mscal >= nscal .AND. nscal /= 0) THEN
  IF (nscal >= 0 .OR. kcount >= icndir) THEN
    mscal = 0
    kcount = 0
!     ------------------------------------------------------------------
!                          SCALE VARIABLES
!     ------------------------------------------------------------------
    DO  i = 1, ndv
      si = scal(i)
      xi = si * x(i)
      sib = si
      IF (nscal > 0) si = ABS(xi)
      IF (si >= 1.0D-10) THEN
        scal(i) = si
        si = 1. / si
        x(i) = xi * si
        IF (nside /= 0) THEN
          vlb(i) = sib * si * vlb(i)
          vub(i) = sib * si * vub(i)
        END IF
      END IF
    END DO
    IF (.NOT.(iprint < 4 .OR. (nscal < 0 .AND. iter > 1))) THEN
      WRITE (6,8600)
      WRITE (6,9900) scal(1:ndv)
    END IF
  END IF
END IF
mscal = mscal + 1
nac = 0
!     ------------------------------------------------------------------
!          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
info = 2
ncal(2) = ncal(2) + 1
IF (nfdg == 1) THEN
  igoto = 2
  GO TO 580
END IF
jgoto = 0

200 CALL cnmn01(jgoto,x,df,g,isc,ic,a,g1,vub,scal,ncal,dx,dx1,  &
                fi,xi,iii,n1,n2,n3)
igoto = 3
IF (jgoto > 0) GO TO 580

210 info = 1
IF (nac >= n3) GO TO 520
IF (nscal /= 0 .AND. nfdg /= 0) THEN
!     ------------------------------------------------------------------
!                              SCALE GRADIENTS
!     ------------------------------------------------------------------
!     SCALE GRADIENT OF OBJECTIVE FUNCTION.
  df(1:ndv) = df(1:ndv) * scal(1:ndv)
  IF (nfdg /= 2 .AND. nac /= 0) THEN
!     SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
    DO  j = 1, ndv
      scj = scal(j)
      a(j,1:nac) = a(j,1:nac) * scj
    END DO
  END IF
END IF

250 IF (iprint >= 3 .AND. ncon /= 0) THEN
!     ------------------------------------------------------------------
!                                   PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
  m1 = 0
  m2 = n3
  IF (nac /= 0) THEN
    DO  i = 1, nac
      j = ic(i)
      IF (j <= ncon) THEN
        gi = g(j)
        c1 = ctam
        IF (isc(j) > 0) c1 = ctbm
        gi = gi - c1
        IF (gi <= 0.) THEN
!     ACTIVE CONSTRAINT.
          m1 = m1 + 1
          ms1(m1) = j
        ELSE
          m2 = m2 + 1
!     VIOLATED CONSTRAINT.
          ms1(m2) = j
        END IF
      END IF
    END DO
  END IF
  m3 = m2 - n3
  WRITE (6,5900) m1
  IF (m1 /= 0) THEN
    WRITE (6,6000)
    WRITE (6,10100) ms1(1:m1)
  END IF
  WRITE (6,6100) m3
  IF (m3 /= 0) THEN
    WRITE (6,6000)
    m3 = n3 + 1
    WRITE (6,10100) ms1(m3:m2)
  END IF
END IF
!     ------------------------------------------------------------------
!            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
!     ------------------------------------------------------------------
IF (nside /= 0) THEN
  mcn1 = ncon
  m1 = 0
  DO  i = 1, ndv
!     LOWER BOUND.
    xi = x(i)
    xid = vlb(i)
    x12 = ABS(xid)
    IF (x12 < 1.) x12 = 1.
    gi = (xid-xi) / x12
    IF (gi >= -1.0D-6) THEN
      m1 = m1 + 1
      ms1(m1) = -i
      nac = nac + 1
      IF (nac >= n3) GO TO 520
      mcn1 = mcn1 + 1
      DO  j = 1, ndv
        a(j,nac) = 0.
      END DO
      a(i,nac) = -1.
      ic(nac) = mcn1
      g(mcn1) = gi
      isc(mcn1) = 1
    END IF
!     UPPER BOUND.
    xid = vub(i)
    x12 = ABS(xid)
    IF (x12 < 1.) x12 = 1.
    gi = (xi-xid) / x12
    IF (gi >= -1.0D-6) THEN
      m1 = m1 + 1
      ms1(m1) = i
      nac = nac + 1
      IF (nac >= n3) GO TO 520
      mcn1 = mcn1 + 1
      DO  j = 1, ndv
        a(j,nac) = 0.
      END DO
      a(i,nac) = 1.
      ic(nac) = mcn1
      g(mcn1) = gi
      isc(mcn1) = 1
    END IF
  END DO
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
  IF (iprint >= 3) THEN
    WRITE (6,6200) m1
    IF (m1 /= 0) THEN
      WRITE (6,6300)
      WRITE (6,10100) (ms1(j),j = 1,m1)
    END IF
  END IF
END IF
!     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
IF (iprint >= 4) THEN
  WRITE (6,8700)
  DO  i = 1, ndv, 6
    m1 = MIN(ndv,i+5)
    WRITE (6,5400) i, df(i:m1)
  END DO
  IF (nac /= 0) THEN
    WRITE (6,8800)
    DO  i = 1, nac
      m1 = ic(i)
      m2 = m1 - ncon
      m3 = 0
      IF (m2 > 0) m3 = ABS(ms1(m2))
      IF (m2 <= 0) WRITE (6,5200) m1
      IF (m2 > 0) WRITE (6,5300) m3
      DO  k = 1, ndv, 6
        m1 = MIN(ndv,k+5)
        WRITE (6,5400) k, (a(j,i),j = k,m1)
      END DO
      WRITE (6,8900)
    END DO
  END IF
END IF
!     ------------------------------------------------------------------
!     ******************  DETERMINE SEARCH DIRECTION *******************
!     ------------------------------------------------------------------
alp = 1.0D+20
IF (nac > 0) GO TO 340
!     ------------------------------------------------------------------
!                        UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION.
!
!  S. N. 575 ADDED ON 2/25/81
!
330 nvc = 0
nfeas = 0
kcount = kcount + 1
!     IF KCOUNT.GT.ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM.
IF (kcount > icndir .OR. iobj == 2) kcount = 1
IF (kcount == 1) jdir = 0
!     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
CALL cnmn02(jdir,slope,dftdf1,df,s)
GO TO 380
!     ------------------------------------------------------------------
!                          CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND USABLE-FEASIBLE DIRECTION.
340 kcount = 0
jdir = 0
phi = 10. * phi
IF (phi > 1000.) phi = 1000.
!
!  THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81
!
!     IF (NFEAS.EQ.1) PHI=5.
!     CALCULATE DIRECTION, S.
CALL cnmn05(g,df,a,s,b,c,slope,phi,isc,ic,ms1,nvc,n1,n2,n3,n4,n5)
!
!  THE FOLLOWING LINE WAS ADDED ON 2/25/81
!
IF (nac == 0) GO TO 330
!
!  THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81
!  REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING
!
!     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
!     SET ABOBJ1=.05 (5 PERCENT).
!     IF (NVC.EQ.0 .AND. NFEAS.GT.1) ABOBJ1=.05
!     IF (NVC.EQ.0) NFEAS=0
IF (iprint >= 3) THEN
  WRITE (6,9000)
  DO  i = 1, nac, 6
    m1 = MIN(nac,i+5)
    WRITE (6,5400) i, (a(ndv1,j),j = i,m1)
  END DO
  WRITE (6,7400) s(ndv1)
END IF
!     ------------------------------------------------------------------
!     ****************** ONE-DIMENSIONAL SEARCH ************************
!     ------------------------------------------------------------------
IF (s(ndv1) < 1.0D-6 .AND. nvc == 0) GO TO 450
!     ------------------------------------------------------------------
!                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
!     ------------------------------------------------------------------
IF (nvc /= 0) THEN
  alp = -1.
  DO  i = 1, nac
    nci = ic(i)
    c1 = g(nci)
    ctc = ctam
    IF (isc(nci) > 0) ctc = ctbm
    IF (c1 > ctc) THEN
      alp1 = DOT_PRODUCT( s(1:ndv), a(1:ndv,i) )
      alp1 = alp1 * a(ndv2,i)
      IF (ABS(alp1) >= 1.0D-20) THEN
        alp1 = -c1 / alp1
        IF (alp1 > alp) alp = alp1
      END IF
    END IF
  END DO
END IF
!     ------------------------------------------------------------------
!                       LIMIT CHANCE TO ABOBJ1*OBJ
!     ------------------------------------------------------------------
380 alp1 = 1.0D+20
si = ABS(obj)
IF (si < .01) si = .01
IF (ABS(slope) > 1.0D-20) alp1 = abobj1 * si / slope
alp1 = ABS(alp1)
IF (nvc > 0) alp1 = 10. * alp1
IF (alp1 < alp) alp = alp1
!     ------------------------------------------------------------------
!                   LIMIT CHANGE IN VARIABLE TO ALPHAX
!     ------------------------------------------------------------------
alp11 = 1.0D+20
DO  i = 1, ndv
  si = ABS(s(i))
  xi = ABS(x(i))
  IF (si >= 1.0D-10 .AND. xi >= 0.1) THEN
    alp1 = alphax * xi / si
    IF (alp1 < alp11) alp11 = alp1
  END IF
END DO
IF (nvc > 0) alp11 = 10. * alp11
IF (alp11 < alp) alp = alp11
IF (alp > 1.0D+20) alp = 1.0D+20
IF (alp <= 1.0D-20) alp = 1.0D-20
IF (iprint >= 3) THEN
  WRITE (6,9100)
  DO  i = 1, ndv, 6
    m1 = MIN(ndv,i+5)
    WRITE (6,5400) i, (s(j),j = i,m1)
  END DO
  WRITE (6,6400) slope, alp
END IF
IF (ncon > 0 .OR. nside > 0) GO TO 420
!     ------------------------------------------------------------------
!           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
jgoto = 0
410 CALL cnmn03(x,s,slope,alp,fff,a1,a2,a3,a4,f1,f2,f3,f4,app,ncal,kount,jgoto)
igoto = 4
IF (jgoto > 0) GO TO 580
jdir = 1
!     PROCEED TO CONVERGENCE CHECK.
GO TO 450
!     ------------------------------------------------------------------
!       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
420 jgoto = 0

430 CALL cnmn06(x,vlb,vub,g,scal,df,s,g1,g2,ctam,ctbm,slope,alp,a2,a3,a4,  &
                f1,f2,f3,cv1,cv2,cv3,cv4,alpca,alpfes,alpln,alpmin,alpnc,  &
                alpsav,alpsid,alptot,isc,ncal,nvc,icount,igood1,  &
                igood2,igood3,igood4,ibest,iii,nlnc,jgoto)
igoto = 5
IF (jgoto > 0) GO TO 580
IF (nac == 0) jdir = 1
!     ------------------------------------------------------------------
!     *******************     UPDATE ALPHAX   **************************
!     ------------------------------------------------------------------
450 IF (alp > 1.0D+19) alp = 0.
!     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I) AND ALHPAX.
alp11 = 0.
DO  i = 1, ndv
  si = ABS(s(i))
  xi = ABS(x(i))
  IF (xi >= 1.0D-10) THEN
    alp1 = alp * si / xi
    IF (alp1 > alp11) alp11 = alp1
  END IF
END DO
alp11 = .5 * (alp11+alphax)
alp12 = 5. * alphax
IF (alp11 > alp12) alp11 = alp12
alphax = alp11
ncobj = ncobj + 1
!     ABSOLUTE CHANGE IN OBJECTIVE.
objd = obj1 - obj
objb = ABS(objd)
IF (objb < 1.0D-10) objb = 0.
IF (nac == 0 .OR. objb > 0.) ncobj = 0
IF (ncobj > 1) ncobj = 0
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
IF (iprint >= 3) THEN
  WRITE (6,9200) alp
END IF
IF (iprint >= 2) THEN
  IF (objb <= 0.) THEN
    IF (iprint == 2) WRITE (6,9300) iter, obj
    IF (iprint > 2) WRITE (6,9400) obj
  ELSE
    IF (iprint /= 2) THEN
      WRITE (6,9500) obj
    ELSE
      WRITE (6,9600) iter, obj
    END IF
  END IF
  WRITE (6,9800)
  DO  i = 1, ndv
    ff1 = 1.
    IF (nscal /= 0) ff1 = scal(i)
    g1(i) = ff1 * x(i)
  END DO
  DO  i = 1, ndv, 6
    m1 = MIN(ndv,i+5)
    WRITE (6,5400) i, (g1(j),j = i,m1)
  END DO
  IF (ncon /= 0) THEN
    WRITE (6,10000)
    DO  i = 1, ncon, 6
      m1 = MIN(ncon,i+5)
      WRITE (6,5400) i, (g(j),j = i,m1)
    END DO
  END IF
END IF
!
!  THE FOLLOWING CODE WAS ADDED ON 3/5/81
!
!  IT HAD NOT BEEN REPORTED AS A FIX TO MAOB
!  BUT WAS SENT TO JEFF STROUD A YEAR AGO
!  SEE OTHER COMMENTS IN CONMIN SUBROUTINE FOR DELETIONS OF CODE
!  ON 3/5/81 PERTAINING TO THIS FIX
!
!
!                   CHECK FEASIBILITY
!
IF (ncon > 0) THEN
  nfeasct = 10
!  added by slp 11/17/94
  DO  i = 1, ncon
    c1 = ctam
    IF (isc(i) > 0) c1 = ctbm
    IF (g(i) > c1) THEN
      nfeas = nfeas + 1
      GO TO 510
    END IF
  END DO
  IF (nfeas > 0) abobj1 = .05
!cc
  nfeas = 0
  phi = 5.
  510 IF (nfeas >= nfeasct) GO TO 520
END IF
!
!  END OF INSERTED FIX
!
!     ------------------------------------------------------------------
!                          CHECK CONVERGENCE
!     ------------------------------------------------------------------
!     STOP IF ITER EQUALS ITMAX.
IF (iter < itmax) THEN
!     ------------------------------------------------------------------
!                     ABSOLUTE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
  objb = ABS(objd)
  kobj = kobj + 1
  IF (objb >= dabfun .OR. nfeas > 0) kobj = 0
!     ------------------------------------------------------------------
!                     RELATIVE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
  IF (ABS(obj1) > 1.0D-10) objd = objd / ABS(obj1)
  abobj1 = .5 * (ABS(abobj)+ABS(objd))
  abobj = ABS(objd)
  iobj = iobj + 1
  IF (nvc > 0 .OR. objd >= delfun) iobj = 0
  IF (iobj < itrm .AND. kobj < itrm) THEN
    obj1 = obj
!     ------------------------------------------------------------------
!           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
!     ------------------------------------------------------------------
    IF (iobj < 1 .OR. nac == 0) GO TO 130
    ct = dct * ct
    ctl = ctl * dctl
    IF (ABS(ct) < ctmin) ct = -ctmin
    IF (ABS(ctl) < ctlmin) ctl = -ctlmin
    GO TO 130
  END IF
END IF

520 IF (nac >= n3) WRITE (6,10200)
!     ------------------------------------------------------------------
!     ****************  FINAL FUNCTION INFORMATION  ********************
!     ------------------------------------------------------------------
IF (nscal /= 0) THEN
!     UN-SCALE THE DESIGN VARIABLES.
  DO  i = 1, ndv
    xi = scal(i)
    IF (nside /= 0) THEN
      vlb(i) = xi * vlb(i)
      vub(i) = xi * vub(i)
    END IF
    x(i) = xi * x(i)
  END DO
END IF
!     ------------------------------------------------------------------
!                           PRINT FINAL RESULTS
!     ------------------------------------------------------------------
IF (iprint /= 0 .AND. nac < n3) THEN
  WRITE (6,10300)
  WRITE (6,9500) obj
  WRITE (6,9800)
  DO  i = 1, ndv, 6
    m1 = MIN(ndv,i+5)
    WRITE (6,5400) i, (x(j),j = i,m1)
  END DO
  IF (ncon /= 0) THEN
    WRITE (6,10000)
    DO  i = 1, ncon, 6
      m1 = MIN(ncon,i+5)
      WRITE (6,5400) i, (g(j),j = i,m1)
    END DO
!     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
    nac = 0
    nvc = 0
    DO  i = 1, ncon
      cta = ctam
      IF (isc(i) > 0) cta = ctbm
      gi = g(i)
      IF (gi <= cta) THEN
        IF (gi < ct .AND. isc(i) == 0) CYCLE
        IF (gi < ctl .AND. isc(i) > 0) CYCLE
        nac = nac + 1
        ic(nac) = i
      ELSE
        nvc = nvc + 1
        ms1(nvc) = i
      END IF
    END DO
    WRITE (6,5900) nac
    IF (nac /= 0) THEN
      WRITE (6,6000)
      WRITE (6,10100) ic(1:nac)
    END IF
    WRITE (6,6100) nvc
    IF (nvc /= 0) THEN
      WRITE (6,6000)
      WRITE (6,10100) ms1(1:nvc)
    END IF
  END IF
  IF (nside /= 0) THEN
!     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
    nac = 0
    DO  i = 1, ndv
      xi = x(i)
      xid = vlb(i)
      x12 = ABS(xid)
      IF (x12 < 1.) x12 = 1.
      gi = (xid-xi) / x12
      IF (gi >= -1.0D-6) THEN
        nac = nac + 1
        ms1(nac) = -i
      END IF
      xid = vub(i)
      x12 = ABS(xid)
      IF (x12 < 1.) x12 = 1.
      gi = (xi-xid) / x12
      IF (gi >= -1.0D-6) THEN
        nac = nac + 1
        ms1(nac) = i
      END IF
    END DO
    WRITE (6,6200) nac
    IF (nac /= 0) THEN
      WRITE (6,6300)
      WRITE (6,10100) (ms1(j),j = 1,nac)
    END IF
  END IF
  WRITE (6,6800)
  IF (iter >= itmax) WRITE (6,6900)
  IF (nfeas >= nfeasct) WRITE (6,7000)
  IF (iobj >= itrm) WRITE (6,7100) itrm
  IF (kobj >= itrm) WRITE (6,7200) itrm
  WRITE (6,7300) iter
  WRITE (6,10400) ncal(1)
  IF (ncon > 0) WRITE (6,10500) ncal(1)
  IF (nfdg /= 0) WRITE (6,10600) ncal(2)
  IF (ncon > 0 .AND. nfdg == 1) WRITE (6,10700) ncal(2)
END IF
!     ------------------------------------------------------------------
!                   RE-SET BASIC PARAMETERS TO INPUT VALUES
!     ------------------------------------------------------------------
itrm = idm1
itmax = idm2
icndir = idm3
delfun = dm1
dabfun = dm2
ct = dm3
ctmin = dm4
ctl = dm5
ctlmin = dm6
theta = dm7
phi = dm8
fdch = dm9
fdchm = dm10
abobj1 = dm11
alphax = dm12
igoto = 0

580 IF (nscal == 0 .OR. igoto == 0) RETURN
!     UN-SCALE VARIABLES.
DO  i = 1, ndv
  c(i) = x(i)
  x(i) = x(i) * scal(i)
END DO
RETURN
!     ------------------------------------------------------------------
!                                FORMATS
!     ------------------------------------------------------------------
!
!
5000 FORMAT (//t6,  &
    'A COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR OBJECTIVE IS SPECIFIED'// &
    t11, 'LINOBJ =', i5/ t11, 'NCON   =', i5/  &
    t11, 'NSIDE  =',i5// t6, 'CONTROL RETURNED TO CALLING PROGRAM')
5100 FORMAT (//t6,   &
    'CONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -1.0E+40'/  &
    t6, 'SOLUTION APPEARS TO BE UNBOUNDED'/ t6, 'OPTIMIZATION IS TERMINATED')
5200 FORMAT (t6, 'CONSTRAINT NUMBER', i5)
5300 FORMAT (t6, 'SIDE CONSTRAINT ON VARIABLE', i5)
5400 FORMAT (t4, i5, ')  ', 6E13.5)
5500 FORMAT (/t6, 'LINEAR CONSTRAINT IDENTIFIERS (ISC)'/  &
    t6, 'NON-ZERO INDICATES LINEAR CONSTRAINT')
5600 FORMAT (t4, i5, ')  ', 15I5)
5700 FORMAT (/t6, 'ALL CONSTRAINTS ARE LINEAR')
5800 FORMAT (/t6, 'ALL CONSTRAINTS ARE NON-LINEAR')
5900 FORMAT (/t6, 'THERE ARE',i5,' ACTIVE CONSTRAINTS')
6000 FORMAT (t6, 'CONSTRAINT NUMBERS ARE')
6100 FORMAT (/t6, 'THERE ARE', i5, ' VIOLATED CONSTRAINTS')
6200 FORMAT (/t6, 'THERE ARE', i5, ' ACTIVE SIDE CONSTRAINTS')
6300 FORMAT (t6, 'DECISION VARIABLES AT LOWER OR UPPER BOUNDS',  &
    ' (MINUS INDICATES LOWER BOUND)')
6400 FORMAT (/t6, 'ONE-DIMENSIONAL SEARCH'/ t6, 'INITIAL SLOPE =', e12.4,  &
    '  PROPOSED ALPHA =', e12.4)
6500 FORMAT (//t6, '* * CONMIN DETECTS VLB(I).GT.VUB(I)'/  &
    t6, 'FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =', i5)
6600 FORMAT (//t6, '* * CONMIN DETECTS INITIAL X(I).LT.VLB(I)'/ t6,   &
    'X(I) =', e12.4, '  VLB(I) =', e12.4/ t6,   &
    'X(I) IS SET EQUAL TO VLB(I) FOR I =',i5)
6700 FORMAT (//t6, '* * CONMIN DETECTS INITIAL X(I).GT.VUB(I)'/t6,   &
    'X(I) =', e12.4, '  VUB(I) =', e12.4/ t6,   &
    'X(I) IS SET EQUAL TO VUB(I) FOR I =',i5)
6800 FORMAT (/t6, 'TERMINATION CRITERION')
6900 FORMAT (t11, 'ITER EQUALS ITMAX')
7000 FORMAT (t11,  &
    'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A FEASIBLE DESIGN')
7100 FORMAT (t11, 'ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR', i3,  &
    ' ITERATIONS')
7200 FORMAT (t11,'ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR',i3,  &
    ' ITERATIONS')
7300 FORMAT (/t6, 'NUMBER OF ITERATIONS =', i5)
7400 FORMAT (/t6, 'CONSTRAINT PARAMETER, BETA =', e14.5)
7500 FORMAT (// t13, 27('* ')/ t13, '*', t65, '*'/ t13,'*',t34,  &
    'C O N M I N', t65, '*'/ t13,'*', t65, '*'/ t13,'*',t29,  &
    ' FORTRAN PROGRAM FOR ', t65, '*'/ t13,'*',t65,'*'/ t13,'*',t23,  &
    'CONSTRAINED FUNCTION MINIMIZATION', t65, '*'/ t13, '*',t65, '*'/  &
    t13, 27('* '))
7600 FORMAT (//t6, 'CONSTRAINED FUNCTION MINIMIZATION'//   &
    t6, 'CONTROL PARAMETERS')
7700 FORMAT (/t6, 'IPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSCAL NFDG'/ &
    8I8//t6, 'LINOBJ  ITRM     N1      N2      N3      N4      N5'/ 8I8)
7800 FORMAT (/t10, 'FDCH', t26, 'FDCHM', t42, 'ALPHAX', t58, 'ABOBJ1'/  &
    ' ', 4('  ', e14.5))
7900 FORMAT (/t10, 'CT', t26, 'CTMIN', t42, 'CTL', t58, 'CTLMIN'/  &
    ' ', 4('  ', e14.5)//  &
    t10, 'THETA', t26, 'PHI', t42, 'DELFUN', t58, 'DABFUN'/  &
    ' ',4('  ', e14.5))
8000 FORMAT (/t6, 'LOWER BOUNDS ON DECISION VARIABLES (VLB)')
8100 FORMAT (/t6, 'UPPER BOUNDS ON DECISION VARIABLES (VUB)')
8200 FORMAT (//t6, 'UNCONSTRAINED FUNCTION MINIMIZATION'//t6,   &
    'CONTROL PARAMETERS')
8300 FORMAT (/t6, 'SCALING VECTOR (SCAL)')
8400 FORMAT (//t6, 'BEGIN ITERATION NUMBER',i5)
8500 FORMAT (/t6, 'CT =', e14.5, '     CTL =', e14.5, '     PHI =', e14.5)
8600 FORMAT (/t6, 'NEW SCALING VECTOR (SCAL)')
8700 FORMAT (/t6, 'GRADIENT OF OBJ')
8800 FORMAT (/t6, 'GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS')
8900 FORMAT (' ')
9000 FORMAT (/t6, 'PUSH-OFF FACTORS, (THETA(I), I=1,NAC)')
9100 FORMAT (/t6, 'SEARCH DIRECTION (S-VECTOR)')
9200 FORMAT (/t6, 'CALCULATED ALPHA =', e14.5)
9300 FORMAT (//t6, 'ITER =', i5, '     OBJ =', e14.5, '     NO CHANGE IN OBJ')
9400 FORMAT (/t6, 'OBJ =', e15.6, '     NO CHANGE ON OBJ')
9500 FORMAT (/t6, 'OBJ =', e15.6)
9600 FORMAT (//t6, 'ITER =', i5, '     OBJ =',e14.5)
9700 FORMAT (//t6, 'INITIAL FUNCTION INFORMATION'// t6, 'OBJ =', e15.6)
9800 FORMAT (/t6, 'DECISION VARIABLES (X-VECTOR)')
9900 FORMAT (t4, 7E13.4)
10000 FORMAT (/t6, 'CONSTRAINT VALUES (G-VECTOR)')
10100 FORMAT (t6,15I5)
10200 FORMAT (/t6,  'THE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEEDS N3-1.'/  &
    t6, 'DIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC IS INSUFFICIENT'/  &
    t6, 'OPTIMIZATION TERMINATED AND CONTROL RETURNED TO MAIN PROGRAM.')
10300 FORMAT (///'    FINAL OPTIMIZATION INFORMATION')
10400 FORMAT (/t6, 'OBJECTIVE FUNCTION WAS EVALUATED        ', i5, '  TIMES')
10500 FORMAT (/t6, 'CONSTRAINT FUNCTIONS WERE EVALUATED', i10, '  TIMES')
10600 FORMAT (/t6, 'GRADIENT OF OBJECTIVE WAS CALCULATED', i9, '  TIMES')
10700 FORMAT (/t6, 'GRADIENTS OF CONSTRAINTS WERE CALCULATED', i5, '  TIMES')
END SUBROUTINE conmin



!----- CNMN01
SUBROUTINE cnmn01(jgoto,x,df,g,isc,ic,a,g1,vub,scal,ncal,dx,  &
                  dx1,fi,xi,iii,n1,n2,n3)

! N.B. Arguments VLB, C & N4 have been removed.

!  ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE.
!  BY G. N. VANDERPLAATS                         JUNE, 1972.
!  NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.

USE common_cnmn1
INTEGER, INTENT(IN OUT)   :: jgoto
REAL (dp), INTENT(IN)     :: vub(:), scal(:)
INTEGER, INTENT(IN)       :: n1, n2, n3
REAL (dp), INTENT(IN OUT) :: x(:), df(:), g(n2), a(n1,n3), g1(n2)
INTEGER, INTENT(IN OUT)   :: ic(n3), ncal(2)
INTEGER, INTENT(IN)       :: isc(n2)
REAL (dp), INTENT(OUT)    :: dx, dx1, fi, xi
INTEGER, INTENT(OUT)      :: iii

! Local variables
REAL (dp)  :: fdch1, x1
INTEGER    :: i, i1, inf, j

IF (jgoto /= 1) THEN
  IF (jgoto == 2) GO TO 40
  infog = 0
  inf = info
  nac = 0
  IF (linobj == 0 .OR. iter <= 1) THEN
!     ------------------------------------------------------------------
!                    GRADIENT OF LINEAR OBJECTIVE
!     ------------------------------------------------------------------
    IF (nfdg == 2) jgoto = 1
    IF (nfdg == 2) RETURN
  END IF
END IF
jgoto = 0
IF (nfdg == 2 .AND. ncon == 0) RETURN
IF (ncon /= 0) THEN
!     ------------------------------------------------------------------
!       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
!     ------------------------------------------------------------------
  DO  i = 1, ncon
    IF (g(i) >= ct) THEN
      IF (isc(i) <= 0 .OR. g(i) >= ctl) THEN
        nac = nac + 1
        IF (nac >= n3) RETURN
        ic(nac) = i
      END IF
    END IF
  END DO
  IF (nfdg == 2 .AND. nac == 0) RETURN
  IF (linobj > 0 .AND. iter > 1 .AND. nac == 0) RETURN
!     ------------------------------------------------------------------
!                  STORE VALUES OF CONSTRAINTS IN G1
!     ------------------------------------------------------------------
  g1(1:ncon) = g(1:ncon)
END IF
jgoto = 0
IF (nac == 0 .AND. nfdg == 2) RETURN
!     ------------------------------------------------------------------
!                            CALCULATE GRADIENTS
!     ------------------------------------------------------------------
infog = 1
info = 1
fi = obj
iii = 0

30 iii = iii + 1
xi = x(iii)
dx = fdch * xi
dx = ABS(dx)
fdch1 = fdchm
IF (nscal /= 0) fdch1 = fdchm / scal(iii)
IF (dx < fdch1) dx = fdch1
x1 = xi + dx
IF (nside /= 0) THEN
  IF (x1 > vub(iii)) dx = -dx
END IF
dx1 = 1. / dx
x(iii) = xi + dx
ncal(1) = ncal(1) + 1
!     ------------------------------------------------------------------
!                         FUNCTION EVALUATION
!     ------------------------------------------------------------------
jgoto = 2
RETURN

40 x(iii) = xi
IF (nfdg == 0) df(iii) = dx1 * (obj-fi)
IF (nac /= 0) THEN
!     ------------------------------------------------------------------
!             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
  DO  j = 1, nac
    i1 = ic(j)
    a(iii,j) = dx1 * (g(i1)-g1(i1))
  END DO
END IF
IF (iii < ndv) GO TO 30
infog = 0
info = inf
jgoto = 0
obj = fi
IF (ncon == 0) RETURN
!     ------------------------------------------------------------------
!             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
!     ------------------------------------------------------------------
g(1:ncon) = g1(1:ncon)
RETURN
END SUBROUTINE cnmn01



!----- CNMN02

SUBROUTINE cnmn02(ncalc,slope,dftdf1,df,s)

! N.B. Argument N1 has been removed.

!  ROUTINE TO DETERMINE CONJUGATE DIRECTION VECTOR OR DIRECTION
!  OF STEEPEST DESCENT FOR UNCONSTRAINED FUNCTION MINIMIZATION.
!  BY G. N. VANDERPLAATS                       APRIL, 1972.
!  NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!  NCALC = CALCULATION CONTROL.
!      NCALC = 0,     S = STEEPEST DESCENT.
!      NCALC = 1,     S = CONJUGATE DIRECTION.
!  CONJUGATE DIRECTION IS FOUND BY FLETCHER-REEVES ALGORITHM.

USE common_cnmn1
INTEGER, INTENT(IN OUT)    :: ncalc
REAL (dp), INTENT(OUT)     :: slope, dftdf1
REAL (dp), INTENT(IN)      :: df(:)
REAL (dp), INTENT(IN OUT)  :: s(:)

! Local variables
INTEGER    :: i
REAL (dp)  :: beta, dfi, dftdf, s1, s2, si
!     ------------------------------------------------------------------
!                   CALCULATE NORM OF GRADIENT VECTOR
!     ------------------------------------------------------------------
dftdf = SUM( df(1:ndv)**2 )
!     ------------------------------------------------------------------
!     **********                FIND DIRECTION S              **********
!     ------------------------------------------------------------------
IF (ncalc == 1) THEN
  IF (dftdf1 >= 1.0D-20) THEN
!     ------------------------------------------------------------------
!                 FIND FLETCHER-REEVES CONJUGATE DIRECTION
!     ------------------------------------------------------------------
    beta = dftdf / dftdf1
    slope = 0.
    DO  i = 1, ndv
      dfi = df(i)
      si = beta * s(i) - dfi
      slope = slope + si * dfi
      s(i) = si
    END DO
    GO TO 40
  END IF
END IF
ncalc = 0
!     ------------------------------------------------------------------
!                  CALCULATE DIRECTION OF STEEPEST DESCENT
!     ------------------------------------------------------------------
s(1:ndv) = -df(1:ndv)
slope = -dftdf
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS VALUE OF UNITY
!     ------------------------------------------------------------------
40 s1 = 0.
DO  i = 1, ndv
  s2 = ABS(s(i))
  IF (s2 > s1) s1 = s2
END DO
IF (s1 < 1.0D-20) s1 = 1.0D-20
s1 = 1. / s1
dftdf1 = dftdf * s1
s(1:ndv) = s1 * s(1:ndv)
slope = s1 * slope
RETURN
END SUBROUTINE cnmn02



!----- CNMN03
SUBROUTINE cnmn03(x,s,slope,alp,fff,a1,a2,a3,a4,f1,f2,f3,f4,app,  &
                  ncal,kount,jgoto)

! N.B. Argument N1 has been removed.

!  ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH IN UNCONSTRAINED
!  MINIMIZATION USING 2-POINT QUADRATIC INTERPOLATION, 3-POINT
!  CUBIC INTERPOLATION AND 4-POINT CUBIC INTERPOLATION.
!  BY G. N. VANDERPLAATS                         APRIL, 1972.
!  NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!  ALP = PROPOSED MOVE PARAMETER.
!  SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF.
!  SLOPE MUST BE NEGATIVE.
!  OBJ = INITIAL FUNCTION VALUE.

USE common_cnmn1
REAL (dp), INTENT(IN OUT)  :: x(:), s(:), slope, alp, fff, a1, a2, a3, a4, &
                              f1, f2, f3, f4
REAL (dp), INTENT(IN OUT)  :: app
INTEGER, INTENT(IN OUT)    :: ncal(2)
INTEGER, INTENT(OUT)       :: kount
INTEGER, INTENT(IN OUT)    :: jgoto

! Local variables
REAL (dp)  :: aa, ab, ab2, ab3, ap, ff, zro= 0.0_dp
INTEGER    :: i, ii

IF (jgoto /= 0) THEN
  SELECT CASE ( jgoto )
    CASE (    1)
      GO TO 30
    CASE (    2)
      GO TO 50
    CASE (    3)
      GO TO 80
    CASE (    4)
      GO TO 110
    CASE (    5)
      GO TO 150
    CASE (    6)
      GO TO 190
    CASE (    7)
      GO TO 230
  END SELECT
END IF
!     ------------------------------------------------------------------
!                     INITIAL INFORMATION  (ALPHA=0)
!     ------------------------------------------------------------------
IF (slope >= 0.) THEN
  alp = 0.
  RETURN
END IF
IF (iprint > 4) WRITE (6,5000)
fff = obj
a1 = 0.
f1 = obj
a2 = alp
a3 = 0.
f3 = 0.
ap = a2
kount = 0
!     ------------------------------------------------------------------
!            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
!     ------------------------------------------------------------------
10 kount = kount + 1
DO  i = 1, ndv
  x(i) = x(i) + ap * s(i)
END DO
IF (iprint > 4) WRITE (6,5100) ap
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 1
RETURN

30 f2 = obj
IF (iprint > 4) WRITE (6,5300) f2
IF (f2 < f1) GO TO 90
!     ------------------------------------------------------------------
!                     CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
IF (kount <= 5) THEN
  ff = 2. * ABS(f1)
  IF (f2 < ff) GO TO 60
  ff = 5. * ABS(f1)
  IF (f2 >= ff) THEN
    a2 = .5 * a2
    ap = -a2
    alp = a2
    GO TO 10
  END IF
END IF
f3 = f2
a3 = a2
a2 = .5 * a2
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
ap = a2 - alp
alp = a2
DO  i = 1, ndv
  x(i) = x(i) + ap * s(i)
END DO
IF (iprint > 4) WRITE (6,5100) a2
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 2
RETURN

50 f2 = obj
IF (iprint > 4) WRITE (6,5300) f2
!     PROCEED TO CUBIC INTERPOLATION.
GO TO 130
!     ------------------------------------------------------------------
!     **********        2-POINT QUADRATIC INTERPOLATION       **********
!     ------------------------------------------------------------------
60 ii = 1
CALL cnmn04(ii,app,zro,a1,f1,slope,a2,f2,zro,zro,zro,zro)
IF (app < zro .OR. app > a2) GO TO 90
f3 = f2
a3 = a2
a2 = app
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
ap = a2 - alp
alp = a2
DO  i = 1, ndv
  x(i) = x(i) + ap * s(i)
END DO
IF (iprint > 4) WRITE (6,5100) a2
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 3
RETURN

80 f2 = obj
IF (iprint > 4) WRITE (6,5300) f2
GO TO 120

90 a3 = 2. * a2
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
ap = a3 - alp
alp = a3
DO  i = 1, ndv
  x(i) = x(i) + ap * s(i)
END DO
IF (iprint > 4) WRITE (6,5100) a3
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 4
RETURN

110 f3 = obj
IF (iprint > 4) WRITE (6,5300) f3

120 IF (f3 < f2) GO TO 170

!     ------------------------------------------------------------------
!     **********       3-POINT CUBIC INTERPOLATION      **********
!     ------------------------------------------------------------------
130 ii = 3
CALL cnmn04(ii,app,zro,a1,f1,slope,a2,f2,a3,f3,zro,zro)
IF (app < zro .OR. app > a3) GO TO 170
!     ------------------------------------------------------------------
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
!     ------------------------------------------------------------------
ap = app - alp
alp = app
x(1:ndv) = x(1:ndv) + ap * s(1:ndv)
IF (iprint > 4) WRITE (6,5100) alp
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 5
RETURN

150 IF (iprint > 4) WRITE (6,5300) obj
!     ------------------------------------------------------------------
!                         CHECK CONVERGENCE
!     ------------------------------------------------------------------
aa = 1. - app / a2
ab2 = ABS(f2)
ab3 = ABS(obj)
ab = ab2
IF (ab3 > ab) ab = ab3
IF (ab < 1.0D-15) ab = 1.0D-15
ab = (ab2-ab3) / ab
IF (ABS(ab) < 1.0D-15 .AND. ABS(aa) < .001) GO TO 260
a4 = a3
f4 = f3
a3 = app
f3 = obj
IF (a3 > a2) GO TO 200
a3 = a2
f3 = f2
a2 = app
f2 = obj
GO TO 200

!     ------------------------------------------------------------------
!     **********        4-POINT CUBIC INTERPOLATION       **********
!     ------------------------------------------------------------------
170 a4 = 2. * a3
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
ap = a4 - alp
alp = a4
x(1:ndv) = x(1:ndv) + ap * s(1:ndv)
IF (iprint > 4) WRITE (6,5100) alp
IF (iprint > 4) WRITE (6,5200) (x(i),i = 1,ndv)
ncal(1) = ncal(1) + 1
jgoto = 6
RETURN

190 f4 = obj
IF (iprint > 4) WRITE (6,5300) f4
IF (f4 <= f3) THEN
  a1 = a2
  f1 = f2
  a2 = a3
  f2 = f3
  a3 = a4
  f3 = f4
  GO TO 170
END IF

200 ii = 4
CALL cnmn04(ii,app,a1,a1,f1,slope,a2,f2,a3,f3,a4,f4)
IF (app <= a1) THEN
  ap = a1 - alp
  alp = a1
  obj = f1
  DO  i = 1, ndv
    x(i) = x(i) + ap * s(i)
  END DO
  GO TO 240
END IF
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
ap = app - alp
alp = app
x(1:ndv) = x(1:ndv) + ap * s(1:ndv)
IF (iprint > 4) WRITE (6,5100) alp
IF (iprint > 4) WRITE (6,5200) x(1:ndv)
ncal(1) = ncal(1) + 1
jgoto = 7
RETURN

230 IF (iprint > 4) WRITE (6,5300) obj

!     ------------------------------------------------------------------
!                    CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
240 IF (obj <= f2 .AND. obj <= f3) THEN
  IF (obj <= f1) GO TO 260
  ap = a1 - alp
  alp = a1
  obj = f1
ELSE
  IF (f2 >= f3) THEN
    obj = f3
    ap = a3 - alp
    alp = a3
  ELSE
    obj = f2
    ap = a2 - alp
    alp = a2
  END IF
END IF
!     ------------------------------------------------------------------
!                       UPDATE DESIGN VECTOR
!     ------------------------------------------------------------------
x(1:ndv) = x(1:ndv) + ap * s(1:ndv)

!     ------------------------------------------------------------------
!                     CHECK FOR MULTIPLE MINIMA
!     ------------------------------------------------------------------
260 IF (obj > fff) THEN
!     INITIAL FUNCTION IS MINIMUM.
  x(1:ndv) = x(1:ndv) - alp * s(1:ndv)
  alp = 0.
  obj = fff
END IF
jgoto = 0
RETURN
!     ------------------------------------------------------------------
!                                 FORMATS
!     ------------------------------------------------------------------
!
5000 FORMAT (///t6,  &
     '* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
5100 FORMAT (/t6, 'ALPHA =', e14.5/ t6, 'X-VECTOR')
5200 FORMAT (t6, 6E13.5)
5300 FORMAT (/t6, 'OBJ =', e14.5)
END SUBROUTINE cnmn03



!----- CNMN04
SUBROUTINE cnmn04(ii,xbar,eps,x1,y1,slope,x2,y2,x3,y3,x4,y4)

!  ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A MINIMUM
!  OF A ONE-DIMENSIONAL REAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!  BY G. N. VANDERPLAATS                          APRIL, 1972.
!  NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!
!  II = CALCULATION CONTROL.
!       1:  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2 AND Y2.
!       2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3 AND Y3.
!       3:  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2,
!           X3 AND Y3.
!       4:  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3,
!           Y3, X4 AND Y4.
!  EPS MAY BE NEGATIVE.
!  IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS
!  ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR INDICATOR.
!  IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!  INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED,
!  AND II WILL BE CHANGED ACCORDINGLY.

INTEGER, INTENT(IN OUT)  :: ii
REAL (dp), INTENT(OUT)   :: xbar
REAL (dp), INTENT(IN)    :: eps
REAL (dp), INTENT(IN)    :: x1, y1, slope, x2, y2, x3, y3, x4, y4

! Local variables
REAL (dp)  :: aa, bac, bb, cc, dnom, dx, q1, q2, q3, q4, q5, q6, qq,  &
              x11, x111, x21, x22, x222, x31, x32, x33, x41, x42, x44, xbar1
INTEGER    :: nslop

xbar1 = eps - 1.
xbar = xbar1
x21 = x2 - x1
IF (ABS(x21) < 1.0D-20) RETURN
nslop = MOD(ii,2)
SELECT CASE ( ii )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO 20
  CASE (    3)
    GO TO 30
  CASE (    4)
    GO TO 40
END SELECT
!     ------------------------------------------------------------------
!                 II=1: 2-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
10 ii = 1
dx = x1 - x2
IF (ABS(dx) < 1.0D-20) RETURN
aa = (slope + (y2-y1)/dx) / dx
IF (aa < 1.0D-20) RETURN
bb = slope - 2. * aa * x1
xbar = -.5 * bb / aa
IF (xbar < eps) xbar = xbar1
RETURN
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
20 ii = 2
x21 = x2 - x1
x31 = x3 - x1
x32 = x3 - x2
qq = x21 * x31 * x32
IF (ABS(qq) < 1.0D-20) RETURN
aa = (y1*x32-y2*x31+y3*x21) / qq
IF (aa >= 1.0D-20) THEN
  bb = (y2-y1) / x21 - aa * (x1+x2)
  xbar = -.5 * bb / aa
  IF (xbar < eps) xbar = xbar1
  RETURN
END IF
IF (nslop == 0) RETURN
GO TO 10
!     ------------------------------------------------------------------
!                   II=3: 3-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
30 ii = 3
x21 = x2 - x1
x31 = x3 - x1
x32 = x3 - x2
qq = x21 * x31 * x32
IF (ABS(qq) < 1.0D-20) RETURN
x11 = x1 * x1
dnom = x2 * x2 * x31 - x11 * x32 - x3 * x3 * x21
IF (ABS(dnom) < 1.0D-20) GO TO 20
aa = ((x31*x31*(y2-y1) - x21*x21*(y3-y1))/(x31*x21) - slope*x32) / dnom
IF (ABS(aa) < 1.0D-20) GO TO 20
bb = ((y2-y1)/x21 - slope - aa*(x2*x2 + x1*x2 - 2.*x11)) / x21
cc = slope - 3. * aa * x11 - 2. * bb * x1
bac = bb * bb - 3. * aa * cc
IF (bac < 0.) GO TO 20
bac = SQRT(bac)
xbar = (bac-bb) / (3.*aa)
IF (xbar < eps) xbar = eps
RETURN
!     ------------------------------------------------------------------
!                    II=4: 4-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
40 x21 = x2 - x1
x31 = x3 - x1
x41 = x4 - x1
x32 = x3 - x2
x42 = x4 - x2
x11 = x1 * x1
x22 = x2 * x2
x33 = x3 * x3
x44 = x4 * x4
x111 = x1 * x11
x222 = x2 * x22
q2 = x31 * x21 * x32
IF (ABS(q2) < 1.0D-30) RETURN
q1 = x111 * x32 - x222 * x31 + x3 * x33 * x21
q4 = x111 * x42 - x222 * x41 + x4 * x44 * x21
q5 = x41 * x21 * x42
dnom = q2 * q4 - q1 * q5
IF (ABS(dnom) >= 1.0D-30) THEN
  q3 = y3 * x21 - y2 * x31 + y1 * x32
  q6 = y4 * x21 - y2 * x41 + y1 * x42
  aa = (q2*q6 - q3*q5) / dnom
  bb = (q3-q1*aa) / q2
  cc = (y2-y1-aa*(x222-x111)) / x21 - bb * (x1+x2)
  bac = bb * bb - 3. * aa * cc
  IF (ABS(aa) >= 1.0D-20 .AND. bac >= 0.) THEN
    bac = SQRT(bac)
    xbar = (bac-bb) / (3.*aa)
    IF (xbar < eps) xbar = xbar1
    RETURN
  END IF
END IF
IF (nslop == 1) GO TO 30
GO TO 20
END SUBROUTINE cnmn04



!----- CNMN05
SUBROUTINE cnmn05(g,df,a,s,b,c,slope,phi,isc,ic,ms1,nvc,n1,n2,n3,n4,n5)

!  ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF
!  FEASIBLE DIRECTIONS.
!  BY G. N. VANDERPLAATS                            MAY, 1972.
!  NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!  NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1.
!  IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE
!  FIND MODIFIED DIRECTION.

USE common_cnmn1
INTEGER, INTENT(IN)        :: n1, n2, n3, n4, n5
REAL (dp), INTENT(IN OUT)  :: df(:), g(n2), a(n1,n3), s(:), c(n4), b(n3,n3)
REAL (dp), INTENT(IN OUT)  :: slope, phi
INTEGER, INTENT(IN OUT)    :: isc(n2), ic(n3), ms1(n5), nvc

! Local variables
REAL (dp)  :: a1, c1, ct1, ct2, cta, ctam, ctb, ctbm, ctc, ctd, gg,  &
              s1, sg, thmax, tht
INTEGER    :: i, j, j1, k, nac1, nci, ncj, ndb, ndv1, ndv2, ner
!     ------------------------------------------------------------------
!     ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  ***
!     ------------------------------------------------------------------
ndv1 = ndv + 1
ndv2 = ndv + 2
nac1 = nac + 1
nvc = 0
thmax = 0.
cta = ABS(ct)
ct1 = 1. / cta
ctam = ABS(ctmin)
ctb = ABS(ctl)
ct2 = 1. / ctb
ctbm = ABS(ctlmin)
a1 = 1.
DO  i = 1, nac
!     CALCULATE THETA
  nci = ic(i)
  ncj = 1
  IF (nci <= ncon) ncj = isc(nci)
  c1 = g(nci)
  ctd = ct1
  ctc = ctam
  IF (ncj > 0) THEN
    ctc = ctbm
    ctd = ct2
  END IF
  IF (c1 > ctc) nvc = nvc + 1
  tht = 0.
  gg = 1. + ctd * c1
  IF (ncj == 0 .OR. c1 > ctc) tht = theta * gg * gg
  IF (tht > 50.) tht = 50.
  IF (tht > thmax) thmax = tht
  a(ndv1,i) = tht
!     ------------------------------------------------------------------
!                    NORMALIZE GRADIENTS OF CONSTRAINTS
!     ------------------------------------------------------------------
  a(ndv2,i) = 1.
  IF (nci <= ncon) THEN
    a1 = 0.
    DO  j = 1, ndv
      a1 = a1 + a(j,i) ** 2
    END DO
    IF (a1 < 1.0D-20) a1 = 1.0D-20
    a1 = SQRT(a1)
    a(ndv2,i) = a1
    a1 = 1. / a1
    DO  j = 1, ndv
      a(j,i) = a1 * a(j,i)
    END DO
  END IF
END DO
!     ------------------------------------------------------------------
!     CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV.
!     ------------------------------------------------------------------
i = 0
40 i = i + 1
50 IF (a(ndv2,i) <= 1.0D-6) THEN
!     ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE.
  IF (iprint >= 2) WRITE (6,5000) ic(i)
!     REDUCE NAC BY ONE.
  nac = nac - 1
!     SHIFT COLUMNS OF A AND ROWS OF IC IF I.LE.NAC.
  IF (i > nac) GO TO 80
!     SHIFT.
  DO  j = i, nac
    j1 = j + 1
    ic(j) = ic(j1)
    DO  k = 1, ndv2
      a(k,j) = a(k,j1)
    END DO
  END DO
  IF (i <= nac) GO TO 50
END IF
IF (i < nac) GO TO 40

80 IF (nac <= 0) RETURN
nac1 = nac + 1
!     DETERMINE IF CONSTRAINTS ARE VIOLATED.
nvc = 0
DO  i = 1, nac
  nci = ic(i)
  ncj = 1
  IF (nci <= ncon) ncj = isc(nci)
  ctc = ctam
  IF (ncj > 0) ctc = ctbm
  IF (g(nci) > ctc) nvc = nvc + 1
END DO
!     ------------------------------------------------------------------
!     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
!     COLUMN OF A
!     ------------------------------------------------------------------
a1 = 0.
DO  i = 1, ndv
  a1 = a1 + df(i) ** 2
END DO
IF (a1 < 1.0D-20) a1 = 1.0D-20
a1 = SQRT(a1)
a1 = 1. / a1
DO  i = 1, ndv
  a(i,nac1) = a1 * df(i)
END DO
!     BUILD C VECTOR.
IF (nvc <= 0) THEN
!     ------------------------------------------------------------------
!                 BUILD C FOR CLASSICAL METHOD
!     ------------------------------------------------------------------
  ndb = nac1
  a(ndv1,ndb) = 1.
  DO  i = 1, ndb
    c(i) = -a(ndv1,i)
  END DO
ELSE
!     ------------------------------------------------------------------
!                   BUILD C FOR MODIFIED METHOD
!     ------------------------------------------------------------------
  ndb = nac
  a(ndv1,nac1) = -phi
!     ------------------------------------------------------------------
!           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
!     ------------------------------------------------------------------
  IF (thmax > 0.00001) thmax = 1. / thmax
  DO  i = 1, ndb
    a(ndv1,i) = a(ndv1,i) * thmax
  END DO
  DO  i = 1, ndb
    c(i) = 0.
    DO  j = 1, ndv1
      c(i) = c(i) + a(j,i) * a(j,nac1)
    END DO
  END DO
END IF
!     ------------------------------------------------------------------
!                      BUILD B MATRIX
!     ------------------------------------------------------------------
DO  i = 1, ndb
  DO  j = 1, ndb
    b(i,j) = 0.
    DO  k = 1, ndv1
      b(i,j) = b(i,j) - a(k,i) * a(k,j)
    END DO
  END DO
END DO
!     ------------------------------------------------------------------
!                    SOLVE SPECIAL L. P. PROBLEM
!     ------------------------------------------------------------------
CALL cnmn08(ndb,ner,c,ms1,b,n3,n4)
IF (iprint > 1 .AND. ner > 0) WRITE (6,5200)
!     CALCULATE RESULTING DIRECTION VECTOR, S.
slope = 0.
!     ------------------------------------------------------------------
!                  USABLE-FEASIBLE DIRECTION
!     ------------------------------------------------------------------
DO  i = 1, ndv
  s1 = 0.
  IF (nvc > 0) s1 = -a(i,nac1)
  DO  j = 1, ndb
    s1 = s1 - a(i,j) * c(j)
  END DO
  slope = slope + s1 * df(i)
  s(i) = s1
END DO
s(ndv1) = 1.
IF (nvc > 0) s(ndv1) = -a(ndv1,nac1)
DO  j = 1, ndb
  s(ndv1) = s(ndv1) - a(ndv1,j) * c(j)
END DO
!     ------------------------------------------------------------------
!     CHECK TO INSURE THE S-VECTOR IS FEASIBLE.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
DO  j = 1, nac
!     S DOT DEL(G).
  sg = DOT_PRODUCT( s(1:ndv), a(1:ndv,j) )
!     IF(SG.GT.0.) GO TO 176
!
!  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS
!
  IF (sg > 1.0D-04) GO TO 240
!     FEASIBLE FOR THIS CONSTRAINT.  CONTINUE.
END DO
GO TO 250

!     S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM.
240 IF (iprint >= 2) WRITE (6,5100)
s(ndv1) = 0.
nvc = 0
RETURN
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS OF UNITY
!     ------------------------------------------------------------------
250 s1 = 0.
DO  i = 1, ndv
  a1 = ABS(s(i))
  IF (a1 > s1) s1 = a1
END DO
!     IF (S1.LT.1.0E-10) RETURN
!
!  E-10 CHANGED TO E-04 ON 1/12/81
!
IF (s1 < 1.0D-04) RETURN
s1 = 1. / s1
DO  i = 1, ndv
  s(i) = s1 * s(i)
END DO
slope = s1 * slope
s(ndv1) = s1 * s(ndv1)
RETURN
5000 FORMAT (t6,'** CONSTRAINT',i5,' HAS ZERO GRADIENT'/t6,   &
    'DELETED FROM ACTIVE SET')
5100 FORMAT (t6,'** CALCULATED S-VECTOR IS NOT FEASIBLE'/t6,   &
    'BETA IS SET TO ZERO')
!     ------------------------------------------------------------------
!                           FORMATS
!     ------------------------------------------------------------------
!
!
5200 FORMAT (//t6, '* * DIRECTION FINDING PROCESS DID NOT CONVERGE'/t6,   &
    '* * S-VECTOR MAY NOT BE VALID')
END SUBROUTINE cnmn05



!----- CNMN06
SUBROUTINE cnmn06(x,vlb,vub,g,scal,df,s,g1,g2,ctam,ctbm,slope,alp,a2,a3,a4, &
                  f1,f2,f3,cv1,cv2,cv3,cv4,alpca,alpfes,alpln,alpmin,alpnc, &
                  alpsav,alpsid,alptot,isc,ncal,nvc,icount,igood1,  &
                  igood2,igood3,igood4,ibest,iii,nlnc,jgoto)

! N.B. Arguments N1 & N2 have been removed.

!  ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED
!  FUNCTION MINIMIZATION.
!  BY G. N. VANDERPLAATS                           AUG., 1974.
!  NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!  OBJ = INITIAL AND FINAL FUNCTION VALUE.
!  ALP = MOVE PARAMETER.
!  SLOPE = INITIAL SLOPE.
!
!  ALPSID = MOVE TO SIDE CONSTRAINT.
!  ALPFES = MOVE TO FEASIBLE REGION.
!  ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT.
!  ALPLN = MOVE TO LINEAR CONSTRAINT.
!  ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
!  ALPMIN = MOVE TO MINIMIZE FUNCTION.
!  ALPTOT = TOTAL MOVE PARAMETER.

USE common_cnmn1
REAL (dp), INTENT(IN OUT)  :: x(:), vlb(:), vub(:), g(:), scal(:), df(:), &
                              s(:), g1(:), g2(:), ctam, ctbm, slope, alp, &
                              a2, a3, a4, f1, f2, f3, cv1, cv2, cv3, cv4, &
                              alpca, alpfes, alpln, alpmin, alpnc, alpsav, &
                              alpsid, alptot
INTEGER, INTENT(IN OUT)    :: isc(:), ncal(2), nvc, icount, igood1, igood2, &
                              igood3, igood4, ibest, iii, nlnc, jgoto

! Local variables
REAL (dp)  :: alpa, alpb, c1, c2, c3, cc, f4, gi, si, xi, xi1, xi2, zro = 0.0_dp
INTEGER    :: i, ii, jbest, ksid, nvc1

IF (jgoto /= 0) THEN
  SELECT CASE ( jgoto )
    CASE (    1)
      GO TO 70
    CASE (    2)
      GO TO 140
    CASE (    3)
      GO TO 230
  END SELECT
END IF
IF (iprint >= 5) WRITE (6,5100)
alpsav = alp
icount = 0
alptot = 0.
!     TOLERANCES.
ctam = ABS(ctmin)
ctbm = ABS(ctlmin)
!     PROPOSED MOVE.
!     ------------------------------------------------------------------
!     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  *****
!     ------------------------------------------------------------------
10 a2 = alpsav
icount = icount + 1
alpsid = 1.0D+20
!     INITIAL ALPHA AND OBJ.
alp = 0.
f1 = obj
ksid = 0
IF (nside /= 0) THEN
!     ------------------------------------------------------------------
!     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
!     SIDE CONSTRAINTS
!     ------------------------------------------------------------------
  DO  i = 1, ndv
    si = s(i)
    IF (ABS(si) <= 1.0D-20) THEN
!     ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
      s(i) = 0.
      slope = slope - si * df(i)
    ELSE
      xi = x(i)
      si = 1. / si
      IF (si <= 0.) THEN
!     LOWER BOUND.
        xi2 = vlb(i)
        xi1 = ABS(xi2)
        IF (xi1 < 1.) xi1 = 1.
!     CONSTRAINT VALUE.
        gi = (xi2-xi) / xi1
        IF (gi > -1.0D-6) GO TO 20
!     PROPOSED MOVE TO LOWER BOUND.
        alpa = (xi2-xi) * si
        IF (alpa < alpsid) alpsid = alpa
        CYCLE
      END IF
!     UPPER BOUND.
      xi2 = vub(i)
      xi1 = ABS(xi2)
      IF (xi1 < 1.) xi1 = 1.
!     CONSTRAINT VALUE.
      gi = (xi-xi2) / xi1
      IF (gi <= -1.0D-6) THEN
!     PROPOSED MOVE TO UPPER BOUND.
        alpa = (xi2-xi) * si
        IF (alpa < alpsid) alpsid = alpa
        CYCLE
      END IF

!     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0.
      20 slope = slope - s(i) * df(i)
      s(i) = 0.
      ksid = ksid + 1
    END IF
  END DO
!     ALPSID IS UPPER BOUND ON ALPHA.
  IF (a2 > alpsid) a2 = alpsid
END IF
!     ------------------------------------------------------------------
!               CHECK ILL-CONDITIONING
!     ------------------------------------------------------------------
IF (ksid == ndv .OR. icount > 10) GO TO 340
IF (nvc == 0 .AND. slope > 0.) GO TO 340
alpfes = -1.
alpmin = -1.
alpln = 1.1 * alpsid
alpnc = alpsid
alpca = alpsid
IF (ncon /= 0) THEN
!     STORE CONSTRAINT VALUES IN G1.
  DO  i = 1, ncon
    g1(i) = g(i)
  END DO
END IF
!     ------------------------------------------------------------------
!                  MOVE A DISTANCE A2*S
!     ------------------------------------------------------------------
alptot = alptot + a2
DO  i = 1, ndv
  x(i) = x(i) + a2 * s(i)
END DO
IF (iprint >= 5) THEN
  WRITE (6,5200) a2
  IF (nscal /= 0) THEN
    DO  i = 1, ndv
      g(i) = scal(i) * x(i)
    END DO
    WRITE (6,5300) (g(i),i = 1,ndv)
  ELSE
    WRITE (6,5300) (x(i),i = 1,ndv)
  END IF
END IF
!     ------------------------------------------------------------------
!                   UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
ncal(1) = ncal(1) + 1
jgoto = 1
RETURN

70 f2 = obj
IF (iprint >= 5) WRITE (6,5400) f2
IF (iprint >= 5 .AND. ncon /= 0) THEN
  WRITE (6,5500)
  WRITE (6,5300) (g(i),i = 1,ncon)
END IF
!     ------------------------------------------------------------------
!               IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2
!     ------------------------------------------------------------------
!     IGOOD = 0 IS ACCAPTABLE.
!     CV = MAXIMUM CONSTRAINT VIOLATION.
igood1 = 0
igood2 = 0
cv1 = 0.
cv2 = 0.
nvc1 = 0
IF (ncon /= 0) THEN
  DO  i = 1, ncon
    cc = ctam
    IF (isc(i) > 0) cc = ctbm
    c1 = g1(i) - cc
    c2 = g(i) - cc
    IF (c2 > 0.) nvc1 = nvc1 + 1
    IF (c1 > cv1) cv1 = c1
    IF (c2 > cv2) cv2 = c2
  END DO
  IF (cv1 > 0.) igood1 = 1
  IF (cv2 > 0.) igood2 = 1
END IF
alp = a2
obj = f2
!     ------------------------------------------------------------------
!     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT
!     VIOLATIONS RETURN
!     ------------------------------------------------------------------
IF (nvc1 < nvc .AND. nvc1 > 0) GO TO 340
!     ------------------------------------------------------------------
!             IDENTIFY BEST OF DESIGNS F1 ANF F2
!     ------------------------------------------------------------------
!     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
!     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM
!     CONSTRAINT VIOLATION.
IF (igood1 /= 0 .OR. igood2 /= 0) THEN
!     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
  ibest = 1
  IF (cv1 >= cv2) ibest = 2
ELSE
!     NO CONSTRAINT VIOLATION.  PICK MINIMUM F.
  ibest = 1
  IF (f2 <= f1) ibest = 2
END IF
ii = 1
!     ------------------------------------------------------------------
!     IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
IF (cv2 > cv1) THEN
  alpln = a2
  alpnc = a2
  alpca = a2
END IF
IF (ncon /= 0) THEN
!     ------------------------------------------------------------------
!     *****                 2 - POINT INTERPOLATION                *****
!     ------------------------------------------------------------------
  iii = 0
  90 iii = iii + 1
  c1 = g1(iii)
  c2 = g(iii)
  IF (isc(iii) /= 0) THEN
!     ------------------------------------------------------------------
!                        LINEAR CONSTRAINT
!     ------------------------------------------------------------------
    IF (c1 >= 1.0D-5 .AND. c1 <= ctbm) GO TO 100
    CALL cnmn07(ii,alp,zro,zro,c1,a2,c2,zro,zro)
    IF (alp <= 0.) GO TO 100
    IF (c1 > ctbm .AND. alp > alpfes) alpfes = alp
    IF (c1 < ctl .AND. alp < alpln) alpln = alp
  ELSE
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
    IF (c1 < 1.0D-5 .OR. c1 > ctam) THEN
      CALL cnmn07(ii,alp,zro,zro,c1,a2,c2,zro,zro)
      IF (alp > 0.) THEN
        IF (c1 > ctam .AND. alp > alpfes) alpfes = alp
        IF (c1 < ct .AND. alp < alpnc) alpnc = alp
      END IF
    END IF
  END IF

  100 IF (iii < ncon) GO TO 90
END IF
IF (linobj <= 0 .AND. slope < 0.) THEN
!     CALCULATE ALPHA TO MINIMIZE FUNCTION.
  CALL cnmn04(ii,alpmin,zro,zro,f1,slope,a2,f2,zro,zro,zro,zro)
END IF
!     ------------------------------------------------------------------
!                         PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
a3 = alpfes
!     MOVE TO MINIMIZE FUNCTION.
IF (alpmin > a3) a3 = alpmin
!     IF A3.LE.0, SET A3 = ALPSID.
IF (a3 <= 0.) a3 = alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
IF (a3 > alpnc) a3 = alpnc
IF (a3 > alpln) a3 = alpln
!     MAKE A3 NON-ZERO.
IF (a3 <= 1.0D-20) a3 = 1.0D-20
!     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
alpb = 1. - a2 / a3
alpa = 1. - alpsid / a3
jbest = 0
IF (ABS(alpb) < 1.0D-10 .AND. ABS(alpa) < 1.0D-10) jbest = 1
IF (jbest == 1 .AND. ibest == 2) GO TO 10
!     SIDE CONSTRAINT CHECK NOT SATISFIED.
IF (ncon /= 0) THEN
!     STORE CONSTRAINT VALUES IN G2.
  DO  i = 1, ncon
    g2(i) = g(i)
  END DO
END IF
!     IF A3=A2, SET A3=.9*A2.
IF (ABS(alpb) < 1.0D-10) a3 = .9 * a2
!     MOVE AT LEAST .01*A2.
IF (a3 < (.01*a2)) a3 = .01 * a2
!     LIMIT MOVE TO 5.*A2.
IF (a3 > (5.*a2)) a3 = 5. * a2
!     LIMIT MOVE TO ALPSID.
IF (a3 > alpsid) a3 = alpsid
!     MOVE A DISTANCE A3*S.
alp = a3 - a2
alptot = alptot + alp
DO  i = 1, ndv
  x(i) = x(i) + alp * s(i)
END DO
IF (iprint >= 5) THEN
  WRITE (6,5600)
  WRITE (6,5200) a3
  IF (nscal /= 0) THEN
    g(1:ndv) = scal(1:ndv) * x(1:ndv)
    WRITE (6,5300) g(1:ndv)
  ELSE
    WRITE (6,5300) x(1:ndv)
  END IF
END IF
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
ncal(1) = ncal(1) + 1
jgoto = 2
RETURN

140 f3 = obj
IF (iprint >= 5) WRITE (6,5400) f3
IF (iprint >= 5 .AND. ncon /= 0) THEN
  WRITE (6,5500)
  WRITE (6,5300) (g(i),i = 1,ncon)
END IF
!     ------------------------------------------------------------------
!       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN
!     ------------------------------------------------------------------
cv3 = 0.
igood3 = 0
nvc1 = 0
IF (ncon /= 0) THEN
  DO  i = 1, ncon
    cc = ctam
    IF (isc(i) > 0) cc = ctbm
    c1 = g(i) - cc
    IF (c1 > cv3) cv3 = c1
    IF (c1 > 0.) nvc1 = nvc1 + 1
  END DO
  IF (cv3 > 0.) igood3 = 1
END IF
!     DETERMINE BEST DESIGN.
IF (ibest /= 2) THEN
!     CHOOSE BETWEEN F1 AND F3.
  IF (igood1 /= 0 .OR. igood3 /= 0) THEN
    IF (cv1 >= cv3) ibest = 3
    GO TO 160
  END IF
  IF (f3 <= f1) ibest = 3
ELSE
!     CHOOSE BETWEEN F2 AND F3.
  IF (igood2 /= 0 .OR. igood3 /= 0) THEN
    IF (cv2 >= cv3) ibest = 3
  ELSE
    IF (f3 <= f2) ibest = 3
  END IF
END IF

160 alp = a3
obj = f3
!     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
IF (nvc1 < nvc) GO TO 340
!     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
IF (linobj /= 0 .AND. nlnc == ncon) GO TO 340
!     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
alpb = 1. - alpln / a3
IF (ABS(alpb) < 1.0D-20 .AND. ibest == 3 .AND. igood3 == 0) GO TO 340
!     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT MODIFICATION.
alpa = 1. - alpsid / a3
IF (ABS(alpa) < 1.0D-20 .AND. ibest == 3) GO TO 10
!     ------------------------------------------------------------------
!     **********            3 - POINT INTERPOLATION            *********
!     ------------------------------------------------------------------
alpnc = alpsid
alpca = alpsid
alpfes = -1.
alpmin = -1.
!     ------------------------------------------------------------------
!     IF A3 IS GREATER THAN A2 AND CV3 IS GREATER THAN CV2, SET
!     MOVE LIMITS TO A3.  PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
IF (a3 > a2 .AND. cv3 > cv2) THEN
  alpln = a3
  alpnc = a3
  alpca = a3
END IF
IF (ncon /= 0) THEN
  iii = 0
  170 iii = iii + 1
  c1 = g1(iii)
  c2 = g2(iii)
  c3 = g(iii)
  IF (isc(iii) /= 0) THEN
!     ------------------------------------------------------------------
!     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
!     ------------------------------------------------------------------
    IF (c1 <= ctbm) GO TO 190
    ii = 1
    CALL cnmn07(ii,alp,zro,zro,c1,a3,c3,zro,zro)
    IF (alp > alpfes) alpfes = alp
  ELSE
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
    ii = 2
    CALL cnmn07(ii,alp,zro,zro,c1,a2,c2,a3,c3)
    IF (alp > zro) THEN
      IF (c1 < ct .OR. c1 > 0.) THEN
        IF (c1 > ctam .OR. c1 < 0.) GO TO 180
      END IF
!     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER.
      alpa = alp
      CALL cnmn07(ii,alp,alpa,zro,c1,a2,c2,a3,c3)
      IF (alp < alpca .AND. alp >= alpa) alpca = alp
      GO TO 190

      180 IF (alp > alpfes .AND. c1 > ctam) alpfes = alp
      IF (alp < alpnc .AND. c1 < 0.) alpnc = alp
    END IF
  END IF

  190 IF (iii < ncon) GO TO 170
END IF
IF (linobj <= 0 .AND. slope <= 0.) THEN
!     ------------------------------------------------------------------
!              CALCULATE ALPHA TO MINIMIZE FUNCTION
!     ------------------------------------------------------------------
  ii = 3
  IF (a2 > a3 .AND. (igood2 == 0 .AND. ibest == 2)) ii = 2
  CALL cnmn04(ii,alpmin,zro,zro,f1,slope,a2,f2,a3,f3,zro,zro)
END IF
!     ------------------------------------------------------------------
!                       PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
a4 = alpfes
!     MOVE TO MINIMIZE FUNCTION.
IF (alpmin > a4) a4 = alpmin
!     IF A4.LE.0, SET A4 = ALPSID.
IF (a4 <= 0.) a4 = alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
IF (a4 > alpln) a4 = alpln
IF (a4 > alpnc) a4 = alpnc
!     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
IF (a4 > alpca) a4 = alpca
!     LIMIT A4 TO 5.*A3.
IF (a4 > (5.*a3)) a4 = 5. * a3
!     UPDATE DESIGN.
IF (ibest == 3 .AND. ncon /= 0) THEN
!     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT.
  DO  i = 1, ncon
    g2(i) = g(i)
  END DO
END IF
!     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3.
alp = a4 - a3
IF (igood1 == 0 .AND. igood3 == 1 .AND. ABS(alp) < 1.0D-20) a4 = .9 * a3
!     ------------------------------------------------------------------
!                   MOVE A DISTANCE A4*S
!     ------------------------------------------------------------------
alp = a4 - a3
alptot = alptot + alp
DO  i = 1, ndv
  x(i) = x(i) + alp * s(i)
END DO
IF (iprint >= 5) THEN
  WRITE (6,5000)
  WRITE (6,5200) a4
  IF (nscal /= 0) THEN
    g(1:ndv) = scal(1:ndv) * x(1:ndv)
    WRITE (6,5300) g(1:ndv)
  ELSE
    WRITE (6,5300) x(1:ndv)
  END IF
END IF
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
ncal(1) = ncal(1) + 1
jgoto = 3
RETURN

230 f4 = obj
IF (iprint >= 5) WRITE (6,5400) f4
IF (iprint >= 5 .AND. ncon /= 0) THEN
  WRITE (6,5500)
  WRITE (6,5300) g(1:ncon)
END IF
!     DETERMINE ACCAPTABILITY OF F4.
igood4 = 0
cv4 = 0.
IF (ncon /= 0) THEN
  DO  i = 1, ncon
    cc = ctam
    IF (isc(i) > 0) cc = ctbm
    c1 = g(i) - cc
    IF (c1 > cv4) cv4 = c1
  END DO
  IF (cv4 > 0.) igood4 = 1
END IF
alp = a4
obj = f4
!     ------------------------------------------------------------------
!                     DETERMINE BEST DESIGN
!     ------------------------------------------------------------------
SELECT CASE ( ibest )
  CASE (    1)
    !     CHOOSE BETWEEN F1 AND F4.
    IF (igood1 /= 0 .OR. igood4 /= 0) THEN
      IF (cv1 > cv4) GO TO 340
    ELSE
      IF (f4 <= f1) GO TO 340
    END IF
    !     F1 IS BEST.
    alptot = alptot - a4
    obj = f1
    x(1:ndv) = x(1:ndv) - a4 * s(1:ndv)
    IF (ncon == 0) GO TO 340
    g(1:ncon) = g1(1:ncon)

  CASE (    2)
    !     CHOOSE BETWEEN F2 AND F4.
    IF (igood2 /= 0 .OR. igood4 /= 0) THEN
      IF (cv2 > cv4) GO TO 340
    ELSE
      IF (f4 <= f2) GO TO 340
    END IF
    !     F2 IS BEST.
    obj = f2
    a2 = a4 - a2
    alptot = alptot - a2
    x(1:ndv) = x(1:ndv) - a2 * s(1:ndv)
    IF (ncon == 0) GO TO 340
    g(1:ncon) = g2(1:ncon)

  CASE (    3)
    !     CHOOSE BETWEEN F3 AND F4.
    IF (igood3 /= 0 .OR. igood4 /= 0) THEN
      IF (cv3 > cv4) GO TO 340
    ELSE
      IF (f4 <= f3) GO TO 340
    END IF
    !     F3 IS BEST.
    obj = f3
    a3 = a4 - a3
    alptot = alptot - a3
    x(1:ndv) = x(1:ndv) - a3 * s(1:ndv)
    IF (ncon /= 0) THEN
      g(1:ncon) = g2(1:ncon)
    END IF

END SELECT

340 alp = alptot
IF (iprint >= 5) WRITE (6,5700)
jgoto = 0
RETURN
!     ------------------------------------------------------------------
!                                  FORMATS
!     ------------------------------------------------------------------
!
!
5000 FORMAT (/t6, 'THREE-POINT INTERPOLATION')
5100 FORMAT (/// '* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
5200 FORMAT (//t6, 'PROPOSED DESIGN'/ t6, 'ALPHA =', e12.5/ t6, 'X-VECTOR')
5300 FORMAT (' ', 8E12.4)
5400 FORMAT (/t6, 'OBJ =', e13.5)
5500 FORMAT (/t6, 'CONSTRAINT VALUES')
5600 FORMAT (/t6, 'TWO-POINT INTERPOLATION')
5700 FORMAT (/t6, '* * * END OF ONE-DIMENSIONAL SEARCH')
END SUBROUTINE cnmn06



!----- CNMN07
SUBROUTINE cnmn07(ii,xbar,eps,x1,y1,x2,y2,x3,y3)

!  ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A REAL ZERO
!  OF A ONE-DIMENSIONAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!  BY G. N. VANDERPLAATS                          APRIL, 1972.
!  NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!  II = CALCULATION CONTROL.
!       1:  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2.
!       2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3 AND Y3.
!  EPS MAY BE NEGATIVE.
!  IF REQUIRED ZERO ON Y DOES NOT EXITS, OR THE FUNCTION IS
!  ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR INDICATOR.
!  IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!  INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED AND
!  II WILL BE CHANGED ACCORDINGLY.

INTEGER, INTENT(IN OUT)    :: ii
REAL (dp), INTENT(OUT)     :: xbar
REAL (dp), INTENT(IN)      :: eps, x1, y1, x2, y2, x3, y3

! Local variables
REAL (dp)  :: aa, bac, bb, cc, dy, qq, x21, x31, x32, xb2, xbar1, yy
INTEGER    :: jj

xbar1 = eps - 1.
xbar = xbar1
jj = 0
x21 = x2 - x1
IF (ABS(x21) < 1.0D-20) RETURN
IF (ii == 2) GO TO 20
!
!     ------------------------------------------------------------------
!                  II=1: 2-POINT LINEAR INTERPOLATION
!     ------------------------------------------------------------------
10 ii = 1
yy = y1 * y2
IF (jj /= 0 .AND. yy >= 0.) THEN
!     INTERPOLATE BETWEEN X2 AND X3.
  dy = y3 - y2
  IF (ABS(dy) >= 1.0D-20) THEN
    xbar = x2 + y2 * (x2-x3) / dy
    IF (xbar < eps) xbar = xbar1
    RETURN
  END IF
END IF
dy = y2 - y1

!     INTERPOLATE BETWEEN X1 AND X2.
IF (ABS(dy) < 1.0D-20) RETURN
xbar = x1 + y1 * (x1-x2) / dy
IF (xbar < eps) xbar = xbar1
RETURN
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
20 jj = 1
x31 = x3 - x1
x32 = x3 - x2
qq = x21 * x31 * x32
IF (ABS(qq) < 1.0D-20) RETURN
aa = (y1*x32-y2*x31+y3*x21) / qq
IF (ABS(aa) < 1.0D-20) GO TO 10
bb = (y2-y1) / x21 - aa * (x1+x2)
cc = y1 - x1 * (aa*x1+bb)
bac = bb * bb - 4. * aa * cc
IF (bac < 0.) GO TO 10
bac = SQRT(bac)
aa = .5 / aa
xbar = aa * (bac-bb)
xb2 = -aa * (bac+bb)
IF (xbar < eps) xbar = xb2
IF (xb2 < xbar .AND. xb2 > eps) xbar = xb2
IF (xbar < eps) xbar = xbar1
RETURN
END SUBROUTINE cnmn07



!----- CNMN08
SUBROUTINE cnmn08(ndb,ner,c,ms1,b,n3,n4)

! N.B. Argument N5 has been removed.

!  ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE
!  TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS.
!  BY G. N. VANDERPLAATS                             APRIL, 1972.
!  NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!  REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS',
!  G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS
!  AND STRUCTURES, VOL 3, PP 739-755, 1973.
!  FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR
!  U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE
!  U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0.
!  NER = ERROR FLAG.  IF NER.NE.0 ON RETURN, PROCESS HAS NOT
!  CONVERGED IN 5*NDB ITERATIONS.
!  VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES.
!  ------------------------------------------------------------------
!  CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1
!  ------------------------------------------------------------------

INTEGER, INTENT(IN)        :: ndb, n3, n4
INTEGER, INTENT(OUT)       :: ner, ms1(:)
REAL (dp), INTENT(IN OUT)  :: c(n4), b(n3,n3)

! Local variables
INTEGER    :: i, ichk, iter1, j, jj, kk, m2, nmax
REAL (dp)  :: bb, bb1, bi, c1, cb, cbmax, cbmin, eps

ner = 1
m2 = 2 * ndb
!     CALCULATE CBMIN AND EPS AND INITIALIZE MS1.
eps = -1.0D+10
cbmin = 0.
DO  i = 1, ndb
  bi = b(i,i)
  cbmax = 0.
  IF (bi < -1.0D-6) cbmax = c(i) / bi
  IF (bi > eps) eps = bi
  IF (cbmax > cbmin) cbmin = cbmax
  ms1(i) = 0
END DO
eps = .0001 * eps
!     IF (EPS.LT.-1.0E-10) EPS=-1.0E-10
!
!  E-10 CHANGED TO E-03 ON 1/12/81
!
IF (eps < -1.0D-03) eps = -1.0D-03
IF (eps > -.0001) eps = -.0001
cbmin = cbmin * 1.0D-6
!     IF (CBMIN.LT.1.0D-10) CBMIN=1.0D-10
!
!  E-10 CHANGED TO E-05 ON 1/12/81
!
IF (cbmin < 1.0D-05) cbmin = 1.0D-05
iter1 = 0
nmax = 5 * ndb
!     ------------------------------------------------------------------
!     **********             BEGIN NEW ITERATION              **********
!     ------------------------------------------------------------------
20 iter1 = iter1 + 1
IF (iter1 > nmax) RETURN
!     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
cbmax = .9 * cbmin
ichk = 0
DO  i = 1, ndb
  c1 = c(i)
  bi = b(i,i)
!     IF (BI.GT.EPS .OR. C1.GT.0.) GO TO 30
  IF (bi <= eps .AND. c1 <= -1.0D-05) THEN
!
!  0. CHANGED TO -1.0E-05 ON 1/12/81
!
    cb = c1 / bi
    IF (cb > cbmax) THEN
      ichk = i
      cbmax = cb
    END IF
  END IF
END DO
IF (cbmax >= cbmin) THEN
  IF (ichk /= 0) THEN
!     UPDATE VECTOR MS1.
    jj = ichk
    IF (ms1(jj) == 0) jj = ichk + ndb
    kk = jj + ndb
    IF (kk > m2) kk = jj - ndb
    ms1(kk) = ichk
    ms1(jj) = 0
!     ------------------------------------------------------------------
!                     PIVOT OF B(ICHK,ICHK)
!     ------------------------------------------------------------------
    bb = 1. / b(ichk,ichk)
    DO  j = 1, ndb
      b(ichk,j) = bb * b(ichk,j)
    END DO
    c(ichk) = cbmax
    b(ichk,ichk) = bb
!     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
!     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.
    DO  i = 1, ndb
      IF (i /= ichk) THEN
        bb1 = b(i,ichk)
        b(i,ichk) = 0.
        DO  j = 1, ndb
          b(i,j) = b(i,j) - bb1 * b(ichk,j)
        END DO
        c(i) = c(i) - bb1 * cbmax
      END IF
    END DO
    GO TO 20
  END IF
END IF
ner = 0
!     ------------------------------------------------------------------
!     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR
!     TEMPORARY STORAGE
!     ------------------------------------------------------------------
b(1:ndb,1) = c(1:ndb)
DO  i = 1, ndb
  c(i) = 0.
  j = ms1(i)
  IF (j > 0) c(i) = b(j,1)
  IF (c(i) < 0.) c(i) = 0.
END DO
RETURN
END SUBROUTINE cnmn08

END MODULE Constrained_minimization

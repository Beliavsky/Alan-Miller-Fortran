MODULE probn
! Replaces two commons in program MDR

! COMMON /probn/ nprob
! COMMON /h3r/ r

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, SAVE    :: nprob
REAL (dp), SAVE  :: r
END MODULE probn



PROGRAM mdr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-17  Time: 10:28:53

!  MASTER DRIVER FOR MORE, GARBOW AND HILLSTROM PROBLEMS

!  AS GIVEN NOW NTEST RUNS FROM 1 TO 2 THUS TESTING ONLY NEWTON'S
!  BASIC PROGRAM.  IT CAN BE RUN FROM 1 TO 9 TO TEST DIFFERENT ASPECTS.

!  PROGRAM RUNS QUICKLY EXCEPT FOR WATSON PROBLEM (#6 AND 7).

!  SUFFIX COMM ON OUTPUT FILES INDICATES FOR COMMUNICATION TO
!  COMPUTERS AND CHEM. ENG.

USE Equation_Solver
USE probn, ONLY: nprob, r
IMPLICIT NONE

INTEGER, PARAMETER  :: mgll = 10, n = 40, nunit = 10
REAL (dp)  :: jac(n,n), lam0, mstpf, nsttol
INTEGER    :: acptcr, i, icount, itsclf, itsclx, jactyp, jupdm, maxexp,  &
              maxit, maxns, maxqns, method, minqns, narmij, ndim, niejev,  &
              njacch, njetot, njupd, ntest, output, qnupdm, stopcr, supprs, &
              trmcod, trupdm
REAL (dp)  :: a(n,n), alpha, boundl(n), boundu(n), confac, delf(n), delfac, &
              delta, epsmch, etafac, fcnnew, fdtolj, fsave(n), ftol,   &
              ftrack(0:mgll-1), fvec(n), fvecc(n), h(n,n), hh, hhpi(n),  &
              omega, plee(n,n), ratiof, rdiag(n), s(n), sbar(n), scalef(n),  &
              scalex(n), sigma, sn(n), ssdhat(n), strack(0:mgll-1), stptol,  &
              ti, vhat(n), xc(n), xplus(n), xsave(n)
LOGICAL    :: absnew, cauchy, deuflh, geoms, linesr, newton, overch
CHARACTER (LEN=6)   :: help
CHARACTER (LEN=25)  :: NAME
INTEGER    :: nfe(6), isum(6)

! EXTERNAL fcn, jacob
! COMMON /nnes_1/ matsup
! COMMON /nnes_2/ wrnsup
! COMMON /nnes_3/ bypass
! COMMON /nnes_4/ nfetot
! COMMON /probn/ nprob
! COMMON /h3r/ r

INTERFACE
  SUBROUTINE fcn (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fcn

  SUBROUTINE jacob (overfl, n, jac, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: jac(n,n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE jacob
END INTERFACE

DO  ntest = 1, 10
  isum(1:6) = 0
  IF (ntest == 1) THEN
    OPEN (UNIT=nunit, FILE='T1.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 2) THEN
    OPEN (UNIT=nunit, FILE='T2.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 3) THEN
    OPEN (UNIT=nunit, FILE='T3.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 4) THEN
    OPEN (UNIT=nunit, FILE='T4.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 5) THEN
    OPEN (UNIT=nunit, FILE='T5.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 6) THEN
    OPEN (UNIT=nunit, FILE='T6.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 7) THEN
    OPEN (UNIT=nunit, FILE='T7.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 8) THEN
    OPEN (UNIT=nunit, FILE='T8.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 9) THEN
    OPEN (UNIT=nunit, FILE='T9.OUT', STATUS='UNKNOWN')
  ELSE IF (ntest == 10) THEN
    OPEN (UNIT=nunit, FILE='T10.OUT', STATUS='UNKNOWN')
  END IF

  WRITE (nunit,5000) ntest
  WRITE(*, 5000) ntest
  WRITE (nunit,5100)
  WRITE (nunit,5200)
  DO  nprob = 20, 20
    WRITE(*, 5300) nprob
    icount = 1
    DO  njupd = 0, 2
      DO  method = 1, 2
        IF (nprob == 1) THEN
          
!                    ROSENBROCK BANANA
          
          NAME = 'Rosenbrock Banana'
          ndim = 2
          xc(1) = -1.2_dp
          xc(2) = 1.0_dp
        ELSE IF (nprob == 2) THEN
          
!                    POWELL SINGULAR
          
          NAME = 'Powell Singular '
          ndim = 4
          xc(1) = 3.0_dp
          xc(2) = -1.0_dp
          xc(3) = 0.0_dp
          xc(4) = 1.0_dp
          
        ELSE IF (nprob == 3) THEN
          
!                    POWELL BADLY-SCALED
          
          NAME = 'Powell Badly-Scaled'
          ndim = 2
          xc(1) = 0.0_dp
          xc(2) = 1.0_dp
          
        ELSE IF (nprob == 4) THEN
          
!                    WOOD
          
          NAME = 'Wood '
          ndim = 4
          xc(1) = -3.0_dp
          xc(2) = -1.0_dp
          xc(3) = -3.0_dp
          xc(4) = -1.0_dp
          
        ELSE IF (nprob == 5) THEN
          
!                    HELICAL VALLEY
          
          NAME = 'Helical Valley '
          ndim = 3
          xc(1) = -1.0_dp
          xc(2) = 0.0_dp
          xc(3) = 0.0_dp
          
        ELSE IF (nprob == 6 .OR. nprob == 7) THEN
          
!                    WATSON
          
          IF (nprob == 6) THEN
            NAME = 'Watson 6-d'
            ndim = 6
          ELSE
            NAME = 'Watson 9-d'
            ndim = 9
          END IF
          DO  i = 1, ndim
            xc(i) = 0.0_dp
          END DO
          
        ELSE IF (nprob == 8 .OR. nprob == 9 .OR. nprob == 10) THEN
          
!                    CHEBYQUAD
          
          IF (nprob == 8) THEN
            NAME = 'Chebyquad 5-d '
            ndim = 5
          ELSE IF (nprob == 9) THEN
            NAME = 'Chebyquad 6-d '
            ndim = 6
          ELSE
            NAME = 'Chebyquad 7-d '
            ndim = 7
          END IF
          DO  i = 1, ndim
            xc(i) = i / (ndim + 1.0_dp)
          END DO
          
        ELSE IF (nprob == 11 .OR. nprob == 12 .OR. nprob == 13) THEN
          
!                    BROWN ALMOST-LINEAR
          
          IF (nprob == 11) THEN
            NAME = 'Brown Almost-Linear 10-d '
            ndim = 10
          ELSE IF (nprob == 12) THEN
            NAME = 'Brown Almost-Linear 30-d '
            ndim = 30
          ELSE
            NAME = 'Brown Almost-Linear 40-d '
            ndim = 40
          END IF
          DO  i = 1, ndim
            xc(i) = 0.5_dp
          END DO
          
        ELSE IF (nprob == 14 .OR. nprob == 15) THEN
          
!                    DISCRETE BOUNDARY VALUE
!                    DISCRETE INTEGRAL
          
          IF (nprob == 14) THEN
            NAME = 'Discrete Boundary Value '
          ELSE
            NAME = 'Discrete Integral '
          END IF
          ndim = 10
          hh = 1.0_dp / (ndim + 1.0_dp)
          DO  i = 1, ndim
            ti = i * hh
            xc(i) = ti * (ti - 1.0_dp)
          END DO
          
        ELSE IF (nprob == 16) THEN
          
!                    TRIGONOMETRIC
          
          NAME = 'Trigonometric '
          ndim = 10
          xc(1:ndim) = 0.1_dp
          
        ELSE IF (nprob == 17) THEN
          
!                    VARIABLY DIMENSIONED
          
          NAME = 'Variably Dimensioned '
          ndim = 10
          DO  i = 1, ndim
            xc(i) = 1.0_dp - 0.1_dp * i
          END DO
          
        ELSE IF (nprob == 18 .OR. nprob == 19) THEN
          
!                    BROYDEN TRIDIAGONAL
!                    BROYDEN BANDED
          
          IF (nprob == 18) THEN
            NAME = 'Broyden Tridiagonal '
          ELSE
            NAME = 'Broyden Banded '
          END IF
          ndim = 10
          xc(1:ndim) = -1.0_dp
          
        ELSE IF (nprob == 20 .OR. nprob == 21) THEN
          
!                    HIEBERT #1
          
          NAME = 'Hiebert #1 '
          ndim = 2
          IF (nprob == 20) THEN
            xc(1) = 0.0_dp
            xc(2) = 0.0_dp
          ELSE
            xc(1) = 10.0_dp
            xc(2) = 10.0_dp
          END IF
          
        ELSE IF (nprob >= 22.AND.nprob <= 25) THEN
          
!                    HIEBERT #2
          
          NAME = 'Hiebert #2 '
          ndim = 6
          IF (nprob /= 24) THEN
            DO  i = 1, ndim
              IF (nprob == 22) xc(i) = 0.0_dp
              IF (nprob == 23) xc(i) = 1.0_dp
              IF (nprob == 25) xc(i) = 10.0_dp
            END DO
          ELSE
            xc(1) = 1.0D-04
            xc(2) = 1.0D-03
            xc(3) = 0.0_dp
            xc(4) = 1.0D-04
            xc(5) = 55.0_dp
            xc(6) = 1.0D-04
          END IF
          
        ELSE IF (nprob >= 26) THEN
          
!                    HIEBERT #3
          
          NAME = 'Hiebert #3 '
          ndim = 10
          xc(1:ndim) = 0.0_dp
          IF (nprob == 26) THEN
            xc(1) = 1.0_dp
            xc(2) = 1.0_dp
            xc(3) = 10.0_dp
            xc(4) = 1.0_dp
            xc(5) = 1.0_dp
            xc(6) = 1.0_dp
            r = 10.0_dp
          ELSE IF (nprob == 27) THEN
            xc(1) = 2.0_dp
            xc(2) = 2.0_dp
            xc(3) = 10.0_dp
            xc(4) = 1.0_dp
            xc(5) = 1.0_dp
            xc(6) = 2.0_dp
            r = 10.0_dp
          ELSE IF (nprob == 28) THEN
            xc(1) = 2.0_dp
            xc(2) = 5.0_dp
            xc(3) = 40.0_dp
            xc(4) = 1.0_dp
            xc(10) = 5.0_dp
            r = 40.0_dp
          ELSE
            xc(1) = 1.0_dp
            xc(2) = 1.0_dp
            xc(3) = 20.0_dp
            xc(4) = 1.0_dp
            xc(10) = 1.0_dp
            r = 40.0_dp
          END IF
        END IF
        CALL setup(absnew, cauchy, deuflh, geoms, linesr, newton, overch,  &
                   acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit,  &
                   maxns, maxqns, minqns, n, narmij, niejev, njacch, output, &
                   qnupdm, stopcr, supprs, trupdm, alpha, confac, delta,  &
                   delfac, epsmch, etafac, fdtolj, ftol, lam0, mstpf,   &
                   nsttol, omega, ratiof, sigma, stptol, boundl, boundu,  &
                   scalef, scalex, help)
!     BUG FIX 7 MAR 1998
        jupdm = njupd
        output = 5
        WRITE (10, 5400) njupd, method
        IF (nprob == 20) THEN
          stopcr = 2
        ELSE IF (nprob >= 22 .AND. nprob <= 25) THEN
          nsttol = 1.0D-16
          stptol = 1.0D-16
          boundl(3) = 0.0_dp
          boundl(4) = 0.0_dp
          boundl(6) = 0.0_dp
          bypass = .true.
        ELSE IF (nprob >= 26) THEN
          boundl(1) = 0.0_dp
          boundl(2) = 0.0_dp
          boundl(3) = 0.0_dp
          boundl(4) = 0.0_dp
        END IF
        IF (method == 1) THEN
          linesr = .true.
        ELSE
          linesr = .false.
        END IF
        IF (ntest == 1) THEN
          newton = .true.
          GO TO 110
        ELSE IF (ntest == 2) THEN
          GO TO 110
        ELSE IF (ntest == 3) THEN
          itsclf = 1
          GO TO 110
        ELSE IF (ntest == 4) THEN
          acptcr = 1
          narmij = maxit
          deuflh = .false.
          GO TO 110
        ELSE IF (ntest == 5) THEN
          trupdm = 1
          deuflh = .false.
          GO TO 110
        ELSE IF (ntest == 6) THEN
          geoms = .false.
          GO TO 110
        ELSE IF (ntest == 7) THEN
          acptcr = 1
          GO TO 110
        ELSE IF (ntest == 8) THEN
          narmij = maxit
          GO TO 110
        ELSE IF (ntest == 9) THEN
          qnupdm = 0
        ELSE IF (ntest == 10) THEN
          narmij = maxit
          acptcr = 1
          GO TO 110
        END IF

        110 CALL nnes(absnew, cauchy, deuflh, geoms, linesr, newton, overch, &
                      acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit,  &
                      maxns, maxqns, mgll, minqns, ndim, narmij, niejev,   &
                      njacch, njetot, nunit, output, qnupdm, stopcr, supprs, &
                      trmcod, trupdm, alpha, confac, delta, delfac, epsmch,  &
                      etafac, fcnnew, fdtolj, ftol, lam0, mstpf, nsttol,   &
                      omega, ratiof, sigma, stptol, a, boundl, boundu, delf, &
                      fsave, ftrack, fvec, fvecc, h, hhpi, jac, plee, rdiag, &
                      s, sbar, scalef, scalex, sn, ssdhat, strack, vhat, xc, &
                      xplus, xsave, help, fcn, jacob)
        nfe(icount) = nfetot
        IF (trmcod > 0) isum(icount) = isum(icount) + nfetot
        icount = icount + 1
      END DO
    END DO

    WRITE(*, '(a, i4)') ' ICOUNT =', icount
    WRITE (nunit,5500) NAME, nfe(1:6)
  END DO
  WRITE (nunit,5600)
  WRITE (nunit,5700) isum(1:6)
END DO
STOP

5000 FORMAT (/t2, 'NTEST: ', i3/)
5100 FORMAT (/t30, 'NEWTON', t46, 'BROYDEN', t63, 'LEE AND LEE'/)
5200 FORMAT (t29, 'LS    TR         LS    TR           LS    TR'/)
5300 FORMAT (t2, 'NPROB: ', i3)
5400 FORMAT ('   njupd = ', i4, '   method = ', i4)
5500 FORMAT (t2, a, i4, '   ', i3, '        ', i3, '   ', i3, '          ', &
             i3, '   ', i3)
5600 FORMAT (t28, '---   ---        ---   ---          ---   ---')
5700 FORMAT (t27, i4, '  ', i4, '       ', i4, '  ', i4, '         ', i4,  &
             '  ', i4)
END PROGRAM mdr



SUBROUTINE fcn(overfl, n, fvec, xc)

USE probn
USE Equation_Solver, ONLY: nfetot
IMPLICIT NONE

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: fvec(n)
REAL (dp), INTENT(IN)   :: xc(n)

! COMMON /nnes_4/ nfetot
! COMMON /probn/ nprob
! COMMON /h3r/ r

REAL (dp)  :: t(0:9,9)
REAL (dp)  :: dfidxj, fi, h, k1, k2, k3, pi, prod, sum, sum2,  &
              t1, theta, ti, tj, tn, tot
INTEGER    :: i, j, k, r1, r2

overfl = .false.
nfetot = nfetot + 1
IF (nprob == 1) THEN
  
!        Rosenbrock Banana Function
  
  fvec(1) = 10.0_dp * (xc(2)-xc(1)**2)
  fvec(2) = 1.0_dp - xc(1)
ELSE IF (nprob == 2) THEN
  
!        Powell Singular Function
  
  fvec(1) = xc(1) + 10.0_dp * xc(2)
  fvec(2) = SQRT(5.0_dp) * (xc(3)-xc(4))
  fvec(3) = (xc(2)-2.0_dp*xc(3)) ** 2
  fvec(4) = SQRT(10.0_dp) * (xc(1)-xc(4)) ** 2
ELSE IF (nprob == 3) THEN
  
!        Powell Badly-Scaled Function
  
  fvec(1) = 1.0D04 * xc(1) * xc(2) - 1.0_dp
  fvec(2) = EXP(-xc(1)) + EXP(-xc(2)) - 1.0001_dp
ELSE IF (nprob == 4) THEN
  
!        Wood Function
  
  
  fvec(1) = -200.0_dp * xc(1) * (xc(2)-xc(1)**2) - (1.0_dp-xc(1))
  fvec(2) = 100.0_dp * (xc(2)-xc(1)**2) + 10.0_dp * (xc(2)+xc(4)-  &
            2.0_dp) + 0.1_dp * (xc(2)-xc(4))
  fvec(3) = -180.0_dp * xc(3) * (xc(4)-xc(3)**2) - (1.0_dp-xc(3))
  fvec(4) = 90.0_dp * (xc(4)-xc(3)**2) + 10.0_dp * (xc(2)+xc(4)-  &
      2.0_dp) - 0.1_dp * (xc(2)-xc(4))
ELSE IF (nprob == 5) THEN
  
!        Helical Valley Function
  
  pi = 3.141592653589793238_dp
  IF (xc(1) > 0.0_dp) THEN
    theta = (1.0_dp/(2.0_dp*pi)) * ATAN(xc(2)/xc(1))
  ELSE
    theta = (1.0_dp/(2.0_dp*pi)) * ATAN(xc(2)/xc(1)) + 0.5_dp
  END IF
  fvec(1) = 10.0_dp * (xc(3) - 10.0_dp*theta)
  fvec(2) = 10.0_dp * (SQRT((xc(1)**2 + xc(2)**2)) - 1.0_dp)
  fvec(3) = xc(3)
ELSE IF (nprob == 6 .OR. nprob == 7) THEN
  
!        Watson Function
  
  DO  j = 1, n
    fvec(j) = 0.0_dp
    DO  i = 1, 29
      ti = i / 29.0_dp
      sum = 0.0_dp
      DO  k = 1, n
        sum = sum + xc(k) * ti ** (k-1)
      END DO
      dfidxj = (j - 1) * ti ** (j-2) - 2.0_dp * sum * ti ** (j-1)
      sum = 0.0_dp
      DO  k = 2, n
        sum = sum + (k - 1) * xc(k) * ti ** (k-2)
      END DO
      sum2 = 0.0_dp
      DO  k = 1, n
        sum2 = sum2 + xc(k) * ti ** (k-1)
      END DO
      sum2 = sum2 ** 2 + 1.0_dp
      fi = sum - sum2
      fvec(j) = fvec(j) + dfidxj * fi
    END DO
  END DO
  fvec(1) = fvec(1) + xc(1) - 2.0_dp * xc(1) * (xc(2) - xc(1)**2 - 1.0_dp)
  fvec(2) = fvec(2) + xc(2) - xc(1) ** 2 - 1.0_dp
ELSE IF (nprob >= 8 .AND. nprob <= 10) THEN
  
!        Chebyquad Function
  
  DO  k = 1, n
    DO  j = 1, n
      t(0,j) = 1.0_dp
      t(1,j) = 2.0_dp * xc(j) - 1.0_dp
      DO  i = 2, n
        t(i,j) = (4.0_dp*xc(j) - 2.0_dp) * t(i-1,j) - t(i-2,j)
      END DO
    END DO
    fvec(k) = 0.0_dp
    DO  i = 1, n
      fvec(k) = fvec(k) + t(k,i)
    END DO
    fvec(k) = fvec(k) / n
    IF (MOD(k,2) == 0) THEN
      fvec(k) = fvec(k) + 1.0_dp / (k**2 - 1)
    END IF
  END DO
ELSE IF (nprob >= 11 .AND. nprob <= 13) THEN
  
!        Brown Almost-Linear Function
  
  sum = 0.0_dp
  prod = 1.0_dp
  DO  j = 1, n
    sum = sum + xc(j)
    prod = prod * xc(j)
  END DO
  sum = sum - (n + 1.0_dp)
  DO  i = 2, n
    fvec(i) = xc(i) + sum
  END DO
  fvec(1) = prod - 1.0_dp
ELSE IF (nprob == 14) THEN
  
!        Discrete Boundary Value Function
  
  h = 1.0_dp / (n + 1)
  DO  i = 2, n - 1
    ti = i * h
    fvec(i) = 2.0_dp * xc(i) - xc(i-1) - xc(i+1) + ((h**2)*(xc(i) +  &
              ti + 1.0_dp)**3) / 2.0_dp
  END DO
  t1 = h
  fvec(1) = 2.0_dp * xc(1) - xc(1+1) + ((h**2)*(xc(1) + t1 + 1.0_dp)**3) / 2.0_dp
  tn = n * h
  fvec(n) = 2.0_dp * xc(n) - xc(n-1) + ((h**2)*(xc(n) + tn + 1.0_dp)**3) / 2.0_dp
ELSE IF (nprob == 15) THEN
  
!        Discrete Integral Equation
  
  h = 1.0_dp / (n + 1.0_dp)
  DO  i = 1, n
    ti = i * h
    sum = 0.0_dp
    DO  j = 1, i
      tj = j * h
      sum = sum + tj * (xc(j) + tj + 1.0_dp) ** 3
    END DO
    sum = sum * (1.0_dp-ti)
    sum2 = 0.0_dp
    DO  j = i + 1, n
      tj = j * h
      sum2 = sum2 + (1.0_dp - tj) * (xc(j) + tj + 1.0_dp) ** 3
    END DO
    sum2 = sum2 * ti
    fvec(i) = xc(i) + (h/2.0_dp) * (sum + sum2)
  END DO
ELSE IF (nprob == 16) THEN
  
!        Trigonometric Function
  
  DO  i = 1, n
    sum = 0.0_dp
    DO  j = 1, n
      sum = sum + COS(xc(j))
    END DO
    fvec(i) = n - sum + i * (1.0_dp - COS(xc(i))) - SIN(xc(i))
  END DO
ELSE IF (nprob == 17) THEN
  
!        Variably Dimensioned Function
  
  sum = 0.0_dp
  DO  j = 1, n
    sum = sum + j * (xc(j) - 1.0_dp)
  END DO
  DO  i = 1, n
    fvec(i) = (xc(i) - 1.0_dp) + i * sum + 2.0_dp * i * sum ** 3
  END DO
ELSE IF (nprob == 18) THEN
  
!        Broyden Tridiagonal Function
  
  DO  i = 2, n - 1
    fvec(i) = (3.0_dp - 2.0_dp*xc(i)) * xc(i) - xc(i-1) - 2.0_dp * xc(i+1) + 1.0_dp
  END DO
  fvec(1) = (3.0_dp - 2.0_dp*xc(1)) * xc(1) - 2.0_dp * xc(2) + 1.0_dp
  fvec(n) = (3.0_dp - 2.0_dp*xc(n)) * xc(n) - xc(n-1) + 1.0_dp
ELSE IF (nprob == 19) THEN
  
!        Broyden Banded Function
  
  k1 = 2.0_dp
  k2 = 5.0_dp
  k3 = 1.0_dp
  r1 = 5
  r2 = 1
  
  DO  i = 1, n
    sum = 0.0_dp
    DO  j = MAX(1,i-r1), MIN(n,i+r2)
      IF (j /= i) sum = sum + xc(j) * (1.0_dp + xc(j))
    END DO
    fvec(i) = (k1 + k2*xc(i)**2) * xc(i) + 1.0_dp - k3 * sum
  END DO
ELSE IF (nprob == 20 .OR. nprob == 21) THEN
  
!        Hiebert's 1st Chemical Engineering Problem
  
  fvec(1) = xc(2) - 10.0_dp
  fvec(2) = xc(1) * xc(2) - 50000.0_dp
ELSE IF (nprob >= 22 .AND. nprob <= 25) THEN
  
!        Hiebert's 2nd Chemical Engineering Problem
  
  fvec(1) = xc(1) + xc(2) + xc(4) - 0.001_dp
  fvec(2) = xc(5) + xc(6) - 55.0_dp
  fvec(3) = xc(1) + xc(2) + xc(3) + 2.0_dp * xc(5) + xc(6) - 110.001_dp
  fvec(4) = xc(1) - 0.1_dp * xc(2)
  fvec(5) = xc(1) - 1.0D04 * xc(3) * xc(4)
  fvec(6) = xc(5) - 5.5D15 * xc(3) * xc(6)
ELSE IF (nprob >= 26) THEN
  
!        Hiebert's 3rd Chemical Reaction Problem (uncorrected)
  
  tot = 0.0_dp
  DO  i = 1, n
    tot = tot + xc(i)
  END DO
  IF (tot < 0.0_dp) THEN
    overfl = .true.
    RETURN
  END IF
  fvec(1) = xc(1) + xc(4) - 3.0_dp
  fvec(2) = 2.0_dp * xc(1) + xc(2) + xc(4) + xc(7) + xc(8) +  &
            xc(9) + 2.0_dp * xc(10) - r
  fvec(3) = 2.0_dp * xc(2) + 2.0_dp * xc(5) + xc(6) + xc(7) - 8.0_dp
  fvec(4) = 2.0_dp * xc(3) + xc(9) - 4.0_dp * r
  fvec(5) = xc(1) * xc(5) - 1.93D-01 * xc(2) * xc(4)
  fvec(6) = xc(6) * SQRT(xc(2)) - 2.597D-03 * SQRT(xc(2)*xc(4)*tot)
  fvec(7) = xc(7) * SQRT(xc(4)) - 3.448D-03 * SQRT(xc(1)*xc(2)*tot)
  fvec(8) = xc(8) * xc(4) - 1.799D-05 * xc(2) * tot
  fvec(9) = xc(9) * xc(4) - 2.155D-04 * xc(1) * SQRT(xc(3)*tot)
  fvec(10) = xc(10) * xc(4) ** 2 - 3.846D-05 * xc(4) ** 2 * tot
END IF
RETURN
END SUBROUTINE fcn



SUBROUTINE jacob(overfl, n, jac, xc)

USE probn
IMPLICIT NONE

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: jac(n,n)
REAL (dp), INTENT(IN)   :: xc(n)


overfl = .false.
jac = 0.0_dp
!*******************************************************************
!     *
!     INSERT ANALYTICAL JACOBIAN IF DESIRED                        *
!     *
!*******************************************************************
RETURN
END SUBROUTINE jacob

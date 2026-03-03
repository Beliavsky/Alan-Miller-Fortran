MODULE Nonlin_Conf_Regions

CONTAINS


SUBROUTINE halprn(param, np, ncols, nobs, np1, np2, p1lo, p1hi, p2lo, p2hi,  &
                  conf, nconf, f, grid, dim1, ngrid1, ngrid2, deriv,    &
                  lindep, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-04-01  Time: 15:54:58

!   ALGORITHM AS290  APPL. STATIST. (1994) VOL. 43, NO. 1

!   Generate a rectangular 2-D grid of variance ratios from which to plot
!   confidence regions for two parameters using Halperin's method.
!   If the model is linear in all parameters other than the two selected,
!   the confidence regions are exact; otherwise they are approximate and the
!   user should test the sensitivity of the confidence regions to variation
!   in the other parameters.

!   Auxiliary routines required: INCLUD, SS, TOLSET and SING, from AS 274
!   and a routine DERIV supplied by the user to calculate first derivatives.

! Arguments:
! PARAM(NP)  Parameter values
! NP         Number of parameters
! NCOLS      Number of first derivatives, usually the same as NP
! NOBS       Number of observations (needed to calculate deg. of freedom
! NP1, NP2   The numbers (positions) of the two parameters for which
!            the confidence regions are needed
! P1LO, P1HI Upper and lower limits for parameter number NP1
! P2LO, P2HI Same for the other parameter
! CONF(NCONF) The confidence levels in percent
! NCONF      The number of confidence levels
! F          (Output) F ratios corresponding to the confidence levels
! GRID(DIM1, NGRID2) (Output) 2-D array of variance ratios on a regular
!            grid of values of the two parameters
! DIM1       1st dimension of array GRID in the calling program
! NGRID1, NGRID2 The numbers of grid points required
! DERIV      The name of the user's routine to calculate 1st derivatives
! LINDEP(NCOLS) (Output) Logical array.   LINDEP(I) = .TRUE. if the Ith
!            column of derivatives is linearly related to previous columns
! IFAULT     (Output) Error indicator
!            = 0 if no errors detected
!            = 1 if NCOLS < 2
!            = 2 if NOBS < NCOLS
!            = 4 if NP1 outside the range 1 to NP
!            = 8 ditto for NP2
!            = 16 if NP1 = NP2
!            = 32 if NGRID1 < 2 or NGRID2 < 2 or NGRID1 > DIM1
!            = -I if the Ith confidence level is not between 0 and 100
! N.B. IFAULT may equal a sum of some of the above values if there is more
!      than one error.

USE lsq, ONLY: includ, rss, ss, sserr, startup, tolset, sing
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN OUT)  :: param(:)
INTEGER, INTENT(IN)        :: np
INTEGER, INTENT(IN)        :: ncols
INTEGER, INTENT(IN)        :: nobs
INTEGER, INTENT(IN)        :: np1
INTEGER, INTENT(IN)        :: np2
REAL (dp), INTENT(IN)      :: p1lo
REAL (dp), INTENT(IN)      :: p1hi
REAL (dp), INTENT(IN)      :: p2lo
REAL (dp), INTENT(IN)      :: p2hi
REAL (dp), INTENT(IN)      :: conf(:)
INTEGER, INTENT(IN)        :: nconf
REAL (dp), INTENT(OUT)     :: f(nconf)
INTEGER, INTENT(IN)        :: dim1
INTEGER, INTENT(IN)        :: ngrid1
INTEGER, INTENT(IN)        :: ngrid2
REAL (dp), INTENT(OUT)     :: grid(dim1, ngrid2)
LOGICAL, INTENT(OUT)       :: lindep(ncols)
INTEGER, INTENT(OUT)       :: ifault

! EXTERNAL deriv
INTERFACE
  SUBROUTINE deriv(param, np, ncols, iobs, grad, resid, wt)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)     :: param(:)
    INTEGER, INTENT(IN)       :: np, ncols, iobs
    REAL (dp), INTENT(OUT)    :: grad(:)
    REAL (dp), INTENT(OUT)    :: resid
    REAL (dp), INTENT(IN)     :: wt
  END SUBROUTINE deriv
END INTERFACE

!     Local variables

INTEGER    :: i, ndf, i1, i2, iobs
REAL (dp)  :: q, resid, step1, step2, p1save, p2save, ssq1
REAL (dp), ALLOCATABLE  :: grad(:)
REAL (dp), PARAMETER    :: zero = 0.0_dp, one = 1.0_dp, hundrd = 100.0_dp,  &
                           half = 0.5_dp, two = 2.0_dp, wt = 1.0_dp

!     Some checks

ifault = 0
IF (ncols < 2) ifault = 1
IF (nobs <= ncols) ifault = ifault + 2
IF (np1 < 1 .OR. np1 > np) ifault = ifault + 4
IF (np2 < 1 .OR. np2 > np) ifault = ifault + 8
IF (np1 == np2) ifault = ifault + 16
IF (ngrid1 < 2 .OR. ngrid2 < 2 .OR. ngrid1 > dim1) ifault = ifault + 32
IF (ifault > 0) RETURN

!     Calculate step sizes.

step1 = (p1hi - p1lo) / (ngrid1 - 1)
step2 = (p2hi - p2lo) / (ngrid2 - 1)
p1save = param(np1)
p2save = param(np2)
!----------------------------------------------------------------------

!     Start of cycle through the grid.

param(np1) = p1lo
ALLOCATE( grad(ncols) )
DO  i1 = 1, ngrid1
  param(np2) = p2lo
  DO  i2 = 1, ngrid2
    
!     Initialize orthogonal reduction.
    
    CALL startup(ncols, .FALSE.)
    DO  iobs = 1, nobs
      CALL deriv(param, np, ncols, iobs, grad, resid, wt)
      
!     Re-order the derivatives so that those for parameters NP1 & NP2
!     are at the end.
      
      IF (np1 < ncols-1) THEN
        IF (np2 /= ncols-1) THEN
          q = grad(ncols-1)
          grad(ncols-1) = grad(np1)
          grad(np1) = q
        ELSE
          q = grad(ncols)
          grad(ncols) = grad(np1)
          grad(np1) = q
        END IF
      END IF
      IF (np2 < ncols-1) THEN
        q = grad(ncols)
        grad(ncols) = grad(np2)
        grad(np2) = q
      END IF
      
!     Update QR factorization.
      
      CALL includ(wt, grad, resid)
    END DO
    
!----------------------------------------------------------------------
    
!     Test for singularities
    
    CALL tolset( 1.0D-8 )
    CALL sing( lindep, ifault)
    
!     If IFAULT < 0, it contains the rank deficiency.
    
    ndf = nobs - ncols - ifault
    
!     Calculate variance ratio.
    
    CALL ss()
    IF (ncols > 2) THEN
      ssq1 = rss(ncols-2)
    ELSE
      ssq1 = rss(2)
    END IF
    grid(i1, i2) = half * (ssq1 - sserr) / (sserr / ndf)
    param(np2) = param(np2) + step2
  END DO
  param(np1) = param(np1) + step1
END DO

param(np1) = p1save
param(np2) = p2save

!     Convert confidence levels in % to F-values, if NCONF > 0.

DO  i = 1, nconf
  q = one - conf(i) / hundrd
  IF (q <= zero .OR. q > one) THEN
    ifault = -i
    RETURN
  END IF
  f(i) = half * ndf * (q**(-two/ndf) - one)
END DO

RETURN
END SUBROUTINE halprn

END MODULE Nonlin_Conf_Regions

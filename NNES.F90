MODULE Equation_Solver

! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-09  Time: 15:20:33

! COPYRIGHT 1991, BY R.S. BAIN
! Originally wrapped by bain@luther.che.wisc.edu on Tue Mar 30 16:27:48 1993

! As of this writing, we do not know how to contact Rod Bain, the
! author of NNES.
!
! -- David M. Gay (dmg@bell-labs.com)
!    Bell Labs, Murray Hill
!    8 February 1999

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         ten = 10.0_dp

! COMMON /nnes_1/ matsup
LOGICAL, SAVE  :: matsup

! COMMON /nnes_2/ wrnsup
LOGICAL, SAVE  :: wrnsup

! COMMON /nnes_3/ bypass
LOGICAL, SAVE  :: bypass

! COMMON /nnes_4/ nfetot
INTEGER, SAVE  :: nfetot

! COMMON /nnes_5/ smallb,bigb,smalls,bigs,bigr
REAL (dp), SAVE  :: smallb, bigb, smalls, bigs, bigr

! COMMON /nnes_6/ itnum,nfunc
INTEGER, SAVE  :: itnum, nfunc

! COMMON /ill/ newstm
! It appears that NEWSTM was never set in the F77 version.


CONTAINS


SUBROUTINE abmul (nradec, nraact, ncbdec, ncbact, ncdec, ncact,  &
                  amat, bmat, cmat, arow)
 
!
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION AB=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN, AROW, BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX

!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCBDEC IS THE NUMBER OF COLUMNS OF B DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION
!
!    MODIFICATION OF MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncbdec
INTEGER, INTENT(IN)     :: ncbact
INTEGER, INTENT(IN)     :: ncdec
INTEGER, INTENT(IN)     :: ncact
REAL (dp), INTENT(IN)   :: amat(nradec,ncdec)
REAL (dp), INTENT(IN)   :: bmat(ncdec,ncbdec)
REAL (dp), INTENT(OUT)  :: cmat(nradec,ncbdec)
REAL (dp), INTENT(OUT)  :: arow(ncdec)

! Local variables

INTEGER    :: i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, &
              ncc32r
REAL (dp)  :: SUM

!
!       FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
ncc32=ncact/32
ncc32r=ncact - 32*ncc32
ncc16=ncc32r/16
ncc16r=ncc32r - 16*ncc16
ncc8=ncc16r/8
ncc8r=ncc16r - 8*ncc8
ncc4=ncc8r/4
ncc4r=ncc8r - 4*ncc4
!
!       REASSIGN ROWS TO VECTOR AROW.
!
DO  i=1,nraact
  k=0
  IF (ncc32 > 0) THEN
    DO  kk=1,ncc32
      k=k+32
      arow(k-31)=amat(i,k-31)
      arow(k-30)=amat(i,k-30)
      arow(k-29)=amat(i,k-29)
      arow(k-28)=amat(i,k-28)
      arow(k-27)=amat(i,k-27)
      arow(k-26)=amat(i,k-26)
      arow(k-25)=amat(i,k-25)
      arow(k-24)=amat(i,k-24)
      arow(k-23)=amat(i,k-23)
      arow(k-22)=amat(i,k-22)
      arow(k-21)=amat(i,k-21)
      arow(k-20)=amat(i,k-20)
      arow(k-19)=amat(i,k-19)
      arow(k-18)=amat(i,k-18)
      arow(k-17)=amat(i,k-17)
      arow(k-16)=amat(i,k-16)
      arow(k-15)=amat(i,k-15)
      arow(k-14)=amat(i,k-14)
      arow(k-13)=amat(i,k-13)
      arow(k-12)=amat(i,k-12)
      arow(k-11)=amat(i,k-11)
      arow(k-10)=amat(i,k-10)
      arow(k-9)=amat(i,k-9)
      arow(k-8)=amat(i,k-8)
      arow(k-7)=amat(i,k-7)
      arow(k-6)=amat(i,k-6)
      arow(k-5)=amat(i,k-5)
      arow(k-4)=amat(i,k-4)
      arow(k-3)=amat(i,k-3)
      arow(k-2)=amat(i,k-2)
      arow(k-1)=amat(i,k-1)
      arow(k)=amat(i,k)
    END DO
  END IF
  IF (ncc16 > 0) THEN
    DO  kk=1,ncc16
      k=k+16
      arow(k-15)=amat(i,k-15)
      arow(k-14)=amat(i,k-14)
      arow(k-13)=amat(i,k-13)
      arow(k-12)=amat(i,k-12)
      arow(k-11)=amat(i,k-11)
      arow(k-10)=amat(i,k-10)
      arow(k-9)=amat(i,k-9)
      arow(k-8)=amat(i,k-8)
      arow(k-7)=amat(i,k-7)
      arow(k-6)=amat(i,k-6)
      arow(k-5)=amat(i,k-5)
      arow(k-4)=amat(i,k-4)
      arow(k-3)=amat(i,k-3)
      arow(k-2)=amat(i,k-2)
      arow(k-1)=amat(i,k-1)
      arow(k)=amat(i,k)
    END DO
  END IF
  IF (ncc8 > 0) THEN
    DO  kk=1,ncc8
      k=k+8
      arow(k-7)=amat(i,k-7)
      arow(k-6)=amat(i,k-6)
      arow(k-5)=amat(i,k-5)
      arow(k-4)=amat(i,k-4)
      arow(k-3)=amat(i,k-3)
      arow(k-2)=amat(i,k-2)
      arow(k-1)=amat(i,k-1)
      arow(k)=amat(i,k)
    END DO
  END IF
  IF (ncc4 > 0) THEN
    DO  kk=1,ncc4
      k=k+4
      arow(k-3)=amat(i,k-3)
      arow(k-2)=amat(i,k-2)
      arow(k-1)=amat(i,k-1)
      arow(k)=amat(i,k)
    END DO
  END IF
  IF (ncc4r > 0) THEN
    DO  kk=1,ncc4r
      k=k+1
      arow(k)=amat(i,k)
    END DO
  END IF
!
!          FIND ENTRY FOR MATRIX C USING COLUMN VECTOR AROW.
!
  DO  j=1,ncbact
    sum=zero
    k=0
    IF (ncc32 > 0) THEN
      DO  kk=1,ncc32
        k=k+32
        sum=sum + arow(k-31)*bmat(k-31,j) + arow(k-30)*bmat(k-30,j) +  &
            arow(k-29)*bmat(k-29,j) + arow(k-28)*bmat(k-28,j) +  &
            arow(k-27)*bmat(k-27,j) + arow(k-26)*bmat(k-26,j) +  &
            arow(k-25)*bmat(k-25,j) + arow(k-24)*bmat(k-24,j)
        sum=sum + arow(k-23)*bmat(k-23,j) + arow(k-22)*bmat(k-22,j) +   &
            arow(k-21)*bmat(k-21,j) + arow(k-20)*bmat(k-20,j) +  &
            arow(k-19)*bmat(k-19,j) + arow(k-18)*bmat(k-18,j) +  &
            arow(k-17)*bmat(k-17,j) + arow(k-16)*bmat(k-16,j)
        sum=sum + arow(k-15)*bmat(k-15,j) + arow(k-14)*bmat(k-14,j) +   &
            arow(k-13)*bmat(k-13,j) + arow(k-12)*bmat(k-12,j) +   &
            arow(k-11)*bmat(k-11,j) + arow(k-10)*bmat(k-10,j) +   &
            arow(k-9)*bmat(k-9,j) + arow(k-8)*bmat(k-8,j)
        sum=sum  +  arow(k-7)*bmat(k-7,j) + arow(k-6)*bmat(k-6,j) +   &
            arow(k-5)*bmat(k-5,j) + arow(k-4)*bmat(k-4,j) +   &
            arow(k-3)*bmat(k-3,j) + arow(k-2)*bmat(k-2,j) +   &
            arow(k-1)*bmat(k-1,j) + arow(k)*bmat(k,j)
      END DO
    END IF
    IF (ncc16 > 0) THEN
      DO  kk=1,ncc16
        k=k+16
        sum=sum + arow(k-15)*bmat(k-15,j) + arow(k-14)*bmat(k-14,j) +   &
            arow(k-13)*bmat(k-13,j) + arow(k-12)*bmat(k-12,j) +   &
            arow(k-11)*bmat(k-11,j) + arow(k-10)*bmat(k-10,j) +   &
            arow(k-9)*bmat(k-9,j) + arow(k-8)*bmat(k-8,j)
        sum=sum + arow(k-7)*bmat(k-7,j) + arow(k-6)*bmat(k-6,j) +   &
            arow(k-5)*bmat(k-5,j) + arow(k-4)*bmat(k-4,j) +   &
            arow(k-3)*bmat(k-3,j) + arow(k-2)*bmat(k-2,j) +   &
            arow(k-1)*bmat(k-1,j) + arow(k)*bmat(k,j)
      END DO
    END IF
    IF (ncc8 > 0) THEN
      DO  kk=1,ncc8
        k=k+8
        sum=sum + arow(k-7)*bmat(k-7,j) + arow(k-6)*bmat(k-6,j) +   &
            arow(k-5)*bmat(k-5,j) + arow(k-4)*bmat(k-4,j) +   &
            arow(k-3)*bmat(k-3,j) + arow(k-2)*bmat(k-2,j) +   &
            arow(k-1)*bmat(k-1,j) + arow(k)*bmat(k,j)
      END DO
    END IF
    IF (ncc4 > 0) THEN
      DO  kk=1,ncc4
        k=k+4
        sum=sum + arow(k-3)*bmat(k-3,j) + arow(k-2)*bmat(k-2,j) +   &
            arow(k-1)*bmat(k-1,j) + arow(k)*bmat(k,j)
      END DO
    END IF
    IF (ncc4r > 0) THEN
      DO  kk=1,ncc4r
        k=k+1
        sum=sum + arow(k)*bmat(k,j)
      END DO
    END IF
    cmat(i,j)=sum
  END DO
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ABMUL.
!
END SUBROUTINE abmul



SUBROUTINE ascalf (n, epsmch, fvecc, jac, scalef)
!
!    FEB. 13, 1991
!
!    THIS SUBROUTINE ESTABLISHES SCALING FACTORS FOR THE RESIDUAL VECTOR
!    IF FUNCTION ADAPTIVE SCALING IS CHOSEN USING INTEGER VARIABLE ITSCLF.
!
!    NOTE: IN QUASI-NEWTON METHODS THE SCALING FACTORS ARE
!          UPDATED ONLY WHEN THE JACOBIAN IS EVALUATED EXPLICITLY.
!
!    SCALING FACTORS ARE DETERMINED FROM THE INFINITY NORMS OF THE ROWS
!    OF THE JACOBIAN AND THE VALUES OF THE CURRENT FUNCTION VECTOR, FVECC.
!
!    A MINIMUM TOLERANCE ON THE SCALING FACTOR IS THE SQUARE
!    ROOT OF THE MACHINE PRECISION, SQRTEP.
!

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(OUT)     :: scalef(n)

! Local variables

REAL (dp)  :: amax, sqrtep
INTEGER    :: i, j

!
sqrtep=SQRT(epsmch)
!
!       I COUNTS THE ROWS.
!
DO  i=1,n
  amax=zero
!
!          FIND MAXIMUM ENTRY IN ROW I.
!
  DO  j=1,n
    amax=MAX(amax, ABS(jac(i,j)))
  END DO
!
  amax=MAX(amax, fvecc(i))
!
!          SET SCALING FACTOR TO A DEFAULT OF ONE IF ITH ROW IS ZEROS.
!
  IF (amax == zero) amax=one
  amax=MAX(amax, sqrtep)
  scalef(i)=one/amax
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ASCALF.
!
END SUBROUTINE ascalf



SUBROUTINE ascalx (n, epsmch, jac, scalex)
!
!    FEB. 13, 1991
!
!    THIS SUBROUTINE ESTABLISHES SCALING FACTORS FOR THE COMPONENET VECTOR
!    IF ADAPTIVE SCALING IS CHOSEN USING INTEGER ITSCLX.
!
!    NOTE: IN QUASI-NEWTON METHODS THE SCALING FACTORS ARE
!          UPDATED ONLY WHEN THE JACOBIAN IS EVALUATED EXPLICITLY.
!
!    SCALING FACTORS ARE DETERMINED FROM THE INFINITY NORMS
!    OF THE COLUMNS OF THE JACOBIAN.
!
!    A MINIMUM TOLERANCE ON THE SCALING FACTOR IS THE SQUARE
!    ROOT OF THE MACHINE PRECISION, SQRTEP.
!

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: epsmch
REAL (dp), INTENT(IN)   :: jac(n,n)
REAL (dp), INTENT(OUT)  :: scalex(n)

! Local variables

REAL (dp)  :: amax, sqrtep
INTEGER    :: i, j

!
sqrtep=SQRT(epsmch)
!
!       J COUNTS COLUMNS.
!
DO  j=1,n
  amax=zero
!
!          FIND MAXIMUM ENTRY IN JTH COLUMN.
!
  DO  i=1,n
    amax=MAX(amax, ABS(jac(i,j)))
  END DO
!
!          IF A COLUMN IS ALL ZEROS SET AMAX TO ONE.
!
  IF (amax == zero) amax=one
  scalex(j)=MAX(amax, sqrtep)
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ASCALX.
!
END SUBROUTINE ascalx



SUBROUTINE atamul (nradec, ncadec, nraact, ncaact, nrbdec, ncbdec, amat, bmat)
!
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^A=B
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NRADEC IS NUMBER OF ROWS IN A DECLARED
!    NCADEC IS NUMBER OF COLUMNS IN A DECLARED
!    NRAACT IS THE LIMIT FOR THE 1ST INDEX IN A
!    NCAACT IS THE LIMIT FOR THE 2ND INDEX IN A
!    NRBDEC IS NUMBER OF ROWS IN B DECLARED
!    NCBDEC IS NUMBER OF COLUMNS IN B DECLARED
!
!    MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES
!    MACKINNON, QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: ncadec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncaact
INTEGER, INTENT(IN)     :: nrbdec
INTEGER, INTENT(IN)     :: ncbdec
REAL (dp), INTENT(IN)   :: amat(nradec,ncadec)
REAL (dp), INTENT(OUT)  :: bmat(nrbdec,ncbdec)

! Local variables

REAL (dp)  :: SUM
INTEGER    :: i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r

!
!       FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
ncc32=nraact/32
ncc32r=nraact - 32*ncc32
ncc16=ncc32r/16
ncc16r=ncc32r - 16*ncc16
ncc8=ncc16r/8
ncc8r=ncc16r - 8*ncc8
ncc4=ncc8r/4
ncc4r=ncc8r - 4*ncc4
!
!       FIND ENTRY IN MATRIX B.
!
DO  i=1,ncaact
  DO  j=i,ncaact
    sum=zero
    k=0
    IF (ncc32 > 0) THEN
      DO  kk=1,ncc32
        k=k+32
        sum=sum + amat(k-31,i)*amat(k-31,j) + amat(k-30,i)*amat(k-30,j) +  &
            amat(k-29,i)*amat(k-29,j) + amat(k-28,i)*amat(k-28,j) +   &
            amat(k-27,i)*amat(k-27,j) + amat(k-26,i)*amat(k-26,j) +   &
            amat(k-25,i)*amat(k-25,j) + amat(k-24,i)*amat(k-24,j)
        sum=sum + amat(k-23,i)*amat(k-23,j) + amat(k-22,i)*amat(k-22,j) +  &
            amat(k-21,i)*amat(k-21,j) + amat(k-20,i)*amat(k-20,j) +   &
            amat(k-19,i)*amat(k-19,j) + amat(k-18,i)*amat(k-18,j) +   &
            amat(k-17,i)*amat(k-17,j) + amat(k-16,i)*amat(k-16,j)
        sum=sum + amat(k-15,i)*amat(k-15,j) + amat(k-14,i)*amat(k-14,j) +  &
            amat(k-13,i)*amat(k-13,j) + amat(k-12,i)*amat(k-12,j) +   &
            amat(k-11,i)*amat(k-11,j) + amat(k-10,i)*amat(k-10,j) +   &
            amat(k-9,i)*amat(k-9,j) + amat(k-8,i)*amat(k-8,j)
        sum=sum + amat(k-7,i)*amat(k-7,j) + amat(k-6,i)*amat(k-6,j) +   &
            amat(k-5,i)*amat(k-5,j) + amat(k-4,i)*amat(k-4,j) +   &
            amat(k-3,i)*amat(k-3,j) + amat(k-2,i)*amat(k-2,j) +   &
            amat(k-1,i)*amat(k-1,j) + amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc16 > 0) THEN
      DO  kk=1,ncc16
        k=k+16
        sum=sum + amat(k-15,i)*amat(k-15,j) + amat(k-14,i)*amat(k-14,j) +  &
            amat(k-13,i)*amat(k-13,j) + amat(k-12,i)*amat(k-12,j) +   &
            amat(k-11,i)*amat(k-11,j) + amat(k-10,i)*amat(k-10,j) +   &
            amat(k-9,i)*amat(k-9,j) + amat(k-8,i)*amat(k-8,j)
        sum=sum + amat(k-7,i)*amat(k-7,j) + amat(k-6,i)*amat(k-6,j) +   &
            amat(k-5,i)*amat(k-5,j) + amat(k-4,i)*amat(k-4,j) +   &
            amat(k-3,i)*amat(k-3,j) + amat(k-2,i)*amat(k-2,j) +   &
            amat(k-1,i)*amat(k-1,j) + amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc8 > 0) THEN
      DO  kk=1,ncc8
        k=k+8
        sum=sum + amat(k-7,i)*amat(k-7,j) + amat(k-6,i)*amat(k-6,j) +   &
            amat(k-5,i)*amat(k-5,j) + amat(k-4,i)*amat(k-4,j) +   &
            amat(k-3,i)*amat(k-3,j) + amat(k-2,i)*amat(k-2,j) +   &
            amat(k-1,i)*amat(k-1,j) + amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc4 > 0) THEN
      DO  kk=1,ncc4
        k=k+4
        sum=sum + amat(k-3,i)*amat(k-3,j) + amat(k-2,i)*amat(k-2,j) +   &
            amat(k-1,i)*amat(k-1,j) + amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc4r > 0) THEN
      DO  kk=1,ncc4r
        k=k+1
        sum=sum + amat(k,i)*amat(k,j)
      END DO
    END IF
    bmat(i,j)=sum
    IF (i /= j) bmat(j,i)=bmat(i,j)
  END DO
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ATAMUL.
!
END SUBROUTINE atamul



SUBROUTINE ataov (overfl, maxexp, n, nunit, output, a, b, scalef)
!
!    SEPT. 8, 1991
!
!    THIS SUBROUTINE FINDS THE PRODUCT OF THE TRANSPOSE OF THE MATRIX A
!    AND MATRIX A.  EACH ENTRY IS CHECKED BEFORE BEING ACCEPTED.
!    IF IT WOULD CAUSE AN OVERFLOW 10**MAXEXP IS INSERTED IN ITS PLACE.
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(IN)   :: a(n,n)
REAL (dp), INTENT(OUT)  :: b(n,n)
REAL (dp), INTENT(IN)   :: scalef(n)

! Local variables

REAL (dp)  :: eps, SUM
INTEGER    :: i, j, k
!
eps=ten**(-maxexp)
overfl=.false.
!
DO  i=1,n
  DO  j=i+1,n
    sum=zero
    DO  k=1,n
      IF (LOG10(ABS(a(k,i))+eps) + LOG10(ABS(a(k,j))+eps) +   &
                           two*LOG10(scalef(k)) > maxexp) THEN
        overfl=.true.
        b(i,j)=SIGN(ten**maxexp,a(k,i))*SIGN(one,a(k,j))
        IF (output > 2 .AND. (.NOT.wrnsup)) THEN
          WRITE (nunit,10)
          10 FORMAT ('   *', t74, '*')
          WRITE (nunit,20) b(i,j)
          20 FORMAT ('   *    WARNING: COMPONENT IN',   &
                     ' MATRIX-MATRIX PRODUCT SET TO ', g12.3, t74, '*')
        END IF
        GO TO 40
      END IF
      sum=sum + a(k,i)*a(k,j)*scalef(k)*scalef(k)
    END DO
    b(i,j)=sum
    b(j,i)=sum
  40 END DO

  sum=zero
  DO  k=1,n
    IF (two*(LOG10(ABS(a(k,i))+eps) + LOG10(scalef(k))) > maxexp) THEN
      overfl=.true.
      b(i,i)=ten**maxexp
      IF (output > 2 .AND. (.NOT.wrnsup)) THEN
        WRITE (nunit,10)
        WRITE (nunit,20) b(i,i)
      END IF
      GO TO 70
    END IF
    sum=sum + a(k,i)*a(k,i)*scalef(k)*scalef(k)
  END DO
  b(i,i)=sum
70 END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ATAOV.
!
END SUBROUTINE ataov



SUBROUTINE atbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact,  &
                   amat, bmat, cmat)
!
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^B=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NCADEC IS 2ND DIM. OF AMAT; NCAACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E.   NCADEC IS NUMBER OF COLUMNS OF A DECLARED
!           NCBDEC IS NUMBER OF COLUMNS OF B DECLARED
!           NCDEC  IS THE NUMBER OF ROWS IN BOTH A AND B DECLARED
!
!    MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES
!    MACKINNON, QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: ncadec
INTEGER, INTENT(IN)     :: ncaact
INTEGER, INTENT(IN)     :: ncbdec
INTEGER, INTENT(IN)     :: ncbact
INTEGER, INTENT(IN)     :: ncdec
INTEGER, INTENT(IN)     :: ncact
REAL (dp), INTENT(IN)   :: amat(ncdec,ncadec)
REAL (dp), INTENT(IN)   :: bmat(ncdec,ncbdec)
REAL (dp), INTENT(OUT)  :: cmat(ncadec,ncbdec)

! Local variables

INTEGER    :: i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, &
              ncc32r
REAL (dp)  :: SUM

!
!       FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
ncc32=ncact/32
ncc32r=ncact - 32*ncc32
ncc16=ncc32r/16
ncc16r=ncc32r - 16*ncc16
ncc8=ncc16r/8
ncc8r=ncc16r - 8*ncc8
ncc4=ncc8r/4
ncc4r=ncc8r - 4*ncc4
!
!       FIND ENTRY IN MATRIX C.
!
DO  i=1,ncaact
  DO  j=1,ncbact
    sum=zero
    k=0
    IF (ncc32 > 0) THEN
      DO  kk=1,ncc32
        k=k+32
        sum=sum + amat(k-31,i)*bmat(k-31,j) + amat(k-30,i)*bmat(k-30,j) +  &
            amat(k-29,i)*bmat(k-29,j) + amat(k-28,i)*bmat(k-28,j) + &
            amat(k-27,i)*bmat(k-27,j) + amat(k-26,i)*bmat(k-26,j) + &
            amat(k-25,i)*bmat(k-25,j) + amat(k-24,i)*bmat(k-24,j)
        sum=sum + amat(k-23,i)*bmat(k-23,j) + amat(k-22,i)*bmat(k-22,j) + &
            amat(k-21,i)*bmat(k-21,j) + amat(k-20,i)*bmat(k-20,j) +  &
            amat(k-19,i)*bmat(k-19,j) + amat(k-18,i)*bmat(k-18,j) +  &
            amat(k-17,i)*bmat(k-17,j) + amat(k-16,i)*bmat(k-16,j)
        sum=sum + amat(k-15,i)*bmat(k-15,j) + amat(k-14,i)*bmat(k-14,j) +   &
            amat(k-13,i)*bmat(k-13,j) + amat(k-12,i)*bmat(k-12,j) +  &
            amat(k-11,i)*bmat(k-11,j) + amat(k-10,i)*bmat(k-10,j) +   &
            amat(k-9,i)*bmat(k-9,j) + amat(k-8,i)*bmat(k-8,j)
        sum=sum + amat(k-7,i)*bmat(k-7,j) + amat(k-6,i)*bmat(k-6,j) +   &
            amat(k-5,i)*bmat(k-5,j) + amat(k-4,i)*bmat(k-4,j) +   &
            amat(k-3,i)*bmat(k-3,j) + amat(k-2,i)*bmat(k-2,j) +   &
            amat(k-1,i)*bmat(k-1,j) + amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc16 > 0) THEN
      DO  kk=1,ncc16
        k=k+16
        sum=sum + amat(k-15,i)*bmat(k-15,j) + amat(k-14,i)*bmat(k-14,j) +  &
            amat(k-13,i)*bmat(k-13,j) + amat(k-12,i)*bmat(k-12,j) +  &
            amat(k-11,i)*bmat(k-11,j) + amat(k-10,i)*bmat(k-10,j) +  &
            amat(k-9,i)*bmat(k-9,j) + amat(k-8,i)*bmat(k-8,j)
        sum=sum + amat(k-7,i)*bmat(k-7,j) + amat(k-6,i)*bmat(k-6,j) +  &
            amat(k-5,i)*bmat(k-5,j) + amat(k-4,i)*bmat(k-4,j) +   &
            amat(k-3,i)*bmat(k-3,j) + amat(k-2,i)*bmat(k-2,j) +   &
            amat(k-1,i)*bmat(k-1,j) + amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc8 > 0) THEN
      DO  kk=1,ncc8
        k=k+8
        sum=sum + amat(k-7,i)*bmat(k-7,j) + amat(k-6,i)*bmat(k-6,j) +  &
            amat(k-5,i)*bmat(k-5,j) + amat(k-4,i)*bmat(k-4,j) +  &
            amat(k-3,i)*bmat(k-3,j) + amat(k-2,i)*bmat(k-2,j) +  &
            amat(k-1,i)*bmat(k-1,j) + amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc4 > 0) THEN
      DO  kk=1,ncc4
        k=k+4
        sum=sum + amat(k-3,i)*bmat(k-3,j) + amat(k-2,i)*bmat(k-2,j) +  &
            amat(k-1,i)*bmat(k-1,j) + amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc4r > 0) THEN
      DO  kk=1,ncc4r
        k=k+1
        sum=sum + amat(k,i)*bmat(k,j)
      END DO
    END IF
    cmat(i,j)=sum
  END DO
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ATBMUL.
!
END SUBROUTINE atbmul



SUBROUTINE atvov (overfl, maxexp, n, nunit, output, amat, bvec, cvec)
!
!    FEB. 8 ,1991
!
!    THIS SUBROUTINE FINDS THE PRODUCT OF THE TRANSPOSE OF THE MATRIX A
!    AND THE VECTOR B WHERE EACH ENTRY IS CHECKED TO PREVENT OVERFLOWS.
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(IN)   :: amat(n,n)
REAL (dp), INTENT(IN)   :: bvec(n)
REAL (dp), INTENT(OUT)  :: cvec(n)

! Local variables

REAL (dp)  :: eps, SUM
INTEGER    :: i, j
!
eps=ten**(-maxexp)
overfl=.false.
!
DO  i=1,n
  sum=zero
  DO  j=1,n
    IF (LOG10(ABS(amat(j,i))+eps) + LOG10(ABS(bvec(j))+eps) > maxexp) THEN
      overfl=.true.
      cvec(i)=SIGN(ten**maxexp, amat(j,i))*SIGN(one,bvec(j))
      IF (output > 2 .AND. (.NOT.wrnsup)) THEN
        WRITE (nunit,10)
        10 FORMAT ('   *', t74, '*')
        WRITE (nunit,20) cvec(i)
        20 FORMAT ('   *    WARNING: COMPONENT IN',  &
                   ' MATRIX-VECTOR PRODUCT SET TO ', g12.3, t74, '*')
      END IF
      GO TO 40
    END IF
    sum=sum + amat(j,i)*bvec(j)
  END DO
  cvec(i)=sum
40 END DO
RETURN
!
!       LAST CARD OF SUBROUTINE ATVOV.
!
END SUBROUTINE atvov



SUBROUTINE avmul (nradec, nraact, ncdec, ncact, amat, bvec, cvec)
!
!    FEB. 8, 1991
!
!    MATRIX-VECTOR MULTIPLICATION AB=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BVEC; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION (COLUMNS OF A AND ROWS OF B)
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncdec
INTEGER, INTENT(IN)     :: ncact
REAL (dp), INTENT(IN)   :: amat(nradec,ncdec)
REAL (dp), INTENT(IN)   :: bvec(ncdec)
REAL (dp), INTENT(OUT)  :: cvec(nradec)

! Local variables

INTEGER    :: i, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, &
              ncc32r
REAL (dp)  :: SUM

!
!       FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
ncc32=ncact/32
ncc32r=ncact - 32*ncc32
ncc16=ncc32r/16
ncc16r=ncc32r - 16*ncc16
ncc8=ncc16r/8
ncc8r=ncc16r - 8*ncc8
ncc4=ncc8r/4
ncc4r=ncc8r - 4*ncc4
DO  i=1,nraact
!
!          FIND ENTRY FOR VECTOR C.
!
  sum=zero
  k=0
  IF (ncc32 > 0) THEN
    DO  kk=1,ncc32
      k=k+32
      sum=sum + amat(i,k-31)*bvec(k-31) + amat(i,k-30)*bvec(k-30) +   &
          amat(i,k-29)*bvec(k-29) + amat(i,k-28)*bvec(k-28) +   &
          amat(i,k-27)*bvec(k-27) + amat(i,k-26)*bvec(k-26) +   &
          amat(i,k-25)*bvec(k-25) + amat(i,k-24)*bvec(k-24)
      sum=sum + amat(i,k-23)*bvec(k-23) + amat(i,k-22)*bvec(k-22) +   &
          amat(i,k-21)*bvec(k-21) + amat(i,k-20)*bvec(k-20) +   &
          amat(i,k-19)*bvec(k-19) + amat(i,k-18)*bvec(k-18) +   &
          amat(i,k-17)*bvec(k-17) + amat(i,k-16)*bvec(k-16)
      sum=sum + amat(i,k-15)*bvec(k-15) + amat(i,k-14)*bvec(k-14) +   &
          amat(i,k-13)*bvec(k-13) + amat(i,k-12)*bvec(k-12) +   &
          amat(i,k-11)*bvec(k-11) + amat(i,k-10)*bvec(k-10) +   &
          amat(i,k-9)*bvec(k-9) + amat(i,k-8)*bvec(k-8)
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) +  &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc16 > 0) THEN
    DO  kk=1,ncc16
      k=k+16
      sum=sum + amat(i,k-15)*bvec(k-15) + amat(i,k-14)*bvec(k-14) +   &
          amat(i,k-13)*bvec(k-13) + amat(i,k-12)*bvec(k-12) +   &
          amat(i,k-11)*bvec(k-11) + amat(i,k-10)*bvec(k-10) +   &
          amat(i,k-9)*bvec(k-9) + amat(i,k-8)*bvec(k-8)
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) +  &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc8 > 0) THEN
    DO  kk=1,ncc8
      k=k+8
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) +  &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc4 > 0) THEN
    DO  kk=1,ncc4
      k=k+4
      sum=sum + amat(i,k-3)*bvec(k-3) + amat(i,k-2)*bvec(k-2) +   &
          amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc4r > 0) THEN
    DO  kk=1,ncc4r
      k=k+1
      sum=sum + amat(i,k)*bvec(k)
    END DO
  END IF
  cvec(i)=sum
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE AVMUL.
!
END SUBROUTINE avmul



SUBROUTINE bakdif (overfl, j, n, deltaj, tempj, fvec, fvecj1,  &
                   jacfdm, xc, fvecev)
!
!       FEB. 6, 1991
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: j
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: deltaj
REAL (dp), INTENT(IN)   :: tempj
REAL (dp), INTENT(IN)   :: fvec(n)
REAL (dp), INTENT(OUT)  :: fvecj1(n)
REAL (dp), INTENT(OUT)  :: jacfdm(n,n)
REAL (dp), INTENT(IN)   :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

deltaj=tempj - xc(j)
CALL fvecev (overfl, n, fvecj1, xc)
IF (.NOT.overfl) THEN
  jacfdm(1:n,j)=(fvec(1:n) - fvecj1(1:n))/deltaj
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE BAKDIF.
!
END SUBROUTINE bakdif



SUBROUTINE bnddif (overfl, j, n, epsmch, boundl, boundu, fvecc,  &
                   fvecj1, jacfdm, xc, fvecev)

! N.B. Argument WV3 has been removed.

!
!   FEB. 15, 1991
!
!   FINITE DIFFERENCE CALCULATION WHEN THE BOUNDS FOR COMPONENT J ARE SO CLOSE
!   THAT NEITHER A FORWARD NOR BACKWARD DIFFERENCE CAN BE PERFORMED.
!

LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: j
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN)      :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: fvecc(n)
REAL (dp), INTENT(OUT)     :: fvecj1(n)
REAL (dp), INTENT(OUT)     :: jacfdm(n,n)
REAL (dp), INTENT(OUT)     :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

INTEGER    :: i
REAL (dp)  :: eps3q, wv3(n)
!
eps3q=epsmch**0.75
!
!       STORE CURRENT
wv3(1:n)=fvecc(1:n)
xc(j)=boundu(j)
CALL fvecev (overfl, n, fvecj1, xc)
IF (.NOT.overfl) THEN
  xc(j)=boundl(j)
  CALL fvecev (overfl, n, fvecc, xc)
  IF (.NOT.overfl) THEN
    DO  i=1,n
!
!                ENSURE THAT THE JACOBIAN CALCULATION ISN'T JUST NOISE.
!
      IF (fvecj1(i)-fvecc(i) > eps3q) THEN
        jacfdm(i,j)=(fvecj1(i)-fvecc(i))/(boundu(j)-boundl(j))
      ELSE
        jacfdm(i,j)=zero
      END IF
    END DO
  END IF
END IF
fvecc(1:n)=wv3(1:n)

RETURN
!
!       LAST CARD OF SUBROUTINE BNDDIF.
!
END SUBROUTINE bnddif



SUBROUTINE broyfa (overch, overfl, sclfch, sclxch, maxexp, n, nunit, output, &
                   epsmch, a, delf, fvec, fvecc, jac, rdiag, s, scalef,   &
                   scalex, t, w, xc, xplus)
!
!    FEB. 23, 1992
!
!    THE BROYDEN QUASI-NEWTON METHOD IS APPLIED TO THE FACTORED
!    FORM OF THE JACOBIAN.
!
!    NOTE: T AND W ARE TEMPORARY WORKING VECTORS ONLY.
!
!    THE UPDATE OCCURS ONLY IF A SIGNIFICANT CHANGE IN THE JACOBIAN
!    WOULD RESULT,  I.E., NOT ALL THE VALUES IN VECTOR W ARE LESS THAN
!    THE THRESHOLD IN MAGNITUDE.  IF THE VECTOR W IS ESSENTIALLY ZERO THEN
!    THE LOGICAL VARIABLE SKIPUP REMAINS SET AT TRUE.
!

LOGICAL, INTENT(IN OUT)    :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(IN)        :: sclfch
LOGICAL, INTENT(IN)        :: sclxch
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: rdiag(n)
REAL (dp), INTENT(OUT)     :: s(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: t(n)
REAL (dp), INTENT(OUT)     :: w(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN)      :: xplus(n)

! Local variables

REAL (dp)  :: denom, eps, SUM, sqrtep
LOGICAL    :: skipup
INTEGER    :: i

!
overfl=.false.
eps=ten**(-maxexp)
sqrtep=SQRT(epsmch)
!
DO  i=1,n
  a(i,i)=rdiag(i)
  s(i)=xplus(i) - xc(i)
END DO
!
!       R IS NOW IN THE UPPER TRIANGLE OF A.
!
skipup=.true.
!
!       THE BROYDEN UPDATE IS CONDENSED INTO THE FORM
!
!       A(NEW) = A(OLD) + T S^
!
!       THE PRODUCT A*S IS FORMED IN TWO STAGES AS R IS IN THE UPPER
!       TRIANGLE OF MATRIX A AND Q^ IS IN JAC.
!
!       FIRST MULTIPLY R*S (A IS CONSIDERED UPPER TRIANGULAR)
!
CALL uvmul (n, n, n, n, a, s, t)
!
!       NOTE: THIS T IS TEMPORARY - IT IS THE T FROM BELOW WHICH
!             IS SENT TO SUBROUTINE QRUPDA.
!
DO  i=1,n
  CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
               sum, jac(1,i), t)
  w(i)=scalef(i)*(fvec(i)-fvecc(i)) - sum
!
!          TEST TO ENSURE VECTOR W IS NONZERO.  ANY VALUE GREATER
!          THAN THE THRESHOLD WILL SET SKIPUP TO FALSE.
!
  IF (ABS(w(i)) > sqrtep*scalef(i)*(ABS(fvec(i)) + ABS(fvecc(i)))) THEN
    skipup=.false.
  ELSE
    w(i)=zero
  END IF
END DO
!
!       IF W(I)=0 FOR ALL I THEN THE UPDATE IS SKIPPED.
!
IF (.NOT.skipup) THEN
!
!          T=Q^W; Q^ IS IN JAC.
!
  CALL avmul (n, n, n, n, jac, w, t)
  IF (sclxch) THEN
    w(1:n)=s(1:n)*scalex(1:n)
  ELSE
    CALL matcop (n, n, 1, 1, n, 1, s, w)
  END IF
  CALL twonrm (overfl, maxexp, n, epsmch, denom, w)
!
!          IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN.
!
  IF (overfl .OR. LOG10(denom+eps) > maxexp/2) THEN
    IF (output > 3) THEN
      WRITE (nunit,40)
      40 FORMAT ('   *', t74, '*')
      WRITE (nunit,50)
      50 FORMAT ('   *    WARNING: JACOBIAN NOT UPDATED',   &
                 ' TO AVOID OVERFLOW IN DENOMINATOR OF', t74, '*')
      WRITE (nunit,60)
      60 FORMAT ('   *    BROYDEN UPDATE', t74, '*')
    END IF
    RETURN
  ELSE
    denom=denom*denom
  END IF
!
!       IF DENOM IS ZERO AVOID DIVIDE BY ZERO AND CONTINUE WITH SAME JACOBIAN.
!
  IF (denom == zero) RETURN
!
!       THE SCALED VERSION OF S REPLACES THE ORIGINAL BEFORE
!       BEING SENT TO QRUPDA.
!
  s(1:n)=s(1:n)*scalex(1:n)*scalex(1:n)/denom
!
!       UPDATE THE QR DECOMPOSITION USING A SERIES OF GIVENS ROTATIONS.
!
  CALL qrupda (overfl, maxexp, n, epsmch, a, jac, t, s)
!
!          RESET RDIAG AS DIAGONAL OF CURRENT R WHICH IS IN
!          THE UPPER TRIANGE OF A.
!
  DO  i=1,n
    rdiag(i)=a(i,i)
  END DO
END IF
!
!       UPDATE THE GRADIENT VECTOR, DELF.  THE NEW Q^ IS IN JAC.
!
!       DELF = (QR)^F = R^Q^F = R^JAC F
!
IF (sclfch) THEN
  w(1:n)=fvec(1:n)*scalef(1:n)
ELSE
  CALL matcop (n, n, 1, 1, n, 1, fvec, w)
END IF
CALL avmul (n, n, n, n, jac, w, t)
CALL utbmul (n, n, 1, 1, n, n, a, t, delf)
RETURN
!
!       LAST CARD OF SUBROUTINE BROYFA.
!
END SUBROUTINE broyfa



SUBROUTINE broyun (overfl, maxexp, n, nunit, output, epsmch, fvec,  &
                   fvecc, jac, scalex, xc, xplus)

! N.B. Argument WV1 has been removed.

!
!    FEB. 23, 1992
!
!    UPDATE THE JACOBIAN USING BROYDEN'S METHOD.
!

LOGICAL, INTENT(IN OUT)    :: overfl
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN)      :: xplus(n)

! Local variables

INTEGER    :: i, j
REAL (dp)  :: denom, eps, sqrtep, SUM, tempi, wv1(n)
!
eps=ten**(-maxexp)
sqrtep=SQRT(epsmch)
!
wv1(1:n)=(xplus(1:n) - xc(1:n))*scalex(1:n)
CALL twonrm (overfl, maxexp, n, epsmch, denom, wv1)
!
!       IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN.
!
IF (overfl .OR. LOG10(denom+eps) > maxexp/2) THEN
  IF (output > 3) THEN
    WRITE (nunit,20)
    20 FORMAT ('   *', t74, '*')
    WRITE (nunit,30)
    30 FORMAT ('   *    WARNING: JACOBIAN NOT UPDATED',   &
               ' TO AVOID OVERFLOW IN DENOMINATOR OF', t74, '*')
    WRITE (nunit,40)
    40 FORMAT ('   *    BROYDEN UPDATE', t74, '*')
  END IF
  RETURN
ELSE
  denom=denom*denom
END IF
!
!       IF DENOM IS ZERO, AVOID OVERFLOW, CONTINUE WITH SAME JACOBIAN.
!
IF (denom == zero) RETURN
!
!       UPDATE JACOBIAN BY ROWS.
!
DO  i=1,n
  sum = DOT_PRODUCT( jac(i,1:n), (xplus(1:n) - xc(1:n)) )
  tempi=fvec(i) - fvecc(i) - sum
!
!          CHECK TO ENSURE THAT SOME MEANINGFUL CHANGE IS BEING MADE
!          TO THE APPROXIMATE JACOBIAN; IF NOT, SKIP UPDATING ROW I.
!
  IF (ABS(tempi) >= sqrtep*(ABS(fvec(i)) + ABS(fvecc(i)))) THEN
    tempi=tempi/denom
    DO  j=1,n
      jac(i,j)=jac(i,j) + tempi*(xplus(j)-xc(j))*scalex(j)*scalex(j)
    END DO
  END IF
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE BROYUN.
!
END SUBROUTINE broyun



SUBROUTINE cholde (n, maxadd, maxffl, sqrtep, h, l)
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS THE CHOLESKY DECOMPOSITION OF THE MATRIX, H,
!    AND RETURNS IT IN THE LOWER TRIANGLE OF MATRIX, L.
!

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: maxadd
REAL (dp), INTENT(IN OUT)  :: maxffl
REAL (dp), INTENT(IN)      :: sqrtep
REAL (dp), INTENT(IN)      :: h(n,n)
REAL (dp), INTENT(IN OUT)  :: l(n,n)

! Local variables

REAL (dp)  :: minl, minl2, minljj, SUM
INTEGER    :: i, j

!
minl=SQRT(sqrtep)*maxffl
!
!       MAXFFL EQUALS 0 WHEN THE MATRIX IS KNOWN TO BE POSITIVE DEFINITE.
!
IF (maxffl == zero) THEN
!
!          FIND SQUARE ROOT OF LARGEST MAGNITUDE DIAGONAL ELEMENT
!          AND SET MINL2.
!
  DO  i=1,n
    maxffl=MAX(maxffl, ABS(h(i,i)))
  END DO
  maxffl=SQRT(maxffl)
  minl2=sqrtep*maxffl
END IF
!
!       MAXADD CONTAINS THE MAXIMUM AMOUNT WHICH IS IMPLICITLY ADDED
!       TO ANY DIAGONAL ELEMENT OF MATRIX H.
!
maxadd=zero
DO  j=1,n
  sum = DOT_PRODUCT( l(j,1:j-1), l(j,1:j-1) )
  l(j,j)=h(j,j)-sum
  minljj=zero
  DO  i=j+1,n
    sum = DOT_PRODUCT( l(i,1:j-1), l(j,1:j-1) )
    l(i,j)=h(j,i) - sum
    minljj=MAX(minljj, ABS(l(i,j)))
  END DO
  minljj=MAX(minljj/maxffl, minl)
  IF (l(j,j) > minljj*minljj) THEN
!
!             NORMAL CHOLESKY DECOMPOSITION.
!
    l(j,j)=SQRT(l(j,j))
  ELSE
!
!             IMPLICITLY PERTURB DIAGONAL OF H.
!
    IF (minljj < minl2) THEN
      minljj=minl2
    END IF
    maxadd=MAX(maxadd, minljj*minljj-l(j,j))
    l(j,j)=minljj
  END IF
  l(j+1:n,j)=l(j+1:n,j)/l(j,j)
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE CHOLDE.
!
END SUBROUTINE cholde



SUBROUTINE chsolv (overch, overfl, maxexp, n, nunit, output, l, rhs, s)

! N.B. Argument WV2 has been removed.
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE USES FORWARD/BACKWARD SUBSTITUTION TO SOLVE THE
!    SYSTEM OF LINEAR EQUATIONS:
!
!         (LL^)S=RHS
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(IN OUT)    :: overfl
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: l(n,n)
REAL (dp), INTENT(IN)      :: rhs(n)
REAL (dp), INTENT(IN OUT)  :: s(n)

! Local array

REAL (dp)  :: wv2(n)
!
CALL lsolv (overch, overfl, maxexp, n, nunit, output, l, wv2, rhs)
CALL ltsolv (overch, overfl, maxexp, n, nunit, output, l, s, wv2)
!
RETURN
!
!       LAST CARD OF SUBROUTINE CHSOLV.
!
END SUBROUTINE chsolv



SUBROUTINE condno (overch, overfl, maxexp, n, nunit, output, connum, a, rdiag)

! N.B. Arguments P, PM & Q have been removed.
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE ESTIMATES THE CONDITION NUMBER OF A QR-DECOMPOSED
!    MATRIX USING THE METHOD OF CLINE, MOLER, STEWART AND WILKINSON
!    (SIAM J. N.A. 16 P368 (1979) ).
!
!    IF A POTENTIAL OVERFLOW IS DETECTED AT ANY POINT THEN A CONDITION NUMBER
!    EQUIVALENT TO THAT OF A SINGULAR MATRIX IS ASSIGNED BY THE CALLING
!    SUBROUTINE.
!

LOGICAL, INTENT(IN)     :: overch
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(OUT)  :: connum
REAL (dp), INTENT(IN)   :: a(n,n)
REAL (dp), INTENT(IN)   :: rdiag(n)

! Local variables

INTEGER    :: i, j
REAL (dp)  :: eps, qm, qnorm, qp, temp, tempm
REAL (dp)  :: p(n), pm(n), q(n)
!
overfl=.false.
eps=ten**(-maxexp)
!
connum=ABS(rdiag(1))
DO  j=2,n
  temp=zero
  DO  i=1,j-1
    IF (overch) THEN
      IF (ABS(a(i,j)) > ten**(maxexp-1)) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    temp=temp + ABS(a(i,j))
  END DO
  temp=temp + ABS(rdiag(j))
  connum=MAX(connum,temp)
END DO
q(1)=one/rdiag(1)
DO  i=2,n
  IF (overch) THEN
    IF (LOG10(ABS(q(1))+eps) + LOG10(ABS(a(1,i)) + eps) > maxexp) THEN
      overfl=.true.
      RETURN
    END IF
  END IF
  p(i)=a(1,i)*q(1)
END DO
DO  j=2,n
  IF (overch) THEN
    IF (LOG10(ABS(p(j))+eps) - LOG10(ABS(rdiag(j)) + eps) > maxexp) THEN
      overfl=.true.
      RETURN
    END IF
  END IF
  qp=(one-p(j))/rdiag(j)
  qm=(-one-p(j))/rdiag(j)
  temp=ABS(qp)
  tempm=ABS(qm)
  DO  i=j+1,n
    IF (overch) THEN
      IF (LOG10(ABS(a(j,i))+eps) + LOG10(ABS(qm)+eps) > maxexp) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    pm(i)=p(i) + a(j,i)*qm
    IF (overch) THEN
      IF (LOG10(ABS(pm(i))+eps) - LOG10(ABS(rdiag(i))+eps) > maxexp) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    tempm=tempm + (ABS(pm(i))/ABS(rdiag(i)))
    IF (overch) THEN
      IF (tempm > ten**(maxexp-1)) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    IF (overch) THEN
      IF (LOG10(ABS(a(j,i))+eps) + LOG10(ABS(qp)+eps) > maxexp) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    p(i)=p(i) + a(j,i)*qp
    IF (overch) THEN
      IF (LOG10(ABS(p(i))+eps) - LOG10(ABS(rdiag(i))+eps) > maxexp) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
    temp=temp + (ABS(p(i))/ABS(rdiag(i)))
    IF (overch) THEN
      IF (temp > ten**(maxexp-1)) THEN
        overfl=.true.
        RETURN
      END IF
    END IF
  END DO
  IF (temp >= tempm) THEN
    q(j)=qp
  ELSE
    q(j)=qm
    p(j+1:n)=pm(j+1:n)
  END IF
END DO
qnorm=zero
DO  j=1,n
  qnorm=qnorm + ABS(q(j))
  IF (overch) THEN
    IF (qnorm > ten**(maxexp-1)) THEN
      overfl=.true.
      RETURN
    END IF
  END IF
END DO
IF (LOG10(connum) - LOG10(qnorm) > maxexp) THEN
  overfl=.true.
  RETURN
END IF
connum=connum/qnorm
CALL rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, q)
IF (overfl) RETURN
qnorm=zero
DO  j=1,n
  qnorm=qnorm + ABS(q(j))
  IF (overch) THEN
    IF (qnorm > ten**(maxexp-1)) THEN
      overfl=.true.
      RETURN
    END IF
  END IF
END DO
connum=connum*qnorm
RETURN
!
!       LAST CARD OF SUBROUTINE CONDNO.
!
END SUBROUTINE condno



SUBROUTINE delcau (cauchy, overch, overfl, itnum, maxexp, n, nunit, output,  &
                   beta, caulen, delta, epsmch, maxstp, newlen, sqrtz,  &
                   a, delf, scalex)

! N.B. Argument WV1 has been removed.

!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE ESTABLISHES AN INITIAL TRUST REGION, DELTA,
!    IF ONE IS NOT SPECIFIED BY THE USER AND FINDS THE LENGTH OF
!    THE SCALED CAUCHY STEP, CAULEN, AT EACH STEP IF THE DOUBLE
!    DOGLEG OPTION IS BEING USED.
!
!    THE USER HAS TWO CHOICES FOR THE INITIAL TRUST REGION:
!       1)  BASED ON THE SCALED CAUCHY STEP
!       2)  BASED ON THE SCALED NEWTON STEP
!

LOGICAL, INTENT(IN)        :: cauchy
LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: itnum
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(OUT)     :: beta
REAL (dp), INTENT(OUT)     :: caulen
REAL (dp), INTENT(IN OUT)  :: delta
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN OUT)  :: maxstp
REAL (dp), INTENT(IN OUT)  :: newlen
REAL (dp), INTENT(OUT)     :: sqrtz
REAL (dp), INTENT(IN)      :: a(n,n)
REAL (dp), INTENT(IN)      :: delf(n)
REAL (dp), INTENT(IN)      :: scalex(n)

! Local variables

REAL (dp), PARAMETER  :: three = 3.0_dp
REAL (dp)  :: eps, temp, wv1(n)
INTEGER    :: i, j
!
overfl=.false.
eps=ten**(-maxexp)
!
!       IF DELTA IS NOT GIVEN EVALUATE IT USING EITHER THE CAUCHY
!       STEP OR THE NEWTON STEP AS SPECIFIED BY THE USER.
!
!       THE SCALED CAUCHY LENGTH, CAULEN, IS REQUIRED IN TWO CASES.
!       1) WHEN SELECTED AS THE CRITERION FOR THE INITIAL DELTA
!       2) IN THE DOUBLE DOGLEG STEP REGARDLESS OF (1)
!
IF (output > 3) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,10)
  WRITE (nunit,20)
  20 FORMAT ('   *    DETERMINATION OF SCALED CAUCHY STEP LENGTH, CAULEN',  &
             t74, '*')
END IF
!
!       FIND FACTOR WHICH GIVES CAUCHY POINT WHEN MULTPLYING
!       STEEPEST DESCENT DIRECTION, DELF.
!
!       CAULEN = ZETA**1.5/BETA
!              =  SQRTZ**3/BETA
!
wv1(1:n)=delf(1:n)/scalex(1:n)
CALL twonrm (overfl, maxexp, n, epsmch, sqrtz, wv1)
IF (overfl) THEN
  caulen=ten**maxexp
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,40) sqrtz
    40 FORMAT ('   *       ZETA SET TO ', g11.3, ' TO AVOID OVERFLOW', t74, '*')
    WRITE (nunit,50) caulen
    50 FORMAT ('   *       SCALED CAUCHY LENGTH, CAULEN SET TO',  &
               g9.2, ' TO AVOID OVERFLOW', t74, '*')
    IF (itnum == 1) THEN
      WRITE (nunit,60)
      60 FORMAT ('   *       THE PROBLEM SHOULD BE RESCALED',  &
                 ' OR A NEW STARTING POINT CHOSEN', t74, '*')
      WRITE (nunit,70)
      70 FORMAT ('   *       EXECUTION CONTINUES WITH',  &
                 ' SUBSTITUTIONS AS LISTED ABOVE', t74, '*')
    END IF
  END IF
ELSE
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,80) sqrtz
    80 FORMAT ('   *       SQUARE ROOT OF ZETA, SQRTZ: ', g12.3, t74, ' *')
  END IF
END IF
!
!       NOTE: THE LOWER TRIANGLE OF MATRIX A NOW CONTAINS THE
!             TRANSPOSE OF R WHERE A=QR.
!
beta=zero
DO  i=1,n
  temp=zero
  DO  j=i,n
    IF (overch) THEN
      IF (LOG10(ABS(a(j,i))+eps) + LOG10(ABS(delf(j))+eps)  > maxexp) THEN
        caulen=SQRT(epsmch)
        GO TO 120
      END IF
    END IF
    temp=temp + a(j,i)*delf(j)/(scalex(j)*scalex(j))
  END DO
  beta=beta + temp*temp
END DO
IF (output > 4) THEN
  WRITE (nunit,10)
  WRITE (nunit,110) beta
  110 FORMAT ('   *       BETA: ', g11.3, '      NOTE: ',  &
              'CAULEN = ZETA**1.5/BETA', t74, '*')
  WRITE (nunit,10)
END IF
!
!       AVOID OVERFLOWS IN FINDING CAULEN.
!
IF (three*LOG10(sqrtz+eps) - LOG10(beta+eps) < maxexp   &
       .AND. (.NOT.overfl) .AND. beta /= zero) THEN
!
!          NORMAL DETERMINATION.
!
  caulen=sqrtz*sqrtz*sqrtz/beta
!
!          THIS STEP AVOIDS DIVIDE BY ZERO IN DOGLEG IN THE (RARE) CASE
!          WHERE DELF(I)=0 FOR ALL I BUT THE POINT IS NOT YET A SOLUTION -
!          MOST LIKELY A BAD STARTING ESTIMATE.
!
  caulen=MAX(caulen,ten**(-maxexp))
!
ELSE
!
!          SUBSTITUTION TO AVOID OVERFLOW.
!
  caulen=ten**maxexp
END IF

120 IF (output > 3) THEN
  WRITE (nunit,10)
  WRITE (nunit,130) caulen
  130 FORMAT ('   *       SCALED CAUCHY LENGTH, CAULEN: ', g12.3, t74, ' *')
END IF
!
!       ESTABLISH INITIAL TRUST REGION IF NEEDED.
!
IF (delta < zero) THEN
!
!          USE DISTANCE TO CAUCHY POINT OR LENGTH OF NEWTON STEP.
!
  IF (cauchy) THEN
    delta=MIN(caulen,maxstp)
  ELSE
    delta=MIN(newlen,maxstp)
  END IF
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,140) delta
    140 FORMAT ('   *       INITIAL TRUST REGION SIZE, DELTA: ',  &
                g12.3, t74, '*')
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE DELCAU.
!
END SUBROUTINE delcau



SUBROUTINE deufls (abort, deuflh, geoms, overch, overfl, qnfail, qrsing,  &
                   restrt, sclfch, sclxch, acpcod, acptcr, contyp, itnum,  &
                   jupdm, maxexp, maxlin, n, nfunc, nunit, output, qnupdm, &
                   stopcr, alpha, confac, delfts, epsmch, fcnmax, fcnnew,  &
                   fcnold, lambda, newmax, sbrnrm, sigma, a, astore, boundl, &
                   boundu, delf, fvec, hhpi, jac, rdiag, rhs, s, sbar,   &
                   scalef, scalex, sn, xc, xplus, fvecev)

! N.B. Argument WV2 has been removed.

!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE CONDUCTS A LINE SEARCH IN THE NEWTON DIRECTION
!    IF NO CONSTRAINTS ARE VIOLATED.  IF THE FIRST TRIAL IS A FAILURE
!    THERE ARE TWO TYPES OF LINE SEARCH AVAILABLE.
!      1)  REDUCE THE RELAXATION FACTOR, LAMBDA, TO SIGMA*LAMBDA
!          WHERE SIGMA IS USER-SPECIFIED (GEOMETRIC LINE SEARCH)
!      2)  AT THE FIRST STEP MINIMIZE A QUADRATIC THROUGH THE OBJECTIVE
!          FUNCTION AT THE CURRENT POINT AND TRIAL ESTIMATE (WHICH MUST BE A
!          FAILURE) WITH INITIAL SLOPE DELFTS.  AT SUBSEQUENT STEPS MINIMIZE
!          A CUBIC THROUGH THE OBJECTIVE FUNCTION AT THE TWO MOST RECENT
!          FAILURES AND THE CURRENT POINT, AGAIN USING THE INITIAL SLOPE,
!          DELFTS.
!
!    CONVIO  INDICATES A CONSTRAINT VIOLATION BY ONE OR MORE COMPONENTS
!    FRSTST  INDICATES THE FIRST STEP IN THE LINE SEARCH.
!
!    RATIO   RATIO OF PROPOSED STEP LENGTH IN (I)TH DIRECTION
!            TO DISTANCE FROM (I)TH COMPONENT TO BOUNDARY VIOLATED
!    RATIOM  MINIMUM OF RATIOS FOR ALL CONSTAINTS VIOLATED
!

LOGICAL, INTENT(OUT)       :: abort
LOGICAL, INTENT(IN)        :: deuflh
LOGICAL, INTENT(IN)        :: geoms
LOGICAL, INTENT(IN OUT)    :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(IN OUT)    :: qnfail
LOGICAL, INTENT(IN)        :: qrsing
LOGICAL, INTENT(IN OUT)    :: restrt
LOGICAL, INTENT(IN)        :: sclfch
LOGICAL, INTENT(IN OUT)    :: sclxch
INTEGER, INTENT(OUT)       :: acpcod
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(IN)        :: contyp
INTEGER, INTENT(IN)        :: itnum
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: maxlin
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN OUT)    :: nfunc
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: qnupdm
INTEGER, INTENT(IN OUT)    :: stopcr
REAL (dp), INTENT(IN)      :: alpha
REAL (dp), INTENT(IN)      :: confac
REAL (dp), INTENT(IN)      :: delfts
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN)      :: fcnmax
REAL (dp), INTENT(OUT)     :: fcnnew
REAL (dp), INTENT(IN)      :: fcnold
REAL (dp), INTENT(IN OUT)  :: lambda
REAL (dp), INTENT(IN)      :: newmax
REAL (dp), INTENT(OUT)     :: sbrnrm
REAL (dp), INTENT(IN)      :: sigma
REAL (dp), INTENT(IN)      :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: astore(n,n)
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN)      :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(OUT)     :: fvec(n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(IN)      :: jac(n,n)
REAL (dp), INTENT(IN)      :: rdiag(n)
REAL (dp), INTENT(OUT)     :: rhs(n)
REAL (dp), INTENT(OUT)     :: s(n)
REAL (dp), INTENT(OUT)     :: sbar(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: sn(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN OUT)  :: xplus(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

REAL (dp)  :: acubic, bcubic, disc, dlftsm, eps, factor, fplpre, lampre,  &
              lamtmp, ratio, ratiom, wv2(n)
INTEGER    :: i, j, k
LOGICAL    :: convio, frstst

REAL (dp), PARAMETER  :: point1 = 0.1_dp, point5 = 0.5_dp, three = 3.0_dp

!
frstst=.true.
overfl=.false.
abort = .FALSE.     ! Added by AJM
eps=ten**(-maxexp)
DO  k=1,maxlin
  ratiom=one
!
  convio=.false.
!
!          FIND TRIAL POINT AND CHECK IF CONSTRAINTS VIOLATED (IF
!          CONTYP IS NOT EQUAL TO ZERO).
!
  DO  i=1,n
!
!          NOTE: WV2 MARKS VIOLATIONS.  WV2(I) CHANGES TO 1 FOR LOWER BOUND
!                VIOLATIONS AND TO 2 FOR UPPER BOUND VIOLATIONS.
!                CONSTRAINT VIOLATIONS CAN ONLY OCCUR AT THE FIRST STEP.
!
    wv2(i)=-one
    xplus(i)=xc(i) + lambda*sn(i)
    IF (contyp > 0 .AND. frstst) THEN
      IF (xplus(i) < boundl(i)) THEN
        convio=.true.
        wv2(i)=one
      ELSE IF (xplus(i) > boundu(i)) THEN
        convio=.true.
        wv2(i)=two
      ELSE
        wv2(i)=-one
      END IF
    END IF
  END DO
!
!       IF CONSTRAINTS ARE VIOLATED FIRST REDUCE THE STEP SIZES FOR THE
!       VIOLATING COMPONENTS TO OBTAIN A FEASIBLE POINT.
!       IF THE DIRECTION TO THIS MODIFIED POINT IS NOT A DESCENT DIRECTION
!       OR IF THE MODIFIED STEP DOES NOT LEAD TO AN ACCEPTABLE POINT THEN
!       RETURN TO THE NEWTON DIRECTION AND START A LINE SEARCH AT A FEASIBLE
!       POINT WHERE THE COMPONENT WHICH HAS THE SMALLEST VALUE OF RATIO
!       (DEFINED BELOW) IS TAKEN TO "CONFAC" OF THE DISTANCE TO THE BOUNDARY.
!       DEFAULT VALUE OF CONFAC IS 0.95.
!
  IF (convio) THEN
    IF (output > 3) THEN
      WRITE (nunit,20)
      20 FORMAT ('   *', t74, '*')
      WRITE (nunit,30) k
      30 FORMAT ('   *          LINE SEARCH STEP:', i3, t74, '*')
      WRITE (nunit,20)
      WRITE (nunit,40) lambda
      40 FORMAT ('   *          LAMBDA FOR ATTEMPTED STEP: ', g12.3, t74, '*')
      WRITE (nunit,50)
      50 FORMAT ('   *          CONSTRAINT VIOLATED', t74, '*')
      WRITE (nunit,20)
      WRITE (nunit,60)
      60 FORMAT ('   *          TRIAL ESTIMATES (VIOLATIONS MARKED)', t74, '*')
      WRITE (nunit,20)
      DO  i=1,n
        IF (wv2(i) > zero) THEN
          WRITE (nunit,70) i,xplus(i)
          70 FORMAT ('   *', t17, 'XPLUS(', i3, ') = ', g12.3, '  *', t74, '*')
        ELSE
          WRITE (nunit,80) i,xplus(i)
          80 FORMAT ('   *', t17, 'XPLUS(', i3, ') = ', g12.3, t74, '*')
        END IF
      END DO
    END IF
    DO  i=1,n
      IF (wv2(i) > zero) THEN
!
!                   FIND RATIO FOR THIS VIOLATING COMPONENT.
!
        IF (wv2(i) == one) THEN
          ratio=-(xc(i)-boundl(i))/(xplus(i)-xc(i))
        ELSE IF (wv2(i) == two) THEN
          ratio=(boundu(i)-xc(i))/(xplus(i)-xc(i))
        END IF
!
!                   NOTE: THIS LINE IS FOR OUTPUT PURPOSES ONLY.
!
        wv2(i)=ratio
!
        ratiom=MIN(ratiom,ratio)
        IF (ratio > point5) THEN
          s(i)=confac*ratio*lambda*sn(i)
        ELSE
!
!                   WITHIN BUFFER ZONE - ONLY TAKE 1/2
!                   OF THE STEP YOU WOULD TAKE OTHERWISE.
!
          s(i)=point5*confac*ratio*lambda*sn(i)
        END IF
!
!                ESTABLISH MODIFIED TRIAL POINT.
!
        xplus(i)=xc(i) + s(i)
      ELSE
!
!                FOR NONVIOLATORS XPLUS REMAINS UNCHANGED BUT THE COMPONENT
!                OF S IS LOADED TO CHECK THE DIRECTIONAL DERIVATIVE.
!
        s(i)=lambda*sn(i)
      END IF
    END DO
    IF (output > 3) THEN
      WRITE (nunit,20)
      WRITE (nunit,110)
      110 FORMAT ('   *       NEW S AND XPLUS VECTORS',   &
                  ' (WITH RATIOS FOR VIOLATIONS)', t74, '*')
      WRITE (nunit,20)
      WRITE (nunit,120)
      120 FORMAT ('   *       NOTE: RATIOS ARE RATIO OF',   &
                  ' LENGTH TO BOUNDARY FROM CURRENT', t74, '*')
      WRITE (nunit,130)
      130 FORMAT ('   *', t16, 'X VECTOR TO MAGNITUDE OF',   &
                  ' CORRESPONDING PROPOSED STEP', t74, '*')
      WRITE (nunit,20)
      DO  i=1,n
        IF (wv2(i) < zero) THEN
          WRITE (nunit,140) i,s(i),i,xplus(i)
          140 FORMAT ('   *       S(', i3, ') = ', g12.3, '    XPLUS(', i3, &
                      ') = ', g12.3, t74, '*')
        ELSE
          WRITE (nunit,150) i,s(i),i,xplus(i),wv2(i)
          150 FORMAT ('   *       S(', i3, ') = ', g12.3, '    XPLUS(', i3,  &
                      ') = ', g12.3, ' ', g11.3, t74, '*')
        END IF
      END DO
      WRITE (nunit,20)
      WRITE (nunit,170) ratiom
      170 FORMAT ('   *       MINIMUM OF RATIOS, RATIOM: ', g12.3, t74, '*')
    END IF
!
!             CHECK DIRECTIONAL DERIVATIVE FOR MODIFIED POINT, DLFTSM.
!
    CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
                 dlftsm, delf, s)
    overfl=.false.
    IF (output > 3) THEN
      WRITE (nunit,20)
      WRITE (nunit,180) dlftsm
      180 FORMAT ('   *       INNER PRODUCT OF DELF AND S FOR',   &
                  ' MODIFIED S: ', g12.3, t74, '*')
    END IF
!
!             IF INNER PRODUCT IS POSITIVE RETURN TO NEWTON DIRECTION
!             AND CONDUCT A LINE SEARCH WITHIN THE FEASIBLE REGION.
!
    IF (dlftsm > zero) THEN
      IF (output > 3) THEN
        WRITE (nunit,190)
        190 FORMAT ('   *       DELFTS > 0  START LINE',   &
                    ' SEARCH AT LAMBDA = CONFAC*LAMBDA*RATIOM', t74, '*')
        WRITE (nunit,200)
        200 FORMAT ('   *       NOTE: NO TRIAL POINT WAS',   &
                    ' EVALUATED AT THIS STEP OF LINE SEARCH', t74, '*')
      END IF
!
!             THE STARTING POINT IS SET AT JUST INSIDE THE MOST
!             VIOLATED BOUNDARY.
!
      lambda=confac*ratiom*lambda
!
!                LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE.
!
      CYCLE
    END IF
!
  END IF
!
!       NO CONSTRAINTS VIOLATED - EVALUATE RESIDUAL VECTOR AT NEW POINT.
!
  CALL fvecev (overfl, n, fvec, xplus)
  nfunc=nfunc + 1
!
!          CHECK FOR OVERFLOW IN FUNCTION VECTOR EVALUATION.
!          IF SO, REDUCE STEP LENGTH AND CONTINUE LINE SEARCH.
!
  IF (overfl) THEN
!
    IF (output > 3) THEN
      WRITE (nunit,210)
      210 FORMAT ('   *       OVERFLOW IN FUNCTION VECTOR',   &
                  ' - STEP LENGTH REDUCED', t74, '*')
!
!             FORCE STEP TO BE WITHIN CONSTRAINTS - DON'T CALL
!             THIS THE FIRST STEP, I.E. FRSTST STAYS AT TRUE.
!
    END IF
    IF (convio) THEN
      lambda=ratiom*confac*lambda
    ELSE
      lambda=sigma*lambda
    END IF
!
!             LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE.
!
    CYCLE
  END IF
!
!          EVALUATE (POSSIBLY SCALED) OBJECTIVE FUNCTION AT NEW POINT.
!
  CALL fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef)
  IF (overfl) THEN
!
    IF (output > 3) THEN
      WRITE (nunit,220)
      220 FORMAT ('   *       OVERFLOW IN OBJECTIVE FUNCTION',   &
                  ' - STEP LENGTH REDUCED', t74, '*')
!
!                FORCE STEP TO BE WITHIN CONSTRAINTS - DON'T CALL
!                THIS THE FIRST STEP, I.E. FRSTST STAYS AT TRUE.
!
    END IF
    IF (convio) THEN
      lambda=ratiom*confac*lambda
    ELSE
      lambda=sigma*lambda
    END IF
    CYCLE
  END IF
!
!       IF DEUFLHARD'S METHOD IS BEING USED FOR EITHER RELAXATION FACTOR
!       INITIALIZATION OR THE SECOND ACCEPTANCE CRITERION THEN EVALUATE SBAR.
!       EVALUATION METHOD DEPENDS UPON WHETHER THE JACOBIAN WAS PERTURBED
!       IN THE SOLUTION OF THE LINEAR SYSTEM.
!       LOGICAL VARIABLE QRSING IS TRUE IF PERTURBATION TOOK PLACE.
!
  IF (deuflh .OR. acptcr == 12) THEN
    IF (qrsing) THEN
!
!             FORM -J^F AS RIGHT HAND SIDE - METHOD DEPENDS ON WHETHER
!             QNUPDM EQUALS 0 OR 1 IF A QUASI-NEWTON UPDATE IS BEING USED.
!             IF JUPDM IS 0 THEN THE NEWTON STEP HAS BEEN FOUND IN SUBROUTINE
!             NSTPUN.
!
      IF (jupdm == 0 .OR. qnupdm == 0) THEN
!
!                   UNSCALED JACOBIAN IN MATRIX JAC.
!
        DO  i=1,n
          IF (sclfch) THEN
            wv2(i)=-fvec(i)*scalef(i)*scalef(i)
          ELSE
            wv2(i)=-fvec(i)
          END IF
        END DO
        CALL atbmul (n, n, 1, 1, n, n, jac, wv2, rhs)
      ELSE
!
!                   R IN UPPER TRIANGLE OF A PLUS RDIAG AND Q^ IN JAC
!                   - FROM QR DECOMPOSITION OF SCALED JACOBIAN.
!
        DO  i=1,n
          wv2(i)=zero
          DO  j=1,n
            wv2(i)=wv2(i) - jac(i,j)*fvec(j)*scalef(j)
          END DO
        END DO
        rhs(1)=rdiag(1)*wv2(1)
        DO  j=2,n
          rhs(j)=zero
          DO  i=1,j-1
            rhs(j)=rhs(j) + a(i,j)*wv2(i)
          END DO
          rhs(j)=rhs(j) + rdiag(j)*wv2(j)
        END DO
      END IF
      CALL chsolv (overch, overfl, maxexp, n, nunit, output, a, rhs, sbar)
    ELSE
!
!             RIGHT HAND SIDE IS -FVEC.
!
      IF (qnupdm == 0 .OR. jupdm == 0) THEN
!
!                QR DECOMPOSITION OF SCALED JACOBIAN STORED IN ASTORE.
!
        sbar(1:n)=-fvec(1:n)*scalef(1:n)
        CALL qrsolv (overch, overfl, maxexp, n, nunit, output,  &
                     astore, hhpi, rdiag, sbar)
      ELSE
!
!                   SET UP RIGHT HAND SIDE - MULTIPLY -FVEC BY Q^
!                   (STORED IN JAC).
!
        wv2(1:n)=-fvec(1:n)*scalef(1:n)
        CALL avmul (n, n, n, n, jac, wv2, sbar)
        CALL rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sbar)
      END IF
    END IF
!
!             NORM OF SCALED SBAR IS NEEDED FOR SECOND ACCEPTANCE TEST.
!
    IF (acptcr == 12) THEN
      wv2(1:n)=scalex(1:n)*sbar(1:n)
      CALL twonrm (overfl, maxexp, n, epsmch, sbrnrm, wv2)
    END IF
  END IF
!
  IF (output > 3) THEN
    WRITE (nunit,20)
    WRITE (nunit,30) k
    WRITE (nunit,20)
    WRITE (nunit,40) lambda
    WRITE (nunit,20)
    IF (.NOT.convio) THEN
      WRITE (nunit,310)
      310 FORMAT ('   *          NEW COMPONENT/FCN VECTORS',   &
                  '  (XPLUS(I) = XC(I) + LAMBDA*SN(I))', t74, '*')
    ELSE
      WRITE (nunit,320)
      320 FORMAT ('   *          NEW FUNCTION VECTORS',   &
                  ' AT MODIFIED POINT', t74, '*')
    END IF
    WRITE (nunit,20)
    DO  i=1,n
      WRITE (nunit,330) i,xplus(i),i,fvec(i)
      330 FORMAT ('   *          XPLUS(', i3, ') = ', g12.3, '     FVEC(',  &
                  i3, ') = ', g12.3, t74, '*')
    END DO
    WRITE (nunit,20)
    IF (.NOT.sclfch) THEN
      WRITE (nunit,350) fcnnew
      350 FORMAT ('   *          OBJECTIVE FUNCTION VALUE',   &
                  ' AT XPLUS: .........', g12.3, t74, '*')
    ELSE
      WRITE (10,360) fcnnew
      360 FORMAT ('   *          SCALED OBJECTIVE FUNCTION VALUE',   &
                  ' AT XPLUS: ..', g12.3, t74, '*')
    END IF
    WRITE (nunit,370) fcnmax+alpha*lambda*delfts
    370 FORMAT ('   *          FCNMAX + ALPHA*LAMBDA*DELFTS:',   &
                ' ................', g12.3, t74, '*')
    IF (deuflh .OR. acptcr == 12) THEN
      IF (itnum > 0) THEN
        IF (.NOT.sclxch) THEN
          WRITE (nunit,20)
          WRITE (nunit,380)
          380 FORMAT ('   *          DEUFLHARD SBAR VECTOR', t74, '*')
          WRITE (nunit,20)
          DO  i=1,n
            WRITE (nunit,390) i,sbar(i)
            390 FORMAT ('   *          SBAR(', I3, ') = ', g12.3, t74, '*')
          END DO
        ELSE
          WRITE (nunit,20)
          WRITE (nunit,410)
          410 FORMAT ('   *          DEUFLHARD SBAR VECTOR', t49,  &
                      'IN SCALE DX UNITS', t74, '*')
          WRITE (nunit,20)
          DO  i=1,n
            WRITE (nunit,420) i,sbar(i),i,scalex(i)*sbar(i)
            420 FORMAT ('   *          SBAR(', I3, ') = ', g12.3,  &
                        '        SBAR(', i3, ') = ', g12.3, t74, '*')
          END DO
        END IF
      END IF
    END IF
    IF (acptcr == 12) THEN
      WRITE (nunit,20)
      IF (.NOT.sclxch) THEN
        WRITE (nunit,440) sbrnrm
        440 FORMAT ('   *          VALUE OF SBRNRM',   &
                    ' AT XPLUS: ..................', g12.3, t74, '*')
      ELSE
        WRITE (nunit,450) sbrnrm
        450 FORMAT ('   *          VALUE OF SCALED SBRNRM',   &
                    ' AT XPLUS: ...........', g12.3, t74, '*')
      END IF
      WRITE (nunit,460) newmax
      460 FORMAT ('   *          NEWMAX: ..............',   &
                  '......................', g12.3, t74, '*')
    END IF
  END IF
!
!          CHECK FOR ACCEPTABLE STEP.
!
  IF (fcnnew < fcnmax + alpha*lambda*delfts) THEN
    acpcod=1
!
!             NOTE: STEP WILL BE ACCEPTED REGARDLESS OF NEXT TEST.
!                   THIS SECTION IS FOR BOOKKEEPING ONLY.
!
    IF (acptcr == 12) THEN
      IF (sbrnrm < newmax) THEN
        acpcod=12
      END IF
    END IF
!
    RETURN
  END IF
  IF (acptcr == 12 .AND. sbrnrm < newmax) THEN
    acpcod=2
    RETURN
  END IF
!
!          FAILURE OF STEP ACCEPTANCE TEST.
!
  IF (convio) THEN
    lambda=confac*ratiom*lambda
!
!             LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE.
!
    CYCLE
  END IF
  IF (lambda == zero) THEN
    IF (output > 0) THEN
      WRITE (nunit,20)
      WRITE (nunit,470)
      470 FORMAT ('   *       LAMBDA IS 0.0:  NO PROGRESS',   &
                  ' POSSIBLE - CHECK BOUNDS OR START', t74, '*')
    END IF
    abort=.true.
    RETURN
  END IF
  IF (geoms) THEN
!
!             GEOMETRIC LINE SEARCH
!
    lambda=sigma*lambda
!
  ELSE
!
    IF (frstst) THEN
      frstst=.false.
!
!                FIND MINIMUM OF QUADRATIC AT FIRST STEP.
!
      lamtmp=-(lambda*lambda)*delfts/(two*(fcnnew-fcnold-lambda* delfts))
      IF (output > 4) THEN
        WRITE (nunit,20)
        WRITE (nunit,480) lamtmp
        480 FORMAT ('   *', t16, 'TEMPORARY LAMBDA FROM',   &
                    ' QUADRATIC MODEL:  ', g11.3, t74, '*')
      END IF
!
    ELSE
!
!                FIND MINIMUM OF CUBIC AT SUBSEQUENT STEPS.
!
      factor=one/(lambda-lampre)
      IF (lambda*lambda == zero) THEN
        lambda=sigma*lambda
!
!                   NOTE: IF THIS LAMBDA**2 WAS ZERO ANY SUBSEQUENT
!                         LAMBDA**2 WILL ALSO BE ZERO.
!
        CYCLE
      END IF
      acubic=factor*((one/lambda*lambda)*(fcnnew-fcnold-lambda*  &
          delfts)-((one/lampre*lambda)*(fplpre-fcnold-lampre*delfts)) )
      bcubic=factor*((-lampre/lambda*lambda)*(fcnnew-fcnold-  &
          lambda*delfts)+((lambda/lampre*lambda)*(fplpre-fcnold- lampre*delfts)))
      IF (two*LOG10(ABS(bcubic)+eps) > DBLE(maxexp)) THEN
        lamtmp=sigma*lambda
        IF (output > 4) THEN
          WRITE (nunit,20)
          WRITE (nunit,490)
          490 FORMAT ('   *', t16, 'POTENTIAL OVERFLOW IN CALCULATING',  &
                      ' TRIAL LAMBDA FROM', t74, '*')
          WRITE (nunit,500)
          500 FORMAT ('   *', t16, 'CUBIC MODEL - LAMBDA',   &
                      ' SET TO SIGMA*LAMBDA', t74, '*')
        END IF
      ELSE
        disc=bcubic*bcubic - three*acubic*delfts
        IF (ABS(acubic) <= epsmch) THEN
          lamtmp=-delfts/(two*bcubic)
        ELSE
          IF (disc < zero) THEN
            lamtmp=sigma*lambda
          ELSE
            lamtmp=(-bcubic + SQRT(disc))/(three*acubic)
          END IF
        END IF
        IF (output > 4) THEN
          WRITE (nunit,20)
          WRITE (nunit,510) lamtmp
          510 FORMAT ('   *', t16, 'TEMPORARY LAMBDA FROM',   &
                      ' CUBIC MODEL : .....', g11.3, t74, '*')
        END IF
      END IF
      IF (lamtmp > sigma*lambda) THEN
        lamtmp=sigma*lambda
        IF (output > 4) THEN
          WRITE (nunit,20)
          WRITE (nunit,520)
          520 FORMAT ('   *', t16, 'LAMTMP TOO LARGE - REDUCED',   &
                      ' TO SIGMA*LAMBDA', t74, '*')
        END IF
      END IF
    END IF
    lampre=lambda
    fplpre=fcnnew
    IF (lamtmp < point1*lambda) THEN
      IF (output > 4) THEN
        WRITE (nunit,20)
        WRITE (nunit,530)
        530 FORMAT ('   *', t16, 'LAMTMP TOO SMALL - INCREASED',   &
                    ' TO 0.1*LAMBDA', t74, '*')
      END IF
      lambda=point1*lambda
    ELSE
      IF (output > 4 .AND. lamtmp /= sigma*lambda) THEN
        WRITE (nunit,20)
        WRITE (nunit,540)
        540 FORMAT ('   *', t16, 'LAMTMP WITHIN LIMITS - ',   &
                    'LAMBDA SET TO LAMTMP', t74, '*')
      END IF
      lambda=lamtmp
    END IF
  END IF
END DO
!
!       FAILURE IN LINE SEARCH
!
acpcod=0
!
!       IF A QUASI-NEWTON STEP HAS FAILED IN THE LINE SEARCH THEN SET QNFAIL
!       TO TRUE ANS RETURN TO SUBROUTINE NNES.  THIS WILL CAUSE THE JACOBIAN
!       TO BE RE-EVALUATED EXPLICITLY AND A LINE SEARCH IN THE NEW DIRECTION
!       CONDUCTED.
!
IF (.NOT.restrt) THEN
  qnfail=.true.
  RETURN
END IF
!
!       FALL THROUGH MAIN LOOP - WARNING GIVEN.
!
IF (output > 2 .AND. (.NOT.wrnsup)) THEN
  WRITE (nunit,20)
  WRITE (nunit,570) maxlin
  570 FORMAT ('   *       WARNING: ', I3, ' CYCLES COMPLETED',   &
              ' IN LINE SEARCH WITHOUT SUCCESS', t74, '*')
END IF
IF (stopcr == 2 .OR. stopcr == 3) THEN
  stopcr=12
  WRITE (nunit,20)
  WRITE (nunit,580)
  580 FORMAT ('   *       STOPPING CRITERION RESET FROM ',   &
              '2 TO 12 TO AVOID HANGING', t74, '*')
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE DEUFLS.
!
END SUBROUTINE deufls



SUBROUTINE dogleg (frstdg, newtkn, overch, overfl, maxexp, n, notrst, nunit, &
                   output, beta, caulen, delta, etafac, newlen, sqrtz,  &
                   delf, s, scalex, sn, ssdhat, vhat)
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS A TRUST REGION STEP USING THE
!    (DOUBLE) DOGLEG METHOD.
!

LOGICAL, INTENT(IN OUT)    :: frstdg
LOGICAL, INTENT(OUT)       :: newtkn
LOGICAL, INTENT(IN OUT)    :: overch
LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN OUT)    :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN OUT)    :: notrst
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN OUT)  :: beta
REAL (dp), INTENT(IN)      :: caulen
REAL (dp), INTENT(IN OUT)  :: delta
REAL (dp), INTENT(IN)      :: etafac
REAL (dp), INTENT(IN)      :: newlen
REAL (dp), INTENT(IN)      :: sqrtz
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(OUT)     :: s(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: sn(n)
REAL (dp), INTENT(IN OUT)  :: ssdhat(n)
REAL (dp), INTENT(IN OUT)  :: vhat(n)

REAL (dp)  :: delfts, eta, factor, gamma, lambda, temp, tempv, zeta
INTEGER    :: i

!
overfl=.false.
eta = 1.0_dp    ! Added by AJM
!
IF (output > 3) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,20) notrst,delta
  20 FORMAT ('   *    TRUST REGION STEP:', i6,   &
             '  TRUST REGION LENGTH, DELTA:', g11.3, t74, '*')
  WRITE (nunit,10)
  WRITE (nunit,30) newlen
  30 FORMAT ('   *       LENGTH OF NEWTON STEP, NEWLEN: ', g11.3, t74, '*')
END IF
!
!       CHECK FOR NEWTON STEP WITHIN TRUST REGION - IF SO USE NEWTON STEP.
!
IF (newlen <= delta) THEN
  s(1:n)=sn(1:n)
  newtkn=.true.
  temp=delta
  delta=newlen
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,50)
    50 FORMAT ('   *       NEWTON STEP WITHIN ACCEPTABLE RANGE',   &
               ' ( <= THAN DELTA)', t74, '*')
    IF (temp == delta) THEN
      WRITE (nunit,60) delta
      60 FORMAT ('   *       DELTA STAYS AT LENGTH OF NEWTON',   &
                 ' STEP: ', g11.3, t74, '*')
    ELSE
      WRITE (nunit,70)
      70 FORMAT ('   *       DELTA SET TO LENGTH OF NEWTON STEP', t74, '*')
    END IF
    WRITE (nunit,10)
    WRITE (nunit,80)
    80 FORMAT ('   *       FULL NEWTON STEP ATTEMPTED', t74, '*')
  END IF
  RETURN
ELSE
!
!          NEWTON STEP NOT WITHIN TRUST REGION - APPLY (DOUBLE)
!          DOGLEG PROCEDURE.  (IF ETAFAC EQUALS 1.0 THEN THE SINGLE
!          DOGLEG PROCEDURE IS BEING APPLIED).
!
  newtkn=.false.
  IF (frstdg) THEN
!
!             SPECIAL SECTION FOR FIRST DOGLEG STEP - CALCULATES
!             CAUCHY POINT (MINIMIZER OF MODEL FUNCTION IN STEEPEST
!             DESCENT DIRECTION OF THE OBJECTIVE FUNCTION).
!
    frstdg=.false.
    IF (output > 4) THEN
      WRITE (nunit,10)
      IF (etafac == one) THEN
        WRITE (nunit,90)
        90 FORMAT ('   *       FIRST SINGLE DOGLEG STEP', t74, '*')
      ELSE
        WRITE (nunit,100)
        100 FORMAT ('   *       FIRST DOUBLE DOGLEG STEP', t74, '*')
      END IF
      WRITE (nunit,10)
      WRITE (nunit,110)
      110 FORMAT ('   *          SCALED CAUCHY STEP', t74, '*')
      WRITE (nunit,10)
    END IF
!
!             NOTE: BETA AND SQRTZ WERE CALCULATED IN SUBROUTINE DELCAU.
!
    zeta=sqrtz*sqrtz
!
!             FIND STEP TO CAUCHY POINT.
!
    factor=-(zeta/beta)
    DO  i=1,n
      ssdhat(i)=factor*(delf(i)/scalex(i))
      IF (output > 4) THEN
        WRITE (nunit,120) i,ssdhat(i)
        120 FORMAT ('   *', t16, 'SSDHAT(', i3, ') = ', g12.3, t74, '*')
      END IF
    END DO
    CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
                 delfts, delf, sn)
    overfl=.false.
!
!             PROTECT AGAINST (RARE) CASE WHEN CALCULATED DIRECTIONAL
!             DERIVATIVE EQUALS ZERO.
!
    IF (delfts /= zero) THEN
!
!                STANDARD EXECUTION.
!
      gamma=(zeta/ABS(delfts))*(zeta/beta)
      eta=etafac + (one-etafac)*gamma
    ELSE
      IF (output > 1 .AND. (.NOT.wrnsup)) THEN
        WRITE (nunit,10)
        WRITE (nunit,140)
        140 FORMAT ('   *    WARNING: DELFTS=0; ETA SET',   &
                    ' TO 1.0 TO AVOID DIVISION BY ZERO', t74, '*')
      END IF
      eta=one
    END IF
    IF (output > 4) THEN
      WRITE (nunit,10)
      WRITE (nunit,150) eta
      150 FORMAT ('   *          ETA = ', g11.3, t74, '*')
      WRITE (nunit,10)
      WRITE (nunit,160)
      160 FORMAT ('   *        VHAT VECTOR     VHAT(I) = ETA*',   &
                  'SN(I)*SCALEX(I) - SSDHAT(I)', t74, '*')
      WRITE (nunit,10)
    END IF
    DO  i=1,n
      vhat(i)=eta*scalex(i)*sn(i) - ssdhat(i)
      IF (output > 4) THEN
        WRITE (nunit,170) i,vhat(i)
        170 FORMAT ('   *', t16, 'VHAT(', i3, ') = ', g12.3, t74, '*')
      END IF
    END DO
  END IF
!
!          ETA*NEWLEN <= DELTA MEANS TAKE STEP IN NEWTON DIRECTION
!          TO TRUST REGION BOUNDARY.
!
  IF (eta*newlen <= delta) THEN
    IF (output > 4) THEN
      WRITE (nunit,10)
      WRITE (nunit,190)
      190 FORMAT ('   *          ETA*NEWLEN <= DELTA     S(I)',   &
                  ' = (DELTA/NEWLEN)*SN(I)', t74, '*')
    END IF
    s(1:n)=(delta/newlen)*sn(1:n)
  ELSE
!
!             DISTANCE TO CAUCHY POINT >= DELTA MEANS TAKE STEP IN
!             STEEPEST DESCENT DIRECTION TO TRUST REGION BOUNDARY.
!
    IF (caulen >= delta) THEN
      IF (output > 4) THEN
        WRITE (nunit,10)
        WRITE (nunit,210)
        210 FORMAT ('   *          CAULEN >= DELTA   S(I)',   &
                    '= (DELTA/CAULEN)*(SSDHAT(I)/SCALEX(I))', t74, '*')
      END IF
      s(1:n) = (delta/caulen)*(ssdhat(1:n)/scalex(1:n))
    ELSE
!
!                TAKE (DOUBLE) DOGLEG STEP.
!
      CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
                   temp, ssdhat, vhat)
      CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
                   tempv, vhat, vhat)
      overfl=.false.
      lambda=(-temp + SQRT(temp*temp - tempv*(caulen*caulen - delta*delta)))/tempv
      IF (output > 4) THEN
        WRITE (nunit,10)
        WRITE (nunit,230)
        230 FORMAT ('   *          S(I) = (SSDHAT(I)+LAMBDA*VHAT(I))',   &
                    '/SCALEX(I)', t74, '*')
        WRITE (nunit,10)
        WRITE (nunit,240) lambda
        240 FORMAT ('   *          WHERE LAMBDA = ', g12.3, t74, '*')
      END IF
      s(1:n)=(ssdhat(1:n) + lambda*vhat(1:n))/scalex(1:n)
    END IF
  END IF
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,260)
    260 FORMAT ('   *          REVISED STEP FROM SUBROUTINE DOGLEG', t74, '*')
    WRITE (nunit,10)
    DO  i=1,n
      WRITE (nunit,270) i,s(i)
      270 FORMAT ('   *', t16, 'S(', I3, ') = ', g12.3, t74, '*')
    END DO
  END IF
END IF
RETURN
END SUBROUTINE dogleg



SUBROUTINE fcnevl (overfl, maxexp, n, nunit, output, epsmch,  &
                   fcnnew, fvec, scalef)

! N.B. Argument WV1 has been removed.

!
!    FEB. 23, 1992
!
!    THE OBJECTIVE FUNCTION, FCNNEW, DEFINED BY:
!
!       FCNNEW: = 1/2(SCALEF*FVEC^SCALEF*FVEC)
!
!    IS EVALUATED BY THIS SUBROUTINE.
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(IN)   :: epsmch
REAL (dp), INTENT(OUT)  :: fcnnew
REAL (dp), INTENT(IN)   :: fvec(n)
REAL (dp), INTENT(IN)   :: scalef(n)

! Local variables

REAL (dp)  :: eps, wv1(n)
!
overfl=.false.
eps=ten**(-maxexp)
!
wv1(1:n)=fvec(1:n)*scalef(1:n)
CALL twonrm (overfl, maxexp, n, epsmch, fcnnew, wv1)
!
!       IF AN OVERFLOW WOULD OCCUR SUBSTITUTE A LARGE VALUE FOR FCNNEW.
!
IF (overfl .OR. two*LOG10(fcnnew+eps) > maxexp) THEN
  overfl=.true.
  fcnnew=ten**maxexp
  IF (output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,20)
    20 FORMAT ('   *', t74, '*')
    WRITE (nunit,30) fcnnew
    30 FORMAT ('   *    WARNING: TO AVOID OVERFLOW',   &
               ' OBJECTIVE FUNCTION SET TO: ', g11.3, t74, '*')
  END IF
  RETURN
END IF
fcnnew=fcnnew*fcnnew/two
RETURN
END SUBROUTINE fcnevl



SUBROUTINE fordif (overfl, j, n, deltaj, fvec, fvecj1, jacfdm, xc, fvecev)
!
!       FEB. 6, 1991
!

LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: j
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: deltaj
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(OUT)     :: fvecj1(n)
REAL (dp), INTENT(IN OUT)  :: jacfdm(n,n)
REAL (dp), INTENT(IN)      :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

CALL fvecev (overfl, n, fvecj1, xc)
IF (.NOT.overfl) THEN
  jacfdm(1:n,j)=(fvecj1(1:n) - fvec(1:n))/deltaj
END IF

RETURN
!
!       LAST CARD OF SUBROUTINE FORDIF.
!
END SUBROUTINE fordif



SUBROUTINE gradf (overch, overfl, restrt, sclfch, sclxch, jupdm, maxexp,  &
                  n, nunit, output, qnupdm, delf, fvecc, jac, scalef, scalex)

! N.B. Argument WV1 has been removed.

!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE COMPUTES THE GRADIENT OF THE FUNCTION
!
!           F=1/2{SCALEF*FVECC)^(SCALEF*FVECC}
!
!    WHICH IS USED AS THE OBJECTIVE FUNCTION FOR MINIMIZATION.
!
!    NOTE: WHEN THE FACTORED FORM OF THE JACOBIAN IS UPDATED IN QUASI-NEWTON
!          METHODS THE GRADIENT IS UPDATED AS WELL IN THE SAME SUBROUTINE -
!          IT IS PRINTED HERE THOUGH.  IN THESE CASES QNUPDM > 0.
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(IN OUT)    :: overfl
LOGICAL, INTENT(IN OUT)    :: restrt
LOGICAL, INTENT(IN)        :: sclfch
LOGICAL, INTENT(IN OUT)    :: sclxch
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: qnupdm
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN OUT)  :: scalex(n)

! Local variables

INTEGER    :: i
REAL (dp)  :: wv1(n)
!
IF (restrt .OR. jupdm == 0 .OR. qnupdm == 0) THEN
!
!          GRADIENT NOT ALREADY UPDATED:    FIND DELF = J^F.
!
  DO  i=1,n
    IF (sclfch) THEN
      wv1(i)=fvecc(i)*scalef(i)*scalef(i)
    ELSE
      wv1(i)=fvecc(i)
    END IF
  END DO
  IF (overch) THEN
!   check each entry individually
    CALL atvov (overfl, maxexp, n, nunit, output, jac, wv1, delf)
  ELSE
    CALL atbmul (n, n, 1, 1, n, n, jac, wv1, delf)
  END IF
!     END IF
!
!       PRINT GRADIENT VECTOR, DELF.
!
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      WRITE (nunit,20)
      20 FORMAT ('   *', t74, '*')
      IF (.NOT.sclfch) THEN
        WRITE (nunit,30)
        30 FORMAT ('   *    GRADIENT OF OBJECTIVE FUNCTION', t74, '*')
      ELSE
        WRITE (nunit,40)
        40 FORMAT ('   *    GRADIENT OF SCALED OBJECTIVE FUNCTION', t74, '*')
      END IF
      WRITE (nunit,20)
      DO  i=1,n
        WRITE (nunit,50) i,delf(i)
        50 FORMAT ('   *', t9, 'DELF(', i3, ') = ', g12.3, t74, '*')
      END DO
    ELSE
      WRITE (nunit,20)
      WRITE (nunit,70)
      70 FORMAT ('   *    GRADIENT OF OBJECTIVE FUNCTION',   &
                 '         IN SCALED X UNITS', t74, '*')
      WRITE (nunit,20)
      DO  i=1,n
        WRITE (nunit,80) i,delf(i),i,scalef(i)*scalef(i)*delf(i)/ scalex(i)
        80 FORMAT ('   *', t9, 'DELF(', i3, ') = ', g12.3,   &
                   '         DELF(', i3, ') = ', g12.3, t74, '*')
      END DO
    END IF
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE GRADF.
!
END SUBROUTINE gradf



SUBROUTINE initch (instop, linesr, newton, overfl, sclfch, sclxch, acptcr,  &
                   contyp, jactyp, jupdm, maxexp, n, nunit, output, qnupdm,  &
                   stopcr, trupdm, epsmch, fcnold, ftol, boundl, boundu,   &
                   fvecc, scalef, scalex, xc, fvecev)

! N.B. Argument WV1 has been removed.

!
!    AUG. 27, 1991
!
!    THIS SUBROUTINE FIRST CHECKS TO SEE IF N IS WITHIN THE ACCEPTABLE RANGE.
!
!    THE SECOND CHECK IS TO SEE IF THE INITIAL ESTIMATE IS
!    ALREADY A SOLUTION BY THE FUNCTION VALUE CRITERION, FTOL.
!
!    THE THIRD CHECK IS MADE TO SEE IF THE NEWTON OPTION IS BEING
!    USED WITH THE LINE SEARCH.  IF NOT A WARNING IS GIVEN AND
!    THE LINE SEARCH OPTION IS INVOKED.
!
!    THE FOURTH CHECK IS TO ENSURE APPLICABILITY OF SELECTED
!    VALUES FOR INTEGER CONSTANTS.
!
!    THE FIFTH CHECK IS TO WARN THE USER IF INITIAL ESTIMATES ARE NOT WITHIN
!    THE RANGES SET BY THE BOUNDL AND BOUNDU VECTORS.
!    CONTYP IS CHANGED FROM 0 TO 1 IF ANY BOUND HAS BEEN SET BY THE USER
!
!    THE SIXTH CHECK ENSURES BOUNDL(I) < BOUNDU(I) FOR ALL I.
!

LOGICAL, INTENT(OUT)       :: instop
LOGICAL, INTENT(IN OUT)    :: linesr
LOGICAL, INTENT(IN OUT)    :: newton
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(OUT)       :: sclfch
LOGICAL, INTENT(OUT)       :: sclxch
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(OUT)       :: contyp
INTEGER, INTENT(IN)        :: jactyp
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: qnupdm
INTEGER, INTENT(IN)        :: stopcr
INTEGER, INTENT(IN)        :: trupdm
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN OUT)  :: fcnold
REAL (dp), INTENT(IN)      :: ftol
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN)      :: boundu(n)
REAL (dp), INTENT(OUT)     :: fvecc(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

LOGICAL    :: frster
REAL (dp)  :: temp1, temp2
INTEGER    :: i
!
instop=.false.
temp1=-ten**maxexp
temp2=ten**maxexp
!
!       CHECK FOR N IN RANGE.
!
IF (n <= 0) THEN
  instop=.true.
  WRITE (nunit,10)
  10 FORMAT (t3,72('*'))
  WRITE (nunit,20)
  20 FORMAT ('   *', t74, '*')
  WRITE (nunit,30)
  30 FORMAT ('   *  N IS OUT OF RANGE - RESET TO POSITIVE INTEGER', t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
!
!       CHECK FOR SCALING FACTORS POSITIVE.
!
frster=.true.
sclfch=.false.
sclxch=.false.
DO  i=1,n
  IF (scalef(i) <= zero) THEN
    IF (frster) THEN
      instop=.true.
      frster=.false.
      WRITE (nunit,10)
    END IF
    WRITE (nunit,20)
    WRITE (nunit,40) i,scalef(i)
    40 FORMAT ('   *       SCALEF(', i3, ') = ', g12.3,   &
               '    SHOULD BE POSITIVE', t74, '*')
  END IF
  IF (scalef(i) /= one) sclfch=.true.
  IF (scalex(i) <= zero) THEN
    IF (frster) THEN
      instop=.true.
      frster=.false.
      WRITE (nunit,10)
    END IF
    WRITE (nunit,20)
    WRITE (nunit,50) i,scalex(i)
    50 FORMAT ('   *       SCALEX(', i3, ') = ', g12.3,   &
               '    SHOULD BE POSITIVE', t74, '*')
  END IF
  IF (scalex(i) /= one) sclxch=.true.
END DO
IF (.NOT.frster) THEN
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
!
!       EVALUATE INITIAL RESIDUAL VECTOR AND OBJECTIVE FUNCTION AND
!       CHECK TO SEE IF THE INITIAL GUESS IS ALREADY A SOLUTION.
!
CALL fvecev (overfl, n, fvecc, xc)
!
!       NOTE: NUMBER OF LINE SEARCH FUNCTION EVALUATIONS, NFUNC,
!             INITIALIZED AT 1 WHICH REPRESENTS THIS EVALUATION.
!
IF (overfl) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,70)
  70 FORMAT ('   *       OVERFLOW IN INITIAL FUNCTION VECTOR EVALUATION',  &
             t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,10)
  instop=.true.
  RETURN
END IF
CALL fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnold, fvecc, scalef)
IF (overfl) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,80)
  80 FORMAT ('   *       OVERFLOW IN INITIAL OBJECTIVE FUNCTION EVALUATION', &
             t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,10)
  instop=.true.
  RETURN
END IF
!
!       CHECK FOR SOLUTION USING SECOND STOPPING CRITERION.
!
DO  i=1,n
  IF (ABS(fvecc(i)) > ftol) GO TO 130
END DO
instop=.true.
WRITE (nunit,10)
WRITE (nunit,20)
WRITE (nunit,100)
100 FORMAT ('   *  WARNING: THIS IS ALREADY A SOLUTION',   &
            ' BY THE CRITERIA OF THE SOLVER', t74, '*')
WRITE (nunit,20)
!
!    IF THE PROBLEM IS BADLY SCALED THE OBJECTIVE FUNCTION MAY MEET THE
!    TOLERANCE ALTHOUGH THE INITIAL ESTIMATE IS NOT THE SOLUTION.
!
WRITE (nunit,110)
110 FORMAT ('   *  THIS MAY POSSIBLY BE ALLEVIATED BY ',   &
            'RESCALING THE PROBLEM IF THE', t74, '*')
WRITE (nunit,120)
120 FORMAT ('   *  INITIAL ESTIMATE IS KNOWN NOT TO BE A SOLUTION', t74, '*')
WRITE (nunit,20)
WRITE (nunit,10)
!
!       CHECK FOR NEWTON'S METHOD REQUESTED BUT LINE SEARCH NOT BEING USED.
!
130 IF (newton .AND. (.NOT.linesr)) THEN
  linesr=.true.
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,140)
  140 FORMAT ('   *  WARNING: INCOMPATIBLE OPTIONS',   &
              ': NEWTON=.TRUE. AND LINESR=.FALSE.', t74, '*')
  WRITE (nunit,150)
  150 FORMAT ('   *  LINESR SET TO .TRUE.; EXECUTION OF NEWTON ',   &
              'METHOD CONTINUING', t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
!
!       CHECK INTEGER CONSTANTS.
!
IF (acptcr /= 1 .AND. acptcr /= 12) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,160) acptcr
  160 FORMAT ('   *  ACPTCR NOT AN ACCEPTABLE VALUE: ',i5, t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
IF (jactyp < 0 .OR. jactyp > 3) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,170) jactyp
  170 FORMAT ('   *  JACTYP:', i5, ' - NOT IN PROPER RANGE', t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
IF (stopcr /= 1 .AND. stopcr /= 12 .AND. stopcr /= 2 .AND. stopcr /= 3) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,180) stopcr
  180 FORMAT ('   *  STOPCR NOT AN ACCEPTABLE VALUE: ', i5, t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
IF (qnupdm < 0 .OR. qnupdm > 1) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,190) qnupdm
  190 FORMAT ('   *  QNUPDM:', i5, ' - NOT IN PROPER RANGE', t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
IF (trupdm < 0 .OR. trupdm > 1) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,200) trupdm
  200 FORMAT ('   *  TRUPDM:', i5, ' - NOT IN PROPER RANGE', t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
IF (jupdm < 0 .OR. jupdm > 2) THEN
  WRITE (nunit,10)
  WRITE (nunit,20)
  WRITE (nunit,210) jupdm
  210 FORMAT ('   *  JUPDM:', i5, ' - NOT IN PROPER RANGE', t74, '*')
  instop=.true.
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
!
!       CHECK FOR INITIAL ESTIMATES NOT WITHIN SPECIFIED BOUNDS AND
!       SET CONTYP TO 1 => AT LEAST ONE BOUND IS IN EFFECT.
!
contyp=0
DO  i=1,n
  IF ((boundl(i) /= temp1 .OR. boundu(i) /= temp2)) THEN
    contyp=1
    EXIT
  END IF
END DO

frster=.true.
IF (contyp /= 0) THEN
  DO  i=1,n
!
!             CHECK FOR INITIAL ESTIMATES OUT OF RANGE AND LOWER
!             BOUND GREATER THAN OR EQUAL TO THE UPPER BOUND.
!
    IF (xc(i) < boundl(i) .OR. xc(i) > boundu(i)) THEN
      IF (frster) THEN
        instop=.true.
        frster=.false.
        WRITE (nunit,10)
        WRITE (nunit,20)
        WRITE (nunit,240)
        240 FORMAT ('   *       COMPONENTS MUST BE WITHIN BOUNDS', t74, '*')
        WRITE (nunit,20)
        WRITE (nunit,250)
        250 FORMAT ('   *', t11, 'NO.', t22, 'XC', t39, 'BOUNDL', t55,   &
                    'BOUNDU', t74, '*')
        WRITE (nunit,20)
      END IF
      WRITE (nunit,260) i,xc(i),boundl(i),boundu(i)
      260 FORMAT ('   *', t10, i3, '   ', g12.3, t39, g12.3, '    ', g12.3, &
                  t74, '*')
    END IF
  END DO
  IF (.NOT.frster) THEN
    WRITE (nunit,20)
    WRITE (nunit,10)
  END IF
!
  frster=.true.
  DO  i=1,n
    IF (boundl(i) >= boundu(i)) THEN
      IF (frster) THEN
        frster=.false.
        WRITE (nunit,10)
        WRITE (nunit,20)
        WRITE (nunit,280)
        280 FORMAT ('   *       LOWER BOUND MUST BE LESS THAN',   &
                    ' UPPER BOUND - VIOLATIONS LISTED', t74, '*')
        WRITE (nunit,20)
      END IF
      WRITE (nunit,290) i,boundl(i),i,boundu(i)
      290 FORMAT ('   *       BOUNDL(', i3, ') = ', g12.3, '    BOUNDU(',  &
                  i3, ') = ', g12.3, t74, '*')
    END IF
  END DO
  IF (.NOT.frster) THEN
    WRITE (nunit,20)
    WRITE (nunit,10)
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE INITCH.
!
END SUBROUTINE initch



SUBROUTINE innerp (overch, overfl, maxexp, ldima, ldimb, n, nunit,  &
                   output, dtpro, a, b)
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE FINDS THE INNER PRODUCT OF TWO VECTORS, A AND B.
!    IF OVERCH IS FALSE, UNROLLED LOOPS ARE USED.
!
!    LDIMA IS THE DIMENSION OF A
!    LDIMB IS THE DIMENSION OF B
!    N IS THE DEPTH INTO A AND B THE INNER PRODUCT IS DESIRED.
!    (USUALLY LDIMA=LDIMB=N)
!

LOGICAL, INTENT(IN)     :: overch
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: ldima
INTEGER, INTENT(IN)     :: ldimb
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(OUT)  :: dtpro
REAL (dp), INTENT(IN)   :: a(ldima)
REAL (dp), INTENT(IN)   :: b(ldimb)

! Local variables

REAL (dp)  :: eps
INTEGER    :: i, k, kk, ng4, ng4r, ng8, ng8r, ng16, ng16r, ng32, ng32r
!
eps=ten**(-maxexp)
overfl=.false.
!
dtpro=zero
IF (overch) THEN
  DO  i=1,n
    IF (LOG10(ABS(a(i))+eps) + LOG10(ABS(b(i))+eps) > maxexp) THEN
      overfl=.true.
      dtpro=SIGN(ten**maxexp,a(i))*SIGN(one,b(i))
      IF (output > 2 .AND. (.NOT.wrnsup)) THEN
        WRITE (nunit,10)
        10 FORMAT ('   *', t74, '*')
        WRITE (nunit,20) dtpro
        20 FORMAT ('   *    WARNING: TO AVOID OVERFLOW,',   &
                   ' INNER PRODUCT SET TO ', g12.3, t74, '*')
      END IF
      RETURN
    END IF
    dtpro=dtpro + a(i)*b(i)
  END DO
ELSE
!
!          SET NUMBER OF GROUPS OF EACH SIZE.
!
  ng32=n/32
  ng32r=n - 32*ng32
  ng16=ng32r/16
  ng16r=ng32r - 16*ng16
  ng8=ng16r/8
  ng8r=ng16r - 8*ng8
  ng4=ng8r/4
  ng4r=ng8r - 4*ng4
!
!          FIND INNER PRODUCT.
!
  k=0
  IF (ng32 > 0) THEN
    DO  kk=1,ng32
      k=k+32
      dtpro=dtpro + a(k-31)*b(k-31) + a(k-30)*b(k-30) + a(k-29)*b(k-29) +   &
            a(k-28)*b(k-28) + a(k-27)*b(k-27) + a(k-26)*b(k-26) +   &
            a(k-25)*b(k-25) + a(k-24)*b(k-24)
      dtpro=dtpro + a(k-23)*b(k-23) + a(k-22)*b(k-22) + a(k-21)*b(k-21) +   &
            a(k-20)*b(k-20) + a(k-19)*b(k-19) + a(k-18)*b(k-18) +   &
            a(k-17)*b(k-17) + a(k-16)*b(k-16)
      dtpro=dtpro + a(k-15)*b(k-15) + a(k-14)*b(k-14) + a(k-13)*b(k-13) +   &
            a(k-12)*b(k-12) + a(k-11)*b(k-11) + a(k-10)*b(k-10) +   &
            a(k-9)*b(k-9) + a(k-8)*b(k-8)
      dtpro=dtpro + a(k-7)*b(k-7) + a(k-6)*b(k-6) + a(k-5)*b(k-5) + a(k-4)  &
          *b(k-4) + a(k-3)*b(k-3) + a(k-2)*b(k-2) + a(k-1)*b(k-1) + a(k)*b(k)
    END DO
  END IF
  IF (ng16 > 0) THEN
    DO  kk=1,ng16
      k=k+16
      dtpro=dtpro + a(k-15)*b(k-15) + a(k-14)*b(k-14) + a(k-13)*b(k-13) +   &
            a(k-12)*b(k-12) + a(k-11)*b(k-11) + a(k-10)*b(k-10) +   &
            a(k-9)*b(k-9) + a(k-8)*b(k-8)
      dtpro=dtpro + a(k-7)*b(k-7) + a(k-6)*b(k-6) + a(k-5)*b(k-5) +   &
            a(k-4)*b(k-4) + a(k-3)*b(k-3) + a(k-2)*b(k-2) + a(k-1)*b(k-1) +  &
            a(k)*b(k)
    END DO
  END IF
  IF (ng8 > 0) THEN
    DO  kk=1,ng8
      k=k+8
      dtpro=dtpro + a(k-7)*b(k-7) + a(k-6)*b(k-6) + a(k-5)*b(k-5) +   &
            a(k-4)*b(k-4) + a(k-3)*b(k-3) + a(k-2)*b(k-2) + a(k-1)*b(k-1) + &
            a(k)*b(k)
    END DO
  END IF
  IF (ng4 > 0) THEN
    DO  kk=1,ng4
      k=k+4
      dtpro=dtpro + a(k-3)*b(k-3) + a(k-2)*b(k-2) + a(k-1)*b(k-1) + a(k)*b(k)
    END DO
  END IF
  IF (ng4r > 0) THEN
    DO  kk=1,ng4r
      k=k+1
      dtpro=dtpro + a(k)*b(k)
    END DO
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE INNERP.
!
END SUBROUTINE innerp



SUBROUTINE jaccd (n, nunit, output, epsmch, fvecj1, fvecj2,  &
                  jacfdm, scalex, xc, fvecev)
!
!    FEB. 11, 1991
!
!    THIS SUBROUTINE EVALUATES THE JACOBIAN USING CENTRAL DIFFERENCES.
!
!    FVECJ1 AND FVECJ2 ARE TEMPORARY VECTORS TO HOLD THE
!    RESIDUAL VECTORS FOR THE CENTRAL DIFFERENCE CALCULATION.
!

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: fvecj1(n)
REAL (dp), INTENT(OUT)     :: fvecj2(n)
REAL (dp), INTENT(OUT)     :: jacfdm(n,n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

LOGICAL    :: overfl
INTEGER    :: j
REAL (dp)  :: curtep, deltaj, tempj
!
overfl=.false.
curtep=epsmch**0.33
!
DO  j=1,n
  deltaj=(curtep)*SIGN((MAX(ABS(xc(j)),one/scalex(j))),xc(j))
  tempj=xc(j)
  xc(j)=xc(j) + deltaj
!
!          NOTE: THIS STEP IS FOR FLOATING POINT ACCURACY ONLY.
!
  deltaj=xc(j) - tempj
!
  CALL fvecev (overfl, n, fvecj1, xc)
  xc(j)=tempj-deltaj
  CALL fvecev (overfl, n, fvecj2, xc)
  IF (overfl .AND. output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,10)
    10 FORMAT ('   *', t74, '*')
    WRITE (nunit,20)
    20 FORMAT ('   *    WARNING: OVERFLOW IN FUNCTION VECTOR IN "JACCD"', &
               t74, '*')
  END IF
  jacfdm(1:n,j)=(fvecj1(1:n) - fvecj2(1:n))/(two*deltaj)
  xc(j)=tempj
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE JACCD.
!
END SUBROUTINE jaccd



SUBROUTINE jacfd (jactyp, n, nunit, output, epsmch, boundl,  &
                  boundu, fvecc, fvecj1, jacfdm, scalex, xc, fvecev)

! N.B. Argument WV3 has been removed.

!
!   FEB. 15, 1991
!
!   THIS SUBROUTINE EVALUATES THE JACOBIAN USING ONE-SIDED FINITE DIFFERENCES.
!
!   JACTYP "1" SIGNIFIES FORWARD DIFFERENCES
!   JACTYP "2" SIGNIFIES BACKWARD DIFFERENCES
!
!   FVECJ1 IS A TEMPORARY VECTOR WHICH STORES THE RESIDUAL
!   VECTOR FOR THE FINITE DIFFERENCE CALCULATION.
!

INTEGER, INTENT(IN)        :: jactyp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN)      :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: fvecj1(n)
REAL (dp), INTENT(OUT)     :: jacfdm(n,n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: xc(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev

  SUBROUTINE jacev (overfl, n, jac, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: jac(n,n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE jacev
END INTERFACE

! Local variables

REAL (dp)  :: deltaj, sqrtep, tempj
LOGICAL    :: overfl
INTEGER    :: j

!
sqrtep=SQRT(epsmch)
!
!       FINITE-DIFFERENCE CALCULATION BY COLUMNS.
!
DO  j=1,n
!
!          DELTAJ IS THE STEP SIZE - IT IS ALWAYS POSITIVE.
!
  deltaj=sqrtep*MAX(ABS(xc(j)), one/scalex(j))
  tempj=xc(j)
!  temporary storage of xc(j)
  IF (jactyp == 1) THEN
    IF (xc(j)+deltaj <= boundu(j)) THEN
!
!                STEP WITHIN BOUNDS - COMPLETE FORWARD DIFFERENCE.
!
      xc(j)=xc(j) + deltaj
      deltaj=xc(j) - tempj
      CALL fordif (overfl, j, n, deltaj, fvecc, fvecj1, jacfdm, xc, fvecev)
    ELSE
!
!                STEP WOULD VIOLATE BOUNDU - TRY BACKWARD DIFFERENCE.
!
      IF (xc(j)-deltaj >= boundl(j)) THEN
        xc(j)=xc(j) - deltaj
        CALL bakdif (overfl, j, n, deltaj, tempj, fvecc, fvecj1,  &
                     jacfdm, xc, fvecev)
      ELSE
!
!             STEP WOULD ALSO VIOLATE BOUNDL - IF THE DIFFERENCE IN THE
!             BOUNDS, (BOUNDU-BOUNDL), IS GREATER THAN DELTAJ CALCULATE
!             THE FUNCTION VECTOR AT EACH BOUND AND USE THIS DIFFERENCE -
!             THIS REQUIRES ONE EXTRA FUNCTION EVALUATION.
!             THE CURRENT FVECC IS STORED IN WV3, THEN REPLACED.
!
        IF (boundu(j)-boundl(j) >= deltaj) THEN
          CALL bnddif (overfl, j, n, epsmch, boundl, boundu,  &
                       fvecc, fvecj1, jacfdm, xc, fvecev)
          IF (output > 2 .AND. (.NOT.wrnsup) .AND. (.NOT.overfl)) THEN
            WRITE (nunit,10)
            10 FORMAT ('   *', t74, '*')
            WRITE (nunit,20)
            20 FORMAT ('   *    WARNING: BOUNDS TOO CLOSE',   &
                       ' FOR 1-SIDED FINITE-DIFFERENCES', t74, '*')
            WRITE (nunit,30) j
            30 FORMAT ('   *', t16, 'LOWER AND UPPER BOUNDS',   &
                       ' USED FOR JACOBIAN COLUMN: ', i3, t74, '*')
            WRITE (nunit,40)
            40 FORMAT ('   *', t16, 'THIS REQUIRED ONE EXTRA',   &
                       ' FUNCTION EVALUATION', t74, '*')
          END IF
        ELSE
!
!                      BOUNDS ARE EXTREMELY CLOSE (BUT NOT EQUAL OR
!                      THE PROGRAM WOULD HAVE STOPPED IN INITCH).
!
          IF (output > 2 .AND. (.NOT.wrnsup) .AND. (.NOT.overfl)) THEN
            WRITE (nunit,10)
            WRITE (nunit,20)
            WRITE (nunit,50) j
            50 FORMAT ('   *', t16, 'BOUNDS ARE EXTREMELY',   &
                       ' CLOSE FOR COMPONENT: ', i3, t74, '*')
            WRITE (nunit,60)
            60 FORMAT ('   *', t16, 'FINITE DIFFERENCE',   &
                       ' JACOBIAN IS UNRELIABLE', t74, '*')
          END IF
          CALL bnddif (overfl, j, n, epsmch, boundl, boundu,  &
                       fvecc, fvecj1, jacfdm, xc, fvecev)
        END IF
      END IF
    END IF
  ELSE
    IF (xc(j)-deltaj >= boundl(j)) THEN
      xc(j)=xc(j) - deltaj
      CALL bakdif (overfl, j, n, deltaj, tempj, fvecc, fvecj1,  &
                   jacfdm, xc, fvecev)
    ELSE
      IF (xc(j)+deltaj <= boundu(j)) THEN
        xc(j)=xc(j) + deltaj
        CALL fordif (overfl, j, n, deltaj, fvecc, fvecj1, jacfdm, xc, fvecev)
      ELSE
        IF (boundu(j)-boundl(j) >= deltaj) THEN
          CALL bnddif (overfl, j, n, epsmch, boundl, boundu,  &
                       fvecc, fvecj1, jacfdm, xc, fvecev)
          IF (output > 2 .AND. (.NOT.wrnsup) .AND. (.NOT.overfl)) THEN
            WRITE (nunit,10)
            WRITE (nunit,20)
            WRITE (nunit,30) j
            WRITE (nunit,40)
          END IF
        ELSE
          CALL bnddif (overfl, j, n, epsmch, boundl, boundu,  &
                       fvecc, fvecj1, jacfdm, xc, fvecev)
          jacfdm(1:n,j)=zero
          IF (output > 2 .AND. (.NOT.wrnsup) .AND. (.NOT.overfl)) THEN
            WRITE (nunit,10)
            WRITE (nunit,20)
            WRITE (nunit,50) j
            WRITE (nunit,60)
          END IF
        END IF
      END IF
    END IF
  END IF
  IF (overfl .AND. output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,10)
    WRITE (nunit,80)
    80 FORMAT ('   *    WARNING: OVERFLOW IN FUNCTION',   &
               ' VECTOR IN SUBROUTINE JACFD', t74, '*')
    WRITE (nunit,50) j
    jacfdm(1:n,j)=zero
  END IF
  xc(j)=tempj
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE JACFD.
!
END SUBROUTINE JACFD



SUBROUTINE jacobi (checkj, jacerr, overfl, jactyp, n, nunit, output, epsmch, &
                   fdtolj, boundl, boundu, fvecc, fvecj1, fvecj2, jac,   &
                   jacfdm, scalex, xc, fvecev, jacev)

! N.B. Argument WV3 has been removed.

!
!    APR. 13, 1991
!
!    THIS SUBROUTINE EVALUATES THE JACOBIAN.  IF CHECKJ IS TRUE
!    THEN THE ANALYTICAL JACOBIAN IS CHECKED NUMERICALLY.
!
!    JACEV IS A USER-SUPPLIED ANALYTICAL JACOBIAN USED ONLY IF JACTYP=0.
!    THE JACOBIAN NAME MAY BE CHANGED BY USING THE EXTERNAL STATEMENT
!    IN THE MAIN DRIVER.
!
!    JACFD ESTIMATES THE JACOBIAN USING FINITE DIFFERENCES:
!    FORWARD IF JACTYP=1 OR BACKWARD IF JACTYP=2.
!
!    JACCD ESTIMATES THE JACOBIAN USING CENTRAL DIFFERENCES.
!
!    IF THE ANALYTICAL JACOBIAN IS CHECKED THE FINITE DIFFERENCE
!    JACOBIAN IS STORED IN "JACFDM" AND THEN COMPARED.
!
!    FRSTER  INDICATES FIRST ERROR - USED ONLY TO SET BORDERS FOR OUTPUT
!    JACERR  FLAG TO INDICATE TO THE CALLING PROGRAM AN ERROR
!            IN THE ANALYTICAL JACOBIAN
!

LOGICAL, INTENT(IN)        :: checkj
LOGICAL, INTENT(OUT)       :: jacerr
LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: jactyp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: fdtolj
REAL (dp), INTENT(IN OUT)  :: boundl(n)
REAL (dp), INTENT(IN OUT)  :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: fvecj1(n)
REAL (dp), INTENT(IN OUT)  :: fvecj2(n)
REAL (dp), INTENT(OUT)     :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: jacfdm(n,n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: xc(n)

! EXTERNAL fvecev,jacev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev

  SUBROUTINE jacev (overfl, n, jac, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: jac(n,n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE jacev
END INTERFACE

! Local variables

LOGICAL    :: frster
INTEGER    :: i, j
!
frster=.true.
jacerr=.false.
overfl=.false.
!
IF (jactyp == 0) THEN
  CALL jacev (overfl, n, jac, xc)
ELSE IF (jactyp == 1 .OR. jactyp == 2) THEN
  CALL jacfd (jactyp, n, nunit, output, epsmch, boundl, boundu,  &
              fvecc, fvecj1, jac, scalex, xc, fvecev)
ELSE
  CALL jaccd (n, nunit, output, epsmch, fvecj1, fvecj2, jac,  &
              scalex, xc, fvecev)
END IF
!
IF (jactyp == 0 .AND. checkj) THEN
!
!          NOTE: JACTYP=0 SENT TO JACFD PRODUCES A FORWARD
!                DIFFERENCE ESTIMATE OF THE JACOBIAN.
!
  CALL jacfd (jactyp, n, nunit, output, epsmch, boundl, boundu,  &
              fvecc, fvecj1, jacfdm, scalex, xc, fvecev)
  DO  j=1,n
    DO  i=1,n
      IF (ABS((jacfdm(i,j)-jac(i,j))/MAX(ABS(jac(i,j)),  &
            ABS(jacfdm(i,j)),one)) > fdtolj) THEN
        jacerr=.true.
        IF (output >= 0) THEN
          IF (frster) THEN
            frster=.false.
            WRITE (nunit,10)
            10 FORMAT (t3, 72('*'))
          END IF
          WRITE (nunit,20)
          20 FORMAT ('   *', t74, '*')
          WRITE (nunit,30) i,j
          30 FORMAT ('   *    CHECK JACOBIAN TERM (', i3, ',', i3, ')', t74, '*')
          WRITE (nunit,20)
          WRITE (nunit,40) jac(i,j)
          40 FORMAT ('   *    ANALYTICAL DERIVATIVE IS ', g12.3, t74, '*')
          WRITE (nunit,50) jacfdm(i,j)
          50 FORMAT ('   *     NUMERICAL DERIVATIVE IS ', g12.3, t74, '*')
          WRITE (nunit,20)
          WRITE (nunit,10)
        END IF
      END IF
    END DO
  END DO
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE JACOBI.
!
END SUBROUTINE jacobi



SUBROUTINE jacrot (overfl, i, maxexp, n, arot, brot, epsmch, a, jac)
!
!    FEB. 11, 1991
!
!    JACOBI (OR GIVENS) ROTATION.
!

LOGICAL, INTENT(IN OUT)    :: overfl
INTEGER, INTENT(IN)        :: i
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: arot
REAL (dp), INTENT(IN)      :: brot
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)

! Local variables

REAL (dp)  :: c, denom, hold(2), s, w, y
INTEGER    :: j, ldhold

IF (arot == zero) THEN
  c=zero
  s=-SIGN(one,brot)
ELSE
  hold(1)=arot
  hold(2)=brot
  ldhold=2
  CALL twonrm (overfl, maxexp, ldhold, epsmch, denom, hold)
  c=arot/denom
  s=-brot/denom
END IF
DO  j=i,n
  y=a(i,j)
  w=a(i+1,j)
  a(i,j)=c*y - s*w
  a(i+1,j)=s*y + c*w
END DO
DO  j=1,n
  y=jac(i,j)
  w=jac(i+1,j)
  jac(i,j)=c*y - s*w
  jac(i+1,j)=s*y + c*w
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE JACROT.
!
END SUBROUTINE jacrot



SUBROUTINE line (abort, absnew, deuflh, geoms, newton, overch, overfl,   &
                 qnfail, qrsing, restrt, sclfch, sclxch, acpcod, acptcr,  &
                 contyp, isejac, itnum, jupdm, maxexp, maxlin, mgll, mnew, &
                 n, narmij, nfunc, nunit, output, qnupdm, stopcr, trmcod,  &
                 alpha, confac, epsmch, fcnmax, fcnnew, fcnold, lam0,   &
                 maxstp, newlen, sbrnrm, sigma, a, boundl, boundu, delf,   &
                 ftrack, fvec, h, hhpi, jac, rdiag, rhs, s, sbar, scalef,  &
                 scalex, sn, strack, xc, xplus, fvecev)

! N.B. Argument WV2 has been removed.

!
!    SEPT. 9, 1991
!

LOGICAL, INTENT(OUT)       :: abort
LOGICAL, INTENT(IN OUT)    :: absnew
LOGICAL, INTENT(IN)        :: deuflh
LOGICAL, INTENT(IN OUT)    :: geoms
LOGICAL, INTENT(IN)        :: newton
LOGICAL, INTENT(IN OUT)    :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(OUT)       :: qnfail
LOGICAL, INTENT(IN OUT)    :: qrsing
LOGICAL, INTENT(IN OUT)    :: restrt
LOGICAL, INTENT(IN OUT)    :: sclfch
LOGICAL, INTENT(IN OUT)    :: sclxch
INTEGER, INTENT(IN OUT)    :: acpcod
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(IN OUT)    :: contyp
INTEGER, INTENT(IN)        :: isejac
INTEGER, INTENT(IN OUT)    :: itnum
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN OUT)    :: maxlin
INTEGER, INTENT(IN)        :: mgll
INTEGER, INTENT(IN)        :: mnew
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: narmij
INTEGER, INTENT(IN OUT)    :: nfunc
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: qnupdm
INTEGER, INTENT(IN OUT)    :: stopcr
INTEGER, INTENT(IN OUT)    :: trmcod
REAL (dp), INTENT(IN OUT)  :: alpha
REAL (dp), INTENT(IN OUT)  :: confac
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(OUT)     :: fcnmax
REAL (dp), INTENT(OUT)     :: fcnnew
REAL (dp), INTENT(IN)      :: fcnold
REAL (dp), INTENT(IN)      :: lam0
REAL (dp), INTENT(IN)      :: maxstp
REAL (dp), INTENT(IN)      :: newlen
REAL (dp), INTENT(IN OUT)  :: sbrnrm
REAL (dp), INTENT(IN OUT)  :: sigma
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN)      :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(IN)      :: ftrack(0:mgll-1)
REAL (dp), INTENT(IN OUT)  :: fvec(n)
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: rdiag(n)
REAL (dp), INTENT(IN OUT)  :: rhs(n)
REAL (dp), INTENT(IN OUT)  :: s(n)
REAL (dp), INTENT(IN OUT)  :: sbar(n)
REAL (dp), INTENT(IN OUT)  :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: sn(n)
REAL (dp), INTENT(IN)      :: strack(0:mgll-1)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(OUT)     :: xplus(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

REAL (dp), SAVE  :: delfts, lambda, mu, newmax, norm, nrmpre
LOGICAL, SAVE    :: convio
INTEGER          :: i, j
REAL (dp)        :: wv2(n)

!
convio=.false.
overfl=.false.
abort = .FALSE.   ! Added by AJM
qnfail = .FALSE.  ! Added by AJM
IF (newton .OR. absnew) THEN
  IF (newton) THEN
!
!          FIND NEXT ITERATE FOR PURE NEWTON'S METHOD.
!
    xplus(1:n)=xc(1:n) + sn(1:n)
  ELSE
!
!       FIND NEXT ITERATE FOR "ABSOLUTE" NEWTON'S METHOD.
!       IF COMPONENT I WOULD BE OUTSIDE ITS BOUND THEN TAKE ABSOLUTE VALUE OF
!       THE VIOLATION AND GO THIS DISTANCE INTO THE FEASIBLE REGION.
!       ENSURE THAT THIS REFLECTION OFF ONE BOUND DOES NOT VIOLATE THE OTHER.
!
    DO  i=1,n
      wv2(i)=zero
      IF (sn(i) >= zero) THEN
        IF (xc(i)+sn(i) > boundu(i)) THEN
          convio=.true.
          wv2(i)=two
          xplus(i)=MAX(two*boundu(i)-xc(i)-sn(i), boundl(i))
        ELSE
          xplus(i)=xc(i) + sn(i)
        END IF
      ELSE
        IF (xc(i)+sn(i) < boundl(i)) THEN
          convio=.true.
          wv2(i)=-two
          xplus(i)=MIN(two*boundl(i)-xc(i)-sn(i), boundu(i))
        ELSE
          xplus(i)=xc(i) + sn(i)
        END IF
      END IF
    END DO
  END IF
  IF (convio .AND. output > 4) THEN
    WRITE (nunit,30)
    30 FORMAT ('   *', t74, '*')
    WRITE (nunit,40)
    40 FORMAT ('   *       CONSTRAINT VIOLATORS IN ABSOLUTE',   &
               ' NEWTON''S METHOD', t74, '*'/ '   *', t74, '*')
    WRITE (nunit,50)
    50 FORMAT ('   *       COMPONENT  PROPOSED POINT  VIOLATED ',  &
               'BOUND  FEASIBLE VALUE', t74, '*'/ '   *', t74, '*')
    DO  i=1,n
      IF (wv2(i) > zero) THEN
        WRITE (nunit,60) i,xc(i)+sn(i),boundu(i),xplus(i)
      ELSE IF (wv2(i) < zero) THEN
        WRITE (nunit,60) i,xc(i)+sn(i),boundl(i),xplus(i)
      END IF
      60 FORMAT ('   *', t10, i6, t21, g12.3, '    ', g12.3, '    ', g12.3, t74, '* ')
    END DO
  END IF
  CALL fvecev (overfl, n, fvec, xplus)
  nfunc=nfunc + 1
  IF (overfl) THEN
    overfl=.false.
    fcnnew=10.0**maxexp
    IF (output > 2) THEN
      WRITE (nunit,30)
      WRITE (nunit,80)
      80 FORMAT ('   *       POTENTIAL OVERFLOW DETECTED',   &
                 ' IN FUNCTION EVALUATION', t74, '*')
    END IF
    GO TO 90
  END IF
  CALL fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef)
!
!       RETURN FROM PURE NEWTON'S METHOD - OTHERWISE CONDUCT LINE SEARCH.
!
  90 RETURN
END IF
IF (output > 3) THEN
  WRITE (nunit,30)
  WRITE (nunit,30)
  WRITE (nunit,100)
  100 FORMAT ('   *    SUMMARY OF LINE SEARCH', t74, '*')
  WRITE (nunit,30)
END IF
!
!       SHORTEN NEWTON STEP IF LENGTH IS GREATER THAN MAXSTP.
!
IF (newlen > maxstp) THEN
  IF (output > 3) THEN
    WRITE (nunit,30)
    WRITE (nunit,110)
    110 FORMAT ('   *       LENGTH OF NEWTON STEP SHORTENED TO MAXSTP',  &
                t74, '*')
  END IF
  sn(1:n)=sn(1:n)*maxstp/newlen
END IF
!
!       CHECK DIRECTIONAL DERIVATIVE (MAGNITUDE AND SIGN).
!
CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output, delfts, delf, sn)
IF (overfl) THEN
  overfl=.false.
  IF (output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,30)
    WRITE (nunit,130) delfts
    130 FORMAT ('   *    WARNING: DIRECTIONAL DERIVATIVE, DELFTS, SET TO ', &
                g12.3, t74, '*')
  END IF
END IF
!
!       REVERSE SEARCH DIRECTION IF DIRECTIONAL DERIVATIVE IS
!       POSITIVE.
!
IF (delfts > zero) THEN
  sn(1:n)=-sn(1:n)
  delfts=-delfts
  IF (output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,30)
    WRITE (nunit,150)
    150 FORMAT ('   *    WARNING: DIRECTIONAL DERIVATIVE IS',  &
                ' POSITIVE: DIRECTION REVERSED', t74, '*')
  END IF
END IF
!
!       OUTPUT INFORMATION.
!
IF (output > 3) THEN
  WRITE (nunit,160) delfts
  160 FORMAT ('   *       INNER PRODUCT OF DELF AND SN, DELFTS: .......',  &
              g12.3, t74, '*')
  IF (.NOT.sclxch) THEN
    WRITE (nunit,170) newlen
    170 FORMAT ('   *       LENGTH OF NEWTON STEP, NEWLEN: ..............', &
                g12.3, t74, '*')
  ELSE
    WRITE (nunit,180) newlen
    180 FORMAT ('   *       LENGTH OF SCALED NEWTON STEP, NEWLEN: .......', &
                g12.3, t74, '*')
  END IF
  WRITE (nunit,190) maxstp
  190 FORMAT ('   *       MAXIMUM STEP LENGTH ALLOWED, MAXSTP: ........',  &
              g12.3, t74, '*')
END IF
!
!       ESTABLISH INITIAL RELAXATION FACTOR.
!
IF (deuflh) THEN
!
!       AT FIRST STEP IN DAMPED NEWTON OR AFTER EXPLICIT JACOBIAN EVALUATION
!       IN QUASI-NEWTON OR IF THE STEP SIZE IS WITHIN STOPPING TOLERANCE BUT
!       STOPCR=3 THE LINE SEARCH IS STARTED AT LAMBDA=1.
!
  IF (isejac == 1 .OR. (trmcod == 1 .AND. stopcr == 3)) THEN
    lambda=lam0
  ELSE
    wv2(1:n)=(sbar(1:n) - sn(1:n))*scalex(1:n)
    CALL twonrm (overfl, maxexp, n, epsmch, norm, wv2)
!
!             PREVENT DIVIDE BY ZERO IF NORM IS ZERO (UNDERFLOWS).
!
    IF (norm < epsmch) THEN
!
!                START LINE SEARCH AT LAMBDA=LAM0, USE DUMMY MU.
!
      mu=ten
    ELSE
      mu=nrmpre*lambda/norm
    END IF
    IF (output > 4) THEN
      WRITE (nunit,30)
      WRITE (nunit,210) mu
      210 FORMAT ('   *       DEUFLHARD TEST RATIO, MU:  ', g11.3, t74, '*')
    END IF
!
!       SET INITIAL LAMBDA DEPENDING ON MU.   THIS IS A MODIFICATION OF
!       DEUFLHARD'S METHOD WHERE THE CUTOFF VALUE WOULD BE 0.7 FOR LAM0=1.0.
!
    IF (mu > lam0/ten) THEN
      lambda=lam0
    ELSE
      lambda=lam0/ten
    END IF
  END IF
ELSE
  lambda=lam0
END IF
!
!       STORE LENGTH OF NEWTON STEP.  IF NEWTON STEP LENGTH WAS
!       GREATER THAN MAXSTP IT WAS SHORTENED TO MAXSTP.
!
nrmpre=MIN(maxstp,newlen)
!
!       ESTABLISH FCNMAX AND NEWMAX FOR NONMONOTONIC LINE SEARCH.
!
newmax=newlen
fcnmax=fcnold
IF (isejac > narmij) THEN
  IF (isejac < narmij+mgll) THEN
    DO  j=1,mnew
      fcnmax=MAX(fcnmax,ftrack(j-1))
      newmax=MAX(newmax,strack(j-1))
    END DO
  ELSE
    DO  j=0,mnew
      fcnmax=MAX(fcnmax, ftrack(j))
      newmax=MAX(newmax, strack(j))
    END DO
  END IF
END IF
!
IF (output > 3) THEN
  WRITE (nunit,30)
  WRITE (nunit,30)
  IF (.NOT.sclxch) THEN
    WRITE (nunit,240)
    240 FORMAT ('   *       LINE SEARCH', t74, '*')
  ELSE
    WRITE (nunit,250)
    250 FORMAT ('   *       LINE SEARCH (X''S GIVEN IN UNSCALED UNITS)',  &
                t74, '*')
  END IF
END IF
!
!       CONDUCT LINE SEARCH.
!
CALL deufls (abort, deuflh, geoms, overch, overfl, qnfail, qrsing, restrt,  &
             sclfch, sclxch, acpcod, acptcr, contyp, itnum, jupdm, maxexp,  &
             maxlin, n, nfunc, nunit, output, qnupdm, stopcr, alpha, confac, &
             delfts, epsmch, fcnmax, fcnnew, fcnold, lambda, newmax, sbrnrm, &
             sigma, a, h, boundl, boundu, delf, fvec, hhpi, jac, rdiag, rhs, &
             s, sbar, scalef, scalex, sn, xc, xplus, fvecev)
RETURN
!
!       LAST CARD OF SUBROUTINE LINE.
!
END SUBROUTINE line



SUBROUTINE llfa (overch, overfl, sclfch, sclxch, isejac, maxexp, n, nunit,  &
                 output, epsmch, omega, a, delf, fvec, fvecc, jac, plee,  &
                 rdiag, s, scalef, scalex, xc, xplus)

! N.B. Arguments T, W & WV3 have been removed.
!
!    FEB. 23, 1991
!
!    THE LEE AND LEE QUASI-NEWTON METHOD IS APPLIED TO
!    THE FACTORED FORM OF THE JACOBIAN.
!
!    NOTE: T AND W ARE TEMPORARY WORKING VECTORS ONLY.
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(IN)        :: sclfch
LOGICAL, INTENT(IN)        :: sclxch
INTEGER, INTENT(IN)        :: isejac
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN)      :: omega
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: plee(n,n)
REAL (dp), INTENT(IN OUT)  :: rdiag(n)
REAL (dp), INTENT(OUT)     :: s(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN)      :: xplus(n)

! Local variables

REAL (dp)  :: denom, sqrtep, SUM, t(n), w(n), wv3(n)
LOGICAL    :: skipup
INTEGER    :: i, j
!
overfl=.false.
sqrtep=SQRT(epsmch)
skipup=.true.
!
DO  i=1,n
  a(i,i)=rdiag(i)
  s(i)=xplus(i) - xc(i)
END DO
!
!       R IS IN THE UPPER TRIANGLE OF A.
!
!       T=RS
!
CALL uvmul (n, n, n, n, a, s, t)
!
!       FORM PART OF NUMERATOR AND CHECK TO SEE IF A SIGNIFICANT
!       CHANGE WOULD BE MADE TO THE JACOBIAN.
!
DO  i=1,n
  CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output,  &
               sum, jac(1,i), t)
  w(i)=scalef(i)*(fvec(i) - fvecc(i)) - sum
!
!          TEST TO ENSURE VECTOR W IS NONZERO.  IF W(I)=0 FOR
!          ALL I THEN THE UPDATE IS SKIPPED - SKIPUP IS TRUE.
!
  IF (ABS(w(i)) > sqrtep*scalef(i)*(ABS(fvec(i)) + ABS(fvecc(i)))) THEN
    skipup=.false.
! update to be performed
  ELSE
    w(i)=zero
  END IF
END DO

IF (.NOT.skipup) THEN
!
!          T=Q^W   Q^ IS STORED IN JAC.
!
  CALL avmul (n, n, n, n, jac, w, t)
!
!          FIND DENOMINATOR; FORM W=S^P (P IS SYMMETRIC SO PS IS FOUND).
!
  CALL avmul (n, n, n, n, plee, s, w)
  IF (sclxch) THEN
!
!             SCALE W TO FIND DENOMINATOR.
!
    wv3(1:n)=w(1:n)*scalex(1:n)*scalex(1:n)
  ELSE
    CALL matcop (n, n, 1, 1, n, 1, w, wv3)
  END IF
  CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output, denom, wv3, s)
!
!          IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN.
!
  IF (overfl) THEN
    IF (output > 3) THEN
      WRITE (nunit,40)
      40 FORMAT ('   *', t74, '*')
      WRITE (nunit,50)
      50 FORMAT ('   *    WARNING: JACOBIAN NOT UPDATED',   &
                 ' TO AVOID OVERFLOW IN DENOMINATOR OF', t74, '*'/  &
                 '   *    LEE AND LEE UPDATE', t74, '*')
    END IF
    RETURN
  END IF
!
!          IF DENOM IS ZERO THE SOLVER IS PROBABLY NEAR SOLUTION -
!          AVOID OVERFLOW AND CONTINUE WITH SAME JACOBIAN.
!
  IF (denom == zero) RETURN
!
!          THE SCALED VERSION OF S IS TAKEN TO THE UPDATE.
!
  w(1:n)=w(1:n)*scalex(1:n)*scalex(1:n)/denom
!
!      UPDATE THE QR DECOMPOSITION USING A SERIES OF GIVENS (JACOBI) ROTATIONS.
!
  CALL qrupda (overfl, maxexp, n, epsmch, a, jac, t, w)
!
!      RESET RDIAG AS DIAGONAL OF CURRENT R WHICH IS IN THE UPPER TRIANGLE OF A.
!
  DO  i=1,n
    rdiag(i)=a(i,i)
  END DO
!
!          UPDATE P MATRIX
!
  denom=omega**(isejac+2) + denom
  plee(1,1)=plee(1,1) - wv3(1)*wv3(1)/denom
  DO  j=2,n
    DO  i=1,j-1
      plee(i,j)=plee(i,j) - wv3(i)*wv3(j)/denom
      plee(j,i)=plee(i,j)
    END DO
    plee(j,j)=plee(j,j) - wv3(j)*wv3(j)/denom
  END DO
END IF
!
!       UPDATE THE GRADIENT VECTOR, DELF.
!
!       DELF = (QR)^F = R^Q^F = R^JAC F
!
IF (sclfch) THEN
  w(1:n)=fvec(1:n)*scalef(1:n)
ELSE
  CALL matcop (n, n, 1, 1, n, 1, fvec, w)
END IF
CALL avmul (n, n, n, n, jac, w, t)
CALL utbmul (n, n, 1, 1, n, n, a, t, delf)
RETURN
!
!       LAST CARD OF SUBROUTINE LLFA.
!
END SUBROUTINE llfa



SUBROUTINE llun (overch, overfl, isejac, maxexp, n, nunit, output, epsmch,  &
                 omega, fvec, fvecc, jac, plee, s, scalex, xc, xplus)

! N.B. Argument WV1 has been removed.
!
!    FEB. 13, 1991
!
!    UPDATE THE JACOBIAN USING THE LEE AND LEE METHOD.
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(IN OUT)    :: overfl
INTEGER, INTENT(IN)        :: isejac
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN)      :: omega
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: plee(n,n)
REAL (dp), INTENT(OUT)     :: s(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN)      :: xplus(n)

! Local variables

REAL (dp)  :: denom, sqrtep, SUM, tempi, wv1(n)
INTEGER    :: i, j
!
sqrtep=SQRT(epsmch)
!
s(1:n)=(xplus(1:n) - xc(1:n))*scalex(1:n)
DO  i=1,n
  wv1(i) = DOT_PRODUCT( s(1:n), plee(1:n,i) )
END DO
CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output, denom, wv1, s)
!
!       IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN.
!
IF (overfl) THEN
  IF (output > 3) THEN
    WRITE (nunit,40)
    40 FORMAT ('   *', t74, '*')
    WRITE (nunit,50)
    50 FORMAT ('   *    WARNING: JACOBIAN NOT UPDATED',   &
               ' TO AVOID OVERFLOW IN DENOMINATOR OF', t74, '*')
    WRITE (nunit,60)
    60 FORMAT ('   *    LEE AND LEE UPDATE', t74, '*')
  END IF
  RETURN
END IF
!
!       IF DENOM IS ZERO THE SOLVER MUST BE VERY NEAR SOLUTION -
!       AVOID OVERFLOW AND CONTINUE WITH SAME JACOBIAN.
!
IF (denom == zero) RETURN
DO  i=1,n
  sum = DOT_PRODUCT( jac(i,1:n), (xplus(1:n) - xc(1:n)) )
  tempi=fvec(i) - fvecc(i) - sum
  IF (ABS(tempi) >= sqrtep*(ABS(fvec(i)) + ABS(fvecc(i)))) THEN
    tempi=tempi/denom
    DO  j=1,n
      jac(i,j)=jac(i,j) + tempi*wv1(j)*scalex(j)
    END DO
  END IF
END DO
!
!       UPDATE P MATRIX.
!
denom=omega**(isejac+2) + denom
plee(1,1)=plee(1,1) - wv1(1)*wv1(1)/denom
DO  j=2,n
  DO  i=1,j-1
    plee(i,j)=plee(i,j) - wv1(i)*wv1(j)/(denom*scalex(i)*scalex(j))
    plee(j,i)=plee(i,j)
  END DO
  plee(j,j)=plee(j,j) - wv1(j)*wv1(j)/denom
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE LLUN.
!
END SUBROUTINE llun



SUBROUTINE lsolv (overch, overfl, maxexp, n, nunit, output, l, b, rhs)
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES:
!
!           LB=RHS
!
!           WHERE    L     IS TAKEN FROM THE CHOLESKY DECOMPOSITION
!                    RHS   IS A GIVEN RIGHT HAND SIDE WHICH IS NOT OVERWRITTEN
!                    B     IS THE SOLUTION VECTOR
!
!    FRSTER IS USED FOR OUTPUT PURPOSES ONLY.
!

LOGICAL, INTENT(IN)     :: overch
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp, n, nunit, output
REAL (dp), INTENT(IN)   :: l(n,n)
REAL (dp), INTENT(OUT)  :: b(n)
REAL (dp), INTENT(IN)   :: rhs(n)

! Local variables

INTEGER    :: i, j, jstar
REAL (dp)  :: eps, maxlog, SUM, tmplog
LOGICAL    :: frster

!
frster=.true.
overfl=.false.
eps=ten**(-maxexp)
!
IF (overch) THEN
  IF (LOG10(ABS(rhs(1))+eps) - LOG10(ABS(l(1,1))+eps) > maxexp) THEN
    overfl=.true.
    b(1)=SIGN(ten**maxexp, rhs(1))*SIGN(one,l(1,1))
    IF (output > 3) THEN
      WRITE (nunit,10)
      10 FORMAT ('   *', t74, '*')
      WRITE (nunit,20) 1,b(1)
      20 FORMAT ('   *    WARNING: COMPONENT ', i3, ' SET TO ',  &
                 g11.3, t74, '*')
    END IF
    GO TO 30
  END IF
END IF
b(1)=rhs(1)/l(1,1)

30 DO  i=2,n
  IF (overch) THEN
!
!             CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW.
!
    maxlog=LOG10(ABS(rhs(i))+eps) - LOG10(ABS(l(i,i))+eps)
    jstar=0
    DO  j=1,i-1
      tmplog=LOG10(ABS(l(i,j))+eps) + LOG10(ABS(b(j))+eps) -  &
             LOG10(ABS(l(i,i))+eps)
      IF (tmplog > maxlog) THEN
        jstar=j
        maxlog=tmplog
      END IF
    END DO
!
!             IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
!             TERM WITH CORRECT SIGN.
!
    IF (maxlog > maxexp) THEN
      overfl=.true.
      IF (jstar == 0) THEN
        b(i)=SIGN(ten**maxexp, rhs(i))*SIGN(one,l(i,i))
      ELSE
        b(i)=-SIGN(ten**maxexp, l(i,jstar))*SIGN(one,b(jstar))*SIGN(one,l(i,i))
      END IF
      IF (frster) THEN
        frster=.false.
        WRITE (nunit,10)
      END IF
      IF (output > 3) WRITE (nunit,20) i,b(i)
      CYCLE
    END IF
  END IF
!
!          SUM FOR EACH TERM, ORDERING OPERATIONS TO MINIMIZE
!          POSSIBILITY OF OVERFLOW.
!
  sum=zero
  DO  j=1,i-1
    sum=sum + (MIN(ABS(l(i,j)),ABS(b(j)))/l(i,i))*(MAX(ABS(l(i,j)),  &
               ABS(b(j))))*SIGN(one,l(i,j))*SIGN(one,b(j))
  END DO
  b(i)=rhs(i)/l(i,i) - sum
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE LSOLV.
!
END SUBROUTINE lsolv



SUBROUTINE ltsolv (overch, overfl, maxexp, n, nunit, output, l, y, b)
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES:
!
!           L^Y=B
!
!           WHERE    L  IS TAKEN FROM THE CHOLESKY DECOMPOSITION
!                    B  IS A GIVEN RIGHT HAND SIDE WHICH IS NOT OVERWRITTEN
!                    Y  IS THE SOLUTION VECTOR
!
!    FRSTER IS USED FOR OUTPUT PURPOSES ONLY.
!

LOGICAL, INTENT(IN)     :: overch
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp, n, nunit, output
REAL (dp), INTENT(IN)   :: l(n,n)
REAL (dp), INTENT(OUT)  :: y(n)
REAL (dp), INTENT(IN)   :: b(n)

! Local variables

REAL (dp)  :: eps, maxlog, SUM, tmplog
LOGICAL    :: frster
INTEGER    :: i, j, jstar
!
frster=.true.
overfl=.false.
eps=ten**(-maxexp)
!
IF (overch) THEN
  IF (LOG10(ABS(b(n))+eps) - LOG10(ABS(l(n,n))+eps) > maxexp) THEN
    overfl=.true.
    y(n)=SIGN(ten**maxexp,b(n))*SIGN(one,l(n,n))
    IF (output > 3) THEN
      frster=.false.
      WRITE (nunit,10)
      10 FORMAT ('   *', t74, '*')
      WRITE (nunit,20) n,y(n)
      20 FORMAT ('   *    WARNING: COMPONENT ', i3, ' SET TO ',  &
                 g11.3, t74, '*')
    END IF
    GO TO 30
  END IF
END IF
y(n)=b(n)/l(n,n)

30 DO  i=n-1,1,-1
  IF (overch) THEN
!
!             CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW.
!
    maxlog=LOG10(ABS(b(i))+eps) - LOG10(ABS(l(i,i))+eps)
    jstar=0
    DO  j=i+1,n
      tmplog=LOG10(ABS(l(j,i))+eps) + LOG10(ABS(y(j))+eps) -  &
             LOG10(ABS(l(i,i))+eps)
      IF (tmplog > maxlog) THEN
        jstar=j
        maxlog=tmplog
      END IF
    END DO
!
!             IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
!             TERM WITH CORRECT SIGN.
!
    IF (maxlog > maxexp) THEN
      overfl=.true.
      IF (jstar == 0) THEN
        y(i)=SIGN(ten**maxexp,b(i))*SIGN(one,l(i,i))
      ELSE
        y(i)=-SIGN(ten**maxexp,l(jstar,i))*SIGN(one,y(jstar))*  &
            SIGN(one,l(i,i))
      END IF
      IF (frster) THEN
        frster=.false.
        WRITE (nunit,10)
      END IF
      IF (output > 3) WRITE (nunit,20) i,y(i)
      CYCLE
    END IF
  END IF
!
!          SUM FOR EACH TERM ORDERING OPERATIONS TO MINIMIZE
!          POSSIBILITY OF OVERFLOW.
!
  sum=zero
  DO  j=i+1,n
    sum=sum + (MIN(ABS(l(j,i)),ABS(y(j)))/l(i,i))*(MAX(ABS(l(j,i)),  &
               ABS(y(j))))*SIGN(one,l(j,i))*SIGN(one,y(j))
  END DO
  y(i)=b(i)/l(i,i) - sum
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE LTSOLV.
!
END SUBROUTINE ltsolv



SUBROUTINE matcop (nradec, nraact, ncadec, ncaact, nrbdec, ncbdec, amat, bmat)
!
!    SEPT. 15, 1991
!
!    COPY A CONTINGUOUS RECTANGULAR PORTION OF ONE MATRIX
!    INTO ANOTHER (ELEMENT (1,1) MUST BE INCLUDED).
!
!    NRADEC IS 1ST DIMENSION OF AMAT, NRAACT IS LIMIT OF 1ST INDEX
!    NCADEC IS 2ND DIMENSION OF AMAT, NCAACT IS LIMIT OF 2ND INDEX
!    NRBDEC IS 1ST DIMENSION OF BMAT
!    NCBDEC IS 2ND DIMENSION OF BMAT
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncadec
INTEGER, INTENT(IN)     :: ncaact
INTEGER, INTENT(IN)     :: nrbdec
INTEGER, INTENT(IN)     :: ncbdec
REAL (dp), INTENT(IN)   :: amat(nradec,ncadec)
REAL (dp), INTENT(OUT)  :: bmat(nrbdec,ncbdec)

! Local variables

INTEGER    :: j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r
!
!       FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
ncc32=nraact/32
ncc32r=nraact - 32*ncc32
ncc16=ncc32r/16
ncc16r=ncc32r - 16*ncc16
ncc8=ncc16r/8
ncc8r=ncc16r - 8*ncc8
ncc4=ncc8r/4
ncc4r=ncc8r - 4*ncc4
DO  j=1,ncaact
!
!          COPY ENTRIES INTO MATRIX B BY COLUMN.
!
  k=0
  IF (ncc32 > 0) THEN
    DO  kk=1,ncc32
      k=k+32
      bmat(k-31,j)=amat(k-31,j)
      bmat(k-30,j)=amat(k-30,j)
      bmat(k-29,j)=amat(k-29,j)
      bmat(k-28,j)=amat(k-28,j)
      bmat(k-27,j)=amat(k-27,j)
      bmat(k-26,j)=amat(k-26,j)
      bmat(k-25,j)=amat(k-25,j)
      bmat(k-24,j)=amat(k-24,j)
      bmat(k-23,j)=amat(k-23,j)
      bmat(k-22,j)=amat(k-22,j)
      bmat(k-21,j)=amat(k-21,j)
      bmat(k-20,j)=amat(k-20,j)
      bmat(k-19,j)=amat(k-19,j)
      bmat(k-18,j)=amat(k-18,j)
      bmat(k-17,j)=amat(k-17,j)
      bmat(k-16,j)=amat(k-16,j)
      bmat(k-15,j)=amat(k-15,j)
      bmat(k-14,j)=amat(k-14,j)
      bmat(k-13,j)=amat(k-13,j)
      bmat(k-12,j)=amat(k-12,j)
      bmat(k-11,j)=amat(k-11,j)
      bmat(k-10,j)=amat(k-10,j)
      bmat(k-9,j)=amat(k-9,j)
      bmat(k-8,j)=amat(k-8,j)
      bmat(k-7,j)=amat(k-7,j)
      bmat(k-6,j)=amat(k-6,j)
      bmat(k-5,j)=amat(k-5,j)
      bmat(k-4,j)=amat(k-4,j)
      bmat(k-3,j)=amat(k-3,j)
      bmat(k-2,j)=amat(k-2,j)
      bmat(k-1,j)=amat(k-1,j)
      bmat(k,j)=amat(k,j)
    END DO
  END IF
  IF (ncc16 > 0) THEN
    DO  kk=1,ncc16
      k=k+16
      bmat(k-15,j)=amat(k-15,j)
      bmat(k-14,j)=amat(k-14,j)
      bmat(k-13,j)=amat(k-13,j)
      bmat(k-12,j)=amat(k-12,j)
      bmat(k-11,j)=amat(k-11,j)
      bmat(k-10,j)=amat(k-10,j)
      bmat(k-9,j)=amat(k-9,j)
      bmat(k-8,j)=amat(k-8,j)
      bmat(k-7,j)=amat(k-7,j)
      bmat(k-6,j)=amat(k-6,j)
      bmat(k-5,j)=amat(k-5,j)
      bmat(k-4,j)=amat(k-4,j)
      bmat(k-3,j)=amat(k-3,j)
      bmat(k-2,j)=amat(k-2,j)
      bmat(k-1,j)=amat(k-1,j)
      bmat(k,j)=amat(k,j)
    END DO
  END IF
  IF (ncc8 > 0) THEN
    DO  kk=1,ncc8
      k=k+8
      bmat(k-7,j)=amat(k-7,j)
      bmat(k-6,j)=amat(k-6,j)
      bmat(k-5,j)=amat(k-5,j)
      bmat(k-4,j)=amat(k-4,j)
      bmat(k-3,j)=amat(k-3,j)
      bmat(k-2,j)=amat(k-2,j)
      bmat(k-1,j)=amat(k-1,j)
      bmat(k,j)=amat(k,j)
    END DO
  END IF
  IF (ncc4 > 0) THEN
    DO  kk=1,ncc4
      k=k+4
      bmat(k-3,j)=amat(k-3,j)
      bmat(k-2,j)=amat(k-2,j)
      bmat(k-1,j)=amat(k-1,j)
      bmat(k,j)=amat(k,j)
    END DO
  END IF
  IF (ncc4r > 0) THEN
    DO  kk=1,ncc4r
      k=k+1
      bmat(k,j)=amat(k,j)
    END DO
  END IF
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE MATCOP.
!
END SUBROUTINE matcop



SUBROUTINE matprt (nrowa, ncola, nrowpr, ncolpr, nunit, a)
!
!    FEB. 6, 1991
!
!    THIS SUBROUTINE PRINTS RECTANGULAR BLOCKS STARTING WITH ELEMENT A(1,1)
!    OF SIZE NROWPR BY NCOLPR FOR MATRIX A (WHICH HAS DECLARED SIZE NROWA BY
!    NCOLA).  THE MATRIX IS PRINTED AS A BLOCK FOR SIZES UP TO 5X5 OR BY
!    COLUMNS IF IT IS LARGER.
!
!    NROWA IS THE NUMBER OF DECLARED ROWS IN THE MATRIX
!    NCOLA IS THE NUMBER OF DECLARED COLUMNS IN THE MATRIX
!
!    NROWPR IS THE NUMBER OF ROWS TO BE PRINTED
!    NCOLPR IS THE NUMBER OF COLUMNS TO BE PRINTED
!
!    IF MATRIX PRINTING IS TO BE SUPPRESSED THEN LOGICAL
!    VARIABLE MATSUP MUST BE SET TO TRUE BEFORE THE CALL TO NNES.
!
INTEGER, INTENT(IN)    :: nrowa
INTEGER, INTENT(IN)    :: ncola
INTEGER, INTENT(IN)    :: nrowpr
INTEGER, INTENT(IN)    :: ncolpr
INTEGER, INTENT(IN)    :: nunit
REAL (dp), INTENT(IN)  :: a(nrowa,ncola)

! Local variables

INTEGER  :: i, j, k, limit
!
IF (matsup) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,20)
  20 FORMAT ('   *       MATRIX PRINTING SUPPRESSED', t74, '*')
  RETURN
END IF
!
!       FOR NCOLPR <= 5 WRITE MATRIX AS A WHOLE.
!
WRITE (nunit,10)
IF (ncolpr <= 5) THEN
  WRITE (nunit,30) (k,k=1,ncolpr)
  30 FORMAT (t74, '*'/ '   *  ', 5(i12:))
  WRITE (nunit,10)
  DO  i=1,nrowpr
    WRITE (nunit,40) i,(a(i,k),k=1,ncolpr)
    40 FORMAT (t74, '*'/ '   *   ', i3, 5(g12.3))
  END DO
ELSE
!
!          LIMIT IS THE NUMBER OF GROUPS OF 5 COLUMNS.
!
  limit=ncolpr/5
!
!          WRITE COMPLETE BLOCKS FIRST (LEFTOVERS LATER).
!
  DO  j=1,limit
    WRITE (nunit,60) (k,k=1+(j-1)*5,5+(j-1)*5)
    60 FORMAT ('   *  ', 5I12, t74, '*')
    WRITE (nunit,10)
    DO  i=1,nrowpr
      WRITE (nunit,70) i,(a(i,k),k=1+(j-1)*5,5+(j-1)*5)
      70 FORMAT ('   *   ', i3, 5g12.3, t74, '*')
    END DO
    WRITE (nunit,10)
  END DO
!
!          WRITE REMAINING ELEMENTS.
!
  WRITE (nunit,100) (k,k=5*limit+1,ncolpr)
  100 FORMAT ('   *', 4(i12), t74, '*')
  WRITE (nunit,10)
  DO  i=1,nrowpr
    WRITE (nunit,110) i,(a(i,k),k=5*limit+1,ncolpr)
    110 FORMAT ('   *   ', i3, 4g12.3, t74, '*')
  END DO
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE MATPRT.
!
END SUBROUTINE matprt



SUBROUTINE maxst (overfl, maxexp, n, nunit, output, epsmch,  &
                  maxstp, mstpf, scalex, xc)

! N.B. Argument WV1 has been removed.
!
!    FEB. 11, 1991
!
!    THIS SUBROUTINE ESTABLISHES A MAXIMUM STEP LENGTH BASED ON THE 2-NORMS
!    OF THE INITIAL ESTIMATES AND THE SCALING FACTORS MULTIPLIED BY A
!    FACTOR MSTPF.
!
!         MAXSTP=MSTPF*MAX{ NORM1 , NORM2 }
!
!              WHERE   MSTPF  USER-CHOSEN FACTOR (DEFAULT: 1000)
!                      NORM1  2-NORM OF SCALED STARTING ESTIMATES
!                      NORM2  2-NORM OF COMPONENT SCALING FACTORS
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nunit
INTEGER, INTENT(IN)     :: output
REAL (dp), INTENT(IN)   :: epsmch
REAL (dp), INTENT(OUT)  :: maxstp
REAL (dp), INTENT(IN)   :: mstpf
REAL (dp), INTENT(IN)   :: scalex(n)
REAL (dp), INTENT(IN)   :: xc(n)

! Local variables

REAL (dp)  :: norm1, norm2, wv1(n)
!
overfl=.false.
!
wv1(1:n)=scalex(1:n)*xc(1:n)
CALL twonrm (overfl, maxexp, n, epsmch, norm1, wv1)
IF (overfl) THEN
  maxstp=ten**maxexp
  IF (output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,20)
    20 FORMAT ('   *', t74, '*')
    WRITE (nunit,30) norm1
    30 FORMAT ('   *    WARNING: NORM OF SCALED INITIAL ESTIMATE SET TO ', &
               g12.3, t74, '*')
    WRITE (nunit,20)
    WRITE (nunit,40) maxstp
    40 FORMAT ('   *       MAXIMUM STEP SIZE, MAXSTP, SET TO ',  &
               g12.3, t74, '*')
  END IF
  RETURN
END IF
CALL twonrm (overfl, maxexp, n, epsmch, norm2, scalex)
IF (overfl) THEN
  maxstp=ten**maxexp
  IF (output > 2 .AND. (.NOT.wrnsup)) THEN
    WRITE (nunit,20)
    WRITE (nunit,50) norm2
    50 FORMAT ('   *    WARNING: NORM OF SCALING FACTORS SET TO ', g12.3, &
               t74, '*')
    WRITE (nunit,20)
    WRITE (nunit,40) maxstp
  END IF
  RETURN
END IF
maxstp=mstpf*MAX(norm1,norm2)
RETURN
!
!       LAST CARD OF SUBROUTINE MAXST.
!
END SUBROUTINE maxst



SUBROUTINE nersl (newton, restrt, sclfch, sclxch, acpcod, jupdm,  &
                  n, nunit, output, fcnnew, fvec, xplus)
!
!    SEPT. 2, 1991
!
!    THE RESULTS OF EACH ITERATION ARE PRINTED.
!

LOGICAL, INTENT(IN)    :: newton
LOGICAL, INTENT(IN)    :: restrt
LOGICAL, INTENT(IN)    :: sclfch
LOGICAL, INTENT(IN)    :: sclxch
INTEGER, INTENT(IN)    :: acpcod
INTEGER, INTENT(IN)    :: jupdm
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN)    :: nunit
INTEGER, INTENT(IN)    :: output
REAL (dp), INTENT(IN)  :: fcnnew
REAL (dp), INTENT(IN)  :: fvec(n)
REAL (dp), INTENT(IN)  :: xplus(n)

! Local variable

INTEGER  :: i
!
WRITE (nunit,10)
10 FORMAT ('   *', t74, '*')
WRITE (nunit,10)
IF (.NOT.sclxch) THEN
  WRITE (nunit,20)
  20 FORMAT ('   *    SUMMARY OF ITERATION RESULTS', t74, '*')
ELSE
  WRITE (nunit,30)
  30 FORMAT ('   *    SUMMARY OF ITERATION RESULTS (X''S GIVEN IN ',  &
             'UNSCALED UNITS)', t74, '*')
END IF
WRITE (nunit,10)
WRITE (nunit,40)
40 FORMAT ('   *        UPDATED ESTIMATES', t45, 'UPDATED FUNCTION VALUES', &
           t74, '*')
WRITE (nunit,10)
DO  i=1,n
  WRITE (nunit,50) i,xplus(i),i,fvec(i)
  50 FORMAT ('   *', t9, 'X(', i3, ') = ', g12.3, t45, 'F(', i3, ') = ',  &
             g12.3, t74, '*')
END DO
WRITE (nunit,10)
IF (.NOT.sclfch) THEN
  WRITE (nunit,70) fcnnew
  70 FORMAT ('   *      OBJECTIVE FUNCTION VALUE: ', g12.3, t74, '*')
ELSE
  WRITE (nunit,80) fcnnew
  80 FORMAT ('   *      SCALED OBJECTIVE FUNCTION VALUE: ', g12.3, t74, '*')
END IF
WRITE (nunit,10)
IF (output > 3 .AND. (.NOT.newton)) THEN
  WRITE (nunit,90) acpcod
  90 FORMAT ('   *      STEP ACCEPTANCE CODE, ACPCOD:',i9, t74, '*')
END IF
IF (restrt .AND. output >= 3 .AND. jupdm /= 0) THEN
  IF (output > 3) WRITE (nunit,10)
  WRITE (nunit,100)
  100 FORMAT ('   *      NOTE: JACOBIAN EVALUATED EXPLICITLY AT THIS STEP',  &
              t74, '*')
  WRITE (nunit,10)
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE NERSL.
!
END SUBROUTINE nersl



SUBROUTINE nestop (absnew, linesr, newton, sclfch, sclxch, acptcr, itnum,  &
                   n, nac1, nac2, nac12, nfunc, njetot, nunit, output,  &
                   stopcr, trmcod, fcnnew, ftol, nsttol, stpmax, stptol,  &
                   fvec, scalef, scalex, xc, xplus)
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE CHECKS TO SEE IF THE CONVERGENCE CRITERIA HAVE BEEN MET.
!

LOGICAL, INTENT(IN OUT)    :: absnew
LOGICAL, INTENT(IN)        :: linesr
LOGICAL, INTENT(IN OUT)    :: newton
LOGICAL, INTENT(IN)        :: sclfch
LOGICAL, INTENT(IN OUT)    :: sclxch
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(IN OUT)    :: itnum
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN OUT)    :: nac1
INTEGER, INTENT(IN OUT)    :: nac2
INTEGER, INTENT(IN OUT)    :: nac12
INTEGER, INTENT(IN OUT)    :: nfunc
INTEGER, INTENT(IN OUT)    :: njetot
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: stopcr
INTEGER, INTENT(IN OUT)    :: trmcod
REAL (dp), INTENT(IN OUT)  :: fcnnew
REAL (dp), INTENT(IN)      :: ftol
REAL (dp), INTENT(IN)      :: nsttol
REAL (dp), INTENT(IN OUT)  :: stpmax
REAL (dp), INTENT(IN)      :: stptol
REAL (dp), INTENT(IN)      :: fvec(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN OUT)  :: xplus(n)

! Local variables

REAL (dp)  :: MAX1, max2, ratio1, SUM
INTEGER    :: i
!
IF (output > 3) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,20)
  20 FORMAT ('   *    CONVERGENCE TESTING', t74, '*')
  WRITE (nunit,10)
END IF
!
!       IF THE NEWTON STEP WAS WITHIN TOLERANCE THEN TRMCOD IS 1.
!
IF (trmcod == 1) THEN
  IF (output > 3) THEN
    WRITE (nunit,10)
    IF (.NOT.sclxch) THEN
      WRITE (nunit,30) stpmax
      30 FORMAT ('   *', t9, 'MAXIMUM NEWTON STEP LENGTH STPMAX:',  &
                 g12.3, t74, '*')
    ELSE
      WRITE (nunit,40) stpmax
      40 FORMAT ('   *', t9, 'MAXIMUM SCALED NEWTON STEP LENGTH STPMAX:',  &
                 g12.3, t74, '*')
    END IF
    WRITE (nunit,50) nsttol
    50 FORMAT ('   *      FIRST CONVERGENCE CRITERION MET;  NSTTOL IS:',  &
               g12.3, t74, '*')
    WRITE (nunit,10)
!
!          SKIP CHECKING OTHER STEP SIZE CRITERION AS TRMCOD IS ALREADY 1.
!
    GO TO 100
  END IF
END IF
!
!       IF THE NEWTON STEP WAS NOT WITHIN TOLERANCE THEN, IF
!       STOPCR IS NOT EQUAL TO 2, THE SECOND STEP SIZE STOPPING
!       CRITERION MUST BE CHECKED.
!
IF (stopcr /= 2 .AND. trmcod /= 1) THEN
  MAX1=zero
  DO  i=1,n
    ratio1=(ABS(xplus(i)-xc(i)))/MAX(ABS(xplus(i)), one/scalex(i))
    MAX1=MAX(MAX1,ratio1)
    IF (output > 4) THEN
      WRITE (nunit,60) i,ratio1
      60 FORMAT ('   *', t9, 'RELATIVE STEP SIZE (', i3, ') = ', g12.3,  &
                 t74, '*')
    END IF
  END DO
  IF (output > 4) WRITE (nunit,10)
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      WRITE (nunit,80) MAX1,stptol
      80 FORMAT ('   *', t9, 'MAXIMUM STEP SIZE:', g12.3, '   STPTOL:',  &
                 g11.3, t74, '*')
    ELSE
      WRITE (nunit,90) MAX1,stptol
      90 FORMAT ('   *', t9, 'MAXIMUM RELATIVE STEP SIZE:', g12.3,   &
                 '   STPTOL:', g11.3, t74, '*')
    END IF
    WRITE (nunit,10)
  END IF
  IF (MAX1 < stptol) THEN
    trmcod=1
  END IF
END IF
!
!       NOTE: CONTINUATION AT 101 MEANS THAT TRMCOD WAS 1 ON ENTRY
!             SO THE STEP SIZE CRITERION ABOVE DID NOT NEED TO BE CHECKED.
!
!
!       THE SECOND STOPPING CRITERION IS CHECKED IF NEEDED.
!
100 IF (stopcr == 2 .OR. stopcr == 12 .OR. (stopcr == 3 .AND. trmcod == 1)) THEN
  max2=zero
  DO  i=1,n
    max2=MAX(max2,scalef(i)*ABS(fvec(i)))
    IF (output > 4) THEN
      IF (.NOT.sclfch) THEN
        WRITE (nunit,110) i,ABS(fvec(i))
        110 FORMAT ('   *', t9, 'ABSOLUTE FUNCTION VECTOR (', i3, ') = ',  &
             g12.3, t74, '*')
      ELSE
        WRITE (nunit,120) i,scalef(i)*ABS(fvec(i))
        120 FORMAT ('   *', t9, 'ABSOLUTE SCALED FUNCTION VECTOR (',  &
            i3,') = ', g12.3, t74, '*')
      END IF
    END IF
  END DO
  IF (output > 4) WRITE (nunit,10)
  IF (output > 3) THEN
    IF (.NOT.sclfch) THEN
      WRITE (nunit,140) max2,ftol
      140 FORMAT ('   *', t9, 'MAXIMUM ABSOLUTE FUNCTION:', g12.3,   &
                  '     FTOL:', g11.3, t74, '*')
    ELSE
      WRITE (nunit,150) max2,ftol
      150 FORMAT ('   *', t9, 'MAX ABSOLUTE SCALED FUNCTION:',  &
                  g12.3, '     FTOL:', g11.3, t74, '*')
    END IF
    WRITE (nunit,10)
  END IF
  IF (max2 < ftol) THEN
    IF (stopcr == 3 .AND. trmcod == 1) THEN
!
!                BOTH NEEDED STOPPING CRITERIA HAVE BEEN MET.
!
      trmcod=3
    ELSE IF (stopcr == 12 .AND. trmcod == 1) THEN
!
!                BOTH STOPPING CRITERIA HAVE BEEN MET ALTHOUGH
!                EITHER ONE WOULD BE SATISFACTORY.
!
      trmcod=12
    ELSE IF (stopcr == 2 .OR. stopcr == 12) THEN
      trmcod=2
    END IF
  ELSE IF (stopcr == 3) THEN
!
!             ONLY THE FIRST STOPPING CRITERION WAS MET - TRMCOD
!             MUST BE RESET FROM 1 BACK TO 0.
!
    trmcod=0
  END IF
END IF
!
!       PRINT FINAL RESULTS IF CONVERGENCE REACHED.
!
IF (trmcod > 0) THEN
  IF (output > 0) THEN
    IF (output == 1) WRITE (nunit,160)
    160 FORMAT (t3, 72('*'))
    WRITE (nunit,10)
    WRITE (nunit,170) trmcod
    170 FORMAT ('   *', t9, 'CONVERGENCE REACHED; TERMINATION CODE: ',  &
                '...............', i6, t74, '*')
!
!             ITERATION RESULTS NOT PRINTED IN NERSL.
!
    WRITE (nunit,10)
    WRITE (nunit,10)
    IF (sclfch .OR. sclxch) THEN
      WRITE (nunit,180)
      180 FORMAT ('   *', t32, 'UNSCALED RESULTS', t74, '*')
      WRITE (nunit,10)
    END IF
    WRITE (nunit,190)
    190 FORMAT ('   *          FINAL ESTIMATES', t32, 'FINAL FUNCTION VALUES', &
                t74, '*')
    WRITE (nunit,10)
    DO  i=1,n
      WRITE (nunit,200) i,xplus(i),i,fvec(i)
      200 FORMAT ('   *', t9, 'X(', i3, ') = ', g14.5, t45, 'F(', i3,   &
                  ') = ', g14.5, t74, '*')
    END DO
    WRITE (nunit,10)
    IF (sclfch) THEN
!
!                NEED UNSCALED OBJECTIVE FUNCTION.
!
      sum = DOT_PRODUCT( fvec(1:n), fvec(1:n) )
      WRITE (nunit,230) sum/two
    ELSE
      WRITE (nunit,230) fcnnew
    END IF
    230 FORMAT ('   *', t9, 'FINAL OBJECTIVE FUNCTION VALUE:', g12.3, t74, '*')
    IF (sclfch .OR. sclxch) THEN
      WRITE (nunit,10)
      WRITE (nunit,10)
      WRITE (nunit,240)
      240 FORMAT ('   *', t33, 'SCALED RESULTS', t74, '*'/ '   *', t74, '*'/ &
                  '   *          FINAL ESTIMATES', t48,   &
                  'FINAL FUNCTION VALUES ', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,200) i,scalex(i)*xplus(i),i,scalef(i)*fvec(i)
      END DO
      WRITE (nunit,10)
      WRITE (nunit,230) fcnnew
    END IF
  END IF
ELSE
  IF (output > 3) THEN
    WRITE (nunit,260)
    260 FORMAT ('   *', t9, 'CONVERGENCE NOT REACHED', t74, '*')
    WRITE (nunit,10)
  END IF
  RETURN
END IF
!
!       TERMINATION HAS BEEN REACHED.
!
IF (output > 0) THEN
  WRITE (nunit,10)
  WRITE (nunit,10)
  WRITE (nunit,270) itnum
  270 FORMAT ('   *', t9, 'TOTAL NUMBER OF ITERATIONS',   &
              ': ..........................', i6, t74, '*')
  IF (.NOT.newton .AND. (.NOT.absnew)) THEN
    IF (linesr) THEN
      WRITE (nunit,280) nfunc
      280 FORMAT ('   *', t9, 'TOTAL NUMBER OF LINE SEARCH ',   &
                  'FUNCTION EVALUATIONS: ....', i6, t74, '*')
    ELSE
      WRITE (nunit,290) nfunc
      290 FORMAT ('   *', t9, 'TOTAL NUMBER OF TRUST REGION',   &
                  ' FUNCTION EVALUATIONS: ...', i6, t74, '*')
    END IF
  END IF
  WRITE (nunit,300) njetot
  300 FORMAT ('   *', t9, 'TOTAL NUMBER OF EXPLICIT JACOBIAN',   &
              ' EVALUATIONS: .......', i6, t74, '*')
  WRITE (nunit,310) nfetot
  310 FORMAT ('   *', t9, 'TOTAL NUMBER OF FUNCTION',   &
              ' EVALUATIONS: ................', i6, t74, '*')
  WRITE (nunit,10)
  IF (.NOT.newton .AND. (.NOT.absnew) .AND. acptcr /= 1 .AND. output > 2) THEN
    WRITE (nunit,320) nac1
    320 FORMAT ('   *', t9, 'NUMBER OF STEPS ACCEPTED',   &
                ' BY FUNCTION VALUE ONLY: .....', i6, t74, '*')
    WRITE (nunit,330) nac2
    330 FORMAT ('   *', t9, 'NUMBER OF STEPS ACCEPTED',   &
                ' BY STEP SIZE VALUE ONLY: ....',i6, t74, '*')
    WRITE (nunit,340) nac12
    340 FORMAT ('   *', t9, 'NUMBER OF STEPS ACCEPTED',   &
                ' BY EITHER CRITERION: ........', i6, t74, '*')
    WRITE (nunit,10)
  END IF
  WRITE (nunit,160)
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE NESTOP.
!
END SUBROUTINE nestop



SUBROUTINE nnes (absnew, cauchy, deuflh, geoms, linesr, newton, overch,   &
                 acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, &
                 maxqns, mgll, minqns, n, narmij, niejev, njacch, njetot,  &
                 nunit, output, qnupdm, stopcr, supprs, trmcod, trupdm, alpha, &
                 confac, delta, delfac, epsmch, etafac, fcnnew, fdtolj, ftol, &
                 lam0, mstpf, nsttol, omega, ratiof, sigma, stptol, a, boundl, &
                 boundu, delf, fsave, ftrack, fvec, fvecc, h, hhpi, jac, plee, &
                 rdiag, s, sbar, scalef, scalex, sn, ssdhat, strack, vhat,  &
                 xc, xplus, xsave, help, fvecev, jacev)

! N.B. Arguments WV1, WV2, WV3 & WV4 have been removed.

!
!    FEB. 28, 1992
!

LOGICAL, INTENT(IN OUT)    :: absnew
LOGICAL, INTENT(IN OUT)    :: cauchy
LOGICAL, INTENT(IN OUT)    :: deuflh
LOGICAL, INTENT(IN OUT)    :: geoms
LOGICAL, INTENT(IN OUT)    :: linesr
LOGICAL, INTENT(IN OUT)    :: newton
LOGICAL, INTENT(IN OUT)    :: overch
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(IN OUT)    :: itsclf
INTEGER, INTENT(IN OUT)    :: itsclx
INTEGER, INTENT(IN)        :: jactyp
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN OUT)    :: maxexp
INTEGER, INTENT(IN)        :: maxit
INTEGER, INTENT(IN)        :: maxns
INTEGER, INTENT(IN)        :: maxqns
INTEGER, INTENT(IN)        :: mgll
INTEGER, INTENT(IN)        :: minqns
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: narmij
INTEGER, INTENT(IN)        :: niejev
INTEGER, INTENT(IN)        :: njacch
INTEGER, INTENT(OUT)       :: njetot
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN OUT)    :: output
INTEGER, INTENT(IN)        :: qnupdm
INTEGER, INTENT(IN OUT)    :: stopcr
INTEGER, INTENT(IN)        :: supprs
INTEGER, INTENT(OUT)       :: trmcod
INTEGER, INTENT(IN OUT)    :: trupdm
REAL (dp), INTENT(IN OUT)  :: alpha
REAL (dp), INTENT(IN OUT)  :: confac
REAL (dp), INTENT(IN OUT)  :: delta
REAL (dp), INTENT(IN OUT)  :: delfac
REAL (dp), INTENT(IN OUT)  :: epsmch
REAL (dp), INTENT(IN OUT)  :: etafac
REAL (dp), INTENT(IN OUT)  :: fcnnew
REAL (dp), INTENT(IN OUT)  :: fdtolj
REAL (dp), INTENT(IN OUT)  :: ftol
REAL (dp), INTENT(IN OUT)  :: lam0
REAL (dp), INTENT(IN OUT)  :: mstpf
REAL (dp), INTENT(IN)      :: nsttol
REAL (dp), INTENT(IN OUT)  :: omega
REAL (dp), INTENT(IN)      :: ratiof
REAL (dp), INTENT(IN OUT)  :: sigma
REAL (dp), INTENT(IN OUT)  :: stptol
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: boundl(n)
REAL (dp), INTENT(IN OUT)  :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(OUT)     :: fsave(n)
REAL (dp), INTENT(OUT)     :: ftrack(0:mgll-1)
REAL (dp), INTENT(IN OUT)  :: fvec(n)
REAL (dp), INTENT(OUT)     :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(OUT)     :: jac(n,n)
REAL (dp), INTENT(OUT)     :: plee(n,n)
REAL (dp), INTENT(IN OUT)  :: rdiag(n)
REAL (dp), INTENT(IN OUT)  :: s(n)
REAL (dp), INTENT(IN OUT)  :: sbar(n)
REAL (dp), INTENT(IN OUT)  :: scalef(n)
REAL (dp), INTENT(IN OUT)  :: scalex(n)
REAL (dp), INTENT(IN OUT)  :: sn(n)
REAL (dp), INTENT(IN OUT)  :: ssdhat(n)
REAL (dp), INTENT(OUT)     :: strack(0:mgll-1)
REAL (dp), INTENT(IN OUT)  :: vhat(n)
REAL (dp), INTENT(IN OUT)  :: xc(n)
REAL (dp), INTENT(IN OUT)  :: xplus(n)
REAL (dp), INTENT(OUT)     :: xsave(n)
CHARACTER (LEN=6), INTENT(IN)  :: help

! EXTERNAL fvecev,jacev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev

  SUBROUTINE jacev (overfl, n, jac, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: jac(n,n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE jacev
END INTERFACE

REAL (dp)  :: beta, caulen, delstr, delta0, fcnmax, fcnmin, fcnold, fcnpre, &
              fstore, maxstp, newlen, newmax, powtau, rellen, sbrnrm, sqrtz,  &
              stpmax, xnorm, wv1(n), wv2(n), wv4(n)
INTEGER    :: acpcod, acpstr, contyp, countr, i, isejac, itemp, itstr, j,  &
              maxlin, maxtrs, mnew, mold, nac1, nac2, nac12, nfail, nfestr,   &
              nfetot, nosubt, notrst, outtmp, retcod
LOGICAL    :: abort, checkj, frstdg, instop, jacerr, newtkn, overfl,  &
              qnfail, qrsing, restrt, savest, sclfch, sclxch
INTEGER    :: newstm = 0   ! Added by AJM

REAL (dp), PARAMETER  :: pt01 = 0.01_dp

!
!       PRINT HELP IF REQUESTED.
!
IF (help(1:4) /= 'NONE') THEN
  CALL olhelp (nunit, help)
  RETURN
END IF

qnfail=.false.
restrt=.true.
savest=.false.
countr=0
isejac=0
mnew=0
nac1=0
nac2=0
nac12=0
nfunc=1
njetot=0
outtmp=output
trmcod=0
retcod = 0    ! Added by AJM
!
!    ESTABLISH INITIAL FUNCTION VALUE AND CHECK FOR STARTING ESTIMATE WHICH
!    IS A SOLUTION.  ALSO, CHECK FOR INCOMPATIBILITIES IN INPUT PARAMETERS.
!
CALL initch (instop, linesr, newton, overfl, sclfch, sclxch, acptcr, contyp, &
             jactyp, jupdm, maxexp, n, nunit, output, qnupdm, stopcr,   &
             trupdm, epsmch, fcnold, ftol, boundl, boundu, fvecc, scalef,   &
             scalex, xc, fvecev)
!
!    IF A FATAL ERROR IS DETECTED IN INITCH RETURN TO MAIN PROGRAM.
!    NOTE: SOME INCOMPATIBILITIES ARE CORRECTED WITHIN INITCH AND EXECUTION
!    CONTINUES.   WARNINGS ARE GENERATED WITHIN INITCH.
!
IF (instop) RETURN
!
!       ESTABLISH MAXIMUM STEP LENGTH ALLOWED (USUALLY THIS IS MUCH
!       LARGER THAN ACTUAL STEP SIZES - IT IS ONLY TO PREVENT
!       EXCESSIVELY LARGE STEPS).  THE FACTOR MSTPF CONTROLS THE
!       MAGNITUDE OF MAXSTP AND IS SET BY THE USER (DEFAULT=1000).
!
CALL maxst (overfl, maxexp, n, nunit, output, epsmch, maxstp,  &
            mstpf, scalex, xc)
!
!       WRITE TITLE AND RECORD PARAMETERS.
!
CALL title (cauchy, deuflh, geoms, linesr, newton, overch, acptcr, contyp,  &
            itsclf, itsclx, jactyp, jupdm, maxit, maxns, maxqns, mgll,   &
            minqns, n, narmij, niejev, njacch, nunit, output, qnupdm, stopcr, &
            trupdm, alpha, confac, delfac, delta, epsmch, etafac, fcnold,  &
            ftol, lam0, maxstp, mstpf, nsttol, omega, ratiof, sigma, stptol,  &
            boundl, boundu, fvecc, scalef, scalex, xc)
!
!    INITIALIZE FTRACK AND STRACK VECTORS (FTRACK STORES (TRACKS) THE FUNCTION
!    VALUES FOR THE NONMONOTONIC COMPARISON AND STRACK, SIMILARLY, STORES THE
!    LENGTH OF THE NEWTON STEPS - WHICH ARE USED IN CONJUNCTION WITH
!    DEUFLHARD'S SECOND ACCEPTANCE CRITERION).
!
DO  j=0,mgll-1
  ftrack(j)=zero
  strack(j)=zero
END DO
!
!       MAIN ITERATIVE LOOP - MAXIT IS SPECIFIED BY USER.
!
!       ITNUM COUNTS OVERALL ITERATIONS.
!       ISEJAC COUNTS ITERATIONS SINCE LAST EXPLICIT JACOBIAN
!       EVALUATION IF A QUASI-NEWTON METHOD IS BEING USED.
!
DO  itnum=1,maxit
!
!          SUPPRESS OUTPUT IF DESIRED - USED IF DETAILED OUTPUT IS DESIRED
!          FOR LATER ITERATIONS ONLY (DEFAULT VALUE FOR SUPPRS IS 0 -
!          I.E. NO SUPPRESSION).
!
  IF (itnum < supprs) THEN
    output=3
  ELSE
    output=outtmp
  END IF
!
  IF (output > 2) THEN
    WRITE (nunit,20)
    20 FORMAT (t3,72('*'))
    WRITE (nunit,30)
    30 FORMAT ('   *', t74, '*')
    WRITE (nunit,40) itnum
    40 FORMAT ('   *  ITERATION NUMBER: ',i5, t74, '*')
  END IF
!
!          UPDATE ITERATION COUNTER, ISEJAC, FOR QUASI-NEWTON METHODS
!          ONLY (I.E. JUPDM > 0).
!
  IF (jupdm > 0 .AND. restrt .AND. (.NOT.newton)) THEN
    isejac=1
  ELSE
    isejac=isejac + 1
  END IF
  IF (output > 4 .AND. jupdm > 0 .AND. (.NOT.newton)) THEN
    WRITE (nunit,30)
    IF (restrt) THEN
      IF (itnum > niejev) THEN
        WRITE (nunit,50)
        50 FORMAT ('   *    RESTRT IS TRUE, ISEJAC SET TO 1', t74, '*')
      END IF
    ELSE
      WRITE (nunit,60) isejac
      60 FORMAT ('   *    # OF ITERATIONS SINCE EXPLICIT JACOBIAN, ',  &
                 'ISEJAC, INCREASED TO', i4, t74, '*')
    END IF
  END IF
!
!       WHEN AN EXPLICIT JACOBIAN IS BEING USED IN QUASI-NEWTON METHODS THEN
!       MAXNS STEPS ARE ALLOWED IN THE LINE SEARCH.   FOR STEPS BASED ON A
!       QUASI-NEWTON APPROXIMATION MAXQNS ARE ALLOWED.  SIMILARLY MAXNS AND
!       MAXQNS STEPS ARE ALLOWED IN TRUST REGION METHODS RESPECTIVELY.
!       THIS IS AN ATTEMPT TO AVOID AN EXCESSIVE NUMBER OF FUNCTION
!       EVALUATIONS IN A DIRECTION WHICH WOULD NOT LEAD TO A SIGNIFICANT
!       REDUCTION.
!
  IF (.NOT.newton) THEN
    IF (linesr) THEN
      IF (jupdm > 0) THEN
        IF (isejac == 1) THEN
!
!                      JACOBIAN UPDATED EXPLICITLY.
!
          maxlin=maxns
        ELSE
!
!                      QUASI-NEWTON UPDATE.
!
          maxlin=maxqns
        END IF
      ELSE
        maxlin=maxns
      END IF
    ELSE
      IF (jupdm > 0) THEN
        IF (isejac == 1) THEN
!
!                      JACOBIAN UPDATED EXPLICITLY.
!
          maxtrs=maxns
        ELSE
!
!                      QUASI-NEWTON UPDATE.
!
          maxtrs=maxqns
        END IF
      ELSE
        maxtrs=maxns
      END IF
    END IF
    IF (jupdm > 0 .AND. output > 4) THEN
      WRITE (nunit,30)
      IF (linesr) THEN
        WRITE (nunit,70) maxlin
        70 FORMAT ('   *    MAXLIN SET TO:',i5, t74, '*')
      ELSE
        WRITE (nunit,80) maxtrs
        80 FORMAT ('   *    MAXTRS SET TO:',i5, t74, '*')
      END IF
    END IF
  END IF
!
!       ESTABLISH WHETHER JACOBIAN IS TO BE CHECKED NUMERICALLY -
!       NJACCH ESTABLISHES THE NUMBER OF ITERATIONS FOR WHICH JACOBIAN
!       CHECKING IS DESIRED.  IF CHECKJ IS TRUE A FORWARD DIFFERENCE NUMERICAL
!       APPROXIMATION OF THE JACOBIAN IS COMPARED TO THE ANALYTICAL VERSION.
!       STORE THE NUMBER OF FUNCTION EVALUATIONS SO THAT THESE "EXTRA" ARE NOT
!       INCLUDED IN OVERALL STATISTICS.
!
  IF (jactyp == 0) THEN
    IF (itnum > njacch) THEN
      checkj=.false.
    ELSE
      checkj=.true.
      nfestr=nfetot
    END IF
  ELSE
    checkj=.false.     ! Added by AJM
  END IF
!
!       EVALUATE JACOBIAN AT FIRST STEP OR IF NO QUASI-NEWTON
!       UPDATE IS BEING USED (RESTRT IS FALSE ONLY IN QUASI-
!       NEWTON METHODS WHEN THE QUASI-NEWTON UPDATE IS BEING USED).
!
  IF (restrt .OR. jupdm == 0) THEN
!
!       IF MORE THAN ONE DAMPED NEWTON STEP HAS BEEN REQUESTED AT THE START,
!       IDENTIFY THIS AS THE REASON FOR EXPLICIT JACOBIAN EVALUATION.
!
    IF (jupdm > 0 .AND. (itnum <= niejev .AND. itnum > 1)  .AND. output > 4) THEN
      WRITE (nunit,30)
      WRITE (nunit,90)
      90 FORMAT ('   *    AS ITNUM <= NIEJEV JACOBIAN EVALUATED EXPLICITLY', &
                 t74, '*')
    END IF
!
!             OTHERWISE JUPDM IS 0 OR RESTRT IS TRUE AND ITNUM IS > NIEJEV.
!
    IF (jupdm > 0 .AND. itnum > 1 .AND. output > 4 .AND. itnum > niejev) THEN
      WRITE (nunit,30)
      WRITE (nunit,100)
      100 FORMAT ('   *    RESTRT IS TRUE - JACOBIAN EVALUATED EXPLICITLY',  &
                  t74, '*')
    END IF
!
!          NOTE: MATRIX H IS USED HERE TO HOLD THE FINITE DIFFERENCE
!                ESTIMATION USED IN CHECKING THE ANALYTICAL JACOBIAN IF CHECKJ
!                IS TRUE.   VECTORS WV1 AND, FOR CENTRAL DIFFERENCES, WV2
!                TEMPORARILY HOLD THE FINITE DIFFERENCE FUNCTION EVALUATIONS.
!
    CALL jacobi (checkj, jacerr, overfl, jactyp, n, nunit, output, epsmch,  &
                 fdtolj, boundl, boundu, fvecc, wv1, wv2, jac, h, scalex,  &
                 xc, fvecev, jacev)
    njetot=njetot + 1
!
!          RETURN IF ANALYTICAL AND NUMERICAL JACOBIANS DON'T AGREE
!          (APPLICABLE ONLY IF CHECKJ IS TRUE AND A DISCREPANCY IS FOUND
!          WITHIN SUBROUTINE JACOBI).  A WARNING IS GIVEN FROM WITHIN JACOBI.
!
    IF (jacerr) RETURN
!
!          RESET TOTAL NUMBER OF FUNCTION EVALUATIONS TO NEGLECT
!          THOSE USED IN CHECKING ANALYTICAL JACOBIAN.
!
    IF (checkj) nfetot=nfestr
!
    IF (jupdm > 0) THEN
!
!          FCNMIN IS THE MINIMUM OF THE OBJECTIVE FUNCTION FOUND SINCE
!          THE LAST EXPLICIT JACOBIAN EVALUATION.  IT IS USED TO ESTABLISH
!          WHICH STEP THE PROGRAM RETURNS TO WHEN A QUASI-NEWTON STEP FAILS.
!
      fcnmin=fcnold
!
!          POWTAU IS THE TAU FROM POWELL'S TRUST REGION UPDATING SCHEME
!          (USED WHEN JUPDM=1).  IT IS RESET TO 1.0 AT EVERY EXPLICIT
!          JACOBIAN EVALUATION IN QUASI-NEWTON METHODS.
!
      powtau=one
!
!          AT EVERY EXPLICIT JACOBIAN EVALUATION EXCEPT POSSIBLY THE FIRST,
!          IN QUASI-NEWTON METHODS,  A NEW TRUST REGION IS CALCULATED
!          INTERNALLY USING EITHER THE NEWTON STEP OR THE CAUCHY STEP AS
!          SPECIFIED BY THE USER IN LOGICAL VARIABLE "CAUCHY" IN
!          SUBROUTINE DELCAU.  THIS IS FORCED BY SETTING DELTA TO A
!          NEGATIVE NUMBER.
!
      IF (itnum > 1) delta=-one
!
!             RESET COUNTER FOR NUMBER OF FAILURES OF RATIO TEST (IS
!             FCNNEW/FCNOLD < RATIOF?) SINCE LAST EXPLICIT JACOBIAN UPDATE.
!
      nfail=0
      IF (output > 4 .AND. itnum > niejev) THEN
        WRITE (nunit,30)
        WRITE (nunit,110)
        110 FORMAT ('   *    NUMBER OF FAILURES OF RATIO TEST, NFAIL, ',  &
                    'SET BACK TO 0', t74, '*')
      END IF
!
!             ESTABLISH A NEW MAXIMUM STEP LENGTH ALLOWED.
!
      IF (itnum > 1) CALL maxst (overfl, maxexp, n, nunit, output, epsmch,  &
                                 maxstp, mstpf, scalex, xc)
!
!             SET "P" MATRIX TO IDENTITY FOR LEE AND LEE UPDATE (JUPDM=2).
!
      IF (jupdm == 2 .AND. itnum >= niejev) THEN
        DO  j=1,n
          plee(1:n,j)=zero
          plee(j,j)=one
        END DO
      END IF
    END IF
  END IF
  IF ((restrt .OR. qnupdm == 0) .AND. output > 4 .AND. (.NOT.matsup)) THEN
!
!             WRITE JACOBIAN MATRIX.
!
    WRITE (nunit,30)
    WRITE (nunit,140)
    140 FORMAT ('   *    JACOBIAN MATRIX', t74, '*')
    CALL matprt (n, n, n, n, nunit, jac)
  END IF
!
!       ESTABLISH SCALING MATRICES IF DESIRED (ITSCLF=0 => NO ADAPTIVE
!       SCALING, WHILE ITSCLF > 0 MEANS ADAPTIVE SCALING STARTS AFTER THE
!       (ITSCLF)TH ITERATION).  SIMILARLY FOR COMPONENT SCALING ... .
!       NOTE: SCALING FACTORS ARE UPDATED ONLY WHEN THE JACOBIAN IS UPDATED
!       EXPLICITLY IN QUASI-NEWTON METHODS.
!
!       FUNCTION SCALING.
!
  IF (restrt .AND. itsclf > 0 .AND. itsclf <= itnum) THEN
!
    CALL ascalf (n, epsmch, fvecc, jac, scalef)
!
    IF (output > 4) THEN
      WRITE (nunit,30)
      WRITE (nunit,150)
      150 FORMAT ('   *    FUNCTION SCALING MATRIX', t74, '*')
      WRITE (nunit,30)
      DO  i=1,n
        WRITE (nunit,160) i,scalef(i)
        160 FORMAT ('   *       SCALEF(', I3, ') = ', g12.3, t74, '*')
      END DO
    END IF
!
!          RECALCULATE OBJECTIVE FUNCTION WITH NEW SCALING FACTORS.
!          THIS AVOIDS PREMATURE FAILURES IF THE CHANGE IS SCALING FACTORS
!          WOULD MAKE THE PREVIOUS OBJECTIVE FUNCTION VALUE SMALLER.
!
    CALL fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnold, fvecc,  &
                 scalef)
  END IF
!
!          COMPONENT SCALING.
!
  IF (restrt .AND. itsclx > 0 .AND. itsclx <= itnum) THEN
!
    CALL ascalx (n, epsmch, jac, scalex)
!
    IF (output > 4) THEN
      WRITE (nunit,30)
      WRITE (nunit,180)
      180 FORMAT ('   *       COMPONENT SCALING MATRIX', t74, '*')
      WRITE (nunit,30)
      DO  i=1,n
        WRITE (nunit,190) i,scalex(i)
        190 FORMAT ('   *          SCALEX(', i3, ') = ', g12.3, t74, '*')
      END DO
    END IF
  END IF
!
!          FIND GRADIENT OF 1/2 FVECC^FVECC (NOT USED IF DECOMPOSED
!          MATRIX IS UPDATED IN WHICH CASE THE GRADIENT IS FOUND
!          WITHIN THAT SUBROUTINE - CALL IS MADE FOR OUTPUT ONLY).
!
!
  CALL gradf (overch, overfl, restrt, sclfch, sclxch, jupdm, maxexp, n,  &
              nunit, output, qnupdm, delf, fvecc, jac, scalef, scalex)
!
!          FIND NEWTON STEP USING QR DECOMPOSITION.
!
!          IF JUPDM = 0 OR QNUPDM = 0 THEN THE UNFACTORED FORM
!          FOR CALCULATING THE NEWTON STEP IS USED.
!
  IF (jupdm == 0 .OR. qnupdm == 0) THEN
!
!             NEWTON STEP - UNFACTORED FORM BEING UPDATED.
!
    CALL nstpun (abort, linesr, overch, overfl, qrsing, sclfch, sclxch,  &
                 itnum, maxexp, n, nunit, output, epsmch, a, delf, fvecc, h, &
                 hhpi, jac, rdiag, scalef, scalex, sn)
  ELSE
!
!             NEWTON STEP - FACTORED FORM BEING UPDATED.
!
    CALL nstpfa (abort, linesr, overch, overfl, qrsing, restrt, sclfch,   &
                 sclxch, itnum, maxexp, n, newstm, nunit, output, epsmch, a,  &
                 delf, fvecc, h, hhpi, jac, rdiag, scalef, scalex, sn)
!
  END IF
!
!          RUN IS ABORTED IF THE JACOBIAN BECOMES ESSENTIALLY
!          ALL ZEROS.  NOTE: 203 FOLLOWS END OF MAIN 200 LOOP.
!
  IF (abort) EXIT
!
!          CHECK FOR CONVERGENCE ON LENGTH OF NEWTON STEP IF
!          STOPCR = 1, 12 OR 3.
!
  IF (stopcr /= 2) THEN
    stpmax=zero
    DO  i=1,n
      stpmax=MAX(stpmax, ABS(sn(i))/MAX(xc(i), scalex(i)))
      wv1(i)=scalex(i)*xc(i)
    END DO
    CALL twonrm (overfl, maxexp, n, epsmch, xnorm, wv1)
    IF (stpmax <= nsttol*(one+xnorm)) THEN
      trmcod=1
!
!                IF STOPCR=3 THEN OBJECTIVE FUNCTION VALUE MUST BE
!                DETERMINED AS WELL - OTHERWISE A SOLUTION HAS BEEN FOUND.
!
      IF (stopcr /= 3) THEN
        GO TO 330
!
!                   NOTE: STATEMENT 202 PRECEDES CONVERGENCE
!                         CHECKING SUBROUTINE.
!
      END IF
    END IF
  END IF
!
!          FIND LENGTH OF (SCALED) NEWTON STEP, NEWLEN.
!
  wv1(1:n)=scalex(1:n)*sn(1:n)
  CALL twonrm (overfl, maxexp, n, epsmch, newlen, wv1)
!
!          FOR ITERATIONS AFTER THE ARMIJO STEPS HAVE BEEN COMPLETED
!          AT THE BEGINNING (IN OTHER WORDS THE MONOTONIC STEPS)
!          STORE THE FUNCTION AND, POSSIBLY, THE NEWTON STEP LENGTHS
!          IN THE FTRACK AND STRACK VECTORS, RESPECTIVELY.
!
  IF (isejac >= narmij) THEN
    IF (isejac == 1) THEN
      strack(0)=newlen
!
!             NEWMAX IS USED TO KEEP A BOUND ON THE ENTRIES IN
!             THE STRACK VECTOR.
!
      newmax=newlen
    ELSE
      strack(countr)=MIN(newmax,newlen)
    END IF
!
!          THE OBJECTIVE FUNCTION VALUE IS STORED EVEN IF IT IS
!          GREATER THAN ANY PRECEEDING FUNCTION VALUE.
!
    ftrack(countr)=fcnold
!
!          WRITE FTRACK AND STRACK VECTORS IF DESIRED.  SINCE ONLY THE LAST
!          MGLL VALUES ARE NEEDED THE COUNTER CIRCULATES THROUGH THE VECTOR
!          CAUSING ONLY THE MGLL MOST RECENT VALUES TO BE KEPT.
!          NOTE: THESE VECTORS ARE NOT APPLICABLE IF NEWTON'S METHOD IS BEING
!                USED.
!
    IF (.NOT.newton .AND. output > 4) THEN
      WRITE (nunit,30)
      WRITE (nunit,30)
!
!             IF ONLY THE FUNCTION VALUE ACCEPTANCE TEST IS BEING USED,
!             THUS ACPTCR=1, THEN ONLY THE FTRACK VECTOR IS APPLICABLE.
!
      IF (acptcr == 1) THEN
        WRITE (nunit,230) countr
        230 FORMAT ('   *    CURRENT FTRACK VECTOR; LATEST CHANGE: ELEMENT', &
                    i4, t74, '*')
        WRITE (nunit,30)
        DO  j=0,mgll-1
          WRITE (nunit,240) j,ftrack(j)
          240 FORMAT ('   *       FTRACK(', i3, ') = ', g11.3, t74, '*')
        END DO
!
      ELSE
!
!                BOTH THE FUNCTION VALUE AND THE STEP SIZE
!                ACCEPTANCE TESTS ARE BEING USED, ACPTCR=12.
!
        WRITE (nunit,260) countr
        260 FORMAT ('   *    CURRENT FTRACK AND STRACK VECTORS;',  &
                    '  LATEST CHANGE:  ELEMENT', i4, t74, '*')
        WRITE (nunit,30)
        DO  j=0,mgll-1
          WRITE (nunit,270) j,ftrack(j),j,strack(j)
          270 FORMAT ('   *       FTRACK(', i3, ') = ', g11.3, '  STRACK(',  &
                      i3, ') = ', g11.3, t74, '*')
        END DO
      END IF
    END IF
!
!          UPDATE COUNTING INTEGER, COUNTR. RECYCLE IF COUNTR
!          HAS REACHED MGLL-1.
!
    IF (countr == mgll-1) THEN
      countr=0
    ELSE
      countr=countr + 1
    END IF
  END IF
!
!          RESET STEP ACCEPTANCE CODE AND DELSTR.
!
  acpcod=0
  IF (.NOT.linesr) delstr=zero
!
!          RESET QNFAIL TO FALSE TO AVOID PREMATURE STOPPING.
!
  qnfail=.false.
!
  IF (linesr) THEN
!
!             THE MAIN LINE SEARCH IS CALLED IN SUBROUTINE LINE.
!
    CALL line (abort, absnew, deuflh, geoms, newton, overch, overfl, qnfail, &
               qrsing, restrt, sclfch, sclxch, acpcod, acptcr, contyp,   &
               isejac, itnum, jupdm, maxexp, maxlin, mgll, mnew, n, narmij,  &
               nfunc, nunit, output, qnupdm, stopcr, trmcod, alpha, confac,  &
               epsmch, fcnmax, fcnnew, fcnold, lam0, maxstp, newlen, sbrnrm, &
               sigma, a, boundl, boundu, delf, ftrack, fvec, h, hhpi, jac,   &
               rdiag, wv1, s, sbar, scalef, scalex, sn, strack, xc, xplus,   &
               fvecev)
    IF (abort) THEN
      IF (output > 0) THEN
        WRITE (nunit,30)
        WRITE (nunit,20)
      END IF
      RETURN
    END IF
!
!             NOTE: 201 PRECEDES PRINTING OF ITERATION RESULTS.
!
    IF (newton) GO TO 320
!
  ELSE
!
!             TRUST REGION METHOD
!
!             ESTABLISH INITIAL TRUST REGION SIZE, DELTA, AND/OR
!             FIND LENGTH OF SCALED DESCENT STEP, CAULEN.
!
    CALL delcau (cauchy, overch, overfl, isejac, maxexp, n, nunit, output,  &
                 beta, caulen, delta, epsmch, maxstp, newlen, sqrtz, a, delf, &
                 scalex)
!
    frstdg=.true.
!
    IF (output > 3) THEN
      WRITE (nunit,30)
      WRITE (nunit,30)
      WRITE (nunit,290)
      290 FORMAT ('   *    SUMMARY OF TRUST REGION METHOD USING (DOUBLE) ',  &
                  'DOGLEG STEP', t74, '*')
    END IF
!
!          MAIN INTERNAL LOOP FOR TRUST REGION METHOD.
!
!          THE TRUST REGION SIZE IS STORED FOR COMPARISON
!          LATER TO SET THE PARAMETER POWTAU USED IN POWELL'S
!          TRUST REGION UPDATING SCHEME (QUASI-NEWTON, TRUPDM=1).
!
    delta0=delta
!
    DO  notrst=1,maxtrs
!
      CALL dogleg (frstdg, newtkn, overch, overfl, maxexp, n, notrst, nunit, &
                   output, beta, caulen, delta, etafac, newlen, sqrtz, delf, &
                   s, scalex, sn, ssdhat, vhat)
!
!                NOTE: WV1 AND WV4 HOLD THE COMPONENT AND RESIDUAL VECTOR
!                      RESPECTIVELY FOR A TRIAL POINT WHICH HAS BEEN FOUND
!                      TO BE ACCEPTABLE WHILE THE TRUST REGION IS EXPANDED
!                      AND A NEW TRIAL POINT TESTED.
!                      WV2 AND WV3 ARE WORK VECTORS.
!                      H IS CALLED ASTORE INSIDE TRSTUP.
!
      CALL trstup (geoms, newtkn, overch, overfl, qrsing, sclfch, sclxch,  &
                   acpcod, acpstr, acptcr, contyp, isejac, jupdm, maxexp,  &
                   mgll, mnew, n, narmij, nfunc, notrst, nunit, output,   &
                   qnupdm, retcod, trupdm, alpha, confac, delfac, delstr,  &
                   delta, epsmch, fcnmax, fcnnew, fcnold, fcnpre, maxstp,  &
                   newlen, newmax, powtau, rellen, stptol, a, h, boundl,   &
                   boundu, delf, wv1, ftrack, fvec, fvecc, hhpi, jac, rdiag, &
                   wv2, s, sbar, scalef, scalex, strack, xc, wv4, xplus, fvecev)
!
      IF (output > 4 .OR. (retcod == 7 .AND. output > 2)) CALL  &
          rcdprt (nunit, retcod, delta, rellen, stptol)
!
!                IF NO PROGRESS WAS BEING MADE (RETCOD=7) IN A QUASI-NEWTON
!                STEP RETRY WITH AN EXPLICIT JACOBIAN EVALUATION.
!
      IF (retcod == 7 .AND. (.NOT.restrt)) qnfail=.true.
!
!                RETURN CODE LESS THAN 8 EXITS FROM TRUST REGION LOOP.
!
      IF (retcod < 8) GO TO 310
!
    END DO
!
!          IF NO SUCCESSFUL STEP FOUND IN A QUASI-NEWTON STEP
!          RETRY WITN AN EXPLICIT JACOBIAN EVALUATION.
!
    IF (.NOT.restrt) qnfail=.true.
!
  END IF
!
  310 IF (.NOT.linesr) THEN
    IF (delta < delta0) powtau=one
    delta=MAX(delta,1.0D-10)
  END IF
!
!          IF RETCOD=7 AND STOPCR=2 RESET STOPPING CRITERION TO
!          AVOID HANGING IN TRUST REGION METHOD. (RETCOD=7 MEANS
!          THE RELATIVE STEP LENGTH WAS LESS THAN STPTOL).
!
  IF (.NOT.linesr .AND. retcod == 7 .AND. stopcr == 2 .AND. (.NOT.qnfail )) THEN
    stopcr=12
  END IF
!
!          RETAIN NUMBER OF STEPS ACCEPTED BY EACH CRITERION
!          FOR PERFORMANCE EVALUATION.
!
  IF (.NOT.newton) THEN
    IF (acpcod == 1) THEN
      nac1=nac1+1
    ELSE IF (acpcod == 2) THEN
      nac2=nac2+1
    ELSE IF (acpcod == 12) THEN
      nac12=nac12+1
    END IF
  END IF
!
!          PRINT RESULTS OF ITERATION.
!
  320 IF (output > 2) THEN
!
    CALL nersl (newton, restrt, sclfch, sclxch, acpcod, jupdm, n,  &
                nunit, output, fcnnew, fvec, xplus)
!
  END IF
!
!       CHECK FOR CONVERGENCE.  STATEMENT 202 IS USED IF THE
!       STEP SIZE OF THE NEWTON STEP IS FOUND TO BE WITHIN
!       THE SPECIFIED TOLERANCE AND STOPCR IS 1 OR 12.
!
  330 IF (.NOT.qnfail) THEN
!
!          IF QNFAIL IS TRUE THE QUASI-NEWTON SEARCH FAILED TO
!          FIND A SATISFACTORY STEP - SINCE THE JACOBIAN IS TO
!          BE RE-EVALUATED AVOID PREMATURE STOPPAGES IN NESTOP.
!
    CALL nestop (absnew, linesr, newton, sclfch, sclxch, acptcr, itnum, n,  &
                 nac1, nac2, nac12, nfunc, njetot, nunit, output, stopcr,   &
                 trmcod, fcnnew, ftol, nsttol, stpmax, stptol, fvec, scalef, &
                 scalex, xc, xplus)
  END IF
!
!          IF THE TERMINATION CODE, TRMCOD, IS GREATER THAN 0 THEN
!          CONVERGENCE HAS BEEN REACHED.
!
  IF (trmcod > 0) RETURN
!
!          QUASI-NEWTON UPDATING - JUPDM > 0.
!
  IF (jupdm > 0) THEN
!
!             QNFAIL MEANS A FAILURE IN THE QUASI-NEWTON SEARCH.
!             RE-EVALUATE JACOBIAN AND TRY A DAMPED NEWTON STEP.
!             MAXLIN IS CHANGED FROM MAXQNS TO MAXNS OR MAXTRS IS
!             CHANGED, SIMILARLY, AT THE START OF THE LOOP.
!
    IF (qnfail) THEN
      IF (output > 4) THEN
        WRITE (nunit,30)
        WRITE (nunit,340)
        340 FORMAT ('   *       FAILURE IN QUASI-NEWTON',   &
                    ' SEARCH: QNFAIL IS TRUE', t74, '*')
      END IF
!
    ELSE IF (.NOT.newton) THEN
!
!                QNFAIL IS FALSE - THE STEP HAS BEEN ACCEPTED.
!
      IF (output > 4) THEN
        WRITE (nunit,30)
        WRITE (nunit,350) fcnnew,fcnold,fcnnew/fcnold
        350 FORMAT ('   *       FCNNEW= ', g11.3, '  FCNOLD= ',  &
                    g11.3, '  RATIO= ', g11.3, t74, '*')
      END IF
      IF (fcnnew/fcnold > ratiof) THEN
!
!                   STEP ACCEPTED BUT NOT A SIGNIFICANT IMPROVEMENT.
!
        nfail=nfail + 1
        IF (output > 4) THEN
          WRITE (nunit,30)
          WRITE (nunit,360) nfail
          360 FORMAT ('   *       RATIO > RATIOF SO NFAIL INCREASED TO: ',  &
                      i5, t74, '*')
          WRITE (nunit,30)
        END IF
      ELSE
!
!                   STEP ACCEPTED WITH A SIGNIFICANT IMPROVEMENT.
!
        IF (fcnnew/fcnold > pt01) THEN
          nosubt=0
        ELSE
          nosubt=1
        END IF
!
!                   ITEMPT IS USED LOCALLY FOR OUTPUT CONTROL.
!
        itemp=nfail
        nfail=MAX(nfail-nosubt,0)
        IF (output > 4) THEN
          WRITE (nunit,30)
          IF (itemp == nfail) THEN
            WRITE (nunit,370) nfail
            370 FORMAT ('   *       NFAIL STAYS AT: ',i5, t74, '*')
          ELSE
            WRITE (nunit,380) nfail
            380 FORMAT ('   *       NFAIL CHANGED TO: ',i5, t74, '*')
          END IF
        END IF
      END IF
!
!                SAVE THE RESULTS FOR RESTART IF A FAILURE IN THE
!                QUASI-NEWTON METHOD OCCURS - ESSENTIALLY THIS
!                FINDS THE BEST POINT SO FAR.
!
      IF (isejac == 1 .OR. nfail == 1 .OR. (nfail <= minqns .AND. fcnnew  &
            /fcnmin < one)) THEN
        savest=.true.
        fcnmin=fcnnew
      END IF
      IF (savest) THEN
        savest=.false.
        fstore=fcnnew
        IF (output > 4) THEN
          WRITE (nunit,30)
          WRITE (nunit,390)
          390 FORMAT ('   *       STEP IS SAVED', t74, '*')
          WRITE (nunit,30)
          WRITE (nunit,400)
          400 FORMAT ('   *          SAVED COMPONENT AND',   &
                      ' FUNCTION VALUES', t74, '*')
          WRITE (nunit,30)
        END IF
        DO  i=1,n
          xsave(i)=xplus(i)
          fsave(i)=fvec(i)
          IF (output > 4) THEN
            WRITE (nunit,410) i,xsave(i),i,fsave(i)
            410 FORMAT ('   *       XSAVE(', i3, ') = ', g12.3,   &
                        '      FSAVE (', i3, ') = ', g12.3, t74, '*')
          END IF
        END DO
        IF (output > 4) WRITE (nunit,30)
        itstr=itnum
      END IF
    END IF
!
!          NOTE: IF QNFAIL IS TRUE THEN NFAIL CANNOT HAVE INCREASED IMPLYING
!                THAT NFAIL CANNOT NOW BE GREATER THAN MINQNS.
!
    IF (qnfail .OR. nfail > minqns) THEN
!
!                RESTART FROM BEST POINT FOUND SO FAR.
!
      restrt=.true.
      IF (output > 4) THEN
        WRITE (nunit,30)
        WRITE (nunit,30)
        WRITE (nunit,430)
        430 FORMAT ('   *       RESTRT IS TRUE', t74, '*')
      END IF
      DO  j=0,mgll-1
        ftrack(j)=zero
      END DO
      IF (acptcr == 12) THEN
        DO  j=0,mgll-1
          strack(j)=zero
        END DO
      END IF
      countr=0
      isejac=0
      mnew=0
      trmcod=0
      fcnold=fstore
      IF (output > 4) THEN
        WRITE (nunit,30)
        WRITE (nunit,460) itstr
        460 FORMAT ('   *       RETURN TO ITERATION:', i5, '   WHERE',  &
                    t74, '*')
        WRITE (nunit,30)
      END IF
      DO  i=1,n
        xc(i)=xsave(i)
        fvecc(i)=fsave(i)
        IF (output > 4) THEN
          WRITE (nunit,470) i,xc(i),i,fvecc(i)
          470 FORMAT ('   *       XC(', i3, ') = ', g12.3, '   FVECC(',  &
                      i3, ') = ', g12.3, t74, '*')
        END IF
      END DO
      IF (output > 4) WRITE (nunit,30)
    ELSE
      IF (itnum >= niejev) restrt=.false.
    END IF
  END IF
!
!          UPDATE JACOBIAN IF DESIRED:
!                QNUPDM = 0  => ACTUAL JACOBIAN BEING UPDATED
!                QNUPDM = 1  => FACTORED JACOBIAN BEING UPDATED
!
  IF (.NOT.restrt) THEN
!
    IF (qnupdm == 0) THEN
!
      IF (jupdm == 1) THEN
!
!                   USE BROYDEN UPDATE.
!
        CALL broyun (overfl, maxexp, n, nunit, output, epsmch,  &
                     fvec, fvecc, jac, scalex, xc, xplus)

      ELSE IF (jupdm == 2) THEN
!
!                   USE LEE AND LEE UPDATE.
!
        CALL llun (overch, overfl, isejac, maxexp, n, nunit, output, epsmch, &
                   omega, fvec, fvecc, jac, plee, s, scalex, xc, xplus)
!
      END IF
!
    ELSE
!
!                THE FACTORED FORM OF THE JACOBIAN IS UPDATED.
!
      IF (jupdm == 1) THEN
!
        CALL broyfa (overch, overfl, sclfch, sclxch, maxexp, n, nunit,   &
                     output, epsmch, a, delf, fvec, fvecc, jac, rdiag,  &
                     s, scalef, scalex, wv1, wv2, xc, xplus)
!
      ELSE IF (jupdm == 2) THEN
!
        CALL llfa (overch, overfl, sclfch, sclxch, isejac, maxexp, n, nunit,  &
                   output, epsmch, omega, a, delf, fvec, fvecc, jac, plee,  &
                   rdiag, s, scalef, scalex, xc, xplus)
!
      END IF
!
    END IF
!
  END IF
!
!          UPDATE CURRENT VALUES - RESET TRMCOD TO ZERO.
!
!          UPDATE M "VECTOR" (ACTUALLY ONLY THE LATEST VALUE IS NEEDED).
!
  mold=mnew
  IF (isejac < narmij) THEN
    mnew=0
  ELSE
    mnew=MIN(mold+1,mgll-1)
  END IF
  IF (jupdm == 0 .OR. (jupdm > 0 .AND. (.NOT.restrt))) THEN
!
    CALL update (mnew, mold, n, trmcod, fcnnew, fcnold, fvec,  &
                 fvecc, xc, xplus)
  END IF
!
END DO
!
IF (output > 0) THEN
  WRITE (nunit,20)
  WRITE (nunit,30)
  WRITE (nunit,510) itnum-1
  510 FORMAT ('   *  NO SOLUTION FOUND AFTER', i6, ' ITERATION(S)', t74, '*')
  WRITE (nunit,30)
  WRITE (nunit,520)
  520 FORMAT ('   *', t12, 'FINAL ESTIMATES', t42, 'FINAL FUNCTION VALUES',  &
              t74, '*')
  WRITE (nunit,30)
  DO  i=1,n
    WRITE (nunit,530) i,xplus(i),i,fvec(i)
    530 FORMAT ('   *', t9, 'X(', i3, ') = ', g12.3, t42, 'F(', i3, ') = ',  &
                g12.3, t74, '*')
  END DO
  WRITE (nunit,30)
  WRITE (nunit,550) fcnnew
  550 FORMAT ('   *  FINAL OBJECTIVE FUNCTION VALUE = ', g12.3,  &
      t74, '*')
  WRITE (nunit,30)
  WRITE (nunit,20)
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE NNES.
!
END SUBROUTINE nnes



SUBROUTINE nstpfa (abort, linesr, overch, overfl, qrsing, restrt, sclfch,  &
                   sclxch, itnum, maxexp, n, newstm, nunit, output, epsmch, &
                   a, delf, fvecc, h, hhpi, jac, rdiag, scalef, scalex, sn)

! N.B. Arguments WV1, WV2 & WV3 have been removed.
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS THE NEWTON STEP.
!
!    IF THE JACOBIAN IS DETECTED AS SINGULAR OR IF THE ESTIMATED CONDITION
!    NUMBER IS TOO HIGH (GREATER THAN EPSMCH**(-2/3)) THEN H:=J^J IS FORMED
!    AND THE DIAGONAL IS PERTURBED BY ADDING SQRT(N*EPSMCH)*H1NORM*SCALEX(I)**2
!    TO THE CORRESPONDING ELEMENT.  A CHOLESKY DECOMPOSITION IS PERFORMED ON
!    THIS MODIFIED MATRIX PRODUCING A PSEUDO-NEWTON STEP.
!
!    IF THE CONDITION NUMBER IS SMALL THEN THE NEWTON STEP, SN,
!    IS FOUND DIRECTLY BY BACK SUBSTITUTION.
!
!    ABORT    IF THE 1-NORM OF MATRIX H BECOMES TOO SMALL
!             ALTHOUGH BUT NOT AT A SOLUTION
!    BYPASS   ALLOWS BYPASSING OF THE SPECIAL TREATMENT FOR
!             BADLY CONDITIONED JACOBIANS
!    QRSING   INDICATES SINGULAR JACOBIAN DETECTED
!

LOGICAL, INTENT(OUT)       :: abort
LOGICAL, INTENT(IN OUT)    :: linesr
LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(OUT)       :: qrsing
LOGICAL, INTENT(IN)        :: restrt
LOGICAL, INTENT(IN OUT)    :: sclfch
LOGICAL, INTENT(IN)        :: sclxch
INTEGER, INTENT(IN)        :: itnum
INTEGER, INTENT(IN OUT)    :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: newstm
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN)      :: delf(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: rdiag(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(OUT)     :: sn(n)

! Local variables

LOGICAL    :: pertrb
REAL (dp)  :: connum, h1norm, maxadd, maxffl, scalfi, scalxj, sqrtep, SUM,  &
              wv1(n)
INTEGER    :: i, j
!
abort=.false.
overfl=.false.
sqrtep=SQRT(epsmch)
!
IF (restrt) THEN
  IF (output > 3) THEN
    WRITE (nunit,10)
    10 FORMAT ('   *', t74, '*')
    WRITE (nunit,10)
    WRITE (nunit,20)
    20 FORMAT ('   *    SOLUTION OF LINEAR SYSTEM FOR NEWTON STEP, SN',  &
               t74, '*')
  END IF
!
!       STORE (POSSIBLY SCALED) JACOBIAN IN MATRIX A.
!
  IF (.NOT.sclfch) THEN
    CALL matcop (n, n, n, n, n, n, jac, a)
  ELSE
    DO  i=1,n
      IF (scalef(i) /= one) THEN
        scalfi=scalef(i)
        a(i,1:n)=jac(i,1:n)*scalfi
      ELSE
        a(i,1:n)=jac(i,1:n)
      END IF
    END DO
  END IF
!
!          SCALED JACOBIAN IS PRINTED ONLY IF AT LEAST ONE SCALING
!          FACTOR IS NOT 1.0.
!
  IF (output > 4 .AND. sclfch) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,60)
    60 FORMAT ('   *       SCALED JACOBIAN MATRIX', t74, '*')
    CALL matprt (n, n, n, n, nunit, a)
  END IF
!
!          QR DECOMPOSITION OF (POSSIBLY SCALED) JACOBIAN.
!
  CALL qrdcom (qrsing, n, epsmch, a, hhpi, rdiag)
!
  IF (output > 4 .AND. n > 1 .AND. (.NOT.matsup)) THEN
    WRITE (nunit,10)
    IF (.NOT.sclfch) THEN
      WRITE (nunit,70)
      70 FORMAT ('   *       QR DECOMPOSITION OF JACOBIAN MATRIX', t74, '*')
    ELSE
      WRITE (nunit,80)
      80 FORMAT ('   *       QR DECOMPOSITION OF SCALED JACOBIAN MATRIX',  &
                 t74, '*')
    END IF
    CALL matprt (n, n, n, n, nunit, a)
    WRITE (nunit,10)
    WRITE (nunit,90)
    90 FORMAT ('   *', t15, 'DIAGONAL OF R', t43, 'PI FACTORS FROM QR ',  &
               'DECOMPOSITION', t74, '*')
    WRITE (nunit,10)
    DO  i=1,n-1
      WRITE (nunit,100) i,rdiag(i),i,hhpi(i)
      100 FORMAT ('   *       RDIAG(', i3, ') = ', g11.3, '        HHPI(',  &
                  i3, ') = ', g12.3, t74, '*')
    END DO
    WRITE (nunit,120) n,rdiag(n)
    120 FORMAT ('   *       RDIAG(', i3, ') = ', g11.3, t74, '*')
    WRITE (nunit,10)
    IF (itnum == 1) THEN
      WRITE (nunit,130)
      130 FORMAT ('   *       NOTE: R IS IN STRICT UPPER',   &
                  ' TRIANGLE OF MATRIX A PLUS RDIAG', t74, '*')
      WRITE (nunit,10)
      WRITE (nunit,140)
      140 FORMAT ('   *', t16, 'THE COLUMNS OF THE LOWER',   &
                  ' TRIANGLE OF MATRIX A PLUS', t74, '*')
      WRITE (nunit,150)
      150 FORMAT ('   *', t16, 'THE ELEMENTS OF VECTOR HHPI',   &
                  ' FORM THE HOUSEHOLDER', t74, '*')
      WRITE (nunit,160)
      160 FORMAT ('   *', t16, 'MATRICES WHICH, WHEN',   &
                  ' MULTIPLIED TOGETHER, FORM Q', t74, '*')
    END IF
  END IF
!
!       FORM THE ACTUAL Q^ MATRIX FROM THE HOUSEHOLDER TRANSFORMATIONS STORED
!       IN THE LOWER TRIANGLE OF A AND THE FACTORS IN HHPI: STORE IT IN JAC.
!
  CALL qform (n, a, hhpi, jac)
!
!          COMPLETE THE UPPER TRIANGULAR R MATRIX BY REPLACING THE
!          DIAGONAL OF A.  THE QR DECOMPOSITION IS NOW AVAILABLE.
!
  DO  i=1,n
    a(i,i)=rdiag(i)
  END DO
!
ELSE
!
!          USING UPDATED FACTORED FORM OF JACOBIAN - CHECK FOR SINGULARITY.
!
  qrsing=.false.
  DO  i=1,n
    IF (a(i,i) == zero) qrsing=.true.
  END DO
!
END IF
!
!       ESTIMATE CONDITION NUMBER IF JACOBIAN IS NOT SINGULAR.
!
IF (.NOT.bypass .AND. (.NOT.qrsing) .AND. n > 1) THEN
  IF (sclxch) THEN
!
!             SET UP FOR CONDITION NUMBER ESTIMATION - SCALE R WRT X'S.
!
    DO  j=1,n
      IF (scalex(j) /= one) THEN
        scalxj=scalex(j)
        rdiag(j)=rdiag(j)/scalxj
        a(1:j-1,j)=a(1:j-1,j)/scalxj
      END IF
    END DO
  END IF
!
  CALL condno (overch, overfl, maxexp, n, nunit, output, connum, a, rdiag)
!
!          UNSCALE R IF IT WAS SCALED BEFORE THE CALL TO CONDNO.
!
  IF (sclxch) THEN
    DO  j=1,n
      IF (scalex(j) /= one) THEN
        scalxj=scalex(j)
        rdiag(j)=rdiag(j)*scalxj
        a(1:j-1,j)=a(1:j-1,j)*scalxj
      END IF
    END DO
  END IF
!
!          IF OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR ASSIGN
!          QRSING AS TRUE SO THAT THE JACOBIAN WILL BE PERTURBED.
!
  IF (overfl) qrsing=.true.
!
!          NOTE: OVERFL SWITCHED TO FALSE BEFORE FORMATION OF H.
!
ELSE
  IF (n == 1) THEN
    connum=one
  ELSE
    connum=zero
  END IF
END IF
!
!       MATRIX H=JAC^JAC IS FORMED IN TWO CASES:
!          1)  THE (SCALED) JACOBIAN IS SINGULAR
!          2)  THE CONDITION NUMBER IS TOO HIGH AND THE OPTION TO BYPASS
!              THE PERTURBATION OF THE JACOBIAN IS NOT BEING USED
!          3)  REQUESTED BY THE USER THROUGH NEWSTM.
!
IF (qrsing .OR. ((.NOT.bypass) .AND. connum > one/sqrtep**1.333)  &
       .OR. newstm == 77) THEN
!
!          FORM H:=(DF*JAC)^(DF*JAC) WHERE DF=DIAG(SCALEF).  USE
!          PREVIOUSLY COMPUTED QR DECOMPOSITION OF (SCALED) JACOBIAN
!          WHERE R IS STORED IN THE UPPER TRIANGLE OF A AND RDIAG.
!
  overfl=.false.
  IF (overch) THEN
    CALL ataov (overfl, maxexp, n, nunit, output, jac, h, scalef)
  ELSE
    CALL rtrmul (n, a, h, rdiag)
  END IF
  IF (output > 3) THEN
    WRITE (nunit,10)
    IF (qrsing .AND. (.NOT.overfl)) THEN
      WRITE (nunit,230)
      230 FORMAT ('   *       SINGULAR JACOBIAN DETECTED:',   &
                  ' JACOBIAN PERTURBED', t74, '*')
    ELSE
      IF (overfl) THEN
        WRITE (nunit,240)
        240 FORMAT ('   *       POTENTIAL OVERFLOW DETECTED',   &
                    ' IN CONDITION NUMBER ESTIMATOR', t74, '*')
        WRITE (nunit,250)
        250 FORMAT ('   *       MATRIX "ASSIGNED" AS ',   &
                    'SINGULAR AND JACOBIAN PERTURBED', t74, '*')
      ELSE
        WRITE (nunit,260) connum
        260 FORMAT ('   *       CONDITION NUMBER TOO HIGH: ', g12.3,  &
                    ' JACOBIAN PERTURBED', t74, '*')
      END IF
    END IF
  END IF
  overfl=.false.
  IF (newstm /= 77) THEN
!
!             FIND 1-NORM OF H MATRIX AND PERTURB DIAGONAL.
!
    CALL onenrm (abort, pertrb, n, nunit, output, epsmch, h1norm, h, scalex)
  END IF
!
!          CHOLESKY DECOMPOSITION OF H MATRIX - MAXFFL=0 INDICATES
!          THAT H IS KNOWN TO BE POSITIVE DEFINITE.
!
  maxffl=zero
  CALL cholde (n, maxadd, maxffl, sqrtep, h, a)
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,270)
    270 FORMAT ('   *     CHOLESKY DECOMPOSITION OF H MATRIX', t74, '*' )
    CALL matprt (n, n, n, n, nunit, a)
  END IF
!
!          FIND NEWTON STEP FROM CHOLESKY DECOMPOSITION.  IF THE
!          DIAGONAL HAS BEEN PERTURBED THEN THIS IS NOT THE ACTUAL
!          NEWTON STEP BUT ONLY AN APPROXIMATION THEREOF.
!
  wv1(1:n) = -delf(1:n)
  CALL chsolv (overch, overfl, maxexp, n, nunit, output, a, wv1, sn)
  overfl=.false.
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      IF (.NOT.pertrb) THEN
        WRITE (nunit,290)
        290 FORMAT ('   *     NEWTON STEP FROM CHOLESKY DECOMPOSITION',  &
                    t74, '*')
      ELSE
        WRITE (nunit,300)
        300 FORMAT ('   *     APPROXIMATE NEWTON STEP FROM',   &
                    ' PERTURBED JACOBIAN', t74, '*')
      END IF
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,310) i,sn(i)
        310 FORMAT ('   *       SN(', i3, ') = ', g12.3, t74, '*')
      END DO
    ELSE
      WRITE (nunit,10)
      IF (.NOT.pertrb) THEN
        WRITE (nunit,330)
        330 FORMAT ('   *     NEWTON STEP FROM CHOLESKY DECOMPOSITION',   &
                    '    IN SCALED UNITS', t74, '*')
      ELSE
        WRITE (nunit,340)
        340 FORMAT ('   *     APPROXIMATE NEWTON STEP FROM',   &
                    ' PERTURBED JACOBIAN   IN SCALED UNITS', t74, '*')
      END IF
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,350) i,sn(i),i,scalex(i)*sn(i)
        350 FORMAT ('   *       SN(', i3, ') = ', g12.3, t48, 'SN(', i3,   &
                    ') = ', g12.3, t74, '*')
      END DO
    END IF
  END IF
!
!          SET QRSING TO TRUE SO THAT THE CORRECT MATRIX FACTORIZATION
!          IS USED IN THE BACK-CALCULATION OF SBAR FOR DEUFLHARD
!          RELAXATION FACTOR INITIALIZATION.
!
  qrsing=.true.
ELSE
  IF (output > 3 .AND. n > 1) THEN
    IF (.NOT.bypass .AND. connum <= one/sqrtep**1.33) THEN
      WRITE (nunit,10)
      WRITE (nunit,370) connum
      370 FORMAT ('   *       CONDITION NUMBER ACCEPTABLE, ', g9.2,  &
                  ', JACOBIAN NOT PERTURBED', t74, '*')
    END IF
    IF (bypass .AND. connum > one/sqrtep**1.33) THEN
      WRITE (nunit,10)
      WRITE (nunit,380) connum
      380 FORMAT ('   *       CONDITION NUMBER HIGH, ', g9.2,   &
                  ', JACOBIAN NOT PERTURBED AS', t74, '*')
      WRITE (nunit,390)
      390 FORMAT ('   *       BYPASS IS TRUE', t74, '*')
    END IF
  END IF
  DO  i=1,n
    sum=zero
    DO  j=1,n
      sum=sum - jac(i,j)*scalef(j)*fvecc(j)
    END DO
    sn(i)=sum
  END DO
  CALL rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sn)
  overfl=.false.
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      WRITE (nunit,10)
      WRITE (nunit,420)
      420 FORMAT ('   *       NEWTON STEP FROM QR DECOMPOSITION ', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,310) i,sn(i)
      END DO
    ELSE
      WRITE (nunit,10)
      WRITE (nunit,440)
      440 FORMAT ('   *       NEWTON STEP FROM QR DECOMPOSITION ',  &
                  '       IN SCALED UNITS', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,350) i,sn(i),i,scalex(i)*sn(i)
      END DO
    END IF
  END IF
!
!       TRANSFORM MATRICES FOR SUBSEQUENT CALCULATIONS IN TRUST REGION METHOD.
!
  IF (.NOT.linesr) THEN
    DO  i=2,n
      DO  j=1,i-1
        a(i,j)=a(j,i)
      END DO
    END DO
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE NSTPFA.
!
END SUBROUTINE nstpfa



SUBROUTINE nstpun (abort, linesr, overch, overfl, qrsing, sclfch, sclxch,  &
                   itnum, maxexp, n, nunit, output, epsmch, a, delf, fvecc, &
                   h, hhpi, jac, rdiag, scalef, scalex, sn)

! N.B. Arguments WV1, WV2 & WV3 have been removed.
!
!   FEB. 23, 1992
!
!   THIS SUBROUTINE FINDS THE NEWTON STEP.
!
!   IF THE JACOBIAN IS DETECTED AS SINGULAR OR IF THE ESTIMATED CONDITION
!   NUMBER IS TOO HIGH (GREATER THAN EPSMCH**(-2/3)) THEN H:=J^J IS FORMED
!   AND THE DIAGONAL IS PERTURBED BY ADDING SQRT(N*EPSMCH)*H1NORM*SCALEX(I)**2
!   TO THE CORRESPONDING ELEMENT.   A CHOLESKY DECOMPOSITION IS PERFORMED ON
!   THIS MODIFIED MATRIX PRODUCING A PSEUDO-NEWTON STEP.
!   NOTE: THIS PROCEDURE MAY BE BE BYPASSED FOR ILL-CONDITIONED JACOBIANS
!   BY SETTING THE LOGICAL VARIABLE BYPASS TO TRUE IN THE DRIVER.
!
!   ABORT    IF THE 1-NORM OF MATRIX H BECOMES TOO SMALL
!            ALTHOUGH BUT NOT AT A SOLUTION
!   BYPASS   ALLOWS BYPASSING OF THE SPECIAL TREATMENT FOR
!            BADLY CONDITIONED JACOBIANS
!   QRSING   INDICATES SINGULAR JACOBIAN DETECTED
!

LOGICAL, INTENT(OUT)       :: abort
LOGICAL, INTENT(IN OUT)    :: linesr
LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(OUT)       :: qrsing
LOGICAL, INTENT(IN OUT)    :: sclfch
LOGICAL, INTENT(IN)        :: sclxch
INTEGER, INTENT(IN)        :: itnum
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(OUT)     :: a(n,n)
REAL (dp), INTENT(IN)      :: delf(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(IN)      :: jac(n,n)
REAL (dp), INTENT(OUT)     :: rdiag(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(OUT)     :: sn(n)

! Local variables

REAL (dp)  :: connum, h1norm, maxadd, maxffl, scalfi, scalxj, sqrtep
REAL (dp)  :: wv1(n)
LOGICAL    :: pertrb
INTEGER    :: i, j
INTEGER    :: newstm = 0   ! Added by AJM
!
abort=.false.
overfl=.false.
pertrb=.false.
sqrtep=SQRT(epsmch)
!
IF (output > 3) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,10)
  WRITE (nunit,20)
  20 FORMAT ('   *    SOLUTION OF LINEAR SYSTEM FOR NEWTON STEP, SN', t74, '*')
END IF
!
!       STORE (POSSIBLY SCALED) JACOBIAN IN MATRIX A.
!
IF (.NOT.sclfch) THEN
  CALL matcop (n, n, n, n, n, n, jac, a)
ELSE
  DO  i=1,n
    IF (scalef(i) /= one) THEN
      scalfi=scalef(i)
      a(i,1:n)=jac(i,1:n)*scalfi
    ELSE
      a(i,1:n)=jac(i,1:n)
    END IF
  END DO
END IF
!
!       SCALED JACOBIAN IS PRINTED ONLY IF AT LEAST ONE SCALING
!       FACTOR IS NOT 1.0.
!
IF (output > 4 .AND. sclfch) THEN
  WRITE (nunit,10)
  WRITE (nunit,10)
  WRITE (nunit,60)
  60 FORMAT ('   *       SCALED JACOBIAN MATRIX', t74, '*')
  CALL matprt (n, n, n, n, nunit, a)
END IF
!
!       QR DECOMPOSITION OF (POSSIBLY SCALED) JACOBIAN.
!
CALL qrdcom (qrsing, n, epsmch, a, hhpi, rdiag)
!
!    SAVE MATRIX A FOR BACK SUBSTITUTION TO CHECK DEUFLHARDS'S SECOND
!    STEP ACCEPTANCE CRITERION IN LINE SEARCH OR TRUST REGION METHOD.
!
CALL matcop (n, n, n, n, n, n, a, h)
!
IF (output > 4 .AND. n > 1 .AND. (.NOT.matsup)) THEN
  WRITE (nunit,10)
  IF (.NOT.sclfch) THEN
    WRITE (nunit,70)
    70 FORMAT ('   *       QR DECOMPOSITION OF JACOBIAN MATRIX', t74, '*')
  ELSE
    WRITE (nunit,80)
    80 FORMAT ('   *       QR DECOMPOSITION OF SCALED JACOBIAN MATRIX',  &
               t74, '*')
  END IF
  CALL matprt (n, n, n, n, nunit, a)
  WRITE (nunit,10)
  WRITE (nunit,90)
  90 FORMAT ('   *', t16, 'DIAGONAL OF R          PI FACTORS',   &
             ' FROM QR DECOMPOSITION', t74, '*')
  WRITE (nunit,10)
  DO  i=1,n-1
    WRITE (nunit,100) i,rdiag(i),i,hhpi(i)
    100 FORMAT ('   *       RDIAG(', i3, ') = ', g11.3, '        HHPI(',  &
                i3, ') = ', g12.3, t74, '*')
  END DO
  WRITE (nunit,120) n,rdiag(n)
  120 FORMAT ('   *       RDIAG(', i3, ') = ', g11.3, t74, '*')
  WRITE (nunit,10)
  IF (itnum == 1) THEN
    WRITE (nunit,130)
    130 FORMAT ('   *       NOTE: R IS IN STRICT UPPER TRIANGLE',   &
                ' OF MATRIX A PLUS RDIAG', t74, '*')
    WRITE (nunit,10)
    WRITE (nunit,140)
    140 FORMAT ('   *', t16, 'THE COLUMNS OF THE LOWER TRIANGLE',   &
                ' OF MATRIX A PLUS', t74, '*'/   &
                '   *', t16, 'THE ELEMENTS OF VECTOR HHPI FORM THE ',  &
                'HOUSEHOLDER', t74, '*')
    WRITE (nunit,150)
    150 FORMAT ('   *', t16, 'MATRICES WHICH, WHEN MULTIPLIED',   &
                ' TOGETHER, FORM Q', t74, '*')
  END IF
END IF
!
!       ESTIMATE CONDITION NUMBER IF (SCALED) JACOBIAN IS NOT SINGULAR.
!
IF (.NOT.bypass .AND. (.NOT.qrsing) .AND. n > 1) THEN
  IF (sclxch) THEN
!
!             SET UP FOR CONDITION NUMBER ESTIMATOR - SCALE R WRT X'S.
!
    DO  j=1,n
      IF (scalex(j) /= one) THEN
        scalxj=scalex(j)
        rdiag(j)=rdiag(j)/scalxj
        a(1:j-1,j)=a(1:j-1,j)/scalxj
      END IF
    END DO
  END IF
!
  CALL condno (overch, overfl, maxexp, n, nunit, output, connum, a, rdiag)
!
!          IF OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR ASSIGN
!          QRSING AS TRUE.
!
  IF (overfl) qrsing=.true.
!
!          NOTE: OVERFL SWITCHED TO FALSE BEFORE FORMATION OF H LATER.
!
ELSE
!
!          ASSIGN DUMMY TO CONNUM FOR SINGULAR JACOBIAN UNLESS N=1.
!
  IF (n == 1) THEN
    connum=one
  ELSE
    connum=zero
  END IF
END IF
!
!    MATRIX H=JAC^JAC IS FORMED IN THREE CASES:
!       1)  THE (SCALED) JACOBIAN IS SINGULAR
!       2)  THE CONDITION NUMBER IS TOO HIGH AND THE OPTION TO BYPASS THE
!           PERTURBATION OF THE JACOBIAN IS NOT BEING USED.
!       3)  REQUESTED BY USER THROUGH NEWSTM.
!
IF (qrsing .OR. ((.NOT.bypass) .AND. connum > one/sqrtep**1.333)  &
       .OR. newstm == 77) THEN
!
!          FORM H:=(DF*JAC)^(DF*JAC) WHERE DF=DIAG(SCALEF).  USE
!          PREVIOUSLY COMPUTED QR DECOMPOSITION OF (SCALED) JACOBIAN
!          WHERE R IS STORED IN THE UPPER TRIANGLE OF A AND RDIAG.
!
  IF (overch) THEN
    overfl=.false.
    CALL ataov (overfl, maxexp, n, nunit, output, jac, h, scalef)
  ELSE
    CALL rtrmul (n, a, h, rdiag)
  END IF
  IF (output > 3) THEN
    WRITE (nunit,10)
    IF (qrsing .AND. (.NOT.overfl)) THEN
      WRITE (nunit,180)
      180 FORMAT ('   *       SINGULAR JACOBIAN DETECTED:',   &
                  ' JACOBIAN PERTURBED', t74, '*')
    ELSE
!
!                NOTE: IF OVERFL IS TRUE THEN QRSING MUST BE TRUE.
!
      IF (overfl) THEN
        WRITE (nunit,190)
        190 FORMAT ('   *       POTENTIAL OVERFLOW DETECTED',   &
                    ' IN CONDITION NUMBER ESTIMATOR', t74, '*')
        WRITE (nunit,200)
        200 FORMAT ('   *       MATRIX "ASSIGNED" AS ',   &
                    'SINGULAR AND JACOBIAN PERTURBED', t74, '*')
      ELSE
        WRITE (nunit,210) connum
        210 FORMAT ('   *       CONDITION NUMBER TOO HIGH: ', g12.3,  &
                    ' JACOBIAN PERTURBED', t74, '*')
      END IF
    END IF
  END IF
  overfl=.false.
  IF (newstm /= 77) THEN
!
!             FIND 1-NORM OF H MATRIX AND PERTURB DIAGONAL.
!
    CALL onenrm (abort, pertrb, n, nunit, output, epsmch, h1norm, h, scalex)
    IF (abort) RETURN
  END IF
!
!          CHOLESKY DECOMPOSITION OF H MATRIX - MAXFFL=0 IMPLIES
!          THAT H IS KNOWN TO BE POSITIVE DEFINITE.
!
  maxffl=zero
  CALL cholde (n, maxadd, maxffl, sqrtep, h, a)
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,220)
    220 FORMAT ('   *    CHOLESKY DECOMPOSITION OF H MATRIX', t74, '*' )
    CALL matprt (n, n, n, n, nunit, a)
  END IF
!
!       FIND NEWTON STEP FROM CHOLESKY DECOMPOSITION.  IF THE DIAGONAL HAS
!       BEEN PERTURBED THEN THIS IS NOT THE ACTUAL NEWTON STEP BUT
!       ONLY AN APPORXIMATION THEREOF.
!
  wv1(1:n)=-delf(1:n)
  CALL chsolv (overch, overfl, maxexp, n, nunit, output, a, wv1, sn)
  overfl=.false.
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      WRITE (nunit,10)
      IF (.NOT.pertrb) THEN
        WRITE (nunit,240)
        240 FORMAT ('   *', t8, 'NEWTON STEP FROM CHOLESKY DECOMPOSITION',  &
                    t74, '*')
      ELSE
        WRITE (nunit,250)
        250 FORMAT ('   *', t8, 'APPROXIMATE NEWTON STEP FROM',   &
                    ' PERTURBED JACOBIAN', t74, '*')
      END IF
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,260) i,sn(i)
        260 FORMAT ('   *       SN(', I3, ') = ', g12.3, t74, '*')
      END DO
    ELSE
      WRITE (nunit,10)
      IF (.NOT.pertrb) THEN
        WRITE (nunit,280)
        280 FORMAT ('   *', t8, 'NEWTON STEP FROM CHOLESKY',   &
                    ' DECOMPOSITION    IN SCALED UNITS', t74, '*')
      ELSE
        WRITE (nunit,290)
        290 FORMAT ('   *', t8, 'APPROXIMATE NEWTON STEP FROM',   &
                    ' PERTURBED JACOBIAN   IN SCALED UNITS', t74, '*')
      END IF
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,300) i,sn(i),i,scalex(i)*sn(i)
        300 FORMAT ('   *       SN(', i3, ') = ', g12.3, t48, 'SN(', i3,  &
                    ')= ', g12.3, t74, '*')
      END DO
    END IF
  END IF
!
!       SET QRSING TO TRUE SO THAT THE CORRECT MATRIX FACTORIZATION IS USED
!       IN THE BACK-CALCULATION OF SBAR FOR DEUFLHARD RELAXATION FACTOR
!       INITIALIZATION IN LINE SEARCH (ONLY MATTERS WHEN JACOBIAN IS ILL-
!       CONDITIONED BUT NOT SINGULAR).
!
  qrsing=.true.
ELSE
  IF (output > 3 .AND. n > 1) THEN
    IF (.NOT.bypass .AND. connum <= one/sqrtep**1.33) THEN
      WRITE (nunit,10)
      WRITE (nunit,320) connum
      320 FORMAT ('   *       CONDITION NUMBER ACCEPTABLE, ', g9.2,  &
                  ', JACOBIAN NOT PERTURBED', t74, '*')
    END IF
    IF (bypass .AND. connum > one/sqrtep**1.33) THEN
      WRITE (nunit,10)
      WRITE (nunit,330) connum
      330 FORMAT ('   *       CONDITION NUMBER HIGH, ', g9.2,   &
                  ', JACOBIAN NOT PERTURBED AS', t74, '*')
      WRITE (nunit,340)
      340 FORMAT ('   *       BYPASS IS TRUE', t74, '*')
    END IF
  END IF
!
!          NOTE: HERE SN STORES THE R.H.S. - IT IS OVERWRITTEN.
!
  sn(1:n) = -fvecc(1:n)*scalef(1:n)
  IF (.NOT.bypass .AND. sclxch) THEN
!
!             R WAS SCALED BEFORE THE CONDITION NUMBER ESTIMATOR -
!             THIS CONVERTS IT BACK TO THE UNSCALED FORM.
!
    DO  j=1,n
      IF (scalex(j) /= one) THEN
        scalxj=scalex(j)
        rdiag(j)=rdiag(j)*scalxj
        a(1:j-1,j)=a(1:j-1,j)*scalxj
      END IF
    END DO
  END IF
!
!          ACCEPTABLE CONDITION NUMBER - USE BACK SUBSTITUTION TO
!          FIND NEWTON STEP FROM QR DECOMPOSITION.
!
  CALL qrsolv (overch, overfl, maxexp, n, nunit, output, a, hhpi, rdiag, sn)
  overfl=.false.
!
  IF (output > 3) THEN
    IF (.NOT.sclxch) THEN
      WRITE (nunit,10)
      WRITE (nunit,380)
      380 FORMAT ('   *       NEWTON STEP FROM QR DECOMPOSITION ', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,260) i,sn(i)
      END DO
    ELSE
      WRITE (nunit,10)
      WRITE (nunit,400)
      400 FORMAT ('   *       NEWTON STEP FROM QR DECOMPOSITION ',  &
                  '       IN SCALED UNITS', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,300) i,sn(i),i,scalex(i)*sn(i)
      END DO
    END IF
  END IF
!
!          TRANSFORM MATRICES FOR SUBSEQUENT CALCULATIONS IN TRUST
!          REGION METHODS (A IS STORED ABOVE IN H).
!
  IF (.NOT.linesr) THEN
    DO  i=1,n
      a(i,i)=rdiag(i)
      a(i,1:i-1)=a(1:i-1,i)
    END DO
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE NSTPUN.
!
END SUBROUTINE nstpun



SUBROUTINE olhelp (nunit, help)
!
!    NOV. 4, 1991
!
!    BRIEF DESCRIPTIONS OF INPUT VARIABLES FOR NNES.
!
!    INDIVIDUAL VARIABLES CAN BE ACCESSED BY ASSIGNING THAT
!    VARIABLE'S NAME TO HELP (LEFT-JUSTIFIED PLEASE).
!
!    GROUPS OF VARIABLES CAN BE ACCESSED BY ASSIGNING:
!
!    ALLLOG   ALL LOGICAL VARIABLES
!    ALLINT   ALL INTEGER VARIABLES
!    ALLREA   ALL REAL VARIABLES
!    ALLMAT   ALL MATRICES (INCLUDING VECTORS)
!
!    OR DESCRIPTIONS OF ALL VARIABLES CAN BE ACCESSED BY SETTING HELP TO ALL.
!
!    ERR STAYS TRUE IF AN INCORRECT INPUT WAS GIVEN BY THE USER.
!

INTEGER, INTENT(IN)            :: nunit
CHARACTER (LEN=6), INTENT(IN)  :: help

! Local variables

LOGICAL  :: ERR, prall, prall1, prall2, prall3, prall4
!
!    PRALL  BECOMES TRUE IF ALL VARIABLES ARE TO BE PRINTED.
!    PRALL1 BECOMES TRUE IF ALL LOGICALS  ARE TO BE PRINTED.
!    PRALL2    "      "   "  "  INTEGERS   "   "  "    "   .
!    PRALL3    "      "   "  "  REALS      "   "  "    "   .
!    PRALL4    "      "   "  "  MATRICES   "   "  "    "   .
!
ERR=.true.
prall=.false.
prall1=.false.
prall2=.false.
prall3=.false.
prall4=.false.
IF (help(1:3) == 'ALL' .AND. help(4:6) /= 'LOG' .AND. help(4:6) /= 'INT'  &
     .AND. help(4:6) /= 'REA' .AND. help(4:6) /= 'MAT') prall=.true.
IF (help == 'ALLLOG') prall1=.true.
IF (help == 'ALLINT') prall2=.true.
IF (help == 'ALLREA') prall3=.true.
IF (help == 'ALLMAT') prall4=.true.
!
!       LOGICAL VARIABLES.
!
WRITE (nunit,10)
10 FORMAT (t3,72('*'))
WRITE (nunit,20)
20 FORMAT ('   *', t74, '*')
WRITE (nunit,30)
30 FORMAT ('   *', t27, 'ON-LINE HELP FOR NNES', t74, '*')
WRITE (nunit,20)
WRITE (nunit,10)
WRITE (nunit,20)
IF (prall .OR. prall1) THEN
  WRITE (nunit,20)
  WRITE (nunit,40)
  40 FORMAT ('   ***  LOGICAL VARIABLES ***', t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,20)
END IF
IF (prall .OR. prall1 .OR. help == 'ABSNEW') THEN
  WRITE (nunit,50)
  50 FORMAT ('   ** ABSNEW **', t74, '*'/  '   *', t74, '*'/   &
             '   *  DESIGNATES WHETHER ABSOLUTE NEWTON''S METHOD IS',  &
             ' TO BE USED.', t74, '*'/  '   *', t74, '*'/   &
             '   *  DEFAULT VALUE: FALSE', t74, '*'/  '   *', t74, '*'/   &
             '   *  CROSS-REFERENCE: LINESR, NEWTON', t74, '*')
  WRITE (nunit,20)
  ERR=.false.
  IF (help == 'ABSNEW') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'CAUCHY') THEN
  WRITE (nunit,60)
  60 FORMAT ('   ** CAUCHY **', t74, '*'/  '   *', t74, '*'/   &
             '   *  THE CAUCHY POINT IS THE POINT AT THE MINIMUM OF THE',   &
             ' QUADRATIC ', t74, '*'/   &
             '   *  MODEL IN THE STEEPEST DESCENT DIRECTION (THE DISTANCE ',  &
             'FROM THE', t74, '*'/   &
             '   *  CURRENT POINT TO THE CAUCHY POINT IS ALWAYS < THE LENGTH', &
             t74, '*'/ '   * OF THE NEWTON STEP).', t74, '*'/  &
             '   *', t74, '*'/   &
             '   *  CAUCHY:  TRUE =>  INITIAL TRUST REGION IS DISTANCE TO ', &
             'CAUCHY POINT', t74, '*'/   &
             '   *', t14, 'FALSE => INITIAL TRUST REGION IS LENGTH OF ',   &
             'NEWTON STEP', t74, '*'/  '   *', t74, '*'/   &
             '   *  DEFAULT VALUE: FALSE (MORE AMBITIOUS; TRUE IS ',   &
             'CONSERVATIVE)', t74, '*'/  '   *', t74, '*'/   &
             '   *  NOT USED IN LINE SEARCH METHODS', t74, '*'/   &
             '   *', t74, '*'/   &
             '   *  CROSS-REFERENCE: LINESR,DELTA', t74, ' *')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'CAUCHY') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'DEUFLH') THEN
  WRITE (nunit,70)
  70 FORMAT ('   ** DEUFLH **', t74, '*'/  '   *', t74, '*'/   &
             '   *  INITIALIZATION OF THE RELAXATION FACTOR FOR LINE',   &
             ' SEARCHES USING', t74, '*'/   &
             '   *  A MODIFICATION OF THE DEUFLHARD METHOD. THE INITIAL ', &
             'LAMBDA IN ANY', t74, '*'/   &
             '   *  LINE SEARCH IS SET TO 1.0 USUALLY, UNLESS',   &
             ' IT IS FOUND INTERNALLY,', t74, '*'/   &
             '   *  USING DEUFLHARD''S METHOD, THAT A SMALL VALUE ',  &
             'OF LAMBDA IS MORE', t74, '*'/   &
             '   *  LIKELY WHEREUPON LAMBDA IS INITIALIZED TO 0.1.', t74, '*'/ &
             '   *', t74, '*'/   &
             '   *  DEUFLH:  TRUE =>  DEUFLHARD INITIALIZATION USED', t74, '*'/ &
             '   *', t14, 'FALSE => INITIAL RELAXATION FACTOR IS',   &
             ' ALWAYS LAM0', t74, '*'/ '   *', t22, ' (LAM0 IS USUALLY 1.0)', &
             t74, '*'/ '   *', t74, '*'/   &
             '   *  DEFAULT VALUE: TRUE', t74, '*'/  '   *', t74, '*'/   &
             '   *  NOT USED IN TRUST REGION METHODS', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'DEUFLH') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help(1:5) == 'GEOMS') THEN
  WRITE (nunit,80)
  80 FORMAT ('   ** GEOMS **', t74, '*'/  '   *', t74, '*'/   &
             '   *  TWO CHOICES ARE AVAILABLE WHEN REDUCING THE RELAXATION', &
             ' FACTOR IN', t74, '*'/   &
             '   *  LINE SEARCHES OR THE TRUST REGION SIZE IN TRUST REGION', &
             'METHODS', t74, '*'/  '   *', t74, '*'/   &
             '   *  GEOMS:   TRUE =>  GEOMETRIC REDUCTION; FACTOR SIGMA ',  &
             'FOR L.S .', t74, '*'/   &
             '   *', t51, 'DELFAC FOR T.R.', t74, '*'/   &
             '   *', t23, 'I.E., LAMBDA(NEW) = 0.5*LAMBDA(OLD) FOR LINE',   &
             t74, '*'/   &
             '   *', t23, 'SEARCHES AND DELTA(NEW) = 0.5*DELTA(OLD) FOR',   &
             t74, '*'/   &
             '   *', t23,' TRUST REGION METHODS.', t74, '*'/   &
             '   *', t14, 'FALSE => SAFEGUARDED POLYNOMIAL INTERPOLATION',  &
             t74, '*'/  '   *',t74, '*'/   &
             '   *  DEFAULT VALUE: TRUE', t74, '*'/  '   *', t74, '*'/   &
             '   *  CROSS-REFERENCES: DELFAC, SIGMA', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'GEOMS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'LINESR') THEN
  WRITE (nunit,90)
  90 FORMAT ('   ** LINESR **', t74, '*'/  '   *', t74, '*'/   &
             '   *  DISTINGUISHES BETWEEN MAJOR SOLUTION METHODS.', t74, '*'/ &
             t3, '*', t74, '*'/   &
             '   *  LINESR:  TRUE =>  LINE SEARCH METHOD', t74, '*'/   &
             '   *', t14, 'FALSE => TRUST REGION METHOD', t74, '*'/   &
             '   *', t74, '*'/  '   *  DEFAULT VALUE: TRUE', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'LINESR') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'MATSUP') THEN
  WRITE (nunit,100)
  100 FORMAT ('   ** MATSUP **', t74, '*'/  '   *', t74, '*'/   &
              '   *  SUPPRESSES MATRIX PRINTING IN DETAILED OUTPUT.', t74, '*'/ &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: FALSE', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: OUTPUT', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'MATSUP') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'OVERCH') THEN
  WRITE (nunit,110)
  110 FORMAT ('   ** OVERCH **', t74, '*'/  '   *', t74, '*'/   &
              '   *  CHECKS FOR POTENTIAL OVERFLOWS AT KEY',   &
              ' LOCATIONS; INSERTS ', t74, '*'/   &
              '   *  (+OR-)10**MAXEXP AS AN APPROXIMATION IF AN OVERFLOW ',  &
              'IS IMMINENT.', t74, '*'/  '   *', t74, '*'/   &
              '   *  OVERCH:  TRUE => OVERFLOW CHECKING', t74, '*'/   &
              '   *', t14, 'FALSE => NORMAL EXECUTION', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: FALSE ', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: MAXEXP', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'OVERCH' .OR. help == 'ALLLOG') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall1 .OR. help == 'WRNSUP') THEN
  WRITE (nunit,120)
  120 FORMAT ('   ** WRNSUP **', t74, '*'/  '   *', t74, '*'/   &
              '   *  SUPPRESSES PRINTING WARNINGS; USED WHEN KNOWN',   &
              ' WARNINGS', t74, '*'/   &
              '   *  WILL CLUTTER OUTPUT.', t74, '*'/  '   *', t74, '*'/  &
              '   *  DEFAULT VALUE: FALSE', t74, '*'/  '   *', t74, '*'/  &
              '   *  CROSS-REFERENCE: OUTPUT', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'WRNSUP') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2) THEN
  WRITE (nunit,20)
  WRITE (nunit,130)
  130 FORMAT ('   ***  INTEGER VARIABLES ***', t74, '*')
  WRITE (nunit,20)
END IF
IF (prall .OR. prall2 .OR. help == 'ACPTCR') THEN
  WRITE (nunit,20)
  WRITE (nunit,140)
  140 FORMAT ('   ** ACPTCR **', t74, '*'/  '   *', t74, '*'/   &
              '   *  A STEP CAN BE ACCEPTED BY THE STANDARD OBJECTIVE',   &
              ' FUNCTION, WHICH', t74, '*'/   &
              '   *  IS 1/2 {SUM OF SQUARES OF RESIDUALS}, ALONE OR ',  &
              'BY DEUFLHARD''S ', t74, '*'/   &
              '   *  "NATURAL" OBJECTIVE FUNCTION AS WELL.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  ACPTCR:  1  => USE ONLY STANDARD OBJECTIVE FUNCTION',  &
              t74, '*'/   &
              t3, '*', t14, '12 => ACCEPT STEP BASED ON EITHER CRITERION',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 12', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'ACPTCR') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'ITSCLF') THEN
  WRITE (nunit,150)
  150 FORMAT ('   ** ITSCLF **', t74, '*'/  '   *', t74, '*'/   &
              '   *  ITERATION AT WHICH ADAPTIVE SCALING OF THE',   &
              ' FUNCTIONS BEGINS.', t74, '*'/  '   *', t74, '*'/   &
              '   *  ITSCLF:  0  => NO ADAPTIVE FUNCTION SCALING', t74, '*'/ &
              '   *', t14, 'K  => ADAPTIVE FUNCTION SCALING BEGINS AT',   &
              ' ITERATION K', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 0', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: SCALEF', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'ITSCLF') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'ITSCLX') THEN
  WRITE (nunit,160)
  160 FORMAT ('   ** ITSCLX **', t74, '*'/  '   *', t74, '*'/   &
              '   *  ITERATION AT WHICH ADAPTIVE SCALING OF THE',   &
              ' VARIABLES BEGINS.', t74, '*'/  '   *', t74, '*'/    &
              '   *  ITSCLX:  0  => NO ADAPTIVE VARIABLE SCALING', t74, '*'/ &
              '   *', t14, 'K  => ADAPTIVE VARIABLE SCALING BEGINS',   &
              ' AT ITERATION K', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 0', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: SCALEX', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'ITSCLX') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'JACTYP') THEN
  WRITE (nunit,170)
  170 FORMAT ('   ** JACTYP **', t74, '*'/  '   *', t74, '*'/   &
              '   *  JACTYP DETERMINES HOW THE JACOBIAN IS TO BE',   &
              ' EVALUATED EXPLICITLY', t74, '*'/   &
              '   *  (VERSUS BEING UPDATED VIA A QUASI-NEWTON METHOD). ',  &
              'LOWER AND UPPER', t74, '*'/   &
              '   *  BOUNDS ARE CHECKED TO PREVENT VIOLATIONS.', t74, '*'/  &
              '   *', t74, '*'/   &
              '   *  JACTYP:  0  => ANALYTICAL JACOBIAN (DECLARE IN ',   &
              'EXTERNAL STATEMENT)', t74, '*'/   &
              '   *', t14, '1  => FORWARD DIFFERENCES', t74, '*'/   &
              '   *', t14, '2  => BACKWARD DIFFERENCES', t74, '*'/   &
              '   *', t14, '3  => CENTRAL DIFFERENCES', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 1 ', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: JUPDM', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'JACTYP') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help(1:5) == 'JUPDM') THEN
  WRITE (nunit,180)
  180 FORMAT ('   ** JUPDM **', t74, '*'/  '   *', t74, '*'/   &
              '   *  JUPDM DETERMINES WHETHER THE JACOBIAN IS TO BE',   &
              ' EVALUATED', t74, '*'/   &
              '   *  EXPLICITLY OR TO BE UPDATED VIA A QUASI-NEWTON METHOD.', &
              t74, '*'/ '   *', t74, '*'/   &
              '   *  JUPDM:   0 => JACOBIAN EVALUATED EXPLICITLY (SEE ',   &
              'JACTYP FOR', t74, '*'/   &
              '   *', t20, 'DIFFERENCING OPTIONS)', t74, '*'/   &
              '   *', t14, '1  => BROYDEN UPDATE', t74, '*'/   &
              '   *', t14, '2  => LEE AND LEE UPDATE', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 0 ', t74, '*'/   &
              '   *', t74, '*'/   '   *  CROSS-REFERENCE: JACTYP', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'JUPDM') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'MAXEXP') THEN
  WRITE (nunit,190)
  190 FORMAT ('   ** MAXEXP **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MAXIMUM EXPONENT ALLOWED, BASE 10, DETERMINED IN',   &
              ' SUBROUTINE SETUP', t74, '*'/    &
              '   *  BY SUBROUTINE MACHAR, E.G., 38 FOR THE VAX, ',   &
              '308 FOR IBM PC''S.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: DETERMINED INTERNALLY', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: OVERCH', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'MAXEXP') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help(1:5) == 'MAXIT') THEN
  WRITE (nunit,200)
  200 FORMAT ('   ** MAXIT **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MAXIMUM NUMBER OF ITERATIONS ALLOWED', t74, '*'/   &
              '   *', t74, '*'/ '   *  DEFAULT VALUE: 100', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'MAXIT') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help(1:5) == 'MAXNS') THEN
  WRITE (nunit,210)
  210 FORMAT ('   ** MAXNS **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MAXIMUM NUMBER OF LINE SEARCH STEPS OR TRUST REGION',  &
              ' REDUCTIONS', t74, '*'/   &
              '   *  ALLOWED WHEN THE JACOBIAN HAS BEEN ',  &
              'CALCULATED EXPLICITLY, I.E.,', t74, '*'/   &
              '   *  BY FINITE DIFFERENCES OR USING A USER-SUPPLIED ',  &
              'JACOBIAN.  THIS', t74, '*'/   &
              '   *  IS USUALLY > MAXQNS, THE NUMBER OF LINE SEARCH STEPS', &
              t74, '*'/   &
              '   *  OR TRUST REGION REDUCTIONS ALLOWED AFTER THE ',  &
              'JACOBIAN HAS BEEN', t74, '*'/   &
              '   *  UPDATED USING A QUASI-NEWTON METHOD.', t74, '*'/   &
              '   *', t74, '*'/ '   *  DEFAULT VALUE: 50', t74, '*'/   &
              '   *', t74, '*'/ '   *  CROSS-REFERENCE: MAXQNS', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'MAXNS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'MAXQNS') THEN
  WRITE (nunit,220)
  220 FORMAT ('   ** MAXQNS **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MAXIMUM NUMBER OF LINE SEARCH STEPS OR TRUST REGION',  &
              ' REDUCTIONS', t74, '*'/   &
              '   *  ALLOWED WHEN THE JACOBIAN HAS BEEN UPDATED BY A ',   &
              'QUASI-NEWTON', t74, '*'/   &
              '   *  METHOD; THIS IS DESIGNED TO PREVENT EXCESSIVE LINE',   &
              ' SEARCH STEPS', t74, '*'/   &
              '   *  IN A POOR SEARCH DIRECTION.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 10', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: MAXNS', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'MAXQNS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help(1:4) == 'MGLL') THEN
  WRITE (nunit,230)
  230 FORMAT ('   ** MGLL **', t74, '*'/  '   *', t74, '*'/   &
              '   *  NUMBER OF PREVIOUS MERIT FUNCTION(S) VALUES USED',   &
              ' FOR THE', t74, '*'/   &
              '   *  NONMONOTONIC STEP ACCEPTANCE CRITERIA; ',   &
              'THE SUM OF SQUARES', t74, '*'/   &
              '   *  AND POSSIBLY DEUFLHARD''S CRITERION IS(ARE)',   &
              ' COMPARED TO THE', t74, '*'/   &
              '   *  GREATEST OF THE MOST RECENT "MGLL" VALUES STORED ',  &
              'IN THE FTRACK', t74, '*'/   &
              '   *  AND STRACK VECTORS RESPECTIVELY.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  E.G. LINE SEARCHES, THE ACCEPTANCE CRITERION FOR A ',  &
              'FUNCTION, F, IS', t74, '*'/  '   *', t74, '*'/   &
              '   *   F(NEW) <= FMAX + LAMBDA*ALPHA*DELF^S   WHERE', t74, '*'/ &
              '   *', t74, '*'/  &
              '   *   F      IS THE OBJECTIVE FUNCTION', t74, '*'/   &
              '   *   LAMBDA IS THE RELAXATION FACTOR', t74, '*'/   &
              '   *   ALPHA IS THE ARMIJO CONSTANT', t74, '*'/   &
              '   *   DELF   IS THE GRADIENT OF F', t74, '*'/   &
              '   *   S      IS THE PROPOSED STEP', t74, '*'/   &
              '   *   FMAX   IS GIVEN BY', t74, '*'/  '   *', t74, '*'/   &
              '   *   MAX(F(K-J)), 0 <= J <= M(K)   WHERE', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *   M(K) = MIN[M(K-1)+1,MGLL]', t74, '*'/  &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: GIVEN BY USER AS A PARAMETER, ',   &
              '10 IS TYPICAL', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: NARMIJ,ALPHA', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:4) == 'MGLL') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'MINQNS') THEN
  WRITE (nunit,240)
  240 FORMAT ('   ** MINQNS **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MINIMUM NUMBER OF QUASI-NEWTON STEPS WHICH MUST',   &
              ' BE TAKEN ', t74, '*'/   &
              '   *  BETWEEN EXPLICIT JACOBIAN EVALUATIONS.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 6', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: JACTYP,JUPDM,RATIOF ', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'MINQNS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'NARMIJ') THEN
  WRITE (nunit,250)
  250 FORMAT ('   ** NARMIJ **', t74, '*'/  '   *', t74, '*'/   &
              '   *  WHEN THE JACOBIAN IS EVALUATED EXPLICITLY AT',   &
              ' EACH ITERATION, IT', t74, '*'/   &
              '   *  IS THE NUMBER OF STEPS WHICH MUST SATISFY THE ',   &
              'ARMIJO CRITERION AT', t74, '*'/   &
              '   *  THE START OF THE PROBLEM, I.E., STRICT DESCENT',   &
              ' STEPS AT THE ', t74, '*'/   &
              '   *  START. IN QUASI-NEWTON METHODS IT IS THE NUMBER ',  &
              'OF STEPS WHICH', t74, '*'/   &
              '   *  MUST SATISFY THE ARMIJO CRITERION AFTER EACH EXPLICIT', &
              ' JACOBIAN', t74, '*'/   &
              '   *  EVALUATION. THE ARMIJO CRITERION TESTS WHETHER', t74, '*'/ &
              '   *', t74, '*'/   &
              '   *   F(NEW) <= F(OLD) + LAMBDA*ALPHA*DELF^S   WHERE', t74, '*'/ &
              '*', t74, '*'/   &
              '   *   F      IS THE OBJECTIVE FUNCTION', t74, '*'/   &
              '   *   LAMBDA IS THE RELAXATION FACTOR', t74, '*'/   &
              '   *   ALPHA  IS THE ARMIJO CONSTANT', t74, '*'/   &
              '   *   DELF   IS THE GRADIENT OF F', t74, '*'/   &
              '   *   S      IS THE PROPOSED STEP', t74, '*'/   &
              '   *', t74, '*'/ '   *  DEFAULT VALUE: 1', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: MGLL,ALPHA', t74, '* ')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'NARMIJ') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'NFETOT') THEN
  WRITE (nunit,260)
  260 FORMAT ('   ** NFETOT **', t74, '*'/  '   *', t74, '*'/   &
              '   *  TOTAL NUMBER OF FUNCTION EVALUATIONS INCLUDING',   &
              ' THOSE REQUIRED', t74, '*'/   &
              '   *  FOR THE FINITE-DIFFERENCE CALCULATION OF JACOBIANS.', &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE : N/A (OUTPUT)', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'NFETOT') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'NINITN') THEN
  WRITE (nunit,270)
  270 FORMAT ('   ** NINITN **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MAY BE USED TO DELAY ONSET OF LINE SEARCHES OR',   &
              ' TRUST REGION STEPS.', t74, '*'/   &
              '   *  IT IS THE NUMBER OF INITIAL NEWTON STEPS BEFORE ',   &
              'LINE SEARCHES OR', t74, '*'/   &
              '   *  TRUST REGION STEPS BEGIN.', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 0', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'NINITN') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'NJACCH') THEN
  WRITE (nunit,280)
  280 FORMAT ('   ** NJACCH **', t74, '*'/  '   *', t74, '*'/   &
              '   *  THE NUMBER TIMES AN ANALYTICAL JACOBIAN IS CHECKED',  &
              ' USING FINITE', t74, '*'/   &
              '   *  DIFFERENCES AT THE START OF THE PROBLEM. THE ',   &
              'TOLERANCE USED FOR', t74, '*'/   &
              '   *  COMPARISON IS GIVEN BY FDTOLJ.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE : 0', t74, '*'/ '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: FDTOLJ', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'NJACCH') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help(1:5) == 'NUNIT') THEN
  WRITE (nunit,290)
  290 FORMAT ('   ** NUNIT **', t74, '*'/  '   *', t74, '*'/   &
              '   *  NUMBER OF THE LOGICAL UNIT FOR OUTPUT.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: GIVEN BY USER', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'NUNIT') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'OUTPUT') THEN
  WRITE (nunit,300)
  300 FORMAT ('   ** OUTPUT **', t74, '*'/  '   *', t74, '*'/   &
              '   *  OUTPUT DETERMINES THE DETAIL OF THE OUTPUT.', t74, '*'/ &
              '   *', t74, '*'/   &
              '   *  OUTPUT:  0  => NO OUTPUT', t74, '*'/   &
              '   *', t14, '1  => ANSWER ONLY', t74, '*'/   &
              '   *', t14, '2  => ECHO INPUT PLUS ANSWER', t74, '*'/   &
              '   *', t14, '3  => SUMMARY OF EACH ITERATION', t74, '*'/  &
              '   *', t14, '4  => DETAILED DESCRIPTION OF EACH',   &
              ' ITERATION', t74, '*'/   &
              '   *', t14, '5  => VERY DETAILED DESCRIPTION OF',  &
              ' EACH ITERATION', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 2', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'OUTPUT') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'QNUPDM') THEN
  WRITE (nunit,310)
  310 FORMAT ('   ** QNUPDM **', t74, '*'/  '   *', t74, '*'/   &
              '   *  QNUPDM DETERMINES HOW THE QUASI-NEWTON',   &
              ' UPDATE IS DONE.', t74, '*'/  '   *', t74, '*'/   &
              '   *  QNUPDM:  0  => UPDATE UNFACTORED JACOBIAN', t74, '*'/  &
              '   *', t14, '1  => UPDATE QR DECOMPOSITION OF JACOBIAN',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  QNUPDM=1 IS FASTER. USE THIS UNLESS YOU WANT TO SEE THE', &
              t74, '*'/   &
              '   *  JACOBIAN AT EACH ITERATION.', t74, '*'/  '   *', t74, '*'/ &
              '   *  DEFAULT VALUE: 1', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'QNUPDM') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'STOPCR') THEN
  WRITE (nunit,320)
  320 FORMAT ('   ** STOPCR **', t74, '*'/  '   *', t74, '*'/   &
              '   *  STOPCR DETERMINES THE STOPPING CRITERIA USED.', t74, '*'/ &
              t3, '*', t74, '*'/   &
              '   *  STOPCR: 1  => STOP BASED ON STEP SIZE ONLY', t74, '*'/   &
              '   *', t14, '2  => STOP BASED ON OBJECTIVE FUNCTION VALUE ONLY', &
              t74, '*'/   &
              '   *', t14, '12 => STOP BASED ON EITHER STEP SIZE OR',   &
              ' FUNCTION VALUE', t74, '*'/   &
              '   *', t14, '3  => STOP BASED ON BOTH BEING SATISFIED', t74, '*'/ &
              '   *', t74, '*'/   &
              '   *  SEE FTOL, NSTTOL AND STPTOL FOR DETAILS OF THE CRITERIA.',  &
              ' THERE ARE', t74, '*'/    &
              '   *  TWO STEP SIZE CRITERIA: THE ONE BASED ON THE FULL NEWTON ',  &
              'STEP IS', t74, '*'/   &
              '   *  GOVERNED BY NSTTOL, AND THE OTHER BASED ON THE',   &
              ' STEP SIZE AFTER', t74, '*'/   &
              '   *  THE LINE SEARCH OR TRUST REGION REDUCTION IS GOVERNED ', &
              'BY STPTOL.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 12', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: TRMCOD,FTOL,NSTTOL,STPTOL', t74, '*'/ &
              '   *', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'STOPCR') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'SUPPRS') THEN
  WRITE (nunit,330)
  330 FORMAT ('   ** SUPPRS **', t74, '*'/  '   *', t74, '*'/   &
              '   *  SUPPRESS DETAILED OUTPUT FOR "SUPPRS" ITERATIONS;',   &
              ' USED PRIMARILY', t74, '*'/   &
              '   *  TO SEE DETAILED OUTPUT BEFORE A FAILURE IN A LARGE ',  &
              'PROBLEM.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 0', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: OUTPUT', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'SUPPRS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'TRMCOD') THEN
  WRITE (nunit,340)
  340 FORMAT ('   ** TRMCOD **', t74, '*'/  '   *', t74, '*'/   &
              '   *  TRMCOD TELLS WHICH STOPPING CRITERIA WERE MET.',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  TRMCOD: 1 => STOP BASED ON STEP SIZE ONLY', t74, '*'/  &
              '   *', t14, '2  => STOP BASED ON OBJECTIVE FUNCTION VALUE ONLY', &
              t74, '*'/  '   *', t14, '12 => STOP BASED ON EITHER ',  &
              'STEP SIZE OR FUNCTION VALUE', t74, '*'/   &
              '   *', t14, '3  => STOP BASED ON BOTH BEING SATISFIED', t74, '*'/ &
              '   *', t74, '*'/   &
              '   *  SEE FTOL, NSTTOL AND STPTOL FOR DETAILS OF THE CRITERIA.', &
              ' THERE ARE', t74, '*'/   &
              '   *  TWO STEP SIZE CRITERIA: THE ONE BASED ON THE FULL ',  &
              'NEWTON STEP IS', t74, '*'/   &
              '   *  GOVERNED BY NSTTOL, AND THE OTHER BASED ON THE',   &
              ' STEP SIZE AFTER', t74, '*'/   &
              '   *  THE LINE SEARCH OR TRUST REGION REDUCTION IS GOVERNED ', &
              'BY STPTOL.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: N/A (OUTPUT)', t74, '*'/  '   *', t74, '*'/ &
              '   *  CROSS-REFERENCE: STOPCR,FTOL,NSTTOL,STPTOL', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'TRMCOD') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall2 .OR. help == 'TRUPDM') THEN
  WRITE (nunit,350)
  350 FORMAT ('   ** TRUPDM **', t74, '*'/  '   *', t74, '*'/   &
              '   *  TRUPDM DETERMINES THE TRUST REGION UPDATING',   &
              ' METHOD.', t74, '*'/  '   *', t74, '*'/   &
              '   *  TRUPDM:  0  => POWELL''S SCHEME', t74, '*'/   &
              '   *', t15, '1 => DENNIS AND SCHNABEL''S SCHEME', t74, '*'/  &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 0', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'TRUPDM' .OR. help == 'ALLINT') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3) THEN
  WRITE (nunit,20)
  WRITE (nunit,360)
  360 FORMAT ('   ***  REAL VARIABLES ***', t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,20)
END IF
IF (prall .OR. prall3 .OR. help(1:5) == 'ALPHA') THEN
  WRITE (nunit,370)
  370 FORMAT ('   ** ALPHA **', t74, '*'/  '   *', t74, '*'/   &
              '   *  ARMIJO CONSTANT. FOR MONOTONIC LINE SEARCHES,',   &
              ' A STEP IS ACCEPTED IF', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *   F(NEW) <= F(OLD) + LAMBDA*ALPHA*(DIRECTIONAL DERIVATIVE)', &
              '   WHERE', t74, '*'/  '   *', t74, '*'/   &
              '   *  F      IS THE OBJECTIVE FUNCTION', t74, '*'/   &
              '   *  LAMBDA IS THE RELAXATION FACTOR IN THE LINE SEARCH',  &
              t74, '*'/   &
              '   *  THE DERIVATIVE IS IN THE LINE SEARCH DIRECTION.',   &
              t74, '*'/  '   *', t74, '*'/    &
              '   *  SIMILARLY, FOR TRUST REGION METHODS, THE CRITERION IS', &
              t74, '*'/  '   *', t74, '*'/   &
              '   *   F(NEW) <= F(OLD) + ALPHA*(DIRECTIONAL DERIVATIVE)   ',  &
              'WHERE' , t74, '*'/  '   *', t74, '*'/   &
              '   *  THE DERIVATIVE IS IN THE TRUST REGION STEP DIRECTION.',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  THE CRITERION IS DIFFERENT FOR NONMONOTONIC SEARCHES, ',  &
              'SEE MGLL.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 1.0D-04', t74, '*'/ '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: MGLL ,NARMIJ', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'ALPHA') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'CONFAC') THEN
  WRITE (nunit,380)
  380 FORMAT ('   ** CONFAC **', t74, '*'/  '   *', t74, '*'/   &
              '   *  CONSTRAINT FACTOR WHICH GIVES THE FRACTION OF',   &
              ' THE DISTANCE', t74, '*'/   &
              '   *  TO THE FIRST VIOLATED CONSTRAINT AT WHICH A LINE ',  &
              'SEARCH WOULD', t74, '*'/   &
              '   *  START OR A TRUST REGION LIMIT SET.', t74, '*'/   &
              '   *', t74, '*'/ '   *  DEFAULT VALUE: 0.95', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'CONFAC') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:5) == 'DELTA') THEN
  WRITE (nunit,390)
  390 FORMAT ('   ** DELTA **', t74, '*'/  '   *', t74, '*'/   &
              '   *  DELTA IS THE TRUST REGION RADIUS. IF IT IS',   &
              ' NEGATIVE ON', t74, '*'/   &
              '   *  INPUT, THE INITIAL TRUST REGION SIZE IS CALCULATED ',  &
              'INTERNALLY;', t74, '*'/   &
              '   *  SEE CAUCHY. A POSITIVE ENTRY SETS THE INITIAL',   &
              ' TRUST REGION.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: -1.0', t74, '*'/  '   *', t74, '*'/   &
              t3, '*  CROSS-REFERENCE: CAUCHY,LINESR', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'DELTA') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'DELFAC') THEN
  WRITE (nunit,400)
  400 FORMAT ('   ** DELFAC **', t74, '*'/  '   *', t74, '*'/   &
              '   *  DELFAC IS THE FACTOR BY WHICH THE TRUST REGION',   &
              ' RADIUS IS CHANGED,', t74, '*'/   &
              '   *  BOTH WHEN INCREASED AND WHEN DECREASED, IF DELTA ',  &
              'IS NOT BEING', t74, '*'/   &
              '   *  BEING INCREASED TO THE LENGTH OF THE NEWTON STEP.',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 2.0D0', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: TRUPDM,DELTA', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'DELFAC') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'EPSMCH') THEN
  WRITE (nunit,410)
  410 FORMAT ('   ** EPSMCH **', t74, '*'/  '   *', t74, '*'/   &
              '   *  MACHINE PRECISION; AN ESTIMATE OF THE',   &
              ' SMALLEST FLOATING POINT', t74, '*'/   &
              '   *  NUMBER SUCH THAT 1.0+X > 1.0.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: CALCULATED INTERNALLY', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'EPSMCH') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'ETAFAC') THEN
  WRITE (nunit,420)
  420 FORMAT ('   ** ETAFAC **', t74, '*'/  '   *', t74, '*'/   &
              '   *  FACTOR USED IN DETERMINING THE SHAPE',   &
              ' OF THE DOUBLE DOGLEG STEP', t74, '*'/   &
              '   *  IN TRUST REGION METHODS; ETAFAC=1 CORRESPONDS TO ',  &
              'SINGLE DOGLEG.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 0.2', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'ETAFAC') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'FCNNEW') THEN
  WRITE (nunit,430)
  430 FORMAT ('   ** FCNNEW **', t74, '*'/  '   *', t74, '*'/   &
              '   *  ON RETURN, FCNNEW HOLDS THE FINAL VALUE OF THE ',  &
              'SUM-OF-SQUARES', t74, '*'/   &
              '   *  OBJECTIVE FUNCTION.', t74, '*'/ '   *' , t74, '*'/   &
              '   *  DEFAULT VALUE: N/A (OUTPUT)', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'FCNNEW') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'FDTOLJ') THEN
  WRITE (nunit,440)
  440 FORMAT ('   ** FDTOLJ **', t74, '*'/   '   *', t74, '*'/   &
              '   *  TOLERANCE FOR A FINITE-DIFFERENCE CHECK OF',   &
              ' AN ANALYTICAL JACOBIAN.', t74, '*'/   &
              '   *  IF JACFD(I) IS THE FINITE-DIFFERENCE APPROXIMATION ',  &
              'AND JACAN(I)', t74, '*'/   &
              '   *  IS THE ANALYTICAL DERIVATIVE, WHEN', t74, '*'/   &
              '   *', t74, '*'/    &
              '   *   ABS(JACFD(I)-JACAN(I))/MAX(ABS(JACFD(I)),1) ',   &
              '>= FDTOLJ', t74, '*'/ '   *', t74, '*'/   &
              '   *  THEN A FAILURE IS DECLARED.', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 1.0D-06', t74, '*'/  &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: NJACCH', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'FDTOLJ') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:4) == 'FTOL') THEN
  WRITE (nunit,450)
  450 FORMAT ('   ** FTOL **', t74, '*'/  '   *', t74, '*'/   &
              '   *  STOPPING TOLERANCE FOR MAX-NORM OF SCALED',   &
              ' FUNCTION VECTOR; IF', t74, '*'/  '   *', t74, '*'/   &
              '   *   MAX(SCALEF(I)*ABS(FVEC(I)) ) I=1,...,N <  FTOL ,STOP.', &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: EPSMCH**(1/3)', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: STOPCR,TRMCOD', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:4) == 'FTOL') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:4) == 'LAM0') THEN
  WRITE (nunit,460)
  460 FORMAT ('   ** LAM0 **', t74, '*'/  '   *', t74, '*'/   &
              '   *  USED TO SET THE INITIAL RELAXATION FACTOR',   &
              ' IN LINE SEARCHES', t74, '*'/   &
              '   *  TO A VALUE LESS THAN 1.0 FOR EXTREMELY NONLINEAR PROBLEMS', &
              t74, '*'/  '   *  THIS OVERRIDES DEUFLHARD INITIALIZATION.',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 1.0', t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: DEUFLH', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:4) == 'LAM0') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:5) == 'MSTPF') THEN
  WRITE (nunit,470)
  470 FORMAT ('   ** MSTPF **', t74, '*'/  '   *', t74, '*'/   &
              '   *  FACTOR USED TO SET THE MAXIMUM STEP SIZE',   &
              ' ALLOWED. USUALLY THE', t74, '*'/   &
              '   *  MAXIMUM STEP SIZE IS SET BY NNES IS MUCH TOO LARGE',  &
              ' TO HAVE ANY', t74, '*'/   &
              '   *  EFFECT, BUT IN SOME CASES THE USER MAY',   &
              ' NEED TO RESTRICT POSSIBLY', t74, '*'/   &
              '   *  FATAL STEPS. THE MAXIMUM STEP IS SET BY:', t74, '*'/ &
              ' *', t74, '*'/   &
              '   *  MAXSTP = MSTPF*MAX(NORM1,NORM2) WHERE', t74, '*'/    &
              '   *', t74, '*'/   &
              '   *   NORM1 = 2-NORM OF SCALED STARTING ESTIMATES', t74, '*'/ &
              '   *   NORM2 = 2-NORM OF COMPONENT SCALING FACTORS', t74, '*'/ &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 1000.', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'MSTPF') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'NSTTOL') THEN
  WRITE (nunit,480)
  480 FORMAT ('   ** NSTTOL **', t74, '*'/  '   *', t74, '*'/   &
              '   *  STOPPING TOLERANCE FOR FULL NEWTON STEP; STOP',   &
              ' IF, FOR ALL I,', t74, '*'/  '   *', t74, '*'/   &
              '   *   MAX(ABS[SN(I)]/MAX(ABS[X(I)],1/SCALEX(I))',   &
              ' < NSTTOL*(1+NORM(DX(I)))', t74, '*'/  '   *', t74, '*'/   &
              '   *   SN(I)     IS THE I(TH) COMPONENT OF THE NEWTON STEP', &
              t74, '*'/   &
              '   *   SCALEX(I) IS THE COMPONENT SCALING FACTOR', t74, '*'/  &
              '   *   DX(I)     IS SCALEX(I)*X(I)', t74, '*'/   &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: EPSMCH**(2/3)',   &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: STOPCR,TRMCOD,SCALEX', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'NSTTOL') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:5) == 'OMEGA') THEN
  WRITE (nunit,490)
  490 FORMAT ('   ** OMEGA **', t74, '*'/  '   *', t74, '*'/   &
              '   *  FACTOR IN THE LEE AND LEE QUASI-NEWTON UPDATES.',  &
              t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 0.1', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: JUPDM', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'OMEGA') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'RATIOF') THEN
  WRITE (nunit,500)
  500 FORMAT ('   ** RATIOF **', t74, '*'/  '   *', t74, '*'/   &
              '   *  FACTOR USED IN QUASI-NEWTON METHODS TO DECIDE',   &
              ' WHETHER AN', t74, '*'/   &
              '   *  EXPLICIT JACOBIAN UPDATE SHOULD BE DONE. IF', t74, '*'/ &
              '   *', t74, '*'/   &
              '   *   F(NEW) > RATIOF*F(OLD)   THEN:  NFAIL=NFAIL+1', t74, '*'/ &
              '   *', t74, '*'/  '   *  ELSE IF', t74, '*'/ '   *', t74, '*'/   &
              '   *   F(NEW) < 0.01*F(OLD)   THEN : NFAIL=NFAIL-1', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  IF NFAIL > MINQNS, NNES RESTARTS AT THE BEST POINT FOUND', &
              ' SO FAR', t74, '*'/   &
              '   *  AS IF THAT POINT WERE A NEW INITIAL ESTIMATE.', t74, '*'/ &
              '   *', t74, '*'/  '   *  DEFAULT VALUE: 7.0D-01', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: MINQNS', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'RATIOF') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help(1:5) == 'SIGMA') THEN
  WRITE (nunit,510)
  510 FORMAT ('   ** SIGMA **', t74, '*'/ '   *', t74, '*'/   &
              '   *  FACTOR USED IN GEOMETRIC STEP REDUCTIONS TO',   &
              ' DECREASE THE', t74, '*'/   &
              '   *  RELAXATION FACTOR IN LINE SEARCHES:',t74, '*'/   &
              '   *', t74, '*'/   &
              '   *   LAMBDA(NEW) = SIGMA*LAMBDA(OLD)', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 5.0D-01', t74, '*'/   &
              '   *', t74, '*'/  '   *  CROSS-REFERENCE: GEOMS', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:5) == 'SIGMA') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall3 .OR. help == 'STPTOL') THEN
  WRITE (nunit,520)
  520 FORMAT ('   ** STPTOL **', t74, '*'/  '   *', t74, '*'/   &
              '   *  STOPPING TOLERANCE FOR STEP AFTER REDUCTION;',   &
              ' STOP IF, FOR ALL I,', t74, '*'/  '   *', t74, '*'/   &
              '   *   MAX(ABS[S(I)]/MAX(ABS[X(I)],1/SCALEX(I))',   &
              ' < STPTOL', t74, '*'/  '   *', t74, '*'/   &
              '   *   S(I)     IS THE I(TH) COMPONENT OF THE REDUCED STEP', &
              t74, '*'/   &
              '   *   SCALEX(I) IS THE COMPONENT SCALING FACTOR',  &
              t74, '*'/   '   *', t74, '*'/    &
              '   *  DEFAULT VALUE: EPSMCH**( 2/3)', t74, '*'/   &
              '   *', t74, '*'/    &
              '   *  CROSS-REFERENCE: STOPCR, TRMCOD,SCALEX', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (prall3 .OR. help == 'STPTOL') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4) THEN
  WRITE (nunit,20)
  WRITE (nunit,530)
  530 FORMAT ('   ***  REAL VECTORS ***', t74, '*')
  WRITE (nunit,20)
  WRITE (nunit,20)
END IF
IF (prall .OR. prall4 .OR. help == 'BOUNDL') THEN
  WRITE (nunit,540)
  540 FORMAT ('   ** BOUNDL **', t74, '*'/  '   *', t74, '*'/   &
              '   *  LOWER BOUNDS FOR THE COMPONENTS.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: -10**MAXEXP', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  CROSS-REFERENCE: MAXEXP,BOUNDU', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'BOUNDL') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4 .OR. help == 'BOUNDU') THEN
  WRITE (nunit,550)
  550 FORMAT ('   ** BOUNDU **', t74, '*'/  '   *', t74, '*'/   &
              '   *  UPPER BOUNDS FOR THE COMPONENTS.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 10**MAXEXP', t74, '*'/  '   *', t74, '*'/ &
              '   *  CROSS-REFERENCE: MAXEXP,BOUNDL', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'BOUNDU') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4 .OR. help == 'SCALEF') THEN
  WRITE (nunit,560)
  560 FORMAT ('   ** SCALEF **', t74, '*'/ '   *', t74, '*'/  &
              '   *  SCALING FACTORS FOR THE FUNCTIONS. THESE ARE',  &
              ' INVERSELY', t74, '*'/   &
              '   *  PROPORTIONAL TO TYPICAL VALUES FOR EACH OF ',  &
              'THE FUNCTIONS.', t74, '*'/  '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: 1.0', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'SCALEF') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4 .OR. help == 'SCALEX') THEN
  WRITE (nunit,570)
  570 FORMAT ('   ** SCALEX **', t74, '*'/   '   *', t74, '*'/   &
              '   *  SCALING FACTORS FOR THE COMPONENTS.  THESE ARE',   &
              ' INVERSELY', t74, '*'/   &
              '   *  PROPORTIONAL TO TYPICAL',   &
              ' VALUES FOR EACH OF THE COMPONENTS.', t74, '*'/   &
              '   *', t74, '*'/ '   *  DEFAULT VALUE: 1.0', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help == 'SCALEX') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4 .OR. help(1:2) == 'XC') THEN
  WRITE (nunit,580)
  580 FORMAT ('   ** XC **', t74, '*'/ '   *', t74, '*'/   &
              '   *  CONTAINS INITIAL ESTIMATE ON ENTRY.', t74, '*'/   &
              '   *', t74, '*'/   &
              '   *  DEFAULT VALUE: GIVEN BY USER', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (help(1:2) == 'XC') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
IF (prall .OR. prall4 .OR. help(1:5) == 'XPLUS') THEN
  WRITE (nunit,590)
  590 FORMAT ('   ** XPLUS **', t74, '*'/  '   *', t74, '*'/   &
              '   *  LATEST ESTIMATE ON RETURN.', t74, '*'/  '   *', t74, '*'/ &
              '   *  DEFAULT VALUE: N/A (OUTPUT)', t74, '*')
  ERR=.false.
  WRITE (nunit,20)
  IF (prall .OR. prall4 .OR. help(1:5) == 'XPLUS') THEN
    WRITE (nunit,10)
    RETURN
  END IF
END IF
!
!       CHECK FOR INCORRECT INPUT INTO HELP FACILITY.
!
IF (ERR) THEN
  WRITE (nunit,20)
  WRITE (nunit,600) help
  600 FORMAT ('   *  INCORRECT INPUT TO HELP FACILITY: ',a, t74, '*' )
  WRITE (nunit,20)
  WRITE (nunit,10)
END IF
!
!       LAST CARD OF SUBROUTINE OLHELP.
!
RETURN
END SUBROUTINE olhelp



SUBROUTINE onenrm (abort, pertrb, n, nunit, output, epsmch, h1norm, h, scalex)
!
!    FEB. 23, 1992
!
!    FIND 1-NORM OF H MATRIX IF PERTURBATION IS DESIRED AND PERTURB DIAGONAL.
!

LOGICAL, INTENT(OUT)       :: abort
LOGICAL, INTENT(OUT)       :: pertrb
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(OUT)     :: h1norm
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN)      :: scalex(n)

! Local variables

REAL (dp)  :: sqrtep, temp
INTEGER    :: i, j
!
sqrtep=SQRT(epsmch)
IF (output > 4) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,20)
  20 FORMAT ('   *       DIAGONAL OF MATRIX H ',  &
             '(= JAC^JAC) BEFORE BEING PERTURBED', t74, '*')
  WRITE (nunit,10)
  DO  i=1,n
    WRITE (nunit,30) i,i,h(i,i)
    30 FORMAT ('   *          H(', i3, ',', i3, ') = ', g12.3, t74, '*')
  END DO
END IF
h1norm=zero
DO  j=1,n
  h1norm=h1norm + ABS(h(1,j))/scalex(j)
END DO
h1norm=h1norm/scalex(1)
DO  i=2,n
  temp=zero
  DO  j=1,i
    temp=temp + ABS(h(j,i))/scalex(j)
  END DO
  DO  j=i+1,n
    temp=temp + ABS(h(i,j))/scalex(j)
  END DO
  h1norm=MAX(h1norm,temp/scalex(i))
END DO
IF (output > 4) THEN
  WRITE (nunit,10)
  WRITE (nunit,90) h1norm
  90 FORMAT ('   *       1-NORM OF MATRIX H: ', g11.3, t74, '*')
END IF
IF (h1norm < epsmch) THEN
  IF (output > 0) THEN
    WRITE (nunit,100)
    100 FORMAT (t3, 72('*'))
    WRITE (nunit,10)
    WRITE (nunit,110)
    110 FORMAT ('   *    PROGRAM FAILS AS 1-NORM OF JACOBIAN IS TOO SMALL',  &
                t74, '*')
    WRITE (nunit,10)
    WRITE (nunit,100)
  END IF
  abort=.true.
  RETURN
ELSE
!
!          PERTURB DIAGONAL OF MATRIX H - USE THIS TO FIND "SN".
!
  pertrb=.true.
  DO  i=1,n
    h(i,i)=h(i,i) + SQRT(DBLE(n))*sqrtep*h1norm*scalex(i)*scalex(i)
  END DO
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,130)
    130 FORMAT ('   *    PERTURBED H MATRIX', t74, '*')
    CALL matprt (n, n, n, n, nunit, h)
  END IF
  abort = .FALSE.    ! Added by AJM
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE ONENRM.
!
END SUBROUTINE onenrm



SUBROUTINE qform (n, a, hhpi, jac)
!
!    FEB. 14, 1991
!
!    FORM Q^  FROM THE HOUSEHOLDER MATRICES STORED IN
!    MATRICES A AND HHPI AND STORE IT IN JAC.
!

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: a(n,n)
REAL (dp), INTENT(IN)   :: hhpi(n)
REAL (dp), INTENT(OUT)  :: jac(n,n)

! Local variables

REAL (dp)  :: tau
INTEGER    :: j, k
!
jac = zero
DO  j=1,n
  jac(j,j) = one
END DO

DO  k=1,n-1
  IF (hhpi(k) /= zero) THEN
    DO  j=1,n
      tau = DOT_PRODUCT( a(k:n,k), jac(k:n,j) )
      tau = tau/hhpi(k)
      jac(k:n,j) = jac(k:n,j) - tau*a(k:n,k)
    END DO
  END IF
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE QFORM.
!
END SUBROUTINE qform



SUBROUTINE qmin (nunit, output, delfts, delta, deltaf, stplen)
!
!    FEB. 9, 1991
!
!    SET THE NEW TRUST REGION SIZE, DELTA, BASED ON A QUADRATIC
!    MINIMIZATION WHERE DELTA IS THE INDEPENDENT VARIABLE.
!
!    DELTAF IS THE DIFFERENCE IN THE SUM-OF-SQUARES OBJECTIVE
!    FUNCTION VALUE AND DELFTS IS THE DIRECTIONAL DERIVATIVE IN
!    THE DIRECTION OF THE CURRENT STEP, S, WHICH HAS STEP LENGTH STPLEN.
!

INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: delfts
REAL (dp), INTENT(IN OUT)  :: delta
REAL (dp), INTENT(IN)      :: deltaf
REAL (dp), INTENT(IN)      :: stplen

! Local variables

REAL (dp)  :: deltmp
REAL (dp), PARAMETER  :: point1 = 0.1_dp

!
IF (deltaf-delfts /= zero) THEN
!
!      CALCULATE DELTA WHERE MINIMUM WOULD OCCUR - DELTMP.
!      THIS IS PROVISIONAL AS IT MUST BE WITHIN CERTAIN LIMITS TO BE ACCEPTED.
!
  deltmp=-delfts*stplen/(two*(deltaf-delfts))
  IF (output > 4) THEN
    WRITE (nunit,10)
    10 FORMAT ('   *', t74, '*')
    WRITE (nunit,20) deltmp
    20 FORMAT ('   *       TEMPORARY DELTA FROM QUADRATIC',   &
               ' MINIMIZATION: ', g12.3, t74, '*')
    WRITE (nunit,30) delta
    30 FORMAT ('   *', t33, 'VERSUS CURRENT DELTA: ', g12.3, t74, '*')
  END IF
!
!          REDUCE DELTA DEPENDING ON THE MAGNITUDE OF DELTMP.
!          IT MUST BE WITHIN [.1DELTA,.5DELTA] TO BE ACCEPTED -
!          OTHERWISE THE NEAREST ENDPOINT OF THE INTERVAL IS USED.
!
  IF (deltmp < point1*delta) THEN
    delta=point1*delta
    IF (output > 4) THEN
      WRITE (nunit,10)
      WRITE (nunit,40)
      40 FORMAT ('   *       NEW DELTA SET TO 0.1 CURRENT DELTA', t74, '*')
    END IF
  ELSE IF (deltmp > delta/two) THEN
    delta=delta/two
    IF (output > 4) THEN
      WRITE (nunit,10)
      WRITE (nunit,50)
      50 FORMAT ('   *       NEW DELTA SET TO 0.5 CURRENT DELTA', t74, '*')
    END IF
  ELSE
    delta=deltmp
    IF (output > 4) THEN
      WRITE (nunit,10)
      WRITE (nunit,60)
      60 FORMAT ('   *       NEW DELTA SET TO DELTMP', t74, '*')
    END IF
  END IF
ELSE
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,70)
    70 FORMAT ('   *       TO AVOID OVERFLOW NEW DELTA',   &
               ' SET TO 0.5 CURRENT DELTA', t74, '*')
  END IF
  delta=delta/two
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE QMIN.
!
END SUBROUTINE qmin



SUBROUTINE qrdcom (qrsing, n, epsmch, a, hhpi, rdiag)
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE COMPUTES THE QR DECOMPOSITION OF THE MATRIX A.
!    THE DECOMPOSITION IS COMPLETED EVEN IF A SINGULARITY IS DETECTED
!    (WHEREUPON QRSING IS SET TO TRUE).
!

LOGICAL, INTENT(OUT)       :: qrsing
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(OUT)     :: hhpi(n)
REAL (dp), INTENT(OUT)     :: rdiag(n)

! Local variables

REAL (dp)  :: eta, sigma, tau
INTEGER    :: i, j, k
!
qrsing=.false.
!
DO  k=1,n-1
  eta=zero
  DO  i=k,n
    eta=MAX(eta, ABS(a(i,k)))
  END DO
  IF (eta < epsmch) THEN
    qrsing=.true.
    hhpi(k)=zero
    rdiag(k)=zero
  ELSE
    a(k:n,k)=a(k:n,k)/eta
    sigma = SUM( a(k:n,k)**2 )
    sigma=SIGN(SQRT(sigma),a(k,k))
    a(k,k)=a(k,k) + sigma
    hhpi(k)=sigma*a(k,k)
    rdiag(k)=-eta*sigma
    DO  j=k+1,n
      tau = DOT_PRODUCT( a(k:n,k), a(k:n,j) )
      tau=tau/hhpi(k)
      a(k:n,j) = a(k:n,j) - tau*a(k:n,k)
    END DO
  END IF
END DO
rdiag(n)=a(n,n)
IF (ABS(rdiag(n)) < epsmch) qrsing=.true.
RETURN
!
!       LAST CARD OF SUBROUTINE QRDCOM.
!
END SUBROUTINE qrdcom



SUBROUTINE qrsolv(overch, overfl, maxexp, n, nunit, output, a, hhpi, rdiag, b)
!
!    FEB. 2, 1991
!
!    THIS SUBROUTINE SOLVES
!
!         (QR)X=B
!
!         WHERE  Q AND R ARE OBTAINED FROM THE QR DECOMPOSITION
!                B IS A GIVEN RIGHT HAND SIDE WHICH IS OVERWRITTEN
!
!                R IS CONTAINED IN THE STRICT UPPER TRIANGLE OF
!                  MATRIX A AND THE VECTOR RDIAG
!                Q IS "CONTAINED" IN THE LOWER TRIANGLE OF MATRIX A
!
!    FRSTOV  INDICATES FIRST OVERFLOW - USED ONLY TO SET BORDER FOR OUTPUT
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: a(n,n)
REAL (dp), INTENT(IN)      :: hhpi(n)
REAL (dp), INTENT(IN)      :: rdiag(n)
REAL (dp), INTENT(IN OUT)  :: b(n)

! Local variables

LOGICAL    :: frstov
REAL (dp)  :: eps, tau
INTEGER    :: i, j

!
eps=ten**(-maxexp)
frstov=.true.
overfl=.false.
!
!    MULTIPLY RIGHT HAND SIDE BY Q^ THEN SOLVE USING R STORED IN MATRIX A.
!
DO  j=1,n-1
  tau=zero
  DO  i=j,n
    IF (overch) THEN
      IF (LOG10(ABS(a(i,j))+eps) + LOG10(ABS(b(i))+eps) -  &
            LOG10(hhpi(j)+eps) > maxexp) THEN
        overfl=.true.
        tau=SIGN(ten**maxexp, a(i,j))*SIGN(one,b(j))
        CYCLE
      END IF
    END IF
    tau=tau + a(i,j)*b(i)/hhpi(j)
  END DO

  DO  i=j,n
    IF (overch) THEN
      IF (LOG10(ABS(tau)+eps) + LOG10(ABS(a(i,j))+eps) > maxexp) THEN
        overfl=.true.
        b(i)=-SIGN(ten**maxexp,tau)*SIGN(one,a(i,j))
        IF (output > 2 .AND. (.NOT.wrnsup)) THEN
          IF (frstov) THEN
            WRITE (nunit,30)
            30 FORMAT ('   *', t74, '*')
          END IF
          WRITE (nunit,40) i,b(i)
          40 FORMAT ('   *    WARNING: COMPONENT ', i3, ' SET TO ',  &
                     g11.3,' IN QRSOLV BEFORE RSOLV', t74, '*')
        END IF
        CYCLE
      END IF
    END IF
    b(i)=b(i) - tau*a(i,j)
  END DO
END DO
CALL rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, b)
RETURN
!
!       LAST CARD OF SUBROUTINE QRSOLV.
!
END SUBROUTINE qrsolv



SUBROUTINE qrupda (overfl, maxexp, n, epsmch, a, jac, u, v)
!
!    FEB. 12, 1991
!
!    UPDATE QR DECOMPOSITION USING A SERIES OF GIVENS ROTATIONS.
!

LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: jac(n,n)
REAL (dp), INTENT(IN OUT)  :: u(n)
REAL (dp), INTENT(IN)      :: v(n)

! Local variables

REAL (dp)  :: eucnrm, hold(2)
INTEGER    :: i, k, l, ldhold
!
!    REPLACE SUBDIAGONAL WITH ZEROS SO THAT WHEN R IS MULTIPLIED BY GIVENS
!    (JACOBI) ROTATIONS THE SUBDIAGONAL ELEMENTS DO NOT AFFECT THE OUTCOME.
!
DO  i=2,n
  a(i,i-1)=zero
END DO
!
!       FIND LARGEST K FOR WHICH U(K) DOES NOT EQUAL ZERO.
!
k=n
DO  l=1,n
  IF (u(k) == zero) THEN
    IF (k > 1) THEN
      k=k-1
    ELSE
      EXIT
    END IF
  ELSE
    EXIT
  END IF
END DO
!
!       MULTIPLY UV^ BY A SERIES OF ROTATIONS SO THAT ALL BUT THE
!       TOP ROW IS MADE ZERO (THEORETICALLY THIS IS WHAT HAPPENS
!       ALTHOUGH THIS MATRIX ISN'T ACTUALLY FORMED).
!
DO  i=k-1,1,-1
  CALL jacrot (overfl, i, maxexp, n, u(i), u(i+1), epsmch, a, jac)
  IF (u(i) == zero) THEN
!
!             THIS STEP JUST AVOIDS ADDING ZERO.
!
    u(i)=ABS(u(i+1))
  ELSE
    hold(1)=u(i)
    hold(2)=u(i+1)
    ldhold=2
    CALL twonrm (overfl, maxexp, ldhold, epsmch, eucnrm, hold)
    u(i)=eucnrm
  END IF
END DO
!
!       ADD THE TOP ROW TO THE TOP ROW OF A - THIS FORMS THE
!       UPPER HESSENBERG MATRIX.
!
a(1,1:n)=a(1,1:n) + u(1)*v(1:n)
!
!       FORM THE UPPER TRIANGULAR R MATRIX BY A SERIES OF ROTATIONS
!       TO ZERO OUT THE SUBDIAGONALS.
!
DO  i=1,k-1
  CALL jacrot (overfl, i, maxexp, n, a(i,i), a(i+1,i), epsmch, a, jac)
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE QRUPDA.
!
END SUBROUTINE qrupda



SUBROUTINE rcdprt (nunit, retcod, delta, rellen, stptol)
!
!    FEB. 14, 1991
!
!    DESCRIBE MEANING OF RETURN CODES, RETCOD, FROM TRUST REGION UPDATING.
!

INTEGER, INTENT(IN)    :: nunit
INTEGER, INTENT(IN)    :: retcod
REAL (dp), INTENT(IN)  :: delta
REAL (dp), INTENT(IN)  :: rellen
REAL (dp), INTENT(IN)  :: stptol

!
WRITE (nunit,10)
10 FORMAT ('   *', t74, '*')
WRITE (nunit,20) retcod
20 FORMAT ('   *       RETCOD, FROM TRUST REGION UPDATING:', i5, t74, '*' )
WRITE (nunit,10)
IF (retcod == 1) THEN
  WRITE (nunit,30)
  30 FORMAT ('   *       PROMISING STEP FOUND; DELTA',   &
             ' HAS BEEN INCREASED TO NEWLEN BUT', t74, '*')
  WRITE (nunit,40)
  40 FORMAT ('   *       BECAUSE OF OVERFLOWS IN THE FUNCTION',   &
             ' VECTOR(S) IN SUBSEQUENT', t74, '*')
  WRITE (nunit,50)
  50 FORMAT ('   *       STEP(S) THE PROJECTED DELTA IS LESS',   &
             ' THAN THAT AT THE ALREADY SUCCESSFUL', t74, '*')
  WRITE (nunit,60)
  60 FORMAT ('   *       STEP - RETURN TO SUCCESSFUL STEP',   &
             ' AND ACCEPT AS NEW POINT', t74, '*')
ELSE IF (retcod == 2) THEN
  WRITE (nunit,30)
  WRITE (nunit,70)
  70 FORMAT ('   *       BECAUSE OF OVERFLOWS IN THE OBJECTIVE',   &
             ' FUNCTION IN SUBSEQUENT', t74, '*')
  WRITE (nunit,50)
  WRITE (nunit,60)
ELSE IF (retcod == 3) THEN
  WRITE (nunit,30)
  WRITE (nunit,80)
  80 FORMAT ('   *       BECAUSE OF SUBSEQUENT FAILURES IN',   &
             ' THE STEP ACCEPTANCE TEST(S)', t74, '*')
  WRITE (nunit,90)
  90 FORMAT ('   *       THE PROJECTED DELTA IS LESS',   &
             ' THAN THAT AT THE ALREADY', t74, '*')
  WRITE (nunit,100)
  100 FORMAT ('   *       SUCCESSFUL STEP - RETURN TO',   &
              ' SUCCESSFUL STEP AND ACCEPT', t74, '*')
ELSE IF (retcod == 4) THEN
  WRITE (nunit,110)
  110 FORMAT ('   *       STEP ACCEPTED BY STEP SIZE CRITERION',   &
              ' ONLY - DELTA REDUCED', t74, '*')
ELSE IF (retcod == 5) THEN
  WRITE (nunit,120)
  120 FORMAT ('   *       STEP ACCEPTED - NEW FUNCTION VALUE > PREVIOUS =>', &
              t74, '*')
  WRITE (nunit,130)
  130 FORMAT ('   *       REDUCE TRUST REGION', t74, '*')
ELSE IF (retcod == 6) THEN
  WRITE (nunit,140)
  140 FORMAT ('   *       STEP ACCEPTED - DELTA CHANGED AS DETAILED ABOVE', &
              t74, '*')
ELSE IF (retcod == 7) THEN
  WRITE (nunit,150)
  150 FORMAT ('   *       NO PROGRESS MADE: RELATIVE STEP',   &
              ' SIZE IS TOO SMALL', t74, '*')
  WRITE (nunit,160) rellen,stptol
  160 FORMAT ('   *       REL. STEP SIZE, RELLEN = ', g12.3,', STPTOL = ', &
              g12.3, t74, '*')
ELSE IF (retcod == 8) THEN
  WRITE (nunit,170)
  170 FORMAT ('   *       POINT MODIFIED BY CONSTRAINTS',   &
              ' NOT A DESCENT DIRECTION', t74, '*')
  WRITE (nunit,180)
  180 FORMAT ('   *       DELTA REDUCED TO CONFAC*RATIOM*DELTA', t74, '*' )
ELSE IF (retcod == 9) THEN
  WRITE (nunit,190)
  190 FORMAT ('   *       OVERFLOW DETECTED IN FUNCTION VECTOR',   &
              ' - DELTA REDUCED', t74, '*')
ELSE IF (retcod == 10) THEN
  WRITE (nunit,200)
  200 FORMAT ('   *       OVERFLOW IN OBJECTIVE FUNCTION',   &
              ' - DELTA REDUCED', t74, '*')
ELSE IF (retcod == 11) THEN
  WRITE (nunit,210)
  210 FORMAT ('   *       STEP NOT ACCEPTED - REDUCE TRUST',   &
              ' REGION SIZE BY MINIMIZATION', t74, '*')
  WRITE (nunit,220)
  220 FORMAT ('   *       OF QUADRATIC MODEL IN STEP DIRECTION', t74, '*')
ELSE
  WRITE (nunit,230)
  230 FORMAT ('   *       PROMISING STEP - INCREASE DELTA TO',   &
              ' NEWLEN AND TRY A NEW STEP', t74, '*')
END IF
WRITE (nunit,10)
WRITE (nunit,240) delta
240 FORMAT ('   *       DELTA ON RETURN FROM TRUST REGION',   &
            ' UPDATING: ', g11.3, t74, '*')
RETURN
!
!       LAST CARD OF SUBROUTINE RCDPRT.
!
END SUBROUTINE rcdprt



SUBROUTINE rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, b)
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES, BY BACKWARDS SUBSTITUTION,
!
!           RX=B
!
!           WHERE    R IS TAKEN FROM THE QR DECOMPOSITION AND IS STORED IN
!                      THE STRICT UPPER TRIANGLE OF MATRIX A AND THE VECTOR,
!                      RDIAG
!                    B IS A GIVEN RIGHT HAND SIDE WHICH IS OVERWRITTEN
!
!    FRSTOV  INDICATES FIRST OVERFLOW - USED ONLY TO SET BORDERS FOR OUTPUT
!

LOGICAL, INTENT(IN)        :: overch
LOGICAL, INTENT(OUT)       :: overfl
INTEGER, INTENT(IN)        :: maxexp
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
REAL (dp), INTENT(IN)      :: a(n,n)
REAL (dp), INTENT(IN)      :: rdiag(n)
REAL (dp), INTENT(IN OUT)  :: b(n)

! Local variables

REAL (dp)  :: eps, maxlog, SUM, tmplog
LOGICAL    :: frstov
INTEGER    :: i, j, jstar
!
frstov=.true.
overfl=.false.
eps=ten**(-maxexp)
!
IF (overch) THEN
  IF (LOG10(ABS(b(n))+eps) - LOG10(ABS(rdiag(n))+eps) > maxexp) THEN
    overfl=.true.
    b(n)=SIGN(ten**maxexp,b(n))*SIGN(one,rdiag(n))
    IF (output > 2 .AND. (.NOT.wrnsup)) THEN
      frstov=.false.
      WRITE (nunit,10)
      10 FORMAT ('   *', t74, '*')
      WRITE (nunit,20) n,b(n)
      20 FORMAT ('   *    WARNING: COMPONENT ', i3, ' SET TO ',  &
                 g12.3, t74, '*')
    END IF
    GO TO 30
  END IF
END IF
b(n)=b(n)/rdiag(n)

30 DO  i=n-1,1,-1
  IF (overch) THEN
!
!         CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW.
!
    maxlog=LOG10(ABS(b(i))+eps) - LOG10(ABS(rdiag(i))+eps)
    jstar=0
    DO  j=i+1,n
      tmplog=LOG10(ABS(a(i,j))+eps) + LOG10(ABS(b(j))+eps) -  &
             LOG10(ABS(rdiag(i))+eps)
      IF (tmplog > maxlog) THEN
        jstar=j
        maxlog=tmplog
      END IF
    END DO
!
!            IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
!            TERM WITH CORRECT SIGN.
!
    IF (maxlog > maxexp) THEN
      overfl=.true.
      IF (jstar == 0) THEN
        b(i)=SIGN(ten**maxexp,b(i))*SIGN(one,rdiag(i))
      ELSE
        b(i)=-SIGN(ten**maxexp,a(i,jstar))*SIGN(one,b(jstar))*  &
              SIGN(one,rdiag(i))
      END IF
      IF (frstov) THEN
        frstov=.false.
        WRITE (nunit,10)
      END IF
      IF (output > 2 .AND. (.NOT.wrnsup)) WRITE (nunit,20) i,b(i)
      CYCLE
    END IF
  END IF
!
!         SUM FOR EACH TERM ORDERING OPERATIONS TO MINIMIZE
!         POSSIBILITY OF OVERFLOW.
!
  sum=zero
  DO  j=i+1,n
    sum=sum + (MIN(ABS(a(i,j)),ABS(b(j)))/rdiag(i))*(MAX(ABS(a(i,j))  &
               ,ABS(b(j))))*SIGN(one,a(i,j))*SIGN(one,b(j))
  END DO
  b(i)=b(i)/rdiag(i)-sum
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE RSOLV.
!
END SUBROUTINE rsolv



SUBROUTINE rtrmul (n, a, h, rdiag)

! N.B. Argument WV1 has been removed.

!
!    SEPT. 4, 1991
!
!    FIND R^R FOR QR-DECOMPOSED JACOBIAN.
!
!    R IS STORED IN STRICT UPPER TRIANGLE OF A AND RDIAG.
!

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: h(n,n)
REAL (dp), INTENT(IN)      :: rdiag(n)

! Local variables

INTEGER    :: i
REAL (dp)  :: wv1(n)
!
!       TEMPORARILY REPLACE DIAGONAL OF R IN A (A IS RESTORED LATER).
!
DO  i=1,n
  wv1(i)=a(i,i)
  a(i,i)=rdiag(i)
END DO
CALL utumul (n, n, n, n, n, n, a, h)
DO  i=1,n
  a(i,i)=wv1(i)
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE RTRMUL.
!
END SUBROUTINE rtrmul



SUBROUTINE setup (absnew, cauchy, deuflh, geoms, linesr, newton, overch,  &
                  acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit,  &
                  maxns, maxqns, minqns, n, narmij, niejev, njacch, output, &
                  qnupdm, stopcr, supprs, trupdm, alpha, confac, delta,   &
                  delfac, epsmch, etafac, fdtolj, ftol, lam0, mstpf,   &
                  nsttol, omega, ratiof, sigma, stptol, boundl, boundu,   &
                  scalef, scalex, help)
!
!    DEC. 7, 1991
!
!    SUBROUTINE SETUP ASSIGNS DEFAULT VALUES TO ALL REQUISITE PARAMETERS.
!

LOGICAL, INTENT(OUT)    :: absnew
LOGICAL, INTENT(OUT)    :: cauchy
LOGICAL, INTENT(OUT)    :: deuflh
LOGICAL, INTENT(OUT)    :: geoms
LOGICAL, INTENT(OUT)    :: linesr
LOGICAL, INTENT(OUT)    :: newton
LOGICAL, INTENT(OUT)    :: overch
INTEGER, INTENT(OUT)    :: acptcr
INTEGER, INTENT(OUT)    :: itsclf
INTEGER, INTENT(OUT)    :: itsclx
INTEGER, INTENT(OUT)    :: jactyp
INTEGER, INTENT(OUT)    :: jupdm
INTEGER, INTENT(OUT)    :: maxexp
INTEGER, INTENT(OUT)    :: maxit
INTEGER, INTENT(OUT)    :: maxns
INTEGER, INTENT(OUT)    :: maxqns
INTEGER, INTENT(OUT)    :: minqns
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(OUT)    :: narmij
INTEGER, INTENT(OUT)    :: niejev
INTEGER, INTENT(OUT)    :: njacch
INTEGER, INTENT(OUT)    :: output
INTEGER, INTENT(OUT)    :: qnupdm
INTEGER, INTENT(OUT)    :: stopcr
INTEGER, INTENT(OUT)    :: supprs
INTEGER, INTENT(OUT)    :: trupdm
REAL (dp), INTENT(OUT)  :: alpha
REAL (dp), INTENT(OUT)  :: confac
REAL (dp), INTENT(OUT)  :: delta
REAL (dp), INTENT(OUT)  :: delfac
REAL (dp), INTENT(OUT)  :: epsmch
REAL (dp), INTENT(OUT)  :: etafac
REAL (dp), INTENT(OUT)  :: fdtolj
REAL (dp), INTENT(OUT)  :: ftol
REAL (dp), INTENT(OUT)  :: lam0
REAL (dp), INTENT(OUT)  :: mstpf
REAL (dp), INTENT(OUT)  :: nsttol
REAL (dp), INTENT(OUT)  :: omega
REAL (dp), INTENT(OUT)  :: ratiof
REAL (dp), INTENT(OUT)  :: sigma
REAL (dp), INTENT(OUT)  :: stptol
REAL (dp), INTENT(OUT)  :: boundl(n)
REAL (dp), INTENT(OUT)  :: boundu(n)
REAL (dp), INTENT(OUT)  :: scalef(n)
REAL (dp), INTENT(OUT)  :: scalex(n)
CHARACTER (LEN=6), INTENT(OUT)  :: help

! INTEGER   :: i1mach
! REAL (dp) :: d1mach

! Local variables

REAL (dp)  :: temp, xmax
INTEGER    :: i, ibeta, it, maxebb, minebb
!
!       LOGICAL VALUES.
!
absnew=.false.
bypass=.false.
cauchy=.false.
deuflh=.true.
geoms=.true.
linesr=.true.
matsup=.false.
newton=.false.
overch=.false.
wrnsup=.false.
!
!       INTEGER VALUES.
!
acptcr=12
itsclf=0
itsclx=0
jactyp=1
jupdm=0
maxit=100    ! Reduced from 250 by AJM
maxns=50
maxqns=10
minqns=7
narmij=1
nfetot=0
niejev=1
njacch=1
output=2
qnupdm=1
stopcr=12
supprs=0
trupdm=0
!
!       REAL VALUES.
!
alpha=1.0D-04
confac=0.95D0
delta=-1.0D0
delfac=2.0D0
etafac=0.2D0
lam0=1.0D0
mstpf=1.0D3
omega=0.1D0
ratiof=0.70D0
sigma=0.5D0
!
!       CHARACTER VARIABLE.
!
help(1:4)='NONE'
!
!       NOTE: NOTATIONAL CHANGES IN CALLING PROGRAM FROM MACHAR
!             1)  EPSMCH DENOTES MACHINE EPSILON
!             2)  MINEBB DENOTES MINIMUM EXPONENT BASE BETA
!             3)  MAXEBB DENOTES MAXIMUM EXPONENT BASE BETA
!
it     = DIGITS(zero)
ibeta  = RADIX(zero)
minebb = MINEXPONENT(zero)
maxebb = MAXEXPONENT(zero)
epsmch = EPSILON(zero)
xmax   = HUGE(zero)
maxexp = INT(DBLE(maxebb)*LOG(DBLE(ibeta))/LOG(ten))
!
!       VALUES FOR TWO-NORM CALCULATIONS.
!
smallb=DBLE(ibeta)**((minebb+1)/2)
bigb=DBLE(ibeta)**((maxebb-it+1)/2)
smalls=DBLE(ibeta)**((minebb-1)/2)
bigs=DBLE(ibeta)**((maxebb+it-1)/2)
bigr=xmax
!
!       SET STOPPING CRITERIA PARAMETERS.
!
fdtolj=1.0D-06
ftol=epsmch**0.333
nsttol=ftol*ftol
stptol=nsttol
!
!       VECTOR VALUES.
!
temp=-10.0D0**maxexp
DO  i=1,n
  boundl(i)=temp
  boundu(i)=-temp
  scalef(i)=1.0D0
  scalex(i)=1.0D0
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE SETUP.
!
END SUBROUTINE setup



SUBROUTINE title (cauchy, deuflh, geoms, linesr, newton, overch, acptcr,  &
                  contyp, itsclf, itsclx, jactyp, jupdm, maxit, maxns,  &
                  maxqns, mgll, minqns, n, narmij, ninitn, njacch, nunit,  &
                  output, qnupdm, stopcr, trupdm, alpha, confac, delfac,  &
                  delta, epsmch, etafac, fcnold, ftol, lam0, maxstp, mstpf, &
                  nsttol, omega, ratiof, sigma, stptol, boundl, boundu,   &
                  fvecc, scalef, scalex, xc)
!
!    APR. 13, 1991
!

LOGICAL, INTENT(IN)    :: cauchy
LOGICAL, INTENT(IN)    :: deuflh
LOGICAL, INTENT(IN)    :: geoms
LOGICAL, INTENT(IN)    :: linesr
LOGICAL, INTENT(IN)    :: newton
LOGICAL, INTENT(IN)    :: overch
INTEGER, INTENT(IN)    :: acptcr
INTEGER, INTENT(IN)    :: contyp
INTEGER, INTENT(IN)    :: itsclf
INTEGER, INTENT(IN)    :: itsclx
INTEGER, INTENT(IN)    :: jactyp
INTEGER, INTENT(IN)    :: jupdm
INTEGER, INTENT(IN)    :: maxit
INTEGER, INTENT(IN)    :: maxns
INTEGER, INTENT(IN)    :: maxqns
INTEGER, INTENT(IN)    :: mgll
INTEGER, INTENT(IN)    :: minqns
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN)    :: narmij
INTEGER, INTENT(IN)    :: ninitn
INTEGER, INTENT(IN)    :: njacch
INTEGER, INTENT(IN)    :: nunit
INTEGER, INTENT(IN)    :: output
INTEGER, INTENT(IN)    :: qnupdm
INTEGER, INTENT(IN)    :: stopcr
INTEGER, INTENT(IN)    :: trupdm
REAL (dp), INTENT(IN)  :: alpha
REAL (dp), INTENT(IN)  :: confac
REAL (dp), INTENT(IN)  :: delfac
REAL (dp), INTENT(IN)  :: delta
REAL (dp), INTENT(IN)  :: epsmch
REAL (dp), INTENT(IN)  :: etafac
REAL (dp), INTENT(IN)  :: fcnold
REAL (dp), INTENT(IN)  :: ftol
REAL (dp), INTENT(IN)  :: lam0
REAL (dp), INTENT(IN)  :: maxstp
REAL (dp), INTENT(IN)  :: mstpf
REAL (dp), INTENT(IN)  :: nsttol
REAL (dp), INTENT(IN)  :: omega
REAL (dp), INTENT(IN)  :: ratiof
REAL (dp), INTENT(IN)  :: sigma
REAL (dp), INTENT(IN)  :: stptol
REAL (dp), INTENT(IN)  :: boundl(n)
REAL (dp), INTENT(IN)  :: boundu(n)
REAL (dp), INTENT(IN)  :: fvecc(n)
REAL (dp), INTENT(IN)  :: scalef(n)
REAL (dp), INTENT(IN)  :: scalex(n)
REAL (dp), INTENT(IN)  :: xc(n)

! Local variables

INTEGER  :: i
!
IF (output < 2) RETURN
WRITE (nunit,10)
10 FORMAT (////t3, 72('*'))
WRITE (nunit,20)
20 FORMAT (t3, 72('*'))
WRITE (nunit,20)
WRITE (nunit,30)
30 FORMAT ('   *', t74, '*')
WRITE (nunit,40)
40 FORMAT ('   *', t36, 'NNES', t74, '*')
WRITE (nunit,30)
WRITE (nunit,50)
50 FORMAT ('   *', t12, 'NONMONOTONIC NONLINEAR EQUATION SOLVER',   &
           '  VERSION 1.05', t74, '*')
WRITE (nunit,30)
WRITE (nunit,60)
60 FORMAT ('   *', t24, 'COPYRIGHT 1991, BY R.S. BAIN', t74, '*')
WRITE (nunit,30)
WRITE (nunit,20)
WRITE (nunit,20)
WRITE (nunit,20)
WRITE (nunit,70)
70 FORMAT (///)
IF (output < 3) GO TO 780
WRITE (nunit,80)
80 FORMAT ('1', t3, 72('*'))
WRITE (nunit,20)
WRITE (nunit,30)
IF (newton) THEN
  IF (jupdm == 0) THEN
    WRITE (nunit,90)
    90 FORMAT ('   *  METHOD: NEWTON (NO LINE SEARCH)', t74, '*')
  ELSE IF (jupdm == 1) THEN
    WRITE (nunit,100)
    100 FORMAT ('   *  METHOD: QUASI-NEWTON (NO LINE SEARCH)',   &
                ' USING BROYDEN UPDATE', t74, '*')
  ELSE IF (jupdm == 2) THEN
    WRITE (nunit,110)
    110 FORMAT ('   *  METHOD: QUASI-NEWTON (NO LINE SEARCH)',   &
                ' USING LEE AND LEE UPDATE', t74, '*')
  END IF
  WRITE (nunit,30)
  IF (overch) THEN
    WRITE (nunit,120)
    120 FORMAT ('   *  OVERLOW CHECKING IN USE', t74, '*')
  ELSE
    WRITE (nunit,130)
    130 FORMAT ('   *  OVERFLOW CHECKING NOT IN USE', t74, '*')
  END IF
  IF (jactyp == 0) THEN
    IF (njacch > 0) THEN
      WRITE (nunit,140) njacch
      140 FORMAT ('   *  ANALYTICAL JACOBIAN USED,',    &
                  ' CHECKED NUMERICALLY, NJACCH: ', i5, t74, '*')
    ELSE
      WRITE (nunit,150)
      150 FORMAT ('   *  ANALYTICAL JACOBIAN USED; NOT',   &
                  ' CHECKED', t74, '*')
    END IF
  ELSE IF (jactyp == 1) THEN
    WRITE (nunit,160)
    160 FORMAT ('   *  JACOBIAN ESTIMATED USING FORWARD',   &
                ' DIFFERENCES', t74, '*')
  ELSE IF (jactyp == 2) THEN
    WRITE (nunit,170)
    170 FORMAT ('   *  JACOBIAN ESTIMATED USING BACKWARD',   &
                ' DIFFERENCES', t74, '*')
  ELSE
    WRITE (nunit,180)
    180 FORMAT ('   *  JACOBIAN ESTIMATED USING CENTRAL',   &
                ' DIFFERENCES', t74, '*')
  END IF
  WRITE (nunit,30)
  WRITE (nunit,20)
ELSE
  IF (linesr) THEN
    WRITE (nunit,30)
    IF (deuflh) THEN
      WRITE (nunit,190)
      190 FORMAT ('   *  DEUFLHARD RELAXATION FACTOR ',   &
                  'INITIALIZATION IN EFFECT', t74, '*')
    ELSE
      WRITE (nunit,200)
      200 FORMAT ('   *  DEUFLHARD RELAXATION FACTOR ',   &
                  'INITIALIZATION NOT IN EFFECT', t74, '*')
    END IF
  ELSE
    IF (etafac == one) THEN
      WRITE (nunit,210)
      210 FORMAT ('   *  METHOD: TRUST REGION USING',   &
                  ' SINGLE DOGLEG STEPS', t74, '*')
    ELSE
      WRITE (nunit,220)
      220 FORMAT ('   *  METHOD: TRUST REGION USING',   &
                  ' DOUBLE DOGLEG STEPS', t74, '*')
    END IF
    WRITE (nunit,30)
    IF (cauchy) THEN
      WRITE (nunit,230)
      230 FORMAT ('   *  INITIAL STEP CONSTRAINED BY',   &
                  ' SCALED CAUCHY STEP', t74, '*')
    ELSE
      WRITE (nunit,240)
      240 FORMAT ('   *  INITIAL STEP CONSTRAINED BY',   &
                  ' SCALED NEWTON STEP', t74, '*')
    END IF
  END IF
  IF (geoms) THEN
    WRITE (nunit,250)
    250 FORMAT ('   *  METHOD: GEOMETRIC SEARCH', t74, '*')
  ELSE
    WRITE (nunit,260)
    260 FORMAT ('   *  METHOD: SEARCH BASED ON',   &
                ' SUCCESSIVE MINIMIZATIONS', t74, '*')
  END IF
  IF (overch) THEN
    WRITE (nunit,120)
  ELSE
    WRITE (nunit,130)
  END IF
  IF (jupdm == 0) THEN
    WRITE (nunit,270)
    270 FORMAT ('   *  NO QUASI-NEWTON UPDATE USED', t74, '*')
  END IF
  IF (jupdm == 1) THEN
    IF (qnupdm == 0) THEN
      WRITE (nunit,280)
      280 FORMAT ('   *  BROYDEN QUASI-NEWTON UPDATE',   &
                  ' OF UNFACTORED JACOBIAN', t74, '*')
    ELSE
      WRITE (nunit,290)
      290 FORMAT ('   *  BROYDEN QUASI-NEWTON UPDATE',   &
                  ' OF FACTORED JACOBIAN', t74, '*')
    END IF
  ELSE IF (jupdm == 2) THEN
    IF (qnupdm == 0) THEN
      WRITE (nunit,300)
      300 FORMAT ('   *  LEE AND LEE QUASI-NEWTON UPDATE',   &
                  ' OF UNFACTORED JACOBIAN', t74, '*')
    ELSE
      WRITE (nunit,310)
      310 FORMAT ('   *  LEE AND LEE QUASI-NEWTON UPDATE',   &
                  ' OF FACTORED JACOBIAN', t74, '*')
    END IF
  END IF
  IF (jactyp == 0) THEN
    IF (njacch > 0) THEN
      WRITE (nunit,140)
    ELSE
      WRITE (nunit,150)
    END IF
  ELSE IF (jactyp == 1) THEN
    WRITE (nunit,160)
  ELSE IF (jactyp == 2) THEN
    WRITE (nunit,170)
  ELSE
    WRITE (nunit,180)
  END IF
  IF (.NOT.linesr) THEN
    WRITE (nunit,30)
    IF (trupdm == 0 .AND. jupdm > 0) THEN
      WRITE (nunit,320)
      320 FORMAT ('   *  TRUST REGION UPDATED USING',   &
                  ' POWELL STRATEGY', t74, '*')
    ELSE
      WRITE (nunit,330)
      330 FORMAT ('   *  TRUST REGION UPDATED USING',   &
                  ' DENNIS AND SCHNABEL STRATEGY', t74, '*')
    END IF
  END IF
  WRITE (nunit,30)
  WRITE (nunit,20)
  WRITE (nunit,30)
  IF (itsclf /= 0) THEN
    WRITE (nunit,340) itsclf
    340 FORMAT ('   *  ADAPTIVE FUNCTION SCALING STARTED AT',   &
                ' ITERATION: ..........', i6, t74, '*')
    WRITE (nunit,30)
  END IF
  IF (itsclx /= 0) THEN
    WRITE (nunit,350) itsclx
    350 FORMAT ('   *  ADAPTIVE VARIABLE SCALING STARTED AT',   &
                ' ITERATION: ..........', i6, t74, '*')
    WRITE (nunit,30)
  END IF
  IF (linesr) THEN
    IF (jupdm == 0) THEN
      WRITE (nunit,360) maxns
      360 FORMAT ('   *  MAXIMUM NUMBER OF STEPS IN LINE',   &
                  ' SEARCH, MAXNS: ...........', i6, t74, '*')
    ELSE
      WRITE (nunit,370) maxns
      370 FORMAT ('   *  MAXIMUM NUMBER OF NEWTON LINE',   &
                  ' SEARCH STEPS, MAXNS: .......', i6, t74, '*')
      WRITE (nunit,380) maxqns
      380 FORMAT ('   *  MAXIMUM NUMBER OF QUASI-NEWTON',   &
                  ' LINE SEARCH STEPS, MAXQNS: ', i6, t74, '*')
    END IF
  ELSE
    IF (jupdm == 0) THEN
      WRITE (nunit,390) maxns
      390 FORMAT ('   *  MAXIMUM NUMBER OF TRUST REGION',   &
                  ' UPDATES, MAXNS: ...........', i6, t74, '*')
    ELSE
      WRITE (nunit,400) maxns
      400 FORMAT ('   *  MAXIMUM NO. OF NEWTON TRUST',   &
                  ' REGION UPDATES, MAXNS: .......', i6, t74, '*')
      WRITE (nunit,410) maxqns
      410 FORMAT ('   *  MAXIMUM NO. OF QUASI-NEWTON',   &
                  ' TRUST REGION UPDATES, MAXQNS: ', i6, t74, '*')
    END IF
  END IF
  IF (narmij < maxit) THEN
    WRITE (nunit,30)
    WRITE (nunit,420) mgll
    420 FORMAT ('   *  NUMBER OF OBJECTIVE FUNCTION',   &
                ' VALUES COMPARED, MGLL: ......', i6, t74, '*')
  END IF
  IF (jupdm > 0) THEN
    IF (narmij == maxit) WRITE (nunit,30)
    WRITE (nunit,430) minqns
    430 FORMAT ('   *  MINIMUM NUMBER OF STEPS BETWEEN',   &
                ' JACOBIAN UPDATES, MINQNS: ', i6, t74, '*')
    WRITE (nunit,440) ninitn
    440 FORMAT ('   *  NUMBER OF NON-QUASI-NEWTON',   &
                ' STEPS AT START, NINITN: .......', i6, t74, '*')
  END IF
  WRITE (nunit,450) narmij
  450 FORMAT ('   *  NUMBER OF ARMIJO STEPS AT START, NARMIJ:',  &
              ' .................', i6, t74, '*')
END IF
WRITE (nunit,30)
IF (stopcr == 3) THEN
  WRITE (nunit,460) stopcr
  460 FORMAT ('   *  FUNCTION AND STEP SIZE STOPPING CRITERIA, ',  &
              'STOPCR: ........', i6, t74, '*')
ELSE IF (stopcr == 12) THEN
  WRITE (nunit,470) stopcr
  470 FORMAT ('   *  FUNCTION OR STEP SIZE STOPPING',   &
              ' CRITERIA, STOPCR: .........', i6, t74, '*')
ELSE IF (stopcr == 1) THEN
  WRITE (nunit,480) stopcr
  480 FORMAT ('   *  STEP SIZE STOPPING CRITERION,',   &
              ' STOPCR: ....................', i6, t74, '*')
ELSE
  WRITE (nunit,490) stopcr
  490 FORMAT ('   *  FUNCTION STOPPING CRITERION,',   &
              ' STOPCR: .....................', i6, t74, '*')
END IF
IF (.NOT.newton) THEN
  WRITE (nunit,30)
  IF (acptcr == 12) THEN
    WRITE (nunit,500) acptcr
    500 FORMAT ('   *  FUNCTION AND STEP SIZE ACCEPTANCE CRITERIA,',  &
                ' ACPTCR: ......', i6, t74, '*')
  ELSE IF (acptcr == 2) THEN
    WRITE (nunit,510) acptcr
    510 FORMAT ('   *  STEP SIZE ACCEPTANCE CRITERION, ',   &
                'ACPTCR: ..................', i6, t74, '*')
  ELSE
    WRITE (nunit,520) acptcr
    520 FORMAT ('   *  FUNCTION ACCEPTANCE CRITERION, ',   &
                'ACPTCR: ...................', i6, t74, '*')
  END IF
  IF (contyp /= 0) THEN
    WRITE (nunit,30)
    WRITE (nunit,530) contyp
    530 FORMAT ('   *  CONSTRAINTS IN USE, CONTYP: ',   &
                '..............................', i6, t74, '*')
  END IF
END IF
WRITE (nunit,30)
WRITE (nunit,20)
WRITE (nunit,30)
WRITE (nunit,540) epsmch
540 FORMAT ('   *  ESTIMATED MACHINE EPSILON, EPSMCH:',   &
            ' ...................', g10.3, t74, '*')
WRITE (nunit,30)
WRITE (nunit,550) mstpf
550 FORMAT ('   *  FACTOR TO ESTABLISH MAXIMUM STEP SIZE',   &
            ', MSTPF : ........', g10.3, t74, '*')
WRITE (nunit,560) maxstp
560 FORMAT ('   *  CALCULATED MAXIMUM STEP SIZE, MAXSTP:',  &
            ' ................', g10.3, t74, '*')
IF (.NOT.linesr) THEN
  IF (delta < zero) THEN
    WRITE (nunit,30)
    WRITE (nunit,570)
    570 FORMAT ('   *  INITIAL TRUST REGION NOT PROVIDED', t74, '*')
  ELSE
    WRITE (nunit,30)
    WRITE (nunit,580) delta
    580 FORMAT ('   *  INITIAL TRUST REGION SIZE, DELTA:',   &
                ' ....................', g10.3, t74, '*')
  END IF
  IF (etafac < one) THEN
    WRITE (nunit,590) etafac
    590 FORMAT ('   *  FACTOR TO SET DIRECTION OF TRUST REGION ',  &
                'STEP, ETAFAC: ... ', f6.4, t74, '*')
  END IF
  WRITE (nunit,600) delfac
  600 FORMAT ('   *  TRUST REGION UPDATING FACTOR, DELFAC: ',  &
              '................', g10.3, t74, '*')
END IF
IF (.NOT.newton) THEN
  WRITE (nunit,30)
  WRITE (nunit,610) alpha
  610 FORMAT ('   *  FACTOR IN OBJECTIVE FUNCTION',  &
              ' COMPARISON, ALPHA: ......', g10.3, t74, '*')
  IF (linesr .AND. (.NOT.newton)) THEN
    WRITE (nunit,30)
    WRITE (nunit,620) sigma
    620 FORMAT ('   *  REDUCTION FACTOR FOR RELAXATION FACTOR, ',  &
                'SIGMA: .......', g10.3, t74, '*')
  END IF
  IF (jupdm /= 0) THEN
    WRITE (nunit,30)
    WRITE (nunit,630) ratiof
    630 FORMAT ('   *  REDUCTION REQUIRED IN OBJ. FUNCTION',  &
                ' FOR QN STEP, RATIOF:  ', f6.4, t74, '*')
  END IF
  IF (jupdm == 2) THEN
    WRITE (nunit,30)
    WRITE (nunit,640) omega
    640 FORMAT ('   *  FACTOR IN LEE AND LEE UPDATE, OMEGA:',  &
                ' .....................', f6.4, t74, '*')
  END IF
END IF
WRITE (nunit,30)
IF (stopcr /= 2) THEN
  WRITE (nunit,650) stptol
  650 FORMAT ('   *  STOPPING TOLERANCE FOR STEP SIZE, STPTOL: ',  &
              '............', g10.3, t74, '*')
  WRITE (nunit,660) nsttol
  660 FORMAT ('   *  STOPPING TOLERANCE FOR NEWTON STEP, NSTTOL: ',  &
              '..........', g10.3, t74, '*')
END IF
IF (stopcr /= 1) THEN
  WRITE (nunit,670) ftol
  670 FORMAT ('   *  STOPPING TOLERANCE FOR OBJECTIVE FUNCTION, ',  &
              'FTOL: .....', g10.3, t74, '*')
END IF
IF (linesr .AND. (.NOT.newton) .AND. lam0 < one) THEN
  WRITE (nunit,30)
  WRITE (nunit,680) lam0
  680 FORMAT ('   *  INITIAL LAMBDA IN LINE SEARCH,',   &
              ' LAM0: .................', g10.3, t74, '*')
END IF
IF (contyp > 0) THEN
  WRITE (nunit,30)
  WRITE (nunit,690) confac
  690 FORMAT ('   *  FACTOR TO ENSURE STEP WITHIN CONSTRAINTS, ',   &
              'CONFAC: ........', f6.4, t74, '*')
END IF
WRITE (nunit,30)
WRITE (nunit,20)
WRITE (nunit,30)
WRITE (nunit,700)
700 FORMAT ('   *  SCALING FACTORS', t74, '*')
WRITE (nunit,30)
WRITE (nunit,710)
710 FORMAT ('   *      COMPONENT VALUES', t50, 'FUNCTION VALUES', t74, '*' )
WRITE (nunit,30)
DO  i=1,n
  WRITE (nunit,720) i,scalex(i),i,scalef(i)
  720 FORMAT ('   *  SCALEX(', i3, ') = ', g10.3, t45, 'SCALEF(', i3,  &
              ') = ', g10.3, t74, '*')
END DO
IF (contyp > 0) THEN
  WRITE (nunit,30)
  WRITE (nunit,20)
  WRITE (nunit,30)
  WRITE (nunit,740)
  740 FORMAT ('   *  LOWER AND UPPER BOUNDS', t74, '*')
  WRITE (nunit,30)
  WRITE (nunit,750)
  750 FORMAT ('   *', t11, 'LOWER BOUNDS', t49, 'UPPER BOUNDS', t74, '*')
  WRITE (nunit,30)
  DO  i=1,n
    WRITE (nunit,760) i,boundl(i),i,boundu(i)
    760 FORMAT ('   *  BOUNDL(', i3, ') = ', g10.3, t45, 'BOUNDU(',  &
                i3, ') = ', g10.3, t74, '*')
  END DO
END IF
WRITE (nunit,30)

780 IF (output == 2) WRITE (nunit,20)
WRITE (nunit,20)
WRITE (nunit,30)
WRITE (nunit,790)
790 FORMAT ('   *    INITIAL ESTIMATES', t41, 'INITIAL FUNCTION VALUES',  &
            t74, '*')
WRITE (nunit,30)
DO  i=1,n
  WRITE (nunit,800) i,xc(i),i,fvecc(i)
  800 FORMAT ('   *  X(', i3, ') = ', g12.3, t42, 'F(', i3, ') = ',  &
              g12.3, t74, '*')
END DO
WRITE (nunit,30)
WRITE (nunit,820) fcnold
820 FORMAT ('   *  INITIAL OBJECTIVE FUNCTION VALUE = ', g10.3, t74, '*')
WRITE (nunit,30)
WRITE (nunit,20)
WRITE (nunit,20)
IF (output < 3) RETURN
WRITE (nunit,830)
830 FORMAT (//t27, 23('*'))
WRITE (nunit,840)
840 FORMAT (t27, 23('*'))
WRITE (nunit,850)
850 FORMAT (t27, '*', t49, '*')
WRITE (nunit,860)
860 FORMAT (t27, '*  UPDATED ESTIMATES  *')
WRITE (nunit,850)
WRITE (nunit,840)
WRITE (nunit,840)
WRITE (nunit,870)
870 FORMAT (//)
RETURN
!
!       LAST CARD OF SUBROUTINE TITLE.
!
END SUBROUTINE title



SUBROUTINE trstup (geoms, newtkn, overch, overfl, qrsing, sclfch, sclxch,  &
                   acpcod, acpstr, acptcr, contyp, isejac, jupdm, maxexp,  &
                   mgll, mnew, n, narmij, nfunc, notrst, nunit, output,   &
                   qnupdm, retcod, trupdm, alpha, confac, delfac, delstr,  &
                   delta, epsmch, fcnmax, fcnnew, fcnold, fcnpre, maxstp,  &
                   newlen, newmax, powtau, rellen, stptol, a, astore,   &
                   boundl, boundu, delf, fplpre, ftrack, fvec, fvecc, hhpi, &
                   jac, rdiag, rhs, s, sbar, scalef, scalex, strack,   &
                   xc, xplpre, xplus, fvecev)

! N.B. Argument WV3 has been removed.

!
!    FEB. 28, 1992
!
!    THIS SUBROUTINE CHECKS FOR ACCEPTANCE OF A TRUST REGION
!    STEP GENERATED BY THE DOUBLE DOGLEG METHOD.
!

LOGICAL, INTENT(IN)        :: geoms
LOGICAL, INTENT(IN OUT)    :: newtkn
LOGICAL, INTENT(IN OUT)    :: overch
LOGICAL, INTENT(OUT)       :: overfl
LOGICAL, INTENT(IN)        :: qrsing
LOGICAL, INTENT(IN OUT)    :: sclfch
LOGICAL, INTENT(IN OUT)    :: sclxch
INTEGER, INTENT(IN OUT)    :: acpcod
INTEGER, INTENT(IN OUT)    :: acpstr
INTEGER, INTENT(IN)        :: acptcr
INTEGER, INTENT(IN)        :: contyp
INTEGER, INTENT(IN)        :: isejac
INTEGER, INTENT(IN)        :: jupdm
INTEGER, INTENT(IN OUT)    :: maxexp
INTEGER, INTENT(IN)        :: mgll
INTEGER, INTENT(IN)        :: mnew
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: narmij
INTEGER, INTENT(IN OUT)    :: nfunc
INTEGER, INTENT(IN)        :: notrst
INTEGER, INTENT(IN)        :: nunit
INTEGER, INTENT(IN)        :: output
INTEGER, INTENT(IN)        :: qnupdm
INTEGER, INTENT(OUT)       :: retcod
INTEGER, INTENT(IN OUT)    :: trupdm
REAL (dp), INTENT(IN)      :: alpha
REAL (dp), INTENT(IN)      :: confac
REAL (dp), INTENT(IN)      :: delfac
REAL (dp), INTENT(IN OUT)  :: delstr
REAL (dp), INTENT(IN OUT)  :: delta
REAL (dp), INTENT(IN)      :: epsmch
REAL (dp), INTENT(IN OUT)  :: fcnmax
REAL (dp), INTENT(OUT)     :: fcnnew
REAL (dp), INTENT(IN)      :: fcnold
REAL (dp), INTENT(IN OUT)  :: fcnpre
REAL (dp), INTENT(IN OUT)  :: maxstp
REAL (dp), INTENT(IN)      :: newlen
REAL (dp), INTENT(IN OUT)  :: newmax
REAL (dp), INTENT(IN OUT)  :: powtau
REAL (dp), INTENT(OUT)     :: rellen
REAL (dp), INTENT(IN)      :: stptol
REAL (dp), INTENT(IN)      :: a(n,n)
REAL (dp), INTENT(IN OUT)  :: astore(n,n)
REAL (dp), INTENT(IN)      :: boundl(n)
REAL (dp), INTENT(IN OUT)  :: boundu(n)
REAL (dp), INTENT(IN OUT)  :: delf(n)
REAL (dp), INTENT(IN OUT)  :: fplpre(n)
REAL (dp), INTENT(IN)      :: ftrack(0:mgll-1)
REAL (dp), INTENT(OUT)     :: fvec(n)
REAL (dp), INTENT(IN)      :: fvecc(n)
REAL (dp), INTENT(IN OUT)  :: hhpi(n)
REAL (dp), INTENT(IN)      :: jac(n,n)
REAL (dp), INTENT(IN)      :: rdiag(n)
REAL (dp), INTENT(OUT)     :: rhs(n)
REAL (dp), INTENT(IN OUT)  :: s(n)
REAL (dp), INTENT(OUT)     :: sbar(n)
REAL (dp), INTENT(IN)      :: scalef(n)
REAL (dp), INTENT(IN)      :: scalex(n)
REAL (dp), INTENT(IN)      :: strack(0:mgll-1)
REAL (dp), INTENT(IN)      :: xc(n)
REAL (dp), INTENT(IN OUT)  :: xplpre(n)
REAL (dp), INTENT(OUT)     :: xplus(n)

! EXTERNAL fvecev

INTERFACE
  SUBROUTINE fvecev (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fvecev
END INTERFACE

! Local variables

LOGICAL    :: convio
INTEGER    :: i, j, k
REAL (dp)  :: delfpr, delfts, deltaf, dlftsm, dmult, powlam, powmu,  &
              ratio, ratiom, sbrnrm, sp, ss, stplen, SUM, wv3(n)

REAL (dp), PARAMETER  :: pt5 = 0.5_dp, threeq = 0.75_dp, onept1 = 1.1_dp
!
!       NOTE: ACCEPTANCE CODE, ACPCOD, IS 0 ON ENTRANCE TO TRSTUP
!
convio=.false.
overfl=.false.
retcod = 0    ! Added by AJM
!
IF (output > 3) THEN
  WRITE (nunit,10)
  10 FORMAT ('   *', t74, '*')
  WRITE (nunit,10)
  IF (.NOT.sclfch .AND. (.NOT.sclxch)) THEN
    WRITE (nunit,20)
    20 FORMAT ('   *    TRUST REGION UPDATING', t74, '*')
  ELSE
    WRITE (nunit,30)
    30 FORMAT ('   *    TRUST REGION UPDATING (ALL X''S',  &
               ' AND F''S IN UNSCALED UNITS)', t74, '*')
  END IF
  WRITE (nunit,10)
END IF
!
!    CHECK TO MAKE SURE "S" IS A DESCENT DIRECTION - FIND DIRECTIONAL
!    DERIVATIVE AT CURRENT XC USING S GENERATED BY DOGLEG SUBROUTINE.
!
CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output, delfts, delf, s)
IF (output > 3) THEN
  WRITE (nunit,10)
  WRITE (nunit,40) delfts
  40 FORMAT ('   *       INNER PRODUCT OF DELF AND S, DELFTS: ',  &
             '........', g13.4, t74, '*')
END IF
IF (delfts > zero) THEN
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,50)
    50 FORMAT ('   *       DIRECTIONAL DERIVATIVE POSITIVE',  &
               '; SEARCH DIRECTION REVERSED', t74, '*')
  END IF
  s(1:n) = -s(1:n)
END IF
!
!       FIND MAXIMUM OBJECTIVE FUNCTION VALUE AND MAXIMIUM STEP
!       LENGTH FOR NONMONOTONIC SEARCH.  THIS HAS TO BE DONE ONLY
!       ONCE DURING EACH ITERATION (WHERE NOTRST=1).
!
IF (notrst == 1) THEN
  newmax=newlen
  fcnmax=fcnold
  IF (isejac > narmij) THEN
    IF (isejac < narmij+mgll) THEN
      DO  j=1,mnew
        fcnmax=MAX(fcnmax, ftrack(j-1))
        newmax=MAX(newmax, strack(j-1))
      END DO
    ELSE
      DO  j=0,mnew
        fcnmax=MAX(fcnmax, ftrack(j))
        newmax=MAX(newmax, strack(j))
      END DO
    END IF
  END IF
END IF
!
!       TEST TRIAL POINT - FIND XPLUS AND TEST FOR CONSTRAINT
!       VIOLATIONS IF CONTYP DOES NOT EQUAL 0.
!
DO  i=1,n
  wv3(i)=-one
!
!          WV3 IS A MARKER FOR "VIOLATORS" - IT CHANGES TO 1 OR 2.
!
  xplus(i)=xc(i) + s(i)
  IF (contyp > 0) THEN
    IF (xplus(i) < boundl(i)) THEN
      convio=.true.
      wv3(i)=one
    ELSE IF (xplus(i) > boundu(i)) THEN
      convio=.true.
      wv3(i)=two
    END IF
  END IF
END DO
!
!       IF CONSTRAINT IS VIOLATED ...
!
IF (convio) THEN
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,100)
    100 FORMAT ('   *       CONSTRAINT VIOLATED', t74, '*'/  &
                '   *', t74, '*'/  &
                '   *          TRIAL ESTIMATES (VIOLATIONS MARKED)', t74, '*')
    WRITE (nunit,10)
    DO  i=1,n
      IF (wv3(i) > zero) THEN
        WRITE (nunit,110) i,xplus(i)
        110 FORMAT ('   *', t17, 'XPLUS(', i3, ') = ', g12.3, '  *', t74, '*')
      ELSE
        WRITE (nunit,120) i,xplus(i)
        120 FORMAT ('   *', t17, 'XPLUS(', I3, ') = ', g12.3, t74, '*')
      END IF
    END DO
  END IF
!
!       FIND STEP WITHIN CONSTRAINED REGION.
!
!       FIND THE RATIO OF THE DISTANCE FROM THE (I)TH COMPONENT TO ITS
!       CONSTRAINT TO THE LENGTH OF THE PROPOSED STEP, XPLUS(I)-XC(I).
!       MULTIPLY THIS BY CONFAC (DEFAULT 0.95) TO ENSURE THE NEW STEP STAYS
!       WITHIN THE ACCEPTABLE REGION UNLESS XC IS CLOSE TO THE BOUNDARY
!       (RATIO <= 1/2).   IN SUCH CASES A FACTOR OF 0.5*CONFAC IS USED.
!
!       NOTE: ONLY THE VIOLATING COMPONENTS ARE REDUCED.
!
  ratiom=one
!
!          RATIOM STORES THE MINIMUM VALUE OF RATIO.
!
  DO  i=1,n
    IF (wv3(i) == one) THEN
      ratio=(boundl(i)-xc(i))/s(i)
    ELSE IF (wv3(i) == two) THEN
      ratio=(boundu(i)-xc(i))/s(i)
    END IF
    IF (wv3(i) > zero) THEN
!
!                NOTE: RATIO IS STORED IN WV3 FOR OUTPUT ONLY.
!
      wv3(i)=ratio
!
      ratiom=MIN(ratiom,ratio)
      IF (ratio > pt5) THEN
        xplus(i)=xc(i) + confac*ratio*s(i)
      ELSE
!
!                   WITHIN BUFFER ZONE.
!
        xplus(i)=xc(i) + confac*ratio*s(i)/two
      END IF
      s(i)=xplus(i) - xc(i)
    END IF
  END DO
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,10)
    WRITE (nunit,150)
    150 FORMAT ('   *       NEW S AND XPLUS VECTORS',  &
                ' (WITH RATIOS FOR VIOLATIONS)', t74, '*')
    WRITE (nunit,10)
    WRITE (nunit,160)
    160 FORMAT ('   *       NOTE: RATIOS ARE RATIO OF',  &
                ' LENGTH TO BOUNDARY FROM CURRENT', t74, '*')
    WRITE (nunit,170)
    170 FORMAT ('   *       X VECTOR TO MAGNITUDE OF',  &
                ' CORRESPONDING PROPOSED STEP', t74, '*')
    WRITE (nunit,10)
    DO  i=1,n
      IF (wv3(i) < zero) THEN
        WRITE (nunit,180) i,s(i),i,xplus(i)
        180 FORMAT ('   *       S(', i3, ') = ', g12.3, '    XPLUS(', i3,  &
                    ') = ', g12.3, t74, '*')
      ELSE
        WRITE (nunit,190) i,s(i),i,xplus(i),wv3(i)
        190 FORMAT ('   *       S(', i3, ') = ', g12.3, '    XPLUS(', i3,  &
                    ') = ', g12.3, ' ', g11.3, t74, '*')
      END IF
    END DO
    WRITE (nunit,10)
    WRITE (nunit,210) ratiom
    210 FORMAT ('   *       MINIMUM OF RATIOS, RATIOM: ', g12.3, t74, '* ')
  END IF
!
!       THE NEW POINT, XPLUS, IS NOT NECESSARILY IN A DESCENT DIRECTION.
!       CHECK DIRECTIONAL DERIVATIVE FOR MODIFIED STEP, DLFTSM.
!
  CALL innerp (overch, overfl, maxexp, n, n, n, nunit, output, dlftsm, delf, s)
  IF (output > 3) THEN
    WRITE (nunit,10)
    WRITE (nunit,220) dlftsm
    220 FORMAT ('   *       INNER PRODUCT OF DELF AND MODIFIED S',  &
                ', DL FTSM: ', g12.3, t74, '*')
  END IF
!
!       IF DLFTSM IS POSITIVE REDUCE TRUST REGION.  IF NOT, TEST NEW POINT.
!
  IF (dlftsm > zero) THEN
    delta=confac*ratiom*delta
    retcod=8
    RETURN
  END IF
END IF
!
!       CONSTRAINTS NOT (OR NO LONGER) VIOLATED - TEST NEW POINT.
!
CALL fvecev (overfl, n, fvec, xplus)
nfunc=nfunc + 1
!
!       IF OVERFLOW AT NEW POINT REDUCE TRUST REGION AND RETURN.
!
IF (overfl) THEN
!
!          IF THE OVERFLOW COMES AS A RESULT OF INCREASING DELTA WITHIN THE
!          CURRENT ITERATION (IMPLYING DELSTR IS POSITIVE) AND DIVIDING DELTA
!          BY 10 WOULD PRODUCE A DELTA WHICH IS SMALLER THAN THAT AT THE
!          STORED POINT, THEN USE STORED POINT AS THE UPDATED RESULT.
!
  IF (delstr > delta/ten) THEN
    CALL matcop (n, n, 1, 1, n, 1, xplpre, xplus)
    CALL matcop (n, n, 1, 1, n, 1, fplpre, fvec)
    acpcod=acpstr
    delta=delstr
    fcnnew=fcnpre
    retcod=1
  ELSE
    delta=delta/ten
    retcod=9
  END IF
  RETURN
END IF
!
!       NO OVERFLOW IN RESIDUAL VECTOR.
!
IF (output > 3) THEN
  WRITE (nunit,10)
  WRITE (nunit,230)
  230 FORMAT ('   *', t15, 'TRIAL ESTIMATES', t47, 'FUNCTION VALUES', t74, '*')
  WRITE (nunit,10)
  DO  i=1,n
    WRITE (nunit,240) i,xplus(i),i,fvec(i)
    240 FORMAT ('   *       XPLUS(', i3, ') = ', g12.3, '         FVEC(',  &
                i3, ') = ', g12.3, t74, '*')
  END DO
END IF
!
!       IF NO OVERFLOW WITHIN RESIDUAL VECTOR FIND OBJECTIVE FUNCTION.
!
CALL fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef)
!
!       IF OVERFLOW IN OBJECTIVE FUNCTION EVALUATION REDUCE
!       TRUST REGION AND RETURN.
!
IF (overfl) THEN
!
!       IF THE OVERFLOW COMES AS A RESULT OF INCREASING DELTA WITHIN THE
!       CURRENT ITERATION (SO THAT DELSTR IS POSITIVE) AND DIVIDING DELTA
!       BY 10 WOULD PRODUCE A DELTA WHICH IS SMALLER THAN THAT AT THE STORED
!       POINT THEN USE STORED POINT AS THE UPDATED RESULT.
!
  IF (delstr > delta/ten) THEN
    CALL matcop (n, n, 1, 1, n, 1, xplpre, xplus)
    CALL matcop (n, n, 1, 1, n, 1, fplpre, fvec)
    acpcod=acpstr
    delta=delstr
    fcnnew=fcnpre
    retcod=2
  ELSE
    delta=delta/ten
    retcod=10
  END IF
  RETURN
ELSE
!
!       NO OVERFLOW AT TRIAL POINT - COMPARE OBJECTIVE FUNCTION TO FCNMAX.
!
  IF (output > 3) THEN
    WRITE (nunit,10)
    IF (.NOT.sclfch) THEN
      WRITE (nunit,260) fcnnew
      260 FORMAT ('   *       OBJECTIVE FUNCTION AT XPLUS,',  &
                  ' FCNNEW: .........', g12.4, t74, '*')
    ELSE
      WRITE (nunit,270) fcnnew
      270 FORMAT ('   *       SCALED OBJECTIVE FUNCTION AT XPLUS',  &
                  ', FCNNEW: ..', g12.4, t74, '*')
    END IF
    WRITE (nunit,280) fcnmax + alpha*delfts
    280 FORMAT ('   *       COMPARE TO FCNMAX+ALPHA*DELFTS: ', 14('.'),  &
                g12.4, t74, '*')
  END IF
END IF
!
!       IF ACPTCR=12 CHECK SECOND DEUFLHARD STEP ACCEPTANCE TEST
!       BY FINDING 2-NORM OF SBAR.  THERE ARE FOUR POSSIBILITIES
!       DEPENDING ON WHETHER THE JACOBIAN IS OR IS NOT SINGULAR
!       AND WHETHER QNUPDM IS 0 OR 1.
!
IF (acptcr == 12) THEN
  IF (qrsing) THEN
!
!             FORM -J^F AS RIGHT HAND SIDE - METHOD DEPENDS ON
!             WHETHER QNUPDM EQUALS 0 OR 1 (UNFACTORED OR FACTORED).
!
    IF (qnupdm == 0) THEN
!
!                UNSCALED JACOBIAN IS IN MATRIX JAC.
!
      wv3(1:n) = -fvec(1:n)*scalef(1:n)*scalef(1:n)
      CALL atbmul (n, n, 1, 1, n, n, jac, wv3, rhs)
!
    ELSE
!
!                R IN UPPER TRIANGLE OF A PLUS RDIAG AND Q^ IN JAC
!                FROM QR DECOMPOSITION OF SCALED JACOBIAN.
!
      DO  i=1,n
        sum=zero
        DO  j=1,n
          sum=sum - jac(i,j)*fvec(j)*scalef(j)
        END DO
        wv3(i)=sum
      END DO
      rhs(1)=rdiag(1)*wv3(1)
      DO  j=2,n
        sum=DOT_PRODUCT( a(1:j-1,j), wv3(1:j-1) )
        rhs(j)=sum + rdiag(j)*wv3(j)
      END DO
    END IF
    CALL chsolv (overch, overfl, maxexp, n, nunit, output, a, rhs, sbar)
  ELSE
!
!             RIGHT HAND SIDE IS -FVEC.
!
    IF (qnupdm == 0 .OR. jupdm == 0) THEN
!
!                QR DECOMPOSITION OF SCALED JACOBIAN STORED IN ASTORE.
!
      sbar(1:n) = -fvec(1:n)*scalef(1:n)
      CALL qrsolv (overch, overfl, maxexp, n, nunit, output,  &
                   astore, hhpi, rdiag, sbar)
    ELSE
!
!                SET UP RIGHT HAND SIDE - MULTIPLY -FVEC BY Q^
!                (STORED IN JAC).  RHS IS A WORK VECTOR ONLY HERE.
!
      wv3(1:n) = -fvec(1:n)*scalef(1:n)
      CALL avmul (n, n, n, n, jac, wv3, sbar)
      CALL rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sbar)
    END IF
  END IF
!
!          NORM OF (SCALED) SBAR IS NEEDED FOR SECOND ACCEPTANCE TEST.
!
  wv3(1:n)=scalex(1:n)*sbar(1:n)
  CALL twonrm (overfl, maxexp, n, epsmch, sbrnrm, wv3)
!
  IF (output > 4) THEN
    WRITE (nunit,10)
    IF (.NOT.sclxch) THEN
      WRITE (nunit,370)
      370 FORMAT ('   *          DEUFLHARD SBAR VECTOR', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,380) i,sbar(i)
        380 FORMAT ('   *          SBAR(', i3, ') = ', g12.3, t74, '*')
      END DO
    ELSE
      WRITE (nunit,400)
      400 FORMAT ('   *          DEUFLHARD SBAR VECTOR              ',  &
                  'IN SCALED X UNITS', t74, '*')
      WRITE (nunit,10)
      DO  i=1,n
        WRITE (nunit,410) i,sbar(i),i,scalex(i)*sbar(i)
        410 FORMAT ('   *          SBAR(', I3, ') = ', g12.3, '        SBAR(', &
                    i3, ') = ', g12.3, t74, '*')
      END DO
    END IF
  END IF
  IF (output > 3) THEN
    WRITE (nunit,10)
    IF (.NOT.sclxch) THEN
      WRITE (nunit,430) sbrnrm
      430 FORMAT ('   *          VALUE OF SBRNRM',  &
                  ' AT XPLUS: .................', g12.4, t74, '*')
    ELSE
      WRITE (nunit,440) sbrnrm
      440 FORMAT ('   *          VALUE OF SCALED SBRNRM',  &
                  ' AT XPLUS: ..........', g12.4, t74, '*')
    END IF
    WRITE (nunit,450) newmax
    450 FORMAT ('   *          NEWMAX: ', 35('.'), g12.4, t74, '*')
  END IF
  IF (sbrnrm < newmax) acpcod=2
!
!          FUNCTION VALUE ACCEPTANCE IS ALSO CHECKED REGARDLESS
!          OF WHETHER SECOND STEP ACCEPTANCE CRITERION WAS MET.
!
END IF
!
!       ESTABLISH DELTAF FOR USE IN COMPARISON TO PREDICTED
!       CHANGE IN OBJECTIVE FUNCTION, DELFPR, LATER.
!
deltaf=fcnnew - fcnold
IF (fcnnew >= fcnmax + alpha*delfts) THEN
!
!          FAILURE OF FIRST STEP ACCEPTANCE TEST. TEST LENGTH OF
!          STEP TO ENSURE PROGRESS IS STILL BEING MADE.
!
  rellen=zero
  DO  i=1,n
    rellen=MAX(rellen, ABS(s(i))/MAX((ABS(xplus(i))),one/scalex(i)) )
  END DO
  IF (rellen < stptol) THEN
!
!             NO PROGRESS BEING MADE - RETCOD = 7 STOPS PROGRAM.
!
    CALL matcop (n, n, 1, 1, n, 1, xc, xplus)
    retcod=7
    RETURN
  ELSE
!
!          FAILURE OF STEP BY OBJECTIVE FUNCTION CRITERION.
!          ESTABLISH A NEW DELTA FROM EITHER SIMPLE DIVISION
!          BY DELFAC OR BY FINDING THE MINIMUM OF A QUADRATIC MODEL.
!
    IF (geoms) THEN
      delta=delta/delfac
    ELSE
!
!                FIRST FIND LENGTH OF TRUST REGION STEP.
!
      wv3(1:n)=s(1:n)*scalex(1:n)
      CALL twonrm (overfl, maxexp, n, epsmch, stplen, wv3)
!
      CALL qmin (nunit, output, delfts, delta, deltaf, stplen)
!
    END IF
    IF (delta < delstr) THEN
!
!                IF DELTA HAS BEEN INCREASED AT THIS ITERATION AND THE
!                DELTA FROM QMIN IS LESS THAN THE DELTA AT THE PREVIOUSLY
!                ACCEPTED (STORED) POINT THEN RETURN TO THAT POINT AND
!                ACCEPT IT AS THE UPDATED ITERATE.
!
      CALL matcop (n, n, 1, 1, n, 1, xplpre, xplus)
      CALL matcop (n, n, 1, 1, n, 1, fplpre, fvec)
      acpcod=acpstr
      delta=delstr
      fcnnew=fcnpre
      retcod=3
      RETURN
    END IF
!
!             IF THE SECOND ACCEPTANCE TEST HAS BEEN PASSED RETURN
!             WITH NEW TRUST REGION AND CONTINUE ON TO NEXT ITERATION;
!             OTHERWISE TRY A NEW STEP WITH REDUCED DELTA.
!
    IF (acpcod == 2) THEN
      retcod=4
    ELSE
!
!                FAILURE OF FIRST STEP ACCEPTANCE TEST.
!
      retcod=11
    END IF
    RETURN
  END IF
ELSE
!
!          OBJECTIVE FUNCTION MEETS FIRST ACCEPTANCE CRITERION.
!          IN NONMONOTONIC SEARCHES IT MAY BE GREATER THAN THE
!          PREVIOUS OBJECTIVE FUNCTION VALUE - CONSIDER THIS CASE FIRST.
!
  IF (deltaf >= alpha*delfts) THEN
!
!             AN ACCEPTABLE STEP HAS BEEN FOUND FOR THE NONMONOTONIC SEARCH
!             BUT THE OBJECTIVE FUNCTION VALUE IS NOT A "DECREASE" FROM THE
!             PREVIOUS ITERATION (ACTUALLY IT MIGHT BE BETWEEN ZERO AND
!             ALPHA*DELFTS).  ACCEPT STEP BUT REDUCE DELTA.
!
    delta=delta/delfac
    retcod=5
    IF (acpcod == 2) THEN
      acpcod=12
    ELSE
      acpcod=1
    END IF
    RETURN
  END IF
!
!       COMPARE DELTAF TO DELTAF PREDICTED, DELFPR, TO DETERMINE NEXT
!       TRUST REGION SIZE.  NOTE: DELTAF MUST BE LESS THAN ALPHA*DELFTS
!       (IN ESSENCE NEGATIVE) TO HAVE REACHED THIS POINT IN TRSTUP.
!       R IS IN UPPER TRIANGLE OF MATRIX A SO THE FOLLOWING CODE FINDS:
!
!       DELFPR = DELF^S + 1/2 S^J^JS = DELF^S + 1/2 S^R^RS
!
  CALL uvmul (n, n, n, n, a, s, wv3)
  delfpr = delfts + DOT_PRODUCT( wv3(1:n), wv3(1:n) ) / two
  IF (output > 4) THEN
    WRITE (nunit,10)
    WRITE (nunit,490) delfpr
    490 FORMAT ('   *       PREDICTED CHANGE IN OBJECTIVE FUNCTION, DELFPR:', &
                g12.3, t74, '*')
    WRITE (nunit,500) deltaf
    500 FORMAT ('   *          ACTUAL CHANGE IN OBJECTIVE FUNCTION, DELTAF:', &
                g12.3, t74, '*')
  END IF
  IF (retcod <= 6 .AND. (ABS(delfpr-deltaf) <= ABS(deltaf)/  &
        ten .OR. deltaf <= delfts) .AND. (.NOT.newtkn) .AND. (.NOT.convio)  &
         .AND. delstr == zero) THEN
    IF (MIN(newlen,maxstp)/delta > onept1) THEN
!
!                PROMISING STEP - INCREASE TRUST REGION.
!
!                STORE CURRENT POINT.
!
      CALL matcop (n, n, 1, 1, n, 1, xplus, xplpre)
      CALL matcop (n, n, 1, 1, n, 1, fvec, fplpre)
      delstr=delta
      fcnpre=fcnnew
!
!                IF NONMONOTONIC STEPS ARE BEING USED EXPAND TRUST
!                REGION TO NEWLEN, OTHERWISE EXPAND BY DELFAC.
!
      IF (isejac > narmij) THEN
        delta=MIN(newlen,maxstp)
      ELSE
        delta=MIN(delfac*delta,maxstp)
      END IF
      retcod=12
      IF (acpcod == 2) THEN
        acpstr=12
      ELSE
        acpstr=1
      END IF
      acpcod=0
    ELSE
      retcod=0
      IF (acpcod == 2) THEN
        acpcod=12
      ELSE
        acpcod=1
      END IF
    END IF
    RETURN
  ELSE
!
!             CHANGE TRUST REGION SIZE DEPENDING ON DELTAF AND DELFPR.
!
    retcod=6
    IF (acpcod == 2) THEN
      acpcod=12
    ELSE
      acpcod=1
    END IF
    IF (deltaf >= delfpr/ten) THEN
      delta=delta/delfac
      IF (output > 3) THEN
        WRITE (nunit,10)
        WRITE (nunit,510)
        510 FORMAT ('   *    CHANGE IN F, DELTAF, IS > .1 DELFPR - REDUCE DELTA',  &
                    t74, '*')
      END IF
    ELSE IF (trupdm == 0 .AND. jupdm > 0) THEN
!
!                POWELL'S UPDATING SCHEME - FIND JAC S FIRST.
!
      IF (qnupdm == 0) THEN
!
!                   UNSCALED JACOBIAN IN JAC.
!
        rhs(1:n) = s(1:n)*scalef(1:n)
        CALL avmul (n, n, n, n, jac, rhs, wv3)
      ELSE
!
!                   MULTIPLY BY R FIRST.
!
        CALL uvmul (n, n, n, n, a, s, rhs)
!
!                   THEN Q (IN JAC^)
!
        CALL atbmul (n, n, 1, 1, n, n, jac, rhs, wv3)
      END IF
      dmult=delfpr/ten - deltaf
      sp=zero
      ss=zero
      DO  k=1,n
        wv3(k)=wv3(k) + fvecc(k)
        sp=sp + ABS(fvec(k)*(fvec(k) - wv3(k)))
        ss=ss + (fvec(k)-wv3(k))*(fvec(k) - wv3(k))
      END DO
      IF (sp + SQRT(sp*sp + dmult*ss) < epsmch) THEN
        powlam=ten
      ELSE
        powlam=one + dmult/(sp + SQRT(sp*sp + dmult*ss))
      END IF
      powlam=SQRT(powlam)
      powmu=MIN(delfac, powlam, powtau)
      powtau=powlam/powmu
      IF (output > 3) THEN
        WRITE (nunit,10)
        WRITE (nunit,540)
        540 FORMAT ('   *       FACTORS IN POWELL UPDATING SCHEME', t74, '*')
        WRITE (nunit,10)
        WRITE (nunit,550) powlam,powmu,powtau
        550 FORMAT ('   *       LAMBDA: ', g12.3, '    MU: ', g12.3,  &
                    '    TAU: ', g12.3, t74, '*')
        WRITE (nunit,560)
        560 FORMAT ('   *       DELTA IS MINIMUM OF MU*DELTA AND MAXSTP',  &
                    t74, '*')
      END IF
        delta=MIN(powmu*delta, maxstp)
    ELSE
      IF (deltaf < threeq*delfpr) THEN
        delta=MIN(delfac*delta, maxstp)
        IF (output > 3) THEN
          WRITE (nunit,10)
          WRITE (nunit,570)
          570 FORMAT ('   *       CHANGE IN F, DELTAF, IS LESS',  &
                      ' THAN .75 X PREDICTED', t74, '*')
          WRITE (nunit,580) delta
          580 FORMAT ('   *       DELTA INCREASED TO: ', g12.3, t74, '*')
        END IF
      ELSE
        IF (output > 3) THEN
          WRITE (nunit,10)
          WRITE (nunit,590)
          590 FORMAT ('   *       DELTAF BETWEEN 0.1 AND 0.75',  &
                      ' DELFPR - LEAVE DELTA UNCHANGED', t74, '*')
        END IF
      END IF
    END IF
  END IF
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE TRSTUP.
!
END SUBROUTINE trstup



SUBROUTINE twonrm (overfl, maxexp, n, epsmch, eucnrm, v)
!
!    FEB. 23 ,1992
!
!    THIS SUBROUTINE EVALUATES THE EUCLIDEAN NORM OF A VECTOR, V.
!    IT FOLLOWS THE ALGORITHM OF J.L. BLUE IN ACM TOMS V4 15 (1978)
!    BUT USES SLIGHTLY DIFFERENT CUTS.   THE CONSTANTS IN COMMON BLOCK
!    NNES_5 ARE CALCULATED IN THE SUBROUTINE MACHAR OR ARE PROVIDED
!    BY THE USER IN THE DRIVER.
!

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: maxexp
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: epsmch
REAL (dp), INTENT(OUT)  :: eucnrm
REAL (dp), INTENT(IN)   :: v(n)

! Local variables

REAL (dp)  :: abig, absvi, amed, asmall, sqrtep, ymax, ymin
INTEGER    :: i
!
overfl=.false.
sqrtep=SQRT(epsmch)
!
asmall=zero
amed=zero
abig=zero
DO  i=1,n
!
!          ACCUMULATE SUMS OF SQUARES IN ONE OF THREE ACCULULATORS,
!          ABIG, AMED AND ASMALL, DEPENDING ON THEIR SIZES.
!
  absvi=ABS(v(i))
!
!          THIS COMPARISON RESTRICTS THE MAXIMUM VALUE OF AMED TO BE
!          B/N => CANNOT SUM SO THAT AMED OVERFLOWS.
!
  IF (absvi > bigb/DBLE(n)) THEN
!
!             THIS DIVISOR OF 10N RESTRICTS ABIG FROM (PATHALOGICALLY)
!             OVERFLOWING FROM SUMMATION.
!
    abig=abig + ((v(i)/(ten*n*bigs)))**2
!
  ELSE IF (absvi < smalls) THEN
!
    asmall=asmall + ((v(i)/smalls))**2
!
  ELSE
!
    amed=amed + v(i)*v(i)
!
  END IF
!
END DO
IF (abig > zero) THEN
!
!          IF OVERFLOW WOULD OCCUR ASSIGN BIGR AS NORM AND SIGNAL TO
!          CALLING SUBROUTINE VIA OVERFL.
!
  IF (LOG10(abig)/two + LOG10(bigs) + 1 + LOG10(DBLE(n)) > maxexp) THEN
    eucnrm=bigr
    overfl=.true.
    RETURN
  END IF
!
!          IF AMED IS POSITIVE IT COULD CONTRIBUTE TO THE NORM -
!          DETERMINATION IS DELAYED UNTIL LATER TO SAVE REPEATING CODE.
!
  IF (amed > zero) THEN
    ymin=MIN(SQRT(amed),ten*n*bigs*SQRT(abig))
    ymax=MAX(SQRT(amed),ten*n*bigs*SQRT(abig))
  ELSE
!
!             AMED DOESN'T CONTRIBUTE AND ASMALL WON'T MATTER IF
!             ABIG IS NONZERO - FIND NORM USING ABIG AND RETURN.
!
    eucnrm=ten*n*bigs*SQRT(abig)
    RETURN
  END IF
ELSE IF (asmall > zero) THEN
  IF (amed > zero) THEN
    ymin=MIN(SQRT(amed),smalls*SQRT(asmall))
    ymax=MAX(SQRT(amed),smalls*SQRT(asmall))
  ELSE
!
!             ABIG AND AMED ARE ZERO SO USE ASMALL ONLY.
!
    eucnrm=smalls*SQRT(asmall)
    RETURN
  END IF
ELSE
  eucnrm=SQRT(amed)
  RETURN
END IF
IF (ymin < sqrtep*ymax) THEN
!
!          SMALLER PORTION DOES NOT CONTRIBUTE TO NORM.
!
  eucnrm=ymax
ELSE
!
!          SMALLER PORTION CONTRIBUTES TO NORM.
!
  eucnrm=ymax*SQRT((one+ymin/ymax)*(one+ymin/ymax))
END IF
RETURN
!
!       LAST CARD OF SUBROUTINE TWONRM.
!
END SUBROUTINE twonrm



SUBROUTINE update (mnew, mold, n, trmcod, fcnnew, fcnold, fvec,  &
                   fvecc, xc, xplus)
!
!    FEB. 9, 1991
!
!    THIS SUBROUTINE RESETS CURRENT ESTIMATES OF SOLUTION AND UPDATES THE
!    OBJECTIVE FUNCTION VALUE, M (USED TO SET HOW MANY PREVIOUS VALUES TO
!    LOOK AT IN THE NON- MONOTONIC COMPARISONS) AND THE TERMINATION CODE,
!    TRMCOD.
!

INTEGER, INTENT(IN)     :: mnew
INTEGER, INTENT(OUT)    :: mold
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(OUT)    :: trmcod
REAL (dp), INTENT(IN)   :: fcnnew
REAL (dp), INTENT(OUT)  :: fcnold
REAL (dp), INTENT(IN)   :: fvec(n)
REAL (dp), INTENT(OUT)  :: fvecc(n)
REAL (dp), INTENT(OUT)  :: xc(n)
REAL (dp), INTENT(IN)   :: xplus(n)

! Local variable

INTEGER  :: i

fcnold=fcnnew
mold=mnew
trmcod=0
DO  i=1,n
  fvecc(i)=fvec(i)
  xc(i)=xplus(i)
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE UPDATE.
!
END SUBROUTINE update



SUBROUTINE utbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact,  &
                   amat, bmat, cmat)
!
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^B=C    WHERE A IS UPPER TRIANGULAR
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NCADEC IS 2ND DIM. OF AMAT; NCAACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E.   NCADEC IS NUMBER OF COLUMNS OF A DECLARED
!           NCBDEC IS NUMBER OF COLUMNS OF B DECLARED
!           NCDEC IS THE NUMBER OF ROWS IN BOTH A AND B DECLARED
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: ncadec
INTEGER, INTENT(IN)     :: ncaact
INTEGER, INTENT(IN)     :: ncbdec
INTEGER, INTENT(IN)     :: ncbact
INTEGER, INTENT(IN)     :: ncdec
INTEGER, INTENT(IN)     :: ncact
REAL (dp), INTENT(IN)   :: amat(ncdec,ncadec)
REAL (dp), INTENT(IN)   :: bmat(ncdec,ncbdec)
REAL (dp), INTENT(OUT)  :: cmat(ncadec, ncbdec)

! Local variables

INTEGER    :: i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r, nend
REAL (dp)  :: SUM

DO  i=1,ncaact
!
!          FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
  nend=MIN(i,ncact)
!
!          THIS ADJUSTMENT IS REQUIRED WHEN NCACT IS LESS THAN NCAACT.
!
  ncc32=nend/32
  ncc32r=nend - 32*ncc32
  ncc16=ncc32r/16
  ncc16r=ncc32r - 16*ncc16
  ncc8=ncc16r/8
  ncc8r=ncc16r - 8*ncc8
  ncc4=ncc8r/4
  ncc4r=ncc8r - 4*ncc4
!
!          FIND ENTRY IN MATRIX C.
!
  DO  j=1,ncbact
    sum=zero
    k=0
    IF (ncc32 > 0) THEN
      DO  kk=1,ncc32
        k=k+32
        sum=sum + amat(k-31,i)*bmat(k-31,j)+amat(k-30,i)*bmat(k-30,  &
            j)+amat(k-29,i)*bmat(k-29,j)+amat(k-28,i)*bmat(k-28,j)+  &
            amat(k-27,i)*bmat(k-27,j)+amat(k-26,i)*bmat(k-26,j)+  &
            amat(k-25,i)*bmat(k-25,j)+amat(k-24,i)*bmat(k-24,j)
        sum=sum + amat(k-23,i)*bmat(k-23,j)+amat(k-22,i)*bmat(k-22,  &
            j)+amat(k-21,i)*bmat(k-21,j)+amat(k-20,i)*bmat(k-20,j)+  &
            amat(k-19,i)*bmat(k-19,j)+amat(k-18,i)*bmat(k-18,j)+  &
            amat(k-17,i)*bmat(k-17,j)+amat(k-16,i)*bmat(k-16,j)
        sum=sum + amat(k-15,i)*bmat(k-15,j)+amat(k-14,i)*bmat(k-14,  &
            j)+amat(k-13,i)*bmat(k-13,j)+amat(k-12,i)*bmat(k-12,j)+  &
            amat(k-11,i)*bmat(k-11,j)+amat(k-10,i)*bmat(k-10,j)+  &
            amat(k-9,i)*bmat(k-9,j)+amat(k-8,i)*bmat(k-8,j)
        sum=sum + amat(k-7,i)*bmat(k-7,j)+amat(k-6,i)*bmat(k-6,j)+  &
            amat(k-5,i)*bmat(k-5,j)+amat(k-4,i)*bmat(k-4,j)+amat(k-3,  &
            i)*bmat(k-3,j)+amat(k-2,i)*bmat(k-2,j)+amat(k-1,i)*  &
            bmat(k-1,j)+amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc16 > 0) THEN
      DO  kk=1,ncc16
        k=k+16
        sum=sum + amat(k-15,i)*bmat(k-15,j)+amat(k-14,i)*bmat(k-14,  &
            j)+amat(k-13,i)*bmat(k-13,j)+amat(k-12,i)*bmat(k-12,j)+  &
            amat(k-11,i)*bmat(k-11,j)+amat(k-10,i)*bmat(k-10,j)+  &
            amat(k-9,i)*bmat(k-9,j)+amat(k-8,i)*bmat(k-8,j)
        sum=sum + amat(k-7,i)*bmat(k-7,j)+amat(k-6,i)*bmat(k-6,j)+  &
            amat(k-5,i)*bmat(k-5,j)+amat(k-4,i)*bmat(k-4,j)+amat(k-3,  &
            i)*bmat(k-3,j)+amat(k-2,i)*bmat(k-2,j)+amat(k-1,i)*  &
            bmat(k-1,j)+amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc8 > 0) THEN
      DO  kk=1,ncc8
        k=k+8
        sum=sum + amat(k-7,i)*bmat(k-7,j)+amat(k-6,i)*bmat(k-6,j)+  &
            amat(k-5,i)*bmat(k-5,j)+amat(k-4,i)*bmat(k-4,j)+amat(k-3,  &
            i)*bmat(k-3,j)+amat(k-2,i)*bmat(k-2,j)+amat(k-1,i)*  &
            bmat(k-1,j)+amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc4 > 0) THEN
      DO  kk=1,ncc4
        k=k+4
        sum=sum + amat(k-3,i)*bmat(k-3,j)+amat(k-2,i)*bmat(k-2,j)+  &
            amat(k-1,i)*bmat(k-1,j)+amat(k,i)*bmat(k,j)
      END DO
    END IF
    IF (ncc4r > 0) THEN
      DO  kk=1,ncc4r
        k=k+1
        sum=sum + amat(k,i)*bmat(k,j)
      END DO
    END IF
    cmat(i,j)=sum
  END DO
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE UTBMUL.
!
END SUBROUTINE utbmul



SUBROUTINE utumul (nradec, ncadec, nraact, ncaact, nrbdec, ncbdec, amat, bmat)
!
!   FEB. 8, 1991
!
!   MATRIX MULTIPLICATION:   A^A=B   WHERE A IS UPPER TRIANGULAR
!
!   VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!   NRADEC IS NUMBER OF ROWS IN A DECLARED
!   NCADEC IS NUMBER OF COLUMNS IN A DECLARED
!   NRAACT IS THE LIMIT FOR THE 1ST INDEX IN A
!   NCAACT IS THE LIMIT FOR THE 2ND INDEX IN A
!   NRBDEC IS NUMBER OF ROWS IN B DECLARED
!   NCBDEC IS NUMBER OF COLUMNS IN B DECLARED
!
!   MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!   QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: ncadec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncaact
INTEGER, INTENT(IN)     :: nrbdec
INTEGER, INTENT(IN)     :: ncbdec
REAL (dp), INTENT(IN)   :: amat(nradec,ncadec)
REAL (dp), INTENT(OUT)  :: bmat(nrbdec,ncbdec)

! Local variables

INTEGER    :: i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r, nend
REAL (dp)  :: SUM

!
!       FIND ENTRY IN MATRIX B.
!
DO  i=1,ncaact
!
!          FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
  nend=MIN(i,nraact)
!
  ncc32=nend/32
  ncc32r=nend - 32*ncc32
  ncc16=ncc32r/16
  ncc16r=ncc32r - 16*ncc16
  ncc8=ncc16r/8
  ncc8r=ncc16r - 8*ncc8
  ncc4=ncc8r/4
  ncc4r=ncc8r - 4*ncc4
  DO  j=i,ncaact
    sum=zero
    k=0
    IF (ncc32 > 0) THEN
      DO  kk=1,ncc32
        k=k+32
        sum=sum + amat(k-31,i)*amat(k-31,j)+amat(k-30,i)*amat(k-30,  &
            j)+amat(k-29,i)*amat(k-29,j)+amat(k-28,i)*amat(k-28,j)+  &
            amat(k-27,i)*amat(k-27,j)+amat(k-26,i)*amat(k-26,j)+  &
            amat(k-25,i)*amat(k-25,j)+amat(k-24,i)*amat(k-24,j)
        sum=sum + amat(k-23,i)*amat(k-23,j)+amat(k-22,i)*amat(k-22,  &
            j)+amat(k-21,i)*amat(k-21,j)+amat(k-20,i)*amat(k-20,j)+  &
            amat(k-19,i)*amat(k-19,j)+amat(k-18,i)*amat(k-18,j)+  &
            amat(k-17,i)*amat(k-17,j)+amat(k-16,i)*amat(k-16,j)
        sum=sum + amat(k-15,i)*amat(k-15,j)+amat(k-14,i)*amat(k-14,  &
            j)+amat(k-13,i)*amat(k-13,j)+amat(k-12,i)*amat(k-12,j)+  &
            amat(k-11,i)*amat(k-11,j)+amat(k-10,i)*amat(k-10,j)+  &
            amat(k-9,i)*amat(k-9,j)+amat(k-8,i)*amat(k-8,j)
        sum=sum + amat(k-7,i)*amat(k-7,j)+amat(k-6,i)*amat(k-6,j)+  &
            amat(k-5,i)*amat(k-5,j)+amat(k-4,i)*amat(k-4,j)+amat(k-3,  &
            i)*amat(k-3,j)+amat(k-2,i)*amat(k-2,j)+amat(k-1,i)*  &
            amat(k-1,j)+amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc16 > 0) THEN
      DO  kk=1,ncc16
        k=k+16
        sum=sum + amat(k-15,i)*amat(k-15,j)+amat(k-14,i)*amat(k-14,  &
            j)+amat(k-13,i)*amat(k-13,j)+amat(k-12,i)*amat(k-12,j)+  &
            amat(k-11,i)*amat(k-11,j)+amat(k-10,i)*amat(k-10,j)+  &
            amat(k-9,i)*amat(k-9,j)+amat(k-8,i)*amat(k-8,j)
        sum=sum + amat(k-7,i)*amat(k-7,j)+amat(k-6,i)*amat(k-6,j)+  &
            amat(k-5,i)*amat(k-5,j)+amat(k-4,i)*amat(k-4,j)+amat(k-3,  &
            i)*amat(k-3,j)+amat(k-2,i)*amat(k-2,j)+amat(k-1,i)*  &
            amat(k-1,j)+amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc8 > 0) THEN
      DO  kk=1,ncc8
        k=k+8
        sum=sum + amat(k-7,i)*amat(k-7,j)+amat(k-6,i)*amat(k-6,j)+  &
            amat(k-5,i)*amat(k-5,j)+amat(k-4,i)*amat(k-4,j)+amat(k-3,  &
            i)*amat(k-3,j)+amat(k-2,i)*amat(k-2,j)+amat(k-1,i)*  &
            amat(k-1,j)+amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc4 > 0) THEN
      DO  kk=1,ncc4
        k=k+4
        sum=sum + amat(k-3,i)*amat(k-3,j)+amat(k-2,i)*amat(k-2,j)+  &
            amat(k-1,i)*amat(k-1,j)+amat(k,i)*amat(k,j)
      END DO
    END IF
    IF (ncc4r > 0) THEN
      DO  kk=1,ncc4r
        k=k+1
        sum=sum + amat(k,i)*amat(k,j)
      END DO
    END IF
    bmat(i,j)=sum
    IF (i /= j) bmat(j,i)=bmat(i,j)
  END DO
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE UTUMUL.
!
END SUBROUTINE utumul



SUBROUTINE uvmul (nradec, nraact, ncdec, ncact, amat, bvec, cvec)
!
!    FEB. 8, 1991
!
!    MATRIX-VECTOR MULTIPLICATION:  AB=C  WHERE A IS UPPER TRIANGULAR
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BVEC; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION (COLUMNS OF A AND ROWS OF B)
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!

INTEGER, INTENT(IN)     :: nradec
INTEGER, INTENT(IN)     :: nraact
INTEGER, INTENT(IN)     :: ncdec
INTEGER, INTENT(IN)     :: ncact
REAL (dp), INTENT(IN)   :: amat(nradec,ncdec)
REAL (dp), INTENT(IN)   :: bvec(ncdec)
REAL (dp), INTENT(OUT)  :: cvec(nradec)

! Local variables

INTEGER    :: i, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r
REAL (dp)  :: SUM

DO  i=1,nraact
!
!          FIND NUMBER OF GROUPS OF SIZE 32, 16 ...
!
  ncc32=(ncact - (i-1))/32
  ncc32r=(ncact - (i-1))-32*ncc32
  ncc16=ncc32r/16
  ncc16r=ncc32r - 16*ncc16
  ncc8=ncc16r/8
  ncc8r=ncc16r - 8*ncc8
  ncc4=ncc8r/4
  ncc4r=ncc8r - 4*ncc4
!
!          FIND ENTRY FOR VECTOR C.
!
  sum=zero
  k=i-1
  IF (ncc32 > 0) THEN
    DO  kk=1,ncc32
      k=k+32
      sum=sum + amat(i,k-31)*bvec(k-31) + amat(i,k-30)*bvec(k-30) +   &
          amat(i,k-29)*bvec(k-29) + amat(i,k-28)*bvec(k-28) +   &
          amat(i,k-27)*bvec(k-27) + amat(i,k-26)*bvec(k-26) +   &
          amat(i,k-25)*bvec(k-25) + amat(i,k-24)*bvec(k-24)
      sum=sum + amat(i,k-23)*bvec(k-23) + amat(i,k-22)*bvec(k-22) +   &
          amat(i,k-21)*bvec(k-21) + amat(i,k-20)*bvec(k-20) +   &
          amat(i,k-19)*bvec(k-19) + amat(i,k-18)*bvec(k-18) +   &
          amat(i,k-17)*bvec(k-17) + amat(i,k-16)*bvec(k-16)
      sum=sum + amat(i,k-15)*bvec(k-15) + amat(i,k-14)*bvec(k-14) +   &
          amat(i,k-13)*bvec(k-13) + amat(i,k-12)*bvec(k-12) +   &
          amat(i,k-11)*bvec(k-11) + amat(i,k-10)*bvec(k-10) +   &
          amat(i,k-9)*bvec(k-9) + amat(i,k-8)*bvec(k-8)
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) + &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc16 > 0) THEN
    DO  kk=1,ncc16
      k=k+16
      sum=sum + amat(i,k-15)*bvec(k-15) + amat(i,k-14)*bvec(k-14) +   &
          amat(i,k-13)*bvec(k-13) + amat(i,k-12)*bvec(k-12) +   &
          amat(i,k-11)*bvec(k-11) + amat(i,k-10)*bvec(k-10) +   &
          amat(i,k-9)*bvec(k-9) + amat(i,k-8)*bvec(k-8)
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) + &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc8 > 0) THEN
    DO  kk=1,ncc8
      k=k+8
      sum=sum + amat(i,k-7)*bvec(k-7) + amat(i,k-6)*bvec(k-6) +   &
          amat(i,k-5)*bvec(k-5) + amat(i,k-4)*bvec(k-4) + amat(i,k-3)*bvec(k-3) + &
          amat(i,k-2)*bvec(k-2) + amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc4 > 0) THEN
    DO  kk=1,ncc4
      k=k+4
      sum=sum + amat(i,k-3)*bvec(k-3) + amat(i,k-2)*bvec(k-2) +   &
                amat(i,k-1)*bvec(k-1) + amat(i,k)*bvec(k)
    END DO
  END IF
  IF (ncc4r > 0) THEN
    DO  kk=1,ncc4r
      k=k+1
      sum=sum + amat(i,k)*bvec(k)
    END DO
  END IF
  cvec(i)=sum
END DO
RETURN
!
!       LAST CARD OF SUBROUTINE UVMUL.
!
END SUBROUTINE uvmul

END MODULE Equation_Solver

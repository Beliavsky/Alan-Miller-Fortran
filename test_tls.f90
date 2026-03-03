PROGRAM test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-04-27  Time: 23:20:16

USE Total_LS
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, PARAMETER  :: in = 5, iout = 6
INTEGER, PARAMETER  :: ldc = 30, ldx = 10
INTEGER    :: m, n, l, nl, rank, ierr, iwarn
LOGICAL    :: inul(11)
REAL (dp)  :: c(ldc,11), x(ldx,10), q(22), theta, tol1, tol2
INTEGER    :: i, k
! INTRINSIC MIN

READ (in,*) m, n, l, rank, theta
WRITE (iout,5800) m, n, l, rank, theta
nl = n + l
WRITE (iout,5400)
DO  i = 1, m
  READ (in,*) c(i,1:nl)
END DO
WRITE (iout,5100)
DO  i = 1, m
  WRITE (iout,5200) c(i,1:nl)
END DO
tol1 = 1.0D-08
tol2 = 1.0D-10
CALL ptls(c,ldc,m,n,l,rank,theta,x,ldx,q,inul,tol1,tol2,ierr,iwarn)
WRITE (iout,5300) ierr, iwarn
WRITE (iout,6300) rank, theta
k = MIN(m,nl)
WRITE (iout,5500) q(1:k)
WRITE (iout,5600) q(k+2:2*k)
WRITE (iout,5900)
k = 0
DO  i = 1, nl
  IF (inul(i)) THEN
    k = k + 1
    WRITE (iout,6000) k, i
    WRITE (iout,5200) c(1:nl,i)
  END IF
END DO
WRITE (iout,6100)
DO  i = 1, l
  WRITE (iout,6200) i, x(1:n,i)
END DO
STOP

5100 FORMAT (/'      C(*,I)        C(*,I+1)       C(*,I+2)       C(*,I+3)' )
5200 FORMAT (' ',4g15.6)
5300 FORMAT (/' IERR = ',i3,'   IWARN = ',i3)
5400 FORMAT (/' PARTIAL TOTAL LEAST SQUARES TEST PROGRAM'/' ',40('-')/)
5500 FORMAT (/' DIAGONAL OF THE PARTIALLY DIAGONALIZED BIDIAGONAL='/  &
             (' ',4g15.6))
5600 FORMAT (/' SUPERDIAGONAL OF THE PARTIALLY DIAGONALIZED ',  &
             'BIDIAGONAL='/ (' ',4g15.6))
5800 FORMAT (/' M =',i3,'  N =',i3,'  L =',i3,'  RANK =',i3, '  THETA =',g12.5)
5900 FORMAT (/' BASIS OF THE COMPUTED RIGHT SINGULAR SUBSPACE :'/ ' ',45( '-'))
6000 FORMAT (' THE ',i2,'-TH BASE VECTOR V(*,',i2,') =')
6100 FORMAT (/' TLS SOLUTION :'/ ' ',12('*'))
6200 FORMAT (/' X(*,',i2,') = '/ (' ',4g15.6))
6300 FORMAT (/' COMPUTED RANK =',i3,'     COMPUTED BOUND THETA = ',g12.5 )
END PROGRAM test

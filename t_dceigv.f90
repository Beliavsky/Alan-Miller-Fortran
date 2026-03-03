PROGRAM test_dceigv

! Generate a random complex matrix and calculate its
! eigenvalues & vectors.

USE complex_eigen
USE constants_NSWC
IMPLICIT NONE

INTEGER                   :: i, ierr, j, n
REAL (dp)                 :: e
REAL (dp), ALLOCATABLE    :: ar(:,:), ai(:,:), wr(:), wi(:), zr(:,:), zi(:,:)
COMPLEX (dp), ALLOCATABLE :: a(:,:), lhs(:), rhs(:)
CHARACTER (LEN=1)         :: key

e = EXP(1.0_dp)
WRITE(*, '(a)', ADVANCE='NO') ' Enter size of matrix: '
READ(*, *) n

ALLOCATE( ar(n,n), ai(n,n), wr(n), wi(n), zr(n,n), zi(n,n), a(n,n),  &
          lhs(n), rhs(n) )

! Generate random matrices ar & ai

CALL RANDOM_NUMBER( ar )
CALL RANDOM_NUMBER( ai )
ar = e * ar
ai = e * ai
a = CMPLX( ar, ai, KIND=dp )

CALL dceigv(.TRUE., ar, ai, n, wr, wi, zr, zi, ierr)

IF (ierr == 0) THEN
  WRITE(*, '(a)') ' Eigenvalues - real parts'
  WRITE(*, '(" ", 6g13.5)') wr
  WRITE(*, '(a)') ' Eigenvalues - imaginary parts'
  WRITE(*, '(" ", 6g13.5)') wi

  DO j = 1, n
    WRITE(*, '(a, i6)') ' Eigenvector number', j
    DO i = 1, n-1, 2
      WRITE(*, '( 2(" (", g13.5, ", ", g13.5, ")     ") )')  &
            zr(i,j), zi(i,j), zr(i+1,j), zi(i+1,j)
    END DO
    IF (n /= 2*(n/2)) THEN
      WRITE(*, '( " (", g13.5, ", ", g13.5, ")     " )') zr(n,j), zi(n,j)
    END IF
  END DO

  WRITE(*, *) 'Press return to continue'
  READ(*, '(a)') key
  IF (ICHAR(key) == 255) STOP

! Test that A.(eigenvector) = (eigenvalue).(eigenvector)

  DO j = 1, n
    lhs = MATMUL( a, CMPLX( zr(:,j), zi(:,j), KIND=dp ) )
    rhs = CMPLX( wr(j), wi(j), KIND=dp ) * CMPLX( zr(:,j), zi(:,j), KIND=dp )
    WRITE(*, '(i5, "  A.v ="/ (" ", 6g13.5) )') j, lhs
    WRITE(*, '(i5, "  lambda.v ="/ (" ", 6g13.5) )') j, rhs
  END DO

ELSE
  WRITE(*, *) 'ierr =', ierr
END IF

STOP
END PROGRAM test_dceigv

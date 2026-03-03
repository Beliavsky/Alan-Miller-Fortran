MODULE Sparse_Gauss

! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-02  Time: 09:35:06

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


SUBROUTINE nspiv(n, ia, ja, a, b, MAX, r, c, ic, x, ierr)

! N.B. Arguments ITEMP & RTEMP have been removed.

!  NSPIV USES SPARSE GAUSSIAN ELIMINATION WITH COLUMN INTERCHANGES TO SOLVE
!  THE LINEAR SYSTEM A X = B.  THE ELIMINATION PHASE PERFORMS ROW OPERATIONS
!  ON A AND B TO OBTAIN A UNIT UPPER TRIANGULAR MATRIX U AND A VECTOR Y.
!  THE SOLUTION PHASE SOLVES U X = Y.

!  INPUT ARGUMENTS---

!  N      INTEGER NUMBER OF EQUATIONS AND UNKNOWNS

!  IA     INTEGER ARRAY OF N+1 ENTRIES CONTAINING ROW POINTERS TO A
!         (SEE MATRIX STORAGE DESCRIPTION BELOW)

!  JA     INTEGER ARRAY WITH ONE ENTRY PER NONZERO IN A, CONTAINING COLUMN
!         NUMBERS OF THE NONZEROES OF A.  (SEE MATRIX STORAGE DESCRIPTION
!         BELOW)

!  A      REAL ARRAY WITH ONE ENTRY PER NONZERO IN A, CONTAINING THE
!         ACTUAL NONZEROES.  (SEE MATRIX STORAGE DESCRIPTION BELOW)

!  B      REAL ARRAY OF N ENTRIES CONTAINING RIGHT HAND SIDE DATA

!  MAX    INTEGER NUMBER SPECIFYING MAXIMUM NUMBER OF OFF-DIAGONAL
!         NONZERO ENTRIES OF U WHICH MAY BE STORED

!  R      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         ROWS OF A (I.E., THE ELIMINATION ORDER FOR THE EQUATIONS)

!  C      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         COLUMNS OF A.  C IS ALSO AN OUTPUT ARGUMENT

!  IC     INTEGER ARRAY OF N ENTRIES WHICH IS THE INVERSE OF C
!         (I.E., IC(C(I)) = I).  IC IS ALSO AN OUTPUT ARGUMENT

!  ITEMP  INTEGER ARRAY OF 2*N + MAX + 2 ENTRIES, FOR INTERNAL USE

!  RTEMP  REAL ARRAY OF N + MAX ENTRIES FOR INTERNAL USE


!  OUTPUT ARGUMENTS---

!  C      INTEGER ARRAY OF N ENTRIES SPECIFYING THE ORDER OF THE
!         COLUMNS OF U.  C IS ALSO AN INPUT ARGUMENT

!  IC     INTEGER ARRAY OF N ENTRIES WHICH IS THE INVERSE OF C
!         (I.E., IC(C(I)) = I).  IC IS ALSO AN INPUT ARGUMENT

!  X      REAL ARRAY OF N ENTRIES CONTAINING THE SOLUTION VECTOR

!  IERR   INTEGER NUMBER WHICH INDICATES ERROR CONDITIONS OR THE ACTUAL
!         NUMBER OF OFF-DIAGONAL ENTRIES IN U (FOR SUCCESSFUL COMPLETION)

!         IERR VALUES ARE---

!         0 < IERR              SUCCESSFUL COMPLETION.  U HAS IERR
!                               OFF-DIAGONAL NONZERO ENTRIES

!         IERR = 0              ERROR.  N = 0

!         -N <= IERR < 0        ERROR.  ROW NUMBER ABS(IERR) OF A IS NULL

!         -2*N <= IERR < -N     ERROR.  ROW NUMBER ABS(IERR+N) HAS A
!                               DUPLICATE ENTRY

!         -3*N <= IERR < -2*N   ERROR.  ROW NUMBER ABS(IERR+2*N)
!                               HAS A ZERO PIVOT

!         -4*N <= IERR < -3*N   ERROR.  ROW NUMBER ABS(IERR+3*N)
!                               EXCEEDS STORAGE

!  STORAGE OF SPARSE MATRICES---

!  THE SPARSE MATRIX A IS STORED USING THREE ARRAYS IA, JA, AND A.
!  THE ARRAY A CONTAINS THE NONZEROES OF THE MATRIX ROW-BY-ROW, NOT
!  NECESSARILY IN ORDER OF INCREASING COLUMN NUMBER.  THE ARRAY JA CONTAINS
!  THE COLUMN NUMBERS CORRESPONDING TO THE NONZEROES STORED IN THE ARRAY A
!  (I.E., IF THE NONZERO STORED IN A(K) IS IN COLUMN J, THEN JA(K) = J).
!  THE ARRAY IA CONTAINS POINTERS TO THE ROWS OF NONZEROES/COLUMN INDICES IN
!  THE ARRAY A/JA (I.E., A(IA(I))/JA(IA(I)) IS THE FIRST ENTRY FOR ROW I IN
!  THE ARRAY A/JA).
!  IA(N+1) IS SET SO THAT IA(N+1) - IA(1) = THE NUMBER OF NONZEROES IN A

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: ia(:)
INTEGER, INTENT(IN)      :: ja(:)
REAL (dp), INTENT(IN)    :: a(:)
REAL (dp), INTENT(IN)    :: b(:)
INTEGER, INTENT(IN)      :: MAX
INTEGER, INTENT(IN)      :: r(:)
INTEGER, INTENT(IN OUT)  :: c(:)
INTEGER, INTENT(IN OUT)  :: ic(:)
REAL (dp), INTENT(OUT)   :: x(:)
INTEGER, INTENT(OUT)     :: ierr

! Local variables

REAL (dp)  :: y(n), u(MAX)
INTEGER    :: p(n+1), iu(n+1), ju(MAX)
REAL (dp)  :: dk, lki, xpv, xpvmax, yk
INTEGER    :: i, j, jaj, jmin, jmax, juj, juptr, k, maxc, maxcl, nzcnt
INTEGER    :: ck, pk, ppk, pv, v, vi, vj, vk
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

IF (n /= 0) THEN
  
!  INITIALIZE WORK STORAGE AND POINTERS TO JU
  
  x(1:n) = zero
  iu(1) = 1
  juptr = 0
  
!  PERFORM SYMBOLIC AND NUMERIC FACTORIZATION ROW BY ROW
!  VK (VI,VJ) IS THE GRAPH VERTEX FOR ROW K (I,J) OF U
  
  DO  k = 1, n
    
!  INITIALIZE LINKED LIST AND FREE STORAGE FOR THIS ROW
!  THE R(K)-TH ROW OF A BECOMES THE K-TH ROW OF U.
    
    p(n+1) = n + 1
    vk = r(k)
    
!  SET UP ADJACENCY LIST FOR VK, ORDERED IN CURRENT COLUMN ORDER OF U.
!  THE LOOP INDEX GOES DOWNWARD TO EXPLOIT ANY COLUMNS
!  FROM A IN CORRECT RELATIVE ORDER
    
    jmin = ia(vk)
    jmax = ia(vk+1) - 1
    IF (jmin > jmax) GO TO 170
    j = jmax
    20 jaj = ja(j)
    vj = ic(jaj)
    
!  STORE A(K,J) IN WORK VECTOR
    
    x(vj) = a(j)
!  THIS CODE INSERTS VJ INTO ADJACENCY LIST OF VK
    ppk = n + 1
    30 pk = ppk
    ppk = p(pk)
    IF (ppk-vj < 0.0) THEN
      GO TO 30
    ELSE IF (ppk-vj == 0.0) THEN
      GO TO 180
    END IF
    p(vj) = ppk
    p(pk) = vj
    j = j - 1
    IF (j >= jmin) GO TO 20
    
!  THE FOLLOWING CODE COMPUTES THE K-TH ROW OF U
    
    vi = n + 1
    yk = b(vk)
    50 vi = p(vi)
    IF (vi < k) THEN
      
!  VI LT VK -- PROCESS THE L(K,I) ELEMENT AND MERGE THE
!  ADJACENCY OF VI WITH THE ORDERED ADJACENCY OF VK
      
      lki = -x(vi)
      x(vi) = zero
      
!  ADJUST RIGHT HAND SIDE TO REFLECT ELIMINATION
      
      yk = yk + lki * y(vi)
      ppk = vi
      jmin = iu(vi)
      jmax = iu(vi+1) - 1
      IF (jmin > jmax) GO TO 50
      DO  j = jmin, jmax
        juj = ju(j)
        vj = ic(juj)
        
!  IF VJ IS ALREADY IN THE ADJACENCY OF VK,
!  SKIP THE INSERTION
        
        IF (x(vj) == zero) THEN
          
!  INSERT VJ IN ADJACENCY LIST OF VK.
!  RESET PPK TO VI IF WE HAVE PASSED THE CORRECT INSERTION SPOT.
!  (THIS HAPPENS WHEN THE ADJACENCY OF
!  VI IS NOT IN CURRENT COLUMN ORDER DUE TO PIVOTING.)
          
          IF (vj-ppk < 0.0) THEN
            GO TO 60
          ELSE IF (vj-ppk == 0.0) THEN
            GO TO 90
          ELSE
            GO TO 70
          END IF
          60 ppk = vi
          70 pk = ppk
          ppk = p(pk)
          IF (ppk-vj < 0.0) THEN
            GO TO 70
          ELSE IF (ppk-vj == 0.0) THEN
            GO TO 90
          END IF
          p(vj) = ppk
          p(pk) = vj
          ppk = vj
        END IF
        
!  COMPUTE L(K,J) = L(K,J) - L(K,I)*U(I,J) FOR L(K,I) NONZERO
!  COMPUTE U*(K,J) = U*(K,J) - L(K,I)*U(I,J) FOR U(K,J) NONZERO
!  (U*(K,J) = U(K,J)*D(K,K))
        
        90 x(vj) = x(vj) + lki * u(j)
      END DO
      GO TO 50
    END IF
    
!  PIVOT--INTERCHANGE LARGEST ENTRY OF K-TH ROW OF U WITH THE DIAGONAL ENTRY.
    
!  FIND LARGEST ENTRY, COUNTING OFF-DIAGONAL NONZEROES
    
    IF (vi > n) GO TO 190
    xpvmax = ABS(x(vi))
    maxc = vi
    nzcnt = 0
    pv = vi
    110 v = pv
    pv = p(pv)
    IF (pv <= n) THEN
      nzcnt = nzcnt + 1
      xpv = ABS(x(pv))
      IF (xpv <= xpvmax) GO TO 110
      xpvmax = xpv
      maxc = pv
      maxcl = v
      GO TO 110
    END IF
    IF (xpvmax == zero) GO TO 190
    
!  IF VI = K, THEN THERE IS AN ENTRY FOR DIAGONAL WHICH MUST BE DELETED.
!  OTHERWISE, DELETE THE ENTRY WHICH WILL BECOME THE DIAGONAL ENTRY
    
    IF (vi /= k) THEN
      IF (vi /= maxc) THEN
        p(maxcl) = p(maxc)
        GO TO 120
      END IF
    END IF
    vi = p(vi)
    
!  COMPUTE D(K) = 1/L(K,K) AND PERFORM INTERCHANGE.
    
    120 dk = one / x(maxc)
    x(maxc) = x(k)
    i = c(k)
    c(k) = c(maxc)
    c(maxc) = i
    ck = c(k)
    ic(ck) = k
    ic(i) = maxc
    x(k) = zero
    
!  UPDATE RIGHT HAND SIDE.
    
    y(k) = yk * dk
    
!  COMPUTE VALUE FOR IU(K+1) AND CHECK FOR STORAGE OVERFLOW
    
    iu(k+1) = iu(k) + nzcnt
    IF (iu(k+1) > MAX+1) GO TO 200
    
!  MOVE COLUMN INDICES FROM LINKED LIST TO JU.
!  COLUMNS ARE STORED IN CURRENT ORDER WITH ORIGINAL
!  COLUMN NUMBER (C(J)) STORED FOR CURRENT COLUMN J
    
    IF (vi <= n) THEN
      j = vi
      130 juptr = juptr + 1
      ju(juptr) = c(j)
      u(juptr) = x(j) * dk
      x(j) = zero
      j = p(j)
      IF (j <= n) GO TO 130
    END IF
  END DO
  
!  BACKSOLVE U X = Y, AND REORDER X TO CORRESPOND WITH A
  
  k = n
  DO  i = 1, n
    yk = y(k)
    jmin = iu(k)
    jmax = iu(k+1) - 1
    IF (jmin <= jmax) THEN
      DO  j = jmin, jmax
        juj = ju(j)
        juj = ic(juj)
        yk = yk - u(j) * y(juj)
      END DO
    END IF
    y(k) = yk
    ck = c(k)
    x(ck) = yk
    k = k - 1
  END DO
  
!  RETURN WITH IERR = NUMBER OF OFF-DIAGONAL NONZEROES IN U
  
  ierr = iu(n+1) - iu(1)
  RETURN
END IF

!  ERROR RETURNS

!  N = 0

ierr = 0
RETURN

!  ROW K OF A IS NULL

170 ierr = -k
RETURN

!  ROW K OF A HAS A DUPLICATE ENTRY

180 ierr = -(n+k)
RETURN

!  ZERO PIVOT IN ROW K

190 ierr = -(2*n+k)
RETURN

!  STORAGE FOR U EXCEEDED ON ROW K

200 ierr = -(3*n+k)
RETURN
END SUBROUTINE nspiv



SUBROUTINE preord(n, ia, r, c, ic)

!  PREORD ORDERS THE ROWS OF A BY INCREASING NUMBER OF NONZEROES.
!  THE ROW PERMUTATION IS RETURNED IN R.  C IS SET TO THE IDENTITY.

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: ia(:)
INTEGER, INTENT(OUT)  :: r(:)
INTEGER, INTENT(OUT)  :: c(:)
INTEGER, INTENT(OUT)  :: ic(:)

INTEGER  :: i, j, k, kdeg

DO  i = 1, n
  r(i) = i
  c(i) = i
  ic(i) = i
END DO
c(1:n) = 0

DO  k = 1, n
  kdeg = ia(k+1) - ia(k)
  IF (kdeg == 0) kdeg = kdeg + 1
  ic(k) = c(kdeg)
  c(kdeg) = k
END DO

i = 0
DO  j = 1, n
  IF (c(j) /= 0) THEN
    k = c(j)
    DO
      i = i + 1
      r(i) = k
      k = ic(k)
      IF (k <= 0) EXIT
    END DO
  END IF
END DO
DO  i = 1, n
  c(i) = i
  ic(i) = i
END DO
RETURN
END SUBROUTINE preord



SUBROUTINE reschk(n, ia, ja, a, b, x)

!  RESCHK COMPUTES THE MAX-NORM AND 2-NORM OF THE RESIDUAL.
!  REAL (dp) IS USED FOR THE COMPUTATION.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ia(:)
INTEGER, INTENT(IN)        :: ja(:)
REAL (dp), INTENT(IN OUT)  :: a(:)
REAL (dp), INTENT(IN OUT)  :: b(:)
REAL (dp), INTENT(IN OUT)  :: x(:)

REAL (dp)  :: resid, residm, rowsum
INTEGER    :: i, j, jaj, jmin, jmax

resid = 0.
residm = 0.
DO  i = 1, n
  rowsum = b(i)
  jmin = ia(i)
  jmax = ia(i+1) - 1
  DO  j = jmin, jmax
    jaj = ja(j)
    rowsum = rowsum - a(j) * x(jaj)
  END DO
  IF (ABS(rowsum) > residm) residm = ABS(rowsum)
  resid = resid + rowsum ** 2
END DO
resid = SQRT(resid)
WRITE (6,5000) resid
WRITE (6,5100) residm
RETURN

5000 FORMAT (' 2-NORM OF RESIDUAL   = ', g14.7)
5100 FORMAT (' MAX NORM OF RESIDUAL = ', g14.7)
END SUBROUTINE reschk

END MODULE Sparse_Gauss

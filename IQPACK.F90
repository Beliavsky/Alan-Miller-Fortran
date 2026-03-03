MODULE iqpack

! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-18  Time: 22:26:02

!  ALGORITHM 655, COLLECTED ALGORITHMS FROM ACM.
!  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!  VOL. 13, NO. 4, P. 399.

!            IQPACK - FORTRAN SUBROUTINES FOR THE WEIGHTS OF
!                       INTERPOLATORY QUADRATURES

!  FOR A DETAILED DESCRIPTION OF THESE ROUTINES SEE THE PAPER
!  WITH THE ABOVE TITLE -

!  GIVEN A SET OF DISTINCT KNOTS, T, AND THEIR MULTIPLICITIES MLT,
!  THIS PACKAGE COMPUTES THE WEIGHTS D    OF THE INTERPOLATORY
!                                     J,I
!  QUADRATURE FORMULA

!                                         (I)
!              SUM       SUM        D    F   (T(J)),
!            J=1,NT   I=0,MLT(J)-1   J,I

!         (I)
!  WHERE F    IS THE I-TH DERIVATIVE OF F, WHERE THE QUADRATURE
!  IS TO APPROXIMATE

!                  INTEGRAL  F(T)W(T) DT,
!                   ÕA,Bå

!  AND WHERE W(T) IS A WEIGHT FUNCTION. FOR CERTAIN CLASSICAL WEIGHT FUNCTIONS,
!  LISTED BELOW, NO OTHER INFORMATION IS NEEDED.  HOWEVER THE PACKAGE CAN
!  COMPUTE THE QUADRATURE WEIGHTS CORRESPONDING TO ANY W(T) FOR WHICH THE
!  ZERO-TH MOMENT AND THE (TRIDIAGONAL SYMMETRIC) JACOBI MATRIX ASSOCIATED
!  WITH THE POLYNOMIALS ORTHOGONAL ON ÕA,Bå WITH RESPECT TO W(T), ARE KNOWN.
!  (A UTILITY ROUTINE IS SUPPLIED TO PROVIDE THIS INFORMATION FOR CLASSICAL
!  WEIGHT FUNCTIONS).  KNOTS AND WEIGHTS OF GAUSS QUADRATURES WITH NO MULTIPLE
!  KNOTS CAN ALSO BE COMPUTED.

!  THE PACKAGE IS AN IMPLEMENTATION OF THE METHOD DESCRIBED IN

!  "CALCULATION OF THE WEIGHTS OF INTERPOLATORY QUADRATURES",
!  J. KAUTSKY AND S. ELHAY,  NUMER MATH 40 (1982) 407-422,

!  TOGETHER WITH VARIOUS UTILITY ROUTINES. WEIGHTS TO SOME OR ALL THE
!  KNOTS CAN BE COMPUTED.

!        TABLE OF CLASSICAL WEIGHT FUNCTIONS

!   KIND  INTERVAL         WEIGHT FUNCTION                NAME
!     1    (A,B)                ONE                     LEGENDRE
!     2    (A,B)        ((B-X)*(X-A))**(-HALF)          CHEBYSHEV
!     3    (A,B)        ((B-X)*(X-A))**ALPHA            GEGENBAUER
!     4    (A,B)      (B-X)**ALPHA*(X-A)**BETA           JACOBI
!     5   (A,INF)     (X-A)**ALPHA*EXP(-B*(X-A))       GEN LAGUERRE
!     6  (-INF,INF)  ABS(X-A)**ALPHA*EXP(-B*(X-A)**2)  GEN HERMITE
!     7    (A,B)        ABS(X-(A+B)/TWO)**ALFA          EXPONENTIAL
!     8   (A,INF)      (X-A)**ALPHA*(B+X)**BETA          RATIONAL

!  THE VALUES B=1 AND
!                      A=-1 FOR WEIGHT FUNCTIONS 1,2,3,4,7
!                      A= 0  FOR WEIGHT FUNCTIONS 5,6,8
!  WILL BE REFERRED TO AS THE DEFAULT VALUES.

!  WE ALSO DEFINE DEL AS
!                      (A+B)/2 FOR WEIGHT FUNCTIONS 1,2,3,4,7
!                        A     FOR WEIGHT FUNCTIONS 5,6,8

!   IQPACK INDEX
!   ------------

!   LEGEND
!   ------
!   GENERALLY I = THIS QUANTITY IS INPUT TO THIS ROUTINE
!             O = THIS QUANTITY IS OUTPUT FROM THIS ROUTINE
!   KNOTS -   M = MULTIPLE KNOTS ALLOWED
!             S = ONLY SIMPLE KNOTS ALLOWED
!   WEIGHTS - C = COMPUTED
!   QF -      I = ANY INTERPOLATORY QUADRATURE FORMULA
!             G = GAUSSIAN QUADRATURE FORMULA
!   EVAL -    Y = THE QUADRATURE SUM IS FORMED
!             N = THE QUADRATURE SUM IS NOT FORMED
!   PRINT -   Y = THE KNOTS AND WEIGHTS OF THE QUADRATURE FORMULA ARE
!                 OPTIONALLY PRINTED AND A CHECK OF THE MOMENTS IS
!                 OPTIONALLY PRINTED
!             N = NO PRINTING POSSIBLE
!   A,B -     A = ANY VALID VALUES OF THE WEIGHT FUNCTION PARAMETERS A,B
!                 ALLOWED
!             D = ONLY THE DEFAULT VALUES OF A,B ALLOWED

!   USER ROUTINES
!   -------------
!         NAME       KNOTS  WEIGHTS  QF  EVAL  PRINT  A,B
!         ------
!         CEGQF      SO     C        G   Y     N      A
!         CEGQFS     SO     C        G   Y     N      D
!         CGQF       SO     OC       G   N     Y      A
!         CGQFS      SO     OC       G   N     Y      D
!         CDGQF      SO     OC       G   N     N      D
!         SGQF       SO     OC       G   N     N      -
!         CLIQF      SI     OC       I   N     Y      A
!         CLIQFS     SI     OC       I   N     Y      D
!         CEIQF      MI     C        I   Y     N      A
!         CEIQFS     MI     C        I   Y     N      D
!         CIQF       MI     CO       I   N     Y      A
!         CIQFS      MI     CO       I   N     Y      D
!         EIQF       MI     I        I   Y     N      -
!         EIQFS      SI     I        I   Y     N      -
!         CAWIQ      MI     C        I   N     N      D

!   UTILITY AND AUXILIARY ROUTINES
!   ------------------------------
!         CLASS      COMPUTE THE ZERO-TH MOMENT AND JACOBI MATRIX FOR
!                    A CLASSICAL WEIGHT FUNCTION
!         WM         COMPUTE THE MOMENTS OF A CLASSICAL WEIGHT FUNCTION
!         PARCHK     CHECK THAT THE PARAMETER VALUES ARE VALID FOR THIS
!                    WEIGHT FUNCTION
!         CHKQFS     CHECK AND OPTIONALLY PRINT A MOMENTS CHECK OF A QF
!                    AND OPTIONALLY PRINT THE KNOTS AND WEIGHTS.  DEFAULT
!                    VALUES OF A,B ONLY
!         CHKQF      CHECK AND OPTIONALLY PRINT A MOMENTS CHECK OF A QF
!                    AND OPTIONALLY PRINT THE KNOTS AND WEIGHTS.  ANY
!                    VALID VALUES OF A,B ALLOWED
!         SCT        SCALE THE KNOTS OF A QF FOR ANY VALID A,B TO THOSE
!                    FOR THE DEFAULT VALUES OF A,B
!         SCQF       SCALE A CLASSICAL WEIGHT FUNCTION QF WITH DEFAULT
!                    VALUES FOR A,B TO THOSE FOR ANY VALID A,B
!         SCMM       SCALE THE MOMENTS OF A CLASSICAL WEIGHT FUNCTION
!                    WITH DEFAULT VALUES FOR A,B TO THOSE FOR ANY VALID A,B
!         WTFN       COMPUTE THE VALUES OF A CLASSICAL WEIGHT FUNCTION
!                    AT GIVEN POINTS
!         CWIQD      FIND ALL THE WEIGHTS TO 1 MULTIPLE KNOT OF A QF
!         IMTQLX     ORTHOGONALLY DIAGONALIZE A JACOBI MATRIX

!     THE FOLLOWING IS A LIST OF PARAMETERS USED THROUGHOUT THE PACKAGE
!     WHICH ALWAYS HAVE THE SAME MEANING.

! NT    NUMBER OF DISTINCT KNOTS.  MUST BE >=1.
! T     KNOT ARRAY.
! MLT   MULTIPLICITY ARRAY.  T(J) HAS MULTIPLICITY MLT(J).
! NWTS  DIMENSION OF THE ARRAY CONTAINING THE WEIGHTS.
! WTS   ARRAY CONTAINING THE WEIGHTS.
! NDX   FLAGS AND POINTERS ARRAY.  THE PACKAGE HAS BEEN DESIGNED TO
!       (1) TREAT ALL OR ONLY SOME OF THE KNOTS SUPPLIED AS INCLUDED IN
!           THE QUADRATURE,
!       (2) COMPUTE THE WEIGHTS FOR ALL OR ONLY SOME OF THE KNOTS
!           INCLUDED IN THE QUADRATURE,
!       (3) TO PACK THE WEIGHTS IN THE OUTPUT ARRAY IN VARIOUS (POSSIBLY
!           FOUR) DIFFERENT WAYS.
!       NDX INDICATES THE STATUS OF EACH KNOT AND POINTS TO THE LOCATION
!       OF THAT KNOT IN THE WTS ARRAY. ITS USE IS DESCRIBED IN CAWIQ.
!       IN MOST STRAIGHTFORWARD APPLICATIONS THE USER WILL ONLY NEED TO
!       DIMENSION THE ARRAY. THE PACKAGE WILL DO THE REST.
! KEY   WEIGHTS ARRAY STRUCTURE FLAG. WILL USUALLY BE SET 1.  USE
!       DESCRIBED IN CAWIQ.
! KIND  AN INTEGER 0 <= KIND <= 8 SPECIFYING WHICH WEIGHT FUNCTION IS TO BE
!       USED.  KIND=0 INDICATES THAT THE WEIGHT FUNCTION IS OF A TYPE NOT
!       LISTED IN THE TABLE BELOW OF CLASSICAL WEIGHT FUNCTIONS.  FOR KIND=0
!       THE USER MUST SUPPLY THE JACOBI MATRIX AND ANY MOMENTS WHICH ARE
!       REQUIRED.
! ALPHA
! BETA
! A
! B     THE WEIGHT FUNCTION AND/OR INTERVAL PARAMETERS.  ANY ONE MAY BE
!       REPLACED BY A DUMMY VARIABLE IF THE WEIGHT FUNCTION IS
!       INDEPENDENT OF IT.
! NWF   AN INTEGER SPECIFYING THE DIMENSION OF THE WORKFIELD WF.
!       MINIMUM VALUES FOR NWF ARE GIVEN IN THE DESCRIPTION OF EACH
!       ROUTINE THAT USES A WORKFIELD.
! WF    FLOATING POINT WORKFIELD ARRAY TO BE SUPPLIED BY THE USER.
! NIWF  AN INTEGER SPECIFYING THE DIMENSION OF IWF
! IWF   INTEGER TYPE WORKFIELD ARRAY TO BE SUPPLIED BY THE USER.
! QFSUM VARIABLE RETURNING THE VALUE OF THE QUADRATURE SUM.
! F     A USER SUPPLIED FUNCTION INVOKED BY A STATEMENT LIKE Y=F(X,I).
!       IT RETURNS THE VALUE OF THE I-TH DERIVATIVE OF F AT X (ZERO-TH
!       DERIVATIVE=FUNCTION).  THE FUNCTION SHOULD BE CAPABLE OF RETURNING
!       DERIVATIVES OF ALL ORDERS UP TO MMAX-1 WHERE MMAX IS THE MAXIMUM
!       MULTIPLICITY USED AT THE KNOTS.  THE ACTUAL PARAMETER USED IN THE CALL
!       TO ROUTINE EIQF, EIQFS, CEIQF AND CEIQFS MUST BE DECLARED IN AN
!       EXTERNAL STATEMENT IN THE CALLING PROGRAM
! LO    INTEGER VARIABLE USED TO CONTROL OUTPUT. IF LO IS SET TO ZERO
!       THEN THERE WILL BE NO OUTPUT PRINTED.  IF LO IS NON-ZERO THEN
!       ABS(LO) WILL BE THE LOGICAL UNIT NUMBER TO WHICH ALL OUTPUT
!       IS DIRECTED.  WHEN LO IS NEGATIVE WEIGHTS ONLY WILL BE PRINTED
!       AND WHEN LO IS POSITIVE THE WEIGHTS AND A CHECK OF THE MOMENTS
!       WILL BE PRINTED.  IN SOME ROUTINES LO.EQ.0 WILL CAUSE A MOMENTS
!       CHECK TO BE COMPUTED EVEN THOUGH THERE IS NO PRINT WHILE IN
!       OTHERS LO.EQ.0 WILL CAUSE ONLY THE WEIGHTS TO BE COMPUTED.  SEE
!       INDIVIDUAL ROUTINES FOR DETAILS.

!  THROUGHOUT THE COMMENTS IN THIS PACKAGE
!       N...IS THE NUMBER OF KNOTS COUNTED ACCORDING TO THEIR MULTIPLICITIES,
!       MMAX...MAXIMUM OF THE MLT(J)
!       RMAX...MAXIMUM OF 2*MMAX AND N+1
!       NSTAR...INTEGER PART OF (N+1)/2

!  ERROR CONDITIONS ARE INDICATED BY THE VARIABLE IER BEING
!  RETURNED WITH A NON-ZERO VALUE.

!  IER =   1       ALPHA > -1 FALSE
!          2       FOR KIND < 8 BETA > -1 IS FALSE
!          3       FOR KIND = 8 NEED BETA < (ALPHA+BETA+2*N) < 0
!                      TO COMPUTE N ELEMENTS OF THE JACOBI MATRIX.
!          4       UNKNOWN WEIGHT FUNCTION. CANNOT GENERATE JACOBI MATRIX
!          5       GAMMA FUNCTION AND MACHINE PARAMETERS ARE NOT
!                      MATCHED IN ACCURACY
!          6       ZERO LENGTH INTERVAL (KIND=1,2,3,4,7)
!          7       B <= 0 FOR KIND=5,6
!          8       A+B <= 0 FOR KIND=8
!          9       NOT ENOUGH  INTEGER WORKFIELD.  NIWF=2*NT WILL DO
!         10       DIMENSION OF WEIGHTS ARRAY TOO SMALL
!         11       JACOBI MATRIX NOT DIAGONALIZED SUCCESSFULLY
!         12       SIZE OF JACOBI MATRIX TOO SMALL FOR NUMBER OF WEIGHTS
!         13       ZERO-TH MOMENT OF WEIGHTS FUNCTION IS NOT > 0
!         14       KNOTS NOT DISTINCT
!         15       SOME KNOT HAS MULTIPLICITY < 1
!         16       POINTERS FOR WGHTS ARRAY CONTRADICTORY
!         17       0 < ABS(KEY) < 5 FALSE (SEE CAWIQ OR EIQF)
!         18       NUMBER OF KNOTS < 1
!       -K, K > 0  AT LEAST K LOCATIONS ARE REQUIRED IN THE FLOATING-POINT
!                  WORKFIELD IN ORDER TO COMPLETE THE CURRENT TASK.

!             SUBROUTINES AND THEIR CALL SEQUENCES

!   CALL CEGQFS(NT, KIND, ALPHA, BETA, F, QFSUM, IER)
!   CALL CEGQF(NT, KIND, ALPHA, BETA, A, B, F, QFSUM, IER)
!   CALL CGQF(NT, T, WTS, KIND, ALPHA, BETA, A, B, LO, IER)
!   CALL CGQFS(NT, T, WTS, KIND, ALPHA, BETA, LO, IER)
!   CALL CDGQF(NT, T, WTS, KIND, ALPHA, BETA, IER)
!   CALL SGQF(NT, T, WTS, AJ, BJ, ZEMU, IER)
!   CALL CLIQFS(NT, T, WTS, KIND, ALPHA, BETA, LO, NWF, WF, NIWF, IWF, IER)
!   CALL CLIQF(NT, T, WTS, KIND, ALPHA, BETA, A, B, LO, NWF, WF, NIWF, IWF,
!              IER)
!   CALL CEIQFS(NT, T, MLT, KIND, ALPHA, BETA, F, QFSUM, NWF, WF, NIWF, IWF,
!               IER)
!   CALL CEIQF(NT, T, MLT, KIND, ALPHA, BETA, A, B, F, QFSUM, NWF, WF, NIWF,
!              IWF, IER)
!   CALL CIQFS(NT, T, MLT, NWTS, WTS, NDX, KEY, KIND, ALPHA, BETA, LO, NWF,
!              WF, IER)
!   CALL CIQF(NT, T, MLT, NWTS, WTS, NDX, KEY, KIND, ALPHA, BETA, A, B, LO,
!             NWF, WF, IER)
!   CALL EIQF(NT, T, MLT, WTS, NDX, KEY, F, QFSUM, IER)
!   CALL EIQFS(NT, T, WTS, F, QFSUM, IER)
!   CALL CAWIQ(NT, T, MLT, NWTS, WTS, NDX, KEY, NST, AJ, BJ, JDF, ZEMU, NWF,
!              WF, IER)
!   CALL CWIQD(M, NM, L, V, XK, NSTAR, PHI, A, WF, Y, R, Z, D)
!   CALL CLASS(KIND, M, ALPHA, BETA, BJ, AJ, ZEMU, IER)
!   CALL WM(W, M, KIND, ALPHA, BETA, IER)
!   CALL PARCHK(KIND, M, ALPHA, BETA, IER)
!   CALL CHKQFS(T, WTS, MLT, NT, NDX, KEY, W, MOP, MEX, KIND, ALPHA,
!               BETA, LO, E, ER, QM, IER)
!   CALL CHKQF(T, WTS, MLT, NT, NDX, KEY, WF, MOP, MEX, KIND, ALPHA,
!              BETA, LO, E, ER, QM, NWF, A, B, IER)
!   CALL SCT(NT, T, ST, KIND, A, B, IER)
!   CALL SCQF(NT, T, MLT, WTS, NDX, SWTS, ST, KIND, ALPHA, BETA, A, B, IER)
!   CALL SCMM(W, M, KIND, ALPHA, BETA, A, B, IER)
!   CALL WTFN(T, W, NT, KIND, ALPHA, BETA, IER)
!   CALL IMTQLX(N, D, E, Z, IER)

!----------------------------------------------------------------------

!   IN THE DESCRIPTIONS OF THE ROUTINES BELOW ALL THE INPUT AND OUTPUT
!   PARAMETERS ARE INDICATED BY THE SINGLE LETTER I OR O ALIGNED TO EACH
!   VARIABLE IN THE CALLING SEQUENCE.  A * INDICATES THAT THE VARIABLE IS
!   SOMETIMES SET ON INPUT AND SOMETIMES SET BY THE ROUTINE.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /ctrlr/ prec(10)
REAL (dp), SAVE  :: prec(1)

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp


CONTAINS


SUBROUTINE cegqf(nt, kind, alpha, beta, a, b, f, qfsum, ier)

! N.B. Arguments NWF, WF, NIWF & IWF have been removed.

!     ROUTINE TO:
!     1.    COMPUTE ALL THE KNOTS AND WEIGHTS OF CLASSICAL WEIGHT FUNCTION
!           GAUSS QUADRATURE FORMULA WITH ALL SIMPLE KNOTS FOR ANY VALID
!           VALUES OF A AND B
!     2.    EVALUATE THE QUADRATURE SUM

!     INPUT AND OUTPUT VARIABLES -

!                       I  I    I     I    I I I O
!      SUBROUTINE CEGQF(NT,KIND,ALPHA,BETA,A,B,F,QFSUM
!     1,NWF,WF,NIWF,IWF,IER)
!       I   O  I    O   O

!     THE USER SUPPLIES A FUNCTION F, WHICH MUST BE DECLARED IN AN EXTERNAL
!     STATEMENT IN THE CALLING PROGRAM, AND WHICH RETURNS VALUES OF F.

!     NEED NWF >= 2*NT
!         NIWF >= 2*NT

!     FUNCTIONS AND SUBROUTINES REFERENCED - CGQF EIQFS F

INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp)  :: wf(2*nt)
INTEGER    :: lu, na, nb

ier = 0

!     SET WORKFIELD FOR WEIGHTS AND KNOTS
lu = 0
na = 1
nb = na + nt
CALL cgqf(nt, wf(nb:), wf(na:), kind, alpha, beta, a, b, lu, ier)
IF (ier /= 0) RETURN

!     EVALUATE THE QUADRATURE SUM
CALL eiqfs(nt, wf(nb:), wf(na:), f, qfsum, ier)

RETURN
END SUBROUTINE cegqf



SUBROUTINE cegqfs(nt, kind, alpha, beta, f, qfsum, ier)

! N.B. Arguments NWF, WF, NIWF & IWF have been removed.

!  ROUTINE TO:
!  1.    COMPUTE ALL THE KNOTS AND WEIGHTS OF CLASSICAL WEIGHT
!        FUNCTION GAUSS QUADRATURE FORMULA WITH ALL SIMPLE KNOTS
!        FOR THE DEFAULT VALUES OF A AND B
!  2.    EVALUATE THE QUADRATURE SUM

!  INPUT AND OUTPUT VARIABLES -

!                     I  I    I     I    I O
!   SUBROUTINE CEGQFS(NT,KIND,ALPHA,BETA,F,QFSUM
! 1,NWF,WF,NIWF,IWF,IER)
!   I   O  I    O   O

!  F MUST BE DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.

!  NEED NWF >= 2*NT
!      NIWF >= 2*NT

!  FUNCTIONS AND SUBROUTINES REFERENCED - CGQFS EIQFS F

INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp)  :: wf(2*nt)
INTEGER    :: lu, na, nb

ier = 0

!     ASSIGN WORKSPACE FOR KNOTS AND WEIGHTS
lu = 0
na = 1
nb = na + nt
CALL cgqfs(nt, wf(nb:), wf(na:), kind, alpha, beta, lu, ier)
IF (ier /= 0) RETURN

!     EVALUATE THE QUADRATURE SUM
CALL eiqfs(nt, wf(nb:), wf(na:), f, qfsum, ier)

RETURN
END SUBROUTINE cegqfs



SUBROUTINE cgqf(nt, t, wts, kind, alpha, beta, a, b, lo, ier)

! N.B. Arguments NWF, WF, NIWF & IWF have been removed.

!     ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF A GAUSS QF WITH
!     1. A CLASSICAL WEIGHT FUNCTION WITH ANY VALID A,B
!     2. ONLY SIMPLE KNOTS
!     3. OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO > 0...COMPUTE AND PRINT KNOTS AND WEIGHTS. PRINT MOMENTS CHECK
!     LO .EQ. 0...COMPUTE KNOTS AND WEIGHTS. PRINT NOTHING
!     LO < 0...COMPUTE AND PRINT KNOTS AND WEIGHTS. NO MOMENTS CHECK.

!     INPUT AND OUTPUT VARIABLES -
!                     I  O O   I    I     I    I I I
!     SUBROUTINE CGQF(NT,T,WTS,KIND,ALPHA,BETA,A,B,LO,
!     I   O  I    O   O
!    1NWF,WF,NIWF,IWF,IER)

!     NEED NWF>= (9*NT+13) IF LO > 0
!                (2*NT)    IF LO .EQ. 0
!                (3*NT+4)  IF LO < 0
!     IWF...DIMENSION MUST BE >= 2*NT

!     USE ROUTINE EIQFS TO EVALUATE A QUADRATURE COMPUTED BY CGQF.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CDGQF CHKQF SCQF

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(OUT)  :: t(:)
REAL (dp), INTENT(OUT)  :: wts(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(IN)     :: lo
INTEGER, INTENT(OUT)    :: ier

REAL (dp), ALLOCATABLE  :: wf(:)
INTEGER  :: iwf(2*nt)
INTEGER  :: i, key, lex, m, mex, mmex, mop, nai, nbi, NE, ner, nqm, nw

!     CHECK THERE IS ENOUGH WORKFIELD AND ASSIGN WORKFIELD
ier = 0
key = 1
mop = 2 * nt
m = mop + 1
mex = m + 2
mmex = MAX(mex,1)
lex = mop
IF (lo /= 0) lex = mex + nt + 1
IF (lo <= 0) mex = 0
NE = 1
ner = NE + mex
nqm = ner + mex
nw = nqm + mex
lex = lex + 3 * mex
IF (ALLOCATED(wf)) DEALLOCATE(wf)
ALLOCATE( wf(lex) )

!     COMPUTE THE GAUSS QF FOR DEFAULT VALUES OF A,B
CALL cdgqf(nt, t, wts, kind, alpha, beta, ier)

!     EXIT IF ERROR
IF (ier /= 0) RETURN

!     PREPARE TO SCALE QF TO OTHER WEIGHT FUNCTION WITH VALID A,B
!     SET UP INTEGER WORK FIELDS
nai = 1
nbi = nai + nt
DO  i = 1, nt
  iwf(nai+i-1) = 1
  iwf(nbi+i-1) = i
END DO

!     IWF(NAI) IS THE MLT ARRAY.  ALL KNOTS MULT=1
!     IWF(NBI) IS THE NDX ARRAY.  NDX(I)=I
!     SCALE THE QUADRATURE
CALL scqf(nt, t, iwf(nai:), wts, iwf(nbi:), wts, t, kind, alpha, beta,  &
          a, b, ier)

!     EXIT IF ERROR OR IF NO PRINT REQUIRED
IF (ier /= 0.OR.lo == 0) RETURN

CALL chkqf(t, wts, iwf(nai:), nt, iwf(nbi:), key, wf(nw:), mop, mmex, kind, &
           alpha, beta, lo, wf(NE:), wf(ner:), wf(nqm:), lex-nw, a, b, ier)

RETURN
END SUBROUTINE cgqf



SUBROUTINE cgqfs(nt, t, wts, kind, alpha, beta, lo, ier)

! N.B. Arguments NWF, WF, NIWF & IWF have been removed.

!     ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF A GAUSS QF WITH
!     1. A CLASSICAL WEIGHT FUNCTION WITH DEFAULT VALUES FOR A,B
!     2. ONLY SIMPLE KNOTS
!     3. OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO>0...COMPUTE AND PRINT KNOTS AND WEIGHTS. PRINT MOMENTS CHECK
!     LO.EQ.0...COMPUTE KNOTS AND WEIGHTS. PRINT NOTHING
!     LO<0...COMPUTE AND PRINT KNOTS AND WEIGHTS. NO MOMENTS CHECK.

!      INPUT AND OUTPUT VARIABLES -
!                       I  O O   I    I     I    I
!      SUBROUTINE CGQFS(NT,T,WTS,KIND,ALPHA,BETA,LO,
!     1NWF,WF,NIWF,IWF,IER)
!      I   O   I   O   O

!     NEED NWF>= (9*NT+13) IF LO > 0
!                (2*NT)    IF LO .EQ. 0
!                (3*NT+4)  IF LO < 0
!     IWF...DIMENSION MUST BE >= 2*NT

!     USE ROUTINE EIQFS TO EVALUATE A QUADRATURE COMPUTED BY CGQFS.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CDGQF CHKQFS

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(OUT)  :: t(:)
REAL (dp), INTENT(OUT)  :: wts(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(IN)     :: lo
INTEGER, INTENT(OUT)    :: ier

INTEGER  :: iwf(2*nt)
INTEGER  :: i, key, lex, m, mex, mmex, mop, nai, nbi, NE, ner, nqm, nw
REAL (dp), ALLOCATABLE  :: wf(:)

!     CHECK THERE IS ENOUGH WORKFIELD AND ASSIGN WORKFIELD
ier = 0
key = 1
mop = 2 * nt
m = mop + 1
mex = m + 2
mmex = MAX(mex,1)
lex = mop
IF (lo /= 0) lex = mex + nt + 1
IF (lo <= 0) mex = 0
NE = 1
ner = NE + mex
nqm = ner + mex
nw = nqm + mex
lex = lex + 3 * mex
IF (ALLOCATED(wf)) DEALLOCATE(wf)
ALLOCATE( wf(lex) )

!     COMPUTE THE GAUSS QF
CALL cdgqf(nt, t, wts, kind, alpha, beta, ier)

!     EXIT IF ERROR OR IF NO PRINT REQUIRED
IF (ier /= 0 .OR. lo == 0) RETURN

!     SET UP INTEGER WORK FIELDS
nai = 1
nbi = nai + nt
DO  i = 1, nt
  iwf(nai+i-1) = 1
  iwf(nbi+i-1) = i
END DO

!     IWF(NAI) IS THE MLT ARRAY. ALL KNOTS MULT=1
!     IWF(NBI) IS THE NDX ARRAY. NDX(I)=I

CALL chkqfs(t, wts, iwf(nai:), nt, iwf(nbi:), key, wf(nw:), mop, mmex,  &
            kind, alpha, beta, lo, wf(NE:), wf(ner:), wf(nqm:), ier)
RETURN
END SUBROUTINE cgqfs



SUBROUTINE cdgqf(nt, t, wts, kind, alpha, beta, ier)

! N.B. Arguments NWF & WF have been removed.

!  ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF A GAUSS QF WITH
!  1. A CLASSICAL WEIGHT FUNCTION WITH DEFAULT VALUES FOR A,B
!  2. ONLY SIMPLE KNOTS
!     NO MOMENTS CHECK OR PRINTING DONE.

!   INPUT AND OUTPUT VARIABLES -
!                    I  O O   I    I     I
!   SUBROUTINE CDGQF(NT,T,WTS,KIND,ALPHA,BETA,
!  1NWF,WF,IER)
!   I   O  O

!  NWF... MUST BE >= 2*NT

!  USE ROUTINE EIQFS TO EVALUATE A QUADRATURE COMPUTED BY CGQFS.

!  FUNCTIONS AND SUBROUTINES REFERENCED - CLASS PARCHK SGQF


INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(OUT)  :: t(:)
REAL (dp), INTENT(OUT)  :: wts(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(OUT)    :: ier

REAL (dp)  :: wf(2*nt)
REAL (dp)  :: zemu
INTEGER    :: na, nb

CALL parchk(kind, 2*nt, alpha, beta, ier)

!     SET UP ARRAYS FOR DIAGONAL AND SUB-DIAGONAL OF JACOBI MATRIX
na = 1
nb = na + nt
IF (ier /= 0) RETURN

!     GET JACOBI MATRIX AND ZERO-TH MOMENT
CALL class(kind, nt, alpha, beta, wf(nb:), wf(na:), zemu, ier)
IF (ier /= 0) RETURN
CALL sgqf(nt, t, wts, wf(na:), wf(nb:), zemu, ier)

RETURN
END SUBROUTINE cdgqf



SUBROUTINE sgqf(nt, t, wts, aj, bj, zemu, ier)
!     ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF A GAUSS QUADRATURE
!     FORMULA (WITH SIMPLE KNOTS) FROM THE JACOBI MATRIX AND THE ZERO-TH
!     MOMENT OF THE WEIGHT FUNCTION, USING THE GOLUB-WELSCH TECHNIQUE

!      INPUT AND OUTPUT VARIABLES -
!                      I  O O   I  I  I    O
!      SUBROUTINE SGQF(NT,T,WTS,AJ,BJ,ZEMU,IER)

!     INPUT PARAMETERS
!     AJ...DIAGONAL OF JACOBI MATRIX
!     BJ...SUB-DIAGONAL OF JACOBI MATRIX ( IN BJ(1)..BJ(NT-1) )
!     ZEMU...ZERO-TH MOMENT OF WEIGHT FUNCTION

!     OUTPUT PARAMETERS
!     AT OUTPUT T AND WTS CONTAIN THE KNOTS AND WEIGHTS
!     THE ARRAY BJ IS OVERWRITTEN DURING EXECUTION

!     FUNCTIONS AND SUBROUTINES REFERENCED - IMTQLX MACHEP SQRT

INTEGER, INTENT(IN)        :: nt
REAL (dp), INTENT(OUT)     :: t(:)
REAL (dp), INTENT(OUT)     :: wts(:)
REAL (dp), INTENT(IN)      :: aj(:)
REAL (dp), INTENT(IN OUT)  :: bj(:)
REAL (dp), INTENT(IN)      :: zemu
INTEGER, INTENT(OUT)       :: ier

INTEGER   :: i

ier = 0

!     COMPUTE MACHINE EPSILON FOR IMTQLX
prec(1) = EPSILON(one)

!     EXIT IF ZERO-TH MOMENT NOT POSITIVE
IF (zemu <= zero) ier = 13
IF (ier /= 0) RETURN

!     SET UP VECTORS FOR IMTQLX
DO  i = 1, nt
  t(i) = aj(i)
  wts(i) = zero
END DO
wts(1) = SQRT(zemu)

!     DIAGONALIZE JACOBI MATRIX
CALL imtqlx(nt, t, bj, wts, ier)

!     CHECK FOR ERROR  RETURN FROM IMTQLX
IF (ier /= 0) THEN
  ier = 11
  RETURN
END IF

wts(1:nt) = wts(1:nt) ** 2

RETURN
END SUBROUTINE sgqf



SUBROUTINE cliqfs(nt, t, wts, kind, alpha, beta, lo, nwf, wf, niwf, iwf, ier)

!     ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF AN INTERPOLATORY
!     QF WITH
!     1. A CLASSICAL WEIGHT FUNCTION WITH DEFAULT VALUES FOR A,B
!     2. ONLY SIMPLE KNOTS
!     3. OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO > 0...COMPUTE WEIGHTS. PRINT WEIGHTS. PRINT MOMENTS CHECK.
!     LO .EQ. 0...COMPUTE WEIGHTS. PRINT NOTHING.
!     LO < 0...COMPUTE WEIGHTS. PRINT WEIGHTS.

!      INPUT AND OUTPUT VARIABLES -
!                        I  I O   I    I     I
!      SUBROUTINE CLIQFS(NT,T,WTS,KIND,ALPHA,BETA,
!     1LO,NWF,WF,NIWF,IWF,IER)
!      I  I   O  I    O   O

!     NEED NWF  >= (5*N+9)/2  IF LO <= 0
!                  (9*N+25)/2 IF LO > 0
!          NIWF >= 2*NT

!     USE ROUTINE EIQFS TO EVALUATE A QUADRATURE COMPUTED BY CLIQFS.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CIQFS

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(OUT)  :: wts(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(IN)     :: lo
INTEGER, INTENT(IN)     :: nwf
REAL (dp), INTENT(OUT)  :: wf(:)
INTEGER, INTENT(IN)     :: niwf
INTEGER, INTENT(OUT)    :: iwf(:)
INTEGER, INTENT(OUT)    :: ier

INTEGER :: key, na, nb

ier = 0
IF (niwf < 2*nt) THEN
  ier = 9
  RETURN
END IF
key = 1

!     SET UP WORKFIELD FOR MLT,NDX
na = 1
nb = na + nt
iwf(1:nt) = 1
CALL ciqfs(nt, t, iwf(na:), nt, wts, iwf(nb:), key, kind, alpha, beta, lo, &
           nwf, wf, ier)
RETURN
END SUBROUTINE cliqfs



SUBROUTINE cliqf(nt, t, wts, kind, alpha, beta, a, b, lo, nwf, wf, niwf,  &
                 iwf, ier)

!     ROUTINE TO COMPUTE ALL THE KNOTS AND WEIGHTS OF AN INTERPOLATORY QF WITH
!     1. ONLY SIMPLE KNOTS AND
!     2. A CLASSICAL WEIGHT FUNCTION WITH ANY VALID A,B
!     3. OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO > 0...COMPUTE WEIGHTS. PRINT WEIGHTS. PRINT MOMENTS CHECK.
!     LO .EQ. 0...COMPUTE WEIGHTS. PRINT NOTHING.
!     LO < 0...COMPUTE WEIGHTS. PRINT WEIGHTS.

!      INPUT AND OUTPUT VARIABLES -
!                       I  I O   I    I     I    I I
!      SUBROUTINE CLIQF(NT,T,WTS,KIND,ALPHA,BETA,A,B,
!     1LO,NWF,WF,NIWF,IWF,IER)
!      I  I   O  I    O   O

!     NEED NWF  >= (5*N+9)/2  IF LO <= 0
!                  (9*N+25)/2 IF LO > 0
!          NIWF >= 2*NT

!     USE ROUTINE EIQFS TO EVALUATE A QUADRATURE COMPUTED BY CLIQF.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CIQF

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(OUT)  :: wts(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(IN)     :: lo
INTEGER, INTENT(IN)     :: nwf
REAL (dp), INTENT(OUT)  :: wf(:)
INTEGER, INTENT(IN)     :: niwf
INTEGER, INTENT(OUT)    :: iwf(:)
INTEGER, INTENT(OUT)    :: ier

INTEGER :: key, na, nb

ier = 0
IF (niwf < 2*nt) THEN
  ier = 9
  RETURN
END IF
key = 1

!     SET UP WORKFIELD FOR MLT,NDX
na = 1
nb = na + nt
iwf(1:nt) = 1
CALL ciqf(nt, t, iwf(na:), nt, wts, iwf(nb:), key, kind, alpha, beta, a, b, &
          lo, nwf, wf, ier)
RETURN
END SUBROUTINE cliqf



SUBROUTINE ceiqfs(nt, t, mlt, kind, alpha, beta, f, qfsum, nwf, wf, niwf,  &
                  iwf, ier)
!  ROUTINE TO:
!  1.    COMPUTE AN INTERPOLATORY QF FOR CLASSICAL
!        WEIGHT FUNCTION WITH DEFAULT VALUES FOR A,B
!  2.    EVALUATE THE QUADRATURE SUM

!   INPUT AND OUTPUT VARIABLES -
!                     I  I I   I    I     I    I O
!   SUBROUTINE CEIQFS(NT,T,MLT,KIND,ALPHA,BETA,F,QFSUM
!  1,NWF,WF,NIWF,IWF,IER)
!    I   O  I    O   O

!  NEED NWF >= NSTAR+RMAX+NT+3*(N+1)
!      NIWF >= NT

!  FUNCTION F, MUST BE DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.

!  FUNCTIONS AND SUBROUTINES REFERENCED - CIQFS EIQF F

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
INTEGER, INTENT(IN)     :: mlt(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(IN)     :: nwf
REAL (dp), INTENT(OUT)  :: wf(:)
INTEGER, INTENT(IN)     :: niwf
INTEGER, INTENT(OUT)    :: iwf(:)
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER :: j, key, l, lex, lu, m, mtm, n, na, nst, nw

ier = 0
IF (niwf < nt) THEN
  ier = 9
  RETURN
END IF
lu = 0
n = 0
mtm = mlt(1)
DO  j = 1, nt
  mtm = MAX(mtm,mlt(j))
  n = n + mlt(j)
END DO
m = n + 1
nst = m / 2
l = MIN(2*mtm,m)
lex = nst + 3 * m + l + nt
IF (nwf < lex) THEN
  ier = -lex
  RETURN
END IF

!     INDICES FOR WTS,NDX,WF (RESP)
na = 1
nw = na + n
key = 1
CALL ciqfs(nt, t, mlt, n, wf(na:), iwf, key, kind, alpha, beta, lu, nwf-nw,  &
           wf(nw:), ier)
IF (ier /= 0) RETURN
CALL eiqf(nt, t, mlt, wf(na:), iwf, key, f, qfsum, ier)
RETURN
END SUBROUTINE ceiqfs



SUBROUTINE ceiqf(nt, t, mlt, kind, alpha, beta, a, b, f, qfsum, nwf, wf,  &
                 niwf, iwf, ier)
!  ROUTINE TO:
!  1.    COMPUTE AN INTERPOLATORY QF WITH CLASSICAL
!        WEIGHT FUNCTION WITH ANY VALID A,B
!  2.    EVALUATE THE QUADRATURE SUM

!   INPUT AND OUTPUT VARIABLES -
!                    I  I I   I    I     I    I I I O
!   SUBROUTINE CEIQF(NT,T,MLT,KIND,ALPHA,BETA,A,B,F,QFSUM
!  1,NWF,WF,NIWF,IWF,IER)
!    I   O  I    O   O

!  NEED NWF  >= NSTAR+RMAX+NT+3*(N+1)
!       NIWF >= NT

!  FUNCTION F, MUST BE DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.

!  FUNCTIONS AND SUBROUTINES REFERENCED - CIQF EIQF F

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
INTEGER, INTENT(IN)     :: mlt(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(IN)     :: nwf
REAL (dp), INTENT(OUT)  :: wf(:)
INTEGER, INTENT(IN)     :: niwf
INTEGER, INTENT(OUT)    :: iwf(:)
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER :: j, key, l, lex, lu, m, mtm, n, na, nst, nw

ier = 0
IF (niwf < nt) THEN
  ier = 9
  RETURN
END IF
lu = 0
n = 0
mtm = mlt(1)
DO  j = 1, nt
  mtm = MAX(mtm,mlt(j))
  n = n + mlt(j)
END DO
m = n + 1
nst = m / 2
l = MIN(2*mtm,m)
lex = nst + 3 * m + l + nt
IF (nwf < lex) THEN
  ier = -lex
  RETURN
END IF

!     INDICES FOR WTS,WF
na = 1
nw = na + n
key = 1
CALL ciqf(nt, t, mlt, n, wf(na:), iwf, key, kind, alpha, beta, a, b, lu,  &
          nwf-nw, wf(nw:), ier)
IF (ier /= 0) RETURN
CALL eiqf(nt, t, mlt, wf(na:), iwf, key, f, qfsum, ier)
RETURN
END SUBROUTINE ceiqf



SUBROUTINE ciqfs(nt, t, mlt, nwts, wts, ndx, key, kind, alpha, beta, lo, nwf, &
                 wf, ier)
!     ROUTINE TO COMPUTE SOME OR ALL THE WEIGHTS OF A
!     QF FOR A CLASSICAL WEIGHT FUNCTION WITH DEFAULT VALUES OF A,B
!     AND A GIVEN SET OF KNOTS AND MULTIPLICITIES. THE WEIGHTS MAY BE
!     PACKED INTO THE OUTPUT ARRAY WTS ACCORDING TO A USER-DEFINED
!     PATTERN OR SEQUENTIALLY. THE ROUTINE WILL ALSO
!     OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO > 0...COMPUTE WEIGHTS. PRINT WEIGHTS. PRINT MOMENTS CHECK.
!     LO .EQ. 0...COMPUTE WEIGHTS. PRINT NOTHING. NO MOMENTS CHECK.
!     LO < 0...COMPUTE WEIGHTS. PRINT WEIGHTS. NO MOMENTS CHECK.

!      INPUT AND OUTPUT VARIABLES -
!                       I  I I   I    O   *   I   I    I     I
!      SUBROUTINE CIQFS(NT,T,MLT,NWTS,WTS,NDX,KEY,KIND,ALPHA,BETA
!     1,LO,NWF,WF,IER)
!       I  I   O  O

!     NEED NWF >= NSTAR + RMAX + 2*(N+1) IF NO MOMENTS CHECK REQUIRED
!                 NSTAR + RMAX + 2*(2*N+5) IF A MOMENTS CHECK REQUIRED
!     KEY...AN INTEGER VARIABLE INDICATING THE STRUCTURE OF THE WTS
!           ARRAY. IT WILL NORMALLY BE SET TO 1.  FOR MORE DETAILS SEE
!           THE COMMENTS IN CAWIQ.
!     NDX...AN INTEGER ARRAY OF DIMENSION NT USED TO INDEX THE OUTPUT
!           ARRAY WTS.  IF KEY=1 NDX NEED NOT BE PRESET.  FOR MORE
!           DETAILS SEE THE COMMENTS IN CAWIQ.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CAWIQ CHKQFS CLASS

INTEGER, INTENT(IN)      :: nt
REAL (dp), INTENT(IN)    :: t(:)
INTEGER, INTENT(IN)      :: mlt(:)
INTEGER, INTENT(IN)      :: nwts
REAL (dp), INTENT(OUT)   :: wts(:)
INTEGER, INTENT(IN OUT)  :: ndx(:)
INTEGER, INTENT(IN)      :: key
INTEGER, INTENT(IN)      :: kind
REAL (dp), INTENT(IN)    :: alpha
REAL (dp), INTENT(IN)    :: beta
INTEGER, INTENT(IN)      :: lo
INTEGER, INTENT(IN)      :: nwf
REAL (dp), INTENT(OUT)   :: wf(:)
INTEGER, INTENT(OUT)     :: ier

REAL (dp) :: zemu
INTEGER   :: ifl, j, jdf, k, l, lex, m, mex, mmex
INTEGER   :: mtm, n, na, nb, nd, NE, nf, nst, nw

ier = 0
jdf = 0
n = 0
mtm = mlt(1)
l = ABS(key)
DO  j = 1, nt
  IF (l /= 1) THEN
    IF (ABS(ndx(j)) == 0) CYCLE
  END IF
  mtm = MAX(mtm,mlt(j))
  n = n + mlt(j)
END DO

!     N KNOTS WHEN COUNTED ACCORDING TO MULT
IF (nwts < n) THEN
  ier = 10
  RETURN
END IF
m = n + 1
mex = 2 + m
nst = m / 2
ifl = 1
IF (lo <= 0) ifl = 0
l = MIN(2*mtm,m)
k = MAX(m,3*(mex)*ifl)
lex = nst + m + l + k
IF (nwf < lex) THEN
  ier = -lex
  RETURN
END IF

!     SET UP WORK FIELD INDICES FOR CLASS AND CAWIQ
na = 1
nb = na + nst
nw = nb + nst

!     GET JACOBI MATRIX
CALL class(kind, nst, alpha, beta, wf(nb:), wf(na:), zemu, ier)
IF (ier /= 0) RETURN

!     CALL WEIGHTS ROUTINE
CALL cawiq(nt, t, mlt, n, wts, ndx, key, nst, wf(na:), wf(nb:), jdf, zemu,  &
           nwf-nw, wf(nw:), ier)

!     RETURN IF NO PRINTING OR CHECKING REQUIRED
IF (ier /= 0 .OR. lo == 0) RETURN
mmex = mex * ifl
nd = 1
NE = nd + mmex
nf = NE + mmex
nw = nf + mmex

!     CALL CHECKING ROUTINE
CALL chkqfs(t, wts, mlt, nt, ndx, key, wf(nw:), m-1, mex, kind, alpha,  &
            beta, lo, wf(nd:), wf(NE:), wf(nf:), ier)
RETURN
END SUBROUTINE ciqfs



SUBROUTINE ciqf(nt, t, mlt, nwts, wts, ndx, key, kind, alpha, beta, a, b, lo, &
                nwf, wf, ier)
!     ROUTINE TO COMPUTE SOME OR ALL THE WEIGHTS OF A
!     QF FOR A CLASSICAL WEIGHT FUNCTION WITH ANY VALID A,B AND
!     A GIVEN SET OF KNOTS AND MULTIPLICITIES. THE WEIGHTS MAY BE
!     PACKED INTO THE OUTPUT ARRAY WTS ACCORDING TO A USER-DEFINED
!     PATTERN OR SEQUENTIALLY. THE ROUTINE WILL ALSO
!     OPTIONALLY PRINT KNOTS AND WEIGHTS AND A CHECK OF THE MOMENTS

!     LO>0...COMPUTE WEIGHTS. PRINT WEIGHTS. PRINT MOMENTS CHECK.
!     LO.EQ.0...COMPUTE WEIGHTS. PRINT NOTHING. NO MOMENTS CHECK.
!     LO<0...COMPUTE WEIGHTS. PRINT WEIGHTS. NO MOMENTS CHECK.

!      INPUT AND OUTPUT VARIABLES -
!                      I  I I   I    O   *   I   I    I     I    I I
!      SUBROUTINE CIQF(NT,T,MLT,NWTS,WTS,NDX,KEY,KIND,ALPHA,BETA,A,B
!     1,LO,NWF,WF,IER)
!       I  I   O  O

!     NEED NWF>=NSTAR+RMAX+2*(N+1) IF NO MOMENTS CHECK REQUIRED
!                 NSTAR+RMAX+5*N+NT+13 IF A MOMENTS CHECK IS REQUIRED
!     KEY...AN INTEGER VARIABLE INDICATING THE STRUCTURE OF THE WTS
!           ARRAY. IT WILL NORMALLY BE SET TO 1. THIS WILL CAUSE THE
!           WEIGHTS TO BE PACKED SEQUENTIALLY IN ARRAY WTS.
!           FOR MORE DETAILS SEE THE COMMENTS IN CAWIQ.
!     NDX...AN INTEGER ARRAY OF DIMENSION NT USED TO INDEX THE OUTPUT
!           ARRAY WTS. IF KEY=1 NDX NEED NOT BE PRESET. FOR MORE
!           DETAILS SEE THE COMMENTS IN CAWIQ.

!     FUNCTIONS AND SUBROUTINES REFERENCED - CHKQF CIQFS SCQF SCT

INTEGER, INTENT(IN)      :: nt
REAL (dp), INTENT(IN)    :: t(:)
INTEGER, INTENT(IN)      :: mlt(:)
INTEGER, INTENT(IN)      :: nwts
REAL (dp), INTENT(OUT)   :: wts(:)
INTEGER, INTENT(IN OUT)  :: ndx(:)
INTEGER, INTENT(IN)      :: key
INTEGER, INTENT(IN)      :: kind
REAL (dp), INTENT(IN)    :: alpha
REAL (dp), INTENT(IN)    :: beta
REAL (dp), INTENT(IN)    :: a
REAL (dp), INTENT(IN)    :: b
INTEGER, INTENT(IN)      :: lo
INTEGER, INTENT(IN)      :: nwf
REAL (dp), INTENT(OUT)   :: wf(:)
INTEGER, INTENT(OUT)     :: ier

INTEGER :: ifl, j, k, l, lex, lu, m, mex, mmex, mtm, nd, NE, nf, nst, nw

ier = 0
m = 1
mtm = 1
l = ABS(key)
DO  j = 1, nt
  IF (l /= 1) THEN
    IF (ABS(ndx(j)) == 0) CYCLE
  END IF
  mtm = MAX(mtm,mlt(j))
  m = m + mlt(j)
END DO
IF (nwts+1 < m) THEN
  ier = 10
  RETURN
END IF
nst = m / 2
mex = 2 + m
ifl = 1
IF (lo <= 0) ifl = 0
l = MIN(2*mtm,m)
k = MAX(m,3*mex*ifl)
lex = nst + l + k + m + (mex+nt) * ifl
IF (nwf < lex) THEN
  ier = -lex
  RETURN
END IF
nst = 1
nf = nst + nt

!     SCALE THE KNOTS TO DEFAULT A,B
CALL sct(nt, t, wf(nst:), kind, a, b, ier)
IF (ier /= 0) RETURN
lu = 0
CALL ciqfs(nt, wf(nst:), mlt, nwts, wts, ndx, key, kind, alpha, beta, lu, &
           nwf-nf+1, wf(nf:), ier)
IF (ier /= 0) RETURN

!     DON'T SCALE USER'S KNOTS - ONLY SCALE WEIGHTS
CALL scqf(nt, wf(nst:), mlt, wts, ndx, wts, wf(nst:), kind, alpha,  &
          beta, a, b, ier)
IF (ier /= 0.OR.lo == 0) RETURN
mmex = mex * ifl
nd = 1
NE = nd + mmex
nf = NE + mmex
nw = nf + mmex
CALL chkqf(t, wts, mlt, nt, ndx, key, wf(nw:), m-1, mex, kind, alpha, &
           beta, lo, wf(nd:), wf(NE:), wf(nf:), nwf-nw, a, b, ier)
RETURN
END SUBROUTINE ciqf



SUBROUTINE eiqf(nt, t, mlt, wts, ndx, key, f, qfsum, ier)

! N.B. Argument NWTS has been removed.

!  ROUTINE TO EVALUATE AN INTERPOLATORY QF WITH KNOTS, WEIGHTS
!  AND INTEGRAND SUPPLIED.
!  ALL KNOTS FOR WHICH NDX(I).NE.0 ARE USED.
!   INPUT AND OUTPUT VARIABLES -
!                   I  I I   I   I    I   I   I O     O
!   SUBROUTINE EIQF(NT,T,MLT,WTS,NWTS,NDX,KEY,F,QFSUM,IER)

!  **************************************************************
!  *
!  *  F.......A FUNCTION WITH 2 PARAMETERS F(X,I)
!  *  TO BE SUPPLIED BY THE USER.  IT MUST RETURN THE I-TH
!  *  DERIVATIVE OF F, THE FUNCTION BEING INTEGRATED, AT X.
!  *  I MUST RANGE FROM 0 (FOR THE FUNCTION VALUES) UP TO
!  *  (THE MAXIMUM VALUE IN MLT)-1. THIS FUNCTION WILL ONLY
!  *  BE CALLED WITH F AND ITS DERIVATIVES AT THE KNOTS SUPPLIED
!  *  SO IT CAN GENERATE THE VALUES FOR F BY TABLE LOOKUP.
!  *  THIS CAN BE ACHIEVED BY REPLACING THE LINE
!  *           QFSUM = QFSUM + WTS(L+I-1)*F(T(J),I-1)/P
!  *  WITH
!  *           QFSUM = QFSUM + WTS(L+I-1)*F(T,J,I-1)/P
!  *  WHERE THE FUNCTION F HAS THE KNOTS ARRAY T AS A PARAMETER
!  *  AND THE REQUIRED KNOT IS INDICATED BY THE INDEX J. F IS
!  *  CALLED ONLY FROM THIS ROUTINE AND EIQFS.
!  *
!  **************************************************************
!  FUNCTIONS AND SUBROUTINES REFERENCED -  F

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
INTEGER, INTENT(IN)     :: mlt(:)
REAL (dp), INTENT(IN)   :: wts(:)
INTEGER, INTENT(IN)     :: ndx(:)
INTEGER, INTENT(IN)     :: key
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp) :: p
INTEGER   :: i, j, l

ier = 0
l = ABS(key)
IF (l < 1.OR.l > 4) THEN
  ier = 17
  RETURN
END IF
qfsum = zero
DO  j = 1, nt
  l = ABS(ndx(j))
  IF (l /= 0) THEN
    p = one
    DO  i = 1, mlt(j)
      qfsum = qfsum + wts(l+i-1) * f(t(j),i-1) / p
      IF (key <= 0) THEN
        p = p * i
      END IF
    END DO
  END IF
END DO
RETURN
END SUBROUTINE eiqf



SUBROUTINE eiqfs(nt, t, wts, f, qfsum, ier)

!     ROUTINE TO EVALUATE AN INTERPOLATORY QF WITH ALL KNOTS SIMPLE AND
!     ALL KNOTS INCLUDED IN THE QUADRATURE. THIS ROUTINE WILL BE USED
!     TYPICALLY AFTER CLIQF OR CLIQFS HAVE BEEN CALLED.

!      INPUT AND OUTPUT VARIABLES -
!                       I  I I   I O     O
!      SUBROUTINE EIQFS(NT,T,WTS,F,QFSUM,IER)

!     **************************************************************
!     *
!     *  F.......A FUNCTION WITH 2 PARAMETERS F(X,I) TO BE SUPPLIED
!     *  BY THE USER.  F MUST RETURN THE VALUE OF F,
!     *  THE FUNCTION BEING INTEGRATED, AT X.
!     *  I MAY BE A DUMMY VARIABLE BUT IS INCLUDED TO MAKE THIS
!     *  DEFINITION CONFORM WITH THAT FOR F ELSEWHERE. THIS FUNCTION
!     *  WILL ONLY BE CALLED WITH F AND ITS DERIVATIVES AT THE KNOTS
!     *  SUPPLIED SO IT CAN GENERATE THE VALUES FOR F BY TABLE LOOKUP.
!     *  THIS CAN BE ACHIEVED BY REPLACING THE LINE
!     *               QFSUM = QFSUM + WTS(J)*F(T(J),0)
!     *  WITH
!     *              QFSUM = QFSUM + WTS(J)*F(T,J,0)
!     *  WHERE THE FUNCTION F HAS THE KNOTS ARRAY T AS A PARAMETER
!     *  AND THE REQUIRED KNOT IS INDICATED BY THE INDEX J. F IS
!     *  CALLED ONLY FROM THIS ROUTINE AND EIQF.
!     *
!     **************************************************************
!     FUNCTIONS AND SUBROUTINES REFERENCED -  F

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(IN)   :: wts(:)
REAL (dp), INTENT(OUT)  :: qfsum
INTEGER, INTENT(OUT)    :: ier

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER  :: j

ier = 0
qfsum = zero
DO  j = 1, nt
  qfsum = qfsum + wts(j) * f(t(j),0)
END DO
RETURN
END SUBROUTINE eiqfs



SUBROUTINE cawiq(nt, t, mlt, nwts, wts, ndx, key, nst, aj, bj, jdf, zemu,  &
                 nwf, wf, ier)

!  THIS ROUTINE, GIVEN A SET OF DISTINCT KNOTS, T, THEIR MULTIPLICITIES MLT,
!  THE JACOBI MATRIX ASSOCIATED WITH THE POLYNOMIALS ORTHOGONAL WITH RESPECT
!  TO THE WEIGHT FUNCTION W(X), AND THE ZERO-TH MOMENT OF W(X) COMPUTES THE
!  WEIGHTS OF THE QUADRATURE

!                                    (I)
!      SUM          SUM         D   F   (T(J))
!    J=1,NT     I=0,MLT(J)-1     J,I

!    WHICH IS TO APPROXIMATE

!    INTEGRAL  F(T)W(T) DT
!     ÕA,Bå

!  THE ROUTINE MAKES VARIOUS CHECKS, AS INDICATED BELOW, SETS UP VARIOUS
!  VECTORS AND, IF NECESSARY, CALLS FOR THE DIAGONALIZATION OF THE JACOBI
!  MATRIX THAT IS ASSOCIATED WITH THE POLYNOMIALS ORTHOGONAL WITH RESPECT
!  TO W(X) ON ÕA,Bå.  THEN FOR EACH KNOT, THE WEIGHTS OF WHICH ARE REQUIRED,
!  IT CALLS THE ROUTINE CWIQD WHICH COMPUTES ALL THE WEIGHTS FOR ANY GIVEN KNOT.

!   INPUT AND OUTPUT VARIABLES -
!                    I  I I   I    O   *   I
!   SUBROUTINE CAWIQ(NT,T,MLT,NWTS,WTS,NDX,KEY
!  1,NST,AJ,BJ,JDF,ZEMU,NWF,WF,IER)
!    I   I  I  I   I    I   O  O

!  PARAMETERS MARKED WITH A * MAY BE CHANGED BY THE SUBROUTINE

! *NDX     THIS ARRAY ASSOCIATES WITH EACH DISTINCT KNOT T(J),
!          AN INTEGER NDX(J) WHICH IS SUCH THAT THE WEIGHT TO THE
!          I-TH DERIV VALUE OF F AT THE J-TH KNOT, IS STORED IN

!                   WTS(ABS(NDX(J))+I) J=1,2,...,NT,
!                                      I=0,1,2,...,MLT(J)-1
!          ALSO:
!          NDX     > 0 MEANS WEIGHTS ARE WANTED FOR THIS KNOT
!                  < 0 MEANS WEIGHTS NOT WANTED FOR THIS KNOT BUT THE
!                            KNOT IS TO BE INCLUDED IN THE QUADRATURE
!                  = 0 MEANS IGNORE THIS KNOT COMPLETELY
!  KEY     SWITCH INDICATING STRUCTURE OF OUTPUT ARRAYS WTS AND NDX.

!          ABS(KEY)=    1      SET UP POINTERS IN NDX FOR ALL KNOTS
!                              IN T ARRAY (ROUTINE CAWIQ DOES THIS).
!                              THE CONTENTS OF NDX ARE NOT TESTED
!                              ON INPUT AND WEIGHTS ARE PACKED
!                              SEQUENTIALLY IN WTS AS INDICATED
!                              ABOVE
!                       2      SET UP POINTERS ONLY FOR KNOTS WHICH
!                              HAVE NDX.NE.0 ON INPUT. ALL KNOTS
!                              WHICH HAVE A NON-ZERO FLAG ARE
!                              ALLOCATED SPACE IN WTS
!                       3      SET UP POINTERS ONLY FOR KNOTS WHICH
!                              HAVE NDX>0 ON INPUT. SPACE IN WTS
!                              ALLOCATED ONLY FOR KNOTS WITH
!                              NDX > 0
!                       4      NDX ASSUMED TO BE PRESET AS POINTER
!                              ARRAY ON INPUT

!           AND KEY>0 FOR WEIGHTS WTS(J) REQUIRED IN STANDARD FORM
!               KEY<0 FOR J]WTS(J) REQUIRED

!  NST     DIMENSION OF JACOBI MATRIX. NST SHOULD BE BETWEEN (N+1)/2
!          AND N. THE USUAL CHOICE WILL BE (N+1)/2
! *JDF     FLAG TO INDICATE WHETHER JACOBI MATRIX NEEDS DIAGONALIZING OR NOT
!                JDF=      0      DIAGONALIZATION REQUIRED
!                          1      DIAGONALIZATION NOT REQUIRED

! *AJ,BJ   IT IS ASSUMED ON INPUT THAT
!          IF JDF = 0 THEN       AJ CONTAINS THE DIAGONAL OF THE JACOBI
!                                MATRIX AND BJ CONTAINS THE SUBDIAGONAL.

!                                NOTE THAT BJ(NST) IS ALSO
!                                DEFINED BUT NOT USED.

!          IF JDF = 1            AJ CONTAINS THE EIGENVALUES OF
!                                THE JACOBI MATRIX AND
!                                BJ CONTAINS THE SQUARES OF THE
!                                ELEMENTS OF THE 1ST ROW OF U THE
!                                ORTHOGONAL MATRIX DIAGONALIZING THE
!                                                       T
!                                JACOBI MATRIX AS  U D U .
!  ZEMU    ZERO-TH MOMENT OF THE WEIGHT FUNCTION W(X)
!  NWF     DIMENSION OF WORK FIELD.  MUST HAVE NWF >= NST + RMAX + N
! *IER     ERROR FLAG: CODED AS FOLLOWS
!          10        NWTS TOO SMALL
!          11        JACOBI MATRIX NOT DIAGONALIZED SUCCESSFULLY
!          12        NST TOO SMALL FOR N
!          13        ZEMU > 0 FALSE
!          14        KNOTS NOT DISTINCT
!          15        MLT(J) < 1 FOR SOME J
!          16        POINTERS FOR WTS CONTRADICTORY
!          17        0 < ABS(KEY) < 5 FALSE
!          18        NT < 1

!  FUNCTIONS AND SUBROUTINES REFERENCED - CWIQD IMTQLX MACHEP SIGN

INTEGER, INTENT(IN)        :: nt
REAL (dp), INTENT(IN)      :: t(:)
INTEGER, INTENT(IN)        :: mlt(:)
INTEGER, INTENT(IN)        :: nwts
REAL (dp), INTENT(OUT)     :: wts(:)
INTEGER, INTENT(IN OUT)    :: ndx(:)
INTEGER, INTENT(IN)        :: key
INTEGER, INTENT(IN)        :: nst
REAL (dp), INTENT(IN OUT)  :: aj(:)
REAL (dp), INTENT(IN OUT)  :: bj(:)
INTEGER, INTENT(IN OUT)    :: jdf
REAL (dp), INTENT(IN)      :: zemu
INTEGER, INTENT(IN)        :: nwf
REAL (dp), INTENT(OUT)     :: wf(:)
INTEGER, INTENT(OUT)       :: ier

REAL (dp) :: p, tmp
INTEGER   :: i, ierrx, ip, ipm, j, jj, jp, k, l, m, mnm, mtj, mtm
INTEGER   :: n, nk, nm, nmax, nr, nw, ny, nz

!     COMPUTE MACHINE EPSILON
prec(1) = EPSILON(one)

!     EXIT IF NT < 1
ier = 18
IF (nt < 1) RETURN
ier = 0

!     CHECK FOR INDISTINCT KNOTS
IF (nt /= 1) THEN
  k = nt - 1
  DO  i = 1, k
    tmp = t(i)
    l = i + 1
    DO  j = l, nt
      IF (ABS(tmp-t(j)) <= prec(1)) THEN
        ier = 14
        RETURN
      END IF
    END DO
  END DO
END IF

!     CHECK MULTIPLICITIES,
!     SET UP VARIOUS USEFUL PARAMETERS AND
!     SET UP OR CHECK POINTERS TO WTS ARRAY
l = ABS(key)
IF (l < 1 .OR. l > 4) THEN
  ier = 17
  RETURN
END IF
k = 1
SELECT CASE ( l )
  CASE (    1)
    DO  i = 1, nt
      ndx(i) = k
      IF (mlt(i) < 1) GO TO 70
      k = k + mlt(i)
    END DO
    n = k - 1
    GO TO 120

  CASE (  2:3)
    n = 0
    DO  i = 1, nt
      IF (ndx(i) /= 0) THEN
        IF (mlt(i) < 1) GO TO 70
        n = n + mlt(i)
        IF (ndx(i) >= 0 .OR. l /= 3) THEN
          ndx(i) = SIGN(k, ndx(i))
          k = k + mlt(i)
        END IF
      END IF
    END DO
    IF (k <= nwts+1) GO TO 120
    ier = 10
    RETURN

  CASE (    4)
    DO  i = 1, nt
      ip = ABS(ndx(i))
      IF (ip /= 0) THEN
        ipm = ip + mlt(i)
        IF (ipm > nwts) GO TO 110
        IF (i == nt) GO TO 120
        l = i + 1
        DO  j = l, nt
          jp = ABS(ndx(j))
          IF (jp /= 0) THEN
            IF (jp <= ipm .AND. ip <= jp+mlt(j)) GO TO 110
          END IF
        END DO
      END IF
    END DO
    GO TO 120

END SELECT

70 ier = 15
RETURN

110 ier = 16
RETURN

!     GET MAXIMUM MULTIPLICITY TO SEE IF ENOUGH STORE IS AVAILABLE
120 mtm = 1
DO  i = 1, nt
  IF (ndx(i) > 0) mtm = MAX(mtm,mlt(i))
END DO

!     SET UP WORK FIELDS AND TEST SOME PARAMETERS
IF (nst < (n+1)/2) ier = 12
nmax = n + nst + MIN(2*mtm,n+1)
IF (zemu <= zero) ier = 13
IF (nmax > nwf) ier = -nmax
IF (ier /= 0) RETURN

!     TREAT A QF WITH 1 SIMPLE KNOT FIRST.
IF (n <= 1) THEN
  DO  i = 1, nt
    IF (ndx(i) > 0) GO TO 150
  END DO

!     WEIGHT NOT WANTED, DO NOTHING.
  RETURN
  150 wts(ABS(ndx(i))) = zemu
  RETURN
END IF

!     SKIP DIAGONALIZATION IF ALREADY DONE
IF (jdf == 0) THEN
  nw = 1
  
!     SET UNIT VECTOR IN WORK FIELD TO GET BACK 1ST ROW OF Q
  DO  i = 1, nst
    k = nw + i - 1
    wf(k) = zero
  END DO
  wf(nw) = one
  ierrx = 0
  
!     DIAGONALIZE JACOBI MATRIX
  CALL imtqlx(nst, aj, bj, wf(nw:), ierrx)
  
!     CHECK FOR ERROR
  IF (ierrx /= 0) THEN
    ier = 11
    RETURN
  END IF
  
!     SIGNAL JACOBI MATRIX NOW DIAGONALIZED SUCCESSFULLY AND SAVE
!     SQUARES OF 1ST ROW OF U IN SUBDIAGONAL ARRAY
  
  jdf = 1
  DO  i = 1, nst
    l = i + nw - 1
    bj(i) = wf(l) ** 2
  END DO
END IF

!     FIND ALL THE WEIGHTS FOR EACH KNOT FLAGGED
DO  i = 1, nt
  IF (ndx(i) > 0) THEN
    m = mlt(i)
    nm = n - m
    mnm = MAX(nm,1)
    l = MIN(m,nm+1)
    
!        SET UP INDICES FOR WORK FIELDS
    nk = nw + nst
    ny = nk + mnm
    nr = ny + m
    nz = nr + l
    
!     SET UP K-HAT MATRIX (FOR CWIQD) WITH KNOTS ACCORDING TO
!     THEIR MULTIPLICITIES
    k = nk
    DO  j = 1, nt
      IF (ndx(j) /= 0) THEN
        IF (j /= i) THEN
          mtj = mlt(j)
          DO  jj = 1, mtj
            wf(k) = t(j)
            k = k + 1
          END DO
        END IF
      END IF
    END DO
    
!        SET UP RIGHT PRINCIPAL VECTOR ARRAY FOR WEIGHTS ROUTINE
    wf(nr) = one / zemu
    DO  j = 2, l
      k = nr + j - 1
      wf(k) = zero
    END DO
    
!        PICK UP POINTER FOR THE LOCATION OF THE WEIGHTS TO BE OUTPUT
    k = ndx(i)
    
!        FIND ALL THE WEIGHTS FOR THIS KNOT
    
    CALL cwiqd(m, mnm, l, t(i), wf(nk:), nst, aj, bj, wf(nw:), wf(ny:),  &
               wf(nr:), wf(nz:), wts(k:))
    IF (key >= 0) THEN
      
!        DIVIDE BY FACTORIALS FOR WEIGHTS IN STANDARD FORM
      IF (m >= 2) THEN
        tmp = one
        ip = m - 1
        DO  j = 2, ip
          p = j
          tmp = tmp * p
          l = k + j
          wts(l) = wts(l) / tmp
        END DO
      END IF
    END IF
  END IF
END DO

RETURN
END SUBROUTINE cawiq



SUBROUTINE cwiqd(m, nm, l, v, xk, nstar, phi, a, wf, y, r, z, d)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-18  Time: 22:28:56
 
!  ROUTINE TO COMPUTE ALL THE WEIGHTS TO A GIVEN KNOT.
!  FOR DETAILS SEE:
!  KAUTSKY AND ELHAY "CALCULATION OF THE WEIGHTS OF INTERPOLATORY
!  QUADRATURES", NUMER MATH 40 (1982) 407-422.

!  VARIABLES NAMES USED CORRESPOND CLOSELY WITH THOSE IN THE ABOVE
!  MENTIONED PAPER
!   INPUT AND OUTPUT VARIABLES -
!                    I I  I I I  I     I   I O  O O O O
!   SUBROUTINE CWIQD(M,NM,L,V,XK,NSTAR,PHI,A,WF,Y,R,Z,D)

!         M      MULTIPLICITY OF THE KNOT IN QUESTION
!         NM     MAX(N-M,1). N=NUMBER OF KNOTS USED COUNTED
!                ACCORDING TO MULTIPLICITY
!         L      MIN(M,N-M+1)
!         V      THE KNOT IN QUESTION
!         XK     ALL BUT THE LAST M ENTRIES IN THE DIAGONAL OF K-HAT
!         NSTAR  DIMENSION OF THE JACOBI MATRIX
!         PHI    THE EIGENVALUES OF THE JACOBI MATRIX J
!         A      THE SQUARE OF THE 1ST ROW OF THE ORTHOGONAL MATRIX
!                DIAGONALIZING J
!         WF     WORK FIELD
!         Y      Y-HAT
!         R      VECTOR USED TO COMPUTE THE RIGHT PRINCIPAL VECTORS
!         Z      VECTOR USED TO COMPUTE THE LEFT PRINCIPAL VECTORS
!         D      OUTPUT ARRAY FOR THE WEIGHTS
!  OTHER VARIABLES ARE FOR TEMPORARY USE ONLY

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: nm
INTEGER, INTENT(IN)     :: l
REAL (dp), INTENT(IN)   :: v
REAL (dp), INTENT(IN)   :: xk(:)
INTEGER, INTENT(IN)     :: nstar
REAL (dp), INTENT(IN)   :: phi(:)
REAL (dp), INTENT(IN)   :: a(:)
REAL (dp), INTENT(OUT)  :: wf(:)
REAL (dp), INTENT(OUT)  :: y(:)
REAL (dp), INTENT(OUT)  :: r(:)
REAL (dp), INTENT(OUT)  :: z(:)
REAL (dp), INTENT(OUT)  :: d(:)

REAL (dp) :: sum, tmp
INTEGER   :: i, j, jr, k, last, minil

!     COMPUTE PRODUCTS REQUIRED FOR Y-HAT
DO  j = 1, nstar
  wf(j) = a(j)
  IF (nm >= 1) THEN
    DO  i = 1, nm
      wf(j) = wf(j) * (phi(j)-xk(i))
    END DO
  END IF
END DO

!     COMPUTE Y-HAT
DO  i = 1, m
  sum = zero
  DO  j = 1, nstar
    sum = sum + wf(j)
    wf(j) = wf(j) * (phi(j)-v)
  END DO
  y(i) = sum
END DO

!     IF N=1 THE RIGHT PRINCIPAL VECTOR IS ALREADY IN R.
IF (nm /= 0) THEN

!     OTHERWISE COMPUTE THE R-PRINCIPAL VECTOR OF GRADE M-1
  DO  i = 1, nm
    tmp = v - xk(i)
    IF (l /= 1) THEN
      last = MIN(l,i+1)
      DO  jr = 2, last
        j = last - jr + 2
        r(j) = tmp * r(j) + r(j-1)
      END DO
    END IF
    r(1) = tmp * r(1)
  END DO
END IF

!     COMPUTE LEFT PRINCIPAL VECTOR(S) AND WEIGHT FOR HIGHEST DERIV
!     THE FOLLOWING STATEMENT CONTAINS THE ONLY DIVISION IN THIS
!     ROUTINE. ANY TEST FOR OVERFLOW SHOULD BE MADE AFTER IT.
z(1) = one / r(1)
d(m) = y(m) * z(1)
IF (m == 1) RETURN

!     COMPUTE L-PRINCIPAL VECTOR
DO  i = 2, m
  sum = zero
  IF (l /= 1) THEN
    minil = MIN(i,l)
    DO  j = 2, minil
      k = i - j + 1
      sum = sum + r(j) * z(k)
    END DO
  END IF
  z(i) = -sum * z(1)
END DO

!     ACCUMULATE WEIGHTS
DO  i = 2, m
  sum = zero
  DO  j = 1, i
    k = m - i + j
    sum = sum + z(j) * y(k)
  END DO
  k = m - i + 1
  d(k) = sum
END DO
RETURN
END SUBROUTINE cwiqd



SUBROUTINE class(kind, m, alpha, beta, bj, aj, zemu, ier)

!  ROUTINE TO COMPUTE THE DIAGONAL (AJ) AND SUB-DIAGONAL (BJ) ELEMENTS OF THE
!  ORDER M (TRIDIAGONAL SYMMETRIC) JACOBI MATRIX ASSOCIATED WITH THE
!  POLYNOMIALS ORTHOGONAL WITH RESPECT TO THE WEIGHT FUNCTION SPECIFIED BY KIND.
!  FOR WEIGHT FUNCTIONS 1-7 M ELEMENTS ARE DEFINED IN BJ EVEN THOUGH ONLY M-1
!  ARE NEEDED.  FOR WEIGHT FUNCTION 8, BJ(M) IS SET TO ZERO.
!  THE ZERO-TH MOMENT OF THE WEIGHT FUNCTION IS RETURNED IN ZEMU.

!   INPUT AND OUTPUT VARIABLES -
!                    I    I I     I    O  O  O    O
!   SUBROUTINE CLASS(KIND,M,ALPHA,BETA,BJ,AJ,ZEMU,IER)

!  ERROR CODES ARE:
!  IER=1,2,3 PARAMETER RANGES ARE WRONG
!  IER=4 WEIGHT FUNCTION OF UNKNOWN TYPE. CANNOT GENERATE JACOBI MATRIX
!  IER=5 GAMMA FUNCTION DOES NOT MATCH MACHINE PARAMETERS IN ACCURACY

!  FUNCTIONS AND SUBROUTINES REFERENCED -  DGAMMA MACHEP SQRT PARCHK

INTEGER, INTENT(IN)     :: kind
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(OUT)  :: bj(:)
REAL (dp), INTENT(OUT)  :: aj(:)
REAL (dp), INTENT(OUT)  :: zemu
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: a2b2, ab, aba, abi, abj, abti, apone, temp

INTEGER :: i

REAL (dp), PARAMETER :: three = 3.0_dp, four = 4.0_dp
REAL (dp), PARAMETER :: pi = 3.14159265358979323846264338327950_dp

temp = EPSILON(one)
CALL parchk(kind, 2*m-1, alpha, beta, ier)
IF (ABS(dgamma(half)**2 - pi) > 5.0D2*temp) ier = 5
IF (ier /= 0) RETURN

!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (    1)
    ab = zero
    GO TO 30

  CASE (    2)
    zemu = pi
    DO  i = 1, m
      aj(i) = zero
      bj(i) = half
    END DO
    bj(1) = SQRT(half)
    RETURN

  CASE (    3)
    ab = alpha * two
    zemu = two ** (ab+one) * dgamma(alpha+one) ** 2 / dgamma(ab+two)
    aj(1) = zero
    bj(1) = one / (two*alpha+three)
    DO  i = 2, m
      aj(i) = zero
      bj(i) = i * (i+ab) / (four*(i+alpha)**2-one)
    END DO
    GO TO 180

  CASE (    4)
    ab = alpha + beta
    abi = two + ab
    zemu = two ** (ab+one) * dgamma(alpha+one) * dgamma(beta+one) / dgamma(abi)
    aj(1) = (beta-alpha) / abi
    bj(1) = four * (one+alpha) * (one+beta) / ((abi+one)*abi*abi)
    a2b2 = beta * beta - alpha * alpha
    DO  i = 2, m
      abi = two * i + ab
      aj(i) = a2b2 / ((abi-two)*abi)
      abi = abi ** 2
      bj(i) = four * i * (i+alpha) * (i+beta) * (i+ab) / ((abi-one)* abi)
    END DO
    GO TO 180

  CASE (    5)
    zemu = dgamma(alpha+one)
    DO  i = 1, m
      aj(i) = two * i - one + alpha
      bj(i) = i * (i+alpha)
    END DO
    GO TO 180

  CASE (    6)
    zemu = dgamma((alpha+one)/two)
    DO  i = 1, m
      aj(i) = zero
      bj(i) = (i+alpha*MOD(i,2)) / two
    END DO
    GO TO 180

  CASE (    7)
    ab = alpha
    GO TO 30

  CASE (    8)
    ab = alpha + beta
    zemu = dgamma(alpha+one) * dgamma(-(ab+one)) / dgamma(-beta)
    apone = alpha + one
    aba = ab * apone
    aj(1) = -apone / (ab+two)
    bj(1) = -aj(1) * (beta+one) / (ab+two) / (ab+three)
    DO  i = 2, m
      abti = ab + two * i
      aj(i) = aba + two * (ab+i) * (i-1)
      aj(i) = -aj(i) / abti / (abti-two)
    END DO
    DO  i = 2, m - 1
      abti = ab + two * i
      bj(i) = i * (alpha+i) / (abti-one) * (beta+i) /  &
              (abti**2) * (ab+i) / (abti+one)
    END DO
    bj(m) = zero
    GO TO 180

END SELECT

30 zemu = two / (ab+one)
DO  i = 1, m
  aj(i) = zero
  abi = i + ab * MOD(i,2)
  abj = 2 * i + ab
  bj(i) = abi * abi / (abj*abj-one)
END DO
GO TO 180

180 DO  i = 1, m
  bj(i) = SQRT(bj(i))
END DO

RETURN
END SUBROUTINE class



SUBROUTINE wm(w, m, kind, alpha, beta, ier)

!     ROUTINE TO EVALUATE THE FIRST M MOMENTS OF CLASSICAL WEIGHT FUNCTIONS

!      INPUT AND OUTPUT VARIABLES -
!                    O I I    I     I    O
!      SUBROUTINE WM(W,M,KIND,ALPHA,BETA,IER)

!     FUNCTIONS AND SUBROUTINES REFERENCED -  DGAMMA SQRT PARCHK

REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: als, sum, tmpa, tmpb, trm
INTEGER   :: i, ja, jb, k

REAL (dp), PARAMETER :: three = 3.0_dp
REAL (dp), PARAMETER :: pi = 3.14159265358979323846264338327950_dp

CALL parchk(kind, m, alpha, beta, ier)
IF (ier /= 0) RETURN
DO  k = 2, m, 2
  w(k) = zero
END DO

!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (    1)
    als = zero
    GO TO 40

  CASE (    2)
    w(1) = pi
    DO  k = 3, m, 2
      w(k) = w(k-2) * (k - two) / (k - one)
    END DO
    RETURN

  CASE (    3)
    w(1) = SQRT(pi) * dgamma(alpha + one) / dgamma(alpha + three/two)
    DO  k = 3, m, 2
      w(k) = w(k-2) * (k - two) / (two*alpha + k)
    END DO
    RETURN

  CASE (    4)
    als = alpha + beta + one
    w(1) = two ** als * dgamma(alpha+one) / dgamma(als+one) * dgamma(beta+one)
    DO  k = 2, m
      sum = zero
      trm = one
      DO  i = 0, (k-2) / 2
        tmpa = trm
        DO  ja = 1, 2 * i
          tmpa = tmpa * (alpha+ja) / (als+ja)
        END DO
        DO  jb = 1, k - 2 * i - 1
          tmpa = tmpa * (beta+jb) / (als+2*i+jb)
        END DO
        tmpa = tmpa / (2*i+one) * (2*i*(beta+alpha) + beta - (k-one)*alpha)  &
               / (beta+k-2*i-one)
        sum = sum + tmpa
        trm = trm * (k-2*i-one) / (2*i+one) * (k-2*i-two) / (2*i+two)
      END DO
      IF (MOD(k,2) /= 0) THEN
        tmpb = one
        DO  i = 1, k - 1
          tmpb = tmpb * (alpha+i) / (als+i)
        END DO
        sum = sum + tmpb
      END IF
      w(k) = sum * w(1)
    END DO
    RETURN

  CASE (    5)
    w(1) = dgamma(alpha+one)
    DO  k = 2, m
      w(k) = (alpha+k-one) * w(k-1)
    END DO
    RETURN

  CASE (    6)
    w(1) = dgamma((alpha+one)/two)
    DO  k = 3, m, 2
      w(k) = w(k-2) * (alpha+k-two) / two
    END DO
    RETURN

  CASE (    7)
    als = alpha
    GO TO 40

  CASE (    8)
    w(1) = dgamma(alpha+one) * dgamma(-alpha-beta-one) / dgamma(-beta)
    DO  k = 2, m
      w(k) = -w(k-1) * (alpha+k-one) / (alpha+beta+k)
    END DO
    RETURN

END SELECT

40 DO  k = 1, m, 2
  w(k) = two / (k + als)
END DO

RETURN
END SUBROUTINE wm



SUBROUTINE parchk(kind, m, alpha, beta, ier)
!  ROUTINE TO CHECK RANGES OF PARAMETERS ALPHA, BETA FOR CLASSICAL WEIGHT
!  FUNCTIONS.   M IS THE ORDER OF THE JACOBI MATRIX REQUIRED AND IS
!  CONSTRAINED BY ALPHA AND BETA FOR THE RATIONAL WEIGHT FUNCTION
!  (SEE BELOW).  M CAN BE REPLACED BY A DUMMY FOR OTHER WEIGHT FUNCTIONS.

!   INPUT AND OUTPUT VARIABLES -
!                     I    I I     I    O
!   SUBROUTINE PARCHK(KIND,M,ALPHA,BETA,IER)

!  FUNCTIONS AND SUBROUTINES REFERENCED - WM

INTEGER, INTENT(IN)    :: kind
INTEGER, INTENT(IN)    :: m
REAL (dp), INTENT(IN)  :: alpha
REAL (dp), INTENT(IN)  :: beta
INTEGER, INTENT(OUT)   :: ier

REAL (dp) :: tmp

!  CONSTRAINTS ON ALPHA,BETA:-
!         1      ALPHA>-1
!         2      FOR KIND<8 NEED BETA>-1
!         3      FOR KIND.EQ.8 NEED BETA < (ALPHA+BETA+2*M) < 0 TO
!                COMPUTE M ELEMENTS OF THE JACOBI MATRIX.
!   INPUT:
!   KIND...1-8 FOR CLASSICAL WEIGHT FUNCTION, 0 FOR UNKNOWN)
!   ALPHA,BETA...AS IN CLASS
!   M...ORDER OF HIGHEST MOMENT TO BE CALCULATED
!   OUTPUT:
!   IER...ERROR INDICATOR - CODED AS FOLLOWS
!         1...ALPHA <= -1
!         2...BETA <= -1
!         3...ALPHA,BETA COMBINATION WRONG FOR RATIONAL WEIGHT
!             FUNCTION
!         4...KIND = 0. PARAMETERS CANNOT BE CHECKED AND JACOBI MATRIX
!             IS NOT OF CLASSICAL TYPE

ier = 0
IF (kind <= 0) ier = 4

!     CHECK GEGENBAUER,JACOBI,LAGUERRE,HERMITE,EXPONENTIAL
IF (kind >= 3 .AND. (alpha <= -one)) ier = 1

!     CHECK BETA FOR JACOBI
IF (kind == 4 .AND. beta <= -one) ier = 2

!     CHECK RANGE FOR RATIONAL
IF (kind < 8) RETURN
tmp = alpha + beta + m + one
IF (tmp >= zero .OR. tmp <= beta) ier = 3

RETURN
END SUBROUTINE parchk



SUBROUTINE chkqfs(t, wts, mlt, nt, ndx, key, w, mop, mex, kind, alpha,  &
                  beta, lo, e, er, qm, ier)

! N.B. Argument NWTS has been removed.

!   ROUTINE TO CHECK THE POLYNOMIAL ACCURACY OF A QUADRATURE FORMULA.
!   IT WILL OPTIONALLY PRINT WEIGHTS, AND RESULTS OF A MOMENTS TEST.

!   INPUT AND OUTPUT VARIABLES -
!                     I I   I   I  I    I   I   * I   I
!   SUBROUTINE CHKQFS(T,WTS,MLT,NT,NWTS,NDX,KEY,W,MOP,MEX,
!  1                 KIND,ALPHA,BETA,LO,E,ER,QM,IER)
!                    I    I     I    I  O O  O  O

!   T...ARRAY OF DISTINCT KNOTS
!   W...MOMENTS ARRAY OF LENGTH MEX
!   MOP...THE EXPECTED ORDER OF PRECISION OF THE QF
!   MEX...THE TOTAL NUMBER (>1) OF MOMENTS REQUIRED TO BE TESTED
!         SET MEX=1 AND LO < 0 FOR NO MOMENTS CHECK
!   KIND...KIND OF CLASSICAL FORMULA.
!          KIND=0 MEANS UNKNOWN WEIGHT FUNCTION.
!          THE FIRST MEX MOMENTS MUST BE SET UP
!          IN ARRAY W BY THE USER FOR THIS CASE.
!   LO...PRINTING (IF ANY) IS DONE ON UNIT ABS(LO).  LO IS CODED
!       AS FOLLOWS:-
!       LO > 0 MEANS PRINT WEIGHTS AND MOMENT TESTS
!       LO .EQ. 0 MEANS PRINT NOTHING. COMPUTE MOMENT TEST
!       LO < 0 MEANS PRINT WEIGHTS ONLY.  DON'T COMPUTE MOMENT TESTS
!   E,ER...ABSOLUTE AND RELATIVE ERRORS OF THE QF APPLIED TO (X-DEL)**N
!   QM...VALUES OF THE QF APPLIED TO (X-DEL)**N
!   IER...ERROR INDICATOR. ANY ERROR RETURN COMES FROM WM.

!  FUNCTIONS AND SUBROUTINES REFERENCED - WM

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: wts(:)
INTEGER, INTENT(IN)        :: mlt(:)
INTEGER, INTENT(IN)        :: nt
INTEGER, INTENT(IN)        :: ndx(:)
INTEGER, INTENT(IN)        :: key
REAL (dp), INTENT(IN OUT)  :: w(:)
INTEGER, INTENT(IN)        :: mop
INTEGER, INTENT(IN)        :: mex
INTEGER, INTENT(IN)        :: kind
REAL (dp), INTENT(IN)      :: alpha
REAL (dp), INTENT(IN)      :: beta
INTEGER, INTENT(IN)        :: lo
REAL (dp), INTENT(OUT)     :: e(:)
REAL (dp), INTENT(OUT)     :: er(:)
REAL (dp), INTENT(OUT)     :: qm(:)
INTEGER, INTENT(OUT)       :: ier

REAL (dp) :: ek, emn, emx, erest, ern, erx, px, tmp, tmpx
INTEGER   :: i, j, jl, k, kindp, kjl, l, lu, m, mx

CHARACTER (LEN=53) :: txt1(10) =  &
 (/ '          INTERPOLATORY QUADRATURE FORMULA           ',  &
    ' TYPE  INTERVAL     WEIGHT FUNCTION        NAME      ',  &
    '   1    (-1,1)           ONE             LEGENDRE    ',  &
    '   2    (-1,1)     (1-X**2)**(-HALF)     CHEBYSHEV   ',  &
    '   3    (-1,1)      (1-X**2)**ALPHA      GEGENBAUER  ',  &
    '   4    (-1,1)  (1-X)**ALPHA*(1+X)**BETA  JACOBI     ',  &
    '   5   (0,INF)      X**ALPHA*EXP(-X)     GEN LAGUERRE',  &
    '   6  (-INF,INF) ABS(X)**ALFA*EXP(-X**2) GEN HERMITE ',  &
    '   7    (-1,1)       ABS(X)**ALPHA       EXPONENTIAL ',  &
    '   8   (0,INF)    X**ALPHA*(1+X)**BETA    RATIONAL   ' /)

lu = ABS(lo)

!     KIND MAY BE SET TO -1 TO ALLOW PRINTING OF MOMENTS ONLY
!     THIS FEATURE IS ONLY USED INTERNALLY (BY CHKQF)
kindp = MAX(0,kind)
IF (lo /= 0.AND.kind /= -1) THEN
  
  IF (kindp /= 0) THEN
    WRITE (lu,5000) txt1(1), txt1(2), txt1(kindp+2)
    IF (kindp >= 3) WRITE (lu,5100) alpha
    IF (kindp == 4 .OR. kindp == 8) WRITE (lu,5200) beta
  END IF
  IF (kind /= -1) WRITE (lu,5300) prec(1)
  
  WRITE (lu,5600)
  DO  i = 1, nt
    k = ABS(ndx(i))
    IF (k /= 0) THEN
      WRITE (lu,5400) i, t(i), mlt(i), wts(k)
      DO  j = k + 1, k + mlt(i) - 1
        WRITE (lu,5500) wts(j)
      END DO
    END IF
  END DO
END IF
IF (lo < 0) RETURN
ier = 0
IF (kindp /= 0) CALL wm(w, mex, kindp, alpha, beta, ier)
IF (ier /= 0) RETURN
qm(1:mex) = zero
erest = zero
DO  k = 1, nt
  tmp = one
  l = ABS(ndx(k))
  IF (l /= 0) THEN
    erest = erest + ABS(wts(l))
    DO  j = 1, mex
      qm(j) = qm(j) + tmp * wts(l)
      tmpx = tmp
      px = one
      DO  jl = 2, MIN(mlt(k),mex-j+1)
        kjl = j + jl - 1
        tmpx = tmpx * (kjl-1)
        qm(kjl) = qm(kjl) + tmpx * wts(l+jl-1) / px
        IF (key <= 0) THEN
          px = px * jl
        END IF
      END DO
      tmp = tmp * t(k)
    END DO
  END IF
END DO
DO  k = 1, mex
  e(k) = w(k) - qm(k)
  er(k) = e(k) / (ABS(w(k))+one)
END DO

!     FOR SOME STRANGE WEIGHT FUNCTIONS W(1) MAY VANISH
erest = erest / (ABS(w(1)) + one)

!     EXIT IF USER DOES NOT WANT PRINTED OUTPUT
IF (lo == 0) RETURN
emx = ABS(e(1))
emn = emx
erx = ABS(er(1))
ern = erx
m = mop + 1
mx = MIN(mop,mex)
DO  k = 2, mx
  emx = MAX(ABS(e(k)),emx)
  emn = MIN(ABS(e(k)),emn)
  erx = MAX(ABS(er(k)),erx)
  ern = MIN(ABS(er(k)),ern)
END DO
WRITE (lu,5700) mop, emn, ern, emx, erx, erest
IF (mex >= m) THEN
  ek = e(m)
  DO  j = 1, mop
    ek = ek / j
  END DO
  WRITE (lu,5800) mop, e(m), ek
END IF
WRITE (lu,5900)
WRITE (lu,6000) (j,w(j),qm(j),e(j),er(j),j = 1,mx)
WRITE (lu,6000)
IF (mex >= m) WRITE (lu,6000) (j,w(j),qm(j),e(j),er(j),j = m,mex)
RETURN

5000 FORMAT ('1', (t9, a72/))
5100 FORMAT (/'   PARAMETER(S) ALPHA   ', f12.5)
5200 FORMAT ('                BETA    ', f12.5)
5300 FORMAT (/'   MACHINE PRECISION   ', g13.1)
5400 FORMAT (2(i4, g26.17))
5500 FORMAT (t38, g26.17)
5600 FORMAT (/t12, 'KNOTS               MULT                WEIGHTS'/)
5700 FORMAT (//' COMPARISON OF MOMENTS'//  &
             ' ORDER OF PRECISION', i4//  &
             '  ERRORS :    ABSOLUTE    RELATIVE'/  &
             '---------+-------------------------'/ &
             ' MINIMUM :', 2g12.3/  &
             ' MAXIMUM :', 2g12.3// &
             ' WEIGHTS RATIO          ', g13.3)
5800 FORMAT (' ERROR FOR ', i3, '-TH POWER ', g13.3/,  &
             ' ERROR CONSTANT         ', g13.3,/)
5900 FORMAT (/'   MOMENTS: '/  &
             t13, 'TRUE            FROM Q.F.         ERROR     RELATIVE'/)
6000 FORMAT (i4, 2g19.10, 2g12.3)
END SUBROUTINE chkqfs



SUBROUTINE chkqf(t, wts, mlt, nt, ndx, key, wf, mop, mex, kind, alpha,  &
                 beta, lo, e, er, qm, nwf, a, b, ier)

! N.B. Argument NWTS has been removed.

!  ROUTINE TO COMPUTE AND PRINT THE MOMENTS OF A QF FOR
!  A CLASICAL WEIGHT FUNCTION WITH ANY VALID A,B
!  NO CHECK CAN BE MADE FOR NON-CLASSICAL WEIGHT FUNCTIONS

!   INPUT AND OUTPUT VARIABLES -
!                    I I   I   I  I    I   I   O  I   I   I
!   SUBROUTINE CHKQF(T,WTS,MLT,NT,NWTS,NDX,KEY,WF,MOP,MEX,KIND,
!  1                  ALPHA,BETA,LO,E,ER,QM,NWF,A,B,IER)
!                     I     I    I  O O  O  I   I I O

!   NWF...SIZE OF WORKFIELD ARRAY. MUST BE >= MEX+NT
!   MOP...THE EXPECTED ORDER OF PRECISION OF THE QF.
!   MEX...THE TOTAL NUMBER (>1) OF MOMENTS REQUIRED TO BE TESTED
!         SET MEX=1 AND LO<0 FOR NO MOMENTS CHECK
!   LO...PRINTING (IF ANY) IS DONE ON UNIT ABS(LO). LO IS CODED
!       AS FOLLOWS:-
!       LO>0 MEANS PRINT WEIGHTS AND MOMENT TESTS
!       LO.EQ.0 MEANS PRINT NOTHING. COMPUTE MOMENT TEST
!       LO<0 MEANS PRINT WEIGHTS ONLY. DON'T COMPUTE MOMENT TESTS
!   E,ER...ABSOLUTE AND RELATIVE ERRORS OF THE QF APPLIED
!          TO (X-DEL)**N
!   QM...VALUES OF THE QF APPLIED TO (X-DEL)**N
!   IER...ERROR INDICATOR. ANY ERROR RETURN COMES FROM WM.

!  IER CODES -
!  1-4       ERROR RETURN FROM PARCHK: ALPHA, BETA WRONG
!                  SEE ROUTINE PARCHECK
!    6       ZERO LENGTH INTERVAL (KIND=1,2,3,4,7)
!    7       B<=0 FOR (KIND=5,6)
!    8       A+B<=0 (KIND=8)

!  FUNCTIONS AND SUBROUTINES REFERENCED - CHKQFS PARCHK SCMM

REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(IN)   :: wts(:)
INTEGER, INTENT(IN)     :: mlt(:)
INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: ndx(:)
INTEGER, INTENT(IN)     :: key
REAL (dp), INTENT(OUT)  :: wf(:)
INTEGER, INTENT(IN)     :: mop
INTEGER, INTENT(IN)     :: mex
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(IN)     :: lo
REAL (dp), INTENT(OUT)  :: e(:)
REAL (dp), INTENT(OUT)  :: er(:)
REAL (dp), INTENT(OUT)  :: qm(:)
INTEGER, INTENT(IN)     :: nwf
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: tmp
INTEGER   :: i, izero, lex, lu, na, neg

CHARACTER (LEN=62) :: txt2(10) =  &
  (/ '          INTERPOLATORY QUADRATURE FORMULA                    ',  &
     '  TYPE  INTERVAL       WEIGHT FUNCTION               NAME     ',  &
     '    1    (A,B)              ONE                    LEGENDRE   ',  &
     '    2    (A,B)      ((B-X)*(X-A))**(-HALF)         CHEBYSHEV  ',  &
     '    3    (A,B)      ((B-X)*(X-A))**ALPHA           GEGENBAUER ',  &
     '    4    (A,B)    (B-X)**ALPHA*(X-A)**BETA          JACOBI    ',  &
     '    5   (A,INF)   (X-A)**ALPHA*EXP(-B*(X-A))      GEN LAGUERRE',  &
     '    6  (-INF,INF) ABS(X-A)**ALFA*EXP(-B*(X-A)**2) GEN HERMITE ',  &
     '    7    (A,B)      ABS(X-(A+B)/TWO)**ALPHA        EXPONENTIAL',  &
     '    8   (A,INF)    (X-A)**ALPHA*(B+X)**BETA         RATIONAL  ' /)

CALL parchk(kind, mex, alpha, beta, ier)
IF (ier /= 0) RETURN
IF (lo /= 0) THEN

!     PRINT WEIGHTS
  izero = 0
  lu = ABS(lo)
  WRITE (lu,5000) txt2(1), txt2(2), txt2(kind+2)
  WRITE (lu,5300) a
  WRITE (lu,5400) b
  IF (kind >= 3) WRITE (lu,5100) alpha
  IF (kind == 4 .OR. kind == 8) WRITE (lu,5200) beta
  CALL chkqfs(t, wts, mlt, nt, ndx, key, wf, mop, mex, izero, alpha,  &
              beta, -lu, e, er, qm, ier)
  IF (ier /= 0 .OR. lo < 0) RETURN
END IF
lex = mex + nt
IF (nwf < lex) THEN
  ier = -lex
  RETURN
END IF

CALL scmm(wf, mex, kind, alpha, beta, a, b, ier)
IF (ier /= 0) RETURN
na = mex + 1
tmp = (b+a) / two
IF (kind == 5 .OR. kind == 6 .OR. kind == 8) tmp = a
DO  i = 1, nt
  wf(na+i-1) = t(i) - tmp
END DO
neg = -1

!     CHECK MOMENTS
CALL chkqfs(wf(na:), wts, mlt, nt, ndx, key, wf, mop, mex, neg, alpha,  &
            beta, lo, e, er, qm, ier)
RETURN

5000 FORMAT ('1',(t9, a72/))
5100 FORMAT ('                  ALPHA      ', f12.5)
5200 FORMAT ('                  BETA       ', f12.5)
5300 FORMAT ('     PARAMETERS   A          ', f12.5)
5400 FORMAT ('                  B          ', f12.5)
END SUBROUTINE chkqf



SUBROUTINE sct(nt, t, st, kind, a, b, ier)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-18  Time: 22:31:11

! ROUTINE TO SCALE DISTINCT KNOTS FOR ANY VALID A,B TO THOSE FOR
! THE DEFAULT VALUES OF A,B. ARRAYS T AND ST MAY COINCIDE.
! ALL KNOTS IN THE T ARRAY ARE SCALED AND ARE OUTPUT IN ST.

!  INPUT AND OUTPUT VARIABLES -
!                 I  I O  I    I I O
!  SUBROUTINE SCT(NT,T,ST,KIND,A,B,IER)

! FUNCTIONS AND SUBROUTINES REFERENCED - MACHEP

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(OUT)  :: st(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: bma, shft, slp, tmp

ier = 0
IF (kind <= 0 .OR. kind > 8) THEN
  ier = 4
  RETURN
END IF

!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (1:4, 7)
    tmp = EPSILON(one)
    bma = b - a
    IF (bma > tmp) THEN
      slp = two / bma
      shft = -(a+b) / bma
    ELSE
      ier = 6
      RETURN
    END IF

  CASE (    5)
    IF (b < zero) THEN
      ier = 7
      RETURN
    END IF
    slp = b
    shft = -a * b

  CASE (    6)
    IF (b < zero) THEN
      ier = 7
      RETURN
    END IF
    slp = SQRT(b)
    shft = -a * slp

  CASE (    8)
    slp = one / (a+b)
    IF (slp <= zero) THEN
      ier = 8
      RETURN
    END IF
    shft = -a * slp

END SELECT

st(1:nt) = shft + slp * t(1:nt)

RETURN
END SUBROUTINE sct



SUBROUTINE scqf(nt, t, mlt, wts, ndx, swts, st, kind, alpha, beta, a, b, ier)

! N.B. Argument NWTS has been removed.

!  ROUTINE TO SCALE WEIGHTS AND KNOTS FOR CLASSICAL WEIGHT FUNCTION
!  WITH DEFAULT VALUES FOR A AND B TO THOSE FOR ANY VALID A,B

!   INPUT AND OUTPUT VARIABLES -
!                   I  I I   I   I    I   O    O
!   SUBROUTINE SCQF(NT,T,MLT,WTS,NWTS,NDX,SWTS,ST,
!  1                  KIND,ALPHA,BETA,A,B,IER)
!                     I    I     I    I I O

!  THE ARRAYS WTS AND SWTS MAY COINCIDE
!  THE ARRAYS T AND ST MAY COINCIDE
!  IER CODES -
!  1-4       ERROR RETURN FROM PARCHK: ALPHA, BETA WRONG
!                  SEE ROUTINE PARCHECK
!    6       ZERO LENGTH INTERVAL (KIND=1,2,3,4,7)
!    7       B<=0 FOR (KIND=5,6)
!    8       A+B<=0 (KIND=8)


!  FUNCTIONS AND SUBROUTINES REFERENCED - MACHEP SQRT PARCHK

INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(IN)   :: t(:)
INTEGER, INTENT(IN)     :: mlt(:)
REAL (dp), INTENT(IN)   :: wts(:)
INTEGER, INTENT(IN)     :: ndx(:)
REAL (dp), INTENT(OUT)  :: swts(:)
REAL (dp), INTENT(OUT)  :: st(:)
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: al, be, p, shft, slp, temp, tmp
INTEGER   :: i, k, l

temp = EPSILON(one)
CALL parchk(kind, 1, alpha, beta, ier)
IF (ier /= 0) RETURN

!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (    1)
    al = zero
    be = zero
    GO TO 60

  CASE (    2)
    al = -half
    be = -half
    GO TO 60

  CASE (    3)
    al = alpha
    be = alpha
    GO TO 60

  CASE (    4)
    al = alpha
    be = beta
    GO TO 60

  CASE (    5)
    IF (b <= zero) THEN
      ier = 7
      RETURN
    END IF
    shft = a
    slp = one / b
    al = alpha
    be = zero
    GO TO 100

  CASE (    6)
    IF (b <= zero) THEN
      ier = 7
      RETURN
    END IF
    shft = a
    slp = one / SQRT(b)
    al = alpha
    be = zero
    GO TO 100

  CASE (    7)
    al = alpha
    be = zero
    GO TO 60

  CASE (    8)
    IF (a+b <= zero) THEN
      ier = 8
      RETURN
    END IF
    shft = a
    slp = a + b
    al = alpha
    be = beta
    GO TO 100

END SELECT

60 IF ((b-a) <= temp) THEN
  ier = 6
  RETURN
END IF
shft = (a+b) / two
slp = (b-a) / two

100 p = slp ** (al+be+one)
DO  k = 1, nt
  st(k) = shft + slp * t(k)
  l = ABS(ndx(k))
  IF (l /= 0) THEN
    tmp = p
    DO  i = l, l + mlt(k) - 1
      swts(i) = wts(i) * tmp
      tmp = tmp * slp
    END DO
  END IF
END DO

RETURN
END SUBROUTINE scqf



SUBROUTINE scmm(w, m, kind, alpha, beta, a, b, ier)

!  ROUTINE TO COMPUTE MOMENTS OF CLASSICAL WEIGHT FUNCTION WITH
!  DEFAULT VALUES FOR A,B AND SCALE THEM TO THOSE FOR ANY VALID A,B

!   INPUT AND OUTPUT VARIABLES -
!                   O I I    I     I    I I O
!   SUBROUTINE SCMM(W,M,KIND,ALPHA,BETA,A,B,IER)

!  MOMENTS ARE OUTPUT TO W
!  FUNCTIONS AND SUBROUTINES REFERENCED - MACHEP SQRT WM

REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(OUT)    :: ier

REAL (dp) :: al, be, p, q, temp, tmp
INTEGER   :: i

temp = EPSILON(one)

!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (    1)
    al = zero
    be = zero
    GO TO 60

  CASE (    2)
    al = -half
    be = -half
    GO TO 60

  CASE (    3)
    al = alpha
    be = alpha
    GO TO 60

  CASE (    4)
    al = alpha
    be = beta
    GO TO 60

  CASE (    5)
    IF (b <= zero) THEN
      ier = 7
      RETURN
    END IF
    q = one / b
    p = q ** (alpha+one)
    GO TO 100

  CASE (    6)
    IF (b <= zero) THEN
      ier = 7
      RETURN
    END IF
    q = one / SQRT(b)
    p = q ** (alpha+one)
    GO TO 100

  CASE (    7)
    al = alpha
    be = zero
    GO TO 60

  CASE (    8)
    IF (a+b <= zero) THEN
      ier = 8
      RETURN
    END IF
    q = a + b
    p = q ** (alpha+beta+one)
    GO TO 100

END SELECT

60 IF ((b-a) <= temp) THEN
  ier = 6
  RETURN
END IF
q = (b-a) / two
p = q ** (al+be+one)

100 CALL wm(w, m, kind, alpha, beta, ier)
IF (ier /= 0) THEN
  RETURN
END IF
tmp = p
DO  i = 1, m
  w(i) = w(i) * tmp
  tmp = tmp * q
END DO

RETURN
END SUBROUTINE scmm



SUBROUTINE wtfn(t, w, nt, kind, alpha, beta, ier)

!  ROUTINE TO EVALUATE THE CLASSICAL WEIGHT FUNCTIONS AT THE POINTS
!  GIVEN IN ARRAY T. THE INPUT, T, AND OUTPUT, W, ARRAYS MAY BE THE SAME.

!   INPUT AND OUTPUT VARIABLES -
!                   I O I  I    I     I    O
!   SUBROUTINE WTFN(T,W,NT,KIND,ALPHA,BETA,IER)

!  *******WARNING*******
!  NO CHECK IS MADE
!     (1) THAT THE WEIGHT FUNCTION IS DEFINED FOR THE POINTS IN T.
!     (2) THAT THE POINTS ARE IN THE APPROPRIATE INTERVAL.
!  HOWEVER THE PARAMETERS ALPHA AND BETA ARE CHECKED FOR THE
!  CLASSICAL WEIGHT FUNCTIONS.

!  FUNCTIONS AND SUBROUTINES REFERENCED -  EXP SQRT PARCHK

REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: nt
INTEGER, INTENT(IN)     :: kind
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: beta
INTEGER, INTENT(OUT)    :: ier

INTEGER :: k

CALL parchk(kind, 1, alpha, beta, ier)
IF (ier /= 0) RETURN
!           LEG,CHEB,GEG,JAC,LAG,HERM,EXP,RAT
SELECT CASE ( kind )
  CASE (    1)
    GO TO 30

  CASE (    2)
    DO  k = 1, nt
      w(k) = one / SQRT((one-t(k))*(one+t(k)))
    END DO
    RETURN

  CASE (    3)
    IF (alpha == zero) GO TO 30
    DO  k = 1, nt
      w(k) = ((one-t(k))*(one+t(k))) ** alpha
    END DO
    RETURN

  CASE (    4)
    w(1:nt) = one
    IF (alpha /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * (one-t(k)) ** alpha
      END DO
    END IF
    IF (beta /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * (one+t(k)) ** beta
      END DO
    END IF
    RETURN

  CASE (    5)
    DO  k = 1, nt
      w(k) = EXP(-t(k))
    END DO
    IF (alpha /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * t(k) ** alpha
      END DO
    END IF
    RETURN

  CASE (    6)
    DO  k = 1, nt
      w(k) = EXP(-t(k)**2)
    END DO
    IF (alpha /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * ABS(t(k)) ** alpha
      END DO
    END IF
    RETURN

  CASE (    7)
    IF (alpha /= zero) THEN
      DO  k = 1, nt
        w(k) = ABS(t(k)) ** alpha
      END DO
      RETURN
    END IF

  CASE (    8)
    w(1:nt) = one
    IF (alpha /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * t(k) ** alpha
      END DO
    END IF
    IF (beta /= zero) THEN
      DO  k = 1, nt
        w(k) = w(k) * (one + t(k)) ** beta
      END DO
    END IF

END SELECT

30 w(1:nt) = one

RETURN
END SUBROUTINE wtfn



SUBROUTINE imtqlx(n, d, e, z, ier)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-18  Time: 22:40:19
 
!  THIS ROUTINE IS A SLIGHTLY MODIFIED VERSION OF THE EISPACK
!  ROUTINE TO PERFORM THE IMPLICIT QL ALGORITHM ON A SYMMETRIC
!  TRIDIAGONAL MATRIX. THE AUTHORS THANK THE AUTHORS OF EISPACK
!  FOR PERMISSION TO USE THIS ROUTINE. FOR DETAILS SEE
!  MARTIN AND WILKINSON: THE IMPLICIT QL ALGORITHM, NUMER MATH
!  12, 277-383 (1968). IT HAS BEEN MODIFIED TO PRODUCE THE
!           T
!  PRODUCT Q Z, WHERE Z IS AN INPUT VECTOR AND Q IS THE
!  ORTHOGONAL MATRIX DIAGONALIZING THE INPUT MATRIX.  THE CHANGES
!  CONSIST (ESSENTIALY) OF APPLYING THE ORTHOGONAL TRANSFORMATIONS
!  DIRECTLY TO Z AS THEY ARE GENERATED.  SEE REFERENCES TO Z NEAR
!  STATEMENT 60.

!  FUNCTIONS AND SUBROUTINES REFERENCED - SIGN SQRT

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: d(:)
REAL (dp), INTENT(OUT)     :: e(:)
REAL (dp), INTENT(IN OUT)  :: z(:)
INTEGER, INTENT(OUT)       :: ier

REAL (dp) :: b, c, f, g, p, r, s
INTEGER   :: i, ii, j, k, l, m, mml

INTEGER, PARAMETER  :: itn = 30

ier = 0
IF (n /= 1) THEN
  e(n) = zero
  DO  l = 1, n
    j = 0
    10 DO  m = l, n
      IF (m == n) EXIT
      IF (ABS(e(m)) <= prec(1)*(ABS(d(m)) + ABS(d(m+1)))) EXIT
    END DO

    p = d(l)
    IF (m /= l) THEN
      IF (j == itn) GO TO 80
      j = j + 1
      g = (d(l+1)-p) / (two*e(l))
      r = SQRT(g*g+one)
      g = d(m) - p + e(l) / (g + SIGN(r,g))
      s = one
      c = one
      p = zero
      mml = m - l
      DO  ii = 1, mml
        i = m - ii
        f = s * e(i)
        b = c * e(i)
        IF (ABS(f) >= ABS(g)) THEN
          c = g / f
          r = SQRT(c*c+one)
          e(i+1) = f * r
          s = one / r
          c = c * s
        ELSE
          s = f / g
          r = SQRT(s*s+one)
          e(i+1) = g * r
          c = one / r
          s = s * c
        END IF
        g = d(i+1) - p
        r = (d(i)-g) * s + two * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f
      END DO
      d(l) = d(l) - p
      e(l) = g
      e(m) = zero
      GO TO 10
    END IF
  END DO
  
  DO  ii = 2, n
    i = ii - 1
    k = i
    p = d(i)
    DO  j = ii, n
      IF (d(j) < p) THEN
        k = j
        p = d(j)
      END IF
    END DO
    IF (k /= i) THEN
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    END IF
  END DO
  GO TO 90
  
  80 ier = l
END IF

90 RETURN
END SUBROUTINE imtqlx



FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
             -1259.139216722289_dp, 771.3234287757674_dp, &
             -176.6150291498386_dp, 12.50734324009056_dp, &
             -0.1385710331296526_dp, 0.9934937113930748D-05, &
              0.1659470187408462D-06 /), zero = 0._dp, one = 1._dp,  &
              lnsqrt2pi = 0.9189385332046727_dp, half = 0.5_dp,  &
              sixpt5 = 6.5_dp, seven = 7._dp, tmp
INTEGER   :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma



FUNCTION dgamma(x) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = EXP( lngamma(x) )

RETURN
END FUNCTION dgamma

END MODULE iqpack

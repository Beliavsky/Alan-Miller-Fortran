PROGRAM wrtest
USE AiryFunction
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! Code converted using TO_F90 by Alan Miller
! Date: 2002-11-04  Time: 21:28:55
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC TEST PROGRAM OF AIZ, BIZ: TEST OF WRONSKIAN RELATIONS.
!CC CHECK THE FOLLOWING RELATIONS FROM
!CC   ABRAMOWITZ & STEGUN, "HANDBOOK OF MATHEMATICAL FUNCTIONS":
!CC  1.  EQ. 10.4.10
!CC  2.  EQ. 10.4.11
!CC  3.  EQ. 10.4.12
!CC
!CC  NOTE THAT CHECKS 1,2 MAKE SENSE FOR |ARG(Z)| <= PI/3
!CC  WHILE CHECK 3 IS VALID FOR -PI/3 <= ARG(Z) <= PI
!CC
!CC  4.  ALSO WE TEST ANOTHER WRONSKIAN RELATION FOR BIZ
!CC      WHICH IS VALID IN THE SECTOR PI/3 <= ARG(Z) <= PI:
!CC
!CC      W[BI,AI(ZEXP(-I2PI/3))] = EXP(I2PI/3))/2/PI
!CC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  WE TEST THE RELATIONS FOR Z=R*EXP(I*PHI) WHERE
!CC               0 <= R <= 100,  0 <= PHI <= PI
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
REAL (dp) :: pi, sqrt3, x, y, c, s, air, aii, pair, paii, bir, bii,  &
             pbir, pbii, u, v, wraibi, ar, ai, par, pai, ar1, ai1,  &
             apr, api, par1, pai1, wrm1, wr01, wraib2, r, ph, w1max,  &
             w2max, w3max, w4max
INTEGER   :: i, j, ifac, ierr1, ierr2, ierrb

pi = 3.1415926535897932385D0
sqrt3 = 1.7320508075688772935D0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C THIS IS A WRONSKIAN TEST FOR THE SCALED FUNCTIONS       C
ifac = 2
!C ONE CAN CHOOSE A WRONSKIAN TEST FOR THE UNSCALED        C
!C FUNCTIONS BY UNCOMMENTING THE FOLLOWING LINE:           C
!C        IFAC=1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

OPEN (UNIT=10,FILE='res1',STATUS='UNKNOWN')
WRITE (10,*) ' '
WRITE (10,*) 'TEST OF WRONSKIAN RELATIONS FOR AIZ, BIZ'
WRITE (10,*) 'CHECK THE FOLLOWING RELATION:'
WRITE (10,*) ' '
WRITE (10,*) 'CHECK1. ABRAMOWITZ & STEGUN, 10.4.10'
WRITE (10,*) ' '
WRITE (10,5000) 'X', 'Y', 'CHECK1'
OPEN (UNIT=11,FILE='res2',STATUS='UNKNOWN')
WRITE (11,*) ' '
WRITE (11,*) 'TESTS OF WRONSKIAN RELATIONS FOR AIZ'
WRITE (11,*) 'CHECK THE FOLLOWING RELATION:'
WRITE (11,*) ' '
WRITE (11,*) 'CHECK2. ABRAMOWITZ & STEGUN, 10.4.11'
WRITE (11,*) ' '
WRITE (11,5000) 'X', 'Y', 'CHECK2'
OPEN (UNIT=12,FILE='res3',STATUS='UNKNOWN')
WRITE (12,*) ' '
WRITE (12,*) 'FURTHER TESTS OF WRONSKIAN RELATIONS FOR AIZ'
WRITE (12,*) 'CHECK THE FOLLOWING RELATION:'
WRITE (12,*) 'CHECK3. ABRAMOWITZ & STEGUN, 10.4.12'
WRITE (12,*) ' '
WRITE (12,5100) 'X', 'Y', 'CHECK3'
OPEN (UNIT=14,FILE='res4',STATUS='UNKNOWN')
WRITE (14,*) ' '
WRITE (14,*) 'ANOTHER TEST OF WRONSKIAN RELATION FOR BIZ:'
WRITE (14,*) ' W[BI,AI(ZEXP(-I2PI/3))]=EXP(I2PI/3))/2/PI '
WRITE (14,*) 'WHICH CAN BE USED IN THE SECTOR PI > |ARG(Z)| > PI/3'
WRITE (14,*) ' '
WRITE (14,5100) 'X', 'Y', 'CHECK4'
w1max = 0.d0
w2max = 0.d0
w3max = 0.d0
w4max = 0.d0
DO  i = 1, 100
  DO  j = 1, 180, 2
    ierr1 = 0
    ierr2 = 0
    ierrb = 0
    r = i
    ph = j * pi / 180.d0
    x = r * COS(ph)
    y = r * SIN(ph)
    IF (ph <= pi/3.d0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC    CHECK A & S 10.4.10
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL aiz(1, ifac, x, y, air, aii, ierr1)
      CALL aiz(2, ifac, x, y, pair, paii, ierr1)
      CALL biz(1, ifac, x, y, bir, bii, ierrb)
      CALL biz(2, ifac, x, y, pbir, pbii, ierrb)
      IF (ierr1 == 0 .AND. ierrb == 0) THEN
        u = pi * ((air*pbir-aii*pbii) - (pair*bir-paii*bii))
        v = pi * ((aii*pbir+air*pbii) - (paii*bir+pair*bii))
        wraibi = ABS(SQRT(u*u + v*v) - 1.d0)
        IF (wraibi > w1max) THEN
          w1max = wraibi
        END IF
        WRITE (10,5200) x, y, wraibi
      END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC    CHECK A & S 10.4.11
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL aiz(1, ifac, x, y, ar, ai, ierr1)
      CALL aiz(2, ifac, x, y, par, pai, ierr1)
      c = -0.5D0
      s = 0.5D0 * sqrt3
      u = x * c - y * s
      v = x * s + y * c
      CALL aiz(1, ifac, u, v, ar1, ai1, ierr2)
      CALL aiz(2, ifac, u, v, apr, api, ierr2)
      IF ((ierr1 == 0) .AND. (ierr2 == 0)) THEN
        u = apr * c - api * s
        v = apr * s + api * c
        par1 = u
        pai1 = v
        u = pi * ((ar*par1-ai*pai1)-(par*ar1-pai*ai1)) - sqrt3 / 4.d0
        v = pi * ((ar*pai1+ai*par1)-(par*ai1+pai*ar1)) + 1.d0 / 4.d0
        wrm1 = SQRT(u*u + v*v)
        IF (wrm1 > w2max) THEN
          w2max = wrm1
        END IF
        WRITE (11,5200) x, y, wrm1
      END IF
    ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   ANOTHER CHECK FOR THE BS
!CC   IN THE SECTOR ARG(Z) > PI/3
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      c = -0.5D0
      s = -0.5D0 * sqrt3
      u = x * c - y * s
      v = x * s + y * c
      CALL aiz(1, ifac, u, v, air, aii, ierr1)
      CALL aiz(2, ifac, u, v, pair, paii, ierr1)
      CALL biz(1, ifac, x, y, bir, bii, ierrb)
      CALL biz(2, ifac, x, y, pbir, pbii, ierrb)
      IF ((ierr1 == 0) .AND. (ierrb == 0)) THEN
        u = pi * ((bir*pair-bii*paii) - sqrt3*(bii*pair+bir*paii) +  &
            2.d0*(pbir*air-pbii*aii)) - 0.5D0
        v = pi * ((bii*pair+bir*paii) + sqrt3*(bir*pair-bii*paii) +  &
            2.d0*(pbii*air+pbir*aii)) + sqrt3 * 0.5D0
        wraib2 = SQRT(u*u + v*v)
        IF (wraib2 > w3max) THEN
          w3max = wraib2
        END IF
        WRITE (14,5200) x, y, wraib2
      END IF
    END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC    CHECK A & S 10.4.12
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    CALL aiz(1, ifac, x, y, ar, ai, ierr1)
    CALL aiz(2, ifac, x, y, par, pai, ierr1)
    c = -0.5D0
    s = -0.5D0 * sqrt3
    u = x * c - y * s
    v = x * s + y * c
    CALL aiz(1, ifac, u, v, ar1, ai1, ierr2)
    CALL aiz(2, ifac, u, v, apr, api, ierr2)
    IF ((ierr1 == 0) .AND. (ierr2 == 0)) THEN
      u = apr * c - api * s
      v = apr * s + api * c
      par1 = u
      pai1 = v
      u = pi * ((ar*par1-ai*pai1) - (par*ar1-pai*ai1)) - sqrt3 / 4.d0
      v = pi * ((ar*pai1+ai*par1) - (par*ai1+pai*ar1)) - 1.d0 / 4.d0
      wr01 = SQRT(u*u + v*v)
      IF (wr01 > w4max) THEN
        w4max = wr01
      END IF
      WRITE (12,5200) x, y, wr01
    END IF
  END DO
END DO
WRITE (10,*) ' '
WRITE (10,*) 'MAXIMUM VALUE OF THE WRONSKIAN CHECK =', w1max
WRITE (11,*) ' '
WRITE (11,*) 'MAXIMUM VALUE OF THE WRONSKIAN CHECK =', w2max
WRITE (14,*) ' '
WRITE (14,*) 'MAXIMUM VALUE OF THE WRONSKIAN CHECK =', w3max
WRITE (12,*) ' '
WRITE (12,*) 'MAXIMUM VALUE OF THE WRONSKIAN CHECK =', w4max
STOP

5000 FORMAT (t7, a1, t26, a1, t39, a6, t56, a6)
5100 FORMAT (t7, a1, t26, a1, t39, a6)
5200 FORMAT (g16.9, ' ', g16.9, ' ', g16.9)
END PROGRAM wrtest

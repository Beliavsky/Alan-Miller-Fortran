MODULE Adventure_Subs

! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-18  Time: 13:11:36

USE Adventure_Common
IMPLICIT NONE


CONTAINS


!***   A5TOA1

SUBROUTINE a5toa1 (a, b, c, qq, leng)

!  A AND B CONTAIN A 1- TO 12-CHARACTER WORD IN A6 FORMAT, C CONTAINS AN
!  WORD AND/OR PUNCTUATION.  THEY ARE UNPACKED TO ONE CHARACTER PER WORD
!  ARRAY "CHARS".
!  THE INDEX OF THE LAST NON-BLANK CHAR IN CHARS IS RETURNED IN LENG.

CHARACTER (LEN=6), INTENT(IN)   :: a
CHARACTER (LEN=6), INTENT(IN)   :: b
CHARACTER (LEN=6), INTENT(IN)   :: c
CHARACTER (LEN=1), INTENT(OUT)  :: qq(20)
INTEGER, INTENT(OUT)            :: leng

! IMPLICIT INTEGER (a-z)
INTEGER  :: i, jj

DO  jj=1,6
  qq(jj)=a(jj:jj)
  qq(jj+6)=b(jj:jj)
END DO
DO  i=1,12
  IF (qq(i) == ' ') GO TO 30
END DO
leng=12
GO TO 40

30 leng=i-1
40 DO  i=1,6
  IF (c(i:i) /= ' ') THEN
    leng=leng+1
    qq(leng)=c(i:i)
  END IF
END DO
DO  i=leng+1,20
  qq(i)=' '
END DO
DO  i=1,20
  IF (qq(i) == '_') qq(i)=' '
END DO

RETURN
END SUBROUTINE a5toa1



!***  AJAR   .TRUE. IF OBJ IS CONTAINER AND IS OPEN
!  THE NEXT LOGICAL FUNCTIONS DESCRIBE ATTRIBUTES OF OBJECTS.
!  (AJAR, HINGED, OPAQUE, PRINTD, TREASR, VESSEL, WEARNG)

FUNCTION ajar(obj) RESULT(lval)

!  AJAR(OBJ)    = TRUE IF OBJECT IS AN OPEN OR UNHINGED CONTAINER.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: bitset,hinged,vessel
! COMMON /bitcom/ openbt,unlkbt,burnbt,wearbt
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),openbt) .OR. (vessel(obj) .AND..NOT.hinged(obj))
RETURN
END FUNCTION ajar



!***  AT     .TRUE. IF AT OBJ

FUNCTION at(obj) RESULT(lval)

!  AT(OBJ)      = TRUE IF ON EITHER SIDE OF TWO-PLACED OBJECT

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

lval=.false.
IF (obj < 1 .OR. obj > maxobj) RETURN
lval=place(obj) == loc .OR. fixed(obj) == loc
RETURN
END FUNCTION at



!***  ATHAND .TRUE. IF OBJ READILY AVAILABLE

FUNCTION athand(obj) RESULT(lval)

!  ATHAND(OBJ)  = TRUE IF OBJ IS READILY REACHABLE.
!                 IT CAN BE LYING HERE, IN HAND OR IN OPEN CONTAINER.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! LOGICAL  :: toting,ajar,enclsd,holdng,aaa

INTEGER  :: contnr

lval=.false.
IF (place(obj) == loc .OR. holdng(obj)) THEN
  lval=.true.
  RETURN
END IF
IF (.NOT.enclsd(obj)) RETURN
contnr=-place(obj)
lval=(ajar(contnr) .AND. (place(contnr) == loc .OR. (toting(obj)  &
     .AND. holdng(contnr))))

RETURN
END FUNCTION athand




!***   BITS

FUNCTION bits(shift) RESULT(ival)

INTEGER, INTENT(IN)  :: shift
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
! ival=(2**shift)
ival = 0
ival = IBSET(ival,shift)
RETURN
END FUNCTION bits




!***   BITOFF

SUBROUTINE bitoff (obj, bit)

!  TURNS OFF (SETS=0) A BIT IN OBJCON.

INTEGER, INTENT(IN)  :: obj
INTEGER, INTENT(IN)  :: bit

! IMPLICIT INTEGER (a-z)
! COMMON /concom/ loccon(250),objcon(150)

!      OBJCON(OBJ)=IAND(OBJCON(OBJ),INOT(BITS(BIT)))
! THE FOLLOWING SHOULD BE EQUIVALENT TO THE ABOVE
objcon(obj)=IOR(objcon(obj),(bits(bit))) - bits(bit)
RETURN
END SUBROUTINE bitoff




!***   BITON

SUBROUTINE biton (obj, bit)

!  TURNS ON (SETS=1) A BIT IN OBJCON.

INTEGER, INTENT(IN)  :: obj
INTEGER, INTENT(IN)  :: bit

! IMPLICIT INTEGER (a-z)
! COMMON /concom/ loccon(250),objcon(150)

objcon(obj)=IOR(objcon(obj),bits(bit))
RETURN
END SUBROUTINE biton





!***   BITSET
!  MISCELLANEOUS LOGICAL FUNCTIONS (BITSET, PCT)
!  ALSO, SUBROUTINES FOR TURNING BITS ON AND OFF (BITON, BITOFF).

FUNCTION bitset(word,n) RESULT(lval)

!  BITSET(COND,L,N) = TRUE IF COND(L) HAS BIT N SET

INTEGER, INTENT(IN)  :: word, n
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! bitset=IAND(word,2**n) /= 0
lval = BTEST(word,n)
RETURN
END FUNCTION bitset




!***  BLIND  .TRUE. IF YOU CAN'T SEE AT THIS LOC
!  LOCATION ATTRIBUTES.  (BLIND, DARK, FORCED, INSIDE, OUTSID, PORTAL)

FUNCTION blind() RESULT(lval)

!  TRUE IF ADVENTURER IS "BLIND" AT THIS LOC, (DARKNESS OR GLARE)

LOGICAL  :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /concom/ loccon(250),objcon(150)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150), points(150)
! LOGICAL  :: dark,athand
INTEGER  :: lamp = 2

lval=dark() .OR. (loc == 200 .AND. athand(lamp) .AND. prop(lamp) == 1)

RETURN
END FUNCTION blind




!***   BUG

SUBROUTINE bug (num)

INTEGER, INTENT(IN)  :: num

! IMPLICIT INTEGER (a-z)

!  THE FOLLOWING CONDITIONS ARE CURRENTLY CONSIDERED FATAL BUGS.  NUMBER
!  ARE DETECTED WHILE READING THE DATABASE; THE OTHERS OCCUR AT "RUN TIM
!       0       MESSAGE LINE > 70 CHARACTERS
!       1       NULL LINE IN MESSAGE
!       2       TOO MANY WORDS OF MESSAGES
!       3       TOO MANY TRAVEL OPTIONS
!       4       TOO MANY VOCABULARY WORDS
!       5       REQUIRED VOCABULARY WORD NOT FOUND
!       6       TOO MANY RTEXT OR MTEXT MESSAGES
!       7       TOO MANY HINTS
!       8       LOCATION HAS COND BIT BEING SET TWICE
!       9       INVALID SECTION NUMBER IN DATABASE
!       10      OUT OF ORDER LOCS OR RSPEAK ENTRIES.
!       11      ILLEGAL MOTION WORD IN TRAVEL TABLE
!       12      ** UNUSED **.
!       13      UNKNOWN OR ILLEGAL WORD IN ADJECTIVE TABLE.
!       14      ILLEGAL WORD IN PREP/OBJ TABLE
!       15      TOO MANY ENTRIES IN PREP/OBJ TABLE
!       16      OBJECT HAS CONDITION BIT SET TWICE
!       17      OBJECT NUMBER TOO LARGE
!       18      TOO MANY ENTRIES IN ADJECTIVE/NOUN TABLE.
!       20      SPECIAL TRAVEL (500>L>300) EXCEEDS GOTO LIST
!       21      RAN OFF END OF VOCABULARY TABLE
!       22      VERB CLASS (N/1000) NOT BETWEEN 1 AND 3
!       23      INTRANSITIVE ACTION VERB EXCEEDS GOTO LIST
!       24      TRANSITIVE ACTION VERB EXCEEDS GOTO LIST
!       25      CONDITIONAL TRAVEL ENTRY WITH NO ALTERNATIVE
!       26      LOCATION HAS NO TRAVEL ENTRIES
!       27      HINT NUMBER EXCEEDS GOTO LIST
!       28      INVALID MONTH RETURNED BY DATE FUNCTION
!       29      ACTION VERB 'LEAVE' HAS NO OBJECT.
!       30      PREPOSITION FOUND IN UNEXPECTED TABLE
!       31      RECEIVED AN UNEXPECTED WORD TERMINATOR FROM A1TOA5
!       32      TRYING TO PUT A CONTAINER INTO ITSELF (TRICKY!)
!       33      UNKNOWN WORD CLASS IN GETWDS
!       35      TRYING TO CARRY A NON-EXISTENT OBJECT

WRITE (*,10) num
10 FORMAT (' FATAL ERROR, SEE SOURCE CODE FOR INTERPRETATION.'/  &
           ' PROBABLE CAUSE: ERRONEOUS INFO IN DATABASE OR BAD ASAVE.DAT'/  &
           ' ERROR CODE =', i2/)
STOP
END SUBROUTINE bug



!***  BURDEN .. RETURNS WEIGHT OF ITEMS BEING CARRIED

FUNCTION burden(obj) RESULT(ival)

!  IF OBJ=0, BURDEN CALCULATES THE TOTAL WEIGHT OF THE ADVENTURER'S BURD
!       INCLUDING EVERYTHING IN ALL CONTAINERS (EXCEPT THE BOAT) THAT HE
!       CARRYING.
!  IF OBJ#0 AND OBJ IS A CONTAINER, CALCULATE THE WEIGHT OF EVERYTHING I
!       THE CONTAINER (INCLUDING THE CONTAINER ITSELF).  SINCE DONKEY FO
!       ISN'T RECURSIVE, WE WILL ONLY CALCULATE WEIGHTS OF CONTAINED CON
!       ONE LEVEL DOWN.  THE ONLY SERIOUS CONTAINED CONTAINER WOULD BE T
!       THE ONLY THINGS WE'LL MISS WILL BE FILLED VS EMPTY BOTTLE OR CAG
!  IF OBJ#0 AND ISN'T A CONTAINER, RETURN ITS WEIGHT.

INTEGER, INTENT(IN)  :: obj
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150), points(150)
! COMMON /hldcom/ holder(150),hlink(150)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! LOGICAL :: toting,wearng
INTEGER  :: boat = 48, temp

ival=0
IF (obj /= 0) GO TO 20
DO  i=1,maxobj
  IF (.NOT.toting(i) .OR. place(i) == -boat) CYCLE
  ival=ival + weight(i)
END DO
RETURN

20 ival=weight(obj)
IF (obj == boat) RETURN
temp=holder(obj)
30 IF (temp == 0) RETURN
ival=ival + weight(temp)
temp=hlink(temp)
GO TO 30

END FUNCTION burden



!***   CARRY

SUBROUTINE carry (object, where)

!  START TOTING AN OBJECT, REMOVING IT FROM THE LIST OF THINGS AT ITS FO
!  LOCATION.  IF OBJECT>MAXOBJ (MOVING "FIXED" SECOND LOC),
!  DON'T CHANGE PLACE.

INTEGER, INTENT(IN)  :: object
INTEGER, INTENT(IN)  :: where

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

INTEGER  :: temp

IF (object > maxobj) GO TO 10
IF (place(object) == -1) RETURN
place(object)=-1
10 IF (atloc(where) /= object) GO TO 20
atloc(where)=link(object)
RETURN

20 temp=atloc(where)
30 IF (link(temp) == object) GO TO 40
temp=link(temp)
IF (temp /= 0) GO TO 30
CALL bug (35)

40 link(temp)=link(object)
RETURN
END SUBROUTINE carry



!***   CLASS

FUNCTION class(word) RESULT(ival)

!  RETURNS WORD CLASS NUMBER (1=MOTION VERB; 2=NOUN; 3=ACTION VERB;
!  4=MISCELLANEOUS WORD; 5=PREPOSITION; 6=ADJECTIVE; 7=CONJUNCTION).

INTEGER, INTENT(IN)  :: word
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)

ival=word/1000 + 1
IF (word < 0) ival=-1
RETURN
END FUNCTION class



!***   CLRLIN

SUBROUTINE clrlin()

!  CLEARS OUT ALL CURRENT SYNTAX ARGS IN PREPARATION FOR A NEW INPUT LIN

! IMPLICIT INTEGER (a-z)
CHARACTER (LEN=6)  :: allzero
! COMMON /wrdcom/ verbs(45),vrbx,objs(45),objx,iobjs(15),iobx,prep,words(45)
! COMMON /sv3com/ dtk(9),atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)

DO  i=1,6
  allzero(i:i)=CHAR(0)
END DO
DO  i=1,45
  objs(i)=0
  verbs(i)=0
  DO  j=1,2
    vtxt(i,j)=allzero
  END DO
END DO

DO  i=1,15
  iobjs(i)=0
  DO  j=1,2
    iotxt(i,j)=allzero
    otxt(i,j)=allzero
  END DO
END DO

vrbx=0
objx=0
iobx=0
prep=0
RETURN

END SUBROUTINE clrlin



!***   CONFUZ

FUNCTION confuz() RESULT(ival)

!  GENERATES SOME VARIANT OF "DON'T UNDERSTAND THAT" MESSAGE.

INTEGER  :: ival

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: pct

ival=60
IF (pct(50)) ival=61
IF (pct(33)) ival=13
IF (pct(25)) ival=347
IF (pct(20)) ival=195
RETURN
END FUNCTION confuz




!***  DARK   .TRUE. IF THERE IS NO LIGHT HERE

FUNCTION dark() RESULT(lval)

!  TRUE IF LOCATION "LOC" IS DARK

LOGICAL  :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /concom/ loccon(250),objcon(150)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150),points(150)
! LOGICAL  :: athand
INTEGER  :: lamp = 2

lval=MOD(loccon(loc),2) == 0 .AND. (prop(lamp)  == 0 .OR. .NOT.athand(lamp))
RETURN
END FUNCTION dark



!***   DEAD  .TRUE. IF OBJ IS NOW DEAD

FUNCTION dead(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),10)
RETURN
END FUNCTION dead



!***   DROP

SUBROUTINE drop (object, where)

!  PLACE AN OBJECT AT A GIVEN LOC, PREFIXING IT ONTO THE ATLOC LIST.

INTEGER, INTENT(IN)  :: object, where

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

IF (object > maxobj) GO TO 10
place(object)=where
GO TO 20

10 fixed(object-maxobj)=where
20 IF (where <= 0) RETURN
link(object)=atloc(where)
atloc(where)=object
RETURN
END SUBROUTINE drop




!***   DSTROY

SUBROUTINE dstroy (object)

!  PERMANENTLY ELIMINATE "OBJECT" BY MOVING TO A NON-EXISTENT LOCATION.

INTEGER, INTENT(IN)  :: object

! IMPLICIT INTEGER (a-z)

CALL move (object, 0)
RETURN
END SUBROUTINE dstroy




!***   EDIBLE  .TRUE. IF OBJ CAN BE EATEN

FUNCTION edible(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),7)
RETURN
END FUNCTION edible



!***  ENCLSD .TURE. IF OBJ INSIDE SOMETHING

FUNCTION enclsd(object) RESULT(lval)

!  ENCLSD(OBJ) = TRUE IF THE OBJ IS IN A CONTAINER

INTEGER, INTENT(IN)  :: object
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

lval=.false.
IF (object < 1 .OR. object > maxobj) RETURN
lval=place(object) < -1
RETURN
END FUNCTION enclsd




!***   FORCED

FUNCTION forced(loc) RESULT(lval)

!  A FORCED LOCATION IS ONE FROM WHICH HE IS IMMEDIATELY BOUNCED TO ANOT
!  NORMAL USE IS FOR DEATH (FORCE TO LOC ZERO) AND FOR DESCRIPTIONS OF
!  JOURNEY FROM ONE PLACE TO ANOTHER.

INTEGER, INTENT(IN)  :: loc
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /concom/ loccon(250),objcon(150)

lval=loccon(loc) == 2

RETURN
END FUNCTION forced



!***   GETLIN

SUBROUTINE getlin()

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: blklin
! COMMON /blkcom/ blklin
CHARACTER (LEN=6) :: txt2(35,2)
! COMMON /utxcom/ wdx
CHARACTER (LEN=1) :: chrs(150),chr2(150),chrx(70)
! COMMON /sv3com/ dtk(9),atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)

INTEGER  :: i, indx, indx2, j

DO  i=1,35
  DO  j=1,2
    txt2(i,j)='      '
    txt(i,j)='      '
  END DO
END DO

20 IF (blklin) WRITE (*,*) '  '
WRITE (*, 30, ADVANCE='NO')
!1266 FORMAT(' >')
!??????????????????????????????? THE FOLLOWING WORKS ON MANY COMPUTERS:
30 FORMAT (' > ')
READ (*,40) chrx
40 FORMAT (70A1)
DO  i=1,70
  chrx(i)=CHAR(IAND(ICHAR(chrx(i)),127))
  IF (chrx(i) < ' ') chrx(i)=' '
  IF (chrx(i) >= '`') chrx(i)=CHAR(ICHAR(chrx(i))-32)
END DO

DO  i=1,70
  IF (chrx(i) /= ' ') GO TO 70
END DO
GO TO 20

70 indx=1
DO  i=1,70
  IF (chrx(i) == '.'.OR.chrx(i) == ';'.OR.chrx(i) == ',') THEN
    chr2(indx)=' '
    indx=indx+1
    chr2(indx)='A'
    indx=indx+1
    chr2(indx)='N'
    indx=indx+1
    chr2(indx)='D'
    indx=indx+1
    chr2(indx)=' '
    indx=indx+1
  ELSE
    chr2(indx)=chrx(i)
    indx=indx+1
  END IF
END DO
chr2(indx)='.'

DO  indx2=1,indx
  IF (chr2(indx2) /= ' ') EXIT
END DO

j=1
DO  i=indx2,indx
  IF (i /= indx2 .AND. chrs(j-1) == ' ' .AND. chr2(i) == ' ') CYCLE
  chrs(j)=chr2(i)
  j=j+1
END DO
IF (chrs(1) == '.') GO TO 20

wdx=1
j=1
DO  i=1,100
  IF (chrs(i) == '.') EXIT
  IF (chrs(i) /= ' ') THEN
    IF (j <= 6) txt2(wdx,1)(j:j)=chrs(i)
    IF (j > 6 .AND. j <= 12) txt2(wdx,2)(j-6:j-6)=chrs(i)
    j=j+1
  ELSE
    j=1
    wdx=wdx+1
  END IF
END DO

txt(1,1)=txt2(1,1)
txt(1,2)=txt2(1,2)
j=1
DO  i=2,35
  IF (txt(j,1) /= 'AND   ' .OR. txt2(i,1) /= 'AND   ') THEN
    j=j+1
    txt(j,1)=txt2(i,1)
    txt(j,2)=txt2(i,2)
  ELSE
    wdx=wdx-1
  END IF
END DO
RETURN

!          WRITE(*,12345)(TXT(IQQ,1),IQQ = 1,35)
!12345  FORMAT(' ',5A6)
END SUBROUTINE getlin



!***   GETOBJ

SUBROUTINE getobj (obj)

!  ANALYSE AN OBJECT WORD.  SEE IF THE THING IS HERE, WHETHER WE'VE GOT
!  YET, AND SO ON.  OBJECT MUST BE HERE UNLESS VERB IS "FIND" OR "INVENT
!  (AND NO NEW VERB YET TO BE ANALYSED).  WATER, OIL AND WINE ARE ALSO
!  FUNNY, SINCE THEY ARE NEVER ACTUALLY DROPPED AT ANY LOCATION, BUT MIG
!  BE HERE INSIDE THE BOTTLE OR AS A FEATURE OF THE LOCATION.
!
!  HAS THREE POSSIBLE RETURN VALUES FOR 'OBJ':
!       VAL > 0 :: A POSITIVE OBJECT NUMBER
!       VAL = 0 :: OBJECT NOT FOUND HERE.  ERROR MESSAGE PRINTED.
!       VAL < 0 :: OBJECT WORD REALLY SOMETHING ELSE.  RETURN NEGATIVE
!                       VALUE OF SUBSTITUTED WORD.

INTEGER, INTENT(IN OUT)  :: obj

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: at
CHARACTER (LEN=1) :: zapp(20)
! COMMON /dwfcom/ dwarf,knife,knfloc,dflag,dseen(6),dloc(6),odloc(6) ,dwfmax
! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /mnecom/ back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
!     bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
!     lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
!     troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
!     iwest,phuce(2,4),tk(20)
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150),points(150)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! CHARACTER (LEN=6) :: txt
! COMMON /utxcom/ wdx
! CHARACTER (LEN=6) :: vtxt,otxt,iotxt,dtk,atab
! COMMON /wrdcom/ verbs(45),vrbx,objs(45),objx,iobjs(15),iobx,prep,words(45)
! LOGICAL :: athand,blind,here,holdng,plural
! COMMON /sv3com/ dtk(9),atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)


IF (holdng(obj)) RETURN
IF (blind()) GO TO 100
IF (fixed(obj) == loc .OR. athand(obj)) GO TO 120
IF (.NOT.here(obj)) GO TO 10
k=335
IF (plural(obj)) k=373
obj=0
CALL rspeak (k)
RETURN

10 IF (obj /= grate) GO TO 20
IF (loc == 1 .OR. loc == 4 .OR. loc == 7) obj=-dprssn
IF (loc > 9 .AND. loc < 15) obj=-entrnc
IF (obj == grate) GO TO 100
RETURN

20 IF (obj /= dwarf) GO TO 40
l1=dwfmax-1
DO  i=1,l1
  IF (dloc(i) == loc .AND. dflag >= 2) GO TO 120
END DO
GO TO 100

40 IF (obj == liqloc(loc) .OR. (athand(bottle) .AND. liq(bottle) == obj)  &
       .OR. (athand(cask) .AND. liq(cask) == obj)) GO TO 120
IF (obj /= plant .OR. .NOT.at(plant2) .OR. prop(plant2) == 0) GO TO 50
obj=plant2
GO TO 120

50 IF (obj /= rocks .OR. .NOT.at(carvng)) GO TO 60
obj=carvng
GO TO 120

60 IF (obj /= rod .OR. .NOT.athand(rod2)) GO TO 70
obj=rod2
GO TO 120

70 IF (obj /= door .OR. .NOT.(at(safe) .OR. at(tdoor) .OR. at(tdoor2)  &
       .OR. at(pdoor))) GO TO 80
obj=tdoor
IF (at(tdoor2)) obj=tdoor2
IF (at(pdoor)) obj=pdoor
IF (at(safe)) obj=safe
GO TO 120

80 IF (obj /= book .OR. .NOT.athand(book2)) GO TO 90
obj=book2
GO TO 120

90 IF (verbs(vrbx) == find .OR. verbs(vrbx) == invent) GO TO 120

!  IT ISN'T HERE.  TELL HIM & RETURN.
100 obj=0
CALL a5toa1 (txt(wdx,1), txt(wdx,2), '_here.', zapp, k)
WRITE (*,110) zapp(1:k)
110 FORMAT (/' I see no ', 20A1)

120 RETURN

END SUBROUTINE getobj




!***   GETWDS

SUBROUTINE getwds()

!  WHEN CALLED, CHECKS IF PREVIOUS WORDS VECTOR HAS BEEN EXHAUSTED.
!  IF NOT, BRANCH AROUND THE CODE WHICH READS IN A NEW LINE.  IF VECTOR
!  SUCK UP A LINE FROM THE TTY, THEN CHECK EACH WORD FOR INTELLIGIBILITY
!  IF THE WORD IS VALID, ITS NUMBER GETS STUCK INTO THE WORDS VECTOR.
!  THEN EACH WORD IS PARSED BY THE APPROPRIATE CODE.  THE LABELS BELOW A
!  100 TIMES THE WORD CLASS.
!
!  THE FOLLOWING VECTORS ARE USED:
!       TXT(WDX,2)              HOLD THE RAW TEXT FROM GETLIN
!       WORDS(WDX)              LIST OF WORD NUMBERS, CONVERTED FROM TXT
!       VTXT(VRBX,2)            HOLD THE TEXT FOR VERB VRBX.
!       VERBS(VRBX)             IS THE LIST OF VALIDATED VERB NUMBERS.
!       OTXT(OBJX,2)            HOLDS THE TEXT OF THE OBJECT OBJX.
!       OBJX(OBJX)              IS THE LIST OF VALIDATED OBJECT NUMBERS.
!       IOTXT(IOBX,2)           HOLDS THE TEXT FOR PREP'S IOBJ.

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: blind,living,pflag,hinged,at,athand,toting,k1
! LOGICAL  :: killed
CHARACTER (LEN=6) :: word1,word2,dkk,dk
CHARACTER (LEN=1) :: zapp(20)
LOGICAL  :: k1, pflag

! COMMON /izwiz/ iswiz
! COMMON /adjcom/ adjkey(50),adjtab(150),adjsiz
! COMMON /diecom/ numdie,maxdie,turns,killed
! COMMON /dwfcom/ dwarf,knife,knfloc,dflag,dseen(6),dloc(6),odloc(6) ,dwfmax
! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /mnecom/ back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
!     bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
!     lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
!     troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
!     iwest,phuce(2,4),tk(20)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! COMMON /prpcom/ vkey(60),ptab(300),vkysiz,ptbsiz
! CHARACTER (LEN=6) :: txt,dtk(9),atab
! COMMON /utxcom/ wdx
! CHARACTER (LEN=6) :: vtxt,otxt,iotxt
! COMMON /wrdcom/ verbs(45),vrbx,objs(45),objx,iobjs(15),iobx,prep,words(45)
! COMMON /sv3com/ dtk,atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)

!            MOTION  NOUN  ACTION  MISC   PREP   ADJ    CONJ
INTEGER  :: classd=1, classn=2, classa=3, classm=4, classc=7
!
!  A FEW MORE ACTION VERBS NOT IN COMMON.
INTEGER  :: drop=2, feed=21, nothng=5, light=7, drink=15, score=24
INTEGER  :: pick=41, put=42, get=44
!
!  AND A MOTION VERB:
INTEGER  :: enter=3
!
!  AND TO GET/DROP EVERYTHING IN SIGHT:
INTEGER  :: all=109
!
!  TAKDIR IS A LIST OF MOTION VERBS WHICH ARE ACCEPTABLE AFTER 'TAKE'.
INTEGER  :: takdir(20) = (/ 2,6,9,10,11,13,14,17,23,25,33,34,36,37,39,  &
                            78,79,80,89,-1 /)
!
!  IF WORDS(WDX+1) HAS SOMETHING IN IT, WE ARE STILL PROCESSING OLD INPU
!  IF WORDS(1) = -2, SOMEONE ELSE HAS CALLED GETLIN (E.G., KILL DRAGON).

INTEGER  :: i, vrbkey, wclass, word

pour=13
IF (words(wdx+1) > 0) GO TO 30

!  IF WORDS(1) HAS BEEN SET TO -2, SOMEONE ELSE HAS ALREADY READ IN
!  THE NEW LINE, PRESUMABLY TO CHECK FOR SOME NON-STANDARD WORD.
!  (THIS HAPPENS WHEN KILLING DRAGON WITH BARE HANDS.)
10 IF (words(1) /= -2) CALL getlin()
wdx=0
DO  i=1,35
  words(i)=0
  IF (txt(i,1) /= '      ') words(i)=vocabx(txt(i,1), -1)
END DO

!  THE FIRST WORD OF EACH CLAUSE GETS SPECIAL CHECKING, MOSTLY LOOKING
!  FOR IDIOMS.

!  SPECIAL HANDLING FOR 'SAY' AND 'CALL'.  WIN IF SAYING/CALLING
!  MAGIC WORDS.  NARKY IF UTTERING ANYTHING ELSE.  IF NO OBJ, PASS
!  ON FOR HIGHER LEVEL PARSING.
30 pflag=.false.
wdx=wdx+1
word=words(wdx)
IF (word < 0.0) THEN
  GO TO   520
ELSE IF (word == 0.0) THEN
  GO TO   490
END IF

kk=class(word)
IF (kk == -1) GO TO 130
IF (kk == classa .OR. kk == classd .OR. kk == classm) CALL clrlin()
k=val(word)
IF (kk /= classa .OR. (k /= say .AND. k /= yell)) GO TO 60

!  'SAY' OR 'CALL'.  IF NO NEXT WORD, PASS ON TO HIGHER POWERS.
!  IF OBJECT IS MAGIC WORD ('SAY XYZZY'), FLUSH 'SAY' & TAKE NEXT WORD.
IF (words(wdx+1) == 0 .OR. class(words(wdx+1)) == classc) GO TO 160
wdx=wdx+1
IF (k == say) CALL a5toa1 (txt(wdx,1), txt(wdx,2), '".    ', zapp, k)
IF (k == yell) CALL a5toa1 (txt(wdx,1), txt(wdx,2), '"!!!!!', zapp, k)
word=words(wdx)
IF (word == 62 .OR. word == 65 .OR. word == 71 .OR. word == 82 .OR.  &
    word == 2025) GO TO 160
WRITE (*,50) (zapp(i),i=1,k)
50 FORMAT (/' Okay, "',20A1)
GO TO 560

!  SPECIAL STUFF FOR 'ENTER'.  CAN'T GO INTO WATER.
!  'ENTER BOAT' MEANS 'TAKE BOAT'.
60 word1=txt(wdx,1)
word2=txt(wdx+1,1)
IF (word1 /= 'ENTER ') GO TO 70
IF (word2 == '      ') GO TO 130
spk=43
IF (liqloc(loc) == water) spk=70
IF (word2 == 'STREAM' .OR. word2 == 'WATER ' .OR. word2 == 'RESERV'  &
     .OR. word2 == 'OCEAN ' .OR. word2 == 'SEA   ' .OR. word2 == 'POOL  ') GO TO 500
IF (word2 /= 'BOAT  ' .AND. word2 /= 'ROWBOA') GO TO 160
word=take+2000
GO TO 160

!  'LEAVE' IS A MOTION VERB, UNLESS LEAVING AN OBJECT,
!  E.G., 'LEAVE BOAT' OR 'LEAVE BOTTLE'.  BUT MAKE SURE TO LEAVE ('DROP'
!  ONLY TOTABLE OBJECTS.
70 kk=words(wdx+1)
IF (word1 /= 'LEAVE ' .OR. class(kk) /= classn) GO TO 80
IF (hinged(val(kk)) .OR. fixed(val(kk)) /= 0) GO TO 160
word=leave+2000
GO TO 160

!  IF 'LIGHT LAMP', LIGHT MUST BE TAKEN AS AN ACTION VERB, NOT A NOUN.
80 IF (word1 /= 'LIGHT' .OR. words(wdx+1) /= (lamp+1000)) GO TO 90
word=light+2000
GO TO 160

!  'WATER PLANT' BECOMES 'POUR WATER', IF WE ARE AT PLANT.
!  'OIL DOOR' BECOMES 'POUR OIL', ETC., ETC.
90 IF ((word1 /= 'WATER ' .AND. word1 /= 'OIL   ') .OR. (word2 /= 'PLANT'  &
       .AND. word2 /= 'DOOR  ' .AND. word2 /= 'SWORD ' .AND.  &
       word2 /= 'ANVIL ')) GO TO 110
IF (.NOT.at(vocabx(word2,classn))) GO TO 100
words(wdx+1)=words(wdx)
txt(wdx+1,1)=word1
txt(wdx+1,2)=txt(wdx,2)
100 word=pour+2000
GO TO 160

!  CHECK FILLING OR EMPTYING A CONTAINER.
110 IF ((word1 /= 'EMPTY ') .OR. (class(words(wdx+1)) /= classn)) GO TO 130
kk=val(words(wdx+1))
!      IF(KK.NE.SACK .AND. KK.NE.SAFE .AND. KK.NE.BOAT .AND. KK.NE.CHEST)
!    1  GOTO 91
!  *** UNFINISHED CODE HERE ***
!    ALL THAT ACTUALLY HAPPENS IS OFF ERROR MESSAGES. THE STOOGE
!    SIMPLY CAN'T SAY 'EMPTY SACK OR 'TAKE ALL FROM SACK' ETC
GO TO 130


!  THIS IS THE 'INNER' LOOP.  DISPATCHING OF ALL WORDS IN A CLAUSE AFTER
!  THE FIRST COMES THRU HERE.

120 wdx=wdx+1
word=words(wdx)
130 IF (word < 0.0) THEN
  GO TO   520
ELSE IF (word == 0.0) THEN
  GO TO   580
ELSE
  GO TO   150
END IF

140 wclass=wclass+1
word=vocabx(txt(wdx,1),-(wclass+1))
IF (word == -1) GO TO 490
words(wdx)=word
150 IF (class(word) /= classn) GO TO 160

!  IT'S NOT THE FIRST: MAKE SURE HE INCLUDED A COMMA OR 'AND'.
!  DIFFERENTIATE BETWEEN DIR & INDIR OBJECTS.
!  CHECK FOR SPECIAL CASE OF MULTIPLE OBJECTS: 'FEED BEAR HONEY' OR
!  'THROW TROLL NUGGET'.
k=objx
IF (pflag) k=iobx
IF (k == 0 .OR. class(words(wdx-1)) == classc) GO TO 160
kk=val(verbs(vrbx))
IF (.NOT.living(objs(objx)) .OR. (kk /= throw .AND. kk /= feed)) GO TO 490
iobx=iobx+1
iobjs(iobx)=objs(objx)
objs(objx)=0
objx=objx-1

160 wclass=class(word)
SELECT CASE ( wclass )
  CASE (    1)
    GO TO 170
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 280
  CASE (    4)
    GO TO 320
  CASE (    5)
    GO TO 330
  CASE (    6)
    GO TO 400
  CASE (    7)
    GO TO 450
END SELECT
!            MWD OBJ AVB MVB PRP ADJ CNJ
CALL bug (33)
!  MOTION VERB.
!  A MOTION VERB IS EITHER A DIRECTION ('WEST') OR A MOTION ('JUMP').
!  MULTIPLE MOTIONS MUST BE SEPARATED BY COMMAS OR AND'S.  THERE ARE
!  SOME IDIOMATIC USES WHICH MUST BE SCANNED FOR, SUCH AS 'TAKE BRIDGE',
!  WHICH BECOMES 'BRIDGE' AND 'GO WEST', WHICH BECOMES 'WEST', AND 'LEAV
!  IS DIFFERENT FROM JUST 'LEAVE'.
!
!  IF ORIGINAL VERB WAS 'GO', FLUSH IT & REPLACE WITH THIS ONE.
!  I.E., 'GO WEST' BECOMES 'WEST'.
!
!  CHECK TAKDIR(20) LIST FOR VALID OBJECT MOTION VERBS FOR 'TAKE'.
!  IF FOUND, THROW AWAY 'TAKE' AND USE THE MOTION VERB.
!
!  SINCE THE ORIGINAL VERB IS AN ACTION VERB, CHECK THIS WORD IN THE
!  NOUN TABLE.  MAYBE IT IS AN OBJECT SYNONYMOUS WITH A VERB ('ROCKS').
!
!  IF IT ISN'T A VALID MOTION-OBJECT OF 'TAKE' OR 'GO', NOR AN OBJECT,
!  CHECK THE PREP TABLE.  IF FOUND, HAND IT TO THE PREPOSITION ANALYZER.

170 IF (vrbx == 0) GO TO 200
k=verbs(vrbx)
IF (class(k) > classa) GO TO 490
IF (class(k) /= classa) GO TO 190
IF (val(k) == GO) GO TO 200

IF (val(k) /= take) GO TO 140
kk=val(word)
DO  i=1,20
  IF (takdir(i) == kk) GO TO 200
END DO
GO TO 140

!  IF ORIGINAL MOTION VERB WAS CRAWL, JUMP OR CLIMB, IGNORE CURRENT WORD
!  I.E., 'CLIMB UP' OR 'JUMP OVER' BECOME 'CLIMB' & 'JUMP' ONLY.
190 IF (k == 17 .OR. k == 39 .OR. k == 56) GO TO 120
!            'CRAWL'     'JUMP'    'CLIMB'

200 verbs(1)=word
vrbx=1
IF (txt(wdx,1) /= 'WEST  ') GO TO 120
iwest=iwest+1
IF (iwest == 10) CALL rspeak (17)
k=val(word)
IF (k == EXIT .OR. k == enter) GO TO 560
GO TO 120
!  ANALYZE OBJECT.
!  IF PFLAG IS TRUE, THEN WE ARE PROCESSING A SET OF INDIRECT (PREP)
!  OBJECTS, NOT DIRECT OBJS.

210 IF (pflag) GO TO 340
IF (vrbx /= 0) GO TO 220
k=vocabx(txt(wdx,1),-(classa+1))
IF (k == -1) GO TO 220
word=k
GO TO 280

220 word=val(word)
IF (word == all) GO TO 250
CALL getobj (word)
IF (word < 0.0) THEN
  GO TO   230
ELSE IF (word == 0.0) THEN
  GO TO   560
ELSE
  GO TO   240
END IF

!  IT WASN'T REALLY AN OBJECT.  GO SEE WHAT IT WAS.
230 word=-word
GO TO 160

!  IT WAS REALLY AN OBJECT & IT IS HERE.
240 objx=objx+1
objs(objx)=word
otxt(objx,1)=txt(wdx,1)
otxt(objx,2)=txt(wdx,2)
GO TO 120

!  TAKE EVERYTHING NOT BATTENED DOWN.

250 kk=val(verbs(vrbx))
k1=.false.
IF (kk == drop .OR. kk == put .OR. kk == leave) GO TO 260
k1=.true.
IF (kk /= take .AND. kk /= pick .AND. kk /= get) GO TO 490
spk=357
IF (blind()) GO TO 500
260 DO  i=1,maxobj
  IF (.NOT.athand(i) .OR. fixed(i) /= 0) CYCLE
  IF (i >= water .AND. i <= wine+1) CYCLE
  IF ((k1 .AND. toting(i)) .OR. (.NOT.k1 .AND. .NOT.toting(i))) CYCLE
  objx=objx+1
  objs(objx)=i
!       OTXT(OBJX,1)=NTXT(I,1)
!       OTXT(OBJX,2)=NTXT(I,2)
  otxt(objx,1)='BUG???'
  otxt(objx,2)='      '
  IF (objx == 44) GO TO 120
END DO
GO TO 120

!  ACTION VERB.
280 IF (vrbx == 0) GO TO 300
IF (val(verbs(vrbx)) /= take) GO TO 290
k=val(word)
IF (k == drink .OR. k == invent .OR. k == score .OR. k == nothng .OR.  &
    k == look) GO TO 310
IF (k /= GO) GO TO 490
dk=txt(wdx,1)
IF (dk == 'WALK  ' .OR. dk == 'RUN   ' .OR. dk == 'HIKE  ') GO TO 310
GO TO 490

290 IF (objx /= 0 .OR. class(words(wdx-1)) /= classc) GO TO 490
300 vrbx=vrbx+1
310 verbs(vrbx)=word
vtxt(vrbx,1)=txt(wdx,1)
vtxt(vrbx,2)=txt(wdx,2)
GO TO 120


!  MISCELLANEOUS WORDS/VERBS.
320 IF (vrbx /= 0) GO TO 490
verbs(1)=word
vrbx=1
GO TO 120
!  ANALYZE A PREPOSITION AND ITS OBJECT.  CHECK THAT PREP
!  IS VALID FOR THIS VERB, AND THEN CHECK THAT THE OBJECT IS VALID
!  FOR THIS PREPOSITION.  IF FIRST CHECK FAILS, SYNTAX IS MESSED
!  UP; IF SECOND PART FAILS, IT MAY MERELY BE AN IMPOSSIBLE ACT.

330 IF (class(verbs(vrbx)) /= classa .OR. iobx /= 0) GO TO 490
IF (pflag) GO TO 340
vrbkey=vkey(val(verbs(vrbx)))
IF (vrbkey == 0) GO TO 490
prep=val(word)
pflag=.true.
wdx=wdx+1
word=words(wdx)
IF (word == 0) GO TO 360
SELECT CASE ( class(word) )
  CASE (    1)
    GO TO 490
  CASE (    2)
    GO TO 340
  CASE (    3)
    GO TO 490
  CASE (    4)
    GO TO 490
  CASE (    5)
    GO TO 490
  CASE (    6)
    GO TO 400
  CASE (    7)
    GO TO 360
END SELECT
GO TO 520

340 word=val(word)
IF (word == all) GO TO 360
CALL getobj (word)
IF (word < 0.0) THEN
  GO TO   390
ELSE IF (word == 0.0) THEN
  GO TO   560
END IF

iobx=iobx+1
iobjs(iobx)=word
iotxt(iobx,1)=txt(wdx,1)
iotxt(iobx,2)=txt(wdx,2)
360 kk=ABS(ptab(vrbkey)/1000)
IF (kk /= prep) GO TO 370

!  PREP IS VALID WITH THIS VERB.  NOW CHECK OBJECT OF PREP.
IF (word == 0 .OR. class(word) == classc) GO TO 380

!  AN OBJ FOLLOWS THE PREP.  SEE IF IT'S PLAUSIBLE.
kk=ABS((MOD(ptab(vrbkey),0001000)))
IF (kk == word .AND. kk == all) GO TO 250
IF (kk == word .OR. kk == 999) GO TO 120
370 vrbkey=vrbkey+1
IF (ptab(vrbkey-1) >= 0) GO TO 360
GO TO 390

!  NO OBJ FOLLOWS PREP.  CHECK SPECIAL CASES.
380 pflag=.false.
wdx=wdx-1
dk=txt(wdx,1)
dkk=vtxt(vrbx,1)
IF ((dk /= 'ON    ' .AND. dk /= 'OFF   ') .AND. (dkk /= 'TURN  '  &
     .OR. objs(objx) /= lamp) .AND. (dkk /= 'TAKE  ' .AND. dkk /= 'PUT   ')) GO TO 390
IF ((dk == 'UP    ' .AND. dkk /= 'PICK  ') .OR. (dk == 'DOWN  '  &
     .AND. (dkk /= 'PUT   ' .AND. verbs(vrbx) /= throw))) GO TO 390
wdx=wdx+1
word=words(wdx)
IF (word == 0) GO TO 580
IF (class(word) /= classc) GO TO 490
GO TO 130

!  YOU CAN'T DO THAT!!
390 spk=noway()
GO TO 500
!  ADJECTIVE HANDLER.
!  SCARF THE NEXT WORD, MAKE SURE IT IS A VALID OBJECT FOR THIS ADJ.
!  THEN CALL GETOBJ TO SEE IF IT IS REALLY THERE, THEN LINK INTO OBJ CODE.

400 adj=val(word)
wdx=wdx+1
word=words(wdx)
IF (word < 0.0) THEN
  GO TO   520
ELSE IF (word == 0.0) THEN
  GO TO   430
END IF

IF (class(word) == classc) GO TO 430
IF (class(word) /= classn) word=vocabx(txt(wdx,1),-(classn+1))
IF (word == -1 .OR. class(word) /= classn .OR. val(word) == all) GO TO 490
words(wdx)=word
kk=val(word)
k=adjkey(adj)
420 IF (kk == ABS(adjtab(k))) GO TO 150
IF (adjtab(k) < 0) GO TO 490
k=k+1
GO TO 420

430 CALL a5toa1 (txt(wdx-1,1), txt(wdx-1,2), '_WHAT?', zapp, k)
WRITE (*,440) zapp(1:k)
440 FORMAT (' ', 20A1)
GO TO 10


!  ANALYZE A CONJUNCTION.  MAY BE A COMMA OR AN EXPLICIT "AND".
!  LOOK AHEAD AT NEXT WORD.  IF IT IS AN ACTION VERB AND NO OBJECT
!  HAS YET BEEN SPECIFIED, PUT IT INTO THE VERB STACK.  IF IT IS
!  AN OBJECT, ADD IT TO THE PILE.
!  ELSE, BUMP BACK THE WORD POINTER, ASSUME END OF CLAUSE, AND RETURN.
450 wdx=wdx+1
word=words(wdx)
IF (word < 0.0) THEN
  GO TO   520
ELSE IF (word == 0.0) THEN
  GO TO   490
END IF

SELECT CASE ( class(word) )
  CASE (    1)
    GO TO 480
  CASE (    2)
    GO TO 150
  CASE (    3)
    GO TO 470
  CASE (    4)
    GO TO 480
  CASE (    5)
    GO TO 490
  CASE (    6)
    GO TO 150
  CASE (    7)
    GO TO 490
END SELECT

!  A NEW ACTION VERB FOLLOWS.  IF NO PREVIOUS VERB HAS BEEN TYPED,
!  HE LOSES.  IF PREVIOUS VERB IS NOT AN ACTION VERB, HE LOSES.
!  IF AN OBJ/IOBJ WAS SPECIFIED FOR PREV ACT VERB, HE LOSES.  ONLY
!  VALID SYNTAX IS: 'GET AND OPEN CAGE'.
470 IF (vrbx /= 0 .AND. class(verbs(vrbx))  &
     == classa .AND. objx == 0 .AND. iobx == 0) GO TO 150
480 wdx=wdx-1
GO TO 580
!  GEE, I DON'T UNDERSTAND.  FLUSH REST OF CURRENT CLAUSE, UP TO
!  EOL OR CONJUNCTION & CONTINUE.
490 spk=confuz()
500 CALL rspeak (spk)
CALL clrlin()
GO TO 10

!  AN IRREGULAR WORD WAS TYPED IN BY USER.  CHECK FOR WIZARDRY.

520 IF (pct(25)) GO TO 540
CALL a5toa1 (txt(wdx,1), txt(wdx,2), '.     ', zapp, k)
WRITE (*,530) (zapp(i),i=1,k)
530 FORMAT (/' I don''T UNDERSTAND THE WORD ', 20A1)
CALL clrlin()
GO TO 10

540 CALL a5toa1 (txt(wdx,1), txt(wdx,2), '?     ', zapp, k)
WRITE (*,550) (zapp(i),i=1,k)
550 FORMAT (/' Mumble?  ', 20A1)
CALL clrlin()
GO TO 10


!  SCAN TO CONJ OR END OF LINE.
560 CALL clrlin()
pflag=.false.
570 wdx=wdx+1
IF (words(wdx) == 0) GO TO 10
IF (class(words(wdx)) == classc) GO TO 120
GO TO 570

!  END OF CLAUSE.  WE APPEAR TO HAVE REACHED THE END OF A SENTENCE.
!  IT WAS TERMINATED EITHER BY CRLF OR A CONJUNCTION.  IF A CONJ,
!  THE CONJ ANALYZER CLAIMS  THAT THE NEXT WORDS ARE NOT PART OF
!  THIS CLAUSE.  DECIDE WHETHER OR NOT WE HAVE ENOUGH TO WORK WITH.
580 pflag=.false.
IF (verbs(1) /= 0) GO TO 610
IF (objs(1) == 0) GO TO 490
IF (objs(2) /= 0) GO TO 600
CALL a5toa1 (otxt(1,1), otxt(1,2), '?     ', zapp, k)
WRITE (*,590) (zapp(i),i=1,k)
590 FORMAT (/' What do you want to do with the ',20A1)
GO TO 10

600 WRITE (*,*) ' What do you want to do with them'
GO TO 10

610 IF (objx > 1 .AND. iobx > 1) GO TO 490
RETURN

END SUBROUTINE getwds



!***  HERE   .TRUE. IF OBJ AT THIS LOCATION

FUNCTION here(obj) RESULT(lval)

!  HERE(OBJ)    = TRUE IF THE OBJ IS AT "LOC" (OR IS BEING CARRIED)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! LOGICAL :: toting

lval=.false.
IF (obj < 1 .OR. obj > maxobj) RETURN
lval=place(obj) == loc .OR. toting(obj)
RETURN
END FUNCTION here



!***  HINGED .TRUE. IF OBJ CAN BE OPENED

FUNCTION hinged(obj) RESULT(lval)

!  HINGED(OBJ)  = TRUE IF OBJECT CAN BE OPENED/SHUT.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),1)

RETURN
END FUNCTION hinged




!***  HOLDNG .TRUE. IF HOLDING OBJ

FUNCTION holdng(obj) RESULT(lval)

!  HOLDNG(OBJ)  = TRUE IF THE OBJ IS BEING CARRIED IN HAND.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

lval=.false.
IF (obj < 1 .OR. obj > maxobj) RETURN
lval=place(obj) == -1
RETURN
END FUNCTION holdng




!***   INSERT

SUBROUTINE insert (object, contnr)

INTEGER, INTENT(IN)  :: object
INTEGER, INTENT(IN)  :: contnr

! IMPLICIT INTEGER (a-z)
! COMMON /hldcom/ holder(150),hlink(150)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

INTEGER  :: temp

IF (contnr == object) CALL bug (32)
CALL carry (object, loc)

temp=holder(contnr)
holder(contnr)=object
hlink(object)=temp
place(object)=-contnr
RETURN

END SUBROUTINE insert



!***   INSIDE .TRUE. IF LOCATION IS WELL WITHIN THE CAVE


FUNCTION inside(loc) RESULT(lval)

!  INSIDE(LOC)  = TRUE IF LOCATION IS WELL WITHIN THE CAVE

INTEGER, INTENT(IN)  :: loc
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: outsid, cave_entry

lval=.NOT.outsid(loc) .AND. .NOT.cave_entry(loc)
RETURN
END FUNCTION inside




!***   JUGGLE

SUBROUTINE juggle (object)

!  JUGGLE AN OBJECT BY PICKING IT UP AND PUTTING IT DOWN AGAIN, THE PURP
!  BEING TO GET THE OBJECT TO THE FRONT OF THE CHAIN OF THINGS AT ITS LO

INTEGER, INTENT(IN OUT)  :: object

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

i=place(object)
j=fixed(object)
CALL move (object, i)
CALL move (object+maxobj, j)
RETURN
END SUBROUTINE juggle



!***   LIQ

FUNCTION liq(obj) RESULT(ival)

INTEGER, INTENT(IN)  :: obj
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150), points(150)

INTEGER  :: iq

!       LIQ=LIQ2(MAX0(PROP(OBJ),-1-PROP(OBJ)))
ival=0
IF (obj /= bottle .AND. obj /= cask) RETURN
iq=MAX(prop(obj)+1, -1-(prop(obj)+1))
IF (iq <= 0) RETURN
ival=liqtyp(iq)
RETURN
END FUNCTION liq




!***   LIQ2
!  NON-LOGICAL (ILLOGICAL?) FUNCTIONS (CLASS,LIQ,LIQ2,LIQLOC,VAL)

FUNCTION liq2(pbotl) RESULT(ival)

INTEGER, INTENT(IN)  :: pbotl
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)

ival=(1-pbotl)*water + (pbotl/2)*(water+oil) + (pbotl/4)*(water+wine-2*oil)
RETURN
END FUNCTION liq2




!***   LIQLOC

FUNCTION liqloc(loc) RESULT(ival)

INTEGER, INTENT(IN)  :: loc
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)
! INTEGER :: wrd(2)
! COMMON /concom/ loccon(250),objcon(150)
! EQUIVALENCE (loccon,wrd)
!      CALL TOOCT(LOCCON(LOC))
!      CALL TOOCT(WRD(LOC*2))

ival=liq2((MOD(loccon(loc)/8,2)*(MOD(loccon(loc)/2*2,16)-9)+1))

RETURN
END FUNCTION liqloc



!***  LIVING .TRUE. IF OBJ IS LIVING, BEAR FOR EXAMPLE

FUNCTION living(obj) RESULT(lval)

!  LIVING(OBJ)  = TRUE IF OBJ IS SOME SORT OF CRITTER

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),9)

RETURN
END FUNCTION living




!***   LOCKED  .TRUE. IF LOCKABLE OBJ IS LOCKED

FUNCTION locked(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),4)
RETURN
END FUNCTION locked



!***   LOCKS  .TRUE. IF YOU CAN LOCK THIS OBJ

FUNCTION locks(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),3)
RETURN
END FUNCTION locks




!***   LOOKIN

SUBROUTINE lookin (contnr)

!  LIST CONTENTS IF OBJ IS A CONTAINER AND IS OPEN OR TRANSPARENT.
!  SAVE INITIAL VALUE OF BLKLIN THRU SUBROUTINE.

INTEGER, INTENT(IN)  :: contnr

! IMPLICIT INTEGER (a-z)
! COMMON /blkcom/ blklin
! COMMON /hldcom/ holder(150),hlink(150)
! LOGICAL  :: vessel,ajar,opaque,blklin,bsave
INTEGER  :: loop, temp
LOGICAL  :: bsave

IF (.NOT.vessel(contnr) .OR. (.NOT.ajar(contnr) .AND. opaque(contnr))) RETURN
temp=holder(contnr)
loop=0
bsave=blklin
10 IF (temp == 0) RETURN
blklin=.false.
IF (loop == 0) CALL rspeak (360)
CALL tnoua()
CALL pspeak (temp, -1)
blklin=bsave
temp=hlink(temp)
loop=-1
GO TO 10

END SUBROUTINE lookin



SUBROUTINE move (object, where)

!  PLACE ANY OBJECT ANYWHERE BY PICKING IT UP AND DROPPING IT.  MAY ALRE
!  TOTING, IN WHICH CASE THE CARRY IS A NO-OP.  MUSTN'T PICK UP OBJECTS
!  ARE NOT AT ANY LOC, SINCE CARRY WANTS TO REMOVE OBJECTS FROM ATLOC CH

INTEGER, INTENT(IN)  :: object
INTEGER, INTENT(IN)  :: where

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! LOGICAL :: enclsd

INTEGER  :: from

IF (object > maxobj) THEN
  from=fixed(object-maxobj)
ELSE
  IF (enclsd(object)) CALL remove (object)
  from=place(object)
END IF
IF (from > 0 .AND. from <= maxobj*2) CALL carry (object, from)
CALL drop (object, where)
RETURN
END SUBROUTINE move




!***   NOWAY

FUNCTION noway() RESULT(ival)

!  GENERATE'S SOME VARIANT OF "CAN'T DO THAT" MESSAGE.

INTEGER  :: ival

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: pct

ival=14
IF (pct(50)) ival=110
IF (pct(33)) ival=147
IF (pct(25)) ival=250
IF (pct(20)) ival=262
IF (pct(17)) ival=25
IF (pct(14)) ival=345
IF (pct(12)) ival=346
RETURN
END FUNCTION noway




!***  OPAQUE .TRUE. IF OBJ IS NON-TRANSPARENT CONTAINER

FUNCTION opaque(obj) RESULT(lval)

!  OPAQUE(OBJ)  = TRUE IF OBJECT IS NOT TRANSPARENT.  E.G., BAG & CHEST
!                 WICKER CAGE & GLASS BOTTLE ARE TRANSPARENT.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),6)

RETURN
END FUNCTION opaque



!***   OUTSID .TRUE. IF LOCATION IS OUTSIDE THE CAVE

FUNCTION outsid(loc) RESULT(lval)

!  OUTSID(LOC)  = TRUE IF LOCATION IS OUTSIDE THE CAVE

INTEGER, INTENT(IN)  :: loc
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(loccon(loc),6)

RETURN
END FUNCTION outsid




!***   PCT

FUNCTION pct(n) RESULT(lval)

!  PCT(N)       = TRUE N% OF THE TIME (N INTEGER FROM 0 TO 100)

INTEGER, INTENT(IN)  :: n
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
lval=ranz(100) < n
RETURN
END FUNCTION pct




!***  PLURAL .TRUE. IF OBJ IS MULTIPLE OBJS

FUNCTION plural(obj) RESULT(lval)

!  PLURAL(OBJ)  = TRUE IF OBJECT IS A "BUNCH" OF THINGS (COINS, SHOES).

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL  :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),13)

RETURN
END FUNCTION plural




!***   cave_entry .TRUE. IF LOCATION IS IN CAVE ENTRANCE
! N.B. Name changed from PORTAL as there is a variable of that name in
!      COMMON /savcom/

FUNCTION cave_entry(loc) RESULT(lval)

!  cave_entry(LOC)  = TRUE IS LOCATION IS IN CAVE "ENTRANCE"

INTEGER, INTENT(IN)  :: loc
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(loccon(loc),5)

RETURN
END FUNCTION cave_entry




!***  PRINTD .TRUE. IF OBJ CAN BE READ

FUNCTION printd(obj) RESULT(lval)

!  PRINTD(OBJ)  = TRUE IF OBJECT CAN BE READ.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),8)

RETURN
END FUNCTION printd




!***   PSPEAK

SUBROUTINE pspeak (msg, skip)

!  FIND THE SKIP+1ST MESSAGE FROM MSG AND PRINT IT.  MSG SHOULD BE THE I
!  THE INVENTORY MESSAGE FOR OBJECT.  (INVEN+N+1 MESSAGE IS PROP=N MESSA

INTEGER, INTENT(IN)  :: msg
INTEGER, INTENT(IN)  :: skip

! IMPLICIT INTEGER (a-z)
! COMMON /txtcom/ lines(25000),rtext(450),ptext(150)

INTEGER  :: m

m=ptext(msg)
IF (skip < 0) GO TO 30
DO  i=0,skip
  10 m=ABS(lines(m))
  IF (lines(m) >= 0) GO TO 10
END DO
30 CALL speak (m)
RETURN
END SUBROUTINE pspeak




!***   PUT

FUNCTION put(object,where,pval) RESULT(ival)

!  PUT IS THE SAME AS MOVE, EXCEPT IT RETURNS A VALUE USED TO SET UP THE
!  NEGATED PROP VALUES FOR THE REPOSITORY OBJECTS.

INTEGER, INTENT(IN)  :: object
INTEGER, INTENT(IN)  :: where
INTEGER, INTENT(IN)  :: pval
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)

CALL move (object, where)
ival=-1-pval
RETURN
END FUNCTION put




!***   RANZ
!  UTILITY ROUTINES (SHIFT, RAN, DATIME, CIAO, BUG, LOG)

FUNCTION ranz(range) RESULT(ival)

INTEGER, INTENT(IN)  :: range
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
INTEGER, SAVE  :: seed = 12345

seed=seed*69069+1
i=16384/range
j=ABS(seed)/i
ival=MOD(j,range)
RETURN
END FUNCTION ranz



!***   RATING

SUBROUTINE rating (score, bonus, gaveup, scorng, closng, closed, hntmax)

!  CALCULATE WHAT THE PLAYER'S SCORE WOULD BE IF HE QUIT NOW.
!  THIS MAY BE THE END OF THE GAME, OR HE MAY JUST BE WONDERING
!  HOW HE IS DOING.

INTEGER, INTENT(OUT)  :: score
INTEGER, INTENT(IN)   :: bonus
LOGICAL, INTENT(IN)   :: gaveup
LOGICAL, INTENT(IN)   :: scorng
LOGICAL, INTENT(IN)   :: closng
LOGICAL, INTENT(IN)   :: closed
INTEGER, INTENT(IN)   :: hntmax

! IMPLICIT INTEGER (a-z)
! LOGICAL :: treasr,hinted,killed
! COMMON /mnecom/ back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
!     bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
!     lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
!     troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
!     iwest,phuce(2,4),tk(20)
! COMMON /diecom/ numdie,maxdie,turns,killed
! COMMON /dwfcom/ dwarf,knife,knfloc,dflag,dseen(6),dloc(6),odloc(6),dwfmax
! COMMON /hntcom/ hintlc(20),hinted(20),hints(20,4),hntsiz,hntmin
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150),points(150)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! COMMON /sv2com/ anvil,batter,bees,billbd,bird,brush,cage,cakes,  &
!     chain,chest,chloc,chloc2,clam,cloak,clsses,coins,crown,daltlc,dog,  &
!     dragon,eggs,fissur,flower,gatloc,grail,hive,honey,horn,jewels,  &
!     keys,lyre,magzin,mirror,mushrm,mxscor,nugget,oyster,pearl,phone,  &
!     pillow,pole,poster,prepat,prepdn,prepfr,prepin,prepof,prepon,  &
!     pyram,radium,ring,rug,sapphi,shield,shoes,shut,slugs,snake,sphere,  &
!     steps,sticks,sword,tablet,tridnt,unlock,vase,wall,wall2,wear,wumpus, &
!     y2,yank

INTEGER  :: ascore, jturns, qk(20)

!  THE PRESENT SCORING ALGORITHM IS AS FOLLOWS:
!  (TREASURE POINTS ARE EXPLAINED IN A FOLLOWING COMMENT)
!     OBJECTIVE:          POINTS:        PRESENT TOTAL POSSIBLE:
!  GETTING WELL INTO CAVE   25                    25
!  TOTAL POSSIBLE FOR TREASURES (+MAG)           426
!  SURVIVING             (MAX-NUM)*10             30
!  NOT QUITTING              4                     4
!  REACHING "CLOSNG"        20                    20
!  "CLOSED": QUIT/KILLED    10
!            KLUTZED        20
!            WRONG WAY      25
!            SUCCESS        30                    30
!  ROUND OUT THE TOTAL      16                    16
!                                       TOTAL:   551
!  (POINTS CAN ALSO BE DEDUCTED FOR USING HINTS.)

score=0
mxscor=0

!  FIRST TALLY UP THE TREASURES.  MUST BE IN BUILDING AND NOT BROKEN.
!  GIVE THE POOR GUY PARTIAL SCORE JUST FOR FINDING EACH TREASURE.
!  GETS FULL SCORE, QK(3), FOR OBJ IF:
!       OBJ IS AT LOC QK(1), AND
!       OBJ HAS PROP VALUE OF QK(2)
!
!               WEIGHT          TOTAL POSSIBLE
!  MAGAZINE     1 (ABSOLUTE)            1
!
!  ALL THE FOLLOWING ARE MULTIPLIED BY 5 (RANGE 5-25):
!  BOOK         2
!  CASK         3 (WITH WINE ONLY)
!  CHAIN        4 (MUST ENTER VIA STYX)
!  CHEST        5
!  CLOAK        3
!  CLOVER       1
!  COINS        5
!  CROWN        2
!  CRYSTAL-BALL 2
!  DIAMONDS     2
!  EGGS         3
!  EMERALD      3
!  GRAIL        2
!  HORN         2
!  JEWELS       1
!  LYRE         1
!  NUGGET       2
!  PEARL        4
!  PYRAMID      4
!  RADIUM       4
!  RING         4
!  RUG          3
!  SAPPHIRE     1
!  SHOES        3
!  SPICES       1
!  SWORD        4
!  TRIDENT      2
!  VASE         2
!  DROPLET      5
!  TREE         5
!       TOTAL: 85 * 5 = 425 + 1 ==> 426

DO  obj=1,maxobj
  IF (points(obj) == 0) CYCLE
  qk(3)=ABS(points(obj))/1000000
  qk(2)=(ABS(points(obj)) - qk(3)*1000000)/1000
  qk(1)=ABS(points(obj)) - qk(3)*1000000-qk(2)*1000
  IF (points(obj) < 0) qk(1)=-qk(1)
  k=0
  IF (.NOT.treasr(obj)) GO TO 10
  k=qk(3)*2
  IF (prop(obj) >= 0) score=score + k
  qk(3)=qk(3)*5
  10 IF (place(obj) == qk(1) .AND. prop(obj) == qk(2) .AND. (place(obj)  &
         /= -chest .OR. place(chest) == 3) .AND. (place(obj) /= -  &
         shield .OR. place(shield) == -safe)) score=score + qk(3) - k
  mxscor=mxscor + qk(3)
END DO
!  NOW LOOK AT HOW HE FINISHED AND HOW FAR HE GOT.  MAXDIE AND NUMDIE TE
!  HOW WELL HE SURVIVED.  GAVEUP SAYS WHETHER HE EXITED VIA QUIT.  DFLAG
!  TELL US IF HE EVER GOT SUITABLY DEEP INTO THE CAVE.  CLOSNG STILL IND
!  WHETHER HE REACHED THE ENDGAME.  AND IF HE GOT AS FAR AS "CAVE CLOSED
!  (INDICATED BY "CLOSED"), THEN BONUS IS ZERO FOR MUNDANE EXITS OR 133,
!  135 IF HE BLEW IT (SO TO SPEAK).

ascore=(maxdie-numdie)*10
mxscor=mxscor + maxdie*10
IF (.NOT.(scorng .OR. gaveup)) ascore=ascore + 4
mxscor=mxscor + 4
IF (dflag /= 0) ascore=ascore + 25
mxscor=mxscor + 25
IF (closng) ascore=ascore + 20
mxscor=mxscor+20
IF (.NOT.closed) GO TO 30
IF (bonus == 0) ascore=ascore + 10
IF (bonus == 135) ascore=ascore + 20
IF (bonus == 134) ascore=ascore + 25
IF (bonus == 133) ascore=ascore + 30
30 mxscor=mxscor + 30
!  ROUND IT OFF.

ascore=ascore + 16
mxscor=mxscor + 16

!  DEDUCT POINTS FOR HINTS.  HINTS < HNTMIN ARE SPECIAL; SEE DATABASE DE

DO  i=1,hntmax
  IF (hinted(i)) score=score - hints(i,2)
END DO
jturns=turns/100
IF (jturns == 0) ascore=0
IF (jturns == 1) ascore=ascore/3
IF (jturns == 2) ascore=(ascore*2)/3
score=score + ascore
IF (score < 0) score=0
RETURN
END SUBROUTINE rating



!***   REMOVE

SUBROUTINE remove (object)

INTEGER, INTENT(IN)  :: object

! IMPLICIT INTEGER (a-z)
! COMMON /hldcom/ holder(150),hlink(150)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj

INTEGER  :: contnr, temp

contnr=-place(object)
place(object)=-1

IF (holder(contnr) /= object) GO TO 10
holder(contnr)=hlink(object)
RETURN

10 temp=holder(contnr)
20 IF (hlink(temp) == object) GO TO 30
temp=hlink(temp)
GO TO 20

30 hlink(temp)=hlink(object)
RETURN
END SUBROUTINE remove




!***   RSPEAK

SUBROUTINE rspeak (i)

!  PRINT THE I-TH "RANDOM" MESSAGE (SECTION 6 OF DATABASE).

INTEGER, INTENT(IN)  :: i

! IMPLICIT INTEGER (a-z)
! COMMON /txtcom/ lines(25000),rtext(450),ptext(150)

INTEGER  :: m

IF (i <= 0) RETURN
m=rtext(i)
CALL speak (m)
RETURN
END SUBROUTINE rspeak




!***   SMALL  .TRUE. IF IT FITS IN SACK OR SMALL CONTAINER

FUNCTION small(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),5)
RETURN
END FUNCTION small



!***   SPEAK

SUBROUTINE speak (n)

!  PRINT THE MESSAGE WHICH STARTS AT LINES(N).  PRECEDE IT WITH A BLANK
!  UNLESS BLKLIN IS FALSE.

INTEGER, INTENT(IN)  :: n

! IMPLICIT INTEGER (a-z)
! LOGICAL :: blklin
! COMMON /tnoux/ indent
! COMMON /txtcom/ lines(25000),rtext(450),ptext(150)
! COMMON /blkcom/ blklin
INTEGER  :: oline(30)
!     DATA ZCLYD/'CLYD'/,ZLS/'<$$<'/
!        ZCLYD = 0
!       ZCLYD = 37+256*(9+256*(20+256*15))

INTEGER  :: xclyd, zclyd, zls

zclyd=ICHAR('c') + 256*(ICHAR('L') + 256*(ICHAR('y') + 256*ICHAR('D')))
zls=60 + 256*(36 + 256*(36 + 256*60))
xclyd=IEOR(zclyd,zls)
IF (n == 0) GO TO 60
IF (lines(n+1) == xclyd) GO TO 60
k=n
10 l=ABS(lines(k))-k-1
DO  i=1,l
  oline(i)=IEOR(lines(k+i),zclyd)
END DO
IF (indent == 0) WRITE (*,40) (oline(i),i=1,l)
IF (indent == 1) WRITE (*,30) (oline(i),i=1,l)
30 FORMAT (t7, 19A4)
40 FORMAT (' ', 19A4)
k=k+l+1
IF (lines(k) >= 0) GO TO 10

60 indent=0
RETURN
END SUBROUTINE speak



!       A HORRIBLE KLUDGE

SUBROUTINE tnoua()
! INTEGER  :: indent
! COMMON /tnoux/ indent
indent=1
RETURN
END SUBROUTINE tnoua



!***  TOTING .TRUE. IF OBJ SOMEWHERE ON PERSON

FUNCTION toting(obj) RESULT(lval)

!  TOTING(OBJ)  = TRUE IF THE OBJ IS BEING CARRIED (IN HAND OR
!                 CONTAINER).  OBJ MAY NOT BE REACHABLE.  SEE
!                 ALSO: ENCLSD, ATHAND, HOLDNG.

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! LOGICAL :: enclsd,aaa,bbb,ccc

INTEGER  :: contnr, outer, outer2

lval=.false.
IF (holdng(obj)) THEN
  lval=.true.
  RETURN
END IF
contnr=-place(obj)
IF (contnr <= 0) RETURN
IF (holdng(contnr)) THEN
  lval=.true.
  RETURN
END IF
outer=-place(contnr)
IF (outer <= 0) RETURN
IF (holdng(outer)) THEN
  lval=.true.
  RETURN
END IF
outer2=-place(outer)
IF (outer2 <= 0) RETURN
IF (holdng(outer2)) THEN
  lval=.true.
  RETURN
END IF
RETURN

END FUNCTION toting



!***   TRAVL
!  FIGURE OUT THE NEW LOCATION
!
!  GIVEN THE CURRENT LOCATION IN "LOC", AND A MOTION VERB NUMBER IN "K",
!  THE NEW LOCATION IN "NEWLOC".  THE CURRENT LOC IS SAVED IN "OLDLOC" I
!  HE WANTS TO RETREAT.  THE CURRENT OLDLOC IS SAVED IN OLDLC2, IN CASE
!  DIES.  (IF HE DOES, NEWLOC WILL BE LIMBO, AND OLDLOC WILL BE WHAT KIL
!  HIM, SO WE NEED OLDLC2, WHICH IS THE LAST PLACE HE WAS SAFE.)

SUBROUTINE travl (k, bcross, tally2)

INTEGER, INTENT(IN OUT)  :: k
INTEGER, INTENT(OUT)     :: bcross
INTEGER, INTENT(OUT)     :: tally2

! IMPLICIT INTEGER (a-z)

! LOGICAL  :: killed
! LOGICAL  :: inside,outsid,here,at,forced
INTEGER  :: k2, kalflg, new1

! COMMON /diecom/ numdie,maxdie,turns,killed
! COMMON /trvcom/ travel(1600)
! COMMON /ltxcom/ ltext(250),stext(250),key(250),abb(250),locsiz
! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150), points(150)
! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
! COMMON /mnecom/ back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
!     bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
!     lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
!     troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
!     iwest,phuce(2,4),tk(20)

killed=.false.
kk=key(loc)
newloc=loc
IF (kk == 0) CALL bug (26)
IF (k == null) RETURN
IF (k == back) GO TO 80
IF (k == cave) GO TO 160
oldlc2=oldloc
oldloc=loc

10 ll=ABS(travel(kk))
IF (MOD(ll,0001000) == 1 .OR. MOD(ll,0001000) == k) GO TO 20
IF (travel(kk) < 0000000) GO TO 170
kk=kk+1
GO TO 10

20 ll=ll/0001000
30 newloc=ll/0001000
k=MOD(newloc,100)
IF (newloc <= maxloc) GO TO 50
IF (prop(k) /= newloc/100-3) GO TO 70
40 IF (travel(kk) < 0) CALL bug (25)
kk=kk+1
new1=ABS(travel(kk))/1000
IF (new1 == ll) GO TO 40
ll=new1
GO TO 30

50 IF (newloc <= 100) GO TO 60
IF (toting(k) .OR. (newloc > 200 .AND. at(k))) GO TO 70
GO TO 40

60 IF (newloc /= 0 .AND. .NOT.pct(newloc)) GO TO 40
70 newloc=MOD(ll,1000)
IF (newloc <= maxloc) GO TO 120
IF (newloc <= 500) GO TO 180
CALL rspeak (newloc-500)
newloc=loc
RETURN
!  HANDLE "GO BACK".  LOOK FOR VERB WHICH GOES FROM LOC TO OLDLOC, OR TO
!  IF OLDLOC HAS FORCED-MOTION.  K2 SAVES ENTRY -> FORCED LOC -> PREVIOU

80 k=oldloc
IF (forced(k)) k=oldlc2
oldlc2=oldloc
oldloc=loc
k2=0
IF (k /= loc) GO TO 90
CALL rspeak (91)
GO TO 120

90 ll=MOD((ABS(travel(kk))/1000),1000)
IF (ll == k) GO TO 150
IF (ll > maxloc) GO TO 100
j=key(ll)
IF (forced(ll) .AND. MOD((ABS(travel(j))/1000),1000) == k) k2=kk
100 IF (travel(kk) < 000000) GO TO 110
kk=kk+1
GO TO 90

110 kk=k2
IF (kk /= 0) GO TO 150
CALL rspeak (140)
120 IF (newloc < 242 .OR. newloc > 247) RETURN
IF (newloc /= 242) GO TO 130
kalflg=0
RETURN
130 IF (newloc /= oldloc+1) GO TO 140
kalflg=kalflg+1
RETURN
140 kalflg=-10
RETURN

150 k=MOD(ABS(travel(kk)),1000)
kk=key(loc)
GO TO 10

!  CAVE.  DIFFERENT MESSAGES DEPENDING ON WHETHER ABOVE GROUND.

160 IF (outsid(loc)) CALL rspeak (57)
IF (.NOT.outsid(loc)) CALL rspeak (58)
RETURN

!  NON-APPLICABLE MOTION.  VARIOUS MESSAGES DEPENDING ON WORD GIVEN.

170 spk=12
IF (k >= 43 .AND. k <= 50) spk=9
IF (k == 29 .OR. k == 30) spk=9
IF (k == 7 .OR. k == 36 .OR. k == 37) spk=10
IF (k == 11 .OR. k == 19) spk=11
IF (k == 62 .OR. k == 65 .OR. k == 82) spk=42
IF (k == 17) spk=80
CALL rspeak (spk)
RETURN
!  SPECIAL MOTIONS COME HERE.  LABELLING CONVENTION: STATEMENT NUMBERS N
!  (XX=00-99) ARE USED FOR SPECIAL CASE NUMBER NNN (NNN=301-500).

180 newloc=newloc-maxloc
SELECT CASE ( newloc )
  CASE (    1)
    GO TO 200
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 220
  CASE (    4)
    GO TO 240
  CASE (    5)
    GO TO 260
  CASE (    6)
    GO TO 290
  CASE (    7)
    GO TO 320
END SELECT
!             ALCOV PLOVR TROLL PHUCE BOOTH BRDGE
WRITE (*,190) newloc
190 FORMAT ('BUG IN TRAVEL TABLES.  NEWLOC= ',i5)
CALL bug (20)

!  TRAVEL 301.  PLOVER-ALCOVE PASSAGE.  CAN CARRY ONLY EMERALD.  NOTE: T
!  TABLE MUST INCLUDE "USELESS" ENTRIES GOING THROUGH PASSAGE, WHICH CAN
!  BE USED FOR ACTUAL MOTION, BUT CAN BE SPOTTED BY "GO BACK".

200 newloc=99+100-loc
kk=burden(0)
IF (kk == 0 .OR. (kk == burden(emrald) .AND. holdng(emrald))) RETURN
newloc=loc
CALL rspeak (117)
RETURN

!  TRAVEL 302.  PLOVER TRANSPORT.  DROP THE EMERALD (ONLY USE SPECIAL TR
!  TOTING IT), SO HE'S FORCED TO USE THE PLOVER-PASSAGE TO GET IT OUT.
!  DROPPED IT, GO BACK AND PRETEND HE WASN'T CARRYING IT AFTER ALL.

210 IF (enclsd(emrald)) CALL remove (emrald)
CALL drop (emrald, loc)
GO TO 40

!  TRAVEL 303.  TROLL BRIDGE.  MUST BE DONE ONLY AS SPECIAL MOTION SO TH
!  DWARVES WON'T WANDER ACROSS AND ENCOUNTER THE BEAR.  (THEY WON'T FOLL
!  PLAYER THERE BECAUSE THAT REGION IS FORBIDDEN TO THE PIRATE.)  IF
!  PROP(TROLL)=1, HE'S CROSSED SINCE PAYING, SO STEP OUT AND BLOCK HIM.
!  (STANDARD TRAVEL ENTRIES CHECK FOR PROP(TROLL)=0.)  SPECIAL STUFF FOR

220 IF (prop(troll) /= 1) GO TO 230
CALL pspeak (troll, 1)
prop(troll)=0
CALL move (troll2, 0)
CALL move (troll2+maxobj, 0)
CALL move (troll, plac(troll))
CALL move (troll+maxobj, fixd(troll))
CALL juggle (chasm)
newloc=loc
RETURN

230 newloc=plac(troll)+fixd(troll)-loc
IF (prop(troll) == 0) prop(troll)=1
IF (.NOT.holdng(bear)) RETURN
CALL rspeak (162)
prop(chasm)=1
prop(troll)=2
CALL drop (bear, newloc)
fixed(bear)=-1
prop(bear)=3
IF (prop(spices) < 0) tally2=tally2+1
oldlc2=newloc
killed=.true.
RETURN

!  TRAVEL 304.  GROWING OR SHRINKING IN AREA OF TINY DOOR.  EACH TIME
!  HE DOES THIS, EVERYTHING MUST BE MOVED TO THE NEW LOC.  PRESUMABLY,
!  ALL HIS POSSESIONS ARE SHRUNK OR STRECHED ALONG WITH HIM.
!  PHUCE(2,4) IS AN ARRAY CONTAINING FOUR PAIRS OF "HERE" (K) AND
!  "THERE" (KK) LOCATIONS.

240 k=phuce(1,loc-161+1)
newloc=phuce(2,loc-161+1)
DO  obj=1,maxobj
  IF (obj == boat) CYCLE
  IF (place(obj) == k .AND. (fixed(obj) == 0 .OR. fixed(obj) == -1))  &
      CALL move (obj, newloc)
END DO
RETURN

!  TRAVEL #5.  PHONE BOOTH IN ROTUNDA.
!  TRYING TO SHOVE PAST GNOME, TO GET INTO PHONE BOOTH.

260 IF ((prop(booth) == 0 .AND. pct(55)) .OR. abb(loc) == 1) GO TO 270
newloc=189
IF (prop(booth) /= 1) RETURN
CALL rspeak (253)
GO TO 280

270 CALL rspeak (263)
prop(booth)=1
CALL move (gnome, 188)
280 newloc=loc
RETURN

!  TRAVEL #6.  COLLAPSING CLAY BRIDGE.  HE CAN CROSS WITH THREE (OR FEWE
!  THINGS.  IF MORE, OR IF CARRYING OBVIOUSLY HEAVY THINGS, HE MAY END U
!  IN THE DRINK.

290 newloc=235
IF (loc == 235) newloc=190
bcross=bcross+1
kk=burden(0)
IF (kk <= 4) RETURN
k=MAX(((kk+bcross)**2)/10,10)
IF (pct(k)) GO TO 300
CALL rspeak (318)
RETURN

300 CALL rspeak (319)
newloc=236
IF (holdng(lamp)) CALL move (lamp, 236)
IF (toting(axe) .AND. enclsd(axe)) CALL remove (axe)
IF (holdng(axe)) CALL move (axe, 208)
DO  obj=1,maxobj
  IF (toting(obj)) CALL dstroy (obj)
END DO
prop(chasm2)=1
RETURN

!      THE KALEIDOSCOPE CODE IS HERE
320 IF (kalflg /= 5) GO TO 330
newloc=248
oldloc=247
RETURN

330 newloc=242+ranz(5)
oldloc=newloc-1
CALL rspeak (406)
kalflg=-10
IF (newloc == 242) kalflg=0
RETURN

!  END OF SPECIALS.

END SUBROUTINE travl



!***  TREASR .TRUE. IF OBJ IS VALUABLE FOR POINTS

FUNCTION treasr(obj) RESULT(lval)

!  TREASR(OBJ)  = TRUE IF OBJECT IS A TREASURE

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),14)

RETURN
END FUNCTION treasr




!***   VAL

FUNCTION val(word) RESULT(ival)

!  RETURNS THE 'VALUE' OF A WORD, MODULO 1000.

INTEGER, INTENT(IN)  :: word
INTEGER              :: ival

! IMPLICIT INTEGER (a-z)
ival=MOD(word,1000)
RETURN
END FUNCTION val




!***  VESSEL .TRUE. IF OBJ CAN HOLD A LIQUID

FUNCTION vessel(obj) RESULT(lval)

!  VESSEL(OBJ)  = TRUE IF OBJECT IS A CONTAINER

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),15)

RETURN
END FUNCTION vessel




!***   VOCABX

FUNCTION vocabx(id,init) RESULT(ival)

!  LOOK UP ID IN THE VOCABULARY (ATAB) AND RETURN ITS "DEFINITION" (KTAB
!  -1 IF NOT FOUND.  IF INIT IS POSITIVE, THIS IS AN INITIALISATION CALL
!  UP A KEYWORD VARIABLE, AND NOT FINDING IT CONSTITUTES A BUG.  IT ALSO
!  THAT ONLY KTAB VALUES WHICH TAKEN OVER 1000 EQUAL INIT MAY BE CONSIDE
!  (THUS "STEPS", WHICH IS A MOTION VERB AS WELL AS AN OBJECT, MAY BE LO
!  AS AN OBJECT.)  AND IT ALSO MEANS THE KTAB VALUE IS TAKEN MOD 1000.

CHARACTER (LEN=6), INTENT(IN)  :: id
INTEGER, INTENT(IN)            :: init
INTEGER                        :: ival

! IMPLICIT INTEGER (a-z)
! CHARACTER (LEN=6) :: atab,dtk,vtxt,otxt,iotxt,txt
! COMMON /voccom/ ktab(600),tabsiz
! COMMON /sv3com/ dtk(9),atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)

INTEGER  :: wdclas

!       HASH=ID.XOR.'PHROG'             (DONE BY CALLER)
wdclas=init
IF (init < 0) wdclas=-init-1
DO  i=1,tabsiz
  IF (ktab(i) == -1) GO TO 20
  IF (atab(i) /= id) CYCLE
  IF (class(ktab(i)) >= wdclas) GO TO 40
END DO
CALL bug (21)

20 ival=-1
IF (init < 0) RETURN
WRITE (*,30) id
30 FORMAT (' VOCAB ERROR: CAN''T FIND WORD ''',a5,''' IN TABLE.')
CALL bug (5)

40 ival=ktab(i)
IF (init >= 0) ival=MOD(ival,1000)
RETURN
END FUNCTION vocabx



!***  WEARNG .TRUE. IF WEARING OBJ

FUNCTION wearng(obj) RESULT(lval)

!  WEARNG(OBJ)  = TRUE IF THE OBJ IS BEING WORN

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! COMMON /bitcom/ openbt,unlkbt,burnbt,wearbt
! COMMON /concom/ loccon(250),objcon(150)
! LOGICAL  :: bitset

lval=bitset(objcon(obj),wearbt)
RETURN
END FUNCTION wearng




!***   WORN  .TRUE. IF OBJ IS BEING WORN

FUNCTION worn(obj) RESULT(lval)

INTEGER, INTENT(IN)  :: obj
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
! LOGICAL :: bitset
! COMMON /concom/ loccon(250),objcon(150)

lval=bitset(objcon(obj),11)
RETURN
END FUNCTION worn




!***   YES

FUNCTION yes(x,y,z) RESULT(lval)

!  PRINT MESSAGE X, WAIT FOR YES/NO ANSWER.

INTEGER, INTENT(IN)  :: x
INTEGER, INTENT(IN)  :: y
INTEGER, INTENT(IN)  :: z
LOGICAL              :: lval

! IMPLICIT INTEGER (a-z)
CHARACTER (LEN=6)  :: reply

10 IF (x /= 0) CALL rspeak (x)
WRITE (*, 20, ADVANCE='NO')
!338    FORMAT(/,' >')
!??????????????? THE FOLLOWING WORKS BETTER ON VAXES, ETC:
20 FORMAT (/' >')
READ (*,30) reply
30 FORMAT (a6)
DO  i=1,6
  IF (reply(i:i) >= 'a' .AND. reply(i:i) <= 'z') reply(i:i)=  &
      CHAR(ICHAR(reply(i:i))-32)
END DO
IF (reply == 'YES   ' .OR. reply == 'Y     ') GO TO 50
IF (reply == 'NO    ' .OR. reply == 'N     ') GO TO 60

WRITE (*,*) ' Please answer the question.'
GO TO 10

50 lval=.true.
IF (y /= 0) CALL rspeak (y)
RETURN

60 lval=.false.
IF (z /= 0) CALL rspeak (z)
RETURN
END FUNCTION yes

END MODULE Adventure_Subs

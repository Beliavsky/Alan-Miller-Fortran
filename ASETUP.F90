PROGRAM setup
!  THIS IS THE FIRST LINE OF THE SETUP OF ADVENTURE.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-18  Time: 13:11:22
 
! THE COMMENTS ARE IN THE MAIN PROGRAM

USE Adventure_Common
USE Adventure_Subs
IMPLICIT NONE

! IMPLICIT INTEGER (a-z)
CHARACTER (LEN= 6) :: djj,dk,dkk
CHARACTER (LEN=72) :: textline
CHARACTER (LEN=16) :: dl
CHARACTER (LEN= 2) :: sp1, sp2

! LOGICAL  :: ajar,at,athand,bitset,blind,closed,closng,dark,dead,  &
!             edible,enclsd,forced,gaveup,here,hinged,holdng,inside,lmwarn,  &
!             locks,outsid,opaque,panic,pct,plural,portal,printd,scorng,small, &
!             toting,treasr,locked,vessel,wearng,worn,yea,yesm

!       DATA LINSIZ/25000/,TRVSIZ/1600/,TABSIZ/600/,LOCSIZ/250/,
!     1  VRBSIZ/60/,RTXSIZ/450/,CLSMAX/12/,HNTSIZ/20/,
!     2  MAXOBJ/150/,MAXLOC/300/,HNTMIN/7/,PTBSIZ/300/,ADJSIZ/50/,
!     3  VKYSIZ/60/,BLKLIN/.TRUE./,DWFMAX/6/,ISWIZ/.FALSE./

! DATA linuse/0/
INTEGER, SAVE  :: linuse = 0

!  PHUCE CONSISTS OF FOUR PAIRS OF ORIGIN/DESTINATION LOCATIONS FROM/TO
!  WHICH ONE IS TRANSPORTED ON UTTERING THE ELFIN CURSE AT THE TINY DOOR.
!  HE CAN GO FROM BIG TO SMALL OR SMALL TO BIG, ON EITHER SIDE OF THE DOOR.

! DATA phuce2/158,160,160,158,167,166,166,167/
INTEGER, SAVE  :: phuce2(2,4) = RESHAPE(  &
                  (/ 158,160,160,158,167,166,166,167 /), (/ 2, 4 /) )

!       DATA DEADBT,OPENBT,LOCKBT,BURNBT,WEARBT /10,2,4,6,12/
!  STATEMENT FUNCTIONS
!
!
!  AJAR(OBJ)    = TRUE IF THE OBJECT IS OPEN
!  AT(OBJ)      = TRUE IF ON EITHER SIDE OF TWO-PLACED OBJECT
!  ATHAND(OBJ)  = TRUE IF OBJECT IS HERE AND NOT IN CLOSED CONTAINER.
!  BITSET(COND,L,N) = TRUE IF COND(L) HAS BIT N SET (BIT 0 IS UNITS BIT)
!  BLIND(DUMMY) = TRUE IF HERO CAN'T SEE (TOO DARK OR GLAREY)
!  DARK(DUMMY)  = TRUE IF LOCATION "LOC" IS DARK
!  DEAD(OBJ)    = TRUE IF CRITTER IS KILLED (OR IN ENCHANTED SLEEP)
!  FORCED(LOC)  = TRUE IF LOC MOVES WITHOUT ASKING FOR INPUT (COND=2)
!  HERE(OBJ)    = TRUE IF THE OBJ IS AT "LOC" (OR IS BEING CARRIED)
!  HINGED(OBJ)  = TRUE IF OBJECT CAN BE OPENED/SHUT.
!  INSIDE(LOC)  = TRUE IF LOCATION IS WELL WITHIN THE CAVE
!  LIQ(DUMMY)   = OBJECT NUMBER OF LIQUID IN BOTTLE
!  LIQLOC(LOC)  = OBJECT NUMBER OF LIQUID (IF ANY) AT LOC
!  LIVING(OBJ)  = TRUE IF OBJ IS SOME SORT OF CRITTER
!  LOCKED(OBJ)  = TRUE IF OBJECT IS LOCKED. (NEED NOT HAVE A LOCK,
!                  E.G., RUSTY DOOR)
!  LOCKS(OBJ)   = TRUE IF OBJECT HAS A LOCK.
!  OPAQUE(OBJ)  = TRUE IF CONTAINER IS NOT TRANSPARENT (SACK, CHEST)
!                 TRANSPARENT OBJS: BOTTLE(GLASS), CAGE(WICKER)
!  OUTSID(LOC)  = TRUE IF LOCATION IS OUTSIDE THE CAVE
!  PCT(N)       = TRUE N% OF THE TIME (N INTEGER FROM 0 TO 100)
!  PLURAL(OBJ)  = TRUE IF IT IS A PLURAL OBJ (SHOES, KEYS, ETC.)
!  PORTAL(LOC)  = TRUE IS LOCATION IS IN CAVE "ENTRANCE"
!  PRINTD(OBJ)  = TRUE IF OBJECT CAN BE READ.
!  SMALL(OBJ)   = TRUE IF OBJ FITS INTO SACK
!  TOTING(OBJ)  = TRUE IF THE OBJ IS BEING CARRIED
!  TREASR(OBJ)  = TRUE IF OBJECT IS A TREASURE
!  VESSEL(OBJ)  = TRUE IF OBJECT IS A CONTAINER
!  WEARNG(OBJ)  = TRUE IF OBJECT IS BEING WORN
!  WORN(OBJ)    = TRUE IF THE OBJECT CAN BE WORN
!
!  CLOSED SAYS WHETHER WE'RE ALL THE WAY CLOSED
!  CLOSNG SAYS WHETHER ITS CLOSING TIME YET
!  DEMO IS TRUE IF THIS IS A PRIME-TIME DEMONSTRATION GAME
!  GAVEUP SAYS WHETHER HE EXITED VIA "QUIT"
!  LMWARN SAYS WHETHER HE'S BEEN WARNED ABOUT LAMP GOING DIM
!  PANIC SAYS WHETHER HE'S FOUND OUT HE'S TRAPPED IN THE CAVE
!  SCORNG INDICATES TO THE RATING ROUTINE WHETHER WE'RE DOING A "SCORE"
!  WZDARK SAYS WHETHER THE LOC HE'S LEAVING WAS DARK
!  YEA IS RANDOM YES/NO REPLY

INTEGER  :: allspace, j12, jl1, jl2, jnew, zclyd


!  CLEAR OUT THE VARIOUS TEXT-POINTER ARRAYS.  ALL TEXT IS STORED IN ARR
!  LINES; EACH LINE IS PRECEDED BY A WORD POINTING TO THE NEXT POINTER (
!  THE WORD FOLLOWING THE END OF THE LINE).  THE POINTER IS NEGATIVE IF
!  FIRST LINE OF A MESSAGE.  THE TEXT-POINTER ARRAYS CONTAIN INDICES OF
!  POINTER-WORDS IN LINES.  STEXT(N) IS SHORT DESCRIPTION OF LOCATION N.
!  LTEXT(N) IS LONG DESCRIPTION.  PTEXT(N) POINTS TO MESSAGE FOR PROP(
!  SUCCESSIVE PROP MESSAGES ARE FOUND BY CHASING POINTERS.  RTEXT CONTAI
!  SECTION 6'S STUFF.  CTEXT(N) POINTS TO A PLAYER-CLASS MESSAGE.  MTEXT
!  SECTION 12.  WE ALSO CLEAR COND.  SEE DESCRIPTION OF SECTION 9 FOR DE

deadbt=10
openbt=2
lockbt=4
burnbt=6
wearbt=12
phuce(1:2,1:4) = phuce2(1:2,1:4)

!       DATA LINSIZ/25000/,TRVSIZ/1600/,TABSIZ/600/,LOCSIZ/250/,
!     1  VRBSIZ/60/,RTXSIZ/450/,CLSMAX/12/,HNTSIZ/20/,
!     2  MAXOBJ/150/,MAXLOC/300/,HNTMIN/7/,PTBSIZ/300/,ADJSIZ/50/,
!     3  VKYSIZ/60/,BLKLIN/.TRUE./,DWFMAX/6/,ISWIZ/.FALSE./
linsiz=25000
trvsiz=1600
tabsiz=600
locsiz=250
vrbsiz=60
rtxsiz=450
clsmax=12
hntsiz=20
maxobj=150
maxloc=300
hntmin=7
ptbsiz=300
adjsiz=50
vkysiz=60
blklin=.true.
dwfmax=6
iswiz=.false.

zclyd=ICHAR('c') + 256*(ICHAR('L') + 256*(ICHAR('y') + 256*ICHAR('D')))
allspace=32 + 256*(32 + 256*(32+256*32))
DO  i=1,maxloc
  IF (i <= maxobj) ptext(i)=0
  IF (i <= rtxsiz) rtext(i)=0
  IF (i <= clsmax) ctext(i)=0
  IF (i > locsiz) CYCLE
  stext(i)=0
  ltext(i)=0
  loccon(i)=0
END DO
DO  i=1,maxobj
  points(i)=0
  objcon(i)=0
END DO
OPEN (UNIT=2, FILE='ADVDAT', STATUS='OLD')
key(1:250)=0
adjkey(1:50)=0
DO  i=1,150
  weight(i)=0
  plac(i)=0
  fixd(i)=0
  adjtab(i)=0
END DO
linuse=1
trvs=1
ptbs=1
atbs=1
clsses=1


!  START NEW DATA SECTION.  SECT IS THE SECTION NUMBER.

80 READ (2, *) sect
WRITE (*,100) sect
100 FORMAT (' READING TABLE ', i2, '...')
oldloc=-1
SELECT CASE ( sect+1 )
  CASE (    1)
    GO TO 660
  CASE (2, 3, 6, 7, 11)
    GO TO 120
  CASE (    4)
    GO TO 330
  CASE (    5)
    GO TO 260
  CASE (    8)
    GO TO 360
  CASE (    9)
    GO TO 380
  CASE (   10)
    GO TO 400
  CASE (   12)
    GO TO 430
  CASE (13, 14)
    GO TO 110
  CASE (   15)
    GO TO 460
  CASE (   16)
    GO TO 550
  CASE (   17)
    GO TO 580
END SELECT

110 CALL bug (9)

!  SECTIONS 1, 2, 5, 6, 10  READ MESSAGES AND SET UP POINTERS.

120 jl1=linuse+1
jl2=linuse+18

READ (2,130) loc,(lines(j12),j12=jl1,jl2),kk
130 FORMAT (i8, 19A4)
IF (kk == allspace) GO TO 150
WRITE(*, 140) loc
140 FORMAT (' LINE FOR LOCN ', i4, ' TOO LONG.')
CALL bug (0)

150 IF (loc == -1) GO TO 80
DO  k=1,18
  jkk=linuse+19-k
  IF (lines(jkk) /= allspace) GO TO 170
END DO
IF (loc == 0) GO TO 120
!  ABOVE KLUGE IS TO AVOID F40 BUG IF CRLF BROKEN ACROSS RECORD BOUNDARY
CALL bug (1)

170 jl1=linuse+1
DO  jk1=jl1,jkk
  lines(jk1)=IEOR(lines(jk1),zclyd)
END DO
lines(linuse)=jkk+1
IF (loc == oldloc) GO TO 250
IF (loc > oldloc .OR. sect == 5) GO TO 200
WRITE(*, 190) loc,sect
190 FORMAT (/' LINE ',i3,' OUT OF ORDER IN SECTION ',i2)
CALL bug (10)

200 lines(linuse)=-lines(linuse)
IF (sect == 10) GO TO 240
IF (sect == 6) GO TO 230
IF (sect == 5) GO TO 220
IF (sect == 1) GO TO 210

stext(loc)=linuse
GO TO 250

210 ltext(loc)=linuse
GO TO 250

220 IF (loc > 0 .AND. loc <= maxobj) ptext(loc)=linuse
GO TO 250

230 IF (loc > rtxsiz) CALL bug (6)
rtext(loc)=linuse
GO TO 250

240 ctext(clsses)=linuse
cval(clsses)=loc
clsses=clsses+1
GO TO 250

250 linuse=jkk+1
lines(linuse)=-1
oldloc=loc
IF (linuse+18 > linsiz) CALL bug (2)
GO TO 120

!  THE STUFF FOR SECTION 4 IS ENCODED HERE.  EACH "FROM-LOCATION" GETS A
!  CONTIGUOUS SECTION OF THE "TRAVEL" ARRAY.  EACH ENTRY IN TRAVEL IS
!  NEWLOC*1000 + KEYWORD (FROM SECTION 3, MOTION VERBS), AND IS NEGATED
!  THIS IS THE LAST ENTRY FOR THIS LOCATION.  KEY(N) IS THE INDEX IN TRA
!  OF THE FIRST OPTION AT LOCATION N.

! 260 READ (2,270) loc,jnew, dtk(1:8)
! 270 FORMAT (2I8,8(a6, '  '))
260 READ (2, '(2i8, a)') loc, jnew, textline
IF (loc == 0) GO TO 260
!  ABOVE KLUGE IS TO AVOID AFOREMENTIONED F40 BUG
IF (loc == -1) GO TO 80

!  Strip elements of dtk from textline
l = (LEN_TRIM(textline) + 7)/8
DO i = 1, l
  READ(textline(8*i-7:8*i), '(a)') dtk(i)
END DO
dtk(l+1:8) = '        '

IF (loc >= oldloc) GO TO 280
WRITE(*, 190) loc, sect
CALL bug (10)

280 IF (key(loc) /= 0) GO TO 290
key(loc)=trvs
GO TO 300

290 travel(trvs-1)=-travel(trvs-1)
300 DO  l=1,8
  IF (dtk(l) == '      ') EXIT
  k=vocabx(dtk(l),-1)
  travel(trvs)=jnew*1000 + k
  trvs=trvs+1
  IF (trvs == trvsiz) CALL bug (3)
END DO
travel(trvs-1)=-travel(trvs-1)
GO TO 260

!  HERE WE READ IN THE VOCABULARY.  KTAB(N) IS THE WORD NUMBER, ATAB(N)
!  THE CORRESPONDING WORD.  THE -1 AT THE END OF SECTION 3 IS LEFT IN KT
!  AS AN END-MARKER.  THE WORDS ARE GIVEN A MINIMAL HASH TO MAKE READING
!  CORE-IMAGE HARDER.

330 DO  tabndx=1,tabsiz
  READ (2,340) ktab(tabndx),atab(tabndx)
  340 FORMAT (i8,a6)
!       IF(KTAB(TABNDX).EQ.0)GOTO 1043
!  ABOVE KLUGE IS TO AVOID AFOREMENTIONED F40 BUG
  IF (ktab(tabndx) == -1) GO TO 80
END DO
CALL bug (4)

!  READ IN THE INITIAL LOCATIONS FOR EACH OBJECT.  ALSO THE IMMOVABILITY
!  PLAC CONTAINS INITIAL LOCATIONS OF OBJECTS.  FIXD IS -1 FOR IMMOVABLE
!  OBJECTS (INCLUDING THE SNAKE), OR = SECOND LOC FOR TWO-PLACED OBJECTS
!  WEIGHT CONTAINS THE HEAVINESS OF EACH OBJ, ON A SCALE OF 1-10.

360 READ (2,370) obj,j,k,kk,dl,tk(1:3)
370 FORMAT (4I8,a16,3I8)
IF (obj == -1) GO TO 80
IF (obj <= 0 .OR. obj > maxobj) CALL bug (17)
plac(obj)=j
fixd(obj)=k
weight(obj)=kk
!  READ DEFAULT OBJECT NAMES.  FOR USE IN 'TAKE ALL' COMMANDS.
!  READ POINT VALUES FOR TREASURES.
k=1
IF (tk(3) < 0) k=-1
points(obj)=(tk(1)*1000000*k) + (tk(2)*1000*k) + tk(3)
GO TO 360

!  READ DEFAULT MESSAGE NUMBERS FOR ACTION VERBS, STORE IN ACTSPK.

380 READ (2,390) verb,j
390 FORMAT (10I8)
IF (verb == -1) GO TO 80
actspk(verb)=j
GO TO 380

!  READ INFO ABOUT AVAILABLE LIQUIDS AND OTHER CONDITIONS, STORE IN COND

400 READ (2,390) k, tk(1:9)
IF (k == -1) GO TO 80
DO  i=1,9
  loc=tk(i)
  IF (loc == 0) GO TO 400
  IF (k > 7) GO TO 410
  IF (bitset(loccon(loc),k)) CALL bug (8)
  loccon(loc)=loccon(loc) + bits(k)
  CYCLE
  410 loccon(loc)=ior(loccon(loc),256*(k-7))
END DO
tk(1:9) = 0
GO TO 400

!  READ DATA FOR HINTS.

430 hntmax=0
440 READ (2,390) k, tk(1:9)
IF (k == -1) GO TO 80
IF (k == 0) GO TO 440
IF (k < 0 .OR. k  > hntsiz) CALL bug (7)
DO  i=1,4
  hints(k,i)=tk(i)
END DO
hntmax=MAX(hntmax,k)
GO TO 440

!  SECTION 14 IS THE PREPOSITION TABLE.
460 READ (2,470) dk,sp1,dkk,sp2,textline
470 FORMAT (2(a6, a2), a)
!       IF(DK.EQ.0)GOTO 1100
!  ABOVE KLUGE IS TO AVOID AFOREMENTIONED F40 BUG
IF (dk == '-1    ') GO TO 80

l = (LEN_TRIM(textline) + 7)/8
DO i = 1, l
  READ(textline(8*i-7:8*i), '(a)') dtk(i)
END DO
DO i = l+1, 8
  dtk(i) = '      '
END DO

verb=val(vocabx(dk,-4))
djj=dk
IF (verb == -1) GO TO 480
prep=val(vocabx(dkk,-6))
djj=dkk
IF (prep == -1) GO TO 480

IF (vkey(verb) /= 0) GO TO 490
vkey(verb)=ptbs
GO TO 500

480 WRITE(*, 510) djj
CALL bug (14)
GO TO 460


490 ptab(ptbs-1)=-ptab(ptbs-1)
500 DO  l=1,8
  IF (dtk(l) == '      ') EXIT
  k=val(vocabx(dtk(l),-3))
  IF (k /= -1) GO TO 520
  k=999
  IF (l == 1 .AND. dtk(1) == 'ANY   ') GO TO 520
  djj=dtk(l)
  WRITE(*, 510) djj
  510 FORMAT (/' UNRECOGNIZED WORD "', a6, '" IN PREP/OBJ TABLE.')
  CALL bug (14)
  
  520 ptab(ptbs)=prep*1000 + k
  ptbs=ptbs+1
  IF (ptbs == ptbsiz) CALL bug (15)
END DO
ptab(ptbs-1)=-ptab(ptbs-1)
GO TO 460

!  READ CONDITION BITS FOR OBJECTS.  KK IS THE BIT; TK(I), THE OBJ LIST.

550 READ (2,390) ikk, tk(1:9)
IF (ikk == -1) GO TO 80
DO  i=1,9
  obj=tk(i)
  IF (obj == 0) GO TO 550
  IF (obj <= 0 .OR. obj > maxobj) CALL bug (17)
  IF (.NOT.bitset(objcon(obj),ikk)) GO TO 570
  WRITE (6,560) obj,ikk
  560 FORMAT ('BIT SET TWICE OBJ=', i5, '  BIT= ', i5)
  CALL bug (16)
  570 objcon(obj)=ior(objcon(obj),bits(ikk))
END DO
GO TO 550

!  SECTION 17 IS THE ADJECTIVE TABLE.
580 READ (2, '(a)') textline
READ(textline(1:8), '(a)') dk
!       IF(DK.EQ.0)GOTO 1140
!  ABOVE KLUGE IS TO AVOID AFOREMENTIONED F40 BUG
IF (dk == '-1    ') GO TO 80

l = (LEN_TRIM(textline) + 7)/8
DO i = 1, l-1
  READ(textline(8*i+1:8*i+8), '(a)') dtk(i)
END DO
dtk(l:9) = '      '

adj=vocabx(dk,-7)
djj=dk
IF (adj == -1 .OR. class(adj) /= 6) GO TO 640
adj=val(adj)
IF (adjkey(adj) /= 0) GO TO 600
adjkey(adj)=atbs
GO TO 610

600 adjtab(atbs-1)=-adjtab(atbs-1)
610 DO  l=1,9
  IF (dtk(l) == '      ') EXIT
  djj=dtk(l)
  k=vocabx(dtk(l),-3)
  IF (k == -1 .OR. class(k) /= 2) GO TO 640
  adjtab(atbs)=val(k)
  atbs=atbs+1
  IF (atbs == maxobj) CALL bug (18)
END DO
adjtab(atbs-1)=-adjtab(atbs-1)
GO TO 580

640 WRITE(*, 650) djj
650 FORMAT (/' UNRECOGNIZED WORD "',a6,'" IN ADJECTIVE TABLE.')
CALL bug (13)

!  EVERYTHING IS READ!  NOW FINISH CONSTRUCTING INTERNAL DATA FORMAT.
!  DEFINE SOME HANDY MNEMONICS.  THESE CORRESPOND TO OBJECT NUMBERS.
!  (INCLUDES TREASURES.)

660 CLOSE (UNIT=2)
anvil=vocabx('ANVIL ',2)
axe=vocabx('AXE   ',2)
batter=vocabx('BATTER',2)
bear=vocabx('BEAR  ',2)
bees=vocabx('BEES  ',2)
billbd=vocabx('BILLBO',2)
bird=vocabx('BIRD  ',2)
boat=vocabx('BOAT  ',2)
book=vocabx('BOOK  ',2)
book2=book+1
booth=vocabx('BOOTH ',2)
bottle=vocabx('BOTTLE',2)
brush=vocabx('BRUSH ',2)
cage=vocabx('CAGE  ',2)
cakes=vocabx('CAKES ',2)
carvng=vocabx('CARVIN',2)
cask=vocabx('CASK  ',2)
chain=vocabx('CHAIN ',2)
chasm=vocabx('CHASM ',2)
chasm2=chasm+1
chest=vocabx('CHEST ',2)
clam=vocabx('CLAM  ',2)
cloak=vocabx('CLOAK ',2)
coins=vocabx('COINS ',2)
crown=vocabx('CROWN ',2)
dog=vocabx('DOG   ',2)
door=vocabx('DOOR  ',2)
dragon=vocabx('DRAGON',2)
dwarf=vocabx('DWARF ',2)
eggs=vocabx('EGGS  ',2)
emrald=vocabx('EMERAL',2)
fissur=vocabx('FISSUR',2)
flower=vocabx('FLOWER',2)
food=vocabx('FOOD  ',2)
gnome=vocabx('GNOME ',2)
grail=vocabx('GRAIL ',2)
grate=vocabx('GRATE ',2)
hive=vocabx('HIVE  ',2)
honey=vocabx('HONEY ',2)
horn=vocabx('HORN  ',2)
jewels=vocabx('JEWELS',2)
keys=vocabx('KEYS  ',2)
knife=vocabx('KNIFE ',2)
lamp=vocabx('LAMP  ',2)
lyre=vocabx('LYRE  ',2)
magzin=vocabx('MAGAZI',2)
messag=vocabx('MESSAG',2)
mirror=vocabx('MIRROR',2)
mushrm=vocabx('MUSHRO',2)
nugget=vocabx('NUGGET',2)
oil=vocabx('OIL   ',2)
!       OIL2=OIL+1
oyster=vocabx('OYSTER',2)
pearl=vocabx('PEARL ',2)
phone=vocabx('PHONE ',2)
pillow=vocabx('PILLOW',2)
plant=vocabx('PLANT ',2)
plant2=plant+1
pole=vocabx('POLE  ',2)
poster=vocabx('POSTER',2)
pyram=vocabx('PYRAMI',2)
radium=vocabx('RADIUM',2)
ring=vocabx('RING  ',2)
rocks=vocabx('ROCKS ',2)
rod=vocabx('ROD   ',2)
rod2=rod+1
rug=vocabx('RUG   ',2)
safe=vocabx('SAFE  ',2)
sapphi=vocabx('SAPPHI',2)
shield=vocabx('TUBE  ',2)
shoes=vocabx('SHOES ',2)
skey=vocabx('KEY   ',2)
slugs=vocabx('SLUGS ',2)
snake=vocabx('SNAKE ',2)
spices=vocabx('SPICES',2)
sphere=vocabx('SPHERE',2)
steps=vocabx('STEPS ',2)
sticks=vocabx('STICKS',2)
sword=vocabx('SWORD ',2)
tablet=vocabx('TABLET',2)
tdoor=door+1
tdoor2=tdoor+1
pdoor=tdoor2+1
tridnt=vocabx('TRIDEN',2)
troll=vocabx('TROLL ',2)
troll2=troll+1
vase=vocabx('VASE  ',2)
vend=vocabx('MACHIN',2)
wall=vocabx('WALL  ',2)
wall2=wall+1
water=vocabx('WATER ',2)
!       WATER2=WATER+1
wine=vocabx('WINE  ',2)
!       WINE2=WINE+1
wumpus=vocabx('WUMPUS',2)

!  THESE ARE MOTION-VERB NUMBERS.

back=vocabx('BACK  ',1)
cave=vocabx('CAVE  ',1)
dprssn=vocabx('DEPRES',1)
entrnc=vocabx('ENTRAN',1)
EXIT=vocabx('EXIT  ',1)
null=vocabx('NULL  ',1)

!  AND SOME ACTION VERBS.

find=vocabx('FIND  ',3)
GO=vocabx('GO    ',3)
hit=vocabx('HIT   ',3)
look=vocabx('LOOK  ',3)
yell=vocabx('CALL  ',3)
invent=vocabx('INVENT',3)
leave=vocabx('LEAVE ',3)
lock=vocabx('LOCK  ',3)
say=vocabx('SAY   ',3)
shut=vocabx('CLOSE ',3)
take=vocabx('TAKE  ',3)
throw=vocabx('THROW ',3)
unlock=vocabx('UNLOCK',3)
wear=vocabx('WEAR  ',3)
yank=vocabx('YANK  ',3)

!  AND A FEW PREPOSITIONS.  PREFIX 'PREP' TO DISTINGUISH THEM FROM FUNCT

prepat=vocabx('AT    ',5)
prepdn=vocabx('DOWN  ',5)
prepfr=vocabx('FROM  ',5)
prepin=vocabx('IN    ',5)
prepof=vocabx('OFF   ',5)
prepon=vocabx('ON    ',5)

!  A POPULAR LOCATION IS:

y2=33

!  HAVING READ IN THE DATABASE, CERTAIN THINGS ARE NOW CONSTRUCTED.  PRO
!  SET TO ZERO.  WE FINISH SETTING UP COND BY CHECKING FOR FORCED-MOTION
!  ENTRIES.  THE PLAC AND FIXD ARRAYS ARE USED TO SET UP ATLOC(N) AS THE
!  OBJECT AT LOCATION N, AND LINK(OBJ) AS THE NEXT OBJECT AT THE SAME LO
!  AS OBJ.  (OBJ>MAXOBJ INDICATES THAT FIXED(OBJ-MAXOBJ)=LOC; LINK(OBJ)
!  THE CORRECT LINK TO USE.)  ABB IS ZEROED; IT CONTROLS WHETHER THE ABB
!  DESCRIPTION IS PRINTED.  COUNTS MOD 5 UNLESS "LOOK" IS USED.

loc=1
DO  i=1,maxobj
  place(i)=0
  prop(i)=0
  holder(i)=0
  hlink(i)=0
  link(i)=0
  link(i+maxobj)=0
END DO

DO  i=1,locsiz
  abb(i)=0
  IF (ltext(i) == 0 .OR. key(i) == 0) GO TO 680
  k=key(i)
  IF (MOD(ABS(travel(k)),0001000) == 1) loccon(i)=2
  680 atloc(i)=0
END DO

!  SET UP THE ATLOC AND LINK ARRAYS AS DESCRIBED ABOVE.  WE'LL USE THE D
!  SUBROUTINE, WHICH PREFACES NEW OBJECTS ON THE LISTS.  SINCE WE WANT T
!  IN THE OTHER ORDER, WE'LL RUN THE LOOP BACKWARDS.  IF THE OBJECT IS I
!  LOCS, WE DROP IT TWICE.  THIS ALSO SETS UP "PLACE" AND "FIXED" AS COP
!  "PLAC" AND "FIXD".  ALSO, SINCE TWO-PLACED OBJECTS ARE TYPICALLY BEST
!  DESCRIBED LAST, WE'LL DROP THEM FIRST.

DO  i=1,maxobj
  k=maxobj+1-i
  IF (fixd(k) <= 0) CYCLE
  CALL drop (k+maxobj, fixd(k))
  CALL drop (k, plac(k))
END DO

DO  i=1,maxobj
  k=maxobj+1-i
  fixed(k)=fixd(k)
  IF (plac(k) /= 0 .AND. fixd(k) <= 0) CALL drop (k, plac(k))
END DO

!  MAKE SURE ALL THE RIGHT THINGS GET CLOSED AND LOCKED, ETC., BEFORE
!  WE GET STARTED.

!  TREASURES, AS NOTED EARLIER, ARE OBJECTS WITH BITSET(14) IN OBJCON.
!  THEIR PROPS ARE INITIALLY -1, AND ARE SET TO 0 THE FIRST TIME THEY AR
!  DESCRIBED.  TALLY KEEPS TRACK OF HOW MANY ARE NOT YET FOUND, SO WE KN
!  WHEN TO CLOSE THE CAVE.  TALLY2 COUNTS HOW MANY CAN NEVER BE FOUND (E
!  LOST BIRD OR BRIDGE).

tally=0
tally2=0
DO  i=1,maxobj
  IF (.NOT.treasr(i)) GO TO 710
  IF (ptext(i) /= 0) prop(i)=-1
  710 tally=tally-prop(i)
END DO


!  CLEAR THE HINT STUFF.  HINTLC(I) IS HOW LONG HE'S BEEN AT LOC WITH CO
!  I.  HINTED(I) IS TRUE IFF HINT I HAS BEEN USED.

DO  i=1,hntmax
  hinted(i)=.false.
  hintlc(i)=0
END DO

!  INITIALISE THE DWARVES.  DLOC IS LOC OF DWARVES, HARD-WIRED IN.  ODLO
!  PRIOR LOC OF EACH DWARF, INITIALLY GARBAGE.  DALTLC IS ALTERNATE INIT
!  FOR DWARF, IN CASE ONE OF THEM STARTS OUT ON TOP OF THE ADVENTURER.
!  OF THE 5 INITIAL LOCS ARE ADJACENT.)  DSEEN IS TRUE IF DWARF HAS SEEN
!  DFLAG CONTROLS THE LEVEL OF ACTIVATION OF ALL THIS:
!       0       NO DWARF STUFF YET (WAIT UNTIL REACHES HALL OF MISTS)
!       1       REACHED HALL OF MISTS, BUT HASN'T MET FIRST DWARF
!       2       MET FIRST DWARF, OTHERS START MOVING, NO KNIVES THROWN Y
!       3       A KNIFE HAS BEEN THROWN (FIRST SET ALWAYS MISSES)
!       3+      DWARVES ARE MAD (INCREASES THEIR ACCURACY)
!  SIXTH DWARF IS SPECIAL (THE PIRATE).  HE ALWAYS STARTS AT HIS CHEST'S
!  EVENTUAL LOCATION INSIDE THE MAZE.  THIS LOC IS SAVED IN CHLOC FOR RE
!  THE DEAD END IN THE OTHER MAZE HAS ITS LOC STORED IN CHLOC2.

chloc=114
chloc2=140
dseen(1:dwfmax)=.false.
dflag=0
dloc(1)=plac(snake)
dloc(2)=plac(booth)
dloc(3)=y2
dloc(4)=44
dloc(5)=plac(clam)
!       DLOC(6)=PLAC(VEND)
dloc(dwfmax)=chloc
daltlc=plac(nugget)

!  OTHER RANDOM FLAGS AND COUNTERS, AS FOLLOWS:
!       ABBNUM  HOW OFTEN WE SHOULD PRINT NON-ABBREVIATED DESCRIPTIONS
!       BCROSS  NUMBER OF TIMES COLLAPSING BRIDGE HAS BEEN TRAVERSED.
!       BONUS   USED TO DETERMINE AMOUNT OF BONUS IF HE REACHES CLOSING
!       CHASE   TELLS HOW CLOSE THE WUMPUS IS TO GOBBLING HIM UP
!       CLOCK1  NUMBER OF TURNS FROM FINDING LAST TREASURE TILL CLOSING
!       CLOCK2  NUMBER OF TURNS FROM FIRST WARNING TILL BLINDING FLASH
!       CLOCK3  NUMBER OF TURNS IN REPOSITORY TILL PHONE RINGS.
!               AFTER TICKING TO 0, TICKS 7 TIMES TO WAKE DWARVES.
!       COMBO   CURRENT PROGRESS IN GIVING SAFE'S COMBINATION
!       DETAIL  HOW OFTEN WE'VE SAID "NOT ALLOWED TO GIVE MORE DETAIL"
!       DKILL   NUMBER OF DWARVES KILLED (UNUSED IN SCORING, NEEDED FOR
!       FOOBAR  CURRENT PROGRESS IN SAYING "FEE FIE FOE FOO".
!       HEALTH  PERCENTAGE OF MAXIMUM (100) FITNESS
!       IWEST   HOW MANY TIMES HE'S SAID "WEST" INSTEAD OF "W"
!       KNFLOC  0 IF NO KNIFE HERE, LOC IF KNIFE HERE, -1 AFTER CAVEAT
!       LIMIT   LIFETIME OF LAMP (NOT SET HERE)
!       MAXDIE  NUMBER OF REINCARNATION MESSAGES AVAILABLE (UP TO 5)
!       NUMDIE  NUMBER OF TIMES KILLED SO FAR
!       TERSE   IF TRUE, NEVER PRINT LONG LOCATION DESCRIPTIONS
!       TURNS   TALLIES HOW MANY COMMANDS HE'S GIVEN (IGNORES YES/NO)
!       WASTE   TELLS HOW LONG HE HAS USED LAMP IN LIGHTED AREA.

!       LOGICALS WERE EXPLAINED EARLIER

abbnum=5
bcross=0
bonus=0
clock1=30
clock2=50
clock3=20 + ranz(20)
chase=0
closed=.false.
closng=.false.
combo=0
detail=0
dkill=0
foobar=0
gaveup=.false.
health=100
iwest=0
knfloc=0
lmwarn=.false.
DO  i=0,4
  IF (rtext(2*i+81) /= 0) maxdie=i+1
END DO
numdie=0
panic=.false.
scorng=.false.
terse=.false.
turns=0
waste=0

!  SETUP THE LIQUIDS ACCORDING TO CONTAINER PROP VALUES
liqtyp(1)=water
liqtyp(2)=0
liqtyp(3)=oil
liqtyp(4)=0
liqtyp(5)=wine

prop(pole)=1
prop(skey)=1
place(water)=-1
CALL insert (water, bottle)
place(book)=-1
CALL insert (book, safe)

!  AND CLEAR OUT ANY LEFTOVER WORD VECTORS...
CALL clrlin()
words(1:35)=0
wdx=0

!  IF SETUP=1, REPORT ON AMOUNT OF ARRAYS ACTUALLY USED, TO PERMIT REDUC


jj=0
DO  k=1,vkysiz
  IF (vkey(k) /= 0) jj=jj+1
END DO

DO  k=1,locsiz
  kk=locsiz+1-k
  IF (ltext(kk) /= 0) GO TO 780
END DO

780 ll=0
obj=0
DO  k=1,maxobj
  IF (treasr(k)) ll=ll+1
  IF (ptext(k) /= 0) obj=obj+1
END DO

DO  k=1,tabndx
  IF (ktab(k)/1000 == 2) verb=ktab(k) - 2000
END DO

DO  k=1,rtxsiz
  j=rtxsiz+1-k
  IF (rtext(j) /= 0) EXIT
END DO


CALL rating (score, 0, .false., .false., .false., .false., hntmax)
k=maxobj
WRITE(*, 830) linuse,linsiz,trvs,trvsiz,tabndx,tabsiz,kk,locsiz,obj,k,ll,   &
              verb,vrbsiz,j,rtxsiz,clsses,clsmax,hntmax,hntsiz,ptbs,ptbsiz, &
              jj,verb,mxscor
830 FORMAT (' TABLE SPACE USED:'/' ',i6,' OF ',i6,' WORDS OF MESSAGES'/  &
            ' ',i6,' OF ',i6,' TRAVEL OPTIONS'/  &
            ' ',i6,' OF ',i6,' VOCABULARY WORDS'/  &
            ' ',i6,' OF ',i6,' LOCATIONS'/  &
            ' ',i6,' OF ',i6,' OBJECTS OF WHICH ',i2,' ARE TREASURES.'/  &
            ' ',i6,' OF ',i6,' ACTION VERBS'/  &
            ' ' ,i6,' OF ',i6,' RTEXT MESSAGES'/  &
            ' ',i6,' OF ',i6,' CLASS MESSAGES' /  &
            ' ',i6,' OF ',i6,' HINTS'/  &
            ' ',i6,' OF ',i6,' VERB/PREP/OBJ OPTIONS'/  &
            ' ',i6,' OF ',i6,' VERBS TAKE PREPOSITIONS'//  &
            ' MAXIMUM SCORE FOR THIS VERSION IS ',i4,' POINTS.'/)
OPEN (UNIT=16, FILE='ADVTXT', STATUS='UNKNOWN', FORM='UNFORMATTED',  &
      ACCESS='SEQUENTIAL')
WRITE (16) iswiz,adjkey,adjtab,adjsiz,openbt,lockbt,burnbt,wearbt
WRITE (16) blklin,loccon,objcon,numdie,maxdie,turns,killed
WRITE (16) dwarf,knife,knfloc,dflag,dseen,dloc,odloc,dwfmax
WRITE (16) holder,hlink,hintlc,hinted,hints,hntsiz,hntmin
WRITE (16) bottle,cask,water,oil,wine,liqtyp
WRITE (16) loc,oldloc,oldlc2,newloc,maxloc
WRITE (16) ltext,stext,key,abb,locsiz
WRITE (16) back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,bear,  &
    boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,lamp,  &
    pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,troll2,  &
    emrald,spices,find,yell,invent,leave,pour,say,take,throw,iwest,phuce,tk

WRITE (16) plac,fixd,weight,prop,points
WRITE (16) atloc,link,place,fixed,maxobj
WRITE (16) vkey,ptab,vkysiz,ptbsiz,travel
WRITE (16) lines,rtext,ptext,wdx,ktab,tabsiz
WRITE (16) verbs,vrbx,objs,objx,iobjs,iobx,prep,words

WRITE (16) abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,  &
    clock2,clock3,closed,closng,clsmax,combo,deadbt,detail,dkill,  &
    dtotal,dwarfn,flg239,foo,foobar,food,gaveup,health,hint,hit,  &
    hntmax,i,ikk,iloc,iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,  &
    lmwarn,lock,logout,messag,obj,panic,portal,ptbs,rdflag,retn,  &
    rtxsiz,score,scorng,sect,skey,sloc,spk,start,stick,tabndx,tally,  &
    tally2,terse,trvs,trvsiz,vend,verb,vrbsiz,waste,wkday,wkend,  &
    wzdark,yea,actspk,ctext,cval,hname
WRITE (16) anvil,batter,bees,billbd,bird,brush,cage,cakes,chain,  &
    chest,chloc,chloc2,clam,cloak,clsses,coins,crown,daltlc,dog,  &
    dragon,eggs,fissur,flower,gatloc,grail,hive,honey,horn,jewels,  &
    keys,lyre,magzin,mirror,mushrm,mxscor,nugget,oyster,pearl,phone,  &
    pillow,pole,poster,prepat,prepdn,prepfr,prepin,prepof,prepon,  &
    pyram,radium,ring,rug,sapphi,shield,shoes,shut,slugs,snake,sphere,  &
    steps,sticks,sword,tablet,tridnt,unlock,vase,wall,wall2,wear,wumpus,  &
    y2,yank
WRITE (16) dtk,atab,vtxt,otxt,iotxt,txt
CLOSE (16)

STOP
END PROGRAM setup

PROGRAM Adventure

!  THIS IS THE FIRST LINE OF ADVENTURE, MAIN MODULE.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-18  Time: 13:11:29

!  CURRENT LIMITS:
!      20000 WORDS OF MESSAGE TEXT (LINES, LINSIZ).
!       1600 TRAVEL OPTIONS (TRAVEL, TRVSIZ).
!       600 VOCABULARY WORDS (KTAB, ATAB, TABSIZ).
!       250 LOCATIONS (LTEXT, STEXT, KEY, LOCCON, ABB, ATLOC,
!                       LOCSIZ, MAXLOC).
!       150 OBJECTS (PLAC, PLACE, FIXD, FIXED, LINK (TWICE), PTEXT, POINTS,
!                       HOLDER, HLINK, OBJCON, PROP, WEIGHT, MAXOBJ).
!        60 "ACTION" VERBS (ACTSPK, VRBSIZ, VKEY).
!       300 VERB/PREP/OBJ COMBINATIONS (PTAB, PTBSIZ).
!        50 ADJECTIVES (ADJKEY, ADJSIZ)
!       150 MODIFIED NOUNS (ADJTAB, MAXOBJ)
!       450 RANDOM MESSAGES (RTEXT, RTXSIZ).
!        12 DIFFERENT PLAYER CLASSIFICATIONS (CTEXT, CVAL, CLSMAX).
!        20 HINTS, LESS 3 (HINTLC, HINTED, HINTS, HNTSIZ, HNTMIN).

!  THERE ARE ALSO LIMITS WHICH CANNOT BE EXCEEDED DUE TO THE STRUCTURE OF
!  THE DATABASE.  (E.G., THE VOCABULARY USES N/1000 TO DETERMINE WORD TYPE,
!  SO THERE CAN'T BE MORE THAN 1000 WORDS OF ANY CLASS.)  THESE UPPER
!  LIMITS ARE:
!     1000 NON-SYNONYMOUS VOCABULARY WORDS OF EACH CLASS.  CURRENTLY
!     DEFINED CLASSES ARE:
!             1. MOTION/DIRECTION WORDS (EAST, WEST, UP, JUMP, ETC.)
!             2. NOUNS/OBJECTS (LAMP, KEYS, TROLL, ETC.)
!             3. ACTION VERBS, TRANSITIVE & INTRANSITIVE (TAKE, DROP, KILL, IN
!             4. MISCELLANEOUS WORDS; MOSTLY THINGS OR ACTIONS WHICH GENERATE
!                FIXED REPLIES (FEE-FIE-FOE-FOO, TREE, CAVE, CURSES)
!             5. PREPOSITIONS, DUPLICATES MANY WORDS IN SECTION 1.
!             6. ADJECTIVES
!             7. CONJUNCTIONS
!     450 LOCATIONS
!     100 OBJECTS WHICH CAN BE USED IN TRAVEL TABLE (PLUS 900 MORE, WHICH CAN


!  DESCRIPTION OF THE DATABASE FORMAT
!
!
!  THE DATA FILE CONTAINS SEVERAL SECTIONS.  EACH BEGINS WITH A LINE CONTAINING
!  A NUMBER IDENTIFYING THE SECTION, AND ENDS WITH A LINE CONTAINING "-1".
!
!  SECTION 1: LONG FORM DESCRIPTIONS.  EACH LINE CONTAINS A LOCATION NUMBER,
!       A TAB, AND A LINE OF TEXT.  THE SET OF (NECESSARILY ADJACENT) LINES
!       WHOSE NUMBERS ARE X FORM THE LONG DESCRIPTION OF LOCATION X.
!
!  SECTION 2: SHORT FORM DESCRIPTIONS.  SAME FORMAT AS LONG FORM.  NOT ALL
!       PLACES HAVE SHORT DESCRIPTIONS.
!
!  SECTION 3: VOCABULARY.  EACH LINE CONTAINS A NUMBER (N), A TAB, AND A
!       FIVE-LETTER WORD.  CALL M=N/1000.  IF M=0, THEN THE WORD IS A MOTION
!       VERB FOR USE IN TRAVELLING (SEE SECTION 4).  ELSE, IF M=1, THE WORD IS
!       AN OBJECT.  ELSE, IF M=2, THE WORD IS AN ACTION VERB (SUCH AS "CARRY"
!       OR "ATTACK").  ELSE, IF M=3, THE WORD IS A SPECIAL CASE VERB (SUCH AS
!       "DIG") AND N MOD 1000 IS AN INDEX INTO SECTION 6.  OBJECTS FROM 50 TO
!       (CURRENTLY, ANYWAY) 79 ARE CONSIDERED TREASURES (FOR PIRATE, CLOSEOUT).
!
!  SECTION 4: TRAVEL TABLE.  EACH LINE CONTAINS A LOCATION NUMBER (X), A SECOND
!       LOCATION NUMBER (Y), AND A LIST OF MOTION VERBS (SEE SECTION 3).
!       EACH MOTION REPRESENTS A VERB WHICH WILL GO TO Y IF CURRENTLY AT X.
!       Y, IN TURN, IS INTERPRETED AS FOLLOWS.  LET M=Y/1000, N=Y MOD 1000.
!               IF N<=MAXLOC    IT IS THE LOCATION TO GO TO.
!               IF MAXLOC<N<=500   N-MAXLOC IS USED IN A COMPUTED GOTO
!                                       TO A SECTION OF SPECIAL CODE.
!               IF N>500        MESSAGE N-500 FROM SECTION 6 IS PRINTED,
!                                       AND HE STAYS WHEREVER HE IS.
!       MEANWHILE, M SPECIFIES THE CONDITIONS ON THE MOTION.
!               IF M=0          IT'S UNCONDITIONAL.
!               IF 0<M<100      IT IS DONE WITH M% PROBABILITY.
!               IF M=100        UNCONDITIONAL, BUT FORBIDDEN TO DWARVES.
!               IF 100<M<=200   HE MUST BE CARRYING OBJECT M-100.
!               IF 200<M<=300   MUST BE CARRYING OR IN SAME ROOM AS M-200.
!               IF 300<M<=400   PROP(M MOD 100) MUST *NOT* BE 0.
!               IF 400<M<=500   PROP(M MOD 100) MUST *NOT* BE 1.
!               IF 500<M<=600   PROP(M MOD 100) MUST *NOT* BE 2, ETC.
!       IF THE CONDITION (IF ANY) IS NOT MET, THEN THE NEXT *DIFFERENT*
!       "DESTINATION" VALUE IS USED (UNLESS IT FAILS TO MEET *ITS* CONDITIONS,
!       IN WHICH CASE THE NEXT IS FOUND, ETC.).  TYPICALLY, THE NEXT DEST WILL
!       BE FOR ONE OF THE SAME VERBS, SO THAT ITS ONLY USE IS AS THE ALTERNATE
!       DESTINATION FOR THOSE VERBS.  FOR INSTANCE:
!               15      110022  29      31      34      35      23      43
!               15      14      29
!       THIS SAYS THAT, FROM LOC 15, ANY OF THE VERBS 29, 31, ETC., WILL TAKE
!       HIM TO 22 IF HE'S CARRYING OBJECT 10, AND OTHERWISE WILL GO TO 14.
!               11      303008  49
!               11      9       50
!       THIS SAYS THAT, FROM 11, 49 TAKES HIM TO 8 UNLESS PROP(3)=0, IN WHICH
!       CASE HE GOES TO 9.  VERB 50 TAKES HIM TO 9 REGARDLESS OF PROP(3).
!       (SEE DESCRIPTION FOR SECTION 14 FOR A SCHEMATIC OF TABLES.)
!
!  SECTION 5: OBJECT DESCRIPTIONS.  EACH LINE CONTAINS A NUMBER (N), A TAB.
!       AND A MESSAGE.  IF N IS FROM 1 TO MAXOBJ, THE MESSAGE IS THE "INVENTORY"
!       MESSAGE FOR OBJECT N.  OTHERWISE, N SHOULD BE 0000, 1000, 2000, ETC., AN
!       THE MESSAGE SHOULD BE THE DESCRIPTION OF THE PRECEDING OBJECT WHEN ITS
!       PROP VALUE IS N/1000.  THE N/1000 IS USED ONLY TO DISTINGUISH MULTIPLE
!       MESSAGES FROM MULTI-LINE MESSAGES; THE PROP INFO ACTUALLY REQUIRES ALL
!       MESSAGES FOR AN OBJECT TO BE PRESENT AND CONSECUTIVE.  PROPERTIES WHICH
!       PRODUCE NO MESSAGE SHOULD BE GIVEN THE MESSAGE "<$$<".  NOTE THAT
!       OBJECTS WITH N>100 CANNOT BE USED FOR CONDITIONAL MOTIONS IN
!       TRAVEL TABLE.
!
!  SECTION 6: ARBITRARY MESSAGES.  SAME FORMAT AS SECTIONS 1, 2, AND 5, EXCEPT
!       THE NUMBERS BEAR NO RELATION TO ANYTHING (EXCEPT FOR SPECIAL VERBS
!       IN SECTION 3).
!
!  SECTION 7: CONTAINS LOTS OF OBJECT INFO:
!       (1) OBJECT LOCATIONS & WEIGHTS.  EACH LINE CONTAINS AN OBJECT NUMBER
!       AND ITS INITIAL LOCATION (ZERO (OR OMITTED) IF NONE).  IF THE OBJECT IS
!       IMMOVABLE, THE LOCATION IS FOLLOWED BY A "-1".  IF IT HAS TWO LOCATIONS
!       (E.G. THE GRATE) THE FIRST LOCATION IS FOLLOWED WITH THE SECOND, AND
!       THE OBJECT IS ASSUMED TO BE IMMOVABLE.  IF THE OBJECT IS MOVABLE, IT
!       HAS A THIRD NUMBER WHICH IS ITS RELATIVE WEIGHT.
!       (2) DEFAULT OBJECT NAMES. TEN CHARACTER MAXIMUM, TO PERMIT PRINTING
!       OBJECT NAMES WHEN HERO SAYS 'TAKE ALL' OR SUCH.
!       (3) POINTS.  OBJECT NUMBER, A SCORE
!       VALUE, A PROP VALUE, AND A LOCATION NUMBER.  THE SCORE IS
!       FOR LEAVING IT AT THE NAMED LOCATION WITH THE RIGHT PROP VALUE.
!       SCORES ARE ASSIGNED ON A DIFFICULTY SCALE OF 1-5 (1 IS
!       FOR EASY TREASURES; 5 IS FOR REAL HARD ONES), AND ARE MULTIPLIED
!       BY SOME APPROPRIATE FACTOR AT SCORING TIME.
!
!  SECTION 8: ACTION DEFAULTS.  EACH LINE CONTAINS AN "ACTION-VERB" NUMBER AND
!       THE INDEX (IN SECTION 6) OF THE DEFAULT MESSAGE FOR THE VERB.
!
!  SECTION 9: LIQUID ASSETS, ETC.  EACH LINE CONTAINS A NUMBER (N) AND UP TO 20
!       LOCATION NUMBERS.
!       FOR THE LOW BYTE BIT N (WHERE 0 IS THE UNITS BIT) IS SET IN LOCCON(LOC)
!       FOR EACH LOC GIVEN.  THE COND BITS CURRENTLY ASSIGNED ARE:
!               0       LIGHT
!               1       IF BIT 3 IS ON: ON FOR OIL, OFF FOR WATER
!               2       IF BIT 3 IS ON: ON FOR WINE, OFF FOR WATER & OIL
!               3       LIQUID ASSET, SEE BITS 1 & 2
!               4       PIRATE DOESN'T GO HERE UNLESS FOLLOWING PLAYER
!               5       ALL LOCATIONS IN EITHER 'PORTAL';  I.E., NOT OUTSIDE
!                       CAVE, BUT NOT FAR IN
!               6       ALL LOCATIONS OUTSIDE THE CAVE
!       THE HIGH BYTE IS USED TO INDICATE AREAS OF INTEREST TO "HINT" ROUTINES:
!     THE NUMBER 256*(BIT-7) INDICATES THE HINT
!               1       LOST IN MAZE
!               2       PONDERING DARK ROOM
!               3       AT WITT'S END
!               4       TRYING TO EXTRACT SWORD
!               5       TRYING TO GO UP SLIDE
!               6       TRYING TO GET INTO CAVE VIA SEA ENTRANCE
!               7      TRYING TO FIND CAVE (ANY ENTRANCE)
!               8      TRYING TO CATCH BIRD
!               9       TRYING TO GET OVER THE RAINBOW
!               10      STYMIED BY DOG
!               11      TRYING TO DEAL WITH SNAKE
!               12      TRYING TO GET TO THE CASTLE
!     ONLY ONE HINT PER LOC IS ALLOWED
!       LOCCON(LOC) IS SET TO 2, OVERRIDING ALL OTHER BITS, IF LOC HAS FORCED
!       MOTION.
!
!  SECTION 10: CLASS MESSAGES.  EACH LINE CONTAINS A NUMBER (N), A TAB, AND A
!       MESSAGE DESCRIBING A CLASSIFICATION OF PLAYER.  THE SCORING SECTION
!       SELECTS THE APPROPRIATE MESSAGE, WHERE EACH MESSAGE IS CONSIDERED TO
!       APPLY TO PLAYERS WHOSE SCORES ARE HIGHER THAN THE PREVIOUS N BUT NOT
!       HIGHER THAN THIS N.  NOTE THAT THESE SCORES PROBABLY CHANGE WITH EVERY
!       MODIFICATION (AND PARTICULARLY EXPANSION) OF THE PROGRAM.
!
!  SECTION 11: HINTS.  EACH LINE CONTAINS A HINT NUMBER (CORRESPONDING TO A
!       COND BIT, SEE SECTION 9), THE NUMBER OF TURNS HE MUST BE AT THE RIGHT
!       LOC(S) BEFORE TRIGGERING THE HINT, THE POINTS DEDUCTED FOR TAKING THE
!       HINT, THE MESSAGE NUMBER (SECTION 6) OF THE QUESTION, AND THE MESSAGE
!       NUMBER OF THE HINT.  THESE VALUES ARE STASHED IN THE "HINTS" ARRAY.
!       HNTMAX IS SET TO THE MAX HINT NUMBER (<= HNTSIZ).  NUMBERS 1-6 ARE
!       UNUSABLE SINCE COND BITS ARE OTHERWISE ASSIGNED, SO 2 IS USED TO
!       REMEMBER IF HE'S READ THE CLUE IN THE REPOSITORY, AND 3 IS USED TO
!       REMEMBER WHETHER HE ASKED FOR INSTRUCTIONS (GETS MORE TURNS, BUT LOSES
!       POINTS).  HNTMIN IS SET TO THE NUMBER OF THE FIRST USABLE HINT.
!
!  SECTION 14: PREPOSITION TABLE.  EACH LINE CONTAINS A VERB, A PREPOSITION
!       AND VALID OBJECTS FOR THAT VERB/PREP COMBINATION.  THEY ARE CONVERTED
!       INTO TWO TABLES SIMILAR IN FORMAT TO THE TRAVEL TABLE.  THE FIRST
!       TABLE, VKEY, HAS ONE ENTRY PER VERB.  A ZERO ENTRY INDICATES NO
!       PREPOSITION IS VALID WITH THAT VERB.  A NON-ZERO ENTRY POINTS TO THE
!       BEGINNING OF THE PREP/OBJ LIST FOR THAT VERB.  THE POSITION OF THE
!       ENTRY IN VKEY CORRESPONDS TO THE VERB NUMBER.  THE PREP/OBJ LIST,
!       PTAB, FOR A GIVEN VERB CONSISTS OF A SERIES OF ONE-WORD ENTRIES
!       DELIMITED BY A NEGATIVE ENTRY.  EACH WORD CONTAINS THE PREPOSITION
!       NUMBER TIMES 1000 PLUS THE NUMBER OF A VALID OBJECT.  A SCHEMATIC
!       FOLLOWS.  ENTRIES IN VKEY ARE REPRESENTED BY V1, V2, ... VN.
!       ENTRIES IN PTAB ARE REPRESENTED BY P1, P2, ..., PN (PREPOSITIONS),
!       AND OB1, OB2,..., OBN (OBJECTS).
!
!           I   VKEY(I)         PTAB(J) J
!           -   -------         ------- -
!           1.  V1==========>>  P1,OB1  1.
!                               P1,OB2  2.
!                               P1,OB3  3.
!                               P1,OB4  4.
!                               P1,OB5  5.
!                               P2,OB1  6.
!                               P2,OB2  7.
!                               P3,OB1  8.
!                               P3,OB2  9.
!                               P3,OB3  10.
!                              -P3,OB4  11.
!           2.  V2==========>>  P1,OB1  12.
!                               P1,OB2  13.
!                               P2,OB1  14.
!                              -P3,OB1  15.
!           3.  0  (THE VERB CORRESPONDING TO THIS POSITION TAKES NO PREPOSITION
!           4.  V3==========>>  P1,OB1  16.
!                              -P2,OB1  17.
!
! SECTION 15: OBJECT ATTRIBUTES.  EACH LINE CONTAINS A BIT NUMBER AND
!       UP TO 20 OBJECT NUMBERS.  BIT N (WHERE ZERO IS THE UNITS BIT) IS SET
!       IN OBJCON(OBJ) FOR EACH OBJECT GIVEN.  THE BITS CURRENTLY ASSIGNED
!       ARE:
!               1       THE OBJ CAN BE OPENED/CLOSED. (DOORS, GRATE, CLAM, ETC.)
!               2       THE OBJ IS CURRENTLY OPEN.
!               3       IT HAS A LOCK. (GRATE, CHAIN, ELFIN DOOR, ETC.)
!               4       IT IS CURRENTLY LOCKED.
!               5       FLAMMABLE.  (IT WILL BURN IF IGNITED.)
!               6       IT IS CURRENTLY BURNING.
!               7       EDIBLE.  FOOD, MUSHROOMS, ETC.
!               8       PRINTED MATERIAL, ANYTHING READABLE
!               9       A LIVING BEASTIE, E.G., DWARF, DOG, WUMPUS, ETC.
!               10      DEAD (KILLED) BEASTIE: WUMPUS, DRAGON, SLEEPING DOG
!               11      CAN BE WORN: CROWN, SHOES, CLOAK, JEWELS
!               12      IS CURRENTLY BEING WORN
!               13      REQUIRES PLURAL RESPONSES (SHOES, COINS, ETC.)
!               14      TREASURE.
!               15      CONTAINER.
!               16      OBJECT IS "SMALL". (CAN FIT INTO SACK OR CHEST)
!               17      CONTAINER IS OPAQUE -- CONTENTS ARE NOT VISIBLE UNLESS
!                       CONTAINER IS OPEN.
!
!  SECTION 16: ADJECTIVE/NOUN LIST.  EACH VALID ADJECTIVE IS FOLLOWED BY
!       ALL NOUNS WHICH IT MAY MODIFY.
!
!  SECTION 0: END OF DATABASE.

USE Adventure_Common
USE Adventure_Subs
IMPLICIT NONE

! IMPLICIT INTEGER (a-z)
CHARACTER (LEN=1) :: zapp(20)
CHARACTER (LEN=2) :: kk2c
! The following line added by M.O. August 18, 1990 (var for name of save file)
CHARACTER (LEN=80) :: filnam


! LOGICAL :: ajar,at,athand,bitset,blind,closed,closng,dark,dead,  &
!    edible,enclsd,forced,gaveup,here,hinged,holdng,inside,lmwarn,  &
!    locks,outsid,opaque,panic,pct,plural,portal,printd,scorng,small,  &
!    toting,treasr,locked,vessel,wearng,worn,yea,yesm

!       DATA LINSIZ/25000/,TRVSIZ/1600/,TABSIZ/600/,LOCSIZ/250/,
!     1  VRBSIZ/60/,RTXSIZ/450/,CLSMAX/12/,HNTSIZ/20/,
!     2  MAXOBJ/150/,MAXLOC/300/,HNTMIN/7/,PTBSIZ/300/,ADJSIZ/50/,
!     3  VKYSIZ/60/,BLKLIN/.TRUE./,DWFMAX/6/,ISWIZ/.FALSE./

INTEGER  :: kqqq

!  PHUCE CONSISTS OF FOUR PAIRS OF ORIGIN/DESTINATION LOCATIONS FROM/TO
!  WHICH ONE IS TRANSPORTED ON UTTERING THE ELFIN CURSE AT THE TINY DOOR.
!  HE CAN GO FROM BIG TO SMALL OR SMALL TO BIG, ON EITHER SIDE OF THE DOOR.


!  STATEMENT FUNCTIONS
!
!  AJAR(OBJ     = TRUE IF THE OBJECT IS OPEN
!  AT(OBJ)      = TRUE IF ON EITHER SIDE OF TWO-PLACED OBJECT
!  ATHAND(OBJ)  = TRUE IF OBJECT IS HERE AND NOT IN CLOSED CONTAINER.
!  BITSET(COND,L,N) = TRUE IF COND(L) HAS BIT N SET (BIT 0 IS UNITS BIT)
!  BLIND()      = TRUE IF HERO CAN'T SEE (TOO DARK OR GLAREY)
!  DARK()       = TRUE IF LOCATION "LOC" IS DARK
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



!  CLEAR OUT THE VARIOUS TEXT-POINTER ARRAYS.  ALL TEXT IS STORED IN ARRAY
!  LINES; EACH LINE IS PRECEDED BY A WORD POINTING TO THE NEXT POINTER (I.E.
!  THE WORD FOLLOWING THE END OF THE LINE).  THE POINTER IS NEGATIVE IF THIS IS
!  FIRST LINE OF A MESSAGE.  THE TEXT-POINTER ARRAYS CONTAIN INDICES OF
!  POINTER-WORDS IN LINES.  STEXT(N) IS SHORT DESCRIPTION OF LOCATION N.
!  LTEXT(N) IS LONG DESCRIPTION.  PTEXT(N) POINTS TO MESSAGE FOR PROP( PROPN)=0
!  SUCCESSIVE PROP MESSAGES ARE FOUND BY CHASING POINTERS.  RTEXT CONTAINS
!  SECTION 6'S STUFF.  CTEXT(N) POINTS TO A PLAYER-CLASS MESSAGE.  MTEXT IS FOR
!  SECTION 12.  WE ALSO CLEAR COND.  SEE DESCRIPTION OF SECTION 9 FOR DETAILS.


!     ****** CALL THE ROUTINE TO RESTORE FORM DISK ALL THE COMMONS

!        DEADBT = 10
!        OPENBT = 2
!        LOCKBT = 4
!        BURNBT = 6
!        WEARBT = 12
!       DATA PHUCE/158,160,160,158,167,166,166,167/

OPEN (UNIT=16, FILE='ADVTXT', STATUS='OLD', FORM='UNFORMATTED',  &
      ACCESS='SEQUENTIAL')
READ (16) iswiz,adjkey,adjtab,adjsiz,openbt,lockbt,burnbt,wearbt
READ (16) blklin,loccon,objcon,numdie,maxdie,turns
READ (16) dwarf,knife,knfloc,dflag,dseen,dloc,odloc,dwfmax
READ (16) holder,hlink,hintlc,hinted,hints,hntsiz,hntmin
READ (16) bottle,cask,water,oil,wine,liqtyp
READ (16) loc,oldloc,oldlc2,newloc,maxloc
READ (16) ltext,stext,key,abb,locsiz
READ (16) back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,bear,boat,  &
    book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,lamp,pdoor,  &
    plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,troll2,emrald,  &
    spices,find,yell,invent,leave,pour,say,take,throw,iwest,phuce,tk

READ (16) plac,fixd,weight,prop,points
READ (16) atloc,link,place,fixed,maxobj
READ (16) vkey,ptab,vkysiz,ptbsiz,travel
READ (16) lines,rtext,ptext,wdx,ktab,tabsiz
READ (16) verbs,vrbx,objs,objx,iobjs,iobx,prep,words

READ (16) abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,clock2,  &
    clock3,closed,closng,clsmax,combo,deadbt,detail,dkill,dtotal,  &
    dwarfn,flg239,foo,foobar,food,gaveup,health,hint,hit,hntmax,i,ikk,  &
    iloc,iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,lmwarn,lock,  &
    logout,messag,obj,panic,portal,ptbs,rdflag,retn,rtxsiz,score,  &
    scorng,sect,skey,sloc,spk,start,stick,tabndx,tally,tally2,terse,  &
    trvs,trvsiz,vend,verb,vrbsiz,waste,wkday,wkend,wzdark,yea,actspk,  &
    ctext,cval,hname
READ (16) anvil,batter,bees,billbd,bird,brush,cage,cakes,chain,  &
    chest,chloc,chloc2,clam,cloak,clsses,coins,crown,daltlc,dog,  &
    dragon,eggs,fissur,flower,gatloc,grail,hive,honey,horn,jewels,  &
    keys,lyre,magzin,mirror,mushrm,mxscor,nugget,oyster,pearl,phone,  &
    pillow,pole,poster,prepat,prepdn,prepfr,prepin,prepof,prepon,  &
    pyram,radium,ring,rug,sapphi,shield,shoes,shut,slugs,snake,sphere,  &
    steps,sticks,sword,tablet,tridnt,unlock,vase,wall,wall2,wear,wumpus,y2,yank
READ (16) dtk,atab,vtxt,otxt,iotxt,txt
CLOSE (16)
loc=1
indent=0
!  FINALLY, SINCE WE'RE CLEARLY SETTING THINGS UP FOR THE FIRST TIME...

!  START-UP, DWARF STUFF

i=ranz(1)
CALL rspeak (325)
hinted(3)=yes(65,1,0)
newloc=1
limit=650
flg239=0
IF (hinted(3)) limit=400

!  CAN'T LEAVE CAVE ONCE IT'S CLOSING (EXCEPT BY MAIN OFFICE).

10 IF (.NOT.outsid(newloc) .OR. newloc == 0 .OR. .NOT.closng) GO TO 20
CALL rspeak (130)
newloc=loc
IF (.NOT.panic) clock2=15
panic=.true.

!  SEE IF A DWARF HAS SEEN HIM AND HAS COME FROM WHERE HE WANTS TO GO.  IF SO
!  THE DWARF'S BLOCKING HIS WAY.  IF COMING FROM PLACE FORBIDDEN TO PIRATE
!  (DWARVES ROOTED IN PLACE) LET HIM GET OUT (AND ATTACKED).

20 IF (newloc == loc) GO TO 40
iloc=loc
IF (forced(iloc)) GO TO 40
IF (IAND(loccon(iloc),16) /= 0) GO TO 40
l1=dwfmax-1
DO  i=1,l1
  IF (odloc(i) /= newloc .OR. .NOT.dseen(i)) CYCLE
  newloc=loc
  CALL rspeak (2)
  EXIT
END DO
40 loc=newloc

!  DWARF STUFF.  SEE EARLIER COMMENTS FOR DESCRIPTION OF VARIABLES.  REMEMBER
!  SIXTH DWARF IS PIRATE AND IS THUS VERY DIFFERENT EXCEPT FOR MOTION RULES.

!  FIRST OFF, DON'T LET THE DWARVES FOLLOW HIM INTO A PIT OR A WALL.  ACTIVATE
!  THE WHOLE MESS THE FIRST TIME HE GETS AS FAR AS THE HALL OF MISTS (LOC 15).
!  IF NEWLOC IS FORBIDDEN TO PIRATE (IN PARTICULAR, IF IT'S BEYOND THE TROLL
!  BRIDGE), BYPASS DWARF STUFF.  THAT WAY PIRATE CAN'T STEAL RETURN TOLL, AND
!  DWARVES CAN'T MEET THE BEAR.  ALSO MEANS DWARVES WON'T FOLLOW HIM INTO DEAD
!  END IN MAZE, BUT C'EST LA VIE.  THEY'LL WAIT FOR HIM OUTSIDE THE DEAD END.

IF (loc == 0 .OR. forced(loc) .OR. IAND(loccon(newloc),16) /= 0) GO TO 280
IF (dflag /= 0) GO TO 50
IF (inside(loc)) dflag=1
GO TO 280

!  WHEN WE ENCOUNTER THE FIRST DWARF, WE KILL 0, 1, OR 2 OF THE DWFMAX DWARVES.
!  IF ANY OF THE SURVIVORS IS AT LOC, REPLACE HIM WITH THE ALTERNATE.

50 IF (dflag /= 1) GO TO 80
IF (.NOT.inside(loc) .OR. pct(95)) GO TO 280
dflag=2
DO  i=1,2
  j=1+ranz(dwfmax-1)
  IF (pct(50)) dloc(j)=0
END DO
l1=dwfmax-1
DO  i=1,l1
  IF (dloc(i) == loc) dloc(i)=daltlc
  odloc(i)=dloc(i)
END DO
CALL rspeak (3)
CALL drop (axe, loc)
GO TO 280

!  THINGS ARE IN FULL SWING.  MOVE EACH DWARF AT RANDOM, EXCEPT IF HE'S SEEN US
!  HE STICKS WITH US.  DWARVES NEVER GO TO LOCS WHICH ARE OUTSIDE OR IN
!  EITHER OF THE TWO PORTAL AREAS.  IF WANDERING AT RANDOM, THEY
!  DON'T BACK UP UNLESS THERE'S NO ALTERNATIVE.  IF THEY DON'T HAVE TO
!  MOVE, THEY ATTACK.  AND, OF COURSE, DEAD DWARVES DON'T DO MUCH OF ANYTHING.

80 dtotal=0
attack=0
stick=0
DO  i=1,dwfmax
  IF (dloc(i) == 0) CYCLE
  j=1
  kk=key(dloc(i))
  IF (kk == 0) GO TO 110
  90 newloc=MOD(ABS(travel(kk))/1000,1000)
  IF (newloc > maxloc .OR. newloc == odloc(i)  &
       .OR. .NOT.inside(newloc) .OR. (j > 1 .AND. newloc == tk(j-1))  &
       .OR. j >= 20 .OR. newloc == dloc(i) .OR. forced(newloc)  &
       .OR. (i == dwfmax .AND. IAND(loccon(newloc),16) /= 0)  &
       .OR. ABS(travel(kk))/1000 == 100) GO TO 100
  tk(j)=newloc
  j=j+1
  100 kk=kk+1
  IF (travel(kk-1) >= 0) GO TO 90
  110 tk(j)=odloc(i)
  IF (j >= 2) j=j-1
  j=1+ranz(j)
  odloc(i)=dloc(i)
  dloc(i)=tk(j)
  dseen(i)=(dseen(i) .AND. inside(loc)) .OR. (dloc(i)  == loc .OR. odloc(i) == loc)
  IF (.NOT.dseen(i)) CYCLE
  dloc(i)=loc
  IF (i /= dwfmax) GO TO 170
  
!  THE PIRATE'S SPOTTED HIM.  HE LEAVES HIM ALONE ONCE WE'VE FOUND CHEST.
!  K COUNTS IF A TREASURE IS HERE.  IF NOT, AND TALLY=TALLY2 PLUS ONE FOR
!  AN UNSEEN CHEST, LET THE PIRATE BE SPOTTED.
  
  IF (loc == chloc .OR. prop(chest) >= 0) CYCLE
  k=0
  DO  j=1,maxobj
!  PIRATE WON'T TAKE PYRAMID FROM PLOVER ROOM OR DARK ROOM (TOO EASY!).
    IF (.NOT.treasr(j) .OR. (j == cask .AND. liq(cask) /= wine)) GO TO 120
    IF (j == pyram .AND. (loc == plac(pyram) .OR. loc == plac(emrald)))  &
        GO TO 120
    IF (toting(j) .AND. athand(j)) GO TO 130
    120 IF (here(j) .AND. treasr(j)) k=1
  END DO
  IF (tally == tally2+1 .AND. k == 0 .AND. place(chest) == 0 .AND.  &
      athand(lamp) .AND. prop(lamp) == 1) GO TO 160
  IF (odloc(dwfmax) /= dloc(dwfmax) .AND. pct(30)) CALL rspeak (127)
  CYCLE
  
  130 CALL rspeak (128)
!  DON'T STEAL CHEST BACK FROM TROLL!
  IF (place(messag) == 0) CALL move (chest, chloc)
  CALL move (messag, chloc2)
  DO  j=1,maxobj
    IF (.NOT.treasr(j) .OR. (j == pyram .AND. (loc == plac(pyram)  &
         .OR. loc == plac(emrald))) .OR. (j == cask .AND. liq(cask)  &
         /= wine) .OR. (enclsd(j) .AND. .NOT.athand(j))) CYCLE
    IF (at(j) .AND. fixed(j) == 0) CALL carry (j, loc)
    IF (enclsd(j)) CALL remove (j)
    IF (.NOT.holdng(j)) CYCLE
    CALL insert (j, chest)
    IF (.NOT.wearng(j)) CYCLE
    prop(j)=0
    CALL bitoff (j, wearbt)
  END DO
  150 dloc(dwfmax)=chloc
  odloc(dwfmax)=chloc
  dseen(dwfmax)=.false.
  CYCLE
  
  160 CALL rspeak (186)
  CALL move (chest, chloc)
  CALL move (messag, chloc2)
  GO TO 150
  
!  THIS THREATENING LITTLE DWARF IS IN THE ROOM WITH HIM!
  
  170 dtotal=dtotal+1
  IF (odloc(i) /= dloc(i)) CYCLE
  attack=attack+1
  IF (knfloc >= 0) knfloc=loc
  IF (ranz(1000) < 250*(dflag-2)) stick=stick+1
END DO

!  NOW WE KNOW WHAT'S HAPPENING.  LET'S TELL THE POOR SUCKER ABOUT IT.

IF (dtotal == 0) GO TO 280
IF (dtotal == 1) GO TO 200
WRITE (*,190) dtotal
190 FORMAT (/' There are ', i1,  &
             ' threatening little dwarves in the room with you!')
GO TO 210

200 CALL rspeak (4)
210 IF (attack == 0) GO TO 280
IF (dflag == 2) dflag=3
IF (attack == 1) GO TO 270
WRITE (*,220) attack
220 FORMAT (/' ',i1,' of them throw knives at you!')
k=6
230 IF (stick > 1) GO TO 240
CALL rspeak (k+stick)
IF (stick == 0) GO TO 280
GO TO 260

240 WRITE (*,250) stick
250 FORMAT (/' ',i1,' of them get you!')
260 oldlc2=loc
GO TO 2930

270 CALL rspeak (5)
k=52
GO TO 230

!  DESCRIBE THE CURRENT LOCATION AND (MAYBE) GET NEXT COMMAND.

!  PRINT TEXT FOR CURRENT LOC.

280 IF (loc == 0) GO TO 2930
jkk=stext(loc)
IF (verb == look .OR. jkk == 0 .OR. (.NOT.terse .AND. MOD(abb(loc),  &
    abbnum) == 0)) jkk=ltext(loc)
IF ((forced(loc) .OR. .NOT.dark()) .AND. loc /= 200) GO TO 300
IF (loc /= 200 .AND. (dark() .OR. prop(lamp) == 0 .OR. .NOT.athand(lamp))) &
    GO TO 290
IF (prop(lamp) == 0 .OR. .NOT.athand(lamp)) GO TO 310
IF (pct(35)) GO TO 2920
jkk=rtext(294)
GO TO 310

290 IF (wzdark .AND. pct(35)) GO TO 2920
jkk=rtext(16)
300 IF (holdng(bear) .AND. .NOT.dark()) CALL rspeak (141)
310 CALL speak (jkk)
k=1
abb(loc)=abb(loc)+1
IF (.NOT.forced(loc)) GO TO 320
CALL travl (k, bcross, tally2)
IF (killed) GO TO 2930
GO TO 10

320 abb(loc)=abb(loc)-1
IF (loc == y2 .AND. pct(25) .AND. .NOT.closng) CALL rspeak (8)
IF (loc == 147 .AND. abb(loc) == 1) CALL rspeak (216)

!  SEE IF HE IS WASTING HIS BATTERIES OUT IN THE OPEN.
k=0
IF (.NOT.outsid(loc) .OR. prop(lamp) == 0) GO TO 330
k=waste+1
IF (k <= 12) GO TO 330
CALL rspeak (324)
k=0
330 waste=k

!  IF WUMPUS IS CHASING STOOGE, SEE IF WUMPUS GETS HIM.
IF (chase == 0) GO TO 340
chase=chase+1
kk=chase/2
prop(wumpus)=kk
CALL move (wumpus, loc)
IF (kk < 5) GO TO 340
IF (dark()) CALL rspeak (270)
CALL pspeak (wumpus, 5)
GO TO 2930

!  CHECK FOR RADIATION POISONING.
340 k=1
IF (outsid(loc)) k=3
health=MIN(health+k,100)
IF (.NOT.here(radium) .OR. (place(radium) == -shield .AND.  &
    .NOT.ajar(shield))) GO TO 350
health=health-7
IF (health >= 60) GO TO 350
CALL rspeak (391+(60-health)/10)
IF (health <= 0) GO TO 2930

!  PRINT OUT DESCRIPTIONS OF OBJECTS AT THIS LOCATION.  IF NOT CLOSING AND
!  PROPERTY VALUE IS NEGATIVE, TALLY OFF ANOTHER TREASURE.  RUG IS SPECIAL
!  CASE; ONCE SEEN, ITS PROP IS 1 (DRAGON ON IT) TILL DRAGON IS KILLED.
!  SIMILARLY FOR CHAIN; PROP IS INITIALLY 1 (LOCKED TO BEAR).
!  LIKEWISE, FOR SWORD (MUST PROVE ELFIN ROYALTY).

350 IF (oldloc /= 188 .OR. loc == 189 .OR. loc == 188 .OR. prop(booth) /= 1)  &
    GO TO 360
CALL move (gnome, 0)
prop(booth)=0
360 IF (blind()) GO TO 450
abb(loc)=abb(loc)+1
i=atloc(loc)
370 IF (i == 0) GO TO 450
obj=i
IF (obj > maxobj) obj=obj-maxobj
IF (obj == steps .AND. toting(nugget)) GO TO 390
IF (prop(obj) >= 0) GO TO 380
IF (closed) GO TO 390
prop(obj)=0
IF (obj == rug .OR. obj == chain .OR. obj == sword .OR. obj == cask) prop(obj)=1
IF (obj == cloak .OR. obj == ring) prop(obj)=2
tally=tally-1

!  IF REMAINING TREASURES TOO ELUSIVE, ZAP HIS LAMP.
IF (tally == tally2 .AND. tally /= 0) limit=MIN(35,limit)
380 kk=prop(obj)
IF (obj == steps .AND. loc == fixed(steps)) kk=1
CALL pspeak (obj, kk)
CALL lookin (obj)
390 i=link(i)
GO TO 370

!  "I DON'T UNDERSTAND THAT!"
400 spk=confuz()
GO TO 430

!  "YOU CAN'T DO THAT!"  (AN IMPOSSIBLE ACT, E.G., "OPEN SWORD", "FEED BOAT", ETC.
410 spk=noway()
GO TO 430

420 spk=54
430 IF (obj == 0 .OR. (objs(2) == 0 .AND. iobjs(2) == 0)) GO TO 440
CALL pspeak (obj, -1)
CALL tnoua()
blklin=.false.
440 CALL rspeak (spk)
blklin=.true.

450 rdflag=.false.
IF (objx == 0) GO TO 460
objx=objx+1
IF (objs(objx) == 0) objx=0
460 IF (objx > 0 .AND. objs(objx) /= 0) GO TO 470
IF (iobx == 0) GO TO 470
iobx=iobx+1
IF (iobjs(iobx) == 0) iobx=0
IF (iobx /= 0 .AND. objs(1) /= 0) objx=1

470 IF (objx > 0 .OR. iobx > 0) GO TO 480
IF (objs(1) /= 0) objx=1
IF (iobjs(1) /= 0) iobx=1
vrbx=vrbx+1
IF (verbs(vrbx) /= 0) GO TO 480
CALL clrlin()
rdflag=.true.

!  CHECK IF THIS LOC IS ELIGIBLE FOR ANY HINTS.  IF BEEN HERE LONG ENOUGH,
!  BRANCH TO HELP SECTION (ON LATER PAGE).  HINTS ALL COME BACK HERE EVENTUALLY
!  TO FINISH THE LOOP.  IGNORE "HINTS" < HNTMIN (SPECIAL STUFF, SEE DATABASE
!  NOTES).

480 DO  hint=hntmin,hntmax
  IF (hinted(hint)) CYCLE
  IF ((loccon(loc)/256) /= hint-6) hintlc(hint)=-1
  hintlc(hint)=hintlc(hint)+1
  IF (hintlc(hint) >= hints(hint,1)) GO TO 2670
END DO

!  KICK THE RANDOM NUMBER GENERATOR JUST TO ADD VARIETY TO THE CHASE.  ALSO,
!  IF CLOSING TIME, CHECK FOR ANY OBJECTS BEING TOTED WITH PROP < 0 AND SET
!  THE PROP TO -1-PROP.  THIS WAY OBJECTS WON'T BE DESCRIBED UNTIL THEY'VE
!  BEEN PICKED UP AND PUT DOWN SEPARATE FROM THEIR RESPECTIVE PILES.  DON'T
!  TICK CLOCK1 UNLESS WELL INTO CAVE (AND NOT AT Y2).

500 IF (.NOT.closed) GO TO 520
IF (prop(oyster) < 0 .AND. toting(oyster)) CALL pspeak (oyster, 1)
DO  i=1,maxobj
  IF (toting(i) .AND. prop(i) < 0) prop(i)=-1-prop(i)
END DO
520 wzdark=dark()
IF (knfloc > 0 .AND. knfloc /= loc) knfloc=0
i=ranz(1)
IF (.NOT.rdflag) GO TO 530

!  GET A NEW INPUT CLAUSE, OR FINISH GETTING CURRENT ONE.

CALL getwds()
vrbx=1
objx=0
IF (objs(1) /= 0) objx=1
iobx=0
IF (iobjs(1) /= 0) iobx=1
rdflag=.true.

!  EVERY INPUT, CHECK "FOOBAR" FLAG.  IF ZERO, NOTHING'S GOING ON.  IF POS,
!  MAKE NEG.  IF NEG, HE SKIPPED A WORD, SO MAKE IT ZERO.

530 foobar=MIN(0,-foobar)
combo=MIN(0,-combo)
turns=turns+1
IF (turns == 310 .AND. abbnum /= 10000 .AND. .NOT.terse) CALL rspeak (273)

!  BUMP ALL THE RIGHT CLOCKS FOR RECONNING BATTERY LIFE AND CLOSING.

IF (closed) clock3=clock3-1
IF (clock3 == -7) GO TO 2910
IF (clock3 /= 0) GO TO 540
prop(phone)=0
prop(booth)=0
CALL rspeak (284)
540 IF (tally == 0 .AND. inside(loc) .AND. loc /= y2) clock1=clock1-1
IF (clock1 == 0) GO TO 2810
IF (clock1 < 0) clock2=clock2-1
IF (clock2 == 0) GO TO 2830
IF (prop(lamp) == 1) limit=limit-1
IF (limit == 0) GO TO 2860
IF (limit < 0 .AND. outsid(loc)) GO TO 2870
IF (limit <= 40) GO TO 2850

550 verb=val(verbs(vrbx))
obj=0
IF (objx /= 0) obj=objs(objx)
iobj=0
IF (iobx /= 0) iobj=iobjs(iobx)
IF (knfloc /= loc .OR. (obj /= knife .AND. iobj /= knife)) GO TO 560
knfloc=-1
spk=116
GO TO 430

560 SELECT CASE ( class(verbs(vrbx)) )
  CASE (    1)
    GO TO 590
  CASE (    2)
    GO TO 570
  CASE (    3)
    GO TO 620
  CASE (    4)
    GO TO 580
END SELECT
570 CALL bug (22)

580 spk=verb
GO TO 430

!  IT IS A MOTION VERB.  ANALYZE IT & LOOP TO 2, IF NOT DEAD.
590 CALL travl (verb, bcross, tally2)
IF (killed) GO TO 2930
GO TO 10

!  ACTION VERB 'LEAVE' (DROP) HAS NO OBJECT.
600 CALL bug (29)

!  VERB 'SAY' OR 'YELL' SLIPPED THROUGH WITH AN OBJECT.
610 CALL bug (34)

!  ANALYSE A VERB.
620 spk=actspk(verb)
IF (obj /= 0 .OR. iobj /= 0) GO TO 630

!  ANALYSE AN INTRANSITIVE VERB (IE, NO OBJECT GIVEN YET).

! GO TO (680,640,640,950,420,950,1010,1020,640,640,430,1040,1150,  &
!    1240,1290,640,640,1510,640,1540,640,1680,1720,1740,1760,1780,640,  &
!    640,640,1890,1920,680,640,640,1990,2020,600,2050,640,640,640,640,  &
!    640,680,680,680,680,2300,950,950,2450,2470,2560,2580,2590,2600,  &
!    420,2610,2640),verb

! 01-10    TAKE  DROP   SAY  OPEN  NOTH CLOSE    ON   OFF  WAVE  CALM
! 11-20    WALK  KILL  POUR   EAT DRINK   RUB THROW  QUIT  FIND INVEN
! 21-30    FEED  FILL BLAST SCORE   FOO BRIEF  READ BREAK  WAKE SUSPD
! 31-40   RESUM  YANK  WEAR   HIT ANSWR  BLOW LEAVE  YELL  DIAL  PLAY
! 41-50    PICK   PUT  TURN   GET INSRT REMOV  BURN GRIPE  LOCK UNLOK
! 51-60  HEALTH  LOOK COMBO SWEEP TERSE WIZ   MAP   GATE   PIRLOC

SELECT CASE ( verb )
  CASE (1)
    GO TO 680
  CASE (2,3)
    GO TO 640
  CASE (4)
    GO TO 950
  CASE (5)
    GO TO 420
  CASE (6)
    GO TO 950
  CASE (7)
    GO TO 1010
  CASE (8)
    GO TO 1020
  CASE (9,10)
    GO TO 640
  CASE (11)
    GO TO 430
  CASE (12)
    GO TO 1040
  CASE (13)
    GO TO 1150
  CASE (14)
    GO TO 1240
  CASE (15)
    GO TO 1290
  CASE (16,17)
    GO TO 640
  CASE (18)
    GO TO 1510
  CASE (19)
    GO TO 640
  CASE (20)
    GO TO 1540
  CASE (21)
    GO TO 640
  CASE (22)
    GO TO 1680
  CASE (23)
    GO TO 1720
  CASE (24)
    GO TO 1740
  CASE (25)
    GO TO 1760
  CASE (26)
    GO TO 1780
  CASE (27:29)
    GO TO 640
  CASE (30)
    GO TO 1890
  CASE (31)
    GO TO 1920
  CASE (32)
    GO TO 680
  CASE (33,34)
    GO TO 640
  CASE (35)
    GO TO 1990
  CASE (36)
    GO TO 2020
  CASE (37)
    GO TO 600
  CASE (38)
    GO TO 2050
  CASE (39:43)
    GO TO 640
  CASE (44:47)
    GO TO 680
  CASE (48)
    GO TO 2300
  CASE (49,50)
    GO TO 950
  CASE (51)
    GO TO 2450
  CASE (52)
    GO TO 2470
  CASE (53)
    GO TO 2560
  CASE (54)
    GO TO 2580
  CASE (55)
    GO TO 2590
  CASE (56)
    GO TO 2600
  CASE (57)
    GO TO 420
  CASE (58)
    GO TO 2610
  CASE (59)
    GO TO 2640
  CASE DEFAULT
    CALL bug (23)
END SELECT

!  ANALYSE A TRANSITIVE VERB.

! 630 GO TO (700,880,610,980,420,1000,1010,1020,1030,430,430,1040,1160,  &
!    1270,1320,1360,1370,430,1520,1520,1590,1690,1720,430,430,1800,  &
!    1810,1820,1870,430,430,1930,1940,1980,2000,2020,880,610,2060,2080,  &
!    2090,2100,2160,2170,2180,2270,430,400,2310,2390,400,2470,410,2580,  &
!    400,400,400,400,400),verb

! 01-10    TAKE  DROP   SAY  OPEN  NOTH CLOSE    ON   OFF  WAVE  CALM
! 11-20    WALK  KILL  POUR   EAT DRINK   RUB THROW  QUIT  FIND INVEN
! 21-30    FEED  FILL BLAST SCORE   FOO BRIEF  READ BREAK  WAKE SUSPD
! 31-40    HOUR  YANK  WEAR   HIT ANSWR  BLOW LEAVE  YELL  DIAL  PLAY
! 41-50    PICK   PUT  TURN   GET INSRT REMOV  BURN GRIPE  LOCK UNLOK
! 51-60  HEALTH  LOOK COMBO SWEEP TERSE WIZ    MAP   GATE  PIRLOC

630 SELECT CASE( verb )
  CASE (1)
    GO TO 700
  CASE (2)
    GO TO 880
  CASE (3)
    GO TO 610
  CASE (4)
    GO TO 980
  CASE (5)
    GO TO 420
  CASE (6)
    GO TO 1000
  CASE (7)
    GO TO 1010
  CASE (8)
    GO TO 1020
  CASE (9)
    GO TO 1030
  CASE (10,11)
    GO TO 430
  CASE (12)
    GO TO 1040
  CASE (13)
    GO TO 1160
  CASE (14)
    GO TO 1270
  CASE (15)
    GO TO 1320
  CASE (16)
    GO TO 1360
  CASE (17)
    GO TO 1370
  CASE (18)
    GO TO 430
  CASE (19,20)
    GO TO 1520
  CASE (21)
    GO TO 1590
  CASE (22)
    GO TO 1690
  CASE (23)
    GO TO 1720
  CASE (24,25)
    GO TO 430
  CASE (26)
    GO TO 1800
  CASE (27)
    GO TO 1810
  CASE (28)
    GO TO 1820
  CASE (29)
    GO TO 1870
  CASE (30,31)
    GO TO 430
  CASE (32)
    GO TO 1930
  CASE (33)
    GO TO 1940
  CASE (34)
    GO TO 1980
  CASE (35)
    GO TO 2000
  CASE (36)
    GO TO 2020
  CASE (37)
    GO TO 880
  CASE (38)
    GO TO 610
  CASE (39)
    GO TO 2060
  CASE (40)
    GO TO 2080
  CASE (41)
    GO TO 2090
  CASE (42)
    GO TO 2100
  CASE (43)
    GO TO 2160
  CASE (44)
    GO TO 2170
  CASE (45)
    GO TO 2180
  CASE (46)
    GO TO 2270
  CASE (47)
    GO TO 430
  CASE (48)
    GO TO 400
  CASE (49)
    GO TO 2310
  CASE (50)
    GO TO 2390
  CASE (51)
    GO TO 400
  CASE (52)
    GO TO 2470
  CASE (53)
    GO TO 410
  CASE (54)
    GO TO 2580
  CASE (55:59)
    GO TO 400
  CASE DEFAULT
    CALL bug (24)
END SELECT

!  ROUTINES FOR PERFORMING THE VARIOUS ACTION VERBS

!  STATEMENT NUMBERS IN THIS SECTION ARE 10000 FOR INTRANSITIVE VERBS, 20000 FOR
!  TRANSITIVE, PLUS 100 TIMES THE VERB NUMBER.  MANY INTRANSITIVE VERBS USE THE
!  TRANSITIVE CODE, AND SOME VERBS USE CODE FOR OTHER VERBS, AS NOTED BELOW.

!  RANDOM INTRANSITIVE VERBS COME HERE.  CLEAR OBJ JUST IN CASE (SEE "ATTACK")

640 CALL a5toa1 (vtxt(vrbx,1), vtxt(vrbx,2), '_What?', zapp, k)
WRITE (*,650) zapp(1:k)
650 FORMAT (/' ',20A1)
objs(1)=0
objx=0
GO TO 480

660 CALL a5toa1 (vtxt(vrbx,1), vtxt(vrbx,2), '_it?  ', zapp, k)
WRITE (*,670) zapp(1:k)
670 FORMAT (/' Where do you want to ', 20A1)
GO TO 480


!  CONSTRUCT MSG: "I DON'T KNOW HOW TO [VERB] THE [OBJ]", AND VARIANTS.

!  CARRY, NO OBJECT GIVEN YET.  OK IF ONLY ONE OBJECT PRESENT.

680 IF (atloc(loc) == 0 .OR. link(atloc(loc)) /= 0 .OR. blind()) GO TO 640
l1=dwfmax-1
DO  i=1,l1
  IF (dloc(i) == loc .AND. dflag >= 2) GO TO 640
END DO
obj=atloc(loc)
IF (verb == yank) GO TO 1930
IF (verb == wear) GO TO 1940

!  CARRY AN OBJECT.  SPECIAL CASES FOR BIRD AND CAGE (IF BIRD IN CAGE, CAN'T
!  TAKE ONE WITHOUT THE OTHER.  LIQUIDS ALSO SPECIAL, SINCE THEY DEPEND ON
!  STATUS OF BOTTLE.  ALSO VARIOUS SIDE EFFECTS, ETC.
!  "YANK" AND "WEAR" ALSO WEAVE INTO THIS CODE, SINCE THEY ARE MOSTLY
!  JUST RESTRICTED CARRY'S.

700 IF (obj == boat) spk=281
IF (plural(obj)) spk=297
IF (obj == bird .AND. .NOT.closed .AND. athand(bird) .AND. place(bird)  &
     /= loc) GO TO 710
IF (prep /= prepof) GO TO 720
IF (obj /= 0 .AND. iobj /= 0) GO TO 400
IF (obj == 0) obj=iobj
iobj=0
GO TO 880

710 CALL rspeak (407)
GO TO 10

720 IF (holdng(obj)) GO TO 430
! ASSIGN 730 TO retn
retn = 730
GO TO 870

730 IF (prep == prepin) GO TO 2180
IF (prep == prepfr .OR. enclsd(obj)) GO TO 2270

!  THE NEXT LINES ARE FOR 'TAKING' LIQUIDS (WATER, OIL & WINE).
!  IF WE ARE HOLDING A CONTAINER (BOTTLE OR CASK), WE CAN TAKE THE
!  THE LIQUID BY FILLING THE CONTAINER.  IF THERE IS A CONTAINER NEARBY
!  HOLDING THE REQUESTED LIQUID, WE WILL PICK UP THE CONTAINER.

IF (iobj == 0) GO TO 740
spk=313
IF (obj /= cask .AND. obj /= bottle) GO TO 430
k=0
IF (obj == cask) k=1
iobj=iobj+k
IF (liq(obj) == iobj) GO TO 770
spk=302+k
IF (prop(obj) /= 1) GO TO 430
GO TO 780

740 IF (obj /= water .AND. obj /= oil .AND. obj /= wine) GO TO 790
iobj=obj
k=0
obj=bottle
IF (.NOT.here(bottle)) GO TO 760
IF (prop(bottle) /= 1) GO TO 750
IF (.NOT.here(cask) .OR. (here(cask) .AND. prop(cask) == 1)) GO TO 780
obj=0
CALL rspeak (304)
GO TO 480

750 IF (liq(bottle) == iobj) GO TO 770
760 spk=312
IF (.NOT.here(cask)) GO TO 430
obj=cask
k=1
IF (prop(cask) == 1) GO TO 780
IF (liq(cask) == iobj) GO TO 770
spk=315
IF (.NOT.athand(bottle)) spk=303
GO TO 430

770 IF (.NOT.holdng(obj)) GO TO 790
spk=302+k
GO TO 430

780 IF (holdng(obj)) GO TO 1690
GO TO 790
!  *** END OF LIQUID STUFF

!  'WEAR' AND 'YANK' WEAVE IN HERE.

790 spk=343
IF (obj == bear .OR. burden(0)+burden(obj) <= 15) GO TO 800
spk=92
IF (.NOT.wearng(obj)) GO TO 430
prop(obj)=0
CALL bitoff (obj, wearbt)
GO TO 430

!  CLOAK.  BIG TROUBLE AHEAD.  CAN ONLY GET HERE VIA 'YANK'.
800 IF (obj /= cloak .OR. prop(cloak) /= 2) GO TO 810
prop(rocks)=1
prop(cloak)=0
fixed(cloak)=0
CALL carry (cloak, loc)
CALL rspeak (241)
IF (at(wumpus) .AND. prop(wumpus) == 0) GO TO 1870
GO TO 450

!  POSTER: HIDES WALL SAFE.
810 IF (obj /= poster .OR. place(safe) /= 0) GO TO 820
prop(poster)=1
spk=362
!  MOVE SAFE AND WALL CONTAINING SAFE INTO VIEW.
CALL drop (safe, loc)
CALL drop (wall2, loc)
GO TO 860

!  BOAT: NEED THE POLE TO PUSH IT
820 IF (obj /= boat) GO TO 830
spk=218
IF (.NOT.toting(pole) .AND. place(pole) /= -boat) GO TO 430
prop(boat)=1
spk=221
GO TO 860

!  BIRD: GOT TO HAVE CAGE, BUT ROD CAN'T BE AROUND TO TAKE BIRD
830 IF (obj /= bird .OR. prop(bird) /= 0) GO TO 840
spk=26
IF (athand(rod)) GO TO 430
spk=27
IF (.NOT.holdng(cage)) GO TO 430
CALL insert (bird, cage)
CALL bitoff (cage, openbt)
GO TO 420

!  SWORD: IF IN ANVIL, NEEDS CROWN & MUST YANK.
840 IF (obj /= sword .OR. prop(sword) == 0) GO TO 860
IF (iobj /= 0 .AND. iobj /= anvil) GO TO 410
IF (verb == yank) GO TO 850

!  HE WANTS THE SWORD, BUT HASN'T ESTABLISHED HIS ROYAL BLOOD, OR HE
!  HASN'T PULLED HARD ENOUGH.  OR NEITHER.

IF (.NOT.yes(215,0,0)) GO TO 420
850 IF (wearng(crown)) GO TO 860
CALL pspeak (sword, 2)
IF (closed) GO TO 2880
fixed(sword)=-1
prop(sword)=3
GO TO 450

860 CALL carry (obj, loc)
IF (obj == pole .OR. obj == skey .OR. obj == sword .OR. ((obj == cloak  &
    .OR. obj == ring) .AND. .NOT.wearng(obj))) prop(obj)=0
IF (verb /= yank .OR. obj == sword) GO TO 430
spk=204
GO TO 430

!  THIS IS A QUASI-SUBROUTINE, CALLED FROM 'TAKE' AND FROM 'INSERT', WHEN
!  THE ITEM IS NOT CURRENTLY BEING TOTED.  'RETN' IS A VARIABLE DEFINED
!  TO BE THE RETURN ADDRESS.

870 spk=noway()
IF (obj == plant .AND. prop(plant) <= 0) spk=115
IF (obj == bear .AND. prop(bear) == 1) spk=169
IF (obj == chain .AND. prop(bear) /= 0) spk=170
IF (obj == sword .AND. prop(sword) == 5) spk=208
IF (obj == cloak .AND. prop(cloak) == 2) spk=242
IF (obj == axe .AND. prop(axe) == 2) spk=246
IF (obj == phone) spk=251
IF (obj == bees .OR. obj == hive) spk=295
IF (obj == sticks) spk=296
IF (fixed(obj) /= 0) GO TO 430
! GO TO retn
IF (retn == 730) THEN
  GO TO 730
ELSE IF (retn == 2190) THEN
  GO TO 2190
END IF
WRITE(*, *) 'Stopping - illegal value for retn'
STOP

!  DROP/DISCARD OBJECT.  "THROW" ALSO COMES HERE FOR MOST OBJECTS.
!  SPECIAL CASES FOR BIRD (MIGHT ATTACK SNAKE OR DRAGON) AND CAGE (MIGHT
!  CONTAIN BIRD) AND VASE.
!  DROP COINS IN VENDING MACHINE FOR EXTRA BATTERIES.

880 IF (holdng(rod2) .AND. obj == rod .AND. .NOT.holdng(rod)) obj=rod2
IF (plural(obj)) spk=105
k=liq(bottle)
IF (k == obj) obj=bottle
IF (obj /= bottle) k=liq(cask)
IF (obj /= bottle .AND. k == obj) obj=cask
IF (.NOT.toting(obj)) GO TO 430
IF (prep == prepin) GO TO 2180
IF (obj /= bird .OR. .NOT.here(snake)) GO TO 890
CALL rspeak (30)
IF (closed) GO TO 2880
CALL remove (bird)
CALL dstroy (snake)
!  SET SNAKE PROP FOR USE BY TRAVEL OPTIONS
prop(snake)=1
CALL drop (bird, loc)
GO TO 450

890 spk=344
IF (verb == leave) spk=353
IF (verb == throw) spk=352
IF (verb == take) spk=54
IF (obj /= pole .OR. .NOT.holdng(boat)) GO TO 900
spk=280
GO TO 430

900 IF (obj /= bird .OR. .NOT.at(dragon) .OR. prop(dragon) /= 0) GO TO 910
CALL rspeak (154)
CALL remove (bird)
CALL dstroy (bird)
IF (place(snake) == plac(snake)) tally2=tally2+1
GO TO 450

910 IF (obj /= bear .OR. .NOT.at(troll)) GO TO 920
spk=163
CALL dstroy (troll)
CALL dstroy (troll+maxobj)
CALL move (troll2, plac(troll))
CALL move (troll2+maxobj, fixd(troll))
CALL juggle (chasm)
prop(troll)=2
GO TO 940

920 IF (obj /= vase .OR. loc == plac(pillow)) GO TO 930
prop(vase)=2
IF (at(pillow)) prop(vase)=0
CALL pspeak (vase, prop(vase)+1)
IF (prop(vase) /= 0) fixed(vase)=-1
GO TO 940

930 IF (worn(obj) .OR. obj == pole .OR. obj == boat) prop(obj)=0
IF (worn(obj)) CALL bitoff (obj, wearbt)
IF (obj == pole) prop(boat)=0
940 IF (enclsd(obj)) CALL remove (obj)
CALL drop (obj, loc)
GO TO 430

!  OPEN/CLOSE/LOCK/UNLOCK: NO OBJECT GIVEN.
!  ASSUME VARIOUS THINGS IF PRESENT.

950 spk=28
k=0
DO  i=1,maxobj
  IF (.NOT.(here(i) .AND. hinged(i))) CYCLE
  obj=i
  k=k+1
END DO
IF (k > 1) GO TO 640
IF (obj /= 0) GO TO 970
IF (verb == lock .OR. verb == unlock) GO TO 430
GO TO 640

970 IF (verb == lock) GO TO 2310
IF (verb == unlock) GO TO 2390
IF (verb == shut) GO TO 1000

!  OPEN.   SPECIAL STUFF FOR OPENING CLAM/OYSTER.
!  THE FOLLOWING CAN BE OPENED WITHOUT A KEY:
!       CLAM/OYSTER, DOOR, PDOOR, BOTTLE, CASK, CAGE

980 IF (.NOT.hinged(obj)) GO TO 410
spk=253
IF (obj == pdoor .AND. prop(pdoor) == 1) GO TO 430
spk=336
IF (ajar(obj)) GO TO 430
IF (locks(obj) .OR. iobj == keys .OR. iobj == skey) GO TO 2390
spk=337
IF (obj == door) spk=111
IF (locked(obj)) GO TO 430
IF (obj == clam .OR. obj == oyster) GO TO 990
CALL biton (obj, openbt)
GO TO 420

!  CLAM/OYSTER.
990 k=0
IF (obj == oyster) k=1
spk=124+k
IF (holdng(obj)) spk=120+k
IF (.NOT.athand(tridnt)) spk=122+k
IF (iobj /= 0 .AND. iobj /= tridnt) spk=376+k
IF (spk /= 124) GO TO 430
CALL dstroy (clam)
CALL drop (oyster, loc)
CALL drop (pearl, 105)
GO TO 430

!  CLOSE.  SHUT.
!  THE FOLLOWING CAN BE CLOSED WITHOUT KEYS:
!       DOOR, PDOOR, BOTTLE, CASK, CAGE

1000 IF (.NOT.hinged(obj)) GO TO 410
spk=338
IF (.NOT.ajar(obj)) GO TO 430
IF (locks(obj)) GO TO 2310
CALL bitoff (obj, openbt)
GO TO 420

!  LIGHT LAMP

1010 IF (.NOT.athand(lamp)) GO TO 430
spk=184
IF (limit < 0) GO TO 430
spk=321
IF (prop(lamp) == 1) GO TO 430
prop(lamp)=1
k=39
IF (loc == 200) k=108
CALL rspeak (k)
IF (wzdark) GO TO 280
GO TO 450

!  LAMP OFF

1020 IF (.NOT.athand(lamp)) GO TO 430
spk=322
IF (prop(lamp) == 0) GO TO 430
prop(lamp)=0
CALL rspeak (40)
IF (dark()) CALL rspeak (16)
GO TO 450

!  WAVE.  NO EFFECT UNLESS WAVING ROD AT FISSURE.

1030 IF ((.NOT.holdng(obj)) .AND. (obj /= rod .OR. .NOT.holdng(rod2))) spk= 29
IF (obj /= rod .OR. .NOT.at(fissur) .OR. .NOT.holdng(obj) .OR. closng) GO TO 430
IF (iobj /= 0 .AND. iobj /= fissur) GO TO 430
prop(fissur)=1-prop(fissur)
CALL pspeak (fissur, 2-prop(fissur))
IF (chase == 0 .OR. prop(fissur) /= 0) GO TO 450

!  DEMISE OF THE WUMPUS.  CHAMP MUST HAVE JUST CROSSED BRIDGE.

IF ((loc == 17 .AND. oldloc /= 27) .OR. (loc == 27 .AND. oldloc /= 17)) GO TO 450
CALL rspeak (244)
chase=0
CALL drop (ring, 209)
CALL move (wumpus, 209)
prop(wumpus)=6
CALL biton (wumpus, deadbt)
IF (place(axe) /= plac(wumpus)) GO TO 450
fixed(axe)=0
prop(axe)=0
GO TO 450

!  ATTACK.  ASSUME TARGET IF UNAMBIGUOUS.  "THROW" ALSO LINKS HERE.  ATTACKABLE
!  OBJECTS FALL INTO TWO CATEGORIES: ENEMIES (SNAKE, DWARF, ETC.)  AND OTHERS
!  (BIRD, CLAM).  AMBIGUOUS IF TWO ENEMIES, OR IF NO ENEMIES BUT TWO OTHERS.

!  KILL OBJ WITH IOBJ.

1040 l1=dwfmax-1
DO  dwarfn=1,l1
  IF (dloc(dwarfn) == loc .AND. dflag >= 2) GO TO 1060
END DO
dwarfn=0
1060 IF (obj /= 0) GO TO 1070
IF (dwarfn /= 0) obj=dwarf
IF (here(snake)) obj=obj*maxobj+snake
IF (at(dragon) .AND. prop(dragon) == 0) obj=obj*maxobj+dragon
IF (at(troll)) obj=obj*maxobj+troll
IF (here(gnome)) obj=obj*maxobj+gnome
IF (here(bear) .AND. prop(bear) == 0) obj=obj*maxobj+bear
IF (here(wumpus) .AND. prop(wumpus) == 0) obj=obj*maxobj+wumpus
IF (obj > maxobj) GO TO 640
IF (obj /= 0) GO TO 1070
!  CAN'T ATTACK BIRD BY THROWING AXE.
IF (here(bird) .AND. verb /= throw) obj=bird
!  CLAM AND OYSTER BOTH TREATED AS CLAM FOR INTRANSITIVE CASE; NO HARM DONE.
IF (here(clam) .OR. here(oyster)) obj=maxobj*obj+clam
IF (obj > maxobj) GO TO 640
1070 IF (obj /= bird) GO TO 1080
spk=137
IF (closed) GO TO 430
CALL dstroy (bird)
prop(bird)=0
IF (place(snake) == plac(snake)) tally2=tally2+1
spk=45
1080 IF (obj == dwarf) GO TO 1130
IF (obj == 0) spk=44
IF (obj == clam .OR. obj == oyster) spk=150
IF (at(dog) .AND. prop(dog) == 1) spk=291
IF (obj == snake) spk=46
IF (obj == dragon .OR. (obj == wumpus .AND. prop(wumpus) == 6)) spk= 167
IF (obj == troll) spk=157
IF (obj == bear) spk=165+(prop(bear)+1)/2
IF (obj == gnome) spk=320
IF (iobj /= axe .OR. verb == throw .OR. (obj /= dog .AND. obj /= wumpus  &
    .AND.obj /= dragon .AND. obj /= troll)) GO TO 1090
iobj=obj
obj=iobjs(iobx)
spk=110
GO TO 1370

1090 IF (iobj /= 0 .AND. iobj /= axe) GO TO 410
IF (.NOT.(obj /= dragon .OR. prop(dragon) /= 0)) GO TO 1100
IF (obj /= troll .AND. spk == 158) spk=110
GO TO 430

!  FUN STUFF FOR DRAGON.  IF HE INSISTS ON ATTACKING IT, WIN!  SET PROP TO DEAD,
!  MOVE DRAGON TO CENTRAL LOC (STILL FIXED), MOVE RUG THERE (NOT FIXED), AND
!  MOVE HIM THERE, TOO.  THEN DO A NULL MOTION TO GET NEW DESCRIPTION.
!  THERE IS SOME AMOUNT OF PAIN HERE, TO FORCE GETWDS TO DO THE RIGHT THING.

1100 CALL rspeak (49)
CALL getlin()
wdx=0
CALL clrlin()
IF (txt(1,1) == 'Y     ' .OR. txt(1,1) == 'YES   ') GO TO 1110
words(1)=-2
rdflag=.true.
GO TO 480

1110 CALL pspeak (dragon, 1)
CALL biton (dragon, deadbt)
prop(dragon)=2
prop(rug)=0
k=(plac(dragon)+fixd(dragon))/2
CALL move (dragon+maxobj, -1)
CALL move (rug+maxobj, 0)
CALL move (dragon, k)
CALL move (rug, k)
DO  obj=1,maxobj
  IF (place(obj) == plac(dragon) .OR. place(obj) == fixd(dragon))  &
      CALL move (obj, k)
END DO
words(1)=0
loc=k
newloc=k
GO TO 10

!  HE IS ATTACKING A DWARF.  IF USING SOMETHING OTHER THAN AXE OR SWORD,
!  GOODBYE CHARLIE.  IF USING NOTHING, DON'T LET HIM.  IF USING AXE OR
!  SWORD, THE FOLLOWING ODDS PREVAIL (IF I CALCULATED THIS MESS RIGHT!)
!  (THE END OF LINE FIGURE IS THE CULULATIVE PROBABILITY OF THE EVENT):
!  .25 - HERO KILLS DWARF (.25)
!  .75 - HERO MISSES
!       .25 - HERO GETS KNIFE IN (HIS) RIBS.  DIES. (.1875)
!       .75 - HERO CAN'T MAKE A CLEAN THRUST
!               .36 - STANDOFF (.2)
!               .64 - DWARF SLASHES
!                       .61 - DWARF MISSES! (.22)
!                       .39 - DWARF KILLS HERO (.14)
!  ADVENTURER HAS 1/3 CHANCE OF GETTING NAILED, 1/4 CHANCE OF NAILING
!  DWARF.  ALL BY WAY OF ENCOURAGING HIM TO THROW THE AXE.

1130 IF (obj == dwarf .AND. closed) GO TO 2880
spk=49
IF (iobj == 0) GO TO 430
spk=355
IF (iobj /= axe .AND. iobj /= sword) GO TO 1140
IF (pct(25)) GO TO 1430
IF (pct(25)) GO TO 1140
CALL rspeak (354)
IF (pct(36)) GO TO 450
CALL rspeak (356)
spk=52
IF (pct(61)) GO TO 430
spk=53

!  HERO IS GONZO.
1140 CALL rspeak (spk)
oldlc2=loc
GO TO 2930

!  POUR.  IF NO OBJECT, ASSUME LIQ IN CONTAINER, IF HOLDING ONLY ONE.
!  SPECIAL TESTS FOR POURING WATER OR OIL ON PLANT OR RUSTY DOOR.

1150 IF (.NOT.holdng(bottle) .AND. .NOT.holdng(cask)) GO TO 640
k=liq(bottle)
kk=liq(cask)
IF (holdng(bottle) .AND. k /= 0 .AND. holdng(cask) .AND. kk /= 0) GO TO 640
IF (kk /= 0 .AND. holdng(cask)) obj=cask
IF (k /= 0 .AND. holdng(bottle)) obj=bottle
IF (obj == 0) GO TO 640

!  POUR OBJ FROM IOBJ.
1160 spk=78
IF (obj /= bottle .AND. obj /= cask) GO TO 1170
iobj=obj
obj=liq(iobj)
spk=316
IF (obj == 0) GO TO 430
GO TO 1180

1170 IF (obj < water .OR. obj > wine+1) GO TO 430
spk=29
IF (.NOT.holdng(bottle) .AND. .NOT.holdng(cask)) GO TO 430
IF (holdng(bottle) .AND. liq(bottle) == obj) iobj=bottle
IF (holdng(cask) .AND. liq(cask) == obj) iobj=cask
IF (iobj == 0) GO TO 430

1180 spk=335
IF (.NOT.ajar(iobj)) GO TO 430
IF (iobj == cask) obj=obj+1
prop(iobj)=1
CALL remove (obj)
place(obj)=0
spk=77
IF (iobj /= cask) GO TO 1190
obj=obj-1
spk=104
1190 IF (.NOT.(at(plant) .OR. at(door) .OR. at(sword)) .OR. at(sword)  &
     .AND. prop(sword) == 0) GO TO 430

IF (at(door)) GO TO 1200
IF (at(sword)) GO TO 1220
spk=112
IF (obj /= water) GO TO 430
CALL pspeak (plant, prop(plant)+1)
prop(plant)=MOD(prop(plant)+2,6)
prop(plant2)=prop(plant)/2
newloc=loc
GO TO 10

1200 prop(door)=0
IF (obj /= oil) GO TO 1210
prop(door)=1
CALL bitoff (door, lockbt)
CALL biton (door, openbt)
1210 spk=113+prop(door)
GO TO 430

!  IF SWORD IS ALREADY OILY, DON'T LET HIM CLEAN IT.  NO SOAP.

1220 IF (prop(sword) == 5) GO TO 1230
prop(sword)=4
IF (obj /= oil) GO TO 1230
prop(sword)=5
fixed(sword)=-1
1230 spk=206+prop(sword)-4
GO TO 430

!  EAT.  INTRANSITIVE: ASSUME EDIBLE IF PRESENT, ELSE ASK WHAT.  TRANSITIVE:
!  FOOD/MUSHROOMS/CAKES OK, SOME THINGS LOSE APPETITE, REST ARE RIDICULOUS.
!  IF HE HAS MORE THAN ONE EDIBLE, OR NONE, 'EAT' IS AMBIGUOUS WITHOUT
!  AN EXPLICIT OBJECT.

1240 k=0
DO  i=1,maxobj
  IF (.NOT.(here(i) .AND. edible(i))) CYCLE
  k=k+1
  kk=i
END DO
IF (k /= 1) GO TO 640
obj=kk
IF (obj /= food .AND. obj /= honey) GO TO 1270
1260 IF (obj == honey) tally2=tally2+1
CALL dstroy (obj)
spk=72
GO TO 430

!  IF HE ATE THE RIGHT THING AND IS IN THE RIGHT PLACE, MOVE HIM TO
!  THE OTHER PLACE WITH ALL HIS JUNK.  OTHERWISE, NARKY MESSAGE.
1270 IF (obj == food .OR. obj == honey) GO TO 1260
IF (obj == bird .OR. obj == snake .OR. obj == clam .OR. obj == oyster .OR.  &
    obj == flower) spk=301
IF (obj == dwarf .OR. obj == dragon .OR. obj == troll .OR. obj == dog .OR.  &
    obj == wumpus .OR. obj == bear .OR. obj == gnome) spk=250
IF (obj /= mushrm .AND. obj /= cakes) GO TO 430

k=obj-mushrm
ll=229+k
k=159-k
kk=skey
IF (obj == mushrm) kk=tdoor
IF (obj == mushrm .AND. loc /= 158) tally2=tally2+1
CALL dstroy (obj)
spk=228
IF (.NOT.(here(kk) .OR. fixed(kk) == loc)) GO TO 430
CALL rspeak (ll)
!  IF HE HASN'T TAKEN TINY KEY OFF SHELF, DON'T LET HIM GET IT FOR FREE!
DO  obj=1,maxobj
  IF (obj == skey .AND. prop(skey) == 1) CYCLE
  IF (place(obj) == plac(kk) .AND. fixed(obj) == 0) CALL move (obj, k)
END DO
IF (loc == plac(skey) .AND. place(skey) == plac(skey)) tally2= tally2+1
loc=k
newloc=k
GO TO 10
!  DRINK.  IF NO OBJECT, ASSUME WATER OR WINE AND LOOK FOR THEM HERE.
!  IF POTABLE IS IN BOTTLE OR CASK, DRINK THAT.  IF NOT, SEE IF THERE
!  IS SOMETHING DRINKABLE NEARBY (STREAM, LAKE, WINE FOUNTAIN, ETC.),
!  AND DRINK THAT.  IF HE HAS STUFF IN BOTH CONTAINERS, ASK WHICH.

!  DRINK OBJ FROM IOBJ
1290 ll=liqloc(loc)
IF (.NOT.athand(bottle) .AND. .NOT.athand(cask)  &
     .AND. ll /= wine .AND. ll /= water) GO TO 640
k=liq(bottle)
kk=liq(cask)
IF (.NOT.athand(bottle) .OR. k == 0) GO TO 1300
IF (athand(cask) .AND. kk /= 0 .AND. kk /= k) GO TO 640
obj=k
iobj=bottle
GO TO 1330

1300 IF (.NOT.athand(cask) .OR. kk == 0) GO TO 1310
obj=kk
iobj=cask
GO TO 1330

1310 IF (ll == 0) GO TO 640
obj=ll
iobj=-1
GO TO 1330

1320 IF (obj == 0 .AND. (iobj == bottle .OR. iobj == cask)) obj=liq(iobj)
spk=110
IF (obj == oil) spk=301
IF (obj /= water .AND. obj /= wine) GO TO 430
IF (iobj /= 0) GO TO 1330
IF (obj == liqloc(loc)) iobj=-1
IF (athand(cask) .AND. obj == liq(cask)) iobj=cask
IF (athand(bottle) .AND. obj == liq(bottle)) iobj=bottle
1330 spk=73
IF (iobj == -1) GO TO 1340
IF (iobj == cask) obj=obj+1
CALL remove (obj)
place(obj)=0
prop(iobj)=1
spk=74
IF (iobj == cask) spk=299
1340 IF (obj == water .OR. obj == water+1) GO TO 430

!  UH-OH.  HE'S A WINO.  LET HIM REAP THE REWARDS OF INCONTINENCE.
!  HE'LL WANDER AROUND FOR AWHILE, THEN WAKE UP SOMEWHERE OR OTHER,
!  HAVING DROPPED MOST OF HIS STUFF.

CALL rspeak (300)
IF (prop(lamp) == 1) limit=limit-ranz(limit)/2
IF (limit < 10) limit=25
k=0
IF (pct(15)) k=49
IF (k == 0 .AND. pct(15)) k=53
IF (k == 0 .AND. pct(25)) k=132
IF (k == 0) k=175
IF (outsid(loc)) k=5
IF (k == loc) GO TO 450
IF (holdng(axe)) CALL move (axe, k)
IF (holdng(lamp)) CALL move (lamp, k)
DO  j=1,maxobj
  IF (wearng(j)) CALL bitoff (j, wearbt)
  IF (holdng(j)) CALL drop (j, loc)
END DO
loc=k
newloc=k
GO TO 10

!  RUB.  YIELDS VARIOUS SNIDE REMARKS.

1360 IF (obj /= lamp) spk=76
GO TO 430

!  THROW OBJ AT IOBJ.
!  SAME AS DISCARD UNLESS AXE.  THEN SAME AS ATTACK EXCEPT IGNORE BIRD,
!  AND IF DWARF IS PRESENT THEN ONE MIGHT BE KILLED.
!  AXE ALSO SPECIAL FOR DRAGON, BEAR, DOG, WUMPUS AND TROLL.
!  TREASURES SPECIAL FOR TROLL.
!  IF THROWING FOOD AT SOMEONE WHO MIGHT BE HUNGRY, GO FEED HIM.

1370 IF (prep == prepdn) GO TO 2100
IF (holdng(rod2) .AND. obj == rod .AND. .NOT.holdng(rod)) obj=rod2
IF (.NOT.holdng(obj)) GO TO 430
IF (obj == boat .OR. obj == bear) GO TO 410
dwarfn=0
IF (iobj /= 0) GO TO 1410

!  NO INDIRECT OBJ WAS SPECIFIED.  IF A DWARF IS PRESENT, ASSUME IT
!  IS THE IOBJ.  IF NOT, LOOK FOR ANY OTHER LIVING THING.  IF NO LIVING
!  THINGS PRESENT, TREAT 'THROW' AS 'DROP'.

l1=dwfmax-1
DO  dwarfn=1,l1
  IF (dloc(dwarfn) == loc .AND. dflag >= 2) GO TO 1400
END DO
dwarfn=0

!  NO DWARVES PRESENT; FIGURE OUT PLAUSIBLE OBJECT.

k=0
DO  i=1,maxobj
  IF (.NOT.(at(i) .AND. living(i))) CYCLE
  iobj=i
  k=k+1
END DO
IF (k == 0) GO TO 880

!  IT IS A BEASTIE OF SOME SORT.  IS THERE MORE THAN ONE?
!  DON'T KILL THE BIRD BY DEFAULT.

IF (k == 1) GO TO 1400
CALL rspeak (43)
GO TO 480

1400 IF (iobj == bird) GO TO 880
IF (treasr(obj) .AND. at(troll)) iobj=troll

1410 IF (treasr(obj) .AND. iobj == troll) GO TO 1490
IF (obj == sword .OR. obj == bottle) GO TO 1820
IF (dwarfn /= 0) iobj=dwarf
IF (obj == flower .AND. iobj == hive) iobj=bees
IF (edible(obj) .AND. living(iobj)) GO TO 1590
IF (obj /= axe) GO TO 880
spk=152
IF (iobj == dragon .AND. prop(dragon) == 0) GO TO 1440
spk=158
IF (iobj == troll) GO TO 1440
IF (iobj /= dwarf) GO TO 1450
spk=48
IF (ranz(4) == 0) GO TO 1440
IF (dwarfn /= 0) GO TO 1430
l1=dwfmax-1
DO  dwarfn=1,l1
  IF (dloc(dwarfn) == loc .AND. dflag >= 2) EXIT
END DO

!  'ATTACK' WITH AXE OR SWORD LINKS IN HERE.
1430 dseen(dwarfn)=.false.
dloc(dwarfn)=0
spk=47
dkill=dkill+1
IF (dkill == 1) spk=149

1440 CALL rspeak (spk)
CALL drop (axe, loc)
newloc=loc
GO TO 10

!  THIS'LL TEACH HIM TO THROW THE AXE AT THE BEAR!
1450 IF (iobj /= bear .OR. prop(bear) /= 0) GO TO 1460
spk=164
CALL drop (axe, loc)
fixed(axe)=-1
prop(axe)=1
CALL juggle (bear)
GO TO 430

!  OR THE WUMPUS!
1460 IF (iobj /= wumpus .OR. prop(wumpus) == 6) GO TO 1470
IF (prop(wumpus) == 6) GO TO 880
spk=245
prop(axe)=2
IF (prop(wumpus) == 0) GO TO 1480
spk=243
CALL dstroy (axe)
GO TO 430

!  OR THE NICE DOGGIE!
1470 IF (iobj /= dog .OR. prop(dog) == 1) GO TO 1490
spk=248
prop(axe)=3
1480 CALL drop (axe, loc)
fixed(axe)=-1
CALL juggle (iobj)
GO TO 430

!  SNARF A TREASURE FOR THE TROLL.
1490 IF (iobj /= troll) GO TO 1500
prep=0
IF (obj == cask .AND. liq(cask) /= wine) GO TO 880
spk=159
CALL drop (obj, 0)
IF (obj == cask) place(wine+1)=0
CALL move (troll, 0)
CALL move (troll+maxobj, 0)
CALL drop (troll2, plac(troll))
CALL drop (troll2+maxobj, fixd(troll))
CALL juggle (chasm)
GO TO 430

!  THROWING AXE AT NONE OF THE ABOVE.  ASSUME 'ATTACK'.
1500 obj=iobj
iobj=objs(objx)
GO TO 1040

!  QUIT.  INTRANSITIVE ONLY.  VERIFY INTENT AND EXIT IF THAT'S WHAT HE WANTS.

1510 gaveup=yes(22,54,54)
IF (gaveup) GO TO 2970
GO TO 450

!  FIND.  MIGHT BE CARRYING IT, OR IT MIGHT BE HERE.  ELSE GIVE CAVEAT.

1520 IF (at(obj) .OR. (liq(bottle) == obj .AND. at(bottle))  &
     .OR. k == liqloc(loc)) spk=94
l1=dwfmax-1
DO  i=1,l1
  IF (dloc(i) == loc .AND. dflag >= 2 .AND. obj == dwarf) spk=94
END DO
IF (closed) spk=138
IF (athand(obj)) spk=24
GO TO 430


!  INVENTORY.  IF OBJECT, TREAT SAME AS FIND.  ELSE REPORT ON CURRENT BURDEN.
!  THERE ARE SOME FUNNY CASES, LIKE THE WEARABLE THINGS.  ALSO, BOAT
!  AND BEAR, WHICH AREN'T REALLY CARRIED.  LIST OUTER-LEVEL CONTAINERS
!  AND CONTENTS, IF CONTAINER IS OPEN OR TRANSPARENT.

1540 spk=98
DO  i=1,maxobj
  IF (i == bear .OR. i == boat .OR. .NOT.holdng(i)) CYCLE
  IF (wearng(i)) CYCLE
  IF (spk == 98) CALL rspeak (99)
  blklin=.false.
  CALL pspeak (i, -1)
  spk=0
  IF (i /= boat) CALL lookin (i)
END DO

!  TELL HIM WHAT HE IS WEARING.

k=0
DO  i=1,maxobj
  IF (.NOT.wearng(i)) CYCLE
  IF (k == 0) WRITE (*,1560)
  1560 FORMAT (' You are wearing:')
  CALL tnoua()
  CALL pspeak (i, -1)
  k=-1
END DO

IF (.NOT.holdng(boat)) GO TO 1580
CALL rspeak (221)
CALL lookin (boat)
1580 IF (holdng(bear)) spk=141
GO TO 430

!  FEED.  IF BIRD, NO SEED.  SNAKE, DRAGON, TROLL: QUIP.  IF DWARF, MAKE HIM
!  MAD.  BEAR, SPECIAL.

!  CASE 1: FEED CRITTER.   *OR*
!  CASE 2: FEED CRITTER EDIBLE.
!        [** THIS CASE TRANSFORMED BY PARSER INTO CASE 3 **]
!  CASE 3: FEED EDIBLE TO CRITTER.

1590 IF (iobj /= 0 .AND. living(iobj)) GO TO 1620
spk=100
IF (obj == bird) GO TO 430
IF (.NOT.living(obj)) GO TO 410

!  SEE IF THERE IS ANYTHING EDIBLE AROUND HERE.
kk=0
k=0
DO  i=1,maxobj
  IF (.NOT.here(i) .OR. .NOT.edible(i)) CYCLE
  k=k+1
  kk=i
END DO
iobj=obj
obj=kk
IF (k == 1 .OR. dead(iobj)) GO TO 1620
CALL a5toa1 (otxt(objx,1), otxt(objx,2), '?     ', zapp, k)
WRITE (*,1610) zapp(1:k)
1610 FORMAT (/' What do you want to feed the ', 20A1)
objs(1)=0
objx=0
GO TO 480

!  FEED OBJ TO IOBJ.
1620 IF (iobj /= snake .AND. iobj /= dragon .AND. iobj /= troll) GO TO 1630
spk=102
IF (iobj == dragon .AND. prop(dragon) /= 0) spk=noway()
IF (iobj == troll) spk=182
IF (iobj /= snake .OR. closed .OR. obj /= bird) GO TO 430
spk=101
CALL dstroy (bird)
prop(bird)=0
tally2=tally2+1
GO TO 430

!  FEED DWARF?
1630 IF (iobj /= dwarf) GO TO 1640
spk=103
dflag=dflag+1
GO TO 430

!  FEED BEAR?
1640 spk=102
IF (iobj /= bear) GO TO 1650
IF (prop(bear) == 3) spk=noway()
IF (prop(bear) == 1 .OR. prop(bear) == 2) spk=264
IF (obj == food) spk=278
IF (obj /= honey) GO TO 430
prop(bear)=1
fixed(axe)=0
prop(axe)=0
spk=168
CALL dstroy (honey)
GO TO 430

!  FEED DOG?
1650 IF (iobj /= dog) GO TO 1660
IF (prop(dog) == 1) spk=291
IF (obj /= food .OR. prop(dog) == 1) GO TO 430
spk=249
CALL dstroy (food)
GO TO 430

!  FEED WUMPUS?
1660 IF (iobj /= wumpus) GO TO 1670
IF (prop(wumpus) == 6) spk=326
IF (prop(wumpus) == 0) spk=327
IF (obj == food) spk=240
GO TO 430

!  FEED BEES?
1670 IF (iobj /= bees .OR. obj /= flower) GO TO 410
IF (enclsd(flower)) CALL remove (flower)
CALL drop (flower, loc)
fixed(flower)=-1
prop(flower)=1
CALL drop (honey, loc)
CALL juggle (honey)
spk=267
prop(hive)=1
GO TO 430

!  FILL.  BOTTLE/CASK MUST BE EMPTY, AND SOME LIQUID AVAILABLE.
!  (VASE IS NASTY & GRAIL IS CRACKED.)
1680 IF ((.NOT.here(bottle) .AND. .NOT.here(cask)) .OR. (here(bottle)  &
     .AND. here(cask))) GO TO 640
IF (prop(cask) == 1 .AND. here(cask)) obj=cask
IF (prop(bottle) == 1 .AND. here(bottle)) obj=bottle
IF (obj == 0) GO TO 640

!  FILL OBJ WITH IOBJ
1690 spk=313
IF (.NOT.vessel(obj)) GO TO 430
IF (iobj == 0) iobj=liqloc(loc)
IF (obj /= bottle .AND. obj /= cask) GO TO 1700
k=0
IF (obj == cask) k=1
spk=0
IF (iobj == 0) spk=304+k
IF (liq(obj) /= 0) spk=302+k
IF (spk /= 0) GO TO 430
spk=306+k
IF (iobj == oil) spk=308+k
IF (iobj == wine) spk=310+k
prop(obj)=MOD(loccon(loc),8)/2*2
place(iobj+k)=-1
CALL insert (iobj+k, obj)
GO TO 430

!  VASE.  (NASTY).
1700 IF (obj /= vase) GO TO 1710
spk=144
IF (iobj == 0 .OR. .NOT.holdng(vase)) GO TO 430
CALL rspeak (145)
prop(vase)=2
fixed(vase)=-1
GO TO 940

!  GRAIL OR OTHER.
1710 spk=339
IF (obj == grail) spk=298
GO TO 430

!  BLAST.  NO EFFECT UNLESS YOU'VE GOT DYNAMITE, WHICH IS A NEAT TRICK!

1720 IF (closed) GO TO 1730
GO TO 430

1730 bonus=135
IF (place(rod2) == 212 .AND. loc == 116) bonus=133
IF (place(rod2) == 116 .AND. loc /= 116) bonus=134
CALL rspeak (bonus)
GO TO 2970

!  SCORE.

1740 scorng=.true.
CALL rating (score, bonus, gaveup, scorng, closng, closed, hntmax)
scorng=.false.
WRITE (*,1750) score,mxscor,turns
1750 FORMAT (/' If you were to quit now, you would score',i4,' out of',  &
    i4,' using',i4,' turns.')
!       GAVEUP=YES(143,54,54)
!       GOTO 11850
GO TO 450

!  FEE FIE FOE FOO (AND FUM).  ADVANCE TO NEXT STATE IF GIVEN IN PROPER ORDER.
!  LOOK UP WD1 IN SECTION 3 OF VOCAB TO DETERMINE WHICH WORD WE'VE GOT.   LAST
!  WORD ZIPS THE EGGS BACK TO THE GIANT ROOM (UNLESS ALREADY THERE).

1760 k=vocabx(vtxt(vrbx,1),4)
spk=42
IF (foobar == 1-k) GO TO 1770
IF (foobar /= 0) spk=151
GO TO 430

1770 foobar=k
IF (k /= 4) GO TO 420
foobar=0
IF (place(eggs) == plac(eggs) .OR. (toting(eggs)  &
     .AND. loc == plac(eggs))) GO TO 430
!  BRING BACK TROLL IF WE STEAL THE EGGS BACK FROM HIM BEFORE CROSSING.
IF (place(eggs) == 0 .AND. place(troll) == 0 .AND. prop(troll) == 0)  &
    prop(troll)=1
k=2
IF (here(eggs)) k=1
IF (loc == plac(eggs)) k=0
CALL move (eggs, plac(eggs))
CALL pspeak (eggs, k)
GO TO 450

!  BRIEF/UNBRIEF.  INTRANSITIVE ONLY.
!  SUPPRESS LONG DESCRIPTIONS AFTER FIRST TIME.

1780 detail=3
terse=.false.
IF (abbnum == 10000) GO TO 1790
spk=156
abbnum=10000
GO TO 430

1790 abbnum=5
spk=374
GO TO 430

1800 IF (obj /= 0 .AND. living(obj) .AND. iobj /= 0) GO TO 430
GO TO 400

!  READ.  MAGAZINES IN DWARVISH, MESSAGE WE'VE SEEN, AND . . . OYSTER?

1810 IF (blind()) GO TO 430
IF (obj /= 0 .AND. iobj /= 0) GO TO 400
spk=confuz()
IF (obj == 0) obj=iobj
IF (obj == book .OR. obj == book2) spk=142
IF (obj == billbd) spk=361
IF (obj == carvng) spk=372
IF (obj == magzin) spk=190
IF (obj == messag) spk=191
IF (obj == oyster .AND. hinted(2) .AND. holdng(oyster)) spk=194
IF (obj == poster) spk=370
IF (obj == tablet) spk=196
IF (obj /= oyster .OR. hinted(2) .OR. .NOT.holdng(oyster)  &
     .OR. .NOT.closed) GO TO 430
hinted(2)=yes(192,193,54)
GO TO 450

!  BREAK.  WORKS FOR MIRROR IN REPOSITORY AND, OF COURSE, THE
!  VASE AND BOTTLE.  ALSO, THE SWORD IS MORE BRITTLE THAN IT APPEARS.

1820 IF (obj == mirror) spk=148
IF (obj == vase .AND. prop(vase) == 0) GO TO 1830
IF (obj == bottle .AND. prop(bottle) /= 3) GO TO 1850
IF (obj == sword) GO TO 1860
IF (obj /= mirror .OR. .NOT.closed) GO TO 430
CALL rspeak (197)
GO TO 2880

1830 spk=198
prop(vase)=2
1840 IF (enclsd(obj)) CALL remove (obj)
IF (holdng(obj)) CALL drop (obj, loc)
fixed(obj)=-1
GO TO 430

1850 spk=231
k=liq(bottle)
prop(bottle)=3
IF (k == 0) GO TO 1840
CALL remove (k)
place(k)=0
GO TO 1840

!  HE'D BETTER NOT SLING THE SWORD AROUND!
1860 spk=29
IF (.NOT.holdng(sword)) GO TO 430
spk=279
prop(sword)=4
GO TO 1840

!  ERWACHE.  ONLY USE IS TO DISTURB THE DWARVES OR THE WUMPUS.
!  OTHER WUMPUS-WAKERS LINK HERE.

1870 IF (.NOT.at(wumpus)) GO TO 1880
chase=1
prop(wumpus)=1
spk=276
GO TO 430

1880 IF (at(dog) .AND. prop(dog) == 1) spk=291
IF (obj /= dwarf .OR. .NOT.closed) GO TO 430
CALL rspeak (199)
GO TO 2880

1890 spk=201

! The following lines added by M.O. Aug 28, 1990

WRITE (*,1900, ADVANCE='NO') 'Enter name of save file: '
READ (*,1910) filnam
1900 FORMAT (' ',a)
1910 FORMAT (a)

!     OPEN (UNIT=16,FILE='ASAVE',STATUS='UNKNOWN',FORM=
OPEN (UNIT=16,FILE=filnam,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
      ACCESS='SEQUENTIAL')
WRITE (16) iswiz,openbt,lockbt,burnbt,wearbt
WRITE (16) blklin,loccon,objcon,numdie,maxdie,turns,killed
WRITE (16) dwarf,knife,knfloc,dflag,dseen,dloc,odloc,dwfmax
WRITE (16) holder,hlink
WRITE (16) bottle,cask,water,oil,wine,liqtyp
WRITE (16) loc,oldloc,oldlc2,newloc,maxloc
WRITE (16) ltext,stext,key,abb,locsiz

WRITE (16) plac,fixd,weight,prop,points
WRITE (16) atloc,link,place,fixed,maxobj
WRITE (16) verbs,vrbx,objs,objx,iobjs,iobx,prep,words

WRITE (16) abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,  &
    clock2,clock3,closed,closng,clsmax,combo,deadbt,detail,dkill,  &
    dtotal,dwarfn,flg239,foo,foobar,food,gaveup,health,hint,hit,  &
    hntmax,i,ikk,iloc,iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,  &
    lmwarn,lock,logout,messag,obj,panic,portal,ptbs,rdflag,retn,  &
    rtxsiz,score,scorng,sect,skey,sloc,spk,start,stick,tabndx,tally,  &
    tally2,terse,trvs,trvsiz,vend,verb,vrbsiz,waste,wkday,wkend,  &
    wzdark,yea,actspk,ctext,cval,hname
CLOSE (16)

GO TO 450

!       THE ZSAVE FUNCTION MUST BE CHANGED IF ANY NEW VARIABLES
!       WHICH CHANGE DURING EXECUTION ARE ADDED.


!   RESUME   RESUME GAME BROM BACKUP FILE

!     ****************************
!     SEE THE COMMENT ON SUSPEND ABOVE!!!!!!!!!!!


! The following lines added by M.O. Aug 28, 1990

1920 WRITE (*,1900, ADVANCE='NO') 'Enter name of save file: '
READ (*,1910) filnam

!     OPEN (UNIT=16,FILE='ASAVE',STATUS='OLD',FORM=
OPEN (UNIT=16, FILE=filnam, STATUS='OLD', FORM='UNFORMATTED',   &
      ACCESS='SEQUENTIAL')
REWIND (16)
READ (16) iswiz,openbt,lockbt,burnbt,wearbt
READ (16) blklin,loccon,objcon,numdie,maxdie,turns
READ (16) dwarf,knife,knfloc,dflag,dseen,dloc,odloc,dwfmax
READ (16) holder,hlink
READ (16) bottle,cask,water,oil,wine,liqtyp
READ (16) loc,oldloc,oldlc2,newloc,maxloc
READ (16) ltext,stext,key,abb,locsiz

READ (16) plac,fixd,weight,prop,points
READ (16) atloc,link,place,fixed,maxobj
READ (16) verbs,vrbx,objs,objx,iobjs,iobx,prep,words

READ (16) abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,clock2,  &
    clock3,closed,closng,clsmax,combo,deadbt,detail,dkill,dtotal,  &
    dwarfn,flg239,foo,foobar,food,gaveup,health,hint,hit,hntmax,i,ikk,  &
    iloc,iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,lmwarn,lock,  &
    logout,messag,obj,panic,portal,ptbs,rdflag,retn,rtxsiz,score,  &
    scorng,sect,skey,sloc,spk,start,stick,tabndx,tally,tally2,terse,  &
    trvs,trvsiz,vend,verb,vrbsiz,waste,wkday,wkend,wzdark,yea,actspk,  &
    ctext,cval,hname
CLOSE (16)

GO TO 450


!  YANK.  A VARIANT OF 'CARRY'.  IN GENERAL, NOT A GOOD IDEA.
!  AT MOST, IT GETS THE CLOAK OR A COUPLE SNIDE COMMENTS.

1930 IF (toting(obj)) GO TO 880
IF (obj == plant .OR. obj == sword .OR. obj == chain) GO TO 700
spk=205
IF (obj == bear .AND. prop(chain) == 1) GO TO 430
IF (obj == cloak .AND. prop(cloak) == 2) GO TO 790
GO TO 700

!  WEAR: ONLY GOOD FOR JEWELS, RUBY SLIPPERS, CLOAK & CROWN.
!  BUT HE MIGHT TRY THE SWORD.  ANYTHING ELSE IS RIDICULOUS.
!  ANOTHER VARIANT OF 'CARRY'.

1940 spk=209
IF (obj == sword .AND. prop(sword) /= 3) GO TO 430
IF (worn(obj)) GO TO 1970
1950 CALL a5toa1 (otxt(objx,1), otxt(objx,2), '?     ', zapp, k)
WRITE (*,1960) (zapp(i),i=1,k)
1960 FORMAT (/' Just exactly how does one wear a ',20A1)
GO TO 450

1970 spk=242
IF (obj == cloak .AND. prop(cloak) == 2) GO TO 430
spk=210
IF (obj == shoes) spk=227
IF (wearng(obj)) GO TO 430
prop(obj)=1
CALL biton (obj, wearbt)
IF (enclsd(obj)) CALL remove (obj)
IF (holdng(obj)) GO TO 420
GO TO 790

!  HIT.  IF NOT PUNCHING OUT TELEPHONE, ASSUME ATTACK.

1980 IF (at(wumpus) .AND. prop(wumpus) == 0) GO TO 1870
IF (obj /= phone) GO TO 1040
IF (closed) GO TO 2890
spk=256
IF (prop(phone) == 2) GO TO 430
CALL drop (slugs, loc)
spk=257
prop(phone)=2
prop(booth)=2
GO TO 430

!  ANSWER (TELEPHONE).  SMARTASS FOR ANYTHING ELSE.
1990 IF (loc /= 189 .OR. prop(phone) /= 0) GO TO 640
obj=phone
2000 IF (obj /= phone) GO TO 2010
spk=269
IF (prop(phone) /= 0) GO TO 430
IF (closed) GO TO 2900
spk=261
prop(phone)=1
prop(booth)=2
GO TO 430

2010 IF (obj == dwarf .OR. obj == wumpus .OR. obj == snake .OR.  &
         obj == bear .OR. obj == dragon) spk=259
IF (obj == troll) spk=258
IF (obj == bird) spk=260
GO TO 430

!  BLOW.  JOSHUA FIT DE BATTLE OF JERICHO, AND DE WALLS....

2020 IF (obj /= 0 .AND. iobj /= 0) GO TO 400
IF (obj == 0) obj=iobj
iobj=0
IF (obj == 0) spk=268
IF (obj /= horn) GO TO 430
spk=266
IF (outsid(loc)) spk=277
IF (.NOT.at(wumpus)) GO TO 2030
IF (prop(wumpus) /= 0) GO TO 430
CALL rspeak (spk)
GO TO 1870

2030 IF (prop(wall) == 1 .OR. (loc /= 102 .AND. loc /= 194)) GO TO 430
k=196
IF (loc == 194) k=195
CALL rspeak (265)
prop(wall)=1
DO  obj=1,maxobj
  IF (place(obj) == loc .OR. fixed(obj) == loc) CALL move (obj, k)
END DO
newloc=k
GO TO 10

!  CALL.  IF NO PHONE IS HANDY, YELL.

2050 IF (.NOT.here(phone)) GO TO 640
GO TO 2070

!  DIAL.  NO EFFECT UNLESS AT PHONE.

2060 IF (obj /= phone) GO TO 430
2070 IF (closed) GO TO 2900
spk=271
GO TO 430

!  PLAY.  ONLY FOR HORN OR LYRE.
2080 IF (obj /= 0 .AND. iobj /= 0) GO TO 400
IF (obj == 0) obj=iobj
iobj=0
IF (obj == horn) GO TO 2020
IF (obj /= lyre) GO TO 430
spk=287
IF (.NOT.here(dog) .OR. dead(dog)) GO TO 430
prop(dog)=1
CALL biton (dog, deadbt)
fixed(axe)=0
prop(axe)=0
spk=288
GO TO 430

!  PICK/PICK UP.  CAN PICK FLOWERS & MUSHROOMS, BUT MUST 'PICK UP' EVERY EL

2090 IF (obj == 0) obj=iobj
iobj=0
IF (obj == flower .OR. obj == mushrm) GO TO 700
IF (prep /= 0) GO TO 700
GO TO 400
!  PUT DOWN: EQUIVALENT TO DROP.
!  PUT IN: IF LIQUID, MEANS 'FILL'.
!  PUT ON: WEAR OR DROP.

2100 IF (prep /= 0) GO TO 2120
CALL a5toa1 (otxt(objx,1), otxt(objx,2), '?     ', zapp, k)
WRITE (*,2110) (zapp(i),i=1,k)
2110 FORMAT (/' Where do you want to put the ',20A1)
GO TO 480

2120 IF (prep == prepin) GO TO 2180

!  PUT ON: "WEAR" OR "PUT OBJ ON IOBJ".
IF (prep /= prepon) GO TO 2140
IF (obj /= 0) GO TO 2130
obj=iobj
otxt(objx,1)=iotxt(iobx,1)
otxt(objx,2)=iotxt(iobx,2)
iobj=0
2130 IF (worn(obj)) GO TO 1940
IF (iobj == 0) GO TO 1950
GO TO 880

!  PUT DOWN: "DROP"
2140 IF (obj == 0 .OR. iobj == 0) GO TO 2150
GO TO 410

2150 IF (obj == 0) obj=iobj
iobj=0
GO TO 880


!  TURN ON/OFF.

2160 IF (prep == 0) GO TO 400
IF (obj == 0 .AND. iobj == lamp) obj=lamp
IF (obj /= lamp) GO TO 410
IF (prep == prepon) GO TO 1010
GO TO 1020


!  GET (NO PREP): "TAKE"
!  GET IN: "ENTER"
!  GET OUT: "LEAVE"
!****** NEEDS WORK

2170 IF (prep == 0 .OR. prep == prepfr) GO TO 700
IF (obj /= 0) GO TO 400
obj=iobj
iobj=0
prep=0
GO TO 700

!  INSERT/PUT IN.

2180 IF (iobj == 0) GO TO 660
spk=noway()
IF (obj == sword .AND. iobj == anvil .AND. prop(sword) == 0) spk=350
IF (.NOT.vessel(iobj)) GO TO 430
! ASSIGN 2190 TO retn
retn = 2190
GO TO 870

2190 IF (iobj /= bottle .AND. iobj /= cask .AND. iobj /= vase .AND.  &
         iobj /= grail .AND. (obj < water .OR. obj > wine+1)) GO TO 2200
obj=iobj
iobj=objs(objx)
GO TO 1690

2200 spk=252
IF (obj == iobj) GO TO 430
spk=358
IF (.NOT.ajar(iobj)) GO TO 430
IF (iobj /= boat .AND. iobj /= chest) GO TO 2210
IF (iobj == chest .AND. obj == boat) GO TO 410
GO TO 2260

!  BIRD GOES INTO CAGE AND ONLY CAGE.  CAGE HOLDS NOTHING ELSE.
!  BAR VASE & PILLOW FROM SAFE, TO FORCE PUTTING DOWN ON FLOOR.
2210 spk=351
IF (obj == bird .AND. iobj /= cage) GO TO 430
spk=329
IF (iobj == cage .AND. obj /= bird) GO TO 430
IF (obj == bird) GO TO 790
IF (iobj == safe .AND. (obj == vase .OR. obj == pillow)) GO TO 430
IF (iobj == shield .AND. obj /= radium) GO TO 430
IF (iobj /= phone) GO TO 2220
IF (obj /= coins .AND. obj /= slugs) GO TO 410
CALL dstroy (obj)
spk=330
GO TO 430

2220 IF (iobj /= vend) GO TO 2240
IF (obj /= coins .AND. obj /= slugs) GO TO 410
CALL dstroy (obj)
CALL move (batter, loc)
IF (prop(batter) /= 1) GO TO 2230
CALL rspeak (317)
prop(vend)=1
2230 prop(batter)=0
CALL pspeak (batter, 0)
GO TO 450

!  PUT BATTERIES IN LAMP.
!  THERE IS A GLITCH HERE, IN THAT IF HE TRIES TO GET A THIRD SET OF
!  BATTERIES BEFORE THE SECOND SET HAS BEEN INSERTED, THE SECOND SET
!  DISAPPEARS!  ***FIX THIS SOMETIME***
2240 IF (iobj /= lamp) GO TO 2250
IF (obj /= batter .OR. prop(batter) /= 0) GO TO 410
prop(batter)=1
IF (enclsd(batter)) CALL remove (batter)
IF (holdng(batter)) CALL drop (batter, loc)
limit=400
prop(lamp)=1
lmwarn=.false.
spk=188
GO TO 430

2250 IF (.NOT.small(obj)) GO TO 430
2260 IF (wearng(obj)) CALL bitoff (obj, wearbt)
IF (worn(obj)) prop(obj)=0
IF (enclsd(obj)) CALL remove (obj)
CALL insert (obj, iobj)
GO TO 420

!  REMOVE/TAKE FROM.

2270 IF (obj /= ring .OR. prop(ring) /= 2) GO TO 2280
prep=0
iobj=0
GO TO 700

2280 spk=343
IF (iobj /= 0) GO TO 2290
IF (.NOT.enclsd(obj)) spk=340
iobj=-place(obj)
2290 IF (place(obj) /= -iobj) spk=341
IF (.NOT.ajar(iobj)) spk=335
IF (obj == water .OR. obj == oil .OR. obj == wine) spk=342
IF (.NOT.toting(obj) .AND. burden(0)+burden(obj) > 15) spk=92
IF (spk /= 343) GO TO 430
CALL remove (obj)
IF (obj == bird) GO TO 880
GO TO 420


!  BURN


!  GRIPE/COMPLAIN/SUGGEST.

!       A SUGGESTION MECHANISM USED TO GO HERE.  IT IS NOT NEEDED FOR A
!       SINGLE-USER GAME
2300 GO TO 430

!  LOCK.  CHAIN, GRATE, CHEST, ELFIN DOOR
!  HERE ARE THE CURRENT LOCK/UNLOCK MESSAGES & NUMBERS:
!       31      YOU HAVE NO KEYS.
!       32      IT HAS NO LOCK.
!       34      IT'S ALREADY LOCKED.
!       35      THE GRATE IS NOW LOCKED.
!       36      THE GRATE IS NOW UNLOCKED.
!       37      IT WAS ALREADY UNLOCKED.
!       55      YOU CAN'T UNLOCK THE KEYS.
!       171     THE CHAIN IS NOW UNLOCKED.
!       172     THE CHAIN IS NOW LOCKED.
!       173     THERE IS NOTHING HERE TO WHICH THE CHAIN CAN BE LOCKED.
!       224     YOUR KEYS ARE ALL TOO LARGE.
!       234     THE WROUGHT-IRON DOOR IS NOW LOCKED.
!       235     THE TINY DOOR IS NOW LOCKED.
!       236     THE WROUGHT-IRON DOOR IS NOW UNLOCKED.
!       237     THE TINY DOOR IS NOW UNLOCKED.
!       375     YOU DON'T HAVE THE RIGHT KEY.
!       333     THE CHEST IS NOW LOCKED.
!       334     THE CHEST IS NOW UNLOCKED.
!       367     THE SAFE'S DOOR SWINGS SHUT.

2310 IF (hinged(obj)) GO TO 2330
CALL a5toa1 (otxt(objx,1), otxt(objx,2), '.     ', zapp, k)
WRITE (*,2320) (zapp(i),i=1,k)
2320 FORMAT (/' I DON''T KNOW HOW TO LOCK OR UNLOCK THE ',20A1)
GO TO 450

2330 spk=375
IF (.NOT.locks(obj)) spk=32
IF (locked(obj)) spk=34
IF (.NOT.(athand(keys) .OR. athand(skey) .OR. obj == safe)) spk=31
IF (spk /= 375) GO TO 430

!  CHAIN.
IF (obj /= chain) GO TO 2340
IF (.NOT.athand(keys)) GO TO 430
spk=173
IF (loc /= plac(chain)) GO TO 430
spk=172
prop(chain)=2
IF (enclsd(chain)) CALL remove (chain)
IF (holdng(chain)) CALL drop (chain, loc)
fixed(chain)=-1
GO TO 2380

!  CHEST.
2340 IF (obj /= chest) GO TO 2350
IF (.NOT.athand(keys)) GO TO 430
spk=334
GO TO 2380

!  ELFIN DOOR.
2350 IF (obj /= tdoor .AND. obj /= tdoor2) GO TO 2360
spk=224
IF (.NOT.toting(skey)) GO TO 430
prop(tdoor)=0
prop(tdoor2)=0
spk=234+2*prop(tdoor)+(tdoor2-obj)
k=tdoor+(tdoor2-obj)
CALL biton (k, lockbt)
CALL bitoff (k, openbt)
GO TO 2380

!  GRATE.
2360 IF (obj /= grate) GO TO 2370
IF (.NOT.athand(keys)) GO TO 430
prop(grate)=0
spk=35
GO TO 2380

!  SAFE.
2370 prop(safe)=0
spk=367
!       GOTO 24990

2380 CALL biton (obj, lockbt)
CALL bitoff (obj, openbt)
GO TO 430

!  UNLOCK.  CHAIN, GRATE, CHEST, ELFIN DOOR.

2390 spk=55
IF (obj == keys .OR. obj == skey) GO TO 430
IF (.NOT.hinged(obj)) GO TO 2310
spk=375
IF (.NOT.locked(obj)) spk=37
IF (.NOT.locks(obj)) spk=32
IF (obj == safe) spk=342
IF (obj == safe .AND. (iobj == keys .OR. iobj == skey)) spk=368
IF (.NOT.(athand(keys) .OR. athand(skey) .OR. obj == safe)) spk=31
IF (spk /= 375) GO TO 430

!  CHAIN.
IF (obj /= chain) GO TO 2400
IF (.NOT.athand(keys)) GO TO 430
spk=171
IF (prop(bear) == 0) spk=41
IF (spk /= 171) GO TO 430
prop(chain)=0
fixed(chain)=0
IF (prop(bear) /= 3) prop(bear)=2
fixed(bear)=2-prop(bear)
GO TO 2430

!  CHEST.
2400 IF (obj /= chest) GO TO 2410
IF (.NOT.athand(keys)) GO TO 430
spk=333
GO TO 2430

!  ELFIN DOOR.
!  STUFF TO LOCK/UNLOCK TINY DOOR W/SPECIAL KEY.
!  THE DAMN THING IS REALLY AT FOUR PLACES, AND WE WANT THE RIGHT
!  MESSAGES IF HE ONLY HAS 'BIG' KEYS (OR NO KEYS).  ALSO, HE
!  CAN UNLOCK IT EITHER WHILE HE IS BIG OR SMALL.
2410 IF (obj /= tdoor .AND. obj /= tdoor2) GO TO 2420
spk=224
IF (.NOT.athand(skey)) GO TO 430
IF (closng) GO TO 2440
prop(tdoor)=1
prop(tdoor2)=1
spk=234+2*prop(tdoor)+(tdoor2-obj)
k=tdoor+(tdoor2-obj)
CALL bitoff (k, lockbt)
CALL biton (k, openbt)
GO TO 2430

!  GRATE.
2420 IF (.NOT.athand(keys)) GO TO 430
IF (closng) GO TO 2440
prop(grate)=1
spk=36

2430 CALL bitoff (obj, lockbt)
CALL biton (obj, openbt)
GO TO 430

!  CLOSING.  NO EXIT THIS WAY.
2440 spk=130
IF (.NOT.panic) clock2=15
panic=.true.
GO TO 430


!  HEALTH.  GIVE HIM A DIAGNOSIS.

2450 IF (health < 100) WRITE (*,2460) health
2460 FORMAT (/' You''RE HEALTH RATING IS ',i2,' OUT OF A POSSIBLE 100.' )
IF (pct(50)) spk=349
IF (health >= 95) GO TO 430
spk=381+(100-health)/20
GO TO 430


!  LOOK.  CAN'T GIVE MORE DETAIL.  PRETEND IT WASN'T DARK (THOUGH IT MAY "NOW"
!  BE DARK) SO HE WON'T FALL INTO A PIT WHILE STARING INTO THE GLOOM.

2470 IF (obj /= 0) GO TO 400
IF (iobj /= 0) GO TO 2490
2480 IF (detail < 3) CALL rspeak (15)
detail=detail+1
wzdark=.false.
abb(loc)=0
newloc=loc
GO TO 10

!  LOOK INTO SOMETHING (A CONTAINER).
2490 IF (.NOT.vessel(iobj)) GO TO 2510
IF (.NOT.ajar(iobj) .AND. opaque(iobj)) GO TO 430
spk=359
IF (holder(iobj) == 0) GO TO 430
WRITE (*,2500)
2500 FORMAT (' ')
CALL lookin (iobj)
GO TO 450

!  LOOK AT SOMETHING.  IF WRITTEN, READ IT.
2510 IF (.NOT.printd(iobj)) GO TO 2520
obj=iobj
iobj=0
GO TO 1810

2520 IF (iobj /= sphere) GO TO 2480
IF (inside(loc) .AND. (.NOT.athand(sapphi))) GO TO 2530
CALL rspeak (42)
GO TO 450
2530 CALL rspeak (400)
WRITE (*,*) '  '
sloc=place(sapphi)
IF ((MOD(loccon(sloc),2) == 0 .OR. enclsd(sapphi))  &
     .AND. sloc /= 200 .AND. .NOT.(place(lamp) == sloc .AND. prop(lamp)  &
     /= 0)) GO TO 2540
CALL speak (ltext(sloc))
IF (sloc /= 239 .OR. flg239 /= 0) GO TO 2550
CALL rspeak (403)
flg239=1
GO TO 2550
2540 CALL rspeak (401)
2550 WRITE (*,*) '  '
CALL rspeak (402)
GO TO 450

!  COMBO: TRYING TO OPEN SAFE.  SEE COMMENTS FOR FEE FIE FOE FOO.

2560 IF (.NOT.at(safe)) GO TO 400
k=vocabx(vtxt(vrbx,1),4)-10
spk=42
IF (combo == 1-k) GO TO 2570
IF (combo /= 0) spk=366
GO TO 430

2570 combo=k
spk=371
IF (k /= 3) GO TO 430
combo=0
CALL bitoff (safe, lockbt)
CALL biton (safe, openbt)
prop(safe)=1
spk=365
IF (prop(book) >= 0) GO TO 430
tally=tally-1
prop(book)=0
!  IF REMAINING TREASURES TOO ELUSIVE, ZAP HIS LAMP.
!  THIS COPIES SOME CODE LOCATED AROUND LABEL 2000.  MUST BE DONE
!  HERE SINCE BOOK IS CONTAINED IN SAFE & TALLY STUFF ONLY WORKS FOR THINGS
!  DEPOSITED AT A LOC.
IF (tally == tally2 .AND. tally /= 0) limit=MIN(35,limit)
GO TO 430

!  DUST/SWEEP.
2580 IF (.NOT.athand(brush)) spk=342
IF (.NOT.at(carvng) .OR. .NOT.athand(brush) .OR. prop(carvng) == 1) GO TO 430
prop(carvng)=1
CALL rspeak (363)
spk=372
GO TO 430


!  TERSE/UNTERSE.  SUPRESS ALL LONG-FORM DESCRIPTIONS.

2590 terse=.NOT.terse
detail=3
GO TO 420
!  WIZ STUFF

! ???????????? UNCOMMENT THE FOLLOWING LINE TO ACTIVATE WIZARD MODE:
!           ISWIZ = .NOT.ISWIZ
2600 GO TO 420

2610 IF (.NOT.iswiz) GO TO 420
WRITE (*, 2620)
2620 FORMAT (' LOCATION ?')
READ (*,2630) gatloc
2630 FORMAT (i4)
IF (gatloc > 0 .AND. gatloc <= maxloc) loc=gatloc
GO TO 420

2640 IF (iswiz) WRITE (*,2650) (dloc(kqqq),kqqq=1,dwfmax-1)
2650 FORMAT (' THE DWARFS ARE AT LOCATIONS '/ 8I6)
IF (iswiz) WRITE (*,2660) dloc(dwfmax)
2660 FORMAT (' THE PIRATE IS AT LOCATION ',i4)
GO TO 420

!  HINTS

!  COME HERE IF HE'S BEEN LONG ENOUGH AT REQUIRED LOC(S) FOR SOME UNUSED HINT.
!  HINT NUMBER IS IN VARIABLE "HINT".  BRANCH TO QUICK TEST FOR ADDITIONAL
!  CONDITIONS, THEN COME BACK TO DO NEAT STUFF.  GOTO 40010 IF CONDITIONS ARE
!  MET AND WE WANT TO OFFER THE HINT.  GOTO 40020 TO CLEAR HINTLC BACK TO ZERO,
!  40030 TO TAKE NO ACTION YET.

2670 SELECT CASE ( (hint+1-hntmin) )
  CASE (    1)
    GO TO 2720
  CASE (    2)
    GO TO 2730
  CASE (    3)
    GO TO 2690
  CASE (    4)
    GO TO 2740
  CASE (    5)
    GO TO 2690
  CASE (    6)
    GO TO 2750
  CASE (    7)
    GO TO 2760
  CASE (    8)
    GO TO 2770
  CASE (    9)
    GO TO 2780
  CASE (   10)
    GO TO 2790
  CASE (   11)
    GO TO 2800
  CASE (   12)
    GO TO 2690
END SELECT
!             MAZE  DARK  WITT  SWORD SLIDE CAVE1 BIRD  CAVE2 RNBOW
!             SNAKE STYX

WRITE (*,2680) hint
2680 FORMAT (' TRYING TO PRINT HINT # ',i1/)
CALL bug (27)

2690 hintlc(hint)=0
IF (.NOT.yes(hints(hint,3),0,54)) GO TO 500
WRITE (*,2700) hints(hint,2)
2700 FORMAT (/' I am prepared to give you a hint, but it will cost you', i2, &
             ' points.')
hinted(hint)=yes(175,hints(hint,4),54)
IF (hinted(hint) .AND. limit > 30) limit=limit + 30*hints(hint,2)
2710 hintlc(hint)=0
GO TO 500

!  NOW FOR THE QUICK TESTS.  SEE DATABASE DESCRIPTION FOR ONE-LINE NOTES

2720 IF (atloc(loc) == 0 .AND. atloc(oldloc) == 0 .AND. atloc(oldlc2) == 0 &
         .AND. burden(0) > 1) GO TO 2690
GO TO 2710

2730 IF (prop(emrald) /= -1 .AND. prop(pyram) == -1) GO TO 2690
GO TO 2710


2740 IF ((prop(sword) == 1 .OR. prop(sword) == 5) .AND. .NOT.toting(crown))  &
    GO TO 2690
GO TO 2710


2750 IF (prop(grate) == 0 .AND. .NOT.athand(keys)) GO TO 2690
GO TO 2710

2760 IF (here(bird) .AND. athand(rod) .AND. obj == bird) GO TO 2690
GO TO 500

2770 IF (abb(159) == 0) GO TO 2690
GO TO 2710

2780 IF (.NOT.toting(shoes) .OR. abb(205) == 0) GO TO 2690
GO TO 2710

2790 IF (.NOT.athand(lyre) .AND. prop(dog) /= 1) GO TO 2690
GO TO 2710

2800 IF (here(snake) .AND. .NOT.here(bird)) GO TO 2690
GO TO 2710

!  CAVE CLOSING AND SCORING

!  THESE SECTIONS HANDLE THE CLOSING OF THE CAVE.  THE CAVE CLOSES "CLOCK1"
!  TURNS AFTER THE LAST TREASURE HAS BEEN LOCATED (INCLUDING THE PIRATE'S
!  CHEST, WHICH MAY OF COURSE NEVER SHOW UP).  NOTE THAT THE TREASURES NEED NOT
!  HAVE BEEN TAKEN YET, JUST LOCATED.  HENCE CLOCK1 MUST BE LARGE ENOUGH TO GET
!  OUT OF THE CAVE (IT ONLY TICKS WHILE INSIDE THE CAVE).  WHEN IT HITS ZERO,
!  WE BRANCH TO 90000 TO START CLOSING THE CAVE, AND THEN SIT BACK AND WAIT FOR
!  HIM TO TRY TO GET OUT.  IF HE DOESN'T WITHIN CLOCK2 TURNS, WE CLOSE THE
!  CAVE; IF HE DOES TRY, WE ASSUME HE PANICS, AND GIVE HIM A FEW ADDITIONAL
!  TURNS TO GET FRANTIC BEFORE WE CLOSE.  WHEN CLOCK2 HITS ZERO, WE BRANCH TO
!  90000 TO TRANSPORT HIM INTO THE FINAL PUZZLE.  NOTE THAT THE PUZZLE DEPENDS
!  UPON ALL SORTS OF RANDOM THINGS.  FOR INSTANCE, THERE MUST BE NO WATER OR
!  OIL, SINCE THERE ARE BEANSTALKS WHICH WE DON'T WANT TO BE ABLE TO WATER,
!  SINCE THE CODE CAN'T HANDLE IT.  ALSO, WE CAN HAVE NO KEYS, SINCE THERE IS A
!  GRATE (HAVING MOVED THE FIXED OBJECT!) THERE SEPARATING HIM FROM ALL THE
!  TREASURES.  MOST OF THESE PROBLEMS ARISE FROM THE USE OF NEGATIVE PROP
!  NUMBERS TO SUPPRESS THE OBJECT DESCRIPTIONS UNTIL HE'S ACTUALLY MOVED THE
!  OBJECTS.

!  WHEN THE FIRST WARNING COMES, WE LOCK THE GRATE, DESTROY THE BRIDGE, KILL
!  ALL THE DWARVES (AND THE PIRATE), REMOVE THE TROLL AND BEAR (UNLESS DEAD),
!  AND SET "CLOSNG" TO TRUE.  LEAVE THE DRAGON; TOO MUCH TROUBLE TO MOVE IT.
!  FROM NOW UNTIL CLOCK2 RUNS OUT, HE CANNOT UNLOCK THE GRATE, MOVE TO ANY
!  LOCATION OUTSIDE THE CAVE (LOC=BITSET(LOCCON,6)), OR CREATE THE BRIDGE.  NOR
!  RESURRECTED IF HE DIES.  NOTE THAT THE SNAKE IS ALREADY GONE, SINCE HE GOT
!  TO THE TREASURE ACCESSIBLE ONLY VIA THE HALL OF THE MT. KING.  ALSO, HE'S
!  BEEN IN GIANT ROOM (TO GET EGGS), SO WE CAN REFER TO IT.  ALSO ALSO, HE'S
!  GOTTEN THE PEARL, SO WE KNOW THE BIVALVE IS AN OYSTER.  *AND*, THE DWARVES
!  MUST HAVE BEEN ACTIVATED, SINCE WE'VE FOUND CHEST.

2810 prop(grate)=0
prop(fissur)=0
prop(tdoor)=0
prop(tdoor2)=0
DO  i=1,6
  dseen(i)=.false.
  dloc(i)=0
END DO
CALL move (troll, 0)
CALL move (troll+maxobj, 0)
CALL move (troll2, plac(troll))
CALL move (troll2+maxobj, fixd(troll))
CALL juggle (chasm)
IF (prop(bear) /= 3) CALL dstroy (bear)
prop(chain)=0
fixed(chain)=0
prop(axe)=0
fixed(axe)=0
CALL rspeak (129)
clock1=-1
closng=.true.
GO TO 550

!  ONCE HE'S PANICKED, AND CLOCK2 HAS RUN OUT, WE COME HERE TO SET UP THE
!  STORAGE ROOM.  THE ROOM HAS TWO LOCS, HARDWIRED AS 115 (NE) AND 116 (SW).
!  AT THE NE END, WE PLACE EMPTY BOTTLES, A NURSERY OF PLANTS, A BED OF
!  OYSTERS, A PILE OF LAMPS, RODS WITH STARS, SLEEPING DWARVES, PHONE BOOTH AND
!  AT THE SW END WE PLACE GRATE OVER TREASURES, SNAKE PIT, COVEY OF CAGED BIRDS,
!  MORE RODS, AND PILLOWS.  A MIRROR STRETCHES ACROSS ONE WALL.  MANY OF THE
!  OBJECTS COME FROM KNOWN LOCATIONS AND/OR STATES (E.G. THE SNAKE IS KNOWN TO
!  HAVE BEEN DESTROYED AND NEEDN'T BE CARRIED AWAY FROM ITS OLD "PLACE"),
!  MAKING THE VARIOUS OBJECTS BE HANDLED DIFFERENTLY.  WE ALSO DROP ALL OTHER
!  OBJECTS HE MIGHT BE CARRYING (LEST HE HAVE SOME WHICH COULD CAUSE TROUBLE,
!  SUCH AS THE KEYS).  WE DESCRIBE THE FLASH OF LIGHT AND TRUNDLE BACK.
!  THE PHONE MAKES IT IMPOSSIBLE FOR THE WALDO TO FART AROUND IN THE
!  REPOSITORY TOO LONG.  WHEN CLOCK3 TICKS TO ZERO, THE PHONE STARTS
!  RINGING.  WHEN IT HITS -7, THE DWARVES WAKE UP.  IF HE FIDDLES WITH
!  THE PHONE, HE GETS ZONKED IN OTHER WAYS.

2830 prop(bottle)=put(bottle,115,1)
prop(plant)=put(plant,115,0)
prop(oyster)=put(oyster,115,0)
prop(lamp)=put(lamp,115,0)
prop(rod)=put(rod,115,0)
prop(dwarf)=put(dwarf,115,0)
loc=115
oldloc=115
newloc=115

!  LEAVE THE GRATE WITH NORMAL (NON-NEGATIVE PROPERTY).

foo=put(grate,116,0)
prop(snake)=put(snake,116,1)
prop(bird)=put(bird,116,1)
prop(cage)=put(cage,116,0)
prop(rod2)=put(rod2,116,0)
prop(pillow)=put(pillow,116,0)

prop(booth)=put(booth,116,-4)
fixed(booth)=115
prop(phone)=put(phone,212,-4)

prop(mirror)=put(mirror,115,0)
prop(book2)=put(book2,115,0)
fixed(mirror)=116

DO  i=1,maxobj
  IF (toting(i) .AND. enclsd(i)) CALL remove (i)
  IF (toting(i)) CALL dstroy (i)
END DO

CALL rspeak (132)
closed=.true.
GO TO 10

!  ANOTHER WAY WE CAN FORCE AN END TO THINGS IS BY HAVING THE LAMP GIVE OUT.
!  WHEN IT GETS CLOSE, WE COME HERE TO WARN HIM.
!  92000 IS FOR CASES OF LAMP DYING.  92400 IS WHEN IT GOES OUT,
!  AND 92600 IS IF HE'S WANDERED OUTSIDE AND THE LAMP IS USED UP, IN WHICH
!  CASE WE FORCE HIM TO GIVE UP.

2850 IF (lmwarn .OR. .NOT.here(lamp)) GO TO 550
lmwarn=.true.
spk=187
IF (prop(batter) == 1) spk=323
IF (place(batter) == 0) spk=183
IF (prop(vend) == 1) spk=189
CALL rspeak (spk)
GO TO 550

2860 limit=-1
prop(lamp)=0
IF (here(lamp)) CALL rspeak (184)
GO TO 550

2870 CALL rspeak (185)
gaveup=.true.
GO TO 2970


!  OH DEAR, HE'S DISTURBED THE DWARVES.

2880 CALL rspeak (136)
GO TO 2970

!  HIT THE PHONE.  JINGLE, JANGLE, CRASH!
2890 CALL rspeak (282)
GO TO 2880

!  WHOOPS.  DOWN THE CHUTE.
2900 CALL rspeak (283)
GO TO 2970

!  DWARVES CAN'T SLEEP THRU ALL THIS RINGING!
2910 CALL rspeak (254)
GO TO 2970

!  "YOU'RE DEAD, FRED."         "YOU DIED, CLYDE?"
!
!  IF THE CURRENT LOC IS ZERO, IT MEANS THE CLOWN GOT HIMSELF KILLED.  WE'LL
!  ALLOW THIS MAXDIE TIMES.  MAXDIE IS AUTOMATICALLY SET BASED ON THE NUMBER OF
!  SNIDE MESSAGES AVAILABLE.  EACH DEATH RESULTS IN A MESSAGE (81, 83, ETC.)
!  WHICH OFFERS REINCARNATION; IF ACCEPTED, THIS RESULTS IN MESSAGE 82, 84,
!  ETC.  THE LAST TIME, IF HE WANTS ANOTHER CHANCE, HE GETS A SNIDE REMARK AS
!  WE EXIT.  WHEN REINCARNATED, ALL OBJECTS BEING CARRIED GET DROPPED AT OLDLC2
!  (PRESUMABLY THE LAST PLACE PRIOR TO BEING KILLED) WITHOUT CHANGE OF PROPS.
!  THE LOOP RUNS BACKWARDS TO ASSURE THAT THE BIRD IS DROPPED BEFORE THE CAGE.
!  (THIS KLUGE COULD BE CHANGED ONCE WE'RE SURE ALL REFERENCES TO BIRD AND CAGE
!  ARE DONE BY KEYWORDS.)  THE LAMP IS A SPECIAL CASE (IT WOULDN'T DO TO LEAVE
!  IT IN THE CAVE).  IT IS TURNED OFF AND LEFT OUTSIDE THE BUILDING (ONLY IF HE
!  WAS CARRYING IT, OF COURSE).  HE HIMSELF IS LEFT INSIDE THE BUILDING (AND
!  HEAVEN HELP HIM IF HE TRIES TO XYZZY BACK INTO THE CAVE WITHOUT THE LAMP!).
!  OLDLOC IS ZAPPED SO HE CAN'T JUST "RETREAT".

!  THE EASIEST WAY TO GET KILLED IS TO FALL INTO A PIT IN PITCH DARKNESS.

2920 CALL rspeak (23)
oldlc2=loc

!  OKAY, HE'S DEAD.  LET'S GET ON WITH IT.

2930 IF (closng) GO TO 2960
yea=yes(81+numdie*2,82+numdie*2,54)
numdie=numdie+1
IF (numdie == maxdie .OR. .NOT.yea) GO TO 2970
IF (chase == 0) GO TO 2940

!  CHAMP WAS BEING CHASED BY WUMPUS, & DIED ANOTHER WAY.
!  PUT WUMPUS BACK TO SLEEP, IN CASE OUR HERO STILL HASN'T GOT THE CLOAK.
chase=0
prop(wumpus)=0
CALL move (wumpus, 174)
2940 IF (toting(lamp)) prop(lamp)=0
DO  j=1,maxobj
  IF (.NOT.holdng(j)) CYCLE
  loc=oldlc2
  IF (j == lamp) loc=1
  CALL drop (j, loc)
  IF (.NOT.wearng(j)) CYCLE
  prop(j)=0
  CALL bitoff (j, wearbt)
END DO
loc=3
oldloc=loc
wdx=0
words(1)=0
CALL clrlin()
health=100
GO TO 280

!  HE DIED DURING CLOSING TIME.  NO RESURRECTION.  TALLY UP A DEATH AND EXIT.

2960 CALL rspeak (131)
numdie=numdie+1
GO TO 2970

!  IT'S OVER.  TALLY SCORE.

2970 CALL rating (score, bonus, gaveup, scorng, closng, closed, hntmax)

!  THAT SHOULD BE GOOD ENOUGH.  LET'S TELL HIM ALL ABOUT IT.

WRITE (*,2980) score,mxscor,turns
2980 FORMAT (///' You scored',i4,' out of a possible',i4,', using',i5,' turns.')

DO  i=1,clsses
  IF (cval(i) >= score) GO TO 3010
END DO
WRITE (*,3000)
3000 FORMAT (/' You just went off my scale!!'/)
GO TO 3050

3010 CALL speak (ctext(i))
IF (i == clsses-1) GO TO 3030
i = MAX(i,1)
k=cval(i)+1-score
kk2c='s.'
IF (k == 1) kk2c='. '
WRITE (*,3020) k,kk2c
3020 FORMAT (/' To achieve the next higher rating, you need', i3,  &
             ' more point', a2/)
GO TO 3050

3030 WRITE (*,3040)
3040 FORMAT (/' To achieve the next higher rating would be a neat trick!'//  &
             ' CONGRATULATIONS!!'/)

3050 STOP

END PROGRAM Adventure

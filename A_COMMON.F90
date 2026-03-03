MODULE Adventure_Common
! The COMMON blocks are all grouped together here.

IMPLICIT NONE

! COMMON /izwiz/ iswiz
LOGICAL, SAVE  :: iswiz

! COMMON /adjcom/ adjkey(50),adjtab(150),adjsiz
INTEGER, SAVE  :: adjkey(50),adjtab(150),adjsiz

! COMMON /bitcom/ openbt,lockbt,burnbt,wearbt
INTEGER, SAVE  :: openbt,lockbt,burnbt,wearbt

! COMMON /blkcom/ blklin
LOGICAL, SAVE  :: blklin

! COMMON /concom/ loccon(250),objcon(150)
INTEGER, SAVE  :: loccon(250),objcon(150)

! COMMON /diecom/ numdie,maxdie,turns,killed
INTEGER, SAVE  :: numdie,maxdie,turns
LOGICAL, SAVE  :: killed

! COMMON /dwfcom/ dwarf,knife,knfloc,dflag,dseen(6),dloc(6),odloc(6),dwfmax
INTEGER, SAVE  :: dwarf,knife,knfloc,dflag,dloc(6),odloc(6),dwfmax
LOGICAL, SAVE  :: dseen(6)

! COMMON /hldcom/ holder(150),hlink(150)
INTEGER, SAVE  :: holder(150),hlink(150)

! COMMON /hntcom/ hintlc(20),hinted(20),hints(20,4),hntsiz,hntmin
INTEGER, SAVE  :: hintlc(20),hints(20,4),hntsiz,hntmin
LOGICAL, SAVE  :: hinted(20)

! COMMON /liqcom/ bottle,cask,water,oil,wine,liqtyp(5)
INTEGER, SAVE  :: bottle,cask,water,oil,wine,liqtyp(5)

! COMMON /loccom/ loc,oldloc,oldlc2,newloc,maxloc
INTEGER, SAVE  :: loc,oldloc,oldlc2,newloc,maxloc

! COMMON /mnecom/ back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
!    bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
!    lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
!    troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
!    iwest,phuce(2,4),tk(20)
INTEGER, SAVE  :: back,cave,dprssn,entrnc,EXIT,GO,look,null,axe,  &
     bear,boat,book,book2,booth,carvng,chasm,chasm2,door,gnome,grate,  &
     lamp,pdoor,plant,plant2,rocks,rod,rod2,safe,tdoor,tdoor2,troll,  &
     troll2,emrald,spices,find,yell,invent,leave,pour,say,take,throw,  &
     iwest,phuce(2,4),tk(20)

! COMMON /ltxcom/ ltext(250),stext(250),key(250),abb(250),locsiz
INTEGER, SAVE  :: ltext(250),stext(250),key(250),abb(250),locsiz

! COMMON /objcom/ plac(150),fixd(150),weight(150),prop(150),points(150)
INTEGER, SAVE  :: plac(150),fixd(150),weight(150),prop(150),points(150)

! COMMON /placom/ atloc(250),link(300),place(150),fixed(150),maxobj
INTEGER, SAVE  :: atloc(250),link(300),place(150),fixed(150),maxobj

! COMMON /prpcom/ vkey(60),ptab(300),vkysiz,ptbsiz
INTEGER, SAVE  :: vkey(60),ptab(300),vkysiz,ptbsiz

! COMMON /trvcom/ travel(1600)
INTEGER, SAVE  :: travel(1600)

! COMMON /txtcom/ lines(25000),rtext(450),ptext(150)
INTEGER, SAVE  :: lines(25000),rtext(450),ptext(150)

! COMMON /utxcom/ wdx
INTEGER, SAVE  :: wdx

! COMMON /voccom/ ktab(600),tabsiz
INTEGER, SAVE  :: ktab(600),tabsiz

! COMMON /wrdcom/ verbs(45),vrbx,objs(45),objx,iobjs(15),iobx,prep,words(45)
INTEGER, SAVE  :: verbs(45),vrbx,objs(45),objx,iobjs(15),iobx,prep,words(45)


! COMMON /savcom/ abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,  &
!    clock2,clock3,closed,closng,clsmax,combo,deadbt,detail,dkill,  &
!    dtotal,dwarfn,flg239,foo,foobar,food,gaveup,health,hint,hit,  &
!    hntmax,i,ikk,iloc,iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,  &
!    lmwarn,lock,logout,messag,obj,panic,portal,ptbs,rdflag,retn,  &
!    rtxsiz,score,scorng,sect,skey,sloc,spk,start,stick,tabndx,tally,  &
!    tally2,terse,trvs,trvsiz,vend,verb,vrbsiz,waste,wkday,wkend,  &
!    wzdark,yea,actspk,ctext,cval,hname(10)
INTEGER, SAVE  :: abbnum,adj,atbs,attack,bcross,bonus,chase,clock1,clock2, &
                  clock3,clsmax,combo,deadbt,detail,dkill,dtotal,dwarfn,  &
                  flg239,foo,foobar,food,health,hint,hit,hntmax,i,ikk,iloc, &
                  iobj,j,jj,jk1,jkk,k,k1,kk,l,l1,limit,linsiz,ll,lock,   &
                  messag,obj,portal,ptbs,retn,rtxsiz,score,sect,skey,sloc,  &
                  spk,stick,tabndx,tally,tally2,trvs,trvsiz,vend,verb,vrbsiz,  &
                  waste,wkday,wkend,actspk(60),ctext(12),cval(12),hname(10)
LOGICAL, SAVE  :: closed,closng,gaveup,lmwarn,logout,panic,rdflag,scorng,  &
                  start,terse,wzdark,yea

! COMMON /sv2com/ anvil,batter,bees,billbd,bird,brush,cage,cakes,  &
!    chain,chest,chloc,chloc2,clam,cloak,clsses,coins,crown,daltlc,dog,  &
!    dragon,eggs,fissur,flower,gatloc,grail,hive,honey,horn,jewels,  &
!    keys,lyre,magzin,mirror,mushrm,mxscor,nugget,oyster,pearl,phone,  &
!    pillow,pole,poster,prepat,prepdn,prepfr,prepin,prepof,prepon,  &
!    pyram,radium,ring,rug,sapphi,shield,shoes,shut,slugs,snake,sphere,  &
!    steps,sticks,sword,tablet,tridnt,unlock,vase,wall,wall2,wear,wumpus,y2,yank
INTEGER, SAVE  :: anvil,batter,bees,billbd,bird,brush,cage,cakes,  &
    chain,chest,chloc,chloc2,clam,cloak,clsses,coins,crown,daltlc,dog,  &
    dragon,eggs,fissur,flower,gatloc,grail,hive,honey,horn,jewels,  &
    keys,lyre,magzin,mirror,mushrm,mxscor,nugget,oyster,pearl,phone,  &
    pillow,pole,poster,prepat,prepdn,prepfr,prepin,prepof,prepon,  &
    pyram,radium,ring,rug,sapphi,shield,shoes,shut,slugs,snake,sphere,  &
    steps,sticks,sword,tablet,tridnt,unlock,vase,wall,wall2,wear,wumpus,  &
    y2,yank

! COMMON /sv3com/ dtk,atab(600),vtxt(45,2),otxt(45,2),iotxt(15,2),txt(35,2)
CHARACTER (LEN=6), SAVE  :: dtk(9),atab(600),vtxt(45,2),otxt(45,2),  &
                            iotxt(15,2),txt(35,2)

! COMMON /tnoux/ indent
INTEGER, SAVE  :: indent

END MODULE Adventure_Common

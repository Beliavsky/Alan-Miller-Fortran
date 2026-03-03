MODULE Genetic_Algorithm
IMPLICIT NONE

!     Common block to make iseed visible to rninit (and to save
!     it between calls)
! COMMON /rnseed/ iseed
INTEGER, SAVE  :: iseed


CONTAINS


SUBROUTINE pikaia(ff,n,ctrl,x,f,STATUS)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-09  Time: 15:54:13
 
!====================================================================
!  Optimization (maximization) of user-supplied "fitness" function ff
!  over n-dimensional parameter space  x  using a basic genetic algorithm
!  method.

!  Paul Charbonneau & Barry Knapp
!  High Altitude Observatory
!  National Center for Atmospheric Research
!  Boulder CO 80307-3000
!  USA
!  <paulchar@hao.ucar.edu>
!  <knapp@hao.ucar.edu>

!  Web site:
!  http://www.hao.ucar.edu/public/research/si/pikaia/pikaia.html

!  Version 1.0   [ 1995 December 01 ]

!  Genetic algorithms are heuristic search techniques that incorporate in a
!  computational setting, the biological notion of evolution by means of
!  natural selection.  This subroutine implements the three basic operations
!  of selection, crossover, and mutation, operating on "genotypes" encoded as
!  strings.

!  References:

!     Charbonneau, Paul.  "Genetic Algorithms in Astronomy and Astrophysics."
!        Astrophysical J. (Supplement), vol 101, in press (December 1995).

!     Goldberg, David E.  Genetic Algorithms in Search, Optimization,
!        & Machine Learning.  Addison-Wesley, 1989.

!     Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
!        Van Nostrand Reinhold, 1991.
!====================================================================
!  USES: ff, urand, setctl, report, rnkpop, select, encode, decode,
!        cross, mutate, genrep, stdrep, newpop, adjmut

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: ctrl(12)
REAL, INTENT(OUT)     :: x(n)
REAL, INTENT(OUT)     :: f
INTEGER, INTENT(OUT)  :: STATUS

INTERFACE
  FUNCTION ff(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION ff
END INTERFACE

! EXTERNAL ff

!  Input:
!  o Integer  n  is the parameter space dimension, i.e., the number
!    of adjustable parameters.

!  o Function  ff  is a user-supplied scalar function of n variables, which
!    must have the calling sequence f = ff(n,x), where x is a real parameter
!    array of length n.  This function must be written so as to bound all
!    parameters to the interval [0,1]; that is, the user must determine
!    a priori bounds for the parameter space, and ff must use these bounds
!    to perform the appropriate scalings to recover true parameter values in
!    the a priori ranges.

!    By convention, ff should return higher values for more optimal
!    parameter values (i.e., individuals which are more "fit").
!    For example, in fitting a function through data points, ff
!    could return the inverse of chi**2.

!    In most cases initialization code will have to be written
!    (either in a driver or in a separate subroutine) which loads
!    in data values and communicates with ff via one or more labeled
!    common blocks.  An example exercise driver and fitness function
!    are provided in the accompanying file, xpkaia.f.


!  Input/Output:


!  o Array  ctrl  is an array of control flags and parameters, to
!    control the genetic behavior of the algorithm, and also printed
!    output.  A default value will be used for any control variable
!    which is supplied with a value less than zero.  On exit, ctrl
!    contains the actual values used as control variables.  The
!    elements of ctrl and their defaults are:

!       ctrl( 1) - number of individuals in a population (default
!                  is 100)
!       ctrl( 2) - number of generations over which solution is
!                  to evolve (default is 500)
!       ctrl( 3) - number of significant digits (i.e., number of
!                  genes) retained in chromosomal encoding (default
!                  is 6)  (Note: This number is limited by the
!                  machine floating point precision.  Most 32-bit
!                  floating point representations have only 6 full
!                  digits of precision.  To achieve greater preci-
!                  sion this routine could be converted to double
!                  precision, but note that this would also require
!                  a double precision random number generator, which
!                  likely would not have more than 9 digits of
!                  precision if it used 4-byte integers internally.)
!       ctrl( 4) - crossover probability; must be  <= 1.0 (default is 0.85)
!       ctrl( 5) - mutation mode; 1/2=steady/variable (default is 2)
!       ctrl( 6) - initial mutation rate; should be small (default is 0.005)
!                  (Note: the mutation rate is the probability that any one
!                  gene locus will mutate in any one generation.)
!       ctrl( 7) - minimum mutation rate; must be >= 0.0 (default is 0.0005)
!       ctrl( 8) - maximum mutation rate; must be <= 1.0 (default is 0.25)
!       ctrl( 9) - relative fitness differential; range from 0
!                  (none) to 1 (maximum).  (default is 1.)
!       ctrl(10) - reproduction plan; 1/2/3=Full generational
!                  replacement/Steady-state-replace-random/Steady-
!                  state-replace-worst (default is 3)
!       ctrl(11) - elitism flag; 0/1=off/on (default is 0)
!                  (Applies only to reproduction plans 1 and 2)
!       ctrl(12) - printed output 0/1/2=None/Minimal/Verbose (default is 0)


! Output:


!  o Array  x(1:n)  is the "fittest" (optimal) solution found,
!     i.e., the solution which maximizes fitness function ff

!  o Scalar  f  is the value of the fitness function at x

!  o Integer  status  is an indicator of the success or failure
!     of the call to pikaia (0=success; non-zero=failure)


! Constants

INTEGER, PARAMETER :: nmax = 32, pmax = 128, dmax = 6

!  o NMAX is the maximum number of adjustable parameters (n <= NMAX)

!  o PMAX is the maximum population (ctrl(1) <= PMAX)

!  o DMAX is the maximum number of Genes (digits) per Chromosome
!        segement (parameter) (ctrl(3) <= DMAX)


!     Local variables
INTEGER :: np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig, ip1,  &
           ip2, NEW, newtot
REAL :: pcross, pmut, pmutmn, pmutmx, fdif

REAL :: ph(nmax,2), oldph(nmax,pmax), newph(nmax,pmax)

INTEGER :: gn1(nmax*dmax), gn2(nmax*dmax)
INTEGER :: ifit(pmax), jfit(pmax)
REAL :: fitns(pmax)

!     User-supplied uniform random number generator
! REAL :: urand
! EXTERNAL urand

! Function urand should not take any arguments.  If the user wishes to be able
! to initialize urand, so that the same sequence of random numbers can be
! repeated, this capability could be implemented with a separate subroutine,
! and called from the user's driver program.  An example urand function
! (and initialization subroutine) which uses the function ran0 (the "minimal
! standard" random number generator of Park and Miller [Comm. ACM 31, 1192-
! 1201, Oct 1988; Comm. ACM 36 No. 7, 105-110, July 1993]) is provided.


!     Set control variables from input and defaults
CALL setctl(ctrl, n, np, ngen, nd, pcross, pmutmn, pmutmx, pmut, imut, fdif, &
            irep, ielite, ivrb, STATUS)
IF (STATUS /= 0) THEN
  WRITE (*, *) ' Control vector (ctrl) argument(s) invalid'
  RETURN
END IF

!     Make sure locally-dimensioned arrays are big enough
IF (n > nmax .OR. np > pmax .OR. nd > dmax) THEN
  WRITE (*, *) ' Number of parameters, population, or genes too large'
  STATUS = -1
  RETURN
END IF

!     Compute initial (random but bounded) phenotypes
DO  ip = 1, np
  DO  k = 1, n
    oldph(k,ip) = urand()
  END DO
  fitns(ip) = ff(n, oldph(:,ip))
END DO

!     Rank initial population by fitness order
CALL rnkpop(np,fitns,ifit,jfit)

!     Main Generation Loop
DO  ig = 1, ngen
  
!        Main Population Loop
  newtot = 0
  DO  ip = 1, np / 2
    
!           1. pick two parents
    CALL select(np,jfit,fdif,ip1)
    30 CALL select(np,jfit,fdif,ip2)
    IF (ip1 == ip2) GO TO 30
    
!           2. encode parent phenotypes
    CALL encode(n,nd,oldph(1,ip1),gn1)
    CALL encode(n,nd,oldph(1,ip2),gn2)
    
!           3. breed
    CALL cross(n,nd,pcross,gn1,gn2)
    CALL mutate(n,nd,pmut,gn1)
    CALL mutate(n,nd,pmut,gn2)
    
!           4. decode offspring genotypes
    CALL decode(n,nd,gn1,ph(1,1))
    CALL decode(n,nd,gn2,ph(1,2))
    
!           5. insert into population
    IF (irep == 1) THEN
      CALL genrep(nmax,n,np,ip,ph,newph)
    ELSE
      CALL stdrep(ff,nmax,n,np,irep,ielite,ph,oldph,fitns,ifit, jfit,NEW)
      newtot = newtot + NEW
    END IF
    
!        End of Main Population Loop
  END DO
  
!        if running full generational replacement: swap populations
  IF (irep == 1) CALL newpop(ff,ielite,nmax,n,np,oldph,newph,ifit,  &
      jfit,fitns,newtot)
  
!        adjust mutation rate?
  IF (imut == 2) CALL adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
  
!        print generation report to standard output?
  IF (ivrb > 0) CALL report(ivrb,nmax,n,np,nd,oldph,fitns,ifit,pmut,ig,newtot)
  
!     End of Main Generation Loop
END DO

!     Return best phenotype and its fitness
DO  k = 1, n
  x(k) = oldph(k,ifit(np))
END DO
f = fitns(ifit(np))

RETURN
END SUBROUTINE pikaia

!********************************************************************

SUBROUTINE setctl(ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,  &
                  imut,fdif,irep,ielite,ivrb,STATUS)
!===================================================================
!     Set control variables and flags from input and defaults
!===================================================================

!     Input
!     Input/Output
REAL, INTENT(IN OUT)  :: ctrl(12)
INTEGER, INTENT(IN)   :: n

!     Output
INTEGER, INTENT(OUT)  :: np
INTEGER, INTENT(OUT)  :: ngen
INTEGER, INTENT(OUT)  :: nd
REAL, INTENT(OUT)     :: pcross
REAL, INTENT(OUT)     :: pmutmn
REAL, INTENT(OUT)     :: pmutmx
REAL, INTENT(OUT)     :: pmut
INTEGER, INTENT(OUT)  :: imut
REAL, INTENT(OUT)     :: fdif
INTEGER, INTENT(OUT)  :: irep
INTEGER, INTENT(OUT)  :: ielite
INTEGER, INTENT(OUT)  :: ivrb
INTEGER, INTENT(OUT)  :: STATUS


!     Local
INTEGER :: i
REAL, SAVE  :: dfault(12) = (/ 100., 500., 5., .85, 2., .005, .0005, .25,  &
                               1., 1., 1., 0. /)

DO  i = 1, 12
  IF (ctrl(i) < 0.) ctrl(i) = dfault(i)
END DO

np = ctrl(1)
ngen = ctrl(2)
nd = ctrl(3)
pcross = ctrl(4)
imut = ctrl(5)
pmut = ctrl(6)
pmutmn = ctrl(7)
pmutmx = ctrl(8)
fdif = ctrl(9)
irep = ctrl(10)
ielite = ctrl(11)
ivrb = ctrl(12)
STATUS = 0

!     Print a header
IF (ivrb > 0) THEN
  
  WRITE (*,5000) ngen, np, n, nd, pcross, pmut, pmutmn, pmutmx, fdif
  IF (imut == 1) WRITE (*,5100) 'Constant'
  IF (imut == 2) WRITE (*,5100) 'Variable'
  IF (irep == 1) WRITE (*,5200) 'Full generational replacement'
  IF (irep == 2) WRITE (*,5200) 'Steady-state-replace-random'
  IF (irep == 3) WRITE (*,5200) 'Steady-state-replace-worst'
END IF

!     Check some control values
IF (imut /= 1 .AND. imut /= 2) THEN
  WRITE (*,5300)
  STATUS = 5
END IF

IF (fdif > 1.) THEN
  WRITE (*,5400)
  STATUS = 9
END IF

IF (irep /= 1 .AND. irep /= 2 .AND. irep /= 3) THEN
  WRITE (*,5500)
  STATUS = 10
END IF

IF (pcross > 1.0 .OR. pcross < 0.) THEN
  WRITE (*,5600)
  STATUS = 4
END IF

IF (ielite /= 0 .AND. ielite /= 1) THEN
  WRITE (*,5700)
  STATUS = 11
END IF

IF (irep == 1 .AND. imut == 1 .AND. pmut > 0.5 .AND. ielite == 0) THEN
  WRITE (*,5800)
END IF

IF (irep == 1 .AND. imut == 2 .AND. pmutmx > 0.5 .AND. ielite == 0) THEN
  WRITE (*,5900)
END IF

IF (fdif < 0.33 .AND. irep /= 3) THEN
  WRITE (*,6000)
END IF

IF (MOD(np,2) > 0) THEN
  np = np - 1
  WRITE (*,6100) np
END IF

RETURN
5000 FORMAT (/' ', 60('*') /   &
    ' *', t16, 'PIKAIA Genetic Algorithm Report ', t60, '*' / &
    ' ', 60('*') //  &
    '   Number of Generations evolving: ', i4 /  &
    '       Individuals per generation: ', i4 /  &
    '    Number of Chromosome segments: ', i4 /  &
    '    Length of Chromosome segments: ', i4 /  &
    '            Crossover probability: ', f9.4 /   &
    '            Initial mutation rate: ', f9.4 /   &
    '            Minimum mutation rate: ', f9.4 /   &
    '            Maximum mutation rate: ', f9.4 /   &
    '    Relative fitness differential: ', f9.4)
5100 FORMAT ('                    Mutation Mode: '/ a)
5200 FORMAT ('                Reproduction Plan: '/ a)
5300 FORMAT (' ERROR: illegal value for imut (ctrl(5))')
5400 FORMAT (' ERROR: illegal value for fdif (ctrl(9))')
5500 FORMAT (' ERROR: illegal value for irep (ctrl(10))')
5600 FORMAT (' ERROR: illegal value for pcross (ctrl(4))')
5700 FORMAT (' ERROR: illegal value for ielite (ctrl(11))')
5800 FORMAT (' WARNING: dangerously high value for pmut (ctrl(6));' /  &
    ' (Should enforce elitism with ctrl(11)=1.)')
5900 FORMAT (' WARNING: dangerously high value for pmutmx (ctrl(8));' /  &
    ' (Should enforce elitism with ctrl(11)=1.)')
6000 FORMAT (' WARNING: dangerously low value of fdif (ctrl(9))')
6100 FORMAT (' WARNING: decreasing population size (ctrl(1)) to np='/ i4 )
END SUBROUTINE setctl

!********************************************************************

SUBROUTINE report(ivrb, ndim, n, np, nd, oldph, fitns, ifit, pmut, ig, nnew)

!     Write generation report to standard output

!     Input:
INTEGER,  INTENT(IN)  :: ivrb
INTEGER,  INTENT(IN)  :: ndim
INTEGER,  INTENT(IN)  :: n
INTEGER,  INTENT(IN)  :: np
INTEGER,  INTENT(IN)  :: nd
REAL, INTENT(IN)      :: oldph(ndim, np)
REAL, INTENT(IN)      :: fitns(np)
INTEGER, INTENT(IN)   :: ifit(np)
REAL, INTENT(IN)      :: pmut
INTEGER, INTENT(IN)   :: ig
INTEGER, INTENT(IN)   :: nnew

!     Output: none

!     Local
REAL, SAVE  :: bestft = 0.0, pmutpv = 0.0
INTEGER  :: ndpwr, k
LOGICAL  :: rpt

rpt = .false.

IF (pmut /= pmutpv) THEN
  pmutpv = pmut
  rpt = .true.
END IF

IF (fitns(ifit(np)) /= bestft) THEN
  bestft = fitns(ifit(np))
  rpt = .true.
END IF

IF (rpt .OR. ivrb >= 2) THEN
  
!        Power of 10 to make integer genotypes for display
  ndpwr = nint(10.**nd)
  
  WRITE (*, '(/i6, i6, f10.6, 4f10.6)') ig, nnew, pmut,  &
      fitns(ifit(np)), fitns(ifit(np-1)), fitns(ifit(np/2))
  DO  k = 1, n
    WRITE (*, '(22x, 3i10)') nint(ndpwr*oldph(k, ifit(np))),  &
        nint(ndpwr*oldph(k, ifit(np-1))), nint(ndpwr*oldph(k, ifit(np/2)))
  END DO
  
END IF
RETURN
END SUBROUTINE report

!**********************************************************************
!                         GENETICS MODULE
!**********************************************************************

!     ENCODE:    encodes phenotype into genotype
!                called by: PIKAIA

!     DECODE:    decodes genotype into phenotype
!                called by: PIKAIA

!     CROSS:     Breeds two offspring from two parents
!                called by: PIKAIA

!     MUTATE:    Introduces random mutation in a genotype
!                called by: PIKAIA

!     ADJMUT:    Implements variable mutation rate
!                called by: PIKAIA

!**********************************************************************

SUBROUTINE encode(n, nd, ph, gn)
!======================================================================
!     encode phenotype parameters into integer genotype
!     ph(k) are x, y coordinates [ 0 < x, y < 1 ]
!======================================================================


INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: nd
REAL, INTENT(IN OUT)  :: ph(n)
INTEGER, INTENT(OUT)  :: gn(n*nd)

!     Inputs:



!     Output:


!     Local:
INTEGER :: ip, i, j, ii
REAL :: z

z = 10. ** nd
ii = 0
DO  i = 1, n
  ip = INT(ph(i)*z)
  DO  j = nd, 1, -1
    gn(ii+j) = MOD(ip, 10)
    ip = ip / 10
  END DO
  ii = ii + nd
END DO

RETURN
END SUBROUTINE encode

!**********************************************************************

SUBROUTINE decode(n, nd, gn, ph)
!======================================================================
!     decode genotype into phenotype parameters
!     ph(k) are x, y coordinates [ 0 < x, y < 1 ]
!======================================================================


INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(IN)  :: nd
INTEGER, INTENT(IN)  :: gn(n*nd)
REAL, INTENT(OUT)    :: ph(n)

!     Inputs:


!     Output:


!     Local:
INTEGER :: ip, i, j, ii
REAL :: z

z = 10. ** (-nd)
ii = 0
DO  i = 1, n
  ip = 0
  DO  j = 1, nd
    ip = 10 * ip + gn(ii+j)
  END DO
  ph(i) = ip * z
  ii = ii + nd
END DO

RETURN
END SUBROUTINE decode

!**********************************************************************

SUBROUTINE cross(n, nd, pcross, gn1, gn2)
!======================================================================
!     breeds two parent chromosomes into two offspring chromosomes
!     breeding occurs through crossover starting at position ispl
!======================================================================
!     USES: urand

!     Inputs:
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: nd
REAL, INTENT(IN)         :: pcross

!     Input/Output:
INTEGER, INTENT(IN OUT)  :: gn1(n*nd)
INTEGER, INTENT(IN OUT)  :: gn2(n*nd)

!     Local:
INTEGER :: i, ispl, t

!     Function
! REAL :: urand
! EXTERNAL urand


!     Use crossover probability to decide whether a crossover occurs
IF (urand() < pcross) THEN
  
!        Compute crossover point
  ispl = INT(urand()*n*nd) + 1
  
!        Swap genes at ispl and above
  DO  i = ispl, n * nd
    t = gn2(i)
    gn2(i) = gn1(i)
    gn1(i) = t
  END DO
END IF

RETURN
END SUBROUTINE cross


!**********************************************************************

SUBROUTINE mutate(n, nd, pmut, gn)
!======================================================================
!     Mutations occur at rate pmut at all gene loci
!======================================================================
!     USES: urand

!     Input:
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: nd
REAL, INTENT(IN)         :: pmut

!     Input/Output:
INTEGER, INTENT(IN OUT)  :: gn(n*nd)


!     Local:
INTEGER :: i

!     Function:
! REAL :: urand
! EXTERNAL urand

!     Subject each locus to mutation at the rate pmut
DO  i = 1, n * nd
  IF (urand() < pmut) THEN
    gn(i) = INT(urand()*10.)
  END IF
END DO

RETURN
END SUBROUTINE mutate

!**********************************************************************

SUBROUTINE adjmut(np, fitns, ifit, pmutmn, pmutmx, pmut)
!======================================================================
!     dynamical adjustment of mutation rate; criterion is relative
!     difference in absolute fitnesses of best and median individuals
!======================================================================

!     Input:
INTEGER, INTENT(IN)   :: np
REAL, INTENT(IN)      :: fitns(:)
INTEGER, INTENT(IN)   :: ifit(:)
REAL, INTENT(IN)      :: pmutmn
REAL, INTENT(IN)      :: pmutmx

!     Input/Output:
REAL, INTENT(IN OUT)  :: pmut

!     Local:
REAL  :: rdif
REAL, PARAMETER  :: rdiflo = 0.05, rdifhi = 0.25, delta = 1.5

rdif = ABS(fitns(ifit(np)) - fitns(ifit(np/2))) / (fitns(ifit(np)) +  &
       fitns(ifit(np/2)))
IF (rdif <= rdiflo) THEN
  pmut = MIN(pmutmx, pmut*delta)
ELSE IF (rdif >= rdifhi) THEN
  pmut = MAX(pmutmn, pmut/delta)
END IF

RETURN
END SUBROUTINE adjmut

!**********************************************************************
!                       REPRODUCTION MODULE
!**********************************************************************

!  SELECT:   Parent selection by roulette wheel algorithm
!            called by: PIKAIA

!  RNKPOP:   Ranks initial population
!            called by: PIKAIA, NEWPOP

!  GENREP:   Inserts offspring into population, for full
!            generational replacement
!            called by: PIKAIA

!  STDREP:   Inserts offspring into population, for steady-state
!            reproduction
!            called by: PIKAIA
!            calls:     FF

!  NEWPOP:   Replaces old generation with new generation
!            called by: PIKAIA
!            calls:     FF, RNKPOP

!**********************************************************************

SUBROUTINE select(np, jfit, fdif, idad)
!======================================================================
!     Selects a parent from the population, using roulette wheel
!     algorithm with the relative fitnesses of the phenotypes as
!     the "hit" probabilities [see Davis 1991, chap. 1].
!======================================================================
!     USES: urand

!     Input:
INTEGER, INTENT(IN)   :: np
INTEGER, INTENT(IN)   :: jfit(np)
REAL, INTENT(IN)      :: fdif

!     Output:
INTEGER, INTENT(OUT)  :: idad

!     Local:
INTEGER :: np1, i
REAL :: dice, rtfit

!     Function:
! REAL :: urand
! EXTERNAL urand


np1 = np + 1
dice = urand() * np * np1
rtfit = 0.
DO  i = 1, np
  rtfit = rtfit + np1 + fdif * (np1-2*jfit(i))
  IF (rtfit >= dice) THEN
    idad = i
    GO TO 20
  END IF
END DO
!     Assert: loop will never exit by falling through

20 RETURN
END SUBROUTINE select

!**********************************************************************

SUBROUTINE rnkpop(n, arrin, indx, rank)
!======================================================================
!     Calls external sort routine to produce key index and rank order
!     of input array arrin (which is not altered).
!======================================================================
!     USES: rqsort

!     Input
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: arrin(:)

!     Output
INTEGER, INTENT(OUT)  :: indx(:)
INTEGER, INTENT(OUT)  :: rank(:)


!     Local
INTEGER :: i

!     External sort subroutine
! EXTERNAL rqsort


!     Compute the key index
CALL rqsort(n, arrin, indx)

!     ...and the rank order
DO  i = 1, n
  rank(indx(i)) = n - i + 1
END DO
RETURN
END SUBROUTINE rnkpop

!***********************************************************************

SUBROUTINE genrep(ndim, n, np, ip, ph, newph)
!=======================================================================
!     full generational replacement: accumulate offspring into new
!     population array
!=======================================================================

!     Input:
INTEGER, INTENT(IN)  :: ndim
INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(IN)  :: np
INTEGER, INTENT(IN)  :: ip
REAL, INTENT(IN)     :: ph(ndim, 2)

!     Output:
REAL, INTENT(OUT)    :: newph(ndim, np)


!     Local:
INTEGER :: i1, i2, k


!     Insert one offspring pair into new population
i1 = 2 * ip - 1
i2 = i1 + 1
DO  k = 1, n
  newph(k, i1) = ph(k, 1)
  newph(k, i2) = ph(k, 2)
END DO

RETURN
END SUBROUTINE genrep

!**********************************************************************

SUBROUTINE stdrep(ff, ndim, n, np, irep, ielite, ph, oldph, fitns, ifit, jfit, nnew)
!======================================================================
!     steady-state reproduction: insert offspring pair into population
!     only if they are fit enough (replace-random if irep=2 or
!     replace-worst if irep=3).
!======================================================================
!     USES: ff, urand

!     Input:
INTEGER, INTENT(IN)      :: ndim
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: np
INTEGER, INTENT(IN)      :: irep
INTEGER, INTENT(IN)      :: ielite
REAL, INTENT(IN)         :: ph(ndim, 2)

!     Input/Output:
REAL, INTENT(IN OUT)     :: oldph(ndim, np)
REAL, INTENT(IN OUT)     :: fitns(np)
INTEGER, INTENT(IN OUT)  :: ifit(np)
INTEGER, INTENT(IN OUT)  :: jfit(np)

!     Output:
INTEGER, INTENT(OUT)     :: nnew

INTERFACE
  FUNCTION ff(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION ff
END INTERFACE

! EXTERNAL ff

!     Local:
INTEGER :: i, j, k, i1, if1
REAL :: fit

!     External function
! REAL :: urand
! EXTERNAL urand


nnew = 0
loop70:  DO  j = 1, 2
  
!        1. compute offspring fitness (with caller's fitness function)
  fit = ff(n, ph(:, j))
  
!        2. if fit enough, insert in population
  DO  i = np, 1, -1
    IF (fit > fitns(ifit(i))) THEN
      
!              make sure the phenotype is not already in the population
      IF (i < np) THEN
        DO  k = 1, n
          IF (oldph(k, ifit(i+1)) /= ph(k, j)) GO TO 20
        END DO
        CYCLE loop70
      END IF
      
!              offspring is fit enough for insertion, and is unique
      
!              (i) insert phenotype at appropriate place in population
      20 IF (irep == 3) THEN
        i1 = 1
      ELSE IF (ielite == 0 .OR. i == np) THEN
        i1 = INT(urand()*np) + 1
      ELSE
        i1 = INT(urand()*(np-1)) + 1
      END IF
      if1 = ifit(i1)
      fitns(if1) = fit
      DO  k = 1, n
        oldph(k, if1) = ph(k, j)
      END DO
      
!              (ii) shift and update ranking arrays
      IF (i < i1) THEN
        
!                 shift up
        jfit(if1) = np - i
        DO  k = i1 - 1, i + 1, -1
          jfit(ifit(k)) = jfit(ifit(k)) - 1
          ifit(k+1) = ifit(k)
        END DO
        ifit(i+1) = if1
      ELSE
        
!                 shift down
        jfit(if1) = np - i + 1
        DO  k = i1 + 1, i
          jfit(ifit(k)) = jfit(ifit(k)) + 1
          ifit(k-1) = ifit(k)
        END DO
        ifit(i) = if1
      END IF
      nnew = nnew + 1
      CYCLE loop70
    END IF
  END DO
  
END DO loop70

RETURN
END SUBROUTINE stdrep

!**********************************************************************

SUBROUTINE newpop(ff, ielite, ndim, n, np, oldph, newph, ifit, jfit, fitns, nnew)
!======================================================================
!     replaces old population by new; recomputes fitnesses & ranks
!======================================================================
!     USES: ff, rnkpop

!     Input:
INTEGER, INTENT(IN)   :: ielite
INTEGER, INTENT(IN)   :: ndim
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: np

!     Input/Output:
REAL, INTENT(IN OUT)  :: oldph(ndim, np)
REAL, INTENT(IN OUT)  :: newph(ndim, np)

!     Output:
INTEGER, INTENT(OUT)  :: ifit(np)
INTEGER, INTENT(OUT)  :: jfit(np)
REAL, INTENT(OUT)     :: fitns(np)
INTEGER, INTENT(OUT)  :: nnew

INTERFACE
  FUNCTION ff(n, x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL, INTENT(IN)     :: x(:)
    REAL                 :: fn_val
  END FUNCTION ff
END INTERFACE

! EXTERNAL ff

!     Local:
INTEGER :: i, k

nnew = np

!     if using elitism, introduce in new population fittest of old
!     population (if greater than fitness of the individual it is
!     to replace)
IF (ielite == 1 .AND. ff(n, newph(:, 1)) < fitns(ifit(np))) THEN
  DO  k = 1, n
    newph(k, 1) = oldph(k, ifit(np))
  END DO
  nnew = nnew - 1
END IF

!     replace population
DO  i = 1, np
  DO  k = 1, n
    oldph(k, i) = newph(k, i)
  END DO
  
!        get fitness using caller's fitness function
  fitns(i) = ff(n, oldph(:, i))
END DO

!     compute new population fitness rank order
CALL rnkpop(np, fitns, ifit, jfit)

RETURN
END SUBROUTINE newpop

!*********************************************************************

FUNCTION urand() RESULT(fn_val)
!=====================================================================
!  Return the next pseudo-random deviate from a sequence which is
!  uniformly distributed in the interval [0, 1]

!  Uses the function ran0, the "minimal standard" random number
!  generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
!  Comm. ACM 36 No. 7, 105-110, July 1993).
!=====================================================================

!     Input - none

!     Output
REAL  :: fn_val

!     Local
! INTEGER :: iseed
! REAL :: ran0
! EXTERNAL ran0

!     Common block to make iseed visible to rninit (and to save
!     it between calls)
! COMMON /rnseed/ iseed

fn_val = ran0()
RETURN
END FUNCTION urand

!*********************************************************************

SUBROUTINE rninit(seed)
!=====================================================================
!     Initialize random number generator urand with given seed
!=====================================================================

!     Input
INTEGER, INTENT(IN)  :: seed

!     Output - none

!     Local
! INTEGER  :: iseed

!     Common block to communicate with urand
! COMMON /rnseed/ iseed

!     Set the seed value
iseed = seed
IF (iseed <= 0) iseed = 123456
RETURN
END SUBROUTINE rninit

!*********************************************************************

FUNCTION ran0() RESULT(fn_val)
!=====================================================================
!  "Minimal standard" pseudo-random number generator of Park and Miller.
!  Returns a uniform random deviate r s.t. 0 < r < 1.0.
!  Set seed to any non-zero integer value to initialize a sequence, then do
!  not change seed between calls for successive deviates in the sequence.

!  References:
!     Park, S. and Miller, K., "Random Number Generators: Good Ones
!        are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
!     Park, S. and Miller, K., in "Remarks on Choosing and Implementing
!        Random Number Generators", Comm. ACM 36 No. 7, 105-110 (July 1993)
!=====================================================================
! *** Declaration section ***

!     Output:
REAL  :: fn_val

!     Constants:

INTEGER, PARAMETER  :: a = 48271, m = 2147483647, q = 44488, r = 3399

REAL, PARAMETER :: scale = 1./m, eps = 1.2E-7, rnmx = 1. - eps

!     Local:
INTEGER  :: j

! *** Executable section ***

j = iseed / q
iseed = a * (iseed - j*q) - r * j
IF (iseed < 0) iseed = iseed + m
fn_val = MIN(iseed*scale, rnmx)

RETURN
END FUNCTION ran0

!**********************************************************************

SUBROUTINE rqsort(n, a, p)
!======================================================================
!  Return integer array p which indexes array a in increasing order.
!  Array a is not disturbed.  The Quicksort algorithm is used.

!  B. G. Knapp, 86/12/23

!  Reference: N. Wirth, Algorithms and Data Structures/
!  Prentice-Hall, 1986
!======================================================================

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: a(:)
INTEGER, INTENT(OUT)  :: p(:)

!     Constants

INTEGER, PARAMETER :: lgn = 32, q = 11
!        (LGN = log base 2 of maximum n;
!         Q = smallest subfile to use quicksort on)

!     Local:
REAL :: x
INTEGER :: stackl(lgn), stackr(lgn), s, t, l, m, r, i, j

!     Initialize the stack
stackl(1) = 1
stackr(1) = n
s = 1

!     Initialize the pointer array
DO  i = 1, n
  p(i) = i
END DO

20 IF (s > 0) THEN
  l = stackl(s)
  r = stackr(s)
  s = s - 1
  
  30 IF ((r-l) < q) THEN
    
!           Use straight insertion
    DO  i = l + 1, r
      t = p(i)
      x = a(t)
      DO  j = i - 1, l, -1
        IF (a(p(j)) <= x) GO TO 50
        p(j+1) = p(j)
      END DO
      j = l - 1
      50 p(j+1) = t
    END DO
  ELSE
    
!           Use quicksort, with pivot as median of a(l), a(m), a(r)
    m = (l+r) / 2
    t = p(m)
    IF (a(t) < a(p(l))) THEN
      p(m) = p(l)
      p(l) = t
      t = p(m)
    END IF
    IF (a(t) > a(p(r))) THEN
      p(m) = p(r)
      p(r) = t
      t = p(m)
      IF (a(t) < a(p(l))) THEN
        p(m) = p(l)
        p(l) = t
        t = p(m)
      END IF
    END IF
    
!           Partition
    x = a(t)
    i = l + 1
    j = r - 1
    70 IF (i <= j) THEN
      80 IF (a(p(i)) < x) THEN
        i = i + 1
        GO TO 80
      END IF
      90 IF (x < a(p(j))) THEN
        j = j - 1
        GO TO 90
      END IF
      IF (i <= j) THEN
        t = p(i)
        p(i) = p(j)
        p(j) = t
        i = i + 1
        j = j - 1
      END IF
      GO TO 70
    END IF
    
!           Stack the larger subfile
    s = s + 1
    IF (j-l > r-i) THEN
      stackl(s) = l
      stackr(s) = j
      l = i
    ELSE
      stackl(s) = i
      stackr(s) = r
      r = j
    END IF
    GO TO 30
  END IF
  GO TO 20
END IF
RETURN
END SUBROUTINE rqsort

END MODULE Genetic_Algorithm

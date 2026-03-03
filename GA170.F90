!#######################################################################
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-05  Time: 17:25:32

MODULE GA_commons
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! include 'params.f'
INTEGER, PARAMETER  :: indmax=200, nchrmax=30, nparmax=2

! COMMON /ga1/ npopsiz, nowrite
INTEGER, SAVE  :: npopsiz, nowrite

! COMMON /ga2/ nparam, nchrome
INTEGER, SAVE  :: nparam, nchrome

! COMMON /ga3/ parent, iparent
REAL (dp), SAVE  :: parent(nparmax,indmax)
INTEGER, SAVE    :: iparent(nchrmax,indmax)

! COMMON /ga4/ fitness
REAL (dp), SAVE  :: fitness(indmax)

! COMMON /ga5/ g0, g1, ig2
REAL (dp), SAVE  :: g0(nparmax), g1(nparmax)
INTEGER, SAVE    :: ig2(nparmax)

! COMMON /ga6/ parmax, parmin, pardel, nposibl
REAL (dp), SAVE  :: parmax(nparmax), parmin(nparmax), pardel(nparmax)
INTEGER, SAVE    :: nposibl(nparmax)

! COMMON /ga7/ child, ichild
REAL (dp), SAVE  :: child(nparmax,indmax)
INTEGER, SAVE    :: ichild(nchrmax,indmax)

! COMMON /ga8/ nichflg
INTEGER, SAVE  :: nichflg(nparmax)

! COMMON /inputga/ pcross, pmutate, pcreep, maxgen, idum, irestrt,  &
!      itourny, ielite, icreep, iunifrm, iniche, iskip, iend, nchild,  &
!      microga, kountmx
REAL (dp), SAVE  :: pcross, pmutate, pcreep
INTEGER, SAVE    :: maxgen, idum, irestrt, itourny, ielite, icreep,  &
                    iunifrm, iniche, iskip, iend, nchild, microga, kountmx

END MODULE GA_commons



MODULE ga
!
!  This is version 1.7, last updated on 12/11/98.
!
!  Copyright David L. Carroll; this code may not be reproduced for sale
!  or for use in part of another code for sale without the express
!  written permission of David L. Carroll.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  David L. Carroll
!  University of Illinois
!  306 Talbot Lab
!  104 S. Wright St.
!  Urbana, IL  61801
!
!  e-mail: carroll@uiuc.edu
!  Phone:  217-333-4741
!  fax:    217-244-0720
!
!  This genetic algorithm (GA) driver is free for public use.  My only
!  request is that the user reference and/or acknowledge the use of this
!  driver in any papers/reports/articles which have results obtained
!  from the use of this driver.  I would also appreciate a copy of such
!  papers/articles/reports, or at least an e-mail message with the
!  reference so I can get a copy.  Thanks.
!
!  This program is a FORTRAN version of a genetic algorithm driver.
!  This code initializes a random sample of individuals with different
!  parameters to be optimized using the genetic algorithm approach, i.e.
!  evolution via survival of the fittest.  The selection scheme used is
!  tournament selection with a shuffling technique for choosing random
!  pairs for mating.  The routine includes binary coding for the
!  individuals, jump mutation, creep mutation, and the option for
!  single-point or uniform crossover.  Niching (sharing) and an option
!  for the number of children per pair of parents has been added.
!  An option to use a micro-GA is also included.
!
!  For companies wishing to link this GA driver with an existing code,
!  I am available for some consulting work.  Regardless, I suggest
!  altering this code as little as possible to make future updates
!  easier to incorporate.
!
!  Any users new to the GA world are encouraged to read David Goldberg's
!  "Genetic Algorithms in Search, Optimization and Machine Learning,"
!  Addison-Wesley, 1989.
!
!  Other associated files are:  ga.inp
!                               ga.out
!                               ga.res
!                               params.f
!                               ReadMe
!                               ga2.inp (w/ different namelist identifier)
!
!  I have provided a sample subroutine "func", but ultimately
!  the user must supply this subroutine "func" which should be your
!  cost function.  You should be able to run the code with the
!  sample subroutine "func" and the provided ga.inp file and obtain
!  the optimal function value of 1.0000 at generation 187 with the
!  uniform crossover micro-GA enabled (this is 935 function evaluations)
!
!  The code is presently set for a maximum population size of 200,
!  30 chromosomes (binary bits) and 8 parameters.  These values can be
!  changed in params.f as appropriate for your problem.  Correspondingly
!  you will have to change a few 'write' and 'format' statements if you
!  change nchrome and/or nparam.  In particular, if you change nchrome
!  and/or nparam, then you should change the 'format' statement numbers
!  1050, 1075, 1275, and 1500 (see ReadMe file).
!
!  Please feel free to contact me with questions, comments, or errors
!  (hopefully none of latter).
!
!  Disclaimer:  this program is not guaranteed to be free of error
!  (although it is believed to be free of error), therefore it should
!  not be relied on for solving problems where an error could result in
!  injury or loss.  If this code is used for such solutions, it is
!  entirely at the user's risk and the author disclaims all liability.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      real*4 cpu,cpu0,cpu1,tarray(2)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Input variable definitions:
!
!  icreep   = 0 for no creep mutations
!           = 1 for creep mutations; creep mutations are recommended.
!  idum     The initial random number seed for the GA run.  Must equal
!           a negative integer, e.g. idum=-1000.
!  ielite   = 0 for no elitism (best individual not necessarily
!               replicated from one generation to the next).
!           = 1 for elitism to be invoked (best individual replicated
!               into next generation); elitism is recommended.
!  iend         = 0 for normal GA run (this is standard).
!           = number of last population member to be looked at in a set
!             of individuals.  Setting iend-0 is only used for debugging
!             purposes and is commonly used in conjunction with iskip.
!  iniche   = 0 for no niching
!           = 1 for niching; niching is recommended.
!  irestrt  = 0 for a new GA run, or for a single function evaluation
!           = 1 for a restart continuation of a GA run.
!  iskip    = 0 for normal GA run (this is standard).
!           = number in population to look at a specific individual or
!             set of individuals.  Setting iskip-0 is only used for
!             debugging purposes.
!  itourny  No longer used.  The GA is presently set up for only
!           tournament selection.
!  iunifrm  = 0 for single-point crossover
!           = 1 for uniform crossover; uniform crossover is recommended.
!  kountmx  = the maximum value of kount before a new restart file is
!             written; presently set to write every fifth generation.
!             Increasing this value will reduce I/O time requirements
!             and reduce wear and tear on your storage device
!  maxgen   The maximum number of generations to run by the GA.
!           For a single function evaluation, set equal to 1.
!  microga  = 0 for normal conventional GA operation
!           = 1 for micro-GA operation (this will automatically reset
!             some of the other input flags).  I recommend using
!             npopsiz=5 when microga=1.
!  nchild   = 1 for one child per pair of parents (this is what I
!               typically use).
!           = 2 for two children per pair of parents (2 is more common
!               in GA work).
!  nichflg  = array of 1/0 flags for whether or not niching occurs on
!             a particular parameter.  Set to 0 for no niching on
!             a parameter, set to 1 for niching to operate on parameter.
!             The default value is 1, but the implementation of niching
!             is still controlled by the flag iniche.
!  nowrite  = 0 to write detailed mutation and parameter adjustments
!           = 1 to not write detailed mutation and parameter adjustments
!  nparam   Number of parameters (groups of bits) of each individual.
!           Make sure that nparam matches the number of values in the
!           parmin, parmax and nposibl input arrays.
!  npopsiz  The population size of a GA run (typically 100 works well).
!           For a single calculation, set equal to 1.
!  nposibl  = array of integer number of possibilities per parameter.
!             For optimal code efficiency set nposibl=2**n, i.e. 2, 4,
!             8, 16, 32, 64, etc.
!  parmax   = array of the maximum allowed values of the parameters
!  parmin   = array of the minimum allowed values of the parameters
!  pcreep   The creep mutation probability.  Typically set this
!           = (nchrome/nparam)/npopsiz.
!  pcross   The crossover probability.  For single-point crossover, a
!           value of 0.6 or 0.7 is recommended.  For uniform crossover,
!           a value of 0.5 is suggested.
!  pmutate  The jump mutation probability.  Typically set = 1/npopsiz.
!
!
!  For single function evaluations, set npopsiz=1, maxgen=1, & irestrt=0
!
!  My favorite initial choices of GA parameters are:
!     microga=1, npopsiz=5, iunifrm=1, maxgen=200
!     microga=1, npopsiz=5, iunifrm=0, maxgen=200
!  I generally get good performance with both the uniform and single-
!  point crossover micro-GA.
!
!  For those wishing to use the more conventional GA techniques,
!  my old favorite choice of GA parameters was:
!     iunifrm=1, iniche=1, ielite=1, itourny=1, nchild=1
!  For most problems I have dealt with, I get good performance using
!     npopsiz=100, pcross=0.5, pmutate=0.01, pcreep=0.02, maxgen=26
!  or
!     npopsiz= 50, pcross=0.5, pmutate=0.02, pcreep=0.04, maxgen=51
!
!  Any negative integer for idum should work.  I typically arbitrarily
!  choose idum=-10000 or -20000.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Code variable definitions (those not defined above):
!
!  best     = the best fitness of the generation
!  child    = the floating point parameter array of the children
!  cpu      = cpu time of the calculation
!  cpu0,cpu1= cpu times associated with 'etime' timing function
!  creep    = +1 or -1, indicates which direction parameter creeps
!  delta    = del/nparam
!  diffrac  = fraction of total number of bits which are different
!             between the best and the rest of the micro-GA population.
!             Population convergence arbitrarily set as diffrac<0.05.
!  evals    = number of function evaluations
!  fbar     = average fitness of population
!  fitness  = array of fitnesses of the parents
!  fitsum   = sum of the fitnesses of the parents
!  genavg   = array of average fitness values for each generation
!  geni     = generation array
!  genmax   = array of maximum fitness values for each generation
!  g0       = lower bound values of the parameter array to be optimized.
!             The number of parameters in the array should match the
!             dimension set in the above parameter statement.
!  g1       = the increment by which the parameter array is increased
!             from the lower bound values in the g0 array.  The minimum
!             parameter value is g0 and the maximum parameter value
!             equals g0+g1*(2**g2-1), i.e. g1 is the incremental value
!             between min and max.
!  ig2      = array of the number of bits per parameter, i.e. the number
!             of possible values per parameter.  For example, ig2=2 is
!             equivalent to 4 (=2**2) possibilities, ig2=4 is equivalent
!             to 16 (=2**4) possibilities.
!  ig2sum   = sum of the number of possibilities of ig2 array
!  ibest    = binary array of chromosomes of the best individual
!  ichild   = binary array of chromosomes of the children
!  icount   = counter of number of different bits between best
!             individual and other members of micro-GA population
!  icross   = the crossover point in single-point crossover
!  indmax   = maximum # of individuals allowed, i.e. max population size
!  iparent  = binary array of chromosomes of the parents
!  istart   = the generation to be started from
!  jbest    = the member in the population with the best fitness
!  jelite   = a counter which tracks the number of bits of an individual
!             which match those of the best individual
!  jend     = used in conjunction with iend for debugging
!  jstart   = used in conjunction with iskip for debugging
!  kount    = a counter which controls how frequently the restart
!             file is written
!  kelite   = kelite set to unity when jelite=nchrome, indicates that
!             the best parent was replicated amongst the children
!  mate1    = the number of the population member chosen as mate1
!  mate2    = the number of the population member chosen as mate2
!  nchrmax  = maximum # of chromosomes (binary bits) per individual
!  nchrome  = number of chromosomes (binary bits) of each individual
!  ncreep   = # of creep mutations which occurred during reproduction
!  nmutate  = # of jump mutations which occurred during reproduction
!  nparmax  = maximum # of parameters which the chromosomes make up
!  paramav  = the average of each parameter in the population
!  paramsm  = the sum of each parameter in the population
!  parent   = the floating point parameter array of the parents
!  pardel   = array of the difference between parmax and parmin
!  rand     = the value of the current random number
!  npossum  = sum of the number of possible values of all parameters
!  tarray   = time array used with 'etime' timing function
!  time0    = clock time at start of run
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Subroutines:
!  ____________
!
!  code     = Codes floating point value to binary string.
!  crosovr  = Performs crossover (single-point or uniform).
!  decode   = Decodes binary string to floating point value.
!  evalout  = Evaluates the fitness of each individual and outputs
!             generational information to the 'ga.out' file.
!  func     = The function which is being evaluated.
!  gamicro  = Implements the micro-GA technique.
!  input    = Inputs information from the 'ga.inp' file.
!  initial  = Program initialization and inputs information from the
!             'ga.restart' file.
!  mutate   = Performs mutation (jump and/or creep).
!  newgen   = Writes child array back into parent array for new
!             generation; also checks to see if best individual was
!             replicated (elitism).
!  niche    = Performs niching (sharing) on population.
!  possibl  = Checks to see if decoded binary string falls within
!             specified range of parmin and parmax.
!  ran3     = The random number generator.
!  restart  = Writes the 'ga.restart' file.
!  select   = A subroutine of 'selectn'.
!  selectn  = Performs selection; tournament selection is the only
!             option in this version of the code.
!  shuffle  = Shuffles the population randomly for selection.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

USE GA_commons
IMPLICIT NONE


CONTAINS


!#######################################################################

SUBROUTINE initial(istart, npossum, ig2sum)
!
!  This subroutine sets up the program by generating the g0, g1 and
!  ig2 arrays, and counting the number of chromosomes required for the
!  specified input.  The subroutine also initializes the random number
!  generator, parent and iparent arrays (reads the ga.restart file).

INTEGER, INTENT(OUT)  :: istart
INTEGER, INTENT(OUT)  :: npossum
INTEGER, INTENT(OUT)  :: ig2sum

INTEGER    :: i, j, k, l, n2j
REAL (dp)  :: rand
!
!
DO  i = 1, nparam
  g0(i) = parmin(i)
  pardel(i) = parmax(i) - parmin(i)
  g1(i) = pardel(i) / DBLE(nposibl(i)-1)
END DO
DO  i = 1, nparam
  DO  j = 1, 30
    n2j = 2 ** j
    IF (n2j >= nposibl(i)) THEN
      ig2(i) = j
      EXIT
    END IF
    IF (j >= 30) THEN
      WRITE (6,5100)
      WRITE (24,5100)
      CLOSE (24)
      STOP
    END IF
  END DO
END DO
!
!  Count the total number of chromosomes (bits) required
nchrome = 0
npossum = 0
ig2sum = 0
DO  i = 1, nparam
  nchrome = nchrome + ig2(i)
  npossum = npossum + nposibl(i)
  ig2sum = ig2sum + (2**ig2(i))
END DO
IF (nchrome > nchrmax) THEN
  WRITE (6,5000) nchrome
  WRITE (24,5000) nchrome
  CLOSE (24)
  STOP
END IF
!
IF (npossum < ig2sum.AND.microga /= 0) THEN
  WRITE (6,5200)
  WRITE (24,5200)
END IF
!
!  Initialize random number generator
CALL ran3(idum, rand)
!
IF (irestrt == 0) THEN
!  Initialize the random distribution of parameters in the individual
!  parents when irestrt=0.
  istart = 1
  DO  i = 1, npopsiz
    DO  j = 1, nchrome
      idum = 1
      CALL ran3(idum, rand)
      iparent(j,i) = 1
      IF (rand < 0.5D0) iparent(j,i) = 0
    END DO
  END DO
  IF (npossum < ig2sum) CALL possibl(parent, iparent)
ELSE
!  If irestrt.ne.0, read from restart file.
  OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
  REWIND (25)
  READ (25,*) istart, npopsiz
  DO  j = 1, npopsiz
    READ (25,*) k, (iparent(l,j),l = 1,nchrome)
  END DO
  CLOSE (25)
END IF
!
IF (irestrt /= 0) THEN
  l = idum-istart
  CALL ran3(l, rand)
END IF
!
!
RETURN

5000 FORMAT (' ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
5100 FORMAT (' ERROR: You have a parameter with a number of '/  &
             '    possibilities > 2**30!  If you really desire this,'/  &
             '    change the DO loop 7 statement and recompile.'//  &
             '    You may also need to alter the code to work with'/  &
             '    REAL numbers rather than INTEGER numbers; Fortran'/  &
             '    does not like to compute 2**j when j>30.')
5200 FORMAT (' WARNING: for some cases, a considerable performance'/  &
             '    reduction has been observed when running a non-'/  &
             '    optimal number of bits with the micro-GA.'/  &
             '    If possible, use values for nposibl of 2**n,'/  &
             '    e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
END SUBROUTINE initial
!#######################################################################

SUBROUTINE evalout(iiskip, iiend, ibest, fbar, best)
!
!  This subroutine evaluates the population, assigns fitness,
!  establishes the best individual, and outputs information.

INTEGER, INTENT(IN)     :: iiskip
INTEGER, INTENT(IN)     :: iiend
INTEGER, INTENT(OUT)    :: ibest(nchrmax)
REAL (dp), INTENT(OUT)  :: fbar
REAL (dp), INTENT(OUT)  :: best

INTERFACE
  SUBROUTINE func(j, funcval)
    USE GA_commons
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: j
    REAL (dp), INTENT(OUT)  :: funcval
  END SUBROUTINE func
END INTERFACE

!
INTEGER    :: j, jend, jstart, k, kk, n
REAL (dp)  :: fitsum, funcval, paramsm(nparmax), paramav(nparmax)
!
fitsum = 0.0D0
best = -1.0D10
paramsm(1:nparam) = 0.0D0
jstart = 1
jend = npopsiz
IF (iiskip /= 0) jstart = iiskip
IF (iiend /= 0) jend = iiend
DO  j = jstart, jend
  CALL decode(j, parent, iparent)
  IF (iiskip /= 0 .AND. iiend /= 0 .AND. iiskip == iiend) WRITE (6,5000)  &
      j, (iparent(k,j),k = 1,nchrome), (parent(kk,j),kk = 1,nparam ), 0.0
!
!  Call function evaluator, write out individual and fitness, and add
!  to the summation for later averaging.
  CALL func(j, funcval)
  fitness(j) = funcval
  WRITE (24,5000) j, (iparent(k,j),k = 1,nchrome), (parent(kk,j),  &
                  kk = 1,nparam), fitness(j)
  fitsum = fitsum + fitness(j)
  DO  n = 1, nparam
    paramsm(n) = paramsm(n) + parent(n,j)
  END DO
!
!  Check to see if fitness of individual j is the best fitness.
  IF (fitness(j) > best) THEN
    best = fitness(j)
    DO  k = 1, nchrome
      ibest(k) = iparent(k,j)
    END DO
  END IF
END DO
!
!  Compute parameter and fitness averages.
fbar = fitsum / DBLE(npopsiz)
DO  n = 1, nparam
  paramav(n) = paramsm(n) / DBLE(npopsiz)
END DO
!
!  Write output information
IF (npopsiz == 1) THEN
  WRITE (24,5000) 1, (iparent(k,1),k = 1,nchrome), (parent(k,1),k  &
                  = 1,nparam), fitness(1)
  WRITE (24,*) ' Average Values:'
  WRITE (24,5300) (parent(k,1),k = 1,nparam), fbar
ELSE
  WRITE (24,5300) (paramav(k),k = 1,nparam), fbar
END IF
WRITE (6,5100) fbar
WRITE (24,5100) fbar
WRITE (6,5200) best
WRITE (24,5200) best
!
RETURN

5000 FORMAT (i3, ' ', 30I1, 2(' ', f7.4), ' ', f8.5)
5100 FORMAT (' Average Function Value of Generation=', f8.5)
5200 FORMAT (' Maximum Function Value              =', f8.5/)
5300 FORMAT (/' Average Values:', t34, 2(' ', f7.4), ' ', f8.5/)
END SUBROUTINE evalout
!#######################################################################

SUBROUTINE niche()
!
!  Implement "niching" through Goldberg's multidimensional phenotypic
!  sharing scheme with a triangular sharing function.  To find the
!  multidimensional distance from the best individual, normalize all
!  parameter differences.
!
!
!   Variable definitions:
!
!  alpha   = power law exponent for sharing function; typically = 1.0
!  del     = normalized multidimensional distance between ii and all
!            other members of the population
!            (equals the square root of del2)
!  del2    = sum of the squares of the normalized multidimensional
!            distance between member ii and all other members of
!            the population
!  nniche  = number of niched parameters
!  sigshar = normalized distance to be compared with del; in some sense,
!            1/sigshar can be viewed as the number of regions over which
!            the sharing function should focus, e.g. with sigshar=0.1,
!            the sharing function will try to clump in ten distinct
!            regions of the phase space.  A value of sigshar on the
!            order of 0.1 seems to work best.
!  share   = sharing function between individual ii and j
!  sumshar = sum of the sharing functions for individual ii
!
!      alpha=1.0

REAL (dp)  :: del, del2, share, sigshar, sumshar
INTEGER    :: ii, j, jj, k, nniche

sigshar = 0.1D0
nniche = 0
DO  jj = 1, nparam
  nniche = nniche + nichflg(jj)
END DO
IF (nniche == 0) THEN
  WRITE (6,5000)
  WRITE (24,5000)
  CLOSE (24)
  STOP
END IF
DO  ii = 1, npopsiz
  sumshar = 0.0D0
  DO  j = 1, npopsiz
    del2 = 0.0D0
    DO  k = 1, nparam
      IF (nichflg(k) /= 0) THEN
        del2 = del2 + ((parent(k,j)-parent(k,ii))/pardel(k)) ** 2
      END IF
    END DO
    del = SQRT(del2) / nniche
    IF (del < sigshar) THEN
!               share=1.0 - ((del/sigshar)**alpha)
      share = 1.0D0 - (del/sigshar)
    ELSE
      share = 0.0D0
    END IF
    sumshar = sumshar + share / npopsiz
  END DO
  IF (sumshar /= 0.0D0) fitness(ii) = fitness(ii) / sumshar
END DO
!
!
RETURN
5000 FORMAT (' ERROR: iniche=1 and all values in nichflg array = 0'/  &
             '       Do you want to niche or not?')
END SUBROUTINE niche
!#######################################################################

SUBROUTINE selectn(ipick, j, mate1, mate2)
!
!  Subroutine for selection operator.  Presently, tournament selection
!  is the only option available.
!

INTEGER, INTENT(IN OUT)  :: ipick
INTEGER, INTENT(IN)      :: j
INTEGER, INTENT(IN OUT)  :: mate1
INTEGER, INTENT(IN OUT)  :: mate2

INTEGER  :: n
!
!  If tournament selection is chosen (i.e. itourny=1), then
!  implement "tournament" selection for selection of new population.
IF (itourny == 1) THEN
  CALL select(mate1, ipick)
  CALL select(mate2, ipick)
!        write(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
  DO  n = 1, nchrome
    ichild(n,j) = iparent(n, mate1)
    IF (nchild == 2) ichild(n,j+1) = iparent(n,mate2)
  END DO
END IF
!
RETURN
END SUBROUTINE selectn
!#######################################################################

SUBROUTINE crosovr(ncross, j, mate1, mate2)
!
!  Subroutine for crossover between the randomly selected pair.

INTEGER, INTENT(OUT)     :: ncross
INTEGER, INTENT(IN)      :: j
INTEGER, INTENT(IN OUT)  :: mate1
INTEGER, INTENT(IN OUT)  :: mate2

REAL (dp)  :: rand
INTEGER    :: icross, n
!
IF (iunifrm == 0) THEN
!  Single-point crossover at a random chromosome point.
  idum = 1
  CALL ran3(idum, rand)
  IF (rand > pcross) GO TO 30
  ncross = ncross + 1
  idum = 1
  CALL ran3(idum, rand)
  icross = 2 + INT(DBLE(nchrome-1)*rand)
  DO  n = icross, nchrome
    ichild(n,j) = iparent(n,mate2)
    IF (nchild == 2) ichild(n,j+1) = iparent(n,mate1)
  END DO
ELSE
!  Perform uniform crossover between the randomly selected pair.
  DO  n = 1, nchrome
    idum = 1
    CALL ran3(idum, rand)
    IF (rand <= pcross) THEN
      ncross = ncross + 1
      ichild(n,j) = iparent(n,mate2)
      IF (nchild == 2) ichild(n,j+1) = iparent(n,mate1)
    END IF
  END DO
END IF
!
30 RETURN
END SUBROUTINE crosovr
!#######################################################################

SUBROUTINE mutate()
!

REAL (dp)  :: creep, rand
INTEGER    :: j, k, ncreep, nmutate
!
!  This subroutine performs mutations on the children generation.
!  Perform random jump mutation if a random number is less than pmutate.
!  Perform random creep mutation if a different random number is less
!  than pcreep.
nmutate = 0
ncreep = 0
DO  j = 1, npopsiz
  DO  k = 1, nchrome
!  Jump mutation
    idum = 1
    CALL ran3(idum, rand)
    IF (rand <= pmutate) THEN
      nmutate = nmutate + 1
      IF (ichild(k,j) == 0) THEN
        ichild(k,j) = 1
      ELSE
        ichild(k,j) = 0
      END IF
      IF (nowrite == 0) WRITE (6,5100) j, k
      IF (nowrite == 0) WRITE (24,5100) j, k
    END IF
  END DO

!  Creep mutation (one discrete position away).
  IF (icreep /= 0) THEN
    DO  k = 1, nparam
      idum = 1
      CALL ran3(idum, rand)
      IF (rand <= pcreep) THEN
        CALL decode(j, child, ichild)
        ncreep = ncreep + 1
        creep = 1.0D0
        idum = 1
        CALL ran3(idum, rand)
        IF (rand < 0.5D0) creep = -1.0D0
        child(k,j) = child(k,j) + g1(k) * creep
        IF (child(k,j) > parmax(k)) THEN
          child(k,j) = parmax(k) - 1.0D0 * g1(k)
        ELSE IF (child(k,j) < parmin(k)) THEN
          child(k,j) = parmin(k) + 1.0D0 * g1(k)
        END IF
        CALL code(j, k, child, ichild)
        IF (nowrite == 0) WRITE (6,5200) j, k
        IF (nowrite == 0) WRITE (24,5200) j, k
      END IF
    END DO
  END IF
END DO
WRITE (6,5000) nmutate, ncreep
WRITE (24,5000) nmutate, ncreep
!
!
RETURN
5000 FORMAT (/'  Number of Jump Mutations  =', i5/  &
              '  Number of Creep Mutations =', i5)
5100 FORMAT ('*** Jump mutation performed on individual  ', i4,  &
             ', chromosome ', i3, ' ***')
5200 FORMAT ('*** Creep mutation performed on individual ', i4,  &
             ', parameter  ', i3, ' ***')
END SUBROUTINE mutate
!#######################################################################

SUBROUTINE newgen(iielite, npossum, ig2sum, ibest)
!
!  Write child array back into parent array for new generation.  Check
!  to see if the best parent was replicated; if not, and if ielite=1,
!  then reproduce the best parent into a random slot.
!

INTEGER, INTENT(IN)  :: iielite
INTEGER, INTENT(IN)  :: npossum
INTEGER, INTENT(IN)  :: ig2sum
INTEGER, INTENT(IN)  :: ibest(nchrmax)

INTEGER    :: irand, j, jelite, kelite, n
REAL (dp)  :: rand
!
IF (npossum < ig2sum) CALL possibl(child, ichild)
kelite = 0
DO  j = 1, npopsiz
  jelite = 0
  DO  n = 1, nchrome
    iparent(n,j) = ichild(n,j)
    IF (iparent(n,j) == ibest(n)) jelite = jelite + 1
    IF (jelite == nchrome) kelite = 1
  END DO
END DO
IF (iielite /= 0 .AND. kelite == 0) THEN
  idum = 1
  CALL ran3(idum, rand)
  irand = 1D0 + INT(DBLE(npopsiz)*rand)
  DO  n = 1, nchrome
    iparent(n,irand) = ibest(n)
  END DO
  WRITE (24,5000) irand
END IF
!
!
RETURN
5000 FORMAT ('  Elitist Reproduction on Individual ', i4)
END SUBROUTINE newgen
!#######################################################################

SUBROUTINE gamicro(i, npossum, ig2sum, ibest)
!
!  Micro-GA implementation subroutine
!

INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(IN)  :: npossum
INTEGER, INTENT(IN)  :: ig2sum
INTEGER, INTENT(IN)  :: ibest(nchrmax)

INTEGER    :: icount, j, n
REAL (dp)  :: diffrac, rand

!
!  First, check for convergence of micro population.
!  If converged, start a new generation with best individual and fill
!  the remainder of the population with new randomly generated parents.
!
!  Count number of different bits from best member in micro-population
icount = 0
DO  j = 1, npopsiz
  DO  n = 1, nchrome
    IF (iparent(n,j) /= ibest(n)) icount = icount + 1
  END DO
END DO
!
!  If icount less than 5% of number of bits, then consider population
!  to be converged.  Restart with best individual and random others.
diffrac = DBLE(icount) / DBLE((npopsiz-1)*nchrome)
IF (diffrac < 0.05D0) THEN
  DO  n = 1, nchrome
    iparent(n,1) = ibest(n)
  END DO
  DO  j = 2, npopsiz
    DO  n = 1, nchrome
      idum = 1
      CALL ran3(idum, rand)
      iparent(n,j) = 1
      IF (rand < 0.5D0) iparent(n,j) = 0
    END DO
  END DO
  IF (npossum < ig2sum) CALL possibl(parent,iparent)
  WRITE (6,5000) i
  WRITE (24,5000) i
END IF
!
!
RETURN
5000 FORMAT (//'%%%%%%%  Restart micro-population at generation', i5,  &
               '  %%%%%%%')
END SUBROUTINE gamicro
!#######################################################################

SUBROUTINE select(mate, ipick)
!
!  This routine selects the better of two possible parents for mating.
!

INTEGER, INTENT(OUT)     :: mate
INTEGER, INTENT(IN OUT)  :: ipick

INTEGER  :: ifirst, isecond
!
IF (ipick+1 > npopsiz) CALL shuffle(ipick)
ifirst = ipick
isecond = ipick + 1
ipick = ipick + 2
IF (fitness(ifirst) > fitness(isecond)) THEN
  mate = ifirst
ELSE
  mate = isecond
END IF
!     write(3,*)'select', ifirst, isecond, fitness(ifirst), fitness(isecond)
!
RETURN
END SUBROUTINE select
!#######################################################################

SUBROUTINE shuffle(ipick)
!
!  This routine shuffles the parent array and its corresponding fitness
!

INTEGER, INTENT(OUT)  :: ipick

INTEGER    :: iother, itemp, j, n
REAL (dp)  :: rand, temp
!
ipick = 1
DO  j = 1, npopsiz - 1
  n = 1
  CALL ran3(n, rand)
  iother = j + 1 + INT(DBLE(npopsiz-j)*rand)
  DO  n = 1, nchrome
    itemp = iparent(n,iother)
    iparent(n,iother) = iparent(n,j)
    iparent(n,j) = itemp
  END DO
  temp = fitness(iother)
  fitness(iother) = fitness(j)
  fitness(j) = temp
END DO
!
RETURN
END SUBROUTINE shuffle
!#######################################################################

SUBROUTINE decode(i, array, iarray)
!
!  This routine decodes a binary string to a real number.
!

INTEGER, INTENT(IN)     :: i
REAL (dp), INTENT(OUT)  :: array(nparmax,indmax)
INTEGER, INTENT(IN)     :: iarray(nchrmax,indmax)

INTEGER  :: iparam, j, k, l, m
!
l = 1
DO  k = 1, nparam
  iparam = 0
  m = l
  DO  j = m, m + ig2(k) - 1
    l = l + 1
    iparam = iparam + iarray(j,i) * (2**(m+ig2(k)-1-j))
  END DO
  array(k,i) = g0(k) + g1(k) * iparam
END DO
!
RETURN
END SUBROUTINE decode
!#######################################################################

SUBROUTINE code(j, k, array, iarray)
!
!  This routine codes a parameter into a binary string.
!

INTEGER, INTENT(IN)        :: j
INTEGER, INTENT(IN)        :: k
REAL (dp), INTENT(IN OUT)  :: array(nparmax,indmax)
INTEGER, INTENT(OUT)       :: iarray(nchrmax,indmax)

INTEGER  :: i, iparam, istart, m

!
!  First, establish the beginning location of the parameter string of interest.

istart = 1
DO  i = 1, k - 1
  istart = istart + ig2(i)
END DO
!
!  Find the equivalent coded parameter value, and back out the binary
!  string by factors of two.
m = ig2(k) - 1
IF (g1(k) == 0.0D0) RETURN
iparam = NINT((array(k,j)-g0(k))/g1(k))
DO  i = istart, istart + ig2(k) - 1
  iarray(i,j) = 0
  IF (iparam+1 > 2**m) THEN
    iarray(i,j) = 1
    iparam = iparam - 2 ** m
  END IF
  m = m - 1
END DO
!     write(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
!
RETURN
END SUBROUTINE code
!#######################################################################
!

SUBROUTINE possibl(array, iarray)
!
!  This subroutine determines whether or not all parameters are within
!  the specified range of possibility.  If not, the parameter is
!  randomly reassigned within the range.  This subroutine is only
!  necessary when the number of possibilities per parameter is not
!  optimized to be 2**n, i.e. if npossum < ig2sum.
!

REAL (dp), INTENT(OUT)   :: array(nparmax,indmax)
INTEGER, INTENT(IN OUT)  :: iarray(nchrmax,indmax)

INTEGER    :: i, irand, j, n2ig2j
REAL (dp)  :: rand
!
DO  i = 1, npopsiz
  CALL decode(i, array, iarray)
  DO  j = 1, nparam
    n2ig2j = 2 ** ig2(j)
    IF (nposibl(j) /= n2ig2j .AND. array(j,i) > parmax(j)) THEN
      idum = 1
      CALL ran3(idum, rand)
      irand = INT(DBLE(nposibl(j))*rand)
      array(j,i) = g0(j) + irand * g1(j)
      CALL code(i, j, array, iarray)
      IF (nowrite == 0) WRITE (6,5000) i, j
      IF (nowrite == 0) WRITE (24,5000) i, j
    END IF
  END DO
END DO
!
!
RETURN
5000 FORMAT ('*** Parameter adjustment to individual     ',i4,  &
    ', parameter  ',i3,' ***')
END SUBROUTINE possibl
!#######################################################################

SUBROUTINE restart(i, istart, kount)
!
!  This subroutine writes restart information to the ga.restart file.
!

INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(IN)   :: istart
INTEGER, INTENT(OUT)  :: kount

INTEGER  :: j, l

kount = kount + 1
IF (i == maxgen+istart-1 .OR. kount == kountmx) THEN
  OPEN (UNIT=25, FILE='ga.res', STATUS='OLD')
  REWIND (25)
  WRITE (25,*) i + 1, npopsiz
  DO  j = 1, npopsiz
    WRITE (25,5000) j, (iparent(l,j),l = 1,nchrome)
  END DO
  CLOSE (25)
  kount = 0
END IF
!
!
RETURN
5000 FORMAT (i5, '   ', 30I2)
END SUBROUTINE restart
!#######################################################################

SUBROUTINE ran3(iidum, rand)
!
!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', "Numerical Recipes" p. 199.
!

INTEGER, INTENT(IN OUT)  :: iidum
REAL (dp), INTENT(OUT)   :: rand

!      implicit real*4(m)
REAL (dp), PARAMETER :: mbig = 4000000., mseed = 1618033., mz = 0.0,  &
                        fac = 1./mbig
!     parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
!
!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.
INTEGER, SAVE  :: ma(55)
INTEGER, SAVE  :: iff = 0

INTEGER        :: i, ii, k, mj, mk
INTEGER, SAVE  :: inext, inextp

IF (iidum < 0 .OR. iff == 0) THEN
  iff = 1
  mj = mseed - DBLE(ABS(iidum))
  mj = MOD(mj,INT(mbig))
  ma(55) = mj
  mk = 1
  DO  i = 1, 54
    ii = MOD(21*i,55)
    ma(ii) = mk
    mk = mj - mk
    IF (mk < mz) mk = mk + mbig
    mj = ma(ii)
  END DO
  DO  k = 1, 4
    DO  i = 1, 55
      ma(i) = ma(i) - ma(1+MOD(i+30,55))
      IF (ma(i) < mz) ma(i) = ma(i) + mbig
    END DO
  END DO
  inext = 0
  inextp = 31
  iidum = 1
END IF
inext = inext + 1
IF (inext == 56) inext = 1
inextp = inextp + 1
IF (inextp == 56) inextp = 1
mj = ma(inext) - ma(inextp)
IF (mj < mz) mj = mj + mbig
ma(inext) = mj
rand = mj * fac

RETURN
END SUBROUTINE ran3

END MODULE ga



PROGRAM gafortran
! Driver program

USE GA_commons
USE ga
IMPLICIT NONE
!
INTEGER    :: ibest(nchrmax)
REAL (dp)  :: geni(1000000), genavg(1000000), genmax(1000000)


INTEGER    :: i, ig2sum, ipick, istart, j, kount, mate1, mate2, ncross, npossum
REAL (dp)  :: best, evals, fbar
!
!      call etime(tarray)
!      write(6,*) tarray(1),tarray(2)
!      cpu0=tarray(1)
!
!  Call the input subroutine.
!      TIME0=SECNDS(0.0)
CALL INPUT()
!
!  Perform necessary initialization and read the ga.restart file.
CALL initial(istart, npossum, ig2sum)
!
!  $$$$$ Main generational processing loop. $$$$$
DO  i = istart, maxgen + istart - 1
  WRITE (6,5100) i
  WRITE (24,5100) i
  WRITE (24,5000)
!
!  Evaluate the population, assign fitness, establish the best
!  individual, and write output information.
  CALL evalout(iskip, iend, ibest, fbar, best)
  geni(i)   = i
  genavg(i) = fbar
  genmax(i) = best
  IF (npopsiz == 1 .OR. iskip /= 0) THEN
    CLOSE (24)
    STOP
  END IF
!
!  Implement "niching".
  IF (iniche /= 0) CALL niche()
!
!  Enter selection, crossover and mutation loop.
  ncross = 0
  ipick = npopsiz
  DO  j = 1, npopsiz, nchild
!
!  Perform selection.
    CALL selectn(ipick, j, mate1, mate2)
!
!  Now perform crossover between the randomly selected pair.
    CALL crosovr(ncross, j, mate1, mate2)
  END DO
  WRITE (6,5200) ncross
  WRITE (24,5200) ncross
!
!  Now perform random mutations.  If running micro-GA, skip mutation.
  IF (microga == 0) CALL mutate()
!
!  Write child array back into parent array for new generation.  Check
!  to see if the best parent was replicated.
  CALL newgen(ielite, npossum, ig2sum, ibest)
!
!  Implement micro-GA if enabled.
  IF (microga /= 0) CALL gamicro(i, npossum, ig2sum, ibest)
!
!  Write to restart file.
  CALL restart(i, istart, kount)
END DO
!  $$$$$ End of main generational processing loop. $$$$$
! 999  continue
WRITE (24,5300)
DO  i = istart, maxgen + istart - 1
  evals = npopsiz * geni(i)
  WRITE (24,5400) geni(i), evals, genavg(i), genmax(i)
END DO
!      call etime(tarray)
!      write(6,*) tarray(1),tarray(2)
!      cpu1=tarray(1)
!      cpu=(cpu1-cpu0)
!      write(6,1400) cpu,cpu/60.0
!      write(24,1400) cpu,cpu/60.0
CLOSE (24)
!
! 1400 format(2x,'CPU time for all generations=', e12.6, ' sec'/
!     +       2x,'                             ', e12.6, ' min')
!
STOP

5000 FORMAT ('  #      Binary Code                Param1  Param2  Fitness')
5100 FORMAT (//'#################  Generation', i5, '  #################')
5200 FORMAT (/'  Number of Crossovers      =', i5)
5300 FORMAT (//' Summary of Output'/  &
               '  Generation   Evaluations   Avg.Fitness   Best Fitness')
5400 FORMAT (t3, 3(e10.4, '    '), e11.5)

!#######################################################################


CONTAINS


SUBROUTINE INPUT()
!
!  This subroutine inputs information from the ga.inp (gafort.in) file.
!
!
NAMELIST /ga/ irestrt, npopsiz, pmutate, maxgen, idum, pcross,  &
    itourny, ielite, icreep, pcreep, iunifrm, iniche, iskip, iend,  &
    nchild, nparam, parmin, parmax, nposibl, nowrite, nichflg, microga, kountmx
!
kountmx = 5
irestrt = 0
itourny = 0
ielite = 0
iunifrm = 0
iniche = 0
iskip = 0
iend = 0
nchild = 1
nichflg(1:nparam) = 1
microga = 0
!
OPEN (UNIT=24, FILE='ga.out', STATUS='UNKNOWN')
REWIND (24)
OPEN (UNIT=23, FILE='ga.inp', STATUS='OLD')
READ (23, nml=ga)
CLOSE (23)
itourny = 1
!      if (itourny.eq.0) nchild=2
!
!  Check for array sizing errors.
IF (npopsiz > indmax) THEN
  WRITE (6,5000) npopsiz
  WRITE (24,5000) npopsiz
  CLOSE (24)
  STOP
END IF
IF (nparam > nparmax) THEN
  WRITE (6,5100) nparam
  WRITE (24,5100) nparam
  CLOSE (24)
  STOP
END IF
!
!  If using the microga option, reset some input variables
IF (microga /= 0) THEN
  pmutate = 0.0D0
  pcreep = 0.0D0
  itourny = 1
  ielite = 1
  iniche = 0
  nchild = 1
  IF (iunifrm == 0) THEN
    pcross = 1.0D0
  ELSE
    pcross = 0.5D0
  END IF
END IF
!
!
RETURN
5000 FORMAT (' ERROR: npopsiz > indmax.  Set indmax = ', i6)
5100 FORMAT (' ERROR: nparam > nparmax.  Set nparmax = ', i6)
END SUBROUTINE INPUT

END PROGRAM gafortran

!#######################################################################
!

SUBROUTINE func(j, funcval)
!

USE GA_commons
IMPLICIT NONE

INTEGER, INTENT(IN)     :: j
REAL (dp), INTENT(OUT)  :: funcval

!
!      dimension parent2(indmax,nparmax),iparent2(indmax,nchrmax)
!
!  This is an N-dimensional version of the multimodal function with
!  decreasing peaks used by Goldberg and Richardson (1987, see ReadMe
!  file for complete reference).  In N dimensions, this function has
!  (nvalley-1)^nparam peaks, but only one global maximum.  It is a
!  reasonably tough problem for the GA, especially for higher dimensions
!  and larger values of nvalley.
!
INTEGER    :: i, nvalley
REAL (dp)  :: f1, f2, pi

nvalley = 6
pi = 4.0D0 * ATAN(1.d0)
funcval = 1.0D0
DO  i = 1, nparam
  f1 = (SIN(5.1D0*pi*parent(i,j)+0.5D0)) ** nvalley
  f2 = EXP(-4.0D0*LOG(2.0D0)*((parent(i,j)-0.0667D0)**2)/0.64D0)
  funcval = funcval * f1 * f2
END DO
!
!  As mentioned in the ReadMe file, The arrays have been rearranged
!  to enable a more efficient caching of system memory.  If this causes
!  interface problems with existing functions used with previous
!  versions of my code, then you can use some temporary arrays to bridge
!  this version with older versions.  I've named the temporary arrays
!  parent2 and iparent2.  If you want to use these arrays, uncomment the
!  dimension statement above as well as the following do loop lines.
!
!      do 11 i=1,nparam
!         parent2(j,i)=parent(i,j)
! 11   continue
!      do 12 k=1,nchrome
!         iparent2(j,k)=iparent(k,j)
! 12   continue
!
RETURN
END SUBROUTINE func
!#######################################################################

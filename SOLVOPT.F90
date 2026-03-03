MODULE Shor_minimization
IMPLICIT NONE

! This is an include file to compile with solvopt.f(or)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-10  Time: 23:37:59
 
CHARACTER (LEN=80) :: errmes='SolvOpt error:', wrnmes='SolvOpt warning:',  &
          error2='Improper space dimension.',  &
          error32='Function equals infinity at the point.',  &
          error42='Gradient equals infinity at the starting point.', &
          error43='Gradient equals zero at the starting point.',  &
          error52='<FUNC> returns infinite value at the point.',  &
          error62='<GRADC> returns infinite vector at the point.', &
          error63='<GRADC> returns zero vector at an infeasible point.',  &
          error5='Function is unbounded.',   &
          error6='Choose another starting point.',  &
          warn1='Gradient is zero, but stopping criteria are not fulfilled.', &
          warn20='Normal re-setting of a transformation matrix.',  &
          warn21='Re-setting due to the use of a new penalty coefficient.', &
          warn4='Iterations limit exceeded.',   &
          warn31='The function is flat in certain directions.',  &
          warn32='Trying to recover by shifting insensitive variables.',  &
          warn09='Re-run from recorded point.',  &
          warn08='Ravine with a flat bottom is detected.',  &
          termwarn0='SolvOpt: Normal termination.',  &
          termwarn1='SolvOpt: Termination warning:',  &
          appwarn='The above warning may be reasoned'//  &
                  ' by inaccurate gradient approximation',  &
          endwarn1='Premature stop is possible.'//  &
                   ' Try to re-run the routine from the obtained point.',  &
          endwarn2='Result may not provide the optimum.'//  &
                   ' The function apparently has many extremum points.',  &
          endwarn3='Result may be inaccurate in the coordinates.'//  &
                   ' The function is flat at the solution.',  &
          endwarn4='Stopping criteria are not fulfilled.'//  &
                   ' The function is very steep at the solution.'

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


subroutine solvopt(n, x, f, fun, flg, grad, options, flfc, flgc, func, gradc)

! N.B. The 2 last arguments, func & gradc are now optional.
!      Argument flgc has been inserted before argument func.

!-------------------------------------------------------------------------
! The subroutine SOLVOPT performs a modified version of Shor's r-algorithm
! in order to find a local minimum resp. maximum of a nonlinear function
! defined on the n-dimensional Euclidean space
! or
! a local minimum for a nonlinear constrained problem:
! min { f(x): g(x) (<)= 0, g(x) in R(m), x in R(n) }.
! Arguments:
! n       is the space dimension (integer*4),
! x       is the n-vector, the coordinates of the starting point
!         at a call to the subroutine and the optimizer at regular return
!         (REAL (dp)),
! f       returns the optimum function value
!         (REAL (dp)),
! fun     is the entry name of a subroutine which computes the value
!         of the function <fun> at a point x, should be declared as
!         external in a calling routine,
!         synopsis: fun(x,f)
! grad    is the entry name of a subroutine which computes the gradient
!         vector of the function <fun> at a point x, should be declared
!         as external in a calling routine,
!         synopsis: grad(x,g)
! func    is the entry name of a subroutine which computes the MAXIMAL
!         RESIDUAL!!! (a scalar) for a set of constraints at a point x,
!         should be declared as external in a calling routine,
!         synopsis: func(x,fc)
! gradc   is the entry name of a subroutine which computes the gradient
!         vector for a constraint with the MAXIMAL RESIDUAL at a point x,
!         should be declared as external in a calling routine,
!         synopsis: gradc(x,gc)
! flg     (logical) is a flag for the use of a subroutine <grad>:
!         .true. means gradients are calculated by the user-supplied routine.
! flfc    (logical) is a flag for a constrained problem:
!         .true. means the maximal residual for a set of constraints
!         is calculated by <func>.
! flgc    (logical) is a flag for the use of a subroutine <gradc>:
!         .true. means gradients of the constraints are calculated
!         by the user-supplied routine.
! options is a vector of optional parameters (REAL (dp)):
!     options(1)= H, where sign(H)=-1 resp. sign(H)=+1 means minimize resp.
!         maximize <fun> (valid only for an unconstrained problem) and
!         H itself is a factor for the initial trial step size
!         (options(1)=-1.d0 by default),
!     options(2)= relative error for the argument in terms of the infinity-norm
!         (1.d-4 by default),
!     options(3)= relative error for the function value (1.d-6 by default),
!     options(4)= limit for the number of iterations (1.5d4 by default),
!     options(5)= control of the display of intermediate results and error
!         resp. warning messages (default value is 0.d0, i.e., no intermediate
!         output but error and warning messages, see the manual for more),
!     options(6)= maximal admissible residual for a set of constraints
!         (options(6)=1.d-8 by default, see the manual for more),
!    *options(7)= the coefficient of space dilation (2.5d0 by default),
!    *options(8)= lower bound for the stepsize used for the difference
!        approximation of gradients (1.d-11 by default,see the manual for more).
!   (* ... changes should be done with care)
! Returned optional values:
!     options(9),  the number of iterations, if positive,
!         or an abnormal stop code, if negative (see manual for more),
!                -1: allocation error,
!                -2: improper space dimension,
!                -3: <fun> returns an improper value,
!                -4: <grad> returns a zero vector or improper value at the
!                    starting point,
!                -5: <func> returns an improper value,
!                -6: <gradc> returns an improper value,
!                -7: function is unbounded,
!                -8: gradient is zero at the point,
!                    but stopping criteria are not fulfilled,
!                -9: iterations limit exceeded,
!               -11: Premature stop is possible,
!               -12: Result may not provide the true optimum,
!               -13: Function is flat: result may be inaccurate
!                    in view of a point.
!               -14: Function is steep: result may be inaccurate
!                    in view of a function value,
!       options(10), the number of objective function evaluations, and
!       options(11), the number of gradient evaluations.
!       options(12), the number of constraint function evaluations, and
!       options(13), the number of constraint gradient evaluations.
! ____________________________________________________________________________
!

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:), options(13)
REAL (dp), INTENT(OUT)     :: f
LOGICAL, INTENT(IN)        :: flg, flfc, flgc

! external fun,grad,func,gradc

INTERFACE
  SUBROUTINE fun(x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
  END SUBROUTINE fun

  SUBROUTINE grad(x, g)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
  END SUBROUTINE grad

  SUBROUTINE func(x, fc)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fc
  END SUBROUTINE func

  SUBROUTINE gradc(x, gc)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: gc(:)
  END SUBROUTINE gradc
END INTERFACE

OPTIONAL                   :: func, gradc

! Local variables

logical    :: constr, app, appconstr
logical    :: FsbPnt, FsbPnt1, termflag, stopf
logical    :: stopping, dispwarn, Reset, ksm,knan,obj
integer    :: kstore, ajp,ajpp,knorms, k, kcheck, numelem
integer    :: dispdata, ld, mxtc, termx, limxterm, nzero, krerun
integer    :: warnno, kflat, stepvanish, i, j, ni, ii, kd, kj, kc, ip
integer    :: iterlimit, kg, k1, k2, kless, allocerr
REAL (dp)  :: doptions(13)
REAL (dp)  :: nsteps(3), gnorms(10), kk, nx
REAL (dp)  :: ajb,ajs, des, dq, du20, du10, du03
REAL (dp)  :: n_float, cnteps
REAL (dp)  :: low_bound, ZeroGrad, ddx, y
REAL (dp)  :: lowxbound, lowfbound, detfr, detxr, grbnd
REAL (dp)  :: fp, fp1, fc, f1, f2, fm, fopt, frec, fst, fp_rate
REAL (dp)  :: PenCoef, PenCoefNew
REAL (dp)  :: gamma, w, wdef, h1, h, hp
REAL (dp)  :: dx, ng, ngc, nng, ngt, nrmz, ng1, d, dd, laststep
REAL (dp), dimension(:,:), allocatable :: B
REAL (dp), dimension(:), allocatable   :: g, g0, g1, gt, gc, z, x1, xopt,  &
                                          xrec, grec, xx, deltax
integer, dimension(:), allocatable     :: idx
CHARACTER (LEN=100)  :: endwarn

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         three = 3.0_dp, four = 4.0_dp, five = 5.0_dp,  &
                         six = 6.0_dp, seven = 7.0_dp, eight = 8.0_dp,  &
                         nine = 9.0_dp, ten = 1.d1,  hundr = 100.0_dp,  &
                         powerm12 = 1.d-12, infty = 1.d100, epsnorm = 1.d-15, &
                         epsnorm2 = 1.d-30
CHARACTER (LEN=19), PARAMETER  :: allocerrstr = 'Allocation Error = '

! Check the dimension:
if (n < 2) then
    WRITE(*, *) errmes
    WRITE(*, *) error2
  options(9)=-one
  GO TO 999
END IF

n_float=dble(n)

! allocate working arrays:
allocate (B(n,n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (g(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (g0(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (g1(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (gt(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (gc(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (z(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (x1(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (xopt(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (xrec(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (grec(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (xx(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (deltax(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

allocate (idx(n),stat=allocerr)
if (allocerr /= 0) then
   options(9)=-one
   WRITE(*, *) allocerrstr, allocerr
END IF

! store flags:
app=.not.flg
constr=flfc
appconstr=.not.flgc

! Default values for options:
call soptions(doptions)
do i=1,8
   if (options(i) == zero) then
      options(i)=doptions(i)
   ELSE IF (i == 2.or.i == 3.or.i == 6) then
      options(i)=MAX(options(i),powerm12)
      options(i)=MIN(options(i),one)
      if (i == 2)options(i)=MAX(options(i),options(8)*hundr)
   ELSE IF (i == 7) then
      options(7)=MAX(options(i),1.5d0)
   END IF
END DO
               
! WORKING CONSTANTS AND COUNTERS ----{
                        
options(10)=zero    !! counter for function calculations
options(11)=zero    !! counter for gradient calculations
options(12)=zero    !! counter for constraint function calculations
options(13)=zero    !! counter for constraint gradient calculations
iterlimit=INT(options(4))
if (constr) then
  h1=-one           !! NLP: restricted to minimization
  cnteps=options(6)
else
  h1=SIGN(one,options(1))  !! Minimize resp. maximize a function
END IF
k=0                         !! Iteration counter
wdef=one/options(7) - one   !! Default space transf. coeff.

! Gamma control ---{
ajb=one+1.d-1/n_float**2    !! Base I
ajp=20
ajpp=ajp                    !! Start value for the power
ajs=1.15_dp
knorms=0
gnorms(1:10)=zero
!---}

! Display control ---{
if (options(5) <= zero) then
   dispdata=0
   if (options(5) == -one) then
      dispwarn=.false.
   else
      dispwarn=.true.
   END IF
else
   dispdata=NINT(options(5))
   dispwarn=.true.
END IF
ld=dispdata
!---}

! Stepsize control ---{
dq=5.1_dp               !! Step divider (at f_{i+1}>gamma*f_{i})
du20=two
du10=1.5_dp
du03=1.05_dp            !! Step multipliers (at certain steps made)
kstore=3
nsteps(1:kstore)=zero   !! Steps made at the last 'kstore' iterations
if (app) then
  des=6.3_dp            !! Desired number of steps per 1-D search
else
  des=3.3_dp
END IF
mxtc=3                  !! Number of trial cycles (steep wall detect)
!---}

termx=0
limxterm=50             !! Counter and limit for x-criterion

! stepsize for gradient approximation
ddx=MAX(1.d-11, options(8))

low_bound=-one + 1.d-4     !! Lower bound cosine used to detect a ravine
ZeroGrad=n_float*1.d-16    !! Lower bound for a gradient norm
nzero=0                    !! Zero-gradient events counter

! Low bound for the values of variables to take into account
lowxbound=MAX(options(2), 1.d-3)

! Lower bound for function values to be considered as making difference
lowfbound=options(3)**2
krerun=0                 !! Re-run events counter
detfr=options(3)*hundr   !! Relative error for f/f_{record}
detxr=options(2)*ten     !! Relative error for norm(x)/norm(x_{record})
warnno=0                 !! the number of warn.mess. to end with
kflat=0                  !! counter for points of flatness
stepvanish=0             !! counter for vanished steps
stopf=.false.
! ----}  End of setting constants
! ----}  End of the preamble
!--------------------------------------------------------------------
! COMPUTE THE FUNCTION  ( FIRST TIME ) ----{
call fun(x, f)
options(10)=options(10) + one
if (ABS(f) >= infty) then
   if (dispwarn) then
      WRITE(*, *) errmes
      WRITE(*, *) error32
      WRITE(*, *) error6
   END IF
   options(9)=-three
   GO TO 999
END IF
xrec(1:n)=x(1:n)
frec=f     !! record point and function value

! Constrained problem
if (constr) then
    kless=0
    fp=f
    call func(x,fc)
    options(12)=options(12)+one
    if (ABS(fc) >= infty) then
       if (dispwarn) then
          WRITE(*, *) errmes
          WRITE(*, *) error52
          WRITE(*, *) error6
       END IF
       options(9)=-five
       GO TO 999
    END IF
  PenCoef=one          !! first rough approximation
  if (fc <= cnteps) then
   FsbPnt=.true.       !! feasible point
   fc=zero
  else
   FsbPnt=.false.
  END IF
  f=f + PenCoef*fc
END IF
! ----}

! COMPUTE THE GRADIENT ( FIRST TIME ) ----{
if (app) then
  deltax(1:n)=h1*ddx
  obj=.true.
  if (constr) then
     call apprgrdn(n, g, x, fp, fun, deltax, obj)
  else
     call apprgrdn(n, g, x, f, fun, deltax, obj)
  END IF
  options(10)=options(10) + n_float
else
  call grad(x, g)
  options(11)=options(11) + one
END IF
ng=SUM( g(1:n)**2 )
ng=SQRT(ng)
if (ng >= infty) then
   if (dispwarn) then
      WRITE(*, *) errmes
      WRITE(*, *) error42
      WRITE(*, *) error6
   END IF
   options(9)=-four
   GO TO 999
ELSE IF (ng < ZeroGrad) then
   if (dispwarn) then
      WRITE(*, *) errmes
      WRITE(*, *) error43
      WRITE(*, *) error6
   END IF
   options(9)=-four
   GO TO 999
END IF
if (constr) then
  if (.not.FsbPnt) then
    if (appconstr) then
       do j=1,n
         if (x(j) >= zero) then
            deltax(j)=ddx
         else
            deltax(j)=-ddx
         END IF
       END DO
       obj=.false.
       call apprgrdn(n, gc, x, fc, func, deltax, obj)
    else
       call gradc(x, gc)
    END IF
    ngc=SUM( gc(1:n)**2 )
    ngc=SQRT(ngc)
    if (ng >= infty) then
       if (dispwarn) then
          WRITE(*, *) errmes
          WRITE(*, *) error62
          WRITE(*, *) error6
       END IF
       options(9)=-six
       GO TO 999
    ELSE IF (ng < ZeroGrad) then
       if (dispwarn) then
          WRITE(*, *) errmes
          WRITE(*, *) error63
       END IF
       options(9)=-six
       GO TO 999
    END IF
    do i=1,n
      g(i)=g(i) + PenCoef*gc(i)
    END DO
    ng=zero
    do i=1,n
      ng=ng + g(i)*g(i)
      grec(i)=g(i)
    END DO
    ng=SQRT(ng)
  END IF
END IF
grec(1:n)=g(1:n)
nng=ng
! ----}

! INITIAL STEPSIZE
d=zero
do i=1,n
  if (d < ABS(x(i))) d=ABS(x(i))
END DO
h=h1*SQRT(options(2))*d                  !! smallest possible stepsize
if (ABS(options(1)) /= one) then
   h=h1*MAX(ABS(options(1)), ABS(h))     !! user-supplied stepsize
else
   h=h1*MAX(one/LOG(ng+1.1_dp), ABS(h))  !! calculated stepsize
END IF

! RESETTING LOOP ----{
do while (.true.)
  kcheck=0                       !! Set checkpoint counter.
  kg=0                           !! stepsizes stored
  kj=0                           !! ravine jump counter
  do i=1,n
    B(i,1:n)=zero
    B(i,i)=one                   !! re-set transf. matrix to identity
    g1(i)=g(i)
  END DO
  fst=f
  dx=0
! ----}

! MAIN ITERATIONS ----{
   
  do while (.true.)
    k=k + 1
    kcheck=kcheck + 1
    laststep=dx
! ADJUST GAMMA --{
    gamma=one + MAX(ajb**((ajp-kcheck)*n), two*options(3))
    gamma=MIN ( gamma, ajs**MAX(one, LOG10(nng+one)) )
! --}
    ngt=zero
    ng1=zero
    dd=zero
    do i=1,n
      d=DOT_PRODUCT( B(1:n,i), g(1:n) )
      gt(i)=d
      dd=dd + d*g1(i)
      ngt=ngt + d*d
      ng1=ng1 + g1(i)*g1(i)
    END DO
    ngt=SQRT(ngt)
    ng1=SQRT(ng1)
    dd=dd/ngt/ng1

    w=wdef

! JUMPING OVER A RAVINE ----{
    if (dd < low_bound) then
      if (kj == 2) then
        xx(1:n)=x(1:n)
      END IF
      if (kj == 0) kd=4
      kj=kj + 1
      w=-.9_dp              !! use large coef. of space dilation
      h=h*two
      if (kj > 2*kd) then
        kd=kd+1
        warnno=1
        endwarn=endwarn1
        do i=1,n
          if (ABS(x(i)-xx(i)) < epsnorm*ABS(x(i))) then
            if (dispwarn)  then
               WRITE(*, *) wrnmes
               WRITE(*, *) warn08
            END IF
          END IF
        END DO
      END IF
    else
      kj=0
    END IF
! ----}

! DILATION ----{
    nrmz=zero
    do i=1,n
      z(i)=gt(i) - g1(i)
      nrmz=nrmz + z(i)*z(i)
    END DO
    nrmz=SQRT(nrmz)
    if (nrmz > epsnorm*ngt) then
      z(1:n)=z(1:n)/nrmz

! New direction in the transformed space: g1=gt+w*(z*gt')*z and
! new inverse matrix: B = B ( I + (1/alpha -1)zz' )
      d = DOT_PRODUCT( z(1:n), gt(1:n) )
      ng1=zero
      d = d*w
      do i=1,n
        dd=zero
        g1(i)=gt(i) + d*z(i)
        ng1=ng1 + g1(i)*g1(i)
        do j=1,n
           dd=dd + B(i,j)*z(j)
        END DO
        dd=w*dd
        B(i,1:n)=B(i,1:n) + dd*z(1:n)
      END DO
      ng1=SQRT(ng1)
    else
       do i=1,n
         z(i)=zero
         g1(i)=gt(i)  
       END DO
       nrmz=zero
    END IF
    gt(1:n)=g1(1:n)/ng1
    do i=1,n
      d=DOT_PRODUCT( B(i,1:n), gt(1:n) )
      g0(i)=d
    END DO
! ----}

! RESETTING ----{
    if (kcheck > 1) then
       numelem=0
       do i=1,n
          if (ABS(g(i)) > ZeroGrad) then
             numelem=numelem + 1
             idx(numelem)=i
          END IF
       END DO
       if (numelem > 0) then
          grbnd=epsnorm*dble(numelem**2)
          ii=0
          do i=1,numelem
             j=idx(i)
             if (ABS(g1(j)) <= ABS(g(j))*grbnd) ii=ii + 1
          END DO
          if (ii == n .or. nrmz == zero) then
            if (dispwarn) then
              WRITE(*, *) wrnmes
              WRITE(*, *) warn20
            END IF
            if (ABS(fst-f) < ABS(f)*1.d-2) then
               ajp=ajp - 10*n
            else
               ajp=ajpp
            END IF
            h=h1*dx/three
            k=k - 1
            exit
          END IF
       END IF
    END IF
! ----}

! STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH
    xopt(1:n)=x(1:n)
    fopt=f
    k1=0
    k2=0
    ksm=.false.
    kc=0
    knan=.false.
    hp=h
    if (constr) Reset=.false.

! 1-D SEARCH ----{
    do while (.true.)
      x1(1:n)=x(1:n)
      f1=f
      if (constr) then
        FsbPnt1=FsbPnt
        fp1=fp
      END IF

! NEW POINT
      x(1:n)=x(1:n) + hp*g0(1:n)
      ii=0
      do i=1,n
        if (ABS(x(i)-x1(i)) < ABS(x(i))*epsnorm) ii=ii + 1
      END DO

! FUNCTION VALUE
      call fun(x, f)
      options(10)=options(10) + one
      if (h1*f >= infty) then
         if (dispwarn) then
           WRITE(*, *) errmes
           WRITE(*, *) error5
         END IF
         options(9)=-seven
         GO TO 999
      END IF
      if (constr) then
         fp=f
         call func(x, fc)
         options(12)=options(12) + one
         if (ABS(fc) >= infty) then
             if (dispwarn) then
                WRITE(*, *) errmes
                WRITE(*, *) error52
                WRITE(*, *) error6
             END IF
             options(9)=-five
             GO TO 999
         END IF
         if (fc <= cnteps) then
            FsbPnt=.true.
            fc=zero
         else
            FsbPnt=.false.
            fp_rate=fp - fp1
            if (fp_rate < -epsnorm) then
              if (.not.FsbPnt1) then
                d=SUM( (x(1:n) - x1(1:n))**2 )
                d=SQRT(d)
                PenCoefNew=-1.5d1*fp_rate/d
                if (PenCoefNew > 1.2_dp*PenCoef) then
                  PenCoef=PenCoefNew
                  Reset=.true.
                  kless=0
                  f=f + PenCoef*fc
                  exit
                END IF
              END IF
            END IF
         END IF
         f=f + PenCoef*fc
      END IF

      if (ABS(f) >= infty) then
          if (dispwarn) then
            WRITE(*, *) wrnmes
            WRITE(*, *) error32
          END IF
          if (ksm.or.kc >= mxtc) then
             options(9)=-three
             GO TO 999
          else
             k2=k2 + 1
             k1=0
             hp=hp/dq
             x(1:n)=x1(1:n)
             f=f1
             knan=.true.
             if (constr) then
               FsbPnt=FsbPnt1
               fp=fp1
             END IF
          END IF

! STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
      ELSE IF (ii == n) then
         stepvanish=stepvanish + 1
         if (stepvanish >= 5) then
             options(9)=-ten - four
             if (dispwarn) then
                WRITE(*, *) termwarn1
                WRITE(*, *) endwarn4
             END IF
             GO TO 999
         else
             do i=1,n
              x(i)=x1(i)
             END DO
             f=f1
             hp=hp*ten
             ksm=.true.
             if (constr) then
                FsbPnt=FsbPnt1
                fp=fp1
             END IF
         END IF

! USE SMALLER STEP
      ELSE IF (h1*f < h1*gamma**INT(SIGN(one,f1))*f1) then
          if (ksm) exit
          k2=k2 + 1
          k1=0
          hp=hp/dq
          x(1:n)=x1(1:n)
          f=f1
          if (constr) then
             FsbPnt=FsbPnt1
             fp=fp1
          END IF
          if (kc >= mxtc) exit

! 1-D OPTIMIZER IS LEFT BEHIND
      else
          if (h1*f <= h1*f1) exit

! USE LARGER STEP
          k1=k1 + 1
          if (k2 > 0) kc=kc + 1
          k2=0
          if (k1 >= 20) then
              hp=du20*hp
          ELSE IF (k1 >= 10) then
              hp=du10*hp
          ELSE IF (k1 >= 3) then
              hp=du03*hp
          END IF
      END IF
    END DO
! ----}  End of 1-D search

! ADJUST THE TRIAL STEP SIZE ----{
    dx=SUM( (xopt(1:n) - x(1:n))**2 )
    dx=SQRT(dx)
    if (kg < kstore)  kg=kg+1
    if (kg >= 2)  then
       do i=kg,2,-1
         nsteps(i)=nsteps(i-1)
       END DO
    END IF
    d=SUM( g0(1:n)**2 )
    d=SQRT(d)
    nsteps(1)=dx/(ABS(h)*d)
    kk=zero
    d=zero
    do i=1,kg
       dd=dble(kg-i+1)
       d=d + dd
       kk=kk + nsteps(i)*dd
    END DO
    kk=kk/d
    if (kk > des) then
       if (kg == 1) then
          h=h*(kk-des+one)
       else
          h=h*SQRT(kk-des+one)
       END IF
    ELSE IF (kk < des) then
       h=h*SQRT(kk/des)
    END IF

    if (ksm) stepvanish=stepvanish + 1
! ----}

! COMPUTE THE GRADIENT ----{
        if (app) then
          do j=1,n
            if (g0(j) >= zero) then
               deltax(j)=h1*ddx
            else
               deltax(j)=-h1*ddx
            END IF
          END DO
          obj=.true.     
          if (constr) then
             call apprgrdn(n, g, x, fp, fun, deltax, obj)
          else 
             call apprgrdn(n, g, x, f, fun, deltax, obj)
          END IF
          options(10)=options(10) + n_float
        else
          call grad(x, g)
          options(11)=options(11) + one
        END IF
        ng=SUM( g(1:n)**2 )
        ng=SQRT(ng)
        if (ng >= infty) then
          if (dispwarn) then
            WRITE(*, *) errmes
            WRITE(*, *) error42
          END IF
          options(9)=-four
          GO TO 999
        ELSE IF (ng < ZeroGrad) then
          if (dispwarn) then
            WRITE(*, *) wrnmes
            WRITE(*, *) warn1
          END IF
          ng=ZeroGrad
        END IF

! Constraints:
        if (constr) then
          if (.not.FsbPnt) then
            if (ng < 1.d-2*PenCoef) then
               kless=kless+1
               if (kless >= 20) then
                  PenCoef=PenCoef/ten
                  Reset=.true.
                  kless=0
               END IF
            else
               kless=0
            END IF
            if (appconstr) then
               do j=1,n
                 if (x(j) >= zero) then
                    deltax(j)=ddx
                 else
                    deltax(j)=-ddx
                 END IF
               END DO
               obj=.false.
               call apprgrdn(n, gc, x, fc, func, deltax, obj)
               options(12)=options(12) + n_float
            else
               call gradc(x,gc)
               options(13)=options(13) + one
            END IF
            ngc=SUM( gc(1:n)**2 )
            ngc=SQRT(ngc)
            if (ngc >= infty) then
               if (dispwarn) then
                  WRITE(*, *) errmes
                  WRITE(*, *) error62
               END IF
               options(9)=-six
               GO TO 999
            ELSE IF (ngc < ZeroGrad .and. .not.appconstr) then
               if (dispwarn) then
                  WRITE(*, *) errmes
                  WRITE(*, *) error63
               END IF
               options(9)=-six
               GO TO 999
            END IF
            g(1:n)=g(1:n) + PenCoef*gc(1:n)
            ng=SUM( g(1:n)**2 )
            ng=SQRT(ng)
            if (Reset) then
               if (dispwarn) then
                  WRITE(*, *) wrnmes
                  WRITE(*, *) warn21
               END IF
               h=h1*dx/three
               k=k - 1
               nng=ng
               exit
            END IF
          END IF
        END IF
        if (h1*f > h1*frec) then
          frec=f 
          do i=1,n
            xrec(i)=x(i)
            grec(i)=g(i)
          END DO
        END IF
! ----}
       if (ng > ZeroGrad) then
         if (knorms < 10) knorms=knorms + 1
         if (knorms >= 2) then
           do i=knorms,2,-1
            gnorms(i)=gnorms(i-1)
           END DO
         END IF
         gnorms(1)=ng
         nng=one
           do i=1,knorms
             nng=nng*gnorms(i)
           END DO
         nng=nng**(one/dble(knorms))
       END IF

! Norm X:
       nx=SUM( x(1:n)**2 )
       nx=sqrt(nx)
  
! DISPLAY THE CURRENT VALUES ----{
       if (k == ld) then
         WRITE(*, *) 'Iteration # ..... Function Value ..... ',  &
                     'Step Value ..... Gradient Norm'
         WRITE(*, '(5x, i5, 7x, g13.5, 6x, g13.5, 7x, g13.5)') k, f, dx, ng
         ld=k + dispdata
       END IF
!----}

! CHECK THE STOPPING CRITERIA ----{
      termflag=.true.
      if (constr) then
        if (.not.FsbPnt) termflag=.false.
      END IF
      if(kcheck <= 5 .or. kcheck <= 12 .and. ng > one) termflag=.false.
      if(kc >= mxtc .or. knan) termflag=.false.

! ARGUMENT
      if (termflag) then
          ii=0
          stopping=.true.
          do i=1,n
            if (ABS(x(i)) >= lowxbound) then
               ii=ii + 1
               idx(ii)=i
               if (ABS(xopt(i)-x(i)) > options(2)*ABS(x(i))) then
                 stopping=.false.
               END IF
            END IF
          END DO
          if (ii == 0 .or. stopping) then
             stopping=.true.
             termx=termx + 1
             d=SUM( (x(1:n) - xrec(1:n))**2 )
             d=SQRT(d)
! FUNCTION
             if(ABS(f-frec) > detfr*ABS(f) .and.  &
                ABS(f-fopt) <= options(3)*ABS(f) .and.  &
                krerun <= 3 .and. .not. constr) then
                stopping=.false.
                if (ii > 0) then
                  do i=1,ii
                   j=idx(i)
                   if (ABS(xrec(j)-x(j)) > detxr*ABS(x(j))) then
                     stopping=.true.
                     exit
                   END IF
                  END DO
                END IF
                if (stopping) then
                   if (dispwarn) then
                     WRITE(*, *) wrnmes
                     WRITE(*, *) warn09
                   END IF
                   ng=zero
                   do i=1,n
                     x(i)=xrec(i)
                     g(i)=grec(i)
                     ng=ng + g(i)*g(i)
                   END DO
                   ng=SQRT(ng)
                   f=frec
                   krerun=krerun + 1
                   h=h1*MAX(dx, detxr*nx)/dble(krerun)
                   warnno=2
                   endwarn=endwarn2
                   exit
                else
                   h=h*ten
                END IF
             ELSE IF(ABS(f-frec) > options(3)*ABS(f) .and.  &
                     d < options(2)*nx .and. constr) then
                h=h*one
             ELSE IF (ABS(f-fopt) <= options(3)*ABS(f) .or.  &
                      ABS(f) <= lowfbound .or.  &
                      (ABS(f-fopt) <= options(3).and.  &
                      termx >= limxterm )) then
               if (stopf) then
                 if (dx <= laststep) then
                   if (warnno == 1 .and. ng < SQRT(options(3))) then
                      warnno=0
                   END IF
                   if (.not.app) then
                     do i=1,n
                      if (ABS(g(i)) <= epsnorm2) then
                        warnno=3
                        endwarn=endwarn3
                        exit
                      END IF
                     END DO
                   END IF
                   if (warnno /= 0) then
                      options(9)=-dble(warnno)-ten
                      if (dispwarn) then
                        WRITE(*, *) termwarn1
                        WRITE(*, *) endwarn
                        if (app) WRITE(*, *) appwarn
                      END IF
                   else
                      options(9)=dble(k)
                      if (dispwarn) WRITE(*, *) termwarn0
                   END IF
                   GO TO 999
                 END IF
               else
                stopf=.true.
               END IF
             ELSE IF (dx < powerm12*MAX(nx,one) .and. termx >= limxterm ) then
               options(9)=-four-ten
               if (dispwarn) then
                 WRITE(*, *) termwarn1
                 WRITE(*, *) endwarn4
                 if (app) WRITE(*, *) appwarn
                 f=frec
                 x(1:n)=xrec(1:n)
               END IF
               GO TO 999
             END IF
          END IF
      END IF

! ITERATIONS LIMIT
      if(k == iterlimit) then
          options(9)=-nine
          if (dispwarn) then
            WRITE(*, *) wrnmes
            WRITE(*, *) warn4
          END IF
          GO TO 999
      END IF
! ----}
! ZERO GRADIENT ----{
      if (constr) then
        if (ng <= ZeroGrad) then
            if (dispwarn) then
              WRITE(*, *) termwarn1
              WRITE(*, *) warn1
            END IF
            options(9)=-eight
            GO TO 999
        END IF
      else
        if (ng <= ZeroGrad) then
         nzero=nzero + 1
         if (dispwarn) then
           WRITE(*, *) wrnmes
           WRITE(*, *) warn1
         END IF
         if (nzero >= 3) then
           options(9)=-eight
           GO TO 999
         END IF
         g0(1:n)=-h*g0(1:n)/two
         do i=1,10
           x(1:n)=x(1:n) + g0(1:n)
           call fun(x, f)
           options(10)=options(10) + one
           if (ABS(f) >= infty) then
             if (dispwarn) then
               WRITE(*, *) errmes
               WRITE(*, *) error32
             END IF
             options(9)=-three
             GO TO 999
           END IF
           if (app) then
               do j=1,n
                 if (g0(j) >= zero) then
                    deltax(j)=h1*ddx
                 else
                    deltax(j)=-h1*ddx
                 END IF
               END DO
               obj=.true.
               call apprgrdn(n, g, x, f, fun, deltax, obj)
               options(10)=options(10) + n_float
           else
               call grad(x,g)
               options(11)=options(11)+one
           END IF
           ng=SUM( g(1:n)**2 )
           ng=SQRT(ng)
           if (ng >= infty) then
              if (dispwarn) then
                WRITE(*, *) errmes
                WRITE(*, *) error42
              END IF
              options(9)=-four
              GO TO 999
           END IF
           if (ng > ZeroGrad) exit
         END DO
         if (ng <= ZeroGrad) then
            if (dispwarn) then
              WRITE(*, *) termwarn1
              WRITE(*, *) warn1
            END IF
            options(9)=-eight
            GO TO 999
         END IF
         h=h1*dx
         exit
        END IF
      END IF
! ----}
! FUNCTION IS FLAT AT THE POINT ----{
      if (.not.constr .and. ABS(f-fopt) < ABS(fopt)*options(3) .and.  &
          kcheck > 5 .and. ng < one ) then

        ni=0
        do i=1,n
          if (ABS(g(i)) <= epsnorm2) then
            ni=ni + 1
            idx(ni)=i
          END IF
        END DO
        if (ni >= 1 .and. ni <= n/2 .and. kflat <= 3) then
          kflat=kflat+1
          if (dispwarn) then
             WRITE(*, *) wrnmes
             WRITE(*, *) warn31
          END IF
          warnno=1
          endwarn=endwarn1
          x1(1:n)=x(1:n)
          fm=f
          do i=1,ni
            j=idx(i)
            f2=fm
            y=x(j)
            if (y == zero) then
              x1(j)=one
            ELSE IF (ABS(y) < one) then
              x1(j)=sign(one, y)
            else
              x1(j)=y
            END IF
            do ip=1,20
              x1(j)=x1(j)/1.15_dp
              call fun(x1, f1)
              options(10)=options(10)+one
              if (ABS(f1) < infty) then
                if (h1*f1 > h1*fm) then
                  y=x1(j)
                  fm=f1
                ELSE IF (h1*f2 > h1*f1) then
                  exit
                ELSE IF (f2 == f1) then
                  x1(j)=x1(j)/1.5_dp
                END IF
                f2=f1
              END IF
            END DO
            x1(j)=y
          END DO
          if (h1*fm > h1*f) then
            if (app) then
              deltax(1:n)=h1*ddx
              obj=.true.
              call apprgrdn(n, gt, x1, fm, fun, deltax, obj)
              options(10)=options(10) + n_float
            else
              call grad(x1, gt)
              options(11)=options(11) + one
            END IF
            ngt=SUM( gt(1:n)**2 )
            if (ngt > epsnorm2 .and. ngt < infty) then
              if (dispwarn) WRITE(*, *) warn32
              do i=1,n
               x(i)=x1(i)
               g(i)=gt(i)
              END DO
              ng=ngt
              f=fm
              h=h1*dx/three
              options(3)=options(3)/five
              exit
            END IF   !! regular gradient
          END IF   !! a better value has been found
        END IF   !! function is flat
      END IF   !! pre-conditions are fulfilled
! ----}
  END DO   !! iterations
END DO   !! restart

! deallocate working arrays:

999 deallocate (idx, deltax, xx, grec, xrec, xopt, x1, z, gc, gt, g1, g0, g, B)

RETURN
end subroutine solvopt



SUBROUTINE apprgrdn(n, g, x, f, fun, deltax, obj)
 
! Subroutine APPRGRDN performs the finite difference approximation
! of the gradient <g> at a point <x>.
! f      is the calculated function value at a point <x>,
! <fun>  is the name of a subroutine that calculates function values,
! deltax is an array of the relative stepsizes.
! obj    is the flag indicating whether the gradient of the objective
!        function (1) or the constraint function (0) is to be calculated.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: g(:)
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: f
REAL (dp), INTENT(IN)      :: deltax(:)
LOGICAL, INTENT(IN)        :: obj

INTERFACE
  SUBROUTINE fun(x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
  END SUBROUTINE fun
END INTERFACE

! Local variables

REAL (dp)  :: d, y, fi
INTEGER    :: i, j
LOGICAL    :: center
REAL (dp), PARAMETER  :: lowbndobj = 2.d-10, lowbndcnt = 5.d-15,   &
                         one = 1.0_dp, ten = 10.0_dp, half = 0.5_dp

DO i=1,n
  y=x(i)
  d=MAX(lowbndcnt, ABS(y))
  d=deltax(i)*d
  IF (obj) THEN
    IF (ABS(d) < lowbndobj) THEN
      d=lowbndobj*SIGN(one, deltax(i))
      center=.true.
    ELSE
      center=.false.
    END IF
  ELSE
    IF (ABS(d) < lowbndcnt) THEN
      d=lowbndcnt*SIGN(one, deltax(i))
    END IF
  END IF
  x(i)=y + d
  CALL fun(x, fi)
  IF (obj) THEN
    IF (fi == f) THEN
      DO j=1,3
        d=d*ten
        x(i)=y + d
        CALL fun(x, fi)
        IF (fi /= f) EXIT
      END DO
    END IF
  END IF
  g(i)=(fi - f)/d
  IF (obj) THEN
    IF (center) THEN
      x(i)=y - d
      CALL fun(x, fi)
      g(i)=half*(g(i) + (f-fi)/d)
    END IF
  END IF
  x(i)=y
END DO

RETURN
END SUBROUTINE apprgrdn



SUBROUTINE soptions(default)
 
! SOPTIONS returns the default values for the optional parameters
! used by SolvOpt.

REAL (dp), INTENT(OUT)  :: default(13)

default(1) = -1.0_dp
default(2) = 1.d-4
default(3) = 1.d-6
default(4) = 15.d3
default(5) = 0.0_dp
default(6) = 1.d-8
default(7) = 2.5D0
default(8) = 1.d-11
default(9) = 0.0_dp
default(10) = 0.0_dp
default(11) = 0.0_dp
default(12) = 0.0_dp
default(13) = 0.0_dp

RETURN
END SUBROUTINE soptions

END MODULE Shor_minimization

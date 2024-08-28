## Classification of [Alan Miller's Fortran codes](https://jblevins.org/mirror/amiller/)

### Linear Algebra
[as6.f90](https://jblevins.org/mirror/amiller/as6.f90) Cholesky factorization of symmetric positive definite matrix with only lower triangle stored as a vector (2 versions).

[as7.f90](https://jblevins.org/mirror/amiller/as7.f90) Inversion of symmetric positive definite matrix with only lower triangle stored as a vector using AS6 (2 versions). Includes short test program and interfaces.

[as60.f90](https://jblevins.org/mirror/amiller/as60.f90) Calculates the eigenvalues/vectors of a real symmetric matrix.

### Linear Regression
[as110.f90](https://jblevins.org/mirror/amiller/as110.f90) Lp-norm fitting of straight lines, particularly for 1 <= p <= 2, by extension of an algorithm of Schlossmacher.

[as132.f90](https://jblevins.org/mirror/amiller/as132.f90) Minimize the sum of absolute deviations for a simple Y on X linear regression.

[as132.f90](https://jblevins.org/mirror/amiller/as132.f90) Minimize the sum of absolute deviations for a simple Y on X linear regression.

[as135.f90](https://jblevins.org/mirror/amiller/as135.f90) Min-Max (L-infinity) estimates for linear multiple regression, includes a test program.

[as282.f90](https://jblevins.org/mirror/amiller/as282.f90) Robust regression using the least median of squares (LMS) criterion. There is also a test program: [t_as282.f90](https://jblevins.org/mirror/amiller/t_as282.f90)

### Optimization
[as298.f90](https://jblevins.org/mirror/amiller/as298.f90) Hybrid minimization routine using simulated annealing and a user-supplied local minimizer. You should also look at Dr. Tibor Csendes' global optimization package.

[as319.f90](https://jblevins.org/mirror/amiller/as319.f90) Variable metric unconstrained function minimization without derivatives.

[cobyla.f90](https://jblevins.org/mirror/amiller/cobyla.f90) Mike Powell's routine for minimization of non-linear functions with smooth non-linear constraints, using local linear approximations to the constraints.

[conmin.zip](https://jblevins.org/mirror/amiller/conmin.zip) The classic CONMIN package for constrained minimization updated to Fortran 90. Test examples and the manual are included in the ZIP file. N.B. CONMIN is included as just one of the algorithms in TOMS algorithm 734 which can be downloaded from [netlib](http://www.netlib.org).

[ga.zip](https://jblevins.org/mirror/amiller/ga.zip) A package for global minimization using genetic algorithms. Please note the copyright conditions if you want to use this for commercial purposes.

[global.f90](https://jblevins.org/mirror/amiller/global.f90) At Arnold Neumaier's [site](https://arnold-neumaier.at/glopt.html), this is recommended as the most successful of the global optimization packages. There is a sample program fit.f90 and the original documentation global.txt for the f77 version. I have included testfunc.f90 which will eventually contain all of Neumaier's 30 test functions. N.B. Users of local optimization packages usually obtain satisfactory convergence after 10s or sometimes 100s of function evaluations. Global optimization routines usually require many 1000s of function evaluations.

[minim.f90](https://jblevins.org/mirror/amiller/minim.f90) The Nelder-Mead simplex algorithm for unconstrained minimization. It does NOT require or use derivatives. N.B. This is NOT for linear programming! [t_minim.f90](https://jblevins.org/mirror/amiller/t_minim.f90) is a very simple test program for minim.f90 which may help users. This is now ELF90-compatible.

[pikaia.zip](https://jblevins.org/mirror/amiller/pikaia.zip) Another genetic algorithm for global optimization without derivatives. This one comes from the High Altitude Observatory of NCAR. This is a ZIP file, containing a number of interesting examples. There is an excellent manual, but you will have to download it from their web site. N.B. The manual is 120 pages long, and in postscript. Where does the name Pikaia come from? Pikaia gracilens is a flattened worm-like beast about 5cm long which crawled in the mud on the sea floor about 530 million years ago!

[solvopt.zip](https://jblevins.org/mirror/amiller/solvopt.zip) The SolvOpt package minimizes or maximizes nonlinear functions, which may be nonsmooth and may have constraints, using the so-called method of exact penalization. This is a zip file containing several driver programs.

[tensolve.zip](https://jblevins.org/mirror/amiller/tensolve.zip) This contains both TOMS 739 (UNCMIN) for unconstrained minimization, and TOMS 768 (TENSOLVE) for the solution of sets of non-linear equations.

[tn.zip](https://jblevins.org/mirror/amiller/tn.zip) Stephen Nash's truncated-Newton code for the minimization of continuous functions. It can use differences instead of derivatives, and bounds may be imposed on the parameters.

[toms667.f90](https://jblevins.org/mirror/amiller/toms667.f90) Global minimization using a stochastic algorithm, with 37 test problems [test667.f90](https://jblevins.org/mirror/amiller/test667.f90). Many of the test problems require hundreds of thousands of function evaluations. See also global.f90 below in the miscellaneous code section.

[toms778.zip](https://jblevins.org/mirror/amiller/toms778.zip) This subroutine solves bound-constrained optimization problems by using the compact formula of the limited memory BFGS updates.

[toms813.zip](https://jblevins.org/mirror/amiller/toms813.zip) This package is for local minimization, with or without convex constraints. It requires the user to supply first derivatives. This version is in standard Fortran. It uses the so-called spectral projected gradient (SPG) method.

[tron.zip](https://jblevins.org/mirror/amiller/tron.zip) Newton's method for large bound-constrained optimization problems by Chih-Jen Lin & Jorge More', a MINPACK-2 project.

[uobyqa.f90](https://jblevins.org/mirror/amiller/uobyqa.f90) Mike Powell's package for unconstrained minimization when derivatives are not available. There is a test/driver program: [t_uobyqa.f90](https://jblevins.org/mirror/amiller/t_uobyqa.f90) and a file of results from the test program: [test.out](https://jblevins.org/mirror/amiller/test.out) The documentation, which is in gzipped postscript, can be downloaded from the [Optimization Decision Tree](https://plato.asu.edu/guide.html).

### Uniform Random Number Generation
[luxury.f90](https://jblevins.org/mirror/amiller/luxury.f90) Another generator of uniformly distributed random numbers. [luxtst.f90](https://jblevins.org/mirror/amiller/luxtst.f90) A program to test luxury.f90.

[taus88.f90](https://jblevins.org/mirror/amiller/taus88.f90) L'Ecuyer's 1996 Tausworthe random number generator, and [lfsr113.f90](https://jblevins.org/mirror/amiller/lfsr113.f90) L'Ecuyer's 1999 Tausworthe random number generator. The first has a cycle of 2^88 while the second is a little slower but has a cycle of 2^113. Both are translations from C. N.B. These both assume that default integers are 32-bit.

[lfsr258.f90](https://jblevins.org/mirror/amiller/lfsr258.f90) A 64-bit random number generator from Pierre L'Ecuyer with a cycle of about 2^258 or more than 10^77.

[dprand.f90](https://jblevins.org/mirror/amiller/dprand.f90) Nick Maclaren's double precision random number generator, translated into ELF90-compatible form.

[mt19937.f90](https://jblevins.org/mirror/amiller/mt19937.f90) The 'Mersenne Twister' random number generator from Japan with a cycle of length (2^19937 - 1). [mt19937a.f90](https://jblevins.org/mirror/amiller/mt19937a.f90) is a version for compilers which stop when there are integer overflows, as some do when compiler check options are enabled for debugging purposes. Before using this module, you are advised to scan the [Japanese web site](https://www.math.keio.ac.jp/~matumoto/emt.html) particularly under `Initialization'. mt19937.f90 was revised on 5 February 2002; mt19937a.f90 has not been revised.

[rand3.f90](https://jblevins.org/mirror/amiller/rand3.f90) Yet another random number generator; this one is based upon an algorithm by Donald Knuth (1997).

[freq2d.f90](https://jblevins.org/mirror/amiller/freq2d.f90) Pairs of uniform random numbers should be uniformly distributed in the unit square. This is a simple test program which failed using a new (1998) compiler. ELF90-compatible.

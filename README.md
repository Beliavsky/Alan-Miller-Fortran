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

[conmin.zip](https://jblevins.org/mirror/amiller/conmin.zip) The classic CONMIN package for constrained minimization updated to Fortran 90. Test examples and the manual are included in the ZIP file. N.B. CONMIN is included as just one of the algorithms in TOMS algorithm 734 which can be downloaded from [netlib](http://www.netlib.org).

[ga.zip](https://jblevins.org/mirror/amiller/ga.zip) A package for global minimization using genetic algorithms. Please note the copyright conditions if you want to use this for commercial purposes.

[global.f90](https://jblevins.org/mirror/amiller/global.f90) At Arnold Neumaier's [site](https://arnold-neumaier.at/glopt.html), this is recommended as the most successful of the global optimization packages. There is a sample program fit.f90 and the original documentation global.txt for the f77 version. I have included testfunc.f90 which will eventually contain all of Neumaier's 30 test functions. N.B. Users of local optimization packages usually obtain satisfactory convergence after 10s or sometimes 100s of function evaluations. Global optimization routines usually require many 1000s of function evaluations.

[minim.f90](https://jblevins.org/mirror/amiller/minim.f90) The Nelder-Mead simplex algorithm for unconstrained minimization. It does NOT require or use derivatives. N.B. This is NOT for linear programming! [t_minim.f90](https://jblevins.org/mirror/amiller/t_minim.f90) is a very simple test program for minim.f90 which may help users. This is now ELF90-compatible.

[pikaia.zip](https://jblevins.org/mirror/amiller/pikaia.zip) Another genetic algorithm for global optimization without derivatives. This one comes from the High Altitude Observatory of NCAR. This is a ZIP file, containing a number of interesting examples. There is an excellent manual, but you will have to download it from their web site. N.B. The manual is 120 pages long, and in postscript. Where does the name Pikaia come from? Pikaia gracilens is a flattened worm-like beast about 5cm long which crawled in the mud on the sea floor about 530 million years ago!

[solvopt.zip](https://jblevins.org/mirror/amiller/solvopt.zip) The SolvOpt package minimizes or maximizes nonlinear functions, which may be nonsmooth and may have constraints, using the so-called method of exact penalization. This is a zip file containing several driver programs.

[tensolve.zip](https://jblevins.org/mirror/amiller/tensolve.zip) This contains both TOMS 739 (UNCMIN) for unconstrained minimization, and TOMS 768 (TENSOLVE) for the solution of sets of non-linear equations.

[toms778.zip](https://jblevins.org/mirror/amiller/toms778.zip) This subroutine solves bound-constrained optimization problems by using the compact formula of the limited memory BFGS updates.

[tron.zip](https://jblevins.org/mirror/amiller/tron.zip) Newton's method for large bound-constrained optimization problems by Chih-Jen Lin & Jorge More', a MINPACK-2 project.

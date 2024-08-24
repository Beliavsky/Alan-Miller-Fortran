## Classification of [Alan Miller's Fortran codes](https://jblevins.org/mirror/amiller/)

### Optimization
[as298.f90](https://jblevins.org/mirror/amiller/as298.f90) Hybrid minimization routine using simulated annealing and a user-supplied local minimizer. You should also look at Dr. Tibor Csendes' global optimization package.

[conmin.zip](https://jblevins.org/mirror/amiller/conmin.zip) The classic CONMIN package for constrained minimization updated to Fortran 90. Test examples and the manual are included in the ZIP file. N.B. CONMIN is included as just one of the algorithms in TOMS algorithm 734 which can be downloaded from [netlib](http://www.netlib.org).

[global.f90](https://jblevins.org/mirror/amiller/global.f90) At Arnold Neumaier's [site](https://arnold-neumaier.at/glopt.html), this is recommended as the most successful of the global optimization packages. There is a sample program fit.f90 and the original documentation global.txt for the f77 version. I have included testfunc.f90 which will eventually contain all of Neumaier's 30 test functions. N.B. Users of local optimization packages usually obtain satisfactory convergence after 10s or sometimes 100s of function evaluations. Global optimization routines usually require many 1000s of function evaluations.

[minim.f90](https://jblevins.org/mirror/amiller/minim.f90) The Nelder-Mead simplex algorithm for unconstrained minimization. It does NOT require or use derivatives. N.B. This is NOT for linear programming! [t_minim.f90](https://jblevins.org/mirror/amiller/t_minim.f90) is a very simple test program for minim.f90 which may help users. This is now ELF90-compatible.

[tensolve.zip](https://jblevins.org/mirror/amiller/tensolve.zip) This contains both TOMS 739 (UNCMIN) for unconstrained minimization, and TOMS 768 (TENSOLVE) for the solution of sets of non-linear equations.

[toms778.zip](https://jblevins.org/mirror/amiller/toms778.zip) This subroutine solves bound-constrained optimization problems by using the compact formula of the limited memory BFGS updates.

[tron.zip](https://jblevins.org/mirror/amiller/tron.zip) Newton's method for large bound-constrained optimization problems by Chih-Jen Lin & Jorge More', a MINPACK-2 project.

## Classification of [Alan Miller's Fortran codes](https://jblevins.org/mirror/amiller/)

### Optimization
[as298.f90](https://jblevins.org/mirror/amiller/as298.f90) Hybrid minimization routine using simulated annealing and a user-supplied local minimizer. You should also look at Dr. Tibor Csendes' global optimization package.

[conmin.zip](https://jblevins.org/mirror/amiller/conmin.zip) The classic CONMIN package for constrained minimization updated to Fortran 90. Test examples and the manual are included in the ZIP file. N.B. CONMIN is included as just one of the algorithms in TOMS algorithm 734 which can be downloaded from netlib (http://www.netlib.org).

[tensolve.zip](https://jblevins.org/mirror/amiller/tensolve.zip) This contains both TOMS 739 (UNCMIN) for unconstrained minimization, and TOMS 768 (TENSOLVE) for the solution of sets of non-linear equations.

[toms778.zip](https://jblevins.org/mirror/amiller/toms778.zip) This subroutine solves bound-constrained optimization problems by using the compact formula of the limited memory BFGS updates.

[tron.zip](https://jblevins.org/mirror/amiller/tron.zip) Newton's method for large bound-constrained optimization problems by Chih-Jen Lin & Jorge More', a MINPACK-2 project.

2D Poisson Solver

This code provides a MATLAB implementation of a 2D Poisson solver using the multigrid method. The solver can be used to solve the Poisson equation of the form:

∇^2 u = f

where u is the solution, f is a given function, and ∇^2 is the Laplace operator.

Files
The following files are included in this repository:
V_cycle.m: Implements the V-cycle for Multigrid solver of the 2D poisson equation.
 
prolong.m: A helper function that prolongs a coarse grid solution to a fine grid.

restrict.m: A helper function that restricts a fine grid residual to a coarse grid.

Jacrex.m: A helper function that performs an under weighted Jacobi smoother on a given grid together with the function that
computes the residuals for any grid.

test_script.m: A script that demonstrates the usage of  function.

Usage
To solve a 2D Poisson equation, follow these steps:

Define the grid size and the number of levels to use in the multigrid method.

Define the right-hand side function f and the boundary conditions.

Call the V_cycle.m function, using the test_script function passing in the grid size, the number of levels, f, and the boundary conditions.

The V_cycle function returns the solution u after a certain threshold is reached.
See test_script.m for an example. Several right hand sides are included in this script.

An exact solution function  (exact_solution.m) is included to compare the approximate solution with this exact.

Two scripts, fd2poissonsp.m and testFdPoisson.m are used to find a solution to the same plot using Gaussian elimination. This is meant to give us a dimension for comparison purposes.



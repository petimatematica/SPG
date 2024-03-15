# spg.jl 

This file contains the implementation of the Spectral Projected Gradient Method (SPG). Additionally, it contains the GLL and BMR non-monotone linear searches.

## Function spg:
Implements the Spectral Projected Gradient Method (SPG). It takes the following parameters:

- ** x0 (Vector): The initial point.
- ** f (Function): The objective function to be minimized.
- ** gradf (Function): The gradient of the objective function.
- ** proj (Function): The projection function.
- ** tol (Float64): The convergence tolerance.
- ** maxiter (Int): The maximum number of iterations allowed.
- ** lambda_min (Float64): The minimum value for the lambda parameter.
- ** lambda_max (Float64): The maximum value for the lambda parameter.
- ** M (Int): The value of M.
- ** sigma1 (Float64): The value of sigma1.
- ** sigma2 (Float64): The value of sigma2.
- ** gamma (Float64): The value of gamma.
- ** linesearch (Function): The line search function.
The function returns the minimization point, an error indicator, optimization information, elapsed time, and the sequence of points during optimization.

Functions spg1 and spg2:
Implement backtracking routines for the SPG method. They take the same parameters as the spg function, except for the linesearch parameter. Each function returns the next point, the gradient at that point, the direction, a history of function values, the step size, elapsed time, and the number of objective function evaluations.

Function spg1:
This function implements the backtracking routine for the SPG method with a specific update rule for the step size. It adjusts the step size based on a condition involving the function values at the current and the next points.

Function spg2:
Similar to spg1, this function also implements the backtracking routine for the SPG method but with a different update rule for the step size. It calculates the step size based on the projected gradient and a condition involving the function values.

These functions are essential components of the SPG method and are used by the main spg function to perform optimization.

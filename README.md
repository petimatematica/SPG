# spg.jl 

This file contains the implementation of the Spectral Projected Gradient Method (SPG). Additionally, it contains the GLL and BMR nonmonotone linear searches.

## Function spg:
Implements the Spectral Projected Gradient Method (SPG). It takes the following parameters:

- x0 (Vector): The initial point.
- f (Function): The objective function to be minimized.
- gradf (Function): The gradient of the objective function.
- proj (Function): The projection function.
- tol (Float64): The convergence tolerance.
- maxiter (Int): The maximum number of iterations allowed.
- lambda_min (Float64): The minimum value for the lambda parameter.
- lambda_max (Float64): The maximum value for the lambda parameter.
- M (Int): Integer number greater than zero.
- sigma1 (Float64): The value of sigma1.
- sigma2 (Float64): The value of sigma2.
- gamma (Float64): The value of gamma.
- linesearch (Function): The line search function.
The function returns the minimization point, an error indicator, optimization information, elapsed time, and the sequence of points during optimization.

Functions spg1 and spg2:
Implement backtracking routines for the SPG method. They take the same parameters as the spg function, except for the linesearch parameter. Each function returns the next point, the gradient at that point, the direction, a history of function values, the step size, elapsed time, and the number of objective function evaluations.

## Linesearch of Grippo, Lampariello and Lucidi (GLL) for SPG1

### Description:
The GLL backtracking routine is responsible for selecting the step size during each iteration of the SPG optimization algorithm. It iteratively adjusts the step size until it satisfies certain conditions based on the objective function and gradient information.

### Parameters:
- k (Int): Iteration index.
- lambda_k (Float64): Current step size.
- x_k (Vector): Current estimate of the minimizer.
- gradf_x_k (Vector): Gradient of the objective function evaluated at the current estimate.
- f_hist (Vector{Float64}): History of objective function values.
- M (Int): Parameter determining the history size for backtracking.
- sigma1 (Float64): Parameter used in the backtracking process.
- sigma2 (Float64): Parameter used in the backtracking process.
- gamma (Float64): Parameter used in the Armijo condition.

### Output:
- x_plus (Vector): Updated estimate of the minimizer.
- gradf_x_kp1 (Vector): Gradient of the objective function evaluated at the updated estimate.
- s_k (Vector): Step taken in the optimization direction.
- y_k (Vector): Change in gradient.
- f_hist (Vector{Float64}): Updated history of objective function values.
- alpha (Float64): Selected step size.
- et (Float64): Elapsed time.
- evalf (Int): Number of objective function evaluations.

### Remarks:
The GLL backtracking routine adjusts the step size iteratively until it satisfies conditions based on the Armijo condition and other parameters.

It employs a line search strategy to find an appropriate step size that ensures sufficient decrease in the objective function value.

GLL is one of the variants of the backtracking routine used in the SPG optimization algorithm, offering different strategies for selecting the step size during optimization.

## Linesearch of Birgin, Martínez and Raydan (BMR) for SPG2

### Description:
The BMR backtracking routine is responsible for determining the appropriate step size during each iteration of the SPG optimization algorithm. It adjusts the step size iteratively based on the objective function and gradient information until certain conditions are met.

Parameters:
- k (Int): Iteration index.
- lambda_k (Float64): Current step size.
- x_k (Vector): Current estimate of the minimizer.
- gradf_x_k (Vector): Gradient of the objective function evaluated at the current estimate.
- f_hist (Vector{Float64}): History of objective function values.
- M (Int): Parameter determining the history size for backtracking.
- sigma1 (Float64): Parameter used in the backtracking process.
- sigma2 (Float64): Parameter used in the backtracking process.
- gamma (Float64): Parameter used in the Armijo condition.

### Output:
- x_plus (Vector): Updated estimate of the minimizer.
- gradf_x_kp1 (Vector): Gradient of the objective function evaluated at the updated estimate.
- s_k (Vector): Step taken in the optimization direction.
- y_k (Vector): Change in gradient.
- f_hist (Vector{Float64}): Updated history of objective function values.
- alpha (Float64): Selected step size.
- et (Float64): Elapsed time.
- evalf (Int): Number of objective function evaluations.

### Remarks:
The BMR backtracking routine adjusts the step size iteratively until it satisfies conditions based on the Armijo condition and other parameters.

It employs a line search strategy to find an appropriate step size that ensures sufficient decrease in the objective function value.

BMR is another variant of the backtracking routine used in the SPG optimization algorithm, providing an alternative strategy for selecting the step size during optimization.

# main.jl 

## Imports and Inclusions:
The script imports necessary packages such as CUTEst for accessing optimization problems and includes essential Julia files (spg.jl).

## Configuration and Parameters:
Parameters such as problem names, dimensions, and solver settings are initialized. The script iterates over two cases of SPG (SPG1 and SPG2).

## Optimization Cycle:
For each problem, the script initializes the problem using CUTEst, defines the objective function and its gradient, and sets upper and lower bounds. The SPG method is then called with the specified parameters, and optimization is performed. Information about the optimization process is saved to a JLD2 file.

## Performance Measurement:
The number of iterations, CPU time, and function evaluations are recorded for each problem and SPG variant. If an error occurs during optimization, corresponding values are set to infinity.

## Performance Profiles:
Performance profiles are generated based on the collected data. Profiles illustrate the percentage of solved problems as a function of the number of iterations, CPU time ratio, and function evaluations. Profiles are saved as images (performanceprofileiters.png, performanceprofiletime.png, performanceprofileevalf.png).

## Instructions for Use:
Ensure that the necessary Julia files (spg.jl) are in the same directory as the script.
Run the script using the Julia interpreter.
After execution, performance profiles will be generated and saved as images in the current directory.

## Remark:
The script provides insights into the performance of the SPG Method across a range of optimization problems.

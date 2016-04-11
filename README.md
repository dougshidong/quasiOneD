# Quasi-One Dimensional Euler Flow Channel Adjoint-Based Optimization

## Flow Solver

The flow solver is based on a quasi-one dimensional channel.

Finite-volume, explicit time-stepping.

Flux Schemes: Steger-Warming (SW), Modified SW (MSW), Corrected-MSW, Roe

Time-stepping Schemes: Euler Explicit, Runge-Kutta 4th order, Jameson's modifified RK

## Optimization

Gradient-based unconstrained optimization.

Design variables are either the parameters of a sine-shaped channel or the discrete area's themselves (work in progress).

Opimization methods: Steepest-Descent, Newton's method, Quasi-Newton (BFGS)

## Gradient

Finite-difference, adjoint variable method or direct method is used to obtain the gradients.

Finite-difference (first order) requires nDes flow solutions. Approximate derivatives

Direct method requires nDes flow solutions. Exact derivatives.

Adjoint Variable requires 1 flow solution. Exact derivatives

## Hessian

Finite-difference, direct-direct, adjoint-direct, adjoint-adjoit, direct-adjoint. 

Refer to "Direct, adjoint and mixed approaches for the computation of Hessian in airfoil design problems" by Papadimitriou and Giannakoglou.

## External Libraries

Eigen is used for almost all the Linear Algebra computations.

PETSc has recently been implemented to solve for AX = B on a the "petsc" branch.

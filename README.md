# finate-element-method-project

This R script implements a numerical solver for a second-order differential equation using a finite element method. And Galerkin method The solver generates a solution for the differential equation on the domain [a,b] with Dirichlet boundary conditions and plots the resulting function. The primary goal of this solver is to approximate solutions to equations of the form G⋅4πρ.

## Functions Overview
calcXi(h, i)
 - Computes the grid point $ξ_i$ for index i with spacing h.

calcEi(x, h, i)
 - Calculates the value of the basis function $E_i(x)$ for a given x, with grid point spacing h and index i. This function is defined piecewise and returns 0 if x is out of the bounds of $ξ_{i−1}$ and $ξ_{i+1}$.

calcEiPrime(x, h, i)
 - Computes the derivative of the basis function $E'_i(x)$ with respect to x, using finite differences.

calcF(x, G)
- A placeholder function intended to calculate the G⋅4πρ(x) term. 
ρ(x) is indicator funcion on bounds (1,2]

integrateL(a, b, h, i)
- Performs numerical integration of the $E_i(x)$ function using the Gauss quadrature method over the range [a,b] for the index i.

integrateB(a, b, h, i, j)
- Integrates the B matrix 

calcB(h, i, j, a, b)
- Calculates the element $B_{ij}$ of the matrix B by defining appropriate integration bounds for i and j.

calcL(G, h, i)
- Constructs the vector L by iterating through all basis functions and calculating the corresponding values for each index.

createBMatrix(a, b, n)
- Constructs the matrix B for the finite element method by calculating $B_{ij}$ for all i,j, while applying boundary conditions.

calcUroof(a, b, u1, u2, x)
- Computes a linear function that satisfies the Dirichlet boundary conditions. This is used to adjust the solution u by adding the boundary function to the numerical solution.

differentialSolver(n, a, b, u1, u2, G)
- The main function that orchestrates the solution process. It:
  - Defines the grid points,
  - Calculates the vector L and matrix 𝐵
  - Solves the system 𝐵 x 𝑤 = 𝐿
  - Combines the solution with the boundary condition correction,
  - Plots the final function u(x).

The function plots the numerical solution Φ(x) over the specified interval [a,b].

## Usage
To use this solver:
- Define boundaries (a, b),
- define n (number of basis functions)
- define boundary values $u_1$ and $u_2$
- define constant G.
- call differentialSolver

## Potential changes:
As this was suposed to be project that just plots differential equation of gravitational potential, a lot can be changed.

For example `calcF(x, G)` is the easiest change and fix, the boundaries should've been passed as parameters of the function.

`calcF` could be also differed for other function that solves:

$\frac{d^2Φ}{dx^2} = f(x)$

the end.

fin.

Why are you still reading this :v
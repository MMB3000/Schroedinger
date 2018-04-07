--------readme.txt--------

Author: Miguel Bengala
Description: Solves the Schroedinger equation and returns wave function as an output




How to compile this project:
Linux, using g++, with Makefile included in the project




Files included in the project:

MatAnyD.cpp     MatAnyD.h
--> Define a Matrix class and the corresponding operations

GaussianElimination.cpp
--> Implement Gaussian elimination a its subproducts (determinant, inverse matrix, solution of linear system)

Eigen.cpp
--> Implements algorithm to calculate the first eigenvector and eigenvalue and next ones by the deflation method

SchroedingerSolver.cpp     SchroedingerSolver.h
--> Outputs to file the solutions of the Schroedinger Equation for given inputs of physical quantities, and prints the correponding energies

Main.cpp
--> Receives the input of parameters (physical quantities, etc.) needed to solve the Schroedinger equation, calls the solver for the equation
and outputs the execution time

WaveFunction.nb
--> Plots the graphs for the wave function solutions




Format of input:
The input should be given as "InputSchroedinger.txt" and it should contain these values in the following order:
|Parameter k|
|Parameter a|
|Parameter b|
|Planck constant in the desired physical units|
|Parameter m|
|Size of x coordinate|
|Number of points|

Example of input values:
1
1
1
1
1
100
99



Format of output of wave function values:
The output consists of three columns: the first corresponds to the index of the solution, the second to the value of x, and the third to the value of Psi of the given x.





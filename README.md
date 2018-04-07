# Schroedinger

Author: Miguel Bengala  
Description: Solves the Time Independent Schroedinger Equation of the form  
![schroedinger_eq](http://latex.codecogs.com/gif.latex?-%5Cfrac%7B%28h/2%5Cpi%29%5E2%7D%7B2m%7D%5Cfrac%7Bd%5E2%7D%7Bdx%5E2%7D%5CPsi%28x%29&plus;%5Cfrac%7B1%7D%7B2%7DKx%5E2%5CPsi%28x%29%3DE%5CPsi%28x%29)  
with  
![kinetic_eq](http://latex.codecogs.com/gif.latex?K%3D-%5Cfrac%7B%5Chbar%5E2%7D%7B2m%7D%5Cfrac%7Bd%5E2%7D%7Bdx%5E2%7D)   and   ![potential_eq](http://latex.codecogs.com/gif.latex?V%28x%29%3D%5Cfrac%7Bk%7D%7B2%7D%5Cleft%20%28%20%5Cfrac%7Ba%5E4%7D%7Bb%5E2%7D-%5Cfrac%7Ba%5E4%7D%7Bx%5E2&plus;b%5E2%7D%20%5Cright%20%29)  
Returns wave function of the different solutions as an output to a file and corresponding energies to the screen.    




## How to compile and run this project:  
On Linux and using g++, it can be compiled with the Makefile included in the project and then run on the terminal  ```./Schroedinger```    



## Files included in the project:

Main.cpp    
--> Receives the input of parameters (physical quantities, etc.) needed to solve the Schroedinger equation, calls the solver for the equation
and outputs the execution time

SchroedingerSolver.cpp  
--> Outputs to file the solutions of the Schroedinger Equation for given inputs of physical quantities, and prints the correponding energies

SchroedingerSolver.h  
--> Header to the methods of SchroedingerSolver

MatAnyD.cpp  
--> Define a Matrix class and the corresponding operations

MatAnyD.h  
--> Header including the Matrix class, its friend methods (including gaussian elimination and methods for eigenvectors) and other related functions.

GaussianElimination.cpp  
--> Implement Gaussian elimination a its subproducts (determinant, inverse matrix, solution of linear system)

Eigen.cpp  
--> Implements algorithm to calculate the first eigenvector and eigenvalue and next ones by the deflation method

WaveFunction.nb  
--> Plots the graphs for the wave function solutions    




## Format of input:  
The input should be given as file named "InputSchroedinger.txt" and it should contain these values, each on a different line, in the following order:  
```
< Parameter k >  
< Parameter a >  
< Parameter b >  
< Planck constant in the desired physical units >  
< Parameter m >  
< Size of x coordinate >  
< Number of points >
```

### Example of input values:  
```
1  
1  
1  
1  
1  
100  
99    
```


## Format of output of wave function values:  
The output consists of three columns: the first corresponds to the index of the solution, the second to the value of x, and the third to the value of Psi of the given x.





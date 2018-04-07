#ifndef SCHROEDINGERSOLVER_H
#define SCHROEDINGERSOLVER_H

//---------- SchroedingerSolver.h ----------
#include <iostream>

//EnergyMatrix receives all the properties from Matrix base, but is limited to receive values in tridiagonal positions
class EnergyMatrix: public Matrix {
public:
  EnergyMatrix(unsigned _Nrows,unsigned _Ncolumns); //standard constructor of EnergyMatrix
  
  void set(unsigned i, unsigned j, double x);

  EnergyMatrix& operator=( const Matrix & r2 ); 
  //redefines operator = to allow assignments using derived and base class

};

//computes the Hamiltonian matrix from the sum of kinetic and potential energy
//receives as arguments the physical quantities and algorithm parameters obtained from input file or calculated in main()
EnergyMatrix computeHamilton(double Nmax, double hx, double kValue, double aValue, double bValue, double Lbd, double hbar, double m);

//solves schroedinger equation and outputs the wave function to txt file
//receives as arguments the Hamiltonian and the position parameters
void solveSchroedinger(Matrix & H, double Lbd, double h);


#endif

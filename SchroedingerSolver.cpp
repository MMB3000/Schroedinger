//---------- SchroedingerSolver.cpp ----------
//Solves the Schroedinger equation, including the computation of the hamiltonian and the algorithm to
//solve the equation from eigenvectors; defines a derived class EnergyMatrix from base class Matrix
//which restricts Matrix to the tridiagonal form
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include "MatAnyD.h"
#include "SchroedingerSolver.h"


//standard constructor of EnergyMatrix
EnergyMatrix::EnergyMatrix(unsigned _Nrows,unsigned _Ncolumns) : Matrix(_Nrows, _Ncolumns) {
  //besides the standard construction of Matrix, EnergyMatrix is initialized with '0' in all positions
  for(unsigned i=0;i<Nrows;i++) {
    for(unsigned j=0;j<Ncolumns;j++) {
      ptr[i][j]=0;
    }
  }
  
}

void EnergyMatrix::set(unsigned i, unsigned j, double x) {
  //the 'if' prohibits from writing outside the matrix and in non tridiagonal positions
  if((i==j||i-1==j||j-1==i) && (i>=0) && (j>=0) && (i<Nrows) && (j<Ncolumns))
    ptr[i][j] = x;
  else
    std::cout << "Impossible to write in position (" << i << "," << j << ") of energy matrix!" << std::endl;
}

//redefines operator = to allow assignments using derived and base class
EnergyMatrix& EnergyMatrix::operator=( const Matrix & r2 ){

  //check if sizes are compatible, if not change dimensions of the first matrix
  if(Nrows != r2.size_rows() || Ncolumns != r2.size_columns()) {
    std::cout << "size does not match!" << std::endl;

    //delete the pointers with incorrect size
    for(unsigned i=0;i<Nrows;i++) delete[] ptr[i];
    delete[] ptr;

    //allocate new pointers with the desired size
    Nrows = r2.size_rows();
    Ncolumns = r2.size_columns();
    ptr = new double*[Nrows];
    for(unsigned i=0;i<Nrows;i++) ptr[i]= new double[Ncolumns];
  }

  //assign the values
  for(unsigned i = 0; i < Nrows; ++i) {
    for(unsigned j = 0; j < Ncolumns; ++j) {
      ptr[i][j] = r2.get(i,j);
    }
  }
  
  return *this;
}


//computes the Hamiltonian matrix from the sum of kinetic and potential energy
//receives as arguments the physical quantities and algorithm parameters obtained from input file or calculated in main()
EnergyMatrix computeHamilton(double Nmax, double hx, double kValue, double aValue, double bValue, double Lbd, double hbar, double m) {
  
  EnergyMatrix K(Nmax, Nmax); //Kinetic energy
  EnergyMatrix V(Nmax, Nmax); //Potential energy
  EnergyMatrix H(Nmax, Nmax); //Hamiltonian
  
  double Vn, x;
  
  //assigns values to K and V matrices
  for(int i=0; i<Nmax; ++i) {
    x=-Lbd + hx + i*hx;
    //value to be assigned to the diagonal of the potential matrix
    Vn=(kValue/2)*((pow(aValue,4)/pow(bValue,2))-(pow(aValue,4)/(pow(x,2)+pow(bValue,2))));

    //assign values to potential and kinetic matrices
    for(int j=0; j<Nmax; ++j)	{
      if(i==j){
	K.set(i,j, -2);
	V.set(i,j, Vn);
	if(i!=Nmax-1) {
	  K.set(i+1, j, 1);
	  K.set(i, j+1, 1);
	}
      }
    }
  }

  //normalization of matrix K to the desired units
  K=(1/(hx*hx))*K;
  K=(-hbar*hbar/(2*m))*K;

  //final computation of the hamiltonian (matrix H)
  H=V+K;
  
  return H;
}

//solves schroedinger equation and outputs the wave function to txt file
//receives as arguments the Hamiltonian and the position parameters
void solveSchroedinger(Matrix & H, double Lbd, double h){
  
  std::ofstream output("OutputWaves.txt");
  if ( output.fail() ) {
    std::cerr << "Impossible to open file OutputWaves.txt" <<  std::endl;
    exit(1); // this terminates the process
  }
  
  Matrix M(H);
  
  M=Invert(H); //first step of the algorithm to solve equation

  
  int nv=4; //number of desired vectors
  
  int Nx=H.size_rows();
  
  Matrix psi(nv,Nx);
  Matrix energia(nv,1);
  
  Matrix aproxi(Nx, 1);
  
  
  for(int i=0; i<Nx; ++i){
    aproxi.set(i,0,1);
  }
  
  Matrix eigenTable(eigenvec_deflate(M,aproxi,nv)); //computation of the nv eigenvectors
  //(organized in a table - each row corresponding to a solution and first element eigenvalue, then eigenvector)
  
  energia=eigenTable.getcol(0); //gets energy values<=>1/eigenvalues

  for(int i=0; i<nv; ++i)
    {
      energia.set(i,0, 1/energia.get(i,0));
    }

  std::cout << "\nThe values for energy of the solutions are, respectively: " << std::endl;
  std::cout << energia << std::endl; //prints the values for energy in the screen

  //assign the values for psi (correspond to the eigenvectors determined)
  for(unsigned i=0;i<nv;i++) {
    for(unsigned j=0;j<Nx;j++) {
      psi.set(i,j,eigenTable.get(i,j+1));
    }
  }
  
  double x; //position
  
  for(int i=0; i<nv; ++i) {
    for(int j=0; j<Nx; ++j) {
      x=-Lbd + j*h; //position calculated from interval and step
      output << i << " " << x << " " << psi.get(i,j) <<std::endl; //output of psi function to a file
    }
  }
  
  output.close();
}

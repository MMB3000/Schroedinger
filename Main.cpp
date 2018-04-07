//---------- Main.cpp ----------
//--> Receives the input of parameters (physical quantities, etc.)
//needed to solve the Schroedinger equation, calls the solver for the equation
//and outputs the execution time
//
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "MatAnyD.h"
#include "SchroedingerSolver.h"


int main()
{
  //used to measure time of execution
  clock_t t;
  t= clock();


  //Introduction to appear on screen
  std::cout << "**************************************************" << std::endl;
  std::cout << "*                                                *" << std::endl;
  std::cout << "*  Author: Miguel Bengala                        *" << std::endl;
  std::cout << "*                                                *" << std::endl;
  std::cout << "*  Description: Solves the Time Independent      *" << std::endl;
  std::cout << "*  Schroedinger Equation and returns wave        *" << std::endl;
  std::cout << "*  function of the different solutions as an     *" << std::endl;
  std::cout << "*  output and prints corresponding energies      *" << std::endl;      
  std::cout << "*                                                *" << std::endl;
  std::cout << "*  For more information: README.md               *" << std::endl;
  std::cout << "*                                                *" << std::endl;
  std::cout << "**************************************************" << std::endl;
  

  
  double kValue, aValue, bValue, hbar, m, Lbd; //declaration of physical quantities
  int Nmax; //declaration of algorithm related variables
  double hx; //step


  
  //INPUT
  std::ifstream input("InputSchroedinger.txt");
  //check for fail opening the file
  if (input.fail()) {
    std::cerr << "\nImpossible to read file InputSchroedinger.txt\n" << std::endl;
    exit(1); // this terminates the process
  }

  //receive the parameters from the file
  input >>kValue >>aValue >>bValue >>hbar >>m >>Lbd >>Nmax;
  
  input.close();



  //CALCULATIONS (with output)
  hx=2*Lbd/(Nmax+1); //calculate step

  Matrix H(Nmax, Nmax); //Hamiltonian
  
  H = computeHamilton(Nmax,hx,kValue,aValue,bValue,Lbd,hbar,m); //computes hamiltonian needed to solve Sch-Eq
  
  solveSchroedinger(H,Lbd,hx); //solves the equation and outputs psi-function to file and energy to the screen
  

  
  //time of execution
  t=clock()-t;
  std::cout << "It took " << t << " clicks (" << ((float)t)/CLOCKS_PER_SEC << " seconds)" << std::endl;
 
  
  return 0;
}

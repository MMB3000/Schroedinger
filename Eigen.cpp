//---------- Eigen.cpp ----------
//Implements algorithm to calculate the first eigenvector and eigenvalue
//and next ones by the deflation method
#include <iostream>
#include <cmath>
#include "MatAnyD.h"


//namespace to save eigenvalue
namespace eigen {
  double eigenval;
}


//Computes eigenvector of a matrix
Matrix eigenvec( Matrix & M, const Matrix & V) {
	
  double _ni, _nj;
  _ni = M.size_rows();
  _nj = M.size_columns();
  
  eigen::eigenval=1;
  
  double mod0=1;
  
  Matrix res0(V);
  Matrix res1(V);
  
  double A0=0.5, A1=144;
  double mod=1.;
  double vpr0=0.5, vpr1=1;
  //iterate to calculate eigenvector until 1e-12 precision is obtained
  while(std::abs(vpr0-vpr1)>1e-12) {
    
    res1=M*res0;
    
    A1=res1.get(0,0);
    for(int j=0; j<_ni; ++j){
      if(std::abs(A1)<std::abs(res1.get(j,0))) {
	A1=res1.get(j, 0);
	A0=res0.get(j,0);
      }
    }
    
    vpr0=vpr1;
    vpr1 = A1/A0;
    mod = sqrt(dot(res1,res1));
    res1=(1/mod)*res1;
    res0=res1;
  }
  
  eigen::eigenval=vpr1;
  return res0;
}


//returns eigenvalue calculated with method eigenVec()
double eigenvalue(Matrix & A,Matrix & V1) {
  Matrix V(V1);
  Matrix valorfunc(V);
  valorfunc=eigenvec(A,V);
  double eigenva=1.;
  eigenva=eigen::eigenval;
  return eigenva;
}


//computes n next eigenvalues and -vectors using Deflation method
Matrix eigenvec_deflate(Matrix & A,Matrix & V1,unsigned n) {
  Matrix a(V1);
  Matrix K(A);
  unsigned a_Nrows= a.size_rows();
  unsigned tab_Nrows=a_Nrows+1;
  Matrix Eigentable(n,tab_Nrows);
  for(unsigned i=0;i<n;i++) {
    a=eigenvec(K,a);
    
    Eigentable.set(i,0,eigen::eigenval);
    
    for(unsigned j=1;j<tab_Nrows;j++) {
      Eigentable.set(i,j,a.get(j-1,0));
    }
    
    K=K-eigen::eigenval*a*transpose(a);
  }
  
  return Eigentable;
}

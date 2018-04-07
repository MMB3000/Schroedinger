//---------- GaussianElimination.cpp ----------
//Implement Gaussian elimination a its subproducts
//(determinant, inverse matrix, solution of linear system)
#include <iostream>
#include <cmath>
#include "MatAnyD.h"


//namespace used to calculate determinant
namespace determin{
  int coef_det;
}

//applies Gaussian elimination to system A|b, with these matrices as arguments
Matrix applyGaussianElimin(Matrix & A, Matrix & b) {
  unsigned A_Nrows= A.size_rows();
  unsigned A_Ncolumns=A.size_columns();
  unsigned M_Ncolumns=A.size_columns()+b.size_columns();
  unsigned M_Nrows=A.size_rows();
  
  Matrix M(M_Nrows,M_Ncolumns); //augmented matrix of the system

  for(unsigned i=0;i<A_Nrows;i++) {
    for(unsigned j=0;j<A_Ncolumns;j++) {
      //assign the values of the first columns of M to the values of A
      M.set(i,j,A.get(i,j));
    }
    
    for(unsigned j=A_Ncolumns;j<M_Ncolumns;j++) {
      //and the following columns to the values of b
      M.set(i,j,b.get(i,j-A_Ncolumns));
    }
  }
  
  determin::coef_det=1; //coefficient for deteminant formula
  double *temp; //pointer used to switch position of greater pivot
  unsigned maior=0; //index of greater pivot
  double maior_value=0; //value of greater pivot
  for(unsigned k=0;k<A_Nrows;k++) {
    
    //searches for the greater pivot
    maior=k;
    maior_value=A.get(k,k);
    for(unsigned i=k;i<A_Nrows;i++) {
      if(std::abs(A[i][k])>maior_value) {
	maior=i;
	maior_value=A.get(i,k);
      }
    }
    if(((k-maior)%2)!=0) determin::coef_det*=-1; //coefficient changes signal when the line change is odd
      //line change
      temp=M.ptr[k];
      M.ptr[k]=M.ptr[maior];
      M.ptr[maior]=temp;
      
      
      //subtraction of lines
      double pivot;
      for(unsigned g=k+1;g<M_Nrows;g++) {
	pivot=M.get(g,k);
	for(unsigned h=0;h<M_Ncolumns;h++) {
	  M.set(g,h,M.get(g,h)-(pivot/M.get(k,k))*M.get(k,h));
	}
      }      
  }
  
  return M;
}

//applies Gauss-Jordan method
Matrix applyGaussJordan(Matrix & A, Matrix & b) {
  Matrix M(applyGaussianElimin(A,b));
  unsigned A_Nrows= A.size_rows();
  unsigned A_Ncolumns=A.size_columns();
  unsigned M_Ncolumns=A.size_columns()+b.size_columns();
  unsigned M_Nrows=A.size_rows();
  
  //subtraction of lines bottom to top
  double pivot,pivot_diag;
  for(int k=M_Nrows-1;k>0;--k) //goes up line by line, from the last to the second
    {
      pivot_diag=M.get(k,k);
      for(int g=k-1;g>=0;--g) { //goes up line by line, from the above the pivot we are eliminating to the first
	pivot=M.get(g,k);
	for(int h=0;h<M_Ncolumns;h++) { //runs the line from left to right
	  M.set(g,h,M.get(g,h)-(pivot/pivot_diag)*M.get(k,h));
	}
      }      
    }
  
  //pivot normalization
  double result;
  for(unsigned i=0;i<M_Nrows;i++) { //moves line by line from first to the last
      pivot_diag=M.get(i,i);
      //moves column by column (runs the lines) of the augmented matrix, from the first to the last
      for(unsigned j=0;j<M_Ncolumns;j++) {
	result=(1./pivot_diag)*M.get(i,j);
	if(result==-0) {
	  M.set(i,j,0);
	}
	else {
	  M.set(i,j,result);
	}
      }
  }
  
  return M;
}

//solves the system Ax=b
Matrix solveGauss(Matrix & A, Matrix & b) {
  Matrix M(applyGaussJordan(A,b)); //Augmented matrix, transformed by Gauss-Jordan Method
  
  unsigned A_Ncolumns=A.size_columns();
  unsigned b_Nrows= b.size_rows();
  unsigned b_Ncolumns=b.size_columns();
  unsigned M_Ncolumns=A.size_columns()+b.size_columns();
  unsigned M_Nrows=A.size_rows();
  
  Matrix x(b_Nrows,b_Ncolumns); //Solution

  for(unsigned i=0;i<M_Nrows;i++) {
    for(unsigned j=A_Ncolumns;j<M_Ncolumns;j++) {
      x.set(i,j-A_Ncolumns,M.get(i,j));
    }
  }

  return x;
}


//inverts matrix
Matrix invert(Matrix & A) {
  unsigned A_Nrows= A.size_rows();
  Matrix identity(Id(A_Nrows));
  Matrix M(solveGauss(A,identity)); // solves system A|Identity --> corresponds to inverse matrix
  return M;
}

//computes determinant of matrix
double det(Matrix & A) {
  double deter=1;
  Matrix nula(0,0);
  
  Matrix M(applyGaussianElimin(A,nula)); //applies Gaussian elimination to matrix A

  deter*=determin::coef_det;
  unsigned A_Nrows= A.size_rows();

  //multiplies all the non normalized pivots of A by the coefficients of line change
  for(unsigned i=0;i<A.size_rows();i++) deter*=M.get(i,i);
  
  return deter;
}

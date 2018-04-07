//---------- MatAnyD.cpp ----------
//--> Define a Matrix class and the corresponding operations
#include <iostream>
#include <cmath>
#include "MatAnyD.h"


//constructor with both dimensions of matrix
Matrix::Matrix(unsigned _Nrows,unsigned _Ncolumns) : Nrows(_Nrows) , Ncolumns(_Ncolumns) {
  ptr = new double*[Nrows];
  for(unsigned i=0;i<Nrows;i++) ptr[i]= new double[Ncolumns];
}

//copy constructor
Matrix::Matrix(const Matrix& mat) : Nrows(mat.Nrows) , Ncolumns(mat.Ncolumns) {
  ptr = new double*[Nrows];

  for(unsigned i=0;i<Nrows;i++) {
    ptr[i]= new double[Ncolumns];
    for(unsigned j=0;j<Ncolumns;j++) {
      ptr[i][j]=mat.ptr[i][j];
    }
  }
}

//destructor
Matrix::~Matrix() {
  for(unsigned i=0;i<Nrows;i++) delete[] ptr[i];
  delete[] ptr;
}


unsigned Matrix::size_rows() const { return Nrows; }
unsigned Matrix::size_columns() const { return Ncolumns; }


double Matrix::get(unsigned i,unsigned j) const { return ptr[i][j]; }

void Matrix::set(unsigned i, unsigned j,double x) {
  //the 'if' prohibits from writing outside the matrix
  if((i>=0) && (j>=0) && (i<Nrows) && (j<Ncolumns))
    ptr[i][j] = x;
  else
    std::cout << "Impossible to write in position (" << i << "," << j << ") of matrix!" << std::endl;
} 


//matrices output
std::ostream& operator<<( std::ostream& out, const Matrix & r ) {
  for( unsigned i = 0; i < r.size_rows(); ++i ) {
    for(unsigned j=0;j<r.size_columns();++j) out  << r.get(i,j) << "\t";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  return out;
}

//matrices assignment
Matrix& Matrix::operator=( const Matrix & r2 ){

  //check if sizes are compatible, if not change dimensions of the first matrix
  if(Nrows != r2.size_rows() || Ncolumns != r2.size_columns()) {
    std::cout << "size does not match!" << std::endl;
    for(unsigned i=0;i<Nrows;i++) delete[] ptr[i];
    delete[] ptr;
    Nrows = r2.size_rows();
    Ncolumns = r2.size_columns();
    ptr = new double*[Nrows];
    for(unsigned i=0;i<Nrows;i++) ptr[i]= new double[Ncolumns];
  }

  //assign the values
  for(unsigned i = 0; i < Nrows; ++i)
    for(unsigned j = 0; j < Ncolumns; ++j) ptr[i][j] = r2.ptr[i][j];
  return *this;
}

//compound assignment (sum)
Matrix& Matrix::operator+=( const Matrix & r2 ) {
  
  //check if sizes are compatible, if not change dimensions of the first matrix
  if(Nrows != r2.size_rows() || Ncolumns != r2.size_columns()) {
    std::cout << "size does not match!" << std::endl;
    for(unsigned i=0;i<Nrows;i++) delete[] ptr[i];
    delete[] ptr;
    Nrows = r2.size_rows();
    Ncolumns = r2.size_columns();
    ptr = new double*[Nrows];
    for(unsigned i=0;i<Nrows;i++) {
      ptr[i]= new double[Ncolumns];
      for(unsigned j;j<Ncolumns;j++) ptr[i][j]=0; //initialize resized matrix
    }
  }
  
  for(unsigned i = 0; i < Nrows; ++i)
    for(unsigned j = 0; j < Ncolumns; ++j)
      ptr[i][j] = ptr[i][j]+r2.ptr[i][j];

  return *this;
}

//logical test ==
bool Matrix::operator==(const Matrix & r2 ) const {
  bool check=true;
  for(int i=0;i<Nrows;++i) {
    for(unsigned j=0;j<Ncolumns;++j) {
      if(ptr[i][j]!=r2.ptr[i][j]) {
	check=false;
	break;
      }
    }
  }
  
  return check;
}

//sum of two matrices
Matrix Matrix::operator+(const Matrix & r2 ) const {
  unsigned _Nrows= size_rows();
  unsigned _Ncolumns=size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i ) {
      for(unsigned j=0;j<_Ncolumns;++j) {
	ret.set(i,j, get(i,j) + r2.get(i,j));
      }
  }
  return ret;
}

//adds a number to every component of matrix (number to the right)
Matrix Matrix::operator+(const double a ) const { 
  unsigned _Nrows = size_rows();
  unsigned _Ncolumns= size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i )
    for(unsigned j=0;j<_Ncolumns;++j) {
      ret.set(i,j, get(i,j) + a);
    }
  return ret;
}

//adds a number to every component of matrix (number to the left)
Matrix operator+(const double a, const Matrix & r2 ) {
  unsigned _Nrows = r2.size_rows();
  unsigned _Ncolumns= r2.size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i ) {
    for(unsigned j=0;j<_Ncolumns;++j) {
      ret.set(i,j, a+ r2.get(i,j));
    }
  }
  return ret;
}

//difference between two matrices
Matrix Matrix::operator-(const Matrix & r2 ) const {
  unsigned _Nrows= size_rows();
  unsigned _Ncolumns=size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i ) {
    for(unsigned j=0;j<_Ncolumns;++j) {
      ret.set(i,j, get(i,j) - r2.get(i,j));
    }
  }
  return ret;
}


//matrix with symmetrical values
Matrix Matrix::operator-() const {
  unsigned _Nrows = size_rows();
  unsigned _Ncolumns=size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i ) {
    for(unsigned j=0;j<_Ncolumns;++j) {
      ret.set(i,j, -get(i,j));
    }
  }
  return ret;
}

double const* Matrix::operator[](size_t i) const { return &ptr[i][0];}
double* Matrix::operator[](size_t i) { return &ptr[i][0];}

//dot product (only works for column matrices, i.e., vectors)
double dot( const Matrix & r1, const Matrix & r2 ) {
  
  //check if the matrices aren't vectors and have different size
  if(r1.size_rows() != r2.size_rows() || r1.size_columns()!= 1 || r2.size_columns()!= 1) {
    std::cout << "size does not match!" << std::endl;
    return 0;
  }
  else { //if not returns dot product
    unsigned _Nrows = r1.size_rows();
    double res = 0;
    for( unsigned i = 0; i < _Nrows; ++i )
      res +=  r1.get(i,0) * r2.get(i,0);
    return res;
  }
}

//product of matrices
Matrix Matrix::operator*(const Matrix & r2) {
  unsigned _Nrows= size_rows();
  unsigned _Ncolumns=r2.size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  double sum=0.;
  if(size_columns()!=r2.size_rows()) { //check matrix dimensions
    std::cout << "It is not possible to multiply these matrices!!" << std::endl;
    Matrix ret2(r2);
    return ret2; //returns the right matrix in case of error
  }
  else {
    for(int i=0;i<_Nrows;i++) {
      for(int j=0;j<_Ncolumns;j++) {
	for(int k=0;k<size_columns();k++) {
	  sum+=get(i,k)*r2.get(k,j);
	}
	ret.set(i,j, sum);
	sum=0.;
      }
    }
    return ret;
  }
}

//product of number and matrix (number to the left)
Matrix operator*(const double a, const Matrix & r2 ) { 
  unsigned _Nrows = r2.size_rows();
  unsigned _Ncolumns= r2.size_columns();
  Matrix ret(_Nrows,_Ncolumns);
  for( unsigned i = 0; i < _Nrows; ++i ) {
    for(unsigned j=0;j<_Ncolumns;++j)	{
      ret.set(i,j, a*r2.get(i,j));
    }
  }
  return ret;
}

//transpose matrix
Matrix transpose(Matrix & A) {
  unsigned _Nrows = A.size_rows();
  unsigned _Ncolumns = A.size_columns();
  Matrix At(_Ncolumns,_Nrows);
  for(unsigned i=0;i<_Nrows;i++) {
    for(unsigned j=0;j<_Ncolumns;j++) {
      At.set(j,i,A.get(i,j));
    }
  }
  return At;
}

//returns the vector corresponding to a given column of matrix
Matrix Matrix::getcol(int j){
  Matrix Vec(size_rows(), 1);
  
  for(int i=0; i<size_rows(); ++i){
    Vec.set(i,0, get(i,j));
  }
  return Vec;
}

//returns the vector corresponding to a given row of matrix
Matrix Matrix::getlin(int j){
  Matrix Vec(1, size_columns());
  
  for(int i=0; i<size_columns(); ++i){
    Vec.set(0,i, get(j,i));
  }
  return Vec;
}


//generates identity matrix of given size
Matrix Id(int n) {
  Matrix identity(n,n);
  for(unsigned i=0;i<n;i++) {
    for(unsigned j=0;j<n;j++) {
      if(i==j) {
	identity.set(i,j,1);
      }
      else {
	identity.set(i,j,0);
      }
    }
  }
  return identity;
}

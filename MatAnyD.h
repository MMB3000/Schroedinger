#ifndef MATANYD_H
#define MATANYD_H

//---------- MatAnyD.h ----------
#include <iostream>

class Matrix {
 protected:
  unsigned Nrows,Ncolumns;
  double **ptr;
  int rows, cols;
  
 public:
  explicit Matrix(unsigned _Nrows,unsigned _Nlines); //standard constructor
  Matrix(const Matrix& mat); //copy constructor
  ~Matrix(); //Destructor
  unsigned size_rows() const; //gives the number of columns of the matrix
  unsigned size_columns() const; //gives the number of rows of the matrix
  double get(unsigned i,unsigned j) const; //gets element i,j
  void set(unsigned i,unsigned j,double x); //sets element i,j to x
  
  
  Matrix& operator=( const Matrix & r2 );  //matrices assignment
  
  Matrix operator+(const Matrix & r2 ) const; //sum of two matrices
  Matrix operator+(const double a ) const; //adds a number to every component of matrix (number to the right)
  friend Matrix operator+(const double a, const Matrix & r2 ); //adds a number to every component of matrix (number to the left)
  
  
  Matrix operator-(const Matrix & r2 ) const; //difference between two matrices
  Matrix operator-() const; //matrix with symmetrical values
  
  double const* operator[](size_t i) const;
  double* operator[](size_t i);
  
  Matrix& operator+=( const Matrix & r2 ) ; //compound assignment (sum)
  
  bool operator==(const Matrix & r2 ) const; //logical test ==
  
  
  friend std::ostream& operator<<( std::ostream& out, const Matrix & r ); //output operator
  
  friend double dot( const Matrix & r1, const Matrix & r2 ); // dot product
  
  Matrix operator*(const Matrix & r2); //product of matrices

  friend Matrix operator*(const double a, const Matrix & r2 ); // product of number by matrix (number to the left)

  friend Matrix transpose(Matrix & A); //transpose matrix
  
  friend Matrix MEG(Matrix & A, Matrix & b); //applies Gaussian elimination to system A|b, with these matrices as arguments
  
  friend Matrix MEGJ(Matrix & A, Matrix & b); //applies Gauss-Jordan method
  
  friend Matrix MEG_solve(Matrix & A, Matrix & b); //solves the system Ax=b
  
  friend Matrix Invert(Matrix & A); //inverts matrix
  
  friend double det(Matrix & A); //computes determinant of matrix

  friend Matrix eigenvec(Matrix & A,const Matrix & V1); //Computes eigenvector of a matrix
  
  friend double eigenvalue(Matrix & A,Matrix & V1); //returns eigenvalue calculated with method eigenVec()
  
  friend Matrix eigenvec_deflate(Matrix & A,Matrix & V1,unsigned n); //computes n next eigenvalues and -vectors using Deflation method

  Matrix getcol(int j); //returns the vector corresponding to a given column of matrix

  Matrix getlin(int j); //returns the vector corresponding to a given row of matrix


};

Matrix Id(int n); //generates identity matrix of given size

#endif

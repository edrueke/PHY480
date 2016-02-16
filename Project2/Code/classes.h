/*
Elizabeth Drueke
PHY 480
March 4, 2016

Project 2

This is the header file for the matrix and vector classes defined for 
project 2 (much the same as project 1).
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

//Vector class
class thevec{
 public:
  int sz; //Size of vector
  double *point; //Dynamic memory array

  //Constructors
  thevec(int s); //Initializer
  thevec(vector<double> vec); //Initializer from a vector of doubles
  thevec(const thevec &vec); //Copy Constructor
  ~thevec(); //Deconstructor

  //Other functions and operators
  string print(); //Print
  friend thevec operator+(const thevec &vec1,const thevec &vec2); //Addition
  friend thevec operator-(const thevec &vec1,const thevec &vec2); //Subtraction
  friend double operator*(const thevec &vec1,const thevec &vec2); //Dot product
  friend thevec operator*(const thevec &vec1,double fact); //Scalar Multiplication
  double operator[](int i); //Component retrieval
  thevec &operator=(const thevec &vec); //Equality
};

thevec operator+(const thevec &vec1,const thevec &vec2);
thevec operator-(const thevec &vec1,const thevec &vec2); 
double operator*(const thevec &vec1,const thevec &vec2); 
thevec operator*(const thevec &vec1,double fact);

//Matrix class
class themat{
 public:
  int sz; //Size of matrix (sz x sz)
  double **point; //Dynamic memory array

  //Constructors
  themat(int s); //Initializer
  themat(vector<double> vec); //Initializer from a vector of doubles
  themat(const themat &mat); //Copy constructor
  ~themat(); //Deconstructor

  //Other functions and operations
  string print(); //Print
  friend themat operator+(const themat &mat1, const themat &mat2); //Addition
  friend themat operator-(const themat &mat1, const themat &mat2); //Subtraction
  friend themat operator*(const themat &mat1, const themat &mat2); //Matrix-matrix multiplication
  friend thevec operator*(const themat &mat, const thevec &vec); //Matrix-vector multiplication
  friend themat operator*(const themat &mat, double fact); //Scalar Multiplication
  thevec operator[](int i); //Component retrieval
  themat &operator=(const themat &mat); //Equality

  friend class thevec;
};

themat operator+(const themat &mat1, const themat &mat2); //Addition
themat operator-(const themat &mat1, const themat &mat2); //Subtraction
themat operator*(const themat &mat1, const themat &mat2); //Matrix-matrix multiplication
thevec operator*(const themat &mat, const thevec &vec); //Matrix-vector multiplication
themat operator*(const themat &mat, double fact); //Scalar multiplication

//Additional functions - Gaussian elimination
vector<vector<double> > one_forw_reduc(themat matr, thevec vec, int i);
vector<vector<double> > one_back_reduc(themat matr, thevec vec, int i);
thevec gauss_elim(themat matr1,thevec vec);

//Additional functions - LU Decomposition

vector<vector<double> > LU_decomp(themat matr);
thevec LU_decomp_solver(themat &L,themat &U,thevec &vec);

//Additional functions - Convert to string

string to_string(double);
double error_calc(double,double);

//Project 1 - Specific Functions

double function(double &x);
double sol_function(double &x);
vector<vector<double> > LU_decomp_special(int sz);
thevec LU_decomp_solver_special(thevec &vec, themat &L, themat &U);

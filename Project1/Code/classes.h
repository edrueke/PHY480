/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This is the header file for the matrix and vector classes defined for 
project 1.
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
//
//#ifndef __THEVEC_H_INCLUDED__
//#define __THEVEC_H_INCLUDED__
//
//#ifndef __THEMAT_H_INCLUDED__
//#define __THEMAT_H_INCLUDED__
//
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
  double operator[](int i); //Component retrieval
  thevec &operator=(const thevec &vec); //Equality
  //scalar multiplication?
};

thevec operator+(const thevec &vec1,const thevec &vec2);
thevec operator-(const thevec &vec1,const thevec &vec2); 
double operator*(const thevec &vec1,const thevec &vec2); 

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
  thevec operator[](int i); //Component retrieval
  themat &operator=(const themat &mat); //Equality
  //scalar multiplication?

  friend class thevec;
};

themat operator+(const themat &mat1, const themat &mat2); //Addition
themat operator-(const themat &mat1, const themat &mat2); //Subtraction
themat operator*(const themat &mat1, const themat &mat2); //Matrix-matrix multiplication
thevec operator*(const themat &mat, const thevec &vec); //Matrix-vector multiplication

//Additional functions - Gaussian elimination
vector<vector<double> > one_forw_reduc(themat matr, thevec vec);
vector<vector<double> > one_back_reduc(themat matr, thevec vec);
thevec gauss_elim(themat matr1,thevec vec);

//Additional functions - LU Decomposition

vector<vector<double> > LU_decomp(themat matr);
 
//Additional functions - Convert to string

string to_string(double);

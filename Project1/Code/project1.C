/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This script solves the 1d Poisson Equation -u''(x)=f(x) for 0<x<1 subject to
the Dirichlet boundary conditions u(0)=u(1)=0 for f(x)=100e**(-10*x).
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const int dim = 2;

//double first_deriv(double &x, double &h){
//  /*Takes in 
//    and returns the first derivative*/
//  
//  return 0;
//}
//
//double second_deriv(double &x, double &h){
//  /*Takes in
//    and returns the second derivative.*/
//
//  return 0;
//}

void project1(){
  /*The main function of the script.*/

  double *vec;
  vec = new double[dim];

  double **matr;
  matr = new double*[dim];

  for(int i=0;i<dim;i++){
    matr[i] = new double[dim];
  }

  for(int i=0;i<dim;i++){
    vec[i]=1;
  }
  matr[0][0]=1; matr[0][1]=1; matr[1][0]=2; matr[1][1]=-1;

  PrintMatr(matr);
  PrintVector(vec);

  double** matr_product = mat_mat_mult(matr,matr);
  PrintMatr(matr_product);
  double* vec_product = mat_vec_mult(matr,vec);
  PrintVector(vec_product);


  vector<vector<double> > LU_decomposed = LU_decomp(matr);


  double** L = MakeMatr(LU_decomposed.at(0));
  double** U = MakeMatr(LU_decomposed.at(1));

  cout<<"L:"<<endl;
  PrintMatr(L);
  cout<<"U:"<<endl;
  PrintMatr(U);
  cout<<"Gauss elim:"<<endl;
  PrintVector(gauss_elim(matr,vec));

  for(int i=0;i<dim;i++){
    delete [] matr[i];
    delete [] L[i];
    delete [] U[i];
    delete [] matr_product[i];
  }

  delete[] matr;
  delete [] U;
  delete [] L;
  delete [] vec;
  delete [] matr_product;
  delete [] vec_product;
}

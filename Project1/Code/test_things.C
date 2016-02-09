/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This script tests the various functions and classes and sets benchmarks for 
reproducibility.
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
#include "classes.C"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const int dim = 2;

void test_things(){
  /*
    The main function which does the tests.
  */

  //Define the vector and matrix to be used.
  thevec vec = thevec(dim);
  themat matr = themat(dim);

  //Fill the matrix and vector with values.
  for(int i=0;i<dim;i++){
    vec.point[i]=1;
  }
  matr.point[0][0]=1; matr.point[0][1]=1; 
  matr.point[1][0]=2; matr.point[1][1]=-1;

  //Print the matrix and vector.
  cout<<"Matrix:"<<endl;
  cout<<matr.print()<<endl;
  cout<<"Vector:"<<endl;
  cout<<vec.print()<<endl;

  //Check matrix and vector multiplication.
  cout<<"Matrix Squared:"<<endl;
  cout<<(matr*matr).print()<<endl;
  cout<<"Matrix-Vector Multiplication:"<<endl;
  cout<<(matr*vec).print()<<endl;

  //Check LU Decomposition
  vector<vector<double> > LU_decomposed = LU_decomp(matr);
  themat L = themat(LU_decomposed.at(0));
  themat U = themat(LU_decomposed.at(1));

  cout<<"L:"<<endl<<L.print()<<endl;
  cout<<"U:"<<endl<<U.print()<<endl;
  cout<<"L*U:"<<endl<<(L*U).print()<<endl;

  //Check Gaussian Elimination
  cout<<"Gaussian elimination:"<<endl<<(gauss_elim(matr,vec)).print()<<endl;

}

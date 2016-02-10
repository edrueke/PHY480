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

#include "classes.C"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const int dim = 10;

double function(double &x){
  /*
    Takes in an x-value and returns our test function at that point 
    f(x) = 100e^{-10x}.
  */

  return 100*exp(-10*x);
}

vector<vector<double> > LU_decomp_special(){
  /*
    Computes the L and U decomposition of the matrix particular to this 
    problem.
  */

  themat L = themat(dim);
  themat U = themat(dim);
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
	U.point[i][j]=(i+2.0)/(i+1.0);
	L.point[i][j]=1;
      }
      else if((i-1)==j){
	L.point[i][j]=(-1.0*j)/i;
	U.point[i][j]=0;
      }
      else if(i+1==j){
	L.point[i][j]=0;
	U.point[i][j]=-1;
      }
      else{
	L.point[i][j]=0;
	U.point[i][j]=0;
      }
    }
  }

  for(int i=1;i<dim;i++){
    if(i!=(dim-1))
      L.point[i][i-1]=L[i+1][i];
    else
      L.point[i][i-1]=(-1.0*i)/(i+1);
  }

  vector<vector<double> > to_ret;
  vector<double> L_ret, U_ret;

  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      L_ret.push_back(L[i][j]);
      U_ret.push_back(U[i][j]);
    }
  }

  to_ret.push_back(L_ret);
  to_ret.push_back(U_ret);

  return to_ret;
}

thevec LU_decomp_solver_special(thevec &vec){
  /*
    Special solver for this particular problem.
  */

  thevec y = thevec(dim);

  y.point[0]=vec[0];
  for(int i=1;i<dim;i++){
    y.point[i]=vec[i]+(i/(i+1.0))*y[i-1];
  }

  thevec to_ret = thevec(dim);
  to_ret.point[dim-1]=1.0*dim*y[dim-1]/(dim+1);
  for(int i=dim-2;i>-1;i--){
    to_ret.point[i] = ((i+1.0)/(i+2))*(y[i]+to_ret[i+1]);
  }

  return to_ret;
}

void project1(){
  /*The main function of the script.*/

  //Define the matrix.
  themat mat = themat(dim);
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j)
	mat.point[i][j]=2.0;
      else if(i==j+1||i==j-1)
	mat.point[i][j]=-1.0;
      else
	mat.point[i][j]=0;
    }
  }

  //Define the vector. NOTE THIS LETS x_i=0.
  thevec vec = thevec(dim);
  double h = 1.0/(dim);
  for(int i=0;i<dim;i++){
    double xi = (i)*h;
    vec.point[i]=function(xi);
  }

  //Solve with full Gaussian elimination.
  thevec gausselim = gauss_elim(mat,vec);

  //Solve with full LU decomposition.
//  vector<vector<double> > lu_vec_full = LU_decomp(mat);
//  themat L_full = themat(lu_vec_full.at(0));
//  themat U_full = themat(lu_vec_full.at(1));
//  thevec ludecomp_full = LU_decomp_solver(L_full,U_full,vec);
//
  //Solve with special LU decomposition.
  vector<vector<double> > lu_vec = LU_decomp_special();
  themat L = themat(lu_vec.at(0));
  themat U = themat(lu_vec.at(1));
  thevec ludecomp = LU_decomp_solver_special(vec);

  cout<<"A:"<<endl<<mat.print()<<endl;
  cout<<"Gaussian elimination: "<<endl<<gausselim.print()<<endl;
//  cout<<"L*U full:"<<endl<<(L_full*U_full).print()<<endl;
//  cout<<"LU Decomposition: "<<endl<<ludecomp_full.print()<<endl;
  cout<<"L*U special:"<<endl<<(L*U).print()<<endl;
  cout<<"LU Decomposition Special: "<<endl<<ludecomp.print()<<endl;

}

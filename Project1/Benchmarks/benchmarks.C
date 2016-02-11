/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This script defines certain benchmarks which should be reproducible. In 
particular, it tests the classes (thevec, themat) and checks the full 
Gaussian, full LU decomposition, and special LU decomposition for various 
matrices.
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
#include "time.h"
#include "TMultiGraph.h"

#include "classes.C"

using namespace std;

void benchmarks(){
  /*
    Runs the various functions and outputs the results to a file which should
    be reproducible.
  */

  const int dim = 2;

  //Define the matrix
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

  //Define the vector. 
  thevec vec = thevec(dim);
  double h = 1.0/(dim);
  for(int i=0;i<dim;i++){
    double xi = (i)*h;
    vec.point[i]=function(xi);
  }

  //Solve with full Gaussian elimination.
  thevec gausselim = gauss_elim(mat,vec);

  //Solve with full LU decomposition.
  vector<vector<double> > lu_vec_full = LU_decomp(mat);
  themat L_full = themat(lu_vec_full.at(0));
  themat U_full = themat(lu_vec_full.at(1));
  thevec ludecomp_full = LU_decomp_solver(L_full,U_full,vec);

  //Solve with special LU decomposition.
  vector<vector<double> > lu_vec = LU_decomp_special(dim);
  themat L = themat(lu_vec.at(0));
  themat U = themat(lu_vec.at(1));
  thevec ludecomp = LU_decomp_solver_special(vec);

  cout<<"A:"<<endl<<mat.print()<<endl;
  cout<<"Gaussian elimination: "<<endl<<gausselim.print()<<endl;
  cout<<"L*U full:"<<endl<<(L_full*U_full).print()<<endl;
  cout<<"LU Decomposition: "<<endl<<ludecomp_full.print()<<endl;
  cout<<"L*U special:"<<endl<<(L*U).print()<<endl;
  cout<<"LU Decomposition Special: "<<endl<<ludecomp.print()<<endl;

}

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

  //Define the output file.
  ofstream myfile;
  myfile.open("benchmarks.txt");

  vector<int> dims;
  dims.push_back(1); dims.push_back(2); dims.push_back(3); dims.push_back(5);
  dims.push_back(10); dims.push_back(20); dims.push_back(50);

  for(unsigned int ct=0;ct<dims.size();ct++){
    
    int dim = dims.at(ct);

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
    thevec ludecomp = LU_decomp_solver_special(vec,L,U);
    
    //Define and fill an output file for comparison later.
    myfile<<"Dimensions: "<<dim<<endl<<endl;
    myfile<<"A:"<<endl<<mat.print()<<endl<<endl<<endl;
    
    myfile<<"Gaussian elimination: "<<endl<<endl;
    myfile<<"Solution:"<<endl<<gausselim.print()<<endl<<endl<<endl;
    
    myfile<<"Full LU Decomposition:"<<endl<<endl;
    myfile<<"L:"<<endl<<L_full.print()<<endl;
    myfile<<"U:"<<endl<<U_full.print()<<endl;
    myfile<<"L*U:"<<endl<<(L_full*U_full).print()<<endl;
    myfile<<"Solution: "<<endl<<ludecomp_full.print()<<endl<<endl<<endl;
    
    myfile<<"Special LU Decomposition:"<<endl<<endl;
    myfile<<"L:"<<endl<<L.print()<<endl;
    myfile<<"U:"<<endl<<U.print()<<endl;
    myfile<<"L*U:"<<endl<<(L*U).print()<<endl;
    myfile<<"Solution: "<<endl<<ludecomp.print()<<endl<<endl<<endl;
    
  }
  myfile.close();
}

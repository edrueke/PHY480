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
#include "time.h"

#include "classes.C"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const int dim = 2;

double function(double &x){
  /*
    Takes in an x-value and returns our test function at that point 
    f(x) = 100e^{-10x}.
  */

  return 100*exp(-10*x);
}

vector<vector<double> > LU_decomp_special(int sz){
  /*
    Computes the L and U decomposition of the matrix particular to this 
    problem.
  */

  themat L = themat(sz);
  themat U = themat(sz);
  
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
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

  for(int i=1;i<sz;i++){
    if(i!=(sz-1))
      L.point[i][i-1]=L[i+1][i];
    else
      L.point[i][i-1]=(-1.0*i)/(i+1);
  }

  vector<vector<double> > to_ret;
  vector<double> L_ret, U_ret;

  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
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

  int sz = vec.sz;

  thevec y = thevec(sz);

  y.point[0]=vec[0];
  for(int i=1;i<sz;i++){
    y.point[i]=vec[i]+(i/(i+1.0))*y[i-1];
  }

  thevec to_ret = thevec(sz);
  to_ret.point[sz-1]=1.0*sz*y[sz-1]/(sz+1);
  for(int i=sz-2;i>-1;i--){
    to_ret.point[i] = ((i+1.0)/(i+2))*(y[i]+to_ret[i+1]);
  }

  return to_ret;
}

void check_timing(){
  /*
    This function defines a clock and determines the time dependence of the 
    various linear system solvers. It plots this dependence based on the 
    number of parameters and fits with a TGraph.

    NOTE: ANY DIMENSION PAST ABOUT 500 TAKES AN EXCEPTIONAL AMOUNT OF TIME
    FOR THE GAUSSIAN ELIMINATION METHOD.
  */
  
  //Define a vector of the dimensions we want to check.
  vector<int> dimensions;
  dimensions.push_back(50); dimensions.push_back(60); 
  dimensions.push_back(75); dimensions.push_back(100); 
  dimensions.push_back(125); dimensions.push_back(150); 
  dimensions.push_back(175); dimensions.push_back(200); 
  dimensions.push_back(250); dimensions.push_back(300); 
  
  cout<<"A"<<endl;
  
  //Define graphs for use with plotting.
  TGraph *g_gauss = new TGraph(10);
  TGraph *g_ludecompfull = new TGraph(10);
  TGraph *g_ludecompspec = new TGraph(10);

  cout<<"B"<<endl;

  for(unsigned int ct=0;ct<dimensions.size();ct++){
   
    cout<<"C"<<endl;

    int dm = dimensions.at(ct);

    cout<<"n = "<<dm<<endl;
    cout<<"D"<<endl;

    //Define the special matrix
    themat mat = themat(dm);
    for(int i=0;i<dm;i++){
      for(int j=0;j<dm;j++){
	if(i==j)
	  mat.point[i][j]=2.0;
	else if(i==j+1||i==j-1)
	  mat.point[i][j]=-1.0;
	else
	  mat.point[i][j]=0;
      }
    }

    cout<<"E"<<endl;

    //Define the solution vector.
    thevec vec = thevec(dm);
    double h = 1.0/(dm);
    for(int i=0;i<dm;i++){
      double xi = (i)*h;
      vec.point[i]=function(xi);
    }

    cout<<"F"<<endl;

    //Check Gaussian Elimination.
    clock_t start_gauss, finish_gauss;
    cout<<"F1"<<endl;
    start_gauss = clock();
    cout<<"F2"<<endl;
    thevec gausselim = gauss_elim(mat,vec);
    cout<<"F3"<<endl;
    finish_gauss = clock();

    cout<<"G"<<endl;

    g_gauss->SetPoint(ct,dm,1.0*(finish_gauss-start_gauss)/CLOCKS_PER_SEC);
    cout<<1.0*(finish_gauss-start_gauss)/CLOCKS_PER_SEC<<endl;
    cout<<"H"<<endl;

    //Check the full LU decomposition.
    clock_t start_lufull, finish_lufull;
    start_lufull = clock();
    vector<vector<double> > lu_vec_full = LU_decomp(mat);
    themat L_full = themat(lu_vec_full.at(0));
    themat U_full = themat(lu_vec_full.at(1));
    thevec ludecomp_full = LU_decomp_solver(L_full,U_full,vec);
    finish_lufull = clock();

    cout<<"I"<<endl;

    g_ludecompfull->SetPoint(ct,dm,1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC);
    cout<<1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC<<endl;
    cout<<"J"<<endl;

    //Check special LU decomposition.
    clock_t start_luspec, finish_luspec;
    start_luspec = clock();
    vector<vector<double> > lu_vec = LU_decomp_special(dm);
    themat L = themat(lu_vec.at(0));
    themat U = themat(lu_vec.at(1));
    thevec ludecomp = LU_decomp_solver_special(vec);
    finish_luspec = clock();

    cout<<"K"<<endl;

    g_ludecompspec->SetPoint(ct,dm,1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC);
    cout<<1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC<<endl;
    cout<<"L"<<endl;
  }

  cout<<"M"<<endl;

  //Plot the graphs
  g_gauss->SetLineColor(kRed);
  g_gauss->SetMarkerColor(kRed);

  g_ludecompfull->SetLineColor(kGreen);
  g_ludecompfull->SetMarkerColor(kGreen);

  g_ludecompspec->SetLineColor(kBlue);
  g_ludecompspec->SetMarkerColor(kBlue);

  cout<<"N"<<endl;

  TLegend *leg = new TLegend(0.15,0.13,0.43,0.35);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  cout<<"O"<<endl;

  leg->AddEntry(g_gauss,"Full Gaussian Elimination","C");
  leg->AddEntry(g_ludecompfull,"Full LU Decomposition","C");
  leg->AddEntry(g_ludecompspec,"Special LU Decomposition","C");

  cout<<"P"<<endl;

  TCanvas *c_time = new TCanvas("c_time","c_time",800,720);
  c_time->SetBorderMode(0);
  c_time->SetLogy(1);

  cout<<"Q"<<endl;

  c_time->cd();
  g_gauss->Draw();
  g_ludecompfull->Draw("SAME");
  g_ludecompspec->Draw("SAME");
  leg->Draw("SAME");
  
  cout<<"R"<<endl;

  c_time->SaveAs("time_plot.png");
  c_time->SaveAs("time_plot.root");

}

void project1(){
  /*The main function of the script.*/

  //Check the timing and create the timing plot
  check_timing();

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

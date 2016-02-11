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
#include "TMultiGraph.h"

#include "classes.C"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const int dim = 10;

void check_timing(){
  /*
    This function defines a clock and determines the time dependence of the 
    various linear system solvers. It plots this dependence based on the 
    number of parameters and fits with a TGraph.

    NOTE: ANY DIMENSION PAST ABOUT 500 TAKES AN EXCEPTIONAL AMOUNT OF TIME
    FOR THE GAUSSIAN ELIMINATION METHOD.
  */
  
  gStyle->SetOptFit();

  //Define a vector of the dimensions we want to check.
  vector<int> dimensions;
  dimensions.push_back(50); dimensions.push_back(60); 
  dimensions.push_back(75); dimensions.push_back(100); 
  dimensions.push_back(125); dimensions.push_back(150); 
  dimensions.push_back(175); dimensions.push_back(200); 
  dimensions.push_back(250); dimensions.push_back(300); 
    
  //Define graphs for use with plotting.
  TGraph *g_gauss = new TGraph(10);
  TGraph *g_ludecompfull = new TGraph(10);
  TGraph *g_ludecompspec = new TGraph(10);

  for(unsigned int ct=0;ct<dimensions.size();ct++){
   
    int dm = dimensions.at(ct);

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

    //Define the solution vector.
    thevec vec = thevec(dm);
    double h = 1.0/(dm);
    for(int i=0;i<dm;i++){
      double xi = (i)*h;
      vec.point[i]=function(xi);
    }

    //Check Gaussian Elimination.
    clock_t start_gauss, finish_gauss;
    start_gauss = clock();
    thevec gausselim = gauss_elim(mat,vec);
    finish_gauss = clock();

    g_gauss->SetPoint(ct,dm,1.0*(finish_gauss-start_gauss)/CLOCKS_PER_SEC);
    
    //Check the full LU decomposition.
    clock_t start_lufull, finish_lufull;
    start_lufull = clock();
    vector<vector<double> > lu_vec_full = LU_decomp(mat);
    themat L_full = themat(lu_vec_full.at(0));
    themat U_full = themat(lu_vec_full.at(1));
    thevec ludecomp_full = LU_decomp_solver(L_full,U_full,vec);
    finish_lufull = clock();

    g_ludecompfull->SetPoint(ct,dm,1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC);

    //Check special LU decomposition.
    clock_t start_luspec, finish_luspec;
    start_luspec = clock();
    vector<vector<double> > lu_vec = LU_decomp_special(dm);
    themat L = themat(lu_vec.at(0));
    themat U = themat(lu_vec.at(1));
    thevec ludecomp = LU_decomp_solver_special(vec);
    finish_luspec = clock();

    g_ludecompspec->SetPoint(ct,dm,1.0*(finish_lufull-start_lufull)/CLOCKS_PER_SEC);

  }

  //Plot the graphs
  g_gauss->SetLineColor(kRed);
  g_gauss->SetMarkerColor(kRed);
  g_gauss->Fit("pol3","q");

  g_ludecompfull->SetLineColor(kGreen);
  g_ludecompfull->SetMarkerColor(kGreen);
  g_ludecompfull->Fit("pol3","q");

  g_ludecompspec->SetLineColor(kBlue);
  g_ludecompspec->SetMarkerColor(kBlue);
  g_ludecompspec->Fit("pol3","q");

  TLegend *leg = new TLegend(0.50,0.13,0.55,0.28);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  leg->AddEntry(g_gauss,"Full Gaussian Elimination","L");
  leg->AddEntry(g_ludecompfull,"Full LU Decomposition","L");
  leg->AddEntry(g_ludecompspec,"Special LU Decomposition","L");

  TMultiGraph *m_time = new TMultiGraph("m_time","Timing for Various Algorithms");
  m_time->SetTitle("Timing for Various Algorithms;Dimension of Matrix;Time (s)");
  m_time->Add(g_gauss);
  m_time->Add(g_ludecompfull);
  m_time->Add(g_ludecompspec);

  TCanvas *c_time = new TCanvas("c_time","c_time",800,720);
  c_time->SetBorderMode(0);
  c_time->SetLogy(1);

  c_time->cd();
  m_time->Draw("AC*");
  leg->Draw("SAME");
  
  c_time->SaveAs("time_plot.png");
  c_time->SaveAs("time_plot.root");

}

void project1(){
  /*The main function of the script.*/

  //Check the timing and create the timing plot
  //check_timing();

  //Define the matrix. We defined it to be dim-1 in size because we know the 
  //boundary conditions.  This way we are solving for i=1,..,n, dim is n+1 and
  //we can include the boundary conditions at i=0 and i=dim at the end.
  themat mat = themat(dim-1);
  for(int i=0;i<dim-1;i++){
    for(int j=0;j<dim-1;j++){
      if(i==j)
	mat.point[i][j]=2.0;
      else if(i==j+1||i==j-1)
	mat.point[i][j]=-1.0;
      else
	mat.point[i][j]=0;
    }
  }

  //Define the vector. 
  thevec vec = thevec(dim-1);
  double h = 1.0/(dim);
  for(int i=0;i<dim-1;i++){
    double xi = (i+1)*h;
    vec.point[i]=function(xi);
  }

  //Solve with special LU decomposition.
  vector<vector<double> > lu_vec = LU_decomp_special(dim-1);
  themat L = themat(lu_vec.at(0));
  themat U = themat(lu_vec.at(1));
  thevec solution = LU_decomp_solver_special(vec);

  //Multiply by h^2 because we haven't taken that into account yet.
  solution=solution*(pow(h,2));

  //Plot the result
  TGraph *g_solution = new TGraph(dim+1);
  
  for(int i=0;i<dim+1;i++){
    //Dirichlet Boundary Conditions
    if(i==0)
      g_solution->SetPoint(i,0,0);
    else if(i==dim+1)
      g_solution->SetPoint(i,1,0);
    
    //Solution vector
    else
      g_solution->SetPoint(i,i*h,solution[i-1]);
  }

  TCanvas *c_solution = new TCanvas("c_solution","c_solution",800,720);
  c_solution->SetBorderMode(0);

  c_solution->cd();
  g_solution->Draw("AC*");
  
  c_solution->SaveAs("solution_plot.png");
  c_solution->SaveAs("solution_plot.root");


}

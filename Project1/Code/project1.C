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
#include "TColor.h"

#include "classes.C"

using namespace std;

void check_timing_all(){
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
  dimensions.push_back(2); dimensions.push_back(3); dimensions.push_back(5);
  dimensions.push_back(10); dimensions.push_back(30);
  dimensions.push_back(50); dimensions.push_back(60); 
  dimensions.push_back(75); dimensions.push_back(100); 
  dimensions.push_back(125); dimensions.push_back(150); 
  dimensions.push_back(175); dimensions.push_back(200); 
  //dimensions.push_back(250); dimensions.push_back(300); 
    
  //Define graphs for use with plotting.
  TGraph *g_gauss = new TGraph(13);
  TGraph *g_ludecompfull = new TGraph(13);
  TGraph *g_ludecompspec = new TGraph(13);

  for(unsigned int ct=0;ct<dimensions.size();ct++){
   
    int dm = dimensions.at(ct);

    cout<<"dim = "<<dm<<endl;

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
    thevec ludecomp = LU_decomp_solver_special(vec,L,U);
    finish_luspec = clock();

    g_ludecompspec->SetPoint(ct,dm,1.0*(finish_luspec-start_luspec)/CLOCKS_PER_SEC);

  }

  //Plot the graphs
  g_gauss->SetLineColor(kRed);
  g_gauss->SetMarkerColor(kRed);

  g_ludecompfull->SetLineColor(kGreen);
  g_ludecompfull->SetMarkerColor(kGreen);

  g_ludecompspec->SetLineColor(kBlue);
  g_ludecompspec->SetMarkerColor(kBlue);

  TLegend *leg = new TLegend(0.20,0.4,0.5,0.6);
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
  c_time->SetLogy(0);

  c_time->cd();
  m_time->Draw("AC*");
  leg->Draw("SAME");
  
  c_time->SaveAs("plots/time_plot.png");
  c_time->SaveAs("plots/time_plot.root");

}

void check_timing_short(){
  /*
    This function defines a clock and determines the time dependence of the 
    various linear system solvers. It plots this dependence based on the 
    number of parameters and fits with a TGraph.

    NOTE: ONLY INVESTIGATES THE TWO LU DECOMP SOLVERS
  */
  
  gStyle->SetOptFit();

  //Define a vector of the dimensions we want to check.
  vector<int> dimensions;
  dimensions.push_back(2); dimensions.push_back(3); dimensions.push_back(5);
  dimensions.push_back(10); dimensions.push_back(30);
  dimensions.push_back(50); dimensions.push_back(60); 
  dimensions.push_back(75); dimensions.push_back(100); 
  dimensions.push_back(125); dimensions.push_back(150); 
  dimensions.push_back(175); dimensions.push_back(200); 
  dimensions.push_back(250); dimensions.push_back(300); 
  dimensions.push_back(500); 

  //Define graphs for use with plotting.
  TGraph *g_ludecompfull = new TGraph(13);
  TGraph *g_ludecompspec = new TGraph(13);

  for(unsigned int ct=0;ct<dimensions.size();ct++){
   
    int dm = dimensions.at(ct);

    cout<<"dim = "<<dm<<endl;

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
    thevec ludecomp = LU_decomp_solver_special(vec,L,U);
    finish_luspec = clock();

    g_ludecompspec->SetPoint(ct,dm,1.0*(finish_luspec-start_luspec)/CLOCKS_PER_SEC);

  }

  //Plot the graphs
  g_ludecompfull->SetLineColor(kGreen);
  g_ludecompfull->SetMarkerColor(kGreen);

  g_ludecompspec->SetLineColor(kBlue);
  g_ludecompspec->SetMarkerColor(kBlue);

  TLegend *leg = new TLegend(0.20,0.4,0.5,0.6);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  leg->AddEntry(g_ludecompfull,"Full LU Decomposition","L");
  leg->AddEntry(g_ludecompspec,"Special LU Decomposition","L");

  TMultiGraph *m_time = new TMultiGraph("m_time","Timing for Various Algorithms");
  m_time->SetTitle("Timing for Various Algorithms;Dimension of Matrix;Time (s)");
  m_time->Add(g_ludecompfull);
  m_time->Add(g_ludecompspec);

  TCanvas *c_time = new TCanvas("c_time","c_time",800,720);
  c_time->SetBorderMode(0);
  c_time->SetLogy(0);

  c_time->cd();
  m_time->Draw("AC*");
  leg->Draw("SAME");
  
  c_time->SaveAs("plots/time_plot_short.png");
  c_time->SaveAs("plots/time_plot_short.root");

}

TGraph* make_sol_plots(int dim, TF1* func){
  /*
    Makes graphs of the solutions with the correct function solution.
  */

  //Define the matrix. We defined it to be dim-1 in size because we know the 
  //boundary conditions. This way we are solving for i=1,..,n, dim is n+1 and
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
  thevec solution = LU_decomp_solver_special(vec,L,U);
  
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
  
  /*TLegend *leg = new TLegend(0.50,0.13,0.55,0.28);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);
  
  g_solution->SetLineColor(kBlack);
  func->SetLineColor(kRed);
  func->SetLineStyle(10);
  
  string leg_entry = "n = "+to_string(dim+1);
  leg->AddEntry(g_solution,leg_entry.c_str(),"L");
  leg->AddEntry(func,"Accepted","L");

  TCanvas *c_solution = new TCanvas("c_solution","c_solution",800,720);
  c_solution->SetBorderMode(0);
  
  c_solution->cd();
  func->Draw();
  g_solution->Draw("SAME AC*");
  leg->Draw("SAME");

  string name1 = "plots/solution_plot"+to_string(dim+1)+".png";
  string name2 = "plots/solution_plot"+to_string(dim+1)+".root";
  c_solution->SaveAs(name1.c_str());
  c_solution->SaveAs(name2.c_str());*/

  return g_solution;
}

void project1(){
  /*
    The main function of the script.
  */

  //Check the timing and create the timing plot
  //check_timing_all();
  //check_timing_short();

  //Define the dimensions to test
  vector<int> dims;
  dims.push_back(2); dims.push_back(3); dims.push_back(5);
  dims.push_back(10); dims.push_back(50); dims.push_back(100); 
  dims.push_back(500); dims.push_back(1000);
  dims.push_back(10000); dims.push_back(100000); 

  //Define various colors
  vector<int> colors;
  for(int i=2;i<20;i++){
    colors.push_back(i);
  }

  //Define the accepted solution
  TF1 *func = new TF1("func","x/exp(10)-x-exp(-10*x)+1",0,1);

  //Define a multigraph to hold all calculated solutions
  TMultiGraph *m_sol = new TMultiGraph("m_sol","Solutions for Various Numbers of Steps");
  m_sol->SetTitle("Solutions for Various Numbers of Steps;x_{i};u(x)");

  //Define a legend 
  TLegend *leg = new TLegend(0.67,0.55,0.8,0.85);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  func->SetLineColor(1);
  leg->AddEntry(func,"Accepted");

  //Work with errors.
  vector<double> errors;

  //Make a graph for each number of steps listed in dims.
  for(unsigned int my_ct=0;my_ct<dims.size();my_ct++){
    
    //Hold the maximum error.
    double max_error = 0;

    int dim = dims.at(my_ct);
    cout<<"dim = "<<dim<<endl;

    //Make the tgraph
    TGraph* hold = make_sol_plots(dim,func);

    //Calculate errors.
    int n_points = hold->GetN();
    double x[n_points]; double y[n_points];
    for(int i=0;i<n_points;i++){
      hold->GetPoint(i,x[i],y[i]);
      double ans = sol_function(x[i]);
      double err = log10(fabs((y[i]-ans)/ans));
      //cout<<"i = "<<i<<"; x = "<<x[i]<<"; y = "<<y[i]<<"; sol = "<<ans<<"; err = "<<err<<endl;
      if(fabs(err)>fabs(max_error) && i!=0 && i!=n_points-1)
	max_error=err;
    }

    errors.push_back(max_error);

    //Style the tgraph
    hold->SetLineStyle(my_ct%10+1);
    hold->SetLineColor(colors.at(my_ct));
    hold->SetFillColor(0);
    hold->SetMarkerColor(colors.at(my_ct));

    //Add tgraph to legend and multigraph
    string name = "n = "+to_string(dim+1);
    leg->AddEntry(hold,name.c_str());
    m_sol->Add(hold);

  }
  
  //Plot the graphs with the legend
  TCanvas *c_solution = new TCanvas("c_solution","c_solution",800,720);
  c_solution->SetBorderMode(0);
  
  func->SetLineColor(1);
  func->SetLineStyle(10);

  c_solution->cd();
  m_sol->Draw("AC");
  func->Draw("SAME");
  leg->Draw("SAME");
  
  c_solution->SaveAs("plots/solution_plot.png");
  c_solution->SaveAs("plots/solution_plot.root");
  
  //Write out the errors to a file
  ofstream myfile;
  myfile.open("plots/errors.txt");
  myfile<<"Maximum errors per step size"<<endl<<endl;
  myfile<<"Number of Steps | Maximum Error"<<endl<<endl;

  for(unsigned int i=0;i<errors.size();i++){
    myfile<<dims.at(i)+1<<" | "<<errors.at(i)<<endl;
  }

  myfile.close();
}

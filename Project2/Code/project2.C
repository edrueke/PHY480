/*
Elizabeth Drueke
PHY 480
March 4, 2016

Project 2

The solution to the problem of electron(s) trapped in a harmonic oscillator 
well.
*/

#include <iostream>

#include "TStyle.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "time.h"
#include "classes.C"

using namespace std;

void partb(){
  /*
    Investigate the various questions posed in part b) of the project 
    outline:
    -How many steps are needed to get the lowest three eigenvalues with 4
    leading digits? (3, 7, 1, ...)
    -What is the dependency of the eigenvalues of the choice of p_max?
    -How many transformations are needed in the Jacobi algorithm?  How does
    this depend on the matrix dimension?
    -Check the time required against Householder's algorithm.
  */

  //First look at the number of points needed to get the lowest three 
  //eigenvalues with four leading digits for various p_max choices.
  double pmin = 0.0;
  vector<double> p_maxes;
  for(int i=10;i<11;i++){
    p_maxes.push_back(i*1.0);
  }
  vector<double> n_steps;
  for(int i=3;i<6;i++){
    n_steps.push_back(i);
  }
  /*  n_steps.push_back(100);
  n_steps.push_back(200);
  n_steps.push_back(300);
  n_steps.push_back(400);*/

  //Create a graph of the number of rotations required for a given number of steps.
  TGraph *g_stepvrot = new TGraph(p_maxes.size()*n_steps.size());
  int pt = 0;

  //Create a graph of the time required for the Householder and Jacobi algorithms
  TGraph *g_timehous = new TGraph(n_steps.size());
  TGraph *g_timejac = new TGraph(n_steps.size());

  for(int i=0;i<p_maxes.size();i++){
    for(int j=0;j<n_steps.size();j++){
      double pmax = p_maxes.at(i);
      double nstep = n_steps.at(j);
      double h = (pmax-pmin)/nstep;
      
      themat mat = themat(nstep);
      
      thevec ds = thevec(nstep);
      thevec es = thevec(nstep);
      
      /*thevec rho = thevec(nstep+1);
      for(int i=0;i<rho.sz;i++){
	rho.point[i]=pmin+i*h;
      }
      rho.point[rho.sz-1]=0;

      thevec vs = thevec(nstep+1);
      for(int i=0;i<vs.sz;i++)
	vs.point[i]=pow((rho[i]/2),2);

      mat.point[0][0]=2/(h*h)+vs[1];
      mat.point[0][1]=-1/(h*h);
      for(int i=1;i<nstep-2;i++){
	mat.point[i,i-1]=-1/(h*h);
	mat.point[i][i]=2/(h*h)+vs[i+1];
	mat.point[i][i+1]=-1/(h*h);
      }

      mat.point[nstep-2][nstep-2]=2/(h*h)+vs[nstep-1];
      mat.point[nstep-2][nstep-3]=-1/(h*h);*/

      for(int m=1;m<nstep-1;m++){
	double p_m = pmin+m*h;
	double V_m = pow(p_m,2);
	double d_m = 2/(h*h)+V_m;
	double e_m = -1/(h*h);
	
	ds.point[m-1]=d_m;
	es.point[m-1]=e_m;

	for(int k=0;k<mat.sz;k++){
	  for(int l=0;l<mat.sz;l++){
	    if(k==l)
	      mat.point[k][l]=d_m;
	    else if(k==l-1 || k==l+1)
	      mat.point[k][l]=e_m;
	    else
	      mat.point[k][l]=0;
	  }
	}
      }
      
      //How long does the Jacobi algorithm take?
      clock_t start_jac, finish_jac;
      start_jac = clock();

      //Diagonalize using the Jacobi Rotation Algorithm
      themat diag = Jacobi_Method(mat,1e-8,g_stepvrot,&pt);
      //cout<<"Matrix: "<<endl<<mat.print()<<endl<<"diag: "<<endl<<diag.print()<<endl;
      /*diag = Jacobi_Method_step(mat,1e-8);
	cout<<"diag: "<<diag.print()<<endl;
      */finish_jac = clock();
      
      if(i==0)
	g_timejac->SetPoint(j+1,nstep-1,finish_jac-start_jac);
      
      //Diagonalize using Householder's Algorithm from lib.cpp
      themat z = themat(nstep-1);
      for(int c=0;c<nstep-1;c++){
	for(int d=0;d<nstep-1;d++){
	    if(c==d)
	      z.point[c][d]=1;
	    else
	      z.point[c][d]=0;
	}
      }

      //How long does the Householder algorithm take?

      clock_t start_hous, finish_hous;
      start_hous = clock();
      tqli (ds.point, es.point, nstep-1,z.point);
      finish_hous = clock();

      if(i==0)
	g_timehous->SetPoint(j+1,nstep-1,finish_hous-start_hous);

      //ds should now hold the eigenvalues
      thevec eigen_jac = thevec(nstep-1);
      for(int m=0;m<nstep-1;m++){
	if(abs(diag[m][m]-3)<1)
	  cout<<"sucess "<<diag[m][m]<<"; nstep = "<<nstep-1<<endl;
	if(abs(diag[m][m])<1)
	  cout<<"fail "<<diag[m][m]<<"; nstep = "<<nstep-1<<endl;
	eigen_jac.point[m] = diag[m][m];
      }

      cout<<"pmax = "<<pmax<<"; nstep = "<<nstep-1<<"; eigenvalues:"<<endl;
      for(int m=0;m<diag.sz;m++)
	cout<<"     "<<eigen_jac[m]<<"     "<<ds[m]<<endl;
      
    }
  }

  //Format and save the plot of steps vs rotations
  gStyle->SetOptFit();
  gStyle->SetStatX(.5);
  gStyle->SetStatY(.8);

  g_stepvrot->SetLineColor(kBlack);
  g_stepvrot->SetMarkerColor(kBlack);
  g_stepvrot->SetMarkerStyle(3);
  g_stepvrot->SetMarkerSize(1.5);
  g_stepvrot->SetTitle("Step Size vs. Number of Rotations to Diagonalize");

  g_stepvrot->Fit("pol2");

  TCanvas *c_stepvrot = new TCanvas("c_stepvrot","c_stepvrot",800,720);
  c_stepvrot->SetBorderMode(0);

  c_stepvrot->cd();
  g_stepvrot->Draw("ap");

  c_stepvrot->SaveAs("plots/stepvrot.pdf");
  c_stepvrot->SaveAs("plots/stepvrot.png");
  c_stepvrot->SaveAs("plots/stepvrot.root");
  c_stepvrot->Close();

  //Format and save the plot of time difference between 
  //householder and jacobi
  gStyle->SetOptFit(0);
  g_timejac->SetLineColor(kRed);
  g_timejac->SetMarkerColor(kRed);
  g_timejac->SetMarkerStyle(3);
  g_timejac->SetMarkerSize(1.5);

  g_timehous->SetLineColor(kBlue);
  g_timehous->SetMarkerColor(kBlue);
  g_timehous->SetMarkerStyle(3);
  g_timehous->SetMarkerSize(1.5);

  TLegend *leg = new TLegend(0.20,0.4,0.5,0.6);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(g_timejac,"Jacobi Algorithm","L");
  leg->AddEntry(g_timehous,"Householder Algorithm","L");

  TMultiGraph *m_time = new TMultiGraph("m_time",
					"Timing for Various Algorithms");
  m_time->SetTitle("Timing for Various Algorithms;Dimension of Matrix;Time (s)");
  m_time->Add(g_timejac);
  m_time->Add(g_timehous);

  TCanvas *c_time = new TCanvas("c_time","c_time",800,720);
  c_time->SetBorderMode(0);
  c_time->cd();
  m_time->Draw("A*");
  leg->Draw("SAME");

  c_time->SaveAs("plots/algtimes.pdf");
  c_time->SaveAs("plots/algtimes.png");
  c_time->SaveAs("plots/algtimes.root");
  c_time->Close();


  
}

void project2(){
  /*
    The main function which calls other functions and attempts to solve the 
    problem.
  */

  //Define the p_min, p_max, and n_steps
  double p_min=0.0; double p_max=10.0;
  int n_step = 5.0;
  double h = (p_max-p_min)/n_step;

  themat my_mat = themat(n_step-1);

  //pretty sure this should be nstep
  for(int i=0;i<p_max+1;i++){
    
    //p_i = p_min+i*h, i = 0,1,...,n
    double p_i = p_min+i*h;

    //Let v_i=p_i^2
    double V_i = pow(p_i,2);

    //Diagonal matrix element d_i = 2/h^2+v_i
    double d_i = 2/(h*h)+V_i;

    //Off-diagonal matrix element e_i = -1/h^2
    double e_i = -1/(h*h);

    //Define the matrix    
    for(int j=0;j<my_mat.sz;j++){
      for(int k=0;k<my_mat.sz;k++){
	if(j==k)
	  my_mat.point[j][k]=d_i;
	else if(j==k-1 || j==k+1)
	  my_mat.point[j][k]=e_i;
	else
	  my_mat.point[j][k]=0;
      }
    }
  }

  //Diagonalize the matrix
  //themat diag = Jacobi_Method(my_mat,1e-8);
  //cout<<"Matrix: "<<endl<<my_mat.print()<<endl<<"Diag:"<<endl<<diag.print()<<endl; 
  /*diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; 
  diag = Jacobi_Method_step(diag,1e-8);
  cout<<"Diag:"<<endl<<diag.print()<<endl; */

  partb();
}

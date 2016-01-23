/*
Elizabeth Drueke
PHY 480
January 20, 2016

Warm Up Exercise
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

inline void MyStyle(){
  
  gStyle->SetPalette(1);
  gStyle->SetHistLineWidth(3.0);
  gStyle->SetTitleX(0.8);
  gStyle->SetTitleY(0.87);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetStatBorderSize(1.0);
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.7);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineWidth(2);
  gStyle->SetOptStat(0);
}

/*Find the total error due to the numerical approximations.*/
/*This is due to the formulae themselves. The f2c has an error that goes
  as h and the f3c has an error that goes as h^2. We can also look at the 
  precision due to the rounding of the computers. In particular, float goes 
  to approximately 10^{-7} and double to 10^{-15}.*/

/*Compute the derivative of arctangent using equation 1 (double precision).*/
/*Arguments: double h is the step length
  double x is the value at which we are computing the derivative*/
double Derivf2c(double h, double x){
  
  double val1 = atan(x+h);
  double val2 = atan(x);
  return (val1-val2)/h;
  
}

/*Compute the derivative of arctangent using equation 1 (single precision).*/
/*Arguments: float &h is the step length
  float x is the value at which we are computing the derivative*/
float Derivf2cSingle(float &h, float x){
  
  float val1 = atan(x+h);
  float val2 = atan(x);
  return (val1-val2)/h;
}

/*Compute the derivative of arctangent using equation 1 (by ref)*/
/*Arguments: double &h is the step length passed by reference
  double &x is the value at which we are computing the derivative*/
double Derivf2cRef(double &h, double x){
  
  double val1 = atan(x+h);
  double val2 = atan(x);
  return (val1-val2)/h;
  
}

/*Compute the derivative of arctangent using equation 2 (single precition).*/
/*Arguments: float &h is the step length
  float x is the value at which we are computing the derivative*/
float Derivf3cSingle(float &h, float x){

  float val1 = atan(x+h);
  float val2 = atan(x-h);
  return (val1-val2)/(2*h);

}

/*Compute the derivative of arctan using equation 2.*/
/*Arguments: double h is the step length
  double x is the value at which we are computing the derivative*/
double Derivf3c(double h, double x){
  
  double val1 = atan(x+h);
  double val2 = atan(x-h);
  return (val1-val2)/(2*h);
  
}

/*Compute the derivative of arctangent using equation 2 (by ref)*/
/*Arguments: double &h is the step length passed by reference
  double &x is the value at which we are computing the derivative*/
double Derivf3cRef(double &h, double x){
  double val1 = atan(x+h);
  double val2 = atan(x-h);
  return (val1-val2)/(2*h);
}

/*Compute the error epsilon*/
/*Arguments: double calc is the calculated value of the derivative*/
double ComputeEpsilon(double calc){
  
  double exact = 0.333333333333333333333;
  return log10(fabs((calc-exact)/exact));
  
}

/*The main program calls the programs above to complete the warm up 
  exercise*/
int Warmup(){
  
  MyStyle();
  
  /*Create a file to be used for the output.*/
  
  ofstream outfile;
  outfile.open("output.txt");
  
  /*Compute the errors*/
  
  /*Create a dynamic array with the various step lengths.*/
  
  double *steps;
  
  steps = new double[10000];
  
  for(int n=0;n<10000;n++){
    steps[n]=0.00001+.0001*n;
  }
  
  cout<<endl;
  
  /*Order the steps*/
  for(int n=0;n<10000;n++){
    for(int i=0;i<n;i++){
      if(steps[i]>steps[n]){
	double hold = steps[i];
	steps[i]=steps[n];
	steps[n]=hold;
	i=0;
	n=0;
      }
    }
  }
  
  /*Compute the derivative for various step values using equation 1.*/
  
  vector<double> f2c;
  vector<float> f2c_single;
  
  for(int n=0;n<10000;n++){
    double sq = sqrt(2);
    double deriv1 = Derivf2c(steps[n],sq);
    double deriv1Ref = Derivf2cRef(steps[n],sq);
    float step = (float)(steps[n]);
    float sq2 = (float)(sq);
    float deriv1Sin = Derivf2cSingle(step,sq2);
    outfile<<"Derivative using equation 1 with step length "<<steps[n]<<" is "<<
      deriv1<<endl;
    f2c.push_back(deriv1);
    outfile<<"Derivative using equation 1 with ref with step length "<<
      steps[n]<<" is "<<deriv1Ref<<endl;
    outfile<<"Derivative using equation 1 with single precision with step "<<
      "length "<<steps[n]<<" is "<<deriv1Sin<<endl<<endl;
    f2c_single.push_back(deriv1Sin);
    
  }
  
  outfile<<"##############################################################"<<
    endl;
  
  /*Compute the derivative for various step values using equation 2. */
  
  vector<double> f3c;
  vector<float> f3c_single;

  for(int n=0;n<10000;n++){
    double sq = sqrt(2);
    double deriv2 = Derivf3c(steps[n],sq);
    double deriv2Ref = Derivf3cRef(steps[n],sq);
    float step = (float)(steps[n]);
    float sq2 = (float)(sq);
    float deriv2Sin = Derivf3cSingle(step,sq2);
    outfile<<"Derivative using equation 2 with setp length "<<steps[n]<<" is "<<
      deriv2<<endl;
    f3c.push_back(deriv2);
    outfile<<"Derivative using equation 1 with ref with step length "<<
      steps[n]<<" is "<<deriv2Ref<<endl;
    outfile<<"Derivative using equation 2 with single precition with step "<<
      "length "<<steps[n]<<" is "<<deriv2Sin<<endl<<endl;
    f3c_single.push_back(deriv2Sin);
  }
  
  outfile<<"###############################################################"<<
    endl;
  
  /*Compute the epsilon errors (Step 3) and plot them against step size.*/
  
  vector<double> epsilonsf2c;
  vector<double> epsilonsf3c;
  vector<double> epsilonsf2c_single;
  vector<double> epsilonsf3c_single;

  for(int n=0;n<10000;n++){
    epsilonsf2c.push_back(ComputeEpsilon(f2c.at(n)));
    epsilonsf3c.push_back(ComputeEpsilon(f3c.at(n)));
    epsilonsf2c_single.push_back(ComputeEpsilon(f2c_single.at(n)));
    epsilonsf3c_single.push_back(ComputeEpsilon(f3c_single.at(n)));
  }
  
  TGraph *g_f2cerrors = new TGraph(11);
  TGraph *g_f3cerrors = new TGraph(11);
  TGraph *g_f2cerrors_sin = new TGraph(11);
  TGraph *g_f3cerrors_sin = new TGraph(11);
  
  for(int n=0;n<10000;n++){
    g_f2cerrors->SetPoint(n,steps[n],epsilonsf2c.at(n));
    g_f3cerrors->SetPoint(n,steps[n],epsilonsf3c.at(n));
    g_f2cerrors_sin->SetPoint(n,steps[n],epsilonsf2c_single.at(n));
    g_f3cerrors_sin->SetPoint(n,steps[n],epsilonsf3c_single.at(n));
  }
  
  g_f2cerrors->SetLineColor(kRed);
  g_f2cerrors->SetMarkerColor(kRed);
  g_f3cerrors->SetLineColor(kCyan);
  g_f3cerrors->SetMarkerColor(kCyan);
  g_f2cerrors_sin->SetLineColor(kMagenta);
  g_f2cerrors_sin->SetMarkerColor(kMagenta);
  g_f3cerrors_sin->SetLineColor(kBlue);
  g_f3cerrors_sin->SetMarkerColor(kBlue);

  TCanvas *c_errors = new TCanvas("c_errors","c_errors",800,720);
  c_errors->SetBorderMode(0);
  
  TLegend *leg = new TLegend(0.45,0.13,0.70,0.35);
  leg->AddEntry(g_f2cerrors,"f2c","lp");
  leg->AddEntry(g_f3cerrors,"f3c","lp");
  leg->AddEntry(g_f2cerrors_sin,"f2c single","lp");
  leg->AddEntry(g_f3cerrors_sin,"f3c single","lp");

  TMultiGraph *m_errors = new TMultiGraph("m_errors","Epsilon Errors");
  m_errors->SetTitle("Epsilon Errors;Step Size;Epsilon");
  
  m_errors->Add(g_f2cerrors);
  m_errors->Add(g_f3cerrors);
  m_errors->Add(g_f2cerrors_sin);
  m_errors->Add(g_f3cerrors_sin);

  c_errors->cd();
  m_errors->Draw("AL*");
  leg->Draw("SAME");
  c_errors->SetLogx(1);
  c_errors->SaveAs("errors.pdf");
  
  /*Delete allocated memory*/
  
  delete[] steps;
  
  outfile.close();
  
  return 0;
}

/*
Elizabeth Drueke
PHY 480
April 29, 2016
Project 4
This is the main file for project 4.  It solves the problems posed in the 
project statement.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

#include "lattice.C"
#include "classes.C"

using namespace std;

void plots_random(int size,double temp){
  /*
    Look at how many MC cycles are needed to approach the expected analytical
    solution for a lattice of size 2 and temperature 1.
  */

  TGraph* g_meanE = new TGraph();
  TGraph* g_meanAbsM = new TGraph();
  TGraph* g_CV = new TGraph();
  TGraph* g_chi = new TGraph();
  TGraph* g_acc = new TGraph();

  int point = 0;

  for(int i=5;i<10;i++){
    lattice lat = lattice(size,temp,i); 
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=10;i<50;i+=5){
    lattice lat = lattice(size,temp,i);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=50;i<1000;i+=10){
    lattice lat = lattice(size,temp,i);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=1000;i<3000000;i+=100000){
    lattice lat = lattice(size,temp,i);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  lattice lat = lattice(size,temp,3000000);
  string name = "plots/plots_random_pofE_size"+to_string(size)+"_temp"+to_string(temp*kB);
  lat.plot_eprobs(name);

  TMultiGraph *m_meanE = new TMultiGraph("m_meanE","<E>");
  m_meanE->SetTitle("<E>;MC cycles;<E>");
  m_meanE->Add(g_meanE);
  TMultiGraph *m_meanAbsM = new TMultiGraph("m_meanAbsM","<|M|>");
  m_meanAbsM->SetTitle("<|M|>;MC cycles;<|M|>");
  m_meanAbsM->Add(g_meanAbsM);
  TMultiGraph *m_CV = new TMultiGraph("m_CV","CV");
  m_CV->SetTitle("Specific Heat;MC cycles;CV");
  m_CV->Add(g_CV);
  TMultiGraph *m_chi = new TMultiGraph("m_chi","Susceptibility");
  m_chi->SetTitle("Magnetic Susceptibility;MC cycles;#chi");
  m_chi->Add(g_chi);
  TMultiGraph *m_acc = new TMultiGraph("m_acc","Accepted MC Samples");
  m_acc->SetTitle("Accepted MC Events;MC cycles;Number of accepted events");
  m_acc->Add(g_acc);

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);
  can->cd();
  m_meanE->Draw("A*");
  can->SaveAs(("plots/plots_random_meanE_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_random_meanE_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_meanAbsM->Draw("A*");
  can->SaveAs(("plots/plots_random_meanAbsM_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_random_meanAbsM_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_CV->Draw("A*");
  can->SaveAs(("plots/plots_random_CV_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_random_CV_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_chi->Draw("A*");
  can->SaveAs(("plots/plots_random_chi_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_random_chi_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_acc->Draw("A*");
  can->SaveAs(("plots/plots_random_acceptedMC_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_random_acceptedMC_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  
  can->Close();

}
  
void plots_steady(int size,double temp){
  /*
    Look at how many MC cycles are needed to approach the expected analytical
    solution for a lattice of size 2 and temperature 1.
  */

  TGraph* g_meanE = new TGraph();
  TGraph* g_meanAbsM = new TGraph();
  TGraph* g_CV = new TGraph();
  TGraph* g_chi = new TGraph();
  TGraph* g_acc = new TGraph();

  int point = 0;

  for(int i=5;i<10;i++){
    lattice lat = lattice(size,temp,i,1); 
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=10;i<50;i+=5){
    lattice lat = lattice(size,temp,i,1);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=50;i<1000;i+=10){
    lattice lat = lattice(size,temp,i,1);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  for(int i=1000;i<3000000;i+=100000){
    lattice lat = lattice(size,temp,i,1);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    g_acc->SetPoint(point,i,lat.get_accepted());
    point++;
  }

  lattice lat = lattice(size,temp,3000000,1);
  string name = "plots/plots_steady_pofE_size"+to_string(size)+"_temp"+to_string(temp*kB);
  lat.plot_eprobs(name);

  TMultiGraph *m_meanE = new TMultiGraph("m_meanE","<E>");
  m_meanE->SetTitle("<E>;MC cycles;<E>");
  m_meanE->Add(g_meanE);
  TMultiGraph *m_meanAbsM = new TMultiGraph("m_meanAbsM","<|M|>");
  m_meanAbsM->SetTitle("<|M|>;MC cycles;<|M|>");
  m_meanAbsM->Add(g_meanAbsM);
  TMultiGraph *m_CV = new TMultiGraph("m_CV","CV");
  m_CV->SetTitle("Specific Heat;MC cycles;CV");
  m_CV->Add(g_CV);
  TMultiGraph *m_chi = new TMultiGraph("m_chi","Susceptibility");
  m_chi->SetTitle("Magnetic Susceptibility;MC cycles;#chi");
  m_chi->Add(g_chi);
  TMultiGraph *m_acc = new TMultiGraph("m_acc","Accepted MC Samples");
  m_acc->SetTitle("Accepted MC Events;MC cycles;Number of accepted events");
  m_acc->Add(g_acc);

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);
  can->cd();
  m_meanE->Draw("A*");
  can->SaveAs(("plots/plots_steady_meanE_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_steady_meanE_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_meanAbsM->Draw("A*");
  can->SaveAs(("plots/plots_steady_meanAbsM_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_steady_meanAbsM_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_CV->Draw("A*");
  can->SaveAs(("plots/plots_steady_CV_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_steady_CV_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_chi->Draw("A*");
  can->SaveAs(("plots/plots_steady_chi_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_steady_chi_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  m_acc->Draw("A*");
  can->SaveAs(("plots/plots_steady_acceptedMC_size"+to_string(size)+"_temp"+to_string(temp*kB)+".png").c_str());
  can->SaveAs(("plots/plots_steady_acceptedMC_size"+to_string(size)+"_temp"+to_string(temp*kB)+".pdf").c_str());
  
  can->Close();

}

void parte(int latsize){
  /*
    Plot the various stat mech quantities for some lattice size against 
    temperature.
  */

  TGraph *g_expE = new TGraph();
  TGraph *g_expAbsM = new TGraph();
  TGraph *g_CV = new TGraph();
  TGraph *g_chi = new TGraph();

  int ct = 0;
  for(double temp=1.0;temp<=4.0;temp+=0.01){
    lattice lat = lattice(latsize,temp/kB,3000000,1);
    g_expE->SetPoint(ct,temp,lat.get_E());
    g_expAbsM->SetPoint(ct,temp,lat.get_absM());
    g_CV->SetPoint(ct,temp,lat.get_CV());
    g_chi->SetPoint(ct,temp,lat.get_susc_absm());
    ct++;
  }

  TMultiGraph *m_expE = new TMultiGraph("m_expE","<E>");
  m_expE->SetTitle(("<E> for L="+to_string(latsize)+";Temperature (#frac{1}{k_{B}});<E>").c_str());
  m_expE->Add(g_expE);
  TMultiGraph *m_expAbsM = new TMultiGraph("m_expAbsM","<|M|>");
  m_expAbsM->SetTitle(("<|M|> for L="+to_string(latsize)+";Temperature (#frac{1}{k_{B}});<|M|>").c_str());
  m_expAbsM->Add(g_expAbsM);
  TMultiGraph *m_CV = new TMultiGraph("m_CV","C_V");
  m_CV->SetTitle(("C_{V} for L="+to_string(latsize)+";Temperature (#frac{1}{k_{B}});C_{V}").c_str());
  m_CV->Add(g_CV);
  TMultiGraph *m_chi = new TMultiGraph("m_chi","#Chi");
  m_chi->SetTitle(("#chi for L="+to_string(latsize)+";Temperature (#frac{1}{k_{B}});#Chi").c_str());
  m_chi->Add(g_chi);

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);
  can->cd();
  m_expE->Draw("A*");
  can->SaveAs(("plots/expE_latsize"+to_string(latsize)+".png").c_str());
  can->SaveAs(("plots/expE_latsize"+to_string(latsize)+".pdf").c_str());
  m_expAbsM->Draw("A*");
  can->SaveAs(("plots/expAbsM_latsize"+to_string(latsize)+".png").c_str());
  can->SaveAs(("plots/expAbsM_latsize"+to_string(latsize)+".pdf").c_str());
  m_CV->Draw("A*");
  can->SaveAs(("plots/CV_latsize"+to_string(latsize)+".png").c_str());
  can->SaveAs(("plots/CV_latsize"+to_string(latsize)+".pdf").c_str());
  m_chi->Draw("A*");
  can->SaveAs(("plots/susceptibility_latsize"+to_string(latsize)+".png").c_str());
  can->SaveAs(("plots/susceptibility_latsize"+to_string(latsize)+".pdf").c_str());
  can->Close();

}
  
void project4(){
  /*
    Main function calls all other functions in order to answer the questions 
    posed in the problem of the project.
  */

  //Look at lattice of size 2 (part b)
  //plots_random(2,1/kB); 

  //Look at lattice of size 20 (part c)
  //plots_random(20,1.0/kB); plots_steady(20,1.0/kB);
  //plots_random(20,2.4/kB); plots_steady(20,2.4/kB);

  //Look at the dependence of various stat mech quantities on temperature
  /*parte(20);
  parte(40);
  parte(60);
  parte(80);*/

  //Compare different random number generators
  plots_random(2,1/kB); 
}


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

void partb(){
  /*
    Look at how many MC cycles are needed to approach the expected analytical
    solution for a lattice of size 2 and temperature 1.
  */

  TGraph* g_meanE = new TGraph();
  TGraph* g_meanAbsM = new TGraph();
  TGraph* g_CV = new TGraph();
  TGraph* g_chi = new TGraph();

  int point = 0;

  for(int i=5;i<10;i++){
    lattice lat = lattice(2,1,i); 
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    point++;
  }
  for(int i=10;i<50;i+=5){
    lattice lat = lattice(2,1,i);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    point++;
  }
  for(int i=50;i<101;i+=10){
    lattice lat = lattice(2,1,i);
    g_meanE->SetPoint(point,i,lat.get_E());
    g_meanAbsM->SetPoint(point,i,lat.get_absM());
    g_CV->SetPoint(point,i,lat.get_CV());
    g_chi->SetPoint(point,i,lat.get_susc());
    point++;
  }

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

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);
  can->cd();
  m_meanE->Draw("AC*");
  can->SaveAs("plots/partb_meanE.png");
  can->SaveAs("plots/partb_meanE.pdf");
  m_meanAbsM->Draw("AC*");
  can->SaveAs("plots/partb_meanAbsM.png");
  can->SaveAs("plots/partb_meanAbsM.pdf");
  m_CV->Draw("AC*");
  can->SaveAs("plots/partb_CV.png");
  can->SaveAs("plots/partb_CV.pdf");
  m_chi->Draw("AC*");
  can->SaveAs("plots/partb_chi.png");
  can->SaveAs("plots/partb_chi.pdf");
  can->Close();

}
  
void project4(){
  /*
    Main function calls all other functions in order to answer the questions 
    posed in the problem of the project.
  */

  partb(); 

}


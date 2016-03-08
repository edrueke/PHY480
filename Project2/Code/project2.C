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

void solve1e(){
  /*
    Solve and plot solutions for the 1-electron problem.
  */

  double pmin = 0.0;
  vector<double> p_maxes;
  for(int i=0;i<10;i++){
    p_maxes.push_back((i+1)*10.0);
  }
  
  for(int i=0;i<p_maxes.size();i++){
    double pmax = p_maxes.at(i);
    
    //Define various plot
    //for various p_maxes
    TMultiGraph *m_wr001C = new TMultiGraph("m_wr001C",
					    ("Solution for Single Electron Case with #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr001C->SetTitle(("Solution for Single Electron Case with  #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr001C1 = new TMultiGraph("m_wr001C1",
					    ("Solution for Single Electron Case with #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr001C1->SetTitle(("Solution for Single Electron Case with  #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr001C2 = new TMultiGraph("m_wr001C2",
					    ("Solution for Single Electron Case with #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr001C2->SetTitle(("Solution for Single Electron Case with  #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TLegend *legC001 = new TLegend(0.450,0.6,1.0,0.9);
    legC001->SetFillColor(0);
    legC001->SetLineColor(0);
    legC001->SetShadowColor(0);
    legC001->SetTextSize(0.04);
    TLegend *legC0011 = new TLegend(0.450,0.6,1.0,0.9);
    legC0011->SetFillColor(0);
    legC0011->SetLineColor(0);
    legC0011->SetShadowColor(0);
    legC0011->SetTextSize(0.04);
    TLegend *legC0012 = new TLegend(0.450,0.6,1.0,0.9);
    legC0012->SetFillColor(0);
    legC0012->SetLineColor(0);
    legC0012->SetShadowColor(0);
    legC0012->SetTextSize(0.04);
    
    vector<int>nsteps;
    for(int o=10;o<110;o+=10){
      nsteps.push_back(o);
    }
    
    for(int j=0;j<nsteps.size();j++){	
      int nstep = nsteps.at(j);
      
      double h = (pmax-pmin)/nstep;
      
      themat mat = themat(nstep-1);
      
      thevec ds = thevec(nstep-1);
      thevec es = thevec(nstep-1);
      
      for(int m=1;m<nstep;m++){
	//With Coulomb interaction
	double p_m = pmin+m*h;
	double V_m = pow(p_m,2);
	double d_m = 2/(h*h)+V_m;
	double e_m = -1/(h*h);
	
	ds.point[m-1]=d_m;
	es.point[m-1]=e_m;
	
      }
      for(int k=0;k<mat.sz;k++){
	for(int l=0;l<mat.sz;l++){
	  if(k==l){
	    mat.point[k][l]=ds[k];
	    if(k<mat.sz-1){
	      mat.point[k][l+1]=es[k];
	      mat.point[k+1][l]=es[k];
	    }
	  }
	  else if(k!=l-1 && k!=l+1){
	    mat.point[k][l]=0;
	  }
	}
      }
      //Diagonalize using Householder's Algorithm from lib.cpp
      themat z = themat(nstep-1);
      for(int c=0;c<nstep-1;c++){
	for(int d=0;d<nstep-1;d++){
	  if(c==d){
	    z.point[c][d]=1;
	  }
	  else{
	    z.point[c][d]=0;
	  }
	}
      }
      
      tqli (ds.point, es.point, nstep-1, z.point);
      
      //Where are the relevant eigenvectors? (corresponding to 
      //smallest eigenvalue - ground state and 1st 2 excited states)
      int c3 = 0;
      int c31 = 0; int c32 = 0;
      double minc3 = abs(ds[0]);
      double minc31 = abs(ds[0]);
      double minc32 = abs(ds[0]);
      for(int c=0;c<ds.sz;c++){
	if(abs(ds[c])<minc3){
	  c3 = c;
	  minc3 = abs(ds[c]);
	}
      }
      if(c3==0){
	minc31 = abs(ds[1]);
	minc32 = abs(ds[1]);
	c32 = 1;
	c31 = 1;
      }
      for(int c=0;c<ds.sz;c++){
	if(abs(ds[c])<minc31 && c!=c3){
	  c31 = c;
	  minc31 = abs(ds[c]);
	}
      }
      if(c31==1 && c32==1){
	c32 = 2;
	minc32 = abs(ds[2]);
      }
      for(int c=0;c<ds.sz;c++){
	if(abs(ds[c])<minc32 && c!=c3 && c!=c31){
	  c32 = c;
	  minc32 = abs(ds[c]);
	}
      }
     
      //All the plots
      TGraph *g_c3 = new TGraph();
      TGraph *g_c31 = new TGraph();
      TGraph *g_c32 = new TGraph();
      for(int ct = -1;ct<ds.sz+1;ct++){
	if(ct==-1){
	  g_c3->SetPoint(ct+2,0,0);
	  g_c31->SetPoint(ct+2,0,0);
	  g_c32->SetPoint(ct+2,0,0);
	}
	else if(ct==ds.sz){
	  g_c3->SetPoint(ct+2,pmax,0);
	  g_c31->SetPoint(ct+2,pmax,0);
	  g_c32->SetPoint(ct+2,pmax,0);
	}
	else{
	  g_c3->SetPoint(ct+2,pmin+(ct+1)*h,z.point[ct][c3]);
	  g_c31->SetPoint(ct+2,pmin+(ct+1)*h,z.point[ct][c31]);
	  g_c32->SetPoint(ct+2,pmin+(ct+1)*h,z.point[ct][c32]);
	}
      }
      g_c3->SetLineColor(j+1);
      g_c3->SetMarkerColor(j+1);
      g_c3->SetMarkerStyle(3);
      g_c3->SetMarkerSize(1.5);
      g_c31->SetLineColor(j+1);
      g_c31->SetMarkerColor(j+1);
      g_c31->SetMarkerStyle(3);
      g_c31->SetMarkerSize(1.5);
      g_c32->SetLineColor(j+1);
      g_c32->SetMarkerColor(j+1);
      g_c32->SetMarkerStyle(3);
      g_c32->SetMarkerSize(1.5);
      gStyle->SetOptFit(0);
      
      legC001->AddEntry(g_c3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc3)).c_str());
      legC0011->AddEntry(g_c31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc31)).c_str());
      legC0012->AddEntry(g_c32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc32)).c_str());
      m_wr001C->Add(g_c3);
      m_wr001C1->Add(g_c31);
      m_wr001C2->Add(g_c32);
      cout<<"nstep = "<<nstep<<"; pmax = "<<pmax<<endl;
    }
    
    TCanvas *c_wr001NC = new TCanvas("c_wr001NC","c_wr001NC",800,720);
    c_wr001NC->SetBorderMode(0);
    
    c_wr001NC->cd();
    m_wr001C->Draw("ap");
    legC001->Draw("SAME");
    
    c_wr001NC->SaveAs(("plots/singlee"+to_string(pmax)+".pdf").c_str());
    c_wr001NC->SaveAs(("plots/singlee"+to_string(pmax)+".png").c_str());
    c_wr001NC->SaveAs(("plots/singlee"+to_string(pmax)+".root").c_str());
    c_wr001NC->Close();

    TCanvas *c_wr001NC1 = new TCanvas("c_wr001NC1","c_wr001NC1",800,720);
    c_wr001NC1->SetBorderMode(0);
    
    c_wr001NC1->cd();
    m_wr001C1->Draw("ap");
    legC0011->Draw("SAME");
    
    c_wr001NC1->SaveAs(("plots/singlee"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr001NC1->SaveAs(("plots/singlee"+to_string(pmax)+"1ststate.png").c_str());
    c_wr001NC1->SaveAs(("plots/singlee"+to_string(pmax)+"1ststate.root").c_str());
    c_wr001NC1->Close();

    TCanvas *c_wr001NC2 = new TCanvas("c_wr001NC2","c_wr001NC2",800,720);
    c_wr001NC2->SetBorderMode(0);
    
    c_wr001NC2->cd();
    m_wr001C2->Draw("ap");
    legC0012->Draw("SAME");
    
    c_wr001NC2->SaveAs(("plots/singlee"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr001NC2->SaveAs(("plots/singlee"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr001NC2->SaveAs(("plots/singlee"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr001NC2->Close();

  }
}

void partc(){
  /*
    Solve the 2-electron problem.
  */

  //Look at various p_max's.
  double pmin = 0.0;
  vector<double> p_maxes;
  for(int i=0;i<10;i++){
    p_maxes.push_back((i+1)*10.0);
  }
  
  //Look at various w_r's
  vector<double>wrs;
  wrs.push_back(0.01); wrs.push_back(0.5); 
  wrs.push_back(1); wrs.push_back(5);
  
  for(int i=0;i<p_maxes.size();i++){
    double pmax = p_maxes.at(i);
    
    //Define various plots - each wr has its own NC and C plot 
    //for various p_maxes
    TMultiGraph *m_wr001C = new TMultiGraph("m_wr001C",
					    ("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr001C->SetTitle(("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr001NC = new TMultiGraph("m_wr001C",
					     ("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr001NC->SetTitle(("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05C = new TMultiGraph("m_wr05C",
					   ("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr05C->SetTitle(("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05NC = new TMultiGraph("m_wr05C",
					    ("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr05NC->SetTitle(("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1C = new TMultiGraph("m_wr1C",
					  ("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr1C->SetTitle(("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1NC = new TMultiGraph("m_wr1C",
					   ("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr1NC->SetTitle(("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5C = new TMultiGraph("m_wr5C",
					  ("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr5C->SetTitle(("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5NC = new TMultiGraph("m_wr5C",
					   ("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State").c_str());
    m_wr5NC->SetTitle(("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - Ground State;r;#Psi(r)").c_str());

    //1st excited state
    TMultiGraph *m_wr001C1 = new TMultiGraph("m_wr001C1",
					    ("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr001C1->SetTitle(("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr001NC1 = new TMultiGraph("m_wr001C1",
					     ("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr001NC1->SetTitle(("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05C1 = new TMultiGraph("m_wr05C1",
					   ("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr05C1->SetTitle(("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05NC1 = new TMultiGraph("m_wr05C1",
					    ("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr05NC1->SetTitle(("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1C1 = new TMultiGraph("m_wr1C1",
					  ("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr1C1->SetTitle(("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1NC1 = new TMultiGraph("m_wr1C1",
					   ("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr1NC1->SetTitle(("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5C1 = new TMultiGraph("m_wr5C1",
					  ("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr5C1->SetTitle(("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5NC1 = new TMultiGraph("m_wr5C1",
					   ("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State").c_str());
    m_wr5NC1->SetTitle(("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 1st Excited State;r;#Psi(r)").c_str());

    //2nd excited state
    TMultiGraph *m_wr001C2 = new TMultiGraph("m_wr001C2",
					    ("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr001C2->SetTitle(("Solution for w_{r}=0.01 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr001NC2 = new TMultiGraph("m_wr001C2",
					     ("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr001NC2->SetTitle(("Solution for w_{r}=0.01 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05C2 = new TMultiGraph("m_wr05C2",
					   ("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr05C2->SetTitle(("Solution for w_{r}=0.5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr05NC2 = new TMultiGraph("m_wr05C2",
					    ("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr05NC2->SetTitle(("Solution for w_{r}=0.5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1C2 = new TMultiGraph("m_wr1C2",
					  ("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr1C2->SetTitle(("Solution for w_{r}=1 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr1NC2 = new TMultiGraph("m_wr1C2",
					   ("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr1NC2->SetTitle(("Solution for w_{r}=1 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5C2 = new TMultiGraph("m_wr5C2",
					  ("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr5C2->SetTitle(("Solution for w_{r}=5 with Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    TMultiGraph *m_wr5NC2 = new TMultiGraph("m_wr5C2",
					   ("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State").c_str());
    m_wr5NC2->SetTitle(("Solution for w_{r}=5 without Coulomb Repulsion #rho_{max} = "+to_string(pmax)+" - 2nd Excited State;r;#Psi(r)").c_str());
    
    TLegend *legC001 = new TLegend(0.450,0.6,1.0,0.9);
    legC001->SetFillColor(0);
    legC001->SetLineColor(0);
    legC001->SetShadowColor(0);
    legC001->SetTextSize(0.04);
    TLegend *legNC001 = new TLegend(0.450,0.6,1.0,0.9);
    legNC001->SetFillColor(0);
    legNC001->SetLineColor(0);
    legNC001->SetShadowColor(0);
    legNC001->SetTextSize(0.04);
    TLegend *legC05 = new TLegend(0.450,0.6,1.0,0.9);
    legC05->SetFillColor(0);
    legC05->SetLineColor(0);
    legC05->SetShadowColor(0);
    legC05->SetTextSize(0.04);
    TLegend *legNC05 = new TLegend(0.450,0.6,1.0,0.9);
    legNC05->SetFillColor(0);
    legNC05->SetLineColor(0);
    legNC05->SetShadowColor(0);
    legNC05->SetTextSize(0.04);
    TLegend *legC1 = new TLegend(0.450,0.6,1.0,0.9);
    legC1->SetFillColor(0);
    legC1->SetLineColor(0);
    legC1->SetShadowColor(0);
    legC1->SetTextSize(0.04);
    TLegend *legNC1 = new TLegend(0.450,0.6,1.0,0.9);
    legNC1->SetFillColor(0);
    legNC1->SetLineColor(0);
    legNC1->SetShadowColor(0);
    legNC1->SetTextSize(0.04);
    TLegend *legC5 = new TLegend(0.450,0.6,1.0,0.9);
    legC5->SetFillColor(0);
    legC5->SetLineColor(0);
    legC5->SetShadowColor(0);
    legC5->SetTextSize(0.04);
    TLegend *legNC5 = new TLegend(0.450,0.6,1.0,0.9);
    legNC5->SetFillColor(0);
    legNC5->SetLineColor(0);
    legNC5->SetShadowColor(0);
    legNC5->SetTextSize(0.04);

    //1st excited state
    TLegend *legC0011 = new TLegend(0.450,0.6,1.0,0.9);
    legC0011->SetFillColor(0);
    legC0011->SetLineColor(0);
    legC0011->SetShadowColor(0);
    legC0011->SetTextSize(0.04);
    TLegend *legNC0011 = new TLegend(0.450,0.6,1.0,0.9);
    legNC0011->SetFillColor(0);
    legNC0011->SetLineColor(0);
    legNC0011->SetShadowColor(0);
    legNC0011->SetTextSize(0.04);
    TLegend *legC051 = new TLegend(0.450,0.6,1.0,0.9);
    legC051->SetFillColor(0);
    legC051->SetLineColor(0);
    legC051->SetShadowColor(0);
    legC051->SetTextSize(0.04);
    TLegend *legNC051 = new TLegend(0.450,0.6,1.0,0.9);
    legNC051->SetFillColor(0);
    legNC051->SetLineColor(0);
    legNC051->SetShadowColor(0);
    legNC051->SetTextSize(0.04);
    TLegend *legC11 = new TLegend(0.450,0.6,1.0,0.9);
    legC11->SetFillColor(0);
    legC11->SetLineColor(0);
    legC11->SetShadowColor(0);
    legC11->SetTextSize(0.04);
    TLegend *legNC11 = new TLegend(0.450,0.6,1.0,0.9);
    legNC11->SetFillColor(0);
    legNC11->SetLineColor(0);
    legNC11->SetShadowColor(0);
    legNC11->SetTextSize(0.04);
    TLegend *legC51 = new TLegend(0.450,0.6,1.0,0.9);
    legC51->SetFillColor(0);
    legC51->SetLineColor(0);
    legC51->SetShadowColor(0);
    legC51->SetTextSize(0.04);
    TLegend *legNC51 = new TLegend(0.450,0.6,1.0,0.9);
    legNC51->SetFillColor(0);
    legNC51->SetLineColor(0);
    legNC51->SetShadowColor(0);
    legNC51->SetTextSize(0.04);

    //2nd excited state
     TLegend *legC0012 = new TLegend(0.450,0.6,1.0,0.9);
    legC0012->SetFillColor(0);
    legC0012->SetLineColor(0);
    legC0012->SetShadowColor(0);
    legC0012->SetTextSize(0.04);
    TLegend *legNC0012 = new TLegend(0.450,0.6,1.0,0.9);
    legNC0012->SetFillColor(0);
    legNC0012->SetLineColor(0);
    legNC0012->SetShadowColor(0);
    legNC0012->SetTextSize(0.04);
    TLegend *legC052 = new TLegend(0.450,0.6,1.0,0.9);
    legC052->SetFillColor(0);
    legC052->SetLineColor(0);
    legC052->SetShadowColor(0);
    legC052->SetTextSize(0.04);
    TLegend *legNC052 = new TLegend(0.450,0.6,1.0,0.9);
    legNC052->SetFillColor(0);
    legNC052->SetLineColor(0);
    legNC052->SetShadowColor(0);
    legNC052->SetTextSize(0.04);
    TLegend *legC12 = new TLegend(0.450,0.6,1.0,0.9);
    legC12->SetFillColor(0);
    legC12->SetLineColor(0);
    legC12->SetShadowColor(0);
    legC12->SetTextSize(0.04);
    TLegend *legNC12 = new TLegend(0.450,0.6,1.0,0.9);
    legNC12->SetFillColor(0);
    legNC12->SetLineColor(0);
    legNC12->SetShadowColor(0);
    legNC12->SetTextSize(0.04);
    TLegend *legC52 = new TLegend(0.450,0.6,1.0,0.9);
    legC52->SetFillColor(0);
    legC52->SetLineColor(0);
    legC52->SetShadowColor(0);
    legC52->SetTextSize(0.04);
    TLegend *legNC52 = new TLegend(0.450,0.6,1.0,0.9);
    legNC52->SetFillColor(0);
    legNC52->SetLineColor(0);
    legNC52->SetShadowColor(0);
    legNC52->SetTextSize(0.04);
   
    for(int r=0;r<wrs.size();r++){
      double wr = wrs.at(r);
      
      vector<int>nsteps;
      for(int o=10;o<110;o+=10){
	nsteps.push_back(o);
      }
      
      for(int j=0;j<nsteps.size();j++){	
	int nstep = nsteps.at(j);
	
	double h = (pmax-pmin)/nstep;
	
	themat matC = themat(nstep-1);
	themat matNC = themat(nstep-1);
	
	thevec dsC = thevec(nstep-1);
	thevec esC = thevec(nstep-1);
	thevec dsNC = thevec(nstep-1);
	thevec esNC = thevec(nstep-1);
	
	for(int m=1;m<nstep;m++){
	  //With Coulomb interaction
	  double p_m = pmin+m*h;
	  double V_m = pow(wr*p_m,2)+1/p_m;
	  double d_m = 2/(h*h)+V_m;
	  double e_m = -1/(h*h);
	  
	  dsC.point[m-1]=d_m;
	  esC.point[m-1]=e_m;
	  
	  //Without Coulomb interaction
	  V_m = pow(wr*p_m,2);
	  d_m=2/(h*h)+V_m;
	  
	  dsNC.point[m-1]=d_m;
	  esNC.point[m-1]=e_m;
	  
	}
	for(int k=0;k<matC.sz;k++){
	  for(int l=0;l<matC.sz;l++){
	    if(k==l){
	      matC.point[k][l]=dsC[k];
	      matNC.point[k][l]=dsNC[k];
	      if(k<matC.sz-1){
		matC.point[k][l+1]=esC[k];
		matC.point[k+1][l]=esC[k];
		matNC.point[k][l+1]=esNC[k];
		matNC.point[k+1][l]=esNC[k];
	      }
	    }
	    else if(k!=l-1 && k!=l+1){
	      matC.point[k][l]=0;
	      matNC.point[k][l]=0;
	    }
	  }
	}
	
	//Diagonalize using Householder's Algorithm from lib.cpp
	themat zC = themat(nstep-1);
	themat zNC = themat(nstep-1);
	for(int c=0;c<nstep-1;c++){
	  for(int d=0;d<nstep-1;d++){
	    if(c==d){
	      zC.point[c][d]=1;
	      zNC.point[c][d]=1;
	    }
	    else{
	      zC.point[c][d]=0;
	      zNC.point[c][d]=0;
	    }
	  }
	}
	
	tqli (dsC.point, esC.point, nstep-1, zC.point);
	tqli (dsNC.point, esNC.point, nstep-1, zNC.point);
	
	//Where are the relevant eigenvectors? (corresponding to 
	//smallest eigenvalue - ground state)
	int c3 = 0;
	double minc3 = abs(dsC[0]);
	int nc3 = 0;
	double minnc3 = abs(dsNC[0]);
	int c31 = 0;
	double minc31 = abs(dsC[0]);
	int nc31 = 0;
	double minnc31 = abs(dsNC[0]);
	int c32 = 0;
	double minc32 = abs(dsC[0]);
	int nc32 = 0;
	double minnc32 = abs(dsNC[0]);
	
	for(int c=0;c<dsC.sz;c++){
	  if(abs(dsC[c])<minc3){
	    c3 = c;
	    minc3 = abs(dsC[c]);
	  }
	  if(abs(dsNC[c])<minnc3){
	    nc3 = c;
	    minnc3 = abs(dsNC[c]);
	  }
	}
	if(c3==0){
	  c31=1;
	  minc31 = abs(dsC[1]);
	  c32=1;
	  minc32 = abs(dsC[1]);
	}
	if(nc3==0){
	  nc31=1;
	  minnc31 = abs(dsNC[1]);
	  c32=1;
	  minnc32 = abs(dsNC[1]);
	}
	for(int c=0;c<dsC.sz;c++){
	  if(abs(dsC[c])<minc31 && c!=c3){
	    c31 = c;
	    minc31 = abs(dsC[c]);
	  }
	  if(abs(dsNC[c])<minnc31 && c!=nc3){
	    nc31 = c;
	    minnc31 = abs(dsNC[c]);
	  }
	}
	if(c31==1){
	  c32=2;
	  minc32 = abs(dsC[2]);
	}
	if(nc31==1){
	  c32=2;
	  minnc32 = abs(dsNC[2]);
	}
	for(int c=0;c<dsC.sz;c++){
	  if(abs(dsC[c])<minc32 && c!=c3 && c!=c31){
	    c32 = c;
	    minc32 = abs(dsC[c]);
	  }
	  if(abs(dsNC[c])<minnc31 && c!=nc3 && c!=nc31){
	    nc32 = c;
	    minnc32 = abs(dsNC[c]);
	  }
	}
	
	//All the plots
	TGraph *g_c3 = new TGraph();
	TGraph *g_nc3 = new TGraph();
	TGraph *g_c31 = new TGraph();
	TGraph *g_nc31 = new TGraph();
	TGraph *g_c32 = new TGraph();
	TGraph *g_nc32 = new TGraph();
	for(int ct = -1;ct<dsC.sz+1;ct++){
	  if(ct==-1){
	    g_c3->SetPoint(ct+2,0,0);
	    g_nc3->SetPoint(ct+2,0,0);
	    g_c31->SetPoint(ct+2,0,0);
	    g_nc31->SetPoint(ct+2,0,0);
	    g_c32->SetPoint(ct+2,0,0);
	    g_nc32->SetPoint(ct+2,0,0);
	  }
	  else if(ct==dsC.sz){
	    g_c3->SetPoint(ct+2,pmax,0);
	    g_nc3->SetPoint(ct+2,pmax,0);
	    g_c31->SetPoint(ct+2,pmax,0);
	    g_nc31->SetPoint(ct+2,pmax,0);
	    g_c32->SetPoint(ct+2,pmax,0);
	    g_nc32->SetPoint(ct+2,pmax,0);
	  }
	  else{
	    g_c3->SetPoint(ct+2,pmin+(ct+1)*h,zC.point[ct][c3]);
	    g_nc3->SetPoint(ct+2,pmin+(ct+1)*h,zNC[ct][nc3]);
	    g_c31->SetPoint(ct+2,pmin+(ct+1)*h,zC.point[ct][c31]);
	    g_nc31->SetPoint(ct+2,pmin+(ct+1)*h,zNC[ct][nc31]);
	    g_c32->SetPoint(ct+2,pmin+(ct+1)*h,zC.point[ct][c32]);
	    g_nc32->SetPoint(ct+2,pmin+(ct+1)*h,zNC[ct][nc32]);
	  }
	}
	g_c3->SetLineColor(j+1);
	g_c3->SetMarkerColor(j+1);
	g_c3->SetMarkerStyle(3);
	g_c3->SetMarkerSize(1.5);
	g_nc3->SetLineColor(j+1);
	g_nc3->SetMarkerColor(j+1);
	g_nc3->SetMarkerStyle(3);
	g_nc3->SetMarkerSize(1.5);
	g_c31->SetLineColor(j+1);
	g_c31->SetMarkerColor(j+1);
	g_c31->SetMarkerStyle(3);
	g_c31->SetMarkerSize(1.5);
	g_nc31->SetLineColor(j+1);
	g_nc31->SetMarkerColor(j+1);
	g_nc31->SetMarkerStyle(3);
	g_nc31->SetMarkerSize(1.5);
	g_c32->SetLineColor(j+1);
	g_c32->SetMarkerColor(j+1);
	g_c32->SetMarkerStyle(3);
	g_c32->SetMarkerSize(1.5);
	g_nc32->SetLineColor(j+1);
	g_nc32->SetMarkerColor(j+1);
	g_nc32->SetMarkerStyle(3);
	g_nc32->SetMarkerSize(1.5);
	gStyle->SetOptFit(0);
	
	if(wr==0.01){
	  legC001->AddEntry(g_c3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc3)).c_str());
	  legNC001->AddEntry(g_nc3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc3)).c_str());
	  m_wr001C->Add(g_c3);
	  m_wr001NC->Add(g_nc3);
	  legC0011->AddEntry(g_c31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc31)).c_str());
	  legNC0011->AddEntry(g_nc31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc31)).c_str());
	  m_wr001C1->Add(g_c31);
	  m_wr001NC1->Add(g_nc31);
	  legC0012->AddEntry(g_c32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc32)).c_str());
	  legNC0012->AddEntry(g_nc32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc32)).c_str());
	  m_wr001C2->Add(g_c32);
	  m_wr001NC2->Add(g_nc32);
	}
	else if(wr==0.5){
	  legC05->AddEntry(g_c3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc3)).c_str());
	  legNC05->AddEntry(g_nc3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc3)).c_str());
	  m_wr05C->Add(g_c3);
	  m_wr05NC->Add(g_nc3);
	  legC051->AddEntry(g_c31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc31)).c_str());
	  legNC051->AddEntry(g_nc31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc31)).c_str());
	  m_wr05C1->Add(g_c31);
	  m_wr05NC1->Add(g_nc31);
	  legC052->AddEntry(g_c32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc32)).c_str());
	  legNC052->AddEntry(g_nc32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc32)).c_str());
	  m_wr05C2->Add(g_c32);
	  m_wr05NC2->Add(g_nc32);
	}
	else if(wr==1){
	  legC1->AddEntry(g_c3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc3)).c_str());
	  legNC1->AddEntry(g_nc3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc3)).c_str());
	  m_wr1C->Add(g_c3);
	  m_wr1NC->Add(g_nc3);
	  legC11->AddEntry(g_c31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc31)).c_str());
	  legNC11->AddEntry(g_nc31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc31)).c_str());
	  m_wr1C1->Add(g_c31);
	  m_wr1NC1->Add(g_nc31);
	  legC12->AddEntry(g_c32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc32)).c_str());
	  legNC12->AddEntry(g_nc32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc32)).c_str());
	  m_wr1C2->Add(g_c32);
	  m_wr1NC2->Add(g_nc32);
	}
	else{
	  legC5->AddEntry(g_c3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc3)).c_str());
	  legNC5->AddEntry(g_nc3,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc3)).c_str());
	  m_wr5C->Add(g_c3);
	  m_wr5NC->Add(g_nc3);
	  legC51->AddEntry(g_c31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc31)).c_str());
	  legNC51->AddEntry(g_nc31,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc31)).c_str());
	  m_wr5C1->Add(g_c31);
	  m_wr5NC1->Add(g_nc31);
	  legC52->AddEntry(g_c32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minc32)).c_str());
	  legNC52->AddEntry(g_nc32,("nsteps = "+to_string(nstep)+"; #lambda = "+to_string(minnc32)).c_str());
	  m_wr5C2->Add(g_c32);
	  m_wr5NC2->Add(g_nc32);
	}
	
	cout<<"nstep = "<<nstep<<"; pmax = "<<pmax<<"; wr = "<<wr<<endl;
      }
    }

    TCanvas *c_wr001NC = new TCanvas("c_wr001NC","c_wr001NC",800,720);
    c_wr001NC->SetBorderMode(0);
    
    c_wr001NC->cd();
    m_wr001NC->Draw("ap");
    legNC001->Draw("SAME");
    
    c_wr001NC->SaveAs(("plots/wr001NC"+to_string(pmax)+".pdf").c_str());
    c_wr001NC->SaveAs(("plots/wr001NC"+to_string(pmax)+".png").c_str());
    c_wr001NC->SaveAs(("plots/wr001NC"+to_string(pmax)+".root").c_str());
    c_wr001NC->Close();
    
    TCanvas *c_wr05NC = new TCanvas("c_wr05NC","c_wr05NC",800,720);
    c_wr05NC->SetBorderMode(0);
    
    c_wr05NC->cd();
    m_wr05NC->Draw("ap");
    legNC05->Draw("SAME");
    
    c_wr05NC->SaveAs(("plots/wr05NC"+to_string(pmax)+".pdf").c_str());
    c_wr05NC->SaveAs(("plots/wr05NC"+to_string(pmax)+".png").c_str());
    c_wr05NC->SaveAs(("plots/wr05NC"+to_string(pmax)+".root").c_str());
    c_wr05NC->Close();
    
    TCanvas *c_wr1NC = new TCanvas("c_wr1NC","c_wr1NC",800,720);
    c_wr1NC->SetBorderMode(0);
    
    c_wr1NC->cd();
    m_wr1NC->Draw("ap");
    legNC1->Draw("SAME");
    
    c_wr1NC->SaveAs(("plots/wr1NC"+to_string(pmax)+".pdf").c_str());
    c_wr1NC->SaveAs(("plots/wr1NC"+to_string(pmax)+".png").c_str());
    c_wr1NC->SaveAs(("plots/wr1NC"+to_string(pmax)+".root").c_str());
    c_wr1NC->Close();
    
    TCanvas *c_wr5NC = new TCanvas("c_wr5NC","c_wr5NC",800,720);
    c_wr5NC->SetBorderMode(0);
    
    c_wr5NC->cd();
    m_wr5NC->Draw("ap");
    legNC5->Draw("SAME");
    
    c_wr5NC->SaveAs(("plots/wr5NC"+to_string(pmax)+".pdf").c_str());
    c_wr5NC->SaveAs(("plots/wr5NC"+to_string(pmax)+".png").c_str());
    c_wr5NC->SaveAs(("plots/wr5NC"+to_string(pmax)+".root").c_str());
    c_wr5NC->Close();
    
    TCanvas *c_wr001C = new TCanvas("c_wr001C","c_wr001C",800,720);
    c_wr001C->SetBorderMode(0);
    
    c_wr001C->cd();
    m_wr001C->Draw("ap");
    legC001->Draw("SAME");
    
    c_wr001C->SaveAs(("plots/wr001C"+to_string(pmax)+".pdf").c_str());
    c_wr001C->SaveAs(("plots/wr001C"+to_string(pmax)+".png").c_str());
    c_wr001C->SaveAs(("plots/wr001C"+to_string(pmax)+".root").c_str());
    c_wr001C->Close();
    
    TCanvas *c_wr05C = new TCanvas("c_wr05C","c_wr05C",800,720);
    c_wr05C->SetBorderMode(0);
    
    c_wr05C->cd();
    m_wr05C->Draw("ap");
    legC05->Draw("SAME");
    
    c_wr05C->SaveAs(("plots/wr05C"+to_string(pmax)+".pdf").c_str());
    c_wr05C->SaveAs(("plots/wr05C"+to_string(pmax)+".png").c_str());
    c_wr05C->SaveAs(("plots/wr05C"+to_string(pmax)+".root").c_str());
    c_wr05C->Close();
    
    TCanvas *c_wr1C = new TCanvas("c_wr1C","c_wr1C",800,720);
    c_wr1C->SetBorderMode(0);
    
    c_wr1C->cd();
    m_wr1C->Draw("ap");
    legC1->Draw("SAME");
    
    c_wr1C->SaveAs(("plots/wr1C"+to_string(pmax)+".pdf").c_str());
    c_wr1C->SaveAs(("plots/wr1C"+to_string(pmax)+".png").c_str());
    c_wr1C->SaveAs(("plots/wr1C"+to_string(pmax)+".root").c_str());
    c_wr1C->Close();
    
    TCanvas *c_wr5C = new TCanvas("c_wr5C","c_wr5C",800,720);
    c_wr5C->SetBorderMode(0);
    
    c_wr5C->cd();
    m_wr5C->Draw("ap");
    legC5->Draw("SAME");
    
    c_wr5C->SaveAs(("plots/wr5C"+to_string(pmax)+".pdf").c_str());
    c_wr5C->SaveAs(("plots/wr5C"+to_string(pmax)+".png").c_str());
    c_wr5C->SaveAs(("plots/wr5C"+to_string(pmax)+".root").c_str());
    c_wr5C->Close();
    
    //1st excited state
    TCanvas *c_wr001NC1 = new TCanvas("c_wr001NC1","c_wr001NC1",800,720);
    c_wr001NC1->SetBorderMode(0);
    
    c_wr001NC1->cd();
    m_wr001NC1->Draw("ap");
    legNC0011->Draw("SAME");
    
    c_wr001NC1->SaveAs(("plots/wr001NC"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr001NC1->SaveAs(("plots/wr001NC"+to_string(pmax)+"1ststate.png").c_str());
    c_wr001NC1->SaveAs(("plots/wr001NC"+to_string(pmax)+"1ststate.root").c_str());
    c_wr001NC1->Close();
    
    TCanvas *c_wr05NC1 = new TCanvas("c_wr05NC1","c_wr05NC1",800,720);
    c_wr05NC1->SetBorderMode(0);
    
    c_wr05NC1->cd();
    m_wr05NC1->Draw("ap");
    legNC051->Draw("SAME");
    
    c_wr05NC1->SaveAs(("plots/wr05NC"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr05NC1->SaveAs(("plots/wr05NC"+to_string(pmax)+"1ststate.png").c_str());
    c_wr05NC1->SaveAs(("plots/wr05NC"+to_string(pmax)+"1ststate.root").c_str());
    c_wr05NC1->Close();
    
    TCanvas *c_wr1NC1 = new TCanvas("c_wr1NC1","c_wr1NC1",800,720);
    c_wr1NC1->SetBorderMode(0);
    
    c_wr1NC1->cd();
    m_wr1NC1->Draw("ap");
    legNC11->Draw("SAME");
    
    c_wr1NC1->SaveAs(("plots/wr1NC"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr1NC1->SaveAs(("plots/wr1NC"+to_string(pmax)+"1ststate.png").c_str());
    c_wr1NC1->SaveAs(("plots/wr1NC"+to_string(pmax)+"1ststate.root").c_str());
    c_wr1NC1->Close();
    
    TCanvas *c_wr5NC1 = new TCanvas("c_wr5NC1","c_wr5NC1",800,720);
    c_wr5NC1->SetBorderMode(0);
    
    c_wr5NC1->cd();
    m_wr5NC1->Draw("ap");
    legNC51->Draw("SAME");
    
    c_wr5NC1->SaveAs(("plots/wr5NC"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr5NC1->SaveAs(("plots/wr5NC"+to_string(pmax)+"1ststate.png").c_str());
    c_wr5NC1->SaveAs(("plots/wr5NC"+to_string(pmax)+"1ststate.root").c_str());
    c_wr5NC1->Close();

    TCanvas *c_wr001C1 = new TCanvas("c_wr001C1","c_wr001C1",800,720);
    c_wr001C1->SetBorderMode(0);
    
    c_wr001C1->cd();
    m_wr001C1->Draw("ap");
    legC0011->Draw("SAME");
    
    c_wr001C1->SaveAs(("plots/wr001C"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr001C1->SaveAs(("plots/wr001C"+to_string(pmax)+"1ststate.png").c_str());
    c_wr001C1->SaveAs(("plots/wr001C"+to_string(pmax)+"1ststate.root").c_str());
    c_wr001C1->Close();
    
    TCanvas *c_wr05C1 = new TCanvas("c_wr05C1","c_wr05C1",800,720);
    c_wr05C1->SetBorderMode(0);
    
    c_wr05C1->cd();
    m_wr05C1->Draw("ap");
    legC051->Draw("SAME");
    
    c_wr05C1->SaveAs(("plots/wr05C"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr05C1->SaveAs(("plots/wr05C"+to_string(pmax)+"1ststate.png").c_str());
    c_wr05C1->SaveAs(("plots/wr05C"+to_string(pmax)+"1ststate.root").c_str());
    c_wr05C1->Close();
    
    TCanvas *c_wr1C1 = new TCanvas("c_wr1C1","c_wr1C1",800,720);
    c_wr1C1->SetBorderMode(0);
    
    c_wr1C1->cd();
    m_wr1C1->Draw("ap");
    legC11->Draw("SAME");
    
    c_wr1C1->SaveAs(("plots/wr1C"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr1C1->SaveAs(("plots/wr1C"+to_string(pmax)+"1ststate.png").c_str());
    c_wr1C1->SaveAs(("plots/wr1C"+to_string(pmax)+"1ststate.root").c_str());
    c_wr1C1->Close();
    
    TCanvas *c_wr5C1 = new TCanvas("c_wr5C1","c_wr5C1",800,720);
    c_wr5C1->SetBorderMode(0);
    
    c_wr5C1->cd();
    m_wr5C1->Draw("ap");
    legC51->Draw("SAME");
    
    c_wr5C1->SaveAs(("plots/wr5C"+to_string(pmax)+"1ststate.pdf").c_str());
    c_wr5C1->SaveAs(("plots/wr5C"+to_string(pmax)+"1ststate.png").c_str());
    c_wr5C1->SaveAs(("plots/wr5C"+to_string(pmax)+"1ststate.root").c_str());
    c_wr5C1->Close();
    
    //2nd excited state
    TCanvas *c_wr001NC2 = new TCanvas("c_wr001NC2","c_wr001NC2",800,720);
    c_wr001NC2->SetBorderMode(0);
    
    c_wr001NC2->cd();
    m_wr001NC2->Draw("ap");
    legNC0012->Draw("SAME");
    
    c_wr001NC2->SaveAs(("plots/wr001NC"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr001NC2->SaveAs(("plots/wr001NC"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr001NC2->SaveAs(("plots/wr001NC"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr001NC2->Close();
    
    TCanvas *c_wr05NC2 = new TCanvas("c_wr05NC2","c_wr05NC2",800,720);
    c_wr05NC2->SetBorderMode(0);
    
    c_wr05NC2->cd();
    m_wr05NC2->Draw("ap");
    legNC052->Draw("SAME");
    
    c_wr05NC2->SaveAs(("plots/wr05NC"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr05NC2->SaveAs(("plots/wr05NC"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr05NC2->SaveAs(("plots/wr05NC"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr05NC2->Close();
    
    TCanvas *c_wr1NC2 = new TCanvas("c_wr1NC2","c_wr1NC2",800,720);
    c_wr1NC2->SetBorderMode(0);
    
    c_wr1NC2->cd();
    m_wr1NC2->Draw("ap");
    legNC12->Draw("SAME");
    
    c_wr1NC2->SaveAs(("plots/wr1NC"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr1NC2->SaveAs(("plots/wr1NC"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr1NC2->SaveAs(("plots/wr1NC"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr1NC2->Close();
    
    TCanvas *c_wr5NC2 = new TCanvas("c_wr5NC2","c_wr5NC2",800,720);
    c_wr5NC2->SetBorderMode(0);
    
    c_wr5NC2->cd();
    m_wr5NC2->Draw("ap");
    legNC52->Draw("SAME");
    
    c_wr5NC2->SaveAs(("plots/wr5NC"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr5NC2->SaveAs(("plots/wr5NC"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr5NC2->SaveAs(("plots/wr5NC"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr5NC2->Close();
    TCanvas *c_wr001C2 = new TCanvas("c_wr001C2","c_wr001C2",800,720);
    c_wr001C2->SetBorderMode(0);
    
    c_wr001C2->cd();
    m_wr001C2->Draw("ap");
    legC0012->Draw("SAME");
    
    c_wr001C2->SaveAs(("plots/wr001C"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr001C2->SaveAs(("plots/wr001C"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr001C2->SaveAs(("plots/wr001C"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr001C2->Close();
    
    TCanvas *c_wr05C2 = new TCanvas("c_wr05C2","c_wr05C2",800,720);
    c_wr05C2->SetBorderMode(0);
    
    c_wr05C2->cd();
    m_wr05C2->Draw("ap");
    legC052->Draw("SAME");
    
    c_wr05C2->SaveAs(("plots/wr05C"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr05C2->SaveAs(("plots/wr05C"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr05C2->SaveAs(("plots/wr05C"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr05C2->Close();
    
    TCanvas *c_wr1C2 = new TCanvas("c_wr1C2","c_wr1C2",800,720);
    c_wr1C2->SetBorderMode(0);
    
    c_wr1C2->cd();
    m_wr1C2->Draw("ap");
    legC12->Draw("SAME");
    
    c_wr1C2->SaveAs(("plots/wr1C"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr1C2->SaveAs(("plots/wr1C"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr1C2->SaveAs(("plots/wr1C"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr1C2->Close();
    
    TCanvas *c_wr5C2 = new TCanvas("c_wr5C2","c_wr5C2",800,720);
    c_wr5C2->SetBorderMode(0);
    
    c_wr5C2->cd();
    m_wr5C2->Draw("ap");
    legC52->Draw("SAME");
    
    c_wr5C2->SaveAs(("plots/wr5C"+to_string(pmax)+"2ndstate.pdf").c_str());
    c_wr5C2->SaveAs(("plots/wr5C"+to_string(pmax)+"2ndstate.png").c_str());
    c_wr5C2->SaveAs(("plots/wr5C"+to_string(pmax)+"2ndstate.root").c_str());
    c_wr5C2->Close();
    

  }
}

void partb(){
  /*
    Solve the 1-electron problem.
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
  for(int i=0;i<5;i++){
    p_maxes.push_back((i+1)*10.0);
  }

  //Create a graph of the number of rotations required for a given number of steps.
  TGraph *g_stepvrot = new TGraph();//p_maxes.size()*n_steps.size());
  int pt = 0;

  //Create a graph of the time required for the Householder and Jacobi algorithms
  TGraph *g_timehous = new TGraph();//n_steps.size());
  TGraph *g_timejac = new TGraph();//n_steps.size());

  for(int i=0;i<p_maxes.size();i++){
    vector<int> nsteps;
    for(int j=5;j<25;j++)
      nsteps.push_back(j);
    for(int j=25;j<50;j+=5)
      nsteps.push_back(j);
    for(int j=50;j<101;j+=25)
      nsteps.push_back(j);
    
    double pmax = p_maxes.at(i);

    //Create a graph of how many steps were required before the eigenvalues converged.
    /*TGraph *g_jconverge3 = new TGraph();
    TGraph *g_jconverge7 = new TGraph();
    TGraph *g_jconverge11 = new TGraph();
    TGraph *g_hconverge3 = new TGraph();
    TGraph *g_hconverge7 = new TGraph();
    TGraph *g_hconverge11 = new TGraph();*/
    
    for(int mi=0;mi<nsteps.size();mi++){
      int nstep = nsteps.at(mi);
      int j=0;

       double h = (pmax-pmin)/nstep;
      
      themat mat = themat(nstep-1);
      
      thevec ds = thevec(nstep-1);
      thevec es = thevec(nstep-1);
      
      for(int m=1;m<nstep;m++){
	//With Coulomb interaction
	double p_m = pmin+m*h;
	double V_m = pow(p_m,2);
	double d_m = 2/(h*h)+V_m;
	double e_m = -1/(h*h);
	
	ds.point[m-1]=d_m;
	es.point[m-1]=e_m;
	
      }
      for(int k=0;k<mat.sz;k++){
	for(int l=0;l<mat.sz;l++){
	  if(k==l){
	    mat.point[k][l]=ds[k];
	    if(k<mat.sz-1){
	      mat.point[k][l+1]=es[k];
	      mat.point[k+1][l]=es[k];
	    }
	  }
	  else if(k!=l-1 && k!=l+1){
	    mat.point[k][l]=0;
	  }
	}
      }
 
     //How long does the Jacobi algorithm take?
      clock_t start_jac, finish_jac;
      start_jac = clock();
      
      //Diagonalize using the Jacobi Rotation Algorithm
      themat diag = Jacobi_Method(mat,1e-8,g_stepvrot,&pt);
      finish_jac = clock();
      
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

      /*g_jconverge3->SetPoint(i+1,nstep,1);
      g_jconverge3->SetPoint(i+1,nstep,1);
      g_jconverge11->SetPoint(i+1,nstep,1);
      g_hconverge3->SetPoint(i+1,nstep,1);
      g_hconverge7->SetPoint(i+1,nstep,1);
      g_hconverge11->SetPoint(i+1,nstep,1);
     for(int m=0;m<nstep-1;m++){
	if(abs(diag[m][m]-3)<2){
	  cout<<"entered 3j"<<endl;
	  g_jconverge3->RemovePoint(i+1);
	  g_jconverge3->SetPoint(i+1,nstep,abs(diag[m][m]-3));
	}
	if(abs(diag[m][m]-7)<2){
	  cout<<"entered 7j"<<endl;
	  g_jconverge7->RemovePoint(i+1);
	  g_jconverge7->SetPoint(i+1,nstep,abs(diag[m][m]-7));
	}
	if(abs(diag[m][m]-11)<2){
	  cout<<"entered 11j"<<endl;
	  g_jconverge11->RemovePoint(i+1);
	  g_jconverge11->SetPoint(i+1,nstep,abs(diag[m][m]-11));
	}
	eigen_jac.point[m] = diag[m][m];
      }
      
      for(int m=0;m<ds.sz;m++){
	if(abs(ds[m]-3)<2){
	  cout<<"entered 3h"<<endl;
	  g_hconverge3->RemovePoint(i+1);
	  g_hconverge3->SetPoint(i+1,nstep,abs(ds[m]-3));
	}
	else if(abs(ds[m]-7)<2){
	  cout<<"entered 7h"<<endl;
	  g_hconverge7->RemovePoint(i+1);
	  g_hconverge7->SetPoint(i+1,nstep,abs(ds[m]-7));
	}
	else if(abs(ds[m]-11)<2){
	  cout<<"entered 11h"<<endl;
	  g_hconverge11->RemovePoint(i+1);
	  g_hconverge11->SetPoint(i+1,nstep,abs(ds[m]-11));
	}
	}*/
      
      cout<<"pmax = "<<pmax<<"; nstep = "<<nstep<<endl;//<<"; eigenvalues:"<<endl;
      j+=1;
    }

    //Format and save the plot of difference in nstep
    //before jacobi and householder converge
    /*gStyle->SetOptFit(0);
    g_jconverge3->SetLineColor(kRed);
    g_jconverge3->SetMarkerColor(kRed);
    g_jconverge3->SetMarkerStyle(3);
    g_jconverge3->SetMarkerSize(1.5);
    g_jconverge7->SetLineColor(kRed);
    g_jconverge7->SetMarkerColor(kRed);
    g_jconverge7->SetMarkerStyle(3);
    g_jconverge7->SetMarkerSize(1.5);
    g_jconverge11->SetLineColor(kRed);
    g_jconverge11->SetMarkerColor(kRed);
    g_jconverge11->SetMarkerStyle(3);
    g_jconverge11->SetMarkerSize(1.5);
    
    g_hconverge3->SetLineColor(kBlue);
    g_hconverge3->SetMarkerColor(kBlue);
    g_hconverge3->SetMarkerStyle(3);
    g_hconverge3->SetMarkerSize(1.5);
    g_hconverge7->SetLineColor(kBlue);
    g_hconverge7->SetMarkerColor(kBlue);
    g_hconverge7->SetMarkerStyle(3);
    g_hconverge7->SetMarkerSize(1.5);
    g_hconverge11->SetLineColor(kBlue);
    g_hconverge11->SetMarkerColor(kBlue);
    g_hconverge11->SetMarkerStyle(3);
    g_hconverge11->SetMarkerSize(1.5);
    
    TLegend *leg = new TLegend(0.20,0.4,0.5,0.6);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(g_jconverge3,"Jacobi Algorithm","L");
    leg->AddEntry(g_hconverge3,"Householder Algorithm","L");

    TMultiGraph *m_converge3 = new TMultiGraph("m_time",
					       ("Steps Until Convergence to 3 for #rho_{max} = "+to_string(pmax)).c_str());
    m_converge3->SetTitle(("Steps Until Convergence to 3 for #rho_{max} = "+to_string(pmax)+";n_{steps};Difference").c_str());
    m_converge3->Add(g_hconverge3);
    m_converge3->Add(g_jconverge3);
    
    TCanvas *c_converge3 = new TCanvas("c_converge3","c_converge3",800,720);
    c_converge3->SetBorderMode(0);
    c_converge3->cd();
    m_converge3->Draw("A*");
    leg->Draw("SAME");
    
    c_converge3->SaveAs(("plots/algconverge3"+to_string(pmax)+".pdf").c_str());
    c_converge3->SaveAs(("plots/algconverge3"+to_string(pmax)+".png").c_str());
    c_converge3->SaveAs(("plots/algconverge3"+to_string(pmax)+".root").c_str());
    c_converge3->Close();
    
    TMultiGraph *m_converge7 = new TMultiGraph("m_time",
					       ("Steps Until Convergence to 7 for #rho_{max} = "+to_string(pmax)).c_str());
    m_converge7->SetTitle(("Steps Until Convergence to 7 for #rho_{max} = "+to_string(pmax)+";n_{steps};Difference").c_str());
    m_converge7->Add(g_hconverge7);
    m_converge7->Add(g_jconverge7);
    
    TCanvas *c_converge7 = new TCanvas("c_converge7","c_converge7",800,720);
    c_converge7->SetBorderMode(0);
    c_converge7->cd();
    m_converge7->Draw("A*");
    leg->Draw("SAME");
    
    c_converge7->SaveAs(("plots/algconverge7"+to_string(pmax)+".pdf").c_str());
    c_converge7->SaveAs(("plots/algconverge7"+to_string(pmax)+".png").c_str());
    c_converge7->SaveAs(("plots/algconverge7"+to_string(pmax)+".root").c_str());
    c_converge7->Close();
    
    TMultiGraph *m_converge11 = new TMultiGraph("m_time",
						("Steps Until Convergence to 11 for #rho_{max} = "+to_string(pmax)).c_str());
    m_converge11->SetTitle(("Steps Until Convergence to 11 for #rho_{max} = "+to_string(pmax)+";n_{steps};Difference").c_str());
    m_converge11->Add(g_hconverge11);
    m_converge11->Add(g_jconverge11);
    
    TCanvas *c_converge11 = new TCanvas("c_converge11","c_converge11",800,720);
    c_converge11->SetBorderMode(0);
    c_converge11->cd();
    m_converge11->Draw("A*");
    leg->Draw("SAME");
    
    c_converge11->SaveAs(("plots/algconverge11"+to_string(pmax)+".pdf").c_str());
    c_converge11->SaveAs(("plots/algconverge11"+to_string(pmax)+".png").c_str());
    c_converge11->SaveAs(("plots/algconverge11"+to_string(pmax)+".root").c_str());
    c_converge11->Close();*/
    
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

  partb(); //Solve the single electron problem. Focus on efficiency analysis.
  partc(); //Solve the 2-electron problem.
  solve1e(); //Solve the single electron problem.
}

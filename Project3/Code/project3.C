/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the main file for project 3.  It goes through every part listed in 
the project requirements and uses the planets and ode solvers from the 
included files (as well as the vecter and matrix classes from projects 1 and 2)
to solve the problem of the solar system.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

#include "planets.C"
#include "odesolvers.C"
#include "classes.C"

using namespace std;

void parta(){
  /*
    Using the discretized versions of the x- and y- components of the 
    velocity and position of the planet Earth (in the single-planet version
    of the solar system), this function calls on the Verlet and RK4 methods 
    to solve the problem.
  */

  gStyle->SetOptFit();

  //Define the earth as a planet.
  planet earth = planet("earth",6e24,1);
  
  //Compute the expected velocity given a circular orbit.
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);

  //cout<<"v = "<<v<<endl;

  //Define required quantities for Verlet method
  double t0 = 0; //Initial time (years)
  double tf = 1; //Final time (years)
  int nsteps = 100; //# of steps

  //Compute the x position and velocity of the earth over 1 year (1 revolution)
  //Assume the planet starts moving entirely in the x-direction.

  //First using Verlet
  vector<thevec> xcomp_v = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp_v = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos_v = xcomp_v.at(0);
  thevec xvel_v = xcomp_v.at(1);
  thevec ypos_v = ycomp_v.at(0);
  thevec yvel_v = ycomp_v.at(1);

  //cout<<"xpos: "<<xpos.print()<<endl<<"ypos: "<<ypos.print()<<endl;
  //cout<<"xvel: "<<xvel.print()<<endl<<"yvel: "<<yvel.print()<<endl;
  
  //Plot the position and velocity
  TGraph *g_pos_v = new TGraph();
  TGraph *g_vel_vx = new TGraph(nsteps+1);
  TGraph *g_vel_vy = new TGraph(nsteps+1);

  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  for(unsigned int i=0;i<nsteps+1;i++){
    double x = xpos_v[i];
    double y = ypos_v[i];
    double t = t0+i*h;
    g_pos_v->SetPoint(i,x,y);
    g_vel_vx->SetPoint(i,t,xvel_v[i]);
    g_vel_vy->SetPoint(i,t,yvel_v[i]);
  }

  g_pos_v->SetLineColor(kBlue);
  g_pos_v->SetMarkerColor(kBlue);
  g_pos_v->SetMarkerStyle(3);
  g_pos_v->SetMarkerSize(1.5);
  g_vel_vx->SetLineColor(kBlue);
  g_vel_vx->SetMarkerColor(kBlue);
  g_vel_vx->SetMarkerStyle(3);
  g_vel_vx->SetMarkerSize(1.5);
  g_vel_vy->SetLineColor(kRed);
  g_vel_vy->SetMarkerColor(kRed);
  g_vel_vy->SetMarkerStyle(3);
  g_vel_vy->SetMarkerSize(1.5);

  TLegend *leg = new TLegend(0.450,0.6,1.0,0.9);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  leg->AddEntry(g_vel_vx,"v_x");
  leg->AddEntry(g_vel_vy,"v_y");

  TMultiGraph *m_pos_v = new TMultiGraph("m_pos_v","Position of Earth");
  m_pos_v->SetTitle("Position of Earth;x (AU);y(AU)");
  TMultiGraph *m_vel_v = new TMultiGraph("m_vel_v","Velocity of Earth");
  m_vel_v->SetTitle("Velocity of Earth;time (yr);v (AU/yr)");

  m_pos_v->Add(g_pos_v);
  m_vel_v->Add(g_vel_vx);
  m_vel_v->Add(g_vel_vy);

  TCanvas *c_pos_v = new TCanvas("c_pos_v","c_pos_v",800,720);
  TCanvas *c_vel_v = new TCanvas("c_vel_v","c_vel_v",800,720);
  
  c_pos_v->SetBorderMode(0);
  c_pos_v->cd();
  m_pos_v->Draw("AC*");
  c_pos_v->SaveAs("plots/earth_pos_verlet.png");
  c_pos_v->SaveAs("plots/earth_pos_verlet.pdf");
  c_pos_v->Close();
  
  c_vel_v->SetBorderMode(0);
  c_vel_v->cd();
  m_vel_v->Draw("AC*");
  leg->Draw("SAME");
  c_vel_v->SaveAs("plots/earth_vel_verlet.png");
  c_vel_v->SaveAs("plots/earth_vel_verlet.pdf");
  c_vel_v->Close();  

  //Include test of kinetic energy (and other tests?)

  //Now solve with RK4
  
}

void project3(){
  /*
    This is the main function which calls other functions in an attempt to 
    solve the problem posed in project 3.
  */

  parta();

}

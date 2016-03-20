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

  //Define the earth as a planet.
  planet earth = planet("earth",6e24,1);
  
  //Compute the expected velocity given a circular orbit.
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);

  //Define required quantities for Verlet method
  double t0 = 0; //Initial time (years)
  double tf = 1; //Final time (years)
  double nsteps = 10; //# of steps

  //Compute the x position and velocity of the earth over 1 year (1 revolution)
  //Assume the planet starts moving entirely in the x-direction.
  /*
    ERROR: NEED TO TAKE INTO ACCOUNT MULTIPLE COMPONENTS OF ACCELERATION
  */
  vector<thevec> xcomp = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v);
  vector<thevec> ycomp = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0);

  thevec xpos = xcomp.at(0);
  thevec xvel = xcomp.at(1);
  thevec ypos = ycomp.at(0);
  thevec yvel = ycomp.at(1);

  cout<<"xpos: "<<xpos.print()<<endl<<"ypos: "<<ypos.print()<<endl;
  cout<<"xvel: "<<xvel.print()<<endl<<"yvel: "<<yvel.print()<<endl;
  
  //Plot the position and velocity
  TGraph *g_pos = new TGraph(nsteps);
  TGraph *g_velx = new TGraph(nsteps);
  TGraph *g_vely = new TGraph(nsteps);

  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  for(int i=0;i<xpos.sz;i++){
    g_pos->SetPoint(i,xpos[i],ypos[i]);
    g_velx->SetPoint(i,t0+i*h,xvel[i]);
    g_vely->SetPoint(i,t0+i*h,yvel[i]);
  }

  g_pos->SetLineColor(kBlue);
  g_pos->SetMarkerColor(kBlue);
  g_pos->SetMarkerStyle(3);
  g_pos->SetMarkerSize(1.5);
  g_velx->SetLineColor(kBlue);
  g_velx->SetMarkerColor(kBlue);
  g_velx->SetMarkerStyle(3);
  g_velx->SetMarkerSize(1.5);
  g_vely->SetLineColor(kRed);
  g_vely->SetMarkerColor(kRed);
  g_vely->SetMarkerStyle(3);
  g_vely->SetMarkerSize(1.5);

  TLegend *leg = new TLegend(0.450,0.6,1.0,0.9);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);

  leg->AddEntry(g_velx,"v_x");
  leg->AddEntry(g_vely,"v_y");

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Position of Earth");
  m_pos->SetTitle("Position of Earth;x (AU);y(AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocity of Earth");
  m_vel->SetTitle("Velocity of Earth;time (yr);v (AU/yr)");

  m_pos->Add(g_pos);
  m_vel->Add(g_velx);
  m_vel->Add(g_vely);

  TCanvas *c_pos = new TCanvas("c_pos","c_pos",800,720);
  TCanvas *c_vel = new TCanvas("c_vel","c_vel",800,720);
  
  c_pos->cd();
  m_pos->Draw("AH*");
  c_pos->SaveAs("plots/earth_pos.png");
  c_pos->SaveAs("plots/earth_pos.pdf");
  c_pos->Close();
  
  c_vel->cd();
  m_vel->Draw("AH*");
  leg->Draw("SAME");
  c_vel->SaveAs("plots/earth_vel.png");
  c_vel->SaveAs("plots/earth_vel.pdf");
  c_vel->Close();  

}

void project3(){
  /*
    This is the main function which calls other functions in an attempt to 
    solve the problem posed in project 3.
  */

  parta();

}

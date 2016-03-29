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

void benchmarks(){
  /*
    Solve the earth-sun system and check for conservation of kinetic and 
    potential energy and angular momentum.
  */

  //Define the earth as a planet
  planet earth = planet("earth",6e24,1);
  
  //Compute the expected velocity given a circular orbit
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);
  
  //Define required quantities for Verlet
  double t0 = 0;
  double tf = 1;
  int nsteps = 100;

  //Solve the system using Verlet
  vector<thevec> xcomp = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos = xcomp.at(0);
  thevec xvel = xcomp.at(1);
  thevec ypos = ycomp.at(0);
  thevec yvel = ycomp.at(1);

  //Check for consistency
  double max_T, min_T, max_V, min_V, max_l, min_l;
  for(int i=0;i<xpos.sz;i++){
    double dist = sqrt(pow(xpos[i],2)+pow(ypos[i],2));
    double vel = sqrt(pow(xvel[i],2)+pow(yvel[i],2));
    double T = earth.kinetic(vel);
    double V = earth.potential(dist);
    double l = earth.ang_mom(vel,dist);

    if(i==0){
      max_T = T;
      min_T = max_T;
      max_V = V;
      min_V = max_V;
      max_l = l;
      min_l = max_l;
    }

    else{
      if(T < min_T)
	min_T = T;
      if(T > max_T)
	max_T = T;
      if(V < min_V)
	min_V = V;
      if(V > max_V)
	max_V = V;
      if(l < min_l)
	min_l = l;
      if(l > min_l)
	max_l = l;
    }
  }

  //Record results
  ofstream myfile;
  myfile.open("plots/benchmarks.txt");

  myfile<<"Max T: "<<max_T<<endl;
  myfile<<"Min T: "<<min_T<<endl;
  myfile<<"% diff: "<<100*(max_T - min_T)/max_T<<"%"<<endl<<endl;

  myfile<<"Max V: "<<max_V<<endl;
  myfile<<"Min V: "<<min_V<<endl;
  myfile<<"% diff: "<<100*(max_V - min_V)/max_V<<"%"<<endl<<endl;

  myfile<<"Max l: "<<max_l<<endl;
  myfile<<"Min l: "<<min_l<<endl;
  myfile<<"% diff: "<<100*(max_l - min_l)/max_l<<"%"<<endl;

  myfile.close();
}

void check_time_steps(){
  /*
    Solve the earth-sun system for various possible time steps and plot
    the stability.
  */

  gStyle->SetOptFit();

  //Define the earth as a planet
  planet earth = planet("earth",6e24,1);

  //Compute the expected velocity given a circular orbit
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);

  //Define required quantities for Verlet
  double t0 = 0;
  double tf = 1;

  //Define vector of possible numbers of steps
  vector<int> nsteps;
  for(int i=1;i<101;i++){
    if(i<10)
      nsteps.push_back(i);
    else if(i<30 && i%5==0)
      nsteps.push_back(i);
    else if(i%10==0)
      nsteps.push_back(i);
  }

  /*
    ISSUE: need to fix legend characteristics
  */

  //Define certain plotting requirements
  TLegend *leg = new TLegend(0.45,0.6,1.0,0.9);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.01);

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Position of Earth");
  m_pos->SetTitle("Position of Earth;x (AU);y (AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocity of Earth");
  m_vel->SetTitle("Velocity of Earth;v_x (AU/yr);v_y (AU/yr)");

  //Compute the solutions using verlet for various nsteps.
  for(unsigned int i=0;i<nsteps.size();i++){
    vector<thevec> xcomp = Verlet(t0,tf,nsteps.at(i),0,0,earth.acc,v,v,earth.dist_sun);
    vector<thevec> ycomp = Verlet(t0,tf,nsteps.at(i),-1.0*earth.dist_sun, -1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

    thevec xpos = xcomp.at(0);
    thevec xvel = xcomp.at(1);
    thevec ypos = ycomp.at(0);
    thevec yvel = ycomp.at(1);

    //Create the plots
    TGraph *g_pos = new TGraph();
    TGraph *g_vel = new TGraph();

    for(int j=0;j<xpos.sz;j++){
      g_pos->SetPoint(j,xpos[j],ypos[j]);
      g_vel->SetPoint(j,xvel[j],yvel[j]);
    }
    
    g_pos->SetLineColor(i);
    g_vel->SetLineColor(i);

    leg->AddEntry(g_pos,("nstep = "+to_string(nsteps.at(i))).c_str());
    m_pos->Add(g_pos);
    m_vel->Add(g_vel);
  }

  //Plot the graphs and save
  TCanvas *c_pos = new TCanvas("c_pos","c_pos",800,720);
  TCanvas *c_vel = new TCanvas("c_vel","c_vel",800,720);

  c_pos->SetBorderMode(0);
  c_pos->cd();
  m_pos->Draw("AC");
  leg->Draw("SAME");
  c_pos->SaveAs("plots/pos_nstep_check.png");
  c_pos->SaveAs("plots/pos_nstep_check.pdf");
  c_pos->Close();

  c_vel->SetBorderMode(0);
  c_vel->cd();
  m_vel->Draw("AC");
  leg->Draw("SAME");
  c_vel->SaveAs("plots/vel_nstep_check.png");
  c_vel->SaveAs("plots/vel_nstep_check.pdf");
  c_vel->Close();
}

void partb(){
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

  //Plot the position and velocity
  TGraph *g_pos_v = new TGraph();
  TGraph *g_vel_vx = new TGraph(nsteps+1);
  TGraph *g_vel_vy = new TGraph(nsteps+1);

  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  for(int i=0;i<nsteps+1;i++){
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

  partb();

  //Part c
  check_time_steps();
  benchmarks();
}

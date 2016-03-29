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

void partd_precise(){
  /*
    Determine by trial and error a more precise value for the escape velocity 
    than that found in partd() for planet earth.
  */

  //Define the planet
  planet earth = planet("earth",6e24,1);

  //Begin with the expected velocity of a stable circular orbit
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);

  //Vary velocities until the change is less than some delta
  double delta = 1e-10;

  double v0 = 0; double vf = v;
  do{

    v0 = vf;
    v = vf;

    double vstep = 0.1*v0;
    double dist = 0;

    double t0 = 0;
    double tf = 1;
    double nsteps = 100;

    do{

      //Solve the systep using Verlet
      vector<thevec> xcomp = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
      vector<thevec> ycomp = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

      thevec xpos = xcomp.at(0);
      thevec xvel = xcomp.at(1);
      thevec ypos = ycomp.at(0);
      thevec yvel = ycomp.at(1);

      dist = sqrt(pow((xpos[xpos.sz-1]-xpos[xpos.sz-2]),2)+pow((ypos[ypos.sz-1]-ypos[ypos.sz-2]),2));

      v = v+vstep;

    }while(dist < 1);

    vf = v-vstep;
  }while(abs(vf-v0)>delta);

  cout<<endl<<"Final Escape Velocity: "<<vf<<endl;
}

void partd(){
  /*
    Determine by trial and error what the initial velocity of a planet must
    be for it to escape from the sun if it starts at a distance of 1AU.
  */

  //First define a list of several possible masses of the planet.
  vector<double> masses;
  for(int i=1;i<11;i++)
    masses.push_back(i*pow(10,24));

  //Define a plot of escape velocity vs mass
  TGraph *g_escape = new TGraph();

  //Determine the escape velocity by trial and error.
  for(unsigned int i=0;i<masses.size();i++){
    double mass = masses.at(i);

    //Define the planet
    planet plan = planet("plan",mass,1);

    //Start with the expected velocity of a stable circular orbit
    double v = sqrt(4*pow(PI,2)/plan.dist_sun);

    //Check velocities every 10%
    double vstep = 0.01*v;

    //Can only check a finite number of years (say 10)
    double t0=0;
    double tf=10;
    double nsteps=100;

    double dist = 0;

    do{
      //Solve the system using Verlet
      vector<thevec> xcomp = Verlet(t0,tf,nsteps,0,0,plan.acc,v,v,plan.dist_sun);
      vector<thevec> ycomp = Verlet(t0,tf,nsteps,-1.0*plan.dist_sun,-1.0*plan.dist_sun,plan.acc,0,0,plan.dist_sun);
      
      thevec xpos = xcomp.at(0);
      thevec xvel = xcomp.at(1);
      thevec ypos = ycomp.at(0);
      thevec yvel = ycomp.at(0);
      
      //We consider a distance of 1AU from the expected final position to mean 
      //the planet escaped
      dist = sqrt(pow((xpos[xpos.sz-1]-xpos[xpos.sz-2]),2)+pow((ypos[ypos.sz-1]-ypos[ypos.sz-2]),2));
      
      v += vstep;

    }while(dist < 1);

    g_escape->SetPoint(i,mass,v-vstep);
  }

  //Plot results
  TMultiGraph *m_escape = new TMultiGraph("m_escape","Escape Velocities");
  m_escape->SetTitle("Escape Velocities;Planet Mass (kg);Escape Velocity(AU/yr)");

  m_escape->Add(g_escape);

  TCanvas *c_escape = new TCanvas("c_escape","c_escape",800,720);
  c_escape->SetBorderMode(0);
  c_escape->cd();
  m_escape->Draw("AC*");
  c_escape->SaveAs("plots/escape_vel.png");
  c_escape->SaveAs("plots/escape_vel.pdf");
  c_escape->Close();

  /*
    PROBLEM: escape velocity does not seem to depend on mass.  Should figure 
    out why and compute.
  */

}

void benchmarks(){
  /*
    Solve the earth-sun system and check for conservation of kinetic and 
    potential energy and angular momentum.
  */

  /*
    PROBLEM: Should make plot of differences
  */

  //Define the earth as a planet
  planet earth = planet("earth",6e24,1);
  
  //Compute the expected velocity given a circular orbit
  double v = sqrt(4*pow(PI,2)/earth.dist_sun);
  
  //Define required quantities for Verlet
  double t0 = 0;
  double tf = 1;
  int nsteps = 100;
  double h = (tf-t0)/nsteps;
  
  //Solve the system using Verlet
  vector<thevec> xcomp = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos = xcomp.at(0);
  thevec xvel = xcomp.at(1);
  thevec ypos = ycomp.at(0);
  thevec yvel = ycomp.at(1);

  //Define plots for energy and angular momentum
  TGraph *g_T = new TGraph();
  TGraph *g_V = new TGraph();
  TGraph *g_l = new TGraph();
  TGraph *g_tot = new TGraph();

  //Check for consistency
  double max_T, min_T, max_V, min_V, max_l, min_l, min_tot, max_tot;
  for(int i=0;i<xpos.sz;i++){
    double dist = sqrt(pow(xpos[i],2)+pow(ypos[i],2));
    double vel = sqrt(pow(xvel[i],2)+pow(yvel[i],2));
    double T = earth.kinetic(vel);
    double V = earth.potential(dist);
    double l = earth.ang_mom(vel,dist);
    double tot = T+V;

    //Plot energy and angular momentum
    g_T->SetPoint(i,h*i,T);
    g_V->SetPoint(i,h*i,V);
    g_tot->SetPoint(i,h*i,tot);
    g_l->SetPoint(i,h*i,l);

    if(i==0){
      max_T = T;
      min_T = max_T;
      max_V = V;
      min_V = max_V;
      max_l = l;
      min_l = max_l;
      min_tot = tot;
      max_tot = min_tot;
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
      if(tot < min_tot)
	min_tot = tot;
      if(tot > max_tot)
	max_tot = tot;
    }
  }

  //Record results
  ofstream myfile;
  myfile.open("benchmarks/benchmarks.txt");

  myfile<<"Max T: "<<max_T<<endl;
  myfile<<"Min T: "<<min_T<<endl;
  myfile<<"% diff: "<<100*(max_T - min_T)/max_T<<"%"<<endl<<endl;

  myfile<<"Max V: "<<max_V<<endl;
  myfile<<"Min V: "<<min_V<<endl;
  myfile<<"% diff: "<<100*(max_V - min_V)/max_V<<"%"<<endl<<endl;

  myfile<<"Max T+V: "<<max_tot<<endl;
  myfile<<"Min T+V: "<<min_tot<<endl;
  myfile<<"% diff: "<<100*(max_tot - min_tot)/max_tot<<"%"<<endl<<endl;

  myfile<<"Max l: "<<max_l<<endl;
  myfile<<"Min l: "<<min_l<<endl;
  myfile<<"% diff: "<<100*(max_l - min_l)/max_l<<"%"<<endl;

  myfile.close();

  //Plot results
  TMultiGraph *m_T = new TMultiGraph("m_T","Kinetic Energy");
  m_T->SetTitle("Kinetic Energy;time (AU);Energy");
  TMultiGraph *m_V = new TMultiGraph("m_V","Potential Energy");
  m_V->SetTitle("Potential Energy;time (AU);Energy");
  TMultiGraph *m_l = new TMultiGraph("m_l","Angular Momentum");
  m_l->SetTitle("Angular Momentum;time (AU);Angular Momentum");
  TMultiGraph *m_tot = new TMultiGraph("m_tot","Total Energy");
  m_tot->SetTitle("Total Energy;time (AU);Energy");

  m_T->Add(g_T);
  m_V->Add(g_V);
  m_l->Add(g_l);
  m_tot->Add(g_tot);

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);

  can->cd();
  m_T->Draw("AC*");
  can->SaveAs("benchmarks/kinetic.png");
  can->SaveAs("benchmarks/kinetic.pdf");
  m_V->Draw("AC*");
  can->SaveAs("benchmarks/potential.pdf");
  can->SaveAs("benchmarks/potential.png");
  m_tot->Draw("AC*");
  can->SaveAs("benchmarks/total_energy.png");
  can->SaveAs("benchmarks/total_energy.pdf");
  m_l->Draw("AC*");
  can->SaveAs("benchmarks/angular_mom.png");
  can->SaveAs("benchmarks/angular_mom.pdf");
  can->Close();

  /*
    PROBLEM: there seems to be an issue with the velocity calculation for Verlet
  */
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

  //Now solve with RK4
  
}

void project3(){
  /*
    This is the main function which calls other functions in an attempt to 
    solve the problem posed in project 3.
  */

  //partb();

  //Part c
  //check_time_steps();
  benchmarks();

  //partd();
  partd_precise();

}

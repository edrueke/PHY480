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

#include "solar_system.C"
#include "planets.C"
#include "odesolvers.C"
#include "classes.C"

using namespace std;

void parteandf(){
  /*
    Add Jupiter and the other planets and solve the solar system again
    using verlet and rk4 methods.  Initial velocities are obtained from 
    NASA site 
  */

  //km/sec->AU/yr
  double convert = 0.210805;

  planet earth = planet("earth",5.972e24,1,30*convert);
  planet jupiter = planet("jupiter",1898.13e24,5.2,47051*convert/3600);
  planet sun = planet("sun",1.989e30,0);

  solar_system milky_way = solar_system(earth);
  milky_way.Add(jupiter);
  milky_way.Add(sun);

  //First solve with just the earth and jupiter
  milky_way.nsteps = 50;
  milky_way.tf = 1;

  //milky_way.Solve_Verlet();
  //milky_way.Draw_Verlet("jupiter_only");

  //milky_way.Solve_RK4();
  //milky_way.Draw_RK4("jupiter_only");

  //Solve 3-body problem with com origin
  /*milky_way.Set_COM();
  milky_way.Solve_Verlet();
  milky_way.Draw_Verlet("jupiter_only_com");

  milky_way.Solve_RK4();
  milky_way.Draw_RK4("jupiter_only_com");*/

  //Now add other planets
  planet mercury = planet("mercury",3.285e23,0.387,47.4*convert);
  planet venus = planet("venus",4.867e24,0.723,126077*convert/3600);
  planet mars = planet("mars",6.39e23,1.524,86871*convert/3600);
  planet saturn = planet("saturn",5.683e26,9.539,34821*convert/3600);
  planet uranus = planet("uranus",8.681e25,19.2,24607*convert/3600);
  planet neptune = planet("neptune",1.024e26,30.1,19720*convert/3600);
  planet pluto = planet("pluto",1.309e22,40,17096*convert/3600);

  milky_way.Add(mercury);
  milky_way.Add(venus);
  milky_way.Add(mars);
  milky_way.Add(saturn);
  milky_way.Add(uranus);
  milky_way.Add(neptune);
  milky_way.Add(pluto);

  milky_way.Set_O();

  milky_way.Solve_Verlet();
  milky_way.Draw_Verlet("milky_way_full");

  milky_way.Solve_RK4();
  milky_way.Draw_RK4("milky_way_full");
  
  milky_way.Set_COM();

  milky_way.Solve_Verlet();
  milky_way.Draw_Verlet("milky_way_full_com");

  milky_way.Solve_RK4();
  milky_way.Draw_RK4("milky_way_full_com");
  
  
}

void partd_precise(){
  /*
    Determine by trial and error a more precise value for the escape velocity 
    than that found in partd() for planet earth using RK4.
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

      //Solve the systep using RK4
      vector<thevec> xcomp = RK4(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
      vector<thevec> ycomp = RK4(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

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
    be for it to escape from the sun if it starts at a distance of 1AU
    using the more precise RK4 method.
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
      //Solve the system using RK4
      vector<thevec> xcomp = RK4(t0,tf,nsteps,0,0,plan.acc,v,v,plan.dist_sun);
      vector<thevec> ycomp = RK4(t0,tf,nsteps,-1.0*plan.dist_sun,-1.0*plan.dist_sun,plan.acc,0,0,plan.dist_sun);
      
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

}

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
  double h = (tf-t0)/nsteps;
  
  //Solve the system using Verlet
  vector<thevec> xcomp_v = Verlet(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp_v = Verlet(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos_v = xcomp_v.at(0);
  thevec xvel_v = xcomp_v.at(1);
  thevec ypos_v = ycomp_v.at(0);
  thevec yvel_v = ycomp_v.at(1);

  //Solve the system using RK4
  vector<thevec> xcomp_r = RK4(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp_r = RK4(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos_r = xcomp_r.at(0);
  thevec xvel_r = xcomp_r.at(1);
  thevec ypos_r = ycomp_r.at(0);
  thevec yvel_r = ycomp_r.at(1);

  //Define plots for energy and angular momentum
  TGraph *g_T_v = new TGraph();
  TGraph *g_V_v = new TGraph();
  TGraph *g_l_v = new TGraph();
  TGraph *g_tot_v = new TGraph();

  TGraph *g_T_r = new TGraph();
  TGraph *g_V_r = new TGraph();
  TGraph *g_l_r = new TGraph();
  TGraph *g_tot_r = new TGraph();

  g_T_v->SetLineColor(kBlue);
  g_T_v->SetMarkerColor(kBlue);
  g_T_v->SetFillColor(0);
  g_T_v->SetFillStyle(0);

  g_V_v->SetLineColor(kBlue);
  g_V_v->SetMarkerColor(kBlue);
  g_V_v->SetFillColor(0);
  g_V_v->SetFillStyle(0);

  g_l_v->SetLineColor(kBlue);
  g_l_v->SetMarkerColor(kBlue);
  g_l_v->SetFillColor(0);
  g_l_v->SetFillStyle(0);

  g_tot_v->SetLineColor(kBlue);
  g_tot_v->SetMarkerColor(kBlue);
  g_tot_v->SetFillColor(0);
  g_tot_v->SetFillStyle(0);

  g_T_r->SetLineColor(kRed);
  g_T_r->SetMarkerColor(kRed);
  g_T_r->SetFillColor(0);
  g_T_r->SetFillStyle(0);

  g_V_r->SetLineColor(kRed);
  g_V_r->SetMarkerColor(kRed);
  g_V_r->SetFillColor(0);
  g_V_r->SetFillStyle(0);

  g_l_r->SetLineColor(kRed);
  g_l_r->SetMarkerColor(kRed);
  g_l_r->SetFillColor(0);
  g_l_r->SetFillStyle(0);

  g_tot_r->SetLineColor(kRed);
  g_tot_r->SetMarkerColor(kRed);
  g_tot_r->SetFillColor(0);
  g_tot_r->SetFillStyle(0);

  TLegend *leg = new TLegend(0.15,0.6,0.4,0.8);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);

  leg->AddEntry(g_T_v,"Verlet");
  leg->AddEntry(g_T_r,"RK4");

  //Check for consistency
  double max_T_v, min_T_v, max_V_v, min_V_v, max_l_v, min_l_v, min_tot_v, max_tot_v;
  double max_T_r, min_T_r, max_V_r, min_V_r, max_l_r, min_l_r, min_tot_r, max_tot_r;

  for(int i=0;i<xpos_v.sz;i++){
    double dist_v = sqrt(pow(xpos_v[i],2)+pow(ypos_v[i],2));
    double vel_v = sqrt(pow(xvel_v[i],2)+pow(yvel_v[i],2));
    double T_v = earth.kinetic(vel_v);
    double V_v = earth.potential(dist_v);
    double l_v = earth.ang_mom(vel_v,dist_v);
    double tot_v = T_v+V_v;

    double dist_r = sqrt(pow(xpos_r[i],2)+pow(ypos_r[i],2));
    double vel_r = sqrt(pow(xvel_r[i],2)+pow(yvel_r[i],2));
    double T_r = earth.kinetic(vel_r);
    double V_r = earth.potential(dist_r);
    double l_r = earth.ang_mom(vel_r,dist_r);
    double tot_r = T_r+V_r;

    //Plot energy and angular momentum
    g_T_v->SetPoint(i,h*i,T_v);
    g_V_v->SetPoint(i,h*i,V_v);
    g_tot_v->SetPoint(i,h*i,tot_v);
    g_l_v->SetPoint(i,h*i,l_v);

    g_T_r->SetPoint(i,h*i,T_r);
    g_V_r->SetPoint(i,h*i,V_r);
    g_tot_r->SetPoint(i,h*i,tot_r);
    g_l_r->SetPoint(i,h*i,l_r);

    if(i==0){

      max_T_v = T_v;
      min_T_v = max_T_v;
      max_V_v = V_v;
      min_V_v = max_V_v;
      max_l_v = l_v;
      min_l_v = max_l_v;
      min_tot_v = tot_v;
      max_tot_v = min_tot_v;

      max_T_r = T_r;
      min_T_r = max_T_r;
      max_V_r = V_r;
      min_V_r = max_V_r;
      max_l_r = l_r;
      min_l_r = max_l_r;
      min_tot_r = tot_r;
      max_tot_r = min_tot_r;

    }

    else{

      if(T_v < min_T_v)
	min_T_v = T_v;
      if(T_v > max_T_v)
	max_T_v = T_v;
      if(V_v < min_V_v)
	min_V_v = V_v;
      if(V_v > max_V_v)
	max_V_v = V_v;
      if(l_v < min_l_v)
	min_l_v = l_v;
      if(l_v > min_l_v)
	max_l_v = l_v;
      if(tot_v < min_tot_v)
	min_tot_v = tot_v;
      if(tot_v > max_tot_v)
	max_tot_v = tot_v;
      
      if(T_r < min_T_r)
	min_T_r = T_r;
      if(T_r > max_T_r)
	max_T_r = T_r;
      if(V_r < min_V_r)
	min_V_r = V_r;
      if(V_r > max_V_r)
	max_V_r = V_r;
      if(l_r < min_l_r)
	min_l_r = l_r;
      if(l_r > min_l_r)
	max_l_r = l_r;
      if(tot_r < min_tot_r)
	min_tot_r = tot_r;
      if(tot_r > max_tot_r)
	max_tot_r = tot_r;

    }
  }

  //Record results
  ofstream myfile;
  myfile.open("benchmarks/benchmarks.txt");

  myfile<<"Max T_v: "<<max_T_v<<endl;
  myfile<<"Min T_v: "<<min_T_v<<endl;
  myfile<<"% diff: "<<100*(max_T_v - min_T_v)/max_T_v<<"%"<<endl<<endl;

  myfile<<"Max V_v: "<<max_V_v<<endl;
  myfile<<"Min V_v: "<<min_V_v<<endl;
  myfile<<"% diff: "<<100*(max_V_v - min_V_v)/max_V_v<<"%"<<endl<<endl;

  myfile<<"Max T_v+V_v: "<<max_tot_v<<endl;
  myfile<<"Min T_v+V_v: "<<min_tot_v<<endl;
  myfile<<"% diff: "<<100*(max_tot_v - min_tot_v)/max_tot_v<<"%"<<endl<<endl;

  myfile<<"Max l_v: "<<max_l_v<<endl;
  myfile<<"Min l_v: "<<min_l_v<<endl;
  myfile<<"% diff: "<<100*(max_l_v - min_l_v)/max_l_v<<"%"<<endl;

  myfile<<"Max T_r: "<<max_T_r<<endl;
  myfile<<"Min T_r: "<<min_T_r<<endl;
  myfile<<"% diff: "<<100*(max_T_r - min_T_r)/max_T_r<<"%"<<endl<<endl;

  myfile<<"Max V_r: "<<max_V_r<<endl;
  myfile<<"Min V_r: "<<min_V_r<<endl;
  myfile<<"% diff: "<<100*(max_V_r - min_V_r)/max_V_r<<"%"<<endl<<endl;

  myfile<<"Max T_r+V_r: "<<max_tot_r<<endl;
  myfile<<"Min T_r+V_r: "<<min_tot_r<<endl;
  myfile<<"% diff: "<<100*(max_tot_r - min_tot_r)/max_tot_r<<"%"<<endl<<endl;

  myfile<<"Max l_r: "<<max_l_r<<endl;
  myfile<<"Min l_r: "<<min_l_r<<endl;
  myfile<<"% diff: "<<100*(max_l_r - min_l_r)/max_l_r<<"%"<<endl;

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

  m_T->Add(g_T_v);
  m_V->Add(g_V_v);
  m_l->Add(g_l_v);
  m_tot->Add(g_tot_v);
  m_T->Add(g_T_r);
  m_V->Add(g_V_r);
  m_l->Add(g_l_r);
  m_tot->Add(g_tot_r);

  TCanvas *can = new TCanvas("can","can",800,720);
  can->SetBorderMode(0);

  can->cd();
  m_T->Draw("AC*");
  leg->Draw("SAME");
  can->SaveAs("benchmarks/kinetic.png");
  can->SaveAs("benchmarks/kinetic.pdf");
  m_V->Draw("AC*");
  leg->Draw("SAME");  
  can->SaveAs("benchmarks/potential.pdf");
  can->SaveAs("benchmarks/potential.png");
  m_tot->Draw("AC*");
  leg->Draw("SAME");
  can->SaveAs("benchmarks/total_energy.png");
  can->SaveAs("benchmarks/total_energy.pdf");
  m_l->Draw("AC*");
  leg->Draw("SAME");
  can->SaveAs("benchmarks/angular_mom.png");
  can->SaveAs("benchmarks/angular_mom.pdf");
  can->Close();

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
  nsteps.push_back(2); nsteps.push_back(3); nsteps.push_back(5);
  nsteps.push_back(10); nsteps.push_back(20); nsteps.push_back(50);

  //Define certain plotting requirements
  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);

  TMultiGraph *m_pos_v = new TMultiGraph("m_pos_v","Position of Earth");
  m_pos_v->SetTitle("Position of Earth (Verlet);x (AU);y (AU)");
  TMultiGraph *m_vel_v = new TMultiGraph("m_vel_v","Velocity of Earth");
  m_vel_v->SetTitle("Velocity of Earth (Verlet);v_x (AU/yr);v_y (AU/yr)");

  TMultiGraph *m_pos_r = new TMultiGraph("m_pos_r","Position of Earth");
  m_pos_r->SetTitle("Position of Earth (RK4);x (AU);y (AU)");
  TMultiGraph *m_vel_r = new TMultiGraph("m_vel_r","Velocity of Earth");
  m_vel_r->SetTitle("Velocity of Earth (RK4);v_x (AU/yr);v_y (AU/yr)");

  //Compute the solutions for various nsteps.
  for(unsigned int i=0;i<nsteps.size();i++){

    //Using Verlet
    vector<thevec> xcomp_v = Verlet(t0,tf,nsteps.at(i),0,0,earth.acc,v,v,earth.dist_sun);
    vector<thevec> ycomp_v = Verlet(t0,tf,nsteps.at(i),-1.0*earth.dist_sun, -1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

    thevec xpos_v = xcomp_v.at(0);
    thevec xvel_v = xcomp_v.at(1);
    thevec ypos_v = ycomp_v.at(0);
    thevec yvel_v = ycomp_v.at(1);

    //Using RK4
    vector<thevec> xcomp_r = RK4(t0,tf,nsteps.at(i),0,0,earth.acc,v,v,earth.dist_sun);
    vector<thevec> ycomp_r = RK4(t0,tf,nsteps.at(i),-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

    thevec xpos_r = xcomp_r.at(0);
    thevec xvel_r = xcomp_r.at(1);
    thevec ypos_r = ycomp_r.at(0);
    thevec yvel_r = ycomp_r.at(1);

    //Create the plots
    TGraph *g_pos_v = new TGraph();
    TGraph *g_vel_v = new TGraph();
    TGraph *g_pos_r = new TGraph();
    TGraph *g_vel_r = new TGraph();

    for(int j=0;j<xpos_v.sz;j++){
      g_pos_v->SetPoint(j,xpos_v[j],ypos_v[j]);
      g_vel_v->SetPoint(j,xvel_v[j],yvel_v[j]);
      g_pos_r->SetPoint(j,xpos_r[j],ypos_r[j]);
      g_vel_r->SetPoint(j,xvel_r[j],yvel_r[j]);
    }
    
    g_pos_v->SetLineColor(i);
    g_vel_v->SetLineColor(i);
    g_pos_v->SetFillColor(0);
    g_pos_v->SetFillStyle(0);
    g_vel_v->SetFillColor(0);
    g_vel_v->SetFillStyle(0);

    g_pos_r->SetLineColor(i);
    g_pos_r->SetFillColor(0);
    g_pos_r->SetFillStyle(0);
    g_vel_r->SetLineColor(i);
    g_vel_r->SetFillColor(0);
    g_vel_r->SetFillStyle(0);

    leg->AddEntry(g_pos_v,("nstep = "+to_string(nsteps.at(i))).c_str());
    m_pos_v->Add(g_pos_v);
    m_vel_v->Add(g_vel_v);
    m_pos_r->Add(g_pos_r);
    m_vel_r->Add(g_vel_r);
  }

  //Plot the graphs and save
  TCanvas *c_pos_v = new TCanvas("c_pos_v","c_pos_v",800,720);
  TCanvas *c_vel_v = new TCanvas("c_vel_v","c_vel_v",800,720);
  TCanvas *c_pos_r = new TCanvas("c_pos_r","c_pos_r",800,720);
  TCanvas *c_vel_r = new TCanvas("c_vel_r","c_vel_r",800,720);

  c_pos_v->SetBorderMode(0);
  c_pos_v->cd();
  m_pos_v->Draw("AC");
  leg->Draw("SAME");
  c_pos_v->SaveAs("plots/pos_nstep_check_verlet.png");
  c_pos_v->SaveAs("plots/pos_nstep_check_verlet.pdf");
  c_pos_v->Close();

  c_vel_v->SetBorderMode(0);
  c_vel_v->cd();
  m_vel_v->Draw("AC");
  leg->Draw("SAME");
  c_vel_v->SaveAs("plots/vel_nstep_check_verlet.png");
  c_vel_v->SaveAs("plots/vel_nstep_check_verlet.pdf");
  c_vel_v->Close();

  c_pos_r->SetBorderMode(0);
  c_pos_r->cd();
  m_pos_r->Draw("AC");
  leg->Draw("SAME");
  c_pos_r->SaveAs("plots/pos_nstep_check_rk4.png");
  c_pos_r->SaveAs("plots/pos_nstep_check_rk4.pdf");
  c_pos_r->Close();

  c_vel_r->SetBorderMode(0);
  c_vel_r->cd();
  m_vel_r->Draw("AC");
  leg->Draw("SAME");  
  c_vel_r->SaveAs("plots/vel_nstep_check_rk4.png");
  c_vel_r->SaveAs("plots/vel_nstep_check_rk4.pdf");
  c_vel_r->Close();
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

  //Then using RK4
  vector<thevec> xcomp_r = RK4(t0,tf,nsteps,0,0,earth.acc,v,v,earth.dist_sun);
  vector<thevec> ycomp_r = RK4(t0,tf,nsteps,-1.0*earth.dist_sun,-1.0*earth.dist_sun,earth.acc,0,0,earth.dist_sun);

  thevec xpos_r = xcomp_r.at(0);
  thevec xvel_r = xcomp_r.at(1);
  thevec ypos_r = ycomp_r.at(0);
  thevec yvel_r = ycomp_r.at(1);

  //Plot the position and velocity
  TGraph *g_pos_v = new TGraph();
  TGraph *g_vel_v = new TGraph();
  TGraph *g_pos_r = new TGraph();
  TGraph *g_vel_r = new TGraph();

  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  for(int i=0;i<nsteps+1;i++){
    double t = t0+i*h;
    g_pos_v->SetPoint(i,xpos_v[i],ypos_v[i]);
    g_vel_v->SetPoint(i,xvel_v[i],yvel_v[i]);
    g_pos_r->SetPoint(i,xpos_r[i],ypos_r[i]);
    g_vel_r->SetPoint(i,xvel_r[i],yvel_r[i]);
  }

  g_pos_v->SetLineColor(kBlue);
  g_pos_v->SetFillColor(0);
  g_pos_v->SetFillStyle(0);
  g_pos_v->SetMarkerColor(kBlue);
  g_pos_v->SetMarkerStyle(3);
  g_pos_v->SetMarkerSize(1.5);

  g_pos_r->SetLineColor(kRed);
  g_pos_r->SetFillColor(0);
  g_pos_r->SetFillStyle(0);
  g_pos_r->SetMarkerColor(kRed);
  g_pos_r->SetMarkerStyle(3);
  g_pos_r->SetMarkerSize(1.5);

  g_vel_v->SetLineColor(kBlue);
  g_vel_v->SetFillColor(0);
  g_vel_v->SetFillStyle(0);
  g_vel_v->SetMarkerColor(kBlue);
  g_vel_v->SetMarkerStyle(3);
  g_vel_v->SetMarkerSize(1.5);

  g_vel_r->SetLineColor(kRed);
  g_vel_r->SetFillColor(0);
  g_vel_r->SetFillStyle(0);
  g_vel_r->SetMarkerColor(kRed);
  g_vel_r->SetMarkerStyle(3);
  g_vel_r->SetMarkerSize(1.5);

  TLegend *leg = new TLegend(0.350,0.4,0.75,0.7);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);

  leg->AddEntry(g_vel_v,"Verlet");
  leg->AddEntry(g_vel_r,"RK4");

  TMultiGraph *m_pos = new TMultiGraph("m_pos_v","Position of Earth");
  m_pos->SetTitle("Position of Earth;x (AU);y(AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel_v","Velocity of Earth");
  m_vel->SetTitle("Velocity of Earth;v_x (AU/yr);v_y (AU/yr)");

  m_pos->Add(g_pos_v);
  m_pos->Add(g_pos_r);
  m_vel->Add(g_vel_v);
  m_vel->Add(g_vel_r);

  TCanvas *c_pos = new TCanvas("c_pos","c_pos",800,720);
  TCanvas *c_vel = new TCanvas("c_vel","c_vel",800,720);
  
  c_pos->SetBorderMode(0);
  c_pos->cd();
  m_pos->Draw("AC*");
  leg->Draw("SAME");
  c_pos->SaveAs("plots/earth_pos.png");
  c_pos->SaveAs("plots/earth_pos.pdf");
  c_pos->Close();
  
  c_vel->SetBorderMode(0);
  c_vel->cd();
  m_vel->Draw("AC*");
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

  //partb();

  //Part c
  //check_time_steps();
  //benchmarks();

  //partd();
  //partd_precise();

  parteandf();
}

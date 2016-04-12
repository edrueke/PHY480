/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the definition file for the solar system which uses the planet class
and ode solvers to determine the positions and velocities of various planets 
added to it.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

#include "solar_system.h"
#include "planets.h"
#include "classes.h"
#include "odesolvers.h"

using namespace std;

//Solar System Functions

solar_system::solar_system(planet &p){
  /*
    Constructor which initializes the solar system from one planet.
  */

  planets.push_back(&p);
  nsteps = 100;
  originx = 0;
  originy = 0;
    
}

solar_system::solar_system(vector<planet*> ps){
  /*
    Constructor which initializes the solar system from a vector of planets.
  */

  planets = ps;
  nsteps = 100;
  originx = 0;
  originy = 0;

}

solar_system::solar_system(const solar_system &ss){
  /*
    Copy constructor.
  */

  /*
    PROBLEM: PROBABLY WANT COPY OF MAP RATHER THAN ACTUAL MAP
  */

  planets = ss.planets;
  nsteps = 100;
  originx = 0;
  originy = 0;
}

void solar_system::Add(planet &p){
  /*
    Adds a planet to the solar system to be included in the solution.
  */

  planets.push_back(&p);
}

void solar_system::Solve_Verlet(){
  /*
    Solves the solar system after all planets are added.
  */

  double h = 1.0*tf/nsteps;

  double fact = -4*pow(PI,2);
  double msun = 1.989e30; 

  planet *theplanets = new planet[planets.size()];

  for(unsigned int i=0;i<planets.size();i++){
    theplanets[i] = *planets.at(i);
  }
  
  for(unsigned int it=0;it<planets.size();it++){
    
    planet myplan = theplanets[it];

    //Solve initial positions and velocities using verlet
    thevec posx_v = thevec(2);
    thevec posy_v = thevec(2);
    thevec velx_v = thevec(1);
    thevec vely_v = thevec(1);

    posx_v.point[0] = 0;
    posx_v.point[1] = posx_v[0]+h*myplan.v0;

    posy_v.point[0] = -(myplan).dist_sun;
    posy_v.point[1] = posy_v[0];

    velx_v.point[0] = (myplan).v0;
    vely_v.point[0] = 0;

    theplanets[it].positionsx_v.point[0] = posx_v[0];
    theplanets[it].positionsx_v.point[1] = posx_v[1];
    theplanets[it].positionsy_v.point[0] = posy_v[0];
    theplanets[it].positionsy_v.point[1] = posy_v[1];
    theplanets[it].velocitiesx_v.point[0] = velx_v[0];
    theplanets[it].velocitiesy_v.point[0] = vely_v[0];

  }
  
  //Solve the first velocities
  for(unsigned int it = 0;it<planets.size();it++){
    
    planet myplan = theplanets[it];
    
    thevec velx_v = myplan.velocitiesx_v;
    thevec vely_v = myplan.velocitiesy_v;
    thevec posx_v = myplan.positionsx_v;
    thevec posy_v = myplan.positionsy_v;
    
    double r = sqrt(pow(posx_v[1]-originx,2)+pow(posy_v[1]-originy,2));
    double vx = velx_v[0];
    double vy = vely_v[0];
    for(unsigned int m = 0;m<planets.size();m++){
    
      planet plan = theplanets[m];
      double vxpi = (plan.positionsx_v)[1];
      double vypi = (plan.positionsy_v)[1];
      
      if(plan!=myplan){
	double rp = sqrt(pow(posx_v[1] - vxpi,2)+pow(posy_v[1]-vypi,2));
	vx += fact*(plan).mass*(posx_v[1] - vxpi)/(pow(rp,3)*msun);
	vy += fact*(plan).mass*(posy_v[1] - vypi)/(pow(rp,3)*msun);
      }
      
    }
    
    theplanets[it].velocitiesx_v.point[1] = vx;
    theplanets[it].velocitiesy_v.point[1] = vy;
        
  }
  
  //Solve the rest of the system
  
  for(int i=1;i<nsteps;i++){
    
    cout<<"i = "<<i<<endl;

    //Compute positions first
    for(unsigned int it = 0;it<planets.size();it++){
      
      planet myplan = theplanets[it];
      
      double vx_prev = myplan.velocitiesx_v[i-1];
      double vy_prev = myplan.velocitiesy_v[i-1];
      double x_prev = myplan.positionsx_v[i-1];
      double y_prev = myplan.positionsy_v[i-1];
      
      double r_prev = sqrt(pow(x_prev-originx,2)+pow(y_prev-originy,2));
      
      double x = x_prev+h*vx_prev;
      double y = y_prev+h*vy_prev;
      
      for(unsigned int m = 0;m<planets.size();m++){
	
	planet plan = theplanets[m];
	
	if(myplan != plan){
	  
	  double xp_prev = plan.positionsx_v[i-1];
	  double yp_prev = plan.positionsy_v[i-1];
	  
	  double rp_prev = sqrt(pow(x_prev-xp_prev,2)+pow(y_prev-yp_prev,2));
	  
	  x += (pow(h,2)/2)*fact*plan.mass*(x_prev-xp_prev)/(msun*pow(rp_prev,3));
	  y += (pow(h,2)/2)*fact*plan.mass*(y_prev-yp_prev)/(msun*pow(rp_prev,3));
	  
	}
      }
      
      theplanets[it].positionsx_v.point[i] = x;
      theplanets[it].positionsy_v.point[i] = y;
      
    }
    
    //Then compute velocities
    for(unsigned int it = 0;it<planets.size();it++){
      
      planet myplan = theplanets[it];
      
      double vx_prev = myplan.velocitiesx_v[i-1];
      double vy_prev = myplan.velocitiesy_v[i-1];
      double x_prev = myplan.positionsx_v[i-1];
      double y_prev = myplan.positionsy_v[i-1];
      
      double r_prev = sqrt(pow(x_prev-originx,2)+pow(y_prev-originy,2));
      
      double x = myplan.positionsx_v[i];
      double y = myplan.positionsy_v[i];
      double newr = sqrt(pow(x-originx,2)+pow(y-originy,2));
    
      double vx = vx_prev;
      double vy = vy_prev;
      
      for(unsigned int m = 0;m<planets.size();m++){
	
	planet plan = theplanets[m];
	
	if(myplan != plan){
	  
	  double xp_prev = plan.positionsx_v[i-1];
	  double yp_prev = plan.positionsy_v[i-1];
	  
	  double xp_cur = plan.positionsx_v[i];
	  double yp_cur = plan.positionsy_v[i];
	  
	  double rp_cur = sqrt(pow(x-xp_cur,2)+pow(y-yp_cur,2));
	  double rp_prev = sqrt(pow(x_prev-xp_prev,2)+pow(y_prev-yp_prev,2));
	  
	  vx += ((h/2)*fact*plan.mass/msun)*((x-xp_cur)/pow(rp_cur,3)+(x_prev-xp_prev)/pow(rp_prev,3));
	  vy += ((h/2)*fact*plan.mass/msun)*((y-yp_cur)/pow(rp_cur,3)+(y_prev-yp_prev)/pow(rp_prev,3));
	  
	}
      }
      
      theplanets[it].velocitiesx_v.point[i] = vx;
      theplanets[it].velocitiesy_v.point[i] = vy;
    
    }
  }
  
  vector<planet*> temp;
  for(unsigned int i=0;i<planets.size();i++){
    temp.push_back(&theplanets[i]);
  }
  
  planets = temp;
  
}

void solar_system::Solve_RK4(){
  /*
    Solve the solar system using the RK4 algorithm.
  */

  double h = 1.0*tf/nsteps;

  double fact = -4.0*pow(PI,2);
  double msun = 1.989e30;

  planet *theplanets = new planet[planets.size()];

  for(unsigned int i=0;i<planets.size();i++)
    theplanets[i] = *planets.at(i);

  //Set initial positions and velocities
  for(unsigned int it=0;it<planets.size();it++){

    planet myplan = theplanets[it];

    theplanets[it].positionsx_r.point[0] = 0;
    theplanets[it].positionsy_r.point[0] = -myplan.dist_sun;

    theplanets[it].velocitiesx_r.point[0] = myplan.v0;
    theplanets[it].velocitiesy_r.point[0] = 0;

  }

  //Solve using RK4 algorithm
  for(int i = 1;i<nsteps;i++){
    
    cout<<"i = "<<i<<endl;

    for(unsigned int it=0;it<planets.size();it++){
      
      planet myplan = theplanets[it];
      
      thevec posx = myplan.positionsx_r;
      thevec posy = myplan.positionsy_r;
      thevec velx = myplan.velocitiesx_r;
      thevec vely = myplan.velocitiesy_r;
      
      double prev_r = sqrt(pow(posx[i-1]-originx,2)+pow(posy[i-1]-originy,2));

      double k1vx = 0;
      double k1vy = 0;
      double k1x = h*velx[i-1];
      double k1y = h*vely[i-1];

      for(unsigned int n=0;n<planets.size();n++){

	planet plan = theplanets[n];

	if(myplan!=plan){

	  double prev_rp = sqrt(pow(posx[i-1]-plan.positionsx_r[i-1],2)+pow(posy[i-1]-plan.positionsy_r[i-1],2));
	  k1vx+=h*fact*plan.mass*(posx[i-1]-plan.positionsx_r[i-1])/pow(prev_rp,3)/msun;
	  k1vy+=h*fact*plan.mass*(posy[i-1]-plan.positionsy_r[i-1])/pow(prev_rp,3)/msun;
	}
      }
      
      double k2vx = 0;
      double k2vy = 0;
      double k2x = h*(velx[i-1]+k1vx/2);
      double k2y = h*(vely[i-1]+k1vy/2);

      for(unsigned int n=0;n<planets.size();n++){

	planet plan = theplanets[n];

	if(myplan!=plan){

	  double prev_rp = sqrt(pow(posx[i-1]-plan.positionsx_r[i-1],2)+pow(posy[i-1]-plan.positionsy_r[i-1],2));
	  k2vx+=h*fact*plan.mass*(posx[i-1]+k1x/2-plan.positionsx_r[i-1])/pow(prev_rp,3)/msun;
	  k2vy+=h*fact*plan.mass*(posy[i-1]+k1y/2-plan.positionsy_r[i-1])/pow(prev_rp,3)/msun;
	}
      }
      
      double k3vx = 0;
      double k3vy = 0;
      double k3x = h*(velx[i-1]+k2vx/2);
      double k3y = h*(vely[i-1]+k2vy/2);

      for(unsigned int n=0;n<planets.size();n++){

	planet plan = theplanets[n];

	if(myplan!=plan){

	  double prev_rp = sqrt(pow(posx[i-1]-plan.positionsx_r[i-1],2)+pow(posy[i-1]-plan.positionsy_r[i-1],2));
	  k3vx+=h*fact*plan.mass*(posx[i-1]+k2x/2-plan.positionsx_r[i-1])/pow(prev_rp,3)/msun;
	  k3vy+=h*fact*plan.mass*(posy[i-1]+k2y/2-plan.positionsy_r[i-1])/pow(prev_rp,3)/msun;
	}
      }
      
      double k4vx = 0;
      double k4vy = 0;
      double k4x = h*(velx[i-1]+k3vx);
      double k4y = h*(vely[i-1]+k3vy);

      for(unsigned int n=0;n<planets.size();n++){

	planet plan = theplanets[n];

	if(myplan!=plan){

	  double prev_rp = sqrt(pow(posx[i-1]-plan.positionsx_r[i-1],2)+pow(posy[i-1]-plan.positionsy_r[i-1],2));
	  k4vx+=h*fact*plan.mass*(posx[i-1]+k3x-plan.positionsx_r[i-1])/pow(prev_rp,3)/msun;
	  k4vy+=h*fact*plan.mass*(posy[i-1]+k3y-plan.positionsy_r[i-1])/pow(prev_rp,3)/msun;
	}
      }

      theplanets[it].positionsx_r.point[i] = posx[i-1]+(1.0/6)*(k1x+2.0*k2x+2.0*k3x+k4x);
      theplanets[it].positionsy_r.point[i] = posy[i-1]+(1.0/6)*(k1y+2.0*k2y+2.0*k3y+k4y);
      theplanets[it].velocitiesx_r.point[i] = velx[i-1]+(1.0/6)*(k1vx+2.0*k2vx+2.0*k3vx+k4vx);
      theplanets[it].velocitiesy_r.point[i] = vely[i-1]+(1.0/6)*(k1vy+2.0*k2vy+2.0*k3vy+k4vy);

    }

  }
  
  vector<planet*> temp;
  for(unsigned int i=0;i<planets.size();i++)
    temp.push_back(&theplanets[i]);
  
  planets = temp;
  
}

void solar_system::Draw_Verlet(string name){
  /*
    Draw the paths of the planets and save it in the plots/ directory as
    name.
  */

  gStyle->SetOptFit();

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Positions");
  m_pos->SetTitle("Positions (Verlet);x (AU);y (AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocities");
  m_vel->SetTitle("Velocities (Verlet);v_x (AU/yr);v_y (AU/yr)");

  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);

  for(int i=0;i<planets.size();i++){
    
    planet p = *(planets.at(i));

    TGraph *g_pos = new TGraph();
    TGraph *g_vel = new TGraph();

    for(int j=0;j<nsteps;j++){
      g_pos->SetPoint(j,p.positionsx_v[j],p.positionsy_v[j]);
    }

    for(int j=0;j<nsteps;j++){
      g_vel->SetPoint(j,p.velocitiesx_v[j],p.velocitiesy_v[j]);
    }

    if(i<9){
      g_pos->SetLineColor(i+1);
      g_pos->SetMarkerColor(i+1);
      g_pos->SetMarkerSize(0);
      g_pos->SetFillColor(0);
      g_pos->SetFillStyle(0);
      
      g_vel->SetLineColor(i+1);
      g_vel->SetMarkerColor(i+1);
      g_vel->SetMarkerSize(0);
      g_vel->SetFillColor(0);
      g_vel->SetFillStyle(0);
    }
    else{
      g_pos->SetLineColor(i+30);
      g_pos->SetMarkerColor(i+30);
      g_pos->SetMarkerSize(0);
      g_pos->SetFillColor(0);
      g_pos->SetFillStyle(0);
      
      g_vel->SetLineColor(i+30);
      g_vel->SetMarkerColor(i+30);
      g_vel->SetMarkerSize(0);
      g_vel->SetFillColor(0);
      g_vel->SetFillStyle(0);
    }

    leg->AddEntry(g_pos,(p.name).c_str());

    m_pos->Add(g_pos);
    m_vel->Add(g_vel);

  }

  TCanvas *c_pos = new TCanvas("c_pos","c_pos",800,720);
  TCanvas *c_vel = new TCanvas("c_vel","c_vel",800,720);

  leg->SetFillStyle(0);

  c_pos->SetBorderMode(0);
  c_pos->cd();
  m_pos->Draw("AC*");
  leg->Draw("SAME");
  c_pos->SaveAs(("plots/"+name+"_verlet_pos.png").c_str());
  c_pos->SaveAs(("plots/"+name+"_verlet_pos.pdf").c_str());
  c_pos->Close();

  c_vel->SetBorderMode(0);
  c_vel->cd();
  m_vel->Draw("AC*");
  leg->Draw("SAME");
  c_vel->SaveAs(("plots/"+name+"_verlet_vel.png").c_str());
  c_vel->SaveAs(("plots/"+name+"_verlet_vel.pdf").c_str());
  c_vel->Close();

}

void solar_system::Draw_RK4(string name){
  /*
    Draw the paths of the planets and save it in the plots/ directory as
    name.
  */

  gStyle->SetOptFit();

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Positions");
  m_pos->SetTitle("Positions (RK4);x (AU);y (AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocities");
  m_vel->SetTitle("Velocities (RK4);v_x (AU/yr);v_y (AU/yr)");

  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);

  for(int i=0;i<planets.size();i++){
    
    planet p = *(planets.at(i));

    TGraph *g_pos = new TGraph();
    TGraph *g_vel = new TGraph();

    for(int j=0;j<nsteps;j++){
      g_pos->SetPoint(j,p.positionsx_r[j],p.positionsy_r[j]);
    }

    for(int j=0;j<nsteps;j++){
      g_vel->SetPoint(j,p.velocitiesx_r[j],p.velocitiesy_r[j]);
    }

    if(i<9){
      g_pos->SetLineColor(i+1);
      g_pos->SetMarkerColor(i+1);
      g_pos->SetMarkerSize(0);
      g_pos->SetFillColor(0);
      g_pos->SetFillStyle(0);
      
      g_vel->SetLineColor(i+1);
      g_vel->SetMarkerColor(i+1);
      g_vel->SetMarkerSize(0);
      g_vel->SetFillColor(0);
      g_vel->SetFillStyle(0);
    }
    else{
      g_pos->SetLineColor(i+30);
      g_pos->SetMarkerColor(i+30);
      g_pos->SetMarkerSize(0);
      g_pos->SetFillColor(0);
      g_pos->SetFillStyle(0);
      
      g_vel->SetLineColor(i+30);
      g_vel->SetMarkerColor(i+30);
      g_vel->SetMarkerSize(0);
      g_vel->SetFillColor(0);
      g_vel->SetFillStyle(0);
    }

    leg->AddEntry(g_pos,(p.name).c_str());

    m_pos->Add(g_pos);
    m_vel->Add(g_vel);

  }

  TCanvas *c_pos = new TCanvas("c_pos","c_pos",800,720);
  TCanvas *c_vel = new TCanvas("c_vel","c_vel",800,720);

  c_pos->SetBorderMode(0);
  c_pos->cd();
  m_pos->Draw("AC*");
  leg->Draw("SAME");
  c_pos->SaveAs(("plots/"+name+"_rk4_pos.png").c_str());
  c_pos->SaveAs(("plots/"+name+"_rk4_pos.pdf").c_str());
  c_pos->Close();

  c_vel->SetBorderMode(0);
  c_vel->cd();
  m_vel->Draw("AC*");
  leg->Draw("SAME");
  c_vel->SaveAs(("plots/"+name+"_rk4_vel.png").c_str());
  c_vel->SaveAs(("plots/"+name+"_rk4_vel.pdf").c_str());
  c_vel->Close();

}

void solar_system::Set_COM(){
  /*
    Set the origin of the solar system to the center of mass of the system.
  */

  double num = 0;
  double den = 0;

  //com = sum(x_{i}m_{i})/su,(m_{i})

  for(unsigned int i=0;i<planets.size();i++){
    
    planet plan = *planets.at(i);

    num+=plan.dist_sun*plan.mass;
    den+=plan.mass;

  }

  originx = 0;
  originy = -1.0*num/den;

}

void solar_system::Set_O(){
  /*
    Re-zero the origin of the solar system
  */

  originx=0;
  originy=0;
}

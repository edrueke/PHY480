/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the definition file for the solar system which uses the planet class
and ode solvers to determine the positions and velocities of various planets 
added to it.
*/

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
}

solar_system::solar_system(vector<planet*> ps){
  /*
    Constructor which initializes the solar system from a vector of planets.
  */

  planets = ps;

}

solar_system::solar_system(const solar_system &ss){
  /*
    Copy constructor.
  */

  /*
    PROBLEM: PROBABLY WANT COPY OF MAP RATHER THAN ACTUAL MAP
  */

  planets = ss.planets;
}

void solar_system::Add(planet &p){
  /*
    Adds a planet to the solar system to be included in the solution.
  */

  planets.push_back(&p);
}

void solar_system::Solve_Verlet(int nsteps, double tf){
  /*
    Solves the solar system after all planets are added.
  */

  double h = 1.0*tf/nsteps;

  double fact = -1.0*4*pow(PI,2);
  double msun = 1.989e30; 

  for(unsigned int it=0;it<planets.size();it++){

    planet *myplan = planets.at(it);

    //Solve positions and velocities using verlet
    thevec posx_v = thevec(2);
    thevec posy_v = thevec(2);
    thevec velx_v = thevec(1);
    thevec vely_v = thevec(1);

    //Initial positions and velocities
    posx_v.point[0] = 0;
    posx_v.point[1] = posx_v[0]+h*(*myplan).v0;

    posy_v.point[0] = -1.0*(*myplan).dist_sun;
    posy_v.point[1] = posy_v[0]+h*(*myplan).v0;

    velx_v.point[0] = (*myplan).v0;
    vely_v.point[0] = 0;

    (*planets.at(it)).positionsx_v = posx_v;
    (*planets.at(it)).positionsy_v = posy_v;
    (*planets.at(it)).velocitiesx_v = velx_v;
    (*planets.at(it)).velocitiesy_v = vely_v;
  }

  //Solve the first velocities
  for(unsigned int it = 0;it<planets.size();it++){

    planet *myplan = planets.at(it);

    thevec velx_v = (*myplan).velocitiesx_v;
    thevec vely_v = (*myplan).velocitiesy_v;
    thevec posx_v = (*myplan).positionsx_v;
    thevec posy_v = (*myplan).positionsy_v;

    //Solve the first velocity
    double r = sqrt(pow(posx_v[1],2)+pow(posy_v[1],2));
    double vx = velx_v[0]+(h/2)*(fact/pow(r,3));
    double vy = vely_v[0]+(h/2)*(fact/pow(r,3));
    for(unsigned int m = 0;m<planets.size();m++){
      
      planet *plan = planets.at(m);
      double vxpi = ((*plan).positionsx_v)[1];
      double vypi = ((*plan).positionsy_v)[1];

      if(plan!=myplan){
	double rp = sqrt(pow(posx_v[1] - vxpi,2)+pow(posy_v[1]-vypi,2));
	vx += fact*(*plan).mass*(posx_v[1] - vxpi)/(pow(rp,3)*msun);
	vy += fact*(*plan).mass*(posy_v[1] - vypi)/(pow(rp,3)*msun);
      }
    
    }

    ((*(planets.at(it))).velocitiesx_v).Add(vx);
    ((*planets.at(it)).velocitiesy_v).Add(vy);

  }
  
  //Solve the first velocity
  
  for(int i=2;i<nsteps+1;i++){
    
    //Compute positions first
    for(unsigned int it = 0;it<planets.size();it++){
      
      planet *myplan = planets.at(it);
      
      double vx_prev = ((*myplan).velocitiesx_v)[i-1];
      //double vy_prev = ((*myplan).velocitiesy_v)[i-1];
      //double x_prev = ((*myplan).positionsx_v)[i-1];
      //double y_prev = ((*myplan).positionsy_v)[i-1];
      
      /*double r_prev = sqrt(pow(x_prev,2)+pow(y_prev,2));
      
      double x = x_prev+h*vx_prev+(pow(h,2)/2)*fact*x_prev/pow(r_prev,3);
      double y = y_prev+h*vy_prev+(pow(h,2)/2)*fact*y_prev/pow(r_prev,3);
      
      for(unsigned int m = 0;m<planets.size();m++){
	
	planet *plan = planets.at(m);
	
	if(myplan != plan){
	  
	  double xp_prev = ((*plan).positionsx_v)[i-1];
	  double yp_prev = ((*plan).positionsy_v)[i-1];
	  
	  double rp_prev = sqrt(pow(x_prev - xp_prev,2)+pow(y_prev - yp_prev,2));
	  
	  x += (pow(h,2)/2)*fact*(*plan).mass*(x_prev - xp_prev)/(msun*pow(rp_prev,3));
	  y += (pow(h,2)/2)*fact*(*plan).mass*(y_prev - yp_prev)/(msun*pow(rp_prev,3));
	  
	}
      }
      
      ((*planets.at(it)).positionsx_v).Add(x);
      ((*planets.at(it)).positionsy_v).Add(y);*/
      
    }
    
    //Then compute velocities
    /*for(unsigned int it = 0;it<planets.size();it++){
      
      planet *myplan = planets.at(it);
      
      double vx_prev = ((*myplan).velocitiesx_v)[i-1];
      double vy_prev = ((*myplan).velocitiesy_v)[i-1];
      double x_prev = ((*myplan).positionsx_v)[i-1];
      double y_prev = ((*myplan).positionsy_v)[i-1];
      
      double r_prev = sqrt(pow(x_prev,2)+pow(y_prev,2));
      
      double x = ((*myplan).positionsx_v)[i];
      double y = ((*myplan).positionsy_v)[i];
      double newr = sqrt(pow(x,2)+pow(y,2));

      double vx = vx_prev + (h/2)*fact*(x/pow(newr,3)+x_prev/pow(r_prev,3));
      double vy = vy_prev + (h/2)*fact*(y/pow(newr,3)+y_prev/pow(r_prev,3));

      for(unsigned int m = 0;m<planets.size();m++){

	planet *plan = planets.at(m);

	if(myplan != plan){

	  double xp_prev = ((*plan).positionsx_v)[i-1];
	  double yp_prev = ((*plan).positionsy_v)[i-1];

	  double xp_cur = ((*plan).positionsx_v)[i];
	  double yp_cur = ((*plan).positionsy_v)[i];

	  double rp_cur = sqrt(pow(x - xp_cur,2) + pow(y - yp_cur,2));
	  double rp_prev = sqrt(pow(x_prev - xp_prev,2)+pow(y_prev - yp_prev,2));
	  
	  vx += ((h/2)*fact*(*plan).mass/msun)*((x - xp_cur)/pow(rp_cur,3) + (x_prev - xp_prev)/pow(rp_prev,3));
	  vy += ((h/2)*fact*(*plan).mass/msun)*((y - yp_cur)/pow(rp_cur,3) + (y_prev - yp_prev)/pow(rp_prev,3));

	}
      }
      
      ((*planets.at(it)).velocitiesx_v).Add(vx);
      ((*planets.at(it)).velocitiesy_v).Add(vy);
    
      }*/
  }
}

void solar_system::Draw_Verlet(string name){
  /*
    Draw the paths of the planets and save it in the plots/ directory as
    name.
  */

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Positions");
  m_pos->SetTitle("Positions (Verlet);x (AU);y (AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocities");
  m_vel->SetTitle("Velocities (Verlet);v_x (AU/yr);v_y (AU/yr)");

  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);

  for(int i=0;i<planets.size();i++){
    
    planet p = *(planets.at(i));

    TGraph *g_pos = new TGraph();
    TGraph *g_vel = new TGraph();

    for(int j=0;j<(p.positionsx_v).sz;j++){
      g_pos->SetPoint(j,p.positionsx_v[j],p.positionsy_v[j]);
      g_vel->SetPoint(j,p.velocitiesx_v[j],p.velocitiesy_v[j]);
    }

    g_pos->SetLineColor(i+1);
    g_pos->SetMarkerColor(i+1);
    g_pos->SetFillColor(0);
    g_pos->SetFillStyle(0);

    g_vel->SetLineColor(i+1);
    g_vel->SetMarkerColor(i+1);
    g_vel->SetFillColor(0);
    g_vel->SetFillStyle(0);

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

  TMultiGraph *m_pos = new TMultiGraph("m_pos","Positions");
  m_pos->SetTitle("Positions (RK4);x (AU);y (AU)");
  TMultiGraph *m_vel = new TMultiGraph("m_vel","Velocities");
  m_vel->SetTitle("Velocities (RK4);v_x (AU/yr);v_y (AU/yr)");

  TLegend *leg = new TLegend(0.2,0.15,0.4,0.4);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);

  for(int i=0;i<planets.size();i++){
    
    planet p = *(planets.at(i));

    TGraph *g_pos = new TGraph();
    TGraph *g_vel = new TGraph();

    for(int j=0;j<(p.positionsx_r).sz;j++){
      g_pos->SetPoint(j,p.positionsx_r[j],p.positionsy_r[j]);
      g_vel->SetPoint(j,p.velocitiesx_r[j],p.velocitiesy_r[j]);
    }

    g_pos->SetLineColor(i+1);
    g_pos->SetMarkerColor(i+1);
    g_pos->SetFillColor(0);
    g_pos->SetFillStyle(0);

    g_vel->SetLineColor(i+1);
    g_vel->SetMarkerColor(i+1);
    g_vel->SetFillColor(0);
    g_vel->SetFillStyle(0);

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

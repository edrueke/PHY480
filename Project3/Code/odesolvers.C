/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the definition file for the ODE solvers required when solving the 
solar system problems.  They work with the thevec and themat classes from 
classes.C and classes.h developed in Projects 1 and 2.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

#include "classes.h"
#include "odesolvers.h"

using namespace std;

vector<thevec> Verlet(double t0, double tf, int nsteps, double x0, double xf, double a, double v0, double vf,double r){
  /*
    Uses the Verlet Algorithm discussed in Chapter 8.3.2 of the lecture notes
    to solve the differential equation from Newton's second law for the 
    position and velocity.

    NOTE: This method does not take into account gravitational effects of one 
    planet on another, only the effects of the sun on individual planets.
    
    NOTE: May need to change so that we are including planets as a class.

    Args:
    -t0: a double for the initial time (years)
    -tf: a double for the final time (years)
    -nsteps: an integer giving the number of steps we want to use in the 
    discretization
    -x0: a double for the initial position
    -xf: a double for the final position
    -a: a double for the acceleration 
    -v0: double for initial velocity
    -vf: double for final velocity
    -r: constant distance between planet and sun (AU)
  */

  //Compute the step size.
  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  //Define a vector of the inputs to the function.
  thevec times = thevec(nsteps+1);
  times.point[0] = t0;
  for(int i=1;i<nsteps+1;i++){
    times.point[i] = t0+i*h;
  }

  //Define a vector for the positions.
  thevec pos = thevec(nsteps+1);
  pos.point[0]=x0;
  //pos.point[nsteps]=xf;

  //Note that x_1 = x_0+h*v_0+O(h^2)
  pos.point[1]=pos[0]+h*v0;

  //Compute the positions using the algorithm: x_{i+1} = 2x_{i}-x_{i-1}+h^{2}a
  for(int i=2;i<nsteps+1;i++){
    pos.point[i] = 2*pos[i-1]-pos[i-2]-a*pow(h,2)*pos[i-1]/r;
  }

  //Define a vector for the velocity.
  thevec vel = thevec(nsteps+1);
  vel.point[0]=v0;
  vel.point[nsteps]=vf;

  //Compute the velocities using the algorithm: v_{i} = v_{i}+(h/2)*(v'_{i+1}+v'_{i})
  for(int i=1;i<nsteps;i++){
    vel.point[i] = vel[i-1]+(h/2)*(-1.0*a*pos[i]/pow(r,3)-1.0*a*pos[i-1]/pow(r,3));//(pos[i+1]-pos[i-1])/(2*h);
  }

  //Define the vector to return
  vector<thevec> to_ret;
  to_ret.push_back(pos); to_ret.push_back(vel);

  return to_ret;
}

vector<thevec> RK4(double t0, double tf, int nsteps, double x0, double xf, double a, double v0, double vf, double r){
  /*
    Uses the 4th-order Runge-Kutta Algorithm, discussed in Chapter 8.4 of 
    the lecture notes to solve the differential equation from Newton's 
    second law for the position and velocity.

    NOTE: This method does not take into account gravitational effects of one 
    planet on another, only the effects of the sun on individual planets.
    
    NOTE: May need to change so that we are including planets as a class.

    Args:
    -t0: a double for the initial time
    -tf: a double for the final time
    -nsteps: an integer giving the number of steps we want to use in the 
    discretization
    -x0: a double for the initial position
    -xf: a double for the final position
    -a: a double for the acceleration 
    -v0: double for initial velocity
    -vf: double for final velocity
    -r: total radius in the case of a circular orbit (AU)
  */

  //Compute the step size.
  double h = (1.0*tf-1.0*t0)/(1.0*nsteps);

  //Define a vector of the inputs to the function.
  thevec times = thevec(nsteps+1);
  times.point[0] = t0;
  for(int i=1;i<nsteps+1;i++){
    times.point[i] = t0+i*h;
  }

  //Define a vector for the positions.
  thevec pos = thevec(nsteps+1);
  pos.point[0]=x0;
  pos.point[nsteps]=xf;

  //Define a vector for the velocities.
  thevec vel = thevec(nsteps+1);
  vel.point[0]=v0;
  vel.point[nsteps]=vf;

  //Compute the positions and velocities using the algorithm

  for(int i=1;i<nsteps;i++){
    double k1p = h*vel[i-1];
    double k1v = -1.0*h*a*pos[i-1]/pow(r,3);
    
    double k2p = h*(vel[i-1]+k1v/2);//v(ti+h/2,yi+k1/2);
    double k2v = -h*a*(pos[i-1]+k1p/2)/pow(r,3);//x(ti+h/2,yi+k1/2)/pow(r,3);

    double k3p = h*(vel[i-1]+k2v/2);//v(ti+h/2,yi+k2/2);
    double k3v = -h*a*(pos[i-1]+k2p/2)/pow(r,3);//x(ti+h/2,yi+k1/2)/pow(r,3);

    double k4p = h*(vel[i-1]+k3v);//v(ti+h,yi+k3);
    double k4v = -h*a*(pos[i-1]+k3p)/pow(r,3);//x(ti+h,yi+k3)/pow(r,3);

    pos.point[i] = pos[i-1]+(1.0/6)*(k1p+2*k2p+2*k3p+k4p);
    vel.point[i] = vel[i-1]+(1.0/6)*(k1v+2*k2v+2*k3v+k4v);
  }

  vector<thevec> to_ret;
  to_ret.push_back(pos);
  to_ret.push_back(vel);
  return to_ret;
}

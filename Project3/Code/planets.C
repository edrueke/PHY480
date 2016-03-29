/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the definition file for the planet class used to solve the 
solar system problems.  The class works with the thevec and themat classes from 
classes.C and classes.h developed in Projects 1 and 2.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

#include "planets.h"
#include "classes.h"

using namespace std;

const double PI = 3.14159265358979;

//Planet functions

planet::planet(string n, double m, double d){
  /*
    Constructor which initializes the planet (or any other orb-shaped
    celestial object).
  */
  
  mass = m;
  name = n;
  dist_sun = d;
  
  //Calculate the acceleration
  acc = 4*PI*PI/(dist_sun*dist_sun);
  
}

planet::planet(const planet &p){ 
  /*
    Copy constructor
  */
  
  mass = p.mass;
  name = p.name;
  dist_sun = p.dist_sun;
  acc = p.acc;
}

string planet::print(){
  /*
    Returns a string with the planet information.
  */

  string to_ret = "planet: "+name+"; mass: "+to_string(mass)+"kg; distance to sun: "+to_string(dist_sun)+"AU";

  return to_ret;
}

double planet::kinetic(double vel){
  /*
    Takes in the velocity at some point and returns the kinetic energy.
  */

  double T = 0.5*mass*pow(vel,2);
  return T;
}

double planet::potential(double pos){
  /*
    Takes in the position at some point and returns the potential energy.
  */

  double V = mass*4*pow(PI,2)/pow(pos,2);
  return V;
}

double planet::ang_mom(double vel, double pos){
  /*
    Takes in the position and velocity at some point and returns the
    angular momentum.
  */

  double l = mass*vel*pos;
  return l;
}

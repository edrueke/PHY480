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
  
  v0 = 0;
  positionsx_v = thevec(1);
  positionsy_v = thevec(1);
  velocitiesx_v = thevec(1);
  velocitiesy_v = thevec(1);
  positionsx_r = thevec(1);
  positionsy_r = thevec(1);
  velocitiesx_r = thevec(1);
  velocitiesy_r = thevec(1);

}

planet::planet(string n, double m, double d, double v){
  /*
    Constructor which initalizes the planet (including v0)
    
    Necessary for solar system
  */
  
  mass = m;
  name = n;
  dist_sun = d;
  v0 = v;
  
  acc = 4*PI*PI/(dist_sun*dist_sun);

  positionsx_v = thevec(1);
  positionsy_v = thevec(1);
  velocitiesx_v = thevec(1);
  velocitiesy_v = thevec(1);
  positionsx_r = thevec(1);
  positionsy_r = thevec(1);
  velocitiesx_r = thevec(1);
  velocitiesy_r = thevec(1);

}

planet::planet(const planet &p){ 
  /*
    Copy constructor
  */
  
  mass = p.mass;
  name = p.name;
  dist_sun = p.dist_sun;
  acc = p.acc;
  v0 = p.v0;
  positionsx_v = p.positionsx_v;
  positionsy_v = p.positionsy_v;
  velocitiesx_v = p.velocitiesx_v;
  velocitiesy_v = p.velocitiesy_v;
  positionsx_r = p.positionsx_r;
  positionsy_r = p.positionsy_r;
  velocitiesx_r = p.velocitiesx_r;
  velocitiesy_r = p.velocitiesy_r;

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

  double V = mass*4*pow(PI,2)/pos;
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

bool operator==(const planet &p1, const planet &p2){
  /*
    Determine whether two planets are the same.
  */
  
  return(p1.name == p2.name);
}

bool operator!=(const planet &p1, const planet &p2){
  /*
    Determine whether two planets are not the same.
  */

  return !(p1==p2);
}
  

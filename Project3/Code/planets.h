#ifndef PLANETS_H
#define PLANETS_H

/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the header file for the planet class used to solve the 
solar system problems.  The class works with the thevec and themat classes from 
classes.C and classes.h developed in Projects 1 and 2.
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"

#include "classes.h"

using namespace std;

const double PI = 3.14159265358979;

//Planet class
class planet{
 public:
  //Components
  double mass; //mass of planet in kg
  double dist_sun; //distance from sun of planet in AU
  string name; //name of planet
  double acc; //acceleration of planet on a circular orbit
  double v0; //Initial velocity
  thevec positionsx_v; //Verlet positions
  thevec positionsy_v; //Verlet positions
  thevec velocitiesx_v; //Verlet velocities
  thevec velocitiesy_v; //Verlet velocities
  thevec positionsx_r; //RK4 positions
  thevec positionsy_r; //RK4 positions
  thevec velocitiesx_r; //RK4 velocities
  thevec velocitiesy_r; //RK4 velocities
  
  //Constructors
  planet(string n, double m, double d); //Initializer
  planet(string n, double m, double d, double v); //Initializer
  planet(const planet &p); //Copy constructor
  ~planet() {}; //Destructor

  //Other functions and operators
  string print(); //String for cout-ing
  double kinetic(double vel); //Calculates T for some velocity
  double potential(double pos); //Calculates V for some position
  double ang_mom(double vel, double pos); //Calculates l for some x, v

  //Standard operations
  friend bool operator!=(const planet &p1, const planet &p2); //Comparison
  friend bool operator==(const planet &p1, const planet &p2); //Comparison
};

bool operator !=(const planet &p1, const planet &p2); //Comparison
bool operator ==(const planet &p1, const planet &p2); //Comparison
#endif

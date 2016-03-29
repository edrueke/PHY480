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

//Planet class
class planet{
 public:
  //Components
  double mass; //mass of planet in kg
  double dist_sun; //distance from sun of planet in AU
  string name; //name of planet
  double acc; //acceleration of planet on a circular orbit
  //vector for velocities?
  //vector for positions?
  
  //Constructors
  planet(string n, double m, double d); //Initializer
  planet(const planet &p); //Copy constructor
  ~planet() {}; //Destructor

  //Other functions and operators
  string print(); //String for cout-ing
  double kinetic(double vel); //Calculates T for some velocity
  double potential(double pos); //Calculates V for some position
  double ang_mom(double vel, double pos); //Calculates l for some x, v

};

#endif

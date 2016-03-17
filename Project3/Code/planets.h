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

using namespace std;

//Planet class
class planet{
 public:
  //Components
  double mass; //mass of planet
  double dist_sun; //distance from sun of planet
  string name; //name of planet
  //vector for velocities?
  //vector for positions?
  
  //Constructors
  planet(string n, double m, double d);

  //Other functions and operators
  //kinetic energy?
  //potential energy?

};

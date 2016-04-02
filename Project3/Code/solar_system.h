#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H

/*
Elizabeth Drueke
PHY 480
April 1, 2016

This is the header file for the solar system class which inherits objects
from the planet class and uses the odesolvers to compute the velocity and 
positions of the various planets added.
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "classes.h"
#include "planets.h"

using namespace std;

class solar_system{
 public:

  //Components
  vector<planet*> planets; //Vector of planets
  int nsteps; //Number of steps for the calculations
  double tf; //Final time for calculations

  //Constructors
  solar_system(planet &p); //Construct from single planet
  solar_system(vector<planet*> ps); //Construct from list of planets
  solar_system(const solar_system &ss); //Copy constructor
  ~solar_system() {}; //Destructor

  //Other functions and operators.
  void Add(planet &p); //Add a planet to the mix.
  void Solve_Verlet(); //Solve once all planets have been added.
  //void Solve_RK4(int nsteps,double tf); //Solve once all planets have been added.
  void Draw_Verlet(string name); //Draw and save to plots/ directory
  void Draw_RK4(string name); //Draw and save to plots/ directory
};

#endif

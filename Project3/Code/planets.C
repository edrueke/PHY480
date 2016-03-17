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

#include "classes.C"
#include "planets.h"

using namespace std;

//Planet functions

planet::planet(string n, double m, double d){
  /*
    Constructor which initializes the planet.
  */

  mass = m;
  name = n;
  dist_sun = d;
}

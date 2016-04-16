#ifndef PLANETS_H
#define PLANETS_H

/*
Elizabeth Drueke
PHY 480
April 29, 2016
Project 4
This is the header file for class lattice, which simulates a ferromagnetic
lattice and solves for statistical quantities using the Ising model.
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

const double PI = 3.14159265358979;
const double kB = 1.38064852e-23;

//Planet class
class lattice{
 public:
  
  //Components
  int size; //Size n of nxn lattice
  double **spins; //Matrix of spins
  double temp; //Temperature of system (in K)
  double *averages[5]; //Average (E,e^2,M,M^2,abs(M)
  int MCcycles; //# of MC cycles to use in calculations
  double CV; //Specific heat
  double chi; //susceptibility

  //Constructors/Destructors
  lattice(); //Default constructor
  lattice(const lattice &p); //Copy constructor
  ~lattice(); //Destructor
  
  lattice(int sz); //Construct from a lattice size
  lattice(double t); //Construct from a temperature
  lattice(int sz, double t); //Construct from a temperature and lattice size
  //lattice(int sz, double t, int MC); //Construct from a temp, lattice size, and number of MC cycles

  //Important calculations
  void calc_stat_quants(); //Use Metropolis algorithm to calculate important statistical quantities
  double nearest_neighbors(int row, int col); //Calculate deltaE from nearest neighbors
  void calc_spec_heat(); //Calculate the specific heat
  void calc_susceptibility(); //Calculate the susceptibility

};

#endif

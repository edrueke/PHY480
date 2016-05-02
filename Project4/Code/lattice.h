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
#include <map>

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
 private:
  
  //Components
  int size; //Size n of nxn lattice
  double **spins; //Matrix of spins
  double temp; //Temperature of system (in K)
  double *averages; //Average (E,E^2,M,M^2,abs(M)
  int MCcycles; //# of MC cycles to use in calculations
  double CV; //Specific heat
  double chi; //susceptibility
  double chi_absm; //susceptibility from <|M|>
  int accepted; //Accepted events in MC simulation
  map<double,int> e_probs; //Count of number of times an energy appears in the calculation

  //Important calculations
  void calc_stat_quants(); //Use Metropolis algorithm to calculate important statistical quantities

 public:

  //Constructors/Destructors
  lattice(); //Default constructor
  lattice(const lattice &p); //Copy constructor
  ~lattice(); //Destructor
  
  //lattice(int sz, double t, int opt=0); //Construct from a temperature and lattice size
  lattice(int sz, double t, int MC, int opt=0); //Construct from a temp, lattice size, and number of MC cycles

  //Important calculations
  double get_E(); //Get expectation value of E
  double get_M(); //Get expectation value of M
  double get_E2(); //Get expectation value of E^2
  double get_M2(); //Get expectation value of M^2
  double get_absM(); //Get expectation value of abs(M)
  double get_CV(); //Get the specific heat
  double get_susc(); //Get the susceptibility
  double get_susc_absm(); //Get the susceptibility from <|M|>
  int get_accepted(); //Get the accepted # of events in MC simulation
  void set_temp(double t); //Reset the temperature
  void set_MC(int MC); //Reset the number of MC cycles
  void plot_eprobs(string name); //Plot P(E)

};

double nearest_neighbors(int row, int col, int size, double **spins); //Calculate deltaE from nearest neighbors

#endif

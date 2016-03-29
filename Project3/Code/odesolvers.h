#ifndef ODESOLVERS_H
#define ODESOLVERS_H

/*
Elizabeth Drueke
PHY 480
April 1, 2016

Project 3

This is the header file for the ODE solvers required when solving the 
solar system problems.  They work with the thevec and themat classes from 
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

vector<thevec> Verlet(double t0, double tf, int nsteps, double x0, double xf, double a, double v0, double vf, double r);
vector<thevec> RK4(double t0, double tf, int nsteps, double x0, double xf, double a, double v0, double vf, double r);

#endif

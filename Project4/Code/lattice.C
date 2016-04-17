/*
Elizabeth Drueke
PHY 480
April 29, 2016
Project 4
This is the definition file for the lattice class.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"

using namespace std;

//Constructors/Destructors

lattice::lattice(){
  /*
    Default constructor.
  */

  size = 2;
  temp = 1;
  MCcycles = 5;
  accepted = 0;

  averages = new double[5];
  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double*[size];
  for(int i=0;i<size;i++){
    spins[i] = new double[size];
  }

  //Initialize the lattice with random seeds
  srand(time(NULL));
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      int rando = rand()%2;
      if(rando == 0)
	spins[i][j] = -1;
      else
	spins[i][j] = 0;
    }
  }

  (*this).calc_stat_quants();
  for(int i=0;i<5;i++){
    averages[i]=averages[i]/MCcycles;
  }

  CV = (1.0/(kB*pow(temp,2)))*(averages[1]-pow(averages[0],2));
  chi = (1.0/(kB*temp))*(averages[3]-pow(averages[2],2));

}

lattice::lattice(const lattice &p){
  /*
    Copy constructor.
  */

  size = p.size;
  temp = p.temp;
  MCcycles = p.MCcycles;
  CV = p.CV;
  chi = p.chi;
  accepted = p.accepted;

  averages = new double[5];

  for(int i=0;i<5;i++)
    averages[i]=p.averages[i];

  spins = new double*[size];
  for(int i=0;i<size;i++){
    spins[i] = new double[size];
  }

  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      spins[i][j] = p.spins[i][j];
    }
  }

}

lattice::~lattice(){
  /*
    Destructor
  */

  for(int i=0;i<size;i++){
    delete [] spins[i];
  }

  delete [] spins;
  delete [] averages;

}

//lattice::lattice(int sz, double t, int opt=0){
  /*
    Construct a lattice from the temp t and size sz
  */

/*  size = sz;
  temp = t;
  MCcycles = 5;
  accepted = 0;

  averages = new double[5];

  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double *[size];
  
  for(int i=0;i<size;i++){
    spins[i] = new double[size];
  }

  //Initialize the lattice with random seeds
  if(opt==0){
    srand(time(NULL));
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	int rando = rand()%2;
	if(rando == 0)
	  spins[i][j] = -1;
	else
	  spins[i][j] = 0;
      }
    }
  }

  else{
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	spins[i][j]=opt;
      }
    }
  }

  (*this).calc_stat_quants();

  for(int i=0;i<5;i++)
    averages[i]=averages[i]/MCcycles;

  CV = (1.0/(kB*pow(temp,2)))*(averages[1]-pow(averages[0],2));
  chi = (1.0/(kB*temp))*(averages[3]-pow(averages[2],2));

  }*/

lattice::lattice(int sz, double t, int MC, int opt=0){
  /*
    Construct a lattice from the temp t and size sz and number of MC cycles MC
  */

  size = sz;
  temp = t;
  MCcycles = MC;
  accepted = 0;

  averages = new double[5];

  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double *[size];
  
  for(int i=0;i<size;i++){
    spins[i] = new double[size];
  }

  if(opt==0){
    //Initialize the lattice with random seeds
    srand(time(NULL));
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	int rando = rand()%2;
	if(rando == 0)
	  spins[i][j] = -1;
	else
	  spins[i][j] = 0;
      }
    }
  }

  else{
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	spins[i][j]=opt;
      }
    }
  }

  (*this).calc_stat_quants();

  for(int i=0;i<5;i++)
    averages[i]=averages[i]/MCcycles;

  CV = (1.0/(kB*pow(temp,2)))*(averages[1]-pow(averages[0],2));
  chi = (1.0/(kB*temp))*(averages[3]-pow(averages[2],2));

}

void lattice::calc_stat_quants(){
  /*
    Implement the Metropolis algorithm to calculate the expectation values of
    various quantities
  */

  //Change the configuration of the spins by flipping one spin only.
  srand(time(NULL));

  //Initialize energy and magnetization
  double E=0; double M=0;

  int mod = size*size;

  for(int cycle = 0;cycle<MCcycles;cycle++){
    int rando = rand()%mod; //Which spin to change
    int ct = 0;
    int therow = 0; int thecol = 0;
    
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	if(ct==rando){
	  spins[i][j]*=-1;
	  therow = i;
	  thecol = j;
	}
	ct+=1;
      }
    }
    
    //Calculate delta-E using nearest neighbors
    double deltaE = nearest_neighbors(therow,thecol,size,spins);
    
    //If deltaE>0, compare with a random number to determine if we should keep
    //the new system
    bool go = true;
    if(deltaE>0){
      double w = exp(-1.0*deltaE/(kB*temp));
      int myrand = rand()%10; //HERE: What range? 
      if(myrand>w){
	spins[therow][thecol]*=-1;
	deltaE=0;
	go = false;
      }
    }

    if(go){
      accepted++;
      M+= 2*spins[therow][thecol];
      E+= deltaE;
    }  
    
    averages[0]+=E;
    averages[1]+=pow(E,2);
    averages[2]+=M;
    averages[3]+=pow(M,2);
    averages[4]+=abs(M);
  }
}
  
double nearest_neighbors(int row, int col,int size, double **spins){
  /*
    Calculate deltaE from the nearest neighbors
  */
  
  int ups=0; int downs=0;
  
  //First deal with row
  if(row==size-1){
    if(spins[0][col]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row-1][col]==1)
      ups+=1;
    else
      downs+=1;
  }
  else if(row==0){
    if(spins[size-1][col]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row+1][col]==1)
      ups+=1;
    else
      downs+=1;
  }
  else{
    if(spins[row-1][col]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row+1][col]==1)
      ups+=1;
    else
      downs+=1;
  }

  //Then deal with col
  if(col==size-1){
    if(spins[row][0]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row][col-1]==1)
      ups+=1;
    else
      downs+=1;
  }
  else if(col==0){
    if(spins[row][size-1]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row][col+1]==1)
      ups+=1;
    else
      downs+=1;
  }
  else{
    if(spins[row][col-1]==1)
      ups+=1;
    else
      downs+=1;
    if(spins[row][col+1]==1)
      ups+=1;
    else
      downs+=1;
  }

  if(ups+downs!=4)
    cout<<"ERROR: ups+downs = "<<ups+downs<<endl;

  //Based on the number of ups and downs and the spin of the final state, 
  //calculate deltaE.
  double deltaE = (-4*ups+8)*spins[row][col];
  return deltaE;
}

double lattice::get_E(){
  /*
    Return <E>
  */
  
  return averages[0];
}

double lattice::get_M(){
  /*
    Return <M>
  */

  return averages[2];
}

double lattice::get_E2(){
  /*
    Return <E^2>
  */

  return averages[1];
}

double lattice::get_M2(){
  /*
    Return <M^2>
  */

  return averages[3];
}

double lattice::get_absM(){
  /*
    Return <|M|>
  */

  return averages[4];
}

double lattice::get_CV(){
  /*
    Return specific heat
  */

  return CV;
}

double lattice::get_susc(){
  /*
    Return susceptibility
  */

  return chi;
}

void lattice::set_temp(double t){
  /*
    Reset the temperature and recalculate everything
  */

  lattice newlat = lattice(size,t,MCcycles);
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      spins[i][j] = newlat.spins[i][j];
    }
  }
  temp = t;
  for(int i=0;i<5;i++)
    averages[i] = newlat.averages[i];
  CV = newlat.CV;
  chi = newlat.chi;
}

void lattice::set_MC(int MC){
  /*
    Reset the number of MC cycles to use in the calculation
  */

  lattice newlat = lattice(size,temp,MC);

  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      spins[i][j] = newlat.spins[i][j];
    }
  }
  MCcycles = MC;
  for(int i=0;i<5;i++)
    averages[i] = newlat.averages[i];
  CV = newlat.CV;
  chi = newlat.chi;
}

int lattice::get_accepted(){
  /*
    Return the number of accepted events in the MC simulation
  */

  return accepted;
}

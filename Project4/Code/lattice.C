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

  calc_stat_quants();
  for(int i=0;i<5;i++){
    averages[i]=averages[i]/MCcycles;
  }

  calc_spec_heat();
  calc_susceptibility();

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

lattice::lattice(int sz){
  /*
    Construct a lattice of a given size sz.
  */

  size = sz;
  temp = 1;
  MCcycles = 5;

  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double *[size];

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

  calc_stat_quants();
  for(int i=0;i<5;i++)
    averages[i]=averages[i]/MCcycles;

  calc_spec_heat();
  calc_susceptibility();
}

lattice::lattice(double t){
  /*
    Construct a lattice from the temp t
  */

  size = 2;
  temp = t;
  MCcycles = 5;

  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double *[size];
  
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

  calc_stat_quants();
  for(int i=0;i<5;i++)
    averages[i]=averages[i]/MCcycles;

  calc_spec_heat();
  calc_susceptibility();
}

lattice::lattice(int sz, double t){
  /*
    Construct a lattice from the temp t and size sz
  */

  size = sz;
  temp = t;
  MCcycles = 5;

  for(int i=0;i<5;i++)
    averages[i]=0;

  spins = new double *[size];
  
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

  for(int i=0;i<5;i++)
    averages[i]=averages[i]/MCcycles;

  calc_spec_heat();
  calc_susceptibility();
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

  for(int cycle = 0;cycle<MCcycles;cycle++){
    int rando = rand()%(pow(size,2)); //Which spin to change
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
    double deltaE = nearest_neighbors(therow,thecol);
    
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
      M+= 2*spins[therow][thecol];
      E+= deltaE;
    }  
    
    averages[0]+=E;
    averages[1]+=pow(e,2);
    averages[2]+=M;
    averages[3]+=pow(M,2);
    averages[4]+=abs(M);
  }
}
  
double nearest_neighbors(int row, int col){
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
      downss+=1;
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

void calc_spec_heat(){
  /*
    Calculate the set the specific heat from the averages
  */

  CV = (1.0/(kB*pow(temp,2)))*(averages[1]-pow(averages[0],2));

}

void calc_susceptibility(){
  /*
    Calculate the susceptibility from the averages
  */
  
  chi = (1.0/(kB*temp))*(averages[3]-pow(averages[2],2));
}

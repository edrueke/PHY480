/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This script solves the 1d Poisson Equation -u''(x)=f(x) for 0<x<1 subject to
the Dirichlet boundary conditions u(0)=u(1)=0 for f(x)=100e**(-10*x).
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

//dim is a global variable giving the dimension of the nxn matrices
const dim = 2;

double first_deriv(double &x, double &h){
  /*Takes in 
    and returns the first derivative*/
  
  return 0;
}

double second_deriv(double &x, double &h){
  /*Takes in
    and returns the second derivative.*/

  return 0;
}

double mat_mat_mult(double matr1[][dim], double matr2[][dim]){
  /*Takes in 2 matrices and returns their multiplication*/

  double matr[dim][dim];
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      matr[i][j]=0;
      for(int k=0;k<dim;k++){
	matr[i][j]+=matr1[i][k]*matr2[k][j];
      }
    }
  }

  return matr;
}

double mat_vec_mult(double matr1[][dim], double vec[]){
  /*Takes in a matrix and a vector and returns their multiplication.*/

  double vec[dim];
  for(int i=0;i<dim;i++){
    vec[i]=0;
    for(int j=0;j<dim;j++){
      vec[i]+=matr[i][k]*vec[k];
    }
  }

  return 0;
}

//What if there is an arbitrary 0 entry? 
double one_forw_reduc(double matr[][dim], double vec[]){
  /*Takes in a matrix and a vector and performs one forward reduction, whose
    matrix and vector it then returns.*/
  
  //Determine what row needs to be reduced.
  int i=0;
  while(i+1<dim and matr[i+1][i]==0){
    i++;
  }

  //If the matrix is fully reduced (forward), then return the matrix.
  if(i+1==dim)
    return matr,vec;
  
  //Else, reduce by the next row.
  double matr_ret[dim][dim];
  double vec_ret[dim];
  for(int j=0;j<i;j++){
    vec_ret[j]=vec[j];
    for(int k=0;k<dim;k++){
      matr_ret[j][k]=matr[j][k];
    }
  }

  for(int j=i;j<dim;j++){
    double fact = matr[j][i]/matr[i][i];
    vec_ret[j]=vec[j]-vec[i]*fact;
    for(int k=0;k<dim;k++){
      if(matr[j][k]==0)
	matr_ret[j][k]=0;
      else{
	matr_ret[j][k]=matr[j][k]-matr[i]*fact;
      }
    }
  }

  //Can I return this? Doubtful. Probably need to work solely with pointers.
  return matr_ret,vec_ret;
}

//What if there is an arbitrary 0 entry? 
double one_back_reduc(double matr[][dim], double vec[]){
  /*Takes in a matrix and a vector and performs one backward reduction, whose
    matrix and vector it then returns.*/
  
  //Determine what row needs to be reduced.
  int i=dim-1;
  while(i-1>-1 and matr[i-1][i]==0){
    i--;
  }

  //If the matrix is fully reduced (forward), then return the matrix.
  if(i-1==0)
    return matr,vec;
  
  //Else, reduce by the next row.
  double matr_ret[dim][dim];
  double vec_ret[dim];
  for(int j=dim-1;j>i;j--){
    vec_ret[j]=vec[j];
    for(int k=0;k<dim;k++){
      matr_ret[j][k]=matr[j][k];
    }
  }

  for(int j=i;j>0;j--){
    double fact = matr[j][i]/matr[i][i];
    vec_ret[j]=vec[j]-vec[i]*fact;
    for(int k=0;k<dim;k++){
      matr_ret[j][k]=0;
    }
  }

  //Can I return this? Doubtful. Probably need to work solely with pointers.
  return matr_ret,vec_ret;
}

double gauss_elim(double matr1[][dim],double vec[]){
  /*Takes in a matrix and a vector (the solution to Ax=b) and returns the 
    vector x which is the solution to this problem, found by Gaussian 
    elimination.*/

  //First to perform forward elimination.
  double matr_ret[dim][dim] = matr1;
  double vec_ret[dim] = vec;
  for(int i=0;i<dim;i++){
    matr_ret, vec_ret = one_forw_reduc(matr_ret,vec_ret);
  }

  //Then perform backward elimination.
  for(int i=0;i<dim;i++){
    matr_ret, vec_ret = one_back_reduc(matr_ret,vec_ret);
  }

  //Then compute the solution x
  double x[dim];
  for(int i=0;i<dim;i++){
    x[i]=vec_ret[i]/matr_ret[i];
  }

  return vec_ret;
}

//Need to double-check this algorithm.
double LU_decomp(double matr[][dim]){
  /*Takes in a matrix and returns the LU decomposition as a pair of pointers 
    (L,U) (?).*/

  double L[dim][dim];
  double U[dim][dim];

  //Predefined 0's and 1's:


  //First, U_1j = matr_1j
  for(int i=0;i<dim;i++){
    U[1][i]=matr[1][i];
  }

  //Then U_ij = matr_ij - sum(l_ik*u_kj,k,1,i-1)
  for(int j=0;j<dim;j++){
    for(int i=0;i<j;i++){
      U[i][j]=matr[i][j];
      for(int k=0;k<i-2;k++){
	U[i][j]-=L[i][k]*U[k][j];
      }
    }
  }

  //Then U_jj = ajj - sum(l_jk*u_kj,k,1,j-1)
  for(int j=0;j<dim;j++){
    U[j][j]=matr[j][j];
    for(int k=0;k<j-2;k++){
      U[j][j]-=L[j][k]*U[k][j];
    }
  }

  //Finally, for i>j, L_ij = (1/U_ij)*(matr_ij - sum(l_ik*u_kj,k,1,i-1))
  for(int j=0;j<dim;j++){
    for(int i=j+1;i<dim;i++){
      L[i][j]=matr[i][j]/U[j][j];
      for(int k=0;k<i-2;k++){
	L[i][j]-=L[i][k]*U[i][k]/U[j][j];
      }
    }
  }

  return 0;
}

/*
  Next steps: test all the functions.
*/

/*
Elizabeth Drueke
PHY 480
March 4, 2016

Project 2

This was the initial file created to test the classes and attempt to determine
a solution to the problem presented in project 2.
*/

#include <iostream>

#include "classes.C"

using namespace std;

void prelims(){
  /*
    The main function which calls other functions and attempts to solve the 
    problem.
  */

  //Define the p_min, p_max, and n_steps
  double p_min=0.0; double p_max=10.0;
  int n_step = 3.0;
  double h = (p_max-p_min)/n_step;

  themat my_mat = themat(n_step-1);

  for(int i=0;i<p_max+1;i++){
    
    //p_i = p_min+i*h, i = 0,1,...,n
    double p_i = p_min+i*h;

    //Let v_i=p_i^2
    double V_i = pow(p_i,2);

    //Diagonal matrix element d_i = 2/h^2+v_i
    double d_i = 2/(h*h)+V_i;

    //Off-diagonal matrix element e_i = -1/h^2
    double e_i = -1/(h*h);

    //Define the matrix    
    for(int j=0;j<my_mat.sz;j++){
      for(int k=0;k<my_mat.sz;k++){
	if(j==k)
	  my_mat.point[j][k]=d_i;
	else if(j==k-1 || j==k+1)
	  my_mat.point[j][k]=e_i;
	else
	  my_mat.point[j][k]=0;
      }
    }
  }

  //Diagonalize the matrix

  thevec vec = thevec(2);
  themat matr = themat(2);

  //Fill the matrix and vector with values.
  for(int i=0;i<2;i++){
    vec.point[i]=1;
  }
  matr.point[0][0]=1; matr.point[0][1]=1; 
  matr.point[1][0]=2; matr.point[1][1]=-1;

  themat diag = Jacobi_Method(matr,1e-8);
  cout<<"Matrix: "<<endl<<matr.print()<<endl<<"Diag:"<<endl<<diag.print()<<endl;
  


}

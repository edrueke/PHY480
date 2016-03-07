/*
Elizabeth Drueke
PHY 480
March 4, 2016

Project 2

This was the initial file created to test the classes and attempt to determine
a solution to the problem presented in project 2. It defines benchmarks which 
should be met by the code.
*/

#include <iostream>
#include <fstream>

#include "classes.C"

using namespace std;

void prelims(){
  /*
    The main function which calls other functions and attempts to solve the 
    problem.
  */

  //Define an output file.
  ofstream myfile;
  myfile.open("plots/benchmarks.txt");
  myfile<<"Check the Jacobi Rotation Algorithm"<<endl;
  
  //Define and diagonalize matrix 1
  themat mat1 = themat(2);
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      if(i==j)
	mat1.point[i][j]=1;
      else
	mat1.point[i][j]=2;
    }
  }
  
  themat diag1 = Jacobi_Method(mat1,1e-8);
  myfile<<endl<<"Matrix 1:"<<endl<<mat1.print()<<endl;
  myfile<<"Expected Eigenvalues: "<<3<<", "<<-1<<endl<<endl;
  myfile<<"Diagonalized:"<<endl<<diag1.print();

  //Define and diagonalize matrix 2
  themat mat2 = themat(3);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(i==j)
	mat2.point[i][j]=1;
      else
	mat2.point[i][j]=2;
    }
  }
  
  themat diag2 = Jacobi_Method(mat2,1e-8);
  myfile<<endl<<"Matrix 1:"<<endl<<mat2.print()<<endl;
  myfile<<"Expected Eigenvalues: "<<5<<", "<<-1<<", "<<-1<<endl<<endl;
  myfile<<"Diagonalized:"<<endl<<diag2.print();

  //Define and diagonalize matrix 3
  themat mat3 = themat(10);
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      if(i==j)
	mat3.point[i][j]=1;
      else
	mat3.point[i][j]=2;
    }
  }
  
  themat diag3 = Jacobi_Method(mat3,1e-8);
  myfile<<endl<<"Matrix 1:"<<endl<<mat3.print()<<endl;
  myfile<<"Expected Eigenvalues: "<<19<<" and nine "<<-1<<"'s"<<endl<<endl;
  myfile<<"Diagonalized:"<<endl<<diag3.print();

  //Define and diagonalize matrix 4
  themat mat4 = themat(10);
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      if(i==j)
	mat4.point[i][j]=1;
      else if(i<j)
	mat4.point[i][j]=2.0/(i+1);
      else
	mat4.point[i][j]=2.0/(j+1);
    }
  }
  
  themat diag4 = Jacobi_Method(mat4,1e-8);
  myfile<<endl<<"Matrix 1:"<<endl<<mat4.print()<<endl;
  myfile<<"Expected Eigenvalues: 9.8166, -3.36209, 0.777778, 0.742703, 0.694927, 0.617823, 0.579657, -0.420483, 0.4059, 0.147189" <<endl<<endl;
  myfile<<"Diagonalized:"<<endl<<diag4.print();
  myfile<<"Eigenvalues:"<<endl;
  for(int i=0;i<diag4.sz;i++)
    myfile<<"     "<<diag4[i][i];

  myfile.close();
}

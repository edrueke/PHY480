/*
Elizabeth Drueke
PHY 480
February 12, 2016

Project 1

This is the definition file for the matrix and vector classes defined for 
project 1.
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "classes.h"

using namespace std;

//Vector functions

thevec::thevec(int s){
  /*
    Initialize a dynamic memory array of size sz.
  */

  point = new double[s];
  sz = s;
}

thevec::thevec(vector<double> vec){
  /*
    Initialize a dynamic array from a vector of doubles.
  */
  
  sz = vec.size();
  point = new double[sz];

  for(int i=0;i<sz;i++){
    point[i]=vec.at(i);
  }
}

thevec::thevec(const thevec &vec){
  /*
    Copy constructor.
  */

  sz = vec.sz;
  point = new double[sz];
  for(int i=0;i<sz;i++){
    point[i] = vec.point[i];
  }
}

string thevec::print(){
  /*
    Return a string of the components of the vector.
  */

  string toret = "";
  for(int i=0;i<sz;i++){
    toret+=(to_string(point[i])+",");
  }
  return toret;
}

thevec operator+(const thevec &vec1, const thevec &vec2){
  /*
    Addition of two vectors of the same size
  */

  int sz = vec1.sz;
  thevec to_ret = thevec(sz);
  for(int i=0;i<sz;i++){
    to_ret.point[i] = vec1.point[i]+vec2.point[i];
  }

  return to_ret;
}
  
thevec operator-(const thevec &vec1, const thevec &vec2){
  /*
    Subtraction of two vectors of the same size
  */

  int sz = vec1.sz;
  thevec to_ret = thevec(sz);
  for(int i=0;i<sz;i++){
    to_ret.point[i] = vec1.point[i]-vec2.point[i];
  }

  return to_ret;
}

double operator*(const thevec &vec1, const thevec &vec2){
  /*
    Dot product of two vectors of the same size.
  */

  int sz = vec1.sz;
  double ret = 0;
  for(int i=0;i<sz;i++){
    ret+=vec1.point[i]*vec2.point[i];
  }

  return ret;
}

double thevec::operator[](int i){
  /*
    Component retrieval.
  */

  return point[i];
}

thevec &thevec::operator=(const thevec &vec){
  /*
    Sets this equal to some vector.
  */

  for(int i=0;i<sz;i++){
    this->point[i]=vec.point[i];
  }
  
  return *this;
}

thevec::~thevec(){
  /*
    Delete the dynamic memory.
  */

  delete [] point;
}

//Matrix functions

themat::themat(int s){
  /*
    Initialize a matrix of size sz x sz using dynamic memory.
  */

  sz = s;
  point = new double* [sz];
  for(int i=0;i<sz;i++){
    this->point[i] = new double [sz];
  }
}

themat::themat(vector<double> vec){
  /*
    Initialize a matrix from a vector of the components.
  */

  sz = sqrt(vec.size());

  point = new double* [sz];
  for(int i=0;i<sz;i++){
    point[i] = new double[sz];
  }

  int row = 0;
  int column = 0;
  for(unsigned int i=0;i<vec.size();i++){
    if(i%sz==0 && i!=0){
      row+=1;
      column=0;
    }
    point[row][column] = vec.at(i);
    column+=1;
  }
}  

themat::themat(const themat &mat){
  /*
    Copy constructor.
  */

  sz = mat.sz;
  point = new double*[sz];
  for(int i=0;i<sz;i++){
    point[i] = new double[sz];
  }

  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      point[i][j]=mat.point[i][j];
    }
  }
}

string themat::print(){
  /*
    Return a string with the components of the vector.
  */

  string to_ret = "";

  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      to_ret+=(to_string(point[i][j])+",");
    }
    to_ret+="\n";
  }
  return to_ret;
}

themat operator+(const themat &mat1, const themat &mat2){
  /*
    Addition of two matrices of the same size.
  */

  int sz = mat1.sz;
  themat ret = themat(sz);
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      ret.point[i][j] = mat1.point[i][j] + mat2.point[i][j];
    }
  }

  return ret;
}

themat operator-(const themat &mat1, const themat &mat2){
  /*
    Subtraction of two matrices of the same size.
  */

  int sz = mat1.sz;
  themat ret = themat(sz);
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      ret.point[i][j] = mat1.point[i][j]-mat2.point[i][j];
    }
  }

  return ret;
}

themat operator*(const themat &mat1, const themat &mat2){
  /*
    Matrix-matrix multiplication.
  */

  int sz = mat1.sz;
  themat ret = themat(sz);

  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      ret.point[i][j]=0;
      for(int k=0;k<sz;k++){
	ret.point[i][j]+=mat1.point[i][k]*mat2.point[k][j];
      }
    }
  }

  return ret;
}

thevec themat::operator[](int i){
  /*
    Component retrieval.
  */

  thevec toret = thevec(sz);
  for(int j=0;j<sz;j++){
    toret.point[j] = point[i][j];
  }
  return toret;
}

thevec operator*(const themat &mat, const thevec &vec){
  /*
    Mathematically defined matrix-vector multiplication.
  */

  int sz = mat.sz;
  thevec ret = thevec(sz);

  for(int i=0;i<sz;i++){
    ret.point[i]=0;
    for(int j=0;j<sz;j++){
      ret.point[i]+=mat.point[i][j]*vec.point[j];
    }
  }

  return ret;
}

themat &themat::operator=(const themat &mat){
  /*
    Takes in a matrix and sets this equal to that matrix.
  */

  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      point[i][j]=mat.point[i][j];
    }
  }
  return *this;
}

themat::~themat(){
  /*
    Delete the dynamic memory.
  */

  for(int i=0;i<sz;i++){
    delete [] point[i];
  }
  
  delete [] point;
}

//Additional functions - Gaussian elimination

//What if there is an arbitrary 0 entry? 
vector<vector<double> > one_forw_reduc(themat matr, thevec vec, int i){
  /*
    Takes in a matrix and a vector and performs one forward reduction, whose
    matrix and vector it then returns a vector of the vectors to be made into
    the matrix and vector.
  */
  
  //Determine what row needs to be reduced.
  int dim = matr.sz;
  
  //If the matrix is fully reduced (forward), then return the matrix.
  vector<double> vec_to_ret, matr_to_ret;
  vector<vector<double> > to_ret;
  
  //Else, reduce by the next row.
  themat matr_ret = themat(dim);
  thevec vec_ret = thevec(dim);

  //Take in the parts already reduced
  for(int j=0;j<i+1;j++){
    vec_ret.point[j]=vec[j];
    for(int k=0;k<dim;k++){
      matr_ret.point[j][k]=matr[j][k];
    }
  }

  for(int j=i+1;j<dim;j++){
    double fact = matr[j][i]/matr[i][i];
    vec_ret.point[j]=vec[j]-vec[i]*fact;
    for(int k=0;k<dim;k++){
      if(matr.point[j][k]==0)
	matr_ret.point[j][k]=0;
      else{
	matr_ret.point[j][k]=matr[j][k]-matr[i][k]*fact;
      }
    }
  }

  //Create the appropriate return.
  for(int k=0;k<dim;k++){
    vec_to_ret.push_back(vec_ret[k]);
    for(int j=0;j<dim;j++){
      matr_to_ret.push_back(matr_ret[k][j]);
    }
  }

  to_ret.push_back(vec_to_ret);
  to_ret.push_back(matr_to_ret);

  return to_ret;
}

//What if there is an arbitrary 0 entry? 
vector<vector<double> > one_back_reduc(themat matr, thevec vec, int i){
  /*
    Takes in a matrix and a vector and performs one backward reduction, whose
    matrix and vector it then returns a vector of the matrix and vector 
    to be turned into the dynamic memory.
  */
  
  //Determine what row needs to be reduced.
  int dim = vec.sz;

  //If the matrix is fully reduced (forward), then return the matrix.
  vector<double> vec_to_ret, matr_to_ret;
  vector<vector<double> > to_ret;

  //Else, reduce by the next row.
  themat matr_ret = themat(dim);
  thevec vec_ret = thevec(dim);

  //Take in already row-reduced parts
  for(int j=0;j<dim;j++){
    vec_ret.point[j]=vec[j];
    for(int k=0;k<dim;k++){
      matr_ret.point[j][k]=matr[j][k];
    }
  }

  for(int j=i-1;j>-1;j--){
    double fact = matr[j][i]/matr[i][i];
    vec_ret.point[j]=vec[j]-vec[i]*fact;
    for(int k=j+1;k<dim;k++){
      matr_ret.point[j][k]=matr[j][k]-fact*matr[i][k];
    }
  }

  for(int k=0;k<dim;k++){
    vec_to_ret.push_back(vec_ret[k]);
    for(int j=0;j<dim;j++){
      matr_to_ret.push_back(matr_ret[k][j]);
    }
  }

  to_ret.push_back(vec_to_ret);
  to_ret.push_back(matr_to_ret);

  return to_ret;
}

thevec gauss_elim(themat matr1,thevec vec){
  /*
    Takes in a matrix and a vector (the solution to Ax=b) and returns the 
    vector x which is the solution to this problem, found by Gaussian 
    elimination.
  */
  
  //First to perform forward elimination.
  int dim = vec.sz;
  themat matr_ret = themat(dim);
  matr_ret = matr1;
  thevec vec_ret = thevec(dim);
  vec_ret = vec;
  for(int i=0;i<dim;i++){
    vector<vector<double> > returned = one_forw_reduc(matr_ret,vec_ret,i);
    matr_ret = themat(returned.at(1));
    vec_ret = thevec(returned.at(0));
  }

  //Then perform backward elimination.
  for(int i=0;i<dim;i++){
    vector<vector<double> > returned = one_back_reduc(matr_ret,vec_ret,dim-1-i);
    matr_ret = themat(returned.at(1));
    vec_ret = thevec(returned.at(0));
  }
  //Then compute the solution x
  thevec x = thevec(dim);
  for(int i=0;i<dim;i++){
    x.point[i]=1.0*vec_ret[i]/matr_ret[i][i];
  }

  return x;
}

//Additional functions - LU Decomposition

vector<vector<double> > LU_decomp(themat matr){
  /*
    Takes in a matrix and returns the LU decomposition as a vector of vectors
    of doubles to be converted into the L and U matrices.
  */
  
  int dim = matr.sz;
  themat L = themat(dim);
  themat U = themat(dim);

  //Predefined 0's and 1's:
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j)
	L.point[i][j]=1;
      else if(i > j)
	U.point[i][j]=0;
      else
	L.point[i][j]=0;
    }
  }

  //By columns
  for(int j=0;j<dim;j++){
    
    //U_1j = A_1j
    U.point[0][j]=1.0*matr[0][j];

    //For i = 2,..,j-1, U_ij = A_ij - sum(l_ik*u_kj,k,1,i-1)
    for(int i=1;i<j;i++){
      U.point[i][j]=1.0*matr[i][j];
      for(int k=0;k<i;k++){
	U.point[i][j]-=1.0*L[i][k]*U[k][j];
      }
    }

    //Diagonal elements U_jj = A_jj - sum(L_jk*U_kj,k=1,j-1)
    U.point[j][j] = 1.0*matr[j][j];
    for(int k=0;k<j;k++){
      U.point[j][j]-=1.0*L[j][k]*U[k][j];
    }
        
    //For i>j, L_ij = (1/U_jj)*(A_ij - sum(L_ik*U_kj,k,1,j-1)
    for(int i=j+1;i<dim;i++){
      L.point[i][j]=1.0*matr[i][j];
      for(int k=0;k<j;k++){
	L.point[i][j]-=1.0*L[i][k]*U[k][j];
      }
      L.point[i][j]/=1.0*U[j][j];
    }
    
  }

  vector<vector<double> > to_ret;
  vector<double> L_ret, U_ret;

  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      L_ret.push_back(L[i][j]);
      U_ret.push_back(U[i][j]);
    }
  }

  to_ret.push_back(L_ret);
  to_ret.push_back(U_ret);

  return to_ret;
}

thevec LU_decomp_solver(themat &L, themat &U, thevec &vec){
  /*
    Take in an upper- and lower-triangular matrix and compute the solution to
    LUx = vec.
  */

  int dim = L.sz;

  //Define the intermediate vector
  thevec y = thevec(dim);

  //Solve for y
  for(int i=0;i<dim;i++){
    y.point[i]=vec[i];
    for(int j=0;j<i;j++){
      y.point[i]-=1.0*L[i][j]*y[j];
    }
  }

  //Solve for the solution
  thevec to_ret = thevec(dim);
  for(int i=dim-1;i>-1;i--){
    to_ret.point[i]=y[i];
    for(int j=i+1;j<dim;j++){
      to_ret.point[i]-=1.0*U[i][j]*to_ret[j];
    }
    to_ret.point[i]/=U.point[i][i];
  }

  return to_ret;
}

//Additional functions - Convert to string

string to_string(double d){
  
  /*
    Takes in a double and converts it to a string.  This function already 
    exists in C++11, but I don't have that version.
  */

  string to_ret;
  ostringstream convert;
  convert << d;
  to_ret = convert.str();

  return to_ret;
}


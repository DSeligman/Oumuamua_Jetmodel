//Jet Torque Model with time Delay in C++ for OUMUAMUA
//v1.0 11/20/2018
//Darryl Seligman
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

int main()
{
  int nx = 3;
  int ny = 3;
  //set up ICs
  double pi = 3.14159265359;
  double theta = pi/12.5;
  double phi=pi/9.001;
  double psi =-0.;
  
  double a=(9000.*1.5);//#cm                                                                                                                       
  double b=(4000.*1.5);
  double c=(1000.*1.5);

  //#initiatiate torque scaling alpha                                                                                                      
  double alpha= 6.e-4;

  //now set up the matrix R ICs which we will then iterate over
  double * R=NULL;
  R = new double[nx*ny];
  double * Rinv=NULL;
  Rinv = new double[nx*ny];
  
  R[0*ny+0] = cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi);

  R[1*ny+0] = -cos(psi)*cos(theta)*sin(phi)-cos(phi)*sin(psi);

  R[2*ny+0] = sin(phi)*sin(theta);


  R[0*ny+1] = cos(psi)*sin(phi)+cos(phi)*cos(theta)*sin(psi);
  R[1*ny+1] = cos(phi)*cos(psi)*cos(theta)-sin(phi)*sin(psi);
  R[2*ny+1] = -cos(phi)*sin(theta);
  R[0*ny+2] = sin(psi)*sin(theta);
  R[1*ny+2] = cos(psi)*sin(theta);
  R[2*ny+2] = cos(theta);

  
  int s, i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      cout<<R[i*ny+j]<<endl;
  cout<<'!'<<endl;
  
  gsl_matrix_view Rmat = gsl_matrix_view_array(R,nx,ny);
  gsl_matrix_view inv= gsl_matrix_view_array(Rinv,nx,ny);  
  gsl_permutation * p = gsl_permutation_alloc (nx);
  
  cout<<"about to invert"<<endl;
  gsl_linalg_LU_decomp (&Rmat.matrix, p, &s);
  gsl_linalg_LU_invert (&Rmat.matrix, p, &inv.matrix);
  cout<<"inversion done"<<endl;
  /*
  printf("The matrix is\n");
  
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      printf(j==2?"%6.3f\n":"%6.3f ", gsl_matrix_get(&Rmat.matrix,i,j));
  */
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      cout<<R[i*ny+j]<<endl;
  return 0;
}

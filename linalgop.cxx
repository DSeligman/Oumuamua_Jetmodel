
#include "linalgop.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
#include <cmath> 

using namespace std;

void invert( double *Rinv, double *R,int nx, int ny){
  
  // now we need to invert the matrix                                                                                                  
  //compute determinant                                                                                                                
  double det = R[0*ny+0]*(R[1*ny+1]*R[2*ny+2]-R[1*ny+2]*R[2*ny+1])-
    R[0*ny+1]*(R[1*ny+0]*R[2*ny+2]-R[1*ny+2]*R[2*ny+0])+
    R[0*ny+2]*(R[1*ny+0]*R[2*ny+1]-R[1*ny+1]*R[2*ny+0]);

  //build the matrix inverse                                                                                                           
  Rinv[0*ny+0]=1./det*(R[1*ny+1]*R[2*ny+2]-R[1*ny+2]*R[2*ny+1]);
  Rinv[0*ny+1]=1./det*(R[0*ny+2]*R[2*ny+1]-R[0*ny+1]*R[2*ny+2]);
  Rinv[0*ny+2]=1./det*(R[0*ny+1]*R[1*ny+2]-R[1*ny+1]*R[0*ny+2]);

  //2nd row                                                                                                                            
  Rinv[1*ny+0]=1./det*(R[1*ny+2]*R[2*ny+0]-R[2*ny+2]*R[1*ny+0]);
  Rinv[1*ny+1]=1./det*(R[0*ny+0]*R[2*ny+2]-R[0*ny+2]*R[2*ny+0]);
  Rinv[1*ny+2]=1./det*(R[0*ny+2]*R[1*ny+0]-R[0*ny+0]*R[1*ny+2]);

  Rinv[2*ny+0]=1./det*(R[1*ny+0]*R[2*ny+1]-R[1*ny+1]*R[2*ny+0]);
  Rinv[2*ny+1]=1./det*(R[0*ny+1]*R[2*ny+0]-R[0*ny+0]*R[2*ny+1]);
  Rinv[2*ny+2]=1./det*(R[0*ny+0]*R[1*ny+1]-R[1*ny+0]*R[0*ny+1]);
}
void transpose(double *Rtran, double *R,int nx, int ny){
  for (int i = 0; i < 3; i++){
    for ( int j = 0; j < 3; j++){
      Rtran[i*ny+j] = R[j*ny+i];}}}

void matmul(double * X,double * A, double * B,int nx, int ny){
  //solve  X = A * B
  X[0*ny+0] = A[0*ny+0]*B[0*ny+0]+A[0*ny+1]*B[1*ny+0]+A[0*ny+2]*B[2*ny+0];
  X[0*ny+1] = A[0*ny+0]*B[0*ny+1]+A[0*ny+1]*B[1*ny+1]+A[0*ny+2]*B[2*ny+1];
  X[0*ny+2] = A[0*ny+0]*B[0*ny+2]+A[0*ny+1]*B[1*ny+2]+A[0*ny+2]*B[2*ny+2];

  X[1*ny+0] = A[1*ny+0]*B[0*ny+0]+A[1*ny+1]*B[1*ny+0]+A[1*ny+2]*B[2*ny+0];
  X[1*ny+1] = A[1*ny+0]*B[0*ny+1]+A[1*ny+1]*B[1*ny+1]+A[1*ny+2]*B[2*ny+1];
  X[1*ny+2] = A[1*ny+0]*B[0*ny+2]+A[1*ny+1]*B[1*ny+2]+A[1*ny+2]*B[2*ny+2];

  X[2*ny+0] = A[2*ny+0]*B[0*ny+0]+A[2*ny+1]*B[1*ny+0]+A[2*ny+2]*B[2*ny+0];
  X[2*ny+1] = A[2*ny+0]*B[0*ny+1]+A[2*ny+1]*B[1*ny+1]+A[2*ny+2]*B[2*ny+1];
  X[2*ny+2] = A[2*ny+0]*B[0*ny+2]+A[2*ny+1]*B[1*ny+2]+A[2*ny+2]*B[2*ny+2];
   
}

void matvecmul(double * X,double * A, double * B,int nx, int ny){
  //same function as above except B is a vector -- therefore so is X
  //solve  X = A * B  
  X[0] = A[0*ny+0]*B[0]+A[0*ny+1]*B[1]+A[0*ny+2]*B[2];
  X[1] = A[1*ny+0]*B[0]+A[1*ny+1]*B[1]+A[1*ny+2]*B[2];
  X[2] = A[2*ny+0]*B[0]+A[2*ny+1]*B[1]+A[2*ny+2]*B[2];
}
void wrapper(double *dRdt,double* dLydt,double *dLzdt,
	     double *R,double Ly,double Lz, double * Ie,
	     double a, double b, double c, double alpha, int nx, int ny){
  //-------------// 
  //all fo the following code will go in wrapper funtion
  //#get Moment of Inertia Tensor in the inertial frame
  //# I = R X Ie X R^T
  double * Rtran=NULL;
  Rtran = new double[nx*ny];
  transpose(Rtran, R,nx, ny);
  double * IexRt=NULL;
  IexRt = new double[nx*ny];
  matmul(IexRt,Ie, Rtran,nx, ny);
  
  double * I=NULL;
  I = new double[nx*ny];
  matmul(I,R,IexRt,nx,ny);
  
  //invert I
  double * Iinv=NULL;
  Iinv = new double[nx*ny];
  invert( Iinv, I,nx,ny);
  
  // build # vector L angular momentum
  double * L=NULL;
  L = new double[ny];
  L[0] = 0.;
  L[1] = Ly;
  L[2] = Lz;
    
  //#construct omega(n)
  //#Omega = I^-1 X L
  double * omega=NULL;
  omega = new double[ny];
  matvecmul(omega,Iinv,L,nx,ny);
  
  //build omega star
  double * omegastar=NULL;
  omegastar = new double[ny*nx];
  omegastar[0*ny+0] = 0.;omegastar[0*ny+1] = -1.*omega[2]; omegastar[0*ny+2] = omega[1];
  omegastar[1*ny+0] = omega[2];omegastar[1*ny+1] = 0.; omegastar[1*ny+2] = -1.*omega[0];
  omegastar[2*ny+0] = -1.*omega[1];omegastar[2*ny+1] = omega[0]; omegastar[2*ny+2] = 0.;
  //build Rinv
  double * Rinv=NULL;
  Rinv = new double[nx*ny];
  invert( Rinv, R,nx,ny);
  //calculate prime coords
  double Rxx = Rinv[0*ny+0];
  double Ryx = Rinv[1*ny+0];
  double Rzx = Rinv[2*ny+0];
    
    
  double norm=sqrt(a*a*Rxx*Rxx+b*b*Ryx*Ryx+c*c*Rzx*Rzx);
  double xp = -a*a*Rxx/norm;
  double yp = -b*b*Ryx/norm;
  double zp = -c*c*Rzx/norm;
  //build primecords
  double * primecords=NULL;
  primecords = new double[ny];
  double * curlycords=NULL;
  curlycords = new double[ny];
  primecords[0] = xp;
  primecords[1] = yp;
  primecords[2] = zp;
  
  double dotprod =   (Rxx*xp/a*a+Ryx*yp/b*b+Rzx*zp/c*c);
  if(dotprod<0.){
    matvecmul(curlycords,R,primecords,nx,ny);
  }
  else{
    primecords[0] = -xp;
    primecords[1] = -yp;
    primecords[2] = -zp;
    matvecmul(curlycords,R,primecords,nx,ny);}
  //read out the derivatives 
  matmul(dRdt,omegastar, R,nx, ny);
  dLydt[0] = alpha*curlycords[2];
  dLzdt[0] = -alpha*curlycords[1];

  delete Rtran;
  delete IexRt;
  delete I;
  delete Iinv;
  delete L;
  delete omega;
  delete omegastar;
  delete Rinv;
  delete primecords;
  delete curlycords;

  
}

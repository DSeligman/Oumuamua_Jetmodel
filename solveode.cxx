//Jet Torque Model with time Delay in C++ for OUMUAMUA
//v1.0 11/20/2018
//Darryl Seligman
#include "linalgop.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

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
  //  double * Rinv=NULL;
  //  Rinv = new double[nx*ny];
  double * Ie=NULL;
  Ie = new double[nx*ny];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      Ie[i*ny+j]=0.;}}
  
  R[0*ny+0] = cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi);
  R[0*ny+1] = -cos(psi)*cos(theta)*sin(phi)-cos(phi)*sin(psi);
  R[0*ny+2] = sin(phi)*sin(theta);
  R[1*ny+0] = cos(psi)*sin(phi)+cos(phi)*cos(theta)*sin(psi);
  R[1*ny+1] = cos(phi)*cos(psi)*cos(theta)-sin(phi)*sin(psi);
  R[1*ny+2] = -cos(phi)*sin(theta);
  R[2*ny+0] = sin(psi)*sin(theta);
  R[2*ny+1] = cos(psi)*sin(theta);
  R[2*ny+2] = cos(theta);
  //invert( Rinv, R,nx,ny);
  cout<<"R"<<endl;
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      cout<<"i="<<i<<endl;
      cout<<R[i*ny+j]<<endl;}}
  //construct moment of inertia tensor in inertial frame
  double M=1.;
  double Ixx = M/5.*(b*b+c*c);
  double Iyy = M/5.*(a*a+c*c);
  double Izz = M/5.*(a*a+b*b);
  Ie[0*ny+0] = Ixx;
  Ie[1*ny+1] = Iyy;
  Ie[2*ny+2] = Izz;
  cout<<"Ixx "<<Ixx<<endl;
  cout<<"Iyy "<<Iyy<<endl;
  cout<<"Izz "<<Izz<<endl;
  double Ly = 0.;
  double Lz = 0.;
  //all variables we need for RK4
  double * dRdt1=NULL;
  dRdt1 = new double[nx*ny];
  double * dLydt1=NULL;
  dLydt1 = new double[1];
  double * dLzdt1=NULL;
  dLzdt1 = new double[1];
  
  double * dRdt2=NULL;
  dRdt2 = new double[nx*ny];
  double * dLydt2=NULL;
  dLydt2 = new double[1];
  double * dLzdt2=NULL;
  dLzdt2 = new double[1];
  
  double * dRdt3=NULL;
  dRdt3 = new double[nx*ny];
  double * dLydt3=NULL;
  dLydt3 = new double[1];
  double * dLzdt3=NULL;
  dLzdt3 = new double[1];
  
  double * dRdt4=NULL;
  dRdt4 = new double[nx*ny];
  double * dLydt4=NULL;
  dLydt4 = new double[1];
  double * dLzdt4=NULL;
  dLzdt4 = new double[1];

  dLydt1[0]=0.;
  dLzdt1[0]=0.;

  //now we start the time integration
  int bigN=10000;//.e4;
  double endtime=.5;
  double dt =endtime*(60.*60.*24.) /(float)bigN;
  double * t=NULL;
  t = new double[bigN];
  t[0] = 0.;
  for (int i = 1; i < bigN; i++){
    t[i]=i*dt+t[i-1];}
  cout<<"dt"<<endl;
  cout<<dt<<endl;
  ofstream outputFile1("Ly.txt");
  ofstream outputFile2("Lz.txt");

  for (int stepper = 0; stepper < bigN; stepper++){
    wrapper(dRdt1, dLydt1,dLzdt1,R,Ly,Lz,Ie,a,b,c,alpha,nx,ny);
    for(int k=0;k<nx*ny;k++)
      R[k] = R[k]+dt*dRdt1[k];
    Ly=Ly+dt*dLydt1[0];
    Lz = Lz +dt*dLzdt1[0];
    outputFile1 <<Ly<<endl;
    outputFile2 <<Lz<<endl;
  }
  
  //  wrapper(dRdt, dLydt,dLzdt,R,Ly,Lz,Ie,a,b,c,alpha,nx,ny);

  /*
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
  omegastar[1*ny+0] = omega[3];omegastar[1*ny+1] = 0.; omegastar[1*ny+2] = -1.*omega[0];
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
  double * dRdt=NULL;
  dRdt = new double[nx*ny];
  matmul(dRdt,omegastar, R,nx, ny);
  double dLydt = alpha*curlycords[2];
  double dLzdt = -alpha*curlycords[1];
  cout<<"dRdt"<<endl;

  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      cout<<"i="<<i<<endl;
      cout<<dRdt[i*ny+j]<<endl;}}
  cout<<"dLydt"<<endl;
  cout<<dLydt<<endl;
  cout<<"dLzdt"<<endl;
  cout<<dLzdt<<endl;*/
  return 0;
  
}

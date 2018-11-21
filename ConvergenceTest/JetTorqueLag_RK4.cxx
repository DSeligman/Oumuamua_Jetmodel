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
  //now run for longer
  endtime = 30.;
  bigN = 4000000;
  
  double dt =endtime*(60.*60.*24.) /(float)bigN;
  double * t=NULL;
  t = new double[bigN];
  t[0] = 0.;
  for (int i = 1; i < bigN; i++){
    t[i]=i*dt+t[i-1];}
  cout<<"dt"<<endl;
  cout<<dt<<endl;

  //hold derivatives
  double * dLydtar=NULL;
  dLydtar = new double[bigN];
  double * dLzdtar=NULL;
  dLzdtar = new double[bigN];
  //set up the lag
  double tlag= 6000.;
  int lag = (int) tlag/dt;
  cout<<"lag"<<lag<<endl;
  ofstream outputFile1("Ly.txt");
  ofstream outputFile2("Lz.txt");
  double * R1=NULL;
  R1 = new double[nx*ny];
  double * R2=NULL;
  R2 = new double[nx*ny];
  double * R3=NULL;
  R3 = new double[nx*ny];
  double Ly1,Lz1,Ly2,Lz2,Ly3,Lz3;
  for (int stepper = 0; stepper < bigN; stepper++){
    //first step in RK
    wrapper(dRdt1, dLydt1,dLzdt1,R,Ly,Lz,Ie,a,b,c,alpha,nx,ny);
    for(int k=0;k<nx*ny;k++)
      R1[k] = R[k]+dt/2.*dRdt1[k];
    Ly1=Ly+dt/2.*dLydt1[0];
    Lz1 = Lz +dt/2.*dLzdt1[0];
    
    wrapper(dRdt2, dLydt2,dLzdt2,R1,Ly1,Lz1,Ie,a,b,c,alpha,nx,ny);
    for(int k=0;k<nx*ny;k++)
      R2[k] = R[k]+dt/2.*dRdt2[k];
    Ly2=Ly+dt/2.*dLydt2[0];
    Lz2 = Lz +dt/2.*dLzdt2[0];
    
    wrapper(dRdt3, dLydt3,dLzdt3,R2,Ly2,Lz2,Ie,a,b,c,alpha,nx,ny);
    for(int k=0;k<nx*ny;k++)
      R3[k] = R[k]+dt*dRdt3[k];
    Ly3=Ly+dt*dLydt3[0];
    Lz3 = Lz +dt*dLzdt3[0];
    
    wrapper(dRdt4, dLydt4,dLzdt4,R3,Ly3,Lz3,Ie,a,b,c,alpha,nx,ny);
    
    //all steps done update original arrays
    for(int k=0;k<nx*ny;k++)
      R[k] = R[k]+dt*(1./6.*dRdt1[k] + 1./3.*dRdt2[k] +1./3.*dRdt3[k] +1./6.*dRdt4[k]);
    //readout the derivative approx
    dLydtar[stepper] = (1./6.*dLydt1[0] + 1./3.*dLydt2[0] +1./3.*dLydt3[0] +1./6.*dLydt4[0]);
    dLzdtar[stepper] = (1./6.*dLzdt1[0] + 1./3.*dLzdt2[0] +1./3.*dLzdt3[0] +1./6.*dLzdt4[0]);
    // this is for no lag
    /*
    Ly = Ly+dt*(1./6.*dLydt1[0] + 1./3.*dLydt2[0] +1./3.*dLydt3[0] +1./6.*dLydt4[0]);
    Lz = Lz+dt*(1./6.*dLzdt1[0] + 1./3.*dLzdt2[0] +1./3.*dLzdt3[0] +1./6.*dLzdt4[0]);
    */
    //now apply the torque by changing angular momentum at some time in the past
    if(stepper>lag){
      Ly = Ly+dt*dLydtar[stepper-lag];
      Lz = Lz +dt*dLzdtar[stepper-lag];}
    else{
      Ly=Ly+0.;
      Lz = Lz +0.;}
    
    outputFile1 <<Ly<<endl;
    outputFile2 <<Lz<<endl;
  }
  

  return 0;
  
}

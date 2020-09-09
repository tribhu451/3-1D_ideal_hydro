#include <iostream>
#include <TF1.h>
#include "TMath.h"
#include <TRandom.h>
#include <TRandom3.h>
#include <TUUID.h>

using namespace std;

class Nucleus
{
  
 public:
  Nucleus(int _A,int _B,double aRadius, double adlt,double aBETA2,double aBETA4);
  ~Nucleus();
  void generate_nucleons_A_position(double* ,double* ,double* );
  void generate_nucleons_B_position(double* ,double* ,double* );
  inline void Set_Impact_Parameter(double s){b=s;}

  // [Info] nucleus orientation angles
  double etaA;
  double psiA;
  double etaB;
  double psiB;

  // [Info]  returns nucleus orientation angles
  inline double Get_ThetaA(){return etaA;}  
  inline double Get_PhiA(){return psiA;}
  inline double Get_ThetaB(){return etaB;}
  inline double Get_PhiB(){return psiB;}
  

  
 private:
  TRandom3* tr1;
  TRandom3* tr2;
  TRandom3* tr3;
  TRandom3* tr4;
  TRandom3* tr5;
  TRandom3* tr6;
  TRandom3* tr7;
  TRandom3* tr8;
  TRandom3* tr9;    
  TRandom3* trA;
  TRandom3* trB;    
  TRandom3* rsduse;
  TF1* angl;
  
  double X[500];double Y[500];double Z[500];
  int A;             //Nucleus-1
  int B;             //Nucleus-2
  double b;          //Impact Parameter

  // [Info] Woods Saxon parameters 
  double R;       // (radius)       
  double dlt;     // (skin depth)
  double ZETA;
  double BETA2;   // (Nuclear deformation parameter)
  double BETA4;    
  
}; 

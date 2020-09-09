#include <iostream>
#include <TMath.h>

using namespace std;

class nPart_nColl
{
  
public:   
  nPart_nColl(int aa,int ab, double ac,double _npp, double _xhard);
  ~nPart_nColl();  
  void Calculate_nPart_nColl(double* ,double* ,double* ,double*);


  inline double Get_Npart(){return Npart;}
  inline double Get_Ncoll(){return Ncoll;}
  inline void Npart_source(double* _x, double* _y){for(int i=0; i<Npart; i++){_x[i]=Npart_x[i]; _y[i]=Npart_y[i];}}
  inline void Ncoll_source(double* _x, double* _y){for(int i=0; i<Ncoll; i++){_x[i]=Ncoll_x[i]; _y[i]=Ncoll_y[i];}}

  
private:
  double occA[1000];double occB[1000];         //flag during calc of Npart
  double Ncoll_x[2000]; double Ncoll_y[2000];  // x & y co-ordinate of binary collision sources
  double Npart_x[1000]; double Npart_y[1000];  // x & y co-ordinate of participant sources

  int A;             //Nucleus-1
  int B;             //Nucleus-2
  double sigma;      //Energy
  int Ncoll;         // number of binary collision sources
  int Npart;         // number of participants
  double X_hard,npp; // parameters for 2-component glauber
  
};

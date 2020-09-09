#include<fstream>
#include<string>
#include<sstream>
#include "global.h"
#include<iostream>
#include<fstream>
#include "TMath.h"
	

class Fluid;
class EoS;
class nPart_nColl;
class Nucleus;

using namespace std;

class mc_glau
     {

       public:
        mc_glau();
        ~mc_glau();
        void set_initial_gubser(Fluid* f, EoS* eos,double t0, int mode);
        void Set_ic(Fluid* f, EoS* eos, double tau_0);




  
     private:
 


  // [Info] Nuclues Mass number
  int A;
  int B;

  // [Info] collision energy
  double sigma=4.2;

  // Wood Saxon Parameter          
  double Radius=6.37;       
  double dlt=0.54;        
  double BETA2=0.00;
  double BETA4=0.00;

  // [Info] two component Glauber model
  double npp=2.49;
  double X_hard=0.13;

  // [Info] Gaussian smearing 
  double DELTA =0.6;

  // [Info] energy density scaling
  double eps0 ;

  // [Info] Impact parameter range
  double bmin,bmax;

   double theta(double _a);


void set_the_parameters()
{

  string a_; double a;
  
  istringstream* iss;
  char   buff[200];
  
  fstream File0;
  File0.open("input.dat",ios::in);
  if(!File0){cout<<"No input file, exit..."; exit(1);}
  
  int number = 0;
  while(!File0.eof())
    {

      File0.getline(buff,200);
      if (!(*buff) || (*buff == '#')) {number ++; continue;}
      iss = new istringstream(buff);
      *iss >> a_ >> a ;
      if(a_ == "A" )    {A = a;}
      if(a_ == "B" )    {B = a;}
      if(a_ == "sigma" )  {sigma = a;}
      if(a_ == "Radius" )  {Radius = a;}
      if(a_ == "dlt" )  {dlt = a;}
      if(a_ == "BETA2" )  {BETA2 = a;}
      if(a_ == "BETA4" )  {BETA4 = a;}
      if(a_ == "npp" )  {npp = a;}
      if(a_ == "X_hard" )  {X_hard = a;}
      if(a_ == "DELTA" )  {DELTA = a;}       
      if(a_ == "mc_eps0" )  {eps0 = a;}       
      if(a_ == "bmin" )  {bmin = a;}       
      if(a_ == "bmax" )  {bmax = a;}       
      delete iss;
      number++;
    } 

File0.close();
}




     };

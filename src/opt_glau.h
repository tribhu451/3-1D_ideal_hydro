#include <iostream>
#include<TMath.h>
#include<fstream>
#include<string>
#include<sstream>

class Fluid;
class EoS;
class nPart_nColl;
class Nucleus;


using namespace std;

class opt_glau
{

  public:
        opt_glau();
        ~opt_glau();
        double NPart();
        void Set_ic(Fluid* f, EoS* eos, double tau_0);
        

  private:
      double beta2=0.0;
      double beta4=0.0;
      double a = 0.55;
      double R = 6.37;
      double sigma=4.2;
      double A=197.0;
      double B=197.0;
      double n_pp=2.49;
      double X_hard=0.13;

      double mThetaA =0;
      double mPhiA = 0;
      double mThetaB =0;
      double mPhiB = 0;
      double eps0 ;
      double pi = 3.1415927;
      double b;
      double norm_const;

        double Modf_WoodSaxon(const double *x);
        double Norm_Function(double* x,double* p);
        double TA_TB(double* x,double* p);
        double Get_Norm_Constant();
        double npartxy(double x,double y);
        double ncollxy(double x,double y);
        double theta(double _a);


void set_the_parameters()
{

  string a_; double ax;
  
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
      *iss >> a_ >> ax ;
      if(a_ == "A" )    {A = ax;}
      if(a_ == "B" )    {B = ax;}
      if(a_ == "sigma" )  {sigma = ax;}
      if(a_ == "Radius" )  {R = ax;}
      if(a_ == "dlt" )  {a = ax;}
      if(a_ == "BETA2" )  {beta2 = ax;}
      if(a_ == "BETA4" )  {beta4 = ax;}
      if(a_ == "npp" )  {n_pp = ax;}
      if(a_ == "X_hard" )  {X_hard = ax;}
      if(a_ == "opt_eps0" )  {eps0 = ax;}       
      if(a_ == "Impact_parameter" )  {b = ax;}       
      delete iss;
      number++;
    } 

File0.close();
}

 
};

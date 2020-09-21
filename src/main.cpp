#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<chrono>
#include "Fluid.h"
#include "global.h"
#include "eos.h"
#include "Hydro.h"
#include "freeze.h"
#include "opt_glau.h"
#include "mc_glau.h"
#include "out.h"

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::istringstream; 
using std::cin;
using std::fstream;
using std::ios;





// This functions reads the input.dat file and sets the input parameters in the code 
void input(double &xmin, double &xmax, int &nx, double &ymin, double &ymax, int &ny,
            double &etamin, double &etamax, int &neta, double &dtau, double &tau0, double &tauMax, double &Efreeze, int &init_mode)
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
      if(a_ == "nx" )    {nx = a;}
      if(a_ == "ny" )    {ny = a;}
      if(a_ == "neta" )  {neta = a;}
      if(a_ == "xmin" )  {xmin = a;}
      if(a_ == "xmax" )  {xmax = a;}
      if(a_ == "ymin" )  {ymin = a;}
      if(a_ == "ymax" )  {ymax = a;}
      if(a_ == "etamin" )  {etamin = a;}
      if(a_ == "etamax" )  {etamax = a;}
      if(a_ == "dtau" )  {dtau = a;}
      if(a_ == "tau0" )  {tau0 = a;}
      if(a_ == "Efreeze" )  {Efreeze = a;}
      if(a_ == "tauMax" )  {tauMax = a;}
      if(a_ == "init_mode" )  {init_mode = a;}
       
      delete iss;
      number++;
    } 

  File0.close();
}






// MAIN function
int main()
{
  
  // Get starting timepoint 
  auto start = high_resolution_clock::now();  
  ofstream file1;
  
  // [Info] input parameters to be set from "input.dat" file
  // [Info] although some default values are given, it will be replaced by input parameters
  //        according to file "input.dat"
  double xmax=14.,ymax=14.,etamax=10.;
  double xmin=-14.,ymin=-14.,etamin=-10.;
  int nx=141,ny=141,neta=101;
  double tau0=0.25,dtau=0.05,tauMax=13.0;
  double Efreeze=0.24;
  int init_mode = 0;

  //setting input paramters ...
  input(xmin,xmax,nx,ymin,ymax,ny,etamin,etamax,neta,dtau,tau0,tauMax, Efreeze,init_mode);



  cout<<"*********************************************"<<endl;
  cout<<"*                                           *"<<endl;
  cout<<"*              3+1D Ideal Hydro             *"<<endl;
  cout<<"*                                           *"<<endl;
  cout<<"*********************************************\n"<<endl;
  

  // hydro-output during evolution
  output_hydro ou;

  //EoS allocation
  EoS* eos = new EoS();
  cout<<"EoS allocation done ...\n"<<endl;



  //Fluid allocation
  Fluid* f = new Fluid(eos,xmin,xmax,nx,ymin,ymax,ny,etamin,etamax,neta,dtau);
  cout<<"maximum X : "<<f->getmaxX()<<"\t   minimum X : "<<f->getminX()<<"\t   nx : "<<f->Get_nx()<<"\t   dx : "<<f->Get_dx()<<endl;
  cout<<"maximum Y : "<<f->getmaxY()<<"\t   minimum Y : "<<f->getminY()<<"\t   ny : "<<f->Get_ny()<<"\t   dy : "<<f->Get_dy()<<endl;
  cout<<"maximum Z : "<<f->getmaxZ()<<"\t   minimum Z : "<<f->getminZ()<<"\t   neta : "<<f->Get_nz()<<"\t   dz : "<<f->Get_dz()<<endl;
  cout<<"Fluid allocation done ...\n"<<endl;
  cout<<"dt/dx : "<<dtau/f->Get_dx()<<"\tdt/dy : "<<dtau/f->Get_dy()<<"\tdt/dz : "<<dtau/f->Get_dz()<<"\n"<<endl;



  // setting the initial condition ...
  if (init_mode == 0)
     {
      opt_glau* In = new opt_glau();
      In->Set_ic(f,eos, tau0);
     }
  if (init_mode == 1)
     {
      mc_glau* In = new mc_glau();
      In->Set_ic(f,eos, tau0);
     } 
  cout<<"Initialisation is done ... \n"<<endl;
  


  // freezeout hypersurface initialization ...
  evolve* map =new evolve();
  cout<<"\n\n"<<endl;  
  cout<<"          MAP INIT        "<<endl;
  cout<<"  ************************"<<endl;
  map->ini(f,tau0,tauMax,dtau);
  cout<<"  *************************"<<endl;
  cout<<"\n\n"<<endl;  
  
 
  
  //Hydro set up ...
  Hydro* h = new Hydro(eos, f, tau0, dtau);
  cout<<"Hydro started..."<<endl;
  

  // at starting time hypersurface is stored. energy density is also stored at tau = tau0.
  map->put(0,f,tau0,eos);
  ou.output(f,eos,tau0); ou.anisotropy_output(f,eos,tau0);

  // evolution loop ...
  int maxstep = ceil((tauMax - tau0) / dtau);
  for (int istep = 0; istep < maxstep; istep++)
    {
      // decrease timestep automatically, but use fixed dtau for output
      int nSubSteps = 1;
      while (dtau / nSubSteps > 1.0 * (tau0 + dtau * istep) * (etamax - etamin) / (neta - 1))
	   {
            nSubSteps *= 2; 
            cout<<"      small tau      "<<endl;
           h->setDtau(dtau / nSubSteps);
           //cout<<"dtau = "<<dtau / nSubSteps<<endl;
	   }
      for (int j = 0; j < nSubSteps; j++)
	{
	  h->evolve();
	}

      int gg=(istep+1)/map->stepsave;
      if ((gg*map->stepsave)==(istep+1))
	{
	  map->put(gg,f,h->getTau(),eos);
	  cout<<" saving for time="<<h->getTau()<<" evolution step="<<istep+1<<
	    ", "<<gg+1<<" step saved"<<endl;
	}


      //uncomment below line to store the energy density value at tau = 1.00fm/c
      if(abs(h->getTau()-1.00)<0.001){ou.output(f,eos,h->getTau());ou.anisotropy_output(f,eos,h->getTau());}      

    }

  


  // Getting the freezeout hypersurface in .xml file .
  cout<<"\n\n freezeout temperature : "<<eos->temperature(Efreeze,0,0,0)*1000<<" MeV"<<endl;
  map->gethypersurface(f,eos,Efreeze);
  map->hypersurface("hydro_output/",0);
  


  
 // Get ending timepoint 
  auto stop = high_resolution_clock::now();
  auto sdxduration = duration_cast<seconds>(stop - start);
  cout << "Total time taken : "
       << sdxduration.count() << "seconds" << endl;
  
  ofstream aFile;
  aFile.open("hydro_complete_flag.txt");
  aFile<<"Job completed ..."<<endl;
  aFile.close();

  delete f;
  delete eos;
  delete h;
  //delete In;
  delete map;
  
  return 0;
}


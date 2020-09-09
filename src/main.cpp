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

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::istringstream; 
using std::cin;
using std::fstream;
using std::ios;


//This function stores the energy density value at each grid point of the entire fluid at a given time. 
void output(Fluid* f,EoS* eos,double time)
{
  
  cout<<"After time "<<time<<" data is recorded"<<endl;
  char name[200];

  ofstream File1;
  sprintf(name,"output_after/longitudinal_at_tau_%f_fm.dat",time);	
  File1.open(name);
  
  ofstream File2;
  sprintf(name,"output_after/transverse_at_tau_%f_fm.dat",time);	
  File2.open(name);

  ofstream File3;
  sprintf(name,"output_after/temperature_at_tau_%f_fm.dat",time);	
  File3.open(name);
  
  double erg,vx,vy,vz,nb,nq,ns,pressure;
  int nx  = f->Get_nx();
  int ny  = f->Get_ny();
  int nz  = f->Get_nz();
  for(int i=0; i<nx ; i++)
    {
      for(int j=0; j<ny; j++)
        {
          for(int k=0; k<nz; k++)
	    {
	      double x_ = f->Get_x(i);
	      double y_ = f->Get_y(j);
	      double z_ = f->Get_z(k);
	      double rt = TMath::Sqrt(x_*x_+y_*y_);
	      f->Get_cell(i,j,k)->Get_physical_var(eos, time, erg, pressure, nb, nq, ns, vx,vy, vz);          
	      double T = eos->temperature(erg,nb, nq,ns);
	      if (i==nx/2 && j==ny/2)                    File1<<z_<<"\t"<<erg<<"\t"<<0.5*log((1.0+vz)/(1.0-vz))<<endl;
	      if(k == nz/2 )                             File2<<f->Get_x(i)<<"\t"<<f->Get_y(j)
                                                              <<"\t"<<rt<<"\t"<<erg<<"\t"<<nb<<"\t"<<vx<<"\t"<<vy<<endl;
	      if(i==(nx)/2 && j==(ny)/2 && k == (nz)/2)  File3<<time<<"\t"<<1000*T<<"\t"<<erg<<endl;
	    }
	}
    }
  File1.close();
  File2.close();
  File3.close();
}




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
  output(f,eos,tau0);


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

     //if( abs (h->getTau() - 1.00) < 0.001 )  output(f,eos,h->getTau());      //uncomment this line to store the energy density value at tau = 1.00fm/c
     //if( abs (h->getTau() - 5.00) < 0.001 )  output(f,eos,h->getTau());      // same as above but at tau = 5.00fm/c 
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
  


  delete f;
  delete eos;
  delete h;
  //delete In;
  delete map;
  
  return 0;
}


#include<fstream>
#include<string>
#include<sstream>
#include "global.h"
#include "mc_glau.h"
#include "Fluid.h"
#include "eos.h"
#include "Nucleus.h"
#include "nPart_nColl.h"

using namespace std;

mc_glau::mc_glau()
    {
  set_the_parameters();
    }

mc_glau::~mc_glau(){}


void mc_glau::Set_ic(Fluid* f, EoS* eos, double tau_0)
{



  cout<<"      Monte Carlo Glauber Model    "<<endl;
  cout<<"      *************************    "<<endl;
  cout<<"Mass number of Nucleus-A : "<<A<<endl;
  cout<<"Mass number of Nucleus-B : "<<B<<endl;
  cout<<"deformation parameter (Beta-2) : "<<BETA2<<" and (Beta-4) : "<<BETA4<<endl;
  if(X_hard != 1.0){cout<<"two component energy deposition with X_hard"<<endl;}
  cout<<"Energy : "<<sigma<<endl;
  cout<<"Gaussian smearing with std deviation : "<<DELTA<<endl;


  
  TRandom* t1=new TRandom();
  t1->SetSeed(0);
  long kss=t1->GetSeed();
  gRandom->SetSeed(kss);
  TF1* f1= new TF1("f1","x",0.0,25.0);
 

 
  
  Nucleus* NN = new Nucleus(A,B,Radius,dlt,BETA2,BETA4);
  nPart_nColl *PC = new nPart_nColl(A,B, sigma,npp,X_hard);
  Cell* c;


  
  double NPart,NColl;                     // no. of participant and binary collisions
  double b;                               // Impact parameter
  


  double XA[A];double YA[A];double ZA[A];
  double XB[B];double YB[B];double ZB[B];
  double npart_x[500],npart_y[500];
  double ncoll_x[2000],ncoll_y[2000];
  for(int j=0;j<=A;j++){XA[j]=0.0;YA[j]=0.0;ZA[j]=0.0;}
  for(int j=0;j<=B;j++){XB[j]=0.0;YB[j]=0.0;ZB[j]=0.0;}
  for(int j=0;j<500;j++){npart_x[j]=0.0;npart_y[j]=0.0;}
  for(int j=0;j<2000;j++){ncoll_x[j]=0.0;ncoll_y[j]=0.0;}
  


  b=f1->GetRandom(bmin,bmax);                      //Random impact parameter betwen bmin and bmax
  cout<<"b = "<<b<<" (fm)"<<endl;



  NN->Set_Impact_Parameter(b);                    //Impact parameter set to generate Nucleons
  NN->generate_nucleons_A_position(XA,YA,ZA);     //Nucleons of A generated
  //double OTA=NN->Get_ThetaA();                  //Orientation Of Nucleus A (Theta)
  //double OPA=NN->Get_PhiA();                    //---do-----A (Phi)
  NN->generate_nucleons_B_position(XB,YB,ZB);     //Nucleons of B generated
  //double OTB=NN->Get_ThetaB();                  //Orientation Of Nucleus B (Theta)
  //double OPB=NN->Get_PhiB();                    //---do-----B (Phi)
  PC->Calculate_nPart_nColl(XA,YA,XB,YB);         //N_part and N-coll calculated
  NColl=PC->Get_Ncoll();                          // Got the Ncoll
  NPart= PC->Get_Npart();                         //Got the Npart
 

 
  cout<<"No. of participants : "<<NPart<<endl;
  cout<<"No. of binary collisions : "<<NColl<<endl;


  PC->Npart_source( npart_x,  npart_y);
  PC->Ncoll_source( ncoll_x,  ncoll_y);
 

 
  double dx= f->Get_dx();
  double dy= f->Get_dy();


   for(int i=0; i<f->Get_nx(); i++)
      for(int j=0; j<f->Get_ny(); j++)
          for(int k=0; k<f->Get_nz(); k++)
           {

          c = f->Get_cell(i,j,k);         // initialising the 0th z

	  double x_ = f->Get_x(i);
	  double y_ = f->Get_y(j);
          double eta = f->Get_z(k);

          double eps=0;
          for(int ks=0; ks<NPart; ks++)
              {
                double temp1=TMath::Power(x_-npart_x[ks],2)+TMath::Power(y_-npart_y[ks],2);
                double value1=0.5*npp*(1-X_hard);
                eps=eps+((dx*dy)*(value1/(TMath::Sqrt(2*TMath::Pi()*DELTA*DELTA)))*(TMath::Exp((-1.0/(2.0*DELTA*DELTA))*(temp1))));
              }

         for(int ks=0; ks<NColl; ks++)
              {
                double temp1=TMath::Power(x_-ncoll_x[ks],2)+TMath::Power(y_-ncoll_y[ks],2);
                double value1=npp*X_hard;
                eps=eps+((dx*dy)*(value1/(TMath::Sqrt(2*TMath::Pi()*DELTA*DELTA)))*(TMath::Exp((-1.0/(2.0*DELTA*DELTA))*(temp1))));
              }

           if(eps < 0.00001 ) {eps =0;}
           eps = eps*exp(((-(( abs(eta) - 1.5)*( abs(eta)-1.5))/(2*1.3*1.3)))*theta(abs(eta)-1.5));  //https://arxiv.org/pdf/0902.4121.pdf  (eqn_2.12)
           double vx=0; double vy=0; double vz= 0.0; double nb= 0; double ns=0; double nq=0;
	   c->Set_prim_var(eos,tau_0,eps0*eps, nb, nq,  ns,  vx,  vy,  vz);

          }
    
     cout<<"\n";

  delete NN;
  delete PC;
}



double mc_glau::theta(double _a)
{
if(_a > 0){return 1;}else{return 0;}
}

/*
void mc_glau::set_initial_gubser(Fluid* f, EoS* eos, double tau, int mode) // setting initial condition
{
  
  double eps_0 =1.0;
  double n_0 = 0.5;
  double q = 1.0 ;
  Cell* c;
  
  cout<<"After time "<<tau<<" analytical data is recorded"<<endl;
  ofstream File1;
  char name[200];
  sprintf(name,"analytical_dist_%f.dat",tau);	
  File1.open(name);
  
  
  for(int i=0; i<f->Get_nx(); i++)
    {
      for(int j=0; j<f->Get_ny(); j++)
        {
          c = f->Get_cell(i,j,ZINIT);            //Initialising 0th z
	  double x_ = f->Get_x(i);
	  double y_ = f->Get_y(j);
          double rt = TMath::Sqrt(x_*x_+y_*y_);
          double eps = ((eps_0*TMath::Power(2.0*q,8.0/3.0))/(TMath::Power(tau,4.0/3.0)))*(TMath::Power((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*TMath::Power(tau*tau-rt*rt,2)),-4.0/3.0));
          double nd = ((n_0*4*q*q)/tau)*(TMath::Power((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*TMath::Power(tau*tau-rt*rt,2)),-1.0));
	  
          double k_ = TMath::ATanH((2.0*q*q*tau*rt)/(1+q*q*tau*tau+q*q*rt*rt));
          double vx =  (x_/rt)*(TMath::TanH(k_)); 
          double vy =  (y_/rt)*(TMath::TanH(k_)); 
	  
	  if(i == j){File1<<f->Get_x(i)<<"\t"<<f->Get_y(j)<<"\t"<<rt<<"\t"<<eps<<"\t"<<nd<<"\t"<<vx<<"\t"<<vy<<endl;}
	  
	  if (mode==0) {c->Set_prim_var(eos,eps,nd,vx,vy);}
        }
    }
  
  File1.close();
}
*/





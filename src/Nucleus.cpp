#include "Nucleus.h"

using std::cout;
using std::endl;

Nucleus::Nucleus(int _A,int _B,double aRadius, double adlt,double aBETA2,double aBETA4)
{
  A=_A; B=_B; R=aRadius; dlt=adlt; BETA2=aBETA2; BETA4=aBETA4;
  
  tr1 = new TRandom3();
  tr1->SetSeed(0);
  tr2 = new TRandom3();
  tr2->SetSeed(0);
  tr3 = new TRandom3();
  tr3->SetSeed(0);
  tr4 = new TRandom3();
  tr4->SetSeed(0);
  tr5 = new TRandom3();
  tr5->SetSeed(0);
  tr6 = new TRandom3();
  tr6->SetSeed(0);
  tr7 = new TRandom3();
  tr7->SetSeed(0);
  tr8 = new TRandom3();
  tr8->SetSeed(0);
  
  trA = new TRandom3();
  trA->SetSeed(0);
  trB = new TRandom3();
  trB->SetSeed(0);
  
  rsduse=new TRandom3();
  rsduse->SetSeed(0);
  long kft=rsduse->GetSeed();
  gRandom->SetSeed(kft);
  angl= new TF1("angl","TMath::Sin(x)",0.0,TMath::Pi());
  
  tr9 = new TRandom3();
  tr9->SetSeed(0);
}


Nucleus::~Nucleus()
{
  delete trA;
  delete trB;
  delete tr1;
  delete tr2;
  delete tr3;
  delete tr4;
  delete rsduse;
}



void Nucleus::generate_nucleons_A_position(double* X1, double* Y1,double* Z1)
{    
  double CMx=0.0;double CMy=0.0;double CMz=0.0;
  int count=0;
  
  do
    {
      double r=(11.0)*(tr1->Rndm());
      double Theta=(TMath::Pi())*(tr2->Rndm());
      double Phi=((2.0)*TMath::Pi())*(tr3->Rndm());
      double test=tr4->Rndm();
      
      
      double Y20=0.25*TMath::Sqrt(5.0/TMath::Pi())*
	(3*TMath::Cos(Theta)*TMath::Cos(Theta)-1.0);
      double Y40=(3.0/(16.0*TMath::Sqrt(TMath::Pi())))* 
	((35*TMath::Power(TMath::Cos(Theta),4))-
	 (30*TMath::Power(TMath::Cos(Theta),2))+3);
      double RAT= R*(1+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(TMath::Sin(Theta)))/(1.0+(TMath::Exp((r-RAT)/dlt)));
      
      
      if(test < rho )
	{      
	  
	  X[count]= (r*TMath::Sin(Theta)*TMath::Cos(Phi));
	  Y[count]=(r*TMath::Sin(Theta)*TMath::Sin(Phi));
	  Z[count]=(r*TMath::Cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
	}   
    }   
  while(count<A);
  
  CMx=CMx/A;CMy=CMy/A;CMz=CMz/A;
  
  for(int j=0;j<A;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  etaA=0;//angl->GetRandom(0.0,TMath::Pi());  // nucleus orientaton angle (theta)
  psiA=0;//((2.0)*TMath::Pi())*(trA->Rndm()); // nucleus orientation angle (phi)
  
  for(int j=0;j<A;j++)
    {
      X1[j]=(TMath::Cos(psiA)*TMath::Cos(etaA)*X[j])+(-TMath::Sin(psiA)*Y[j])+(-TMath::Cos(psiA)*TMath::Sin(etaA)*Z[j]);
      Y1[j]=(TMath::Sin(psiA)*TMath::Cos(etaA)*X[j])+(TMath::Cos(psiA)*Y[j])+(-TMath::Sin(psiA)*TMath::Sin(etaA)*Z[j]);
      Z1[j]=(TMath::Sin(etaA)*X[j])+(TMath::Cos(etaA)*Z[j]);
    }
  
}




void Nucleus::generate_nucleons_B_position(double* X1, double* Y1,double* Z1)
{
    
  double CMx=0;
  double CMy=0;
  double CMz=0;
  
  int count=0;
  do
    {
      
      double r=(11.0)*(tr1->Rndm());
      double Theta=(TMath::Pi())*(tr2->Rndm());
      double Phi=((2.0)*TMath::Pi())*(tr3->Rndm());
      double test=tr4->Rndm();
     
      double Y20=0.25*TMath::Sqrt(5.0/TMath::Pi())*
	(3*TMath::Cos(Theta)*TMath::Cos(Theta)-1.0);
      double Y40=(3.0/(16.0*TMath::Sqrt(TMath::Pi())))* 
	((35*TMath::Power(TMath::Cos(Theta),4))-
	 (30*TMath::Power(TMath::Cos(Theta),2))+3);
      double RAT= R*(1.0+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(TMath::Sin(Theta)))/(1.0+(TMath::Exp((r-RAT)/dlt)));
       
      if(test < rho)
	{      
	  X[count]= (r*TMath::Sin(Theta)*TMath::Cos(Phi));
	  Y[count]=(r*TMath::Sin(Theta)*TMath::Sin(Phi));
	  Z[count]=(r*TMath::Cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
	}  
    }   
  while(count<B);
  
  CMx=CMx/B;CMy=CMy/B;CMz=CMz/B;
  for(int j=0;j<B;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  
  etaB=0;//angl->GetRandom(0.0,TMath::Pi());   // nucleus orientation angle (theta)
  psiB=0;//((2.0)*TMath::Pi())*(trB->Rndm());  //  nucleus orientation angle (phi)
  
  
  for(int j=0;j<B;j++){
    X1[j]=(TMath::Cos(psiB)*TMath::Cos(etaB)*X[j])+(-TMath::Sin(psiB)*Y[j])+(-TMath::Cos(psiB)*TMath::Sin(etaB)*Z[j]);
    Y1[j]=(TMath::Sin(psiB)*TMath::Cos(etaB)*X[j])+(TMath::Cos(psiB)*Y[j])+(-TMath::Sin(psiB)*TMath::Sin(etaB)*Z[j]);
    Z1[j]=(TMath::Sin(etaB)*X[j])+(TMath::Cos(etaB)*Z[j]);
  }
  
  
  // [Info]  Here the intention is to shift the B-Nucleus randomly in xy plane .But impact paramter will remain same.
  // [Info]  This thing is only done for B-Nucleus. It's not been done for Nucleus-A.  
  double zhi=0;//(2*(TMath::Pi()))*tr9->Rndm();
  for(int j=0;j<B;j++){ X1[j]=X1[j]+((b)*TMath::Cos(zhi)); Y1[j]=Y1[j]+((b)*TMath::Sin(zhi));}
  
}


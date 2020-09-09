#include "nPart_nColl.h"

using namespace std;

nPart_nColl::nPart_nColl(int aa,int ab, double ac, double _npp, double _xhard)
{
  A=aa;B=ab;sigma=ac;npp=_npp; X_hard=_xhard;
}

nPart_nColl::~nPart_nColl(){}

void nPart_nColl::Calculate_nPart_nColl(double* vxA,double* vyA,double* vxB,double* vyB)
{
  
  Ncoll=0;
  Npart=0;
  
  for(int i=0;i<A;i++){occA[i]=0;}
  for(int i=0;i<B;i++){occB[i]=0;}
  
  for (int i=0; i<A; i++)
    {
      for (int j=0; j<B; j++)
	{  
	  double d=TMath::Sqrt( TMath::Power((vxB[j]-vxA[i]),2) + 
				TMath::Power ( (vyB[j]-vyA[i]),2));
	  double D=TMath::Sqrt(sigma/ (TMath::Pi())); 
	  
	  if( d <= D)
	    { 
	      Ncoll_x[Ncoll]=(vxA[i]+vxB[j])/2;
	      Ncoll_y[Ncoll]=(vyA[i]+vyB[j])/2;
	      Ncoll=Ncoll+1;
	      
	      if(occA[i]==0)
		{ 
		  occA[i]=1;Npart_x[Npart]=vxA[i]; Npart_y[Npart]=vyA[i];Npart=Npart+1;
		} 
	      
              
	      if(occB[j]==0)
		{
		  occB[j]=1;Npart_x[Npart]=vxB[j]; Npart_y[Npart]=vyB[j];Npart=Npart+1;
		}
	      
	      
	    }                                                           
	}                                                          
    }                                                         
  

      // [Info]  shifting the energy distributions center to (0,0,0)   
      double xref1=0.0;
      double yref1=0.0;
      double wref1=0.0;
      for(int k=0;k<Npart;k++)
	{ 
	  xref1=xref1+(Npart_x[k]*(0.5*npp*X_hard));
	  yref1=yref1+(Npart_y[k]*(0.5*npp*X_hard));
	  wref1=wref1+(0.5*npp*X_hard);
	}
      
      double xref2=0.0;
      double yref2=0.0;
      double wref2=0.0;
      for(int k=0;k<Ncoll;k++)
	{ 
	  xref2=xref2+(Ncoll_x[k]*(npp*(1-X_hard)));
	  yref2=yref2+(Ncoll_y[k]*(npp*(1-X_hard)));
	  wref2=wref2+(npp*(1-X_hard));
	}
      
      double xAverage=((xref1+xref2)/(wref1+wref2));
      double yAverage=((yref1+yref2)/(wref1+wref2));
      
      // cout<<xAverage<<"  "<<yAverage<<"\n";
      
      for(int k=0;k<Npart;k++)
	{ 
          Npart_x[k] = Npart_x[k]-xAverage;
          Npart_y[k] = Npart_y[k]-yAverage;
	}
     for(int k=0;k<Ncoll;k++)
	{ 
          Ncoll_x[k] = Ncoll_x[k]-xAverage;
          Ncoll_y[k] = Ncoll_y[k]-yAverage;
	}

}

 




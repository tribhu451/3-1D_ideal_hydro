#include "opt_glau.h"
#include<TF1.h>
#include<cmath>
#include<TF2.h>
#include<TMath.h>
#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "global.h"
#include "Fluid.h"
#include "eos.h"
#include "Nucleus.h"
#include "nPart_nColl.h"

using namespace std;

opt_glau::opt_glau(){set_the_parameters(); norm_const = Get_Norm_Constant(); }
opt_glau::~opt_glau(){}


double opt_glau::Modf_WoodSaxon(const double* x){

        double theta=0;
        double phi=0;

        double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
        double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
        double ZP = (- ( (TMath::Sin(theta))*x[0] ) -  ((TMath::Cos(theta))*(TMath::Sin(phi))*x[1])  +  
                                           ((TMath::Cos(theta))*(TMath::Cos(phi))*x[2]));            
        double rsq = ((x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2])) ;
        double RP = R* (        1+
                        (A1*(((3*ZP*ZP)/rsq)-1))+
                        A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3    ) 
               );


        double r = TMath::Sqrt(rsq);
        return 1/(TMath::Exp((r - RP)/a)+1);
}



double opt_glau::Get_Norm_Constant(){


        ROOT::Math::Functor wf(this, &opt_glau::Modf_WoodSaxon,3);
        double min[3] = {-3*R,-3*R,-3*R};
        double max[3] = {3.1*R,3.1*R,3.2*R};
        ROOT::Math::IntegratorMultiDim pd;
        pd.SetFunction(wf);
        double val = pd.Integral(min,max);
        return 1/val;

}




double opt_glau::Norm_Function(double* x,double* p){

        double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
        double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
        double ZP = (- ( (TMath::Sin(p[2]))*(TMath::Cos(p[3]))*p[0] ) -  ((TMath::Sin(p[3]))*p[1])  +  
                                            ((TMath::Cos(p[2]))*(TMath::Cos(p[3]))*x[0]));            
        double rsq = ((x[0]*x[0])+p[0]*p[0]+p[1]*p[1]);
        double RP = R* ( 1+
                        (A1*(((3*ZP*ZP)/rsq)-1))+
                        A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3 ) 

                            );
        double r = TMath::Sqrt(rsq);
        return norm_const/(TMath::Exp((r - RP)/a)+1);

}


double opt_glau::npartxy(double x,double y){


        TF1* f1;
        f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x+(b/2.0)); 
        f1->SetParameter(1,y); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-5);

        f1->SetParameter(0,x-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-5);
        double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
        delete f1;
        return  va;

}


double opt_glau::ncollxy(double x,double y){


        TF1* f1;
        f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x+(b/2.0)); 
        f1->SetParameter(1,y); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-05);
        f1->SetParameter(0,x-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-05);
        delete f1;
        return  A*B*sigma*TA*TB;

}


double opt_glau::TA_TB(double* x,double* p){

        TF1* f1;
        f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-5);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-5);
        double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
        delete f1;
        return  va;

}

double opt_glau::NPart(){
        TF2 *f3;
        f3=new TF2("T_AB",this,&opt_glau::TA_TB,-4*R,4*R,-4*R,4*R,0);
        double TAAB=f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-02);
        double N_Part=TAAB;
        delete f3;
        return N_Part;
}

void opt_glau::Set_ic(Fluid* f, EoS* eos, double tau_0)
{


  cout<<"      Optical Glauber Model    "<<endl;
  cout<<"      *********************    "<<endl;

  cout<<"Mass number of Nucleus-A : "<<A<<endl;
  cout<<"Mass number of Nucleus-B : "<<B<<endl;
  cout<<"deformation parameter (Beta-2) : "<<beta2<<" and (Beta-4) : "<<beta4<<endl;
  cout<<"Impact parameter : "<<b<<" fm"<<endl;
  double npart = NPart();
  cout<<"No. of participants : "<<npart<<endl;
  cout<<"Multiplicity scaling factor eps0 : "<<eps0<<endl;

   double total_deposited = 0.0;  // total deposited energy

  Cell* c;
   for(int i=0; i<f->Get_nx(); i++)
      for(int j=0; j<f->Get_ny(); j++)
        for(int k=0; k<f->Get_nz(); k++)
           {
          c = f->Get_cell(i,j,k);         

	  double x_ = f->Get_x(i);
	  double y_ = f->Get_y(j);
          double eta = f->Get_z(k);
          double eps;
          double nchxy =((1-X_hard)*n_pp*(npartxy(x_, y_)/2.0))+(X_hard*n_pp*ncollxy(x_, y_));
          eps = nchxy*exp(((-(( abs(eta) - 1.5)*( abs(eta)-1.5))/(2*1.3*1.3)))*theta(abs(eta)-1.5));  //https://arxiv.org/pdf/0902.4121.pdf  (eqn_2.12)
            
          total_deposited = total_deposited + eps ;

          double nb= 0; double nq = 0; double ns =0; double vx=0; double vy=0; double vz= 0;
          c->Set_prim_var(eos,tau_0, eps0*eps, nb, nq,  ns,  vx,  vy,  vz);
           }

  
     cout<<"total amount of deposited energy is = "<<total_deposited<<endl;
     cout<<"\n";

}

double opt_glau::theta(double _a)
{
if(_a > 0){return 1;}else{return 0;}
}










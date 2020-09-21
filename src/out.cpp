#include "out.h"
#include "Fluid.h"
#include "eos.h"

//This function stores the energy density value at each grid point of the entire fluid at a given time. 
void output_hydro::output(Fluid* f,EoS* eos,double time)
{
  
  cout<<"After time "<<time<<" data is recorded"<<endl;
  char name[200];

  ofstream File1;
  sprintf(name,"output_after/longitudinal_at_tau_%f_fm.dat",time);	
  File1.open(name);
  
  ofstream File2;
  sprintf(name,"output_after/transverse_at_tau_%f_fm.dat",time);	
  File2.open(name);

  //#####################################
  ofstream File30;
  sprintf(name,"output_after/temperature_at_R_0_tau_%f_fm.dat",time);	
  File30.open(name);

  ofstream File31;
  sprintf(name,"output_after/temperature_at_R_3_tau_%f_fm.dat",time);	
  File31.open(name);


  ofstream File32;
  sprintf(name,"output_after/temperature_at_R_5_tau_%f_fm.dat",time);	
  File32.open(name);
  
  ofstream File40;
  sprintf(name,"output_after/entropy_at_R_0_tau_%f_fm.dat",time);	
  File40.open(name);

  
  ofstream File41;
  sprintf(name,"output_after/entropy_at_R_3_tau_%f_fm.dat",time);	
  File41.open(name);

  
  ofstream File42;
  sprintf(name,"output_after/entropy_at_R_5_tau_%f_fm.dat",time);	
  File42.open(name);

  //#####################################


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
              double S_entropy = eos->entropy(erg,nb,nq,ns);
	      if (i==nx/2 && j==ny/2)  File1<<z_<<"\t"<<erg<<"\t"<<0.5*log((1.0+vz)/(1.0-vz))<<endl;    
	      if(k == nz/2 )                             File2<<f->Get_x(i)<<"\t"<<f->Get_y(j)
                                                              <<"\t"<<rt<<"\t"<<erg<<"\t"<<nb<<"\t"<<vx<<"\t"<<vy<<endl;

              // at (x,y,z) = (0,0,0) :
	      if(i==(nx)/2 && j==(ny)/2 && k == (nz)/2)  
                {
                  File30<<time<<"\t"<<1000*T<<"\t"<<erg<<endl;     // 1000 is multiplied with temp to convert the unit from GeV to MeV
                  File40<<time<<"\t"<<S_entropy<<"\t"<<erg<<endl;
                }

             // at (x,y,z) = ( 3,3 ,0) :
	      if( x_> 0. && y_> 0. && abs(rt-3.0)<1E-10 && k == (nz)/2  )  
                {
                  cout<<"the rt "<<rt<<endl;

                  File31<<time<<"\t"<<1000*T<<"\t"<<erg<<endl;     // 1000 is multiplied with temp to convert the unit from GeV to MeV
                  File41<<time<<"\t"<<S_entropy<<"\t"<<erg<<endl;
                }

             // at (x,y,z) = (5 ,5 ,0) :
	      if( x_> 0. && y_> 0. && abs(rt-5.0)<1E-10 && k == (nz)/2 )  
                {
                  cout<<"the rt "<<rt<<endl;
                  File32<<time<<"\t"<<1000*T<<"\t"<<erg<<endl;     // 1000 is multiplied with temp to convert the unit from GeV to MeV
                  File42<<time<<"\t"<<S_entropy<<"\t"<<erg<<endl;
                }



	    }
	}
    }
  File1.close();
  File2.close();
  File30.close();
  File31.close();
  File32.close();
  File40.close();
  File41.close();
  File42.close();

}

void output_hydro::anisotropy_output(Fluid* f,EoS* eos,double time)
{

  char name[200];
  ofstream File1;
  sprintf(name,"output_after/ex_ep_at_tau_%f_fm.dat",time);	
  File1.open(name);

  double erg,vx,vy,vz,nb,nq,ns,pressure;
  double y2mx2_avg =0.0;     // y^{2} - y^{2} average
  double y2px2_avg =0.0;     // y^{2} + x^{2} average
  double txxmtyy =0.0;       // T^{xx} - T^{yy}
  double txxptyy =0.0;       // T^{xx} + T^{yy}

  int nx  = f->Get_nx();
  int ny  = f->Get_ny();
  int nz  = f->Get_nz();
  for(int i=0; i<nx ; i++)
      for(int j=0; j<ny; j++)
          for(int k=(nz/2); k<(nz/2)+1; k++)  // only at the mid rapidity slice
	    {
	      double x_ = f->Get_x(i);
	      double y_ = f->Get_y(j);
	      double z_ = f->Get_z(k);
	      double rt = TMath::Sqrt(x_*x_+y_*y_);
	      f->Get_cell(i,j,k)->Get_physical_var(eos, time, erg, pressure, nb, nq, ns, vx,vy, vz);          
              y2mx2_avg = y2mx2_avg + ((y_*y_ - x_*x_)*erg);
              y2px2_avg = y2px2_avg + ((y_*y_ + x_*x_)*erg);

              double gamma2 = 1./ (1.-vx*vx-vy*vy-vz*vz);
              txxmtyy = txxmtyy + ((erg+pressure)*gamma2*(vx*vx-vy*vy));
              txxptyy = txxptyy + (((erg+pressure)*gamma2*(vx*vx+vy*vy))+(2*pressure));

            }
            double ex =  y2mx2_avg/ y2px2_avg;
            double ep = txxmtyy/ txxptyy;
           
            File1<<time<<"\t"<<ex<<"\t"<<ep<<endl;
            File1.close();

}






































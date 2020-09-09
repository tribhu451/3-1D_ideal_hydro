#include "global.h"
#include "Fluid.h"
#include "eos.h"


using namespace std;

Fluid::Fluid(EoS* _eos, double _xmin, double _xmax, int _nx, double _ymin, double _ymax, int _ny,double _zmin, double _zmax, int _nz, double _dt)
{
  eos = _eos;
  dt = _dt;
  nx = _nx;
  ny = _ny;
  nz = _nz;
  xmin = _xmin;
  xmax = _xmax;
  ymin = _ymin;
  ymax = _ymax;
  zmin = _zmin;
  zmax = _zmax;
  dx = (xmax - xmin) / (nx - 1);
  dy = (ymax - ymin) / (ny - 1);
  dz = (zmax - zmin) / (nz - 1);
  if(dt > dx/2. || dt > dy/2. || dt > dz/2.){cout<<"[Error]  too big delta_tau"<<endl; exit(1);}
  
  cell = new Cell[nx * ny * nz];

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++)
	{
	  Get_cell(ix, iy, iz)->Set_prev_cell(X_, Get_cell(ix - 1, iy, iz));
	  Get_cell(ix, iy, iz)->Set_next_cell(X_, Get_cell(ix + 1, iy, iz));
	  Get_cell(ix, iy, iz)->Set_prev_cell(Y_, Get_cell(ix, iy - 1, iz));
	  Get_cell(ix, iy, iz)->Set_next_cell(Y_, Get_cell(ix, iy + 1, iz));
	  Get_cell(ix, iy, iz)->Set_prev_cell(Z_, Get_cell(ix, iy, iz - 1));
	  Get_cell(ix, iy, iz)->Set_next_cell(Z_, Get_cell(ix, iy, iz + 1));
	  Get_cell(ix, iy, iz)->Set_pos(ix, iy, iz);
	}
     
}


Fluid::~Fluid(){ delete[] cell; }


void Fluid::getCMFvariablesUmunu(Cell *c, double tau, double &e, double &nb, 
				 double &nq, double &ns, double &vx, double &vy, double &Y){
  double p; double vz;
  c->Get_physical_var(eos, tau, e, p, nb, nq, ns,vx, vy, vz);
  double eta = Get_z(c->Get_iz());
  Y = eta + 0.5*log((1.+vz)/(1.-vz));
  double gf=1./sqrt(1.-vx*vx-vy*vy-vz*vz);
  vx = vx*gf;
  vy = vy*gf;
}
           //check and use



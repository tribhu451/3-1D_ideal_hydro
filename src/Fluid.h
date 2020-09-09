#include<iostream>
#include "Cell.h"
#include "global.h"


using namespace std;

class Fluid
{
  
 private:
  EoS* eos ;
  Cell *cell;            
  int nx, ny,nz;
  double xmin, xmax, ymin, ymax,zmin, zmax;
  double dx, dy,dz, dt; 
  
 public:
  Fluid(EoS* _eps, double _xmin, double _xmax, int _nx, double _ymin, double _ymax, int _ny,double _zmin, double _zmax, int _nz,double _dt);
  ~Fluid();
  double Get_nx(){return nx;}  
  double Get_ny(){return ny;}
  double Get_nz(){return nz;}
  inline double getmaxX() { return xmax; }
  inline double getmaxY() { return ymax; }
  inline double getmaxZ() { return zmax; }
  inline double getminX() { return xmin; }
  inline double getminY() { return ymin; }
  inline double getminZ() { return zmin; }
  double Get_dx() { return dx; } 
  double Get_dy(){return dy;}
  double Get_dz(){return dz;}
  double Get_dt(){return dt;}
  double Get_x(int ix) { return xmin + ix * dx; }  //returns x-coordinate of a cell
  double Get_y(int iy) {return  ymin + iy * dy;}      //returns y-coordinate of a cell
  double Get_z(int iz) {return  zmin + iz * dz;}      //returns eta-coordinate of a cell
  void getCMFvariablesUmunu(Cell *c, double tau, double &e, double &nb, 
			    double &nq, double &ns, double &vx, double &vy, double &Y);
  

 inline Cell *Get_cell(int ix, int iy, int iz) {
  ix = ix > 0 ? ix : 0;
  ix = ix < nx ? ix : nx - 1;
  iy = iy > 0 ? iy : 0;
  iy = iy < ny ? iy : ny - 1;
  iz = iz > 0 ? iz : 0;
  iz = iz < nz ? iz : nz - 1;
  return &cell[ix + nx * iy + nx * ny * iz];
 }
  
};

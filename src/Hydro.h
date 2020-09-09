#include "global.h"

class EoS;
class Fluid;
class Cell;



class Hydro
{

private:
Fluid *f;
EoS *eos;
double dt;
double tau;


public:
Hydro(EoS* _eos, Fluid* _f,double _tau0, double _dt);
~Hydro();
void evolve();
void HLLE_flux(Cell* left, Cell* right, int direction, int mode);
void sourcestep(int mode, int ix, int iy,int iz, double _tau);
inline double getTau() { return tau; }
 void setDtau(double deltaTau);  // change the timestep
 double getDtau() { return dt; }  // return current value of timestep
};


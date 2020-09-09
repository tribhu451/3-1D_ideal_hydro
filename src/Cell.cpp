#include "Cell.h"
#include "Convert.h"


using std::cout;
using std::endl;



Cell::Cell()
{
  for(int i = 0; i<7; i++)   //T^00, T^0x, T^0y, T^0z, nb, nq, ns (7 conserved quantities)
     {
	Q[i] = 0.;
	Qh[i] = 0.;
	Qprev[i] = 0.;
	flux[i] = 0.;
      }
}


double Cell::minmod(double a, double b) {
 if (a * b <= 0.) return 0.;
 //	else return (a*a*b+a*b*b)/(a*a+b*b) ;
 if (fabs(a) > fabs(b))
  return b;
 else
  return a;
}
/*
double Cell::minmod(double dl, double dr)
{ 
if((dl*dr)>0. && abs(dl)<abs(dr)){return dl;}
if((dl*dr)>0. && abs(dl)>abs(dr)){return dr;}
else{return 0.;}
}
*/

 
void Cell::Set_prim_var(EoS* eos,double tau, double _eps,double _nb,
                                 double _nq, double _ns, double _vx, double _vy, double _vz)
{
  double gamma2 = 1.0/(1.0- (_vx*_vx+_vy*_vy+_vz*_vz));
  double gamma = TMath::Sqrt(gamma2);
  double p = eos->pressure(_eps,_nb, _nq, _ns);
  Q[T_] = tau*(((_eps + p)*gamma2) - p);
  Q[X_] = tau*((_eps+p)*gamma2*_vx);
  Q[Y_] = tau*((_eps+p)*gamma2*_vy);
  Q[Z_] = tau*((_eps+p)*gamma2*_vz);
  Q[NB_] = tau*_nb*gamma;
  Q[NQ_] = tau*_nq*gamma;
  Q[NS_] = tau*_ns*gamma;
}



void Cell::Get_physical_var(EoS *eos, double tau, double &_e, double &_p,
                                      double &_nb, double &_nq, double &_ns, double &_vx,
                                      double &_vy, double &_vz)
{
  double _Q[7];
  for(int i =0; i<7 ; i++){_Q[i] = Q[i]/tau; }
  CALC_2_LRF( eos,  _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}



void Cell::Get_left_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{

 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->Get_Q(Qr);
 prev_cell[dir - 1]->Get_Q(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Q[i] - dQ[i]) / tau;
 CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);

}



void Cell::Get_right_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{
 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->Get_Q(Qr);
 prev_cell[dir - 1]->Get_Q(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Q[i] + dQ[i]) / tau;
 CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);

}


void Cell::Get_left_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{

 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->Get_Qh(Qr);
 prev_cell[dir - 1]->Get_Qh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Qh[i] - dQ[i]) / tau;
 CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);




}



void Cell::Get_right_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{
 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->Get_Qh(Qr);
 prev_cell[dir - 1]->Get_Qh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Qh[i] + dQ[i]) / tau;
 CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);

}



void Cell::addFlux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns)
{

  if(std::isinf(Ft) or std::isnan(Ft)) {
   std::cout<<"at pos -> "<<Get_ix()<<"\t"<<Get_iy()<<"\t"<<Get_iz()<<"\t   : "<<Ft<<std::endl;
   std::cout << "Cell::addFlux inf/nan\n"; exit(1);
  }
flux[T_] += Ft;
flux[X_] += Fx;
flux[Y_] += Fy;
flux[Z_] += Fz;
flux[NB_] += Fnb;
flux[NQ_] += Fnq;
flux[NS_] += Fns;
}

void Cell :: updateQtoQhByFlux()
{
   for(int i = 0; i<7; i++) Qh[i] = Q[i] + flux[i];
}

void Cell :: updateByFlux()
{
if(Q[0] + flux[0] < 0.)
	return;
for(int i = 0; i<7; i++) Q[i] += flux[i];
}


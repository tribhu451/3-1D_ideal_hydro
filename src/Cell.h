#pragma once
#include<iostream>
#include "TMath.h"
#include "global.h"

class EoS;

using namespace std;

class Cell
{
  
private:
  int ix, iy, iz;
  double Q[7];
  double Qh[7];
  double Qprev[7];
  double flux[7];

  double minmod(double dl, double dr);
  
  
  public:
  Cell();
  ~Cell(){};
  inline void Set_pos(int _ix, int _iy, int _iz) {ix = _ix;iy = _iy;iz = _iz;} //sets the cell position

  int Get_ix(){return ix;} //returns cell position number along X-axis
  int Get_iy(){return iy;} //returns cell position number along Y-axis
  int Get_iz(){return iz;} //returns cell position number along Y-axis

  void Set_prev_cell(int i, Cell* c){prev_cell[i-1]=c;}  //sets previous cell adress. prev_cell[0] means previous cell along X-axis and prev_cell[1] means previous cell along Y-axis 
  void Set_next_cell(int i, Cell* c){next_cell[i-1]=c;}  //sets next cell adress

  Cell* Get_prev_cell(int i){return prev_cell[i-1];}     //returns previous cell adress
  Cell* Get_next_cell(int i){return next_cell[i-1];}     //returns next cell adress

  void Set_prim_var(EoS* eos,double tau, double _eps,double _nb, double _nq, double _ns, double _vx, double _vy, double _vz); //set the calculational frame variables like E, Mx, My and R.
  void Get_physical_var(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz);

  inline void saveQprev(void){
  for(int i = 0; i<7; i++) Qprev[i] = Q[i];
  }

  inline void clearFlux(void){
  for(int i = 0; i<7; i++) flux[i] = 0.;
  }

  void Get_Q(double* _Q){for(int i=0; i<7; i++){_Q[i]=Q[i];}};
  void Get_Qh(double* _Qh){for(int i=0; i<7; i++){_Qh[i]=Qh[i];}};

  void Get_left_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void Get_right_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void Get_left_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void Get_right_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void addFlux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns);

  void updateQtoQhByFlux();
  void updateByFlux();
  Cell* next_cell[3];
  Cell* prev_cell[3];

};

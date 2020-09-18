#include<iostream>
#include "TMath.h"
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using std::istringstream; 
using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::ios;


using namespace std;

class EoS
{

private:
void is_file_exist(fstream& file);
  fstream infile1_d;
  fstream infile2_d;
  fstream infile3_d;
  fstream infile4_d;
  fstream infile5_d;
  fstream infile6_d;
  fstream infile7_d;

  fstream infile1_t;
  fstream infile2_t;
  fstream infile3_t;
  fstream infile4_t;
  fstream infile5_t;
  fstream infile6_t;
  fstream infile7_t;

  double e0_1,de_1; int ne_1;
  double e0_2,de_2; int ne_2;
  double e0_3,de_3; int ne_3;
  double e0_4,de_4; int ne_4;
  double e0_5,de_5; int ne_5;
  double e0_6,de_6; int ne_6;
  double e0_7,de_7; int ne_7;

  double eps_1[501]={0.0};
  double prs_1[501]={0.0};
  double ent_1[501]={0.0}; 
  double tmp_1[501]={0.0};

  double eps_2[501]={0.0};
  double prs_2[501]={0.0};
  double ent_2[501]={0.0}; 
  double tmp_2[501]={0.0};

  double eps_3[501]={0.0};
  double prs_3[501]={0.0};
  double ent_3[501]={0.0}; 
  double tmp_3[501]={0.0};

  double eps_4[501]={0.0};
  double prs_4[501]={0.0};
  double ent_4[501]={0.0}; 
  double tmp_4[501]={0.0};

  double eps_5[501]={0.0};
  double prs_5[501]={0.0};
  double ent_5[501]={0.0}; 
  double tmp_5[501]={0.0};

  double eps_6[501]={0.0};
  double prs_6[501]={0.0};
  double ent_6[501]={0.0}; 
  double tmp_6[501]={0.0};

  double eps_7[501]={0.0};
  double prs_7[501]={0.0};
  double ent_7[501]={0.0}; 
  double tmp_7[501]={0.0};


  double nb,mub,mus,muq;

  istringstream* iss;
  char   buff[200];


  gsl_interp_accel *acc;

gsl_spline *spline_1; 
gsl_spline *spline_2; 
gsl_spline *spline_3; 
gsl_spline *spline_4;
gsl_spline *spline_5; 
gsl_spline *spline_6;
gsl_spline *spline_7; 

gsl_spline *spline_prs_1;
gsl_spline *spline_prs_2;
gsl_spline *spline_prs_3;
gsl_spline *spline_prs_4;
gsl_spline *spline_prs_5;
gsl_spline *spline_prs_6;
gsl_spline *spline_prs_7;

public:
EoS();
~EoS();
double pressure( double eg,double _nb, double _nq, double _ns);
double temperature( double eg,double _nb, double _nq, double _ns);
double cs(){return 1./TMath::Sqrt(3);};
double cs2(){return 1./3.;};
double cs_2( double eg,double _nb, double _nq, double _ns);
};

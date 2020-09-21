/**********************************************
 * using the functions of this class          *
 * various quantities like  energy density,   *
 * temperature, entropy, momentum anisotropy, *
 * spatial anisotropy etc. will be saved      *
 * as a function of tau.                      *
 **********************************************/


#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<chrono>
#include<cmath>
#include "global.h"

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::istringstream; 
using std::cin;
using std::fstream;
using std::ios;

class Fluid;
class EoS;

class output_hydro
{
public :
void output(Fluid* f,EoS* eos,double time);
void anisotropy_output(Fluid* f,EoS* eos,double time);



};

/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <sstream>
#include <TMath.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "Model_BWA.h"

using namespace TMath;
using namespace std;

extern Configurator* sMainConfig;
extern TString	sModelINI;
extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;

Model_BWA::Model_BWA()
: Model(),
  mThermo(0),
  mRapPRange(0.0), mRapSRange(0.0), mRhoMax(0.0), mTau(0.0), mVt(0.0), mA(0.0), mDelay(0.0)
{
}

Model_BWA::Model_BWA(TRandom2* aRandom)
: Model(aRandom)
{
  if      (sModel == 2) mName = "BWAVT";
  else if (sModel == 3) mName = "BWAVTDelay";
  else if (sModel == 4) mName = "BWAVLinear";
  else if (sModel == 5) mName = "BWAVLinearDelay";
  else if (sModel == 6) mName = "BWAVLinearFormation";
  mThermo = new Thermodynamics();
  ReadParameters();
  Description();
  mHyperCube = mRapSRange * 2.*Pi() * mRhoMax * mRapPRange * 1.0 * 2.*Pi();
}


Model_BWA::~Model_BWA()
{
  delete mThermo;
}

double Model_BWA::GetIntegrand(ParticleType* aPartType)
{
  double dSigmaP, PdotU;
  double Spin, Statistics;
  double Tau, Rho, PhiS, RapS;
  double Mt,  Pt,  PhiP, RapP;
  double Vt;
  double dPt;
  
// Type of statistics Bose-Einstein or Fermi-Dirac
  Spin	= aPartType->GetSpin();
  Statistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );
// Generate spacial components
  Rho	= mRhoMax    * mRandom->Rndm();
  PhiS	= 2.0 * Pi() * mRandom->Rndm();
  RapS	= mRapSRange * mRandom->Rndm() - 0.5 * mRapSRange;
  Tau   = mTau + mA * Rho;
// Generate momentum components
  {
    double Zet = mRandom->Rndm();
    Pt	= Zet / (1.0 - Zet);
    dPt	= 1.0 / ( (1.0 - Zet) * (1.0 - Zet) );
  }
  PhiP	= 2.0 * Pi() * mRandom->Rndm();
  RapP	= mRapPRange * mRandom->Rndm() - 0.5 * mRapPRange;
  Mt	= Hypot(aPartType->GetMass(), Pt);
// Transverse velocity
  Vt	= mVt;
  if((sModel == 4) || (sModel == 5) || (sModel == 6))
    // Linear Transverse velocity profile
    Vt	= (Rho / mRhoMax) / (mVt + Rho / mRhoMax);
// Invariants
  PdotU	= 1.0 / Sqrt(1 - Vt * Vt) * (Mt * CosH(RapS - RapP) - Vt * Pt * Cos(PhiS - PhiP));
  dSigmaP = Tau * Rho * (Mt * CosH(RapS - RapP) - mA * Pt * Cos(PhiS - PhiP));
// Return values
  Xt = Tau * CosH(RapS);
  Xx = Rho * Cos(PhiS);
  Xy = Rho * Sin(PhiS);
  Xz = Tau * SinH(RapS);
  Pe = Mt  * CosH(RapP);
  Px = Pt  * Cos(PhiP);
  Py = Pt  * Sin(PhiP);
  Pz = Mt  * SinH(RapP);
  if((sModel == 3) || (sModel == 5))
    Xt += -mDelay * Log(mRandom->Rndm());
  else if(sModel == 6) {
    double Energy = Hypot(Mt,Pz);
    Xt += -mDelay * Log(mRandom->Rndm());
    Xx += Xt * Px / Energy;
    Xy += Xt * Py / Energy;
    Xz += Xt * Pz / Energy;
  } 
// integrand
  return (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * Pt * dPt * dSigmaP * 1.0 / (Exp( ( PdotU - mThermo->GetChemicalPotential(aPartType) ) / mThermo->GetTemperature() ) + Statistics);
}

void Model_BWA::Description()
{
  ostringstream oss;
  oss << "##################################################"<< endl;
  oss << MODEL_NAME(mName);
  oss << "# - rapidity range         : " <<MODEL_PAR_DESC(mRapPRange,		"[units]");
  oss << "# - spatial rapidity range : " <<MODEL_PAR_DESC(mRapSRange,		"[units]");
  oss << "# - max cylinder size      : " <<MODEL_PAR_DESC(mRhoMax * kHbarC,	"[fm]");
  oss << "# - Blast-Wave time        : " <<MODEL_PAR_DESC(mTau    * kHbarC,	"[fm]");
  oss << "# - transverse velocity    : " <<MODEL_PAR_DESC(mVt,		"[c]");
  oss << "# - parameter A            : " <<MODEL_PAR_DESC(mA,			"[1]");
  if((sModel == 3) || (sModel == 5) || (sModel == 6))
    oss << "# - delay                  : " <<MODEL_PAR_DESC(mDelay  * kHbarC,	"[fm]");
  oss << "# - freeze-out temperature : " <<MODEL_PAR_DESC(mThermo->GetTemperature() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_B   : " <<MODEL_PAR_DESC(mThermo->GetMuB() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_I3  : " <<MODEL_PAR_DESC(mThermo->GetMuI() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_S   : " <<MODEL_PAR_DESC(mThermo->GetMuS() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_C   : " <<MODEL_PAR_DESC(mThermo->GetMuC() * 1000.0,	"[MeV]");
  oss << "# Parameters hash (CRC32)  : " <<MODEL_PAR_DESC(mHash,		"");
  oss << "# Integration samples      : " <<MODEL_PAR_DESC(sIntegrateSample,	"");
  oss << "# Random seed              : " <<MODEL_PAR_DESC((sRandomize ? "yes" : "no"),"");
  oss << "# Generation date          : " <<sTimeStamp<<" #"<<endl;
  oss << "##################################################"<< endl;  
  mDescription = oss.str();
}

void Model_BWA::AddParameterBranch(TTree* aTree)
{
  Model_t_BWA tPar;
  
  tPar.RapPRange = mRapPRange;
  tPar.RapSRange = mRapSRange;
  tPar.RhoMax    = mRhoMax  * kHbarC;
  tPar.Tau       = mTau     * kHbarC;
  tPar.Vt        = mVt;
  tPar.A         = mA;
  tPar.Delay     = mDelay   * kHbarC;
  tPar.Temp      = mThermo->GetTemperature() * 1000.0;
  tPar.MuB       = mThermo->GetMuB() * 1000.0;
  tPar.MuI       = mThermo->GetMuI() * 1000.0;
  tPar.MuS       = mThermo->GetMuS() * 1000.0;
  tPar.MuC       = mThermo->GetMuC() * 1000.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_BWA_)->Fill();
}

void Model_BWA::ReadParameters()
{
  Configurator*	tModelParam;
  Parser*	tParser;
  
  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;
  
  try {
    mRapPRange	 = tModelParam->GetParameter("RapPRange").Atof();		// [1]
    mRapSRange	 = tModelParam->GetParameter("RapSRange").Atof();		// [1]
    mRhoMax	 = tModelParam->GetParameter("RhoMax").Atof()	/ kHbarC;	// [GeV^-1]
    mTau	 = tModelParam->GetParameter("Tau").Atof()	/ kHbarC;	// [GeV^-1]
    mVt		 = tModelParam->GetParameter("VelT").Atof();			// [c]
    mA		 = tModelParam->GetParameter("ParA").Atof();			// [1]
    mDelay	 = tModelParam->GetParameter("Delay").Atof()	/ kHbarC;	// [GeV^-1]
    mThermo->SetTemperature(tModelParam->GetParameter("Temperature").Atof() * 0.001);	// [GeV]
    mThermo->SetChemistry(  tModelParam->GetParameter("MuB").Atof() * 0.001,
			    tModelParam->GetParameter("MuI").Atof() * 0.001,
			    tModelParam->GetParameter("MuS").Atof() * 0.001,
			    tModelParam->GetParameter("MuC").Atof() * 0.001);	// [GeV]
  } catch (TString tError) {
    PRINT_MESSAGE("<Model_BWA::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }

// calculate parameter hash
  ostringstream oss;
  oss << sModel;
  oss << mRapPRange << mRapSRange << mRhoMax << mTau << mVt << mA << mDelay;
  oss << mThermo->GetTemperature() << mThermo->GetMuB() << mThermo->GetMuI() << mThermo->GetMuS() << mThermo->GetMuC();
  CalculateHash(TString(oss.str()));

// create event subdirectory if needed
  try {
    sEventDIR += tModelParam->GetParameter("EventSubDir");
    CreateEventSubDir();
  } catch (TString tError) {
  }
  
  delete tModelParam;
}

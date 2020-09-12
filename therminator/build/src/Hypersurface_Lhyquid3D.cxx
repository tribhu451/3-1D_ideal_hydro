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

#include <TMath.h>
#include "Configurator.h"
#include "Parser.h"
#include "Hypersurface_Library.h"
#include "Hypersurface_Lhyquid3D.h"
#include "THGlobal.h"

using namespace std;
using namespace TMath;

extern void AddLogEntry(const char* aEntry);
extern Configurator* sMainConfig;
extern TString	sModelINI;
extern TString	sHyperXML;
extern TString	sEventDIR;

Hypersurface_Lhyquid3D::Hypersurface_Lhyquid3D()
{
  mThermo = new Thermodynamics();
  ReadParameters();
}

Hypersurface_Lhyquid3D::~Hypersurface_Lhyquid3D()
{
  delete mThermo;
  delete mDistance;
  delete mDistanceDPhi;
  delete mDistanceDZeta;
  delete mDistanceDTheta;
  delete mFluidUx;
  delete mFluidUy;
  delete mFluidRapidity;
}

double Hypersurface_Lhyquid3D::GetDSigmaP(double aMt, double aPt, double aPhiP, double aRapP)
{
  return 
    Dhs * Dhs * Sin(Theta) * (
      Tau / mLambda * (
          dDdZeta * Cos(Zeta)							* (-aMt * Cos(Zeta) * CosH(aRapP - RapS) + aPt * Sin(Zeta) * Cos(PhiS - aPhiP))
        + (Dhs * Sin(Theta) - dDdTheta * Cos(Theta)) * Cos(Zeta) * Sin(Theta)	* ( aMt * Sin(Zeta) * CosH(aRapP - RapS) + aPt * Cos(Zeta) * Cos(PhiS - aPhiP))
        + dDdPhi * aPt * Sin(PhiS - aPhiP)
      ) 
      +   (Dhs * Cos(Theta) + dDdTheta * Sin(Theta)) 		 * Sin(Theta)	*   aMt * Cos(Zeta) * SinH(aRapP - RapS)
    );
}

double Hypersurface_Lhyquid3D::GetPdotU(double aMt, double aPt, double aPhiP, double aRapP)
{
  return Sqrt(1 + Ux * Ux + Uy * Uy) * aMt * CosH(RapF - aRapP) - aPt * (Ux * Cos(aPhiP) + Uy * Sin(aPhiP));
}

double Hypersurface_Lhyquid3D::GetHyperCubeSpatialVolume()
{
  return (mDistance->GetXMax() - mDistance->GetXMin()) * (mDistance->GetYMax() - mDistance->GetYMin()) * (mDistance->GetZMax() - mDistance->GetZMin());
}

void Hypersurface_Lhyquid3D::SetPositionOnHypersurface(TRandom2* aRandom, double* aXt, double* aXx, double* aXy, double* aXz)
{
  Zeta    = mDistance->GetXMin() + (mDistance->GetXMax() - mDistance->GetXMin()) * aRandom->Rndm();
  PhiS    = mDistance->GetYMin() + (mDistance->GetYMax() - mDistance->GetYMin()) * aRandom->Rndm();
  Theta   = mDistance->GetZMin() + (mDistance->GetZMax() - mDistance->GetZMin()) * aRandom->Rndm();
  Dhs     = mDistance      ->Interpolate(Zeta, PhiS, Theta);
  dDdZeta = mDistanceDZeta ->Interpolate(Zeta, PhiS, Theta);
  dDdPhi  = mDistanceDPhi  ->Interpolate(Zeta, PhiS, Theta);
  dDdTheta= mDistanceDTheta->Interpolate(Zeta, PhiS, Theta);
  Ux	  = mFluidUx       ->Interpolate(Zeta, PhiS, Theta);
  Uy	  = mFluidUy       ->Interpolate(Zeta, PhiS, Theta);
  RapF	  = mFluidRapidity ->Interpolate(Zeta, PhiS, Theta);
  Tau     = mTauI + Dhs * Sin(Theta) * Sin(Zeta);
  Rho     = Dhs * Sin(Theta) * Cos(Zeta);
  RapS	  = Dhs * Cos(Theta) / mLambda;
  (*aXt)  = Tau * CosH(RapS);
  (*aXx)  = Rho * Cos(PhiS);
  (*aXy)  = Rho * Sin(PhiS);
  (*aXz)  = Tau * SinH(RapS);
}

void Hypersurface_Lhyquid3D::ReadParameters()
{
  Hypersurface_Library*	tLib;
  Configurator*		tModelParam;
  Parser*		tParser;
  
  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;
  
  try {
    mRapPRange  = tModelParam->GetParameter("RapPRange").Atof();						// [1]
    if(sHyperXML.IsNull()) {
      sHyperXML = sMainConfig->GetParameter("FreezeOutDir");
      sHyperXML.Prepend("./");
      sHyperXML += tModelParam->GetParameter("FreezeFile");
    }
  } catch (TString tError) {
    PRINT_MESSAGE("<Hypersurface_Lhyquid3D::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
 
  char tBuff[2*kFileNameMaxChar];
  sprintf(tBuff,"[input]\t%s",sHyperXML.Data()); 
  AddLogEntry(tBuff);
//###############################################
//		XML file with hypersurface	#
//###############################################

  PRINT_DEBUG_1("Hydro based model input. FreezeFile: " << sHyperXML.Data());
  tLib    = new Hypersurface_Library;
  tParser = new Parser(sHyperXML.Data());
  tParser->ReadXML(tLib);
  delete tParser;
  
  try {
    mTauI	     = tLib->GetXMLTag("PARAMETER","name","Tau_i")->GetXMLContent().Atof()  / kHbarC;		// [GeV^-1]
    mLambda	     = tLib->GetXMLTag("PARAMETER","name","Lambda")->GetXMLContent().Atof() / kHbarC;		// [GeV^-1]
    mThermo->SetTemperature(tLib->GetXMLTag("PARAMETER","name","Temperature")->GetXMLContent().Atof() * 0.001);	// [GeV]
    mThermo->SetChemistry(  tLib->GetXMLTag("PARAMETER","name","Mu_B")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_I")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_S")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_C")->GetXMLContent().Atof() * 0.001);	// [GeV]
    mCollidingEnergy = tLib->GetXMLTag("PARAMETER","name","colliding_energy")->GetXMLContent().Atof();		// [GeV]
    mCentralityMin   = tLib->GetXMLTag("PARAMETER","name","centrality_min")->GetXMLContent().Atof();		// [%]
    mCentralityMax   = tLib->GetXMLTag("PARAMETER","name","centrality_max")->GetXMLContent().Atof();		// [%]
    mImpactParameter = tLib->GetXMLTag("PARAMETER","name","impact_parameter")->GetXMLContent().Atof();		// [fm]
    mTempI	     = tLib->GetXMLTag("PARAMETER","name","temperature_at_center")->GetXMLContent().Atof();	// [MeV]
    sprintf(mDeviceName,     "%s",tLib->GetXMLTag("PARAMETER","name","device")->GetXMLContent().Data());
    sprintf(mCollidingSystem,"%s",tLib->GetXMLTag("PARAMETER","name","colliding_system")->GetXMLContent().Data());
    mDistance	     = tLib->GetXMLTag("VECTOR3D", "name","Distance")->GetXMLVector3D();			// [GeV^-1]
    mFluidUx	     = tLib->GetXMLTag("VECTOR3D", "name","FluidUx" )->GetXMLVector3D();			// [1]
    mFluidUy	     = tLib->GetXMLTag("VECTOR3D", "name","FluidUy" )->GetXMLVector3D();			// [1]
    mFluidRapidity   = tLib->GetXMLTag("VECTOR3D", "name","FluidRap")->GetXMLVector3D();			// [1]
  } catch (int tError) {
    PRINT_MESSAGE("<Hypersurface_Lhyquid3D::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the XML file.");
    exit(_ERROR_LIBRARY_TAG_NOT_FOUND_);
  }

// mDistance derivatives
  try {
    mDistanceDZeta = tLib->GetXMLTag("VECTOR3D","name","DistanceDZeta")->GetXMLVector3D();		// [GeV^-1/rad]
  } catch (int tError) {
    PRINT_DEBUG_1("<Hypersurface_Lhyquid3D::ReadParameters>\tCalculating derivative Distance->DerivativeX()");
    mDistanceDZeta = mDistance->DerivativeX("DistanceDZeta");   
  }
  try {
    mDistanceDPhi = tLib->GetXMLTag("VECTOR3D","name","DistanceDPhi")->GetXMLVector3D();		// [GeV^-1/rad]
  } catch (int tError) {
    PRINT_DEBUG_1("<Hypersurface_Lhyquid3D::ReadParameters>\tCalculating derivative Distance->DerivativeY()");
    mDistanceDPhi = mDistance->DerivativeY("DistanceDPhi");   
  }
  try {
    mDistanceDTheta = tLib->GetXMLTag("VECTOR3D","name","DistanceDTheta")->GetXMLVector3D();		// [GeV^-1/rad]
  } catch (int tError) {
    PRINT_DEBUG_1("<Hypersurface_Lhyquid3D::ReadParameters>\tCalculating derivative Distance->DerivativeZ()");
    mDistanceDTheta = mDistance->DerivativeZ("DistanceDTheta");   
  }
  
// event subdirectory
  try {
    sEventDIR += tModelParam->GetParameter("EventSubDir");
  }
  catch (TString tError) {
    TString tTemp = tModelParam->GetParameter("FreezeFile");
    tTemp.ReplaceAll("/","-");
    tTemp.ReplaceAll(".xml","/");
    sEventDIR += tTemp;
  }
  
  delete tLib; 
  delete tModelParam;
}

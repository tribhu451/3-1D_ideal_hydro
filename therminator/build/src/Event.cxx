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
#include <TDatime.h>
#include "Crc32.h"
#include "Configurator.h"
#include "ParticleDecayer.h"
#include "Event.h"
#include "THGlobal.h"

extern Configurator* sMainConfig;
extern TString	sTimeStamp;
extern int	sRandomize;

using namespace std;

Event::Event()
: mPartDB(0), mInteg(0), mRandom(0), mDistribution(0)
{
  mMultiplicities.clear();
  Reset();
}

Event::Event(ParticleDB* aDB, Integrator* aInteg)
: mPartDB(aDB), mInteg(aInteg), mDistribution(0)
{ 
  mRandom = new TRandom2();
#ifdef _ROOT_4_
  mRandom->SetSeed2(31851, 14327);
#else
  mRandom->SetSeed(31851);
#endif
  mMultiplicities.clear();
  mMultiplicities.resize(mPartDB->GetParticleTypeCount());
  Reset();
  ReadParameters();
}

Event::~Event()
{
  mParticles.clear();
  mMultiplicities.clear();
  delete mRandom;
}

void Event::Reset(int aEventIter)
{
  ostringstream oss;
  Crc32 tEventID;
  
  mParticles.clear();
  Particle::ZeroEID();
  
  oss << sTimeStamp.Data() << "Event: " << aEventIter;
  tEventID.Update(oss.str().data(), oss.str().length());
  tEventID.Finish(); 
  mEventID = tEventID.GetValue();
}

list<Particle>* Event::GetParticleList()
{
  return &mParticles;
}

Integrator* Event::GetIntegrator() const
{
  return mInteg;
}

ParticleDB* Event::GetParticleDB() const
{
  return mPartDB;
}

unsigned int Event::GetEventID() const
{
  return mEventID;
}

void Event::GeneratePrimordials(int aSeed)
{ 
#ifdef _ROOT_4_
  if (aSeed) mRandom->SetSeed2(aSeed, (aSeed*2) % (7*11*23*31));
#else
  if (aSeed) mRandom->SetSeed(aSeed);
#endif

  GenerateMultiplicities();
  for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
    if(! strstr(mPartDB->GetParticleType(tIter)->GetName(),"gam000zer")) // disable primordial photons production
      mInteg->GenerateParticles(mPartDB->GetParticleType(tIter), mMultiplicities[tIter], &mParticles);
}

void Event::DecayParticles(int aSeed)
{
  list<Particle>::iterator tIter;
  ParticleType*    tFatherType;
  ParticleDecayer* tDecayer;
  
  tDecayer = new ParticleDecayer(mPartDB, &mParticles);

  if (sRandomize)
    tDecayer->Randomize();
  else
    tDecayer->SeedSet(aSeed);
  
  tIter = mParticles.begin();
// as new particles are added from decays the end() of the list moves until all particles had decayed
  do {
    tFatherType = tIter->GetParticleType();
    // if not stable or stable but has a decay table with at least one decay channel
    if((tFatherType->GetGamma() >= 0.0) && (tFatherType->GetTable()) && (tFatherType->GetTable()->GetChannelCount() + 1 > 0))
      tDecayer->DecayParticle( &(*tIter) );
    tIter++;
  } while (tIter != mParticles.end());
  delete tDecayer;
}

void Event::GenerateMultiplicities()
{
  if(mDistribution == 0) { // Poisson
    for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
      mMultiplicities[tIter] = mRandom->Poisson(mPartDB->GetParticleType(tIter)->GetMultiplicity());
  } else if(mDistribution == 1) { // Negative Binomial
    for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
      mMultiplicities[tIter] = 0; // HOW?
  }
}

void Event::Randomize()
{
  TDatime tDate;

#ifdef _ROOT_4_
  mRandom->SetSeed2(tDate.Get() / 2 * 3, tDate.Get() / 11 * 9);
#else
  mRandom->SetSeed(tDate.Get() / 2 * 3);
#endif
}

void Event::ReadParameters()
{
  TString tDistribution; 
  try {
    tDistribution	= sMainConfig->GetParameter("MultiplicityDistribution");
    if (tDistribution.Contains("NegativeBinomial"))
      mDistribution = 1;
  }
  catch (TString tError) {
    PRINT_DEBUG_1("<Event::ReadParameters>\tUsing default multiplicity distribution: Poissonian");
  }
}

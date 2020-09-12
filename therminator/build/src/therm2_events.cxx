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

#include <fstream>
#include <TString.h>
#include <TDatime.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "EventGenerator.h"

using namespace std;

Configurator *sMainConfig;
TString	sMainINI;
TString	sModelINI;
TString sHyperXML;
TString	sEventDIR;
TString	sTimeStamp;
int	sModel;
int	sRandomize;
int	sIntegrateSample;
int	sParentPID;

void ReadParameters();
void ReadSHARE(ParticleDB* aPartDB);
void CheckSHARE(ParticleDB* aPartDB);
void MessageIntro();
void MessageHelp();
void MessageVersion();
void CopyINIFile();
void AddLogEntry(const char* aEntry);

int main(int argc, char **argv)
{
  ParticleDB*	  tPartDB;
  EventGenerator* tEventGen;

  sMainINI   = "./events.ini";
  sParentPID = 0;
  sHyperXML  = "";
/*
## ParentPID ##
Applicable if the user runs multiple copies of therm2_events - e.g. a cluster.
A temporary file is created named "event_<tParentPID>.tmp" with two lines:
 - current event storage directory
 - number of event*.root files generated by this run + files already in that directory (previous runs)
That information can be passed to other programs i.e. ROOT figures or HBT in one BASH script
*/
  if (argc > 1) {
    TString tDummy;
    for(int i=1; i<argc;i++) {
      tDummy = argv[i];
      if((tDummy == "-h") || (tDummy == "--help")) {
        MessageHelp();
        return 0;
      } else if((tDummy == "-v") || (tDummy == "--version")) {
        MessageVersion();
        return 0;
      } else if (tDummy.EndsWith(".xml"))
        sHyperXML = tDummy;
      else if (tDummy.EndsWith(".ini"))
	sMainINI  = tDummy;
      else if (tDummy.IsDigit())
	sParentPID = tDummy.Atoi();
    }
  }

  MessageIntro();

  sMainConfig = new Configurator;
  ReadParameters();

  {
    char tBuff[2*kFileNameMaxChar];
    std::sprintf(tBuff,"[input]\t%s\t%i",sMainINI.Data(),sParentPID);
    AddLogEntry(tBuff);
    std::sprintf(tBuff,"[input]\t%s",sModelINI.Data());
    AddLogEntry(tBuff);
  }

  tPartDB     = new ParticleDB();
  ReadSHARE(tPartDB);
  //CheckSHARE(tPartDB);

  tEventGen   = new EventGenerator(tPartDB);
  tEventGen->GenerateEvents();
  tEventGen->SetEventsTemp();

  delete tEventGen;
  delete tPartDB;
  delete sMainConfig;

  return 0;
}

// ##############################################################
// #			---===[ END OF MAIN ]===--		#
// ##############################################################

void ReadParameters()
{
  TString tModel;
  TString tModelINI;
  TDatime tDate;
  Parser* tParser;

  tDate.Set();
  sTimeStamp = tDate.AsSQLString();

  tParser = new Parser(sMainINI);
  tParser->ReadINI(sMainConfig);
  delete tParser;
  try {
    sRandomize = sMainConfig->GetParameter("Randomize").Atoi();
    tModelINI  = sMainConfig->GetParameter("FreezeOutDir"); tModelINI.Prepend("./");
    tModel     = sMainConfig->GetParameter("FreezeOutModel");
    if      (tModel == "KrakowSFO")		{ sModel = 0;  tModelINI += "krakow.ini";	}
    else if (tModel == "BlastWave")		{ sModel = 1;  tModelINI += "blastwave.ini";	}
    else if (tModel == "BWAVT")			{ sModel = 2;  tModelINI += "bwa.ini";		}
    else if (tModel == "BWAVTDelay")		{ sModel = 3;  tModelINI += "bwa.ini";		}
    else if (tModel == "BWAVLinear")		{ sModel = 4;  tModelINI += "bwa.ini";		}
    else if (tModel == "BWAVLinearDelay")	{ sModel = 5;  tModelINI += "bwa.ini";		}
    else if (tModel == "BWAVLinearFormation")	{ sModel = 6;  tModelINI += "bwa.ini";		}
    else if (tModel == "Lhyquid3D")		{ sModel = 10; tModelINI += "lhyquid3d.ini";	}
    else if (tModel == "Lhyquid2DBI")		{ sModel = 11; tModelINI += "lhyquid2dbi.ini";	}
/**********************************************************************************************
 // [1] Associate a number (sModel) and ini file name (sModelINI) to a your Model_* name
 // given in the therminator.ini file. See Integrator.cxx for further instructions.
 else if (tModel == "Example")		{ sModel = 99; sModelINI = "example.ini";	}
 **********************************************************************************************/
    else {
      PRINT_MESSAGE("<therm2_events::ReadParameters>\tUnknown model type: " << tModel);
      PRINT_MESSAGE("\tPlease provide the proper model name in the events.ini file.");
      exit(_ERROR_GENERAL_MODEL_UNKNOWN_);
    }
  } catch (TString tError) {
    PRINT_MESSAGE("<therm2_events::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  // Custom Model ini file
  try {
    sModelINI  = sMainConfig->GetParameter("FreezeOutModelINI");
    PRINT_MESSAGE("<therm2_events::ReadParameters>\tUsing custom Freeze-Out-Model INI file " << sModelINI);
  }
  catch (TString tError) {
    sModelINI  = tModelINI;
  }
}

void ReadSHARE(ParticleDB* aPartDB)
{
  TString tShareDir;
  Parser* tParser;

  try {
    tShareDir = sMainConfig->GetParameter("ShareDir"); tShareDir.Prepend("./");
  } catch (TString tError) {
    PRINT_DEBUG_1("<Parser::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find SHARE input file location.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }

  tParser = new Parser((tShareDir + "particles.data").Data());
  tParser->ReadSHAREParticles(aPartDB);
  delete tParser;

  tParser = new Parser((tShareDir + "decays.data").Data());
  tParser->ReadSHAREDecays(aPartDB);
  delete tParser;
}

void MessageIntro()
{
  PRINT_MESSAGE("  ***********************************************************************"	);
  PRINT_MESSAGE("  *\t\tTHERMINATOR 2 EVENTS version "<<_THERMINATOR2_VERSION_<<"\t\t\t*"	);
  PRINT_MESSAGE("  *\t\t\t\t\t\t\t\t\t*"							);
  PRINT_MESSAGE("  * authors: M.Chojnacki, A.Kisiel, W.Florkowski, W.Broniowski\t\t*"		);
  PRINT_MESSAGE("  * cite as: arXiv:1102.0273\t\t\t\t\t\t*"					);
  PRINT_MESSAGE("  * webpage: http://therminator2.ifj.edu.pl/\t\t\t\t*"				);
  PRINT_MESSAGE("  ***********************************************************************"	);
}

void MessageHelp()
{
  PRINT_MESSAGE("Usage:");
  PRINT_MESSAGE("therm2_events [EVENTS_INI] [PPID] [HYPER_XML]");
  PRINT_MESSAGE("therm2_events [OPTION]");
  PRINT_MESSAGE("  [EVENTS_INI]\t\tmain settings file;\t\tdefault: events.ini");
  PRINT_MESSAGE("  [PPID]\t\tparent's system process ID;\tdefault: 0");
  PRINT_MESSAGE("  [HYPER_XML]\tlocation of the hypersurface XML file;\tdefault:");
  PRINT_MESSAGE("  [OPTION]");
  PRINT_MESSAGE("    -h | --help\t\tthis screen");
  PRINT_MESSAGE("    -v | --version\tversion information");
}

void MessageVersion()
{
  PRINT_MESSAGE("version:\tTHERMINATOR 2 EVENTS version "<<_THERMINATOR2_VERSION_);
  PRINT_MESSAGE("compiled with:\t"<<_CXX_VER_<<", ROOT("<<_ROOT_VER_<<")");
  std::cout <<  "  preprocessor: ";
#ifdef _DEBUG_LEVEL_
  std::cout << "DEBUG="<<_DEBUG_LEVEL_<<" ";
#endif
#ifdef _MODEL_LHYQUID_ONLY_BACK_FLOW_
  std::cout << "BACK_FLOW="<<_MODEL_LHYQUID_ONLY_BACK_FLOW_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_RESCALE_CHANNELS_
  std::cout << "RESCALE_CHANNELS="<<_PARTICLE_DECAYER_RESCALE_CHANNELS_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_
  std::cout << "DISABLE_THREE_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_<<" ";
#endif
#ifdef _PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_
  std::cout << "DISABLE_TWO_BODY_DECAYS="<<_PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_<<" ";
#endif
  std::cout << std::endl;
}

void CopyINIFile()
{
  TString  tINI;
  ifstream ifs;
  ofstream ofs;

  tINI = sMainINI;
  tINI.ReplaceAll("./",sEventDIR);
  ifs.open(sMainINI.Data(), std::ios::binary);
  ofs.open(tINI.Data(), std::ios::binary);
  if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) {
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
  } else
    PRINT_MESSAGE("<therm2_events::CopyINIFile>\tUnable to copy "<<sMainINI<<" to "<<tINI);

  tINI = sModelINI;
  tINI.ReplaceAll(("./" + sMainConfig->GetParameter("FreezeOutDir")),sEventDIR);
  ifs.open(sModelINI.Data(), std::ios::binary);
  ofs.open(tINI.Data(), std::ios::binary);
  if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) {
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
  } else {
    PRINT_MESSAGE("<therm2_events::CopyINIFile>\tUnable to copy "<< sModelINI<<" to "<<tINI);
  }
}

void AddLogEntry(const char* aEntry)
{
  TString tLogName;
  TDatime tDate;
  ofstream tFile;

  tDate.Set();
  try {
    tLogName = sMainConfig->GetParameter("LogFile"); tLogName.Prepend("./");
  }
  catch (TString tError) {
    return;
  }

  tFile.open(tLogName, std::ios_base::app);
  if (static_cast<long>(tFile.tellp()) == 0) {
    tFile << "# THERMINATOR 2 Log File"<<std::endl;
  }
  tFile << '['<<tDate.AsSQLString()<<"]\ttherm2_events\t"<<sParentPID<<'\t';
  tFile << aEntry << std::endl;
  tFile.close();
}

void CheckSHARE(ParticleDB* aPartDB) {
  ParticleType*	tType;
  DecayTable*	tDecTable;
  double	SumBR;

  PRINT_DEBUG_2("<therm2_events::ParserCheck>\tRead "<< aPartDB->GetParticleTypeCount()<<" particle types.");
  for(int tPart=0; tPart<aPartDB->GetParticleTypeCount(); tPart++) {
    tType = aPartDB->GetParticleType(tPart);
    PRINT_DEBUG_2("\tParticle " << tType->GetNumber() << ": " << tType->GetName()
      <<", Mass = "	<<tType->GetMass()
      <<", Gamma = "	<<tType->GetGamma()
      <<", Spin = "	<<tType->GetSpin()
      <<", I  = "	<<tType->GetI()
      <<", I3 = "	<<tType->GetI3()
      <<", BarionN = "	<<tType->GetBarionN()
      <<", StrangeN = "	<<tType->GetStrangeN()
      <<", CharmN = "	<<tType->GetCharmN()
      <<", Charge = "	<<tType->GetCharge()
      <<", MC# = "	<<tType->GetPDGCode()
    );

    SumBR = 0.0;
    if (tDecTable = tType->GetTable()) {
      for (int tChanIndex = 0; tChanIndex < tDecTable->GetChannelCount() + 1; tChanIndex++) {
	if (tDecTable->GetDecayChannel(tChanIndex)->Is3Particle()) {
          PRINT_DEBUG_2("\t\tChannel " << tChanIndex << ": "
	    << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChanIndex)->GetParticle1()))->GetName() << " + "
	    << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChanIndex)->GetParticle2()))->GetName() << " + "
	    << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChanIndex)->GetParticle3()))->GetName()
	    << ", BR = " << tType->GetTable()->GetDecayChannel(tChanIndex)->GetBranchingRatio()
	    << ", Step = " << tType->GetTable()->GetDecayStep(tChanIndex)
	  );
	} else {
	  PRINT_DEBUG_2("\t\tChannel " << tChanIndex << ": "
	    << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChanIndex)->GetParticle1()))->GetName() << " + "
	    << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChanIndex)->GetParticle2()))->GetName()
	    << ", BR = " << tType->GetTable()->GetDecayChannel(tChanIndex)->GetBranchingRatio()
	    << ", Step = " << tType->GetTable()->GetDecayStep(tChanIndex)
	  );
	}
        SumBR += tType->GetTable()->GetDecayChannel(tChanIndex)->GetBranchingRatio();
      }
      PRINT_DEBUG_2("\t\tSum BR = " << SumBR);
    }
  }
}

/*! @file therm2_events.cxx
 * @brief <c><b>THERMINATOR 2</b></c> event generating program.
 *
 * @fn int main(int argc, char **argv)
 * @brief Main program.
 * @param [in] argc number of program arguments
 * @param [in] argv vector with program arguments
 *
 * @fn void ReadParameters()
 * @brief Reads some of the global parameters.
 *
 * The <c><b>THERMINATOR 2</b></c> main settings file <b>therminator.ini</b> options used here:
 * <table>
 *   <tr><th>Keyword</th>		<th>Description</th></tr>
 *   <tr><td>Randomize</td>		<td>all random generators are initiated with system date or have constant seed </td></tr>
 *   <tr><td>FreezeOutDir</td>		<td>directory with freeze-out models setting files</td></tr>
 *   <tr><td>FreezeOutModel</td>	<td>name of the freeze-out model to use</td></tr>
 * </table>
 *
 * @fn void ReadSHARE(ParticleDB* aPartDB)
 * @brief Reads the <c><b>SHARE</b></c> database and builds the ParticleDB.
 *
 * The <c><b>THERMINATOR 2</b></c> main settings file <b>therminator.ini</b> options used here:
 * <table>
 *   <tr><th>Keyword</th>		<th>Description</th></tr>
 *   <tr><td>ShareDir</td>		<td>directory with the <c><b>SHARE</b></c> database files</td></tr>
 * </table>
 * @param [in] aPartDB particles database
 *
 * @fn void CheckSHARE(ParticleDB* aPartDB)
 * @brief Procedure to check if Parser has correctly filled ParticleDB.
 * @param [in] aPartDB particles database
 *
 * @fn void MessageIntro()
 * @brief Prints introduction message.
 *
 * @fn void MessageHelp();
 * @brief Prints help message.
 *
 * @fn void MessageVersion();
 * @brief Prints version message.
 *
 * @fn void CopyINIFile()
 * @brief Copy INI file to sEventDIR
 *
 * @fn void AddLogEntry(const char* aEntry)
 * @brief Adds a log entry to the log file
 * @param [in] aEntry text string with the entry
 *
 * @var Configurator *sMainConfig
 * @brief Global variable with Configurator holding settings form <b>events.ini</b> file.
 *
 * @var TString sEventDIR
 * @brief Global variable with directory to store event files.
 *
 * @var TString	sMainINI
 * @brief Global variable with name of the main ini file [<i>default: therminator.ini</i>]
 *
 * @var TString sModelINI
 * @brief Global variable with name and location of selected model ini file.
 *
 * @var TString sHyperXML
 * @brief Global variable with name and location of selected hypersurface XML file.
 *
 * @var char sTimeStamp[21]
 * @brief Global variable with time-stamp character string. (default format "YYYY-MM-DD hh:mm:ss")
 *
 * @var int sModel
 * @brief Global variable with number associated to a chosen freeze-out Model.
 *
 * @var int sRandomize
 * @brief Global variable that controls the random number generators initiation.
 *
 * @var int sIntegrateSample
 * @brief Global variable with number of Monte-Carlo samples used by Integrator to find the average multiplicity and max of the
 * Cooper-Frye integrand for each ParticleType in the selected Model.
 *
 * @var int sParentPID
 * @brief Parent's system process ID.
*/

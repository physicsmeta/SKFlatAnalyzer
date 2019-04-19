#ifndef HNtypeI_WZ_CR_h
#define HNtypeI_WZ_CR_h

#include "AnalyzerCore.h"

class HNtypeI_WZ_CR : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  double weight_Prefire;

  HNtypeI_WZ_CR();
  ~HNtypeI_WZ_CR();

};



#endif


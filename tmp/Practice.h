#ifndef Practice_h
#define Practice_h

#include "AnalyzerCore.h"

class Practice : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param, Long64_t Nentry); 
  void executeEvent(Long64_t Nentry); 

  TString EleTriggerName, IsoMuTriggerName;
  double ele0ptcut, ele1ptcut, MuTriggerSafePtCut;

  vector<TString> EleIDs, EleIDSFKeys;
  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Electron> AllEles;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  double weight_Prefire;

  Practice();
  ~Practice();

};



#endif


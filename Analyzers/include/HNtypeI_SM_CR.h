#ifndef HNtypeI_SM_CR_h
#define HNtypeI_SM_CR_h

#include "HNAnalyzerCore.h"

class HNtypeI_SM_CR : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  
  vector<TString> MuonTriggers;
  vector<TString> ElectronTriggers;
  double MuonPtCut1;
  double MuonPtCut2;
  double ElectronPtCut1;
  double ElectronPtCut2;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

//  double weight_Prefire;

  HNtypeI_SM_CR();
  ~HNtypeI_SM_CR();

};

#endif

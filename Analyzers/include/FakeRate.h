#ifndef FakeRate_h
#define FakeRate_h

#include "AnalyzerCore.h"

class FakeRate : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;
  
  vector<TString> MuonTriggers;
  vector<TString> ElectronTriggers;
  double MuonPtCut1, MuonPtCut2, MuonPtCut3;
  double MuonPtconeCut1, MuonPtconeCut2, MuonPtconeCut3;
  double ElectronPtCut1, ElectronPtCut2, ElectronPtCut3;
  double ElectronPtconeCut1, ElectronPtconeCut2, ElectronPtconeCut3;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

//  double weight_Prefire;

  FakeRate();
  ~FakeRate();

};

#endif

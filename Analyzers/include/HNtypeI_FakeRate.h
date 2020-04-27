#ifndef HNtypeI_FakeRate_h
#define HNtypeI_FakeRate_h

#include "HNAnalyzerCore.h"

class HNtypeI_FakeRate : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  // Trigger  
  vector<TString> MuonTriggers;
  vector<TString> ElectronTriggers;

  // Lepton ID
  vector<TString> MuonVetoIDs;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronVetoIDs;
  vector<TString> ElectronLooseIDs;
  vector<TString> ElectronTightIDs;

  TString MuonTrig1, MuonTrig2, MuonTrig3;
  TString ElectronTrig1, ElectronTrig2, ElectronTrig3, ElectronTrig4;

  // Lepton pT, pTcone cut
  double MuonPtCut1, MuonPtCut2, MuonPtCut3;
  double MuonPtconeCut1, MuonPtconeCut2, MuonPtconeCut3;
  double ElectronPtCut1, ElectronPtCut2, ElectronPtCut3, ElectronPtCut4;
  double ElectronPtconeCut1, ElectronPtconeCut2, ElectronPtconeCut3, ElectronPtconeCut4;

  // Luminosity
  double MuonLumi1, MuonLumi2, MuonLumi3;
  double ElectronLumi1, ElectronLumi2, ElectronLumi3, ElectronLumi4;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

//  double weight_Prefire;

  HNtypeI_FakeRate();
  ~HNtypeI_FakeRate();

};

#endif

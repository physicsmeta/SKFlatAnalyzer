#ifndef HNtypeI_SR_h
#define HNtypeI_SR_h

#include "AnalyzerCore.h"

class HNtypeI_SR : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst, RunFake, RunCF;
  bool RunNewPDF;
  bool RunXSecSyst;

  
  vector<TString> MuonTriggers;
  vector<TString> ElectronTriggers;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronLooseIDs;
  vector<TString> ElectronTightIDs;
  vector<TString> FakeRateIDs;
  double MuonPtCut1;
  double MuonPtCut2;
  double ElectronPtCut1;
  double ElectronPtCut2;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<FatJet> AllFatJets;

//  double weight_Prefire;

  HNtypeI_SR();
  ~HNtypeI_SR();

};

#endif

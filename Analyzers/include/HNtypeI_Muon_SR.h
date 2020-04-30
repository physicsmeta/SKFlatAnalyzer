#ifndef HNtypeI_Muon_SR_h
#define HNtypeI_Muon_SR_h

#include "HNAnalyzerCore.h"

class HNtypeI_Muon_SR : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst, RunFake, RunCF;
  bool RunNewPDF;
  bool RunXSecSyst;

  
  vector<TString> MuonTriggers;
  vector<TString> ElectronTriggers;
  vector<TString> EMuTriggers;
  vector<TString> MuonVetoIDs;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronVetoIDs;
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

  HNtypeI_Muon_SR();
  ~HNtypeI_Muon_SR();

};

#endif

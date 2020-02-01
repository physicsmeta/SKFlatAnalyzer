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

  // Trigger
  vector<TString> MuonTriggers;
  vector<TString> MuonTriggersH;
  vector<TString> ElectronTriggers;
  vector<TString> EMuTriggers;
  vector<TString> EMuTriggersH;

  // Lepton ID
  vector<TString> MuonVetoIDs;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronVetoIDs;
  vector<TString> ElectronLooseIDs;
  vector<TString> ElectronTightIDs;

  // Fake rate file
  vector<TString> FakeRateIDs;

  // Lepton pT cut
  double MuonPtCut1;
  double MuonPtCut2;
  double ElectronPtCut1;
  double ElectronPtCut2;
  double EMuPtCut1;
  double EMuPtCut2;

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

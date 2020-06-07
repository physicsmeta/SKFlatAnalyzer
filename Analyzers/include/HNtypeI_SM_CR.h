#ifndef HNtypeI_SM_CR_h
#define HNtypeI_SM_CR_h

#include "HNAnalyzerCore.h"

class HNtypeI_SM_CR : public HNAnalyzerCore {

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
//  vector<TString> EMuTriggers;
//  vector<TString> EMuTriggersH;
  vector<TString> Mu8Ele23Triggers;
  vector<TString> Mu23Ele12Triggers;

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
//  double EMuPtCut1;
//  double EMuPtCut2;

//  vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<FatJet> AllFatJets;

//  double weight_Prefire;

  HNtypeI_SM_CR();
  ~HNtypeI_SM_CR();

};

#endif

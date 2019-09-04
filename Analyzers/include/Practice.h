#ifndef Practice_h
#define Practice_h

#include "AnalyzerCore.h"

class Practice : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEventFromParameter(AnalyzerParameter param, Long64_t Nentry); 
  void executeEvent();
  void executeEvent(Long64_t Nentry); 

  TString EleTriggerName;
  double lep0ptcut, lep1ptcut;

  vector<TString> EleIDs, EleIDSFKeys;
  vector<Electron> AllEles;

  double weight_Prefire;

  Practice();
  ~Practice();

};



#endif


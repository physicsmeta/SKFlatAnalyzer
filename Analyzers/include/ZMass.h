#ifndef ZMass_h
#define ZMass_h

#include "AnalyzerCore.h"

class ZMass : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  TString EleTriggerName;
  double lep0ptcut, lep1ptcut;

  vector<TString> EleIDs, EleIDSFKeys;
  vector<Electron> AllEles;

  double weight_Prefire;

  ZMass();
  ~ZMass();

};



#endif


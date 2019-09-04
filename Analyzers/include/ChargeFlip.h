#ifndef ChargeFlip_h
#define ChargeFlip_h

#include "AnalyzerCore.h"

class ChargeFlip : public AnalyzerCore {

public:

  void initializeAnalyzer();
//  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent(Long64_t Nentry);
	double GetCFweight(std::vector<Electron> eles);
	double GetCFweight_SF(std::vector<Electron> eles);

  TString EleTriggerName;
  double lep0ptcut, lep1ptcut;

  vector<TString> EleIDs, EleIDSFKeys;
  vector<Electron> AllEles;

  ChargeFlip();
  ~ChargeFlip();

};



#endif


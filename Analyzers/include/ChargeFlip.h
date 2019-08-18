#ifndef ChargeFlip_h
#define ChargeFlip_h

#include "AnalyzerCore.h"

class ChargeFlip : public AnalyzerCore {

public:

//  void initializeAnalyzer();
//  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  vector<Electron> AllEles;

  ChargeFlip();
  ~ChargeFlip();

};



#endif


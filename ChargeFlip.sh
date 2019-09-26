#!/bin/bash

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightID
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXY
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdZ
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXYdZ
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags ChargeFlipHE &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 100 --userflags ChargeFlipHE &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 20 --userflags ChargeFlipHE &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags ChargeFlipHE,ClosureTest &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags CFrate,ClosureTest &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 20 --userflags ScaleFactor,RunSyst &
python python/SKFlat.py -a ChargeFlip -y 2017 -i DYJets -n 100 --userflags CFrate,ClosureTest &
python python/SKFlat.py -a ChargeFlip -y 2017 -i DoubleEG -n 20 --userflags ScaleFactor,RunSyst &
python python/SKFlat.py -a ChargeFlip -y 2018 -i DYJets_MG -n 100 --userflags CFrate,ClosureTest &
python python/SKFlat.py -a ChargeFlip -y 2018 -i EGamma -n 20 --userflags ScaleFactor,RunSyst &

#!/bin/bash

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightID
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXY
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdZ
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXYdZ
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags ChargeFlipHE &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 100 --userflags ChargeFlipHE &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 20 --userflags ChargeFlipHE &

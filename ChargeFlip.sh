#!/bin/bash

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passDoubleTrigger
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passSingleTrigger
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passMETFilter,passTightID
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passDoubleTrigger,passTightID
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightID
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXY
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdZ
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passTightChargeTightIDdXYdZ
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags ChargeFlipHE

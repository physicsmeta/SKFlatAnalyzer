#!/bin/bash

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passDoubleTrigger
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passSingleTrigger
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passMETFilter,passSingleTrigger,passTightID
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passMETFilter,passDoubleTrigger,passTightID
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 100 --userflags passMETFilter,passDoubleTrigger

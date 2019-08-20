#!/bin/bash

python python/SKFlat.py -a ZMass -y 2016 -i TTLL_powheg -n 20 
python python/SKFlat.py -a ZMass -y 2016 -i WJets_MG -n 20 
python python/SKFlat.py -a ZMass -y 2016 -i WW_pythia -n 10 
python python/SKFlat.py -a ZMass -y 2016 -i WZ_pythia -n 10 
python python/SKFlat.py -a ZMass -y 2016 -i ZZ_pythia -n 10 
python python/SKFlat.py -a ZMass -y 2016 -i DoubleEG -n 20 
python python/SKFlat.py -a ZMass -y 2016 -i DoubleMuon -n 20 
python python/SKFlat.py -a ZMass -y 2016 -i DYJets -n 100 

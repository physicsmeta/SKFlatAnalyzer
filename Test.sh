#!/bin/bash

python python/SKFlat.py -a test -y 2016 -i DYTypeI_SS_EE_M100 -n 50 --nmax 20 &
python python/SKFlat.py -a test_2016H -y 2016 -i DoubleMuon:H -n 50 --userflags RunOS --nmax 20 &


#python python/SKFlat.py -a test -y 2018 -i DYTypeI_SS_EE_M100 -n 50 &
#python python/SKFlat.py -a test -y 2018 -i DYTypeI_SS_MuMu_M100 -n 50 &
#python python/SKFlat.py -a test -y 2018 -i DYJets -n 50 &
#python python/SKFlat.py -a test -y 2018 -i EGamma -n 50 &
#python python/SKFlat.py -a test -y 2018 -i DoubleMuon -n 50 &
#python python/SKFlat.py -a test -y 2018 -i TTLL_powheg -n 50 &
#python python/SKFlat.py -a test -y 2018 -i WZTo3LNu_powheg -n 50 &
#python python/SKFlat.py -a test -y 2018 -i ZZTo4L_powheg -n 50 &

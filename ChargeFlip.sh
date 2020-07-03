#!/bin/bash

python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --skim SkimTree_Dilepton --userflags CFrate &

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 10 --userflags ScaleFactor &


#CFrate, ClosureTest, HalfSampleTest first run (to calculate the Even set CFrate)

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DYJets -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i TTLL_powheg -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i DYJets_MG -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i TTLL_powheg -n 50 --userflags CFrate &
#
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DYJets -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i TTLL_powheg -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i DYJets_MG -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i TTLL_powheg -n 50 --userflags ClosureTest &
#
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DYJets -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i TTLL_powheg -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i DYJets_MG -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i TTLL_powheg -n 50 --userflags HalfSampleTest &

#HalfSampleTest(second run), ScaleFactor, RunSyst using the above results

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DYJets -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i TTLL_powheg -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i DYJets_MG -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i TTLL_powheg -n 50 --userflags HalfSampleTest &

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 10 --userflags ScaleFactor,RunSyst &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DoubleEG -n 10 --userflags ScaleFactor,RunSyst &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i EGamma -n 10 --userflags ScaleFactor,RunSyst &

#ScaleFactor(second run) using the above results

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DoubleEG -n 10 --userflags ScaleFactor &
#python python/SKFlat.py -a ChargeFlip -y 2017 -i DoubleEG -n 10 --userflags ScaleFactor &
#python python/SKFlat.py -a ChargeFlip -y 2018 -i EGamma -n 10 --userflags ScaleFactor &


###################################### etc ###############################################




#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --userflags CFrate &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --userflags ClosureTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --userflags HalfSampleTest &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --userflags HalfSampleTest &

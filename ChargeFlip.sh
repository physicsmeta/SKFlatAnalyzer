#!/bin/bash

python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-100to200 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-200to400 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-400to500 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-500to700 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-700to800 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-800to1000 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-1000to1500 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-1500to2000 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-2000to3000 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-50To100 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --skim SkimTree_Dilepton --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLJ_powheg -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-70to100 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-100to200 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-200to400 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-400to600 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-600to800 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-800to1200 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-1200to2500 -n 50 --userflags CFrate &
python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-2500toInf -n 50 --userflags CFrate &

#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-100to200 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-200to400 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-400to500 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-500to700 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-700to800 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-800to1000 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-1000to1500 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-1500to2000 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_M-2000to3000 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-50To100 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-100To250 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-250To400 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-400To650 -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_Pt-650ToInf -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i TTLL_powheg -n 50 --skim SkimTree_Dilepton --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-70to100 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-100to200 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-200to400 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-400to600 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-600to800 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-800to1200 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-1200to2500 -n 50 --userflags CFrateValidation &
#python python/SKFlat.py -a ChargeFlip -y 2016 -i DYJets_MG_HT-2500toInf -n 50 --userflags CFrateValidation &

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

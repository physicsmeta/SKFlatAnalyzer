#!/bin/bash
############################################
### SR, CR
############################################

### Test ###
#python python/SKFlat.py -a Signal_2016H -y 2016 -i DoubleMuon:H -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i DoubleEG:B_ver2 -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i WZTo3LNu_powheg -n 80 --skim SkimTree_Dilepton &

### MC ###
#python python/SKFlat.py -a Signal -y 2016 -l submitList/Dilepton_SR_2016.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/NoSkim_SR.txt -n 80 &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/Dilepton_SM_CR_2016.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/NoSkim_SM_CR.txt -n 80 &

#python python/SKFlat.py -a Signal -y 2016 -l submitList/JH_Dilepton_2016.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/JH_NoSkim.txt -n 80 &

#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_Signal.txt -n 50 &

#python python/SKFlat.py -a Control -y 2016 -l submitList/Dilepton_SM_CR_2016.txt -n 50 --skim SkimTree_Dilepton --userflags SM &
#python python/SKFlat.py -a Control -y 2016 -l submitList/NoSkim_SM_CR.txt -n 50 --userflags SM &

### CF ###
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 80 --skim SkimTree_Dilepton --userflags RunCF &

### Fake ###
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal_2016H -y 2016 -l submitList/2016_periodH.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &

### Overlap (Fake & prompt)
#python python/SKFlat.py -a Signal -y 2016 -l submitList/Dilepton_SR_2016.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/NoSkim_SR.txt -n 80 --userflags RunFake &

### DATA ###
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal_2016H -y 2016 -l submitList/2016_periodH.txt -n 80 --skim SkimTree_Dilepton &

#python python/SKFlat.py -a Control -y 2016 -l submitList/2016_DoubleMuon_BtoG.txt -n 50 --skim SkimTree_Dilepton --userflags SM &
#python python/SKFlat.py -a Control_2016H -y 2016 -i DoubleMuon:H -n 50 --skim SkimTree_Dilepton --userflags SM &
#python python/SKFlat.py -a Control -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Control -y 2016 -l submitList/2016_MuonEG_BtoG.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Control_2016H -y 2016 -l submitList/2016_periodH.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Control -y 2016 -l submitList/2016_SingleMuon_BtoH.txt -n 50 --userflags SM &

####TEST####
#python python/SKFlat.py -a Control_2016H -y 2016 -l submitList/2016_DoubleMuon_H.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Control_2016H -y 2016 -l submitList/2016_MuonEG_H.txt -n 50 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Control -y 2016 -i DYJets -n 50 --skim SkimTree_Dilepton &
############

#python python/SKFlat.py -a Signal_ee -y 2016 -l submitList/2016_DoubleEG_BtoH.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal_ee -y 2017 -l submitList/2017_DoubleEG.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal_ee -y 2018 -l submitList/2018_EGamma.txt -n 80 --skim SkimTree_Dilepton &

### Checking signal-like events for using CMSShow
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_cmsshow.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i DoubleEG:B_ver2 -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i MuonEG:B_ver2 -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i DoubleMuon:C -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i DoubleEG:B_ver2 -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2016 -i MuonEG:B_ver2 -n 80 --skim SkimTree_Dilepton &

### Checking event yields using HLT_Mu17_Mu8_SameSign_DZ
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_Fake_DoubleMuon.txt -n 80 &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/2016_Fake_DoubleMuon.txt -n 80 --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2016 -l submitList/Dilepton_SR_2016.txt -n 80 &


### MC ###
#python python/SKFlat.py -a Signal -y 2017 -l submitList/Dilepton_SR_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/NoSkim_SR.txt -n 80 &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/Dilepton_SM_CR_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/NoSkim_SM_CR.txt -n 80 &

#python python/SKFlat.py -a Signal -y 2017 -l submitList/JH_Dilepton_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/JH_NoSkim.txt -n 80 &

### CF ###
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_DoubleEG.txt -n 80 --skim SkimTree_Dilepton --userflags RunCF &

### Fake ###
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_DoubleMuon.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_DoubleEG.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_MuonEG.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &

### DATA ###
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_DoubleMuon.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_DoubleEG.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2017 -l submitList/2017_MuonEG.txt -n 80 --skim SkimTree_Dilepton &



### MC ###
#python python/SKFlat.py -a Signal -y 2018 -l submitList/Dilepton_SR_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/NoSkim_SR.txt -n 80 &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/Dilepton_SM_CR_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/NoSkim_SM_CR.txt -n 80 &

#python python/SKFlat.py -a Signal -y 2018 -l submitList/JH_Dilepton_2017.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/JH_NoSkim.txt -n 80 &

### CF ###
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_EGamma.txt -n 80 --skim SkimTree_Dilepton --userflags RunCF &

### Fake ###
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_DoubleMuon.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_EGamma.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_MuonEG.txt -n 80 --skim SkimTree_Dilepton --userflags RunFake &

### DATA ###
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_DoubleMuon.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_EGamma.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a Signal -y 2018 -l submitList/2018_MuonEG.txt -n 80 --skim SkimTree_Dilepton &


############################################
### Signal acceptance check
############################################

#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M700 -n 50 --userflags M700 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M700 -n 50 --userflags M700,pp &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M700 -n 50 --userflags M700,mm &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMupMup_Tch_M700 -n 50 --userflags M700 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMumMum_Tch_M700 -n 50 --userflags M700 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1000 -n 50 --userflags M1000 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1000 -n 50 --userflags M1000,pp &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1000 -n 50 --userflags M1000,mm &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMupMup_Tch_M1000 -n 50 --userflags M1000 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMumMum_Tch_M1000 -n 50 --userflags M1000 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1500 -n 50 --userflags M1500 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1500 -n 50 --userflags M1500,pp &
#python python/SKFlat.py -a Signal_opt -y 2016 -i HNToMuMu_Tch_M1500 -n 50 --userflags M1500,mm &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMupMup_Tch_M1500 -n 50 --userflags M1500 &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMupMup_Tch_M1500 -n 50 --userflags M1500,pp &
#python python/SKFlat.py -a Signal_opt -y 2016 -i Last_HNToMupMup_Tch_M1500 -n 50 --userflags M1500,mm &


############################################
### Cutflow
############################################

### Checking the dilepton mass cut in the dimuon trigger
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2017 -l submitList/2017_DoubleMuon.txt -n 80 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2017 -l submitList/VV_exclusive.txt -n 80 &

### Checking the dilepton mass cut in the dimuon trigger
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2018 -l submitList/2018_DoubleMuon.txt -n 80 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2018 -l submitList/VV_exclusive.txt -n 80 &

#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i DYJets -n 40 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i TTLL_powheg -n 40 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i WZTo3LNu_powheg -n 20 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i ZZTo4L_powheg -n 30 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i MajoranaNeutrinoToMupMup_Schannel_M100 -n 3 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i MajoranaNeutrinoToMumMum_Schannel_M100 -n 3 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -l submitList/Cutflow_MC_2016.txt -n 80 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i DoubleMuon:B_ver2 -n 80 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -l submitList/2016_DoubleMuon.txt -n 80 &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -l submitList/2016_DoubleMuon.txt -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i DoubleMuon:B_ver2 -n 80 --skim SkimTree_Dilepton &
#python python/SKFlat.py -a HNtypeI_Cutflow -y 2016 -i WZTo3LNu_powheg -n 80 --skim SkimTree_Dilepton &

#################################################
### MC validation year by year at presel level.
#################################################

#python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_SS_EE_M100 -n 50 &
#python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_SS_EE_M700 -n 50 &
#python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_SS_EE_M100 -n 50 &
#python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_SS_EE_M700 -n 50 &
#python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_SS_MuMu_M100 -n 50 &
#python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_SS_MuMu_M700 -n 50 &
#python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_SS_MuMu_M100 -n 50 &
#python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_SS_MuMu_M700 -n 50 &
python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_OS_EE_M100 -n 50 &
python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_OS_EE_M700 -n 50 &
python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_OS_EE_M100 -n 50 &
python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_OS_EE_M700 -n 50 &
python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_OS_MuMu_M100 -n 50 &
python python/SKFlat.py -a Presel -y 2017 -i DYTypeI_OS_MuMu_M700 -n 50 &
python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_OS_MuMu_M100 -n 50 &
python python/SKFlat.py -a Presel -y 2018 -i DYTypeI_OS_MuMu_M700 -n 50 &

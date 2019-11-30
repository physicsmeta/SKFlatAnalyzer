import commands as cmd
import os

os.chdir('/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0000')
ls1 = cmd.getoutput('ls')
os.chdir('/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0001')
ls2 = cmd.getoutput('ls')
os.chdir('/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0002')
ls3 = cmd.getoutput('ls')

os.chdir('/data4/Users/jihkim/SKFlatAnalyzer/data/Run2Legacy_v3/2018/Sample/ForSNU')

list_ls1 = ls1.split('\n')
list_ls2 = ls2.split('\n')
list_ls3 = ls3.split('\n')

for i in range(len(list_ls1)):
  list_ls1[i] = '/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0000/'+list_ls1[i]
for i in range(len(list_ls2)):
  list_ls2[i] = '/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0001/'+list_ls2[i]
for i in range(len(list_ls3)):
  list_ls3[i] = '/gv0/DATA/SKFlat/Run2Legacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/190930_135722/0002/'+list_ls3[i]

with open('DYJets.txt','w') as f:
  for i in range(len(list_ls1)):
    f.write(list_ls1[i]+'\n')
  for i in range(len(list_ls2)):
    f.write(list_ls2[i]+'\n')
  for i in range(len(list_ls3)):
    f.write(list_ls3[i]+'\n')

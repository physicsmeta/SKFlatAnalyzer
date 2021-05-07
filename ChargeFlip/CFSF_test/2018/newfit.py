import ROOT as rt
rt.gROOT.LoadMacro('./histFitter.C+')
#rt.gROOT.LoadMacro('./RooCBExGaussShape.cc+')
#rt.gROOT.LoadMacro('./RooCMSShape.cc+')
from ROOT import tnpFitter
import time

channels = ["BB","EE"]

for channel in channels:

  #fileTruth  = rt.TFile("CFSF_test_DYJets_MG_All.root",'read')
  #fileTruth  = rt.TFile("CFSF_test_DYJets_MG.root",'read')
  #fileTruth  = rt.TFile("CFSF_test_DYJets.root",'read')
  fileTruth  = rt.TFile("CFSF_test_DYJets_All.root",'read')

  funcs = [
      "Gaussian::sigResPass(x,meanOS,sigmaOS)",
      "Gaussian::sigResFail(x,meanSS,sigmaSS)",
      "RooCMSShape::bkgPass(x, acmsOS, betaOS, gammaOS, peakOS)",
      "RooCMSShape::bkgFail(x, acmsSS, betaSS, gammaSS, peakSS)",
      ]
 
  if channel == "BB" and fileTruth.GetName() == 'CFSF_test_DYJets_All.root':
    pars = [
        "meanOS[-0.0,-5.0,5.0]","sigmaOS[0.9,0.5,5.0]",
        "meanSS[-0.0,-5.0,5.0]","sigmaSS[0.9,0.5,5.0]",
        "acmsOS[60.,50.,80.]","betaOS[0.05,0.01,0.08]","gammaOS[0.1, -2, 2]","peakOS[90.0]",
        "acmsSS[60.,50.,80.]","betaSS[0.05,0.01,0.08]","gammaSS[0., -2, 0.015]","peakSS[90.0]",
        ]
  else:
    pars = [
        "meanOS[-0.0,-5.0,5.0]","sigmaOS[0.9,0.5,5.0]",
        "meanSS[-0.0,-5.0,5.0]","sigmaSS[0.9,0.5,5.0]",
        "acmsOS[60.,50.,80.]","betaOS[0.05,0.01,0.08]","gammaOS[0.1, -2, 2]","peakOS[90.0]",
        "acmsSS[60.,50.,80.]","betaSS[0.05,0.01,0.08]","gammaSS[0.1, -2, 2]","peakSS[90.0]",
        ]

  
  this_workspace = []
  this_workspace.extend(pars)
  this_workspace.extend(funcs)
  
  infile = rt.TFile("CFSF_test_EGamma.root", "read")
  if channel == "BE":
    hOS = infile.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_OS_CFSFweighted_shifted_0.8%")
  else:
    hOS = infile.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_OS_CFweighted_shifted_0.8%")
  hSS = infile.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_SS_MET0")
  fitter = tnpFitter( hOS, hSS, "myhist_"+channel )
  infile.Close()
  
  fitter.useMinos()
  rootfile = rt.TFile("newfit.root",'update')
  fitter.setOutputFile( rootfile )
  
  if channel == "BE":
    histZLineShapeOS = fileTruth.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_OS_CFSFweighted_shifted_0.8%")
  else:
    histZLineShapeOS = fileTruth.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_OS_CFweighted_shifted_0.8%")
  histZLineShapeSS = fileTruth.Get("HNTightV1/ScaleFactor/"+channel+"_ZMass_SS_MET0")
  fitter.setZLineShapes(histZLineShapeOS,histZLineShapeSS)
  
  fileTruth.Close()
  
  workspace = rt.vector("string")()
  for i in this_workspace:
    workspace.push_back(i)
  fitter.setWorkspace( workspace )
  
  title = "mytitle_"+channel
  fitter.fits(False,title)
  print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  time.sleep(2)
  rootfile.Close()

print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"

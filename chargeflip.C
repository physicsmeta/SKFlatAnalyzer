R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatAnalyzer/lib/libAnalyzerTools.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib/libLHAPDF.so)


void chargeflip(){

  ChargeFlip m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 1000;
  m.MCSample = "DYJets";
  //m.MCSample = "DYJets_Pt-250To400";
  //m.MCSample = "DYJets_Pt-650ToInf";
  m.IsDATA = false;
  m.xsec = 6225.42;
  m.sumW = 80924255;
  //m.xsec = 3.047;
  //m.sumW = 7781193;
  //m.xsec = 0.03636;
  //m.sumW = 422153;
  m.IsFastSim = false;
  m.DataYear = 2016;
  m.Userflags = {
    "CFrate",
    //"bar",
  };
  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_1.root");
  //m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190428_140356/0000/SKFlatNtuple_2016_MC_1.root");
  //m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190428_143538/0000/SKFlatNtuple_2016_MC_1.root");
  m.SetOutfilePath("chargeflip.root");
  //m.SetOutfilePath("chargeflip_DYJets_Pt-250To400.root");
  //m.SetOutfilePath("chargeflip_DYJets_Pt-650ToInf.root");
  m.Init();
  m.initializeAnalyzerTools();
  m.initializeAnalyzer();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}

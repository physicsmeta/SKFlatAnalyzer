R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(/data4/Users/jihkim/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/data4/Users/jihkim/SKFlatAnalyzer/lib/libAnalyzerTools.so)
R__LOAD_LIBRARY(/data4/Users/jihkim/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib/libLHAPDF.so)


void practice(){

  Practice m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 1000;
  m.MCSample = "DYJets";
  m.IsDATA = false;
  m.xsec = 6225.42;
  m.sumW = 80924255;
  m.IsFastSim = false;
  m.DataYear = 2016;
  m.Userflags = {
    "A",
    "B",
  };
  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_1.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_10.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_100.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_101.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_102.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_103.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_104.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_105.root");
//  m.AddFile("/data7/DATA/SKFlat/Run2Legacy_v3/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/190423_002438/0000/SKFlatNtuple_2016_MC_106.root");
  m.SetOutfilePath("practice.root");
  m.Init();
  m.initializeAnalyzerTools();
  m.initializeAnalyzer();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}

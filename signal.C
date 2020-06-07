R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatRunlog/2020_03_14_174250__70143__HNtypeI_SR__Year2016__TAMSA1/lib/libDataFormats.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatRunlog/2020_03_14_174250__70143__HNtypeI_SR__Year2016__TAMSA1/lib/libAnalyzerTools.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatRunlog/2020_03_14_174250__70143__HNtypeI_SR__Year2016__TAMSA1/lib/libGEScaleSyst.so)
R__LOAD_LIBRARY(/data6/Users/jihkim/SKFlatRunlog/2020_03_14_174250__70143__HNtypeI_SR__Year2016__TAMSA1/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib/libLHAPDF.so)


void run_0(){

  Signal m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 1000;
  m.MCSample = "DYJets";
  m.IsDATA = false;
  m.xsec = 6077.22;
  m.sumW = 80924255;
  m.IsFastSim = false;
  m.DataYear = 2016;
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_1.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_10.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_100.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_101.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_102.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_103.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_104.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_105.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_106.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_107.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_108.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_109.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_11.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_110.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_111.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_112.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_113.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/191229_214240/0000/SKFlatNtuple_2016_MC_114.root");
  m.SetOutfilePath("hists.root");
  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}

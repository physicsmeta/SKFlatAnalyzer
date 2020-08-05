#include "HNtypeI_SM_CR.h"

HNtypeI_SM_CR::HNtypeI_SM_CR(){

}

void HNtypeI_SM_CR::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  cout << "[HNtypeI_SM_CR::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  RunFake = HasFlag("RunFake");
  RunCF = HasFlag("RunCF");

  /*MuonTightIDs = {"HNTight", "HNTightV2", "HNTight2016"};
  MuonLooseIDs = {"HNLoose", "HNLoose", "HNLoose2016"};
  MuonVetoIDs  = {"POGLoose", "POGLoose", "HNVeto2016"};
  //MuonIDSFKeys = { "NUM_TightID_DEN_genTracks" };
  ElectronTightIDs = {"HNTight", "HNTightV2", "HNTight2016"};
  ElectronLooseIDs = {"HNLoose", "HNLooseV23", "HNLoose2016"};
  ElectronVetoIDs  = {"passVetoID", "passVetoID", "HNVeto2016"};
  FakeRateIDs = {"HNtypeI_V1", "HNtypeI_V2", "HNtypeI_16"};*/
  MuonTightIDs     = {"HNTight2016", "POGHighPtWithLooseTrkIso", "HNTightV1", "HNTightV2", "HNTightV1"};
  MuonLooseIDs     = {"HNLoose2016", "POGMedium", "HNLooseV1", "HNLooseV1", "HNLooseV2"};
  MuonVetoIDs      = {"HNVeto2016", "POGLoose", "HNVeto", "HNVeto", "HNVeto"};
  ElectronTightIDs = {"HNTight2016", "HEEP_dZ", "HNTightV1", "HNTightV2", "HNMVATight"};
  ElectronLooseIDs = {"HNLoose2016", "passLooseID", "HNLooseV1", "HNLooseV1", "HNMVALoose"};
  ElectronVetoIDs  = {"HNVeto2016", "passVetoID", "HNVeto", "HNVeto", "HNMVAVeto"};
  FakeRateIDs      = {"HNtypeI_16", "HNtypeI_16", "HNtypeI_16", "HNtypeI_16", "HNtypeI_16"};


  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SM_CR.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  MuonTriggers.clear();
  MuonTriggersH.clear();
  ElectronTriggers.clear();
  EMuTriggers.clear();
  EMuTriggersH.clear();

  if(DataYear==2016){                                                                   // Lumi values for trigger weight (/pb)
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                       // 27267.591112919
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                     // 27267.591112919
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                    // 35918.219492947
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                  // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                   // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                 // 35918.219492947
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 27267.591112919
    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
    EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2018){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

  //cout << "[HNtypeI_SM_CR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_SM_CR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_SM_CR::~HNtypeI_SM_CR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_SM_CR::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_SM_CR.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets());

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SM_CR.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){
    TString MuonTightID = MuonTightIDs.at(it_id);
    TString MuonLooseID = MuonLooseIDs.at(it_id); 
    TString MuonVetoID  = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);
    TString FakeRateID = FakeRateIDs.at(it_id);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;

    //param.Name = MuonID+"_"+"Central";

    //==== Muon ID
    param.Muon_Tight_ID       = MuonTightID;
    param.Muon_Loose_ID       = MuonLooseID;
    param.Muon_Veto_ID        = MuonVetoID;
    param.Muon_FR_ID          = FakeRateID;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key         = "AwayJetPt40";  // histname
    //param.Muon_ID_SF_Key      = "NUM_TightID_DEN_genTracks";
    //param.Muon_ISO_SF_Key     = "NUM_TightRelIso_DEN_TightIDandIPCut";
    param.Muon_ID_SF_Key      = "";
    param.Muon_ISO_SF_Key     = "";
    param.Muon_Trigger_SF_Key = "";
    param.Muon_UsePtCone      = true;

    //==== Electron ID
    param.Electron_Tight_ID       = ElectronTightID;
    if(DataYear==2018 && it_id==1) param.Electron_Tight_ID = "HEEP2018_dZ";
    param.Electron_Loose_ID       = ElectronLooseID;
    param.Electron_Veto_ID        = ElectronVetoID;
    param.Electron_FR_ID          = FakeRateID;     // ID name in histmap_Electron.txt
    param.Electron_FR_Key         = "AwayJetPt40";  // histname
    //param.Electron_ID_SF_Key      = "passTightID";
    param.Electron_ID_SF_Key      = "";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone      = true;

    //==== Jet ID
    //param.Jet_ID = "tightLepVeto";
    param.Jet_ID    = "HNTight";
    param.FatJet_ID = "HNTight";

    if(RunFake){
      if(it_id > 0) continue;
    }

    if(it_id != 1) continue;

    executeEventFromParameter(param);

    /*if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }*/
  }
}

void HNtypeI_SM_CR::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> regions = {"DYmm", "DYee", "DYemu", "TTmm", "TTee", "TTemu", "WZ", "ZG", "WG", "ZZ"}; 
  vector<TString> channels3L = {"mmm", "mme", "mee", "eee"};
  vector<TString> channels4L = {"mmmm", "mmee", "eeee"};
  TString IDsuffix = "HN16";
  if(param.Electron_Tight_ID.Contains("HEEP")) IDsuffix = "LRSM";
  if(param.Electron_Tight_ID.Contains("V1")) IDsuffix = "HNV1";
  if(param.Electron_Tight_ID.Contains("V2")) IDsuffix = "HNV2";
  if(param.Electron_Tight_ID.Contains("MVA")) IDsuffix = "MVA";
  TString NjetBin = "", VtxBin = "";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1., weight_noPU = 1., weight_vertex = 1., weight_rho = 1., PUweight_up = 1., PUweight_down = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;
 
  Event ev = GetEvent();

  //========================================================
  //==== No Cut
  //========================================================

  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);

    weight_noPU *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight_noPU *= ev.MCweight();
    weight_noPU *= GetPrefireWeight(0);
  }

  int Nvtx = nPV;
  //if(!IsDATA) Nvtx = nPileUp+1;

  FillHist("nPileUp", nPileUp, weight, 200., 0., 200.);
  //FillHist("nPV", nPV, weight, 200., 0., 200.);
  FillHist("Nvtx", Nvtx, weight, 200., 0., 200.);

  FillHist("nPileUp_noPU", nPileUp, weight_noPU, 200., 0., 200.);
  //FillHist("nPV_noPU", nPV, weight_noPU, 200., 0., 200.);
  FillHist("Nvtx_noPU", Nvtx, weight_noPU, 200., 0., 200.);

  // Cutflow : No Cuts
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, -0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, -0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  // Cutflow : MET filter
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== Trigger
  //========================================================

  //if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers))) return; 
  if(DataYear==2016){
    if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers) || ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersH))) return;
  }
  else{
    if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers) || ev.PassTrigger(EMuTriggers))) return;
  }

  //========================================================
  //==== Copy AllObjects
  //========================================================

  vector<Muon> this_AllMuons;
  if(param.Muon_Tight_ID.Contains("HighPt")) this_AllMuons = UseTunePMuon(AllMuons);
  else this_AllMuons = AllMuons;
  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Jet> this_AllJets = AllJets;
  vector<FatJet> this_AllFatJets = AllFatJets;
  vector<Gen> gens = GetGens();

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  /*if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnUp){
    this_AllMuons = ScaleMuons( this_AllMuons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnDown){
    this_AllMuons = ScaleMuons( this_AllMuons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResUp){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResDown){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnUp){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnDown){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
  }
  else{
    cout << "[HNtypeI_SM_CR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  // Leptons
  TString MuonID = param.Muon_Tight_ID;
  TString ElectronID = param.Electron_Tight_ID;
  if(RunFake){
    MuonID = param.Muon_Loose_ID;
    ElectronID = param.Electron_Loose_ID;
  }

  vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 10., 2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);

  // Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  // Charge flip
  vector<Electron> electrons_beforeShift;
  vector<Electron> electrons_afterShift;
  electrons_beforeShift.clear();
  electrons_afterShift.clear();

  // Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 10., 2.7);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  vector<Jet> jets = JetsPassPileupMVA(jets_lepveto);

  vector<Jet> jets_Pt10to18;
  vector<Jet> jets_Pt20;
  vector<Jet> jets_Pt30;
  jets_Pt10to18.clear();
  jets_Pt20.clear();
  jets_Pt30.clear();

  for(unsigned int i=0; i<jets_lepveto.size(); i++){
    if(jets_lepveto.at(i).Pt()>10. && jets_lepveto.at(i).Pt()<18.) jets_Pt10to18.push_back(jets_lepveto.at(i));
  }

  for(unsigned int i=0; i<jets.size(); i++){
    if(jets.at(i).Pt() > 20.) jets_Pt20.push_back(jets.at(i));
    if(jets.at(i).Pt() > 30.) jets_Pt30.push_back(jets.at(i));
  }
  //vector<Jet> jets = JetsAwayFromFatJet(jets_PUveto, fatjets);  // AK4jets used in SR, CR
  

  std::vector<Lepton*> leptons, leptons_minus, leptons_plus, leptons_veto;

  //========================================================
  //==== Sort in pT-order
  //========================================================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);
  std::sort(jets_Pt10to18.begin(), jets_Pt10to18.end(), PtComparing);
  std::sort(jets_Pt20.begin(), jets_Pt20.end(), PtComparing);
  std::sort(jets_Pt30.begin(), jets_Pt30.end(), PtComparing);

  //========================================================
  //==== B-Tagging
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0, Nbjet_Pt30_loose = 0., Nbjet_Pt30_medium = 0, Nbjet_medium_PUveto = 0;
  JetTagging::Parameters jtp_DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
    //double this_discr = jets_nolepveto.at(ij).GetTaggerResult(JetTagging::DeepCSV);
    //==== No SF
    //if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
    if(jets_nolepveto.at(ij).Pt() > 20.){
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
    }
    if(jets_nolepveto.at(ij).Pt() > 30.){
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_Pt30_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_Pt30_medium++;
    }
  }

  for(unsigned int ij=0; ij<jets_Pt20.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_Pt20.at(ij))) Nbjet_medium_PUveto++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  //ev.SetMET(pfMET_Type1_PhiCor_pt, pfMET_Type1_PhiCor_phi);
  //Particle METv_dxy = ev.GetMETVector();

  if((muons.size()+electrons.size() > 1) && (muons.size()+electrons.size() < 5)){
    METv = UpdateMETMuon(METv, muons);
    METv = UpdateMETElectron(METv, electrons);
    //METv_dxy = UpdateMETMuon(METv_dxy, muons);
    //METv_dxy = UpdateMETElectron(METv_dxy, electrons);
  }

  double MET = METv.Pt();
  //double MET_dxy = METv_dxy.Pt();

  //========================================================
  //==== Define particles, variables
  //========================================================

  double HT = 0., ST = 0., MET2ST = 0.;
  double Mt = 0., Mt3l = 0.;
  double MZ = 91.1876;
  //double MW = 80.379;
  double muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_recosf = 1., ele_idsf = 1.;
  int lepton_veto_size = 0;

  bool passPtCut = false;
  Particle ZCand, Wtemp1, Wtemp2, WCand1, WCand2;
  Particle llj, l1j, l2j,  lljj, l1jj, l2jj, l1J, l2J;
  Particle WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp, Ztemp1, Ztemp2, Ztemp3, Ztemp4, ZCand1, ZCand2, GammaCand, GammaLep1, GammaLep2;
  int ossf_mass10 = 0;
  
  // Set up pTcone when RunFake=true
  double mu_tight_iso = 0.15, el_tight_iso = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(RunFake){
    //if(IDsuffix == "HNV2") mu_tight_iso = 0.1;
    if(IDsuffix == "HN16") mu_tight_iso = 0.07;

    for(unsigned int i=0; i<muons.size(); i++){
      this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso);
      muons.at(i).SetPtCone(this_ptcone_muon);
    }

    for(unsigned int i=0; i<electrons.size(); i++){
      /*el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt();
      if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
      if(IDsuffix == "HNV2"){
        el_tight_iso = std::min(0.08, 0.0287+0.506/electrons.at(i).UncorrPt());
        if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/electrons.at(i).UncorrPt());
      }*/
      if(IDsuffix == "HN16") el_tight_iso = 0.08;
      this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
      electrons.at(i).SetPtCone(this_ptcone_electron);
    }

    // Correct MET when RunFake=true, because pT was replaced by pTcone
    if((muons.size()+electrons.size() > 1) && (muons.size()+electrons.size() < 5)){
      METv = UpdateMETFake(METv, electrons, muons);
      //METv_dxy = UpdateMETFake(METv_dxy, electrons, muons);
    }

    muons = MuonUsePtCone(muons);
    electrons = ElectronUsePtCone(electrons);
    std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
  }

  /*if(muons.size()==2 && electrons.size()==0){
    FillHist("Pt_muon1", muons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_muon2", muons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon1", muons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon2", muons.at(1).PtCone(), weight, 1000, 0., 1000.);
  }
  if(muons.size()==0 && electrons.size()==2){
    FillHist("Pt_electron1", electrons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_electron2", electrons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron1", electrons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron2", electrons.at(1).PtCone(), weight, 1000, 0., 1000.);
  }

  if(electrons.size() > 0) cout << electrons.at(0).PtCone() << endl;*/

  // Define leptons (pT order)
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));
  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  // Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  // leptons (minus, plus charge)
  for(unsigned int i=0; i<muons.size(); i++){
    if(muons.at(i).Charge() < 0) leptons_minus.push_back(&muons.at(i));
    if(muons.at(i).Charge() > 0) leptons_plus.push_back(&muons.at(i));
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    if(electrons.at(i).Charge() < 0) leptons_minus.push_back(&electrons.at(i));
    if(electrons.at(i).Charge() > 0) leptons_plus.push_back(&electrons.at(i));
  }

  lepton_veto_size = leptons_veto.size() - leptons.size();

  // Define HT, ST, MET^2/ST
  MET = METv.Pt();
  //MET_dxy = METv_dxy.Pt();

  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;

  for(unsigned int i=0; i<jets_Pt10to18.size(); i++){
    HT += jets_Pt10to18.at(i).Pt();
  }

  //=====================================================================================
  //==== SM background CR (DYmm, DYee, DYemu, WZ, ZG, WG, ZZ)
  //=====================================================================================

  if(RunCF) return;

  // Period-dependent trigger weight (only for 2016 MC)
  trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
  if(!IsDATA){
    if(DataYear==2016){
      if(ev.PassTrigger(MuonTriggers)){
        if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
        if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
        trigger_lumi = dimu_trig_weight;
      }
      if(ev.PassTrigger(ElectronTriggers)){
        trigger_lumi = ev.GetTriggerLumi("Full"); 
      }
      if(ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersH)){
        if(ev.PassTrigger(EMuTriggers)) emu_trig_weight += 27267.591;
        if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628;
        trigger_lumi = emu_trig_weight;
      }
    }
    else{
      trigger_lumi = ev.GetTriggerLumi("Full");
    }
  }

  //========================================================
  //==== Event selection
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_recosf = 1., ele_idsf = 1.;
    ossf_mass10 = 0;

    if(!IsDATA){
      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);
    }
    // Cutflow : passing dilepton triggers (dimu || diel || emu)
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }

    if(param.Muon_Tight_ID.Contains("HighPt")){
      MuonPtCut1 = 50., MuonPtCut2 = 50.;
      ElectronPtCut1 = 35., ElectronPtCut2 = 35.;
    }
    else{
      MuonPtCut1 = 20., MuonPtCut2 = 10.;
      ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    }

    //========================================================
    //==== DY, ttbar control region
    //========================================================

    if(it_rg < 6){

      // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)

      if(leptons.size() != 2) continue;
      //if(muons.size()==1 && electrons.size()==1) continue;

      // Passing triggers
      if(it_rg==0 || it_rg==3){ if(!ev.PassTrigger(MuonTriggers)) continue; }
      if(it_rg==1 || it_rg==4){ if(!ev.PassTrigger(ElectronTriggers)) continue; }
      if(it_rg==2 || it_rg==5){ if(!ev.PassTrigger(EMuTriggers)) continue; }

      trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
      if(!IsDATA){
        if(DataYear==2016){
          if(it_rg==0 || it_rg==3){
            if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
            if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
            trigger_lumi = dimu_trig_weight;
          }
          if(it_rg==1 || it_rg==4) trigger_lumi = ev.GetTriggerLumi("Full");
          if(it_rg==2 || it_rg==5){
            if(ev.PassTrigger(EMuTriggers)) emu_trig_weight += 27267.591;
            if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628;
            trigger_lumi = emu_trig_weight;
          }
        }
        else{
          trigger_lumi = ev.GetTriggerLumi("Full");
        }
      } 

      // pT > trigger thresholds
      passPtCut = false;

      if(it_rg==0 || it_rg==3){
        if(!(muons.size()==2 && electrons.size()==0)) continue;
        if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      }
      if(it_rg==1 || it_rg==4){
        if(!(muons.size()==0 && electrons.size()==2)) continue;
        if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      }
      if(it_rg==2 || it_rg==5){
        if(!(muons.size()==1 && electrons.size()==1)) continue;
        if(ev.PassTrigger(Mu8Ele23Triggers)){
          if(electrons.at(0).Pt()>EMuPtCut1 && muons.at(0).Pt()>EMuPtCut2) passPtCut = true;
        }
        if(ev.PassTrigger(Mu23Ele12Triggers)){
          if(muons.at(0).Pt()>EMuPtCut1 && electrons.at(0).Pt()>EMuPtCut2) passPtCut = true;;
        }
      }

      if(!passPtCut) continue;

      // Truth matching
      muons_prompt.clear();
      electrons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

      if(it_rg==0 || it_rg==3){
        if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) continue;
      }
      if(it_rg==1 || it_rg==4){
        if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) continue;
      }
      if(it_rg==2 || it_rg==5){
        if(!(muons_prompt.size()==1 && electrons_prompt.size()==1)) continue;
      }

      weight = 1., weight_noPU = 1., weight_vertex = 1., weight_rho = 1., PUweight_up = 1., PUweight_down = 1.;

      if(it_rg < 2){
        weight_vertex = GetVertexWeight(Nvtx, regions.at(it_rg));
        weight_rho = GetRhoWeight(Rho, regions.at(it_rg));
      }

      // weights for MC
      if(!IsDATA){

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        weight_noPU *= weight_norm_1invpb*trigger_lumi;
        weight_noPU *= ev.MCweight();
        weight_noPU *= GetPrefireWeight(0);

        PUweight_up = GetPileUpWeight(nPileUp,1);
        PUweight_down = GetPileUpWeight(nPileUp,-1);

        for(unsigned int i=0; i<muons.size(); i++){
          if(param.Muon_Tight_ID.Contains("HighPt")){
            muon_miniaodP = sqrt( muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz() );
            muon_recosf   = mcCorr->MuonReco_SF("HighPtMuonRecoSF", muons.at(i).Eta(), muon_miniaodP, 0);
            muon_idsf     = mcCorr->MuonID_SF("NUM_HighPtID_DEN_genTracks",  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
            muon_isosf    = mcCorr->MuonISO_SF("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut", muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          }
          else{
            muon_recosf = 1.;
            muon_idsf   = 1.;
            muon_isosf  = 1.;
          }
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int i=0; i<electrons.size(); i++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }

      }
      if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

      ZCand = *leptons.at(0) + *leptons.at(1);

      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      // Cutflow : OS event
      if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) continue;

      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      // Cutflow : veto 3rd leptons using veto ID
      if(lepton_veto_size > 0) continue;

      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      // Cutflow : m(ll) > 50 GeV (Because there is no DY10to50 sample)
      if(!(ZCand.M() > 50.)) continue;

      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);      

      if(it_rg < 3){

        // Cutflow : No b jets
        if(!(Nbjet_medium == 0)) continue;

        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_Pt20_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_Jets_Pt30_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        //FillHist(regions.at(it_rg)+"/MET_dxy_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

        FillHist(regions.at(it_rg)+"/Number_Vertex_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
        FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
        FillHist(regions.at(it_rg)+"/Rho_"+IDsuffix, Rho, weight, 400, 0., 100.);
        FillHist(regions.at(it_rg)+"/HT_"+IDsuffix, HT, weight, 1000, 0., 1000.);

        if(ZCand.M()>80. && ZCand.M()<100.){
          FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Jets_Pt20_Mass80to100_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_Jets_Pt30_Mass80to100_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Loose_Mass80to100_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Medium_Mass80to100_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_FatJets_Mass80to100_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/ZCand_Mass_Mass80to100_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(regions.at(it_rg)+"/ZCand_Pt_Mass80to100_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Pt_Mass80to100_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep2_Pt_Mass80to100_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Eta_Mass80to100_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/Lep2_Eta_Mass80to100_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);

          FillHist(regions.at(it_rg)+"/MET_Mass80to100_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          //FillHist(regions.at(it_rg)+"/MET_dxy_Mass80to100_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET_PUup_Mass80to100_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET_PUdown_Mass80to100_PUdown_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET_vertex_Mass80to100_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET_rho_Mass80to100_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET2ST_Mass80to100_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          FillHist(regions.at(it_rg)+"/Number_Vertex_Mass80to100_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_Mass80to100_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_PUup_Mass80to100_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_PUdown_Mass80to100_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_vertex_Mass80to100_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_rho_Mass80to100_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
 
          FillHist(regions.at(it_rg)+"/Rho_Mass80to100_"+IDsuffix, Rho, weight, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/Rho_noPU_Mass80to100_"+IDsuffix, Rho, weight_noPU, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/Rho_PUup_Mass80to100_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/Rho_PUdown_Mass80to100_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/Rho_vertex_Mass80to100_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/Rho_rho_Mass80to100_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
          
          FillHist(regions.at(it_rg)+"/HT_Mass80to100_"+IDsuffix, HT, weight, 1000, 0., 1000.);

          FillHist(regions.at(it_rg)+"/NvtxRho2D_Mass80to100_"+IDsuffix, Nvtx, Rho, weight, 200, 0., 200., 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/NvtxRho2D_PUup_Mass80to100_"+IDsuffix, Nvtx, Rho, weight_noPU*PUweight_up, 200, 0., 200., 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/NvtxRho2D_PUdown_Mass80to100_"+IDsuffix, Nvtx, Rho, weight_noPU*PUweight_down, 200, 0., 200., 400, 0., 100.);

          if(param.Electron_Tight_ID.Contains("HEEP")){
          
            /*if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/Jet1_Pt20_Pt_Mass80to100_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
            if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/Jet2_Pt20_Pt_Mass80to100_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
            if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/Jet1_Pt30_Pt_Mass80to100_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
            if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/Jet2_Pt30_Pt_Mass80to100_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/


            //==== PV bin
            if(Nvtx <= 10){
              VtxBin = "Vtx1to10";
  
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt20_"+VtxBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt30_"+VtxBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUup_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUdown_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_vertex_"+VtxBin+"_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_rho_"+VtxBin+"_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/VertexBin/MET_dxy_"+VtxBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_noPU_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUup_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUdown_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_rho_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_"+VtxBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUup_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUdown_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_vertex_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_rho_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/HT_"+VtxBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
  
              /*if(jets_Pt20.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt20_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt30_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/
            }
            else if(Nvtx <= 20){
              VtxBin = "Vtx11to20";
  
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt20_"+VtxBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt30_"+VtxBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUup_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUdown_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_vertex_"+VtxBin+"_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_rho_"+VtxBin+"_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/VertexBin/MET_dxy_"+VtxBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_noPU_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUup_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUdown_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_rho_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_"+VtxBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUup_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUdown_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_vertex_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_rho_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/HT_"+VtxBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
  
              /*if(jets_Pt20.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt20_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt30_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/
            }
            else if(Nvtx <= 30){
              VtxBin = "Vtx21to30";
  
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt20_"+VtxBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt30_"+VtxBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUup_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUdown_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_vertex_"+VtxBin+"_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_rho_"+VtxBin+"_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/VertexBin/MET_dxy_"+VtxBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_noPU_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUup_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUdown_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_rho_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_"+VtxBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUup_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUdown_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_vertex_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_rho_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/HT_"+VtxBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
  
              /*if(jets_Pt20.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt20_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt30_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/
            }
            else if(Nvtx <= 40){
              VtxBin = "Vtx31to40";
  
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt20_"+VtxBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt30_"+VtxBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUup_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUdown_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_vertex_"+VtxBin+"_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_rho_"+VtxBin+"_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/VertexBin/MET_dxy_"+VtxBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_noPU_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUup_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUdown_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_rho_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_"+VtxBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUup_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUdown_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_vertex_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_rho_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/HT_"+VtxBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
  
              /*if(jets_Pt20.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt20_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt30_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/
            }
            else{
              VtxBin = "Vtx41toInf";
  
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt20_"+VtxBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Jets_Pt30_"+VtxBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUup_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_up, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_PUdown_"+VtxBin+"_"+IDsuffix, MET, weight_noPU*PUweight_down, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_vertex_"+VtxBin+"_"+IDsuffix, MET, weight*weight_vertex, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/MET_rho_"+VtxBin+"_"+IDsuffix, MET, weight*weight_rho, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/VertexBin/MET_dxy_"+VtxBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_noPU_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUup_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_up, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_PUdown_"+VtxBin+"_"+IDsuffix, Nvtx, weight_noPU*PUweight_down, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_vertex_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_vertex, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Number_Vertex_rho_"+VtxBin+"_"+IDsuffix, Nvtx, weight*weight_rho, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_"+VtxBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUup_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_up, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_PUdown_"+VtxBin+"_"+IDsuffix, Rho, weight_noPU*PUweight_down, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_vertex_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_vertex, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/Rho_rho_"+VtxBin+"_"+IDsuffix, Rho, weight*weight_rho, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/VertexBin/HT_"+VtxBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
  
              /*if(jets_Pt20.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt20_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt20.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt20_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() == 0) FillHist(regions.at(it_rg)+"/VertexBin/MET_Njet0_Pt30_"+VtxBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 0) FillHist(regions.at(it_rg)+"/VertexBin/Jet1_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              if(jets_Pt30.size() > 1) FillHist(regions.at(it_rg)+"/VertexBin/Jet2_Pt30_Pt_"+VtxBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);*/
            }
  
  
            //==== Njet bin when pT(j) > 20 GeV
            if(jets_Pt20.size() == 0){
              NjetBin = "Njet0";
  
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Jets_Pt20_"+NjetBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
            }
            if(jets_Pt20.size() == 1){
              NjetBin = "Njet1";
  
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Jets_Pt20_"+NjetBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet1_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
            }
            if(jets_Pt20.size() == 2){
              NjetBin = "Njet2";
  
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Jets_Pt20_"+NjetBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet1_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet2_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
            }
            if(jets_Pt20.size() == 3){
              NjetBin = "Njet3";
  
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Jets_Pt20_"+NjetBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet1_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet2_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
            }
            if(jets_Pt20.size() >= 4){
              NjetBin = "Njet4";
  
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Jets_Pt20_"+NjetBin+"_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt20/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet1_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt20/Jet2_Pt20_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt20.at(1).Pt(), weight, 1000, 0., 1000.);
            }
  
  
            //===== Njet bin when pT(j) > 30 GeV
            if(jets_Pt30.size() == 0){
              NjetBin = "Njet0";
              
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Jets_Pt30_"+NjetBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
            } 
            if(jets_Pt30.size() == 1){
              NjetBin = "Njet1";
              
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Jets_Pt30_"+NjetBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet1_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
            } 
            if(jets_Pt30.size() == 2){
              NjetBin = "Njet2";
    
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Jets_Pt30_"+NjetBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet1_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet2_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);
            } 
            if(jets_Pt30.size() == 3){
              NjetBin = "Njet3";
            
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Jets_Pt30_"+NjetBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet1_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet2_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);
            }
            if(jets_Pt30.size() >= 4){
              NjetBin = "Njet4";
              
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Jets_Pt30_"+NjetBin+"_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_"+NjetBin+"_"+IDsuffix, MET, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/MET_dxy_"+NjetBin+"_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_"+NjetBin+"_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Number_Vertex_noPU_"+NjetBin+"_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/Rho_"+NjetBin+"_"+IDsuffix, Rho, weight, 400, 0., 100.);
              FillHist(regions.at(it_rg)+"/NjetBinPt30/HT_"+NjetBin+"_"+IDsuffix, HT, weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet1_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(0).Pt(), weight, 1000, 0., 1000.);
              //FillHist(regions.at(it_rg)+"/NjetBinPt30/Jet2_Pt30_Pt_"+NjetBin+"_"+IDsuffix, jets_Pt30.at(1).Pt(), weight, 1000, 0., 1000.);
            }
          }

        }

      }
      else{

        // Cutflow : At least 2 jets, 1 b jets
        if(jets_Pt20.size()>1 && Nbjet_medium_PUveto>0){

          FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Medium_PUveto_"+IDsuffix, Nbjet_medium_PUveto, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          //FillHist(regions.at(it_rg)+"/MET_dxy_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          FillHist(regions.at(it_rg)+"/Number_Vertex_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Rho_"+IDsuffix, Rho, weight, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/HT_"+IDsuffix, HT, weight, 1000, 0., 1000.);

          if(MET > 40.){
            FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/Number_Jets_METgt40_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_BJets_Loose_METgt40_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_BJets_Medium_METgt40_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_BJets_Medium_PUveto_METgt40_"+IDsuffix, Nbjet_medium_PUveto, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_FatJets_METgt40_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/ZCand_Mass_METgt40_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
            FillHist(regions.at(it_rg)+"/ZCand_Pt_METgt40_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep1_Pt_METgt40_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep2_Pt_METgt40_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep1_Eta_METgt40_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/Lep2_Eta_METgt40_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/MET_METgt40_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            //FillHist(regions.at(it_rg)+"/MET_dxy_METgt40_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/MET2ST_METgt40_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

            FillHist(regions.at(it_rg)+"/Number_Vertex_METgt40_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
            FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_METgt40_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
            FillHist(regions.at(it_rg)+"/Rho_METgt40_"+IDsuffix, Rho, weight, 400, 0., 100.);
            FillHist(regions.at(it_rg)+"/HT_METgt40_"+IDsuffix, HT, weight, 1000, 0., 1000.);
          }
        }

        /*if(jets_Pt30.size()>1 && Nbjet_Pt30_medium>0){

          FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/Number_Jets_Pt30_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Pt30_Loose_"+IDsuffix, Nbjet_Pt30_loose, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_BJets_Pt30_Medium_"+IDsuffix, Nbjet_Pt30_medium, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/Number_FatJets_Pt30_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/ZCand_Mass_Pt30_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(regions.at(it_rg)+"/ZCand_Pt_Pt30_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Pt_Pt30_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep2_Pt_Pt30_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/Lep1_Eta_Pt30_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/Lep2_Eta_Pt30_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/MET_Pt30_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET_dxy_Pt30_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/MET2ST_Pt30_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          FillHist(regions.at(it_rg)+"/Number_Vertex_Pt30_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_Pt30_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
          FillHist(regions.at(it_rg)+"/Rho_Pt30_"+IDsuffix, Rho, weight, 400, 0., 100.);
          FillHist(regions.at(it_rg)+"/HT_Pt30_"+IDsuffix, HT, weight, 1000, 0., 1000.);

          if(MET > 40.){
            FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/Number_Jets_Pt30_METgt40_"+IDsuffix, jets_Pt30.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_BJets_Pt30_Loose_METgt40_"+IDsuffix, Nbjet_Pt30_loose, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_BJets_Pt30_Medium_METgt40_"+IDsuffix, Nbjet_Pt30_medium, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/Number_FatJets_Pt30_METgt40_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/ZCand_Mass_Pt30_METgt40_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
            FillHist(regions.at(it_rg)+"/ZCand_Pt_Pt30_METgt40_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep1_Pt_Pt30_METgt40_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep2_Pt_Pt30_METgt40_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/Lep1_Eta_Pt30_METgt40_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/Lep2_Eta_Pt30_METgt40_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/MET_Pt30_METgt40_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/MET_dxy_Pt30_METgt40_"+IDsuffix, MET_dxy, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/MET2ST_Pt30_METgt40_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

            FillHist(regions.at(it_rg)+"/Number_Vertex_Pt30_METgt40_"+IDsuffix, Nvtx, weight, 200, 0., 200.);
            FillHist(regions.at(it_rg)+"/Number_Vertex_noPU_Pt30_METgt40_"+IDsuffix, Nvtx, weight_noPU, 200, 0., 200.);
            FillHist(regions.at(it_rg)+"/Rho_Pt30_METgt40_"+IDsuffix, Rho, weight, 400, 0., 100.);
            FillHist(regions.at(it_rg)+"/HT_Pt30_METgt40_"+IDsuffix, HT, weight, 1000, 0., 1000.);
          }
        }*/

      }

    }

    //========================================================
    //==== WZ, ZG, WG control region
    //========================================================

    if(it_rg > 5){
      if(param.Electron_Tight_ID.Contains("HEEP")) continue;
    }

    if(it_rg>5 && it_rg<9){
      if(leptons.size() != 3) continue;

      // Passing triggers
      if(muons.size() >= 2){ if(!ev.PassTrigger(MuonTriggers)) continue; }
      if(electrons.size() >= 2){ if(!ev.PassTrigger(ElectronTriggers)) continue; }

      trigger_lumi = 1., dimu_trig_weight = 0.;
      if(!IsDATA){
        if(DataYear==2016){
          if(muons.size() >= 2){
            if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
            if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
            trigger_lumi = dimu_trig_weight;
          }
          if(electrons.size() >= 2) trigger_lumi = ev.GetTriggerLumi("Full");
        }
        else{
          trigger_lumi = ev.GetTriggerLumi("Full");
        }
      }

      // pT > trigger thresholds
      passPtCut = false;

      if(muons.size() >= 2){
        if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      }
      if(electrons.size() >= 2){
        if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      }

      if(!passPtCut) continue; 

      // Truth matching
      muons_prompt.clear();
      electrons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

      if(muons.size() != muons_prompt.size()) continue;
      if(electrons.size() != electrons_prompt.size()) continue;
 
      weight = 1.; 
      // weights for MC 
      if(!IsDATA){

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
          if(param.Muon_Tight_ID.Contains("HighPt")){
            muon_miniaodP = sqrt( muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz() );
            muon_recosf   = mcCorr->MuonReco_SF("HighPtMuonRecoSF", muons.at(i).Eta(), muon_miniaodP, 0);
            muon_idsf     = mcCorr->MuonID_SF("NUM_HighPtID_DEN_genTracks",  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
            muon_isosf    = mcCorr->MuonISO_SF("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut", muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          }
          else{
            muon_recosf = 1.;
            muon_idsf   = 1.;
            muon_isosf  = 1.;
          }
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int i=0; i<electrons.size(); i++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }

      }
      if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

      // Cutflow : 3 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      int l1 = -999, l2 = -999, l3 = -999, l4 = -999, wlepWZ = -999, wlepWG = -999;
      // OSSF lepton pair, W-tagged lepton
      if(muons.size()==2 && electrons.size()==1){
        if(muons.at(0).Charge()*muons.at(1).Charge() > 0) continue;
        ZCand = muons.at(0) + muons.at(1);
        WtagLep = electrons.at(0);
        ZtagLep1 = muons.at(0);
        ZtagLep2 = muons.at(1);
        GammaCand = ZCand;
        GammaLep1 = ZtagLep1;
        GammaLep2 = ZtagLep2;
      }

      if(muons.size()==1 && electrons.size()==2){
        if(electrons.at(0).Charge()*electrons.at(1).Charge() > 0) continue;
        ZCand = electrons.at(0) + electrons.at(1);
        WtagLep = muons.at(0);
        ZtagLep1 = electrons.at(0);
        ZtagLep2 = electrons.at(1);
        GammaCand = ZCand;
        GammaLep1 = ZtagLep1;
        GammaLep2 = ZtagLep2;
      }

      if(muons.size()==3 || electrons.size()==3){
        if(fabs(leptons.at(0)->Charge() + leptons.at(1)->Charge() + leptons.at(2)->Charge()) == 1){

          // ZCand, GammaCand
          double tmpMassDiff = 1000000., tmpMass = 100000.; 
          for(int ilep1=0; ilep1<2; ilep1++){
            for(int ilep2=ilep1+1; ilep2<3; ilep2++){
              if(leptons.at(ilep1)->Charge()*leptons.at(ilep2)->Charge() > 0) continue;
              Ztemp = *leptons.at(ilep1) + *leptons.at(ilep2);
              // For WZ, ZG
              if(!(Ztemp.M() > 10.)) ossf_mass10++;
              if(fabs(Ztemp.M() - MZ) < tmpMassDiff){
                tmpMassDiff = fabs(Ztemp.M() - MZ);
                ZCand = Ztemp; l1 = ilep1; l2 = ilep2;
              }
              // For WG
              if(Ztemp.M() < tmpMass){
                tmpMass = Ztemp.M();
                GammaCand = Ztemp; l3 = ilep1; l4 = ilep2;
              }
            }
          }

          ZtagLep1 = *leptons.at(l1);
          ZtagLep2 = *leptons.at(l2);
          GammaLep1 = *leptons.at(l3);
          GammaLep2 = *leptons.at(l4);

          // Set the lepton from W
          for(int ilep3=0; ilep3<3; ilep3++){
            if(fabs(ilep3-l1)>0 && fabs(ilep3-l2)>0) wlepWZ = ilep3;
            if(fabs(ilep3-l3)>0 && fabs(ilep3-l4)>0) wlepWG = ilep3;
          }

          if(it_rg < 8) WtagLep = *leptons.at(wlepWZ);
          else WtagLep = *leptons.at(wlepWG);

        }
        else continue;
      } 

      TriLep = *leptons.at(0) + *leptons.at(1) + *leptons.at(2);
      Mt = MT(WtagLep, METv);
      Mt3l = MT(TriLep, METv);

      // WZ, ZG control region
      if(it_rg < 8){
        if(!(ossf_mass10 == 0)) continue;
      
        // Cutflow : m(ll) > 10 GeV
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(Nbjet_medium == 0)) continue;

        // Cutflow : No b jets
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

        if(it_rg == 6){
          if(!IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET > 50.)) continue;
          if(!(Mt > 20.)) continue;
          if(!(TriLep.M() > MZ + 15.)) continue;
        }
        if(it_rg == 7){
          if(IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET < 50.)) continue;
          if(!IsOnZ(TriLep.M(), 15.)) continue;
        }
      }

      // WG control region 
      if(it_rg == 8){
        if(!(GammaCand.M() < 4.)) continue;

        // Cutflow : m(ll) < 4 GeV
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(Nbjet_medium == 0)) continue;

        // Cutflow : No b jets
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(MET > 30.)) continue;
        if(!(Mt3l > 30.)) continue;
      }

      // Histograms
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.); 
      FillHist(regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
      FillHist(regions.at(it_rg)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
      FillHist(regions.at(it_rg)+"/GammaCand_Mass_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
      FillHist(regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/GammaCand_Pt_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/GammaLep1_Pt_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/GammaLep2_Pt_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

      if(Nbjet_loose == 0){
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
        FillHist(regions.at(it_rg)+"/TriLep_Mass_NoLooseBJet_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
        FillHist(regions.at(it_rg)+"/GammaCand_Mass_NoLooseBJet_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
        FillHist(regions.at(it_rg)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/GammaCand_Pt_NoLooseBJet_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/WtagLep_Pt_NoLooseBJet_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/ZtagLep1_Pt_NoLooseBJet_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/ZtagLep2_Pt_NoLooseBJet_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/GammaLep1_Pt_NoLooseBJet_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/GammaLep2_Pt_NoLooseBJet_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Mt_NoLooseBJet_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      }

      for(unsigned int it_ch=0; it_ch<channels3L.size(); it_ch++){
        if(it_ch == electrons.size()){
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Jets_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaCand_Mass_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaCand_Pt_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaLep1_Pt_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaLep2_Pt_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          if(Nbjet_loose == 0){
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/TriLep_Mass_NoLooseBJet_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaCand_Mass_NoLooseBJet_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaCand_Pt_NoLooseBJet_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/WtagLep_Pt_NoLooseBJet_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZtagLep1_Pt_NoLooseBJet_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/ZtagLep2_Pt_NoLooseBJet_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaLep1_Pt_NoLooseBJet_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/GammaLep2_Pt_NoLooseBJet_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Mt_NoLooseBJet_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          }
        }
      } 

    }

    //========================================================
    //==== ZZ control region
    //========================================================

    if(it_rg == 9){
      if(leptons.size() != 4) continue;
      if((muons.size()==1 && electrons.size()==3) || (muons.size()==3 && electrons.size()==1)) continue;

      // Passing triggers
      if(muons.size() >= 2){ if(!ev.PassTrigger(MuonTriggers)) continue; }
      if(electrons.size() == 4){ if(!ev.PassTrigger(ElectronTriggers)) continue; }

      trigger_lumi = 1., dimu_trig_weight = 0.;
      if(!IsDATA){
        if(DataYear==2016){
          if(muons.size() >= 2){
            if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
            if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
            trigger_lumi = dimu_trig_weight;
          }
          if(electrons.size() == 4) trigger_lumi = ev.GetTriggerLumi("Full");
        }
        else{
          trigger_lumi = ev.GetTriggerLumi("Full");
        }
      }

      // pT > trigger thresholds
      passPtCut = false;

      if(muons.size() >= 2){
        if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      }
      if(electrons.size() == 4){
        if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      }

      if(!passPtCut) continue;

      // Truth matching
      muons_prompt.clear();
      electrons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

      if(muons.size() != muons_prompt.size()) continue;
      if(electrons.size() != electrons_prompt.size()) continue;

      weight = 1.;
      // weights for MC
      if(!IsDATA){

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
          if(param.Muon_Tight_ID.Contains("HighPt")){
            muon_miniaodP = sqrt( muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz() );
            muon_recosf   = mcCorr->MuonReco_SF("HighPtMuonRecoSF", muons.at(i).Eta(), muon_miniaodP, 0);
            muon_idsf     = mcCorr->MuonID_SF("NUM_HighPtID_DEN_genTracks",  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
            muon_isosf    = mcCorr->MuonISO_SF("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut", muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          }
          else{
            muon_recosf = 1.;
            muon_idsf   = 1.;
            muon_isosf  = 1.;
          }
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int i=0; i<electrons.size(); i++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }

      }
      if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

      // Cutflow : 4 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      // OSSF lepton pairs
      if(muons.size()==2 && electrons.size()==2){
        if(muons.at(0).Charge()*muons.at(1).Charge() > 0) continue;
        if(electrons.at(0).Charge()*electrons.at(1).Charge() > 0) continue;
        ZCand1 = muons.at(0) + muons.at(1);
        ZCand2 = electrons.at(0) + electrons.at(1);
      }

      if(muons.size()==4 || electrons.size()==4){
        if(leptons_minus.size() == leptons_plus.size()){
          Ztemp1 = *leptons_minus.at(0) + *leptons_plus.at(0);
          Ztemp2 = *leptons_minus.at(1) + *leptons_plus.at(1);
          Ztemp3 = *leptons_minus.at(0) + *leptons_plus.at(1);
          Ztemp4 = *leptons_minus.at(1) + *leptons_plus.at(0);
          if(!(Ztemp1.M()>10. && Ztemp2.M()>10. && Ztemp3.M()>10. && Ztemp4.M()>10.)) ossf_mass10++;
          ZCand1 = Ztemp1; ZCand2 = Ztemp2;

          if(!(IsOnZ(ZCand1.M(), 15.) && IsOnZ(ZCand2.M(), 15.))){
            ZCand1 = Ztemp3; ZCand2 = Ztemp4;
          }
        }
      }

      if(!(ossf_mass10 == 0)) continue;

      // Cutflow : m(ll) > 10 GeV
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(Nbjet_medium == 0)) continue;

      // Cutflow : No b jets
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(IsOnZ(ZCand1.M(), 15.) && IsOnZ(ZCand2.M(), 15.))) continue;

      // Histograms 
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/ZCand1_Mass_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
      FillHist(regions.at(it_rg)+"/ZCand2_Mass_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
      FillHist(regions.at(it_rg)+"/ZCand1_Pt_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/ZCand2_Pt_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep4_Pt_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep4_Eta_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      if(Nbjet_loose == 0){
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/ZCand1_Mass_NoLooseBJet_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
        FillHist(regions.at(it_rg)+"/ZCand2_Mass_NoLooseBJet_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
        FillHist(regions.at(it_rg)+"/ZCand1_Pt_NoLooseBJet_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/ZCand2_Pt_NoLooseBJet_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep4_Pt_NoLooseBJet_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep4_Eta_NoLooseBJet_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      }

      for(unsigned int it_ch=0; it_ch<channels4L.size(); it_ch++){
        if(it_ch == electrons.size()/2){
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Jets_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand1_Mass_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand2_Mass_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand1_Pt_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand2_Pt_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep4_Pt_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep4_Eta_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          if(Nbjet_loose == 0){
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets_Pt20.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand1_Mass_NoLooseBJet_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand2_Mass_NoLooseBJet_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand1_Pt_NoLooseBJet_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/ZCand2_Pt_NoLooseBJet_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep4_Pt_NoLooseBJet_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Lep4_Eta_NoLooseBJet_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          }
        } 
      }

    } // ZZ Control region

  } // Control Region Loop

}




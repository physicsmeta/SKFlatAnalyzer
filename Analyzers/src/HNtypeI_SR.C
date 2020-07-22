#include "HNtypeI_SR.h"

HNtypeI_SR::HNtypeI_SR(){

}

void HNtypeI_SR::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  cout << "[HNtypeI_SR::initializeAnalyzer] RunSyst = " << RunSyst << endl;
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
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SR.h 
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
    //MuonTriggers.push_back("HLT_Mu17_Mu8_SameSign_DZ_v");
    //MuonTriggersH.push_back("HLT_Mu17_Mu8_SameSign_DZ_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 27267.591112919
    //EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
    EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028
    //EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
    //EMuPtCut1 = 25., EMuPtCut2 = 10.;
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

  //cout << "[HNtypeI_SR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_SR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_SR::~HNtypeI_SR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_SR::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_SR.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  //AllFatJets = GetAllFatJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets());

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SR.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<MuonTightIDs.size(); it_id++){
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

    if(RunFake || RunCF){
      if(it_id > 0) continue;
    }

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

void HNtypeI_SR::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> channels = {"dimu", "diel", "emu"};
  vector<TString> regions = {"fakeCR1", "lowSR1", "lowCR1", "highSR1", "highCR1", "lowSR2", "lowCR2", "highSR2", "highCR2"};
  TString IDsuffix = "HN16";
  if(param.Electron_Tight_ID.Contains("HEEP")) IDsuffix = "LRSM";
  if(param.Electron_Tight_ID.Contains("V1")) IDsuffix = "HNV1";
  if(param.Electron_Tight_ID.Contains("V2")) IDsuffix = "HNV2";
  if(param.Electron_Tight_ID.Contains("MVA")) IDsuffix = "MVA";
  //TString LepCategory = "TT";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;
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
  }

  // Cutflow : No Cuts
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, -0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, -0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, -0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, -0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  // Cutflow : MET filter
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== Trigger
  //========================================================

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
    cout << "[HNtypeI_SR::executeEventFromParameter] Wrong syst" << endl;
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

  // For charge flip
  vector<Electron> electrons_beforeShift;
  vector<Electron> electrons_afterShift;
  electrons_beforeShift.clear();
  electrons_afterShift.clear();

  // Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  vector<Jet> jets_PUveto = JetsPassPileupMVA(jets_lepveto);
  vector<Jet> jets = JetsAwayFromFatJet(jets_PUveto, fatjets);  // AK4jets used in SR, CR

  vector<Jet> jets_WCandLowMass;
  vector<Jet> jets_WCandHighMass;
  FatJet fatjets_WCand;
  jets_WCandLowMass.clear();
  jets_WCandHighMass.clear();

  std::vector<Lepton*> leptons, leptons_veto;

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

  //========================================================
  //==== B-Tagging 
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0;
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
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  if(muons.size()+electrons.size() == 2){
    METv = UpdateMETMuon(METv, muons);
    METv = UpdateMETElectron(METv, electrons);
  }

  double MET = METv.Pt();

  //========================================================
  //==== Define particles, variables
  //========================================================

  double ST = 0., MET2ST = 0.;
  double dRll = 0., dRl2WCand = 0.;
  //double MZ = 91.1876;
  double MW = 80.379;
  double muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_recosf = 1., ele_idsf = 1.;
  int lepton_veto_size = 0;

  bool passPtCut = false;
  Particle ZCand, Ztemp;
  Particle WCand, lljjLow, l1jjLow, l2jjLow, lljjHigh, l1jjHigh, l2jjHigh;  // Low Mass SR1, High Mass SR1
  Particle llj, l1j, l2j;                                                   // Low Mass SR2
  Particle l1J, l2J;                                                        // High Mass SR2
  
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
    if(muons.size()+electrons.size() == 2){
      METv = UpdateMETFake(METv, electrons, muons);
    }

    muons = MuonUsePtCone(muons);
    electrons = ElectronUsePtCone(electrons);
    std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
  }

  // Shift electron energy and MET when RunCF=true
  if(RunCF){
    if(muons.size()==0 && electrons.size()==2){
      electrons_beforeShift.push_back(electrons.at(0));
      electrons_beforeShift.push_back(electrons.at(1));
      electrons = ShiftElectronEnergy(electrons, param, true);
      electrons_afterShift.push_back(electrons.at(0));
      electrons_afterShift.push_back(electrons.at(1));
      METv = UpdateMETElectronCF(METv, electrons_beforeShift, electrons_afterShift);
    }
  }

  // Define leptons
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));

  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  // Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  lepton_veto_size = leptons_veto.size() - leptons.size();

  // Define ST, MET^2/ST
  MET = METv.Pt();

  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;

  //========================================================
  //==== Event selection
  //========================================================

  // Loop for each channel : it_ch (0,1,2) = (mumu, ee, emu)
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    trigger_lumi = 1.;

    if(it_ch==0 || it_ch==2){ if(RunCF) continue; }

    // Cutflow : passing dilepton triggers
    if(it_ch==0){ if(!ev.PassTrigger(MuonTriggers)) continue; }
    if(it_ch==1){ if(!ev.PassTrigger(ElectronTriggers)) continue; }
    if(it_ch==2){ if(!ev.PassTrigger(EMuTriggers)) continue; }

    trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
    if(!IsDATA){
      if(DataYear==2016){
        if(it_ch==0){
          if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
          if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
          trigger_lumi = dimu_trig_weight;
        }
        if(it_ch==1) trigger_lumi = ev.GetTriggerLumi("Full");
        if(it_ch==2){
          if(ev.PassTrigger(EMuTriggers)) emu_trig_weight += 27267.591;
          if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628;
          trigger_lumi = emu_trig_weight;
        }
      }
      else{
        trigger_lumi = ev.GetTriggerLumi("Full");
      }
    }

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      weight = 1.;
      if(!IsDATA){
        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max); 

    //========================================================
    //==== Preselection
    //========================================================

    if(leptons.size() == 2){ 

      // Cutflow : 2 tight leptons (truth-matched, pT > trigger thresholds)

      // pT > trigger thresholds
      passPtCut = false;

      if(it_ch==0){
        if(!(muons.size()==2 && electrons.size()==0)) continue;
        if(leptons.at(0)->Pt()>MuonPtCut1 && leptons.at(1)->Pt()>MuonPtCut2) passPtCut = true;
      }
      if(it_ch==1){
        if(!(muons.size()==0 && electrons.size()==2)) continue;
        if(leptons.at(0)->Pt()>ElectronPtCut1 && leptons.at(1)->Pt()>ElectronPtCut2) passPtCut = true;
      }
      if(it_ch==2){
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

      if(it_ch==0){
        if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) continue;
      }
      if(it_ch==1){
        if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) continue;
      }
      if(it_ch==2){
        if(!(muons_prompt.size()==1 && electrons_prompt.size()==1)) continue;
      }

      // weights for MC
      weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_recosf = 1., ele_idsf = 1.;

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

      // weights for fake, CF
      if(RunFake) weight *= fakeEst->GetWeight(leptons, param);
      if(RunCF) weight *= GetCFweight(leptons, param, true, 0);

      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      ZCand = *leptons.at(0) + *leptons.at(1);
      dRll  = leptons.at(0)->DeltaR(*leptons.at(1));

      // Cutflow : same-sign (oppsite-sign when RunCF=true) 

      if(!RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) continue;
      if(RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) continue;

      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      // Cutflow : veto 3rd leptons using veto ID

      if(lepton_veto_size > 0) continue;

      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      // Cutflow : m(ll) > 10 GeV (|m(ll)-m(Z)| > 10 GeV in dielectron channel)

      if(!(ZCand.M() > 10.)) continue;
      if(it_ch==1 && IsOnZ(ZCand.M(), 10.)) continue;

      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);     

      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"PreNoJetCut/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      // non-prompt CR2 : no jets && same-sign back-to-back 2 leptons
      if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){
       
        // Cutflow : jet requirement for non-prompt CR2 
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5)) continue;
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/fakeCR2/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      }

      // Cutflow : jet requirement
      //if(!(jets.size()>0 || fatjets.size()>0)) continue; 
      if(!(fatjets.size()>0) && !(jets.size()>1 && fatjets.size()==0) && !(jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.)) continue;
      
      FillHist(channels.at(it_ch)+"/Pre/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      // Event selections for each CR
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        jets_WCandLowMass.clear();
        jets_WCandHighMass.clear();

        // This is the number or events at preselection
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

        // non-prompt CR1 : SS 2 leptons with b-tagged jets
        if(it_rg == 0){

          if(!(Nbjet_medium > 0)) continue;

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

        }
      
        // Low mass SR1, CR1 & High mass SR1, CR1
        if(it_rg>=1 && it_rg<5){

          if(!(jets.size()>=2 && fatjets.size()==0)) continue;

          jets_WCandLowMass  = JetsWCandLowMass(*leptons.at(0), *leptons.at(1), jets, MW);
          jets_WCandHighMass = JetsWCandHighMass(jets, MW);

          WCand     = jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
          lljjLow   = *leptons.at(0) + *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
          l1jjLow   = *leptons.at(0) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
          l2jjLow   = *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
          lljjHigh  = *leptons.at(0) + *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
          l1jjHigh  = *leptons.at(0) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
          l2jjHigh  = *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
          dRl2WCand = leptons.at(0)->DeltaR(WCand);

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max); 
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand_Mass_nocut_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCandJet1_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCandJet2_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljjLow_Mass_nocut_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jjLow_Mass_nocut_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jjLow_Mass_nocut_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljjHigh_Mass_nocut_"+IDsuffix, lljjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jjHigh_Mass_nocut_"+IDsuffix, l1jjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jjHigh_Mass_nocut_"+IDsuffix, l2jjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2WCand_DeltaR_nocut_"+IDsuffix, dRl2WCand, weight, 60, 0., 6.);

          // Low mass SR1 
          if(it_rg == 1){
            if(!(Nbjet_medium == 0)) continue;
            if(!(lljjLow.M() < 300.)) continue;
            if(!(MET < 80.)) continue;
          }

          // Low mass CR1
          if(it_rg == 2){
            if(!(lljjLow.M() < 300.)) continue;
            if(!(Nbjet_medium>0 || MET>100.)) continue;
          }

          // High mass SR1
          if(it_rg == 3){
            if(!(Nbjet_medium == 0)) continue;
            if(!(WCand.M() < 150.)) continue;
            if(!(MET2ST < 15.)) continue;
          }

          // High mass CR1
          if(it_rg == 4){
            if(!(WCand.M() < 150.)) continue;
            if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand_Mass_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCandJet1_Pt_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCandJet2_Pt_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljjLow_Mass_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jjLow_Mass_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jjLow_Mass_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljjHigh_Mass_"+IDsuffix, lljjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jjHigh_Mass_"+IDsuffix, l1jjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jjHigh_Mass_"+IDsuffix, l2jjHigh.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2WCand_DeltaR_"+IDsuffix, dRl2WCand, weight, 60, 0., 6.);

          /*if(it_rg==1) cout << "In LowSR1, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;
          if(it_rg==3) cout << "In HighSR1, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;
          if(it_ch==0 && it_rg==3){
            if(run==276283 && event==1252562683){
              cout << "Muon1 : (" << muons.at(0).Pt() << ", " << muons.at(0).Eta() << ", " << muons.at(0).Phi() << ", " << muons.at(0).E() << ")" << endl;
              cout << "Muon2 : (" << muons.at(1).Pt() << ", " << muons.at(1).Eta() << ", " << muons.at(1).Phi() << ", " << muons.at(1).E() << ")" << endl;
              for(unsigned int i=0; i<jets.size(); i++){
                cout << TString::Itoa(i, 10)+"th Jets : (" << jets.at(i).Pt() << ", " << jets.at(i).Eta() << ", " << jets.at(i).Phi() << ", " << jets.at(i).E() << ")" << endl;
              }
            }
          }*/

        }

        // Low mass SR2, CR2
        if(it_rg>=5 && it_rg<7){

          if(!(jets.size()==1 && fatjets.size()==0)) continue;

          llj = *leptons.at(0) + *leptons.at(1) + jets.at(0);
          l1j = *leptons.at(0) + jets.at(0);
          l2j = *leptons.at(1) + jets.at(0);

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/llj_Mass_nocut_"+IDsuffix, llj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1j_Mass_nocut_"+IDsuffix, l1j.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2j_Mass_nocut_"+IDsuffix, l2j.M(), weight, 2000, 0., 2000.);

          // Low mass SR2
          if(it_rg == 5){
            if(!(Nbjet_medium == 0)) continue;
            if(!(llj.M() < 300.)) continue;
            if(!(MET < 80.)) continue;
          }

          // Low mass CR2
          if(it_rg == 6){
            if(!(llj.M() < 300.)) continue;
            if(!(Nbjet_medium>0 || MET>100.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/llj_Mass_"+IDsuffix, llj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1j_Mass_"+IDsuffix, l1j.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2j_Mass_"+IDsuffix, l2j.M(), weight, 2000, 0., 2000.);

          //if(it_rg==5) cout << "In LowSR2, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;

        }

        // High mass SR2, CR2
        if(it_rg >= 7){

          if(!(fatjets.size() > 0)) continue;

          fatjets_WCand = FatJetWCand(fatjets, MW);

          l1J = *leptons.at(0) + fatjets_WCand;
          l2J = *leptons.at(1) + fatjets_WCand;

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Pt_nocut_"+IDsuffix, fatjets_WCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Mass_nocut_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1J_Mass_nocut_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2J_Mass_nocut_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

          // High mass SR2
          if(it_rg == 7){
            if(!(Nbjet_medium == 0)) continue;
            if(!(fatjets_WCand.SDMass() < 150.)) continue;
            if(!(MET2ST < 15.)) continue;
          }

          // High mass CR2
          if(it_rg == 8){
            if(!(fatjets_WCand.SDMass() < 150.)) continue;
            if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Pt_"+IDsuffix, fatjets_WCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Mass_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2J_Mass_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

          //if(it_rg==7) cout << "In HighSR2, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;

        }

      } // Region Loop 
    } // Dilepton selection
  } // Channel Loop

}


#include "HNtypeI_SM_CR_2016H.h"

HNtypeI_SM_CR_2016H::HNtypeI_SM_CR_2016H(){

}

void HNtypeI_SM_CR_2016H::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[HNtypeI_SM_CR_2016H::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  RunFake = HasFlag("RunFake");
  RunCF = HasFlag("RunCF");

/*  MuonTightIDs = {"HNTight", "HNTightV2", "HNTight2016"};
  MuonLooseIDs = {"HNLoose", "HNLoose", "HNLoose2016"};
  MuonVetoIDs  = {"POGLoose", "POGLoose", "HNVeto2016"};
//  MuonIDSFKeys = { "NUM_TightID_DEN_genTracks" };
  ElectronTightIDs = {"HNTight", "HNTightV2", "HNTight2016"};
  ElectronLooseIDs = {"HNLoose", "HNLooseV23", "HNLoose2016"};
  ElectronVetoIDs  = {"passVetoID", "passVetoID", "HNVeto2016"};
  FakeRateIDs = {"HNtypeI_V1", "HNtypeI_V2", "HNtypeI_16"};*/
  MuonTightIDs = {"HNTightV2", "HNTight2016"};
  MuonLooseIDs = {"HNLoose", "HNLoose2016"};
  MuonVetoIDs  = {"POGLoose", "HNVeto2016"};
  ElectronTightIDs = {"HNTightV2", "HNTight2016"};
  ElectronLooseIDs = {"HNLooseV23", "HNLoose2016"};
  ElectronVetoIDs  = {"passVetoID", "HNVeto2016"};
  FakeRateIDs = {"HNtypeI_V2", "HNtypeI_16"};


  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SM_CR_2016H.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  MuonTriggers.clear();
  MuonTriggersH.clear();
  ElectronTriggers.clear();
//  EMuTriggers.clear();
//  EMuTriggersH.clear();

  if(DataYear==2016){                                                                   // Lumi values for trigger weight (/pb)
//    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                       // 27267.591112919   // NOTE : Change for 2016H
//    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                     // 27267.591112919   // NOTE : Change for 2016H
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                    // 35918.219492947
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                  // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                   // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                 // 35918.219492947
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
//    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
//    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 27267.591112919
//    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
//   EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028

    // These are needed for applying lepton pT cuts 
/*    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");*/

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
//    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
//    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
//    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    // These are needed for applying lepton pT cuts 
//    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
//    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
//    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2018){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
//    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
//    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
//    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

//  cout << "[HNtypeI_SM_CR_2016H::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[HNtypeI_SM_CR_2016H::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_SM_CR_2016H::~HNtypeI_SM_CR_2016H(){

  //==== Destructor of this Analyzer

}

void HNtypeI_SM_CR_2016H::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_SM_CR_2016H.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SM_CR_2016H.h
//  weight_Prefire = GetPrefireWeight(0);

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

//  param.Name = MuonID+"_"+"Central";

    // Muon ID
    param.Muon_Tight_ID       = MuonTightID;
    param.Muon_Loose_ID       = MuonLooseID;
    param.Muon_Veto_ID        = MuonVetoID;
    param.Muon_FR_ID          = FakeRateID;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key         = "AwayJetPt40";  // histname
    param.Muon_ID_SF_Key      = "NUM_TightID_DEN_genTracks";
    param.Muon_ISO_SF_Key     = "NUM_TightRelIso_DEN_TightIDandIPCut";
    param.Muon_Trigger_SF_Key = "";
    param.Muon_UsePtCone      = true;

    // Electron ID
    param.Electron_Tight_ID       = ElectronTightID;
    param.Electron_Loose_ID       = ElectronLooseID;
    param.Electron_Veto_ID        = ElectronVetoID;
    param.Electron_FR_ID          = FakeRateID;     // ID name in histmap_Electron.txt
    param.Electron_FR_Key         = "AwayJetPt40";  // histname
    param.Electron_ID_SF_Key      = "passTightID";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone      = true;

    // Jet ID
//    param.Jet_ID = "tightLepVeto";
    param.Jet_ID    = "HNTight";
    param.FatJet_ID = "HNTight";

    executeEventFromParameter(param);

/*    if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }*/
  }
}

void HNtypeI_SM_CR_2016H::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> regions = {"DYmm", "DYee", "WZ", "ZG", "WG", "ZZ"}; 
  vector<TString> channels3L = {"mmm", "mme", "mee", "eee"};
  vector<TString> channels4L = {"mmmm", "mmee", "eeee"};
  TString IDsuffix = "HNV1";
  if(param.Electron_Tight_ID.Contains("V2")) IDsuffix = "HNV2";
  if(param.Electron_Tight_ID.Contains("2016")) IDsuffix = "HN16";
  TString LepCategory = "TT";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0.;
 
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

  if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers))) return; 

  //========================================================
  //==== Copy AllObjects
  //========================================================

  vector<Muon> this_AllMuons = AllMuons;
  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Jet> this_AllJets = AllJets;
  vector<FatJet> this_AllFatJets = AllFatJets;
  vector<Gen> gens = GetGens();

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

/*  if(param.syst_ == AnalyzerParameter::Central){

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
    cout << "[HNtypeI_SM_CR_2016H::executeEventFromParameter] Wrong syst" << endl;
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

//  int lepton_count1 = 0, lepton_count2 = 0, fatjet_count = 0; 

/*  for(unsigned int i=0; i<this_AllFatJets.size(); i++){
    lepton_count1 = 0;
    if(!(this_AllFatJets.at(i).PassID(param.FatJet_ID))) continue;
    if(!(this_AllFatJets.at(i).Pt() > 200.)) continue;
    if(!(fabs(this_AllFatJets.at(i).Eta()) < 2.7)) continue;
    for(unsigned int j=0; j<muons_veto.size(); j++){
      if(this_AllFatJets.at(i).DeltaR(muons_veto.at(j)) < 1.0) lepton_count1++;
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllFatJets.at(i).DeltaR(electrons_veto.at(j)) < 1.0) lepton_count1++;
    } 
    if(lepton_count1 > 0) continue;
    fatjets.push_back(this_AllFatJets.at(i));
  }

  for(unsigned int i=0; i<this_AllJets.size(); i++){
    lepton_count2 = 0, fatjet_count = 0;
    if(!(this_AllJets.at(i).PassID(param.Jet_ID))) continue;
    if(!(this_AllJets.at(i).Pt() > 20.)) continue;
    if(!(fabs(this_AllJets.at(i).Eta()) < 2.7)) continue;
    for(unsigned int j=0; j<muons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(muons_veto.at(j)) < 0.4) lepton_count2++;
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(electrons_veto.at(j)) < 0.4) lepton_count2++;
    }
    for(unsigned int j=0; j<fatjets.size(); j++){
      if(this_AllJets.at(i).DeltaR(fatjets.at(j)) < 0.8) fatjet_count++;
    }
    if(lepton_count2 > 0) continue;
    if(fatjet_count > 0) continue;
    jets.push_back(this_AllJets.at(i));
  }*/

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

  //========================================================
  //==== B-Tagging
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0;
  JetTagging::Parameters jtp_DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
//  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
//    double this_discr = jets_nolepveto.at(ij).GetTaggerResult(JetTagging::DeepCSV);
      //==== No SF
//      if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  METv = UpdateMETMuon(METv, this_AllMuons);
  METv = UpdateMETElectron(METv, this_AllElectrons);

  double MET = METv.Pt();

  //========================================================
  //==== Define particles, variables
  //========================================================

  double Mt = 0.;
  double Mt3l = 0.;
  double ST = 0.;
  double MET2ST = 0.;
  double MZ = 91.1876;
//  double MW = 80.379;
  double muon_idsf = 1.;
  double muon_isosf = 1.;
  double ele_idsf = 1.;
  double ele_recosf = 1.;
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
    if(IDsuffix == "HNV2") mu_tight_iso = 0.1;
    if(IDsuffix == "HN16") mu_tight_iso = 0.07;

    for(unsigned int i=0; i<muons.size(); i++){
      this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso);
      muons.at(i).SetPtCone(this_ptcone_muon);
    }

    for(unsigned int i=0; i<electrons.size(); i++){
      el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt();
      if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
      if(IDsuffix == "HNV2"){
        el_tight_iso = std::min(0.08, 0.0287+0.506/electrons.at(i).UncorrPt());
        if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/electrons.at(i).UncorrPt());
      }
      if(IDsuffix == "HN16") el_tight_iso = 0.08;
      this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
      electrons.at(i).SetPtCone(this_ptcone_electron);
    }

    // Correct MET when RunFake=true, because pT was replaced by pTcone
    if(leptons.size()>1 && leptons.size()<5){
      METv = UpdateMETFake(METv, electrons, muons);
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

  // Define ST, MET^2/ST
  MET = METv.Pt();

  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;

  //=====================================================================================
  //==== SM background CR (DYmm, DYee, WZ, ZG, WG, ZZ)
  //=====================================================================================

  if(RunCF) return;

  // Period-dependent trigger weight (only for 2016 MC)
  trigger_lumi = 1., dimu_trig_weight = 0.;
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
    }
    else{
      trigger_lumi = ev.GetTriggerLumi("Full");
    }
  }

  //========================================================
  //==== Event selection
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
    ossf_mass10 = 0;

    if(!IsDATA){
      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);
    }

    // Cutflow : passing dilepton triggers (dimu || diel)
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }

    //========================================================
    //==== DY control region
    //========================================================

    if(it_rg < 2){
      if(leptons.size() != 2) continue;
      if(muons.size()==1 && electrons.size()==1) continue;

      // Passing triggers
      if(it_rg == 0){ if(!ev.PassTrigger(MuonTriggers)) continue; }
      if(it_rg == 1){ if(!ev.PassTrigger(ElectronTriggers)) continue; }

      trigger_lumi = 1., dimu_trig_weight = 0.;
      if(!IsDATA){
        if(DataYear==2016){
          if(it_rg == 0){
            if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
            if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
            trigger_lumi = dimu_trig_weight;
          }
          if(it_rg == 1) trigger_lumi = ev.GetTriggerLumi("Full");
        }
        else{
          trigger_lumi = ev.GetTriggerLumi("Full");
        }
      } 

      // pT > trigger thresholds
      passPtCut = false;

      if(it_rg == 0){
        if(!(muons.size()==2 && electrons.size()==0)) continue;
        if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      }
      if(it_rg == 1){
        if(!(muons.size()==0 && electrons.size()==2)) continue;
        if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      }

      if(!passPtCut) continue;

      // Truth matching
      muons_prompt.clear();
      electrons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

      if(it_rg == 0){
        if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) continue;
      }
      if(it_rg == 1){
        if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) continue;
      }

      weight = 1.;
      // weights for MC
      if(!IsDATA){

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
          weight *= muon_idsf*muon_isosf;
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      ZCand = *leptons.at(0) + *leptons.at(1);

      // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) continue;

      // Cutflow : OS event
      FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
      FillHist(regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

      if(!(Nbjet_medium == 0)) continue;

      // Cutflow : No b jets
      FillHist(regions.at(it_rg)+"/Number_Events_NoMediumBJet_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_NoMediumBJet_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Jets_NoMediumBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Loose_NoMediumBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Medium_NoMediumBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_FatJets_NoMediumBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/ZCand_Mass_NoMediumBJet_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
      FillHist(regions.at(it_rg)+"/ZCand_Pt_NoMediumBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Pt_NoMediumBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep2_Pt_NoMediumBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Eta_NoMediumBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep2_Eta_NoMediumBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/MET_NoMediumBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/MET2ST_NoMediumBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      if(Nbjet_loose == 0){
        FillHist(regions.at(it_rg)+"/Number_Events_NoLooseBJet_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(regions.at(it_rg)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(regions.at(it_rg)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions.at(it_rg)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(regions.at(it_rg)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      }

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID
      FillHist(regions.at(it_rg)+"/Number_Events_No3rdLep_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Events_unweighted_No3rdLep_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regions.at(it_rg)+"/Number_Jets_No3rdLep_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Loose_No3rdLep_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_BJets_Medium_No3rdLep_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/Number_FatJets_No3rdLep_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regions.at(it_rg)+"/ZCand_Mass_No3rdLep_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
      FillHist(regions.at(it_rg)+"/ZCand_Pt_No3rdLep_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Pt_No3rdLep_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep2_Pt_No3rdLep_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/Lep1_Eta_No3rdLep_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/Lep2_Eta_No3rdLep_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions.at(it_rg)+"/MET_No3rdLep_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regions.at(it_rg)+"/MET2ST_No3rdLep_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        
    }

    //========================================================
    //==== WZ, ZG, WG control region
    //========================================================

    if(it_rg>1 && it_rg<5){
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
          weight *= muon_idsf*muon_isosf;
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

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
          if(it_rg < 4) WtagLep = *leptons.at(wlepWZ);
          else WtagLep = *leptons.at(wlepWG);

        }
        else continue;
      } 

      TriLep = *leptons.at(0) + *leptons.at(1) + *leptons.at(2);
      Mt = MT(WtagLep, METv);
      Mt3l = MT(TriLep, METv);

      // WZ, ZG control region
      if(it_rg < 4){
        if(!(ossf_mass10 == 0)) continue;
      
        // Cutflow : m(ll) > 10 GeV
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(Nbjet_medium == 0)) continue;

        // Cutflow : No b jets
        FillHist(regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

        if(it_rg == 2){
          if(!IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET > 50.)) continue;
          if(!(Mt > 20.)) continue;
          if(!(TriLep.M() > MZ + 15.)) continue;
        }
        if(it_rg == 3){
          if(IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET < 50.)) continue;
          if(!IsOnZ(TriLep.M(), 15.)) continue;
        }
      }

      // WG control region 
      if(it_rg == 4){
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
      FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.); 
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
        FillHist(regions.at(it_rg)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
          FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels3L.at(it_ch)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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

    if(it_rg == 5){
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
          weight *= muon_idsf*muon_isosf;
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

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
      FillHist(regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
        FillHist(regions.at(it_rg)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(regions.at(it_rg)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
          FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(regions.at(it_rg)+"/"+channels4L.at(it_ch)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
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
    }
  } 
}




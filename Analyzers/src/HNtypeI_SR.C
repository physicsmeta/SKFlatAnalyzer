#include "HNtypeI_SR.h"

HNtypeI_SR::HNtypeI_SR(){

}

void HNtypeI_SR::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[HNtypeI_SR::initializeAnalyzer] RunSyst = " << RunSyst << endl;
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
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SR.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  // We can use SkimTree_SMP for dimuon and dielectron channel. (See https://github.com/sansan9401/SKFlatAnalyzer/blob/Run2Legacy_hsseo/Analyzers/src/SkimTree_SMP.C)

  MuonTriggers.clear();
  ElectronTriggers.clear();

  if(DataYear==2016){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
/*    MuonTriggers.push_back("HLT_IsoMu24_v");
    MuonTriggers.push_back("HLT_IsoTkMu24_v");
    ElectronTriggers.push_back("HLT_Ele27_WPTight_Gsf_v");
    MuonPtCut = 26.;
    ElectronPtCut = 29.;*/
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
/*    MuonTriggers.push_back("HLT_IsoMu27_v");
//    ElectronTriggers.push_back("HLT_Ele27_WPTight_Gsf_L1DoubleEG_v");
    ElectronTriggers.push_back("HLT_Ele35_WPTight_Gsf_v");
    MuonPtCut = 29.;
    ElectronPtCut = 37.;*/
  }
  else if(DataYear==2018){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
  }

//  cout << "[HNtypeI_SR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[HNtypeI_SR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== Test btagging code
  //==== add taggers and WP that you want to use in analysis
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);
//  vtaggers.push_back(Jet::CSVv2);  // available for 2016 only

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Loose);
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs
  SetupBTagger(vtaggers, v_wps, true, true);

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
//  AllFatJets = GetAllFatJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets());

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SR.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_EleID=0; it_EleID<ElectronTightIDs.size(); it_EleID++){
    for(unsigned int it_MuonID=0; it_MuonID<MuonTightIDs.size(); it_MuonID++){
      if(it_EleID != it_MuonID) continue; 
      TString MuonTightID = MuonTightIDs.at(it_MuonID);
      TString MuonLooseID = MuonLooseIDs.at(it_MuonID); 
      TString MuonVetoID  = MuonVetoIDs.at(it_MuonID);
      TString ElectronTightID = ElectronTightIDs.at(it_EleID);
      TString ElectronLooseID = ElectronLooseIDs.at(it_EleID);
      TString ElectronVetoID  = ElectronVetoIDs.at(it_EleID);
      TString FakeRateID = FakeRateIDs.at(it_EleID);

      param.Clear();

      param.syst_ = AnalyzerParameter::Central;

//      param.Name = MuonID+"_"+"Central";

      // Muon ID
      param.Muon_Tight_ID = MuonTightID;
      param.Muon_Loose_ID = MuonLooseID;
      param.Muon_Veto_ID  = MuonVetoID;
      param.Muon_FR_ID = FakeRateID;     // ID name in histmap_Muon.txt
      param.Muon_FR_Key = "AwayJetPt40"; // histname
      param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
      param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
      param.Muon_Trigger_SF_Key = "";
      param.Muon_UsePtCone = true;

      // Electron Id
      param.Electron_Tight_ID = ElectronTightID;
      param.Electron_Loose_ID = ElectronLooseID;
      param.Electron_Veto_ID  = ElectronVetoID;
      param.Electron_FR_ID = FakeRateID;     // ID name in histmap_Electron.txt
      param.Electron_FR_Key = "AwayJetPt40"; // histname
      param.Electron_ID_SF_Key = "passTightID";
      param.Electron_Trigger_SF_Key = "";
      
      param.Electron_UsePtCone = true;

      // Jet ID
//      param.Jet_ID = "tightLepVeto";
      param.Jet_ID = "HNTight";
      param.FatJet_ID = "HNTight";

      executeEventFromParameter(param);

/*  if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }*/
    }
  }
}

void HNtypeI_SR::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> channels = {"dimu", "diel"};
  vector<TString> regions = {"fakeCR1", "lowSR1", "lowCR1", "highSR1", "highCR1", "lowSR2", "lowCR2", "highSR2", "highCR2"};
  vector<TString> regionsSM = {"WZ", "ZG", "ZZ"}; // "WG"
  vector<TString> channels3L = {"eee", "eem", "emm", "mmm"};
  vector<TString> channels4L = {"eeee", "eemm", "mmmm"};
  TString IDsuffix = "HNV1";
  if(param.Electron_Tight_ID.Contains("V2")) IDsuffix = "HNV2";
  if(param.Electron_Tight_ID.Contains("2016")) IDsuffix = "HN16";
  TString LepCategory = "TT";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;

  Event ev = GetEvent();

  //=============
  //==== No Cut
  //=============

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
  for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, -0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, -0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================
  //==== MET Filter
  //========================

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
  for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  Particle METv = ev.GetMETVector();

  //==============
  //==== Trigger
  //==============

  if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers))) return;

/*  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    if(it_ch==0 && !ev.PassTrigger(MuonTriggers)) continue;
    if(it_ch==1 && !ev.PassTrigger(ElectronTriggers)) continue;
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){  // For SM CR, apply the trigger selection later
      FillHist("Number_Events_"+channels.at(it_ch)+"_"+regions.at(it_rg)+"_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    }
  }*/

  //======================
  //==== Copy AllObjects
  //======================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;
  vector<FatJet> this_AllFatJets = AllFatJets;
  vector<Gen> gens = GetGens();

//  FillHist("Nfatjet_"+IDsuffix, this_AllFatJets.size(), weight, 5, 0., 5.);

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
    cout << "[HNtypeI_SR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

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
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, "tight", 20., 2.4);
//  vector<FatJet> fatjets = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

//  FillHist("Njet_"+IDsuffix, jets_nolepveto.size(), weight, 8, 0., 8.);

  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  vector<Jet> jets;
  vector<FatJet> fatjets;
  jets.clear();
  fatjets.clear();
  int lepton_count1 = 0, lepton_count2 = 0, fatjet_count = 0; 

  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124
  for(unsigned int i=0; i<this_AllFatJets.size(); i++){
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
  }

//  FillHist("Nfatjet_hn_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
//  FillHist("Njet_hn_"+IDsuffix, jets.size(), weight, 8, 0., 8.); 

  std::vector<Lepton*> leptons, leptons_minus, leptons_plus, leptons_veto;

  //=======================
  //==== Sort in pt-order
  //=======================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);

//  int Nbjet_deepcsv_m=0;
//  int Nbjet_deepcsv_m_noSF=0;
  int Nbjet_loose = 0, Nbjet_medium = 0;

  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) Nbjet_deepcsv_m++; // method for getting btag with SF applied to MC
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,false,0)) Nbjet_deepcsv_m_noSF++; // method for getting btag with no SF applied to MC
//    if(IsBTagged(jets.at(ij), Jet::CSVv2, Jet::Loose,false,0)) Nbjet_csvv2++;
    if(IsBTagged(jets_nolepveto.at(ij), Jet::DeepCSV, Jet::Loose, true, 0)) Nbjet_loose++; // For b-jet veto in SM CR
    if(IsBTagged(jets_nolepveto.at(ij), Jet::DeepCSV, Jet::Medium, true, 0)) Nbjet_medium++; // For b-jet veto in SR
  }

//  FillHist("Nbjet_loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
//  FillHist("Nbjet_medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);

  //===================================
  //==== Set up pTcone, lepton vector
  //===================================

  double MET = ev.GetMETVector().Pt();
  double Mt = 0.;
  double ST = 0.;
  double MET2ST = 0.;
  double MZ = 91.1876;
  double MW = 80.379;
  double muon_idsf = 1.;
  double muon_isosf = 1.;
  double ele_idsf = 1.;
  double ele_recosf = 1.;
  int lepton_veto_size = 0;
  double LeptonPtCut1 = 0., LeptonPtCut2 = 0.;
  Particle ZCand, Wtemp1, Wtemp2, WCand1, WCand2;
  Particle llj, l1j, l2j,  lljj, l1jj, l2jj, l1J, l2J;
  Particle WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp, Ztemp1, Ztemp2, Ztemp3, Ztemp4, ZCand1, ZCand2;
  int ossf_mass10 = 0;
  
  // Set tight_iso cut & calculate pTcone
  double mu_tight_iso = 0.15;
  if(IDsuffix == "HNV2") mu_tight_iso = 0.1;
  if(IDsuffix == "HN16") mu_tight_iso = 0.07;

  double el_tight_iso = 0.;

  // Set pTcone
  for(unsigned int i=0; i<muons.size(); i++){
    double this_ptcone = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso);
    muons.at(i).SetPtCone(this_ptcone);
  }
   
  for(unsigned int i=0; i<electrons.size(); i++){
    el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt();
    if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
    if(IDsuffix == "HNV2"){
      el_tight_iso = std::min(0.08, 0.0287+0.506/electrons.at(i).UncorrPt());
      if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/electrons.at(i).UncorrPt());
    } 
    if(IDsuffix == "HN16") el_tight_iso = 0.08;
    double this_ptcone = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
    electrons.at(i).SetPtCone(this_ptcone);
  }

//  if(electrons.size() > 0) cout << electrons.at(0).PtCone() << endl;

  if(RunCF){
    if(electrons.size() == 2) electrons = ShiftElectronEnergy(electrons, param, true); 
  }

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
  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();
  ST += MET;
  MET2ST = MET*MET/ST;

  //=========================
  //==== Event selections..
  //=========================

  // Loop for each channel (mumu, ee)
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    if(it_ch == 0){ LeptonPtCut1 = MuonPtCut1; LeptonPtCut2 = MuonPtCut2; }
    if(it_ch == 1){ LeptonPtCut1 = ElectronPtCut1; LeptonPtCut2 = ElectronPtCut2; }
    if(it_ch==0 && RunCF) continue;

    // Triggers for each channel
    if(it_ch==0 && !ev.PassTrigger(MuonTriggers)) continue;
    if(it_ch==1 && !ev.PassTrigger(ElectronTriggers)) continue;

    // Cutflow : dilepton triggers
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max); 

    if(leptons.size() == 2){ 
      if(it_ch == 0){ if(!(muons.size()==2 && electrons.size()==0)) continue; }
      if(it_ch == 1){ if(!(muons.size()==0 && electrons.size()==2)) continue; }

      ZCand = *leptons.at(0) + *leptons.at(1);

      weight = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
      // weights for MC
      if(!IsDATA){
        // Gen matching with dR < 0.1
        Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        if(truth_lep1.PID() == 0) continue;
        if(truth_lep2.PID() == 0) continue;
 
        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          weight *= muon_idsf*muon_isosf;
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }

      // weight for fake, CF
      if(RunFake){
        weight = fakeEst->GetWeight(leptons, param);
//        FillHist("weight_FR", weight, 1., 220, -1.1, 1.1);
      }
      if(RunCF) weight = GetCFweight(leptons, param, true, 0);

      /////////////////////////////////////////////////////////
      //// Preselection (triggers have been already applied.)
      /////////////////////////////////////////////////////////

      if(!(leptons.at(0)->Pt()>LeptonPtCut1 && leptons.at(1)->Pt()>LeptonPtCut2)) continue;
      
      // Cutflow : 2 tight leptons (gen-matched, pT larger than trigger thresholds)
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(!RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) continue;
      if(RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) continue;

      // Cutflow : same-sign (oppsite-sign when RunCF=true)
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto 3rd leptons using veto ID
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(ZCand.M() > 10.)) continue;
      if(it_ch==1 && IsOnZ(ZCand.M(), 10.)) continue;

      // Cutflow : m(ll) > 10 GeV, |m(ll)-m(Z)| > 10 GeV for ee 
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);     

      // Lepton categories for the fake background
      if(it_ch == 0){
        if(muons.at(0).PassID(param.Muon_Tight_ID) && !(muons.at(1).PassID(param.Muon_Tight_ID))) LepCategory = "TL";
        if(!(muons.at(0).PassID(param.Muon_Tight_ID)) && muons.at(1).PassID(param.Muon_Tight_ID)) LepCategory = "LT";
        if(!(muons.at(1).PassID(param.Muon_Tight_ID)) && !(muons.at(0).PassID(param.Muon_Tight_ID))) LepCategory = "LL";
      }
      if(it_ch == 1){
        if(electrons.at(0).PassID(param.Electron_Tight_ID) && !(electrons.at(1).PassID(param.Electron_Tight_ID))) LepCategory = "TL";
        if(!(electrons.at(0).PassID(param.Electron_Tight_ID)) && electrons.at(1).PassID(param.Electron_Tight_ID)) LepCategory = "LT";
        if(!(electrons.at(0).PassID(param.Electron_Tight_ID)) && !(electrons.at(1).PassID(param.Electron_Tight_ID))) LepCategory = "LL";
      }

      // non-prompt CR2 : no jets && same-sign back-to-back 2 leptons
      if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){
       
        // Cutflow : jet requirement for non-prompt CR2 
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5)) continue;
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/fakeCR2/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

        if(RunFake){
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/fakeCR2/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
        }
      }

      if(!(jets.size()+fatjets.size() >= 1)) continue; 

      
      FillHist(channels.at(it_ch)+"/"+"Pre/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
      FillHist(channels.at(it_ch)+"/"+"Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/"+"Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/"+"Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/"+"Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      if(RunFake){
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/ZCand_Pt_unweighted_"+IDsuffix, ZCand.Pt(), 1., 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);  
        FillHist(channels.at(it_ch)+"/"+"Pre/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
      }

      // Event selections for each CR
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

        // Cutflow : jet requirement (Number or events at preselection)
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

        // non-prompt CR1 : SS 2 leptons with b-tagged jets
        if(it_rg == 0){
          if(!(Nbjet_medium > 0)) continue;
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_unweighted_"+IDsuffix, ZCand.Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
          }
        }
      
        // Low mass SR1, CR1 && High mass SR1, CR1
        if(it_rg>=1 && it_rg<5){
          if(!(jets.size()>=2 && fatjets.size()==0)) continue;

          // Select two jets that makes m(lljj), m(jj) closest to m(W)
          double tmpMassDiff1 = 10000., tmpMassDiff2 = 10000.;
          int j1 = 0, j2 = 0, j3 = 0, j4 = 0;
          for(unsigned int k=0; k<jets.size(); k++){
            for(unsigned int l=k+1; l<jets.size(); l++){
              Wtemp1 = *leptons.at(0) + *leptons.at(1) + jets.at(k) + jets.at(l);
              Wtemp2 = jets.at(k) + jets.at(l);
              if(fabs(Wtemp1.M() - MW) < tmpMassDiff1){
                tmpMassDiff1 = fabs(Wtemp1.M() - MW); 
                j1 = k; j2 = l;
              }
              if(fabs(Wtemp2.M() - MW) < tmpMassDiff2){
                tmpMassDiff2 = fabs(Wtemp2.M() - MW);
                j3 = k; j4 = l;
              }
            }
          }
          WCand1 = *leptons.at(0) + *leptons.at(1) + jets.at(j1) + jets.at(j2);
          WCand2 = jets.at(j3) + jets.at(j4);
          lljj = *leptons.at(0) + *leptons.at(1) + jets.at(j3) + jets.at(j4);
          l1jj = *leptons.at(0) + jets.at(j3) + jets.at(j4);
          l2jj = *leptons.at(1) + jets.at(j3) + jets.at(j4);

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max); 
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand1_Mass_nocut_"+IDsuffix, WCand1.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand2_Mass_nocut_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljj_Mass_nocut_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jj_Mass_nocut_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jj_Mass_nocut_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand1_Mass_nocut_"+IDsuffix, WCand1.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand2_Mass_nocut_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/lljj_Mass_nocut_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1jj_Mass_nocut_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2jj_Mass_nocut_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand1_Mass_nocut_unweighted_"+IDsuffix, WCand1.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand2_Mass_nocut_unweighted_"+IDsuffix, WCand2.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/lljj_Mass_nocut_unweighted_"+IDsuffix, lljj.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1jj_Mass_nocut_unweighted_"+IDsuffix, l1jj.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2jj_Mass_nocut_unweighted_"+IDsuffix, l2jj.M(), 1., 2000, 0., 2000.);
          }

          // Low mass SR1 
          if(it_rg == 1){
            if(!(Nbjet_medium == 0)) continue;
            if(!(WCand1.M() < 300.)) continue;
            if(!(MET < 80.)) continue;
          }

          // Low mass CR1
          if(it_rg == 2){
            if(!(WCand1.M() < 300.)) continue;
            if(!(Nbjet_medium>0 || MET>100.)) continue;
          }

          // High mass SR1
          if(it_rg == 3){
            if(!(Nbjet_medium == 0)) continue;
            if(!(WCand2.M() < 150.)) continue;
            if(!(MET2ST < 15.)) continue;
          }

          // High mass CR1
          if(it_rg == 4){
            if(!(WCand2.M() < 150.)) continue;
            if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand1_Mass_"+IDsuffix, WCand1.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand2_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand1_Mass_"+IDsuffix, WCand1.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand2_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_unweighted_"+IDsuffix, ZCand.Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand1_Mass_unweighted_"+IDsuffix, WCand1.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/WCand2_Mass_unweighted_"+IDsuffix, WCand2.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/lljj_Mass_unweighted_"+IDsuffix, lljj.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1jj_Mass_unweighted_"+IDsuffix, l1jj.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2jj_Mass_unweighted_"+IDsuffix, l2jj.M(), 1., 2000, 0., 2000.);
          }
        }

        // Low mass SR2, CR2
        if(it_rg>=5 && it_rg<7){
          if(!(jets.size()==1 && fatjets.size()==0)) continue;
          llj = *leptons.at(0) + *leptons.at(1) + jets.at(0);
          l1j = *leptons.at(0) + jets.at(0);
          l2j = *leptons.at(1) + jets.at(0);

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/llj_Mass_nocut_"+IDsuffix, llj.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1j_Mass_nocut_"+IDsuffix, l1j.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2j_Mass_nocut_"+IDsuffix, l2j.M(), weight, 1000, 0., 1000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/llj_Mass_nocut_"+IDsuffix, llj.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1j_Mass_nocut_"+IDsuffix, l1j.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2j_Mass_nocut_"+IDsuffix, l2j.M(), weight, 1000, 0., 1000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/llj_Mass_nocut_unweighted_"+IDsuffix, llj.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1j_Mass_nocut_unweighted_"+IDsuffix, l1j.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2j_Mass_nocut_unweighted_"+IDsuffix, l2j.M(), 1., 1000, 0., 1000.);
          }

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
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/llj_Mass_"+IDsuffix, llj.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1j_Mass_"+IDsuffix, l1j.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2j_Mass_"+IDsuffix, l2j.M(), weight, 1000, 0., 1000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/llj_Mass_"+IDsuffix, llj.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1j_Mass_"+IDsuffix, l1j.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2j_Mass_"+IDsuffix, l2j.M(), weight, 1000, 0., 1000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_unweighted_"+IDsuffix, ZCand.Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/llj_Mass_unweighted_"+IDsuffix, llj.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1j_Mass_unweighted_"+IDsuffix, l1j.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2j_Mass_unweighted_"+IDsuffix, l2j.M(), 1., 1000, 0., 1000.);
          }
        }

        // High mass SR2, CR2
        if(it_rg >= 7){
          if(!(fatjets.size() > 0)) continue;
          double tmpMassDiff3 = 10000.;
          int j5 = 0;
          for(unsigned int k=0; k<fatjets.size(); k++){
            if(fabs(fatjets.at(k).M() - MW) < tmpMassDiff3){
              tmpMassDiff3 = fabs(fatjets.at(k).SDMass() - MW);
              j5 = k;
            }
          }
//          WCand3 = fatjets.at(j5);
          l1J = *leptons.at(0) + fatjets.at(j5);
          l2J = *leptons.at(1) + fatjets.at(j5);

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Mass_nocut_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1J_Mass_nocut_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2J_Mass_nocut_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Fatjet_Mass_nocut_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1J_Mass_nocut_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2J_Mass_nocut_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Fatjet_Mass_nocut_unweighted_"+IDsuffix, fatjets.at(j5).SDMass(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1J_Mass_nocut_unweighted_"+IDsuffix, l1J.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2J_Mass_nocut_unweighted_"+IDsuffix, l2J.M(), 1., 2000, 0., 2000.);
          }

          // High mass SR2
          if(it_rg == 7){
            if(!(Nbjet_medium == 0)) continue;
            if(!(fatjets.at(j5).SDMass() < 150.)) continue;
            if(!(MET2ST < 15.)) continue;
          }

          // High mass CR2
          if(it_rg == 8){
            if(!(fatjets.at(j5).SDMass() < 150.)) continue;
            if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Fatjet_Mass_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2J_Mass_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

          if(RunFake){
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Fatjet_Mass_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2J_Mass_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_Jets_unweighted_"+IDsuffix, jets.size(), 1., 8, 0., 8.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Loose_unweighted_"+IDsuffix, Nbjet_loose, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_BJets_Medium_unweighted_"+IDsuffix, Nbjet_medium, 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Number_FatJets_unweighted_"+IDsuffix, fatjets.size(), 1., 5, 0., 5.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Mass_unweighted_"+IDsuffix, ZCand.M(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/ZCand_Pt_unweighted_"+IDsuffix, ZCand.Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Pt_unweighted_"+IDsuffix, leptons.at(0)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Pt_unweighted_"+IDsuffix, leptons.at(1)->Pt(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep1_Eta_unweighted_"+IDsuffix, leptons.at(0)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Lep2_Eta_unweighted_"+IDsuffix, leptons.at(1)->Eta(), 1., 50, -2.5, 2.5);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET_unweighted_"+IDsuffix, MET, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/MET2ST_unweighted_"+IDsuffix, MET2ST, 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/Fatjet_Mass_unweighted_"+IDsuffix, fatjets.at(j5).SDMass(), 1., 1000, 0., 1000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l1J_Mass_unweighted_"+IDsuffix, l1J.M(), 1., 2000, 0., 2000.);
            FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/"+LepCategory+"/l2J_Mass_unweighted_"+IDsuffix, l2J.M(), 1., 2000, 0., 2000.);
          }
        }
      }  
    }
  }

  //////////////////////////////////////
  //// SM background CR
  //////////////////////////////////////
  if(RunCF) return;

  if(IsDATA){
    if(DataStream.Contains("DoubleMuon") && !ev.PassTrigger(MuonTriggers)) return;
    if(DataYear==2016 || DataYear==2017){
      if(DataStream.Contains("DoubleEG")){
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
    if(DataYear==2018){
      if(DataStream.Contains("EGamma")){
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
  }

/*  weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);
  }

  for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
    FillHist("Number_Events_"+regionsSM.at(it_rg2)+"_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
  }*/ 

  //=========================
  //==== Event selections..
  //========================= 
  
  for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
    weight = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
    ossf_mass10 = 0;

    // Cutflow : passing dilepton triggers (dimu || diel)
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);

    // Requirements after passing triggers
    if(ev.PassTrigger(MuonTriggers)){
      if(muons.size() < 2) continue;
      if(muons.size()>=2 && !(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2)) continue;
    }
    if(!ev.PassTrigger(MuonTriggers) && ev.PassTrigger(ElectronTriggers)){
      if(electrons.size() < 2) continue;
      if(electrons.size()>=2 && !(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2)) continue;
    }

    // WZ, ZG
    if(leptons.size()==3 && it_rg2<2){
   
      // weights for MC 
      if(!IsDATA){
        Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        Gen truth_lep3 = GetGenMatchedLepton(*leptons.at(2), gens);
        if(truth_lep1.PID() == 0) continue;
        if(truth_lep2.PID() == 0) continue;
        if(truth_lep3.PID() == 0) continue;

        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
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

      // Cutflow : 3 tight leptons (gen-matched, pT larger than trigger thresholds)
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      int l1 = -999, l2 = -999, wlep = -999;
      // OSSF lepton pair, W-tagged lepton
      if(muons.size()==2 && muons.at(0).Charge()*muons.at(1).Charge()<0){
        ZCand = muons.at(0) + muons.at(1);
        WtagLep = electrons.at(0);
        ZtagLep1 = muons.at(0);
        ZtagLep2 = muons.at(1);
      }
      else if(electrons.size()==2 && electrons.at(0).Charge()*electrons.at(1).Charge()<0){
        ZCand = electrons.at(0) + electrons.at(1);
        WtagLep = muons.at(0);
        ZtagLep1 = electrons.at(0);
        ZtagLep2 = electrons.at(1);
      }
      else if(muons.size()==3 || electrons.size()==3){
        if(fabs(leptons.at(0)->Charge() + leptons.at(1)->Charge() + leptons.at(2)->Charge()) == 1){
          double tmpMassDiff = 1000000.;
 
          for(int ilep1=0; ilep1<2; ilep1++){
            for(int ilep2=ilep1+1; ilep2<3; ilep2++){
              if(leptons.at(ilep1)->Charge()*leptons.at(ilep2)->Charge()>0) continue;
              Ztemp = *leptons.at(ilep1) + *leptons.at(ilep2);
              if(!(Ztemp.M() > 10.)) ossf_mass10++;
              if(fabs(Ztemp.M() - MZ) < tmpMassDiff) {
                  tmpMassDiff= fabs(Ztemp.M() - MZ);
                  ZCand = Ztemp; l1 = ilep1; l2 = ilep2;
              }
            }
          }
          for(int ilep3=0; ilep3<3; ilep3++){
            if(fabs(ilep3-l1)>0 && fabs(ilep3-l2)>0) { wlep = ilep3; break;}
          }
          WtagLep = *leptons.at(wlep);
          ZtagLep1 = *leptons.at(l1);
          ZtagLep2 = *leptons.at(l2);

        }
      } 
      else continue;

      TriLep = ZCand + WtagLep;
      Mt = MT(WtagLep, METv);

      if(!(ossf_mass10 == 0)) continue;
      
      // Cutflow : m(ll) > 10 GeV
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(Nbjet_loose == 0)) continue;

      // Cutflow : No b jets
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(it_rg2 == 0){
        if(!IsOnZ(ZCand.M(), 15.)) continue;
        if(!(MET > 50.)) continue;
        if(!(Mt > 20.)) continue;
        if(!(TriLep.M() > MZ + 15.)) continue;
      }
      if(it_rg2 == 1){
        if(IsOnZ(ZCand.M(), 15.)) continue;
        if(!(MET < 50.)) continue;
        if(!IsOnZ(TriLep.M(), 15.)) continue;
      }

      // weights for MC
/*      if(!IsDATA){
        Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        if(truth_lep1.PID() == 0) continue;
        if(truth_lep2.PID() == 0) continue;

        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          weight *= muon_idsf*muon_isosf;
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        } 
      }
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);*/

      // Histograms
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.); 
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

      for(unsigned int it_ch2=0; it_ch2<channels3L.size(); it_ch2++){
        if(muons.size() == it_ch2){
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        }
      } 
    }

    // ZZ 
    if(leptons.size()==4 && it_rg2==2){

      // weights for MC
      if(!IsDATA){
        Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        Gen truth_lep3 = GetGenMatchedLepton(*leptons.at(2), gens);
        Gen truth_lep4 = GetGenMatchedLepton(*leptons.at(3), gens);
        if(truth_lep1.PID() == 0) continue;
        if(truth_lep2.PID() == 0) continue;
        if(truth_lep3.PID() == 0) continue;
        if(truth_lep4.PID() == 0) continue;

        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
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

      // Cutflow : 4 tight leptons (gen-matched, pT larger than trigger thresholds)
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      // OSSF lepton pairs
      if(muons.size()==2 && muons.at(0).Charge()*muons.at(1).Charge()<0 && electrons.at(0).Charge()*electrons.at(1).Charge()<0){
        ZCand1 = muons.at(0) + muons.at(1);
        ZCand2 = electrons.at(0) + electrons.at(1);
      }
      else if(muons.size()==4 || electrons.size()==4){
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
      else continue;

      if(!(ossf_mass10 == 0)) continue;

      // Cutflow : m(ll) > 10 GeV
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(Nbjet_loose == 0)) continue;

      // Cutflow : No b jets
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(IsOnZ(ZCand1.M(), 15.) && IsOnZ(ZCand2.M(), 15.))) continue;

      // weights for MC
/*      if(!IsDATA){
        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          weight *= muon_idsf*muon_isosf;
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);*/

      // Histograms 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand1_Mass_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand2_Mass_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand1_Pt_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand2_Pt_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep4_Pt_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep4_Eta_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      for(unsigned int it_ch2=0; it_ch2<channels4L.size(); it_ch2++){
        if(muons.size() == 2*it_ch2){
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 8, 0., 8.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand1_Mass_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand2_Mass_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand1_Pt_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand2_Pt_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep4_Pt_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep4_Eta_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        } 
      }
    }
  } 
}




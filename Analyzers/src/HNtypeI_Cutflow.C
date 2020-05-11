#include "HNtypeI_Cutflow.h"

HNtypeI_Cutflow::HNtypeI_Cutflow(){

}

void HNtypeI_Cutflow::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[HNtypeI_Cutflow::initializeAnalyzer] RunSyst = " << RunSyst << endl;
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
  MuonTightIDs = {"HNTight2016"};
  MuonLooseIDs = {"HNLoose2016"};
  MuonVetoIDs  = {"HNVeto2016"};
  ElectronTightIDs = {"HNTight2016"};
  ElectronLooseIDs = {"HNLoose2016"};
  ElectronVetoIDs  = {"HNVeto2016"};
  FakeRateIDs = {"HNtypeI_16"};


  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_Cutflow.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  // We can use SkimTree_SMP for dimuon and dielectron channel. (See https://github.com/sansan9401/SKFlatAnalyzer/blob/Run2Legacy_hsseo/Analyzers/src/SkimTree_SMP.C)

  MuonTriggers.clear();
  ElectronTriggers.clear();

  if(DataYear==2016){
//    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
//    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
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

//  cout << "[HNtypeI_Cutflow::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[HNtypeI_Cutflow::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_Cutflow::~HNtypeI_Cutflow(){

  //==== Destructor of this Analyzer

}

void HNtypeI_Cutflow::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_Cutflow.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_Cutflow.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_EleID=0; it_EleID<ElectronTightIDs.size(); it_EleID++){
    TString MuonTightID     = MuonTightIDs.at(it_EleID);
    TString MuonLooseID     = MuonLooseIDs.at(it_EleID); 
    TString MuonVetoID      = MuonVetoIDs.at(it_EleID);
    TString ElectronTightID = ElectronTightIDs.at(it_EleID);
    TString ElectronLooseID = ElectronLooseIDs.at(it_EleID);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_EleID);
    TString FakeRateID      = FakeRateIDs.at(it_EleID);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;

//      param.Name = MuonID+"_"+"Central";

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
//    param.Jet_ID    = "HNTight";
    param.Jet_ID    = "tight";
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

void HNtypeI_Cutflow::executeEventFromParameter(AnalyzerParameter param){

  TString IDsuffix = "HN16";
//  if(param.Electron_Tight_ID.Contains("V2")) IDsuffix = "HNV2";
//  if(param.Electron_Tight_ID.Contains("2016")) IDsuffix = "HN16";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;

  Event ev = GetEvent();

  //================================================================
  //==== No Cut
  //================================================================

  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);
  }

  FillHist("Number_Events_"+IDsuffix, -0.5, weight, cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_unweighted_"+IDsuffix, -0.5, 1., cutflow_bin, 0., cutflow_max);

  //================================================================
  //==== MET Filter
  //================================================================

  if(!PassMETFilter()) return;
//  FillHist("Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
//  FillHist("Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);

  //================================================================
  //==== Trigger
  //================================================================

  if(!(ev.PassTrigger(MuonTriggers))) return;

  FillHist("Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);

  //================================================================
  //==== Copy AllObjects
  //================================================================

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
    cout << "[HNtypeI_Cutflow::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //================================================================
  //==== Then, apply ID selections using this_AllXXX
  //================================================================

  TString MuonID = param.Muon_Tight_ID;
  TString ElectronID = param.Electron_Tight_ID;

  vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 10., 2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);  // AK4jets used for b tagging
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  vector<Muon> muons_prompt;
  muons_prompt.clear();

  // Jet selection to avoid double counting due to jets matched geometrically with a lepton
  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
//  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
//  vector<Jet> jets = JetsPassPileupMVA(jets_lepveto); // We use only AK4jets here
  //==== For test jet selection
//  vector<Jet> jets = jets_nolepveto;
//  vector<Jet> jets = GetJets("tight", 20., 2.7);
//  vector<Jet> jets = jets_lepveto;
  vector<Jet> jets = JetsPassPileupMVA(jets_lepveto);
//  vector<Jet> jets = JetsAwayFromFatJet(jets_PUveto, fatjets);  // AK4jets used in SR, CR

  //================================================================
  //==== Sort in pT-order
  //================================================================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);

  //================================================================
  //==== b tagging
  //================================================================

  int Nbjet_loose = 0, Nbjet_medium = 0, Nbjet_lepveto_medium = 0;
  JetTagging::Parameters jtp_DeepCSV_Loose  = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
//  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int i=0; i<jets_nolepveto.size(); i++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(i))) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(i))) Nbjet_medium++;
  }

  for(unsigned int i=0; i<jets.size(); i++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets.at(i))) Nbjet_lepveto_medium++;
  }

  FillHist("Number_Jets_nolepveto_Nocut_"+IDsuffix, jets_nolepveto.size(), weight, 10, 0., 10.); 
  FillHist("Number_Jets_lepveto_Nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
  FillHist("Number_BJets_nolepveto_Nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
  FillHist("Number_BJets_lepveto_Nocut_"+IDsuffix, Nbjet_lepveto_medium, weight, 10, 0., 10.);

  //================================================================
  //==== Set up MET
  //================================================================

  Particle METv = ev.GetMETVector();

  // MET correction
  METv = UpdateMETMuon(METv, this_AllMuons);           // Rochester momentum correction
  METv = UpdateMETElectron(METv, this_AllElectrons);   // electron pT in miniAOD = UncorrPt

  double MET = METv.Pt();

  //================================================================
  //==== Set up particles, variables
  //================================================================

//  double MZ = 91.1876;
  double MW = 80.379;
  Particle ZCand, WCand, Wtemp;

  //================================================================
  //==== Event selections
  //================================================================
  //TODO : Add twiki here!!

  if(muons.size() == 2){
    if(!(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2)) return;

    // Truth matching
    muons_prompt = MuonPromptOnlyHNtypeI(muons, gens); 
    if(muons_prompt.size() < 2) return;

/*    if(!IsDATA){
      Gen truth_lep1 = GetGenMatchedLepton(muons.at(0), gens);
      Gen truth_lep2 = GetGenMatchedLepton(muons.at(1), gens);
      if(truth_lep1.PID() == 0) return;
      if(truth_lep2.PID() == 0) return;
    }*/

    FillHist("Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max); 
    FillHist("Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);

/*    if(!(muons.at(0).Charge() == muons.at(1).Charge())) return;
    FillHist("Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max); 
    FillHist("Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);*/

    if(!(muons_veto.size() == 2)) return;
    FillHist("Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(electrons_veto.size() == 0)) return;
    FillHist("Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

    ZCand = muons.at(0) + muons.at(1);
    if(!(ZCand.M() > 10.)) return;
    FillHist("Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

    FillHist("Njet_"+IDsuffix, jets.size(), 1., 10, 0., 10.);

    if(!(jets.size() >= 2)) return;
    FillHist("Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(MET < 50.)) return;
    FillHist("Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

/*    if(!(fatjets.size()>0) && !(jets.size()>1 && fatjets.size()==0) && !(jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.)) return;
    FillHist("Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(Nbjet_lepveto_medium == 0)) return;
    FillHist("Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(Nbjet_medium == 0)) return;
    FillHist("Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);*/

    double tmpMassDiff = 10000.;
    int j1 = 0, j2 = 0;
    for(unsigned int k=0; k<jets.size(); k++){
      for(unsigned int l=k+1; l<jets.size(); l++){
        Wtemp = jets.at(k) + jets.at(l);
        if(fabs(Wtemp.M() - MW) < tmpMassDiff){
          tmpMassDiff = fabs(Wtemp.M() - MW);
          j1 = k; j2 = l;
        }
      }
    }

    WCand = jets.at(j1) + jets.at(j2);

    if(!(WCand.M() < 200.)) return;
    FillHist("Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

    FillHist("Number_BJets_nolepveto_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
    FillHist("Number_BJets_lepveto_"+IDsuffix, Nbjet_lepveto_medium, weight, 10, 0., 10.);

    if(!(Nbjet_lepveto_medium == 0)) return;
    FillHist("Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(Nbjet_medium == 0)) return;
    FillHist("Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
  }
}

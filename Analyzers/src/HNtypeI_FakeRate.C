#include "HNtypeI_FakeRate.h"

HNtypeI_FakeRate::HNtypeI_FakeRate(){

}

void HNtypeI_FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  MuonTightIDs     = {"HNTight2016", "ISRTight", "HNTightV1", "HNTightV2", "HNTightV1"};
  MuonLooseIDs     = {"HNLoose2016", "ISRLoose", "HNLooseV1", "HNLooseV1", "HNLooseV2"};
  MuonVetoIDs      = {"HNVeto2016", "ISRVeto", "HNVeto", "HNVeto", "HNVeto"};
  ElectronTightIDs = {"HNTight2016", "ISRTight", "HNTightV1", "HNTightV2", "HNMVATight"};
  ElectronLooseIDs = {"HNLoose2016", "ISRLoose", "HNLooseV1", "HNLooseV1", "HNMVALoose"};
  ElectronVetoIDs  = {"HNVeto2016", "ISRVeto", "HNVeto", "HNVeto", "HNMVAVeto"};

  /*MuonTightIDs = {"HNTightV2"};
  MuonLooseIDs = {"HNLoose"};
  ElectronTightIDs = {"HNTight2016"};
  ElectronLooseIDs = {"HNLoose2016"};*/

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_FakeRate.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
  MuonTriggers.clear();
  ElectronTriggers.clear();

  MuonTrig1 = "HLT_Mu3_PFJet40_v";       // DoubleMuon(2016), SingleMuon(2017,2018)
  MuonTrig2 = "HLT_Mu8_TrkIsoVVL_v";     // DoubleMuon
  MuonTrig3 = "HLT_Mu17_TrkIsoVVL_v";    // DoubleMuon

  MuonTriggers.push_back(MuonTrig1);      
  MuonTriggers.push_back(MuonTrig2);    
  MuonTriggers.push_back(MuonTrig3);   
  MuonPtCut1 = 5., MuonPtCut2 = 10., MuonPtCut3 = 20.;
  MuonPtconeCut1 = 5., MuonPtconeCut2 = 20., MuonPtconeCut3 = 30.;

  // DoubleEG (2016), SingleElectron (2017), EGamma (2018)
  ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; 
  ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  ElectronTrig3 = "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v";
  ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

  ElectronTriggers.push_back(ElectronTrig1);
  ElectronTriggers.push_back(ElectronTrig2);
  ElectronTriggers.push_back(ElectronTrig3);
  ElectronTriggers.push_back(ElectronTrig4);
  ElectronPtCut1 = 9.5, ElectronPtCut2 = 15., ElectronPtCut3 = 20., ElectronPtCut4 = 25.;
  ElectronPtconeCut1 = 10., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 40.;

  // luminosity of prescaled triggers
  // Without normalizationfactor
  /*if(DataYear==2016){
    MuonLumi1 = 7.408, MuonLumi2 = 7.801, MuonLumi3 = 216.748;
    ElectronLumi1 = 6.988, ElectronLumi2 = 14.851 , ElectronLumi3 = 62.761, ElectronLumi4 = 62.808;  // Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v : 58.639
  }
  if(DataYear==2017){
    MuonLumi1 = 4.612, MuonLumi2 = 2.903, MuonLumi3 = 65.943;
    ElectronLumi1 = 3.973, ElectronLumi2 = 27.698, ElectronLumi3 = 35.594, ElectronLumi4 = 43.468;
  }
  if(DataYear==2018){
    MuonLumi1 = 2.696, MuonLumi2 = 8.561, MuonLumi3 = 45.781;
    ElectronLumi1 = 6.412, ElectronLumi2 = 38.849, ElectronLumi3 = 38.861, ElectronLumi4 = 38.906;
  }

  // With normalization factor
  if(DataYear==2016){
    MuonLumi1 = 7.408*0.680163, MuonLumi2 = 7.801*1.22923, MuonLumi3 = 216.748*0.912288;
    ElectronLumi1 = 6.988*1.05547, ElectronLumi2 = 14.851*0.972706 , ElectronLumi3 = 62.761*0.929557, ElectronLumi4 = 62.808*0.928212;  // Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v : 58.639
  }
  if(DataYear==2017){
    MuonLumi1 = 4.612*1.01916, MuonLumi2 = 2.903*1.22549, MuonLumi3 = 65.943*0.905105;
    ElectronLumi1 = 3.973*0.962824, ElectronLumi2 = 27.698*0.846783, ElectronLumi3 = 35.594*0.792109, ElectronLumi4 = 43.468*0.767242;
  }
  if(DataYear==2018){
    MuonLumi1 = 2.696*2.05535, MuonLumi2 = 8.561*1.03262, MuonLumi3 = 45.781*0.901847;
    ElectronLumi1 = 6.412*0.89131, ElectronLumi2 = 38.849*0.92211, ElectronLumi3 = 38.861*0.79613, ElectronLumi4 = 38.906*0.791474;
  }*/

  //cout << "[HNtypeI_FakeRate::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_FakeRate::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_FakeRate::~HNtypeI_FakeRate(){

  //==== Destructor of this Analyzer

}

void HNtypeI_FakeRate::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_FakeRate.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_FakeRate.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){

    /*if(it_EleID < 3){
      if(it_EleID != it_MuonID) continue;
    }
    if(it_EleID >= 3){
      if(it_MuonID != 0) continue;  // See NO_MUON_EVENT
    }*/

    //TString MuonID = "HNTight2016";
    //TString MuonIDSFKey = "NUM_TightID_DEN_genTracks";
    TString MuonTightID = MuonTightIDs.at(it_id);
    TString MuonLooseID = MuonLooseIDs.at(it_id);
    TString MuonVetoID  = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);

    param.Clear();

    param.fakesyst_ = AnalyzerParameter::FakeCentral;

    //param.Name = MuonID+"_"+"Central";
    param.Name = "FakeCentral";

    // Muon ID
    param.Muon_Tight_ID = MuonTightID;
    param.Muon_Loose_ID = MuonLooseID;
    param.Muon_Veto_ID  = MuonVetoID;
    param.Muon_ID_SF_Key = "";
    param.Muon_ISO_SF_Key = "";

    // Electron ID
    param.Electron_Tight_ID = ElectronTightID;
    param.Electron_Loose_ID = ElectronLooseID;
    param.Electron_Veto_ID  = ElectronVetoID;
    param.Electron_ID_SF_Key = "";

    // Jet ID
    param.Jet_ID = "HNTight";

    executeEventFromParameter(param);

    /*if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NFakeSyst; it_syst++){
        param.fakesyst_ = AnalyzerParameter::FakeSyst(it_syst);
        //param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        param.Name  = "FakeSyst_"+param.GetFakeSystType();
        executeEventFromParameter(param);
      }
    }*/

  }

}

void HNtypeI_FakeRate::executeEventFromParameter(AnalyzerParameter param){
 
  // ID version
  vector<TString> regions_mu = {"muonFR_2016", "muonDY_2016", "muonWJ_2016"};                               // 2016 ID
  if(param.Muon_Tight_ID.Contains("ISR")) regions_mu = {"muonFR_ISR", "muonDY_ISR", "muonWJ_ISR"};          // ISR ID
  if(param.Muon_Tight_ID.Contains("HNTightV1")) regions_mu = {"muonFR_V1", "muonDY_V1", "muonWJ_V1"};       // HN ID V1 based on POG cut-based ID
  if(param.Muon_Tight_ID.Contains("HNTightV2")) regions_mu = {"muonFR_V2", "muonDY_V2", "muonWJ_V2"};       // HN ID V2 based on POG cut-based ID
  if(param.Muon_Loose_ID.Contains("HNLooseV2")) regions_mu = {"muonFR_MVA", "muonDY_MVA", "muonWJ_MVA"};    // HN ID based on POG cut-based ID (Loose V2)

  vector<TString> regions_el = {"eleFR_2016", "eleDY_2016", "eleWJ_2016"};                                  // 2016 ID
  if(param.Electron_Tight_ID.Contains("ISR")) regions_el = {"eleFR_ISR", "eleDY_ISR", "eleWJ_ISR"};         // ISR ID
  if(param.Electron_Tight_ID.Contains("HNTightV1")) regions_el = {"eleFR_V1", "eleDY_V1", "eleWJ_V1"};      // HN ID V1 based on POG cut-based ID
  if(param.Electron_Tight_ID.Contains("HNTightV2")) regions_el = {"eleFR_V2", "eleDY_V2", "eleWJ_V2"};      // HN ID V2 based on POG cut-based ID
  if(param.Electron_Tight_ID.Contains("HNMVA")) regions_el = {"eleFR_MVA", "eleDY_MVA", "eleWJ_MVA"};       // HN ID based on POG MVA ID

  TString systName = param.Name;

  //========================================================
  //==== Luminosity of prescaled triggers
  //========================================================

  TString IsNorm = true;

  // With normalization factor

  // Obsolete
  /*if(DataYear==2016){
    MuonLumi1 = 7.408*0.680163, MuonLumi2 = 7.801*1.22923, MuonLumi3 = 216.748*0.912288;
    ElectronLumi1 = 6.988*1.05547, ElectronLumi2 = 14.851*0.972706 , ElectronLumi3 = 62.761*0.929557, ElectronLumi4 = 62.808*0.928212;  // Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v : 58.639
  }
  if(DataYear==2017){
    MuonLumi1 = 4.612*1.01916, MuonLumi2 = 2.903*1.22549, MuonLumi3 = 65.943*0.905105;
    ElectronLumi1 = 3.973*0.962824, ElectronLumi2 = 27.698*0.846783, ElectronLumi3 = 35.594*0.792109, ElectronLumi4 = 43.468*0.767242;
  }
  if(DataYear==2018){
    MuonLumi1 = 2.696*2.05535, MuonLumi2 = 8.561*1.03262, MuonLumi3 = 45.781*0.901847;
    ElectronLumi1 = 6.412*0.89131, ElectronLumi2 = 38.849*0.92211, ElectronLumi3 = 38.861*0.79613, ElectronLumi4 = 38.906*0.791474;
  }*/

  // HN 2016 ID
  if(param.Muon_Tight_ID.Contains("2016")){
    if(DataYear==2016){
      MuonLumi1 = 7.408*0.679317, MuonLumi2 = 7.801*1.26668, MuonLumi3 = 216.748*0.940637;
      ElectronLumi1 = 6.988*1.04596, ElectronLumi2 = 14.851*0.968804, ElectronLumi3 = 62.761*0.933706, ElectronLumi4 = 62.808*0.933263;
    }
    if(DataYear==2017){
      MuonLumi1 = 4.612*1.14633, MuonLumi2 = 2.903*1.32612, MuonLumi3 = 65.943*0.981435;
      ElectronLumi1 = 3.973*1.05833, ElectronLumi2 = 27.698*0.913181, ElectronLumi3 = 35.594*0.854491, ElectronLumi4 = 43.468*0.821393;
    }
    if(DataYear==2018){
      MuonLumi1 = 2.696*1.89483, MuonLumi2 = 8.561*1.0494, MuonLumi3 = 45.781*0.915734;
      ElectronLumi1 = 6.412*0.932546, ElectronLumi2 = 38.849*0.928263, ElectronLumi3 = 38.861*0.797685, ElectronLumi4 = 38.906*0.793358;
    }
  }

  // ISR ID
  if(param.Muon_Tight_ID.Contains("ISR")){
    if(DataYear==2016){
      MuonLumi1 = 7.408*0.697318, MuonLumi2 = 7.801*1.25839, MuonLumi3 = 216.748*0.933965;
      ElectronLumi1 = 6.988*1.09595, ElectronLumi2 = 14.851*1.02164, ElectronLumi3 = 62.761*0.975231, ElectronLumi4 = 62.808*0.974035;
    }
    if(DataYear==2017){
      MuonLumi1 = 4.612*1.09987, MuonLumi2 = 2.903*1.3113, MuonLumi3 = 65.943*0.969653;
      ElectronLumi1 = 3.973*1.06166, ElectronLumi2 = 27.698*0.928937, ElectronLumi3 = 35.594*0.86962, ElectronLumi4 = 43.468*0.841484;
    }
    if(DataYear==2018){
      MuonLumi1 = 2.696*1.92987, MuonLumi2 = 8.561*1.08561, MuonLumi3 = 45.781*0.949186;
      ElectronLumi1 = 6.412*0.998446, ElectronLumi2 = 38.849*0.99398, ElectronLumi3 = 38.861*0.85411, ElectronLumi4 = 38.906*0.848167;
    }
  }

  // HN POG ID
  if(param.Muon_Tight_ID.Contains("HNTight")){
    if(DataYear==2016){
      MuonLumi1 = 7.408*0.704597, MuonLumi2 = 7.801*1.24946, MuonLumi3 = 216.748*0.931095;
      ElectronLumi1 = 6.988*1.0641, ElectronLumi2 = 14.851*0.981272, ElectronLumi3 = 62.761*0.941996, ElectronLumi4 = 62.808*0.94131;  
    }
    if(DataYear==2017){
      MuonLumi1 = 4.612*1.11176, MuonLumi2 = 2.903*1.28788, MuonLumi3 = 65.943*0.965029;
      ElectronLumi1 = 3.973*1.01877, ElectronLumi2 = 27.698*0.886662, ElectronLumi3 = 35.594*0.825669, ElectronLumi4 = 43.468*0.795745;
    }
    if(DataYear==2018){
      MuonLumi1 = 2.696*1.95092, MuonLumi2 = 8.561*1.06968, MuonLumi3 = 45.781*0.943577;
      ElectronLumi1 = 6.412*0.922875, ElectronLumi2 = 38.849*0.939071, ElectronLumi3 = 38.861*0.809609, ElectronLumi4 = 38.906*0.804223;
    }
  }

  // Without normalization factors
  if(!IsNorm){
    if(DataYear==2016){
      MuonLumi1 = 7.408, MuonLumi2 = 7.801, MuonLumi3 = 216.748;
      ElectronLumi1 = 6.988, ElectronLumi2 = 14.851 , ElectronLumi3 = 62.761, ElectronLumi4 = 62.808;  // Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v : 58.639
    }
    if(DataYear==2017){
      MuonLumi1 = 4.612, MuonLumi2 = 2.903, MuonLumi3 = 65.943;
      ElectronLumi1 = 3.973, ElectronLumi2 = 27.698, ElectronLumi3 = 35.594, ElectronLumi4 = 43.468;
    }
    if(DataYear==2018){
      MuonLumi1 = 2.696, MuonLumi2 = 8.561, MuonLumi3 = 45.781;
      ElectronLumi1 = 6.412, ElectronLumi2 = 38.849, ElectronLumi3 = 38.861, ElectronLumi4 = 38.906;
    }
  }


  Event ev = GetEvent();
  ev.SetMET(pfMET_Type1_pt, pfMET_Type1_phi);

  //========================================================
  //==== No Cut
  //========================================================

  //JSFillHist(param.Name, "NoCut_"+param.Name, 0., 1., 1, 0., 1.);

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //========================================================
  //==== Trigger
  //========================================================

  //if(! (ev.PassTrigger(MuonTriggers) )) return;

  //========================================================
  //==== Copy AllObjects
  //========================================================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;
  vector<Gen> gens = GetGens();

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  double dphiCut = 2.5;
  //double jetPtCut_syst = 40.;
  double PtRatioCut = 1.;

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
    cout << "[HNtypeI_FakeRate::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  if(param.fakesyst_ == AnalyzerParameter::FakeCentral){

  }

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  //==== Leptons
  vector<Muon> muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  vector<Muon> muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, MuonPtCut1, 2.4);
  vector<Muon> muons_veto  = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Electron> electrons_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  vector<Electron> electrons_loose = SelectElectrons(this_AllElectrons, param.Electron_Loose_ID, ElectronPtCut1, 2.5);
  vector<Electron> electrons_veto  = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);

  //==== Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  //==== Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);
  vector<Jet> jets = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  //vector<Jet> jets_PUveto = JetsPassPileupMVA(jets_lepveto);
  
  vector<Jet> jets_awayFromMuon;
  vector<Jet> jets_awayFromElectron;
  jets_awayFromMuon.clear();
  jets_awayFromElectron.clear();

  //========================================================
  //==== Sort in pt-order
  //========================================================

  std::sort(muons_tight.begin(), muons_tight.end(), PtComparing);
  std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons_tight.begin(), electrons_tight.end(), PtComparing);
  std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);

  //========================================================
  //==== B-Tagging
  //========================================================

  int Nbjet_medium=0;
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv_central = ev.GetMETVector();

  /*if(muons_loose.size()==1 || electrons_loose.size()==1){
    METv = UpdateMETMuon(METv, muons_loose);
    METv = UpdateMETElectron(METv, electrons_loose);
  }*/

  double MET = METv_central.Pt();

  //========================================================    
  //==== Define particles, variables
  //========================================================

  double mu_tight_iso = 0.07;
  if(param.Muon_Tight_ID.Contains("ISR")) mu_tight_iso = 0.15;
  if(param.Muon_Tight_ID.Contains("HNTight")) mu_tight_iso = 0.1;
  double el_tight_iso = 0.;   

  // POG cut-based Medium  
  // barrel : 0.0478+0.506/pT, endcap : 0.0658+0.963/pT
  // POG cut-based Tight
  // barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT

  //double pi = 3.14159265358979323846;
  double MZ = 91.1876;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;
  double jet_emfraction = 0.;

  double trigLumi = 1.; 
  double jetPtCut = 40.;

  double ptcone_mu = 0.;
  double ptcone_el = 0.;
  //double trkiso_Pt = 0.;
  double trkiso_miniaodPt = 0.;
  //double ptcone_mu1 = 0.;
  TString PtConeRange = "";
  Particle ZCand, METv;

  /*Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2);*/

  //========================================================
  //==== Muon Fake Rate Measurement
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions_mu.size(); it_rg++){
    weight = 1.;
    //if(param.Electron_Loose_ID.Contains("V23")) break; // NO_MUON_EVENT : Save muon hists only once!! 

    // Fake rate measurement region
    if(it_rg == 0){

      if(!IsNorm) continue;

      if(!(muons_loose.size()==1 && electrons_loose.size()==0)) continue;
      if(!(jets.size() >= 1)) continue;

      // MET
      METv = UpdateMETMuon(METv_central, muons_loose);
      MET = METv.Pt();

      // Set up pTcone
      ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
      //trkiso_Pt = muons_loose.at(0).TrkIso()/muons_loose.at(0).Pt();
      trkiso_miniaodPt = muons_loose.at(0).TrkIso()/muons_loose.at(0).MiniAODPt();
      //ptcone_mu1 = muons_loose.at(0).Pt()*(1.+std::max(0., muons_loose.at(0).RelIso()-mu_tight_iso));
      //FillHist("PtCone_ratio", ptcone_mu1/ptcone_mu, weight, 20, 0., 2.);

      // Truth matching
      muons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons_loose, gens);
      if(!(muons_prompt.size() == 1)) continue;

      // Weights except trigger luminosity
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      // For checking TrkIso and RelIso
      if(ev.PassTrigger(MuonTrig1)){
        //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_Pt_PassMu3", trkiso_Pt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_PassMu3", trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_RelIso_PassMu3", muons_loose.at(0).RelIso(), weight, 20, 0., 1.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_Pt_PassMu8", trkiso_Pt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_PassMu8", trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_RelIso_PassMu8", muons_loose.at(0).RelIso(), weight, 20, 0., 1.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_Pt_PassMu17", trkiso_Pt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_PassMu17", trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_RelIso_PassMu17", muons_loose.at(0).RelIso(), weight, 20, 0., 1.);
      }

      trigLumi = 1.;
      // One prescaled trigger for each pTcone range, setup lumi
      if(!(ptcone_mu >= MuonPtconeCut1)) continue;
      if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2){
        if(!(muons_loose.at(0).Pt() > MuonPtCut1)) continue;
        if(!ev.PassTrigger(MuonTrig1)) continue;
        if(!IsDATA) trigLumi = MuonLumi1;
        jetPtCut = 50.;
        PtConeRange = "Range0";
      }
      if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut2)) continue;
        if(!ev.PassTrigger(MuonTrig2)) continue;
        if(!IsDATA) trigLumi = MuonLumi2; 
        PtConeRange = "Range1";
      }
      if(ptcone_mu >= MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut3)) continue;
        if(!ev.PassTrigger(MuonTrig3)) continue;
        if(!IsDATA) trigLumi = MuonLumi3; 
        PtConeRange = "Range2";
      }

      weight *= trigLumi;

      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_Eta_NoDijet_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      // Away jet selection
      jets_awayFromMuon.clear();
      jets_awayFromMuon = JetsAwayFromLepton(jets, muons_loose.at(0), dphiCut);
      std::sort(jets_awayFromMuon.begin(), jets_awayFromMuon.end(), PtComparing);

      /*for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        // define dphi between a jet and the loose lepton
        dphi = fabs(jets.at(ijet).Phi() - muons_loose.at(0).Phi());
        if(dphi > pi) dphi = 2.*pi-dphi;
        //dphi = fabs(muons_loose.at(0).DeltaPhi(jets.at(ijet)));
        FillHist(regions_mu.at(it_rg)+"/dphi_"+PtConeRange, dphi, weight, 32, 0., 3.2);

        if(dphi > 2.5) awayjet++; 
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }*/

      if(!(jets_awayFromMuon.size() > 0)) continue;
      if(!(jets_awayFromMuon.at(0).Pt() > jetPtCut)) continue;
 
      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets_awayFromMuon.at(0).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_Pt_NoCut_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_NoCut_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_NoCut_"+PtConeRange, 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_Eta_NoCut_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_TrkIso_Pt_NoCut_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt_NoCut_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_NoCut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;

      // Histograms after applying cuts
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_Pt_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
      FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_"+PtConeRange, 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_Eta_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        //FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_TrkIso_Pt_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        // Loose ID
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_IB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_IB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_IB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){
        // Loose ID
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_OB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_OB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_OB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

      // Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.5){
        // Loose ID
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_EC_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_LooseID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_EC_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Muon_TightID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_EC_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

    }

    // DY control region
    if(it_rg == 1){
      if(!(muons_tight.size()==2 && electrons_tight.size()==0)) continue;

      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      ZCand = muons_tight.at(0) + muons_tight.at(1);

      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep1_Pt_NoCut_Mu3", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep2_Pt_NoCut_Mu3", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_NoCut_Mu3", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep1_Pt_NoCut_Mu8", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep2_Pt_NoCut_Mu8", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_NoCut_Mu8", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep1_Pt_NoCut_Mu17", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep2_Pt_NoCut_Mu17", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_NoCut_Mu17", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }

      // Event selection
      if(!(muons_tight.at(0).Pt()>20. && muons_tight.at(1).Pt()>10.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Inclusive_Mu3", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu3", 0.5, weight*trigLumi, 2, 0., 2.);
        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Mu3", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu3", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Inclusive_Mu8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu8", 0.5, weight*trigLumi, 2, 0., 2.);
        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Mu8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu8", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Inclusive_Mu17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu17", 0.5, weight*trigLumi, 2, 0., 2.);
        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/ZCand_Mass_Mu17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu17", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }

    }

    // W+jets control region
    if(it_rg == 2){
      if(!(muons_tight.size()==1 && electrons_tight.size()==0)) continue;

      // MET
      METv = UpdateMETMuon(METv_central, muons_tight);
      MET = METv.Pt();

      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      Mt = MT(muons_tight.at(0), METv);

      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep_Pt_NoCut_Mu3", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/MET_NoCut_Mu3", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_NoCut_Mu3", Mt, weight*trigLumi, 500, 0., 500.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep_Pt_NoCut_Mu8", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/MET_NoCut_Mu8", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_NoCut_Mu8", Mt, weight*trigLumi, 500, 0., 500.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Lep_Pt_NoCut_Mu17", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/MET_NoCut_Mu17", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_NoCut_Mu17", Mt, weight*trigLumi, 500, 0., 500.);
      }

      // Event selection
      if(!(muons_tight.at(0).Pt() > 20.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;
      
      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_Mu3", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu3", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_Mu8", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Mt_Mu17", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_mu.at(it_rg)+"/Number_Events_Mu17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }

  }

  //========================================================
  //==== Electron Fake Rate Measurement
  //========================================================
 
  for(unsigned int it_rg2=0; it_rg2<regions_el.size(); it_rg2++){
    weight = 1.;
    //if(param.Muon_Tight_ID.Contains("V2") || param.Muon_Tight_ID.Contains("2016")) break;

    // Fake rate measurement region 
    if(it_rg2 == 0){

      if(!IsNorm) continue;

      if(!(muons_loose.size()==0 && electrons_loose.size()==1)) continue;
      if(!(jets.size() >= 1)) continue;

      // MET
      METv = UpdateMETElectron(METv_central, electrons_loose);
      MET = METv.Pt();

      // Set up pTcone
      el_tight_iso = 0.08;
      //el_tight_iso = 0.0287+0.506/electrons_loose.at(0).UncorrPt();
      //if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons_loose.at(0).UncorrPt();
      if(param.Electron_Tight_ID.Contains("ISR")){
        el_tight_iso = 0.0478+0.506/electrons_loose.at(0).UncorrPt();
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0658+0.963/electrons_loose.at(0).UncorrPt();
      }
      if(param.Electron_Tight_ID.Contains("HNTight")){
        el_tight_iso = std::min(0.08, 0.0287+0.506/electrons_loose.at(0).UncorrPt());
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/electrons_loose.at(0).UncorrPt());
      }

      // For the same bin in 2016 analysis
      if(param.Electron_Tight_ID.Contains("2016")) ElectronPtconeCut2 = 23.;
      else ElectronPtconeCut2 = 25.;

      ptcone_el = electrons_loose.at(0).CalcPtCone(electrons_loose.at(0).RelIso(), el_tight_iso);
      
      // Truth matching
      electrons_prompt.clear();
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons_loose, gens);
      if(!(electrons_prompt.size() == 1)) continue;

      // Weights except trigger luminosity
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      trigLumi = 1.;
      // One prescaled trigger for each PtCone range, setup lumi
      if(!(ptcone_el >= ElectronPtconeCut1)) continue;
      if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut1)) continue;
        if(!ev.PassTrigger(ElectronTrig1)) continue;
        if(!IsDATA) trigLumi = ElectronLumi1;
        PtConeRange = "Range0";
      }
      if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut2)) continue;      
        if(!ev.PassTrigger(ElectronTrig2)) continue;
        if(!IsDATA) trigLumi = ElectronLumi2;
        PtConeRange = "Range1";
      }
      if(ptcone_el >= ElectronPtconeCut3){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut3)) continue;
        if(!ev.PassTrigger(ElectronTrig3)) continue;
        if(!IsDATA) trigLumi = ElectronLumi3;
        PtConeRange = "Range2";
      }
      if(ptcone_el >= ElectronPtconeCut4){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut4)) continue;
        if(!ev.PassTrigger(ElectronTrig4)) continue;
        if(!IsDATA) trigLumi = ElectronLumi4;
        PtConeRange = "Range3";
      }

      weight *= trigLumi;

      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_Eta_NoDijet_"+PtConeRange, electrons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);        
      
      // Away jet selection
      jets_awayFromElectron.clear();
      jets_awayFromElectron = JetsAwayFromLepton(jets, electrons_loose.at(0), dphiCut);
      std::sort(jets_awayFromElectron.begin(), jets_awayFromElectron.end(), PtComparing);

      /*for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        //dphi = fabs(jets.at(ijet).Phi() - electrons_loose.at(0).Phi());
        //if(dphi > pi) dphi = 2.*pi-dphi;
        dphi = fabs(electrons_loose.at(0).DeltaPhi(jets.at(ijet)));
        FillHist("dphi_"+PtConeRange+"_"+regions_el.at(it_rg2), dphi, weight, 32, 0., 3.2);

        if(dphi > 2.5) awayjet++;
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }*/

      if(!(jets_awayFromElectron.size() > 0)) continue;
      if(!(jets_awayFromElectron.at(0).Pt() > jetPtCut)) continue;

      Mt = MT(electrons_loose.at(0), METv);
      Pt_ratio = jets_awayFromElectron.at(0).Pt()/electrons_loose.at(0).Pt();
      jet_emfraction = jets_awayFromElectron.at(0).ChargedEmEnergyFraction();

      // Histograms before applying cuts
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_Eta_NoCut_"+PtConeRange, electrons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Jet_ChargedEmEnergyFraction_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_NoCut_"+PtConeRange, 0.5, weight, 2, 0., 2.); 

      if(electrons_tight.size() > 0){
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_Eta_NoCut_"+PtConeRange, electrons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_NoCut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;
      if(!(jet_emfraction < 0.65)) continue;

      // Histograms after applying cuts
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_Eta_"+PtConeRange, electrons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_"+PtConeRange, 0.5, weight, 2, 0., 2.);
      if(electrons_tight.size() > 0){
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_Eta_"+PtConeRange, electrons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(electrons_loose.at(0).Eta()) < 0.8){
        // Loose ID
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_IB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_IB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_IB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(electrons_loose.at(0).Eta()) >= 0.8 && fabs(electrons_loose.at(0).Eta()) < 1.479){
        // Loose ID
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_OB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_OB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_OB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

      // Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(electrons_loose.at(0).Eta()) >= 1.479 && fabs(electrons_loose.at(0).Eta()) < 2.5){
        // Loose ID
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_EC_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){
          if(Nbjet_medium == 0){
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_LooseID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_EC_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){
            if(Nbjet_medium == 0){
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Electron_TightID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_EC_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
          }
        }

      }

    }

    // DY control region
    if(it_rg2 == 1){
      if(!(muons_tight.size()==0 && electrons_tight.size()==2)) continue;

      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      ZCand = electrons_tight.at(0) + electrons_tight.at(1);

      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep1_Pt_NoCut_Ele8", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep2_Pt_NoCut_Ele8", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_NoCut_Ele8", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep1_Pt_NoCut_Ele12", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep2_Pt_NoCut_Ele12", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_NoCut_Ele12", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep1_Pt_NoCut_Ele17", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep2_Pt_NoCut_Ele17", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_NoCut_Ele17", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep1_Pt_NoCut_Ele23", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep2_Pt_NoCut_Ele23", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_NoCut_Ele23", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }

      // Event selection
      if(!(electrons_tight.at(0).Pt()>25. && electrons_tight.at(1).Pt()>15.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Inclusive_Ele8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele8", 0.5, weight*trigLumi, 2, 0., 2.);
        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Ele8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele8", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Inclusive_Ele12", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele12", 0.5, weight*trigLumi, 2, 0., 2.);
        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Ele12", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele12", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Inclusive_Ele17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele17", 0.5, weight*trigLumi, 2, 0., 2.);
        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Ele17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele17", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Inclusive_Ele23", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele23", 0.5, weight*trigLumi, 2, 0., 2.);
        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/ZCand_Mass_Ele23", ZCand.M(), weight*trigLumi, 40, 70., 110.);
          FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele23", 1.5, weight*trigLumi, 2, 0., 2.);
        }
      }

    }

    // W+jets control region
    if(it_rg2 == 2){
      if(!(muons_tight.size()==0 && electrons_tight.size()==1)) continue;

      // MET
      METv = UpdateMETElectron(METv_central, electrons_tight);
      MET = METv.Pt();

      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      Mt = MT(electrons_tight.at(0), METv);

      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep_Pt_NoCut_Ele8", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/MET_NoCut_Ele8", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_NoCut_Ele8", Mt, weight*trigLumi, 500, 0., 500.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep_Pt_NoCut_Ele12", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/MET_NoCut_Ele12", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_NoCut_Ele12", Mt, weight*trigLumi, 500, 0., 500.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep_Pt_NoCut_Ele17", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/MET_NoCut_Ele17", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_NoCut_Ele17", Mt, weight*trigLumi, 500, 0., 500.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Lep_Pt_NoCut_Ele23", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/MET_NoCut_Ele23", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_NoCut_Ele23", Mt, weight*trigLumi, 500, 0., 500.);
      }

      // Event selection
      if(!(electrons_tight.at(0).Pt() > 25.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_Ele8", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_Ele12", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele12", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_Ele17", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Mt_Ele23", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(systName+"/"+regions_el.at(it_rg2)+"/Number_Events_Ele23", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }

  }

}


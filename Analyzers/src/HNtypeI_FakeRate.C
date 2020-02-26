#include "HNtypeI_FakeRate.h"

HNtypeI_FakeRate::HNtypeI_FakeRate(){

}

void HNtypeI_FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  MuonTightIDs = {"HNTight", "HNTightV2", "HNTight2016"};
  MuonLooseIDs = {"HNLoose", "HNLoose", "HNLoose2016"};
  ElectronTightIDs = {"ISRTight", "HNTight2016", "HNTight", "HNTightV2", "HNTightV2", "HNTightV2", "HNMVATight", "HNMVATightV2"};
  ElectronLooseIDs = {"ISRLoose", "HNLoose2016", "HNLoose", "HNLooseV21", "HNLooseV22", "HNLooseV23", "HNMVALoose", "HNMVALooseV2"};
/*  MuonTightIDs = {"HNTightV2"};
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
  }
//  cout << "[HNtypeI_FakeRate::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[HNtypeI_FakeRate::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

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
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_EleID=0; it_EleID<ElectronTightIDs.size(); it_EleID++){
    for(unsigned int it_MuonID=0; it_MuonID<MuonTightIDs.size(); it_MuonID++){    
      if(it_EleID<3 && (it_EleID != it_MuonID)) continue;
      if(it_EleID >= 3){
        if(it_MuonID != 0) continue;
      }

//      TString MuonID = "HNTight2016";
//      TString MuonIDSFKey = "NUM_TightID_DEN_genTracks";
      TString MuonTightID = MuonTightIDs.at(it_MuonID);
      TString MuonLooseID = MuonLooseIDs.at(it_MuonID);
      TString ElectronTightID = ElectronTightIDs.at(it_EleID);
      TString ElectronLooseID = ElectronLooseIDs.at(it_EleID);

      param.Clear();

      param.syst_ = AnalyzerParameter::Central;

//      param.Name = MuonID+"_"+"Central";

      param.Electron_Tight_ID = ElectronTightID;
      param.Electron_Loose_ID = ElectronLooseID;
//      param.Electron_Tight_ID = "HNTight";
//      param.Electron_Loose_ID = "HNLoose";
      param.Electron_Veto_ID = "";
      param.Electron_ID_SF_Key = "";
      param.Muon_Tight_ID = MuonTightID;
      param.Muon_Loose_ID = MuonLooseID;
      param.Muon_Veto_ID = "";
      param.Muon_ID_SF_Key = "";
      param.Muon_ISO_SF_Key = "";
      param.Jet_ID = "tight";

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

void HNtypeI_FakeRate::executeEventFromParameter(AnalyzerParameter param){
 
  // Version : Loose ID

  vector<TString> regions_mu = {"muonV1FR", "muonV1DY", "muonV1Wj"};
  if(param.Muon_Tight_ID.Contains("V2")) regions_mu = {"muonV2FR", "muonV2DY", "muonV2Wj"};
  if(param.Muon_Tight_ID.Contains("2016")) regions_mu = {"muon16FR", "muon16DY", "muon16Wj"};

  vector<TString> regions_el = {"eleV1FR", "eleV1DY", "eleV1Wj"};
  if(param.Electron_Loose_ID.Contains("V21")) regions_el = {"eleV21FR", "eleV21DY", "eleV21Wj"};
  if(param.Electron_Loose_ID.Contains("V22")) regions_el = {"eleV22FR", "eleV22DY", "eleV22Wj"};
  if(param.Electron_Loose_ID.Contains("V23")) regions_el = {"eleV23FR", "eleV23DY", "eleV23Wj"};
  if(param.Electron_Tight_ID.Contains("MVA")){ 
    regions_el = {"mvaV1FR", "mvaV1DY", "mvaV1Wj"};
    if(param.Electron_Tight_ID.Contains("V2")) regions_el = {"mvaV2FR", "mvaV2DY", "mvaV2Wj"};
  }
  if(param.Electron_Tight_ID.Contains("2016")) regions_el = {"ele16FR", "ele16DY", "ele16Wj"};
  if(param.Electron_Tight_ID.Contains("ISR")) regions_el = {"isrFR", "isrDY", "isrWj"};

  //=============
  //==== No Cut
  //=============

//  JSFillHist(param.Name, "NoCut_"+param.Name, 0., 1., 1, 0., 1.);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  //==============
  //==== Trigger
  //==============
//  if(! (ev.PassTrigger(MuonTriggers) )) return;



  //======================
  //==== Copy AllObjects
  //======================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;

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
    cout << "[HNtypeI_FakeRate::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Electron> ele_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  vector<Electron> ele_loose = SelectElectrons(this_AllElectrons, param.Electron_Loose_ID, ElectronPtCut1, 2.5);
//  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);
  vector<Muon> muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  vector<Muon> muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, MuonPtCut1, 2.4);
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);
  vector<Gen> gens = GetGens(); 
//  std::vector<Lepton*> leptons;

  //=======================
  //==== Sort in pt-order
  //=======================

  std::sort(ele_loose.begin(), ele_loose.end(), PtComparing);
  std::sort(ele_tight.begin(), ele_tight.end(), PtComparing);
  std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
  std::sort(muons_tight.begin(), muons_tight.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);

  //==== B-Tagging
  int Nbjet_medium=0;
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV,
                                                                     JetTagging::Medium,
                                                                     JetTagging::incl, JetTagging::comb); 

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
//  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets.at(ij))) Nbjet_medium++;
  }
    
  //=========================
  //==== Event selections..
  //=========================

  double mu_tight_iso = 0.15;
  if(param.Muon_Tight_ID.Contains("V2")) mu_tight_iso = 0.1;
  if(param.Muon_Tight_ID.Contains("2016")) mu_tight_iso = 0.07;
  double el_tight_iso = 0.;   
  // POG cut-based Medium  
  // barrel : 0.0478+0.506/pT, endcap : 0.0658+0.963/pT
  // POG cut-based Tight
  // barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT
//  double pi = 3.14159265358979323846;
  double MZ = 91.1876;
  double MET = ev.GetMETVector().Pt(); 
  double dphi = 0.;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;
  double jet_emfraction = 0.;

  double trigLumi = 1.; 
  double awayjet_ptcut = 40.;
  int awayjet = 0;
  int leadingjet = 0;

  double ptcone_mu = 0.;
  double ptcone_el = 0.;
  double trkiso_Pt = 0.;
  double trkiso_miniaodPt = 0.;
//  double ptcone_mu1 = 0.;
  TString PtConeRange = "";
  Particle ZCand;

/*  Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2); */

  //===================================================================================================
  //==== Muon Fake Rate Measurement
  //===================================================================================================

  for(unsigned int it_rg=0; it_rg<regions_mu.size(); it_rg++){
    weight = 1.;
    if(param.Electron_Tight_ID.Contains("V2") || param.Electron_Tight_ID.Contains("MVA")) break; // Stack muon hists only once!!

    // Fake rate measurement region 
    if(it_rg==0 && muons_loose.size()==1){
      if(jets.size() == 0) continue;
      ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
      trkiso_Pt = muons_loose.at(0).TrkIso()/muons_loose.at(0).Pt();
      trkiso_miniaodPt = muons_loose.at(0).TrkIso()/muons_loose.at(0).MiniAODPt();
//      ptcone_mu1 = muons_loose.at(0).Pt()*(1.+std::max(0., muons_loose.at(0).RelIso()-mu_tight_iso));
//      FillHist("PtCone_ratio", ptcone_mu1/ptcone_mu, weight, 20, 0., 2.);

      if(!IsDATA){
        // Gen matching
        // No matched gen lepton -> PID == 0 (See DataFormats/src/Gen.C)
        Gen truth_mu = GetGenMatchedLepton(muons_loose.at(0), gens);
        if(truth_mu.PID() == 0) continue;

        // Weights except trigger luminosity
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      // For checking TrkIso and RelIso
      if(ev.PassTrigger(MuonTrig1)){
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_Pt_PassMu3", trkiso_Pt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_MiniAODPt_PassMu3", trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_RelIso_PassMu3", muons_loose.at(0).RelIso(), weight, 20, 0., 1.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_Pt_PassMu8", trkiso_Pt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_MiniAODPt_PassMu8", trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_RelIso_PassMu8", muons_loose.at(0).RelIso(), weight, 20, 0., 1.);
      }

      trigLumi = 1.;
      // only 1 prescaled trigger for each PtCone range, setup lumi
      if(ptcone_mu < MuonPtconeCut1) continue;
      if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2){
        if(muons_loose.at(0).Pt() < MuonPtCut1) continue;
        if(!ev.PassTrigger(MuonTrig1)) continue;
        if(!IsDATA) trigLumi = MuonLumi1;
        awayjet_ptcut = 50.;
        PtConeRange = "range0";
      }
      if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3){
        if(muons_loose.at(0).Pt() < MuonPtCut2) continue;
        if(!ev.PassTrigger(MuonTrig2)) continue;
        if(!IsDATA) trigLumi = MuonLumi2; 
        PtConeRange = "range1";
      }
      if(ptcone_mu >= MuonPtconeCut3){
        if(muons_loose.at(0).Pt() < MuonPtCut3) continue;
        if(!ev.PassTrigger(MuonTrig3)) continue;
        if(!IsDATA) trigLumi = MuonLumi3; 
        PtConeRange = "range2";
      }

      awayjet = 0, leadingjet = 0;
      weight *= trigLumi;

      FillHist(regions_mu.at(it_rg)+"/Muon_loose_Eta_nodijet_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        // define dphi between a jet and the loose lepton
//        dphi = fabs(jets.at(ijet).Phi() - muons_loose.at(0).Phi());
//        if(dphi > pi) dphi = 2.*pi-dphi;
        dphi = fabs(muons_loose.at(0).DeltaPhi(jets.at(ijet)));
        FillHist(regions_mu.at(it_rg)+"/dphi_"+PtConeRange, dphi, weight, 32, 0., 3.2);

        if(dphi > 2.5) awayjet++; 
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }

      // away jet selection
      if(awayjet == 0) continue;
      if(jets.at(leadingjet).Pt() < awayjet_ptcut) continue;
 
      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets.at(leadingjet).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist(regions_mu.at(it_rg)+"/MET_nocut_"+PtConeRange, MET, weight, 200, 0., 200.);
      FillHist(regions_mu.at(it_rg)+"/Mt_nocut_"+PtConeRange, Mt, weight, 200, 0., 200.);
      FillHist(regions_mu.at(it_rg)+"/Ptratio_nocut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_PtCone_nocut_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_Eta_nocut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_Pt_nocut_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_MiniAODPt_nocut_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
      FillHist(regions_mu.at(it_rg)+"/Number_Events_nocut_"+PtConeRange, 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_PtCone_nocut_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_Eta_nocut_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_TrkIso_Pt_nocut_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_TrkIso_MiniAODPt_nocut_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_nocut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // additional cuts to reduce prompt contribution
      if(MET > 80.) continue;
      if(Mt > 25.) continue;
      if(Pt_ratio < 1.) continue;

      // Histograms after applying cuts
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_PtCone_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_Pt_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
      FillHist(regions_mu.at(it_rg)+"/Muon_loose_TrkIso_MiniAODPt_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
      FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_"+PtConeRange, 0.5, weight, 2, 0., 2.);
      if(ptcone_mu > 10.) FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone10_"+PtConeRange, 0.5, weight, 2, 0., 2.);
      if(muons_tight.size() > 0){
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_PtCone_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_Eta_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_TrkIso_Pt_"+PtConeRange, trkiso_Pt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Muon_tight_TrkIso_MiniAODPt_"+PtConeRange, trkiso_miniaodPt, weight, 20, 0., 1.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_PtCone_barrel1_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_barrel1_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone10_barrel1_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist(regions_mu.at(it_rg)+"/Muon_tight_PtCone_barrel1_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
          FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_barrel1_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }

      // outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_PtCone_barrel2_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_barrel2_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone10_barrel2_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist(regions_mu.at(it_rg)+"/Muon_tight_PtCone_barrel2_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
          FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_barrel2_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }

      // endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.5){
        FillHist(regions_mu.at(it_rg)+"/Muon_loose_PtCone_endcap_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_endcap_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone10_endcap_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist(regions_mu.at(it_rg)+"/Muon_tight_PtCone_endcap_"+PtConeRange, ptcone_mu, weight, 200, 0., 200.);
          FillHist(regions_mu.at(it_rg)+"/Number_Events_PtCone5_endcap_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }
    }

    // DY control region
    if(it_rg==1 && muons_tight.size()==2){
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
        FillHist(regions_mu.at(it_rg)+"/Lep1_Pt_nocut_Mu3", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Lep2_Pt_nocut_Mu3", muons_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_nocut_Mu3", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(regions_mu.at(it_rg)+"/Lep1_Pt_nocut_Mu8", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Lep2_Pt_nocut_Mu8", muons_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_nocut_Mu8", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(regions_mu.at(it_rg)+"/Lep1_Pt_nocut_Mu17", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Lep2_Pt_nocut_Mu17", muons_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_nocut_Mu17", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }

      // event selection
      if(muons_tight.at(0).Pt() < 20. || muons_tight.at(1).Pt() < 10.) continue;
      if(fabs(ZCand.M() - MZ) > 10.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_Mu3", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu3", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_Mu8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(regions_mu.at(it_rg)+"/ZCand_Mass_Mu17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }

    // W+jets control region
    if(it_rg==2 && muons_tight.size()==1){
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
        FillHist(regions_mu.at(it_rg)+"/Lep_Pt_nocut_Mu3", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/MET_nocut_Mu3", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Mt_nocut_Mu3", Mt, weight*trigLumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(regions_mu.at(it_rg)+"/Lep_Pt_nocut_Mu8", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/MET_nocut_Mu8", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Mt_nocut_Mu8", Mt, weight*trigLumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(regions_mu.at(it_rg)+"/Lep_Pt_nocut_Mu17", muons_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/MET_nocut_Mu17", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Mt_nocut_Mu17", Mt, weight*trigLumi, 200, 0., 200.);
      }

      // event selection
      if(muons_tight.at(0).Pt() < 20.) continue;
      if(MET < 40.) continue;
      if(Mt < 60. || Mt > 100.) continue;
      
      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1;
        FillHist(regions_mu.at(it_rg)+"/Mt_Mu3", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu3", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2;
        FillHist(regions_mu.at(it_rg)+"/Mt_Mu8", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3;
        FillHist(regions_mu.at(it_rg)+"/Mt_Mu17", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_mu.at(it_rg)+"/Number_Events_Mu17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }
  }

  //==================================================================================================
  //==== Electron Fake Rate Measurement
  //==================================================================================================
 
  for(unsigned int it_rg2=0; it_rg2<regions_el.size(); it_rg2++){
    weight = 1.;
//    if(param.Muon_Tight_ID.Contains("V2") || param.Muon_Tight_ID.Contains("2016")) break;

    // Fake rate measurement region 
    if(it_rg2==0 && ele_loose.size()==1){
      if(jets.size() == 0) continue;
      el_tight_iso = 0.0287+0.506/ele_loose.at(0).UncorrPt();
      if(fabs(ele_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/ele_loose.at(0).UncorrPt();
      if(param.Electron_Tight_ID.Contains("V2")){
        el_tight_iso = std::min(0.08, 0.0287+0.506/ele_loose.at(0).UncorrPt());
        if(fabs(ele_loose.at(0).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/ele_loose.at(0).UncorrPt());
      }
      if(param.Electron_Tight_ID.Contains("2016")) el_tight_iso = 0.08;
      if(param.Electron_Tight_ID.Contains("ISR")){
        el_tight_iso = 0.0478+0.506/ele_loose.at(0).UncorrPt();
        if(fabs(ele_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0658+0.963/ele_loose.at(0).UncorrPt();
      }

      // For the same bin in 2016 analysis
      if(param.Electron_Tight_ID.Contains("2016")) ElectronPtconeCut2 = 23.;
      else ElectronPtconeCut2 = 25.;

      ptcone_el = ele_loose.at(0).CalcPtCone(ele_loose.at(0).RelIso(), el_tight_iso);
      
      if(!IsDATA){
        // Gen matching
        // No matched gen lepton -> PID == 0 (See DataFormats/src/Gen.C)
        Gen truth_el = GetGenMatchedLepton(ele_loose.at(0), gens);
        if(truth_el.PID() == 0) continue; 
 
        // weight except trigger luminosity
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      trigLumi = 1.;
      // only 1 prescaled trigger for each PtCone range, setup lumi
      if(ptcone_el < ElectronPtconeCut1) continue;
      if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2){
        if(ele_loose.at(0).Pt() < ElectronPtCut1) continue;
        if(!ev.PassTrigger(ElectronTrig1)) continue;
        if(!IsDATA) trigLumi = ElectronLumi1;
        PtConeRange = "range0";
      }
      if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3){
        if(ele_loose.at(0).Pt() < ElectronPtCut2) continue;      
        if(!ev.PassTrigger(ElectronTrig2)) continue;
        if(!IsDATA) trigLumi = ElectronLumi2;
        PtConeRange = "range1";
      }
      if(ptcone_el >= ElectronPtconeCut3){
        if(ele_loose.at(0).Pt() < ElectronPtCut3) continue;
        if(!ev.PassTrigger(ElectronTrig3)) continue;
        if(!IsDATA) trigLumi = ElectronLumi3;
        PtConeRange = "range2";
      }
      if(ptcone_el >= ElectronPtconeCut4){
        if(ele_loose.at(0).Pt() < ElectronPtCut4) continue;
        if(!ev.PassTrigger(ElectronTrig4)) continue;
        if(!IsDATA) trigLumi = ElectronLumi4;
        PtConeRange = "range3";
      }

      awayjet = 0, leadingjet = 0;
      weight *= trigLumi;

      FillHist(regions_el.at(it_rg2)+"/Electron_loose_Eta_nodijet_"+PtConeRange, ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);        
      
      for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        dphi = fabs(jets.at(ijet).Phi() - ele_loose.at(0).Phi());
        if(dphi > pi) dphi = 2.*pi-dphi;
//        dphi = fabs(ele_loose.at(0).DeltaPhi(jets.at(ijet)));
        FillHist("dphi_"+PtConeRange+"_"+regions_el.at(it_rg2), dphi, weight, 32, 0., 3.2);

        if(dphi > 2.5) awayjet++;
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }

      // away jet selection      
      if(awayjet == 0) continue;
      if(jets.at(leadingjet).Pt() < awayjet_ptcut) continue;

      Mt = MT(ele_loose.at(0), METv);
      Pt_ratio = jets.at(leadingjet).Pt()/ele_loose.at(0).Pt();
      jet_emfraction = jets.at(leadingjet).ChargedEmEnergyFraction();

      // Histograms before applying cuts
      FillHist(regions_el.at(it_rg2)+"/MET_nocut_"+PtConeRange, MET, weight, 200, 0., 200.);
      FillHist(regions_el.at(it_rg2)+"/Mt_nocut_"+PtConeRange, Mt, weight, 200, 0., 200.);
      FillHist(regions_el.at(it_rg2)+"/Ptratio_nocut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(regions_el.at(it_rg2)+"/Electron_loose_PtCone_nocut_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
      FillHist(regions_el.at(it_rg2)+"/Electron_loose_Eta_nocut_"+PtConeRange, ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions_el.at(it_rg2)+"/Jet_ChargedEmEnergyFraction_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);
      FillHist(regions_el.at(it_rg2)+"/Number_Events_nocut_"+PtConeRange, 0.5, weight, 2, 0., 2.); 

      if(ele_tight.size() > 0){
        FillHist(regions_el.at(it_rg2)+"/Electron_tight_PtCone_nocut_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Electron_tight_Eta_nocut_"+PtConeRange, ele_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_nocut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // additional cuts to reduce prompt contribution
      if(MET > 80.) continue;
      if(Mt > 25.) continue;
      if(Pt_ratio < 1.) continue;
      if(jet_emfraction > 0.65) continue;

      // Histograms after applying cuts
      FillHist(regions_el.at(it_rg2)+"/Electron_loose_PtCone_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
      FillHist(regions_el.at(it_rg2)+"/Electron_loose_Eta_"+PtConeRange, ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(regions_el.at(it_rg2)+"/Number_Events_"+PtConeRange, 0.5, weight, 2, 0., 2.);
      if(ele_tight.size() > 0){
        FillHist(regions_el.at(it_rg2)+"/Electron_tight_PtCone_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Electron_tight_Eta_"+PtConeRange, ele_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // inner barrel ( |eta| < 0.8 )
      if(fabs(ele_loose.at(0).Eta()) < 0.8){
        FillHist(regions_el.at(it_rg2)+"/Electron_loose_PtCone_barrel1_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_barrel1_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist(regions_el.at(it_rg2)+"/Electron_tight_PtCone_barrel1_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
          FillHist(regions_el.at(it_rg2)+"/Number_Events_barrel1_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }

      // outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(ele_loose.at(0).Eta()) >= 0.8 && fabs(ele_loose.at(0).Eta()) < 1.479){
        FillHist(regions_el.at(it_rg2)+"/Electron_loose_PtCone_barrel2_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_barrel2_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist(regions_el.at(it_rg2)+"/Electron_tight_PtCone_barrel2_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
          FillHist(regions_el.at(it_rg2)+"/Number_Events_barrel2_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }

      // endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(ele_loose.at(0).Eta()) >= 1.479 && fabs(ele_loose.at(0).Eta()) < 2.5){
        FillHist(regions_el.at(it_rg2)+"/Electron_loose_PtCone_endcap_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_endcap_"+PtConeRange, 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist(regions_el.at(it_rg2)+"/Electron_tight_PtCone_endcap_"+PtConeRange, ptcone_el, weight, 200, 0., 200.);
          FillHist(regions_el.at(it_rg2)+"/Number_Events_endcap_"+PtConeRange, 1.5, weight, 2, 0., 2.);
        }
      }
    }

    // DY control region
    if(it_rg2==1 && ele_tight.size()==2){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      ZCand = ele_tight.at(0) + ele_tight.at(1);
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(regions_el.at(it_rg2)+"/Lep1_Pt_nocut_Ele8", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Lep2_Pt_nocut_Ele8", ele_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_nocut_Ele8", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(regions_el.at(it_rg2)+"/Lep1_Pt_nocut_Ele12", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Lep2_Pt_nocut_Ele12", ele_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_nocut_Ele12", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(regions_el.at(it_rg2)+"/Lep1_Pt_nocut_Ele17", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Lep2_Pt_nocut_Ele17", ele_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_nocut_Ele17", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(regions_el.at(it_rg2)+"/Lep1_Pt_nocut_Ele23", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Lep2_Pt_nocut_Ele23", ele_tight.at(1).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_nocut_Ele23", ZCand.M(), weight*trigLumi, 80, 50., 130.);
      }

      // event selection
      if(ele_tight.at(0).Pt() < 25. || ele_tight.at(1).Pt() < 15.) continue;
      if(fabs(ZCand.M() - MZ) > 10.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_Ele8", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_Ele12", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele12", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_Ele17", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(regions_el.at(it_rg2)+"/ZCand_Mass_Ele23", ZCand.M(), weight*trigLumi, 40, 70., 110.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele23", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }

    // W+jets control region
    if(it_rg2==2 && ele_tight.size()==1){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      Mt = MT(ele_tight.at(0), METv);
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(regions_el.at(it_rg2)+"/Lep_Pt_nocut_Ele8", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/MET_nocut_Ele8", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Mt_nocut_Ele8", Mt, weight*trigLumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(regions_el.at(it_rg2)+"/Lep_Pt_nocut_Ele12", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/MET_nocut_Ele12", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Mt_nocut_Ele12", Mt, weight*trigLumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(regions_el.at(it_rg2)+"/Lep_Pt_nocut_Ele17", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/MET_nocut_Ele17", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Mt_nocut_Ele17", Mt, weight*trigLumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(regions_el.at(it_rg2)+"/Lep_Pt_nocut_Ele23", ele_tight.at(0).Pt(), weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/MET_nocut_Ele23", MET, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Mt_nocut_Ele23", Mt, weight*trigLumi, 200, 0., 200.);
      }

      // event selection
      if(ele_tight.at(0).Pt() < 25.) continue;
      if(MET < 40.) continue;
      if(Mt < 60. || Mt > 100.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1;
        FillHist(regions_el.at(it_rg2)+"/Mt_Ele8", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele8", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2;
        FillHist(regions_el.at(it_rg2)+"/Mt_Ele12", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele12", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3;
        FillHist(regions_el.at(it_rg2)+"/Mt_Ele17", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele17", 0.5, weight*trigLumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4;
        FillHist(regions_el.at(it_rg2)+"/Mt_Ele23", Mt, weight*trigLumi, 200, 0., 200.);
        FillHist(regions_el.at(it_rg2)+"/Number_Events_Ele23", 0.5, weight*trigLumi, 2, 0., 2.);
      }
    }
  }
}


#include "FakeRate.h"

FakeRate::FakeRate(){

}

void FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;

//  MuonIDs = { "POGTight" };
//  MuonIDSFKeys = { "NUM_TightID_DEN_genTracks" };
  ElectronTightIDs = {"HNTight", "HNMVATight", "ISRTight"};
  ElectronLooseIDs = {"HNLoose", "HNMVALoose", "ISRLoose"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/FakeRate.h 
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
  ElectronPtconeCut1 = 10., ElectronPtconeCut2 = 23., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 40.;

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
//  cout << "[FakeRate::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[FakeRate::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== Test btagging code
  //==== add taggers and WP that you want to use in analysis
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);
//  vtaggers.push_back(Jet::CSVv2);

  std::vector<Jet::WP> v_wps;
//  v_wps.push_back(Jet::Medium);
  v_wps.push_back(Jet::Loose);

  //=== list of taggers, WP, setup systematics, use period SFs
//  SetupBTagger(vtaggers,v_wps, true, true);
  SetupBTagger(vtaggers,v_wps, true, true);  

}

FakeRate::~FakeRate(){

  //==== Destructor of this Analyzer

}

void FakeRate::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/FakeRate.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/FakeRate.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_EleID=0; it_EleID<ElectronTightIDs.size(); it_EleID++){

    TString MuonID = "HNTight";
//    TString MuonID = "HNTight2016";
//    TString MuonIDSFKey = "NUM_TightID_DEN_genTracks";
    TString ElectronTightID = ElectronTightIDs.at(it_EleID);
    TString ElectronLooseID = ElectronLooseIDs.at(it_EleID);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;

    param.Name = MuonID+"_"+"Central";

    param.Electron_Tight_ID = ElectronTightID;
    param.Electron_Loose_ID = ElectronLooseID;
//    param.Electron_Tight_ID = "HNTight";
//    param.Electron_Loose_ID = "HNLoose";
    param.Electron_Veto_ID = "";
    param.Electron_ID_SF_Key = "";
    param.Muon_Tight_ID = MuonID;
    param.Muon_Loose_ID = "HNLoose";
//    param.Muon_Loose_ID = "HNLoose2016";
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

void FakeRate::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> regions_mu = {"muonFR", "muonDY", "muonWj"};
  vector<TString> regions_el = {"eleFR", "eleDY", "eleWj"};
  if(param.Electron_Tight_ID.Contains("MVA")) regions_el = {"mvaFR", "mvaDY", "mvaWj"};
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

  if(param.syst_ == AnalyzerParameter::Central){

  }
/*  else if(param.syst_ == AnalyzerParameter::JetResUp){
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
    cout << "[FakeRate::executeEventFromParameter] Wrong syst" << endl;
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
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
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

//  int n_bjet_deepcsv_m=0;
//  int n_bjet_deepcsv_m_noSF=0;
  int n_bjet_deepcsv=0;

  for(unsigned int ij = 0 ; ij < jets.size(); ij++){
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) n_bjet_deepcsv_m++; // method for getting btag with SF applied to MC
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,false,0)) n_bjet_deepcsv_m_noSF++; // method for getting btag with no SF applied to MC
    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Loose, true, 0)) n_bjet_deepcsv++;
  }
  
  //=========================
  //==== Event selections..
  //=========================

  double mu_tight_iso = 0.15;
  double el_tight_iso = 0.;     // barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT
  double pi = 3.14159265358979323846;
  double MZ = 91.1876;
  double MET = ev.GetMETVector().Pt(); 
  double dphi = 0.;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;
  double jet_emfraction = 0.;

  double triglumi = 1.; 
  double awayjet_ptcut = 40.;
  int awayjet = 0;
  int leadingjet = 0;

  double ptcone_mu = 0.;
  double ptcone_el = 0.;
//  double ptcone_mu1 = 0.;
  TString PtConeRange = "";
  Particle ZCand;

/*  Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2); */

/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Muon ///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

  for(unsigned int it_rg=0; it_rg<regions_mu.size(); it_rg++){
    weight = 1.;
    if(param.Electron_Tight_ID.Contains("MVA") || param.Electron_Tight_ID.Contains("ISR")) break; // Stack muon hists only once!!

    // Fake rate measurement region 
    if(muons_loose.size()==1 && it_rg==0){
      if(jets.size() == 0) continue;
      ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
//      ptcone_mu1 = muons_loose.at(0).Pt()*(1.+std::max(0., muons_loose.at(0).RelIso()-mu_tight_iso));
//      FillHist("PtCone_ratio", ptcone_mu1/ptcone_mu, weight, 20, 0., 2.);

      if(!IsDATA){
        // Gen matching
        // No matched gen lepton -> PID == 0 (See DataFormats/src/Gen.C)
        Gen truth_mu = GetGenMatchedLepton(muons_loose.at(0), gens);
        if(truth_mu.PID() == 0) continue;

        // weights except trigger luminosity
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

/*      if(ev.PassTrigger(MuonTrig1) && muons_loose.at(0).Pt() > MuonPtCut1){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi1;
        FillHist("Mu_loose_Eta_Mu3_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2) FillHist("Mu_loose_Eta_Mu3_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(!ev.PassTrigger(MuonTrig2) && !ev.PassTrigger(MuonTrig3)){
          FillHist("Mu_loose_Eta_Mu3only_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
          if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2) FillHist("Mu_loose_Eta_Mu3only_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        }
      }
      if(ev.PassTrigger(MuonTrig2) && muons_loose.at(0).Pt() > MuonPtCut2){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi2;
        FillHist("Mu_loose_Eta_Mu8_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3) FillHist("Mu_loose_Eta_Mu8_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(!ev.PassTrigger(MuonTrig1) && !ev.PassTrigger(MuonTrig3)){
          FillHist("Mu_loose_Eta_Mu8only_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
          if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3) FillHist("Mu_loose_Eta_Mu8only_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        }
      }
      if(ev.PassTrigger(MuonTrig3) && muons_loose.at(0).Pt() > MuonPtCut3){
         triglumi = 1.;
         if(!IsDATA) triglumi = MuonLumi3;
         FillHist("Mu_loose_Eta_Mu17_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
         if(ptcone_mu >= MuonPtconeCut3) FillHist("Mu_loose_Eta_Mu17_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
         if(!ev.PassTrigger(MuonTrig1) && !ev.PassTrigger(MuonTrig2)){
           FillHist("Mu_loose_Eta_Mu17only_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
           if(ptcone_mu >= MuonPtconeCut3) FillHist("Mu_loose_Eta_Mu17only_PtConeCut_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
         }
      }*/

      triglumi = 1.;
      // only 1 prescaled trigger for each PtCone range, setup lumi
      if(ptcone_mu < MuonPtconeCut1) continue;
      if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2){
        if(muons_loose.at(0).Pt() < MuonPtCut1) continue;
//        if(!(ev.PassTrigger(MuonTrig1) && !ev.PassTrigger(MuonTrig2) && !ev.PassTrigger(MuonTrig3))) continue;
        if(!ev.PassTrigger(MuonTrig1)) continue;
        if(!IsDATA) triglumi = MuonLumi1;
        awayjet_ptcut = 50.;
        PtConeRange = "range0";
      }
      if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3){
        if(muons_loose.at(0).Pt() < MuonPtCut2) continue;
//        if(!(!ev.PassTrigger(MuonTrig1) && ev.PassTrigger(MuonTrig2) && !ev.PassTrigger(MuonTrig3))) continue;
        if(!ev.PassTrigger(MuonTrig2)) continue;
        if(!IsDATA) triglumi = MuonLumi2; 
        PtConeRange = "range1";
      }
      if(ptcone_mu >= MuonPtconeCut3){
        if(muons_loose.at(0).Pt() < MuonPtCut3) continue;
//        if(!(!ev.PassTrigger(MuonTrig1) && !ev.PassTrigger(MuonTrig2) && ev.PassTrigger(MuonTrig3))) continue;
        if(!ev.PassTrigger(MuonTrig3)) continue;
        if(!IsDATA) triglumi = MuonLumi3; 
        PtConeRange = "range2";
      }

      awayjet = 0, leadingjet = 0;
      weight *= triglumi;

      FillHist("Mu_loose_Eta_nodijet_"+PtConeRange+"_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        // define dphi between a jet and the loose lepton
        dphi = fabs(jets.at(ijet).Phi() - muons_loose.at(0).Phi());
        if(dphi > pi) dphi = 2.*pi-dphi;
        FillHist("dphi_"+PtConeRange+"_"+regions_mu.at(it_rg), dphi, weight, 32, 0., 3.2);

        if(dphi > 2.5) awayjet++; 
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }

      // away jet selection
      if(awayjet == 0) continue;
      if(jets.at(leadingjet).Pt() < awayjet_ptcut) continue;
 
      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets.at(leadingjet).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist("MET_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), MET, weight, 200, 0., 200.);
      FillHist("Mt_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), Mt, weight, 200, 0., 200.);
      FillHist("Ptratio_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), Pt_ratio, weight, 50, 0., 5.);
      FillHist("Mu_loose_PtCone_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
      FillHist("Mu_loose_Eta_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Number_Events_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist("Mu_tight_PtCone_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Mu_tight_Eta_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist("Number_Events_nocut_"+PtConeRange+"_"+regions_mu.at(it_rg), 1.5, weight, 2, 0., 2.);
      }

      // additional cuts to reduce prompt contribution
      if(MET > 80.) continue;
      if(Mt > 25.) continue;
      if(Pt_ratio < 1.) continue;

      // Histograms after applying cuts
      FillHist("Mu_loose_PtCone_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
      FillHist("Mu_loose_Eta_"+PtConeRange+"_"+regions_mu.at(it_rg), muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Number_Events_PtCone5_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
      if(ptcone_mu > 10.) FillHist("Number_Events_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
      if(muons_tight.size() > 0){
        FillHist("Mu_tight_PtCone_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Mu_tight_Eta_"+PtConeRange+"_"+regions_mu.at(it_rg), muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist("Number_Events_"+PtConeRange+"_"+regions_mu.at(it_rg), 1.5, weight, 2, 0., 2.);
      }

      // inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        FillHist("Mu_loose_PtCone_barrel1_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Number_Events_PtCone5_barrel1_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist("Number_Events_barrel1_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_barrel1_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
          FillHist("Number_Events_barrel1_"+PtConeRange+"_"+regions_mu.at(it_rg), 1.5, weight, 2, 0., 2.);
        }
      }

      // outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){
        FillHist("Mu_loose_PtCone_barrel2_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Number_Events_PtCone5_barrel2_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist("Number_Events_barrel2_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_barrel2_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
          FillHist("Number_Events_barrel2_"+PtConeRange+"_"+regions_mu.at(it_rg), 1.5, weight, 2, 0., 2.);
        }
      }

      // endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.5){
        FillHist("Mu_loose_PtCone_endcap_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Number_Events_PtCone5_endcap_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(ptcone_mu > 10.) FillHist("Number_Events_endcap_"+PtConeRange+"_"+regions_mu.at(it_rg), 0.5, weight, 2, 0., 2.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_endcap_"+PtConeRange+"_"+regions_mu.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
          FillHist("Number_Events_endcap_"+PtConeRange+"_"+regions_mu.at(it_rg), 1.5, weight, 2, 0., 2.);
        }
      }
    }

    // DY control region
    if(muons_tight.size()==2 && it_rg==1){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      ZCand = muons_tight.at(0) + muons_tight.at(1);
      if(ev.PassTrigger(MuonTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi1;
        FillHist("Lep1_Pt_nocut_Mu3_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Mu3_"+regions_mu.at(it_rg), muons_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Mu3_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi2;
        FillHist("Lep1_Pt_nocut_Mu8_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Mu8_"+regions_mu.at(it_rg), muons_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Mu8_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi3;
        FillHist("Lep1_Pt_nocut_Mu17_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Mu17_"+regions_mu.at(it_rg), muons_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Mu17_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }

      // event selection
      if(muons_tight.at(0).Pt() < 20. || muons_tight.at(1).Pt() < 10.) continue;
      if(fabs(ZCand.M() - MZ) > 10.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi1;
        FillHist("ZCand_Mass_Mu3_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Mu3_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi2;
        FillHist("ZCand_Mass_Mu8_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Mu8_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi3;
        FillHist("ZCand_Mass_Mu17_"+regions_mu.at(it_rg), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Mu17_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
    }

    // W+jets control region
    if(muons_tight.size()==1 && it_rg==2){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      Mt = MT(muons_tight.at(0), METv);
      if(ev.PassTrigger(MuonTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi1;
        FillHist("Lep_Pt_nocut_Mu3_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Mu3_"+regions_mu.at(it_rg), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Mu3_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi2;
        FillHist("Lep_Pt_nocut_Mu8_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Mu8_"+regions_mu.at(it_rg), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Mu8_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi3;
        FillHist("Lep_Pt_nocut_Mu17_"+regions_mu.at(it_rg), muons_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Mu17_"+regions_mu.at(it_rg), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Mu17_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
      }

      // event selection
      if(muons_tight.at(0).Pt() < 20.) continue;
      if(MET < 40.) continue;
      if(Mt < 60. || Mt > 100.) continue;
      
      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi1;
        FillHist("Mt_Mu3_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Mu3_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi2;
        FillHist("Mt_Mu8_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Mu8_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = MuonLumi3;
        FillHist("Mt_Mu17_"+regions_mu.at(it_rg), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Mu17_"+regions_mu.at(it_rg), 0.5, weight*triglumi, 2, 0., 2.);
      }
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Electron ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
 
  for(unsigned int it_rg2=0; it_rg2<regions_el.size(); it_rg2++){
    weight = 1.;

    // Fake rate measurement region 
    if(ele_loose.size()==1 && it_rg2==0){
      if(jets.size() == 0) continue;
      el_tight_iso = 0.0287+0.506/ele_loose.at(0).Pt();
      if(fabs(ele_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/ele_loose.at(0).Pt(); 
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

/*      if(ev.PassTrigger(ElectronTrig1) && ele_loose.at(0).Pt() > ElectronPtCut1){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi1;
        FillHist("Ele_loose_Eta_Ele8_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2) FillHist("Ele_loose_Eta_Ele8_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//        if(!ev.PassTrigger(ElectronTrig2) && !ev.PassTrigger(ElectronTrig3)){
          FillHist("Ele_loose_Eta_Ele8only_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//          if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2) FillHist("Ele_loose_Eta_Ele8only_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        }
      }
      if(ev.PassTrigger(ElectronTrig2) && ele_loose.at(0).Pt() > ElectronPtCut2){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi2;
        FillHist("Ele_loose_Eta_Ele12_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3) FillHist("Ele_loose_Eta_Ele12_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//        if(!ev.PassTrigger(ElectronTrig1) && !ev.PassTrigger(ElectronTrig3)){
          FillHist("Ele_loose_Eta_Ele12only_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//          if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3) FillHist("Ele_loose_Eta_Ele12only_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        }
      }
      if(ev.PassTrigger(ElectronTrig3) && ele_loose.at(0).Pt() > ElectronPtCut3){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi3;
        FillHist("Ele_loose_Eta_Ele17_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_el >= ElectronPtconeCut3) FillHist("Ele_loose_Eta_Ele17_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
      }
      if(ev.PassTrigger(ElectronTrig4) && ele_loose.at(0).Pt() > ElectronPtCut4){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi4;
        FillHist("Ele_loose_Eta_Ele23_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
        if(ptcone_el >= ElectronPtconeCut4) FillHist("Ele_loose_Eta_Ele23_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//         if(!ev.PassTrigger(ElectronTrig1) && !ev.PassTrigger(ElectronTrig2)){
           FillHist("Ele_loose_Eta_Ele23only_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
//           if(ptcone_el >= ElectronPtconeCut3) FillHist("Ele_loose_Eta_Ele23only_PtConeCut_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight*triglumi, 50, -2.5, 2.5);
         }
      }*/

      triglumi = 1.;
      // only 1 prescaled trigger for each PtCone range, setup lumi
      if(ptcone_el < ElectronPtconeCut1) continue;
      if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2){
        if(ele_loose.at(0).Pt() < ElectronPtCut1) continue;
        if(!ev.PassTrigger(ElectronTrig1)) continue;
        if(!IsDATA) triglumi = ElectronLumi1;
        PtConeRange = "range0";
      }
      if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3){
        if(ele_loose.at(0).Pt() < ElectronPtCut2) continue;      
        if(!ev.PassTrigger(ElectronTrig2)) continue;
        if(!IsDATA) triglumi = ElectronLumi2;
        PtConeRange = "range1";
      }
      if(ptcone_el >= ElectronPtconeCut3){
        if(ele_loose.at(0).Pt() < ElectronPtCut3) continue;
        if(!ev.PassTrigger(ElectronTrig3)) continue;
        if(!IsDATA) triglumi = ElectronLumi3;
        PtConeRange = "range2";
      }
      if(ptcone_el >= ElectronPtconeCut4){
        if(ele_loose.at(0).Pt() < ElectronPtCut4) continue;
        if(!ev.PassTrigger(ElectronTrig4)) continue;
        if(!IsDATA) triglumi = ElectronLumi4;
        PtConeRange = "range3";
      }

      awayjet = 0, leadingjet = 0;
      weight *= triglumi;

      FillHist("Ele_loose_Eta_nodijet_"+PtConeRange+"_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);        
      
      for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        dphi = fabs(jets.at(ijet).Phi() - ele_loose.at(0).Phi());
        if(dphi > pi) dphi = 2.*pi-dphi;
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
      FillHist("MET_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), MET, weight, 200, 0., 200.);
      FillHist("Mt_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), Mt, weight, 200, 0., 200.);
      FillHist("Ptratio_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), Pt_ratio, weight, 50, 0., 5.);
      FillHist("Ele_loose_PtCone_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
      FillHist("Ele_loose_Eta_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Jet_ChargedEmEnergyFraction_"+PtConeRange+"_"+regions_el.at(it_rg2), jet_emfraction, weight, 100, 0., 1.);
      FillHist("Number_Events_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), 0.5, weight, 2, 0., 2.); 

      if(ele_tight.size() > 0){
        FillHist("Ele_tight_PtCone_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
        FillHist("Ele_tight_Eta_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), ele_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist("Number_Events_nocut_"+PtConeRange+"_"+regions_el.at(it_rg2), 1.5, weight, 2, 0., 2.);
      }

      // additional cuts to reduce prompt contribution
      if(MET > 80.) continue;
      if(Mt > 25.) continue;
      if(Pt_ratio < 1.) continue;
      if(jet_emfraction > 0.65) continue;

      // Histograms after applying cuts
      FillHist("Ele_loose_PtCone_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
      FillHist("Ele_loose_Eta_"+PtConeRange+"_"+regions_el.at(it_rg2), ele_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Number_Events_"+PtConeRange+"_"+regions_el.at(it_rg2), 0.5, weight, 2, 0., 2.);
      if(ele_tight.size() > 0){
        FillHist("Ele_tight_PtCone_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
        FillHist("Ele_tight_Eta_"+PtConeRange+"_"+regions_el.at(it_rg2), ele_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist("Number_Events_"+PtConeRange+"_"+regions_el.at(it_rg2), 1.5, weight, 2, 0., 2.);
      }

      // inner barrel ( |eta| < 0.8 )
      if(fabs(ele_loose.at(0).Eta()) < 0.8){
        FillHist("Ele_loose_PtCone_barrel1_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
        FillHist("Number_Events_barrel1_"+PtConeRange+"_"+regions_el.at(it_rg2), 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist("Ele_tight_PtCone_barrel1_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
          FillHist("Number_Events_barrel1_"+PtConeRange+"_"+regions_el.at(it_rg2), 1.5, weight, 2, 0., 2.);
        }
      }

      // outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(ele_loose.at(0).Eta()) >= 0.8 && fabs(ele_loose.at(0).Eta()) < 1.479){
        FillHist("Ele_loose_PtCone_barrel2_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
        FillHist("Number_Events_barrel2_"+PtConeRange+"_"+regions_el.at(it_rg2), 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist("Ele_tight_PtCone_barrel2_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
          FillHist("Number_Events_barrel2_"+PtConeRange+"_"+regions_el.at(it_rg2), 1.5, weight, 2, 0., 2.);
        }
      }

      // endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(ele_loose.at(0).Eta()) >= 1.479 && fabs(ele_loose.at(0).Eta()) < 2.5){
        FillHist("Ele_loose_PtCone_endcap_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
        FillHist("Number_Events_endcap_"+PtConeRange+"_"+regions_el.at(it_rg2), 0.5, weight, 2, 0., 2.);
        if(ele_tight.size() > 0){
          FillHist("Ele_tight_PtCone_endcap_"+PtConeRange+"_"+regions_el.at(it_rg2), ptcone_el, weight, 200, 0., 200.);
          FillHist("Number_Events_endcap_"+PtConeRange+"_"+regions_el.at(it_rg2), 1.5, weight, 2, 0., 2.);
        }
      }
    }

    // DY control region
    if(ele_tight.size()==2 && it_rg2==1){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      ZCand = ele_tight.at(0) + ele_tight.at(1);
      if(ev.PassTrigger(ElectronTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi1;
        FillHist("Lep1_Pt_nocut_Ele8_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Ele8_"+regions_el.at(it_rg2), ele_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Ele8_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi2;
        FillHist("Lep1_Pt_nocut_Ele12_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Ele12_"+regions_el.at(it_rg2), ele_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Ele12_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi3;
        FillHist("Lep1_Pt_nocut_Ele17_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Ele17_"+regions_el.at(it_rg2), ele_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Ele17_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi4;
        FillHist("Lep1_Pt_nocut_Ele23_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("Lep2_Pt_nocut_Ele23_"+regions_el.at(it_rg2), ele_tight.at(1).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("ZCand_Mass_nocut_Ele23_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 80, 50., 130.);
      }

      // event selection
      if(ele_tight.at(0).Pt() < 25. || ele_tight.at(1).Pt() < 15.) continue;
      if(fabs(ZCand.M() - MZ) > 10.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi1;
        FillHist("ZCand_Mass_Ele8_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Ele8_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi2;
        FillHist("ZCand_Mass_Ele12_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Ele12_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi3;
        FillHist("ZCand_Mass_Ele17_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Ele17_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi4;
        FillHist("ZCand_Mass_Ele23_"+regions_el.at(it_rg2), ZCand.M(), weight*triglumi, 40, 70., 110.);
        FillHist("Number_Events_Ele23_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
    }

    // W+jets control region
    if(ele_tight.size()==1 && it_rg2==2){
      if(!IsDATA){
        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      Mt = MT(ele_tight.at(0), METv);
      if(ev.PassTrigger(ElectronTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi1;
        FillHist("Lep_Pt_nocut_Ele8_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Ele8_"+regions_el.at(it_rg2), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Ele8_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi2;
        FillHist("Lep_Pt_nocut_Ele12_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Ele12_"+regions_el.at(it_rg2), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Ele12_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi3;
        FillHist("Lep_Pt_nocut_Ele17_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Ele17_"+regions_el.at(it_rg2), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Ele17_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi4;
        FillHist("Lep_Pt_nocut_Ele23_"+regions_el.at(it_rg2), ele_tight.at(0).Pt(), weight*triglumi, 200, 0., 200.);
        FillHist("MET_nocut_Ele23_"+regions_el.at(it_rg2), MET, weight*triglumi, 200, 0., 200.);
        FillHist("Mt_nocut_Ele23_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
      }

      // event selection
      if(ele_tight.at(0).Pt() < 25.) continue;
      if(MET < 40.) continue;
      if(Mt < 60. || Mt > 100.) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi1;
        FillHist("Mt_Ele8_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Ele8_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig2)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi2;
        FillHist("Mt_Ele12_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Ele12_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig3)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi3;
        FillHist("Mt_Ele17_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Ele17_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
      if(ev.PassTrigger(ElectronTrig4)){
        triglumi = 1.;
        if(!IsDATA) triglumi = ElectronLumi4;
        FillHist("Mt_Ele23_"+regions_el.at(it_rg2), Mt, weight*triglumi, 200, 0., 200.);
        FillHist("Number_Events_Ele23_"+regions_el.at(it_rg2), 0.5, weight*triglumi, 2, 0., 2.);
      }
    }
  }
}


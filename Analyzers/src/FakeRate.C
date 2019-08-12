#include "FakeRate.h"

FakeRate::FakeRate(){

}

void FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;

//  MuonIDs = { "POGTight" };
//  MuonIDSFKeys = { "NUM_TightID_DEN_genTracks" };

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/FakeRate.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
  MuonTriggers.push_back("HLT_Mu3_PFJet40_v");      // SingleMuon
  MuonTriggers.push_back("HLT_Mu8_TrkIsoVVL_v");    // DoubleMuon
  MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_v");   // DoubleMuon
  MuonPtCut1 = 5., MuonPtCut2 = 10., MuonPtCut3 = 20.;
  MuonPtconeCut1 = 5., MuonPtconeCut2 = 20., MuonPtconeCut3 = 30.;

  ElectronTriggers.push_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  ElectronTriggers.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  ElectronTriggers.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  ElectronPtCut1 = 9.5, ElectronPtCut2 = 15., ElectronPtCut3 = 25.;
  ElectronPtconeCut1 = 10., ElectronPtconeCut2 = 23., ElectronPtconeCut3 = 35.;
 
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
//  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/FakeRate.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

//  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

  TString MuonID = "HNTight";
//  TString MuonIDSFKey = "NUM_TightID_DEN_genTracks";

  param.Clear();

  param.syst_ = AnalyzerParameter::Central;

  param.Name = MuonID+"_"+"Central";

/*  param.Electron_Tight_ID = "";
  param.Electron_Veto_ID = "";
  param.Electron_ID_SF_Key = "";*/
  param.Muon_Tight_ID = MuonID;
  param.Muon_Loose_ID = "HNLoose";
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
//  }
}

void FakeRate::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> regions = {"muonFR", "muonDY", "muonWj"};
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
  if(! (ev.PassTrigger(MuonTriggers) )) return;



  //======================
  //==== Copy AllObjects
  //======================

//  vector<Electron> this_AllElectrons = AllElectrons;
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

//  vector<Electron> electrons = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
//  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);
  vector<Muon> muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  vector<Muon> muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, 5., 2.4);
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
  vector<Gen> gens = GetGens(); 
//  std::vector<Lepton*> leptons;

  //=======================
  //==== Sort in pt-order
  //=======================

//  std::sort(electrons.begin(), electrons.end(), PtComparing);
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
  double pi = 3.14159265358979323846;
  double MZ = 91.1876;
  double MET = ev.GetMETVector().Pt(); 
  double dphi = 0.;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;

  double triggerlumi = 0.; 
  double awayjet_ptcut = 40.;
  int awayjet = 0;
  int leadingjet = 0;

  double ptcone_mu = 0.;
  TString PtConeRange = "";
  Particle ZCand;

/*  Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2); */

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1.;

    if(muons_loose.size()==1 && it_rg==0){
      ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
      // only 1 prescaled trigger for each PtCone range, setup lumi
      if(ptcone_mu < MuonPtconeCut1) continue; 
      if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2){
        if(muons_loose.at(0).Pt() < MuonPtCut1) continue;
        if(!(ev.PassTrigger("HLT_Mu3_PFJet40_v") && !ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v") && !ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v"))) continue;
        triggerlumi = 7.408;                     // private : 7.299
        if(DataYear==2017) triggerlumi = 4.612;  // private : 4.667
        awayjet_ptcut = 50.;
        PtConeRange = "range0";
      }
      else if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3){
        if(muons_loose.at(0).Pt() < MuonPtCut2) continue;
        if(!(!ev.PassTrigger("HLT_Mu3_PFJet40_v") && ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v") && !ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v"))) continue;
        triggerlumi = 7.801;                     // private : 7.618
        if(DataYear==2017) triggerlumi = 2.932;  // calculated privately
        PtConeRange = "range1";
      }
      else{
        if(muons_loose.at(0).Pt() < MuonPtCut3) continue;
        if(!(!ev.PassTrigger("HLT_Mu3_PFJet40_v") && !ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v") && ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v"))) continue;
        triggerlumi = 216.748;                   // private : 210.097
        if(DataYear==2017) triggerlumi = 66.625; // calculated privately
        PtConeRange = "range1";
      }

      if(!IsDATA){
        // Gen matching
        Gen truth_lep = GetGenMatchedLepton(muons_loose.at(0), gens);
        if(truth_lep.PID() == 0) continue; // No matched gen lepton -> PID == 0 (See DataFormats/src/Gen.C)
      
        // weight
        weight *= weight_norm_1invpb*triggerlumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        // define dphi
        dphi = fabs(jets.at(ijet).Eta()-muons_loose.at(0).Eta());
        if(dphi > pi) dphi = 2*pi-dphi;

        if(dphi > 2.5) awayjet++; 
        if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      }

      if(awayjet == 0) continue;
      if(jets.at(leadingjet).Pt() < awayjet_ptcut) continue;
 
      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets.at(leadingjet).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist("MET_nocut_"+PtConeRange+"_"+regions.at(it_rg), MET, weight, 200, 0., 200.);
      FillHist("Mt_nocut_"+PtConeRange+"_"+regions.at(it_rg), Mt, weight, 200, 0., 200.);
      FillHist("Ptratio_nocut_"+PtConeRange+"_"+regions.at(it_rg), Pt_ratio, weight, 30, 0., 3.);
      FillHist("Mu_loose_PtCone_nocut_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
      FillHist("Mu_loose_Eta_nocut_"+PtConeRange+"_"+regions.at(it_rg), muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      if(muons_tight.size() > 0){
        FillHist("Mu_tight_PtCone_nocut_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Mu_tight_Eta_nocut_"+PtConeRange+"_"+regions.at(it_rg), muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
      }

      // cuts on kinematic variables
      if(MET > 80.) continue;
      if(Mt > 25.) continue;
      if(Pt_ratio < 1.) continue;

      // Histograms after applying cuts
      FillHist("Mu_loose_PtCone_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
      FillHist("Mu_loose_Eta_"+PtConeRange+"_"+regions.at(it_rg), muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5); 
      if(muons_tight.size() > 0){
        FillHist("Mu_tight_PtCone_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        FillHist("Mu_tight_Eta_"+PtConeRange+"_"+regions.at(it_rg), muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
      }

      // inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        FillHist("Mu_loose_PtCone_barrel1_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_barrel1_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        }
      }

      // outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){
        FillHist("Mu_loose_PtCone_barrel2_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_barrel2_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        }
      }

      // endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.5){
        FillHist("Mu_loose_PtCone_endcap_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        if(muons_tight.size() > 0){
          FillHist("Mu_tight_PtCone_endcap_"+PtConeRange+"_"+regions.at(it_rg), ptcone_mu, weight, 200, 0., 200.);
        }
      }
    }

    if(muons_tight.size()==2 && it_rg==1){
      ZCand = muons_tight.at(0) + muons_tight.at(1);
      FillHist("Mu1_Pt_"+regions.at(it_rg), muons_tight.at(0).Pt(), weight, 200, 0., 200.);
      FillHist("Mu2_Pt_"+regions.at(it_rg), muons_tight.at(1).Pt(), weight, 200, 0., 200.);
      FillHist("ZCand_Mass_"+regions.at(it_rg), ZCand.M(), weight, 80, 50., 130.);
 
      if(muons_tight.at(0).Pt() < 20. || muons_tight.at(1).Pt() < 10.) continue;
      if(fabs(ZCand.M() - MZ) > 20.) continue;

      if(!IsDATA){
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }

      // weight for each trigger
      if(ev.PassTrigger("HLT_Mu3_PFJet40_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*7.408;
          if(DataYear == 2017) weight *= weight_norm_1invpb*4.612; 
        }
        FillHist("ZCand_Mass_Mu3_"+regions.at(it_rg), ZCand.M(), weight, 40, 70., 110.);
        FillHist("Number_Events_Mu3_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
      if(ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*7.801;
          if(DataYear == 2017) weight *= weight_norm_1invpb*2.932;
        }
        FillHist("ZCand_Mass_Mu8_"+regions.at(it_rg), ZCand.M(), weight, 40, 70., 110.);
        FillHist("Number_Events_Mu8_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
      if(ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*216.748;
          if(DataYear == 2017) weight *= weight_norm_1invpb*66.625;
        }
        FillHist("ZCand_Mass_Mu17_"+regions.at(it_rg), ZCand.M(), weight, 40, 70., 110.);
        FillHist("Number_Events_Mu17_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
    }

    if(muons_tight.size()==1 && it_rg==2){
      Mt = MT(muons_tight.at(0), METv);
      FillHist("Mu1_Pt_"+regions.at(it_rg), muons_tight.at(0).Pt(), weight, 200, 0., 200.);
      FillHist("MET_"+regions.at(it_rg), MET, weight, 200, 0., 200.);
      FillHist("Mt_"+regions.at(it_rg), Mt, weight, 200, 0., 200.);

      if(muons_tight.at(0).Pt() < 20.) continue;
      if(MET < 40.) continue;
      if(Mt < 60. || Mt > 100.) continue;

      if(!IsDATA){
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      }
      
      // weight for each trigger
      if(ev.PassTrigger("HLT_Mu3_PFJet40_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*7.408;
          if(DataYear == 2017) weight *= weight_norm_1invpb*4.612; 
        }
        FillHist("Mt_Mu3_"+regions.at(it_rg), Mt, weight, 200, 0., 200.);
        FillHist("Number_Events_Mu3_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
      if(ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*7.801;
          if(DataYear == 2017) weight *= weight_norm_1invpb*2.932;
        }
        FillHist("Mt_Mu8_"+regions.at(it_rg), Mt, weight, 200, 0., 200.);
        FillHist("Number_Events_Mu8_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
      if(ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v")){
        if(!IsDATA){
          if(DataYear == 2016) weight *= weight_norm_1invpb*216.748;
          if(DataYear == 2017) weight *= weight_norm_1invpb*66.625;
        }
        FillHist("Mt_Mu17_"+regions.at(it_rg), Mt, weight, 200, 0., 200.);
        FillHist("Number_Events_Mu17_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
      }
    }
  }
}


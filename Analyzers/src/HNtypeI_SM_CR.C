#include "HNtypeI_SM_CR.h"

HNtypeI_SM_CR::HNtypeI_SM_CR(){

}

void HNtypeI_SM_CR::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[HNtypeI_SM_CR::initializeAnalyzer] RunSyst = " << RunSyst << endl;

//  MuonIDs = { "POGTight" };
//  MuonIDSFKeys = { "NUM_TightID_DEN_genTracks" };

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SM_CR.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
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

//  cout << "[HNtypeI_SM_CR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[HNtypeI_SM_CR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== Test btagging code
  //==== add taggers and WP that you want to use in analysis
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);
//  vtaggers.push_back(Jet::CSVv2);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);
//  v_wps.push_back(Jet::Loose);

  //=== list of taggers, WP, setup systematics, use period SFs
  SetupBTagger(vtaggers, v_wps, true, true);

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
  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SM_CR.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

//  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

//  TString MuonID = "POGTightWithTightIso";  // POG ID
  TString MuonID = "HNTight2016";               // HN ID

  param.Clear();

  param.syst_ = AnalyzerParameter::Central;

  param.Name = MuonID+"_"+"Central";

  // POG ID
//  param.Electron_Tight_ID = "passTightID";
/*  param.Electron_Tight_ID = "TightWithIPcut";
  param.Electron_Veto_ID = "passVetoID";
  param.Electron_ID_SF_Key = "passTightID";
  param.Electron_Trigger_SF_Key = "";
  param.Muon_Tight_ID = MuonID;
  param.Muon_Veto_ID = "POGLoose";
  param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  param.Muon_Trigger_SF_Key = "";*/

  // HN ID
  param.Electron_Tight_ID = "HNTight2016";
  param.Electron_Veto_ID = "HNVeto2016";
  param.Electron_ID_SF_Key = "";
  param.Electron_Trigger_SF_Key = "";
  param.Muon_Tight_ID = MuonID;
  param.Muon_Veto_ID = "HNVeto2016";
  param.Muon_ID_SF_Key = "";
  param.Muon_ISO_SF_Key = "";
  param.Muon_Trigger_SF_Key = "";

  // Jet ID
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

void HNtypeI_SM_CR::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> regions = {"WZ", "ZG", "ZZ", "DY"};
  vector<TString> channels = {"eee", "eem", "emm", "mmm"};
  TString pair_sign, DYchannel;
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  //=============
  //==== No Cut
  //=============

//  JSFillHist(param.Name, "NoCut_"+param.Name, 0., 1., 1, 0., 1.);
  FillHist("Number_Events_WZ", 0.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZG", 0.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZZ", 0.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_OS_DY", 0.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_SS_DY", 0.5, 1., cutflow_bin, 0., cutflow_max);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;
  FillHist("Number_Events_WZ", 1.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZG", 1.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZZ", 1.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_OS_DY", 1.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_SS_DY", 1.5, 1., cutflow_bin, 0., cutflow_max);

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  //==============
  //==== Trigger
  //==============

  if(!IsDATA && !(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers))) return;
  if(IsDATA){
/*    if(DataStream.Contains("SingleMuon") && !ev.PassTrigger(MuonTriggers)) return;
    else if(DataStream.Contains("SingleElectron")){
      if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
    }*/
    if(DataYear==2016 || DataYear==2017){
      if(DataStream.Contains("DoubleMuon") && !ev.PassTrigger(MuonTriggers)) return;
      else if(DataStream.Contains("DoubleEG")){
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
    if(DataYear==2018){
      if(DataStream.Contains("SignleMuon") && !ev.PassTrigger(MuonTriggers)) return;
      else if(DataStream.Contains("EGamma")){
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
  }
  FillHist("Number_Events_WZ", 2.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZG", 2.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_ZZ", 2.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_OS_DY", 2.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("Number_Events_SS_DY", 2.5, 1., cutflow_bin, 0., cutflow_max);

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
    cout << "[HNtypeI_SM_CR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Electron> electrons = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);
  vector<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);
  std::vector<Lepton*> leptons;
  std::vector<Lepton*> leptons_minus;
  std::vector<Lepton*> leptons_plus;

  //=======================
  //==== Sort in pt-order
  //=======================

  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);

//  int n_bjet_deepcsv_m=0;
//  int n_bjet_deepcsv_m_noSF=0;
  int n_bjet_deepcsv=0;

  for(unsigned int ij=0; ij<jets.size(); ij++){
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) n_bjet_deepcsv_m++; // method for getting btag with SF applied to MC
//    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,false,0)) n_bjet_deepcsv_m_noSF++; // method for getting btag with no SF applied to MC
//    if(IsBTagged(jets.at(ij), Jet::CSVv2, Jet::Loose,false,0)) n_bjet_csvv2++;
    if(fabs(jets.at(ij).Eta()) > 2.4) continue;
    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium, true, 0)) n_bjet_deepcsv++;
  }
  
  //=========================
  //==== Event selections..
  //=========================

  int lepton_size = electrons.size() + muons.size();
  int lepton_veto_size = electrons_veto.size() - electrons.size() + muons_veto.size() - muons.size();
  double MET = ev.GetMETVector().Pt();
  double Mt = 0.;
  double MZ = 91.1876;
  double weight = 1.;
  double muon_idsf = 1.;
  double muon_isosf = 1.;
  double ele_idsf = 1.;
  double ele_recosf = 1.;
  int ossf_mass10 = 0;

  // Sort leptons with pT order, separate leptons with their charge
  if(lepton_size < 5){
    for(unsigned int imu=0; imu<muons.size(); imu++){
      leptons.push_back(&muons.at(imu));
    }
    for(unsigned int iel=0; iel<electrons.size(); iel++){
      leptons.push_back(&electrons.at(iel));
    }
    std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

    for(unsigned int it_mu2=0; it_mu2<muons.size(); it_mu2++){
      if(muons.at(it_mu2).Charge() < 0) leptons_minus.push_back(&muons.at(it_mu2));
      if(muons.at(it_mu2).Charge() > 0) leptons_plus.push_back(&muons.at(it_mu2));
    }
    for(unsigned int it_el2=0; it_el2<electrons.size(); it_el2++){
      if(electrons.at(it_el2).Charge() < 0) leptons_minus.push_back(&electrons.at(it_el2));
      if(electrons.at(it_el2).Charge() > 0) leptons_plus.push_back(&electrons.at(it_el2));
    }    
  }

  Particle ZCand, WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp, ZCand1, ZCand2;

  // Event selections for each CR
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
    ossf_mass10 = 0;
    
    // WZ or ZG CR
    if(lepton_size==3 && it_rg<2){
      FillHist("Number_Events_"+regions.at(it_rg), 3.5, 1., cutflow_bin, 0., cutflow_max);

      // Veto additional leptons
      if(lepton_veto_size>0) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 4.5, 1., cutflow_bin, 0., cutflow_max);

      // Triggers
      if(ev.PassTrigger(MuonTriggers)){
        if(muons.size()<2) continue;
        else if(muons.size()>=2 && (muons.at(0).Pt()<MuonPtCut1 || muons.at(1).Pt()<MuonPtCut2)) continue;
      }
      else if(!ev.PassTrigger(MuonTriggers) && ev.PassTrigger(ElectronTriggers)){
        if(electrons.size()<2) continue;
        else if(electrons.size()>=2 && (electrons.at(0).Pt()<ElectronPtCut1 || electrons.at(1).Pt()<ElectronPtCut2)) continue;
      }

      FillHist("Number_Events_"+regions.at(it_rg), 5.5, 1., cutflow_bin, 0., cutflow_max);

      // sort 3 leptons with pT order
/*      for(unsigned int imu=0; imu<muons.size(); imu++){
        leptons.push_back(&muons.at(imu));
      }
      for(unsigned int iel=0; iel<electrons.size(); iel++){
        leptons.push_back(&electrons.at(iel));
      }
      std::sort(leptons.begin(), leptons.end(), PtComparingPtr);*/

      // Z candidate
//      Particle ZCand, WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp;
      int l1 = -999, l2 = -999, wlep = -999;  
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
//        if(muons.size()>2) leptons = MakeLeptonPointerVector(muons);
//        else leptons = MakeLeptonPointerVector(electrons);
        // No same-sign 3 leptons
        if(fabs(leptons.at(0)->Charge() + leptons.at(1)->Charge() + leptons.at(2)->Charge()) == 1){
          double tmpMassDiff=1000000.; 

          
          //ZCand = *leptons.at(0) + *leptons.at(1);
          //if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) ZCand = *leptons.at(0) + *leptons.at(2);
          // Select OSSF lepton pair
          for(int ilep1=0; ilep1<2; ilep1++){
            for(int ilep2=ilep1+1; ilep2<3; ilep2++){
              // require OS dilepton
              if(leptons.at(ilep1)->Charge()*leptons.at(ilep2)->Charge()>0) continue;
              Ztemp = *leptons.at(ilep1) + *leptons.at(ilep2);
              if(Ztemp.M() < 10) ossf_mass10++;
              if(fabs(Ztemp.M() - MZ) <  tmpMassDiff) {
                  tmpMassDiff= fabs(Ztemp.M() - MZ);
                  ZCand = Ztemp; l1 = ilep1; l2 = ilep2;
              }
            }
          }
          // Set W-tagged lepton
          for(int ilep3=0; ilep3<3; ilep3++){
            if(fabs(ilep3-l1)>0 && fabs(ilep3-l2)>0) { wlep = ilep3; break;}
          }
          WtagLep = *leptons.at(wlep);
          ZtagLep1 = *leptons.at(l1);
          ZtagLep2 = *leptons.at(l2);
        }
      }
      else continue;

      FillHist("Number_Events_"+regions.at(it_rg), 6.5, 1., cutflow_bin, 0., cutflow_max);

      TriLep = ZCand + WtagLep;
      Mt = MT(WtagLep, METv);

      // cuts on kinematic variables
      if(ZCand.M()<10. || ossf_mass10>0) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 7.5, 1., cutflow_bin, 0., cutflow_max);
      if(n_bjet_deepcsv > 0) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 8.5, 1., cutflow_bin, 0., cutflow_max);

      if(it_rg == 0){
        if(!IsOnZ(ZCand.M(), 15.)) continue;
        if(MET < 50.) continue;
        if(Mt < 20.) continue;
        if(TriLep.M() < MZ + 15.) continue;
      }
      if(it_rg == 1){
        if(IsOnZ(ZCand.M(), 15.)) continue;
        if(MET > 50.) continue;
        if(!IsOnZ(TriLep.M(), 15.)) continue;
      }

      FillHist("Number_Events_"+regions.at(it_rg), 9.5, 1., cutflow_bin, 0., cutflow_max);

      // weights for MC
      if(!IsDATA){
      //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
      //==== So, you have to multiply trigger luminosity
      //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"

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
      
      // Histograms
      FillHist("Number_Jets_"+regions.at(it_rg), jets.size(), weight, 8, 0., 8.);
      FillHist("ZCand_Mass_"+regions.at(it_rg), ZCand.M(), weight, 80, 50., 130.);
      FillHist("TriLep_Mass_"+regions.at(it_rg), TriLep.M(), weight, 80, 50., 130.);
      FillHist("WtagLep_Pt_"+regions.at(it_rg), WtagLep.Pt(), weight, 400, 0., 400.);
      FillHist("ZtagLep1_Pt_"+regions.at(it_rg), ZtagLep1.Pt(), weight, 400, 0., 400.);
      FillHist("ZtagLep2_Pt_"+regions.at(it_rg), ZtagLep2.Pt(), weight, 400, 0., 400.);
      FillHist("Lep1_Pt_"+regions.at(it_rg), leptons.at(0)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep2_Pt_"+regions.at(it_rg), leptons.at(1)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep3_Pt_"+regions.at(it_rg), leptons.at(2)->Pt(), weight, 400, 0., 400.);
      FillHist("MET_"+regions.at(it_rg), MET, weight, 400, 0., 400.);
      FillHist("Mt_"+regions.at(it_rg), Mt, weight, 400, 0., 400.);
      FillHist("Number_Events_weighted_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
     
      for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
        if(muons.size() == it_ch){
          FillHist("Number_Jets_"+regions.at(it_rg)+"_"+channels.at(it_ch), jets.size(), weight, 8, 0., 8.);
          FillHist("ZCand_Mass_"+regions.at(it_rg)+"_"+channels.at(it_ch), ZCand.M(), weight, 80, 50., 130.);
          FillHist("TriLep_Mass_"+regions.at(it_rg)+"_"+channels.at(it_ch), TriLep.M(), weight, 80, 50., 130.);
          FillHist("WtagLep_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), WtagLep.Pt(), weight, 400, 0., 400.);
          FillHist("ZtagLep1_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), ZtagLep1.Pt(), weight, 400, 0., 400.);
          FillHist("ZtagLep2_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), ZtagLep2.Pt(), weight, 400, 0., 400.);
          FillHist("Lep1_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), leptons.at(0)->Pt(), weight, 400, 0., 400.);
          FillHist("Lep2_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), leptons.at(1)->Pt(), weight, 400, 0., 400.);
          FillHist("Lep3_Pt_"+regions.at(it_rg)+"_"+channels.at(it_ch), leptons.at(2)->Pt(), weight, 400, 0., 400.);
          FillHist("MET_"+regions.at(it_rg)+"_"+channels.at(it_ch), MET, weight, 400, 0., 400.);
          FillHist("Mt_"+regions.at(it_rg)+"_"+channels.at(it_ch), Mt, weight, 400, 0., 400.);
          FillHist("Number_Events_weighted_"+regions.at(it_rg)+"_"+channels.at(it_ch), 0.5, weight, 2, 0., 2.);
        }
      }
    }

    // ZZ CR
    if(lepton_size==4 && it_rg==2){

      FillHist("Number_Events_"+regions.at(it_rg), 3.5, 1., cutflow_bin, 0., cutflow_max);

      // Veto additional leptons
      if(lepton_veto_size>0) continue;

      FillHist("Number_Events_"+regions.at(it_rg), 4.5, 1., cutflow_bin, 0., cutflow_max);

      // Triggers
      if(ev.PassTrigger(MuonTriggers)){
        if(muons.size()<2) continue;
        else if(muons.size()>=2 && (muons.at(0).Pt()<MuonPtCut1 || muons.at(1).Pt()<MuonPtCut2)) continue;
      }
      else if(!ev.PassTrigger(MuonTriggers) && ev.PassTrigger(ElectronTriggers)){
        if(electrons.size()<2) continue;
        else if(electrons.size()>=2 && (electrons.at(0).Pt()<ElectronPtCut1 || electrons.at(1).Pt()<ElectronPtCut2)) continue;
      }

      FillHist("Number_Events_"+regions.at(it_rg), 5.5, 1., cutflow_bin, 0., cutflow_max);

      // sort 4 leptons with pT order
/*      for(unsigned int imu=0; imu<muons.size(); imu++){
        leptons.push_back(&muons.at(imu));
      }
      for(unsigned int iel=0; iel<electrons.size(); iel++){
        leptons.push_back(&electrons.at(iel));
      }
      std::sort(leptons.begin(), leptons.end(), PtComparingPtr);*/

      // Select OSSF lepton pairs
//      Particle ZCand1, ZCand2;
      if(muons.size()==2 && muons.at(0).Charge()*muons.at(1).Charge()<0 && electrons.at(0).Charge()*electrons.at(1).Charge()<0){
        ZCand1 = muons.at(0) + muons.at(1);
        ZCand2 = electrons.at(0) + electrons.at(1);
      }
      else if(muons.size()==4 || electrons.size()==4){
/*        if(muons.size()==4){
          for(unsigned int it_mu=0; it_mu<4; it_mu++){
            if(muons.at(it_mu).Charge() < 0) leptons_minus.push_back(&muons.at(it_mu));
            if(muons.at(it_mu).Charge() > 0) leptons_plus.push_back(&muons.at(it_mu));
          }
        }
        if(electrons.size()==4){
          for(unsigned int it_el=0; it_el<4; it_el++){
            if(electrons.at(it_el).Charge() < 0) leptons_minus.push_back(&electrons.at(it_el));
            if(electrons.at(it_el).Charge() > 0) leptons_plus.push_back(&electrons.at(it_el));
          }
        }*/

        if(leptons_minus.size() == leptons_plus.size()){
          ZCand1 = *leptons_minus.at(0) + *leptons_plus.at(0);
          ZCand2 = *leptons_minus.at(1) + *leptons_plus.at(1);
          if(!IsOnZ(ZCand1.M(), 15.) || !IsOnZ(ZCand2.M(), 15.)){
            ZCand1 = *leptons_minus.at(0) + *leptons_plus.at(1);
            ZCand2 = *leptons_minus.at(1) + *leptons_plus.at(0);
          }
        }
      }
      else continue;

      FillHist("Number_Events_"+regions.at(it_rg), 6.5, 1., cutflow_bin, 0., cutflow_max);

      // cuts on kinematic variables
      if(ZCand1.M()<10. || ZCand2.M()<10.) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 7.5, 1., cutflow_bin, 0., cutflow_max);

      if(n_bjet_deepcsv > 0) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 8.5, 1., cutflow_bin, 0., cutflow_max);
   
      if(!IsOnZ(ZCand1.M(), 15.) || !IsOnZ(ZCand2.M(), 15.)) continue;
      FillHist("Number_Events_"+regions.at(it_rg), 9.5, 1., cutflow_bin, 0., cutflow_max);

      // weights for MC
      if(!IsDATA){
        weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//         muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          weight *= muon_idsf*muon_isosf;
        }
        for(unsigned int j=0; j<electrons.size(); j++){
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          weight *= ele_recosf*ele_idsf;
        }
      }
 
      // Histograms
      FillHist("Number_Jets_"+regions.at(it_rg), jets.size(), weight, 8, 0., 8.);
      FillHist("ZCand1_Mass_"+regions.at(it_rg), ZCand1.M(), weight, 80, 50., 130.);
      FillHist("ZCand2_Mass_"+regions.at(it_rg), ZCand2.M(), weight, 80, 50., 130.);
      FillHist("Lep1_Pt_"+regions.at(it_rg), leptons.at(0)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep2_Pt_"+regions.at(it_rg), leptons.at(1)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep3_Pt_"+regions.at(it_rg), leptons.at(2)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep4_Pt_"+regions.at(it_rg), leptons.at(3)->Pt(), weight, 400, 0., 400.);
      FillHist("MET_"+regions.at(it_rg), MET, weight, 400, 0., 400.);
      FillHist("Number_Events_weighted_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
    }

    // DY CR
    if(lepton_size==2 && it_rg==3){
      if(muons.size()==1 && electrons.size()==1) continue; 
      FillHist("Number_Events_OS_"+regions.at(it_rg), 3.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist("Number_Events_SS_"+regions.at(it_rg), 3.5, 1., cutflow_bin, 0., cutflow_max);

      // Triggers
      if(ev.PassTrigger(MuonTriggers)){
        if(muons.size()<2) continue;
        else if(muons.size()==2 && (muons.at(0).Pt()<MuonPtCut1 || muons.at(1).Pt()<MuonPtCut2)) continue;
      }
      if(ev.PassTrigger(ElectronTriggers)){
        if(electrons.size()<2) continue;
        else if(electrons.size()==2 && (electrons.at(0).Pt()<ElectronPtCut1 || electrons.at(1).Pt()<ElectronPtCut2)) continue;
      }      

      FillHist("Number_Events_OS_"+regions.at(it_rg), 4.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist("Number_Events_SS_"+regions.at(it_rg), 4.5, 1., cutflow_bin, 0., cutflow_max);

      ZCand = *leptons.at(0) + *leptons.at(1);

      if(!IsDATA){
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

      if(n_bjet_deepcsv > 0) continue;
      FillHist("Number_Events_OS_"+regions.at(it_rg), 5.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist("Number_Events_SS_"+regions.at(it_rg), 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(muons.size()==2){ DYchannel="mm"; }
      else{ DYchannel="ee"; }

      if(leptons.at(0)->Charge()*leptons.at(1)->Charge() < 0){ pair_sign = "OS"; }
      else{ pair_sign = "SS"; }      

      FillHist("Number_VetoLeptons"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), lepton_veto_size, weight, 5, 0., 5.);
      FillHist("Number_Jets_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), jets.size(), weight, 8, 0., 8.);
      FillHist("ZCand_Mass_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), ZCand.M(), weight, 400, 0., 400.);
      FillHist("ZCand_Pt_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), ZCand.Pt(), weight, 400, 0., 400.);
      FillHist("Lep1_Pt_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), leptons.at(0)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep2_Pt_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), leptons.at(1)->Pt(), weight, 400, 0., 400.);
      FillHist("Lep1_Eta_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("Lep2_Eta_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("MET_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), MET, weight, 400, 0., 400.);
      FillHist("Number_Events_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist("Number_Events_weighted_"+pair_sign+"_"+DYchannel+"_"+regions.at(it_rg), 0.5, weight, 2, 0., 2.);
    } 
  }
}




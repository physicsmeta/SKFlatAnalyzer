#include "Presel.h"

Presel::Presel(){

}

void Presel::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
//  RunSyst = HasFlag("RunSyst");
//  cout << "[Presel::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  RunFake = HasFlag("RunFake");
  RunCF = HasFlag("RunCF");

  MuonTightIDs = {"HNTightV1"};
  MuonLooseIDs = {"HNLooseV3"};
  MuonVetoIDs  = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV1"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};
  MuonFRNames      = {"HNRun2"};
  ElectronFRNames  = {"HNRun2"}; //JH : NOTE This is used in fakeEst->ReadHistograms() in m.initializeAnalyzerTools() 

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/Presel.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  MuonTriggers.clear();
  MuonTriggersH.clear();
  ElectronTriggers.clear();
  EMuTriggers.clear();
  EMuTriggersH.clear();

  if(DataYear==2016){                                                                   // Lumi values for trigger weight (/pb)
    MuonTriggersBtoG.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                   // 27267.591112919 
    MuonTriggersBtoG.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                 // 27267.591112919 //JH : NOTE these two are prescaled at 2016H -> https://its.cern.ch/jira/browse/CMSHLT-1002
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                       // 35918.219492947
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                     // 35918.219492947 
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                    // 35918.219492947
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                  // 35918.219492947 
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                   // 8650.628380028
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                 // 8650.628380028
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
    EMuTriggersBtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");      // 27267.591112919
    EMuTriggersBtoG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");     // 27267.591112919
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 35918.219492947
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 35918.219492947
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");       // 35918.219492947
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 35918.219492947
    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
    EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.; 
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2018){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

//  cout << "[Presel::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
//  cout << "[Presel::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps); //JH : NOTE This is used in mcCorr->SetupJetTagging() in m.initializeAnalyzerTools();
}

Presel::~Presel(){

  //==== Destructor of this Analyzer

}

void Presel::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/Presel.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
//  AllFatJets = GetAllFatJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets()); //JH : puppiCorr = new FakeBackgroundEstimator(); in the constructor of AnalyzerCore.C; apply correction to fatjet.SDMass(); the total weight = gen correction * reco correction, from SKFlatAnalyzer/data/Run2Legacy_v4/DataYear/PuppiSoftdropMassCorr/puppiCorr.root

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/Presel.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){
    TString MuonTightID = MuonTightIDs.at(it_id);
    TString MuonLooseID = MuonLooseIDs.at(it_id); 
    TString MuonVetoID  = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);
    TString MuonFRName      = MuonFRNames.at(it_id);
    TString ElectronFRName  = ElectronFRNames.at(it_id);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;

//    param.Name = MuonID+"_"+"Central";

    // Muon ID
    param.Muon_Tight_ID = MuonTightID;
    param.Muon_Loose_ID = MuonLooseID;
    param.Muon_Veto_ID  = MuonVetoID;
    param.Muon_FR_ID = MuonFRName;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key = "AwayJetPt40"; // histname
    //param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
    //param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
    param.Muon_ID_SF_Key = "";
    param.Muon_ISO_SF_Key = "";
    param.Muon_Trigger_SF_Key = "";
    param.Muon_UsePtCone = true;

    // Electron ID
    param.Electron_Tight_ID = ElectronTightID;
    param.Electron_Loose_ID = ElectronLooseID;
    param.Electron_Veto_ID  = ElectronVetoID;
    param.Electron_FR_ID = ElectronFRName;     // ID name in histmap_Electron.txt
    param.Electron_FR_Key = "AwayJetPt40"; // histname
    //param.Electron_ID_SF_Key = "passTightID";
    param.Electron_ID_SF_Key = "";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone = true;

    // Jet ID
//    param.Jet_ID = "tightLepVeto";
    param.Jet_ID = "HNTight";
    //param.FatJet_ID = "HNTight";
    if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45"; //JH : TODO

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

void Presel::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> channels;
  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) channels.push_back("dimu");
    else if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) channels.push_back("diel");
    else if(DataStream.Contains("MuonEG")) channels.push_back("emu");
  }
  else channels = {"dimu", "diel", "emu"};
  TString IDsuffix = "Run2";
  TString LepCategory = "TT";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;
 
  Event ev = GetEvent();

  //=============
  //==== No Cut
  //=============

  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full"); //JH : weight_norm_1invpb = xsec/sumW; Lumi = 35.9, 41.5, 59.7(fb-1) total 137fb-1
    weight *= ev.MCweight(); //JH : gen_weight in MiniAOD
    weight *= GetPrefireWeight(0); //JH : No issue in 2018, otherwise returns L1PrefireReweight_Central in MiniAOD
    weight *= GetPileUpWeight(nPileUp,0); //JH : mcCorr->GetPileUpWeight(N_pileup, syst); mcCorr->GetPileUpWeight2017(N_pileup, syst); NOTE 2018 not yet added.
  } //JH : total weight calculation done.

  // Cutflow : No Cuts
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  // Cutflow : MET filter
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //==============
  //==== Trigger
  //==============

  if(!(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers) || ev.PassTrigger(EMuTriggers))) return; 

  //======================
  //==== Copy AllObjects
  //======================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons;
  if(param.Muon_Tight_ID.Contains("HighPt")) this_AllMuons = UseTunePMuon(AllMuons);
  else this_AllMuons = AllMuons;
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
    cout << "[Presel::executeEventFromParameter] Wrong syst" << endl;
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

  vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 5., 2.4);
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 10., 2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5); //JH : lepton selection done
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7); //JH : to reject bjets
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
    lepton_count1 = 0;
    if(!(this_AllFatJets.at(i).PassID(param.FatJet_ID))) continue; //JH : "HNTight"
    if(!(this_AllFatJets.at(i).Pt() > 200.)) continue;
    if(!(fabs(this_AllFatJets.at(i).Eta()) < 2.7)) continue;
    for(unsigned int j=0; j<muons_veto.size(); j++){
      if(this_AllFatJets.at(i).DeltaR(muons_veto.at(j)) < 0.8) lepton_count1++; //JH : muon cleaning
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllFatJets.at(i).DeltaR(electrons_veto.at(j)) < 0.8) lepton_count1++; //JH : electron cleaning
    } 
    if(lepton_count1 > 0) continue;
    fatjets.push_back(this_AllFatJets.at(i));
  }

  for(unsigned int i=0; i<this_AllJets.size(); i++){
    lepton_count2 = 0, fatjet_count = 0;
    if(!(this_AllJets.at(i).PassID(param.Jet_ID))) continue; //JH :"HNTight"
    if(!(this_AllJets.at(i).Pt() > 20.)) continue;
    if(!(fabs(this_AllJets.at(i).Eta()) < 2.7)) continue;
    for(unsigned int j=0; j<muons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(muons_veto.at(j)) < 0.4) lepton_count2++; //JH : muon cleaning
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(electrons_veto.at(j)) < 0.4) lepton_count2++; //JH : electron cleaning
    }
    for(unsigned int j=0; j<fatjets.size(); j++){
      if(this_AllJets.at(i).DeltaR(fatjets.at(j)) < 1.0) fatjet_count++; //JH : fatjet cleaning
    }
    if(lepton_count2 > 0) continue;
    if(fatjet_count > 0) continue;
    jets.push_back(this_AllJets.at(i));
  }

//JH : jet, fatjet selection done.

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

  //==== B-Tagging 
  int Nbjet_loose = 0, Nbjet_medium = 0;
  JetTagging::Parameters jtp_DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb); //JH : Set b-tagging parameters

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
//  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){ 
//    double this_discr = jets_nolepveto.at(ij).GetTaggerResult(JetTagging::DeepCSV);
      //==== No SF
//      if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++; 
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++; //JH : count Nbjet. NOTE : AN says they used CVSv2 and medium WP.
  } 

//  FillHist("Nbjet_loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
//  FillHist("Nbjet_medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);

  //===================================
  //==== Set up pTcone, lepton vector
  //===================================

  Particle METv = ev.GetMETVector();
  METv = UpdateMETMuon(METv, muons);
  METv = UpdateMETElectron(METv, electrons);
  double MET = METv.Pt(); // JH : MET propagated

  double Mt = 0.;
  double Mt3l = 0.;
  double ST = 0.;
  double MET2ST = 0.;
  double MZ = 91.1876;
  double MW = 80.379;
  double muon_recosf = 1.;
  double muon_idsf = 1.;
  double muon_isosf = 1.;
  double ele_idsf = 1.;
  double ele_recosf = 1.;
  int lepton_veto_size = 0;
  double LeptonPtCut1 = 0., LeptonPtCut2 = 0.;
  Particle ZCand, Wtemp1, Wtemp2, WCand1, WCand2;
  Particle llj, l1j, l2j,  lljj, l1jj, l2jj, l1J, l2J;
  Particle WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp, Ztemp1, Ztemp2, Ztemp3, Ztemp4, ZCand1, ZCand2, GammaCand, GammaLep1, GammaLep2;
  int ossf_mass10 = 0;
  
  // Set tight_iso cut & calculate pTcone
  double mu_tight_iso = 0.05;
  if(IDsuffix == "HN16") mu_tight_iso = 0.07;

  double el_tight_iso = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  // Define leptons (pT order)
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));
  std::sort(leptons.begin(), leptons.end(), PtComparingPtr); //JH : lepton sorting

  // Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  lepton_veto_size = leptons_veto.size() - leptons.size();

  // Define ST, MET^2/ST
  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();
  ST += MET;
  MET2ST = MET*MET/ST;

  //=====================================================================================
  //=====================================================================================
  //==== Preselection, low/high mass signal regions
  //=====================================================================================
  //=====================================================================================

  //=========================
  //==== Event selections..
  //=========================

  // Loop for each channel : it_ch (0,1,2) = (mumu, ee, emu)
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){

    if(channels.at(it_ch)=="dimu"){ LeptonPtCut1 = MuonPtCut1; LeptonPtCut2 = MuonPtCut2; }
    else if(channels.at(it_ch)=="diel"){ LeptonPtCut1 = ElectronPtCut1; LeptonPtCut2 = ElectronPtCut2; }
    else if(channels.at(it_ch)=="emu"){ LeptonPtCut1 = EMuPtCut1; LeptonPtCut2 = EMuPtCut2; }
    if( (channels.at(it_ch).Contains("mu")) && RunCF) continue; //JH : mumu, emu are irrelevant to CF

    // Triggers for each channel
    if(channels.at(it_ch)=="dimu" && !ev.PassTrigger(MuonTriggers)) continue; //JH : NOTE PassTrigger runs for loop and returns true even if a single item in triggers vector is fired;
    if(channels.at(it_ch)=="diel" && !ev.PassTrigger(ElectronTriggers)) continue;
    if(channels.at(it_ch)=="emu" && !ev.PassTrigger(EMuTriggers)) continue; //JH : NOTE logically, I can use the following cut to avoid possible trigger double counting : e.g. (1) pass MuonTrigger (2) pass Electron Trigger && not pass Muon Trigger (3) pass EMu Trigger && not pass Muon Trigger nor Electron Trigger. But then this would not be consistent with the MC trigger lumi weight. So I will just let it as is. the possible (small) double counting happen in data and MC both, so doesn't matter.

    // Period-dependent trigger weight (only for 2016 MC)
    trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
    if(!IsDATA){
      if(DataYear==2016){
        if(channels.at(it_ch)=="dimu"){
          if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
          if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628; //JH : ??? in AN, they used dz with all period. why should we allow (dimu & not dz)?
          trigger_lumi = dimu_trig_weight;
        }
        if(channels.at(it_ch)=="diel") trigger_lumi = ev.GetTriggerLumi("Full");
        if(channels.at(it_ch)=="emu"){
          if(ev.PassTrigger(EMuTriggers)) emu_trig_weight += 27267.591;
          if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628; //JH
          trigger_lumi = emu_trig_weight; 
        }
      }
      else{
        trigger_lumi = ev.GetTriggerLumi("Full");
      }
    }

    weight = 1.;
    if(!IsDATA){
      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);
    } //JH : recalculate total weight for 2016 period dependency.

    // Cutflow : passing dilepton triggers
    FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

    if(leptons.size() == 2){ 
      if(channels.at(it_ch)=="dimu"){ if(!(muons.size()==2 && electrons.size()==0)) continue; }
      if(channels.at(it_ch)=="diel"){ if(!(muons.size()==0 && electrons.size()==2)) continue; }
      if(channels.at(it_ch)=="emu"){ if(!(muons.size()==1 && electrons.size()==1)) continue; }

      ZCand = *leptons.at(0) + *leptons.at(1);

      weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
      // weights for MC
      if(!IsDATA){
        
        // Select prompt only
        if(-4<=GetLeptonType(*leptons.at(0), gens)&&GetLeptonType(*leptons.at(0), gens)<=0) continue;
        if(-4<=GetLeptonType(*leptons.at(1), gens)&&GetLeptonType(*leptons.at(1), gens)<=0) continue;

        weight *= weight_norm_1invpb*trigger_lumi; //JH : trigger_lumi for period dependency
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
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
          weight *= muon_recosf*muon_idsf*muon_isosf; //JH : FIXME no muon id, iso SF applied
        }
        for(unsigned int j=0; j<electrons.size(); j++){
//          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0); //JH : FIXME no electron ID SF applied
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf; //JH : recalculate total weight applying 2016 trigger lumi dependency && lepton SF
        }
      } //JH : Now total weight including 2016 trigger lumi and lepton SF done && lepton gen-matching (for first 2 leptons only) done

      /////////////////////////////////////////////////////////
      //// Preselection (triggers have been already applied.)
      /////////////////////////////////////////////////////////

      if(!(leptons.at(0)->Pt()>LeptonPtCut1 && leptons.at(1)->Pt()>LeptonPtCut2)) continue;
      
      // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      if(!RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) continue;
      if(RunCF && leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) continue;

      // Cutflow : same-sign (oppsite-sign when RunCF=true)
      FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto 3rd leptons using veto ID
      FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(ZCand.M() > 10.)) continue; 
      if(channels.at(it_ch)=="diel" && IsOnZ(ZCand.M(), 10.)) continue; //JH : see p.12 of preapproval -> https://indico.cern.ch/event/694943/contributions/2849972/attachments/1583026/2501796/180115__JaesungKim__JetsX_Meeting__HN_DiLepton_PreApproval.pdf

      // Cutflow : m(ll) > 10 GeV, |m(ll)-m(Z)| > 10 GeV for ee 
      FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

      //if(!( fatjets.size()>0 || (jets.size()>1 && fatjets.size()==0) || (jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.) )) continue;
      if(!( fatjets.size()>0 || jets.size()>0 )) continue; //JH : new preselection?
      
      FillHist(channels.at(it_ch)+"/Pre/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      if(fatjets.size()>0) FillHist(channels.at(it_ch)+"/Pre/FatJet_Mass_"+IDsuffix, fatjets.at(0).M(), weight, 1000, 0., 1000.);
      if(jets.size()>0){
        FillHist(channels.at(it_ch)+"/Pre/Jet1_Mass_"+IDsuffix, jets.at(0).M(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/Pre/Jet1_Pt_"+IDsuffix, jets.at(0).Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/Pre/Jet1_Eta_"+IDsuffix, jets.at(0).Eta(), weight, 54, -2.7, 2.7);
      }
      if(jets.size()>1){
        FillHist(channels.at(it_ch)+"/Pre/DiJet_Mass_"+IDsuffix, (jets.at(0)+jets.at(1)).M(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/Pre/Jet2_Pt_"+IDsuffix, jets.at(1).Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/Pre/Jet2_Eta_"+IDsuffix, jets.at(1).Eta(), weight, 54, -2.7, 2.7);
      }
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      // Cutflow : jet requirement (This is the number or events at preselection)
      FillHist(channels.at(it_ch)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

    } //JH : if lepton size 2.
  } //JH : for loop in ee, mm, em channel 
}

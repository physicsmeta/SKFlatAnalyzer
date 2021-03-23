#include "test.h"

test::test(){

}

void test::initializeAnalyzer(){

  RunSyst = HasFlag("RunSyst");
  RunFake = HasFlag("RunFake");
  RunCF = HasFlag("RunCF");
  RunOS   = HasFlag("RunOS");

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
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/test.h 
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

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps); //JH : NOTE This is used in mcCorr->SetupJetTagging() in m.initializeAnalyzerTools();
}

test::~test(){

  //==== Destructor of this Analyzer

}

void test::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/test.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/test.h
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
//    param.Name = "Central";

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
    param.Electron_ID_SF_Key = "";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone = true;

    // Jet ID
    param.Jet_ID = "HNTight";
    if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45"; //JH : TODO

    executeEventFromParameter(param);

  }
}

void test::executeEventFromParameter(AnalyzerParameter param){

  vector<TString> channels = {"dimu", "diel", "emu"};
  vector<TString> regions = {"fakeCR1", "lowSR1", "lowCR1", "highSR1", "highCR1", "lowSR2", "lowCR2", "highSR2", "highCR2"}; 
  TString IDsuffix = "Run2";
  double cutflow_max = 18.;
  int cutflow_bin = 18;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;
 
  Event ev = GetEvent();
  vector<Gen> gens_signcheck = GetGens(); //JH

  //=============
  //==== No Cut
  //=============

  vector<Muon> truth_muons;
  vector<Electron> truth_eles;
  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full"); //JH : weight_norm_1invpb = xsec/sumW; Lumi = 35.9, 41.5, 59.7(fb-1) total 137fb-1
    weight *= ev.MCweight(); //JH : gen_weight in MiniAOD
    weight *= GetPrefireWeight(0); //JH : No issue in 2018, otherwise returns L1PrefireReweight_Central in MiniAOD
    weight *= GetPileUpWeight(nPileUp,0); //JH : mcCorr->GetPileUpWeight(N_pileup, syst); mcCorr->GetPileUpWeight2017(N_pileup, syst); NOTE 2018 not yet added.

    int IsPP;
    //Count total number of pp / mm
    for(int i=0;i<gens_signcheck.size();i++){
      if(gens_signcheck.at(i).PID()==13&&gens_signcheck.at(i).fromHardProcessFinalState()){
        FillHist("SampleCharge_unweighted",-1,1,3,-1,2);
        FillHist("SampleCharge",-1,weight,3,-1,2);
        IsPP = 0;
        //JHPrintGen(gens_signcheck);
        break;
      }
      if(gens_signcheck.at(i).PID()==-13&&gens_signcheck.at(i).fromHardProcessFinalState()){
        FillHist("SampleCharge_unweighted",1,1,3,-1,2);
        FillHist("SampleCharge",1,weight,3,-1,2);
        IsPP = 1;
        //JHPrintGen(gens_signcheck);
        break;
      }
    }

    //if(HasFlag("pp")){
    //  if(IsPP == 0) return;
    //}
    //if(HasFlag("mm")){
    //  if(IsPP == 1) return;
    //}

    for(int i=0; i<AllMuons.size(); i++){
      Gen truth_mu = GetGenMatchedLepton(AllMuons.at(i), gens_signcheck);
      if(truth_mu.fromHardProcessFinalState()) truth_muons.push_back(AllMuons.at(i));
    }
    for(int i=0; i<AllElectrons.size(); i++){
      Gen truth_el = GetGenMatchedLepton(AllElectrons.at(i), gens_signcheck);
      if(truth_el.fromHardProcessFinalState()) truth_eles.push_back(AllElectrons.at(i));
    }

    //PrintGen(gens_signcheck);

    FillHist("Nocut/truth_muon/Size", truth_muons.size(), weight, 10, 0., 10.);
    for(int i=0; i<truth_muons.size(); i++){
      FillHist("Nocut/truth_muon/Pt", truth_muons.at(i).Pt(), weight, 1000, 0., 1000.); //JH : rc pt
      FillHist("Nocut/truth_muon/Eta", truth_muons.at(i).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Nocut/truth_muon/Phi", truth_muons.at(i).Phi(), weight, 63, -3.15, 3.15);
      FillHist("Nocut/truth_muon/Energy", truth_muons.at(i).Energy(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_muon/Mass", truth_muons.at(i).M(), weight, 1000, 0., 10.); //JH : from muon_mass in SKFlatNtuple.h
      FillHist("Nocut/truth_muon/dXY", truth_muons.at(i).dXY(), weight, 1000, -0.1, 0.1);
      FillHist("Nocut/truth_muon/dXYSig", fabs(truth_muons.at(i).dXY())/truth_muons.at(i).dXYerr(), weight, 100, 0., 10.); 
      FillHist("Nocut/truth_muon/dZ", truth_muons.at(i).dZ(), weight, 1000, -0.5, 0.5);
      FillHist("Nocut/truth_muon/dZSig", fabs(truth_muons.at(i).dZ())/truth_muons.at(i).dZerr(), weight, 100, 0., 10.); 
      FillHist("Nocut/truth_muon/IP3D", truth_muons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
      FillHist("Nocut/truth_muon/IP3DSig", fabs(truth_muons.at(i).IP3D())/truth_muons.at(i).IP3Derr(), weight, 100, 0., 10.);
      FillHist("Nocut/truth_muon/RelIso", truth_muons.at(i).RelIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_muon/TrkIso", truth_muons.at(i).TrkIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_muon/ValidMuonHits", truth_muons.at(i).ValidMuonHits(), weight, 10, 0., 10.);
      FillHist("Nocut/truth_muon/MatchedStations", truth_muons.at(i).MatchedStations(), weight, 10, 0., 10.);
      FillHist("Nocut/truth_muon/PixelHits", truth_muons.at(i).PixelHits(), weight, 10, 0., 10.);
      FillHist("Nocut/truth_muon/TrackerLayers", truth_muons.at(i).TrackerLayers(), weight, 10, 0., 10.);
      FillHist("Nocut/truth_muon/Chi2", truth_muons.at(i).Chi2(), weight, 500, 0., 50.);
    }
    FillHist("Nocut/truth_ele/Size", truth_eles.size(), weight, 10, 0., 10.);
    for(int i=0; i<truth_eles.size(); i++){
      FillHist("Nocut/truth_ele/Pt", truth_eles.at(i).Pt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/Eta", truth_eles.at(i).Eta(), weight, 50, -2.5, 2.5);
      FillHist("Nocut/truth_ele/Phi", truth_eles.at(i).Phi(), weight, 63, -3.15, 3.15);
      FillHist("Nocut/truth_ele/Energy", truth_eles.at(i).Energy(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/Mass", truth_eles.at(i).M(), weight, 1000, 0., 10.); 
      FillHist("Nocut/truth_ele/UncorrPt", truth_eles.at(i).UncorrPt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/UncorrE", truth_eles.at(i).UncorrE(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/scE", truth_eles.at(i).scE(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/scEta", truth_eles.at(i).scEta(), weight, 50, -2.5, 2.5);
      FillHist("Nocut/truth_ele/scPhi", truth_eles.at(i).scPhi(), weight, 63, -3.15, 3.15);
      FillHist("Nocut/truth_ele/dXY", truth_eles.at(i).dXY(), weight, 1000, -0.1, 0.1);
      FillHist("Nocut/truth_ele/dXYSig", fabs(truth_eles.at(i).dXY())/truth_eles.at(i).dXYerr(), weight, 100, 0., 10.); 
      FillHist("Nocut/truth_ele/dZ", truth_eles.at(i).dZ(), weight, 1000, -0.5, 0.5);
      FillHist("Nocut/truth_ele/dZSig", fabs(truth_eles.at(i).dZ())/truth_eles.at(i).dZerr(), weight, 100, 0., 10.); 
      FillHist("Nocut/truth_ele/IP3D", truth_eles.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
      FillHist("Nocut/truth_ele/IP3DSig", fabs(truth_eles.at(i).IP3D())/truth_eles.at(i).IP3Derr(), weight, 100, 0., 10.);
      FillHist("Nocut/truth_ele/RelIso", truth_eles.at(i).RelIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/TrkIso", truth_eles.at(i).TrkIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/PassConversionVeto", truth_eles.at(i).PassConversionVeto(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/NMissingHits", truth_eles.at(i).NMissingHits(), weight, 10, 0., 10.);
      FillHist("Nocut/truth_ele/MVAIso", truth_eles.at(i).MVAIso(), weight, 200, -1., 1.);
      FillHist("Nocut/truth_ele/MVANoIso", truth_eles.at(i).MVANoIso(), weight, 200, -1., 1.);
      FillHist("Nocut/truth_ele/IsGsfCtfScPixChargeConsistent", truth_eles.at(i).IsGsfCtfScPixChargeConsistent(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/Full5x5_sigmaIetaIeta", truth_eles.at(i).Full5x5_sigmaIetaIeta(), weight, 10, 0., 0.1);
      FillHist("Nocut/truth_ele/dEtaSeed", truth_eles.at(i).dEtaSeed(), weight, 20, -0.1, 0.1);
      FillHist("Nocut/truth_ele/dPhiIn", truth_eles.at(i).dPhiIn(), weight, 30, -0.15, 0.15);
      FillHist("Nocut/truth_ele/HoverE", truth_eles.at(i).HoverE(), weight, 50, 0., 0.5);
      FillHist("Nocut/truth_ele/InvEminusInvP", truth_eles.at(i).InvEminusInvP(), weight, 200, -1., 1.);
      FillHist("Nocut/truth_ele/e2x5OverE5x5", truth_eles.at(i).e2x5OverE5x5(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/e1x5OverE5x5", truth_eles.at(i).e1x5OverE5x5(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/dr03EcalRecHitSumEt", truth_eles.at(i).dr03EcalRecHitSumEt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/dr03HcalDepth1TowerSumEt", truth_eles.at(i).dr03HcalDepth1TowerSumEt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/dr03HcalTowerSumEt", truth_eles.at(i).dr03HcalTowerSumEt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/dr03TkSumPt", truth_eles.at(i).dr03TkSumPt(), weight, 1000, 0., 1000.);
      FillHist("Nocut/truth_ele/ecalPFClusterIso", truth_eles.at(i).ecalPFClusterIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/hcalPFClusterIso", truth_eles.at(i).hcalPFClusterIso(), weight, 100, 0., 1.);
      FillHist("Nocut/truth_ele/passVetoID", truth_eles.at(i).passVetoID(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passLooseID", truth_eles.at(i).passLooseID(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passMediumID", truth_eles.at(i).passMediumID(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passTightID", truth_eles.at(i).passTightID(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passMVAID_noIso_WP80", truth_eles.at(i).passMVAID_noIso_WP80(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passMVAID_noIso_WP90", truth_eles.at(i).passMVAID_noIso_WP90(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passMVAID_iso_WP80", truth_eles.at(i).passMVAID_iso_WP80(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/passMVAID_iso_WP90", truth_eles.at(i).passMVAID_iso_WP90(), weight, 2, 0., 2.);
      FillHist("Nocut/truth_ele/Pass_HNTight2016", truth_eles.at(i).Pass_HNTight2016(), weight, 2, 0., 2.);
    }

  } //JH : total weight calculation done.


  // Cutflow : No Cuts
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }

  FillHist("Nocut/AllMuon/Size", AllMuons.size(), weight, 100, 0., 100.);
  for(int i=0; i<AllMuons.size(); i++){
    FillHist("Nocut/AllMuon/Pt", AllMuons.at(i).Pt(), weight, 1000, 0., 1000.); //JH : rc pt
    FillHist("Nocut/AllMuon/Eta", AllMuons.at(i).Eta(), weight, 50, -2.5, 2.5);
    FillHist("Nocut/AllMuon/Phi", AllMuons.at(i).Phi(), weight, 63, -3.15, 3.15);
    FillHist("Nocut/AllMuon/Energy", AllMuons.at(i).Energy(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllMuon/Mass", AllMuons.at(i).M(), weight, 1000, 0., 10.); //JH : from muon_mass in SKFlatNtuple.h
    FillHist("Nocut/AllMuon/dXY", AllMuons.at(i).dXY(), weight, 1000, -0.1, 0.1);
    FillHist("Nocut/AllMuon/dXYSig", fabs(AllMuons.at(i).dXY())/AllMuons.at(i).dXYerr(), weight, 100, 0., 10.); 
    FillHist("Nocut/AllMuon/dZ", AllMuons.at(i).dZ(), weight, 1000, -0.5, 0.5);
    FillHist("Nocut/AllMuon/dZSig", fabs(AllMuons.at(i).dZ())/AllMuons.at(i).dZerr(), weight, 100, 0., 10.); 
    FillHist("Nocut/AllMuon/IP3D", AllMuons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
    FillHist("Nocut/AllMuon/IP3DSig", fabs(AllMuons.at(i).IP3D())/AllMuons.at(i).IP3Derr(), weight, 100, 0., 10.);
    FillHist("Nocut/AllMuon/RelIso", AllMuons.at(i).RelIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllMuon/TrkIso", AllMuons.at(i).TrkIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllMuon/ValidMuonHits", AllMuons.at(i).ValidMuonHits(), weight, 10, 0., 10.);
    FillHist("Nocut/AllMuon/MatchedStations", AllMuons.at(i).MatchedStations(), weight, 10, 0., 10.);
    FillHist("Nocut/AllMuon/PixelHits", AllMuons.at(i).PixelHits(), weight, 10, 0., 10.);
    FillHist("Nocut/AllMuon/TrackerLayers", AllMuons.at(i).TrackerLayers(), weight, 10, 0., 10.);
    FillHist("Nocut/AllMuon/Chi2", AllMuons.at(i).Chi2(), weight, 500, 0., 50.);
  }
  FillHist("Nocut/AllElectron/Size", AllElectrons.size(), weight, 100, 0., 100.);
  for(int i=0; i<AllElectrons.size(); i++){
    FillHist("Nocut/AllElectron/Pt", AllElectrons.at(i).Pt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/Eta", AllElectrons.at(i).Eta(), weight, 50, -2.5, 2.5);
    FillHist("Nocut/AllElectron/Phi", AllElectrons.at(i).Phi(), weight, 63, -3.15, 3.15);
    FillHist("Nocut/AllElectron/Energy", AllElectrons.at(i).Energy(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/Mass", AllElectrons.at(i).M(), weight, 1000, 0., 10.); 
    FillHist("Nocut/AllElectron/UncorrPt", AllElectrons.at(i).UncorrPt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/UncorrE", AllElectrons.at(i).UncorrE(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/scE", AllElectrons.at(i).scE(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/scEta", AllElectrons.at(i).scEta(), weight, 50, -2.5, 2.5);
    FillHist("Nocut/AllElectron/scPhi", AllElectrons.at(i).scPhi(), weight, 63, -3.15, 3.15);
    FillHist("Nocut/AllElectron/dXY", AllElectrons.at(i).dXY(), weight, 1000, -0.1, 0.1);
    FillHist("Nocut/AllElectron/dXYSig", fabs(AllElectrons.at(i).dXY())/AllElectrons.at(i).dXYerr(), weight, 100, 0., 10.); 
    FillHist("Nocut/AllElectron/dZ", AllElectrons.at(i).dZ(), weight, 1000, -0.5, 0.5);
    FillHist("Nocut/AllElectron/dZSig", fabs(AllElectrons.at(i).dZ())/AllElectrons.at(i).dZerr(), weight, 100, 0., 10.); 
    FillHist("Nocut/AllElectron/IP3D", AllElectrons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
    FillHist("Nocut/AllElectron/IP3DSig", fabs(AllElectrons.at(i).IP3D())/AllElectrons.at(i).IP3Derr(), weight, 100, 0., 10.);
    FillHist("Nocut/AllElectron/RelIso", AllElectrons.at(i).RelIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/TrkIso", AllElectrons.at(i).TrkIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/PassConversionVeto", AllElectrons.at(i).PassConversionVeto(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/NMissingHits", AllElectrons.at(i).NMissingHits(), weight, 10, 0., 10.);
    FillHist("Nocut/AllElectron/MVAIso", AllElectrons.at(i).MVAIso(), weight, 200, -1., 1.);
    FillHist("Nocut/AllElectron/MVANoIso", AllElectrons.at(i).MVANoIso(), weight, 200, -1., 1.);
    FillHist("Nocut/AllElectron/IsGsfCtfScPixChargeConsistent", AllElectrons.at(i).IsGsfCtfScPixChargeConsistent(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/Full5x5_sigmaIetaIeta", AllElectrons.at(i).Full5x5_sigmaIetaIeta(), weight, 10, 0., 0.1);
    FillHist("Nocut/AllElectron/dEtaSeed", AllElectrons.at(i).dEtaSeed(), weight, 20, -0.1, 0.1);
    FillHist("Nocut/AllElectron/dPhiIn", AllElectrons.at(i).dPhiIn(), weight, 30, -0.15, 0.15);
    FillHist("Nocut/AllElectron/HoverE", AllElectrons.at(i).HoverE(), weight, 50, 0., 0.5);
    FillHist("Nocut/AllElectron/InvEminusInvP", AllElectrons.at(i).InvEminusInvP(), weight, 200, -1., 1.);
    FillHist("Nocut/AllElectron/e2x5OverE5x5", AllElectrons.at(i).e2x5OverE5x5(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/e1x5OverE5x5", AllElectrons.at(i).e1x5OverE5x5(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/dr03EcalRecHitSumEt", AllElectrons.at(i).dr03EcalRecHitSumEt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/dr03HcalDepth1TowerSumEt", AllElectrons.at(i).dr03HcalDepth1TowerSumEt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/dr03HcalTowerSumEt", AllElectrons.at(i).dr03HcalTowerSumEt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/dr03TkSumPt", AllElectrons.at(i).dr03TkSumPt(), weight, 1000, 0., 1000.);
    FillHist("Nocut/AllElectron/ecalPFClusterIso", AllElectrons.at(i).ecalPFClusterIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/hcalPFClusterIso", AllElectrons.at(i).hcalPFClusterIso(), weight, 100, 0., 1.);
    FillHist("Nocut/AllElectron/passVetoID", AllElectrons.at(i).passVetoID(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passLooseID", AllElectrons.at(i).passLooseID(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passMediumID", AllElectrons.at(i).passMediumID(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passTightID", AllElectrons.at(i).passTightID(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passMVAID_noIso_WP80", AllElectrons.at(i).passMVAID_noIso_WP80(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passMVAID_noIso_WP90", AllElectrons.at(i).passMVAID_noIso_WP90(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passMVAID_iso_WP80", AllElectrons.at(i).passMVAID_iso_WP80(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/passMVAID_iso_WP90", AllElectrons.at(i).passMVAID_iso_WP90(), weight, 2, 0., 2.);
    FillHist("Nocut/AllElectron/Pass_HNTight2016", AllElectrons.at(i).Pass_HNTight2016(), weight, 2, 0., 2.);
  }

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  // Cutflow : MET filter
  for(unsigned int it_ch=0; it_ch<channels.size(); it_ch++){
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
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
    cout << "[test::executeEventFromParameter] Wrong syst" << endl;
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
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5); //JH : lepton selection done
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7); //JH : to reject bjets

  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  vector<Jet> jets;
  vector<FatJet> fatjets;
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
      if(this_AllFatJets.at(i).DeltaR(muons_veto.at(j)) < 1.0) lepton_count1++; //JH : muon cleaning
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllFatJets.at(i).DeltaR(electrons_veto.at(j)) < 1.0) lepton_count1++; //JH : electron cleaning
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
      if(this_AllJets.at(i).DeltaR(fatjets.at(j)) < 0.8) fatjet_count++; //JH : fatjet cleaning
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

  double ST = 0.;
  double MET2ST = 0.;
  double MW = 80.379;
  double muon_recosf = 1., muon_idsf = 1., muon_isosf = 1.;
  double ele_idsf = 1., ele_recosf = 1.;
  double fatjetTau21SF = 1.;
  int lepton_veto_size = 0;
  double LeptonPtCut1 = 0., LeptonPtCut2 = 0.;
  Particle ZCand, Wtemp1, Wtemp2, WCand1, WCand2, WCand3;
  Particle llj, l1j, l2j,  lljj, l1jj, l2jj, l1J, l2J;
  Particle WtagLep, TriLep, ZtagLep1, ZtagLep2, Ztemp, Ztemp1, Ztemp2, Ztemp3, Ztemp4, ZCand1, ZCand2, GammaCand, GammaLep1, GammaLep2;
  
  // Set tight_iso cut & calculate pTcone
  double mu_tight_iso = 0.07, el_tight_iso = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(RunFake){
    for(unsigned int i=0; i<muons.size(); i++){
      this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso);
      muons.at(i).SetPtCone(this_ptcone_muon);
    }

    for(unsigned int i=0; i<electrons.size(); i++){

      if(param.Electron_Tight_ID.Contains("HNTight")){ // POG cut-based tight WP
        el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt();
        if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
      }

      //if(param.Electron_Tight_ID.Contains("HNMVA")) el_tight_iso = 0.08;

      this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
      electrons.at(i).SetPtCone(this_ptcone_electron);

    }
  }

  // Define leptons (pT order)
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));
  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  // Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  lepton_veto_size = leptons_veto.size() - leptons.size();

/*
  // leptons (minus, plus charge)
  for(unsigned int i=0; i<muons.size(); i++){
    if(muons.at(i).Charge() < 0) leptons_minus.push_back(&muons.at(i));
    if(muons.at(i).Charge() > 0) leptons_plus.push_back(&muons.at(i));
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    if(electrons.at(i).Charge() < 0) leptons_minus.push_back(&electrons.at(i));
    if(electrons.at(i).Charge() > 0) leptons_plus.push_back(&electrons.at(i));
  }
*/

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

    if(it_ch==0){ LeptonPtCut1 = MuonPtCut1; LeptonPtCut2 = MuonPtCut2; }
    if(it_ch==1){ LeptonPtCut1 = ElectronPtCut1; LeptonPtCut2 = ElectronPtCut2; }
    if(it_ch==2){ LeptonPtCut1 = EMuPtCut1; LeptonPtCut2 = EMuPtCut2; }
    if((it_ch==0||it_ch==2) && RunCF) continue; //JH : mumu, emu are irrelevant to CF

    // Triggers for each channel
    if(it_ch==0 && !ev.PassTrigger(MuonTriggers)) continue; //JH : NOTE PassTrigger runs for loop and returns true even if a single item in triggers vector is fired;
    if(it_ch==1 && !ev.PassTrigger(ElectronTriggers)) continue;
    if(it_ch==2 && !ev.PassTrigger(EMuTriggers)) continue; //JH : NOTE logically, I can use the following cut to avoid possible trigger double counting : e.g. (1) pass MuonTrigger (2) pass Electron Trigger && not pass Muon Trigger (3) pass EMu Trigger && not pass Muon Trigger nor Electron Trigger. But then this would not be consistent with the MC trigger lumi weight. So I will just let it as is. the possible (small) double counting happen in data and MC both, so doesn't matter.


    // Period-dependent trigger weight (only for 2016 MC)
    trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
    if(!IsDATA){
      if(DataYear==2016){
        if(it_ch==0){
          if(ev.PassTrigger(MuonTriggers)) dimu_trig_weight += 27267.591;
          if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628; //JH
          trigger_lumi = dimu_trig_weight;
        }
        if(it_ch==1) trigger_lumi = ev.GetTriggerLumi("Full");
        if(it_ch==2){
          if(ev.PassTrigger(EMuTriggers)) emu_trig_weight += 27267.591;
          if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628; //JH
          trigger_lumi = emu_trig_weight; 
        }
      }
      else{
        trigger_lumi = ev.GetTriggerLumi("Full");
      }
    }

    // Cutflow : passing dilepton triggers
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      weight = 1.;
      if(!IsDATA){
        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);
      } //JH : recalculate total weight for 2016 period dependency.
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      /*
      if(!IsDATA){
        FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Size", truth_muons.size(), weight, 10, 0., 10.);
        for(int i=0; i<truth_muons.size(); i++){
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Pt", truth_muons.at(i).Pt(), weight, 1000, 0., 1000.); //JH : rc pt
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Eta", truth_muons.at(i).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Phi", truth_muons.at(i).Phi(), weight, 63, -3.15, 3.15);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Energy", truth_muons.at(i).Energy(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Mass", truth_muons.at(i).M(), weight, 1000, 0., 10.); //JH : from muon_mass in SKFlatNtuple.h
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/dXY", truth_muons.at(i).dXY(), weight, 1000, -0.1, 0.1);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/dXYSig", fabs(truth_muons.at(i).dXY())/truth_muons.at(i).dXYerr(), weight, 100, 0., 10.); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/dZ", truth_muons.at(i).dZ(), weight, 1000, -0.5, 0.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/dZSig", fabs(truth_muons.at(i).dZ())/truth_muons.at(i).dZerr(), weight, 100, 0., 10.); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/IP3D", truth_muons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/IP3DSig", fabs(truth_muons.at(i).IP3D())/truth_muons.at(i).IP3Derr(), weight, 100, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/RelIso", truth_muons.at(i).RelIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/TrkIso", truth_muons.at(i).TrkIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/ValidMuonHits", truth_muons.at(i).ValidMuonHits(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/MatchedStations", truth_muons.at(i).MatchedStations(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/PixelHits", truth_muons.at(i).PixelHits(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/TrackerLayers", truth_muons.at(i).TrackerLayers(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_muon/Chi2", truth_muons.at(i).Chi2(), weight, 500, 0., 50.);
        }
        FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Size", truth_eles.size(), weight, 10, 0., 10.);
        for(int i=0; i<truth_eles.size(); i++){
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Pt", truth_eles.at(i).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Eta", truth_eles.at(i).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Phi", truth_eles.at(i).Phi(), weight, 63, -3.15, 3.15);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Energy", truth_eles.at(i).Energy(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Mass", truth_eles.at(i).M(), weight, 1000, 0., 10.); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/UncorrPt", truth_eles.at(i).UncorrPt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/UncorrE", truth_eles.at(i).UncorrE(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/scE", truth_eles.at(i).scE(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/scEta", truth_eles.at(i).scEta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/scPhi", truth_eles.at(i).scPhi(), weight, 63, -3.15, 3.15);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dXY", truth_eles.at(i).dXY(), weight, 1000, -0.1, 0.1);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dXYSig", fabs(truth_eles.at(i).dXY())/truth_eles.at(i).dXYerr(), weight, 100, 0., 10.); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dZ", truth_eles.at(i).dZ(), weight, 1000, -0.5, 0.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dZSig", fabs(truth_eles.at(i).dZ())/truth_eles.at(i).dZerr(), weight, 100, 0., 10.); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/IP3D", truth_eles.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/IP3DSig", fabs(truth_eles.at(i).IP3D())/truth_eles.at(i).IP3Derr(), weight, 100, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/RelIso", truth_eles.at(i).RelIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/TrkIso", truth_eles.at(i).TrkIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/PassConversionVeto", truth_eles.at(i).PassConversionVeto(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/NMissingHits", truth_eles.at(i).NMissingHits(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/MVAIso", truth_eles.at(i).MVAIso(), weight, 200, -1., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/MVANoIso", truth_eles.at(i).MVANoIso(), weight, 200, -1., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/IsGsfCtfScPixChargeConsistent", truth_eles.at(i).IsGsfCtfScPixChargeConsistent(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Full5x5_sigmaIetaIeta", truth_eles.at(i).Full5x5_sigmaIetaIeta(), weight, 10, 0., 0.1);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dEtaSeed", truth_eles.at(i).dEtaSeed(), weight, 20, -0.1, 0.1);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dPhiIn", truth_eles.at(i).dPhiIn(), weight, 30, -0.15, 0.15);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/HoverE", truth_eles.at(i).HoverE(), weight, 50, 0., 0.5);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/InvEminusInvP", truth_eles.at(i).InvEminusInvP(), weight, 200, -1., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/e2x5OverE5x5", truth_eles.at(i).e2x5OverE5x5(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/e1x5OverE5x5", truth_eles.at(i).e1x5OverE5x5(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dr03EcalRecHitSumEt", truth_eles.at(i).dr03EcalRecHitSumEt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dr03HcalDepth1TowerSumEt", truth_eles.at(i).dr03HcalDepth1TowerSumEt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dr03HcalTowerSumEt", truth_eles.at(i).dr03HcalTowerSumEt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/dr03TkSumPt", truth_eles.at(i).dr03TkSumPt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/ecalPFClusterIso", truth_eles.at(i).ecalPFClusterIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/hcalPFClusterIso", truth_eles.at(i).hcalPFClusterIso(), weight, 100, 0., 1.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passVetoID", truth_eles.at(i).passVetoID(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passLooseID", truth_eles.at(i).passLooseID(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passMediumID", truth_eles.at(i).passMediumID(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passTightID", truth_eles.at(i).passTightID(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passMVAID_noIso_WP80", truth_eles.at(i).passMVAID_noIso_WP80(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passMVAID_noIso_WP90", truth_eles.at(i).passMVAID_noIso_WP90(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passMVAID_iso_WP80", truth_eles.at(i).passMVAID_iso_WP80(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/passMVAID_iso_WP90", truth_eles.at(i).passMVAID_iso_WP90(), weight, 2, 0., 2.);
          FillHist(channels.at(it_ch)+"/AfterTrig/truth_ele/Pass_HNTight2016", truth_eles.at(i).Pass_HNTight2016(), weight, 2, 0., 2.);
        }
      }
      FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Size", AllMuons.size(), weight, 100, 0., 100.);
      for(int i=0; i<AllMuons.size(); i++){
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Pt", AllMuons.at(i).Pt(), weight, 1000, 0., 1000.); //JH : rc pt
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Eta", AllMuons.at(i).Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Phi", AllMuons.at(i).Phi(), weight, 63, -3.15, 3.15);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Energy", AllMuons.at(i).Energy(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Mass", AllMuons.at(i).M(), weight, 1000, 0., 10.); //JH : from muon_mass in SKFlatNtuple.h
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/dXY", AllMuons.at(i).dXY(), weight, 1000, -0.1, 0.1);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/dXYSig", fabs(AllMuons.at(i).dXY())/AllMuons.at(i).dXYerr(), weight, 100, 0., 10.); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/dZ", AllMuons.at(i).dZ(), weight, 1000, -0.5, 0.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/dZSig", fabs(AllMuons.at(i).dZ())/AllMuons.at(i).dZerr(), weight, 100, 0., 10.); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/IP3D", AllMuons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/IP3DSig", fabs(AllMuons.at(i).IP3D())/AllMuons.at(i).IP3Derr(), weight, 100, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/RelIso", AllMuons.at(i).RelIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/TrkIso", AllMuons.at(i).TrkIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/ValidMuonHits", AllMuons.at(i).ValidMuonHits(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/MatchedStations", AllMuons.at(i).MatchedStations(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/PixelHits", AllMuons.at(i).PixelHits(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/TrackerLayers", AllMuons.at(i).TrackerLayers(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllMuon/Chi2", AllMuons.at(i).Chi2(), weight, 500, 0., 50.);
      }
      FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Size", AllElectrons.size(), weight, 100, 0., 100.);
      for(int i=0; i<AllElectrons.size(); i++){
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Pt", AllElectrons.at(i).Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Eta", AllElectrons.at(i).Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Phi", AllElectrons.at(i).Phi(), weight, 63, -3.15, 3.15);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Energy", AllElectrons.at(i).Energy(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Mass", AllElectrons.at(i).M(), weight, 1000, 0., 10.); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/UncorrPt", AllElectrons.at(i).UncorrPt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/UncorrE", AllElectrons.at(i).UncorrE(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/scE", AllElectrons.at(i).scE(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/scEta", AllElectrons.at(i).scEta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/scPhi", AllElectrons.at(i).scPhi(), weight, 63, -3.15, 3.15);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dXY", AllElectrons.at(i).dXY(), weight, 1000, -0.1, 0.1);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dXYSig", fabs(AllElectrons.at(i).dXY())/AllElectrons.at(i).dXYerr(), weight, 100, 0., 10.); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dZ", AllElectrons.at(i).dZ(), weight, 1000, -0.5, 0.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dZSig", fabs(AllElectrons.at(i).dZ())/AllElectrons.at(i).dZerr(), weight, 100, 0., 10.); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/IP3D", AllElectrons.at(i).IP3D(), weight, 1000, -0.5, 0.5); 
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/IP3DSig", fabs(AllElectrons.at(i).IP3D())/AllElectrons.at(i).IP3Derr(), weight, 100, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/RelIso", AllElectrons.at(i).RelIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/TrkIso", AllElectrons.at(i).TrkIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/PassConversionVeto", AllElectrons.at(i).PassConversionVeto(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/NMissingHits", AllElectrons.at(i).NMissingHits(), weight, 10, 0., 10.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/MVAIso", AllElectrons.at(i).MVAIso(), weight, 200, -1., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/MVANoIso", AllElectrons.at(i).MVANoIso(), weight, 200, -1., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/IsGsfCtfScPixChargeConsistent", AllElectrons.at(i).IsGsfCtfScPixChargeConsistent(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Full5x5_sigmaIetaIeta", AllElectrons.at(i).Full5x5_sigmaIetaIeta(), weight, 10, 0., 0.1);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dEtaSeed", AllElectrons.at(i).dEtaSeed(), weight, 20, -0.1, 0.1);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dPhiIn", AllElectrons.at(i).dPhiIn(), weight, 30, -0.15, 0.15);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/HoverE", AllElectrons.at(i).HoverE(), weight, 50, 0., 0.5);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/InvEminusInvP", AllElectrons.at(i).InvEminusInvP(), weight, 200, -1., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/e2x5OverE5x5", AllElectrons.at(i).e2x5OverE5x5(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/e1x5OverE5x5", AllElectrons.at(i).e1x5OverE5x5(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dr03EcalRecHitSumEt", AllElectrons.at(i).dr03EcalRecHitSumEt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dr03HcalDepth1TowerSumEt", AllElectrons.at(i).dr03HcalDepth1TowerSumEt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dr03HcalTowerSumEt", AllElectrons.at(i).dr03HcalTowerSumEt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/dr03TkSumPt", AllElectrons.at(i).dr03TkSumPt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/ecalPFClusterIso", AllElectrons.at(i).ecalPFClusterIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/hcalPFClusterIso", AllElectrons.at(i).hcalPFClusterIso(), weight, 100, 0., 1.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passVetoID", AllElectrons.at(i).passVetoID(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passLooseID", AllElectrons.at(i).passLooseID(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passMediumID", AllElectrons.at(i).passMediumID(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passTightID", AllElectrons.at(i).passTightID(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passMVAID_noIso_WP80", AllElectrons.at(i).passMVAID_noIso_WP80(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passMVAID_noIso_WP90", AllElectrons.at(i).passMVAID_noIso_WP90(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passMVAID_iso_WP80", AllElectrons.at(i).passMVAID_iso_WP80(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/passMVAID_iso_WP90", AllElectrons.at(i).passMVAID_iso_WP90(), weight, 2, 0., 2.);
        FillHist(channels.at(it_ch)+"/AfterTrig/AllElectron/Pass_HNTight2016", AllElectrons.at(i).Pass_HNTight2016(), weight, 2, 0., 2.);
      }
      */
    }
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max); 

    if(leptons.size() == 2){ 
      if(it_ch==0){ if(!(muons.size()==2 && electrons.size()==0)) continue; }
      if(it_ch==1){ if(!(muons.size()==0 && electrons.size()==2)) continue; }
      if(it_ch==2){ if(!(muons.size()==1 && electrons.size()==1)) continue; } //JH : require 2 tight leptons
      if(leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) FillHist(channels.at(it_ch)+"/IsSameSign_"+IDsuffix,0,1,4,0,4);
      if(leptons.at(0)->Charge()*leptons.at(1)->Charge()>0){
        FillHist(channels.at(it_ch)+"/IsSameSign_"+IDsuffix,1,1,4,0,4);
        if(leptons.at(0)->Charge()>0) FillHist(channels.at(it_ch)+"/IsSameSign_"+IDsuffix,2,1,4,0,4);
        else if(leptons.at(0)->Charge()<0) FillHist(channels.at(it_ch)+"/IsSameSign_"+IDsuffix,3,1,4,0,4); //JH : Check leading lepton's sign. (probably +)
      }

      ZCand = *leptons.at(0) + *leptons.at(1);

      weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1., fatjetTau21SF = 1.;
      // weights for MC
      if(!IsDATA){
        //// Gen matching with dR < 0.1
        //Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        //Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        //if(truth_lep1.PID() == 0) continue;
        //if(truth_lep2.PID() == 0) continue;
        
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
          else if(param.Muon_Tight_ID.Contains("HNTight")){
            muon_recosf = 1.;
            muon_idsf   = mcCorr->MuonID_SF_HNtypeI(param.Muon_Tight_ID, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0); //JH : id+iso once?
            muon_isosf  = 1.;
            if(RunFake){  // When subtracting prompt contribution from fake contribution, we apply ID SF only for muons passing the tight ID
              if(!muons.at(i).PassID(param.Muon_Tight_ID)) muon_idsf = 1.;
            }
          }
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int j=0; j<electrons.size(); j++){
//          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          }
          else if(param.Electron_Tight_ID.Contains("HNTight")){
            ele_idsf = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
            if(RunFake){  // When subtracting prompt contribution from fake contribution, we apply ID SF only for electrons passing the tight ID
              if(!electrons.at(j).PassID(param.Electron_Tight_ID)) ele_idsf = 1.;
            }
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf; //JH : recalculate total weight applying 2016 trigger lumi dependency && lepton SF
        }

        for(unsigned int i=0; i<fatjets.size(); i++){
          fatjetTau21SF = mcCorr->FatJetWtagSF(param.FatJet_ID, 0); //JH : NOTE NOTE NOTE
          weight *= fatjetTau21SF;
        }

      } //JH : Now total weight including 2016 trigger lumi, lepton SF and Wtag SF done && lepton gen-matching (for 2 tight leptons) done

      if(RunFake) weight *= fakeEst->GetWeight(leptons, param);
      if(RunCF) weight *= GetCFweight(leptons, param, true, 0); //JH : this is from the DATA

      /////////////////////////////////////////////////////////
      //// Preselection (triggers have been already applied.)
      /////////////////////////////////////////////////////////

      if(!(leptons.at(0)->Pt()>LeptonPtCut1 && leptons.at(1)->Pt()>LeptonPtCut2)) continue;
      
      // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(RunCF||RunOS) && leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) continue;
      if( (RunCF||RunOS) && leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) continue;

      // Cutflow : same-sign (oppsite-sign when RunCF or RunOS == true) 
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto 3rd leptons using veto ID
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(ZCand.M() > 10.)) continue; 
      if(it_ch==1 && IsOnZ(ZCand.M(), 10.)) continue; //JH : see p.12 of preapproval -> https://indico.cern.ch/event/694943/contributions/2849972/attachments/1583026/2501796/180115__JaesungKim__JetsX_Meeting__HN_DiLepton_PreApproval.pdf

      // Cutflow : m(ll) > 10 GeV, |m(ll)-m(Z)| > 10 GeV for ee 
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      }
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);     

      // non-prompt CR2 : no jets && same-sign back-to-back 2 leptons
      if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){
       
        // Cutflow : jet requirement for non-prompt CR2 
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5)) continue;
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/fakeCR2/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(channels.at(it_ch)+"/fakeCR2/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        
      }

      //if(!( fatjets.size()>0 || (jets.size()>1 && fatjets.size()==0) || (jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.) )) continue; //JH : the last analysis
      if(!( fatjets.size()>0 || jets.size()>0 )) continue; //JH 
     
      FillHist(channels.at(it_ch)+"/Pre/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(channels.at(it_ch)+"/Pre/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.); 
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Pt_HighEnd_"+IDsuffix, leptons.at(0)->Pt(), weight, 10000, 0., 10000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Pt_HighEnd_"+IDsuffix, leptons.at(1)->Pt(), weight, 10000, 0., 10000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(channels.at(it_ch)+"/Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      // Event selections for each CR
      for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

        // Cutflow : jet requirement (This is the number or events at preselection)
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

        // non-prompt CR1 : preselection with b-tagged jets
        if(it_rg == 0){
          if(!(Nbjet_medium > 0)) continue;
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.); 
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

        }
      
        // Low mass SR1, CR1 && High mass SR1, CR1
        if(it_rg>=1 && it_rg<5){

          if(it_rg < 3) continue; // No low mass for full Run2

          if(!(jets.size()>=2 && fatjets.size()==0)) continue; //JH : at least 2jet and no fatjet

          // Cutflow : jet requirement for SR1, CR1
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);

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
                j3 = k; j4 = l; //JH : this saves (k,l) tuple if that combination gives a smaller difference than the former combination
              }
            }
          }
          WCand1 = *leptons.at(0) + *leptons.at(1) + jets.at(j1) + jets.at(j2); // This is for Low mass
          WCand2 = jets.at(j3) + jets.at(j4); // This is for High mass
          lljj = *leptons.at(0) + *leptons.at(1) + jets.at(j3) + jets.at(j4); // High mass W jets
          l1jj = *leptons.at(0) + jets.at(j3) + jets.at(j4); // High mass W jets
          l2jj = *leptons.at(1) + jets.at(j3) + jets.at(j4);

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

          // Cut flow : selection for SR1, CR1
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Njets_"+IDsuffix, jets.size(), weight, 10, 0, 10);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Pt_HighEnd_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Pt_HighEnd_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet1_Pt_"+IDsuffix, jets.at(0).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet2_Pt_"+IDsuffix, jets.at(1).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet1_Pt_"+IDsuffix, jets.at(j3).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet2_Pt_"+IDsuffix, jets.at(j4).Pt(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet1_Mass_"+IDsuffix, jets.at(0).M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet2_Mass_"+IDsuffix, jets.at(1).M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet1_Mass_"+IDsuffix, jets.at(j3).M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet2_Mass_"+IDsuffix, jets.at(j4).M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet1_Eta_"+IDsuffix, jets.at(0).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Jet2_Eta_"+IDsuffix, jets.at(1).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet1_Eta_"+IDsuffix, jets.at(j3).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WJet2_Eta_"+IDsuffix, jets.at(j4).Eta(), weight, 50, -2.5, 2.5);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand1_Mass_"+IDsuffix, WCand1.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/WCand2_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);

          // Optional : Now start the optimization for High mass SR1
          
          if(it_rg == 3){

            if (HasFlag("M700")){
              if(!(jets.size()<4)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(jets.at(j3).Pt()>25.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(leptons.at(0)->Pt()>110.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(50.<WCand2.M()&&WCand2.M()<120.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 13.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 13.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(lljj.M()>800.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 14.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 14.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(370.<l1jj.M()&&l1jj.M()<885.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 15.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 15.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(MET2ST < 7.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 16.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 16.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/WCand_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Njets_"+IDsuffix, jets.size(), weight, 5, 0, 5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ptj1_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            }
            if (HasFlag("M1000")){
              if(!(jets.size()<4)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(jets.at(j3).Pt()>25.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(leptons.at(0)->Pt()>110.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(50.<WCand2.M()&&WCand2.M()<120.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 13.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 13.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(lljj.M()>800.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 14.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 14.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(370.<l1jj.M()&&l1jj.M()<1230.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 15.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 15.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(MET2ST < 7.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 16.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 16.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/WCand_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Njets_"+IDsuffix, jets.size(), weight, 5, 0, 5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ptj1_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            }
            if (HasFlag("M1500")){
              if(!(jets.size()<4)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(jets.at(j3).Pt()>25.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(leptons.at(0)->Pt()>110.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(50.<WCand2.M()&&WCand2.M()<120.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 13.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 13.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(lljj.M()>800.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 14.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 14.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(370.<l1jj.M()&&l1jj.M()<2220.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 15.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 15.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(MET2ST < 7.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 16.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 16.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/WCand_Mass_"+IDsuffix, WCand2.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/lljj_Mass_"+IDsuffix, lljj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1jj_Mass_"+IDsuffix, l1jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l2jj_Mass_"+IDsuffix, l2jj.M(), weight, 2000, 0., 2000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Njets_"+IDsuffix, jets.size(), weight, 5, 0, 5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ptj1_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
            }

          }

        }

        // Low mass SR2, CR2
        if(it_rg>=5 && it_rg<7){

          continue; // No low mass in full Run2

          if(!(jets.size()==1 && fatjets.size()==0)) continue; //JH : only 1 jet and no fatjet -> this jet is the proxy of WCand1
          llj = *leptons.at(0) + *leptons.at(1) + jets.at(0);
          l1j = *leptons.at(0) + jets.at(0);
          l2j = *leptons.at(1) + jets.at(0);

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

        }

        // High mass SR2, CR2
        if(it_rg >= 7){
          if(!(fatjets.size() > 0)) continue; //JH : at least 1 fatjet

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);

          double tmpMassDiff3 = 10000.;
          int j5 = 0;
          for(unsigned int k=0; k<fatjets.size(); k++){
            if(fabs(fatjets.at(k).M() - MW) < tmpMassDiff3){
              tmpMassDiff3 = fabs(fatjets.at(k).SDMass() - MW);
              j5 = k;
            }
          }
          WCand3 = fatjets.at(j5);
          l1J = *leptons.at(0) + fatjets.at(j5);
          l2J = *leptons.at(1) + fatjets.at(j5);

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

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
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

          //Start the optimization for High mass SR2

          if(it_rg == 7){

            if (HasFlag("M700")){
              if(!(leptons.at(0)->Pt()>140.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(40.<WCand3.M()&&WCand3.M()<130.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(635.<l1J.M()&&l1J.M()<825.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Fatjet_Mass_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
            }
            if (HasFlag("M1000")){
              if(!(leptons.at(0)->Pt()>140.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(40.<WCand3.M()&&WCand3.M()<130.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(900.<l1J.M()&&l1J.M()<1205.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Fatjet_Mass_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
            }
            if (HasFlag("M1500")){
              if(!(leptons.at(0)->Pt()>140.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(40.<WCand3.M()&&WCand3.M()<130.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
              if(!(1330.<l1J.M()&&l1J.M()<1800.)) continue;
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 12.5, weight, cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 12.5, 1., cutflow_bin, 0., cutflow_max);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/MET2ST_"+IDsuffix, MET2ST, weight, 16, 0, 16);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/Fatjet_Mass_"+IDsuffix, fatjets.at(j5).SDMass(), weight, 1000, 0., 1000.);
                FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Final/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
            }

          }


        } //JH : if highmass SR2, CR2

      } //JH : for loop in each region
    } //JH : if lepton size 2.
  } //JH : for loop in ee, mm, em channel 

}

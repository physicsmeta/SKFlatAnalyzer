#include "SampleValidation.h"

SampleValidation::SampleValidation(){

}

void SampleValidation::initializeAnalyzer(){

}

SampleValidation::~SampleValidation(){

  //==== Destructor of this Analyzer

}

void SampleValidation::executeEvent(){

  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  AllFatJets = GetAllFatJets();
  //AllFatJets = puppiCorr->Correct(GetAllFatJets()); //JH : puppiCorr = new FakeBackgroundEstimator(); in the constructor of AnalyzerCore.C; apply correction to fatjet.SDMass(); the total weight = gen correction * reco correction, from SKFlatAnalyzer/data/Run2Legacy_v4/DataYear/PuppiSoftdropMassCorr/puppiCorr.root

  AnalyzerParameter param;

  param.Muon_Tight_ID = "HNTightV1";
  param.Muon_Veto_ID  = "ISRVeto";
  param.Electron_Tight_ID = "HNTightV1";
  param.Electron_Veto_ID  = "ISRVeto";

  param.Jet_ID = "HNTight";
  if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
  else param.FatJet_ID = "HNTight0p45";

  executeEventFromParameter(param);

}

void SampleValidation::executeEventFromParameter(AnalyzerParameter param){

  double weight = 1.;
  double muon_miniaodP = 0.;
 
  Event ev = GetEvent();

  //=============
  //==== No Cut
  //=============

  if(!IsDATA){
    //weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full"); //JH : weight_norm_1invpb = xsec/sumW; Lumi = 35.9, 41.5, 59.7(fb-1) total 137fb-1
    if(ev.MCweight()<0) weight*=-1.;
    else if(ev.MCweight()>0) weight*=1.; //JH : gen_weight in MiniAOD
    //weight *= GetPrefireWeight(0); //JH : No issue in 2018, otherwise returns L1PrefireReweight_Central in MiniAOD
    //weight *= GetPileUpWeight(nPileUp,0); //JH : mcCorr->GetPileUpWeight(N_pileup, syst); mcCorr->GetPileUpWeight2017(N_pileup, syst); NOTE 2018 not yet added.
  } //JH : total weight calculation done.

  // Cutflow : No Cuts

  //========================
  //==== MET Filter
  //========================

  //if(!PassMETFilter()) return;

  //==============
  //==== Trigger
  //==============

  //if(!ev.PassTrigger(ElectronTriggers)) return;

  //======================
  //==== Copy AllObjects
  //======================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;
  vector<FatJet> this_AllFatJets = AllFatJets;
  vector<Gen> gens = GetGens();

  //Gen hard_l, HN_l, last_HN, W, ISR_parton;
  //vector<Gen> hard_partons, N_partons;

  //for(unsigned int i=0; i<gens.size(); i++){

    //FillHist("Pt_HN_mu", HN_mu.Pt(), weight, 1000, 0., 1000.);
    //FillHist("Eta_HN_mu", HN_mu.Eta(), weight, 100, -5., 5.);
    //FillHist("Phi_HN_mu", HN_mu.Phi(), weight, 63, -3.15, 3.15);

  //}


  cout << "================================================" << endl;
  cout << "all muon size: " << this_AllMuons.size() << endl;
  cout << "all electron size: " << this_AllElectrons.size() << endl;
  cout << "all jet size: " << this_AllJets.size() << endl;
  cout << "all fatjet size: " << this_AllFatJets.size() << endl;

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Muon> muons = this_AllMuons;
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 0., 10.);
  vector<Electron> electrons = this_AllElectrons;
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 0., 10.); //JH : lepton selection done

  //vector<Jet> jets_nolepveto = SelectJets(this_AllJets, "tight", 20., 2.7); //JH : to reject bjets
//  vector<FatJet> fatjets = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

//  FillHist("Njet_"+IDsuffix, jets_nolepveto.size(), weight, 8, 0., 8.);


  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  vector<Jet> jets;
  vector<FatJet> fatjets;
  jets.clear();
  fatjets.clear();
  int lepton_count1 = 0, lepton_count2 = 0, fatjet_count = 0; 

  // Fatjet selection
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

  cout << "jet size: " << jets.size() << endl;
  cout << "fatjet size: " << fatjets.size() << endl;

//  FillHist("Nfatjet_hn_"+IDsuffix, fatjets.size(), weight, 5, 0., 5.);
//  FillHist("Njet_hn_"+IDsuffix, jets.size(), weight, 8, 0., 8.); 

  //std::vector<Lepton*> leptons, leptons_minus, leptons_plus, leptons_veto;

  //=======================
  //==== Sort in pt-order
  //=======================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  //std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);

/*

  // muon selection

  Muon HN_mu, hard_mu;
  for(unsigned int i=0; i<muons.size(); i++){
    if (GetGenMatchedLepton(muons.at(i),gens).PID()!=0){
      Gen truth_mu = GetGenMatchedLepton(muons.at(i),gens);
      if (gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID() == 9900012){
        HN_mu = muons[i];
        cout << "HN_mu exists;" << endl;
        cout << "pt: " << HN_mu.Pt() << ", eta: " << HN_mu.Eta() << ", phi: " << HN_mu.Phi() << endl;
        cout << "matched gen Index: " << truth_mu.Index() << ", PID: " << truth_mu.PID() << ", mother's PID: " << gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID() << endl;
        cout << "pt: " << truth_mu.Pt() << ", eta: " << truth_mu.Eta() << ", phi: " << truth_mu.Phi() << endl;
      }
      if ((abs(gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID()) == 24) || (abs(gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID()) <= 4) || abs(gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID()) == 21 ){
        hard_mu = muons[i];
        cout << "hard_mu exists;" << endl;
        cout << "pt: " << hard_mu.Pt() << ", eta: " << hard_mu.Eta() << ", phi: " << hard_mu.Phi() << endl;
        cout << "matched gen Index: " << truth_mu.Index() << ", PID: " << truth_mu.PID() << ", mother's PID: " << gens.at(TrackGenSelfHistory(truth_mu,gens).at(1)).PID() << endl;
        cout << "pt: " << truth_mu.Pt() << ", eta: " << truth_mu.Eta() << ", phi: " << truth_mu.Phi() << endl;
      }
    }
  }

  if(HN_mu.Chi2()!=999.){
    FillHist("Pt_HN_mu", HN_mu.Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_HN_mu", HN_mu.Eta(), weight, 100, -5., 5.);
    FillHist("Phi_HN_mu", HN_mu.Phi(), weight, 63, -3.15, 3.15);
    if(fabs(HN_mu.Eta())<2.4){
      FillHist("fid_Pt_HN_mu", HN_mu.Pt(), weight, 1000, 0., 1000.);
      FillHist("fid_Eta_HN_mu", HN_mu.Eta(), weight, 100, -5., 5.);
      FillHist("fid_Phi_HN_mu", HN_mu.Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(hard_mu.Chi2()!=999.){
    FillHist("Pt_hard_mu", hard_mu.Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_hard_mu", hard_mu.Eta(), weight, 100, -5., 5.);
    FillHist("Phi_hard_mu", hard_mu.Phi(), weight, 63, -3.15, 3.15);
    if(fabs(hard_mu.Eta())<2.4){
      FillHist("fid_Pt_hard_mu", hard_mu.Pt(), weight, 1000, 0., 1000.);
      FillHist("fid_Eta_hard_mu", hard_mu.Eta(), weight, 100, -5., 5.);
      FillHist("fid_Phi_hard_mu", hard_mu.Phi(), weight, 63, -3.15, 3.15);
    }
  }

  if(muons.size()>0){
    FillHist("Pt_mu1", muons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_mu1", muons.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_mu1", muons.at(0).Phi(), weight, 63, -3.15, 3.15);
    if(fabs(muons[0].Eta())<2.4){
      FillHist("fid_Pt_mu1", muons.at(0).Pt(), weight, 1000, 0., 1000.);
      FillHist("fid_Eta_mu1", muons.at(0).Eta(), weight, 100, -5., 5.);
      FillHist("fid_Phi_mu1", muons.at(0).Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(muons.size()>1){
    FillHist("Pt_mu2", muons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_mu2", muons.at(1).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_mu2", muons.at(1).Phi(), weight, 63, -3.15, 3.15);
    if(fabs(muons[1].Eta())<2.4){
      FillHist("fid_Pt_mu2", muons.at(1).Pt(), weight, 1000, 0., 1000.);
      FillHist("fid_Eta_mu2", muons.at(1).Eta(), weight, 100, -5., 5.);
      FillHist("fid_Phi_mu2", muons.at(1).Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(electrons.size()>0){
    FillHist("Pt_e1", electrons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_e1", electrons.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_e1", electrons.at(0).Phi(), weight, 63, -3.15, 3.15);
  }
  if(electrons.size()>1){
    FillHist("Pt_e2", electrons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_e2", electrons.at(1).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_e2", electrons.at(1).Phi(), weight, 63, -3.15, 3.15);
  }
  if(jets.size()>0){
    if(jets.at(0).Pt()>20.&&fabs(jets.at(0).Eta())<2.7){
      FillHist("Pt_j1", jets.at(0).Pt(), weight, 1000, 0., 1000.);
      FillHist("Eta_j1", jets.at(0).Eta(), weight, 100, -5., 5.);
      FillHist("Phi_j1", jets.at(0).Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(jets.size()>1){
    if(jets.at(1).Pt()>20.&&fabs(jets.at(1).Eta())<2.7){
      FillHist("Pt_j2", jets.at(1).Pt(), weight, 1000, 0., 1000.);
      FillHist("Eta_j2", jets.at(1).Eta(), weight, 100, -5., 5.);
      FillHist("Phi_j2", jets.at(1).Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(fatjets.size()>0){
    FillHist("Pt_fatjet", fatjets.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_fatjet", fatjets.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_fatjet", fatjets.at(0).Phi(), weight, 63, -3.15, 3.15);
    FillHist("M_fatjet", fatjets.at(0).M(), weight, 1000, 0., 1000.);
  }

*/

  Gen *last_HN=NULL;
  Gen *hard_l=NULL;
  Gen *HN_l=NULL;
  Gen *W_l=NULL;
  vector<Gen*> hard_Ws, hard_partons, leptons, N_partons;

  for (unsigned int i=0; i<gens.size(); i++){
    if(gens[i].isHardProcess()){
      if(abs(gens[i].PID())==24) hard_Ws.push_back(&gens[i]);
      else if(abs(gens[i].PID())<=4||gens[i].PID()==21) hard_partons.push_back(&gens[i]);
    }
  }
  PrintGen(gens);
  for(unsigned int i=0; i<gens.size(); i++){
    cout << i << "th particle:" << endl;
    if(gens[i].PID()==9900012) last_HN=&gens[i];
    if(last_HN) cout << "^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~this is the last_HN : " << last_HN << endl;
    if(((abs(gens[i].PID())==11)||(abs(gens[i].PID())==13))&&((gens[i].MotherIndex()==hard_partons.at(0)->Index())||(abs(gens[gens[i].MotherIndex()].PID())==24))){
      hard_l=&gens[i]; 
      cout << "^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~this is the hard_l : " << hard_l << endl;
    }
    if(last_HN){
      if((abs(gens[i].PID())==11||abs(gens[i].PID())==13)&&(gens[i].MotherIndex()==last_HN->Index())){
        HN_l=&gens[i];
        cout << "^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~this is the HN_l : " << HN_l << endl;
      }
    }
    if((abs(gens[i].PID())==11||abs(gens[i].PID())==13)&&(abs(gens[gens[i].MotherIndex()].PID())==24)){
      W_l=&gens[i];
      cout << "^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~this is the W_l : " << W_l << endl;
    }
    if(gens[i].isPromptFinalState()){
      if(abs(gens[i].PID())==11||abs(gens[i].PID())==13) leptons.push_back(&gens[i]);
    }    
  }

  for(int i=0;i<hard_partons.size();i++){
    if(abs((gens[hard_partons.at(i)->MotherIndex()]).PID())==24||(gens[hard_partons.at(i)->MotherIndex()]).PID()==9900012) N_partons.push_back(hard_partons.at(i));
  }

  if(hard_Ws.size()>0){
    for(int i=0;i<hard_Ws.size();i++) hard_Ws.at(i) = FindLastCopy(hard_Ws.at(i),gens);
  }
  if(hard_l) hard_l = FindLastCopy(hard_l,gens);
  cout << "hard_l: " << hard_l << endl;
  if(W_l) W_l = FindLastCopy(W_l,gens);
  cout << "W_l: " << W_l << endl;

  FillHist("Gen_Pt_Nlepton", HN_l->Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_Nlepton", HN_l->Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_Nlepton", HN_l->Phi(), weight, 63, -3.15, 3.15); 
  FillHist("Gen_Mass_Nreco", (*HN_l+*N_partons.at(0)+*N_partons.at(1)).M(), weight, 1000, 0., 1000.); 

/*

  vector<Gen> gen_muons, gen_N_partons;
  Gen gen_N, gen_HN_mu, gen_hard_mu;
  for (unsigned int i=0; i<gens.size(); i++){
    if((abs(gens[i].PID())==13)&&gens[i].fromHardProcessFinalState()) gen_muons.push_back(gens[i]);
    if((abs(gens[i].PID())==13)&&gens[i].fromHardProcessFinalState()&&((abs(gens.at(TrackGenSelfHistory(gens[i],gens).at(1)).PID()) == 24) || (abs(gens.at(TrackGenSelfHistory(gens[i],gens).at(1)).PID()) <= 4) || (abs(gens.at(TrackGenSelfHistory(gens[i],gens).at(1)).PID()) == 21))) gen_hard_mu = gens[i];
    if((abs(gens[i].PID())==13)&&gens[i].fromHardProcessFinalState()&&(gens.at(TrackGenSelfHistory(gens[i],gens).at(1)).PID() == 9900012)) gen_HN_mu = gens[i];
    if(gens[i].PID()==9900012) gen_N = gens[i];
    if((abs(gens[i].PID())<=4||gens[i].PID()==21)&&gens[i].isHardProcess()&&(abs((gens[gens[i].MotherIndex()]).PID())==24||(gens[gens[i].MotherIndex()]).PID()==9900012)) gen_N_partons.push_back(gens[i]);
  }

  FillHist("Gen_Pt_HN_mu", gen_HN_mu.Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_HN_mu", gen_HN_mu.Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_HN_mu", gen_HN_mu.Phi(), weight, 63, -3.15, 3.15);
  if(fabs(gen_HN_mu.Eta())<2.4){
    FillHist("fid_Gen_Pt_HN_mu", gen_HN_mu.Pt(), weight, 1000, 0., 1000.);
    FillHist("fid_Gen_Eta_HN_mu", gen_HN_mu.Eta(), weight, 100, -5., 5.);
    FillHist("fid_Gen_Phi_HN_mu", gen_HN_mu.Phi(), weight, 63, -3.15, 3.15);
  }
  FillHist("Gen_Pt_hard_mu", gen_hard_mu.Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_hard_mu", gen_hard_mu.Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_hard_mu", gen_hard_mu.Phi(), weight, 63, -3.15, 3.15);
  if(fabs(gen_hard_mu.Eta())<2.4){
    FillHist("fid_Gen_Pt_hard_mu", gen_hard_mu.Pt(), weight, 1000, 0., 1000.);
    FillHist("fid_Gen_Eta_hard_mu", gen_hard_mu.Eta(), weight, 100, -5., 5.);
    FillHist("fid_Gen_Phi_hard_mu", gen_hard_mu.Phi(), weight, 63, -3.15, 3.15);
  }
  
  std::sort(gen_muons.begin(), gen_muons.end(), PtComparing);
  std::sort(gen_N_partons.begin(), gen_N_partons.end(), PtComparing);
  FillHist("Gen_Pt_mu1", gen_muons.at(0).Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_mu1", gen_muons.at(0).Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_mu1", gen_muons.at(0).Phi(), weight, 63, -3.15, 3.15);
  if(fabs(gen_muons.at(0).Eta())<2.4){
    FillHist("fid_Gen_Pt_mu1", gen_muons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("fid_Gen_Eta_mu1", gen_muons.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("fid_Gen_Phi_mu1", gen_muons.at(0).Phi(), weight, 63, -3.15, 3.15);
  }
  FillHist("Gen_Pt_mu2", gen_muons.at(1).Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_mu2", gen_muons.at(1).Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_mu2", gen_muons.at(1).Phi(), weight, 63, -3.15, 3.15);
  if(fabs(gen_muons.at(1).Eta())<2.4){
    FillHist("fid_Gen_Pt_mu2", gen_muons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("fid_Gen_Eta_mu2", gen_muons.at(1).Eta(), weight, 100, -5., 5.);
    FillHist("fid_Gen_Phi_mu2", gen_muons.at(1).Phi(), weight, 63, -3.15, 3.15);
  }
  if(gen_N_partons.at(0).Pt()>20.&&fabs(gen_N_partons.at(0).Eta())<2.7){
    FillHist("Gen_Pt_j1", gen_N_partons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Gen_Eta_j1", gen_N_partons.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Gen_Phi_j1", gen_N_partons.at(0).Phi(), weight, 63, -3.15, 3.15);
  }
  if(gen_N_partons.at(1).Pt()>20.&&fabs(gen_N_partons.at(1).Eta())<2.7){
    FillHist("Gen_Pt_j2", gen_N_partons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("Gen_Eta_j2", gen_N_partons.at(1).Eta(), weight, 100, -5., 5.);
    FillHist("Gen_Phi_j2", gen_N_partons.at(1).Phi(), weight, 63, -3.15, 3.15);
  }
  FillHist("Gen_Pt_N", gen_N.Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_N", gen_N.Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_N", gen_N.Phi(), weight, 63, -3.15, 3.15);
  FillHist("Gen_M_N", gen_N.M(), weight, 1000, 0., 1000.);

  FillHist("nPileUp", nPileUp, weight, 100, 0., 100.);

  vector<Jet> W_jets1, W_jets2;
  cout << "For each jet:" << endl;
  for(unsigned int i=0; i<jets.size(); i++){
    cout << "pt: " << jets[i].Pt() << ", eta: " << jets[i].Eta() << ", phi: " << jets[i].Phi() << endl;
    if (jets[i].DeltaR(gen_N_partons[0])<0.3) W_jets1.push_back(jets[i]);
    if (jets[i].DeltaR(gen_N_partons[1])<0.3) W_jets2.push_back(jets[i]);
  }
  std::sort(W_jets1.begin(), W_jets1.end(), PtComparing);
  auto last1 = unique(W_jets1.begin(), W_jets1.end());
  W_jets1.erase(last1, W_jets1.end());
  std::sort(W_jets2.begin(), W_jets2.end(), PtComparing);
  auto last2 = unique(W_jets2.begin(), W_jets2.end());
  W_jets2.erase(last2, W_jets2.end());
  if(W_jets1.size()>0){

    cout << "Number of first parton matched jet: " << W_jets1.size() << endl;
    for(unsigned int i=0; i<W_jets1.size(); i++){
      cout << "pt: " << W_jets1[i].Pt() << ", eta: " << W_jets1[i].Eta() << ", phi: " << W_jets1[i].Phi() << endl;
    }
    if(W_jets1.at(0).Pt()>20.&&fabs(W_jets1.at(0).Eta())<2.7){
      FillHist("Pt_Wjet1", W_jets1.at(0).Pt(), weight, 1000, 0., 1000.);
      FillHist("Eta_Wjet1", W_jets1.at(0).Eta(), weight, 100, -5., 5.);
      FillHist("Phi_Wjet1", W_jets1.at(0).Phi(), weight, 63, -3.15, 3.15);
    }
  }
  if(W_jets2.size()>0){

    cout << "Number of second parton matched jet: " << W_jets2.size() << endl;
    for(unsigned int i=0; i<W_jets2.size(); i++){
      cout << "pt: " << W_jets2[i].Pt() << ", eta: " << W_jets2[i].Eta() << ", phi: " << W_jets2[i].Phi() << endl;
    }

    if(W_jets1.size()==0){
      if(W_jets2.at(0).Pt()>20.&&fabs(W_jets2.at(0).Eta())<2.7){
        FillHist("Pt_Wjet2", W_jets2.at(0).Pt(), weight, 1000, 0., 1000.);
        FillHist("Eta_Wjet2", W_jets2.at(0).Eta(), weight, 100, -5., 5.);
        FillHist("Phi_Wjet2", W_jets2.at(0).Phi(), weight, 63, -3.15, 3.15);
      }
    }
    else if(W_jets1.size()>0&&(&W_jets2.at(0)!=&W_jets1.at(0))){
      if(W_jets2.at(0).Pt()>20.&&fabs(W_jets2.at(0).Eta())<2.7){
        FillHist("Pt_Wjet2", W_jets2.at(0).Pt(), weight, 1000, 0., 1000.);
        FillHist("Eta_Wjet2", W_jets2.at(0).Eta(), weight, 100, -5., 5.);
        FillHist("Phi_Wjet2", W_jets2.at(0).Phi(), weight, 63, -3.15, 3.15);
      }
    }
  }

  for(unsigned int i=0; i<gen_N_partons.size(); i++){
    cout << "W parton Index: " << gen_N_partons[i].Index() << ", PID: " << gen_N_partons[i].PID() << endl;
    cout << "pt: " << gen_N_partons[i].Pt() << ", eta: " << gen_N_partons[i].Eta() << ", phi: " << gen_N_partons[i].Phi() << endl;
  }
  cout << "dR: " << gen_N_partons[0].DeltaR(gen_N_partons[1]) << endl;

  PrintGen(gens);

*/

}

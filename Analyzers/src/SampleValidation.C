#include "SampleValidation.h"

SampleValidation::SampleValidation(){

}

void SampleValidation::initializeAnalyzer(){

}

SampleValidation::~SampleValidation(){

  //==== Destructor of this Analyzer

}

void SampleValidation::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/SampleValidation.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  AllFatJets = GetAllFatJets();
  //AllFatJets = puppiCorr->Correct(GetAllFatJets()); //JH : puppiCorr = new FakeBackgroundEstimator(); in the constructor of AnalyzerCore.C; apply correction to fatjet.SDMass(); the total weight = gen correction * reco correction, from SKFlatAnalyzer/data/Run2Legacy_v4/DataYear/PuppiSoftdropMassCorr/puppiCorr.root

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/SampleValidation.h
//  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  param.Muon_Tight_ID = "POGTight";
  param.Muon_Veto_ID  = "POGLoose";
  param.Electron_Tight_ID = "passTightID";
  param.Electron_Veto_ID  = "passVetoID";

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

  cout << "================================================" << endl;
  cout << "all muon size: " << this_AllMuons.size() << endl;
  cout << "all electron size: " << this_AllElectrons.size() << endl;
  cout << "all jet size: " << this_AllJets.size() << endl;
  cout << "all fatjet size: " << this_AllFatJets.size() << endl;

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
    cout << "[SampleValidation::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

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

  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124
  for(unsigned int i=0; i<this_AllFatJets.size(); i++){
    lepton_count1 = 0;
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
    for(unsigned int j=0; j<muons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(muons_veto.at(j)) < 0.4) lepton_count2++; //JH : muon cleaning
    }
    for(unsigned int j=0; j<electrons_veto.size(); j++){
      if(this_AllJets.at(i).DeltaR(electrons_veto.at(j)) < 0.4) lepton_count2++; //JH : electron cleaning
    }
    //for(unsigned int j=0; j<fatjets.size(); j++){
    //  if(this_AllJets.at(i).DeltaR(fatjets.at(j)) < 0.8) fatjet_count++; //JH : fatjet cleaning
    //}
    if(lepton_count2 > 0) continue;
    //if(fatjet_count > 0) continue;
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
  if(electrons.size()>0) std::sort(electrons.begin(), electrons.end(), PtComparing);
  if(electrons_veto.size()>0) std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  //std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  if(fatjets.size()>0) std::sort(fatjets.begin(), fatjets.end(), PtComparing);

  ////==== B-Tagging 
  //int Nbjet_loose = 0, Nbjet_medium = 0;
  //JetTagging::Parameters jtp_DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV,
  //                                                                   JetTagging::Loose,
  //                                                                   JetTagging::incl, JetTagging::comb);
  //JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV,
  //                                                                   JetTagging::Medium,
  //                                                                   JetTagging::incl, JetTagging::comb); //JH : Set b-tagging parameters

  ////==== method 1a)
  ////==== multiply "btagWeight" to the event weight
////  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  ////==== method 2a)
  //for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){ 
////    double this_discr = jets_nolepveto.at(ij).GetTaggerResult(JetTagging::DeepCSV);
  //    //==== No SF
////      if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
  //  if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++; 
  //  if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++; //JH : count Nbjet. NOTE : AN says they used CVSv2 and medium WP.
  //} 

//  FillHist("Nbjet_loose_"+IDsuffix, Nbjet_loose, weight, 5, 0., 5.);
//  FillHist("Nbjet_medium_"+IDsuffix, Nbjet_medium, weight, 5, 0., 5.);

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
  if(jets.size()>1){
    FillHist("Pt_j1", jets.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_j1", jets.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_j1", jets.at(0).Phi(), weight, 63, -3.15, 3.15);
  }
  if(jets.size()>1){
    FillHist("Pt_j2", jets.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_j2", jets.at(1).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_j2", jets.at(1).Phi(), weight, 63, -3.15, 3.15);
  }
  if(fatjets.size()>0){
    FillHist("Pt_fatjet", fatjets.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_fatjet", fatjets.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_fatjet", fatjets.at(0).Phi(), weight, 63, -3.15, 3.15);
    FillHist("M_fatjet", fatjets.at(0).M(), weight, 1000, 0., 1000.);
  }

  //Gen *last_HN, *hard_l, *HN_l, *W_l;
  //vector<Gen*> hard_Ws, hard_partons, leptons, N_partons;

  //for (unsigned int i=0; i<gens.size(); i++){
  //  if(gens[i].isHardProcess()){
  //    if(abs(gens[i].PID())==24) hard_Ws.push_back(&gens[i]);
  //    else if(abs(gens[i].PID())<=4||gens[i].PID()==21) hard_partons.push_back(&gens[i]);
  //    cout << "a" << endl;
  //  }
  //}
	//PrintGen(gens);
  //for(unsigned int i=0; i<gens.size(); i++){
	//	cout << "&gens[i]: " << &gens[i] << endl;
  //  if(gens[i].PID()==9900012) last_HN=&gens[i];
	//	cout << "last_HN: " << last_HN << endl;
  //  cout << i << "th particle:" << endl;
  //  cout << "b" << endl;
	//	cout << "gens[i].PID(): " << gens[i].PID() << endl;
	//	cout << "abs(gens[i].PID()): " << abs(gens[i].PID()) << endl;
	//	cout << "gens[i].MotherIndex(): " << gens[i].MotherIndex() << endl;
	//	cout << "hard_partons.at(0)->Index(): " << hard_partons.at(0)->Index() << endl;
	//	cout << "((abs(gens[i].PID())==11)||(abs(gens[i].PID())==13)): " << ((abs(gens[i].PID())==11)||(abs(gens[i].PID())==13)) << endl;
  //  cout << "(gens[i].MotherIndex()==hard_partons.at(0)->Index()): " << (gens[i].MotherIndex()==hard_partons.at(0)->Index()) << endl;
  //  cout << (((abs(gens[i].PID())==11)||(abs(gens[i].PID())==13))&&(gens[i].MotherIndex()==hard_partons.at(0)->Index())) << endl;
  //  cout <<"bb" << endl;
  //  if(((abs(gens[i].PID())==11)||(abs(gens[i].PID())==13))&&((gens[i].MotherIndex()==hard_partons.at(0)->Index())||(abs(gens[gens[i].MotherIndex()].PID())==24))){
  //    cout << "C" << endl;
  //    hard_l=&gens[i]; 
  //    cout << "D" << endl;
  //  }
	//	if(last_HN->PID()!=0){
	//		cout << "E" << endl;
  //    if((abs(gens[i].PID())==11||abs(gens[i].PID())==13)&&(gens[i].MotherIndex()==last_HN->Index())){
  //      HN_l=&gens[i];
  //      cout << "e" << endl;
  //    }
	//	}
  //  if((abs(gens[i].PID())==11||abs(gens[i].PID())==13)&&(abs(gens[gens[i].MotherIndex()].PID())==24)){
  //    W_l=&gens[i];
  //    cout << "f" << endl;
  //  }
  //  if(gens[i].isPromptFinalState()){
  //    if(abs(gens[i].PID())==11||abs(gens[i].PID())==13) leptons.push_back(&gens[i]);
  //    cout << "g" << endl;
  //  }    
  //}

  //for(int i=0;i<hard_partons.size();i++){
  //  if(abs((gens[hard_partons.at(i)->MotherIndex()]).PID())==24||(gens[hard_partons.at(i)->MotherIndex()]).PID()==9900012) N_partons.push_back(hard_partons.at(i));
  //  cout << "h" << endl;
  //}

  //if(hard_Ws.size()>0){
  //  for(int i=0;i<hard_Ws.size();i++) hard_Ws.at(i) = FindLastCopy(hard_Ws.at(i),gens);
  //  cout << "i" << endl;
  //}
  //cout << "hard_l: " << hard_l << endl;
  //if(hard_l->PID()!=0) hard_l = FindLastCopy(hard_l,gens);
  //cout << "j" << endl;
	//cout << "W_l: " << W_l << endl;
  //if(W_l->PID()!=0) W_l = FindLastCopy(W_l,gens);
  //cout << "k" << endl;

  //FillHist("Gen_Pt_Nlepton", HN_l->Pt(), weight, 1000, 0., 1000.);
  //cout << "l" << endl;
  //FillHist("Gen_Eta_Nlepton", HN_l->Eta(), weight, 100, -5., 5.);
  //FillHist("Gen_Phi_Nlepton", HN_l->Phi(), weight, 63, -3.15, 3.15); 

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
  FillHist("Gen_Pt_j1", gen_N_partons.at(0).Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_j1", gen_N_partons.at(0).Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_j1", gen_N_partons.at(0).Phi(), weight, 63, -3.15, 3.15);
  FillHist("Gen_Pt_j2", gen_N_partons.at(1).Pt(), weight, 1000, 0., 1000.);
  FillHist("Gen_Eta_j2", gen_N_partons.at(1).Eta(), weight, 100, -5., 5.);
  FillHist("Gen_Phi_j2", gen_N_partons.at(1).Phi(), weight, 63, -3.15, 3.15);
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

    FillHist("Pt_Wjet1", W_jets1.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Eta_Wjet1", W_jets1.at(0).Eta(), weight, 100, -5., 5.);
    FillHist("Phi_Wjet1", W_jets1.at(0).Phi(), weight, 63, -3.15, 3.15);
  }
  if(W_jets2.size()>0){

    cout << "Number of second parton matched jet: " << W_jets2.size() << endl;
    for(unsigned int i=0; i<W_jets2.size(); i++){
      cout << "pt: " << W_jets2[i].Pt() << ", eta: " << W_jets2[i].Eta() << ", phi: " << W_jets2[i].Phi() << endl;
    }

    if(W_jets1.size()==0){
      FillHist("Pt_Wjet2", W_jets2.at(0).Pt(), weight, 1000, 0., 1000.);
      FillHist("Eta_Wjet2", W_jets2.at(0).Eta(), weight, 100, -5., 5.);
      FillHist("Phi_Wjet2", W_jets2.at(0).Phi(), weight, 63, -3.15, 3.15);
    }
		else if(W_jets1.size()>0&&(&W_jets2.at(0)!=&W_jets1.at(0))){
      FillHist("Pt_Wjet2", W_jets2.at(0).Pt(), weight, 1000, 0., 1000.);
      FillHist("Eta_Wjet2", W_jets2.at(0).Eta(), weight, 100, -5., 5.);
      FillHist("Phi_Wjet2", W_jets2.at(0).Phi(), weight, 63, -3.15, 3.15);
    }
	}

  for(unsigned int i=0; i<gen_N_partons.size(); i++){
    cout << "W parton Index: " << gen_N_partons[i].Index() << ", PID: " << gen_N_partons[i].PID() << endl;
    cout << "pt: " << gen_N_partons[i].Pt() << ", eta: " << gen_N_partons[i].Eta() << ", phi: " << gen_N_partons[i].Phi() << endl;
  }
	cout << "dR: " << gen_N_partons[0].DeltaR(gen_N_partons[1]) << endl;

  PrintGen(gens);

/*
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
  double mu_tight_iso = 0.15;
  //if(IDsuffix == "HNV2") mu_tight_iso = 0.1;
  if(IDsuffix == "HN16") mu_tight_iso = 0.07;

  double el_tight_iso = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  // Set pTcone
  for(unsigned int i=0; i<muons.size(); i++){
    this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso); //JH : CalcPtCone() in Lepton.h; this returns (i) pt for more tightly isolated leptons than the tight_iso, or (ii) pt + pt*(RelIso-tight_iso) which is the proxy for the mother parton's pt -> used for fake estimation
    muons.at(i).SetPtCone(this_ptcone_muon);
  }
   
  for(unsigned int i=0; i<electrons.size(); i++){
    //el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt(); //JH : TODO electron uses UncorrPt() but I don't understand the meaning yet
    //if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
    //if(IDsuffix == "HNV2"){
    //  el_tight_iso = std::min(0.08, 0.0287+0.506/electrons.at(i).UncorrPt());
    //  if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = std::min(0.08, 0.0445+0.963/electrons.at(i).UncorrPt());
    //} 
    if(IDsuffix == "HN16") el_tight_iso = 0.08;
    this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
    electrons.at(i).SetPtCone(this_ptcone_electron);
  }

  if(muons.size()==2 && electrons.size()==0){
    FillHist("Pt_muon1", muons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_muon2", muons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon1", muons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon2", muons.at(1).PtCone(), weight, 1000, 0., 1000.);
  }
  if(muons.size()==0 && electrons.size()==2){
    FillHist("Pt_electron1", electrons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_electron2", electrons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron1", electrons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron2", electrons.at(1).PtCone(), weight, 1000, 0., 1000.);
  } //JH : Draw lepton pt and ptcone

//  if(electrons.size() > 0) cout << electrons.at(0).PtCone() << endl;

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

  //=====================================================================================
  //=====================================================================================
  //==== Preselection, low/high mass signal regions
  //=====================================================================================
  //=====================================================================================

  //=========================
  //==== Event selections..
  //=========================

  LeptonPtCut1 = ElectronPtCut1; LeptonPtCut2 = ElectronPtCut2;

  // Period-dependent trigger weight (only for 2016 MC)
  trigger_lumi = 1.;
  if(!IsDATA) trigger_lumi = ev.GetTriggerLumi("Full");

  // Cutflow : passing dilepton triggers
  weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb*trigger_lumi;
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);
  } //JH : recalculate total weight for 2016 period dependency.
  FillHist("ee/DY/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
  FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  FillHist("ee/Pre/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
  FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);

*/

/*

  if(leptons.size() == 2){
    if(!(muons.size()==0 && electrons.size()==2)) return;

    ZCand = *leptons.at(0) + *leptons.at(1);

    weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
    // weights for MC
    if(!IsDATA){
      
      // Select prompt only
      if(GetLeptonType(*leptons.at(0), gens)<=0) return;
      if(GetLeptonType(*leptons.at(1), gens)<=0) return;

      weight *= weight_norm_1invpb*trigger_lumi; //JH : trigger_lumi for period dependency
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);

      for(unsigned int i=0; i<muons.size(); i++){
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
        weight *= muon_recosf*muon_idsf*muon_isosf; 
      }
      for(unsigned int j=0; j<electrons.size(); j++){
        ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
        if(param.Electron_Tight_ID.Contains("HEEP")){
          ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
        }
        else ele_idsf = 1.;
        weight *= ele_recosf*ele_idsf;
      }
    } //JH : Now total weight including trigger lumi and lepton SF done && lepton gen-matching (for first 2 leptons only) done

    /////////////////////////////////////////////////////////
    //// Preselection (triggers have been already applied.)
    /////////////////////////////////////////////////////////

    if(!(leptons.at(0)->Pt()>LeptonPtCut1 && leptons.at(1)->Pt()>LeptonPtCut2)) return;
    
    // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)
    FillHist("ee/DY/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

    if(lepton_veto_size > 0) return;

    // Cutflow : veto 3rd leptons using veto ID
    FillHist("ee/DY/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

    if(!(ZCand.M() > 10.)) return; 

    // Cutflow : m(ll) > 10 GeV 
    FillHist("ee/DY/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

    // Now separate OS / SS

    if(leptons.at(0)->Charge()*leptons.at(1)->Charge()<0){
      // Cutflow : opposite sign 
      FillHist("ee/DY/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(Nbjet_medium == 0)) return;

      // Cutflow : No b jets
      FillHist("ee/DY/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist("ee/DY/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

      FillHist("ee/DY/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist("ee/DY/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/DY/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/DY/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/DY/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("ee/DY/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("ee/DY/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist("ee/DY/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
    } // OS (DY CR)

    else if (leptons.at(0)->Charge()*leptons.at(1)->Charge()>0){
      // Cutflow : same sign 
      FillHist("ee/Pre/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(IsOnZ(ZCand.M(), 10.)) return; //JH : see p.12 of preapproval -> https://indico.cern.ch/event/694943/contributions/2849972/attachments/1583026/2501796/180115__JaesungKim__JetsX_Meeting__HN_DiLepton_PreApproval.pdf

      // Cutflow : |m(ll)-m(Z)| > 10 GeV for ee 
      FillHist("ee/Pre/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);

      if(!( fatjets.size()>0 || (jets.size()>1 && fatjets.size()==0) || (jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.) )) return; //JH : jet requirement
     
      // Cutflow : jet requirement (This is the number or events at preselection)
      FillHist("ee/Pre/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist("ee/Pre/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

      FillHist("ee/Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist("ee/Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist("ee/Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("ee/Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist("ee/Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist("ee/Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

    } // SS (preselection)

  } //JH : if lepton size 2.

*/

  //=====================================================================================
  //=====================================================================================
  //==== SM background CR (DYmm, DYee, DYemu, WZ, ZG, WG, ZZ)
  //=====================================================================================
  //=====================================================================================

/*

  if(RunCF) return;

  if(IsDATA){
    //if(DataStream.Contains("DoubleMuon") && !ev.PassTrigger(MuonTriggersH)) return;
    if(DataStream.Contains("DoubleMuon") && !ev.PassTrigger(MuonTriggers)) return;
    if(DataYear==2016 || DataYear==2017){
      if(DataStream.Contains("DoubleEG")){
        //if(ev.PassTrigger(MuonTriggersH) || !ev.PassTrigger(ElectronTriggers)) return;
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
    if(DataYear==2018){
      if(DataStream.Contains("EGamma")){
        //if(ev.PassTrigger(MuonTriggersH) || !ev.PassTrigger(ElectronTriggers)) return;
        if(ev.PassTrigger(MuonTriggers) || !ev.PassTrigger(ElectronTriggers)) return;
      }
    }
    if(DataStream.Contains("MuonEG")){
      //if(ev.PassTrigger(MuonTriggersH) || ev.PassTrigger(ElectronTriggers) || !ev.PassTrigger(EMuTriggersH)) return;
      if(ev.PassTrigger(MuonTriggers) || ev.PassTrigger(ElectronTriggers) || !ev.PassTrigger(EMuTriggers)) return;
    }
  } //JH : I'm not sure DataStream info is necessary.. any possible bias?

  // Period-dependent trigger weight (only for 2016 MC)
  trigger_lumi = 1., dimu_trig_weight = 0., emu_trig_weight = 0.;
  if(!IsDATA){
    if(DataYear==2016){
      if(ev.PassTrigger(MuonTriggers)){ 
        dimu_trig_weight += 27267.591;
        trigger_lumi = dimu_trig_weight;
      }
      if(ev.PassTrigger(MuonTriggersH)){
        dimu_trig_weight += 8650.628;
        trigger_lumi = dimu_trig_weight;
      }
      if(!ev.PassTrigger(MuonTriggers)&&ev.PassTrigger(ElectronTriggers)){
        trigger_lumi = ev.GetTriggerLumi("Full"); 
      }
      if(!ev.PassTrigger(MuonTriggers)&&!ev.PassTrigger(ElectronTriggers)&&ev.PassTrigger(EMuTriggers)){
        emu_trig_weight += 27267.591;
        trigger_lumi = emu_trig_weight;
      }
      if(!ev.PassTrigger(MuonTriggers)&&!ev.PassTrigger(ElectronTriggers)&&ev.PassTrigger(EMuTriggersH)){
        emu_trig_weight += 8650.628;
        trigger_lumi = emu_trig_weight;
      } //JH : Need to check this is alright
    }
    else{
      trigger_lumi = ev.GetTriggerLumi("Full");
    }
  }

  //=========================
  //==== Event selections..
  //========================= 
  
  for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
    weight = 1., muon_recosf = 1., muon_idsf = 1., muon_isosf = 1., ele_idsf = 1., ele_recosf = 1.;
    ossf_mass10 = 0;

    if(!IsDATA){
      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);
    }

    // Cutflow : passing dilepton triggers (dimu || diel || emu)
    for(unsigned int it_rg2=0; it_rg2<regionsSM.size(); it_rg2++){
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }

    // Requirements after passing triggers
    if(ev.PassTrigger(MuonTriggers)){
      if(muons.size() < 2) continue;
      if(muons.size()>=2 && !(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2)) continue;
    }
    else if(ev.PassTrigger(ElectronTriggers)){
      if(electrons.size() < 2) continue;
      if(electrons.size()>=2 && !(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2)) continue;
    }
    else if(ev.PassTrigger(EMuTriggers)){
      if(muons.size()*electrons.size() == 0) continue;
      if(muons.size()*electrons.size()>0 && !(leptons.at(0)->Pt()>EMuPtCut1 && leptons.at(1)->Pt()>EMuPtCut2)) continue;
    } // JH : Note that 2016BtoG trigger is already a superset of 2016H trigger, so no need to set BtoG

    //=====================================
    //==== DY control region
    //=====================================
    if(it_rg2<3 && leptons.size()==2){ //JH : DYmm, DYee, DYemu
      if(it_rg2==0 && muons.size()<2) continue;
      if(it_rg2==1 && electrons.size()<2) continue; 
      if(it_rg2==2 && !(muons.size()==1&&electrons.size()==1)) continue; 

      weight = 1.;
      // weights for MC
      if(!IsDATA){
        //Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        //Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        //if(truth_lep1.PID() == 0) continue;
        //if(truth_lep2.PID() == 0) continue;

        // Select prompt only
        if(GetLeptonType(*leptons.at(0), gens)<=0) continue;
        if(GetLeptonType(*leptons.at(1), gens)<=0) continue;

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
//          weight *= muon_idsf*muon_isosf;
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
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int j=0; j<electrons.size(); j++){
//          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }
      }

      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      ZCand = *leptons.at(0) + *leptons.at(1);

      // Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(ZCand.M() > 10.)) continue;

      // Cutflow : m(ll) > 10 GeV
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) continue;

      // Cutflow : OS event
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
      //FillHist(regionsSM.at(it_rg2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      //FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      //FillHist(regionsSM.at(it_rg2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

      if(!(Nbjet_medium == 0)) continue;

      // Cutflow : No b jets
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_NoMediumBJet_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_NoMediumBJet_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
      //FillHist(regionsSM.at(it_rg2)+"/Number_Jets_NoMediumBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_NoMediumBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_NoMediumBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_NoMediumBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      //FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_NoMediumBJet_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_NoMediumBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_NoMediumBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_NoMediumBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_NoMediumBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      //FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_NoMediumBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      //FillHist(regionsSM.at(it_rg2)+"/MET_NoMediumBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      //FillHist(regionsSM.at(it_rg2)+"/MET2ST_NoMediumBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      if(Nbjet_loose == 0){
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_NoLooseBJet_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
        //FillHist(regionsSM.at(it_rg2)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        //FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        //FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        //FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
        //FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        //FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        //FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        //FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        //FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        //FillHist(regionsSM.at(it_rg2)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        //FillHist(regionsSM.at(it_rg2)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      }

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_No3rdLep_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_No3rdLep_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Jets_No3rdLep_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_No3rdLep_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_No3rdLep_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_No3rdLep_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_No3rdLep_"+IDsuffix, ZCand.M(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_No3rdLep_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_No3rdLep_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_No3rdLep_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_No3rdLep_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_No3rdLep_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/MET_No3rdLep_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET2ST_No3rdLep_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

        
    } //JH : DYmm, DYee, DYemu done


    //=====================================
    //==== WZ, ZG, WG control region
    //=====================================
    if(it_rg2>2 && it_rg2<6 && leptons.size()==3){ //JH : WZ, ZG, WG, 3 tight leptons
  
      weight = 1.; 
      // weights for MC 
      if(!IsDATA){
        //Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        //Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        //Gen truth_lep3 = GetGenMatchedLepton(*leptons.at(2), gens);
        //if(truth_lep1.PID() == 0) continue;
        //if(truth_lep2.PID() == 0) continue;
        //if(truth_lep3.PID() == 0) continue; //JH : require all 3 lepton to be gen-matched

        // Select prompt only
        if(GetLeptonType(*leptons.at(0), gens)<=0) continue;
        if(GetLeptonType(*leptons.at(1), gens)<=0) continue;
        if(GetLeptonType(*leptons.at(2), gens)<=0) continue;

        weight *= weight_norm_1invpb*trigger_lumi;
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
          weight *= muon_recosf*muon_idsf*muon_isosf;
        }

        for(unsigned int j=0; j<electrons.size(); j++){
//          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }
      }

      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      // Cutflow : 3 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue; ///JH : 4th lepton veto

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      int l1 = -999, l2 = -999, l3 = -999, l4 = -999, wlepWZ = -999, wlepWG = -999;
      // OSSF lepton pair, W-tagged lepton
      if(muons.size()==2 && muons.at(0).Charge()*muons.at(1).Charge()<0){ //JH : mme
        ZCand = muons.at(0) + muons.at(1);
        WtagLep = electrons.at(0);
        ZtagLep1 = muons.at(0);
        ZtagLep2 = muons.at(1);
        GammaCand = ZCand;
        GammaLep1 = ZtagLep1;
        GammaLep2 = ZtagLep2;
      }
      else if(electrons.size()==2 && electrons.at(0).Charge()*electrons.at(1).Charge()<0){ //JH : eem
        ZCand = electrons.at(0) + electrons.at(1);
        WtagLep = muons.at(0);
        ZtagLep1 = electrons.at(0);
        ZtagLep2 = electrons.at(1);
        GammaCand = ZCand;
        GammaLep1 = ZtagLep1;
        GammaLep2 = ZtagLep2;
      }
      else if(muons.size()==3 || electrons.size()==3){ //JH : mmm / eee
        if(fabs(leptons.at(0)->Charge() + leptons.at(1)->Charge() + leptons.at(2)->Charge()) == 1){ //JH : 1 OSSF

          // ZCand, GammaCand
          double tmpMassDiff = 1000000., tmpMass = 100000.; 
          for(int ilep1=0; ilep1<2; ilep1++){
            for(int ilep2=ilep1+1; ilep2<3; ilep2++){ //JH : for each pair (01, 02, 12)
              if(leptons.at(ilep1)->Charge()*leptons.at(ilep2)->Charge()>0) continue; //JH : skip same sign
              Ztemp = *leptons.at(ilep1) + *leptons.at(ilep2);
              // For WZ, ZG
              if(!(Ztemp.M() > 10.)) ossf_mass10++; //JH : count m(OSSF) < 10GeV
              if(fabs(Ztemp.M() - MZ) < tmpMassDiff){
                tmpMassDiff = fabs(Ztemp.M() - MZ);
                ZCand = Ztemp; l1 = ilep1; l2 = ilep2; //JH : l1, l2 are the closest to m(Z)
              }
              // For WG
              if(Ztemp.M() < tmpMass){
                tmpMass = Ztemp.M();
                GammaCand = Ztemp; l3 = ilep1; l4 = ilep2; //JH : l3, l4 are the smallest mass
              }
            }
          }

          ZtagLep1 = *leptons.at(l1);
          ZtagLep2 = *leptons.at(l2);
          GammaLep1 = *leptons.at(l3);
          GammaLep2 = *leptons.at(l4);

          // Set the lepton from W
          for(int ilep3=0; ilep3<3; ilep3++){
            if(fabs(ilep3-l1)>0 && fabs(ilep3-l2)>0) wlepWZ = ilep3; //JH : ilep3 != l1 nor l2
            if(fabs(ilep3-l3)>0 && fabs(ilep3-l4)>0) wlepWG = ilep3; //JH : ilep3 != l3 nor l4
          }
          if(it_rg2 < 5) WtagLep = *leptons.at(wlepWZ); //JH : WZ, ZG
          else WtagLep = *leptons.at(wlepWG); //JH : WG
        }
        else continue;
      } 
      else continue;

      TriLep = *leptons.at(0) + *leptons.at(1) + *leptons.at(2);
      Mt = MT(WtagLep, METv);
      Mt3l = MT(TriLep, METv);

      if(it_rg2 < 5){   // WZ, ZG control region
        if(!(ossf_mass10 == 0)) continue; //JH : m(OSSF) > 10GeV
      
        // Cutflow : m(ll) > 10 GeV
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(Nbjet_medium == 0)) continue;

        // Cutflow : No b jets
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

        if(it_rg2 == 3){ //JH : WZ
          if(!IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET > 50.)) continue;
          if(!(Mt > 20.)) continue;
          if(!(TriLep.M() > MZ + 15.)) continue;
        }
        if(it_rg2 == 4){ //JH : ZG
          if(IsOnZ(ZCand.M(), 15.)) continue;
          if(!(MET < 50.)) continue;
          if(!IsOnZ(TriLep.M(), 15.)) continue;
        }
      }
      else{   // WG control region
        if(!(GammaCand.M() < 4.)) continue;

        // Cutflow : m(ll) < 4 GeV
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(Nbjet_medium == 0)) continue;

        // Cutflow : No b jets
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

        if(!(MET > 30.)) continue;
        if(!(Mt3l > 30.)) continue;
      }

      // weights for MC
      if(!IsDATA){
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
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      // Histograms
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.); 
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
      FillHist(regionsSM.at(it_rg2)+"/GammaCand_Mass_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
      FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/GammaCand_Pt_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/GammaLep1_Pt_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/GammaLep2_Pt_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      FillHist(regionsSM.at(it_rg2)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      FillHist(regionsSM.at(it_rg2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.); 

      //if(Nbjet_loose == 0){
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
      //  FillHist(regionsSM.at(it_rg2)+"/TriLep_Mass_NoLooseBJet_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
      //  FillHist(regionsSM.at(it_rg2)+"/GammaCand_Mass_NoLooseBJet_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/GammaCand_Pt_NoLooseBJet_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/WtagLep_Pt_NoLooseBJet_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZtagLep1_Pt_NoLooseBJet_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZtagLep2_Pt_NoLooseBJet_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/GammaLep1_Pt_NoLooseBJet_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/GammaLep2_Pt_NoLooseBJet_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Mt_NoLooseBJet_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      //} //JH : TODO meaning of NoLooseBJet?

      for(unsigned int it_ch2=0; it_ch2<channels3L.size(); it_ch2++){ //JH : for each lepton channels
        if(it_ch2 == electrons.size()){
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/TriLep_Mass_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaCand_Mass_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaCand_Pt_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/WtagLep_Pt_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep1_Pt_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep2_Pt_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaLep1_Pt_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaLep2_Pt_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Pt_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Eta_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Mt_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

          //if(Nbjet_loose == 0){
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Mass_NoLooseBJet_"+IDsuffix, ZCand.M(), weight, 80, 50., 130.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/TriLep_Mass_NoLooseBJet_"+IDsuffix, TriLep.M(), weight, 80, 50., 130.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaCand_Mass_NoLooseBJet_"+IDsuffix, GammaCand.M(), weight, 50, 0., 5.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZCand_Pt_NoLooseBJet_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaCand_Pt_NoLooseBJet_"+IDsuffix, GammaCand.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/WtagLep_Pt_NoLooseBJet_"+IDsuffix, WtagLep.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep1_Pt_NoLooseBJet_"+IDsuffix, ZtagLep1.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/ZtagLep2_Pt_NoLooseBJet_"+IDsuffix, ZtagLep2.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaLep1_Pt_NoLooseBJet_"+IDsuffix, GammaLep1.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/GammaLep2_Pt_NoLooseBJet_"+IDsuffix, GammaLep2.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/Mt_NoLooseBJet_"+IDsuffix, Mt, weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels3L.at(it_ch2)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          //}
        }
      } 

    }


    //=====================================
    //==== ZZ control region
    //=====================================
    if(it_rg2==6 && leptons.size()==4){

      weight = 1.;
      // weights for MC
      if(!IsDATA){
        //Gen truth_lep1 = GetGenMatchedLepton(*leptons.at(0), gens);
        //Gen truth_lep2 = GetGenMatchedLepton(*leptons.at(1), gens);
        //Gen truth_lep3 = GetGenMatchedLepton(*leptons.at(2), gens);
        //Gen truth_lep4 = GetGenMatchedLepton(*leptons.at(3), gens);
        //if(truth_lep1.PID() == 0) continue;
        //if(truth_lep2.PID() == 0) continue;
        //if(truth_lep3.PID() == 0) continue;
        //if(truth_lep4.PID() == 0) continue;

        // Select prompt only
        if(GetLeptonType(*leptons.at(0), gens)<=0) continue;
        if(GetLeptonType(*leptons.at(1), gens)<=0) continue;

        weight *= weight_norm_1invpb*trigger_lumi;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons.size(); i++){
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
          weight *= muon_recosf*muon_idsf*muon_isosf;
//          muon_idsf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
//          muon_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        for(unsigned int j=0; j<electrons.size(); j++){
//          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
//          ele_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(j).scEta(), electrons.at(j).Pt(), 0);
          ele_recosf = mcCorr->ElectronReco_SF(electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          if(param.Electron_Tight_ID.Contains("HEEP")){
            ele_idsf   = mcCorr->ElectronID_SF("HEEP", electrons.at(j).scEta(), electrons.at(j).UncorrPt(), 0);
          }
          else ele_idsf = 1.;
          weight *= ele_recosf*ele_idsf;
        }
      }

      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      // Cutflow : 4 tight leptons (gen-matched, pT > trigger thresholds)
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 2.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);

      if(lepton_veto_size > 0) continue;

      // Cutflow : veto additional leptons using veto ID 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

      // OSSF lepton pairs
      if(muons.size()==2 && muons.at(0).Charge()*muons.at(1).Charge()<0 && electrons.at(0).Charge()*electrons.at(1).Charge()<0){ //JH : 2 OSSF
        ZCand1 = muons.at(0) + muons.at(1);
        ZCand2 = electrons.at(0) + electrons.at(1);
      }
      else if(muons.size()==4 || electrons.size()==4){
        if(leptons_minus.size() == leptons_plus.size()){ //JH : 2 OSSF
          Ztemp1 = *leptons_minus.at(0) + *leptons_plus.at(0);
          Ztemp2 = *leptons_minus.at(1) + *leptons_plus.at(1);
          Ztemp3 = *leptons_minus.at(0) + *leptons_plus.at(1);
          Ztemp4 = *leptons_minus.at(1) + *leptons_plus.at(0);
          if(!(Ztemp1.M()>10. && Ztemp2.M()>10. && Ztemp3.M()>10. && Ztemp4.M()>10.)) ossf_mass10++; //JH : all combination of m(OSSF) > 10GeV
          ZCand1 = Ztemp1; ZCand2 = Ztemp2;

          if(!(IsOnZ(ZCand1.M(), 15.) && IsOnZ(ZCand2.M(), 15.))){ //JH : 1st OSSF combination must be IsOnZ
            ZCand1 = Ztemp3; ZCand2 = Ztemp4; //JH : 2nd OSSF combination
          }
        }
      }
      else continue;

      if(!(ossf_mass10 == 0)) continue;

      // Cutflow : m(ll) > 10 GeV
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(Nbjet_medium == 0)) continue;

      // Cutflow : No b jets
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

      if(!(IsOnZ(ZCand1.M(), 15.) && IsOnZ(ZCand2.M(), 15.))) continue; //JH : 2nd OSSF combination also must be IsOnZ

      // weights for MC
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
      if(RunFake) weight = fakeEst->GetWeight(leptons, param);

      // Histograms 
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(regionsSM.at(it_rg2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
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

      //if(Nbjet_loose == 0){
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand1_Mass_NoLooseBJet_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand2_Mass_NoLooseBJet_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand1_Pt_NoLooseBJet_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/ZCand2_Pt_NoLooseBJet_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep4_Pt_NoLooseBJet_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/Lep4_Eta_NoLooseBJet_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
      //  FillHist(regionsSM.at(it_rg2)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
      //  FillHist(regionsSM.at(it_rg2)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      //}

      for(unsigned int it_ch2=0; it_ch2<channels4L.size(); it_ch2++){ //JH : for each lepton channels
        if(it_ch2 == electrons.size()/2){
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
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

          //if(Nbjet_loose == 0){
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_NoLooseBJet_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Events_unweighted_NoLooseBJet_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_Jets_NoLooseBJet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Loose_NoLooseBJet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_BJets_Medium_NoLooseBJet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Number_FatJets_NoLooseBJet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand1_Mass_NoLooseBJet_"+IDsuffix, ZCand1.M(), weight, 80, 50., 130.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand2_Mass_NoLooseBJet_"+IDsuffix, ZCand2.M(), weight, 80, 50., 130.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand1_Pt_NoLooseBJet_"+IDsuffix, ZCand1.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/ZCand2_Pt_NoLooseBJet_"+IDsuffix, ZCand2.Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep1_Pt_NoLooseBJet_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep2_Pt_NoLooseBJet_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep3_Pt_NoLooseBJet_"+IDsuffix, leptons.at(2)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep4_Pt_NoLooseBJet_"+IDsuffix, leptons.at(3)->Pt(), weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep1_Eta_NoLooseBJet_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep2_Eta_NoLooseBJet_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep3_Eta_NoLooseBJet_"+IDsuffix, leptons.at(2)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/Lep4_Eta_NoLooseBJet_"+IDsuffix, leptons.at(3)->Eta(), weight, 50, -2.5, 2.5);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/MET_NoLooseBJet_"+IDsuffix, MET, weight, 1000, 0., 1000.);
          //  FillHist(regionsSM.at(it_rg2)+"/"+channels4L.at(it_ch2)+"/MET2ST_NoLooseBJet_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
          //}
        } 
      }
      
    }
  } */
}




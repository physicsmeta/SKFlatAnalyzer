#include "Practice.h"

Practice::Practice(){

}

void Practice::initializeAnalyzer(){

  if(DataYear==2016){
    EleIDs = { "passMediumID" }; // PassID() in Electron.C
		EleIDSFKeys = {	"passMediumID" }; // histmap.txt
    EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"; // referred to HS's SKFlatValidation.C
    lep0ptcut = 25.;
		lep1ptcut = 15.;
  }
//	else if(DataYear==2017){
//		EleIDs = {
//			"passMediumID",
//			"passTightID",
//		};
//		EleIDSFKeys = {
//			"passMediumID",
//			"passTightID",
//		};
//		EleTriggerName = "Ele35_WPTight_Gsf";
//		TriggerSafePtCut = 38.;
//	}

}

Practice::~Practice(){

}

void Practice::executeEvent(){

  AllEles = GetAllElectrons();

	weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

	for(unsigned int i=0; i<EleIDs.size(); i++){

		TString EleID = EleIDs.at(i);
		TString EleIDSFKey = EleIDSFKeys.at(i);

		param.Clear();

		param.syst_ = AnalyzerParameter::Central;

		param.Name = EleID+"_"+"syst_Central";

		param.Electron_Tight_ID = EleID;
		param.Electron_ID_SF_Key = EleIDSFKey;

  executeEventFromParameter(param);

  }

}

void Practice::executeEvent(Long64_t Nentry){

  AllEles = GetAllElectrons();

	weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

	for(unsigned int i=0; i<EleIDs.size(); i++){

		TString EleID = EleIDs.at(i);
		TString EleIDSFKey = EleIDSFKeys.at(i);

		param.Clear();

		param.syst_ = AnalyzerParameter::Central;

		param.Name = EleID+"_"+"syst_Central";

		param.Electron_Tight_ID = EleID;
		param.Electron_ID_SF_Key = EleIDSFKey;

  executeEventFromParameter(param, Nentry);

  }

}

void Practice::executeEventFromParameter(AnalyzerParameter param){

  JSFillHist("Practice", "Cutflow", 0, 1, 7, 0, 7);

  /* MET Filter */

  if(!PassMETFilter()) return;

  JSFillHist("Practice", "Cutflow", 1, 1, 7, 0, 7);

  Event ev = GetEvent();
//Particle METv = ev.GetMETVector(); 

	/* Trigger */

	if(! (ev.PassTrigger(EleTriggerName) )) return;

  JSFillHist("Practice", "Cutflow", 2, 1, 7, 0, 7);

	/* copy */

  vector<Electron> this_AllEles = AllEles;

	/* ID selection */

  vector<Electron> eles = SelectElectrons(this_AllEles, param.Electron_Tight_ID, 25., 2.5);

  /* sort */

	std::sort(eles.begin(), eles.end(), PtComparing);

	/* Event selection */

  JSFillHist("Practice", "Cutflow", 3, 1, 7, 0, 7);

	if(eles.size() != 2) return;
	
  JSFillHist("Practice", "Cutflow", 4, 1, 7, 0, 7);

	if( eles.at(0).Pt() < lep0ptcut || eles.at(1).Pt() < lep1ptcut ) return; 

  JSFillHist("Practice", "Cutflow", 5, 1, 7, 0, 7);

	Particle ZCand = eles.at(0) + eles.at(1);
	if( ZCand.M() < 60. ) return;

  JSFillHist("Practice", "Cutflow", 6, 1, 7, 0, 7);

	double weight = 1.;

  if(!IsDATA){

		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= weight_Prefire;

		for(unsigned int i=0; i<eles.size(); i++){
	
			double this_recosf = mcCorr->ElectronReco_SF(eles.at(i).scEta(), eles.at(i).Pt());
			double this_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, eles.at(i).scEta(), eles.at(i).Pt());

			weight *= this_recosf*this_idsf;

		}

	  vector<Gen> gens = GetGens();
    Gen truth_lep = GetGenMatchedLepton(eles.at(0), gens);
		if(truth_lep.PID() == 0) return; // TODO check the meaning

		int truth_lep_Charge;
		if(truth_lep.PID() == 11) truth_lep_Charge = -1;
		else if(truth_lep.PID() == -11) truth_lep_Charge = 1;

    if(MCSample.Contains("DYJets")){
			if(abs(eles.at(0).scEta())<0.8){
		    JSFillHist("ChargeFlip", "EtaRegion1_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion1_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
		  }
			else if(0.8<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<1.4442){
		    JSFillHist("ChargeFlip", "EtaRegion2_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion2_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			}
			else if(1.556<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<2.5){
		    JSFillHist("ChargeFlip", "EtaRegion3_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion3_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			} // TODO Don't need this. delete later.
		}
	}

/* Draw ZMass (weighted if MC) */

  JSFillHist("ZMass", "Total", ZCand.M(), weight, 40, 70., 110.);

  if(eles.at(0).Charge()*eles.at(1).Charge()<0){
		JSFillHist("ZMass", "OS", ZCand.M(), weight, 40, 70., 110.);
	}
	else{
		JSFillHist("ZMass", "SS", ZCand.M(), weight, 40, 70., 110.);
	}


/* See the energy, momentum */

  JSFillHist("Practice", "ele1_pt", eles.at(0).Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_pt", eles.at(1).Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_pt", ZCand.Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_EsinTheta", eles.at(0).E()*TMath::Sin(eles.at(0).Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_EsinTheta", eles.at(1).E()*TMath::Sin(eles.at(1).Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_EsinTheta", ZCand.E()*TMath::Sin(ZCand.Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_px", eles.at(0).Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_px", eles.at(1).Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_px", ZCand.Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_py", eles.at(0).Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_py", eles.at(1).Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_py", ZCand.Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_pz", eles.at(0).Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_pz", eles.at(1).Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_pz", ZCand.Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_m", eles.at(0).M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_m", eles.at(1).M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_m", ZCand.M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_sqrt(E^2-P^2)", sqrt(pow(eles.at(0).E(),2)-pow(eles.at(0).P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_sqrt(E^2-P^2)", sqrt(pow(eles.at(1).E(),2)-pow(eles.at(1).P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_sqrt(E^2-P^2)", sqrt(pow(ZCand.E(),2)-pow(ZCand.P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_E", eles.at(0).E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_E", eles.at(1).E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_E", ZCand.E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_P", eles.at(0).P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_P", eles.at(1).P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_P", ZCand.P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_sqrt(px^2+py^2+pz^2)", sqrt(pow(eles.at(0).Px(),2)+pow(eles.at(0).Py(),2)+pow(eles.at(0).Pz(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_sqrt(px^2+py^2+pz^2)", sqrt(pow(eles.at(1).Px(),2)+pow(eles.at(1).Py(),2)+pow(eles.at(1).Pz(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_sqrt(px^2+py^2+pz^2)", sqrt(pow(ZCand.Px(),2)+pow(ZCand.Py(),2)+pow(ZCand.Pz(),2)), weight, 200, 0., 200.);
		
}


void Practice::executeEventFromParameter(AnalyzerParameter param, Long64_t Nentry){

  JSFillHist("Practice", "Cutflow", 0, 1, 7, 0, 7);

  /* MET Filter */

  if(!PassMETFilter()) return;

  JSFillHist("Practice", "Cutflow", 1, 1, 7, 0, 7);

  Event ev = GetEvent();
//Particle METv = ev.GetMETVector(); 

	/* Trigger */

	if(! (ev.PassTrigger(EleTriggerName) )) return;

  JSFillHist("Practice", "Cutflow", 2, 1, 7, 0, 7);

	/* copy */

  vector<Electron> this_AllEles = AllEles;

	/* ID selection */

  vector<Electron> eles = SelectElectrons(this_AllEles, param.Electron_Tight_ID, 25., 2.5);

  /* sort */

	std::sort(eles.begin(), eles.end(), PtComparing);

	/* Event selection */

  JSFillHist("Practice", "Cutflow", 3, 1, 7, 0, 7);

	if(eles.size() != 2) return;
	
  JSFillHist("Practice", "Cutflow", 4, 1, 7, 0, 7);

	if( eles.at(0).Pt() < lep0ptcut || eles.at(1).Pt() < lep1ptcut ) return; 

  JSFillHist("Practice", "Cutflow", 5, 1, 7, 0, 7);

	Particle ZCand = eles.at(0) + eles.at(1);
	if( ZCand.M() < 60. ) return;

  JSFillHist("Practice", "Cutflow", 6, 1, 7, 0, 7);

	double weight = 1.;

  if(!IsDATA){

		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= weight_Prefire;

		for(unsigned int i=0; i<eles.size(); i++){
	
			double this_recosf = mcCorr->ElectronReco_SF(eles.at(i).scEta(), eles.at(i).Pt());
			double this_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, eles.at(i).scEta(), eles.at(i).Pt());

			weight *= this_recosf*this_idsf;

		}

	  vector<Gen> gens = GetGens();
    Gen truth_lep = GetGenMatchedLepton(eles.at(0), gens);
		if(truth_lep.PID() == 0) return; // TODO check the meaning

		int truth_lep_Charge;
		if(truth_lep.PID() == 11) truth_lep_Charge = -1;
		else if(truth_lep.PID() == -11) truth_lep_Charge = 1;

    if(MCSample.Contains("DYJets")){
			if(abs(eles.at(0).scEta())<0.8){
		    JSFillHist("ChargeFlip", "EtaRegion1_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "In loop " << Nentry << ", !!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion1_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
		  }
			else if(0.8<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<1.4442){
		    JSFillHist("ChargeFlip", "EtaRegion2_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "In loop " << Nentry << ", !!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion2_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			}
			else if(1.556<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<2.5){
		    JSFillHist("ChargeFlip", "EtaRegion3_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "In loop " << Nentry << ", !!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist("ChargeFlip", "EtaRegion3_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			}
		}
    
    if(Nentry%(LogEvery)==0){
      PrintGen(gens);
			cout << "===========================================================" << endl;
    }

	}

/* Draw ZMass (weighted if MC) */

  JSFillHist("ZMass", "Total", ZCand.M(), weight, 40, 70., 110.);

  if(eles.at(0).Charge()*eles.at(1).Charge()<0){
		JSFillHist("ZMass", "OS", ZCand.M(), weight, 40, 70., 110.);
	}
	else{
		JSFillHist("ZMass", "SS", ZCand.M(), weight, 40, 70., 110.);
	}


/* See the energy, momentum */

  JSFillHist("Practice", "ele1_pt", eles.at(0).Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_pt", eles.at(1).Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_pt", ZCand.Pt(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_EsinTheta", eles.at(0).E()*TMath::Sin(eles.at(0).Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_EsinTheta", eles.at(1).E()*TMath::Sin(eles.at(1).Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_EsinTheta", ZCand.E()*TMath::Sin(ZCand.Theta()), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_px", eles.at(0).Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_px", eles.at(1).Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_px", ZCand.Px(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_py", eles.at(0).Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_py", eles.at(1).Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_py", ZCand.Py(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_pz", eles.at(0).Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_pz", eles.at(1).Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_pz", ZCand.Pz(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_m", eles.at(0).M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_m", eles.at(1).M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_m", ZCand.M(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_sqrt(E^2-P^2)", sqrt(pow(eles.at(0).E(),2)-pow(eles.at(0).P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_sqrt(E^2-P^2)", sqrt(pow(eles.at(1).E(),2)-pow(eles.at(1).P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_sqrt(E^2-P^2)", sqrt(pow(ZCand.E(),2)-pow(ZCand.P(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_E", eles.at(0).E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_E", eles.at(1).E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_E", ZCand.E(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_P", eles.at(0).P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_P", eles.at(1).P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_P", ZCand.P(), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele1_sqrt(px^2+py^2+pz^2)", sqrt(pow(eles.at(0).Px(),2)+pow(eles.at(0).Py(),2)+pow(eles.at(0).Pz(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "ele2_sqrt(px^2+py^2+pz^2)", sqrt(pow(eles.at(1).Px(),2)+pow(eles.at(1).Py(),2)+pow(eles.at(1).Pz(),2)), weight, 200, 0., 200.);
  JSFillHist("Practice", "diele_sqrt(px^2+py^2+pz^2)", sqrt(pow(ZCand.Px(),2)+pow(ZCand.Py(),2)+pow(ZCand.Pz(),2)), weight, 200, 0., 200.);

  if(Nentry%(LogEvery)==0){
    cout << "ele1_E = " << eles.at(0).E() << endl;
    cout << "ele1_P = " << eles.at(0).P() << endl;
    cout << "ele1_Px = " << eles.at(0).Px() << endl;
    cout << "ele1_Py = " << eles.at(0).Py() << endl;
    cout << "ele1_Pz = " << eles.at(0).Pz() << endl;
    cout << "ele1_sqrt(E^2-P^2) = " << sqrt(pow(eles.at(0).E(),2)-pow(eles.at(0).P(),2)) << endl;
    cout << "ele1_m = " << eles.at(0).M() << endl;
    cout << "ele2_E = " << eles.at(1).E() << endl;
    cout << "ele2_P = " << eles.at(1).P() << endl;
    cout << "ele2_Px = " << eles.at(1).Px() << endl;
    cout << "ele2_Py = " << eles.at(1).Py() << endl;
    cout << "ele2_Pz = " << eles.at(1).Pz() << endl;
    cout << "ele2_sqrt(E^2-P^2) = " << sqrt(pow(eles.at(1).E(),2)-pow(eles.at(1).P(),2)) << endl;
    cout << "ele2_m = " << eles.at(1).M() << endl;
    cout << "diele_E = " << ZCand.E() << endl;
    cout << "diele_P = " << ZCand.P() << endl;
    cout << "diele_Px = " << ZCand.Px() << endl;
    cout << "diele_Py = " << ZCand.Py() << endl;
    cout << "diele_Pz = " << ZCand.Pz() << endl;
    cout << "diele_sqrt(E^2-P^2) = " << sqrt(pow(ZCand.E(),2)-pow(ZCand.P(),2)) << endl;
    cout << "diele_m = " << ZCand.M() << endl;
  }

}




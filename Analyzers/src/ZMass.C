#include "ZMass.h"

ZMass::ZMass(){

}

void ZMass::initializeAnalyzer(){

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

ZMass::~ZMass(){

}

void ZMass::executeEvent(){

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

void ZMass::executeEventFromParameter(AnalyzerParameter param){

  /* MET Filter */

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
//Particle METv = ev.GetMETVector(); 

	/* Trigger */

	if(! (ev.PassTrigger(EleTriggerName) )) return;

	/* copy */

  vector<Electron> this_AllEles = AllEles;

	/* ID selection */

  vector<Electron> eles = SelectElectrons(this_AllEles, param.Electron_Tight_ID, 25., 2.5);

  /* sort */

	std::sort(eles.begin(), eles.end(), PtComparing);

	/* Event selection */

	if(eles.size() != 2) return;

	if( eles.at(0).Pt() < lep0ptcut || eles.at(1).Pt() < lep1ptcut ) return; 

	Particle ZCand = eles.at(0) + eles.at(1);
	if( ZCand.M() < 60. ) return;

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
		    JSFillHist(param.Name, "ChargeFlip_EtaRegion1_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist(param.Name, "ChargeFlip_EtaRegion1_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
		  }
			else if(0.8<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<1.4442){
		    JSFillHist(param.Name, "ChargeFlip_EtaRegion2_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist(param.Name, "ChargeFlip_EtaRegion2_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			}
			else if(1.556<=abs(eles.at(0).scEta())&&abs(eles.at(0).scEta())<2.5){
		    JSFillHist(param.Name, "ChargeFlip_EtaRegion3_Denom_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    if(truth_lep_Charge*eles.at(0).Charge()<0){
			    cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(0).Charge() << endl;
				  JSFillHist(param.Name, "ChargeFlip_EtaRegion3_Num_"+param.Name, 1/eles.at(0).Pt(), 1., 40, 0., 0.04);
		    }
			} // TODO Don't need this. delete later.
		}
	}

/* Fill histo */

  JSFillHist(param.Name, "ZCand_Mass_"+param.Name, ZCand.M(), weight, 40, 70., 110.);

  if(eles.at(0).Charge()*eles.at(1).Charge()<0){
		JSFillHist(param.Name, "ZCand_Mass_"+param.Name+"_OS", ZCand.M(), weight, 40, 70., 110.);
	}
	else{
		JSFillHist(param.Name, "ZCand_Mass_"+param.Name+"_SS", ZCand.M(), weight, 40, 70., 110.);
	}

}




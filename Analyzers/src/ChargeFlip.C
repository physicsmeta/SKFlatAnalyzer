#include "ChargeFlip.h" // For now, 2016 DYJets only

ChargeFlip::ChargeFlip(){

}

ChargeFlip::~ChargeFlip(){

}

void ChargeFlip::executeEvent(){

	AllEles = GetAllElectrons();

	/* MET Filter */

	if(HasFlag("passMETFilter")){ 

		if(!PassMETFilter()) return;

	}

  Event ev = GetEvent();

	/* Trigger */

	if(HasFlag("passSingleTrigger")) {

		TString EleTriggerName = "HLT_Ele27_WPTight_Gsf_v";

		if(! (ev.PassTrigger(EleTriggerName) )) return;

	}

	else if(HasFlag("passDoubleTrigger")){

		TString EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

		if(! (ev.PassTrigger(EleTriggerName) )) return;

	}

	/* ChargeFlip ID selection */

  vector<Electron> eles;

  if(!HasFlag("passTightID")) eles = SelectChargeFlipElectrons(AllEles, 25., 2.5);

	else if(HasFlag("passTightChargeTightID")) eles = SelectChargeTightElectrons(AllEles, "passTightID", 25., 2.5);

  /* sort */

	std::sort(eles.begin(), eles.end(), PtComparing);

  if(IsDATA){
		cout << "[ChargeFlip::executeEvent()] Put MC samples only!!" << endl;
		exit(EXIT_FAILURE);
	}

  vector<Gen> gens = GetGens();
	
	for(unsigned int i=0; i<eles.size(); i++){
    Gen truth_lep = GetGenMatchedLepton(eles.at(i), gens);
  	if(truth_lep.PID() == 0) return; // TODO check the meaning
  
  	int truth_lep_Charge;
  	if(truth_lep.PID() == 11) truth_lep_Charge = -1;
  	else if(truth_lep.PID() == -11) truth_lep_Charge = 1;
  
    if(MCSample.Contains("DYJets")){
  		if(abs(eles.at(i).scEta())<0.8){
  	    JSFillHist("ChargeFlip", "EtaRegion1_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    if(truth_lep_Charge*eles.at(i).Charge()<0){
  		    cout << "!!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
  			  JSFillHist("ChargeFlip", "EtaRegion1_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    }
  	  }
  		else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
  	    JSFillHist("ChargeFlip", "EtaRegion2_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    if(truth_lep_Charge*eles.at(i).Charge()<0){
  		    cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
  			  JSFillHist("ChargeFlip", "EtaRegion2_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    }
  		}
  		else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
  	    JSFillHist("ChargeFlip", "EtaRegion3_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    if(truth_lep_Charge*eles.at(i).Charge()<0){
  		    cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
  			  JSFillHist("ChargeFlip", "EtaRegion3_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
  	    }
  		}
  	}
		else if(!MCSample.Contains("DYJets")){
		  cout << "[ChargeFlip::executeEvent()] Put DY samples only!!" << endl;
		  exit(EXIT_FAILURE);
							  }

	}
}




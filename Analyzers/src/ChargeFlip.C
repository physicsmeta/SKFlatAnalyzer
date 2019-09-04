#include "ChargeFlip.h" // For now, 2016 DYJets only

ChargeFlip::ChargeFlip(){

}

void ChargeFlip::initializeAnalyzer(){

  if(DataYear==2016){
    EleIDs = { "passMediumID" }; // PassID() in Electron.C
    EleIDSFKeys = { "passMediumID" }; // histmap.txt
    EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"; // referred to HS's SKFlatValidation.C
    lep0ptcut = 25.;
    lep1ptcut = 15.;
  }
//else if(DataYear==2017){
//    EleIDs = {
//      "passMediumID",
//      "passTightID",
//    };
//    EleIDSFKeys = {
//      "passMediumID",
//      "passTightID",
//    };
//    EleTriggerName = "Ele35_WPTight_Gsf";
//    TriggerSafePtCut = 38.;
//  }

}

ChargeFlip::~ChargeFlip(){

}

void ChargeFlip::executeEvent(Long64_t Nentry){

  AllEles = GetAllElectrons();

  Event ev = GetEvent();

  vector<Electron> eles;

  /* Measure CF rate using DY MC */

  if(!IsDATA){

    /* CF ID selection */

    if(HasFlag("ChargeFlipID")) eles = SelectChargeFlipElectrons(AllEles, 25., 2.5);
    else if(HasFlag("passTightChargeTightID")) eles = SelectChargeTightElectrons(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdXY")) eles = SelectChargeTightElectronsDXY(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdZ")) eles = SelectChargeTightElectronsDZ(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdXYdZ")) eles = SelectChargeTightElectronsDXYDZ(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("ChargeFlipHE")) eles = SelectElectronsHE(AllEles, "passTightID", 25., 2.5);
    std::sort(eles.begin(), eles.end(), PtComparing);

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
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            JSFillHist("ChargeFlip", "EtaRegion1_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            JSFillHist("ChargeFlip", "EtaRegion1_IsLeading", 0, 1, 2, 0, 2);
            }
						//
          }
        }
        else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
          JSFillHist("ChargeFlip", "EtaRegion2_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            JSFillHist("ChargeFlip", "EtaRegion2_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            JSFillHist("ChargeFlip", "EtaRegion2_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            JSFillHist("ChargeFlip", "EtaRegion2_IsLeading", 0, 1, 2, 0, 2);
            }
						//
          }
        }
        else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
          JSFillHist("ChargeFlip", "EtaRegion3_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            JSFillHist("ChargeFlip", "EtaRegion3_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            JSFillHist("ChargeFlip", "EtaRegion3_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            JSFillHist("ChargeFlip", "EtaRegion3_IsLeading", 0, 1, 2, 0, 2);
            }
						//
          }
        }
      }
      else if(!MCSample.Contains("DYJets")){
        cout << "[ChargeFlip::executeEvent()] Put DY samples only!!" << endl;
        exit(EXIT_FAILURE);
      }
  	}

    /* MC Closure test start */
    
    if(eles.size() != 2) return;
    if(!PassMETFilter()) return;
    Particle METv = ev.GetMETVector();

    Particle ZCand = eles.at(0)+eles.at(1);
    if(! (70.<=ZCand.M()&&ZCand.M()<110.) ) return;
    double weight = GetCFweight(eles);

    if(Nentry%(LogEvery*10)==0){
    cout << "OS event weighted by: " << weight << endl;
    }
  
    JSFillHist("ClosureTest", "ZMass_total", ZCand.M(), 1., 40, 70., 110.);
  
    if(eles.at(0).Charge()*eles.at(1).Charge()<0){
      JSFillHist("ClosureTest", "ZMass_OS_Unweighted", ZCand.M(), 1., 40, 70., 110.);
      JSFillHist("ClosureTest", "ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
      JSFillHist("ClosureTest", "pt1_OS_CFweighted", eles.at(0).Pt(), weight, 70, 20., 90.);
      JSFillHist("ClosureTest", "pt2_OS_CFweighted", eles.at(1).Pt(), weight, 70, 20., 90.);
      JSFillHist("ClosureTest", "MET_OS_CFweighted", METv.Pt(), weight, 100, 0., 100.);
    }
    else{
      JSFillHist("ClosureTest", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
      JSFillHist("ClosureTest", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
      JSFillHist("ClosureTest", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
      JSFillHist("ClosureTest", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
    }

    /* There's disagreement between OS_CFweighted and SS. */
		/* Now, Let's shift the electrons' energy */
    
    Particle ZCand_tmp;
    Particle METv_tmp;
    double weight_tmp;

    for(int i=0;i<50;i++){
      vector<Electron> eles_tmp = eles; // copy the vector
      for(int j=0;j<2;j++){
        eles_tmp.at(j).SetE(eles_tmp.at(j).E()*(1-0.001*(i+1)));
        eles_tmp.at(j).SetPtEtaPhiE(eles_tmp.at(j).E() * TMath::Sin(eles_tmp.at(j).Theta()), eles_tmp.at(j).Eta(), eles_tmp.at(j).Phi(), eles_tmp.at(j).E());
      }

			ZCand_tmp = eles_tmp.at(0) + eles_tmp.at(1);
      METv_tmp.SetPxPyPzE(METv.Px()-ZCand.Px()+ZCand_tmp.Px(),METv.Py()-ZCand.Py()+ZCand_tmp.Py(),0,METv.E()-ZCand.E()+ZCand_tmp.E());
      weight_tmp = GetCFweight(eles_tmp);

      if(! (70.<=ZCand_tmp.M()&&ZCand_tmp.M()<110.) ) continue;

      if(eles.at(0).Charge()*eles.at(1).Charge()<0){
        JSFillHist("ClosureTest", "ZMass_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
        JSFillHist("ClosureTest", "pt1_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), eles_tmp.at(0).Pt(), weight_tmp, 70, 20., 90.);
        JSFillHist("ClosureTest", "pt2_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), eles_tmp.at(1).Pt(), weight_tmp, 70, 20., 90.);
        JSFillHist("ClosureTest", "MET_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), METv_tmp.Pt(), weight_tmp, 100, 0., 100.);
      }
		}
  }

  /* Scale Factor measurement */

	else{
    
    if(!PassMETFilter()) return;
    if(! (ev.PassTrigger(EleTriggerName) )) return;
    
    /* CF SF ID selection */

    if(HasFlag("ChargeFlipID")) eles = SelectChargeFlipElectrons(AllEles, 25., 2.5);
    else if(HasFlag("passTightChargeTightID")) eles = SelectChargeTightElectrons(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdXY")) eles = SelectChargeTightElectronsDXY(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdZ")) eles = SelectChargeTightElectronsDZ(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("passTightChargeTightIDdXYdZ")) eles = SelectChargeTightElectronsDXYDZ(AllEles, "passTightID", 25., 2.5);
    else if(HasFlag("ChargeFlipHE")) eles = SelectElectronsHE(AllEles, "passTightID", 25., 2.5);
    std::sort(eles.begin(), eles.end(), PtComparing);

    //

    if(eles.size() != 2) return;
		//if(eles.at(0).Pt()<lep0ptcut||eles.at(1).Pt()<lep1ptcut) return; //No need already pt min = 25

    Particle METv = ev.GetMETVector();
    Particle ZCand = eles.at(0)+eles.at(1);
    if(! (70.<=ZCand.M()&&ZCand.M()<110.) ) return;
    double weight = GetCFweight(eles);
    double weight_SF = GetCFweight_SF(eles);

    if(Nentry%(LogEvery*10)==0){
    cout << "OS event weighted by: " << weight << endl;
    }
  
    JSFillHist("ScaleFactor", "ZMass_total", ZCand.M(), 1., 40, 70., 110.);
  
    if(eles.at(0).Charge()*eles.at(1).Charge()<0){
      JSFillHist("ScaleFactor", "ZMass_OS_Unweighted", ZCand.M(), 1., 40, 70., 110.);
      JSFillHist("ScaleFactor", "ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
      JSFillHist("ScaleFactor", "pt1_OS_CFweighted", eles.at(0).Pt(), weight, 70, 20., 90.);
      JSFillHist("ScaleFactor", "pt2_OS_CFweighted", eles.at(1).Pt(), weight, 70, 20., 90.);
      JSFillHist("ScaleFactor", "MET_OS_CFweighted", METv.Pt(), weight, 100, 0., 100.);
    }
    else{
      JSFillHist("ScaleFactor", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
      JSFillHist("ScaleFactor", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
      JSFillHist("ScaleFactor", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
      JSFillHist("ScaleFactor", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
    }

    // BB
    if(abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())<1.4442){
    JSFillHist("ScaleFactor/BB", "ZMass_total", ZCand.M(), 1., 40, 70., 110.);

      if(eles.at(0).Charge()*eles.at(1).Charge()<0){
        JSFillHist("ScaleFactor/BB", "ZMass_OS_Unweighted", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/BB", "ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
        JSFillHist("ScaleFactor/BB", "pt1_OS_CFweighted", eles.at(0).Pt(), weight, 70, 20., 90.);
        JSFillHist("ScaleFactor/BB", "pt2_OS_CFweighted", eles.at(1).Pt(), weight, 70, 20., 90.);
        JSFillHist("ScaleFactor/BB", "MET_OS_CFweighted", METv.Pt(), weight, 100, 0., 100.);
      }
      else{
        JSFillHist("ScaleFactor/BB", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/BB", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/BB", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/BB", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
      }
    }

    // BE
    if((abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())>=1.556)||(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())<1.4442)){
    JSFillHist("ScaleFactor/BE", "ZMass_total", ZCand.M(), 1., 40, 70., 110.);

      if(eles.at(0).Charge()*eles.at(1).Charge()<0){
        JSFillHist("ScaleFactor/BE", "ZMass_OS_Unweighted", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/BE", "ZMass_OS_CFweighted", ZCand.M(), weight_SF, 40, 70., 110.);
        JSFillHist("ScaleFactor/BE", "pt1_OS_CFweighted", eles.at(0).Pt(), weight_SF, 70, 20., 90.);
        JSFillHist("ScaleFactor/BE", "pt2_OS_CFweighted", eles.at(1).Pt(), weight_SF, 70, 20., 90.);
        JSFillHist("ScaleFactor/BE", "MET_OS_CFweighted", METv.Pt(), weight_SF, 100, 0., 100.);
      }
      else{
        JSFillHist("ScaleFactor/BE", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/BE", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/BE", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/BE", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
      }
    }

    // EE
    if(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())>=1.556){
    JSFillHist("ScaleFactor/EE", "ZMass_total", ZCand.M(), 1., 40, 70., 110.);

      if(eles.at(0).Charge()*eles.at(1).Charge()<0){
        JSFillHist("ScaleFactor/EE", "ZMass_OS_Unweighted", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/EE", "ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
        JSFillHist("ScaleFactor/EE", "pt1_OS_CFweighted", eles.at(0).Pt(), weight, 70, 20., 90.);
        JSFillHist("ScaleFactor/EE", "pt2_OS_CFweighted", eles.at(1).Pt(), weight, 70, 20., 90.);
        JSFillHist("ScaleFactor/EE", "MET_OS_CFweighted", METv.Pt(), weight, 100, 0., 100.);
      }
      else{
        JSFillHist("ScaleFactor/EE", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
        JSFillHist("ScaleFactor/EE", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/EE", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
        JSFillHist("ScaleFactor/EE", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
      }
    }

		/* Now let's shift the electrons' energy 1.3% */
    
    Particle ZCand_tmp;
    Particle METv_tmp;
    double weight_tmp, weight_tmp_SF;

    vector<Electron> eles_tmp = eles; // copy the vector
    for(int j=0;j<2;j++){
      eles_tmp.at(j).SetE(eles_tmp.at(j).E()*(1-0.013));
      eles_tmp.at(j).SetPtEtaPhiE(eles_tmp.at(j).E() * TMath::Sin(eles_tmp.at(j).Theta()), eles_tmp.at(j).Eta(), eles_tmp.at(j).Phi(), eles_tmp.at(j).E());
    }

		ZCand_tmp = eles_tmp.at(0) + eles_tmp.at(1);
    METv_tmp.SetPxPyPzE(METv.Px()-ZCand.Px()+ZCand_tmp.Px(),METv.Py()-ZCand.Py()+ZCand_tmp.Py(),0,METv.E()-ZCand.E()+ZCand_tmp.E());
    weight_tmp = GetCFweight(eles_tmp);
    weight_tmp_SF = GetCFweight_SF(eles_tmp);

		if(70.<=ZCand_tmp.M()&&ZCand_tmp.M()<110.){
      if(eles.at(0).Charge()*eles.at(1).Charge()<0){
        JSFillHist("ScaleFactor", "ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
        JSFillHist("ScaleFactor", "pt1_OS_CFweighted_shifted_1.3%", eles_tmp.at(0).Pt(), weight_tmp, 70, 20., 90.);
        JSFillHist("ScaleFactor", "pt2_OS_CFweighted_shifted_1.3%", eles_tmp.at(1).Pt(), weight_tmp, 70, 20., 90.);
        JSFillHist("ScaleFactor", "MET_OS_CFweighted_shifted_1.3%", METv_tmp.Pt(), weight_tmp, 100, 0., 100.);

        // BB
        if(abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())<1.4442){
          JSFillHist("ScaleFactor/BB", "ZMass_OS_Unweighted_shifted_1.3%", ZCand_tmp.M(), 1., 40, 70., 110.);
          JSFillHist("ScaleFactor/BB", "ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
          JSFillHist("ScaleFactor/BB", "pt1_OS_CFweighted_shifted_1.3%", eles_tmp.at(0).Pt(), weight_tmp, 70, 20., 90.);
          JSFillHist("ScaleFactor/BB", "pt2_OS_CFweighted_shifted_1.3%", eles_tmp.at(1).Pt(), weight_tmp, 70, 20., 90.);
          JSFillHist("ScaleFactor/BB", "MET_OS_CFweighted_shifted_1.3%", METv_tmp.Pt(), weight_tmp, 100, 0., 100.);
        }
    
        // BE
        if((abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())>=1.556)||(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())<1.4442)){
          JSFillHist("ScaleFactor/BE", "ZMass_OS_Unweighted_shifted_1.3%", ZCand_tmp.M(), 1., 40, 70., 110.);
          JSFillHist("ScaleFactor/BE", "ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp_SF, 40, 70., 110.);
          JSFillHist("ScaleFactor/BE", "pt1_OS_CFweighted_shifted_1.3%", eles_tmp.at(0).Pt(), weight_tmp_SF, 70, 20., 90.);
          JSFillHist("ScaleFactor/BE", "pt2_OS_CFweighted_shifted_1.3%", eles_tmp.at(1).Pt(), weight_tmp_SF, 70, 20., 90.);
          JSFillHist("ScaleFactor/BE", "MET_OS_CFweighted_shifted_1.3%", METv_tmp.Pt(), weight_tmp_SF, 100, 0., 100.);
        }
    
        // EE
        if(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())>=1.556){
          JSFillHist("ScaleFactor/EE", "ZMass_OS_Unweighted_shifted_1.3%", ZCand_tmp.M(), 1., 40, 70., 110.);
          JSFillHist("ScaleFactor/EE", "ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
          JSFillHist("ScaleFactor/EE", "pt1_OS_CFweighted_shifted_1.3%", eles_tmp.at(0).Pt(), weight_tmp, 70, 20., 90.);
          JSFillHist("ScaleFactor/EE", "pt2_OS_CFweighted_shifted_1.3%", eles_tmp.at(1).Pt(), weight_tmp, 70, 20., 90.);
          JSFillHist("ScaleFactor/EE", "MET_OS_CFweighted_shifted_1.3%", METv_tmp.Pt(), weight_tmp, 100, 0., 100.);
        }
				
      }
    }
  }

}

double ChargeFlip::GetCFweight(std::vector<Electron> eles){
  
  double prob[2];

  for(int i=0;i<2;i++){
    if(abs(eles.at(i).scEta())<0.8){
      if(1/eles.at(i).Pt()<0.021) prob[i] = 1.63348e-04-4.97821e-03/eles.at(i).Pt();
  	  else prob[i] = 5.11382e-05+3.88422e-04/eles.at(i).Pt();
    }
    else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
  	  if(1/eles.at(i).Pt()<0.0155) prob[i] = 1.94800e-03-9.32813e-02/eles.at(i).Pt();
  	  else if(0.0155<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.023) prob[i] = 6.34964e-04-6.98093e-03/eles.at(i).Pt();
  	  else prob[i] = 5.80094e-04-4.85519e-03/eles.at(i).Pt();
  	}
  	else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
  	  if(1/eles.at(i).Pt()<0.0105) prob[i] = 1.30143e-02-7.31072e-01/eles.at(i).Pt();
  	  else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.02) prob[i] = 7.46291e-03-2.07999e-01/eles.at(i).Pt();
  	  else prob[i] = 3.74301e-03-2.23102e-02/eles.at(i).Pt();
    }
  }

  return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

}

double ChargeFlip::GetCFweight_SF(std::vector<Electron> eles){
  
  double prob[2];

  for(int i=0;i<2;i++){
    if(abs(eles.at(i).scEta())<0.8){
      if(1/eles.at(i).Pt()<0.021) prob[i] = (1.63348e-04-4.97821e-03/eles.at(i).Pt())*0.614879;
  	  else prob[i] = (5.11382e-05+3.88422e-04/eles.at(i).Pt())*0.614879;
 	  }
    else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
  	  if(1/eles.at(i).Pt()<0.0155) prob[i] = (1.94800e-03-9.32813e-02/eles.at(i).Pt())*0.614879;
  	  else if(0.0155<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.023) prob[i] = (6.34964e-04-6.98093e-03/eles.at(i).Pt())*0.614879;
  	  else prob[i] = (5.80094e-04-4.85519e-03/eles.at(i).Pt())*0.614879;
  	}
  	else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
  	  if(1/eles.at(i).Pt()<0.0105) prob[i] = (1.30143e-02-7.31072e-01/eles.at(i).Pt())*0.849782;
  	  else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.02) prob[i] = (7.46291e-03-2.07999e-01/eles.at(i).Pt())*0.849782;
  	  else prob[i] = (3.74301e-03-2.23102e-02/eles.at(i).Pt())*0.849782;
    }
  }

  return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

}

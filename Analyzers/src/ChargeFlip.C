#include "ChargeFlip.h" 

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

  /* Measure CF rate using MC */

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

    /* MC Closure test start */
    
    if(!PassMETFilter()) return;
    Particle METv = ev.GetMETVector();

    if(eles.size() == 2){

      double weight = GetCFweight(eles);

      Particle ZCand = eles.at(0)+eles.at(1);

      if(70.<=ZCand.M()&&ZCand.M()<110.){
        if(eles.at(0).Charge()*eles.at(1).Charge()<0){
          JSFillHist("ClosureTest", "ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
        }
        else{
          JSFillHist("ClosureTest", "ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
          JSFillHist("ClosureTest", "pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
          JSFillHist("ClosureTest", "pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
          JSFillHist("ClosureTest", "MET_SS", METv.Pt(), 1., 100, 0., 100.);
        }
      }
    }

    /* Now see what is changed when requiring only prompt electrons in SS events */
		
    vector<Electron> eles_prmt = ElectronPromptOnly(eles, gens); // Get prompt electrons only
 
    if(Nentry%(LogEvery)==0){
      cout << "electrons pt:" << endl;
      for(unsigned int i=0; i<eles.size(); i++){
        cout << eles.at(i).Pt() << endl;
      }
      cout << "prompt electrons pt:" << endl;
      for(unsigned int i=0; i<eles_prmt.size(); i++){
        cout << eles_prmt.at(i).Pt() << endl;
      }
    } // To see how many electrons are cut off 
		
    if(eles_prmt.size() == 2){
      Particle ZCand_prmt = eles_prmt.at(0)+eles_prmt.at(1);
      if(70.<=ZCand_prmt.M()&&ZCand_prmt.M()<110.){
        if(eles_prmt.at(0).Charge()*eles_prmt.at(1).Charge()>0){
          JSFillHist("ClosureTest", "ZMass_prmt_SS", ZCand_prmt.M(), 1., 40, 70., 110.);
          JSFillHist("ClosureTest", "pt1_prmt_SS", eles_prmt.at(0).Pt(), 1., 70, 20., 90.);
          JSFillHist("ClosureTest", "pt2_prmt_SS", eles_prmt.at(1).Pt(), 1., 70, 20., 90.);
          JSFillHist("ClosureTest", "MET_prmt_SS", METv.Pt(), 1., 100, 0., 100.);
        }
      }
    }

    //

    /* There's disagreement between OS_CFweighted and SS(where NO requirement on its promptness?). */
    /* Now, Let's shift the electrons' energy */
    
    if(eles.size() != 2) return;

    Particle ZCand = eles.at(0)+eles.at(1);
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
      METv_tmp.SetPxPyPzE(METv.Px()+ZCand.Px()-ZCand_tmp.Px(),METv.Py()+ZCand.Py()-ZCand_tmp.Py(),0,METv.E()+ZCand.E()-ZCand_tmp.E());
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

    Particle ZCand = eles.at(0)+eles.at(1);
    if(! (70.<=ZCand.M()&&ZCand.M()<110.) ) return;

    // BB
    if(abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())<1.4442){

      if(eles.at(0).Charge()*eles.at(1).Charge()>0){
        JSFillHist("ScaleFactor", "BB_ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
      }

    }

    // BE
    if((abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())>=1.556)||(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())<1.4442)){

      if(eles.at(0).Charge()*eles.at(1).Charge()>0){
        JSFillHist("ScaleFactor", "BE_ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
      }

    }

    // EE
    if(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())>=1.556){

      if(eles.at(0).Charge()*eles.at(1).Charge()>0){
        JSFillHist("ScaleFactor", "EE_ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
      }

    }

    /* Now let's shift the electrons' energy 1.3% */
    
    Particle ZCand_tmp;
    double weight_tmp, weight_tmp_SF;

    vector<Electron> eles_tmp = eles; // copy the vector
    for(int j=0;j<2;j++){
      eles_tmp.at(j).SetE(eles_tmp.at(j).E()*(1-0.013));
      eles_tmp.at(j).SetPtEtaPhiE(eles_tmp.at(j).E() * TMath::Sin(eles_tmp.at(j).Theta()), eles_tmp.at(j).Eta(), eles_tmp.at(j).Phi(), eles_tmp.at(j).E());
    }

    ZCand_tmp = eles_tmp.at(0) + eles_tmp.at(1);
    weight_tmp = GetCFweight(eles_tmp);
    weight_tmp_SF = GetCFweight_SF(eles_tmp);

    if(70.<=ZCand_tmp.M()&&ZCand_tmp.M()<110.){
      if(eles.at(0).Charge()*eles.at(1).Charge()<0){

        // BB
        if(abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())<1.4442){
          JSFillHist("ScaleFactor", "BB_ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
          JSFillHist("ScaleFactor", "BB_ZMass_OS_CFSFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp_SF, 40, 70., 110.);
        }
    
        // BE
        if((abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())>=1.556)||(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())<1.4442)){
          JSFillHist("ScaleFactor", "BE_ZMass_OS_CFSFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp_SF, 40, 70., 110.);
        }
    
        // EE
        if(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())>=1.556){
          JSFillHist("ScaleFactor", "EE_ZMass_OS_CFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
          JSFillHist("ScaleFactor", "EE_ZMass_OS_CFSFweighted_shifted_1.3%", ZCand_tmp.M(), weight_tmp_SF, 40, 70., 110.);
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

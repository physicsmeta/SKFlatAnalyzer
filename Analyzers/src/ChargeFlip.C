#include "ChargeFlip.h" 

ChargeFlip::ChargeFlip(){

}

void ChargeFlip::initializeAnalyzer(){

  if(DataYear==2016){
    EleIDs = { "HNTight2016",
               "HEID",
    }; // PassID() in Electron.C
    EleIDSFKeys = {
      "",
    }; // histmap.txt
    EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"; // script/PDandTriggers/2016/DoubleEG.txt
    lep0ptcut = 25.;
    lep1ptcut = 15.;
  }
  else if(DataYear==2017){
    EleIDs = {
      "HNTight2016",
      "HEID",
    };
    EleIDSFKeys = {
      "",
    };
    EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v";
    lep0ptcut = 25.;
    lep1ptcut = 15.;
  }
  else if(DataYear==2018){
    EleIDs = {
      "HNTight2016",
      "HEID",
    };
    EleIDSFKeys = {
      "",
    };
    EleTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    lep0ptcut = 25.;
    lep1ptcut = 15.;
  }

}

ChargeFlip::~ChargeFlip(){

}

void ChargeFlip::executeEvent(Long64_t Nentry){

  AllEles = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  AnalyzerParameter param;

  for(unsigned int i=0; i<EleIDs.size(); i++){

    TString EleID = EleIDs.at(i);
    //TString EleIDSFKey = EleIDSFKeys.at(i);
  
    param.Clear();
  
    param.CFsyst_ = AnalyzerParameter::CF_Central;
  
    param.Name = EleID;
  
    param.Electron_User_ID = EleID;
  
    executeEventFromParameter(param, Nentry);
  
    if(HasFlag("RunSyst")){
  
      for(int i=1; i<AnalyzerParameter::N_CFSyst; i++){
  
      param.CFsyst_ = AnalyzerParameter::CFSyst(i);
      param.Name = EleID+"_"+"Syst_"+param.GetCFSystType();
      executeEventFromParameter(param, Nentry);
  
      }

    }

  }

}

void ChargeFlip::executeEventFromParameter(AnalyzerParameter param, Long64_t Nentry){

  Event ev = GetEvent();

  if(param.CFsyst_ == AnalyzerParameter::CF_Central){
  
  }
  else if(param.CFsyst_ == AnalyzerParameter::MllRangeUp){
    MllLeft = 65;
    MllRight = 115;
  }
  else if(param.CFsyst_ == AnalyzerParameter::MllRangeDown){
    MllLeft = 75;
    MllRight = 105;
  }
  else if(param.CFsyst_ == AnalyzerParameter::MinPtUp){
    MinPt = 28;
  }
  else if(param.CFsyst_ == AnalyzerParameter::MinPtDown){
    MinPt = 22;
  }
  else if(param.CFsyst_ == AnalyzerParameter::NBinUp){
    NBin = 35;
  }
  else if(param.CFsyst_ == AnalyzerParameter::NBinDown){
    NBin = 45;
  }
  else{
    cout << "[ExampleRun::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  vector<Electron> eles;
  vector<Muon> muons;
  vector<Jet> jets;

  if(!IsDATA){

    /* CF ID selection */

    eles = SelectElectrons(AllEles, param.Electron_User_ID, 25., 2.5);
    //if(HasFlag("passTightChargeTightID")) eles = SelectChargeTightElectrons(AllEles, "passTightID", 25., 2.5);
    //if(HasFlag("passTightChargeTightIDdXY")) eles = SelectChargeTightElectronsDXY(AllEles, "passTightID", 25., 2.5);
    //if(HasFlag("passTightChargeTightIDdZ")) eles = SelectChargeTightElectronsDZ(AllEles, "passTightID", 25., 2.5);
    //if(HasFlag("passTightChargeTightIDdXYdZ")) eles = SelectChargeTightElectronsDXYDZ(AllEles, "passTightID", 25., 2.5);

    std::sort(eles.begin(), eles.end(), PtComparing);

    vector<Gen> gens = GetGens();
  	
    if(HasFlag("CFrate")){

      /* Measure CF rate using MC */

      for(unsigned int i=0; i<eles.size(); i++){
        Gen truth_lep = GetGenMatchedLepton(eles.at(i), gens);
        if(truth_lep.PID() == 0) return; // TODO check the meaning
    
        int truth_lep_Charge;
        if(truth_lep.PID() == 11) truth_lep_Charge = -1;
        else if(truth_lep.PID() == -11) truth_lep_Charge = 1;
      
        if(abs(eles.at(i).scEta())<0.8){
          FillHist(param.Name+"/CFrate/EtaRegion1_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "!!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            FillHist(param.Name+"/CFrate/EtaRegion1_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            FillHist(param.Name+"/CFrate/EtaRegion1_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            FillHist(param.Name+"/CFrate/EtaRegion1_IsLeading", 0, 1, 2, 0, 2);
            }
            //
          }
        }
        else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
          FillHist(param.Name+"/CFrate/EtaRegion2_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            FillHist(param.Name+"/CFrate/EtaRegion2_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            FillHist(param.Name+"/CFrate/EtaRegion2_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            FillHist(param.Name+"/CFrate/EtaRegion2_IsLeading", 0, 1, 2, 0, 2);
            }
            //
          }
        }
        else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
          FillHist(param.Name+"/CFrate/EtaRegion3_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            FillHist(param.Name+"/CFrate/EtaRegion3_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            // SH's ask //
            if(eles.size() == 2&&i == 0){
            FillHist(param.Name+"/CFrate/EtaRegion3_IsLeading", 1, 1, 2, 0, 2);
            }
            else if(eles.size() == 2&&i == 1){
            FillHist(param.Name+"/CFrate/EtaRegion3_IsLeading", 0, 1, 2, 0, 2);
            }
            //
          }
        }
      }
    }

    if(HasFlag("ClosureTest")){

      /* MC Closure test start */
      
      if(!PassMETFilter()) return;
      Particle METv = ev.GetMETVector();
  
      if(eles.size() == 2){
  
        double weight = GetCFweight(eles, param.Electron_User_ID);
  
        Particle ZCand = eles.at(0)+eles.at(1);
  
        if(70.<=ZCand.M()&&ZCand.M()<110.){
          if(eles.at(0).Charge()*eles.at(1).Charge()<0){
            FillHist(param.Name+"/ClosureTest/ZMass_OS_CFweighted", ZCand.M(), weight, 40, 70., 110.);
          }
          else{
            FillHist(param.Name+"/ClosureTest/ZMass_SS", ZCand.M(), 1., 40, 70., 110.);
            FillHist(param.Name+"/ClosureTest/pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
            FillHist(param.Name+"/ClosureTest/pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
            FillHist(param.Name+"/ClosureTest/MET_SS", METv.Pt(), 1., 100, 0., 100.);
          }
        }
      }
  
      /* Now see what is changed when requiring only prompt electrons in SS events */
  		
      vector<Electron> eles_prmt = ElectronPromptOnly(eles, gens); // Get prompt electrons only
  
      if(eles.size() == 0) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 1) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 2) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 3) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 4) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5); // To check the number of electrons
      
      if((eles.size() >= 2)&&(eles.size() == eles_prmt.size())){
        FillHist(param.Name+"/etc/IsCut", 0, 1, 2, 0, 2);
      }
      else if((eles.size() >= 2)&&(eles.size() != eles_prmt.size())){
        FillHist(param.Name+"/etc/IsCut", 1, 1, 2, 0, 2);
        cout << "[Loop " << Nentry << "]" << endl;
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
            FillHist(param.Name+"/ClosureTest/ZMass_prmt_SS", ZCand_prmt.M(), 1., 40, 70., 110.);
            FillHist(param.Name+"/ClosureTest/pt1_prmt_SS", eles_prmt.at(0).Pt(), 1., 70, 20., 90.);
            FillHist(param.Name+"/ClosureTest/pt2_prmt_SS", eles_prmt.at(1).Pt(), 1., 70, 20., 90.);
            FillHist(param.Name+"/ClosureTest/MET_prmt_SS", METv.Pt(), 1., 100, 0., 100.);
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
        weight_tmp = GetCFweight(eles_tmp, param.Electron_User_ID);
  
        if(! (70.<=ZCand_tmp.M()&&ZCand_tmp.M()<110.) ) continue;
  
        if(eles.at(0).Charge()*eles.at(1).Charge()<0){
          FillHist(param.Name+"/ClosureTest/ZMass_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), ZCand_tmp.M(), weight_tmp, 40, 70., 110.);
          FillHist(param.Name+"/ClosureTest/pt1_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), eles_tmp.at(0).Pt(), weight_tmp, 70, 20., 90.);
          FillHist(param.Name+"/ClosureTest/pt2_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), eles_tmp.at(1).Pt(), weight_tmp, 70, 20., 90.);
          FillHist(param.Name+"/ClosureTest/MET_OS_CFweighted_shifted_"+TString::Itoa(i+1,10), METv_tmp.Pt(), weight_tmp, 100, 0., 100.);
        }
      }
    }

    if(HasFlag("HalfSampleTest")){

      /* Start the Half Sample Test : Even -> observed, Odd -> predicted */

      /* Call the needed particles */

      if(!PassMETFilter()) return;
      Particle METv = ev.GetMETVector();
 
      muons = SelectMuons(AllMuons, "POGTight", 20., 2.4); //Just to calculate HT
      jets = SelectJets(AllJets, "tight", 30., 2.4); //Just to calculate HT
      vector<Jet> jets_LeptonVeto = JetsVetoLeptonInside(jets, eles, muons); //To calculate HT

      /* Calculate the needed variables */

      double MET = METv.Pt(); //MET

      double LT = 0;
      for(unsigned int i=0; i<eles.size(); i++){
        LT += eles.at(i).Pt();
      } //LT : Use electron pT sum

      if(muons.size()>0&&eles.size()>0){
        FillHist(param.Name+"/HalfSampleTest/Error_on_LT", 1, 1, 2, 0, 2);
      } //To Check if there is any event having ele and muon both (so that causing LT miscalculation - maybe removed after returning NO gen electron tho)

      double HT = 0;
      for(unsigned int i=0; i<jets_LeptonVeto.size(); i++){
        HT += jets_LeptonVeto.at(i).Pt();
      } //HT

      double ST = MET + LT + HT; //
      double METsquaredOverST = pow(MET,2)/ST; //

      if(Nentry%2==0){ //Even : measure CFrate with pT, eta as before

        for(unsigned int i=0; i<eles.size(); i++){
          Gen truth_lep = GetGenMatchedLepton(eles.at(i), gens);
          if(truth_lep.PID() == 0) return; 
      
          int truth_lep_Charge;
          if(truth_lep.PID() == 11) truth_lep_Charge = -1;
          else if(truth_lep.PID() == -11) truth_lep_Charge = 1;
        
          if(abs(eles.at(i).scEta())<0.8){
            FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion1_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            if(truth_lep_Charge*eles.at(i).Charge()<0){
              cout << "!!EtaRegion1!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
              FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion1_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            }
          }
          else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
            FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion2_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            if(truth_lep_Charge*eles.at(i).Charge()<0){
              cout << "!!EtaRegion2!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
              FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion2_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            }
          }
          else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
            FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion3_Denom", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            if(truth_lep_Charge*eles.at(i).Charge()<0){
              cout << "!!EtaRegion3!! truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
              FillHist(param.Name+"/HalfSampleTest/Even/EtaRegion3_Num", 1/eles.at(i).Pt(), 1., 40, 0., 0.04);
            }
          }
        }

      }

      if(Nentry%2==1){ //Odd : measure CFrate with MET, METsquaredOverST. Also apply the weight from the Even set

        for(unsigned int i=0; i<eles.size(); i++){
          Gen truth_lep = GetGenMatchedLepton(eles.at(i), gens);
          if(truth_lep.PID() == 0) return; 
      
          int truth_lep_Charge;
          if(truth_lep.PID() == 11) truth_lep_Charge = -1;
          else if(truth_lep.PID() == -11) truth_lep_Charge = 1;
        
          FillHist(param.Name+"/HalfSampleTest/Odd/MET_Denom", MET, 1., 100, 0., 100.);
          FillHist(param.Name+"/HalfSampleTest/Odd/METsquaredOverST_Denom", METsquaredOverST, 1., 50, 0., 50.);
          if(truth_lep_Charge*eles.at(i).Charge()<0){
            cout << "truth lepton charge : " << truth_lep_Charge << ", reco lepton charge : " << eles.at(i).Charge() << endl;
            FillHist(param.Name+"/HalfSampleTest/Odd/MET_Num", MET, 1., 100, 0., 100.);
            FillHist(param.Name+"/HalfSampleTest/Odd/METsquaredOverST_Num", METsquaredOverST, 1., 50, 0., 50.);
          }

          double weight = GetHalfSampleWeight(eles.at(i), param.Electron_User_ID);

          FillHist(param.Name+"/HalfSampleTest/Odd/MET_weight", MET, weight, 100, 0., 100.);
          FillHist(param.Name+"/HalfSampleTest/Odd/METsquaredOverST_weight", METsquaredOverST, weight, 50, 0., 50.);

        }

      }

    }

  }

  /* Scale Factor measurement */

  else{
    
    if(HasFlag("ScaleFactor")){

      if(!PassMETFilter()) return;
      if(! (ev.PassTrigger(EleTriggerName) )) return;
      
      Particle METv = ev.GetMETVector();
  
      /* CF SF ID selection */
  
      eles = SelectElectrons(AllEles, param.Electron_User_ID, MinPt, 2.5);
      //if(HasFlag("passTightChargeTightID")) eles = SelectChargeTightElectrons(AllEles, "passTightID", 25., 2.5);
      //if(HasFlag("passTightChargeTightIDdXY")) eles = SelectChargeTightElectronsDXY(AllEles, "passTightID", 25., 2.5);
      //if(HasFlag("passTightChargeTightIDdZ")) eles = SelectChargeTightElectronsDZ(AllEles, "passTightID", 25., 2.5);
      //if(HasFlag("passTightChargeTightIDdXYdZ")) eles = SelectChargeTightElectronsDXYDZ(AllEles, "passTightID", 25., 2.5);

      std::sort(eles.begin(), eles.end(), PtComparing);
  
      //
  		
      if(eles.size() == 0) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 1) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 2) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 3) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5);
      else if(eles.size() == 4) FillHist(param.Name+"/etc/eles_size", eles.size(), 1, 5, 0, 5); // To check the number of electrons
  
      if(eles.size() != 2) return;
      //if(eles.at(0).Pt()<lep0ptcut||eles.at(1).Pt()<lep1ptcut) return; //No need already pt min = 25
  
      Particle ZCand = eles.at(0)+eles.at(1);
      if(! (MllLeft<=ZCand.M()&&ZCand.M()<MllRight) ) return;
  
      // BB
      if(abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())<1.4442){
  
        if(eles.at(0).Charge()*eles.at(1).Charge()>0){
          FillHist(param.Name+"/ScaleFactor/BB_ZMass_SS", ZCand.M(), 1., NBin, MllLeft, MllRight);
          FillHist(param.Name+"/ScaleFactor/BB_pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/BB_pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/BB_MET_SS", METv.Pt(), 1., 100, 0., 100.);
        }
  
      }
  
      // BE
      if((abs(eles.at(0).scEta())<1.4442&&abs(eles.at(1).scEta())>=1.556)||(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())<1.4442)){
  
        if(eles.at(0).Charge()*eles.at(1).Charge()>0){
          FillHist(param.Name+"/ScaleFactor/BE_ZMass_SS", ZCand.M(), 1., NBin, MllLeft, MllRight);
          FillHist(param.Name+"/ScaleFactor/BE_pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/BE_pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/BE_MET_SS", METv.Pt(), 1., 100, 0., 100.);
        }
  
      }
  
      // EE
      if(abs(eles.at(0).scEta())>=1.556&&abs(eles.at(1).scEta())>=1.556){
  
        if(eles.at(0).Charge()*eles.at(1).Charge()>0){
          FillHist(param.Name+"/ScaleFactor/EE_ZMass_SS", ZCand.M(), 1., NBin, MllLeft, MllRight);
          FillHist(param.Name+"/ScaleFactor/EE_pt1_SS", eles.at(0).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/EE_pt2_SS", eles.at(1).Pt(), 1., 70, 20., 90.);
          FillHist(param.Name+"/ScaleFactor/EE_MET_SS", METv.Pt(), 1., 100, 0., 100.);
        }
  
      }
  
      /* Now let's shift the electrons' energy 1.4% */
      
      vector<Electron> eles_tmp = eles; // copy the vector
  		
      for(int j=0;j<2;j++){
        eles_tmp.at(j).SetE(eles_tmp.at(j).E()*(1-0.014));
        eles_tmp.at(j).SetPtEtaPhiE(eles_tmp.at(j).E() * TMath::Sin(eles_tmp.at(j).Theta()), eles_tmp.at(j).Eta(), eles_tmp.at(j).Phi(), eles_tmp.at(j).E());
      }
  
  		Particle ZCand_tmp = eles_tmp.at(0) + eles_tmp.at(1);
      double weight_tmp = GetCFweight(eles_tmp, param.Electron_User_ID);
      double weight_tmp_SF = GetCFweight_SF(eles_tmp, param.Electron_User_ID);
  
      Particle METv_tmp;
      METv_tmp.SetPxPyPzE(METv.Px()+ZCand.Px()-ZCand_tmp.Px(),METv.Py()+ZCand.Py()-ZCand_tmp.Py(),0,METv.E()+ZCand.E()-ZCand_tmp.E());
  
      if(MllLeft<=ZCand_tmp.M()&&ZCand_tmp.M()<MllRight){
        if(eles.at(0).Charge()*eles.at(1).Charge()<0){
  
          // BB
          if(abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())<1.4442){
            FillHist(param.Name+"/ScaleFactor/BB_ZMass_OS_CFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/BB_ZMass_OS_CFSFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp_SF, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/BB_pt1_OS_CFSFweighted_shifted_1.4%", eles_tmp.at(0).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/BB_pt2_OS_CFSFweighted_shifted_1.4%", eles_tmp.at(1).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/BB_MET_OS_CFSFweighted_shifted_1.4%", METv_tmp.Pt(), weight_tmp_SF, 100, 0., 100.);
          }
      
          // BE
          if((abs(eles_tmp.at(0).scEta())<1.4442&&abs(eles_tmp.at(1).scEta())>=1.556)||(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())<1.4442)){
            FillHist(param.Name+"/ScaleFactor/BE_ZMass_OS_CFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/BE_ZMass_OS_CFSFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp_SF, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/BE_pt1_OS_CFSFweighted_shifted_1.4%", eles.at(0).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/BE_pt2_OS_CFSFweighted_shifted_1.4%", eles.at(1).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/BE_MET_OS_CFSFweighted_shifted_1.4%", METv.Pt(), weight_tmp_SF, 100, 0., 100.);
          }
      
          // EE
          if(abs(eles_tmp.at(0).scEta())>=1.556&&abs(eles_tmp.at(1).scEta())>=1.556){
            FillHist(param.Name+"/ScaleFactor/EE_ZMass_OS_CFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/EE_ZMass_OS_CFSFweighted_shifted_1.4%", ZCand_tmp.M(), weight_tmp_SF, NBin, MllLeft, MllRight);
            FillHist(param.Name+"/ScaleFactor/EE_pt1_OS_CFSFweighted_shifted_1.4%", eles.at(0).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/EE_pt2_OS_CFSFweighted_shifted_1.4%", eles.at(1).Pt(), weight_tmp_SF, 70, 20., 90.);
            FillHist(param.Name+"/ScaleFactor/EE_MET_OS_CFSFweighted_shifted_1.4%", METv.Pt(), weight_tmp_SF, 100, 0., 100.);
          }
  				
        }
      }
    }

  }

}

double ChargeFlip::GetCFweight(std::vector<Electron> eles, TString id){
  
  double prob[2];

/* == DY only CF ==

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

*/ 

  if(id == "HEID"){

    for(int i=0;i<2;i++){ //DY+TTLL CF

      if(abs(eles.at(i).scEta())<0.8){
        if(1/eles.at(i).Pt()<0.0075) prob[i] = 4.59433e-04-4.41286e-02/eles.at(i).Pt();
        else if(0.0075<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = 2.01545e-04-8.19670e-03/eles.at(i).Pt();
        else prob[i] = 8.74309e-05-5.72703e-04/eles.at(i).Pt();
      }
      else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
        if(1/eles.at(i).Pt()<0.0055) prob[i] = 3.81412e-03-4.74900e-01/eles.at(i).Pt();
        else if(0.0055<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = 1.60397e-03-6.97824e-02/eles.at(i).Pt();
        else prob[i] = 6.59614e-04-7.13587e-03/eles.at(i).Pt();
      }
      else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
        if(1/eles.at(i).Pt()<0.0105) prob[i] = 1.23549e-02-6.80049e-01/eles.at(i).Pt();
        else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0205) prob[i] = 7.27935e-03-1.88277e-01/eles.at(i).Pt();
        else prob[i] = 4.06341e-03-3.22562e-02/eles.at(i).Pt();
      }

    }
  
    return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

  }

  else if(id == "HNTight2016"){

    for(int i=0;i<2;i++){ //DY+TTLL CF

      if(abs(eles.at(i).scEta())<0.8){
        if(1/eles.at(i).Pt()<0.0075) prob[i] = 4.59433e-04-4.41286e-02/eles.at(i).Pt();
        else if(0.0075<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = 2.01545e-04-8.19670e-03/eles.at(i).Pt();
        else prob[i] = 8.74309e-05-5.72703e-04/eles.at(i).Pt();
      }
      else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
        if(1/eles.at(i).Pt()<0.0055) prob[i] = 3.81412e-03-4.74900e-01/eles.at(i).Pt();
        else if(0.0055<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = 1.60397e-03-6.97824e-02/eles.at(i).Pt();
        else prob[i] = 6.59614e-04-7.13587e-03/eles.at(i).Pt();
      }
      else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
        if(1/eles.at(i).Pt()<0.0105) prob[i] = 1.23549e-02-6.80049e-01/eles.at(i).Pt();
        else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0205) prob[i] = 7.27935e-03-1.88277e-01/eles.at(i).Pt();
        else prob[i] = 4.06341e-03-3.22562e-02/eles.at(i).Pt();
      }

    }
  
    return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

  }

}

double ChargeFlip::GetCFweight_SF(std::vector<Electron> eles, TString id){
  
  double prob[2];

/* == DY only CF ==

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

*/

  if(id == "HEID"){

    for(int i=0;i<2;i++){ //DY+TTLL CF

      if(abs(eles.at(i).scEta())<0.8){
        if(1/eles.at(i).Pt()<0.0075) prob[i] = (4.59433e-04-4.41286e-02/eles.at(i).Pt())*0.585841;
        else if(0.0075<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = (2.01545e-04-8.19670e-03/eles.at(i).Pt())*0.585841;
        else prob[i] = (8.74309e-05-5.72703e-04/eles.at(i).Pt())*0.585841;
      }
      else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
        if(1/eles.at(i).Pt()<0.0055) prob[i] = (3.81412e-03-4.74900e-01/eles.at(i).Pt())*0.585841;
        else if(0.0055<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = (1.60397e-03-6.97824e-02/eles.at(i).Pt())*0.585841;
        else prob[i] = (6.59614e-04-7.13587e-03/eles.at(i).Pt())*0.585841;
      }
      else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
        if(1/eles.at(i).Pt()<0.0105) prob[i] = (1.23549e-02-6.80049e-01/eles.at(i).Pt())*0.831539;
        else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0205) prob[i] = (7.27935e-03-1.88277e-01/eles.at(i).Pt())*0.831539;
        else prob[i] = (4.06341e-03-3.22562e-02/eles.at(i).Pt())*0.831539;
      }

    }

    return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

  }

  else if(id == "HNTight2016"){

    for(int i=0;i<2;i++){ //DY+TTLL CF

      if(abs(eles.at(i).scEta())<0.8){
        if(1/eles.at(i).Pt()<0.0075) prob[i] = (4.59433e-04-4.41286e-02/eles.at(i).Pt())*0.585841;
        else if(0.0075<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = (2.01545e-04-8.19670e-03/eles.at(i).Pt())*0.585841;
        else prob[i] = (8.74309e-05-5.72703e-04/eles.at(i).Pt())*0.585841;
      }
      else if(0.8<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<1.4442){
        if(1/eles.at(i).Pt()<0.0055) prob[i] = (3.81412e-03-4.74900e-01/eles.at(i).Pt())*0.585841;
        else if(0.0055<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0155) prob[i] = (1.60397e-03-6.97824e-02/eles.at(i).Pt())*0.585841;
        else prob[i] = (6.59614e-04-7.13587e-03/eles.at(i).Pt())*0.585841;
      }
      else if(1.556<=abs(eles.at(i).scEta())&&abs(eles.at(i).scEta())<2.5){
        if(1/eles.at(i).Pt()<0.0105) prob[i] = (1.23549e-02-6.80049e-01/eles.at(i).Pt())*0.831539;
        else if(0.0105<=1/eles.at(i).Pt()&&1/eles.at(i).Pt()<0.0205) prob[i] = (7.27935e-03-1.88277e-01/eles.at(i).Pt())*0.831539;
        else prob[i] = (4.06341e-03-3.22562e-02/eles.at(i).Pt())*0.831539;
      }

    }

    return prob[0]/(1.-prob[0])+prob[1]/(1.-prob[1]);

  }

}

double ChargeFlip::GetHalfSampleWeight(const Electron& electron, TString id){;

  if(id == "HEID"){

    if(abs(electron.scEta())<0.8){
      if(1/electron.Pt()<0.0075) return (4.59433e-04-4.41286e-02/electron.Pt())*0.585841;
      else if(0.0075<=1/electron.Pt()&&1/electron.Pt()<0.0155) return (2.01545e-04-8.19670e-03/electron.Pt())*0.585841;
      else return (8.74309e-05-5.72703e-04/electron.Pt())*0.585841;
    }
    else if(0.8<=abs(electron.scEta())&&abs(electron.scEta())<1.4442){
      if(1/electron.Pt()<0.0055) return (3.81412e-03-4.74900e-01/electron.Pt())*0.585841;
      else if(0.0055<=1/electron.Pt()&&1/electron.Pt()<0.0155) return (1.60397e-03-6.97824e-02/electron.Pt())*0.585841;
      else return (6.59614e-04-7.13587e-03/electron.Pt())*0.585841;
    }
    else if(1.556<=abs(electron.scEta())&&abs(electron.scEta())<2.5){
      if(1/electron.Pt()<0.0105) return (1.23549e-02-6.80049e-01/electron.Pt())*0.831539;
      else if(0.0105<=1/electron.Pt()&&1/electron.Pt()<0.0205) return (7.27935e-03-1.88277e-01/electron.Pt())*0.831539;
      else return (4.06341e-03-3.22562e-02/electron.Pt())*0.831539;
    }

  }

  else if(id == "HNTight2016"){

    if(abs(electron.scEta())<0.8){
      if(1/electron.Pt()<0.0075) return (4.59433e-04-4.41286e-02/electron.Pt())*0.585841;
      else if(0.0075<=1/electron.Pt()&&1/electron.Pt()<0.0155) return (2.01545e-04-8.19670e-03/electron.Pt())*0.585841;
      else return (8.74309e-05-5.72703e-04/electron.Pt())*0.585841;
    }
    else if(0.8<=abs(electron.scEta())&&abs(electron.scEta())<1.4442){
      if(1/electron.Pt()<0.0055) return (3.81412e-03-4.74900e-01/electron.Pt())*0.585841;
      else if(0.0055<=1/electron.Pt()&&1/electron.Pt()<0.0155) return (1.60397e-03-6.97824e-02/electron.Pt())*0.585841;
      else return (6.59614e-04-7.13587e-03/electron.Pt())*0.585841;
    }
    else if(1.556<=abs(electron.scEta())&&abs(electron.scEta())<2.5){
      if(1/electron.Pt()<0.0105) return (1.23549e-02-6.80049e-01/electron.Pt())*0.831539;
      else if(0.0105<=1/electron.Pt()&&1/electron.Pt()<0.0205) return (7.27935e-03-1.88277e-01/electron.Pt())*0.831539;
      else return (4.06341e-03-3.22562e-02/electron.Pt())*0.831539;
    }

  }

}

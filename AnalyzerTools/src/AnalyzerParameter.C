#include "AnalyzerParameter.h"

void AnalyzerParameter::Clear(){

  Name = "";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "";
  Electron_Loose_ID = "";
  Electron_Veto_ID = "";
  Electron_User_ID = "";
  Electron_ID_SF_Key = "";
  Electron_FR_ID = "";
  Electron_FR_Key = "";
  Electron_CF_ID = "";
  Electron_CF_Key = "";
  Electron_Tight_RelIso = 999.;
  Electron_Loose_RelIso = 999.;
  Electron_Veto_RelIso = 999.;
  Electron_UseMini = false;
  Electron_UsePtCone = false;
  Electron_MinPt = 10.;

  Muon_Tight_ID = "";

  Muon_Loose_ID = "";
  Muon_Veto_ID = "";
  Muon_RECO_SF_Key = "";
  Muon_ID_SF_Key = "";
  Muon_ISO_SF_Key = "";
  Muon_Trigger_SF_Key = "";
  Muon_FR_ID = "";
  Muon_FR_Key = "";
  Muon_CF_ID = "";
  Muon_CF_Key = "";
  Muon_Tight_RelIso = 999.;
  Muon_Loose_RelIso = 999.;
  Muon_Veto_RelIso = 999.;
  Muon_UseMini = false;
  Muon_UsePtCone = false;
  Muon_UseTuneP = false;
  Muon_MinPt = 10.;

  Jet_ID = "";
  FatJet_ID = "";

  syst_ = Central;
  CFsyst_ = CF_Central;

}

AnalyzerParameter::AnalyzerParameter(){

  Name = "Default";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "passTightID";
  Electron_Loose_ID = "passLooseID";
  Electron_Veto_ID = "passVetoID";
  Electron_User_ID = "HEID"; //JH
  Electron_ID_SF_Key = "passTightID";

  Muon_Tight_ID = "POGTightWithTightIso";
  Muon_Loose_ID = "POGLoose";
  Muon_Veto_ID = "POGLoose";
  Muon_RECO_SF_Key = "";
  Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  Muon_Trigger_SF_Key = "POGTight";

  Jet_ID = "HN";
  FatJet_ID = "HN";

  syst_ = Central; //JH: enum Syst{ Central, JetResUp, JetResDown, ..., NSyst }; Syst syst_; AnalyzerParameter.h
  CFsyst_ = CF_Central; //JH for CF SF syst

}

TString AnalyzerParameter::GetSystType(){

  if(syst_==Syst::Central){
    return "Central";
  }
  else if(syst_==Syst::JetResUp){
    return "JetResUp";
  }
  else if(syst_==Syst::JetResDown){
    return "JetResDown";
  }
  else if(syst_==Syst::JetEnUp){
    return "JetEnUp";
  }
  else if(syst_==Syst::JetEnDown){
    return "JetEnDown";
  }
  else if(syst_==Syst::MuonEnUp){
    return "MuonEnUp";
  }
  else if(syst_==Syst::MuonEnDown){
    return "MuonEnDown";
  }
  else if(syst_==Syst::ElectronResUp){
    return "ElectronResUp";
  }
  else if(syst_==Syst::ElectronResDown){
    return "ElectronResDown";
  }
  else if(syst_==Syst::ElectronEnUp){
    return "ElectronEnUp";
  }
  else if(syst_==Syst::ElectronEnDown){
    return "ElectronEnDown";
  }
  else{
    cout << "[AnalyzerParameter::GetSystType] Wrong Syst" << endl;
    exit(EXIT_FAILURE);
    return "ERROR";
  }

}

TString AnalyzerParameter::GetCFSystType(){

  if(CFsyst_==CFSyst::CF_Central){
    return "Central";
  }
  else if(CFsyst_==CFSyst::MllRangeUp){
    return "MllRangeUp";
  }
  else if(CFsyst_==CFSyst::MllRangeDown){
    return "MllRangeDown";
  }
  else if(CFsyst_==CFSyst::MinPtUp){
    return "MinPtUp";
  }
  else if(CFsyst_==CFSyst::MinPtDown){
    return "MinPtDown";
  }
  else if(CFsyst_==CFSyst::NBinUp){
    return "NBinUp";
  }
  else if(CFsyst_==CFSyst::NBinDown){
    return "NBinDown";
  }
  else{
    cout << "[AnalyzerParameter::GetCFSystType] Wrong Syst" << endl;
    exit(EXIT_FAILURE);
    return "ERROR";
  }

}

AnalyzerParameter::~AnalyzerParameter(){
}

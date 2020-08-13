vector<TString> channels = {"dimu", "diel", "emu"};
vector<TString> regions = {"lowSR1", "lowSR2", "highSR1", "highSR2"};
vector<TString> IDs = {"HN16", "HNV2"};
vector<TString> periods = {"B_ver2", "C", "D", "E", "F", "G", "H"};
vector<TString> weights = {"Number_Events", "Number_Events_unweighted"};

vector<TString> indices = {"SingleTrigger_pt_50_50","dimuTrigger_pt20_10","dimuTrigger_pt50_50"};
vector<TString> samples = {"HNToMuMu_Sch_M70","HNToMuMu_Sch_M90","HNToMuMu_Sch_M100","HNToMuMu_Sch_M200"};
vector<TString> Signal_IDs = {"LRSM", "HNV2"};

ofstream DATA_signal_eff("Cutflow/2016/signal_eff.txt");
ofstream signal_eff("Cutflow/Signal/2016/signal_eff.txt");

void DrawCutflow(TString channel, TString region, TString ID, TString period, int doWeight = 0, TString SaveAs = "n"){ //ex : dimu, highSR1, HNV2, B_ver2

  TString filename = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Signal/2016/DATA/Signal_DoubleMuon_"+period+".root";
  TFile* f1 = new TFile(filename);
  
  for(int i=0; i<IDs.size(); i++){
    gSystem->Exec("mkdir -p Cutflow/2016/"+IDs.at(i));
  }
  
  TString title;
  if(doWeight==0) title = channel+"_"+region+"_"+ID+"_"+period+"_Cutflow";
  else if(doWeight==1) title = channel+"_"+region+"_"+ID+"_"+period+"_Cutflow (weighted)";

  TH1D *h1 = (TH1D*)f1->Get(channel+"/"+region+"/"+weights.at(doWeight)+"_"+ID);
  TH1D *h2 = new TH1D("h2",title,10,0,10);
  for(int i=0; i<11; i++) h2->SetBinContent(i+1,h1->GetBinContent(i)/h1->GetBinContent(0));

  TCanvas *c1 = new TCanvas("c1",title,200,200,900,800);
  c1->SetLogy();

  h2->SetTitle(title);
  h2->SetStats(0);
  h2->GetXaxis()->ChangeLabel(1,15,0.02,12,-1,-1,"Nocut");
  h2->GetXaxis()->ChangeLabel(2,15,0.02,12,-1,-1,"MET filter");
  h2->GetXaxis()->ChangeLabel(3,15,0.02,12,-1,-1,"Dilepton Trigger");
  h2->GetXaxis()->ChangeLabel(4,15,0.02,12,-1,-1,"2 tight leptons");
  h2->GetXaxis()->ChangeLabel(5,15,0.02,12,-1,-1,"same sign");
  h2->GetXaxis()->ChangeLabel(6,15,0.02,12,-1,-1,"3rd lepton veto");
  h2->GetXaxis()->ChangeLabel(7,15,0.02,12,-1,-1,"m(ll) > 10GeV");
  h2->GetXaxis()->ChangeLabel(8,15,0.02,12,-1,-1,"jet preselection");
  h2->GetXaxis()->ChangeLabel(9,15,0.02,12,-1,-1,"Signal selection");
  h2->GetXaxis()->ChangeLabel(10,-1,0,-1,-1,-1,"");
  h2->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,"");
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetTitleOffset(1);
  h2->GetYaxis()->SetTitle("Efficiency");
  h2->GetYaxis()->SetRangeUser(5e-07,5);
  h2->Draw();

  long int Nevent = h1->GetBinContent(0);
  TString Nevent_t = Form("%ld",Nevent);
  TPaveText *pt = new TPaveText(0.75,0.65,0.9,0.9,"NDC");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->AddText("total Nevent : "+Nevent_t);
  for(int i=0; i<9; i++){
    TString i_t = Form("%d",i+1);
    TString item_t = Form("%f",100*h2->GetBinContent(i+1));
    pt->AddText(i_t+"th eff. : "+item_t+"%");
  }
  pt->Draw();
  
  DATA_signal_eff << title << " : " << 100*h2->GetBinContent(9) << " (%)" << "\n"; //to check signal eff

  if(SaveAs == "y"){
    c1->SaveAs("Cutflow/2016/"+ID+"/"+title+".png");
  }

}

void DrawSignalCutflow(TString index, TString sample, TString ID, TString channel = "dimu", TString region = "lowSR1", TString SaveAs = "n"){

  TString filename = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/HNtypeI_SR/2016/"+index+"/HNtypeI_SR_"+sample+".root";
  TFile* f1 = new TFile(filename);
  
  //for(int i=0; i<IDs.size(); i++){
  //  gSystem->Exec("mkdir -p Cutflow/2016/"+IDs.at(i));
  //}
  for(int i=0; i<Signal_IDs.size(); i++){
    gSystem->Exec("mkdir -p Cutflow/Signal/2016/"+index+"/"+Signal_IDs.at(i));
	}

  TString title;
  title = index+"_"+sample+"_"+channel+"_Preselection_"+ID+"_Cutflow";

  TH1D *h1 = (TH1D*)f1->Get(channel+"/"+region+"/Number_Events_"+ID);
  TH1D *h2 = new TH1D("h2",title,10,0,10);
  for(int i=0; i<8; i++) h2->SetBinContent(i+1,h1->GetBinContent(i)/h1->GetBinContent(0));

  TCanvas *c1 = new TCanvas("c1","c1",200,200,900,800);
  c1->SetLogy();

  h2->SetTitle(sample+"_"+channel+"_Preselection_Cutflow");
  h2->SetStats(0);
  h2->GetXaxis()->ChangeLabel(1,15,0.02,12,-1,-1,"Nocut");
  h2->GetXaxis()->ChangeLabel(2,15,0.02,12,-1,-1,"MET filter");
  h2->GetXaxis()->ChangeLabel(3,15,0.02,12,-1,-1,"Dilepton Trigger");
  h2->GetXaxis()->ChangeLabel(4,15,0.02,12,-1,-1,"2 tight leptons");
  h2->GetXaxis()->ChangeLabel(5,15,0.02,12,-1,-1,"same sign");
  h2->GetXaxis()->ChangeLabel(6,15,0.02,12,-1,-1,"3rd lepton veto");
  h2->GetXaxis()->ChangeLabel(7,15,0.02,12,-1,-1,"m(ll) > 10GeV");
  h2->GetXaxis()->ChangeLabel(8,15,0.02,12,-1,-1,"jet selectioin(final)");
  //h2->GetXaxis()->ChangeLabel(9,15,0.02,12,-1,-1,"Signal selection");
  h2->GetXaxis()->ChangeLabel(9,-1,0,-1,-1,-1,"");
  h2->GetXaxis()->ChangeLabel(10,-1,0,-1,-1,-1,"");
  h2->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,"");
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetTitleOffset(1);
  h2->GetYaxis()->SetTitle("Efficiency");
  h2->GetYaxis()->SetRangeUser(1.e-04,5);
  h2->Draw();

  //long int Nevent = h1->GetBinContent(0);
  //TString Nevent_t = Form("%ld",Nevent);
  TPaveText *pt = new TPaveText(0.75,0.65,0.9,0.9,"NDC");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  //pt->AddText("total Nevent : "+Nevent_t);
  for(int i=0; i<8; i++){
    TString i_t = Form("%d",i+1);
    TString item_t = Form("%f",100*h2->GetBinContent(i+1));
    pt->AddText(i_t+"th eff. : "+item_t+"%");
  }
  pt->Draw();
  
  signal_eff << title << " : " << 100*h2->GetBinContent(8) << " (%)" << "\n"; //to check signal eff

  if(SaveAs == "y"){
    c1->SaveAs("Cutflow/Signal/2016/"+ID+"/"+title+".png");
  }

}

void SaveAll_Signal(){

  for(int i=0; i<indices.size(); i++){
    for(int j=0; j<samples.size(); j++){
      for(int k=0; k<Signal_IDs.size(); k++){
        DrawSignalCutflow(indices.at(i), samples.at(j), Signal_IDs.at(k), "dimu", "lowSR1", "y");
      }
    }
  }

}

void SaveAll(){

  for(int j=0; j<regions.size(); j++){
    for(int k=0; k<IDs.size(); k++){
      for(int l=0; l<periods.size()-1; l++){
        DrawCutflow("dimu", regions.at(j), IDs.at(k), periods.at(l), 0, "y");
      }
    }
  }

}

void NoSaveAll(){ //to get calculation results only; Use this with root -l -b options

  for(int j=0; j<regions.size(); j++){
    for(int k=0; k<IDs.size(); k++){
      for(int l=0; l<periods.size()-1; l++){
        DrawCutflow("dimu", regions.at(j), IDs.at(k), periods.at(l), 0, "n");
      }
    }
  }

}


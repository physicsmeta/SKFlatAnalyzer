void DrawCutflow(TString channel, TString region, TString ID, int weight = 0, TString SaveAs = "n"){ //ex : dimu, highSR1, 1, HNV2

  TString filename = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Signal/2016/DATA/Signal_DoubleMuon_B_ver2.root";
  TFile* f1 = new TFile(filename);
  
  //gSystem->Exec("mkdir Cutflow");
  
  vector<TString> channels = {"dimu", "diel", "emu"};
  vector<TString> regions = {"lowSR1", "lowSR2", "highSR1", "highSR2"};
  vector<TString> weights = {"Number_Events", "Number_Events_unweighted"};
  vector<TString> IDs = {"HN16", "HNV2"};
  
  TString title;
  if(weight==0) title = channel+"_"+region+"_"+ID+"_Cutflow";
  else if(weight==1) title = channel+"_"+region+"_"+ID+"_Cutflow (weighted)";

  TH1D *h1 = (TH1D*)f1->Get(channel+"/"+region+"/"+weights.at(weight)+"_"+ID);
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
  h2->Draw();

  if(SaveAs == "y"){
    c1->SaveAs("Cutflow"+title+".pdf");
    //c1->SaveAs("Cutflow"+title+".png");
  }

}

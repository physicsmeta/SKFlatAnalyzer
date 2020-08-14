void plot_kinematics(TString sample, TString var1, TString var2, TString title, double xran1, double xran2, double yran, int rebin, TString SaveAs = "n"){

  TString filename = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/SampleValidation/2016/SampleValidation_"+sample+".root";
  TFile* file = new TFile(filename);
  TH1D* hist1 = (TH1D*)file->Get(var1);
  TH1D* hist2 = (TH1D*)file->Get(var2);
  
  //gSystem->Exec("rootls "+filename);
  gSystem->Exec("mkdir -p Sample_kinematics/"+sample);

  TCanvas* c1 = new TCanvas("c1",title,200,350,700,650);
  c1->cd();
  
  //hist1->SetTitle(title+" comparison");
  hist1->SetStats(0);
  hist1->Rebin(rebin);
  hist1->GetXaxis()->SetRangeUser(xran1,xran2);
  hist1->GetYaxis()->SetRangeUser(0,yran);
  if(var1.Contains("Pt")) hist1->GetXaxis()->SetTitle(title+" [GeV]");
  else hist1->GetXaxis()->SetTitle(title);
  //hist1->SetLineColor(kRed);
  hist1->SetLineWidth(3);
  hist1->Draw("hist");
  hist2->SetStats(0);
  hist2->Rebin(rebin);
  hist2->SetLineColor(kRed);
  hist2->SetLineWidth(3);
  hist2->SetLineStyle(7);
  hist2->Draw("same hist");
  
  TLegend* legend;
  if(var1.Contains("Pt")){
    double ptdiff_1 = 100*(hist2->GetMean()-hist1->GetMean())/hist1->GetMean();
    TString ptdiff_1_t = Form("%f",ptdiff_1);
    legend = new TLegend(0.62,0.75,0.9,0.9);
    legend->AddEntry((TObject*)0,sample,"");
    legend->AddEntry(hist1,var1,"l");
    legend->AddEntry(hist2,var2,"l");
    legend->AddEntry((TObject*)0,"#bf{mean p_{T}} diff : "+ptdiff_1_t+"%","");
    legend->Draw();
  }
  else{
    legend = new TLegend(0.62,0.75,0.9,0.9);
    legend->AddEntry((TObject*)0,sample,"");
    legend->AddEntry(hist1,var1,"l");
    legend->AddEntry(hist2,var2,"l");
    legend->Draw();
  }

  if(SaveAs=="y") c1->SaveAs("Sample_kinematics/"+sample+"/"+var2+".png");

}



void SaveAll(){
  
  vector<TString> samples = {"HNToMuMu_Sch_M70","HNToMuMu_Sch_M90","HNToMuMu_Sch_M100","HNToMuMu_Sch_M200"};
  for(unsigned int i=0; i<samples.size(); i++){
    plot_kinematics(samples.at(i), "fid_Gen_Pt_HN_mu", "fid_Pt_HN_mu", "HN_mu pt", 0, 100, 10000, 1, "y");
    plot_kinematics(samples.at(i), "fid_Gen_Pt_hard_mu", "fid_Pt_hard_mu", "hard_mu pt", 0, 100, 10000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Pt_j1", "Pt_j1", "j1 pt", 0, 100, 10000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Pt_j2", "Pt_j2", "j2 pt", 0, 100, 10000, 1, "y");
    plot_kinematics(samples.at(i), "fid_Gen_Eta_HN_mu", "fid_Eta_HN_mu", "HN_mu eta", -2.4, 2.4, 2000, 1, "y");
    plot_kinematics(samples.at(i), "fid_Gen_Eta_hard_mu", "fid_Eta_hard_mu", "hard_mu eta", -2.4, 2.4, 2000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Eta_j1", "Eta_j1", "j1 eta", -2.7, 2.7, 2000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Eta_j2", "Eta_j2", "j2 eta", -2.7, 2.7, 2000, 1, "y");
    plot_kinematics(samples.at(i), "fid_Gen_Phi_HN_mu", "fid_Phi_HN_mu", "HN_mu eta", -3.1, 3.1, 2000, 1, "y");
    plot_kinematics(samples.at(i), "fid_Gen_Phi_hard_mu", "fid_Phi_hard_mu", "hard_mu eta", -3.1, 3.1, 2000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Phi_j1", "Phi_j1", "j1 eta", -3.1, 3.1, 2000, 1, "y");
    plot_kinematics(samples.at(i), "Gen_Phi_j2", "Phi_j2", "j2 eta", -3.1, 3.1, 2000, 1, "y");
  }
}

void plot_presel(){

ifstream in("list_presel.txt");
string line;

while(getline(in,line)){
  istringstream is(line);
  TString this_line = line;
  if(this_line.Contains("#")) continue;

  TString sample, var, title, SaveAs;
  int xran1, xran2, rebin;
  double yran;
  is >> sample;
  is >> var;
  is >> title;
  is >> xran1;
  is >> xran2;
  is >> yran;
  is >> rebin;
  is >> SaveAs;

  TString filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2016/Presel_"+sample+".root";
  TString filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2017/Presel_"+sample+".root";
  TString filename3 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+".root";
  TFile* file1 = new TFile(filename1);
  TFile* file2 = new TFile(filename2);
  TFile* file3 = new TFile(filename3);
  TH1D* var1;
  TH1D* var2;
  TH1D* var3;
  if(sample.Contains("EE") || sample.Contains("DoubleEG") || sample.Contains("EGamma")){
  var1 = (TH1D*)file1->Get("diel/Pre/"+var+"_Run2");
  var2 = (TH1D*)file2->Get("diel/Pre/"+var+"_Run2");
  var3 = (TH1D*)file3->Get("diel/Pre/"+var+"_Run2");
  }
  else if(sample.Contains("MuMu") || sample.Contains("DoubleMuon")){
  var1 = (TH1D*)file1->Get("dimu/Pre/"+var+"_Run2");
  var2 = (TH1D*)file2->Get("dimu/Pre/"+var+"_Run2");
  var3 = (TH1D*)file3->Get("dimu/Pre/"+var+"_Run2");
  }
  else if(sample.Contains("MuonEG")){
  var1 = (TH1D*)file1->Get("emu/Pre/"+var+"_Run2");
  var2 = (TH1D*)file2->Get("emu/Pre/"+var+"_Run2");
  var3 = (TH1D*)file3->Get("emu/Pre/"+var+"_Run2");
  }
  
  //gSystem->Exec("rootls "+filename1);
  gSystem->Exec("mkdir -p Presel");

  TCanvas* c1 = new TCanvas("c1",title,200,350,700,650);
  c1->cd();
  
  var1->SetTitle(title+" #scale[0.5]{"+sample+"}");
  var1->SetStats(0);
  var1->Rebin(rebin);
  var1->GetXaxis()->SetRangeUser(xran1,xran2);
  var1->GetYaxis()->SetRangeUser(0,yran);
  if(var.Contains("Pt") || var.Contains("Mass")) var1->GetXaxis()->SetTitle(title+" [GeV]");
  else var1->GetXaxis()->SetTitle(title);
  var1->SetMarkerStyle(8);
  var1->SetMarkerSize(0.8);
  var1->SetMarkerColor(kRed);
  var1->SetLineColor(kRed);
  var1->SetLineWidth(2);
  //var1->SetLineStyle(7);
  var1->Draw("hist ep");
  var2->Scale(var1->GetEntries()/var2->GetEntries());
  var2->Rebin(rebin);
  var2->SetMarkerStyle(8);
  var2->SetMarkerSize(0.8);
  var2->SetMarkerColor(kGreen+1);
  var2->SetLineColor(kGreen+1);
  var2->SetLineWidth(2);
  //var2->SetLineStyle(7);
  var2->Draw("same hist ep");
  var3->Scale(var1->GetEntries()/var3->GetEntries());
  var3->Rebin(rebin);
  var3->SetMarkerStyle(8);
  var3->SetMarkerSize(0.8);
  var3->SetMarkerColor(kBlue);
  var3->SetLineColor(kBlue);
  var3->SetLineWidth(2);
  //var3->SetLineStyle(7);
  var3->Draw("same hist ep");
  
  TLegend* legend;
  if(var.Contains("pt")){
    legend = new TLegend(0.62,0.65,0.9,0.9);
    legend->AddEntry(var1,"#scale[0.9]{2016 "+sample+"}","l");
    legend->AddEntry(var2,"#scale[0.9]{2017 "+sample+"}","l");
    legend->AddEntry(var3,"#scale[0.9]{2018 "+sample+"}","l");
    legend->Draw();
  }
  else{
    legend = new TLegend(0.62,0.75,0.9,0.9);
    legend->AddEntry(var1,"#scale[0.9]{2016 "+sample+"}","l");
    legend->AddEntry(var2,"#scale[0.9]{2017 "+sample+"}","l");
    legend->AddEntry(var3,"#scale[0.9]{2018 "+sample+"}","l");
    legend->Draw();
  }

  if(SaveAs=="y") c1->SaveAs("Presel/"+sample+"_"+var+".png");
}

}

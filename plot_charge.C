void plot_charge(){

ifstream in("list_charge.txt");
string line;

while(getline(in,line)){
  istringstream is(line);
  TString this_line = line;
  if(this_line(0,1)=="#"||this_line=="") continue;

  TString sample, channel, var, title, this_xran1, this_xran2, this_yran, SaveAs;
  double xran1, xran2, yran;
  int rebin;
  is >> sample;
  is >> channel;
  is >> var;
  is >> title;
  is >> this_xran1;
  is >> this_xran2;
  is >> this_yran;
  is >> rebin;
  is >> SaveAs;

  TString sample_c = sample;
  if(sample.Contains("SS")) sample_c.ReplaceAll("SS","OS");
  else if(sample.Contains("OS")) sample_c.ReplaceAll("OS","SS");
  TString filename;
  TString filename_c;
  TFile* file;
  TFile* file_c;
  TH1D* hist;
  TH1D* hist_c;

  filename = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+".root";
  filename_c = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample_c+".root";
  file = new TFile(filename);
  file_c = new TFile(filename_c);
  hist = (TH1D*)file->Get(channel+"/Pre/"+var+"_Run2");
  hist_c = (TH1D*)file_c->Get(channel+"/Pre/"+var+"_Run2");
  
  //gSystem->Exec("rootls "+filename);
  //gSystem->Exec("mkdir -p Presel/"+sample);
  gSystem->Exec("mkdir -p HNplots/Presel/charge");

  TCanvas* c1 = new TCanvas("c1",title,200,200,900,800);
  c1->cd();

  TPad* p1 = new TPad("p1", "p1", 0, 0.3, 1, 1);
  p1->SetBottomMargin(0.01);
  p1->Draw();
  p1->cd();
  
  if(sample.Contains("DYTypeI")) hist->SetTitle(title+" #scale[0.5]{: S-channel}");
  else if(sample.Contains("VBFTypeI")) hist->SetTitle(title+" #scale[0.5]{: T-channel}");
  else hist->SetTitle(title+" #scale[0.5]{"+sample+"_"+channel+"}");
  hist->SetStats(0);
  hist->Rebin(rebin);
  if(this_xran1 == "auto") xran1 = hist->GetXaxis()->GetXmin();
  else if(this_xran1 != "auto") xran1 = this_xran1.Atof();
  if(this_xran2 == "auto") xran2 = hist->GetXaxis()->GetXmax();
  else if(this_xran2 != "auto") xran2 = this_xran2.Atof();
  hist->GetXaxis()->SetRangeUser(xran1,xran2);
  if(this_yran == "auto"){
    if(var.Contains("Eta")) yran = hist->GetMaximum()*1.8;
    else yran = hist->GetMaximum()*1.2;
  }
  else if(this_yran != "auto") yran = this_yran.Atof();
  hist->GetYaxis()->SetRangeUser(0,yran);
  hist->SetMarkerStyle(8);
  hist->SetMarkerSize(0.8);
  hist->SetMarkerColor(kRed);
  hist->SetLineColor(kRed);
  hist->SetLineWidth(2);
  //hist->SetLineStyle(7);
  hist->GetXaxis()->SetLabelSize(0);
  hist->Draw("hist ep");
  hist_c->Scale(hist->GetEntries()/hist_c->GetEntries());
  hist_c->Rebin(rebin);
  hist_c->SetMarkerStyle(8);
  hist_c->SetMarkerSize(0.8);
  hist_c->SetMarkerColor(kBlue);
  hist_c->SetLineColor(kBlue);
  hist_c->SetLineWidth(2);
  //hist_c->SetLineStyle(7);
  hist_c->Draw("same hist ep");
  
  TLegend* legend;
  legend = new TLegend(0.7,0.75,0.9,0.9);
  if(sample.Contains("DYTypeI")){
    legend->AddEntry(hist,"SS","l");
    legend->AddEntry(hist_c,"OS","l");
  }
  if(sample.Contains("VBFTypeI")){
    legend->AddEntry(hist,"SS","l");
    legend->AddEntry(hist_c,"OS","l");
  }
  legend->Draw();

  c1->cd();

  TPad* p2 = new TPad("p2", "p2", 0, 0, 1, 0.3);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.2);
  p2->Draw();
  p2->cd();

	TH1D* hist_ratio = (TH1D*)hist->Clone();
  hist_ratio->Divide(hist_c);
  if(var.Contains("Pt") || var.Contains("Mass")) hist_ratio->GetXaxis()->SetTitle("#scale[2.2]{"+title+" [GeV]}");
  else hist_ratio->GetXaxis()->SetTitle("#scale[2.2]{"+title+"}");
  hist_ratio->SetTitle("");
  hist_ratio->SetLineColor(kGreen+1);
  hist_ratio->SetMarkerColor(kGreen+1);
  hist_ratio->GetXaxis()->SetTitleOffset(1.8);
  hist_ratio->GetXaxis()->SetTickLength(0.05);
  hist_ratio->GetXaxis()->SetLabelSize(0.07);
  hist_ratio->GetYaxis()->SetRangeUser(0,2);
  hist_ratio->GetYaxis()->SetLabelSize(0.07);
  hist_ratio->GetYaxis()->SetTitle("#scale[3]{SS/OS}");
  hist_ratio->Draw("hist ep");

  TLine* LineAtOne = new TLine(hist_ratio->GetXaxis()->GetXmin(),1,hist_ratio->GetXaxis()->GetXmax(),1);
  LineAtOne->SetLineStyle(kDashed);
  LineAtOne->SetLineWidth(3);
  LineAtOne->SetLineColor(2);
  LineAtOne->Draw();

  if(SaveAs=="y") c1->SaveAs("HNplots/Presel/charge/"+sample+"_"+channel+"_"+var+".png");
}

}

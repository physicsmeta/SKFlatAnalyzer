void plot_mass(){

ifstream in("list_mass.txt");
string line;

while(getline(in,line)){
  istringstream is(line);
  TString this_line = line;
  if(this_line(0,1)=="#"||this_line=="") continue;

  TString sample, channel, var, title, this_xran1, this_yran, SaveAs;
  double xran1, xran2, yran;
  int rebin;
  is >> sample;
  is >> channel;
  is >> var;
  is >> title;
  is >> this_xran1;
  is >> xran2;
  is >> this_yran;
  is >> rebin;
  is >> SaveAs;

  TString filename1;
  TString filename2;
  TString filename3;
  TFile* file1;
  TFile* file2;
  TFile* file3;
  TH1D* var1;
  TH1D* var2;
  TH1D* var3;
  if(sample.Contains("DYTypeI")){
    filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+"_M100.root";
    filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+"_M500.root";
    filename3 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+"_M1200.root";
    file1 = new TFile(filename1);
    file2 = new TFile(filename2);
    file3 = new TFile(filename3);
    var1 = (TH1D*)file1->Get(channel+"/Pre/"+var+"_Run2");
    var2 = (TH1D*)file2->Get(channel+"/Pre/"+var+"_Run2");
    var3 = (TH1D*)file3->Get(channel+"/Pre/"+var+"_Run2");
  }
  if(sample.Contains("VBFTypeI")){
    filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+"_M500.root";
    filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/Presel/2018/Presel_"+sample+"_M1200.root";
    file1 = new TFile(filename1);
    file2 = new TFile(filename2);
    var1 = (TH1D*)file1->Get(channel+"/Pre/"+var+"_Run2");
    var2 = (TH1D*)file2->Get(channel+"/Pre/"+var+"_Run2");
  }
  
  //gSystem->Exec("rootls "+filename1);
  //gSystem->Exec("mkdir -p Presel/"+sample);
  gSystem->Exec("mkdir -p Presel/mass");

  TCanvas* c1 = new TCanvas("c1",title,200,350,700,650);
  c1->cd();
  
  if(sample.Contains("DYTypeI")) var1->SetTitle(title+" #scale[0.5]{: S-channel}");
  else if(sample.Contains("VBFTypeI")) var1->SetTitle(title+" #scale[0.5]{: T-channel}");
  else var1->SetTitle(title+" #scale[0.5]{"+sample+"_"+channel+"}");
  var1->SetStats(0);
  var1->Rebin(rebin);
  if(this_xran1 == "auto") xran1 = var1->GetXaxis()->GetXmin();
  else if(this_xran1 != "auto") xran1 = this_xran1.Atof();
  var1->GetXaxis()->SetRangeUser(xran1,xran2);
  if(this_yran == "auto"){
    if(var.Contains("Eta")) yran = var1->GetMaximum()*1.8;
    else yran = var1->GetMaximum()*1.2;
  }
  else if(this_yran != "auto") yran = this_yran.Atof();
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
  if(sample.Contains("DYTypeI")){
    var2->SetMarkerColor(kGreen+1);
    var2->SetLineColor(kGreen+1);
  }
  if(sample.Contains("VBFTypeI")){
    var2->SetMarkerColor(kBlue);
    var2->SetLineColor(kBlue);
  }
  var2->SetLineWidth(2);
  //var2->SetLineStyle(7);
  var2->Draw("same hist ep");
  if(sample.Contains("DYTypeI")){
    var3->Scale(var1->GetEntries()/var3->GetEntries());
    var3->Rebin(rebin);
    var3->SetMarkerStyle(8);
    var3->SetMarkerSize(0.8);
    var3->SetMarkerColor(kBlue);
    var3->SetLineColor(kBlue);
    var3->SetLineWidth(2);
    //var3->SetLineStyle(7);
    var3->Draw("same hist ep");
  }
  
  TLegend* legend;
  legend = new TLegend(0.7,0.75,0.9,0.9);
  if(sample.Contains("DYTypeI")){
    legend->AddEntry(var1,"M100","l");
    legend->AddEntry(var2,"M500","l");
    legend->AddEntry(var3,"M1200","l");
  }
  if(sample.Contains("VBFTypeI")){
    legend->AddEntry(var1,"M500","l");
    legend->AddEntry(var2,"M1200","l");
  }
  legend->Draw();

  if(SaveAs=="y") c1->SaveAs("Presel/mass/"+sample+"_"+channel+"_"+var+".png");
}

}

void plot_years(){

ifstream in("list_years.txt");
string line;

while(getline(in,line)){
  istringstream is(line);
  TString this_line = line;
  if(this_line(0,1)=="#"||this_line=="") continue;

  TString type, analyzer, sample, channel, var, title, this_xran1, this_xran2, this_yran, SaveAs, flag;
  double xran1, xran2, yran;
  int rebin;
  is >> type;
  is >> analyzer;
  is >> sample;
  is >> channel;
  is >> var;
  is >> title;
  is >> this_xran1;
  is >> this_xran2;
  is >> this_yran;
  is >> rebin;
  is >> SaveAs;
  is >> flag;

  TString mass = sample(sample.Last('M'),sample.Length());
  TString filename1, filename2, filename3;

  if(type == "MC"){
    filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2016/"+analyzer+"_"+sample+".root";
    filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2017/"+analyzer+"_"+sample+".root";
    filename3 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2018/"+analyzer+"_"+sample+".root";
  }
  else if(type == "DATA"){
    if(channel == "diel"){
      filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2016/DATA/"+analyzer+"_SkimTree_Dilepton_"+sample+".root";
      filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2017/DATA/"+analyzer+"_SkimTree_Dilepton_"+sample+".root";
      filename3 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2018/DATA/"+analyzer+"_SkimTree_Dilepton_EGamma.root";
    }
    else{
      filename1 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2016/DATA/"+analyzer+"_SkimTree_Dilepton_"+sample+".root";
      filename2 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2017/DATA/"+analyzer+"_SkimTree_Dilepton_"+sample+".root";
      filename3 = "/data6/Users/jihkim/SKFlatOutput/Run2Legacy_v4/"+analyzer+"/2018/DATA/"+analyzer+"_SkimTree_Dilepton_"+sample+".root";
    }
  }
  TFile* file1 = new TFile(filename1);
  TFile* file2 = new TFile(filename2);
  TFile* file3 = new TFile(filename3);

  TH1D *var1, *var2, *var3;
  if(analyzer == "Presel"){
    var1 = (TH1D*)file1->Get(channel+"/Pre/"+var+"_Run2");
    var2 = (TH1D*)file2->Get(channel+"/Pre/"+var+"_Run2");
    var3 = (TH1D*)file3->Get(channel+"/Pre/"+var+"_Run2");
  }
  else if(analyzer == "test"){
    var1 = (TH1D*)file1->Get(var);
    var2 = (TH1D*)file2->Get(var);
    var3 = (TH1D*)file3->Get(var);
  }

  //gSystem->Exec("rootls "+filename1);
  //gSystem->Exec("mkdir -p "+analyzer+"/"+sample);
  gSystem->Exec("mkdir -p HNplots/"+analyzer+"/years/"+sample);

  TCanvas* c1 = new TCanvas("c1",title,200,350,700,650);
  c1->cd();
  
  if(sample.Contains("DYTypeI")) var1->SetTitle(title+"#scale[0.5]{ : S-ch, "+channel+", "+mass+"}");
  else if(sample.Contains("VBFTypeI")) var1->SetTitle(title+"#scale[0.5]{ : T-ch, "+channel+", "+mass+"}");
  else var1->SetTitle(title+"#scale[0.5]{ : "+sample+"_"+channel+"}");
  var1->SetStats(0);
  var1->Rebin(rebin);
  if(this_xran1 == "auto") xran1 = var1->GetXaxis()->GetXmin();
  else if(this_xran1 != "auto") xran1 = this_xran1.Atof();
  if(this_xran2 == "auto") xran2 = var1->GetXaxis()->GetXmax();
  else if(this_xran2 != "auto") xran2 = this_xran2.Atof();
  var1->GetXaxis()->SetRangeUser(xran1,xran2);
  if(this_yran == "auto") yran = var1->GetMaximum()*1.2;
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
  legend = new TLegend(0.7,0.75,0.9,0.9);
  legend->AddEntry(var1,"2016","l");
  legend->AddEntry(var2,"2017","l");
  legend->AddEntry(var3,"2018","l");
  legend->Draw();

  if(SaveAs=="y"){
		if(flag != "") c1->SaveAs("HNplots/"+analyzer+"/years/"+sample+"/"+channel+"_"+var+"__"+flag+".png");
		else c1->SaveAs("HNplots/"+analyzer+"/years/"+sample+"/"+channel+"_"+var+".png");
	}
}

}

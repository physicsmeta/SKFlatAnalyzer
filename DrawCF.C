{
TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/ChargeFlipHE__/ChargeFlip_DYJets.root");
//TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/passTightChargeTightIDdXY__/ChargeFlip_DYJets.root");
//TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/passTightChargeTightID__/ChargeFlip_DYJets.root");
//TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/passTightID__/ChargeFlip_DYJets.root");
TH1D* h1 = (TH1D*)f1->Get("ChargeFlip/EtaRegion1_Denom");
TH1D* h2 = (TH1D*)f1->Get("ChargeFlip/EtaRegion1_Num");
TH1D* h3 = (TH1D*)f1->Get("ChargeFlip/EtaRegion2_Denom");
TH1D* h4 = (TH1D*)f1->Get("ChargeFlip/EtaRegion2_Num");
TH1D* h5 = (TH1D*)f1->Get("ChargeFlip/EtaRegion3_Denom");
TH1D* h6 = (TH1D*)f1->Get("ChargeFlip/EtaRegion3_Num");


//double x[40], ex[40];
//for (int i=0; i<40; i++) { x[i] = (2*i+1)*(0.04/80.); ex[i] = 0.04/80.; }
//
//double y_1[40], y_2[40], y_3[40];
//for (int i=0; i<40; i++) { y_1[i] = h2->GetBinContent(i+1)/h1->GetBinContent(i+1); }
//for (int i=0; i<40; i++) { y_2[i] = h4->GetBinContent(i+1)/h3->GetBinContent(i+1); }
//for (int i=0; i<40; i++) { y_3[i] = h6->GetBinContent(i+1)/h5->GetBinContent(i+1); }
//y_3[0]=0; //NaN
//
//double ey_1[40], ey_2[40], ey_3[40];
//for (int i=0; i<40; i++) { ey_1[i] = y_1[i]*(sqrt(pow(h2->GetBinError(i+1)/h2->GetBinContent(i+1),2)+pow(h1->GetBinError(i+1)/h1->GetBinContent(i+1),2))); }
//for (int i=0; i<40; i++) { ey_2[i] = y_2[i]*(sqrt(pow(h4->GetBinError(i+1)/h4->GetBinContent(i+1),2)+pow(h3->GetBinError(i+1)/h3->GetBinContent(i+1),2))); }
//for (int i=0; i<40; i++) { ey_3[i] = y_3[i]*(sqrt(pow(h6->GetBinError(i+1)/h6->GetBinContent(i+1),2)+pow(h5->GetBinError(i+1)/h5->GetBinContent(i+1),2))); } TODO to be deleted


//////////////////////// Let's use vector instead of array, to remove zero bins and nan bins /////////////////////////////////
vector<double> X_1, EX_1, X_2, EX_2, X_3, EX_3;
for (int i=0; i<40; i++) {
  X_1.push_back((2*i+1)*(0.04/80.)); EX_1.push_back(0.04/80.); 
  X_2.push_back((2*i+1)*(0.04/80.)); EX_2.push_back(0.04/80.);
  X_3.push_back((2*i+1)*(0.04/80.)); EX_3.push_back(0.04/80.);
}
vector<double> Y_1, EY_1, Y_2, EY_2, Y_3, EY_3;
for (int i=0; i<40; i++) {
  Y_1.push_back(h2->GetBinContent(i+1)/h1->GetBinContent(i+1));
  EY_1.push_back(Y_1[i]*(sqrt(pow(h2->GetBinError(i+1)/h2->GetBinContent(i+1),2)+pow(h1->GetBinError(i+1)/h1->GetBinContent(i+1),2))));
  Y_2.push_back(h4->GetBinContent(i+1)/h3->GetBinContent(i+1));
  EY_2.push_back(Y_2[i]*(sqrt(pow(h4->GetBinError(i+1)/h4->GetBinContent(i+1),2)+pow(h3->GetBinError(i+1)/h3->GetBinContent(i+1),2))));
  Y_3.push_back(h6->GetBinContent(i+1)/h5->GetBinContent(i+1));
  EY_3.push_back(Y_3[i]*(sqrt(pow(h6->GetBinError(i+1)/h6->GetBinContent(i+1),2)+pow(h5->GetBinError(i+1)/h5->GetBinContent(i+1),2))));
}

for (int i=0; i<X_1.size();) {
  if ( (Y_1.at(i) == 0) || (isnan(Y_1.at(i)) != 0) ) {
		X_1.erase(X_1.begin()+i);
		EX_1.erase(EX_1.begin()+i);
    Y_1.erase(Y_1.begin()+i);
		EY_1.erase(EY_1.begin()+i);
	}
	else i++;
}
for (int i=0; i<X_2.size();) {
  if ( (Y_2.at(i) == 0) || (isnan(Y_2.at(i)) != 0) ) {
		X_2.erase(X_2.begin()+i);
		EX_2.erase(EX_2.begin()+i);
    Y_2.erase(Y_2.begin()+i);
		EY_2.erase(EY_2.begin()+i);
	}
	else i++;
}
for (int i=0; i<X_3.size();) {
  if ( (Y_3.at(i) == 0) || (isnan(Y_3.at(i)) != 0) ) {
		X_3.erase(X_3.begin()+i);
		EX_3.erase(EX_3.begin()+i);
    Y_3.erase(Y_3.begin()+i);
		EY_3.erase(EY_3.begin()+i);
	}
	else i++;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas* c1 = new TCanvas("c1","ChargeFlip_EtaRegion1",200,200,900,800);
TCanvas* c2 = new TCanvas("c2","ChargeFlip_EtaRegion2",250,150,900,800);
TCanvas* c3 = new TCanvas("c3","ChargeFlip_EtaRegion3",300,100,900,800);

c1->cd();
//TGraphErrors* gr1 = new TGraphErrors(40,x,y_1,ex,ey_1);
TGraphErrors* gr1 = new TGraphErrors(X_1.size(),&X_1[0],&Y_1[0],&EX_1[0],&EY_1[0]);
gr1->SetMarkerStyle(20);
//gr1->SetMarkerSize(0.8);
//gr1->SetMarkerColor(kMagenta+2);
gr1->SetLineColor(15);
gr1->SetTitle("ChargeFlip_|#eta|<0.8");
gr1->SetMinimum(0.);
gr1->Draw("ZAP"); //Z : do not draw small horizontal/vertical lines the end of the error bars
gr1->GetXaxis()->SetRangeUser(0.,0.04);
gr1->GetXaxis()->SetTitle("#scale[0.8]{1/p_{T} (GeV^{-1})}");
gr1->GetXaxis()->SetTickLength(0.025);
gr1->GetXaxis()->SetLabelSize(0.025);
gr1->GetYaxis()->SetLabelSize(0.025);

TF1 *gr1_fit1 = new TF1("gr1_fit1","pol1",0,0.021);
TF1 *gr1_fit2 = new TF1("gr1_fit2","pol1",0.021,0.04);

gr1_fit1->SetLineWidth(3);
gr1_fit2->SetLineWidth(3);

gr1_fit1->SetLineColor(4);
gr1_fit2->SetLineColor(4);

cout << "//////////////////// Now fitting on EtaRegion1 ... ////////////////////" << endl;

gr1->Fit(gr1_fit1,"R");

TGraphErrors *gr1_fit1_err = new TGraphErrors(22);
for(int i=0; i<22; i++) gr1_fit1_err->SetPoint(i,0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr1_fit1_err);
gr1_fit1_err->SetFillColor(4);
gr1_fit1_err->SetFillStyle(3001);
gr1_fit1_err->Draw("3 same");

gr1->Fit(gr1_fit2,"R+");

TGraphErrors *gr1_fit2_err = new TGraphErrors(21);
for(int i=0; i<21; i++) gr1_fit2_err->SetPoint(i,0.021+0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr1_fit2_err);
gr1_fit2_err->SetFillColor(4);
gr1_fit2_err->SetFillStyle(3001);
gr1_fit2_err->Draw("3 same");



c2->cd();
//TGraph* gr2 = new TGraph(40,x,y_2);
//TGraphErrors* gr2 = new TGraphErrors(40,x,y_2,ex,ey_2);
TGraphErrors* gr2 = new TGraphErrors(X_2.size(),&X_2[0],&Y_2[0],&EX_2[0],&EY_2[0]);
gr2->SetMarkerStyle(20);
//gr2->SetMarkerSize(0.8);
//gr2->SetMarkerColor(kMagenta+2);
//gr2->SetMarkerColor(kBlue-4);
gr2->SetLineColor(15);
gr2->SetTitle("ChargeFlip_0.8#leq|#eta|<1.4442");
gr2->SetMinimum(0.);
gr2->Draw("ZAP");
gr2->GetXaxis()->SetRangeUser(0.,0.04);
gr2->GetXaxis()->SetTitle("#scale[0.8]{1/p_{T} (GeV^{-1})}");
gr2->GetXaxis()->SetTickLength(0.025);
gr2->GetXaxis()->SetLabelSize(0.025);
gr2->GetYaxis()->SetLabelSize(0.025);

TF1 *gr2_fit1 = new TF1("gr2_fit1","pol1",0,0.0155);
TF1 *gr2_fit2 = new TF1("gr2_fit2","pol1",0.0155,0.023);
TF1 *gr2_fit3 = new TF1("gr2_fit3","pol1",0.023,0.04);

gr2_fit1->SetLineWidth(3);
gr2_fit2->SetLineWidth(3);
gr2_fit3->SetLineWidth(3);

gr2_fit1->SetLineColor(4);
gr2_fit2->SetLineColor(4);
gr2_fit3->SetLineColor(4);

cout << "//////////////////// Now fitting on EtaRegion2 ... ////////////////////" << endl;

gr2->Fit(gr2_fit1,"R");

TGraphErrors *gr2_fit1_err = new TGraphErrors(16);
for(int i=0; i<16; i++) gr2_fit1_err->SetPoint(i,0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr2_fit1_err);
gr2_fit1_err->SetFillColor(4);
gr2_fit1_err->SetFillStyle(3001);
gr2_fit1_err->Draw("3 same");

gr2->Fit(gr2_fit2,"R+");

TGraphErrors *gr2_fit2_err = new TGraphErrors(9);
for(int i=0; i<9; i++) gr2_fit2_err->SetPoint(i,0.015+0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr2_fit2_err);
gr2_fit2_err->SetFillColor(4);
gr2_fit2_err->SetFillStyle(3001);
gr2_fit2_err->Draw("3 same");

gr2->Fit(gr2_fit3,"R+");

TGraphErrors *gr2_fit3_err = new TGraphErrors(19);
for(int i=0; i<19; i++) gr2_fit3_err->SetPoint(i,0.023+0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr2_fit3_err);
gr2_fit3_err->SetFillColor(4);
gr2_fit3_err->SetFillStyle(3001);
gr2_fit3_err->Draw("3 same");



c3->cd();
//TGraph* gr3 = new TGraph(40,x,y_3);
//TGraphErrors* gr3 = new TGraphErrors(40,x,y_3,ex,ey_3);
TGraphErrors* gr3 = new TGraphErrors(X_3.size(),&X_3[0],&Y_3[0],&EX_3[0],&EY_3[0]);
gr3->SetMarkerStyle(20);
//gr3->SetMarkerSize(0.8);
//gr3->SetMarkerColor(kMagenta+2);
//gr3->SetMarkerColor(kGreen+3);
gr3->SetLineColor(15);
gr3->SetTitle("ChargeFlip_1.556#leq|#eta|<2.5");
gr3->Draw("ZAP");
gr3->GetXaxis()->SetRangeUser(0.,0.04);
gr3->GetXaxis()->SetTitle("#scale[0.8]{1/p_{T} (GeV^{-1})}");
gr3->GetXaxis()->SetTickLength(0.025);
gr3->GetXaxis()->SetLabelSize(0.025);
gr3->GetYaxis()->SetLabelSize(0.025);

TF1 *gr3_fit1 = new TF1("gr3_fit1","pol1",0,0.0105);
TF1 *gr3_fit2 = new TF1("gr3_fit2","pol1",0.0105,0.02);
TF1 *gr3_fit3 = new TF1("gr3_fit3","pol1",0.02,0.04);

gr3_fit1->SetLineWidth(3);
gr3_fit2->SetLineWidth(3);
gr3_fit3->SetLineWidth(3);

gr3_fit1->SetLineColor(4);
gr3_fit2->SetLineColor(4);
gr3_fit3->SetLineColor(4);

cout << "//////////////////// Now fitting on EtaRegion3 ... ////////////////////" << endl;

gr3->Fit(gr3_fit1,"R");

TGraphErrors *gr3_fit1_err = new TGraphErrors(12);
for(int i=0; i<12; i++) gr3_fit1_err->SetPoint(i,0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr3_fit1_err);
gr3_fit1_err->SetFillColor(4);
gr3_fit1_err->SetFillStyle(3001);
gr3_fit1_err->Draw("3 same");

gr3->Fit(gr3_fit2,"R+");

TGraphErrors *gr3_fit2_err = new TGraphErrors(10);
for(int i=0; i<10; i++) gr3_fit2_err->SetPoint(i,0.011+0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr3_fit2_err);
gr3_fit2_err->SetFillColor(4);
gr3_fit2_err->SetFillStyle(3001);
gr3_fit2_err->Draw("3 same");

gr3->Fit(gr3_fit3,"R+");

TGraphErrors *gr3_fit3_err = new TGraphErrors(22);
for(int i=0; i<22; i++) gr3_fit3_err->SetPoint(i,0.02+0.001*i,0);
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr3_fit3_err);
gr3_fit3_err->SetFillColor(4);
gr3_fit3_err->SetFillStyle(3001);
gr3_fit3_err->Draw("3 same");


//double par[6];
//gr3_fit1->GetParameters(&par[0]);
//gr3_fit2->GetParameters(&par[2]);
//gr3_fit3->GetParameters(&par[4]);
//
//TF1 *gr3_fit4 = new TF1("gr3_fit4","CFfit",0,0.04,6);
//gr3_fit4->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5]);
//gr3_fit4->SetLineColor(4);
//gr3_fit4->Draw("SAME");

}


double CFfit(double *x, double *par)
{
  double xx = x[0];
  double f = (xx<0.012)*(par[0]+xx*par[1]) + (xx>=0.012 && xx<0.021)*(par[2]+xx*par[3]) + (xx>=0.021 && xx<0.04)*(par[4]+xx*par[5]);
	return f;
}

{
TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/passTightChargeTightID__/ChargeFlip_DYJets.root");
//TFile* f1 = new TFile("/data4/Users/jihkim/SKFlatOutput/Run2Legacy_v3/ChargeFlip/2016/passTightID__/ChargeFlip_DYJets.root");
TH1D* h1 = (TH1D*)f1->Get("ChargeFlip/EtaRegion1_Denom");
TH1D* h2 = (TH1D*)f1->Get("ChargeFlip/EtaRegion1_Num");
TH1D* h3 = (TH1D*)f1->Get("ChargeFlip/EtaRegion2_Denom");
TH1D* h4 = (TH1D*)f1->Get("ChargeFlip/EtaRegion2_Num");
TH1D* h5 = (TH1D*)f1->Get("ChargeFlip/EtaRegion3_Denom");
TH1D* h6 = (TH1D*)f1->Get("ChargeFlip/EtaRegion3_Num");


double x[40], ex[40];
for (int i=0; i<40; i++) { x[i] = (2*i+1)*(0.04/80.); ex[i] = 0.04/80.; }

double y_1[40], y_2[40], y_3[40];
for (int i=0; i<40; i++) { y_1[i] = h2->GetBinContent(i+1)/h1->GetBinContent(i+1); }
for (int i=0; i<40; i++) { y_2[i] = h4->GetBinContent(i+1)/h3->GetBinContent(i+1); }
for (int i=0; i<40; i++) { y_3[i] = h6->GetBinContent(i+1)/h5->GetBinContent(i+1); }
y_3[0]=0; //NaN

double ey_1[40], ey_2[40], ey_3[40];
for (int i=0; i<40; i++) { ey_1[i] = y_1[i]*(sqrt(pow(h2->GetBinError(i+1)/h2->GetBinContent(i+1),2)+pow(h1->GetBinError(i+1)/h1->GetBinContent(i+1),2))); }
for (int i=0; i<40; i++) { ey_2[i] = y_2[i]*(sqrt(pow(h4->GetBinError(i+1)/h4->GetBinContent(i+1),2)+pow(h3->GetBinError(i+1)/h3->GetBinContent(i+1),2))); }
for (int i=0; i<40; i++) { ey_3[i] = y_3[i]*(sqrt(pow(h6->GetBinError(i+1)/h6->GetBinContent(i+1),2)+pow(h5->GetBinError(i+1)/h5->GetBinContent(i+1),2))); }


TCanvas* c1 = new TCanvas("c1","ChargeFlip_EtaRegion1",200,200,900,800);
TCanvas* c2 = new TCanvas("c2","ChargeFlip_EtaRegion2",250,150,900,800);
TCanvas* c3 = new TCanvas("c3","ChargeFlip_EtaRegion3",300,100,900,800);

c1->cd();
TGraphErrors* gr1 = new TGraphErrors(40,x,y_1,ex,ey_1);
gr1->SetMarkerStyle(20);
//gr1->SetMarkerSize(0.8);
//gr1->SetMarkerColor(kMagenta+2);
gr1->SetLineColor(15);
gr1->SetTitle("ChargeFlip_|#eta|<0.8");
gr1->Draw("ZAP"); //Z : do not draw small horizontal/vertical lines the end of the error bars
gr1->GetXaxis()->SetRangeUser(0.,0.04);
gr1->GetXaxis()->SetTitle("#scale[0.8]{1/p_{T} (GeV^{-1})}");
gr1->GetXaxis()->SetTickLength(0.025);
gr1->GetXaxis()->SetLabelSize(0.025);
gr1->GetYaxis()->SetLabelSize(0.025);

c2->cd();
//TGraph* gr2 = new TGraph(40,x,y_2);
TGraphErrors* gr2 = new TGraphErrors(40,x,y_2,ex,ey_2);
gr2->SetMarkerStyle(20);
//gr2->SetMarkerSize(0.8);
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

c3->cd();
//TGraph* gr3 = new TGraph(40,x,y_3);
TGraphErrors* gr3 = new TGraphErrors(40,x,y_3,ex,ey_3);
gr3->SetMarkerStyle(20);
//gr3->SetMarkerSize(0.8);
//gr3->SetMarkerColor(kGreen+3);
gr3->SetLineColor(15);
gr3->SetTitle("ChargeFlip_1.556#leq|#eta|<2.5");
gr3->Draw("ZAP");
gr3->GetXaxis()->SetRangeUser(0.,0.04);
gr3->GetXaxis()->SetTitle("#scale[0.8]{1/p_{T} (GeV^{-1})}");
gr3->GetXaxis()->SetTickLength(0.025);
gr3->GetXaxis()->SetLabelSize(0.025);
gr3->GetYaxis()->SetLabelSize(0.025);
}

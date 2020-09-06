void Acceptance(){

TString workdir = "/data6/Users/jihkim/SKFlatOutput";
TString SKFlatVersion = "Run2Legacy_v4";
TString analyzer = "Signal_opt";
TString year = "2016";
vector<TString> mass = {"M700","M1000","M1500"};
vector<TString> flag = {"M700__","M1000__","M1500__pp__"};
TString file_path = "";

vector<TString> ID = {"HN16"};

int maxBinNumber_temp = 0;

gSystem->Exec("mkdir -p Acceptance/Signal/2016/");
ofstream Mean_Sigma("Acceptance/Signal/2016/Mean_Sigma_comparison.txt");
ofstream acceptance("Acceptance/Signal/2016/acceptance.txt");
ofstream cutflow_eff("Acceptance/Signal/2016/cutflow_eff.txt");

  for(unsigned int i=0; i<flag.size(); i++){
  
    TFile *f_on, *f_off;
    TH1D *h_on_SR1_cutflow, *h_off_SR1_cutflow, *h_on_SR2_cutflow, *h_off_SR2_cutflow, *h_Temp, *h_on_tot, *h_off_tot, *h_on_final, *h_off_final, *h_on_acc, *h_off_acc;
  
    file_path = workdir+"/"+SKFlatVersion+"/"+analyzer+"/"+year+"/"+flag.at(i);
  
    f_on = new TFile(file_path+"/Signal_opt_HNToMuMu_Tch_"+mass.at(i)+".root");
    f_off = new TFile(file_path+"/Signal_opt_Last_HNToMuMu_Tch_"+mass.at(i)+".root");
  
    h_on_SR1_cutflow = (TH1D*)f_on->Get("dimu/highSR1/Number_Events_HN16");
    h_on_SR2_cutflow = (TH1D*)f_on->Get("dimu/highSR2/Number_Events_HN16");
    h_off_SR1_cutflow = (TH1D*)f_off->Get("dimu/highSR1/Number_Events_HN16");
    h_off_SR2_cutflow = (TH1D*)f_off->Get("dimu/highSR2/Number_Events_HN16");

    h_Temp = (TH1D*)h_on_SR1_cutflow->Clone();
    maxBinNumber_temp = h_Temp->GetNbinsX();
    for(int j=0; j<maxBinNumber_temp+2; j++){
      h_Temp->SetBinContent(i, 0.);
      h_Temp->SetBinError(i, 0.);
    }

    h_on_tot = (TH1D*)h_Temp->Clone();
    h_on_final = (TH1D*)h_Temp->Clone();
    h_on_acc = (TH1D*)h_Temp->Clone();
    h_off_tot = (TH1D*)h_Temp->Clone();
    h_off_final = (TH1D*)h_Temp->Clone();
    h_off_acc = (TH1D*)h_Temp->Clone();

    h_on_tot->SetBinContent(1,h_on_SR1_cutflow->GetBinContent(0));
    h_on_tot->SetBinError(1,h_on_SR1_cutflow->GetBinError(0));
    h_on_tot->SetBinContent(2,h_on_SR2_cutflow->GetBinContent(0));
    h_on_tot->SetBinError(2,h_on_SR2_cutflow->GetBinError(0));
  
    h_on_final->SetBinContent(1,h_on_SR1_cutflow->GetBinContent(15));
    h_on_final->SetBinError(1,h_on_SR1_cutflow->GetBinError(15));
    h_on_final->SetBinContent(2,h_on_SR2_cutflow->GetBinContent(11));
    h_on_final->SetBinError(2,h_on_SR2_cutflow->GetBinError(11));
  
    h_off_tot->SetBinContent(1,h_off_SR1_cutflow->GetBinContent(0));
    h_off_tot->SetBinError(1,h_off_SR1_cutflow->GetBinError(0));
    h_off_tot->SetBinContent(2,h_off_SR2_cutflow->GetBinContent(0));
    h_off_tot->SetBinError(2,h_off_SR2_cutflow->GetBinError(0));
  
    h_off_final->SetBinContent(1,h_off_SR1_cutflow->GetBinContent(15));
    h_off_final->SetBinError(1,h_off_SR1_cutflow->GetBinError(15));
    h_off_final->SetBinContent(2,h_off_SR2_cutflow->GetBinContent(11));
    h_off_final->SetBinError(2,h_off_SR2_cutflow->GetBinError(11));

    for(int j=1;j<3;j++){
      h_on_acc->SetBinContent(j,0.);
      h_on_acc->SetBinError(j,0.);
      h_off_acc->SetBinContent(j,0.);
      h_off_acc->SetBinError(j,0.);
    }

    h_on_acc->Divide(h_on_final,h_on_tot,1.,1.,"B");
    h_off_acc->Divide(h_off_final,h_off_tot,1.,1.,"B");

    acceptance << "mass point : " << mass.at(i) << "\n";
    acceptance << "QEDon acceptance (SR1) : " << 100*h_on_acc->GetBinContent(1) << " +- " << 100*h_on_acc->GetBinError(1) << " %" << "\n";
    acceptance << "QEDon acceptance (SR2) : " << 100*h_on_acc->GetBinContent(2) << " +- " << 100*h_on_acc->GetBinError(2) << " %" << "\n";
    acceptance << "QEDoff acceptance (SR1) : " << 100*h_off_acc->GetBinContent(1) << " +- " << 100*h_off_acc->GetBinError(1) << " %" << "\n";
    acceptance << "QEDoff acceptance (SR2) : " << 100*h_off_acc->GetBinContent(2) << " +- " << 100*h_off_acc->GetBinError(2) << " %" << "\n";
    acceptance << "\n";
    cout << "mass point : " << mass.at(i) << endl;
    cout << "QEDon acceptance (SR1) : " << 100*h_on_acc->GetBinContent(1) << " +- " << 100*h_on_acc->GetBinError(1) << " %" << endl;
    cout << "QEDon acceptance (SR2) : " << 100*h_on_acc->GetBinContent(2) << " +- " << 100*h_on_acc->GetBinError(2) << " %" << endl;
    cout << "QEDoff acceptance (SR1) : " << 100*h_off_acc->GetBinContent(1) << " +- " << 100*h_off_acc->GetBinError(1) << " %" << endl;
    cout << "QEDoff acceptance (SR2) : " << 100*h_off_acc->GetBinContent(2) << " +- " << 100*h_off_acc->GetBinError(2) << " %" << endl;

    //Cutflow efficiency

    cutflow_eff << "===================="+mass.at(i)+"====================" << "\n";

    TH1D *h2_on_SR1_cutflow = new TH1D("SR1","cutflow",16,0,16);
    TH1D *h2_off_SR1_cutflow = new TH1D("SR1","cutflow",16,0,16);
    TH1D *h2_on_SR2_cutflow = new TH1D("SR2","cutflow",12,0,12);
    TH1D *h2_off_SR2_cutflow = new TH1D("SR2","cutflow",12,0,12);
    for(int j=0; j<maxBinNumber_temp+1; j++){
      h2_on_SR1_cutflow->SetBinContent(j+1,h_on_SR1_cutflow->GetBinContent(j)/h_on_SR1_cutflow->GetBinContent(0));
      //h2_on_SR1_cutflow->GetXaxis()->SetBinLabel(j+1,"0");
      h2_off_SR1_cutflow->SetBinContent(j+1,h_off_SR1_cutflow->GetBinContent(j)/h_off_SR1_cutflow->GetBinContent(0));
      //h2_off_SR1_cutflow->GetXaxis()->SetBinLabel(j+1,"0");
      h2_on_SR2_cutflow->SetBinContent(j+1,h_on_SR2_cutflow->GetBinContent(j)/h_on_SR2_cutflow->GetBinContent(0));
      //h2_on_SR1_cutflow->GetXaxis()->SetBinLabel(j+1,"0");
      h2_off_SR2_cutflow->SetBinContent(j+1,h_off_SR2_cutflow->GetBinContent(j)/h_off_SR2_cutflow->GetBinContent(0));
      //h2_off_SR1_cutflow->GetXaxis()->SetBinLabel(j+1,"0");
    }
   
    TCanvas *c1 = new TCanvas("c1","c1",200,200,900,800);
    c1->SetLogy();
   
    TCanvas *c2 = new TCanvas("c2","c2",300,300,1000,900);
    c2->SetLogy();

    TString title1 = "highmass_SR1_Cutflow";
    TString title2 = "highmass_SR2_Cutflow";

    c1->cd();

    h2_on_SR1_cutflow->SetTitle(title1);
    h2_on_SR1_cutflow->SetStats(0);
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(1,15,0.02,-1,-1,-1,"Nocut");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(2,15,0.02,-1,-1,-1,"MET filter");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(3,15,0.02,-1,-1,-1,"Dilepton Trigger");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(4,15,0.02,-1,-1,-1,"2 tight leptons");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(5,15,0.02,-1,-1,-1,"same sign");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(6,15,0.02,-1,-1,-1,"3rd lepton veto");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(7,15,0.02,-1,-1,-1,"m(ll) > 10GeV");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(8,15,0.02,-1,-1,-1,"jet selection (pre)");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(9,15,0.02,-1,-1,-1,"signal selection");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(10,15,0.02,-1,-1,-1,"Njet");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(11,15,0.02,-1,-1,-1,"ptj1");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(12,15,0.02,-1,-1,-1,"ptl1");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(13,15,0.02,-1,-1,-1,"M(jj)");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(14,15,0.02,-1,-1,-1,"M(lljj)");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(15,15,0.02,-1,-1,-1,"M(ljj)");
    //h2_on_SR1_cutflow->GetXaxis()->ChangeLabel(16,15,0.02,-1,-1,-1,"MET2ST");
    h2_on_SR1_cutflow->GetYaxis()->SetLabelSize(0.03);
    h2_on_SR1_cutflow->GetYaxis()->SetTitleOffset(1);
    h2_on_SR1_cutflow->GetYaxis()->SetTitle("Efficiency");
    h2_on_SR1_cutflow->GetYaxis()->SetRangeUser(1.e-02,5);
    h2_on_SR1_cutflow->SetLineColor(kRed);
    h2_on_SR1_cutflow->Draw();
   
    h2_off_SR1_cutflow->SetTitle("");
    h2_off_SR1_cutflow->SetStats(0);
    h2_off_SR1_cutflow->SetLineColor(kBlue);
    h2_off_SR1_cutflow->Draw("same");

    cutflow_eff << "//high mass SR1//\n";
    cutflow_eff << "QEDon\t\tQEDoff\n";
    for(int j=0; j<16; j++){
      TString i_t_1 = Form("%d",j+1);
      TString item_t_on_SR1 = Form("%f",100*h2_on_SR1_cutflow->GetBinContent(j+1));
      TString item_t_off_SR1 = Form("%f",100*h2_off_SR1_cutflow->GetBinContent(j+1));
      cutflow_eff << i_t_1+"th eff. : "+item_t_on_SR1+"%\t\t"+item_t_off_SR1+"%\n";
    }

    TLegend* legend1 = new TLegend(0.75,0.8,0.9,0.9);
    legend1->AddEntry(h2_on_SR1_cutflow,"QEDon","lp");
    legend1->AddEntry(h2_off_SR1_cutflow,"QEDoff","lp");
    legend1->Draw();


    c2->cd();

    h2_on_SR2_cutflow->SetTitle(title2);
    h2_on_SR2_cutflow->SetStats(0);
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(1,15,0.02,-1,-1,-1,"Nocut");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(2,15,0.02,-1,-1,-1,"MET filter");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(3,15,0.02,-1,-1,-1,"Dilepton Trigger");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(4,15,0.02,-1,-1,-1,"2 tight leptons");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(5,15,0.02,-1,-1,-1,"same sign");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(6,15,0.02,-1,-1,-1,"3rd lepton veto");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(7,15,0.02,-1,-1,-1,"m(ll) > 10GeV");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(8,15,0.02,-1,-1,-1,"jet selection (pre)");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(9,15,0.02,-1,-1,-1,"signal selection");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(10,15,0.02,-1,-1,-1,"ptl1");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(11,15,0.02,-1,-1,-1,"M(J)");
    //h2_on_SR2_cutflow->GetXaxis()->ChangeLabel(12,15,0.02,-1,-1,-1,"M(lJ)");
    h2_on_SR2_cutflow->GetYaxis()->SetLabelSize(0.03);
    h2_on_SR2_cutflow->GetYaxis()->SetTitleOffset(1);
    h2_on_SR2_cutflow->GetYaxis()->SetTitle("Efficiency");
    h2_on_SR2_cutflow->GetYaxis()->SetRangeUser(1.e-02,5);
    h2_on_SR2_cutflow->SetLineColor(kRed);
    h2_on_SR2_cutflow->Draw();
   
    h2_off_SR2_cutflow->SetTitle("");
    h2_off_SR2_cutflow->SetStats(0);
    h2_off_SR2_cutflow->SetLineColor(kBlue);
    h2_off_SR2_cutflow->Draw("same");

    cutflow_eff << "//high mass SR2//\n";
    cutflow_eff << "QEDon\t\tQEDoff\n";
    for(int i=0; i<12; i++){
      TString i_t_2 = Form("%d",i+1);
      TString item_t_on_SR2 = Form("%f",100*h2_on_SR2_cutflow->GetBinContent(i+1));
      TString item_t_off_SR2 = Form("%f",100*h2_off_SR2_cutflow->GetBinContent(i+1));
      cutflow_eff << i_t_2+"th eff. : "+item_t_on_SR2+"%\t\t"+item_t_off_SR2+"%\n";
    }
    cutflow_eff << "\n";

    TLegend* legend2 = new TLegend(0.75,0.8,0.9,0.9);
    legend2->AddEntry(h2_on_SR2_cutflow,"QEDon","lp");
    legend2->AddEntry(h2_off_SR2_cutflow,"QEDoff","lp");
    legend2->Draw();

    for(unsigned int j=0; j<ID.size(); j++){ 
      gSystem->Exec("mkdir -p Cutflow/Signal/2016/"+flag.at(i)+"/"+ID.at(j));
      c1->SaveAs("Cutflow/Signal/2016/"+flag.at(i)+"/"+ID.at(j)+"/"+title1+".png");
      c2->SaveAs("Cutflow/Signal/2016/"+flag.at(i)+"/"+ID.at(j)+"/"+title2+".png");
    }

    delete c1;
    delete h2_on_SR1_cutflow;
    delete h2_on_SR2_cutflow;
    delete h2_off_SR1_cutflow;
    delete h2_off_SR2_cutflow;


    //Now compare mean and SD

    //vector<TString> SR1_variables = {"Lep1_Pt","Lep1_Eta","Lep2_Pt","Lep2_Eta","MET","MET2ST","Njets","WCand_Mass","ZCand_Mass","ZCand_Pt","l1jj_Mass","l2jj_Mass","lljj_Mass","ptj1"};
    //vector<TString> SR2_variables = {"Lep1_Pt","Lep1_Eta","Lep2_Pt","Lep2_Eta","MET","MET2ST","Fatjet_Mass","ZCand_Mass","ZCand_Pt","l1J_Mass"};
    vector<TString> SR1_variables = {"Lep1_Pt","MET2ST","Njets","WCand_Mass","l1jj_Mass","lljj_Mass","ptj1"};
    vector<TString> SR2_variables = {"Lep1_Pt","MET2ST","Fatjet_Mass","l1J_Mass"};

    Mean_Sigma << "==============================="+mass.at(i)+"==============================" << "\n";
    Mean_Sigma << "//////////high mass SR1//////////" << "\n";
    cout << "=============================="+mass.at(i)+"==============================" << endl;
    cout << "//////////high mass SR1//////////" << endl;
    for(unsigned int j=0; j<SR1_variables.size(); j++){
      TH1D* h_on_SR1_variable = (TH1D*)f_on->Get("dimu/highSR1/Final/"+SR1_variables.at(j)+"_HN16");
      TH1D* h_off_SR1_variable = (TH1D*)f_off->Get("dimu/highSR1/Final/"+SR1_variables.at(j)+"_HN16");
      Mean_Sigma << "#### "+SR1_variables.at(j)+" ####" << "\n";
      Mean_Sigma << "QEDon mean : " << h_on_SR1_variable->GetMean() << ", sigma : " << h_on_SR1_variable->GetStdDev() << "\n";
      Mean_Sigma << "QEDoff mean : " << h_off_SR1_variable->GetMean() << ", sigma : " << h_off_SR1_variable->GetStdDev() << "\n";
      Mean_Sigma << "mean diff : " << 100*(h_on_SR1_variable->GetMean()-h_off_SR1_variable->GetMean())/h_off_SR1_variable->GetMean() << "%" << "\n";
      cout << "#### "+SR1_variables.at(j)+" ####" << endl;
      cout << "QEDon mean : " << h_on_SR1_variable->GetMean() << ", sigma : " << h_on_SR1_variable->GetStdDev() << endl;
      cout << "QEDoff mean : " << h_off_SR1_variable->GetMean() << ", sigma : " << h_off_SR1_variable->GetStdDev() << endl;
      cout << "mean diff : " << 100*(h_on_SR1_variable->GetMean()-h_off_SR1_variable->GetMean())/h_off_SR1_variable->GetMean() << "%" << endl;
    }
    Mean_Sigma << "\n";
    Mean_Sigma << "///////////high mass SR2//////////" << "\n";
    cout << endl;
    cout << "//////////high mass SR2///////////" << endl;
    for(unsigned int j=0; j<SR2_variables.size(); j++){
      TH1D* h_on_SR2_variable = (TH1D*)f_on->Get("dimu/highSR2/Final/"+SR2_variables.at(j)+"_HN16");
      TH1D* h_off_SR2_variable = (TH1D*)f_off->Get("dimu/highSR2/Final/"+SR2_variables.at(j)+"_HN16");
      Mean_Sigma << "#### "+SR2_variables.at(j)+" ####" << "\n";
      Mean_Sigma << "QEDon mean : " << h_on_SR2_variable->GetMean() << ", sigma : " << h_on_SR2_variable->GetStdDev() << "\n";
      Mean_Sigma << "QEDoff mean : " << h_off_SR2_variable->GetMean() << ", sigma : " << h_off_SR2_variable->GetStdDev() << "\n";
      Mean_Sigma << "mean diff : " << 100*(h_on_SR2_variable->GetMean()-h_off_SR2_variable->GetMean())/h_off_SR2_variable->GetMean() << "%" << "\n";
      cout << "#### "+SR2_variables.at(j)+" ####" << endl;
      cout << "QEDon mean : " << h_on_SR2_variable->GetMean() << ", sigma : " << h_on_SR2_variable->GetStdDev() << endl;
      cout << "QEDoff mean : " << h_off_SR2_variable->GetMean() << ", sigma : " << h_off_SR2_variable->GetStdDev() << endl;
      cout << "mean diff : " << 100*(h_on_SR2_variable->GetMean()-h_off_SR2_variable->GetMean())/h_off_SR2_variable->GetMean() << "%" << endl;
    }
    Mean_Sigma << "\n";
    cout << endl;

  }


}

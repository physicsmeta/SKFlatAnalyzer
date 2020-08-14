{
	TFile *file = TFile::Open("/gv0/DATA/SKFlat/Run2Legacy_v4/2016/PrivateMC/HNtypeI/HeavyMajoranaNeutrinoToDiLepton_Schannel_NLO_MuMu_M90/HeavyMajoranaNeutrinoToDiLepton_Schannel_NLO_MuMu_M90_2016_Ntuple.root");
	TTree *tree = (TTree*)file->Get("recoTree/SKFlat");

	vector<string> *HLT_TriggerName;
	vector<string> *HLT_TriggerFilterName;
	vector<double> *HLTObject_pt;
	vector<double> *HLTObject_eta;
	vector<double> *HLTObject_phi;
	vector<string> *HLTObject_FiredFilters;
	vector<string> *HLTObject_FiredPaths;
	tree->SetBranchAddress("HLT_TriggerName", &HLT_TriggerName);
	tree->SetBranchAddress("HLT_TriggerFilterName", &HLT_TriggerFilterName);
	tree->SetBranchAddress("HLTObject_pt", &HLTObject_pt);
	tree->SetBranchAddress("HLTObject_eta", &HLTObject_eta);
	tree->SetBranchAddress("HLTObject_phi", &HLTObject_phi);
	tree->SetBranchAddress("HLTObject_FiredFilters", &HLTObject_FiredFilters);
	tree->SetBranchAddress("HLTObject_FiredPaths", &HLTObject_FiredPaths);

	int Nevents = tree->GetEntries();

	for(int i=0; i<Nevents; i++){
		tree->GetEntry(i);
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << i << "th event!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "**********All passed triggers: " << HLT_TriggerName->size() << "**********" << endl;
    for(unsigned int j=0; j<HLT_TriggerName->size(); j++){
      //cout << HLT_TriggerName->at(j) << endl;
			//TString this_TriggerName = HLT_TriggerName->at(j);
			//int idx = this_TriggerName.Index("HLT_Mu17_Mu8_DZ_v");
			//if (idx>=0&&this_TriggerName!="HLT_Mu17_Mu8_DZ_v7"){
		  //  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << i << "th event!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;		 
			//	cout << this_TriggerName << endl;
			//}
    }
    cout << "**********All (interested) triggers final filter: " << HLT_TriggerFilterName->size() << "**********" << endl;
		if (HLT_TriggerFilterName->size() != 141) return;
    for(unsigned int j=0; j<HLT_TriggerFilterName->size(); j++){
      cout << HLT_TriggerFilterName->at(j)	<< endl;
    }
    cout << "**********All trigger objects: " << HLTObject_pt->size() << "**********" << endl;
    for(unsigned int j=0; j<HLTObject_pt->size(); j++){
      cout << "=====" << j << "th object=====" << endl;
      cout << "pt: " << HLTObject_pt->at(j) << ", eta: " << HLTObject_eta->at(j) << ", phi: " << HLTObject_phi->at(j) << endl;
      cout << "fired " << " filters = " << HLTObject_FiredFilters->at(j) << endl;
      cout << "fired " << " paths = " << HLTObject_FiredPaths->at(j) << endl;
			if (HLTObject_FiredFilters->size() != HLTObject_FiredPaths->size()) return;
    }
	}

}

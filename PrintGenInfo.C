{
	TFile *file = TFile::Open("/data8/Users/jihkim/GeneratorTools/external/CMSSW_10_2_18/src/SKFlatMaker/SKFlatMaker/test/SKFlatNtuple.root");
	TTree *tree = (TTree*)file->Get("recoTree/SKFlat");

  vector<double>  *gen_phi;
  vector<double>  *gen_eta;
  vector<double>  *gen_pt;
  vector<double>  *gen_mass;
  vector<double>  *gen_charge;
  vector<int>     *gen_mother_index;
  vector<int>     *gen_mother2_index;
  vector<int>     *gen_daughter_index;
  vector<int>     *gen_daughter2_index;
  vector<int>     *gen_status;
  vector<int>     *gen_PID;
	tree->SetBranchAddress("gen_mass", &gen_mass);
	tree->SetBranchAddress("gen_charge", &gen_charge);
	tree->SetBranchAddress("gen_pt", &gen_pt);
	tree->SetBranchAddress("gen_eta", &gen_eta);
	tree->SetBranchAddress("gen_phi", &gen_phi);
	tree->SetBranchAddress("gen_mother_index", &gen_mother_index);
	tree->SetBranchAddress("gen_mother2_index", &gen_mother2_index);
	tree->SetBranchAddress("gen_daughter_index", &gen_daughter_index);
	tree->SetBranchAddress("gen_daughter2_index", &gen_daughter2_index);
	tree->SetBranchAddress("gen_status", &gen_status);
	tree->SetBranchAddress("gen_PID", &gen_PID);

	int Nevents = tree->GetEntries();

	for(int i=0; i<Nevents; i++){
		tree->GetEntry(i);
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << i << "th event!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "index\tPID\tStatus\tMIdx\tM2Idx\tDIdx\tD2IDx\tPt\tEta\tPhi\tM" << endl;
    for(int j=0; j<gen_mass->size(); j++){
      cout << j << "\t" << gen_PID->at(j) << "\t" << gen_status->at(j) << "\t" << gen_mother_index->at(j) << "\t" << gen_mother2_index->at(j) << "\t" << gen_daughter_index->at(j) << "\t" << gen_daughter2_index->at(j) << "\t" << gen_pt->at(j) << "\t" << gen_eta->at(j) << "\t" << gen_phi->at(j) << "\t" << gen_mass->at(j) << endl;
    }
	}

}

#include "runMVA.h"

void b1j4(std::string imva)
{   
   // tH Hut vs ttbar
   if( imva == "HutSTallb1j4" || imva == "HutSTeb1j4" || imva == "HutSTmb1j4" )
     {
	std::string outfileName = "TMVA_HutST_all_b1j4";
	if( imva == "HutSTeb1j4" ) outfileName = "TMVA_HutST_e_b1j4";
	else if( imva == "HutSTmb1j4" ) outfileName = "TMVA_HutST_m_b1j4";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;
	
	double HiggsMass_TOPHLEPBB;
	double TopHadMass_TOPTOPLEPHAD;
	double MVA_TOPHLEPBB;
	double MVA_TOPTOPLEPHAD;
	int LepCharge;
	double HiggsEta_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double TopLepMass_TOPTOPLEPHAD;
	double TopLepPt_TOPHLEPBB;
	double TopLepEta_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double TopLepTopHadDr_TOPTOPLEPHAD;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPTOPLEPHAD;
	double TopHadBJetCSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("TopHadMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	factory->AddVariable("LepCharge",'I');
	factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');

	std::string trName = "tr_all_b1j4";
	if( imva == "HutSTeb1j4" ) trName = "tr_e_b1j4";
	else if( imva == "HutSTmb1j4" ) trName = "tr_m_b1j4";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED_36/ST_TH_1L3B.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED_36/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ZZ_TuneCUETP8M1_13TeV-pythia8.root");

	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("LepCharge",&LepCharge);
	tr_sig.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("LepCharge",&LepCharge);
	tr_bkg.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);

	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(20);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = TopHadMass_TOPTOPLEPHAD;
	     vars[2] = MVA_TOPHLEPBB;
	     vars[3] = MVA_TOPTOPLEPHAD;
	     vars[4] = LepCharge;
	     vars[5] = HiggsEta_TOPHLEPBB;
	     vars[6] = TopLepMass_TOPHLEPBB;
	     vars[7] = TopLepMass_TOPTOPLEPHAD;
	     vars[8] = TopLepPt_TOPHLEPBB;
	     vars[9] = TopLepEta_TOPHLEPBB;
	     vars[10] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[11] = TopLepHiggsDr_TOPHLEPBB;
	     vars[12] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[13] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[14] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[15] = TopLepBJetCSVv2_TOPHLEPBB;
	     vars[16] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[17] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[18] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[19] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddSignalTrainingEvent(vars,weight);
	     else
	       factory->AddSignalTestEvent(vars,weight);
	  }
	
	for(int i=0;i<tr_bkg.GetEntries();i++)
	  {
	     tr_bkg.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = TopHadMass_TOPTOPLEPHAD;
	     vars[2] = MVA_TOPHLEPBB;
	     vars[3] = MVA_TOPTOPLEPHAD;
	     vars[4] = LepCharge;
	     vars[5] = HiggsEta_TOPHLEPBB;
	     vars[6] = TopLepMass_TOPHLEPBB;
	     vars[7] = TopLepMass_TOPTOPLEPHAD;
	     vars[8] = TopLepPt_TOPHLEPBB;
	     vars[9] = TopLepEta_TOPHLEPBB;
	     vars[10] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[11] = TopLepHiggsDr_TOPHLEPBB;
	     vars[12] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[13] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[14] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[15] = TopLepBJetCSVv2_TOPHLEPBB;
	     vars[16] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[17] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[18] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[19] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddBackgroundTrainingEvent(vars,weight);
	     else
	       factory->AddBackgroundTestEvent(vars,weight);
	  }

	factory->PrepareTrainingAndTestTree("","",
					    "SplitMode=Random:NormMode=NumEvents:!V");
	
	factory->BookMethod(TMVA::Types::kBDT,"BDT",
			    "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:IgnoreNegWeightsInTraining");
	
	factory->TrainAllMethods();
	
	factory->TestAllMethods();
	
	factory->EvaluateAllMethods();

	outfile->Write();
	outfile->Close();
     }
   
   // tH Hct vs ttbar
   if( imva == "HctSTallb1j4" || imva == "HctSTeb1j4" || imva == "HctSTmb1j4" )
     {
	std::string outfileName = "TMVA_HctST_all_b1j4";
	if( imva == "HctSTeb1j4" ) outfileName = "TMVA_HctST_e_b1j4";
	else if( imva == "HctSTmb1j4" ) outfileName = "TMVA_HctST_m_b1j4";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;

	double HiggsMass_TOPHLEPBB;
	double TopHadMass_TOPTOPLEPHAD;
	double MVA_TOPHLEPBB;
	double MVA_TOPTOPLEPHAD;
	double HiggsEta_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double TopLepMass_TOPTOPLEPHAD;
	double TopLepPt_TOPHLEPBB;
	double TopLepEta_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double TopLepTopHadDr_TOPTOPLEPHAD;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPTOPLEPHAD;
	double TopHadBJetCSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("TopHadMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');
	
	std::string trName = "tr_all_b1j4";
	if( imva == "HctSTeb1j4" ) trName = "tr_e_b1j4";
	else if( imva == "HctSTmb1j4" ) trName = "tr_m_b1j4";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED_36/ST_TH_1L3B_Hct.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED_36/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ZZ_TuneCUETP8M1_13TeV-pythia8.root");

	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
	
	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(19);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = TopHadMass_TOPTOPLEPHAD;
	     vars[2] = MVA_TOPHLEPBB;
	     vars[3] = MVA_TOPTOPLEPHAD;
	     vars[4] = HiggsEta_TOPHLEPBB;
	     vars[5] = TopLepMass_TOPHLEPBB;
	     vars[6] = TopLepMass_TOPTOPLEPHAD;
	     vars[7] = TopLepPt_TOPHLEPBB;
	     vars[8] = TopLepEta_TOPHLEPBB;
	     vars[9] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[10] = TopLepHiggsDr_TOPHLEPBB;
	     vars[11] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[12] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[13] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[14] = TopLepBJetCSVv2_TOPHLEPBB;
	     vars[15] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[16] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[17] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[18] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddSignalTrainingEvent(vars,weight);
	     else
	       factory->AddSignalTestEvent(vars,weight);
	  }
	
	for(int i=0;i<tr_bkg.GetEntries();i++)
	  {
	     tr_bkg.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = TopHadMass_TOPTOPLEPHAD;
	     vars[2] = MVA_TOPHLEPBB;
	     vars[3] = MVA_TOPTOPLEPHAD;
	     vars[4] = HiggsEta_TOPHLEPBB;
	     vars[5] = TopLepMass_TOPHLEPBB;
	     vars[6] = TopLepMass_TOPTOPLEPHAD;
	     vars[7] = TopLepPt_TOPHLEPBB;
	     vars[8] = TopLepEta_TOPHLEPBB;
	     vars[9] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[10] = TopLepHiggsDr_TOPHLEPBB;
	     vars[11] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[12] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[13] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[14] = TopLepBJetCSVv2_TOPHLEPBB;
	     vars[15] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[16] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[17] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[18] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddBackgroundTrainingEvent(vars,weight);
	     else
	       factory->AddBackgroundTestEvent(vars,weight);
	  }		
	
	factory->PrepareTrainingAndTestTree("","",
					    "SplitMode=Random:NormMode=NumEvents:!V");
	
	factory->BookMethod(TMVA::Types::kBDT,"BDT",
			    "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:IgnoreNegWeightsInTraining");
	
	factory->TrainAllMethods();
	
	factory->TestAllMethods();
	
	factory->EvaluateAllMethods();

	outfile->Write();
	outfile->Close();
     }
   
   // ttbar Hut vs ttbar
   else if( imva == "HutTTallb1j4" || imva == "HutTTeb1j4" || imva == "HutTTmb1j4" )
     {
	std::string outfileName = "TMVA_HutTT_all_b1j4";
	if( imva == "HutTTeb1j4" ) outfileName = "TMVA_HutTT_e_b1j4";
	else if( imva == "HutTTmb1j4" ) outfileName = "TMVA_HutTT_m_b1j4";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;

	double HiggsMass_TOPTOPLEPHBB;
	double MVA_TOPTOPLEPHBB;
	double MVA_TOPTOPLEPHAD;
	double TopLepMass_TOPTOPLEPHBB;
	double TopLepMass_TOPTOPLEPHAD;
	double HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	double TopLepHiggsDr_TOPTOPLEPHBB;
	double TopLepTopHadDr_TOPTOPLEPHAD;
	double HiggsBJet1CSVv2_TOPTOPLEPHBB;
	double HiggsBJet2CSVv2_TOPTOPLEPHBB;
	double TopLepBJetCSVv2_TOPTOPLEPHBB;
	double TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	double TopLepBJetCSVv2_TOPTOPLEPHAD;
	double TopHadBJetCSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPTOPLEPHBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');

	std::string trName = "tr_all_b1j4";
	if( imva == "HutTTeb1j4" ) trName = "tr_e_b1j4";
	else if( imva == "HutTTmb1j4" ) trName = "tr_m_b1j4";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED_36/TT_TopLeptonicDecay_TH_1L3B_Eta_Hut.root");
	tr_sig.Add("../histTEST_MERGED_36/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED_36/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ZZ_TuneCUETP8M1_13TeV-pythia8.root");
   
	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
	
	TRandom3 *r = new TRandom3(666);

	std::vector<double> vars(16);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     
	     vars[0] = HiggsMass_TOPTOPLEPHBB;
	     vars[1] = MVA_TOPTOPLEPHBB;
	     vars[2] = MVA_TOPTOPLEPHAD;
	     vars[3] = TopLepMass_TOPTOPLEPHBB;
	     vars[4] = TopLepMass_TOPTOPLEPHAD;
	     vars[5] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	     vars[6] = TopLepHiggsDr_TOPTOPLEPHBB;
	     vars[7] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[8] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
	     vars[9] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
	     vars[10] = TopLepBJetCSVv2_TOPTOPLEPHBB;
	     vars[11] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	     vars[12] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[13] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[14] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[15] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddSignalTrainingEvent(vars,weight);
	     else
	       factory->AddSignalTestEvent(vars,weight);
	  }
	
	for(int i=0;i<tr_bkg.GetEntries();i++)
	  {
	     tr_bkg.GetEntry(i);

	     if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     
	     vars[0] = HiggsMass_TOPTOPLEPHBB;
	     vars[1] = MVA_TOPTOPLEPHBB;
	     vars[2] = MVA_TOPTOPLEPHAD;
	     vars[3] = TopLepMass_TOPTOPLEPHBB;
	     vars[4] = TopLepMass_TOPTOPLEPHAD;
	     vars[5] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	     vars[6] = TopLepHiggsDr_TOPTOPLEPHBB;
	     vars[7] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[8] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
	     vars[9] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
	     vars[10] = TopLepBJetCSVv2_TOPTOPLEPHBB;
	     vars[11] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	     vars[12] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[13] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[14] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[15] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddBackgroundTrainingEvent(vars,weight);
	     else
	       factory->AddBackgroundTestEvent(vars,weight);
	  }

	factory->PrepareTrainingAndTestTree("","",
					    "SplitMode=Random:NormMode=NumEvents:!V" );
	
	factory->BookMethod(TMVA::Types::kBDT,"BDT",
			    "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:IgnoreNegWeightsInTraining" );
	
	factory->TrainAllMethods();
	
	factory->TestAllMethods();
	
	factory->EvaluateAllMethods();

	outfile->Write();
	outfile->Close();
     }

   // ttbar Hct vs ttbar
   else if( imva == "HctTTallb1j4" || imva == "HctTTeb1j4" || imva == "HctTTmb1j4" )
     {
	std::string outfileName = "TMVA_HctTT_all_b1j4";
	if( imva == "HctTTeb1j4" ) outfileName = "TMVA_HctTT_e_b1j4";
	else if( imva == "HctTTmb1j4" ) outfileName = "TMVA_HctTT_m_b1j4";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;

	double HiggsMass_TOPTOPLEPHBB;
	double MVA_TOPTOPLEPHBB;
	double MVA_TOPTOPLEPHAD;
	double TopLepMass_TOPTOPLEPHBB;
	double TopLepMass_TOPTOPLEPHAD;
	double HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	double TopLepHiggsDr_TOPTOPLEPHBB;
	double TopLepTopHadDr_TOPTOPLEPHAD;
	double HiggsBJet1CSVv2_TOPTOPLEPHBB;
	double HiggsBJet2CSVv2_TOPTOPLEPHBB;
	double TopLepBJetCSVv2_TOPTOPLEPHBB;
	double TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	double TopLepBJetCSVv2_TOPTOPLEPHAD;
	double TopHadBJetCSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	double TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPTOPLEPHBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHBB",'D');
	factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');

	std::string trName = "tr_all_b1j4";
	if( imva == "HctTTeb1j4" ) trName = "tr_e_b1j4";
	else if( imva == "HctTTmb1j4" ) trName = "tr_m_b1j4";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED_36/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
	tr_sig.Add("../histTEST_MERGED_36/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED_36/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED_36/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED_36/ZZ_TuneCUETP8M1_13TeV-pythia8.root");
   
	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_sig.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
	tr_bkg.SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
	
	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(16);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     
	     vars[0] = HiggsMass_TOPTOPLEPHBB;
	     vars[1] = MVA_TOPTOPLEPHBB;
	     vars[2] = MVA_TOPTOPLEPHAD;
	     vars[3] = TopLepMass_TOPTOPLEPHBB;
	     vars[4] = TopLepMass_TOPTOPLEPHAD;
	     vars[5] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	     vars[6] = TopLepHiggsDr_TOPTOPLEPHBB;
	     vars[7] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[8] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
	     vars[9] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
	     vars[10] = TopLepBJetCSVv2_TOPTOPLEPHBB;
	     vars[11] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	     vars[12] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[13] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[14] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[15] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddSignalTrainingEvent(vars,weight);
	     else
	       factory->AddSignalTestEvent(vars,weight);
	  }
	
	for(int i=0;i<tr_bkg.GetEntries();i++)
	  {
	     tr_bkg.GetEntry(i);

	     if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	     if( TopLepMass_TOPTOPLEPHAD > 500. ) TopLepMass_TOPTOPLEPHAD = 500.;
	     
	     vars[0] = HiggsMass_TOPTOPLEPHBB;
	     vars[1] = MVA_TOPTOPLEPHBB;
	     vars[2] = MVA_TOPTOPLEPHAD;
	     vars[3] = TopLepMass_TOPTOPLEPHBB;
	     vars[4] = TopLepMass_TOPTOPLEPHAD;
	     vars[5] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
	     vars[6] = TopLepHiggsDr_TOPTOPLEPHBB;
	     vars[7] = TopLepTopHadDr_TOPTOPLEPHAD;
	     vars[8] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
	     vars[9] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
	     vars[10] = TopLepBJetCSVv2_TOPTOPLEPHBB;
	     vars[11] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
	     vars[12] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	     vars[13] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	     vars[14] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	     vars[15] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
	     
	     float rnd = r->Rndm();
	     if( rnd < 0.5 )
	       factory->AddBackgroundTrainingEvent(vars,weight);
	     else
	       factory->AddBackgroundTestEvent(vars,weight);
	  }

	factory->PrepareTrainingAndTestTree("","",
					    "SplitMode=Random:NormMode=NumEvents:!V" );
	
	factory->BookMethod(TMVA::Types::kBDT,"BDT",
			    "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:IgnoreNegWeightsInTraining" );
	
	factory->TrainAllMethods();
	
	factory->TestAllMethods();
	
	factory->EvaluateAllMethods();

	outfile->Write();
	outfile->Close();
     }
}

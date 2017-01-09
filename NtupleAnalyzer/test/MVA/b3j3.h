#include "runMVA.h"

void b3j3(std::string imva)
{   
   // tH Hut vs ttbar
   if( imva == "HutSTallb3j3" || imva == "HutSTeb3j3" || imva == "HutSTmb3j3" )
     {
	std::string outfileName = "TMVA_HutST_all_b3j3";
	if( imva == "HutSTeb3j3" ) outfileName = "TMVA_HutST_e_b3j3";
	else if( imva == "HutSTmb3j3" ) outfileName = "TMVA_HutST_m_b3j3";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;
	
	double HiggsMass_TOPHLEPBB;
	double MVA_TOPHLEPBB;
	int LepCharge;
	double HiggsEta_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double TopLepPt_TOPHLEPBB;
	double TopLepEta_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("LepCharge",'I');
	factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');

	std::string trName = "tr_all_b3j3";
	if( imva == "HutSTeb3j3" ) trName = "tr_e_b3j3";
	else if( imva == "HutSTmb3j3" ) trName = "tr_m_b3j3";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED/ST_TH_1L3B.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ZZ_TuneCUETP8M1_13TeV-pythia8.root");

	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("LepCharge",&LepCharge);
	tr_sig.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("LepCharge",&LepCharge);
	tr_bkg.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(12);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = LepCharge;
	     vars[3] = HiggsEta_TOPHLEPBB;
	     vars[4] = TopLepMass_TOPHLEPBB;
	     vars[5] = TopLepPt_TOPHLEPBB;
	     vars[6] = TopLepEta_TOPHLEPBB;
	     vars[7] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[8] = TopLepHiggsDr_TOPHLEPBB;
	     vars[9] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[10] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[11] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;

	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = LepCharge;
	     vars[3] = HiggsEta_TOPHLEPBB;
	     vars[4] = TopLepMass_TOPHLEPBB;
	     vars[5] = TopLepPt_TOPHLEPBB;
	     vars[6] = TopLepEta_TOPHLEPBB;
	     vars[7] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[8] = TopLepHiggsDr_TOPHLEPBB;
	     vars[9] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[10] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[11] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
   if( imva == "HctSTallb3j3" || imva == "HctSTeb3j3" || imva == "HctSTmb3j3" )
     {
	std::string outfileName = "TMVA_HctST_all_b3j3";
	if( imva == "HctSTeb3j3" ) outfileName = "TMVA_HctST_e_b3j3";
	else if( imva == "HctSTmb3j3" ) outfileName = "TMVA_HctST_m_b3j3";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;
	
	double HiggsMass_TOPHLEPBB;
	double MVA_TOPHLEPBB;
	double HiggsEta_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double TopLepPt_TOPHLEPBB;
	double TopLepEta_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');

	std::string trName = "tr_all_b3j3";
	if( imva == "HctSTeb3j3" ) trName = "tr_e_b3j3";
	else if( imva == "HctSTmb3j3" ) trName = "tr_m_b3j3";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED/ST_TH_1L3B_Hct.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ZZ_TuneCUETP8M1_13TeV-pythia8.root");

	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(11);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = HiggsEta_TOPHLEPBB;
	     vars[3] = TopLepMass_TOPHLEPBB;
	     vars[4] = TopLepPt_TOPHLEPBB;
	     vars[5] = TopLepEta_TOPHLEPBB;
	     vars[6] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[7] = TopLepHiggsDr_TOPHLEPBB;
	     vars[8] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[9] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[10] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
	     if( TopLepPt_TOPHLEPBB > 1000. ) TopLepPt_TOPHLEPBB = 1000.;

	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = HiggsEta_TOPHLEPBB;
	     vars[3] = TopLepMass_TOPHLEPBB;
	     vars[4] = TopLepPt_TOPHLEPBB;
	     vars[5] = TopLepEta_TOPHLEPBB;
	     vars[6] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[7] = TopLepHiggsDr_TOPHLEPBB;
	     vars[8] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[9] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[10] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
   else if( imva == "HutTTallb3j3" || imva == "HutTTeb3j3" || imva == "HutTTmb3j3" )
     {
	std::string outfileName = "TMVA_HutTT_all_b3j3";
	if( imva == "HutTTeb3j3" ) outfileName = "TMVA_HutTT_e_b3j3";
	else if( imva == "HutTTmb3j3" ) outfileName = "TMVA_HutTT_m_b3j3";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;

	double HiggsMass_TOPHLEPBB;
	double MVA_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');

	std::string trName = "tr_all_b3j3";
	if( imva == "HutTTeb3j3" ) trName = "tr_e_b3j3";
	else if( imva == "HutTTmb3j3" ) trName = "tr_m_b3j3";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_Hut.root");
	tr_sig.Add("../histTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ZZ_TuneCUETP8M1_13TeV-pythia8.root");
   
	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	
	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(8);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = TopLepMass_TOPHLEPBB;
	     vars[3] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[4] = TopLepHiggsDr_TOPHLEPBB;
	     vars[5] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[6] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[7] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = TopLepMass_TOPHLEPBB;
	     vars[3] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[4] = TopLepHiggsDr_TOPHLEPBB;
	     vars[5] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[6] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[7] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
   else if( imva == "HctTTallb3j3" || imva == "HctTTeb3j3" || imva == "HctTTmb3j3" )
     {
	std::string outfileName = "TMVA_HctTT_all_b3j3";
	if( imva == "HctTTeb3j3" ) outfileName = "TMVA_HctTT_e_b3j3";
	else if( imva == "HctTTmb3j3" ) outfileName = "TMVA_HctTT_m_b3j3";
	TFile* outfile = TFile::Open((outfileName+".root").c_str(),"RECREATE");
	
	double weight;

	double HiggsMass_TOPHLEPBB;
	double MVA_TOPHLEPBB;
	double TopLepMass_TOPHLEPBB;
	double HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	double TopLepHiggsDr_TOPHLEPBB;
	double HiggsBJet1CSVv2_TOPHLEPBB;
	double HiggsBJet2CSVv2_TOPHLEPBB;
	double TopLepBJetCSVv2_TOPHLEPBB;

	TMVA::Factory *factory = new TMVA::Factory(outfileName.c_str(),outfile,
						   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
   
	factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	factory->AddVariable("MVA_TOPHLEPBB",'D');
	factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');

	std::string trName = "tr_all_b3j3";
	if( imva == "HctTTeb3j3" ) trName = "tr_e_b3j3";
	else if( imva == "HctTTmb3j3" ) trName = "tr_m_b3j3";
	
	TChain tr_sig(trName.c_str());
	tr_sig.Add("../histTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
	tr_sig.Add("../histTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct.root");

	TChain tr_bkg(trName.c_str());
	tr_bkg.Add("../histTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root");
	tr_bkg.Add("../histTEST_MERGED/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WW_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/WZ_TuneCUETP8M1_13TeV-pythia8.root");
	tr_bkg.Add("../histTEST_MERGED/ZZ_TuneCUETP8M1_13TeV-pythia8.root");
   
	tr_sig.SetBranchAddress("weight",&weight);
	tr_sig.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_sig.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);

	tr_bkg.SetBranchAddress("weight",&weight);
	tr_bkg.SetBranchAddress("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB);
	tr_bkg.SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB);
	
	TRandom3 *r = new TRandom3(666);
	
	std::vector<double> vars(8);

	for(int i=0;i<tr_sig.GetEntries();i++)
	  {
	     tr_sig.GetEntry(i);

	     if( HiggsMass_TOPHLEPBB > 500. ) HiggsMass_TOPHLEPBB = 500.;
	     if( TopLepMass_TOPHLEPBB > 500. ) TopLepMass_TOPHLEPBB = 500.;
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = TopLepMass_TOPHLEPBB;
	     vars[3] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[4] = TopLepHiggsDr_TOPHLEPBB;
	     vars[5] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[6] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[7] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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
	     
	     vars[0] = HiggsMass_TOPHLEPBB;
	     vars[1] = MVA_TOPHLEPBB;
	     vars[2] = TopLepMass_TOPHLEPBB;
	     vars[3] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
	     vars[4] = TopLepHiggsDr_TOPHLEPBB;
	     vars[5] = HiggsBJet1CSVv2_TOPHLEPBB;
	     vars[6] = HiggsBJet2CSVv2_TOPHLEPBB;
	     vars[7] = TopLepBJetCSVv2_TOPHLEPBB;
	     
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

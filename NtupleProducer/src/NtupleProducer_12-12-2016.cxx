#include "../include/NtupleProducer.h"

#include "Riostream.h"
#include "TSystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

//#include <Cintex/Cintex.h>

Tree *ntP;
TChain *ch;
Ntuple *nt;
std::vector<int> *evdebug;
int _isdata;

BTagCalibration *calib;
BTagCalibrationReader *reader_iterativefit;

unsigned int idx;

// ------------ csv applying functions -------------
void fillCSVhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors,
                                  int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

// CSV reweighting
TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

int main(int argc, char *argv[])
{
   if( argc < 4 )
     {
	std::cout << "NtupleProducer usage:" << std::endl;
	std::cout << "--file: input filename" << std::endl;
	std::cout << "--tree: TTree name" << std::endl;
	std::cout << "--noe : Number of events" << std::endl;
	std::cout << "--xsec: Cross-section" << std::endl;
	std::cout << "--nmax: Max number of events" << std::endl;
	std::cout << "--outfile: output file" << std::endl;
	std::cout << "--isdata: is data flag" << std::endl;
	std::cout << "--stream: data stream" << std::endl;
	std::cout << "--issig: is signal" << std::endl;
	exit(1);
     }
   
   const char *fname_str = "output.root";
   const char *stream_str = "FlatTree/tree";
   float noe = 1.;
   float xsec = 1.;
   int nmax = -1;
   const char *outfile_str = "666";
   _isdata = 0;
   int dataStream = -1;
   int issig = 0;
   
   for(int i=0;i<argc;i++)
     {
	if( ! strcmp(argv[i],"--file") ) fname_str = argv[i+1];
	if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
	if( ! strcmp(argv[i],"--noe") ) noe = atof(argv[i+1]);
	if( ! strcmp(argv[i],"--xsec") ) xsec = atof(argv[i+1]);
	if( ! strcmp(argv[i],"--nmax") ) nmax = atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--outfile") ) outfile_str = argv[i+1];
	if( ! strcmp(argv[i],"--isdata") ) _isdata = atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--stream") ) dataStream = atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--issig") ) issig = atoi(argv[i+1]);
     }   

   const char *fname  = fname_str;
   const char *stream = stream_str;
   const char *outfile = outfile_str;
   
   std::cout << "--file=" << fname << std::endl;
   std::cout << "--tree=" << stream << std::endl;
   std::cout << "--noe=" << noe << std::endl;
   std::cout << "--xsec=" << xsec << std::endl;
   std::cout << "--nmax=" << nmax << std::endl;
   std::cout << "--outfile=" << outfile << std::endl;
   std::cout << "--isdata=" << _isdata << std::endl;
   std::cout << "--stream=" << dataStream << std::endl;
   std::cout << "--issig=" << issig << std::endl;

   std::string sfoutput(outfile);
   replace(sfoutput, "output", "sf");
   TFile *sff = new TFile((sfoutput+std::string(".root")).c_str(), "recreate");
   TTree *sft = new TTree("sft", "sft");

   std::cout << fname << std::endl;
   Tree tree(0,const_cast<char*>(fname),stream);
   ntP = &tree;
   
   ch = tree.fChain;
   Long64_t nentries = ch->GetEntries();
   ntP->registerInputBranches(ch);
   
   nt = new Ntuple();

   nt->Init(std::string(outfile));
   nt->createVar();
   nt->setBranchAddress();

   Electron el;
   Muon mu;
   Jet jet;
   Event ev;
   Truth truth;

   calib = new BTagCalibration("csvv2","/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/test/CSVv2_ichep.csv");
   reader_iterativefit = new BTagCalibrationReader(BTagEntry::OP_RESHAPING,"central",
						   {"up_jes","down_jes","up_lf","down_lf",
							"up_hfstats1","down_hfstats1",
							"up_hfstats2","down_hfstats2",
							"up_cferr1","down_cferr1",
							"up_cferr2","down_cferr2"});
   reader_iterativefit->load(*calib,BTagEntry::FLAV_B,"iterativefit");
   reader_iterativefit->load(*calib,BTagEntry::FLAV_C,"iterativefit");
   reader_iterativefit->load(*calib,BTagEntry::FLAV_UDSG,"iterativefit");

   std::cout << "BTagCalibration initialized" << std::endl;
   
   evdebug = new std::vector<int>();
//   evdebug->push_back(91408);
   
   std::cout << "Input events = " << nentries << std::endl;

  //std::string inputFileHF = "data/csv_rwt_fit_hf_2015_11_20.root";
  std::string inputFileHF = "tHFCNC/NtupleProducer/data/csv_rwt_fit_hf_76x_2016_02_08.root";
  //std::string inputFileLF = "data/csv_rwt_fit_lf_2015_11_20.root";
  std::string inputFileLF = "tHFCNC/NtupleProducer/data/csv_rwt_fit_lf_76x_2016_02_08.root";

  TFile* f_CSVwgt_HF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileHF).c_str());
  TFile* f_CSVwgt_LF = new TFile ((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileLF).c_str());

  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

  double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, wgt_csv;
   //
  sft->Branch("wgt_csv_hf", &wgt_csv_hf, "wgt_csv_hf/D", 1);
  sft->Branch("wgt_csv_lf", &wgt_csv_lf, "wgt_csv_lf/D", 1);
  sft->Branch("wgt_csv_cf", &wgt_csv_cf, "wgt_csv_cf/D", 1);
  sft->Branch("wgt_csv", &wgt_csv, "wgt_csv/D", 1);
   
   for(Long64_t i=0;i<nentries;i++)
     {
//	std::cout << i << std::endl;
	if( nmax >= 0 && i >= nmax ) break;

	ch->GetEntry(i);

	nt->clearVar();	
	
	// event
	ev.init();
	ev.read(xsec,noe,dataStream,issig);
	
	nt->NtEvent->push_back(ev);
	
	// muons
	for(int j=0;j<ntP->mu_n;j++)
	  {
	     idx = j;
	     
	     mu.init();
	     mu.read();
	     mu.sel();
	
	     if( mu.isLoose() ) nt->NtMuonLoose->push_back(mu);	     
	     if( mu.isTight() ) nt->NtMuonTight->push_back(mu);	     
	  }		

	// electrons
	for(int j=0;j<ntP->el_n;j++)
	  {
	     idx = j;
	     
	     el.init();
	     el.read();
	     el.sel();

	     if( el.isLoose() ) nt->NtElectronLoose->push_back(el);
	     if( el.isTight() ) nt->NtElectronTight->push_back(el);
	  }


    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int> jetFlavors;
    int iSys = 0; //nominal case

    //std::cout << "CIRKOVIC JETS:" << std::endl;
	// jets
	for(int j=0;j<ntP->jet_n;j++)
	  {
	     idx = j;
	     
	     jet.init();
	     jet.read();
	     jet.sel();

         if( jet.isLoose() ) {
             //std::cout << "\t(pt, eta, csv, flavour) = (" << jet.pt() << ", " << jet.eta() << ", " << jet.CSVv2() << ", " << jet.hadronFlavour() << ")" << std::endl;

             jetPts.push_back(jet.pt());
             jetEtas.push_back(jet.eta());
             jetCSVs.push_back(jet.CSVv2());
             jetFlavors.push_back(jet.hadronFlavour());
         }

	     if( jet.isLoose() ) nt->NtJetLoose->push_back(jet);
	     if( jet.isTight() ) nt->NtJetTight->push_back(jet);
	     if( jet.isTight() && jet.isBTag() ) nt->NtBJetTight->push_back(jet);
	  }
    //std::cout << std::endl;

// Loop over jets and fill vectors: jetPts, jetEtas, jetCSVs, and jetFlavors with std::vector's of those quantities

   //double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
   ////double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
   //double wgt_csv = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
   wgt_csv_hf = 1.0;
   wgt_csv_lf = 1.0;
   wgt_csv_cf = 1.0;
   wgt_csv = 1.0;
   // std::cout << "CIRKOVIC: (wgt_csv, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf) = (" << wgt_csv << ", " << wgt_csv_hf << ", " << wgt_csv_lf << ", " << wgt_csv_cf << ")" << std::endl;
   wgt_csv = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

    //double wgt_lumi_csv = wgt_lumi * wgt_csv;
   // std::cout << "CIRKOVIC: (wgt_csv, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf) = (" << wgt_csv << ", " << wgt_csv_hf << ", " << wgt_csv_lf << ", " << wgt_csv_cf << ")" << std::endl;
   // std::cout << std::endl;
    //sft->Fill();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if( !_isdata )
	  {	     
	     // truth
	     truth.init();
	     truth.read();
	     
	     nt->NtTruth->push_back(truth);
	  }	

	int nElecLoose = nt->NtElectronLoose->size();
	int nMuonLoose = nt->NtMuonLoose->size();
	
	int nElecTight = nt->NtElectronTight->size();
	int nMuonTight = nt->NtMuonTight->size();

	int nJetLoose = nt->NtJetLoose->size();
	int nJetTight = nt->NtJetTight->size();
	int nJetBTag = nt->NtBJetTight->size();

	if(
	   nElecTight+nMuonTight == 1 &&
	   nElecLoose+nMuonLoose == 1 &&
	   nJetTight >= 3 
	  )
	  {
	     nt->fill();
         sft->Fill();
	  }
     }
   
   std::cout << "Output events = " << nt->m_tree->GetEntries() << std::endl;

   nt->write();
   nt->close();
   //delete nt;

   sff->cd();
   sft->Write();
   sff->Close();
}


// fill the histograms (done once)
void fillCSVhistos(TFile* fileHF, TFile* fileLF){
  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }
  return;
}


double get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, 
                                  int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){
  int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPts.size()); iJet++ ){

    double csv = jetCSVs[iJet];
    double jetPt = jetPts[iJet];
    double jetAbsEta = fabs(jetEtas[iJet]);
    int flavor = jetFlavors[iJet];

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100) iPt = 4;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}

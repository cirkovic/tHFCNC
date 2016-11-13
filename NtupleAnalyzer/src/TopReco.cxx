#include "../include/TopReco.h"

#include "TRandom3.h"

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/type_traits.hpp>

struct sortFunc
{
   bool operator()(const std::pair<float,int> &left, const std::pair<float,int> &right)
     {
	return left.first > right.first;
     }
};

TopReco::TopReco()
{
   help = new Helper();
   
   _v_ElectronTight = new std::vector<Electron>;
   _v_MuonTight = new std::vector<Muon>;
   _v_JetTight = new std::vector<Jet>;
   _v_BJetTight = new std::vector<Jet>;
   _v_NonBJetTight = new std::vector<Jet>;
   
   kfTopTopLepHbb = new KINFIT::kfit();
   kfTopHLepbb = new KINFIT::kfit();
   kfTopTopLepHad = new KINFIT::kfit();
   
   MVAFullRecoReaderTOPTOPLEPHBB = new TMVA::Reader("!Color:!Silent");
   MVAFullRecoReaderTOPTOPLEPHBB->AddVariable("HiggsRecMCSVv2M",&MVAFullReco_HiggsRecM_TOPTOPLEPHBB);
   MVAFullRecoReaderTOPTOPLEPHBB->AddVariable("TopLepRecMCSVv2M",&MVAFullReco_TopLepRecM_TOPTOPLEPHBB);
   MVAFullRecoReaderTOPTOPLEPHBB->AddVariable("HiggsTopLepRecDrCSVv2M",&MVAFullReco_HiggsTopLepRecDr_TOPTOPLEPHBB);
   MVAFullRecoReaderTOPTOPLEPHBB->AddVariable("TopLepRecPtCSVv2M",&MVAFullReco_TopLepRecPt_TOPTOPLEPHBB);
   std::string weightsFileFullRecoTOPTOPLEPHBB = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoCSVv2M_BDT.weights.xml";
   MVAFullRecoReaderTOPTOPLEPHBB->BookMVA("BDTG method",weightsFileFullRecoTOPTOPLEPHBB.c_str());
   
   MVAPartRecoReaderTOPTOPLEPHBB = new TMVA::Reader("!Color:!Silent");
   MVAPartRecoReaderTOPTOPLEPHBB->AddVariable("HiggsRecMCSVv2M",&MVAPartReco_HiggsRecM_TOPTOPLEPHBB);
   MVAPartRecoReaderTOPTOPLEPHBB->AddVariable("TopLepRecMTCSVv2M",&MVAPartReco_TopLepRecMT_TOPTOPLEPHBB);
   MVAPartRecoReaderTOPTOPLEPHBB->AddVariable("HiggsTopLepRecDphiTCSVv2M",&MVAPartReco_HiggsTopLepRecDphiT_TOPTOPLEPHBB);
   MVAPartRecoReaderTOPTOPLEPHBB->AddVariable("TopLepRecPtTCSVv2M",&MVAPartReco_TopLepRecPtT_TOPTOPLEPHBB);
   std::string weightsFilePartRecoTOPTOPLEPHBB = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoCSVv2M_BDT.weights.xml";
   MVAPartRecoReaderTOPTOPLEPHBB->BookMVA("BDTG method",weightsFilePartRecoTOPTOPLEPHBB.c_str());

   MVAFullRecoReaderTOPHLEPBB = new TMVA::Reader("!Color:!Silent");
   MVAFullRecoReaderTOPHLEPBB->AddVariable("HiggsRecMCSVv2M",&MVAFullReco_HiggsRecM_TOPHLEPBB);
   MVAFullRecoReaderTOPHLEPBB->AddVariable("TopLepRecMCSVv2M",&MVAFullReco_TopLepRecM_TOPHLEPBB);
   MVAFullRecoReaderTOPHLEPBB->AddVariable("HiggsTopLepRecDrCSVv2M",&MVAFullReco_HiggsTopLepRecDr_TOPHLEPBB);
   MVAFullRecoReaderTOPHLEPBB->AddVariable("TopLepRecPtCSVv2M",&MVAFullReco_TopLepRecPt_TOPHLEPBB);
   std::string weightsFileFullRecoTOPHLEPBB = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAFullRecoHutCSVv2M_BDT.weights.xml";
   MVAFullRecoReaderTOPHLEPBB->BookMVA("BDTG method",weightsFileFullRecoTOPHLEPBB.c_str());
   
   MVAPartRecoReaderTOPHLEPBB = new TMVA::Reader("!Color:!Silent");
   MVAPartRecoReaderTOPHLEPBB->AddVariable("HiggsRecMCSVv2M",&MVAPartReco_HiggsRecM_TOPHLEPBB);
   MVAPartRecoReaderTOPHLEPBB->AddVariable("TopLepRecMTCSVv2M",&MVAPartReco_TopLepRecMT_TOPHLEPBB);
   MVAPartRecoReaderTOPHLEPBB->AddVariable("HiggsTopLepRecDphiTCSVv2M",&MVAPartReco_HiggsTopLepRecDphiT_TOPHLEPBB);
   MVAPartRecoReaderTOPHLEPBB->AddVariable("TopLepRecPtTCSVv2M",&MVAPartReco_TopLepRecPtT_TOPHLEPBB);
   std::string weightsFilePartRecoTOPHLEPBB = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAPartRecoHutCSVv2M_BDT.weights.xml";
   MVAPartRecoReaderTOPHLEPBB->BookMVA("BDTG method",weightsFilePartRecoTOPHLEPBB.c_str());

   MVAFullRecoReaderTOPTOPLEPHAD = new TMVA::Reader("!Color:!Silent");
   MVAFullRecoReaderTOPTOPLEPHAD->AddVariable("TopHadRecMCSVv2M",&MVAFullReco_TopHadRecM_TOPTOPLEPHAD);
   MVAFullRecoReaderTOPTOPLEPHAD->AddVariable("TopLepRecMCSVv2M",&MVAFullReco_TopLepRecM_TOPTOPLEPHAD);
   MVAFullRecoReaderTOPTOPLEPHAD->AddVariable("TopLepTopHadRecDrCSVv2M",&MVAFullReco_TopLepTopHadRecDr_TOPTOPLEPHAD);
   MVAFullRecoReaderTOPTOPLEPHAD->AddVariable("TopLepRecPtCSVv2M",&MVAFullReco_TopLepRecPt_TOPTOPLEPHAD);
   std::string weightsFileFullRecoTOPTOPLEPHAD = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoCSVv2M_BDT.weights.xml";
   MVAFullRecoReaderTOPTOPLEPHAD->BookMVA("BDTG method",weightsFileFullRecoTOPTOPLEPHAD.c_str());
   
   MVAPartRecoReaderTOPTOPLEPHAD = new TMVA::Reader("!Color:!Silent");
   MVAPartRecoReaderTOPTOPLEPHAD->AddVariable("TopHadRecMCSVv2M",&MVAPartReco_TopHadRecM_TOPTOPLEPHAD);
   MVAPartRecoReaderTOPTOPLEPHAD->AddVariable("TopLepRecMTCSVv2M",&MVAPartReco_TopLepRecMT_TOPTOPLEPHAD);
   MVAPartRecoReaderTOPTOPLEPHAD->AddVariable("TopLepTopHadRecDphiTCSVv2M",&MVAPartReco_TopLepTopHadRecDphiT_TOPTOPLEPHAD);
   MVAPartRecoReaderTOPTOPLEPHAD->AddVariable("TopLepRecPtTCSVv2M",&MVAPartReco_TopLepRecPtT_TOPTOPLEPHAD);
   std::string weightsFilePartRecoTOPTOPLEPHAD = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoCSVv2M_BDT.weights.xml";
   MVAPartRecoReaderTOPTOPLEPHAD->BookMVA("BDTG method",weightsFilePartRecoTOPTOPLEPHAD.c_str());
}

TopReco::~TopReco()
{
   delete help;
   
   delete _v_ElectronTight;
   delete _v_MuonTight;
   delete _v_JetTight;
   delete _v_BJetTight;
   delete _v_NonBJetTight;
   
   delete kfTopTopLepHbb;
   delete kfTopHLepbb;
   delete kfTopTopLepHad;
   
   delete MVAFullRecoReaderTOPTOPLEPHBB;
   delete MVAPartRecoReaderTOPTOPLEPHBB;

   delete MVAFullRecoReaderTOPHLEPBB;
   delete MVAPartRecoReaderTOPHLEPBB;

   delete MVAFullRecoReaderTOPTOPLEPHAD;
   delete MVAPartRecoReaderTOPTOPLEPHAD;
}

void TopReco::init()
{
   rnd = new TRandom3(666);

   kfTopTopLepHbb->Init(TOPTOPLEPHBB);
   
   std::string pdfFileNameTopTopLepHbb = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
   
   kfTopTopLepHbb->SetPDF("TopWMass",pdfFileNameTopTopLepHbb.c_str(),"TopLepWM_Fit");
   kfTopTopLepHbb->SetPDF("TopMass",pdfFileNameTopTopLepHbb.c_str(),"TopLepRecM_Fit");
   kfTopTopLepHbb->SetPDF("HiggsMass",pdfFileNameTopTopLepHbb.c_str(),"HiggsRecM_Fit");
   kfTopTopLepHbb->SetPDF("TopHadMass",pdfFileNameTopTopLepHbb.c_str(),"TopHadRecM_Fit");
   kfTopTopLepHbb->SetPDF("MetPx",pdfFileNameTopTopLepHbb.c_str(),"dMetPx_Gaus");
   kfTopTopLepHbb->SetPDF("MetPy",pdfFileNameTopTopLepHbb.c_str(),"dMetPy_Gaus");
   kfTopTopLepHbb->SetPDF("BJetPx",pdfFileNameTopTopLepHbb.c_str(),"dBJetPx_Fit");
   kfTopTopLepHbb->SetPDF("BJetPy",pdfFileNameTopTopLepHbb.c_str(),"dBJetPy_Fit");
   kfTopTopLepHbb->SetPDF("BJetPz",pdfFileNameTopTopLepHbb.c_str(),"dBJetPz_Fit");
   kfTopTopLepHbb->SetPDF("BJetE",pdfFileNameTopTopLepHbb.c_str(),"dBJetE_Fit");
   kfTopTopLepHbb->SetPDF("NonBJetPx",pdfFileNameTopTopLepHbb.c_str(),"dNonBJetPx_Fit");
   kfTopTopLepHbb->SetPDF("NonBJetPy",pdfFileNameTopTopLepHbb.c_str(),"dNonBJetPy_Fit");
   kfTopTopLepHbb->SetPDF("NonBJetPz",pdfFileNameTopTopLepHbb.c_str(),"dNonBJetPz_Fit");
   kfTopTopLepHbb->SetPDF("NonBJetE",pdfFileNameTopTopLepHbb.c_str(),"dNonBJetE_Fit");
   kfTopTopLepHbb->SetPDF("ElecPx",pdfFileNameTopTopLepHbb.c_str(),"dElecPx_Fit");
   kfTopTopLepHbb->SetPDF("ElecPy",pdfFileNameTopTopLepHbb.c_str(),"dElecPy_Fit");
   kfTopTopLepHbb->SetPDF("ElecPz",pdfFileNameTopTopLepHbb.c_str(),"dElecPz_Fit");
   kfTopTopLepHbb->SetPDF("ElecE",pdfFileNameTopTopLepHbb.c_str(),"dElecE_Fit");
   kfTopTopLepHbb->SetPDF("MuonPx",pdfFileNameTopTopLepHbb.c_str(),"dMuonPx_Fit");
   kfTopTopLepHbb->SetPDF("MuonPy",pdfFileNameTopTopLepHbb.c_str(),"dMuonPy_Fit");
   kfTopTopLepHbb->SetPDF("MuonPz",pdfFileNameTopTopLepHbb.c_str(),"dMuonPz_Fit");
   kfTopTopLepHbb->SetPDF("MuonE",pdfFileNameTopTopLepHbb.c_str(),"dMuonE_Fit");

   kfTopTopLepHbb->SetNToy(20);

   kfTopHLepbb->Init(TOPHLEPBB);
   
   std::string pdfFileNameTopHLepbb = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";
   
   kfTopHLepbb->SetPDF("TopWMass",pdfFileNameTopHLepbb.c_str(),"TopLepWM_Fit");
   kfTopHLepbb->SetPDF("TopMass",pdfFileNameTopHLepbb.c_str(),"TopLepRecM_Fit");
   kfTopHLepbb->SetPDF("HiggsMass",pdfFileNameTopHLepbb.c_str(),"HiggsRecM_Fit");
   kfTopHLepbb->SetPDF("MetPx",pdfFileNameTopHLepbb.c_str(),"dMetPx_Gaus");
   kfTopHLepbb->SetPDF("MetPy",pdfFileNameTopHLepbb.c_str(),"dMetPy_Gaus");
   kfTopHLepbb->SetPDF("BJetPx",pdfFileNameTopHLepbb.c_str(),"dBJetPx_Fit");
   kfTopHLepbb->SetPDF("BJetPy",pdfFileNameTopHLepbb.c_str(),"dBJetPy_Fit");
   kfTopHLepbb->SetPDF("BJetPz",pdfFileNameTopHLepbb.c_str(),"dBJetPz_Fit");
   kfTopHLepbb->SetPDF("BJetE",pdfFileNameTopHLepbb.c_str(),"dBJetE_Fit");
   kfTopHLepbb->SetPDF("ElecPx",pdfFileNameTopHLepbb.c_str(),"dElecPx_Fit");
   kfTopHLepbb->SetPDF("ElecPy",pdfFileNameTopHLepbb.c_str(),"dElecPy_Fit");
   kfTopHLepbb->SetPDF("ElecPz",pdfFileNameTopHLepbb.c_str(),"dElecPz_Fit");
   kfTopHLepbb->SetPDF("ElecE",pdfFileNameTopHLepbb.c_str(),"dElecE_Fit");
   kfTopHLepbb->SetPDF("MuonPx",pdfFileNameTopHLepbb.c_str(),"dMuonPx_Fit");
   kfTopHLepbb->SetPDF("MuonPy",pdfFileNameTopHLepbb.c_str(),"dMuonPy_Fit");
   kfTopHLepbb->SetPDF("MuonPz",pdfFileNameTopHLepbb.c_str(),"dMuonPz_Fit");
   kfTopHLepbb->SetPDF("MuonE",pdfFileNameTopHLepbb.c_str(),"dMuonE_Fit");

   kfTopHLepbb->SetNToy(20);

   kfTopTopLepHad->Init(TOPTOPLEPHAD);
   
   std::string pdfFileNameTopTopLepHad = "/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
   
   kfTopTopLepHad->SetPDF("TopWMass",pdfFileNameTopTopLepHad.c_str(),"TopLepWM_Fit");
   kfTopTopLepHad->SetPDF("TopMass",pdfFileNameTopTopLepHad.c_str(),"TopLepRecM_Fit");
   kfTopTopLepHad->SetPDF("TopWHadMass",pdfFileNameTopTopLepHad.c_str(),"TopHadWRecM_Fit");
   kfTopTopLepHad->SetPDF("TopHadMass",pdfFileNameTopTopLepHad.c_str(),"TopHadRecM_Fit");
   kfTopTopLepHad->SetPDF("MetPx",pdfFileNameTopTopLepHad.c_str(),"dMetPx_Gaus");
   kfTopTopLepHad->SetPDF("MetPy",pdfFileNameTopTopLepHad.c_str(),"dMetPy_Gaus");
   kfTopTopLepHad->SetPDF("BJetPx",pdfFileNameTopTopLepHad.c_str(),"dBJetPx_Fit");
   kfTopTopLepHad->SetPDF("BJetPy",pdfFileNameTopTopLepHad.c_str(),"dBJetPy_Fit");
   kfTopTopLepHad->SetPDF("BJetPz",pdfFileNameTopTopLepHad.c_str(),"dBJetPz_Fit");
   kfTopTopLepHad->SetPDF("BJetE",pdfFileNameTopTopLepHad.c_str(),"dBJetE_Fit");
   kfTopTopLepHad->SetPDF("NonBJetPx",pdfFileNameTopTopLepHad.c_str(),"dNonBJetPx_Fit");
   kfTopTopLepHad->SetPDF("NonBJetPy",pdfFileNameTopTopLepHad.c_str(),"dNonBJetPy_Fit");
   kfTopTopLepHad->SetPDF("NonBJetPz",pdfFileNameTopTopLepHad.c_str(),"dNonBJetPz_Fit");
   kfTopTopLepHad->SetPDF("NonBJetE",pdfFileNameTopTopLepHad.c_str(),"dNonBJetE_Fit");
   kfTopTopLepHad->SetPDF("ElecPx",pdfFileNameTopTopLepHad.c_str(),"dElecPx_Fit");
   kfTopTopLepHad->SetPDF("ElecPy",pdfFileNameTopTopLepHad.c_str(),"dElecPy_Fit");
   kfTopTopLepHad->SetPDF("ElecPz",pdfFileNameTopTopLepHad.c_str(),"dElecPz_Fit");
   kfTopTopLepHad->SetPDF("ElecE",pdfFileNameTopTopLepHad.c_str(),"dElecE_Fit");
   kfTopTopLepHad->SetPDF("MuonPx",pdfFileNameTopTopLepHad.c_str(),"dMuonPx_Fit");
   kfTopTopLepHad->SetPDF("MuonPy",pdfFileNameTopTopLepHad.c_str(),"dMuonPy_Fit");
   kfTopTopLepHad->SetPDF("MuonPz",pdfFileNameTopTopLepHad.c_str(),"dMuonPz_Fit");
   kfTopTopLepHad->SetPDF("MuonE",pdfFileNameTopTopLepHad.c_str(),"dMuonE_Fit");

   kfTopTopLepHad->SetNToy(20);
   
   std::cout << "TopReco initialisation done" << std::endl;
}

void TopReco::run()
{
   float w = _v_Event->at(0).mc_weight();

   bool isData = _v_Event->at(0).isData();
   
   float xsec = _v_Event->at(0).xsec();
   float noe = _v_Event->at(0).noe();
   float wil = 19200.*xsec/noe;
   if( !isData ) w *= wil;

   _v_ElectronTight->clear();
   _v_MuonTight->clear();
   _v_JetTight->clear();
   _v_BJetTight->clear();
   _v_NonBJetTight->clear();

   for(int i=0;i<_v_Electron->size();i++)
     {
	if( _v_Electron->at(i).isTight() ) _v_ElectronTight->push_back(_v_Electron->at(i));
     }       
   
   for(int i=0;i<_v_Muon->size();i++)
     {
	if( _v_Muon->at(i).isTight() ) _v_MuonTight->push_back(_v_Muon->at(i));
     }      
   
   for(int i=0;i<_v_Jet->size();i++)
     {
	if( _v_Jet->at(i).isTight() )
	  {	     
	     _v_JetTight->push_back(_v_Jet->at(i));
	     if( _v_Jet->at(i).isBTag() ) _v_BJetTight->push_back(_v_Jet->at(i));
	     else _v_NonBJetTight->push_back(_v_Jet->at(i));
	  }	
     }

   int id = _v_Event->at(0).id();
   int run = _v_Event->at(0).run();
   int lumi = _v_Event->at(0).lumi();

   float metpt = _v_Event->at(0).metpt();
   float metphi = _v_Event->at(0).metphi();
   float metpx = metpt*cos(metphi);
   float metpy = metpt*sin(metphi);
   float metE = sqrt(metpx*metpx+metpy*metpy);
   
   int njets = _v_JetTight->size();
   int nbjets = _v_BJetTight->size();
   int nnonbjets = _v_NonBJetTight->size();
   
   int nLep = _v_ElectronTight->size()+_v_MuonTight->size();

   _chi2_TOPTOPLEPHBB = 10E+10;
   _chi2_TOPHLEPBB = 10E+10;
   _chi2_TOPTOPLEPHAD = 10E+10;

   _MVA_TOPTOPLEPHBB = 10E+10;
   _MVA_TOPHLEPBB = 10E+10;
   _MVA_TOPTOPLEPHAD = 10E+10;
   
   _TopLepBJet_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _HiggsBJet1_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _HiggsBJet2_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWLep_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWNu_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLep_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _Higgs_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepW_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHadNonBJet_TOPTOPLEPHBB_p4.SetPtEtaPhiE(0,0,0,0);

   _TopLepBJet_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _HiggsBJet1_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _HiggsBJet2_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWLep_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWNu_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLep_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _Higgs_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepW_TOPHLEPBB_p4.SetPtEtaPhiE(0,0,0,0);

   _TopLepBJet_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHadWNonBJet1_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHadWNonBJet2_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWLep_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepWNu_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLep_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHad_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHadW_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopLepW_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   _TopHadBJet_TOPTOPLEPHAD_p4.SetPtEtaPhiE(0,0,0,0);
   
   // tH ttbar
   if( nLep == 1 &&
       ((nbjets >= 3 && nnonbjets >= 1) || (nbjets == 2 && nnonbjets >= 2)) )
     {		
	kfTopTopLepHbb->Init(TOPTOPLEPHBB);
	
	std::vector<float> BJetPt;
	std::vector<float> BJetEta;
	std::vector<float> BJetPhi;
	std::vector<float> BJetE;
	
	std::vector<float> NonBJetPt;
	std::vector<float> NonBJetEta;
	std::vector<float> NonBJetPhi;
	std::vector<float> NonBJetE;
	
	std::vector<float> ElectronPt;
	std::vector<float> ElectronEta;
	std::vector<float> ElectronPhi;
	std::vector<float> ElectronE;
	
	std::vector<float> MuonPt;
	std::vector<float> MuonEta;
	std::vector<float> MuonPhi;
	std::vector<float> MuonE;
	
	for(int i=0;i<_v_ElectronTight->size();i++)
	  {
	     ElectronPt.push_back(_v_ElectronTight->at(i).pt());
	     ElectronEta.push_back(_v_ElectronTight->at(i).eta());
	     ElectronPhi.push_back(_v_ElectronTight->at(i).phi());
	     ElectronE.push_back(_v_ElectronTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_MuonTight->size();i++)
	  {
	     MuonPt.push_back(_v_MuonTight->at(i).pt());
	     MuonEta.push_back(_v_MuonTight->at(i).eta());
	     MuonPhi.push_back(_v_MuonTight->at(i).phi());
	     MuonE.push_back(_v_MuonTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_BJetTight->size();i++)
	  {
	     BJetPt.push_back(_v_BJetTight->at(i).pt());
	     BJetEta.push_back(_v_BJetTight->at(i).eta());
	     BJetPhi.push_back(_v_BJetTight->at(i).phi());
	     BJetE.push_back(_v_BJetTight->at(i).E());
	  }   
	     
	for(int i=0;i<_v_NonBJetTight->size();i++)
	  {
	     if( nbjets == 2 && i == 0 )
	       {
		  BJetPt.push_back(_v_NonBJetTight->at(i).pt());
		  BJetEta.push_back(_v_NonBJetTight->at(i).eta());
		  BJetPhi.push_back(_v_NonBJetTight->at(i).phi());
		  BJetE.push_back(_v_NonBJetTight->at(i).E());
	       }
	     else
	       {
		  NonBJetPt.push_back(_v_NonBJetTight->at(i).pt());
		  NonBJetEta.push_back(_v_NonBJetTight->at(i).eta());
		  NonBJetPhi.push_back(_v_NonBJetTight->at(i).phi());
		  NonBJetE.push_back(_v_NonBJetTight->at(i).E());
	       }	     
	  }   
	
	kfTopTopLepHbb->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
	kfTopTopLepHbb->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
	kfTopTopLepHbb->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
	kfTopTopLepHbb->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
	kfTopTopLepHbb->SetMet(metpx,metpy);
	
	kfTopTopLepHbb->Run();

	std::vector<std::pair<float,int> > MVADisc;
	
	int nPerm = kfTopTopLepHbb->GetNPerm();
	
	for(int ip=0;ip<nPerm;ip++)
	  {
	     float disc = kfTopTopLepHbb->GetDisc(ip);

	     int idx_ELECTRON_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(ELECTRON_TOPTOPLEPHBB,ip);
	     int idx_MUON_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(MUON_TOPTOPLEPHBB,ip);
	     int idx_BJETLEP_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJETLEP_TOPTOPLEPHBB,ip);
	     int idx_NONBJETHAD_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip);
	     int idx_BJET1_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJET1_TOPTOPLEPHBB,ip);
	     int idx_BJET2_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJET2_TOPTOPLEPHBB,ip);

	     float NuPx = kfTopTopLepHbb->GetNuPx(ip,0);
	     float NuPy = kfTopTopLepHbb->GetNuPy(ip,0);
	     float NuPz = kfTopTopLepHbb->GetNuPz(ip,0);
	     float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	     TLorentzVector *TopLepWNu_p4 = new TLorentzVector();
	     TopLepWNu_p4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
	     
	     TLorentzVector *HiggsBJet1_p4 = new TLorentzVector();
	     HiggsBJet1_p4->SetPtEtaPhiE(BJetPt[idx_BJET1_TOPTOPLEPHBB],
					 BJetEta[idx_BJET1_TOPTOPLEPHBB],
					 BJetPhi[idx_BJET1_TOPTOPLEPHBB],
					 BJetE[idx_BJET1_TOPTOPLEPHBB]);
	     
	     TLorentzVector *HiggsBJet2_p4 = new TLorentzVector();
	     HiggsBJet2_p4->SetPtEtaPhiE(BJetPt[idx_BJET2_TOPTOPLEPHBB],
					BJetEta[idx_BJET2_TOPTOPLEPHBB],
					BJetPhi[idx_BJET2_TOPTOPLEPHBB],
					BJetE[idx_BJET2_TOPTOPLEPHBB]);
	     
	     TLorentzVector *TopLepBJet_p4 = new TLorentzVector();
	     TopLepBJet_p4->SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPTOPLEPHBB],
					 BJetEta[idx_BJETLEP_TOPTOPLEPHBB],
					 BJetPhi[idx_BJETLEP_TOPTOPLEPHBB],
					 BJetE[idx_BJETLEP_TOPTOPLEPHBB]);
	     
	     TLorentzVector *TopHadNonBJet_p4 = new TLorentzVector();
	     TopHadNonBJet_p4->SetPtEtaPhiE(NonBJetPt[idx_NONBJETHAD_TOPTOPLEPHBB],
					    NonBJetEta[idx_NONBJETHAD_TOPTOPLEPHBB],
					    NonBJetPhi[idx_NONBJETHAD_TOPTOPLEPHBB],
					    NonBJetE[idx_NONBJETHAD_TOPTOPLEPHBB]);
	     
	     TLorentzVector *TopLepWLep_p4 = new TLorentzVector();
	     if( idx_ELECTRON_TOPTOPLEPHBB >= 0 )
	       TopLepWLep_p4->SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPTOPLEPHBB],
					   ElectronEta[idx_ELECTRON_TOPTOPLEPHBB],
					   ElectronPhi[idx_ELECTRON_TOPTOPLEPHBB],
					   ElectronE[idx_ELECTRON_TOPTOPLEPHBB]);
	     else
	       TopLepWLep_p4->SetPtEtaPhiE(MuonPt[idx_MUON_TOPTOPLEPHBB],
					   MuonEta[idx_MUON_TOPTOPLEPHBB],
					   MuonPhi[idx_MUON_TOPTOPLEPHBB],
					   MuonE[idx_MUON_TOPTOPLEPHBB]);
	     
	     TLorentzVector Higgs_p4 = *HiggsBJet1_p4+*HiggsBJet2_p4;
	     TLorentzVector TopLepW_p4 = *TopLepWLep_p4+*TopLepWNu_p4;
	     TLorentzVector TopLep_p4 = TopLepW_p4+*TopLepBJet_p4;

	     float VarHiggsRecM = Higgs_p4.M();
	     float VarTopLepRecM = TopLep_p4.M();
	     float VarHiggsTopLepRecDr = Higgs_p4.DeltaR(TopLep_p4);
	     float VarTopLepRecPt = TopLep_p4.Pt();
	     
	     TLorentzVector *HiggsT_p4 = new TLorentzVector();
	     HiggsT_p4->SetPxPyPzE(Higgs_p4.Px(),Higgs_p4.Py(),0.,Higgs_p4.Et());
	     
	     TLorentzVector *TopLepWLepT_p4 = new TLorentzVector();
	     TopLepWLepT_p4->SetPxPyPzE(TopLepWLep_p4->Px(),TopLepWLep_p4->Py(),0.,TopLepWLep_p4->Et());
	     
	     TLorentzVector *TopLepWNuT_p4 = new TLorentzVector();
	     TopLepWNuT_p4->SetPxPyPzE(TopLepWNu_p4->Px(),TopLepWNu_p4->Py(),0.,TopLepWNu_p4->Et());
	     
	     TLorentzVector *TopLepBJetT_p4 = new TLorentzVector();
	     TopLepBJetT_p4->SetPxPyPzE(TopLepBJet_p4->Px(),TopLepBJet_p4->Py(),0.,TopLepBJet_p4->Et());
	     
	     TLorentzVector TopLepT_p4 = *TopLepWLepT_p4+*TopLepWNuT_p4+*TopLepBJetT_p4;
	     
	     float VarTopLepRecMT = TopLepT_p4.M();
	     float VarHiggsTopLepRecDphiT = HiggsT_p4->DeltaPhi(TopLepT_p4);
	     float VarTopLepRecPtT = TopLepT_p4.Pt();
	     
	     if( disc < 10E+8 )
	       {
		  MVAFullReco_HiggsRecM_TOPTOPLEPHBB = VarHiggsRecM;
		  MVAFullReco_TopLepRecM_TOPTOPLEPHBB = VarTopLepRecM;
		  MVAFullReco_HiggsTopLepRecDr_TOPTOPLEPHBB = VarHiggsTopLepRecDr;
		  MVAFullReco_TopLepRecPt_TOPTOPLEPHBB = VarTopLepRecPt;
		  
		  MVADisc.push_back(std::make_pair(MVAFullRecoReaderTOPTOPLEPHBB->EvaluateMVA("BDTG method"),ip));
	       }	
	     else
	       {
		  MVAPartReco_HiggsRecM_TOPTOPLEPHBB = VarHiggsRecM;
		  MVAPartReco_TopLepRecMT_TOPTOPLEPHBB = VarTopLepRecMT;
		  MVAPartReco_HiggsTopLepRecDphiT_TOPTOPLEPHBB = VarHiggsTopLepRecDphiT;
		  MVAPartReco_TopLepRecPtT_TOPTOPLEPHBB = VarTopLepRecPtT;
		  
		  MVADisc.push_back(std::make_pair(MVAPartRecoReaderTOPTOPLEPHBB->EvaluateMVA("BDTG method"),ip));
	       }	     

	     delete TopLepWNu_p4;
	     delete HiggsBJet1_p4;
	     delete HiggsBJet2_p4;
	     delete TopLepBJet_p4;
	     delete TopHadNonBJet_p4;
	     delete TopLepWLep_p4;
	     delete HiggsT_p4;
	     delete TopLepWLepT_p4;
	     delete TopLepWNuT_p4;
	     delete TopLepBJetT_p4;
	  }
	
	std::sort(MVADisc.begin(),MVADisc.end(),sortFunc());
	int IdxPerm = MVADisc.at(0).second;
	
	_chi2_TOPTOPLEPHBB = kfTopTopLepHbb->GetDisc(IdxPerm);
	_MVA_TOPTOPLEPHBB = MVADisc.at(0).first;

	int idx_ELECTRON_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(ELECTRON_TOPTOPLEPHBB,IdxPerm);
	int idx_MUON_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(MUON_TOPTOPLEPHBB,IdxPerm);
	int idx_BJETLEP_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJETLEP_TOPTOPLEPHBB,IdxPerm);
	int idx_NONBJETHAD_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(NONBJETHAD_TOPTOPLEPHBB,IdxPerm);
	int idx_BJET1_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJET1_TOPTOPLEPHBB,IdxPerm);
	int idx_BJET2_TOPTOPLEPHBB = kfTopTopLepHbb->GetIndex(BJET2_TOPTOPLEPHBB,IdxPerm);

	float NuPx = kfTopTopLepHbb->GetNuPx(IdxPerm,0);
	float NuPy = kfTopTopLepHbb->GetNuPy(IdxPerm,0);
	float NuPz = kfTopTopLepHbb->GetNuPz(IdxPerm,0);
	float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	_TopLepWNu_TOPTOPLEPHBB_p4.SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
	
	_HiggsBJet1_TOPTOPLEPHBB_p4.SetPtEtaPhiE(BJetPt[idx_BJET1_TOPTOPLEPHBB],
						 BJetEta[idx_BJET1_TOPTOPLEPHBB],
						 BJetPhi[idx_BJET1_TOPTOPLEPHBB],
						 BJetE[idx_BJET1_TOPTOPLEPHBB]);
	
	_HiggsBJet2_TOPTOPLEPHBB_p4.SetPtEtaPhiE(BJetPt[idx_BJET2_TOPTOPLEPHBB],
						 BJetEta[idx_BJET2_TOPTOPLEPHBB],
						 BJetPhi[idx_BJET2_TOPTOPLEPHBB],
						 BJetE[idx_BJET2_TOPTOPLEPHBB]);
	
	_TopLepBJet_TOPTOPLEPHBB_p4.SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPTOPLEPHBB],
						 BJetEta[idx_BJETLEP_TOPTOPLEPHBB],
						 BJetPhi[idx_BJETLEP_TOPTOPLEPHBB],
						 BJetE[idx_BJETLEP_TOPTOPLEPHBB]);
	
	_TopHadNonBJet_TOPTOPLEPHBB_p4.SetPtEtaPhiE(NonBJetPt[idx_NONBJETHAD_TOPTOPLEPHBB],
						    NonBJetEta[idx_NONBJETHAD_TOPTOPLEPHBB],
						    NonBJetPhi[idx_NONBJETHAD_TOPTOPLEPHBB],
						    NonBJetE[idx_NONBJETHAD_TOPTOPLEPHBB]);
	
	if( idx_ELECTRON_TOPTOPLEPHBB >= 0 )
	  _TopLepWLep_TOPTOPLEPHBB_p4.SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPTOPLEPHBB],
						   ElectronEta[idx_ELECTRON_TOPTOPLEPHBB],
						   ElectronPhi[idx_ELECTRON_TOPTOPLEPHBB],
						   ElectronE[idx_ELECTRON_TOPTOPLEPHBB]);
	else
	  _TopLepWLep_TOPTOPLEPHBB_p4.SetPtEtaPhiE(MuonPt[idx_MUON_TOPTOPLEPHBB],
						   MuonEta[idx_MUON_TOPTOPLEPHBB],
						   MuonPhi[idx_MUON_TOPTOPLEPHBB],
						   MuonE[idx_MUON_TOPTOPLEPHBB]);
		
	_Higgs_TOPTOPLEPHBB_p4 = _HiggsBJet1_TOPTOPLEPHBB_p4+_HiggsBJet2_TOPTOPLEPHBB_p4;
	_TopLepW_TOPTOPLEPHBB_p4 = _TopLepWLep_TOPTOPLEPHBB_p4+_TopLepWNu_TOPTOPLEPHBB_p4;
	_TopLep_TOPTOPLEPHBB_p4 = _TopLepW_TOPTOPLEPHBB_p4+_TopLepBJet_TOPTOPLEPHBB_p4;
     }

   // tH single top
   if( nLep == 1 && 
       ((nbjets >= 3) || (nbjets == 2 && nnonbjets >= 1)) )
     {		
	kfTopHLepbb->Init(TOPHLEPBB);

	std::vector<float> BJetPt;
	std::vector<float> BJetEta;
	std::vector<float> BJetPhi;
	std::vector<float> BJetE;
	
	std::vector<float> NonBJetPt;
	std::vector<float> NonBJetEta;
	std::vector<float> NonBJetPhi;
	std::vector<float> NonBJetE;
	
	std::vector<float> ElectronPt;
	std::vector<float> ElectronEta;
	std::vector<float> ElectronPhi;
	std::vector<float> ElectronE;
	
	std::vector<float> MuonPt;
	std::vector<float> MuonEta;
	std::vector<float> MuonPhi;
	std::vector<float> MuonE;
	
	for(int i=0;i<_v_ElectronTight->size();i++)
	  {
	     ElectronPt.push_back(_v_ElectronTight->at(i).pt());
	     ElectronEta.push_back(_v_ElectronTight->at(i).eta());
	     ElectronPhi.push_back(_v_ElectronTight->at(i).phi());
	     ElectronE.push_back(_v_ElectronTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_MuonTight->size();i++)
	  {
	     MuonPt.push_back(_v_MuonTight->at(i).pt());
	     MuonEta.push_back(_v_MuonTight->at(i).eta());
	     MuonPhi.push_back(_v_MuonTight->at(i).phi());
	     MuonE.push_back(_v_MuonTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_BJetTight->size();i++)
	  {
	     BJetPt.push_back(_v_BJetTight->at(i).pt());
	     BJetEta.push_back(_v_BJetTight->at(i).eta());
	     BJetPhi.push_back(_v_BJetTight->at(i).phi());
	     BJetE.push_back(_v_BJetTight->at(i).E());
	  }   
	     
	for(int i=0;i<_v_NonBJetTight->size();i++)
	  {
	     if( nbjets == 2 && i == 0 )
	       {		  
		  BJetPt.push_back(_v_NonBJetTight->at(i).pt());
		  BJetEta.push_back(_v_NonBJetTight->at(i).eta());
		  BJetPhi.push_back(_v_NonBJetTight->at(i).phi());
		  BJetE.push_back(_v_NonBJetTight->at(i).E());
	       }
	     else
	       {
		  NonBJetPt.push_back(_v_NonBJetTight->at(i).pt());
		  NonBJetEta.push_back(_v_NonBJetTight->at(i).eta());
		  NonBJetPhi.push_back(_v_NonBJetTight->at(i).phi());
		  NonBJetE.push_back(_v_NonBJetTight->at(i).E());
	       }	     
	  }   
	
	kfTopHLepbb->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
	kfTopHLepbb->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
	kfTopHLepbb->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
	kfTopHLepbb->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
	kfTopHLepbb->SetMet(metpx,metpy);
	
	kfTopHLepbb->Run();

	std::vector<std::pair<float,int> > MVADisc;
	
	int nPerm = kfTopHLepbb->GetNPerm();
	
	for(int ip=0;ip<nPerm;ip++)
	  {
	     float disc = kfTopHLepbb->GetDisc(ip);

	     int idx_ELECTRON_TOPHLEPBB = kfTopHLepbb->GetIndex(ELECTRON_TOPHLEPBB,ip);
	     int idx_MUON_TOPHLEPBB = kfTopHLepbb->GetIndex(MUON_TOPHLEPBB,ip);
	     int idx_BJETLEP_TOPHLEPBB = kfTopHLepbb->GetIndex(BJETLEP_TOPHLEPBB,ip);
	     int idx_BJET1_TOPHLEPBB = kfTopHLepbb->GetIndex(BJET1_TOPHLEPBB,ip);
	     int idx_BJET2_TOPHLEPBB = kfTopHLepbb->GetIndex(BJET2_TOPHLEPBB,ip);

	     float NuPx = kfTopHLepbb->GetNuPx(ip,0);
	     float NuPy = kfTopHLepbb->GetNuPy(ip,0);
	     float NuPz = kfTopHLepbb->GetNuPz(ip,0);
	     float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	     TLorentzVector *TopLepWNu_p4 = new TLorentzVector();
	     TopLepWNu_p4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
	     
	     TLorentzVector *HiggsBJet1_p4 = new TLorentzVector();
	     HiggsBJet1_p4->SetPtEtaPhiE(BJetPt[idx_BJET1_TOPHLEPBB],
					 BJetEta[idx_BJET1_TOPHLEPBB],
					 BJetPhi[idx_BJET1_TOPHLEPBB],
					 BJetE[idx_BJET1_TOPHLEPBB]);
	     
	     TLorentzVector *HiggsBJet2_p4 = new TLorentzVector();
	     HiggsBJet2_p4->SetPtEtaPhiE(BJetPt[idx_BJET2_TOPHLEPBB],
					BJetEta[idx_BJET2_TOPHLEPBB],
					BJetPhi[idx_BJET2_TOPHLEPBB],
					BJetE[idx_BJET2_TOPHLEPBB]);
	     
	     TLorentzVector *TopLepBJet_p4 = new TLorentzVector();
	     TopLepBJet_p4->SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPHLEPBB],
					 BJetEta[idx_BJETLEP_TOPHLEPBB],
					 BJetPhi[idx_BJETLEP_TOPHLEPBB],
					 BJetE[idx_BJETLEP_TOPHLEPBB]);
	     
	     TLorentzVector *TopLepWLep_p4 = new TLorentzVector();
	     if( idx_ELECTRON_TOPHLEPBB >= 0 )
	       TopLepWLep_p4->SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPHLEPBB],
					   ElectronEta[idx_ELECTRON_TOPHLEPBB],
					   ElectronPhi[idx_ELECTRON_TOPHLEPBB],
					   ElectronE[idx_ELECTRON_TOPHLEPBB]);
	     else
	       TopLepWLep_p4->SetPtEtaPhiE(MuonPt[idx_MUON_TOPHLEPBB],
					   MuonEta[idx_MUON_TOPHLEPBB],
					   MuonPhi[idx_MUON_TOPHLEPBB],
					   MuonE[idx_MUON_TOPHLEPBB]);
	     
	     TLorentzVector Higgs_p4 = *HiggsBJet1_p4+*HiggsBJet2_p4;
	     TLorentzVector TopLepW_p4 = *TopLepWLep_p4+*TopLepWNu_p4;
	     TLorentzVector TopLep_p4 = TopLepW_p4+*TopLepBJet_p4;

	     float VarHiggsRecM = Higgs_p4.M();
	     float VarTopLepRecM = TopLep_p4.M();
	     float VarHiggsTopLepRecDr = Higgs_p4.DeltaR(TopLep_p4);
	     float VarTopLepRecPt = TopLep_p4.Pt();
	     
	     TLorentzVector *HiggsT_p4 = new TLorentzVector();
	     HiggsT_p4->SetPxPyPzE(Higgs_p4.Px(),Higgs_p4.Py(),0.,Higgs_p4.Et());
	     
	     TLorentzVector *TopLepWLepT_p4 = new TLorentzVector();
	     TopLepWLepT_p4->SetPxPyPzE(TopLepWLep_p4->Px(),TopLepWLep_p4->Py(),0.,TopLepWLep_p4->Et());
	     
	     TLorentzVector *TopLepWNuT_p4 = new TLorentzVector();
	     TopLepWNuT_p4->SetPxPyPzE(TopLepWNu_p4->Px(),TopLepWNu_p4->Py(),0.,TopLepWNu_p4->Et());
	     
	     TLorentzVector *TopLepBJetT_p4 = new TLorentzVector();
	     TopLepBJetT_p4->SetPxPyPzE(TopLepBJet_p4->Px(),TopLepBJet_p4->Py(),0.,TopLepBJet_p4->Et());
	     
	     TLorentzVector TopLepT_p4 = *TopLepWLepT_p4+*TopLepWNuT_p4+*TopLepBJetT_p4;
	     
	     float VarTopLepRecMT = TopLepT_p4.M();
	     float VarHiggsTopLepRecDphiT = HiggsT_p4->DeltaPhi(TopLepT_p4);
	     float VarTopLepRecPtT = TopLepT_p4.Pt();
	     
	     if( disc < 10E+8 )
	       {
		  MVAFullReco_HiggsRecM_TOPHLEPBB = VarHiggsRecM;
		  MVAFullReco_TopLepRecM_TOPHLEPBB = VarTopLepRecM;
		  MVAFullReco_HiggsTopLepRecDr_TOPHLEPBB = VarHiggsTopLepRecDr;
		  MVAFullReco_TopLepRecPt_TOPHLEPBB = VarTopLepRecPt;
		  
		  MVADisc.push_back(std::make_pair(MVAFullRecoReaderTOPHLEPBB->EvaluateMVA("BDTG method"),ip));
	       }	
	     else
	       {
		  MVAPartReco_HiggsRecM_TOPHLEPBB = VarHiggsRecM;
		  MVAPartReco_TopLepRecMT_TOPHLEPBB = VarTopLepRecMT;
		  MVAPartReco_HiggsTopLepRecDphiT_TOPHLEPBB = VarHiggsTopLepRecDphiT;
		  MVAPartReco_TopLepRecPtT_TOPHLEPBB = VarTopLepRecPtT;
		  
		  MVADisc.push_back(std::make_pair(MVAPartRecoReaderTOPHLEPBB->EvaluateMVA("BDTG method"),ip));
	       }	     

	     delete TopLepWNu_p4;
	     delete HiggsBJet1_p4;
	     delete HiggsBJet2_p4;
	     delete TopLepBJet_p4;
	     delete TopLepWLep_p4;
	     delete HiggsT_p4;
	     delete TopLepWLepT_p4;
	     delete TopLepWNuT_p4;
	     delete TopLepBJetT_p4;
	  }
	
	std::sort(MVADisc.begin(),MVADisc.end(),sortFunc());
	int IdxPerm = MVADisc.at(0).second;

	_chi2_TOPHLEPBB = kfTopHLepbb->GetDisc(IdxPerm);
	_MVA_TOPHLEPBB = MVADisc.at(0).first;
	
	int idx_ELECTRON_TOPHLEPBB = kfTopHLepbb->GetIndex(ELECTRON_TOPHLEPBB,IdxPerm);
	int idx_MUON_TOPHLEPBB = kfTopHLepbb->GetIndex(MUON_TOPHLEPBB,IdxPerm);
	int idx_BJETLEP_TOPHLEPBB = kfTopHLepbb->GetIndex(BJETLEP_TOPHLEPBB,IdxPerm);
	int idx_BJET1_TOPHLEPBB = kfTopHLepbb->GetIndex(BJET1_TOPHLEPBB,IdxPerm);
	int idx_BJET2_TOPHLEPBB = kfTopHLepbb->GetIndex(BJET2_TOPHLEPBB,IdxPerm);

	float NuPx = kfTopHLepbb->GetNuPx(IdxPerm,0);
	float NuPy = kfTopHLepbb->GetNuPy(IdxPerm,0);
	float NuPz = kfTopHLepbb->GetNuPz(IdxPerm,0);
	float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	_TopLepWNu_TOPHLEPBB_p4.SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
	
	_HiggsBJet1_TOPHLEPBB_p4.SetPtEtaPhiE(BJetPt[idx_BJET1_TOPHLEPBB],
					      BJetEta[idx_BJET1_TOPHLEPBB],
					      BJetPhi[idx_BJET1_TOPHLEPBB],
					      BJetE[idx_BJET1_TOPHLEPBB]);
	
	_HiggsBJet2_TOPHLEPBB_p4.SetPtEtaPhiE(BJetPt[idx_BJET2_TOPHLEPBB],
					      BJetEta[idx_BJET2_TOPHLEPBB],
					      BJetPhi[idx_BJET2_TOPHLEPBB],
					      BJetE[idx_BJET2_TOPHLEPBB]);
	
	_TopLepBJet_TOPHLEPBB_p4.SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPHLEPBB],
					      BJetEta[idx_BJETLEP_TOPHLEPBB],
					      BJetPhi[idx_BJETLEP_TOPHLEPBB],
					      BJetE[idx_BJETLEP_TOPHLEPBB]);
	
	if( idx_ELECTRON_TOPHLEPBB >= 0 )
	  _TopLepWLep_TOPHLEPBB_p4.SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPHLEPBB],
						ElectronEta[idx_ELECTRON_TOPHLEPBB],
						ElectronPhi[idx_ELECTRON_TOPHLEPBB],
						ElectronE[idx_ELECTRON_TOPHLEPBB]);
	else
	  _TopLepWLep_TOPHLEPBB_p4.SetPtEtaPhiE(MuonPt[idx_MUON_TOPHLEPBB],
						MuonEta[idx_MUON_TOPHLEPBB],
						MuonPhi[idx_MUON_TOPHLEPBB],
						MuonE[idx_MUON_TOPHLEPBB]);
		
	_Higgs_TOPHLEPBB_p4 = _HiggsBJet1_TOPHLEPBB_p4+_HiggsBJet2_TOPHLEPBB_p4;
	_TopLepW_TOPHLEPBB_p4 = _TopLepWLep_TOPHLEPBB_p4+_TopLepWNu_TOPHLEPBB_p4;
	_TopLep_TOPHLEPBB_p4 = _TopLepW_TOPHLEPBB_p4+_TopLepBJet_TOPHLEPBB_p4;
     }

   // ttbar
   if( nLep == 1 &&
       ((nbjets >= 3 && nnonbjets >= 1) || (nbjets == 2 && nnonbjets >= 2))
     )
     {		
	kfTopTopLepHad->Init(TOPTOPLEPHAD);
	
	std::vector<float> BJetPt;
	std::vector<float> BJetEta;
	std::vector<float> BJetPhi;
	std::vector<float> BJetE;
	
	std::vector<float> NonBJetPt;
	std::vector<float> NonBJetEta;
	std::vector<float> NonBJetPhi;
	std::vector<float> NonBJetE;
	
	std::vector<float> ElectronPt;
	std::vector<float> ElectronEta;
	std::vector<float> ElectronPhi;
	std::vector<float> ElectronE;
	
	std::vector<float> MuonPt;
	std::vector<float> MuonEta;
	std::vector<float> MuonPhi;
	std::vector<float> MuonE;
	
	for(int i=0;i<_v_ElectronTight->size();i++)
	  {
	     ElectronPt.push_back(_v_ElectronTight->at(i).pt());
	     ElectronEta.push_back(_v_ElectronTight->at(i).eta());
	     ElectronPhi.push_back(_v_ElectronTight->at(i).phi());
	     ElectronE.push_back(_v_ElectronTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_MuonTight->size();i++)
	  {
	     MuonPt.push_back(_v_MuonTight->at(i).pt());
	     MuonEta.push_back(_v_MuonTight->at(i).eta());
	     MuonPhi.push_back(_v_MuonTight->at(i).phi());
	     MuonE.push_back(_v_MuonTight->at(i).E());
	  }   
	
	for(int i=0;i<_v_BJetTight->size();i++)
	  {
	     if( i < 2 )
	       {		  
		  BJetPt.push_back(_v_BJetTight->at(i).pt());
		  BJetEta.push_back(_v_BJetTight->at(i).eta());
		  BJetPhi.push_back(_v_BJetTight->at(i).phi());
		  BJetE.push_back(_v_BJetTight->at(i).E());
	       }
	     else
	       {
		  NonBJetPt.push_back(_v_BJetTight->at(i).pt());
		  NonBJetEta.push_back(_v_BJetTight->at(i).eta());
		  NonBJetPhi.push_back(_v_BJetTight->at(i).phi());
		  NonBJetE.push_back(_v_BJetTight->at(i).E());
	       }	     
	  }   
	     
	for(int i=0;i<_v_NonBJetTight->size();i++)
	  {
	     NonBJetPt.push_back(_v_NonBJetTight->at(i).pt());
	     NonBJetEta.push_back(_v_NonBJetTight->at(i).eta());
	     NonBJetPhi.push_back(_v_NonBJetTight->at(i).phi());
	     NonBJetE.push_back(_v_NonBJetTight->at(i).E());
	  }   
	
	kfTopTopLepHad->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
	kfTopTopLepHad->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
	kfTopTopLepHad->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
	kfTopTopLepHad->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
	kfTopTopLepHad->SetMet(metpx,metpy);
	
	kfTopTopLepHad->Run();

	std::vector<std::pair<float,int> > MVADisc;
	
	int nPerm = kfTopTopLepHad->GetNPerm();
	
	for(int ip=0;ip<nPerm;ip++)
	  {
	     float disc = kfTopTopLepHad->GetDisc(ip);

	     int idx_ELECTRON_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(ELECTRON_TOPTOPLEPHAD,ip);
	     int idx_MUON_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(MUON_TOPTOPLEPHAD,ip);
	     int idx_BJETLEP_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(BJETLEP_TOPTOPLEPHAD,ip);
	     int idx_BJETHAD_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(BJETHAD_TOPTOPLEPHAD,ip);
	     int idx_NONBJET1_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(NONBJET1_TOPTOPLEPHAD,ip);
	     int idx_NONBJET2_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(NONBJET2_TOPTOPLEPHAD,ip);

	     float NuPx = kfTopTopLepHad->GetNuPx(ip,0);
	     float NuPy = kfTopTopLepHad->GetNuPy(ip,0);
	     float NuPz = kfTopTopLepHad->GetNuPz(ip,0);
	     float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	     TLorentzVector *TopLepWNu_p4 = new TLorentzVector();
	     TopLepWNu_p4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
	     
	     TLorentzVector *TopHadWNonBJet1_p4 = new TLorentzVector();
	     TopHadWNonBJet1_p4->SetPtEtaPhiE(NonBJetPt[idx_NONBJET1_TOPTOPLEPHAD],
					      NonBJetEta[idx_NONBJET1_TOPTOPLEPHAD],
					      NonBJetPhi[idx_NONBJET1_TOPTOPLEPHAD],
					      NonBJetE[idx_NONBJET1_TOPTOPLEPHAD]);
	     
	     TLorentzVector *TopHadWNonBJet2_p4 = new TLorentzVector();
	     TopHadWNonBJet2_p4->SetPtEtaPhiE(NonBJetPt[idx_NONBJET2_TOPTOPLEPHAD],
					      NonBJetEta[idx_NONBJET2_TOPTOPLEPHAD],
					      NonBJetPhi[idx_NONBJET2_TOPTOPLEPHAD],
					      NonBJetE[idx_NONBJET2_TOPTOPLEPHAD]);
	     
	     TLorentzVector *TopLepBJet_p4 = new TLorentzVector();
	     TopLepBJet_p4->SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPTOPLEPHAD],
					 BJetEta[idx_BJETLEP_TOPTOPLEPHAD],
					 BJetPhi[idx_BJETLEP_TOPTOPLEPHAD],
					 BJetE[idx_BJETLEP_TOPTOPLEPHAD]);
	     
	     TLorentzVector *TopHadBJet_p4 = new TLorentzVector();
	     TopHadBJet_p4->SetPtEtaPhiE(BJetPt[idx_BJETHAD_TOPTOPLEPHAD],
					 BJetEta[idx_BJETHAD_TOPTOPLEPHAD],
					 BJetPhi[idx_BJETHAD_TOPTOPLEPHAD],
					 BJetE[idx_BJETHAD_TOPTOPLEPHAD]);
	     
	     TLorentzVector *TopLepWLep_p4 = new TLorentzVector();
	     if( idx_ELECTRON_TOPTOPLEPHAD >= 0 )
	       TopLepWLep_p4->SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPTOPLEPHAD],
					   ElectronEta[idx_ELECTRON_TOPTOPLEPHAD],
					   ElectronPhi[idx_ELECTRON_TOPTOPLEPHAD],
					   ElectronE[idx_ELECTRON_TOPTOPLEPHAD]);
	     else
	       TopLepWLep_p4->SetPtEtaPhiE(MuonPt[idx_MUON_TOPTOPLEPHAD],
					   MuonEta[idx_MUON_TOPTOPLEPHAD],
					   MuonPhi[idx_MUON_TOPTOPLEPHAD],
					   MuonE[idx_MUON_TOPTOPLEPHAD]);
	     
	     TLorentzVector TopHadW_p4 = *TopHadWNonBJet1_p4+*TopHadWNonBJet2_p4;
	     TLorentzVector TopLepW_p4 = *TopLepWLep_p4+*TopLepWNu_p4;
	     TLorentzVector TopLep_p4 = TopLepW_p4+*TopLepBJet_p4;
	     TLorentzVector TopHad_p4 = TopHadW_p4+*TopHadBJet_p4;

	     float VarTopHadRecM = TopHad_p4.M();
	     float VarTopLepRecM = TopLep_p4.M();
	     float VarTopLepTopHadRecDr = TopLep_p4.DeltaR(TopHad_p4);
	     float VarTopLepRecPt = TopLep_p4.Pt();
	     
	     TLorentzVector *TopHadT_p4 = new TLorentzVector();
	     TopHadT_p4->SetPxPyPzE(TopHad_p4.Px(),TopHad_p4.Py(),0.,TopHad_p4.Et());
	     
	     TLorentzVector *TopLepWLepT_p4 = new TLorentzVector();
	     TopLepWLepT_p4->SetPxPyPzE(TopLepWLep_p4->Px(),TopLepWLep_p4->Py(),0.,TopLepWLep_p4->Et());
	     
	     TLorentzVector *TopLepWNuT_p4 = new TLorentzVector();
	     TopLepWNuT_p4->SetPxPyPzE(TopLepWNu_p4->Px(),TopLepWNu_p4->Py(),0.,TopLepWNu_p4->Et());
	     
	     TLorentzVector *TopLepBJetT_p4 = new TLorentzVector();
	     TopLepBJetT_p4->SetPxPyPzE(TopLepBJet_p4->Px(),TopLepBJet_p4->Py(),0.,TopLepBJet_p4->Et());
	     
	     TLorentzVector TopLepT_p4 = *TopLepWLepT_p4+*TopLepWNuT_p4+*TopLepBJetT_p4;
	     
	     float VarTopLepRecMT = TopLepT_p4.M();
	     float VarTopLepTopHadRecDphiT = TopLepT_p4.DeltaPhi(*TopHadT_p4);
	     float VarTopLepRecPtT = TopLepT_p4.Pt();
	     
	     if( disc < 10E+8 )
	       {
		  MVAFullReco_TopHadRecM_TOPTOPLEPHAD = VarTopHadRecM;
		  MVAFullReco_TopLepRecM_TOPTOPLEPHAD = VarTopLepRecM;
		  MVAFullReco_TopLepTopHadRecDr_TOPTOPLEPHAD = VarTopLepTopHadRecDr;
		  MVAFullReco_TopLepRecPt_TOPTOPLEPHAD = VarTopLepRecPt;
		  
		  MVADisc.push_back(std::make_pair(MVAFullRecoReaderTOPTOPLEPHAD->EvaluateMVA("BDTG method"),ip));
	       }	
	     else
	       {
		  MVAPartReco_TopHadRecM_TOPTOPLEPHAD = VarTopHadRecM;
		  MVAPartReco_TopLepRecMT_TOPTOPLEPHAD = VarTopLepRecMT;
		  MVAPartReco_TopLepTopHadRecDphiT_TOPTOPLEPHAD = VarTopLepTopHadRecDphiT;
		  MVAPartReco_TopLepRecPtT_TOPTOPLEPHAD = VarTopLepRecPtT;
		  
		  MVADisc.push_back(std::make_pair(MVAPartRecoReaderTOPTOPLEPHAD->EvaluateMVA("BDTG method"),ip));
	       }	     

	     delete TopLepWNu_p4;
	     delete TopHadWNonBJet1_p4;
	     delete TopHadWNonBJet2_p4;
	     delete TopLepBJet_p4;
	     delete TopHadBJet_p4;
	     delete TopLepWLep_p4;
	     delete TopHadT_p4;
	     delete TopLepWLepT_p4;
	     delete TopLepWNuT_p4;
	     delete TopLepBJetT_p4;
	  }
	
	std::sort(MVADisc.begin(),MVADisc.end(),sortFunc());
	int IdxPerm = MVADisc.at(0).second;
	
	_chi2_TOPTOPLEPHAD = kfTopTopLepHad->GetDisc(IdxPerm);
	_MVA_TOPTOPLEPHAD = MVADisc.at(0).first;

	int idx_ELECTRON_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(ELECTRON_TOPTOPLEPHAD,IdxPerm);
	int idx_MUON_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(MUON_TOPTOPLEPHAD,IdxPerm);
	int idx_BJETLEP_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(BJETLEP_TOPTOPLEPHAD,IdxPerm);
	int idx_BJETHAD_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(BJETHAD_TOPTOPLEPHAD,IdxPerm);
	int idx_NONBJET1_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(NONBJET1_TOPTOPLEPHAD,IdxPerm);
	int idx_NONBJET2_TOPTOPLEPHAD = kfTopTopLepHad->GetIndex(NONBJET2_TOPTOPLEPHAD,IdxPerm);

	float NuPx = kfTopTopLepHad->GetNuPx(IdxPerm,0);
	float NuPy = kfTopTopLepHad->GetNuPy(IdxPerm,0);
	float NuPz = kfTopTopLepHad->GetNuPz(IdxPerm,0);
	float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

	_TopLepWNu_TOPTOPLEPHAD_p4.SetPxPyPzE(NuPx,NuPy,NuPz,NuE);

	_TopHadWNonBJet1_TOPTOPLEPHAD_p4.SetPtEtaPhiE(NonBJetPt[idx_NONBJET1_TOPTOPLEPHAD],
						      NonBJetEta[idx_NONBJET1_TOPTOPLEPHAD],
						      NonBJetPhi[idx_NONBJET1_TOPTOPLEPHAD],
						      NonBJetE[idx_NONBJET1_TOPTOPLEPHAD]);
	
	_TopHadWNonBJet2_TOPTOPLEPHAD_p4.SetPtEtaPhiE(NonBJetPt[idx_NONBJET2_TOPTOPLEPHAD],
						      NonBJetEta[idx_NONBJET2_TOPTOPLEPHAD],
						      NonBJetPhi[idx_NONBJET2_TOPTOPLEPHAD],
						      NonBJetE[idx_NONBJET2_TOPTOPLEPHAD]);
	
	_TopLepBJet_TOPTOPLEPHAD_p4.SetPtEtaPhiE(BJetPt[idx_BJETLEP_TOPTOPLEPHAD],
						 BJetEta[idx_BJETLEP_TOPTOPLEPHAD],
						 BJetPhi[idx_BJETLEP_TOPTOPLEPHAD],
						 BJetE[idx_BJETLEP_TOPTOPLEPHAD]);
	
	_TopHadBJet_TOPTOPLEPHAD_p4.SetPtEtaPhiE(BJetPt[idx_BJETHAD_TOPTOPLEPHAD],
						 BJetEta[idx_BJETHAD_TOPTOPLEPHAD],
						 BJetPhi[idx_BJETHAD_TOPTOPLEPHAD],
						 BJetE[idx_BJETHAD_TOPTOPLEPHAD]);
	
	if( idx_ELECTRON_TOPTOPLEPHAD >= 0 )
	  _TopLepWLep_TOPTOPLEPHAD_p4.SetPtEtaPhiE(ElectronPt[idx_ELECTRON_TOPTOPLEPHAD],
						   ElectronEta[idx_ELECTRON_TOPTOPLEPHAD],
						   ElectronPhi[idx_ELECTRON_TOPTOPLEPHAD],
						   ElectronE[idx_ELECTRON_TOPTOPLEPHAD]);
	else
	  _TopLepWLep_TOPTOPLEPHAD_p4.SetPtEtaPhiE(MuonPt[idx_MUON_TOPTOPLEPHAD],
						   MuonEta[idx_MUON_TOPTOPLEPHAD],
						   MuonPhi[idx_MUON_TOPTOPLEPHAD],
						   MuonE[idx_MUON_TOPTOPLEPHAD]);
		
	_TopHadW_TOPTOPLEPHAD_p4 = _TopHadWNonBJet1_TOPTOPLEPHAD_p4+_TopHadWNonBJet2_TOPTOPLEPHAD_p4;
	_TopLepW_TOPTOPLEPHAD_p4 = _TopLepWLep_TOPTOPLEPHAD_p4+_TopLepWNu_TOPTOPLEPHAD_p4;
	_TopLep_TOPTOPLEPHAD_p4 = _TopLepW_TOPTOPLEPHAD_p4+_TopLepBJet_TOPTOPLEPHAD_p4;
	_TopHad_TOPTOPLEPHAD_p4 = _TopHadW_TOPTOPLEPHAD_p4+_TopHadBJet_TOPTOPLEPHAD_p4;
     }
}

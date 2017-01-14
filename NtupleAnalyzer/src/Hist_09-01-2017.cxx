#include "../include/Hist.h"

#include "TRandom3.h"

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/type_traits.hpp>

Hist::Hist(std::string home)
{
   help = new Helper();
   
   _home = home;
   
//   std::string foutlog = "log/"+std::string(flog);
   std::string foutlog = "log/list.txt";
   _fevc = fopen(foutlog.c_str(),"w");
//   std::string foutlogVal = "log/"+std::string(flog)+".val";
   std::string foutlogVal = "log/list.val";
   _fevcVal.open(foutlogVal.c_str());
      
   _v_ElectronTight = new std::vector<Electron>;
   _v_ElectronLoose = new std::vector<Electron>;
   _v_MuonTight = new std::vector<Muon>;
   _v_MuonLoose = new std::vector<Muon>;
   _v_JetTight = new std::vector<Jet>;
   _v_BJetTight = new std::vector<Jet>;
   _v_NonBJetTight = new std::vector<Jet>;
   
   _v_Lepton = new std::vector<Lepton>;
}

Hist::~Hist()
{
   delete _v_ElectronTight;
   delete _v_ElectronLoose;
   delete _v_MuonTight;
   delete _v_MuonLoose;
   delete _v_JetTight;
   delete _v_BJetTight;
   delete _v_NonBJetTight;
   
   delete _v_Lepton;

   delete help;
   delete _trec;
   delete _mva;
}

void Hist::init()
{
   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);
   rnd = new TRandom3(666);
   
//   std::string foutput = _home+"/"+std::string(flog)+".root";
   std::string foutput = std::string(flog)+".root";
   _fout = new TFile(foutput.c_str(),"RECREATE");

   _trec = new TopReco();
   _trec->init();

   _fout->cd();
   
   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);
   _mva = new ApplyMVA(_home);
   _mva->init();

   //std::cout << "CIRKOVIC EXIT LAST" << std::endl;    exit(0);
   
   //_h_PassSel_all = new TH1D("h_PassSel_all","h_PassSel_all",3,0.,3.);
   //_h_PassSel_e = new TH1D("h_PassSel_e","h_PassSel_e",3,0.,3.);
   //_h_PassSel_m = new TH1D("h_PassSel_m","h_PassSel_m",3,0.,3.);

   hname.clear(); 

   //histname_n = 38;
   histname_n = 62;
   histname[0] = "h_chi2_TOPTOPLEPHBB_";
   histname[1] = "h_chi2_TOPHLEPBB_";
   histname[2] = "h_chi2_TOPTOPLEPHAD_";
   histname[3] = "h_MVA_TOPTOPLEPHBB_";
   histname[4] = "h_MVA_TOPHLEPBB_";
   histname[5] = "h_MVA_TOPTOPLEPHAD_";
   histname[6] = "h_HiggsMass_TOPTOPLEPHBB_";
   histname[7] = "h_HiggsMass_TOPHLEPBB_";
   histname[8] = "h_TopHadMass_TOPTOPLEPHAD_";
   histname[9] = "h_LepCharge_";
   histname[10] = "h_HiggsEta_TOPTOPLEPHBB_";
   histname[11] = "h_HiggsEta_TOPHLEPBB_";
   histname[12] = "h_TopLepMass_TOPTOPLEPHBB_";
   histname[13] = "h_TopLepMass_TOPHLEPBB_";
   histname[14] = "h_TopLepMass_TOPTOPLEPHAD_";
   histname[15] = "h_TopLepPt_TOPTOPLEPHBB_";
   histname[16] = "h_TopLepPt_TOPHLEPBB_";
   histname[17] = "h_TopLepPt_TOPTOPLEPHAD_";
   histname[18] = "h_TopLepEta_TOPTOPLEPHBB_";
   histname[19] = "h_TopLepEta_TOPHLEPBB_";
   histname[20] = "h_TopLepEta_TOPTOPLEPHAD_";
   histname[21] = "h_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_";
   histname[22] = "h_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_";
   histname[23] = "h_TopLepHiggsDr_TOPTOPLEPHBB_";
   histname[24] = "h_TopLepHiggsDr_TOPHLEPBB_";
   histname[25] = "h_TopLepTopHadDr_TOPTOPLEPHAD_";
   /*
   histname[26] = "h_MVAb3j4HutST_";
   histname[27] = "h_MVAb3j4HctST_";
   histname[28] = "h_MVAb3j4HutTT_";
   histname[29] = "h_MVAb3j4HctTT_";
   histname[30] = "h_MVAb3j3HutST_";
   histname[31] = "h_MVAb3j3HctST_";
   histname[32] = "h_MVAb3j3HutTT_";
   histname[33] = "h_MVAb3j3HctTT_";
   histname[34] = "h_MVAb2j4HutST_";
   histname[35] = "h_MVAb2j4HctST_";
   histname[36] = "h_MVAb2j4HutTT_";
   histname[37] = "h_MVAb2j4HctTT_";
   */
   histname[26] = "h_MVAb0j3HutST_";
   histname[27] = "h_MVAb0j3HctST_";
   histname[28] = "h_MVAb0j3HutTT_";
   histname[29] = "h_MVAb0j3HctTT_";
   histname[30] = "h_MVAb1j3HutST_";
   histname[31] = "h_MVAb1j3HctST_";
   histname[32] = "h_MVAb1j3HutTT_";
   histname[33] = "h_MVAb1j3HctTT_";
   histname[34] = "h_MVAb2j3HutST_";
   histname[35] = "h_MVAb2j3HctST_";
   histname[36] = "h_MVAb2j3HutTT_";
   histname[37] = "h_MVAb2j3HctTT_";
   histname[38] = "h_MVAb3j3HutST_";
   histname[39] = "h_MVAb3j3HctST_";
   histname[40] = "h_MVAb3j3HutTT_";
   histname[41] = "h_MVAb3j3HctTT_";
   histname[42] = "h_MVAb0j4HutST_";
   histname[43] = "h_MVAb0j4HctST_";
   histname[44] = "h_MVAb0j4HutTT_";
   histname[45] = "h_MVAb0j4HctTT_";
   histname[46] = "h_MVAb1j4HutST_";
   histname[47] = "h_MVAb1j4HctST_";
   histname[48] = "h_MVAb1j4HutTT_";
   histname[49] = "h_MVAb1j4HctTT_";
   histname[50] = "h_MVAb2j4HutST_";
   histname[51] = "h_MVAb2j4HctST_";
   histname[52] = "h_MVAb2j4HutTT_";
   histname[53] = "h_MVAb2j4HctTT_";
   histname[54] = "h_MVAb3j4HutST_";
   histname[55] = "h_MVAb3j4HctST_";
   histname[56] = "h_MVAb3j4HutTT_";
   histname[57] = "h_MVAb3j4HctTT_";
   histname[58] = "h_MVAb4j4HutST_";
   histname[59] = "h_MVAb4j4HctST_";
   histname[60] = "h_MVAb4j4HutTT_";
   histname[61] = "h_MVAb4j4HctTT_";

   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);

   type_n = 1;
//   type[0] = "nonQCD";
//   type[1] = "QCD";
   type[0] = "ALL";

   /*
   sel_n = 3;
   sel[0] = "b2j4";
   sel[1] = "b3j3";
   sel[2] = "b3j4";
   */

   sel_n = 9;
   sel[0] = "b0j3";
   sel[1] = "b1j3";
   sel[2] = "b2j3";
   sel[3] = "b3j3";
   sel[4] = "b0j4";
   sel[5] = "b1j4";
   sel[6] = "b2j4";
   sel[7] = "b3j4";
   sel[8] = "b4j4";

   chan_n = 3;
   chan[0] = "all";
   chan[1] = "e";
   chan[2] = "m";
   
   sys_low_n = 1;
   sys_low[0]   = "";
//   sys_low[1]   = "_pu_low";

   sys_up_n = sys_low_n-1;
//   sys_up[0]   = "_pu_up";

   sys_n = sys_low_n + sys_up_n;
   
   for(int is1=0;is1<sys_low_n;is1++)
     {
	sys[is1] = sys_low[is1];
     }   
   for(int is2=0;is2<sys_up_n;is2++)
     {
	sys[sys_low_n+is2] = sys_up[is2];
     }   

   for(int ic=0;ic<chan_n;ic++)
     {
	for(int is=0;is<sel_n;is++)
	  {
	     //std::string trName = "tr_"+chan[ic]+"_"+sel[is];
	     std::string trName = "sf_"+chan[ic]+"_"+sel[is];
	     _trout[ic][is] = new TTree(trName.c_str(),trName.c_str());
	     
/*
	     _trout[ic][is]->Branch("weight",&m_weight,"weight/D");
*/
	     _trout[ic][is]->Branch("wgt_csv",&_wgt_csv,"wgt_csv/D");
         /*
	     _trout[ic][is]->Branch("HiggsMass_TOPTOPLEPHBB",&m_HiggsMass_TOPTOPLEPHBB,"HiggsMass_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("HiggsMass_TOPHLEPBB",&m_HiggsMass_TOPHLEPBB,"HiggsMass_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopHadMass_TOPTOPLEPHAD",&m_TopHadMass_TOPTOPLEPHAD,"TopHadMass_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("chi2_TOPTOPLEPHBB",&m_chi2_TOPTOPLEPHBB,"chi2_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("chi2_TOPHLEPBB",&m_chi2_TOPHLEPBB,"chi2_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("chi2_TOPTOPLEPHAD",&m_chi2_TOPTOPLEPHAD,"chi2_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("MVA_TOPTOPLEPHBB",&m_MVA_TOPTOPLEPHBB,"MVA_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("MVA_TOPHLEPBB",&m_MVA_TOPHLEPBB,"MVA_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("MVA_TOPTOPLEPHAD",&m_MVA_TOPTOPLEPHAD,"MVA_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("LepCharge",&m_LepCharge,"LepCharge/I");
	     _trout[ic][is]->Branch("HiggsEta_TOPTOPLEPHBB",&m_HiggsEta_TOPTOPLEPHBB,"HiggsEta_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("HiggsEta_TOPHLEPBB",&m_HiggsEta_TOPHLEPBB,"HiggsEta_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepMass_TOPTOPLEPHBB",&m_TopLepMass_TOPTOPLEPHBB,"TopLepMass_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("TopLepMass_TOPHLEPBB",&m_TopLepMass_TOPHLEPBB,"TopLepMass_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepMass_TOPTOPLEPHAD",&m_TopLepMass_TOPTOPLEPHAD,"TopLepMass_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("TopLepPt_TOPTOPLEPHBB",&m_TopLepPt_TOPTOPLEPHBB,"TopLepPt_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("TopLepPt_TOPHLEPBB",&m_TopLepPt_TOPHLEPBB,"TopLepPt_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepPt_TOPTOPLEPHAD",&m_TopLepPt_TOPTOPLEPHAD,"TopLepPt_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("TopLepEta_TOPTOPLEPHBB",&m_TopLepEta_TOPTOPLEPHBB,"TopLepEta_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("TopLepEta_TOPHLEPBB",&m_TopLepEta_TOPHLEPBB,"TopLepEta_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepEta_TOPTOPLEPHAD",&m_TopLepEta_TOPTOPLEPHAD,"TopLepEta_TOPTOPLEPHAD/D");
	     _trout[ic][is]->Branch("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB,"HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB,"HiggsBJet1HiggsBJet2Dr_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepHiggsDr_TOPTOPLEPHBB",&m_TopLepHiggsDr_TOPTOPLEPHBB,"TopLepHiggsDr_TOPTOPLEPHBB/D");
	     _trout[ic][is]->Branch("TopLepHiggsDr_TOPHLEPBB",&m_TopLepHiggsDr_TOPHLEPBB,"TopLepHiggsDr_TOPHLEPBB/D");
	     _trout[ic][is]->Branch("TopLepTopHadDr_TOPTOPLEPHAD",&m_TopLepTopHadDr_TOPTOPLEPHAD,"TopLepTopHadDr_TOPTOPLEPHAD/D");
         */
         /*
         _trout[ic][is]->Branch("HiggsMass_TOPHLEPBB",&m_HiggsMass_TOPHLEPBB,"HiggsMass_TOPHLEPBB/D");
         _trout[ic][is]->Branch("MVA_TOPHLEPBB",&m_MVA_TOPHLEPBB,"MVA_TOPHLEPBB/D");
         _trout[ic][is]->Branch("HiggsEta_TOPHLEPBB",&m_HiggsEta_TOPHLEPBB,"HiggsEta_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopLepMass_TOPHLEPBB",&m_TopLepMass_TOPHLEPBB,"TopLepMass_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopLepPt_TOPHLEPBB",&m_TopLepPt_TOPHLEPBB,"TopLepPt_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopLepEta_TOPHLEPBB",&m_TopLepEta_TOPHLEPBB,"TopLepEta_TOPHLEPBB/D");
         _trout[ic][is]->Branch("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB,"HiggsBJet1HiggsBJet2Dr_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopLepHiggsDr_TOPHLEPBB",&m_TopLepHiggsDr_TOPHLEPBB,"TopLepHiggsDr_TOPHLEPBB/D");
         _trout[ic][is]->Branch("HiggsBJet1CSVv2_TOPHLEPBB",&m_HiggsBJet1CSVv2_TOPHLEPBB,"HiggsBJet1CSVv2_TOPHLEPBB/D");
         _trout[ic][is]->Branch("HiggsBJet2CSVv2_TOPHLEPBB",&m_HiggsBJet2CSVv2_TOPHLEPBB,"HiggsBJet2CSVv2_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopLepBJetCSVv2_TOPHLEPBB",&m_TopLepBJetCSVv2_TOPHLEPBB,"TopLepBJetCSVv2_TOPHLEPBB/D");
         _trout[ic][is]->Branch("TopHadMass_TOPTOPLEPHAD",&m_TopHadMass_TOPTOPLEPHAD,"TopHadMass_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("MVA_TOPTOPLEPHAD",&m_MVA_TOPTOPLEPHAD,"MVA_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopLepMass_TOPTOPLEPHAD",&m_TopLepMass_TOPTOPLEPHAD,"TopLepMass_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopLepTopHadDr_TOPTOPLEPHAD",&m_TopLepTopHadDr_TOPTOPLEPHAD,"TopLepTopHadDr_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopLepBJetCSVv2_TOPTOPLEPHAD",&m_TopLepBJetCSVv2_TOPTOPLEPHAD,"TopLepBJetCSVv2_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopHadBJetCSVv2_TOPTOPLEPHAD",&m_TopHadBJetCSVv2_TOPTOPLEPHAD,"TopHadBJetCSVv2_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&m_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD,"TopHadWNonBJet1CSVv2_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&m_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD,"TopHadWNonBJet2CSVv2_TOPTOPLEPHAD/D");
         _trout[ic][is]->Branch("HiggsMass_TOPTOPLEPHBB",&m_HiggsMass_TOPTOPLEPHBB,"HiggsMass_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("MVA_TOPTOPLEPHBB",&m_MVA_TOPTOPLEPHBB,"MVA_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("TopLepMass_TOPTOPLEPHBB",&m_TopLepMass_TOPTOPLEPHBB,"TopLepMass_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB,"HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("TopLepHiggsDr_TOPTOPLEPHBB",&m_TopLepHiggsDr_TOPTOPLEPHBB,"TopLepHiggsDr_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("HiggsBJet1CSVv2_TOPTOPLEPHBB",&m_HiggsBJet1CSVv2_TOPTOPLEPHBB,"HiggsBJet1CSVv2_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("HiggsBJet2CSVv2_TOPTOPLEPHBB",&m_HiggsBJet2CSVv2_TOPTOPLEPHBB,"HiggsBJet2CSVv2_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("TopLepBJetCSVv2_TOPTOPLEPHBB",&m_TopLepBJetCSVv2_TOPTOPLEPHBB,"TopLepBJetCSVv2_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&m_TopHadNonBJetCSVv2_TOPTOPLEPHBB,"TopHadNonBJetCSVv2_TOPTOPLEPHBB/D");
         _trout[ic][is]->Branch("LepCharge",&m_LepCharge,"LepCharge/I");

         if (string(sel[is]) == string("b0j3")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab0j3HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab0j3HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab0j3HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab0j3HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b1j3")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab1j3HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab1j3HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab1j3HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab1j3HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b2j3")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab2j3HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab2j3HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab2j3HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab2j3HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b3j3")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab3j3HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab3j3HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab3j3HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab3j3HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b0j4")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab0j4HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab0j4HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab0j4HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab0j4HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b1j4")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab1j4HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab1j4HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab1j4HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab1j4HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b2j4")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab2j4HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab2j4HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab2j4HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab2j4HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b3j4")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab3j4HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab3j4HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab3j4HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab3j4HctTT,"mvaHctTT/D");
         }
         else if (string(sel[is]) == string("b4j4")) {
             _trout[ic][is]->Branch("mvaHutST",&_mvab4j4HutST,"mvaHutST/D");
             _trout[ic][is]->Branch("mvaHctST",&_mvab4j4HctST,"mvaHctST/D");
             _trout[ic][is]->Branch("mvaHutTT",&_mvab4j4HutTT,"mvaHutTT/D");
             _trout[ic][is]->Branch("mvaHctTT",&_mvab4j4HctTT,"mvaHctTT/D");
         }
         _trout[ic][is]->Branch("njets",&njets,"njets/int");
         _trout[ic][is]->Branch("nbjets",&nbjets,"nbjets/int");
        */
	  }
     }
   /*
   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);

   _s_Hist = new std::vector<std::pair<std::vector<std::string>,double*> >();
   _m1d_Hist = new std::map<std::string, TH1D*>();

   std::vector<double*> set_hist;
   set_hist.clear();

   // Ranges
   set_hist.push_back(RANGE::set_chi2);
   set_hist.push_back(RANGE::set_chi2);
   set_hist.push_back(RANGE::set_chi2);
   set_hist.push_back(RANGE::set_bMVA);
   set_hist.push_back(RANGE::set_bMVA);
   set_hist.push_back(RANGE::set_bMVA);
   set_hist.push_back(RANGE::set_H_m);
   set_hist.push_back(RANGE::set_H_m);
   set_hist.push_back(RANGE::set_top_m);
   set_hist.push_back(RANGE::set_l_charge);
   set_hist.push_back(RANGE::set_H_eta);
   set_hist.push_back(RANGE::set_H_eta);
   set_hist.push_back(RANGE::set_top_m);
   set_hist.push_back(RANGE::set_top_m);
   set_hist.push_back(RANGE::set_top_m);
   set_hist.push_back(RANGE::set_top_pt);
   set_hist.push_back(RANGE::set_top_pt);
   set_hist.push_back(RANGE::set_top_pt);
   set_hist.push_back(RANGE::set_top_eta);
   set_hist.push_back(RANGE::set_top_eta);
   set_hist.push_back(RANGE::set_top_eta);
   set_hist.push_back(RANGE::set_Hb1_Hb2_dr);
   set_hist.push_back(RANGE::set_Hb1_Hb2_dr);
   set_hist.push_back(RANGE::set_W_topb_dr);
   set_hist.push_back(RANGE::set_W_topb_dr);
   set_hist.push_back(RANGE::set_W_topb_dr);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);
   set_hist.push_back(RANGE::set_MVA);

   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);
   
   std::string titl;

   for(int s=0;s<sel_n;s++)
     {
	for(int t=0;t<type_n;t++)
	  {
	     for(int c=0;c<chan_n;c++)
	       {
		  for(int h=0;h<histname_n;h++)
		    {
		       std::string hn = histname[h]+chan[c]+"_"+sel[s]+"_"+type[t];
		       hname.push_back(hn);
				 
		       for(int s=0;s<sys_n;s++)
			 {
			    titl = hn+sys[s];
			    std::vector<std::string> names;
			    names.clear();
			    names.push_back(titl);
			    names.push_back(sys[s]);
			    
			    _s_Hist->push_back(std::make_pair(names,set_hist.at(h)));
			 }
		    }		  
	       }
	  }
     }   

   //std::cout << "CIRKOVIC BREAK" << std::endl;    exit(0);
   
   for(int i=0;i<_s_Hist->size();i++)
     {
	TH1D *_h1d = new TH1D(_s_Hist->at(i).first.at(0).c_str(),
			      _s_Hist->at(i).first.at(0).c_str(),
			      _s_Hist->at(i).second[0],
			      _s_Hist->at(i).second[1],
			      _s_Hist->at(i).second[2]);
	
	_h1d->Sumw2();
	
	_m1d_Hist->insert(std::pair<std::string,TH1D*>(_s_Hist->at(i).first.at(0),_h1d));
     }

   for(int ss=0;ss<sel_n;ss++)
     {
	for(int t=0;t<type_n;t++)
	  {
	     for(int c=0;c<chan_n;c++)
	       {
		  for(int s=0;s<sys_n;s++)
		    {
		       for(int h=0;h<histname_n;h++)
			 {
			    histNAMES[c][t][ss][h][s] =
			      histname[h]+chan[c]+"_"+sel[ss]+"_"+type[t]+sys[s];
			 }		       
		    }		  
	       }
	  }
     }

   std::cout << "Hist initialisation done" << std::endl;
   */
}

void Hist::fill()
{
   float w = _v_Event->at(0).mc_weight();
   
   bool isTrigElec = _v_Event->at(0).isTrigElec();
   bool isTrigMuon = _v_Event->at(0).isTrigMuon();
   bool isSignal = _v_Event->at(0).isSignal();
   
   bool isData = _v_Event->at(0).isData();                                                                              

   float xsec = _v_Event->at(0).xsec();
   float noe = _v_Event->at(0).noe();

   float wil = 12900.*xsec/noe;
   if( !isData ) w *= wil;

/*   if( !isSignal )
     {	
	if( !isTrigElec && !isTrigMuon ) return;
     }*/

   bool passTrig = 0;
   if( isData )
     {
	if( !isTrigElec && isTrigMuon ) passTrig = 1;
	else if( isTrigElec && !isTrigMuon ) passTrig = 1;
     }   
   else passTrig = 1;
   
   m_weight = w;

   _v_ElectronTight->clear();
   _v_ElectronLoose->clear();
   _v_MuonTight->clear();
   _v_MuonLoose->clear();
   _v_JetTight->clear();
   _v_BJetTight->clear();
   _v_NonBJetTight->clear();
   
   for(int i=0;i<_v_Electron->size();i++)
     {
	if( fabs(_v_Electron->at(i).eta()) > 2.4 ) continue;

	if( _v_Electron->at(i).isTight() )
	  {		 		 
	     _v_ElectronTight->push_back(_v_Electron->at(i));
	  }	
	else
	  {
	     _v_ElectronLoose->push_back(_v_Electron->at(i));
	  }	
     }       

   for(int i=0;i<_v_Muon->size();i++)
     {
	if( fabs(_v_Muon->at(i).eta()) > 2.4 ) continue;

	if( _v_Muon->at(i).isTight() )
	  {	     
	     _v_MuonTight->push_back(_v_Muon->at(i));
	  }	
	else
	  {
	     _v_MuonLoose->push_back(_v_Muon->at(i));
	  }	
     }      

   std::sort(_v_Jet->begin(),_v_Jet->end(),Jet::sortCSVv2Predicate);
   
   for(int i=0;i<_v_Jet->size();i++)
     {
	if( fabs(_v_Jet->at(i).eta()) > 2.4 ) continue;

	if( _v_Jet->at(i).isTight() )
	  {	     
	     _v_JetTight->push_back(_v_Jet->at(i));
	     
	     //	if( _v_Jet->at(i).CSVv2() >= 0.605 ) _v_BJetTight->push_back(_v_Jet->at(i));
	     //	if( i < 2 && _v_Jet->at(i).CSVv2() >= 0.890 ) _v_BJetTight->push_back(_v_Jet->at(i));
	     if( _v_Jet->at(i).CSVv2() >= 0.800 ) _v_BJetTight->push_back(_v_Jet->at(i));
	     else _v_NonBJetTight->push_back(_v_Jet->at(i));
	     //	if( _v_Jet->at(i).CSVv2() >= 0.970 ) _v_BJetTight->push_back(_v_Jet->at(i));
	  }	
     }      
   /*
   int id = _v_Event->at(0).id();
   int run = _v_Event->at(0).run();
   int lumi = _v_Event->at(0).lumi();
   
   float metpt = _v_Event->at(0).metpt();
   float metphi = _v_Event->at(0).metphi();
   */
   //int njets = _v_JetTight->size();
   njets = _v_JetTight->size();
   //int nbjets = _v_BJetTight->size();
   nbjets = _v_BJetTight->size();
   
   int nElecLoose = _v_ElectronLoose->size();
   int nMuonLoose = _v_MuonLoose->size();

   int nElecTight = _v_ElectronTight->size();
   int nMuonTight = _v_MuonTight->size();
   
   bool is1L = (nElecTight+nMuonTight == 1);
   bool is1L_looseVETO = (is1L && nElecLoose == 0 && nMuonLoose == 0);
   bool is1E = (nElecTight == 1 && is1L_looseVETO);
   bool is1M = (nMuonTight == 1 && is1L_looseVETO);
   
   passSel_all = 0x0;
   /*
   passSel_all |= (is1L_looseVETO && njets >= 3 && nbjets == 2)   << 0;
   passSel_all |= (is1L_looseVETO && njets == 3 && nbjets == 3)   << 1;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 3)   << 2;
   */
   passSel_all |= (is1L_looseVETO && njets == 3 && nbjets == 0)   << 0;
   passSel_all |= (is1L_looseVETO && njets == 3 && nbjets == 1)   << 1;
   passSel_all |= (is1L_looseVETO && njets == 3 && nbjets == 2)   << 2;
   passSel_all |= (is1L_looseVETO && njets == 3 && nbjets == 3)   << 3;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 0)   << 4;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 1)   << 5;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 2)   << 6;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 3)   << 7;
   passSel_all |= (is1L_looseVETO && njets >= 4 && nbjets == 4)   << 8;

   passSel_e = 0x0;
   passSel_m = 0x0;
   //for(int i=0;i<=3;i++)
   for(int i=0;i<=9;i++)
     {	
	passSel_e |= (CHECK_BIT(passSel_all,i) && is1E) << i;
	passSel_m |= (CHECK_BIT(passSel_all,i) && is1M) << i;
     }   
   
   bool chan_pass[chan_n];
   for(int c=0;c<chan_n;c++)
     {	
	chan_pass[c] = 0;
	if( (chan[c] == "all") ||
	    (chan[c] == "e" && is1E) ||
	    (chan[c] == "m" && is1M)
	  )
	  chan_pass[c] = 1;
     }   
   /*
   // Selection crtieria for plots
   //bool pass = (CHECK_BIT(passSel_all,0) || CHECK_BIT(passSel_all,1) || CHECK_BIT(passSel_all,2));
   bool pass = (CHECK_BIT(passSel_all,0) || CHECK_BIT(passSel_all,1) || CHECK_BIT(passSel_all,2) ||
                CHECK_BIT(passSel_all,3) || CHECK_BIT(passSel_all,4) || CHECK_BIT(passSel_all,5) ||
                CHECK_BIT(passSel_all,6) || CHECK_BIT(passSel_all,7) || CHECK_BIT(passSel_all,8));
   if( pass )
     {
	_v_Lepton->clear();

	for(int ie=0;ie<_v_ElectronTight->size();ie++) 
	  {
	     Lepton lep;
	     lep.setLepton(&_v_ElectronTight->at(ie),ie,1);
	     _v_Lepton->push_back(lep);
	  }	
	for(int im=0;im<_v_MuonTight->size();im++)
	  {
	     Lepton lep;
	     lep.setLepton(&_v_MuonTight->at(im),im,0);
	     _v_Lepton->push_back(lep);
	  }
	
	_trec->setElectron(_v_ElectronTight);
	_trec->setMuon(_v_MuonTight);
	_trec->setJet(_v_JetTight);
	_trec->setBJet(_v_BJetTight);
	_trec->setNonBJet(_v_NonBJetTight);
	_trec->setEvent(_v_Event);

	// top reconstruction
	_trec->run();	

    / *
	m_HiggsMass_TOPTOPLEPHBB = _trec->Higgs_TOPTOPLEPHBB_p4().M();
	m_HiggsMass_TOPHLEPBB = _trec->Higgs_TOPHLEPBB_p4().M();
	m_TopHadMass_TOPTOPLEPHAD = _trec->TopHad_TOPTOPLEPHAD_p4().M();
	//m_chi2_TOPTOPLEPHBB = _trec->chi2_TOPTOPLEPHBB();
	//m_chi2_TOPHLEPBB = _trec->chi2_TOPHLEPBB();
	//m_chi2_TOPTOPLEPHAD = _trec->chi2_TOPTOPLEPHAD();
	m_MVA_TOPTOPLEPHBB = _trec->MVA_TOPTOPLEPHBB();
	m_MVA_TOPHLEPBB = _trec->MVA_TOPHLEPBB();
	m_MVA_TOPTOPLEPHAD = _trec->MVA_TOPTOPLEPHAD();   
	m_LepCharge = _v_Lepton->at(0).charge();
	//m_HiggsEta_TOPTOPLEPHBB = _trec->Higgs_TOPTOPLEPHBB_p4().PseudoRapidity();
	m_HiggsEta_TOPHLEPBB = _trec->Higgs_TOPHLEPBB_p4().PseudoRapidity();
	m_TopLepMass_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().M();
	m_TopLepMass_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().M();
	m_TopLepMass_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().M();
	//m_TopLepPt_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().Pt();
	m_TopLepPt_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().Pt();
	//m_TopLepPt_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().Pt();
	//m_TopLepEta_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().PseudoRapidity();
	m_TopLepEta_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().PseudoRapidity();
	//m_TopLepEta_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().PseudoRapidity();
	m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = _trec->HiggsBJet1_TOPTOPLEPHBB_p4().DeltaR(_trec->HiggsBJet2_TOPTOPLEPHBB_p4());
	m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB = _trec->HiggsBJet1_TOPHLEPBB_p4().DeltaR(_trec->HiggsBJet2_TOPHLEPBB_p4());
	m_TopLepHiggsDr_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().DeltaR(_trec->Higgs_TOPTOPLEPHBB_p4());
	m_TopLepHiggsDr_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().DeltaR(_trec->Higgs_TOPHLEPBB_p4());
	m_TopLepTopHadDr_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().DeltaR(_trec->TopHad_TOPTOPLEPHAD_p4());
    * /

    m_HiggsMass_TOPHLEPBB = _trec->Higgs_TOPHLEPBB_p4().M();
    m_MVA_TOPHLEPBB = _trec->MVA_TOPHLEPBB();
    m_HiggsEta_TOPHLEPBB = _trec->Higgs_TOPHLEPBB_p4().PseudoRapidity();
    m_TopLepMass_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().M();
    m_TopLepPt_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().Pt();
    m_TopLepEta_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().PseudoRapidity();
    m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB = _trec->HiggsBJet1_TOPHLEPBB_p4().DeltaR(_trec->HiggsBJet2_TOPHLEPBB_p4());
    m_TopLepHiggsDr_TOPHLEPBB = _trec->TopLep_TOPHLEPBB_p4().DeltaR(_trec->Higgs_TOPHLEPBB_p4());
    //m_HiggsBJet1CSVv2_TOPHLEPBB = 2*(std::rand() / (double) INT_MAX)-1;
    //m_HiggsBJet2CSVv2_TOPHLEPBB = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopLepBJetCSVv2_TOPHLEPBB = 2*(std::rand() / (double) INT_MAX)-1;
    m_TopHadMass_TOPTOPLEPHAD = _trec->TopHad_TOPTOPLEPHAD_p4().M();
    m_MVA_TOPTOPLEPHAD = _trec->MVA_TOPTOPLEPHAD();
    m_TopLepMass_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().M();
    m_TopLepTopHadDr_TOPTOPLEPHAD = _trec->TopLep_TOPTOPLEPHAD_p4().DeltaR(_trec->TopHad_TOPTOPLEPHAD_p4());
    //m_TopLepBJetCSVv2_TOPTOPLEPHAD = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopHadBJetCSVv2_TOPTOPLEPHAD = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD = 2*(std::rand() / (double) INT_MAX)-1;
    m_HiggsMass_TOPTOPLEPHBB = _trec->Higgs_TOPTOPLEPHBB_p4().M();
    m_MVA_TOPTOPLEPHBB = _trec->MVA_TOPTOPLEPHBB();
    m_TopLepMass_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().M();
    m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = _trec->HiggsBJet1_TOPTOPLEPHBB_p4().DeltaR(_trec->HiggsBJet2_TOPTOPLEPHBB_p4());
    m_TopLepHiggsDr_TOPTOPLEPHBB = _trec->TopLep_TOPTOPLEPHBB_p4().DeltaR(_trec->Higgs_TOPTOPLEPHBB_p4());
    //m_HiggsBJet1CSVv2_TOPTOPLEPHBB = 2*(std::rand() / (double) INT_MAX)-1;
    //m_HiggsBJet2CSVv2_TOPTOPLEPHBB = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopLepBJetCSVv2_TOPTOPLEPHBB = 2*(std::rand() / (double) INT_MAX)-1;
    //m_TopHadNonBJetCSVv2_TOPTOPLEPHBB = 2*(std::rand() / (double) INT_MAX)-1;
    m_LepCharge = _v_Lepton->at(0).charge();


   for(int i=0;i<_v_Jet->size();i++)
     {
    if( fabs(_v_Jet->at(i).eta()) > 2.4 ) continue;

    if( _v_Jet->at(i).isTight() )
      {
         //std::cout << "jet compare: (" << _v_Jet->at(i).pt() << ", " << _v_Jet->at(i).eta()  << ", " << _v_Jet->at(i).phi() << ", " << _v_Jet->at(i).E() << ") =? ("
         //          << TopLepBJetCSVv2_TOPTOPLEPHAD_jet.Pt() << ", " << TopLepBJetCSVv2_TOPTOPLEPHAD_jet.Eta() << ", " << TopLepBJetCSVv2_TOPTOPLEPHAD_jet.Phi() << ", " << TopLepBJetCSVv2_TOPTOPLEPHAD_jet.E() << ") "; // << std::endl;
         double MAXDIF = 1e-2;

         if (
            ((_v_Jet->at(i).pt()-_trec->HiggsBJet1_TOPHLEPBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->HiggsBJet1_TOPHLEPBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->HiggsBJet1_TOPHLEPBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->HiggsBJet1_TOPHLEPBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_HiggsBJet1CSVv2_TOPHLEPBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->HiggsBJet2_TOPHLEPBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->HiggsBJet2_TOPHLEPBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->HiggsBJet2_TOPHLEPBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->HiggsBJet2_TOPHLEPBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_HiggsBJet2CSVv2_TOPHLEPBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopLepBJet_TOPHLEPBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopLepBJet_TOPHLEPBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopLepBJet_TOPHLEPBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopLepBJet_TOPHLEPBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopLepBJetCSVv2_TOPHLEPBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopLepBJet_TOPTOPLEPHAD_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopLepBJet_TOPTOPLEPHAD_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopLepBJet_TOPTOPLEPHAD_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopLepBJet_TOPTOPLEPHAD_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopLepBJetCSVv2_TOPTOPLEPHAD = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopHadBJet_TOPTOPLEPHAD_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopHadBJet_TOPTOPLEPHAD_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopHadBJet_TOPTOPLEPHAD_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopHadBJet_TOPTOPLEPHAD_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopHadBJetCSVv2_TOPTOPLEPHAD = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopHadWNonBJet1_TOPTOPLEPHAD_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopHadWNonBJet1_TOPTOPLEPHAD_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopHadWNonBJet1_TOPTOPLEPHAD_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopHadWNonBJet1_TOPTOPLEPHAD_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopHadWNonBJet2_TOPTOPLEPHAD_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopHadWNonBJet2_TOPTOPLEPHAD_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopHadWNonBJet2_TOPTOPLEPHAD_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopHadWNonBJet2_TOPTOPLEPHAD_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD = _v_Jet->at(i).CSVv2();
            
         if (
            ((_v_Jet->at(i).pt()-_trec->HiggsBJet1_TOPTOPLEPHBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->HiggsBJet1_TOPTOPLEPHBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->HiggsBJet1_TOPTOPLEPHBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->HiggsBJet1_TOPTOPLEPHBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_HiggsBJet1CSVv2_TOPTOPLEPHBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->HiggsBJet2_TOPTOPLEPHBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->HiggsBJet2_TOPTOPLEPHBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->HiggsBJet2_TOPTOPLEPHBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->HiggsBJet2_TOPTOPLEPHBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_HiggsBJet2CSVv2_TOPTOPLEPHBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopLepBJet_TOPTOPLEPHBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopLepBJet_TOPTOPLEPHBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopLepBJet_TOPTOPLEPHBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopLepBJet_TOPTOPLEPHBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopLepBJetCSVv2_TOPTOPLEPHBB = _v_Jet->at(i).CSVv2();

         if (
            ((_v_Jet->at(i).pt()-_trec->TopHadNonBJet_TOPTOPLEPHBB_p4().Pt())/_v_Jet->at(i).pt() < MAXDIF) &&
            ((_v_Jet->at(i).eta()-_trec->TopHadNonBJet_TOPTOPLEPHBB_p4().Eta())/_v_Jet->at(i).eta() < MAXDIF) &&
            ((_v_Jet->at(i).phi()-_trec->TopHadNonBJet_TOPTOPLEPHBB_p4().Phi())/_v_Jet->at(i).phi() < MAXDIF) &&
            ((_v_Jet->at(i).E()-_trec->TopHadNonBJet_TOPTOPLEPHBB_p4().E())/_v_Jet->at(i).E() < MAXDIF)
         )
            m_TopHadNonBJetCSVv2_TOPTOPLEPHBB = _v_Jet->at(i).CSVv2();

      }
     }
    //std::cout << std::endl;

    / *
	_mva->setVariable_HiggsMass_TOPHLEPBB(m_HiggsMass_TOPHLEPBB);
	_mva->setVariable_TopHadMass_TOPTOPLEPHAD(m_TopHadMass_TOPTOPLEPHAD);
	_mva->setVariable_MVA_TOPHLEPBB(m_MVA_TOPHLEPBB);
	_mva->setVariable_MVA_TOPTOPLEPHAD(m_MVA_TOPTOPLEPHAD);
	_mva->setVariable_LepCharge(m_LepCharge);
	_mva->setVariable_HiggsEta_TOPHLEPBB(m_HiggsEta_TOPHLEPBB);
	_mva->setVariable_TopLepMass_TOPHLEPBB(m_TopLepMass_TOPHLEPBB);
	_mva->setVariable_TopLepMass_TOPTOPLEPHAD(m_TopLepMass_TOPTOPLEPHAD);
	_mva->setVariable_TopLepPt_TOPHLEPBB(m_TopLepPt_TOPHLEPBB);
	_mva->setVariable_TopLepEta_TOPHLEPBB(m_TopLepEta_TOPHLEPBB);
	_mva->setVariable_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB(m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
	_mva->setVariable_TopLepHiggsDr_TOPHLEPBB(m_TopLepHiggsDr_TOPHLEPBB);
	_mva->setVariable_TopLepTopHadDr_TOPTOPLEPHAD(m_TopLepTopHadDr_TOPTOPLEPHAD);
	
	_mva->setVariable_HiggsMass_TOPTOPLEPHBB(m_HiggsMass_TOPTOPLEPHBB);
	_mva->setVariable_MVA_TOPTOPLEPHBB(m_MVA_TOPTOPLEPHBB);
	//_mva->setVariable_HiggsEta_TOPTOPLEPHBB(m_HiggsEta_TOPTOPLEPHBB);
	_mva->setVariable_TopLepMass_TOPTOPLEPHBB(m_TopLepMass_TOPTOPLEPHBB);
	//_mva->setVariable_TopLepPt_TOPTOPLEPHBB(m_TopLepPt_TOPTOPLEPHBB);
	//_mva->setVariable_TopLepEta_TOPTOPLEPHBB(m_TopLepEta_TOPTOPLEPHBB);
	_mva->setVariable_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB(m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
	_mva->setVariable_TopLepHiggsDr_TOPTOPLEPHBB(m_TopLepHiggsDr_TOPTOPLEPHBB);
    * /
    _mva->setVariable_HiggsMass_TOPHLEPBB(m_HiggsMass_TOPHLEPBB);
    _mva->setVariable_MVA_TOPHLEPBB(m_MVA_TOPHLEPBB);
    _mva->setVariable_HiggsEta_TOPHLEPBB(m_HiggsEta_TOPHLEPBB);
    _mva->setVariable_TopLepMass_TOPHLEPBB(m_TopLepMass_TOPHLEPBB);
    _mva->setVariable_TopLepPt_TOPHLEPBB(m_TopLepPt_TOPHLEPBB);
    _mva->setVariable_TopLepEta_TOPHLEPBB(m_TopLepEta_TOPHLEPBB);
    _mva->setVariable_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB(m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
    _mva->setVariable_TopLepHiggsDr_TOPHLEPBB(m_TopLepHiggsDr_TOPHLEPBB);
    _mva->setVariable_HiggsBJet1CSVv2_TOPHLEPBB(m_HiggsBJet1CSVv2_TOPHLEPBB);
    _mva->setVariable_HiggsBJet2CSVv2_TOPHLEPBB(m_HiggsBJet2CSVv2_TOPHLEPBB);
    _mva->setVariable_TopLepBJetCSVv2_TOPHLEPBB(m_TopLepBJetCSVv2_TOPHLEPBB);
    _mva->setVariable_TopHadMass_TOPTOPLEPHAD(m_TopHadMass_TOPTOPLEPHAD);
    _mva->setVariable_MVA_TOPTOPLEPHAD(m_MVA_TOPTOPLEPHAD);
    _mva->setVariable_TopLepMass_TOPTOPLEPHAD(m_TopLepMass_TOPTOPLEPHAD);
    _mva->setVariable_TopLepTopHadDr_TOPTOPLEPHAD(m_TopLepTopHadDr_TOPTOPLEPHAD);
    _mva->setVariable_TopLepBJetCSVv2_TOPTOPLEPHAD(m_TopLepBJetCSVv2_TOPTOPLEPHAD);
    _mva->setVariable_TopHadBJetCSVv2_TOPTOPLEPHAD(m_TopHadBJetCSVv2_TOPTOPLEPHAD);
    _mva->setVariable_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD(m_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
    _mva->setVariable_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD(m_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
    _mva->setVariable_HiggsMass_TOPTOPLEPHBB(m_HiggsMass_TOPTOPLEPHBB);
    _mva->setVariable_MVA_TOPTOPLEPHBB(m_MVA_TOPTOPLEPHBB);
    _mva->setVariable_TopLepMass_TOPTOPLEPHBB(m_TopLepMass_TOPTOPLEPHBB);
    _mva->setVariable_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB(m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
    _mva->setVariable_TopLepHiggsDr_TOPTOPLEPHBB(m_TopLepHiggsDr_TOPTOPLEPHBB);
    _mva->setVariable_HiggsBJet1CSVv2_TOPTOPLEPHBB(m_HiggsBJet1CSVv2_TOPTOPLEPHBB);
    _mva->setVariable_HiggsBJet2CSVv2_TOPTOPLEPHBB(m_HiggsBJet2CSVv2_TOPTOPLEPHBB);
    _mva->setVariable_TopLepBJetCSVv2_TOPTOPLEPHBB(m_TopLepBJetCSVv2_TOPTOPLEPHBB);
    _mva->setVariable_TopHadNonBJetCSVv2_TOPTOPLEPHBB(m_TopHadNonBJetCSVv2_TOPTOPLEPHBB);
    _mva->setVariable_LepCharge(m_LepCharge);
	
	// Evaluate MVA discriminator
	/ *
	_mvab3j4HutST = _mva->run("b3j4HutST");
	_mvab3j4HctST = _mva->run("b3j4HctST");
	_mvab3j4HutTT = _mva->run("b3j4HutTT");
	_mvab3j4HctTT = _mva->run("b3j4HctTT");

	_mvab3j3HutST = _mva->run("b3j3HutST");
	_mvab3j3HctST = _mva->run("b3j3HctST");
	_mvab3j3HutTT = _mva->run("b3j3HutTT");
	_mvab3j3HctTT = _mva->run("b3j3HctTT");

	_mvab2j4HutST = _mva->run("b2j4HutST");
	_mvab2j4HctST = _mva->run("b2j4HctST");
	_mvab2j4HutTT = _mva->run("b2j4HutTT");
	_mvab2j4HctTT = _mva->run("b2j4HctTT");
    * /
    _mvab0j3HutST = _mva->run("b0j3HutST");
    _mvab0j3HctST = _mva->run("b0j3HctST");
    _mvab0j3HutTT = _mva->run("b0j3HutTT");
    _mvab0j3HctTT = _mva->run("b0j3HctTT");
                           
    _mvab1j3HutST = _mva->run("b1j3HutST");
    _mvab1j3HctST = _mva->run("b1j3HctST");
    _mvab1j3HutTT = _mva->run("b1j3HutTT");
    _mvab1j3HctTT = _mva->run("b1j3HctTT");
                           
    _mvab2j3HutST = _mva->run("b2j3HutST");
    _mvab2j3HctST = _mva->run("b2j3HctST");
    _mvab2j3HutTT = _mva->run("b2j3HutTT");
    _mvab2j3HctTT = _mva->run("b2j3HctTT");
                           
    _mvab3j3HutST = _mva->run("b3j3HutST");
    _mvab3j3HctST = _mva->run("b3j3HctST");
    _mvab3j3HutTT = _mva->run("b3j3HutTT");
    _mvab3j3HctTT = _mva->run("b3j3HctTT");
                           
    _mvab0j4HutST = _mva->run("b0j4HutST");
    _mvab0j4HctST = _mva->run("b0j4HctST");
    _mvab0j4HutTT = _mva->run("b0j4HutTT");
    _mvab0j4HctTT = _mva->run("b0j4HctTT");
                           
    _mvab1j4HutST = _mva->run("b1j4HutST");
    _mvab1j4HctST = _mva->run("b1j4HctST");
    _mvab1j4HutTT = _mva->run("b1j4HutTT");
    _mvab1j4HctTT = _mva->run("b1j4HctTT");
                           
    _mvab2j4HutST = _mva->run("b2j4HutST");
    _mvab2j4HctST = _mva->run("b2j4HctST");
    _mvab2j4HutTT = _mva->run("b2j4HutTT");
    _mvab2j4HctTT = _mva->run("b2j4HctTT");
                           
    _mvab3j4HutST = _mva->run("b3j4HutST");
    _mvab3j4HctST = _mva->run("b3j4HctST");
    _mvab3j4HutTT = _mva->run("b3j4HutTT");
    _mvab3j4HctTT = _mva->run("b3j4HctTT");
                           
    _mvab4j4HutST = _mva->run("b4j4HutST");
    _mvab4j4HctST = _mva->run("b4j4HctST");
    _mvab4j4HutTT = _mva->run("b4j4HutTT");
    _mvab4j4HctTT = _mva->run("b4j4HctTT");

    / *
    if (nbjets == 0 && njets == 3) std::cout << "b0j3" << std::endl;
    else if (nbjets == 1 && njets == 3) std::cout << "b1j3" << std::endl;
    else if (nbjets == 2 && njets == 3) std::cout << "b2j3" << std::endl;
    else if (nbjets == 3 && njets == 3) std::cout << "b3j3" << std::endl;
    else if (nbjets == 0 && njets == 4) std::cout << "b0j4" << std::endl;
    else if (nbjets == 1 && njets == 4) std::cout << "b1j4" << std::endl;
    else if (nbjets == 2 && njets == 4) std::cout << "b2j4" << std::endl;
    else if (nbjets == 3 && njets == 4) std::cout << "b3j4" << std::endl;
    else if (nbjets == 4 && njets == 4) std::cout << "b4j4" << std::endl;
    * /
	*/
	bool sel_pass[sel_n];
	for(int c=0;c<sel_n;c++)
	  {	
	     sel_pass[c] = 0;
	     if( 
        /*
		(sel[c] == "b2j4" && nbjets == 2 && njets == 4) ||
		(sel[c] == "b3j3" && nbjets == 3 && njets == 3) ||
		(sel[c] == "b3j4" && nbjets == 3 && njets >= 4)
        */
        (sel[c] == "b0j3" && nbjets == 0 && njets == 3) ||
        (sel[c] == "b1j3" && nbjets == 1 && njets == 3) ||
        (sel[c] == "b2j3" && nbjets == 2 && njets == 3) ||
        (sel[c] == "b3j3" && nbjets == 3 && njets == 3) ||
        (sel[c] == "b0j4" && nbjets == 0 && njets == 4) ||
        (sel[c] == "b1j4" && nbjets == 1 && njets == 4) ||
        (sel[c] == "b2j4" && nbjets == 2 && njets == 4) ||
        (sel[c] == "b3j4" && nbjets == 3 && njets == 4) ||
        (sel[c] == "b4j4" && nbjets == 4 && njets == 4)
	       )
	       sel_pass[c] = 1;
	  }   
	
	for(int isys=0;isys<sys_n;isys++)
	  {	
         /*
	     double sfj = w;
	     
	     // 1d
	     std::vector<std::string> histNAMESSEL;
	     std::vector<int> histSYS;
	     std::vector<int> histVAR;

	     for(int ih=0;ih<histname_n;ih++)
	       {
		  for(int ic=0;ic<chan_n;ic++)
		    {
		       if( !chan_pass[ic] ) continue;
		       
		       for(int is=0;is<sel_n;is++)
			 {
			    if( !sel_pass[is] ) continue;

			    if( ih == 0 && isys == 0 )
			      {
				 fillTree(ic,is);
			      }				 
				 
			    // [CHAN][TYPE][SEL][VAR][2*(NSYS-1)+1]
			    histNAMESSEL.push_back(histNAMES[ic][0][is][ih][isys]);
			    histSYS.push_back(isys);
			    histVAR.push_back(ih);
			 }		       
		    }		       
	       }
         */
          for(int ic=0;ic<chan_n;ic++)
            {
               if( !chan_pass[ic] ) continue;

               for(int is=0;is<sel_n;is++)
             {
                if( !sel_pass[is] ) continue;

                 fillTree(ic,is);
             }
            }
         /*
	     int nHISTSEL = histNAMESSEL.size();
	     for(int ih=0;ih<nHISTSEL;ih++)
	       {		       
		  TH1D *h = _m1d_Hist->find(histNAMESSEL.at(ih))->second;
		  
		  int hidx = ih;
		  std::string varName = histname[histVAR[hidx]];
		  fillHisto1D(h,sfj,sys[histSYS[hidx]],0,varName);
	       }
         */
	  }
    /*
     }
   
   fillPassSel(_h_PassSel_all,_h_PassSel_e,_h_PassSel_m,w);
   */
}

void Hist::close()
{
   _fout->Write();
   _fout->Close();
//   _fevc.close();
   _fevcVal.close();
   
//   delete rnd;
}

void Hist::fillTree(int ic,int is)
{
   _trout[ic][is]->Fill();
}

bool Hist::printout(bool doPrint)
{
//	if( CHECK_BIT(passSel,1)
//	  )
//	  {	     
/*	     fprintf(_fevc,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %+2d  %6.2f %+4.2f %+4.2f    %6.1f  %+4.2f    %d \n",
		     run, lumi, id,
		     l1id, l1pt, l1eta, l1phi,
		     l2id, l2pt, l2eta, l2phi,
		     metpt, metphi,
		     njets);*/
	     
//	     _trout->Fill();
//	  }
//     }   
   
   return 1;
}

void Hist::fillHisto1D(TH1D *h,float sfj,std::string sys,int ilep,std::string varName)
{
/*   if( strcmp(varName.c_str(),"h_H_m_") == 0 )
     {	
	float H_m = (_resTop) ? _H_p4.M() : -666.;
	h->Fill(H_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_H_pt_") == 0 )
     {	
	float H_pt = (_resTop) ? _H_p4.Pt() : -666.;
	h->Fill(H_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_H_eta_") == 0 )
     {	
	float H_eta = (_resTop) ? _H_p4.PseudoRapidity() : -666.;
	h->Fill(H_eta,sfj);
     }
   else if( strcmp(varName.c_str(),"h_top_m_") == 0 )
     {	
	float top_m = (_resTop) ? _top_p4.M() : -666.;
	h->Fill(top_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_top_pt_") == 0 )
     {	
	float top_pt = (_resTop) ? _top_p4.Pt() : -666.;
	h->Fill(top_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_top_eta_") == 0 )
     {	
	float top_eta = (_resTop) ? _top_p4.PseudoRapidity() : -666.;
	h->Fill(top_eta,sfj);
     }
   else if( strcmp(varName.c_str(),"h_W_m_") == 0 )
     {	
	float W_m = (_resTop) ? _W_p4.M() : -666.;
	h->Fill(W_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_W_pt_") == 0 )
     {	
	float W_pt = (_resTop) ? _W_p4.Pt() : -666.;
	h->Fill(W_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_W_eta_") == 0 )
     {	
	float W_eta = (_resTop) ? _W_p4.PseudoRapidity() : -666.;
	h->Fill(W_eta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_l_pt_") == 0 )
     {	
	float l_pt = (_resTop) ? _l_p4.Pt() : -666.;
	h->Fill(l_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_l_eta_") == 0 )
     {	
	float l_eta = (_resTop) ? _l_p4.PseudoRapidity() : -666.;
	h->Fill(l_eta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_nu_pt_") == 0 )
     {	
	float nu_pt = (_resTop) ? _nu_p4.Pt() : -666.;
	h->Fill(nu_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_nu_eta_") == 0 )
     {	
	float nu_eta = (_resTop) ? _nu_p4.PseudoRapidity() : -666.;
	h->Fill(nu_eta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_Hb_pt_") == 0 )
     {	
	float Hb1_pt = (_resTop) ? _Hb1_p4.Pt() : -666.;
	h->Fill(Hb1_pt,sfj);
	float Hb2_pt = (_resTop) ? _Hb2_p4.Pt() : -666.;
	h->Fill(Hb2_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_Hb_eta_") == 0 )
     {	
	float Hb1_eta = (_resTop) ? _Hb1_p4.PseudoRapidity() : -666.;
	h->Fill(Hb1_eta,sfj);
	float Hb2_eta = (_resTop) ? _Hb2_p4.PseudoRapidity() : -666.;
	h->Fill(Hb2_eta,sfj);
     }
   else if( strcmp(varName.c_str(),"h_topb_pt_") == 0 )
     {	
	float topb_pt = (_resTop) ? _topb_p4.Pt() : -666.;
	h->Fill(topb_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_topb_eta_") == 0 )
     {	
	float topb_eta = (_resTop) ? _topb_p4.PseudoRapidity() : -666.;
	h->Fill(topb_eta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_Hb1_Hb2_dr_") == 0 )
     {	
	float Hb1_Hb2_dr = (_resTop) ? _Hb1_p4.DeltaR(_Hb2_p4) : -666.;
	h->Fill(Hb1_Hb2_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_H_top_dr_") == 0 )
     {	
	float H_top_dr = (_resTop) ? _H_p4.DeltaR(_top_p4) : -666.;
	h->Fill(H_top_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_l_nu_dr_") == 0 )
     {	
	float l_nu_dr = (_resTop) ? _l_p4.DeltaR(_nu_p4) : -666.;
	h->Fill(l_nu_dr,sfj);
     }
   else if( strcmp(varName.c_str(),"h_W_topb_dr_") == 0 )
     {	
	float W_topb_dr = (_resTop) ? _W_p4.DeltaR(_topb_p4) : -666.;
	h->Fill(W_topb_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_Hb1_Hb2_costheta_") == 0 )
     {	
	float Hb1_Hb2_costheta = (_resTop) ? cos(_Hb1_p4.Angle(_Hb2_p4.Vect())) : -666.;
	h->Fill(Hb1_Hb2_costheta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_H_top_costheta_") == 0 )
     {	
	float H_top_costheta = (_resTop) ? cos(_H_p4.Angle(_top_p4.Vect())) : -666.;
	h->Fill(H_top_costheta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_l_nu_costheta_") == 0 )
     {	
	float l_nu_costheta = (_resTop) ? cos(_l_p4.Angle(_nu_p4.Vect())) : -666.;
	h->Fill(l_nu_costheta,sfj);
     }
   else if( strcmp(varName.c_str(),"h_W_topb_costheta_") == 0 )
     {	
	float W_topb_costheta = (_resTop) ? cos(_W_p4.Angle(_topb_p4.Vect())) : -666.;
	h->Fill(W_topb_costheta,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_HT_") == 0 )
     {
	float HT = 0.;
	int njets = _v_JetTight->size();
	for(int ij=0;ij<njets;ij++) 
	  HT += _v_JetTight->at(ij).p4().Pt();
	h->Fill(HT,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_MET_") == 0 )
     {
	float metpt = _v_Event->at(0).metpt();
	h->Fill(metpt,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_njet_") == 0 )
     {
	int njets = _v_JetTight->size();
	h->Fill(njets,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_nu_phi_") == 0 )
     {	
	float nu_phi = (_resTop) ? _nu_p4.Phi() : -666.;
	h->Fill(nu_phi,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_MET_phi_") == 0 )
     {
	float metphi = _v_Event->at(0).metphi();
	h->Fill(metphi,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_H_nu_dr_") == 0 )
     {	
	float H_nu_dr = (_resTop) ? _H_p4.DeltaR(_nu_p4) : -666.;
	h->Fill(H_nu_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_H_l_dr_") == 0 )
     {	
	float H_l_dr = (_resTop) ? _H_p4.DeltaR(_l_p4) : -666.;
	h->Fill(H_l_dr,sfj);
     }*/
   if( strcmp(varName.c_str(),"h_chi2_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->chi2_TOPTOPLEPHBB(),sfj);
   else if( strcmp(varName.c_str(),"h_chi2_TOPHLEPBB_") == 0 ) h->Fill(_trec->chi2_TOPHLEPBB(),sfj);
   else if( strcmp(varName.c_str(),"h_chi2_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->chi2_TOPTOPLEPHAD(),sfj);
   else if( strcmp(varName.c_str(),"h_MVA_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->MVA_TOPTOPLEPHBB(),sfj);
   else if( strcmp(varName.c_str(),"h_MVA_TOPHLEPBB_") == 0 ) h->Fill(_trec->MVA_TOPHLEPBB(),sfj);
   else if( strcmp(varName.c_str(),"h_MVA_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->MVA_TOPTOPLEPHAD(),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsMass_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->Higgs_TOPTOPLEPHBB_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsMass_TOPHLEPBB_") == 0 ) h->Fill(_trec->Higgs_TOPHLEPBB_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_TopHadMass_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->TopHad_TOPTOPLEPHAD_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_LepCharge_") == 0 ) h->Fill(_v_Lepton->at(0).charge(),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsEta_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->Higgs_TOPTOPLEPHBB_p4().PseudoRapidity(),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsEta_TOPHLEPBB_") == 0 ) h->Fill(_trec->Higgs_TOPHLEPBB_p4().PseudoRapidity(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepMass_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHBB_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepMass_TOPHLEPBB_") == 0 ) h->Fill(_trec->TopLep_TOPHLEPBB_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepMass_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHAD_p4().M(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepPt_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHBB_p4().Pt(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepPt_TOPHLEPBB_") == 0 ) h->Fill(_trec->TopLep_TOPHLEPBB_p4().Pt(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepPt_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHAD_p4().Pt(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepEta_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHBB_p4().PseudoRapidity(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepEta_TOPHLEPBB_") == 0 ) h->Fill(_trec->TopLep_TOPHLEPBB_p4().PseudoRapidity(),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepEta_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHAD_p4().PseudoRapidity(),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->HiggsBJet1_TOPTOPLEPHBB_p4().DeltaR(_trec->HiggsBJet2_TOPTOPLEPHBB_p4()),sfj);
   else if( strcmp(varName.c_str(),"h_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_") == 0 ) h->Fill(_trec->HiggsBJet1_TOPHLEPBB_p4().DeltaR(_trec->HiggsBJet2_TOPHLEPBB_p4()),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepHiggsDr_TOPTOPLEPHBB_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHBB_p4().DeltaR(_trec->Higgs_TOPTOPLEPHBB_p4()),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepHiggsDr_TOPHLEPBB_") == 0 ) h->Fill(_trec->TopLep_TOPHLEPBB_p4().DeltaR(_trec->Higgs_TOPHLEPBB_p4()),sfj);
   else if( strcmp(varName.c_str(),"h_TopLepTopHadDr_TOPTOPLEPHAD_") == 0 ) h->Fill(_trec->TopLep_TOPTOPLEPHAD_p4().DeltaR(_trec->TopHad_TOPTOPLEPHAD_p4()),sfj);
   
/*   else if( strcmp(varName.c_str(),"h_l_charge_") == 0 )
     {	
	h->Fill(_v_Lepton->at(0).charge(),sfj);
     }

   else if( strcmp(varName.c_str(),"h_TTbar_topLep_m_") == 0 )
     {	
	float TTbar_top_m = (_resTop) ? _topLep_TTbar_p4.M() : -666.;
	h->Fill(TTbar_top_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_topHad_m_") == 0 )
     {	
	float TTbar_top_m = (_resTop) ? _topLep_TTbar_p4.M() : -666.;
	h->Fill(TTbar_top_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_topLep_pt_") == 0 )
     {	
	float TTbar_top_pt = (_resTop) ? _topLep_TTbar_p4.Pt() : -666.;
	h->Fill(TTbar_top_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_topHad_pt_") == 0 )
     {	
	float TTbar_top_pt = (_resTop) ? _topLep_TTbar_p4.Pt() : -666.;
	h->Fill(TTbar_top_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_topLep_eta_") == 0 )
     {	
	float TTbar_top_eta = (_resTop) ? _topLep_TTbar_p4.PseudoRapidity() : -666.;
	h->Fill(TTbar_top_eta,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_topHad_eta_") == 0 )
     {	
	float TTbar_top_eta = (_resTop) ? _topLep_TTbar_p4.PseudoRapidity() : -666.;
	h->Fill(TTbar_top_eta,sfj);
     }	     

   else if( strcmp(varName.c_str(),"h_TTbar_WLep_m_") == 0 )
     {	
	float TTbar_W_m = (_resTop) ? _WLep_TTbar_p4.M() : -666.;
	h->Fill(TTbar_W_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_WHad_m_") == 0 )
     {	
	float TTbar_W_m = (_resTop) ? _WLep_TTbar_p4.M() : -666.;
	h->Fill(TTbar_W_m,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_WLep_pt_") == 0 )
     {	
	float TTbar_W_pt = (_resTop) ? _WLep_TTbar_p4.Pt() : -666.;
	h->Fill(TTbar_W_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_WHad_pt_") == 0 )
     {	
	float TTbar_W_pt = (_resTop) ? _WLep_TTbar_p4.Pt() : -666.;
	h->Fill(TTbar_W_pt,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_WLep_eta_") == 0 )
     {	
	float TTbar_W_eta = (_resTop) ? _WLep_TTbar_p4.PseudoRapidity() : -666.;
	h->Fill(TTbar_W_eta,sfj);
     }	     
   else if( strcmp(varName.c_str(),"h_TTbar_WHad_eta_") == 0 )
     {	
	float TTbar_W_eta = (_resTop) ? _WLep_TTbar_p4.PseudoRapidity() : -666.;
	h->Fill(TTbar_W_eta,sfj);
     }	     

   else if( strcmp(varName.c_str(),"h_TTbar_tbLep_tWLep_Dr_") == 0 )
     {	
	float W_topb_dr = (_resTop) ? _WLep_TTbar_p4.DeltaR(_topbLep_TTbar_p4) : -666.;
	h->Fill(W_topb_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_TTbar_tbHad_tWHad_Dr_") == 0 )
     {	
	float W_topb_dr = (_resTop) ? _WHad_TTbar_p4.DeltaR(_topbHad_TTbar_p4) : -666.;
	h->Fill(W_topb_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_TTbar_tWj1_tWj2_Dr_") == 0 )
     {	
	float topWj1_topWj2_dr = (_resTop) ? _topWj1_TTbar_p4.DeltaR(_topWj2_TTbar_p4) : -666.;
	h->Fill(topWj1_topWj2_dr,sfj);
     }   
   else if( strcmp(varName.c_str(),"h_TTbar_chi2_") == 0 )
     {	
	h->Fill(_trec->chi2TTbar(),sfj);
     }*/
   /*
   else if( strcmp(varName.c_str(),"h_MVAb3j4HutST_") == 0 ) h->Fill(_mvab3j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HctST_") == 0 ) h->Fill(_mvab3j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HutTT_") == 0 ) h->Fill(_mvab3j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HctTT_") == 0 ) h->Fill(_mvab3j4HctTT,sfj);
   
   else if( strcmp(varName.c_str(),"h_MVAb3j3HutST_") == 0 ) h->Fill(_mvab3j3HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HctST_") == 0 ) h->Fill(_mvab3j3HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HutTT_") == 0 ) h->Fill(_mvab3j3HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HctTT_") == 0 ) h->Fill(_mvab3j3HctTT,sfj);

   else if( strcmp(varName.c_str(),"h_MVAb2j4HutST_") == 0 ) h->Fill(_mvab2j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HctST_") == 0 ) h->Fill(_mvab2j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HutTT_") == 0 ) h->Fill(_mvab2j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HctTT_") == 0 ) h->Fill(_mvab2j4HctTT,sfj);
   */
   else if( strcmp(varName.c_str(),"h_MVAb0j3HutST_") == 0 ) h->Fill(_mvab0j3HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j3HctST_") == 0 ) h->Fill(_mvab0j3HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j3HutTT_") == 0 ) h->Fill(_mvab0j3HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j3HctTT_") == 0 ) h->Fill(_mvab0j3HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb1j3HutST_") == 0 ) h->Fill(_mvab1j3HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j3HctST_") == 0 ) h->Fill(_mvab1j3HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j3HutTT_") == 0 ) h->Fill(_mvab1j3HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j3HctTT_") == 0 ) h->Fill(_mvab1j3HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb2j3HutST_") == 0 ) h->Fill(_mvab2j3HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j3HctST_") == 0 ) h->Fill(_mvab2j3HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j3HutTT_") == 0 ) h->Fill(_mvab2j3HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j3HctTT_") == 0 ) h->Fill(_mvab2j3HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb3j3HutST_") == 0 ) h->Fill(_mvab3j3HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HctST_") == 0 ) h->Fill(_mvab3j3HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HutTT_") == 0 ) h->Fill(_mvab3j3HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j3HctTT_") == 0 ) h->Fill(_mvab3j3HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb0j4HutST_") == 0 ) h->Fill(_mvab0j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j4HctST_") == 0 ) h->Fill(_mvab0j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j4HutTT_") == 0 ) h->Fill(_mvab0j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb0j4HctTT_") == 0 ) h->Fill(_mvab0j4HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb1j4HutST_") == 0 ) h->Fill(_mvab1j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j4HctST_") == 0 ) h->Fill(_mvab1j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j4HutTT_") == 0 ) h->Fill(_mvab1j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb1j4HctTT_") == 0 ) h->Fill(_mvab1j4HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb2j4HutST_") == 0 ) h->Fill(_mvab2j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HctST_") == 0 ) h->Fill(_mvab2j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HutTT_") == 0 ) h->Fill(_mvab2j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb2j4HctTT_") == 0 ) h->Fill(_mvab2j4HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb3j4HutST_") == 0 ) h->Fill(_mvab3j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HctST_") == 0 ) h->Fill(_mvab3j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HutTT_") == 0 ) h->Fill(_mvab3j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb3j4HctTT_") == 0 ) h->Fill(_mvab3j4HctTT,sfj);
                                                                         
   else if( strcmp(varName.c_str(),"h_MVAb4j4HutST_") == 0 ) h->Fill(_mvab4j4HutST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb4j4HctST_") == 0 ) h->Fill(_mvab4j4HctST,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb4j4HutTT_") == 0 ) h->Fill(_mvab4j4HutTT,sfj);
   else if( strcmp(varName.c_str(),"h_MVAb4j4HctTT_") == 0 ) h->Fill(_mvab4j4HctTT,sfj);
}

void Hist::fillPassSel(TH1D *h,TH1D *he,TH1D *hm,float w)
{
   int nSel = h->GetXaxis()->GetNbins();
   for(int i=0;i<nSel;i++)
     {	
	bool pass = CHECK_BIT(passSel_all,i);
	if( pass ) h->Fill(i,w);

	bool pass_e = CHECK_BIT(passSel_e,i);
	if( pass_e ) he->Fill(i,w);

	bool pass_m = CHECK_BIT(passSel_m,i);
	if( pass_m ) hm->Fill(i,w);
     }   
}

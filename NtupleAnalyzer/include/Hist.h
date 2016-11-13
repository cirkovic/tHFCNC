#ifndef HIST_H
#define HIST_H

#include "Analyzer.h"

#include "TopReco.h"
#include "ApplyMVA.h"
#include "Helper.h"

//#include "FakeWeight.h"

class Hist
{
   
 public:
   
   Hist(std::string home);
   virtual ~Hist();
   
   void setElectron(std::vector<Electron> *v)                  {_v_Electron = v;};
   void setMuon(std::vector<Muon> *v)                      {_v_Muon = v;};
   void setJet(std::vector<Jet> *v)                      {_v_Jet = v;};
   void setEvent(std::vector<Event> *v)                      {_v_Event = v;};
   void setTruth(std::vector<Truth> *v)                      {_v_Truth = v;};

//   void setFakeWeight(SKYPLOT::FakeWeight fakeWeight)      {_fakeWeight = fakeWeight;}; 
	
   void init();
   void fill();
   void close();
   bool printout(bool doPrint);
   void fillTree(int ic,int is);

   std::vector<Lepton> filterPt(double Pt1, double Pt2, double Pt3,
				std::vector<Electron>* ntElectron,
				std::vector<Muon>* ntMuon);
   
   unsigned int passSel_all;
   unsigned int passSel_e;
   unsigned int passSel_m;
   
   void fillHisto1D(TH1D *h,float sfj,std::string sys,int ilep,std::string varName);

   void fillPassSel(TH1D *h,TH1D *he,TH1D *hm,float sfj);
   
//   float getVar(std::string sys,int ijet,std::string varName,int ibin);
//   std::vector<float> getVarVec(std::string sys,int ijet,std::string varName,int ibin);
//   std::pair<float,float> getVar2d(std::string sys,int ijet,std::string varName,int ibin);
//   std::vector<float> getVar3d(std::string sys,int ijet,std::string varName,int ibin);
//   float getPt(std::string sys);
//   double DeltaR(double eta1,double phi1,double eta2,double phi2);

//   double getSF(std::vector<std::pair<double,double> > vBin,
//		std::vector<std::pair<double,double> > vSf,
//		double var);
//   double get2DSF(std::vector<std::pair<double,double> > vBinX,
//		  std::vector<std::pair<double,double> > vBinY,
//		  std::vector<std::pair<double,double> > vSf,
//		  double varX,
//		  double varY);
	
   std::string chan[100];
//   std::string char[100];                                                                                           
   std::string sys[100];                                                                                               
   std::string sys_low[100];                                                                                           
   std::string sys_up[100];                                                                                            
   std::string histname[100];
   std::string type[100];                                                                                              
   std::string sel[100];                                                                                               
   
   std::vector<std::string> hname;
	
   // var / ptbin / bins / x1,x2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwBin;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwSf;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwBin_btag;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwSf_btag;
   
   // var[0] / ptbin / bins / v1,v2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinX;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinY;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinZ;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DSf;
   
   // var[0] / ptbin / bins / v1,v2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DBinX;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DBinY;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DSf;

   int sel_n;
   int chan_n;
   int histname_n;
   int type_n;
   int sys_n;                                                                                                          
   int sys_low_n;
   int sys_up_n;
   int n_jets;
	
   bool fillThis;

   // Tree variables

   double m_weight;
   /*
   double m_HiggsMass_TOPTOPLEPHBB;
   double m_HiggsMass_TOPHLEPBB;
   double m_TopHadMass_TOPTOPLEPHAD;
   double m_chi2_TOPTOPLEPHBB;
   double m_chi2_TOPHLEPBB;
   double m_chi2_TOPTOPLEPHAD;
   double m_MVA_TOPTOPLEPHBB;
   double m_MVA_TOPHLEPBB;
   double m_MVA_TOPTOPLEPHAD;
   int m_LepCharge;
   double m_HiggsEta_TOPTOPLEPHBB;
   double m_HiggsEta_TOPHLEPBB;
   double m_TopLepMass_TOPTOPLEPHBB;
   double m_TopLepMass_TOPHLEPBB;
   double m_TopLepMass_TOPTOPLEPHAD;
   double m_TopLepPt_TOPTOPLEPHBB;
   double m_TopLepPt_TOPHLEPBB;
   double m_TopLepPt_TOPTOPLEPHAD;
   double m_TopLepEta_TOPTOPLEPHBB;
   double m_TopLepEta_TOPHLEPBB;
   double m_TopLepEta_TOPTOPLEPHAD;
   double m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
   double m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
   double m_TopLepHiggsDr_TOPTOPLEPHBB;
   double m_TopLepHiggsDr_TOPHLEPBB;
   double m_TopLepTopHadDr_TOPTOPLEPHAD;
   */

   double m_HiggsMass_TOPHLEPBB;
   double m_MVA_TOPHLEPBB;
   double m_HiggsEta_TOPHLEPBB;
   double m_TopLepMass_TOPHLEPBB;
   double m_TopLepPt_TOPHLEPBB;
   double m_TopLepEta_TOPHLEPBB;
   double m_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
   double m_TopLepHiggsDr_TOPHLEPBB;
   double m_HiggsBJet1CSVv2_TOPHLEPBB;
   double m_HiggsBJet2CSVv2_TOPHLEPBB;
   double m_TopLepBJetCSVv2_TOPHLEPBB;
   double m_TopHadMass_TOPTOPLEPHAD;
   double m_MVA_TOPTOPLEPHAD;
   double m_TopLepMass_TOPTOPLEPHAD;
   double m_TopLepTopHadDr_TOPTOPLEPHAD;
   double m_TopLepBJetCSVv2_TOPTOPLEPHAD;
   double m_TopHadBJetCSVv2_TOPTOPLEPHAD;
   double m_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
   double m_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
   double m_HiggsMass_TOPTOPLEPHBB;
   double m_MVA_TOPTOPLEPHBB;
   double m_TopLepMass_TOPTOPLEPHBB;
   double m_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
   double m_TopLepHiggsDr_TOPTOPLEPHBB;
   double m_HiggsBJet1CSVv2_TOPTOPLEPHBB;
   double m_HiggsBJet2CSVv2_TOPTOPLEPHBB;
   double m_TopLepBJetCSVv2_TOPTOPLEPHBB;
   double m_TopHadNonBJetCSVv2_TOPTOPLEPHBB;
   int    m_LepCharge;

   int    njets;
   int    nbjets;

 protected:

   std::vector<Lepton> *_v_Lepton;
   
   std::vector<Electron>             *_v_Electron;
   std::vector<Muon>             *_v_Muon;
   std::vector<Electron>             *_v_ElectronTight;
   std::vector<Electron>             *_v_ElectronLoose;
   std::vector<Muon>             *_v_MuonTight;
   std::vector<Muon>             *_v_MuonLoose;
   
   std::vector<Event>             *_v_Event;
   std::vector<Jet>             *_v_Jet;
   std::vector<Jet>             *_v_JetTight;
   std::vector<Jet>             *_v_BJetTight;
   std::vector<Jet>             *_v_NonBJetTight;
   
   std::vector<Truth>             *_v_Truth;
   
//   SKYPLOT::FakeWeight           _fakeWeight;
  
   std::map<std::string, TH1D*> *_m1d_Hist;
   
   TH1D *_h_PassSel_all;
   TH1D *_h_PassSel_e;
   TH1D *_h_PassSel_m;

   std::vector<std::pair<std::vector<std::string>,double*> > *_s_Hist;
	
 private:

   // [CHAN][TYPE][SEL][VAR][2*(NSYS-1)+1]
   std::string histNAMES[5][5][5][100][2];
   
//   JetCorrectionUncertainty *jesTotal;
   
   double cJER[5];
   double cJER_down[5];
   double cJER_up[5];
	
//	TLorentzVector *v_jet;
//	TLorentzVector *v_jet_sys_jesTotalUp;
//	TLorentzVector *v_jet_sys_jesTotalLow;
//	TLorentzVector *v_jet_sys_jerTotalUp;
//	TLorentzVector *v_jet_sys_jerTotalLow;
//	TLorentzVector *v_mu;
   
   FILE *_fevc;
   std::ofstream _fevcVal;
   
   std::string _home;
   
   TRandom3 *rnd;
   
   TFile *_fout;
   TTree *_trout[10][10];
   
   TopReco *_trec;
   ApplyMVA *_mva;
   
   Helper *help;
   
   /*
   double _mvab3j4HutST;
   double _mvab3j4HctST;
   double _mvab3j4HutTT;
   double _mvab3j4HctTT;

   double _mvab3j3HutST;
   double _mvab3j3HctST;
   double _mvab3j3HutTT;
   double _mvab3j3HctTT;

   double _mvab2j4HutST;
   double _mvab2j4HctST;
   double _mvab2j4HutTT;
   double _mvab2j4HctTT;
   */
   double _mvab0j3HutST;
   double _mvab0j3HctST;
   double _mvab0j3HutTT;
   double _mvab0j3HctTT;
              
   double _mvab1j3HutST;
   double _mvab1j3HctST;
   double _mvab1j3HutTT;
   double _mvab1j3HctTT;
              
   double _mvab2j3HutST;
   double _mvab2j3HctST;
   double _mvab2j3HutTT;
   double _mvab2j3HctTT;
              
   double _mvab3j3HutST;
   double _mvab3j3HctST;
   double _mvab3j3HutTT;
   double _mvab3j3HctTT;
              
   double _mvab0j4HutST;
   double _mvab0j4HctST;
   double _mvab0j4HutTT;
   double _mvab0j4HctTT;
              
   double _mvab1j4HutST;
   double _mvab1j4HctST;
   double _mvab1j4HutTT;
   double _mvab1j4HctTT;
              
   double _mvab2j4HutST;
   double _mvab2j4HctST;
   double _mvab2j4HutTT;
   double _mvab2j4HctTT;
              
   double _mvab3j4HutST;
   double _mvab3j4HctST;
   double _mvab3j4HutTT;
   double _mvab3j4HctTT;
              
   double _mvab4j4HutST;
   double _mvab4j4HctST;
   double _mvab4j4HutTT;
   double _mvab4j4HctTT;
};

#endif

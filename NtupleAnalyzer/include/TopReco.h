#ifndef TOPRECO_H
#define TOPRECO_H

#include "Analyzer.h"

#include "Helper.h"

#include "kinfit.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class TopReco
{
   
 public:

   TopReco();
   virtual ~TopReco();
   
   void setElectron(std::vector<Electron> *v)                {_v_Electron = v;};
   void setMuon(std::vector<Muon> *v)                        {_v_Muon = v;};
   void setJet(std::vector<Jet> *v)                          {_v_Jet = v;};
   void setBJet(std::vector<Jet> *v)                         {_v_BJet = v;};
   void setNonBJet(std::vector<Jet> *v)                      {_v_NonBJet = v;};
   void setEvent(std::vector<Event> *v)                      {_v_Event = v;};
   void setLepton(std::vector<Lepton> *v)                    {_v_Lepton = v;};
	
   void init();
   void run();
   void close();

   TLorentzVector TopLepBJet_TOPTOPLEPHBB_p4()      {return _TopLepBJet_TOPTOPLEPHBB_p4;};
   TLorentzVector HiggsBJet1_TOPTOPLEPHBB_p4()      {return _HiggsBJet1_TOPTOPLEPHBB_p4;};
   TLorentzVector HiggsBJet2_TOPTOPLEPHBB_p4()      {return _HiggsBJet2_TOPTOPLEPHBB_p4;};
   TLorentzVector TopLepWLep_TOPTOPLEPHBB_p4()      {return _TopLepWLep_TOPTOPLEPHBB_p4;};
   TLorentzVector TopLepWNu_TOPTOPLEPHBB_p4()       {return _TopLepWNu_TOPTOPLEPHBB_p4;};
   TLorentzVector TopLep_TOPTOPLEPHBB_p4()          {return _TopLep_TOPTOPLEPHBB_p4;};
   TLorentzVector Higgs_TOPTOPLEPHBB_p4()           {return _Higgs_TOPTOPLEPHBB_p4;};
   TLorentzVector TopLepW_TOPTOPLEPHBB_p4()         {return _TopLepW_TOPTOPLEPHBB_p4;};
   TLorentzVector TopHadNonBJet_TOPTOPLEPHBB_p4()   {return _TopHadNonBJet_TOPTOPLEPHBB_p4;};

   TLorentzVector TopLepBJet_TOPHLEPBB_p4()      {return _TopLepBJet_TOPHLEPBB_p4;};
   TLorentzVector HiggsBJet1_TOPHLEPBB_p4()      {return _HiggsBJet1_TOPHLEPBB_p4;};
   TLorentzVector HiggsBJet2_TOPHLEPBB_p4()      {return _HiggsBJet2_TOPHLEPBB_p4;};
   TLorentzVector TopLepWLep_TOPHLEPBB_p4()      {return _TopLepWLep_TOPHLEPBB_p4;};
   TLorentzVector TopLepWNu_TOPHLEPBB_p4()       {return _TopLepWNu_TOPHLEPBB_p4;};
   TLorentzVector TopLep_TOPHLEPBB_p4()          {return _TopLep_TOPHLEPBB_p4;};
   TLorentzVector Higgs_TOPHLEPBB_p4()           {return _Higgs_TOPHLEPBB_p4;};
   TLorentzVector TopLepW_TOPHLEPBB_p4()         {return _TopLepW_TOPHLEPBB_p4;};

   TLorentzVector TopLepBJet_TOPTOPLEPHAD_p4()        {return _TopLepBJet_TOPTOPLEPHAD_p4;};
   TLorentzVector TopHadWNonBJet1_TOPTOPLEPHAD_p4()   {return _TopHadWNonBJet1_TOPTOPLEPHAD_p4;};
   TLorentzVector TopHadWNonBJet2_TOPTOPLEPHAD_p4()   {return _TopHadWNonBJet2_TOPTOPLEPHAD_p4;};
   TLorentzVector TopLepWLep_TOPTOPLEPHAD_p4()        {return _TopLepWLep_TOPTOPLEPHAD_p4;};
   TLorentzVector TopLepWNu_TOPTOPLEPHAD_p4()         {return _TopLepWNu_TOPTOPLEPHAD_p4;};
   TLorentzVector TopLep_TOPTOPLEPHAD_p4()            {return _TopLep_TOPTOPLEPHAD_p4;};
   TLorentzVector TopHad_TOPTOPLEPHAD_p4()            {return _TopHad_TOPTOPLEPHAD_p4;};
   TLorentzVector TopHadW_TOPTOPLEPHAD_p4()           {return _TopHadW_TOPTOPLEPHAD_p4;};
   TLorentzVector TopLepW_TOPTOPLEPHAD_p4()           {return _TopLepW_TOPTOPLEPHAD_p4;};
   TLorentzVector TopHadBJet_TOPTOPLEPHAD_p4()        {return _TopHadBJet_TOPTOPLEPHAD_p4;};
   
   double chi2_TOPTOPLEPHBB()       {return _chi2_TOPTOPLEPHBB;};
   double chi2_TOPHLEPBB()          {return _chi2_TOPHLEPBB;};
   double chi2_TOPTOPLEPHAD()       {return _chi2_TOPTOPLEPHAD;};

   double MVA_TOPTOPLEPHBB()       {return _MVA_TOPTOPLEPHBB;};
   double MVA_TOPHLEPBB()          {return _MVA_TOPHLEPBB;};
   double MVA_TOPTOPLEPHAD()       {return _MVA_TOPTOPLEPHAD;};
   
 protected:

   std::vector<Lepton>             *_v_Lepton;
   
   std::vector<Electron>           *_v_Electron;
   std::vector<Muon>               *_v_Muon;
   std::vector<Electron>           *_v_ElectronTight;
   std::vector<Muon>               *_v_MuonTight;
   
   std::vector<Event>              *_v_Event;
   std::vector<Jet>                *_v_Jet;
   std::vector<Jet>                *_v_BJet;
   std::vector<Jet>                *_v_NonBJet;
   std::vector<Jet>                *_v_JetTight;
   std::vector<Jet>                *_v_BJetTight;
   std::vector<Jet>                *_v_NonBJetTight;

   Helper *help;
   
 private:

   TLorentzVector _TopLepBJet_TOPTOPLEPHBB_p4;
   TLorentzVector _HiggsBJet1_TOPTOPLEPHBB_p4;
   TLorentzVector _HiggsBJet2_TOPTOPLEPHBB_p4;
   TLorentzVector _TopLepWLep_TOPTOPLEPHBB_p4;
   TLorentzVector _TopLepWNu_TOPTOPLEPHBB_p4;
   TLorentzVector _TopLep_TOPTOPLEPHBB_p4;
   TLorentzVector _Higgs_TOPTOPLEPHBB_p4;
   TLorentzVector _TopLepW_TOPTOPLEPHBB_p4;
   TLorentzVector _TopHadNonBJet_TOPTOPLEPHBB_p4;

   TLorentzVector _TopLepBJet_TOPHLEPBB_p4;
   TLorentzVector _HiggsBJet1_TOPHLEPBB_p4;
   TLorentzVector _HiggsBJet2_TOPHLEPBB_p4;
   TLorentzVector _TopLepWLep_TOPHLEPBB_p4;
   TLorentzVector _TopLepWNu_TOPHLEPBB_p4;
   TLorentzVector _TopLep_TOPHLEPBB_p4;
   TLorentzVector _Higgs_TOPHLEPBB_p4;
   TLorentzVector _TopLepW_TOPHLEPBB_p4;

   TLorentzVector _TopLepBJet_TOPTOPLEPHAD_p4;
   TLorentzVector _TopHadWNonBJet1_TOPTOPLEPHAD_p4;
   TLorentzVector _TopHadWNonBJet2_TOPTOPLEPHAD_p4;
   TLorentzVector _TopLepWLep_TOPTOPLEPHAD_p4;
   TLorentzVector _TopLepWNu_TOPTOPLEPHAD_p4;
   TLorentzVector _TopLep_TOPTOPLEPHAD_p4;
   TLorentzVector _TopHad_TOPTOPLEPHAD_p4;
   TLorentzVector _TopHadW_TOPTOPLEPHAD_p4;
   TLorentzVector _TopLepW_TOPTOPLEPHAD_p4;
   TLorentzVector _TopHadBJet_TOPTOPLEPHAD_p4;
   
   double _chi2_TOPTOPLEPHBB;
   double _chi2_TOPHLEPBB;
   double _chi2_TOPTOPLEPHAD;

   double _MVA_TOPTOPLEPHBB;
   double _MVA_TOPHLEPBB;
   double _MVA_TOPTOPLEPHAD;
   
   TRandom3 *rnd;

   KINFIT::kfit *kfTopTopLepHbb;
   KINFIT::kfit *kfTopHLepbb;
   KINFIT::kfit *kfTopTopLepHad;

   TMVA::Reader* MVAFullRecoReaderTOPTOPLEPHBB;
   TMVA::Reader* MVAPartRecoReaderTOPTOPLEPHBB;

   float MVAFullReco_HiggsRecM_TOPTOPLEPHBB;
   float MVAFullReco_TopLepRecM_TOPTOPLEPHBB;
   float MVAFullReco_HiggsTopLepRecDr_TOPTOPLEPHBB;
   float MVAFullReco_TopLepRecPt_TOPTOPLEPHBB;
   
   float MVAPartReco_HiggsRecM_TOPTOPLEPHBB;
   float MVAPartReco_TopLepRecMT_TOPTOPLEPHBB;
   float MVAPartReco_HiggsTopLepRecDphiT_TOPTOPLEPHBB;
   float MVAPartReco_TopLepRecPtT_TOPTOPLEPHBB;

   TMVA::Reader* MVAFullRecoReaderTOPHLEPBB;
   TMVA::Reader* MVAPartRecoReaderTOPHLEPBB;

   float MVAFullReco_HiggsRecM_TOPHLEPBB;
   float MVAFullReco_TopLepRecM_TOPHLEPBB;
   float MVAFullReco_HiggsTopLepRecDr_TOPHLEPBB;
   float MVAFullReco_TopLepRecPt_TOPHLEPBB;
   
   float MVAPartReco_HiggsRecM_TOPHLEPBB;
   float MVAPartReco_TopLepRecMT_TOPHLEPBB;
   float MVAPartReco_HiggsTopLepRecDphiT_TOPHLEPBB;
   float MVAPartReco_TopLepRecPtT_TOPHLEPBB;

   TMVA::Reader* MVAFullRecoReaderTOPTOPLEPHAD;
   TMVA::Reader* MVAPartRecoReaderTOPTOPLEPHAD;

   float MVAFullReco_TopHadRecM_TOPTOPLEPHAD;
   float MVAFullReco_TopLepRecM_TOPTOPLEPHAD;
   float MVAFullReco_TopLepTopHadRecDr_TOPTOPLEPHAD;
   float MVAFullReco_TopLepRecPt_TOPTOPLEPHAD;
   
   float MVAPartReco_TopHadRecM_TOPTOPLEPHAD;
   float MVAPartReco_TopLepRecMT_TOPTOPLEPHAD;
   float MVAPartReco_TopLepTopHadRecDphiT_TOPTOPLEPHAD;
   float MVAPartReco_TopLepRecPtT_TOPTOPLEPHAD;
};

#endif

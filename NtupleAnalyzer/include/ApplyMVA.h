#ifndef APPLYMVA_H
#define APPLYMVA_H

#include "Analyzer.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <vector>

using namespace std;

class ApplyMVA
{
   
 public:

   ApplyMVA(std::string home);
   virtual ~ApplyMVA();

   /*
   void setVariable_HiggsMass_TOPHLEPBB(float v)       {HiggsMass_TOPHLEPBB = v;};
   void setVariable_TopHadMass_TOPTOPLEPHAD(float v)       {TopHadMass_TOPTOPLEPHAD = v;};
   void setVariable_MVA_TOPHLEPBB(float v)       {MVA_TOPHLEPBB = v;};
   void setVariable_MVA_TOPTOPLEPHAD(float v)       {MVA_TOPTOPLEPHAD = v;};
   void setVariable_LepCharge(float v)       {LepCharge = v;};
   void setVariable_HiggsEta_TOPHLEPBB(float v)       {HiggsEta_TOPHLEPBB = v;};
   void setVariable_TopLepMass_TOPHLEPBB(float v)       {TopLepMass_TOPHLEPBB = v;};
   void setVariable_TopLepMass_TOPTOPLEPHAD(float v)       {TopLepMass_TOPTOPLEPHAD = v;};
   void setVariable_TopLepPt_TOPHLEPBB(float v)       {TopLepPt_TOPHLEPBB = v;};
   void setVariable_TopLepEta_TOPHLEPBB(float v)       {TopLepEta_TOPHLEPBB = v;};
   void setVariable_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB(float v)       {HiggsBJet1HiggsBJet2Dr_TOPHLEPBB = v;};
   void setVariable_TopLepHiggsDr_TOPHLEPBB(float v)       {TopLepHiggsDr_TOPHLEPBB = v;};
   void setVariable_TopLepTopHadDr_TOPTOPLEPHAD(float v)       {TopLepTopHadDr_TOPTOPLEPHAD = v;};
   
   void setVariable_HiggsMass_TOPTOPLEPHBB(float v)       {HiggsMass_TOPTOPLEPHBB = v;};
   void setVariable_MVA_TOPTOPLEPHBB(float v)       {MVA_TOPTOPLEPHBB = v;};
   void setVariable_HiggsEta_TOPTOPLEPHBB(float v)       {HiggsEta_TOPTOPLEPHBB = v;};
   void setVariable_TopLepMass_TOPTOPLEPHBB(float v)       {TopLepMass_TOPTOPLEPHBB = v;};
   void setVariable_TopLepPt_TOPTOPLEPHBB(float v)       {TopLepPt_TOPTOPLEPHBB = v;};
   void setVariable_TopLepEta_TOPTOPLEPHBB(float v)       {TopLepEta_TOPTOPLEPHBB = v;};
   void setVariable_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB(float v)       {HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = v;};
   void setVariable_TopLepHiggsDr_TOPTOPLEPHBB(float v)       {TopLepHiggsDr_TOPTOPLEPHBB = v;};
   */

   void setVariable_HiggsMass_TOPHLEPBB(float v)       {HiggsMass_TOPHLEPBB = v;};
   void setVariable_MVA_TOPHLEPBB(float v)       {MVA_TOPHLEPBB = v;};
   void setVariable_HiggsEta_TOPHLEPBB(float v)       {HiggsEta_TOPHLEPBB = v;};
   void setVariable_TopLepMass_TOPHLEPBB(float v)       {TopLepMass_TOPHLEPBB = v;};
   void setVariable_TopLepPt_TOPHLEPBB(float v)       {TopLepPt_TOPHLEPBB = v;};
   void setVariable_TopLepEta_TOPHLEPBB(float v)       {TopLepEta_TOPHLEPBB = v;};
   void setVariable_HiggsBJet1HiggsBJet2Dr_TOPHLEPBB(float v)       {HiggsBJet1HiggsBJet2Dr_TOPHLEPBB = v;};
   void setVariable_TopLepHiggsDr_TOPHLEPBB(float v)       {TopLepHiggsDr_TOPHLEPBB = v;};
   void setVariable_HiggsBJet1CSVv2_TOPHLEPBB(float v)       {HiggsBJet1CSVv2_TOPHLEPBB = v;};
   void setVariable_HiggsBJet2CSVv2_TOPHLEPBB(float v)       {HiggsBJet2CSVv2_TOPHLEPBB = v;};
   void setVariable_TopLepBJetCSVv2_TOPHLEPBB(float v)       {TopLepBJetCSVv2_TOPHLEPBB = v;};
   void setVariable_TopHadMass_TOPTOPLEPHAD(float v)       {TopHadMass_TOPTOPLEPHAD = v;};
   void setVariable_MVA_TOPTOPLEPHAD(float v)       {MVA_TOPTOPLEPHAD = v;};
   void setVariable_TopLepMass_TOPTOPLEPHAD(float v)       {TopLepMass_TOPTOPLEPHAD = v;};
   void setVariable_TopLepTopHadDr_TOPTOPLEPHAD(float v)       {TopLepTopHadDr_TOPTOPLEPHAD = v;};
   void setVariable_TopLepBJetCSVv2_TOPTOPLEPHAD(float v)       {TopLepBJetCSVv2_TOPTOPLEPHAD = v;};
   void setVariable_TopHadBJetCSVv2_TOPTOPLEPHAD(float v)       {TopHadBJetCSVv2_TOPTOPLEPHAD = v;};
   void setVariable_TopHadWNonBJet1CSVv2_TOPTOPLEPHAD(float v)       {TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = v;};
   void setVariable_TopHadWNonBJet2CSVv2_TOPTOPLEPHAD(float v)       {TopHadWNonBJet2CSVv2_TOPTOPLEPHAD = v;};
   void setVariable_HiggsMass_TOPTOPLEPHBB(float v)       {HiggsMass_TOPTOPLEPHBB = v;};
   void setVariable_MVA_TOPTOPLEPHBB(float v)       {MVA_TOPTOPLEPHBB = v;};
   void setVariable_TopLepMass_TOPTOPLEPHBB(float v)       {TopLepMass_TOPTOPLEPHBB = v;};
   void setVariable_HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB(float v)       {HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = v;};
   void setVariable_TopLepHiggsDr_TOPTOPLEPHBB(float v)       {TopLepHiggsDr_TOPTOPLEPHBB = v;};
   void setVariable_HiggsBJet1CSVv2_TOPTOPLEPHBB(float v)       {HiggsBJet1CSVv2_TOPTOPLEPHBB = v;};
   void setVariable_HiggsBJet2CSVv2_TOPTOPLEPHBB(float v)       {HiggsBJet2CSVv2_TOPTOPLEPHBB = v;};
   void setVariable_TopLepBJetCSVv2_TOPTOPLEPHBB(float v)       {TopLepBJetCSVv2_TOPTOPLEPHBB = v;};
   void setVariable_TopHadNonBJetCSVv2_TOPTOPLEPHBB(float v)       {TopHadNonBJetCSVv2_TOPTOPLEPHBB = v;};
   void setVariable_LepCharge(float v)       {LepCharge = v;};
   
   void init();
   double run(std::string chan);

 protected:
   /*
   TMVA::Reader* b3j4HutSTreader;
   TMVA::Reader* b3j4HctSTreader;
   TMVA::Reader* b3j4HutTTreader;
   TMVA::Reader* b3j4HctTTreader;

   TMVA::Reader* b3j3HutSTreader;
   TMVA::Reader* b3j3HctSTreader;
   TMVA::Reader* b3j3HutTTreader;
   TMVA::Reader* b3j3HctTTreader;

   TMVA::Reader* b2j4HutSTreader;
   TMVA::Reader* b2j4HctSTreader;
   TMVA::Reader* b2j4HutTTreader;
   TMVA::Reader* b2j4HctTTreader;
   */

   TMVA::Reader* b0j3HctSTreader;
   TMVA::Reader* b1j3HctSTreader;
   TMVA::Reader* b2j3HctSTreader;
   TMVA::Reader* b3j3HctSTreader;
   TMVA::Reader* b0j4HctSTreader;
   TMVA::Reader* b1j4HctSTreader;
   TMVA::Reader* b2j4HctSTreader;
   TMVA::Reader* b3j4HctSTreader;
   TMVA::Reader* b4j4HctSTreader;
   TMVA::Reader* b0j3HctTTreader;
   TMVA::Reader* b1j3HctTTreader;
   TMVA::Reader* b2j3HctTTreader;
   TMVA::Reader* b3j3HctTTreader;
   TMVA::Reader* b0j4HctTTreader;
   TMVA::Reader* b1j4HctTTreader;
   TMVA::Reader* b2j4HctTTreader;
   TMVA::Reader* b3j4HctTTreader;
   TMVA::Reader* b4j4HctTTreader;
   TMVA::Reader* b0j3HutSTreader;
   TMVA::Reader* b1j3HutSTreader;
   TMVA::Reader* b2j3HutSTreader;
   TMVA::Reader* b3j3HutSTreader;
   TMVA::Reader* b0j4HutSTreader;
   TMVA::Reader* b1j4HutSTreader;
   TMVA::Reader* b2j4HutSTreader;
   TMVA::Reader* b3j4HutSTreader;
   TMVA::Reader* b4j4HutSTreader;
   TMVA::Reader* b0j3HutTTreader;
   TMVA::Reader* b1j3HutTTreader;
   TMVA::Reader* b2j3HutTTreader;
   TMVA::Reader* b3j3HutTTreader;
   TMVA::Reader* b0j4HutTTreader;
   TMVA::Reader* b1j4HutTTreader;
   TMVA::Reader* b2j4HutTTreader;
   TMVA::Reader* b3j4HutTTreader;
   TMVA::Reader* b4j4HutTTreader;
   
   std::string _home;

   double _disc;

   //vector<string> new_categories;
   //vector<vector<string> > new_variables;
   
 private:

   /*
   float HiggsMass_TOPHLEPBB;
   float TopHadMass_TOPTOPLEPHAD;
   float MVA_TOPHLEPBB;
   float MVA_TOPTOPLEPHAD;
   float LepCharge;
   float HiggsEta_TOPHLEPBB;
   float TopLepMass_TOPHLEPBB;
   float TopLepMass_TOPTOPLEPHAD;
   float TopLepPt_TOPHLEPBB;
   float TopLepEta_TOPHLEPBB;
   float HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
   float TopLepHiggsDr_TOPHLEPBB;
   float TopLepTopHadDr_TOPTOPLEPHAD;

   float HiggsMass_TOPTOPLEPHBB;
   float MVA_TOPTOPLEPHBB;
   float HiggsEta_TOPTOPLEPHBB;
   float TopLepMass_TOPTOPLEPHBB;
   float TopLepPt_TOPTOPLEPHBB;
   float TopLepEta_TOPTOPLEPHBB;
   float HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
   float TopLepHiggsDr_TOPTOPLEPHBB;
   */

   float HiggsBJet1CSVv2_TOPHLEPBB;
   float HiggsBJet1CSVv2_TOPTOPLEPHBB;
   float HiggsBJet1HiggsBJet2Dr_TOPHLEPBB;
   float HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
   float HiggsBJet2CSVv2_TOPHLEPBB;
   float HiggsBJet2CSVv2_TOPTOPLEPHBB;
   float HiggsEta_TOPHLEPBB;
   float HiggsMass_TOPHLEPBB;
   float HiggsMass_TOPTOPLEPHBB;
   float LepCharge;
   float MVA_TOPHLEPBB;
   float MVA_TOPTOPLEPHAD;
   float MVA_TOPTOPLEPHBB;
   float TopHadBJetCSVv2_TOPTOPLEPHAD;
   float TopHadMass_TOPTOPLEPHAD;
   float TopHadNonBJetCSVv2_TOPTOPLEPHBB;
   float TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
   float TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
   float TopLepBJetCSVv2_TOPHLEPBB;
   float TopLepBJetCSVv2_TOPTOPLEPHAD;
   float TopLepBJetCSVv2_TOPTOPLEPHBB;
   float TopLepEta_TOPHLEPBB;
   float TopLepHiggsDr_TOPHLEPBB;
   float TopLepHiggsDr_TOPTOPLEPHBB;
   float TopLepMass_TOPHLEPBB;
   float TopLepMass_TOPTOPLEPHAD;
   float TopLepMass_TOPTOPLEPHBB;
   float TopLepPt_TOPHLEPBB;
   float TopLepTopHadDr_TOPTOPLEPHAD;

};

#endif

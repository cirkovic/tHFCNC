#include "../include/ApplyMVA.h"

ApplyMVA::ApplyMVA(std::string home)
{
   _home = home;
}

ApplyMVA::~ApplyMVA()
{
   /*
   delete b3j4HutSTreader;
   delete b3j4HctSTreader;
   delete b3j4HutTTreader;
   delete b3j4HctTTreader;

   delete b3j3HutSTreader;
   delete b3j3HctSTreader;
   delete b3j3HutTTreader;
   delete b3j3HctTTreader;

   delete b2j4HutSTreader;
   delete b2j4HctSTreader;
   delete b2j4HutTTreader;
   delete b2j4HctTTreader;
   */
   delete b0j3HctSTreader;
   delete b1j3HctSTreader;
   delete b2j3HctSTreader;
   delete b3j3HctSTreader;

   delete b0j4HctSTreader;
   delete b1j4HctSTreader;
   delete b2j4HctSTreader;
   delete b3j4HctSTreader;
   delete b4j4HctSTreader;

   delete b0j3HctTTreader;
   delete b1j3HctTTreader;
   delete b2j3HctTTreader;
   delete b3j3HctTTreader;

   delete b0j4HctTTreader;
   delete b1j4HctTTreader;
   delete b2j4HctTTreader;
   delete b3j4HctTTreader;
   delete b4j4HctTTreader;

   delete b0j3HutSTreader;
   delete b1j3HutSTreader;
   delete b2j3HutSTreader;
   delete b3j3HutSTreader;

   delete b0j4HutSTreader;
   delete b1j4HutSTreader;
   delete b2j4HutSTreader;
   delete b3j4HutSTreader;
   delete b4j4HutSTreader;

   delete b0j3HutTTreader;
   delete b1j3HutTTreader;
   delete b2j3HutTTreader;
   delete b3j3HutTTreader;

   delete b0j4HutTTreader;
   delete b1j4HutTTreader;
   delete b2j4HutTTreader;
   delete b3j4HutTTreader;
   delete b4j4HutTTreader;

}

void ApplyMVA::init()
{

   /*

   b3j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HutSTreader->AddVariable("LepCharge",&LepCharge);
   b3j4HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);

   b3j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HctSTreader->AddVariable("LepCharge",&LepCharge);
   b3j4HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);

   b3j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HutTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("LepCharge",&LepCharge);
   b3j4HutTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);

   b3j4HctTTreader = new TMVA::Reader("!Color:!Silent");

   b3j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HctTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("LepCharge",&LepCharge);
   b3j4HctTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
  
   std::string b3j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b3j4_BDT.weights.xml";
   b3j4HutSTreader->BookMVA("BDTG method",b3j4HutSTweightsFile.c_str());
   std::string b3j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b3j4_BDT.weights.xml";
   b3j4HctSTreader->BookMVA("BDTG method",b3j4HctSTweightsFile.c_str());
   std::string b3j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b3j4_BDT.weights.xml";
   b3j4HutTTreader->BookMVA("BDTG method",b3j4HutTTweightsFile.c_str());
   std::string b3j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b3j4_BDT.weights.xml";
   b3j4HctTTreader->BookMVA("BDTG method",b3j4HctTTweightsFile.c_str());

   b3j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j3HutSTreader->AddVariable("LepCharge",&LepCharge);
   b3j3HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);

   b3j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j3HctSTreader->AddVariable("LepCharge",&LepCharge);
   b3j3HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);

   b3j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j3HutTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j3HutTTreader->AddVariable("LepCharge",&LepCharge);
   b3j3HutTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);

   b3j3HctTTreader = new TMVA::Reader("!Color:!Silent");
 
   b3j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j3HctTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j3HctTTreader->AddVariable("LepCharge",&LepCharge);
   b3j3HctTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   
   std::string b3j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b3j3_BDT.weights.xml";
   b3j3HutSTreader->BookMVA("BDTG method",b3j3HutSTweightsFile.c_str());
   std::string b3j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b3j3_BDT.weights.xml";
   b3j3HctSTreader->BookMVA("BDTG method",b3j3HctSTweightsFile.c_str());
   std::string b3j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b3j3_BDT.weights.xml";
   b3j3HutTTreader->BookMVA("BDTG method",b3j3HutTTweightsFile.c_str());
   std::string b3j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b3j3_BDT.weights.xml";
   b3j3HctTTreader->BookMVA("BDTG method",b3j3HctTTweightsFile.c_str());

   b2j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HutSTreader->AddVariable("LepCharge",&LepCharge);
   b2j4HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);

   b2j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HctSTreader->AddVariable("LepCharge",&LepCharge);
   b2j4HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);

   b2j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HutTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("LepCharge",&LepCharge);
   b2j4HutTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);

   b2j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HctTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("LepCharge",&LepCharge);
   b2j4HctTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   
   std::string b2j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b2j4_BDT.weights.xml";
   b2j4HutSTreader->BookMVA("BDTG method",b2j4HutSTweightsFile.c_str());
   std::string b2j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b2j4_BDT.weights.xml";
   b2j4HctSTreader->BookMVA("BDTG method",b2j4HctSTweightsFile.c_str());
   std::string b2j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b2j4_BDT.weights.xml";
   b2j4HutTTreader->BookMVA("BDTG method",b2j4HutTTweightsFile.c_str());
   std::string b2j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b2j4_BDT.weights.xml";
   b2j4HctTTreader->BookMVA("BDTG method",b2j4HctTTweightsFile.c_str());

   */

   b0j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   b0j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j3HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b0j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b0j3_BDT.weights.xml";
   b0j3HctSTreader->BookMVA("BDTG method",b0j3HctSTweightsFile.c_str());

   //std::cout << "CIRKOVIC EXIT" << std::endl;    exit(0);

   b1j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   b1j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j3HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b1j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b1j3_BDT.weights.xml";
   b1j3HctSTreader->BookMVA("BDTG method",b1j3HctSTweightsFile.c_str());

   b2j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   b2j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j3HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b2j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b2j3_BDT.weights.xml";
   b2j3HctSTreader->BookMVA("BDTG method",b2j3HctSTweightsFile.c_str());

   b3j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   b3j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b3j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b3j3_BDT.weights.xml";
   b3j3HctSTreader->BookMVA("BDTG method",b3j3HctSTweightsFile.c_str());

   b0j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   b0j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b0j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b0j4HctSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b0j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b0j4_BDT.weights.xml";
   b0j4HctSTreader->BookMVA("BDTG method",b0j4HctSTweightsFile.c_str());

   b1j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   b1j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b1j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b1j4HctSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b1j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b1j4_BDT.weights.xml";
   b1j4HctSTreader->BookMVA("BDTG method",b1j4HctSTweightsFile.c_str());

   //std::cout << "CIRKOVIC EXIT" << std::endl;    exit(0);

   b2j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   b2j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b2j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b2j4_BDT.weights.xml";
   b2j4HctSTreader->BookMVA("BDTG method",b2j4HctSTweightsFile.c_str());

   b3j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   b3j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b3j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b3j4_BDT.weights.xml";
   b3j4HctSTreader->BookMVA("BDTG method",b3j4HctSTweightsFile.c_str());

   b4j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   b4j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b4j4HctSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b4j4HctSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b4j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_all_b4j4_BDT.weights.xml";
   b4j4HctSTreader->BookMVA("BDTG method",b4j4HctSTweightsFile.c_str());

   b0j3HctTTreader = new TMVA::Reader("!Color:!Silent");
   b0j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j3HctTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b0j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b0j3_BDT.weights.xml";
   b0j3HctTTreader->BookMVA("BDTG method",b0j3HctTTweightsFile.c_str());

   //std::cout << "CIRKOVIC EXIT" << std::endl;    exit(0);

   b1j3HctTTreader = new TMVA::Reader("!Color:!Silent");
   b1j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j3HctTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b1j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b1j3_BDT.weights.xml";
   b1j3HctTTreader->BookMVA("BDTG method",b1j3HctTTweightsFile.c_str());

   b2j3HctTTreader = new TMVA::Reader("!Color:!Silent");
   b2j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j3HctTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b2j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b2j3_BDT.weights.xml";
   b2j3HctTTreader->BookMVA("BDTG method",b2j3HctTTweightsFile.c_str());

   b3j3HctTTreader = new TMVA::Reader("!Color:!Silent");
   b3j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b3j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b3j3_BDT.weights.xml";
   b3j3HctTTreader->BookMVA("BDTG method",b3j3HctTTweightsFile.c_str());

   b0j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   b0j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b0j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b0j4HctTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b0j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b0j4_BDT.weights.xml";
   b0j4HctTTreader->BookMVA("BDTG method",b0j4HctTTweightsFile.c_str());

   b1j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   b1j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b1j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b1j4HctTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b1j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b1j4_BDT.weights.xml";
   b1j4HctTTreader->BookMVA("BDTG method",b1j4HctTTweightsFile.c_str());

   b2j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   b2j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b2j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b2j4_BDT.weights.xml";
   b2j4HctTTreader->BookMVA("BDTG method",b2j4HctTTweightsFile.c_str());

   //std::cout << "CIRKOVIC EXIT" << std::endl;    exit(0);

   b3j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   b3j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b3j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b3j4_BDT.weights.xml";
   b3j4HctTTreader->BookMVA("BDTG method",b3j4HctTTweightsFile.c_str());

   b4j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   b4j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b4j4HctTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b4j4HctTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b4j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_all_b4j4_BDT.weights.xml";
   b4j4HctTTreader->BookMVA("BDTG method",b4j4HctTTweightsFile.c_str());

   b0j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   b0j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("LepCharge", &LepCharge);
   b0j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j3HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b0j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b0j3_BDT.weights.xml";
   b0j3HutSTreader->BookMVA("BDTG method",b0j3HutSTweightsFile.c_str());

   b1j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   b1j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("LepCharge", &LepCharge);
   b1j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j3HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b1j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b1j3_BDT.weights.xml";
   b1j3HutSTreader->BookMVA("BDTG method",b1j3HutSTweightsFile.c_str());

   b2j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   b2j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("LepCharge", &LepCharge);
   b2j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j3HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b2j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b2j3_BDT.weights.xml";
   b2j3HutSTreader->BookMVA("BDTG method",b2j3HutSTweightsFile.c_str());

   b3j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   b3j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("LepCharge", &LepCharge);
   b3j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b3j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b3j3_BDT.weights.xml";
   b3j3HutSTreader->BookMVA("BDTG method",b3j3HutSTweightsFile.c_str());

   b0j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   b0j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("LepCharge", &LepCharge);
   b0j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b0j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b0j4HutSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b0j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b0j4_BDT.weights.xml";
   b0j4HutSTreader->BookMVA("BDTG method",b0j4HutSTweightsFile.c_str());

   b1j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   b1j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("LepCharge", &LepCharge);
   b1j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b1j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b1j4HutSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b1j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b1j4_BDT.weights.xml";
   b1j4HutSTreader->BookMVA("BDTG method",b1j4HutSTweightsFile.c_str());

   b2j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   b2j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("LepCharge", &LepCharge);
   b2j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b2j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b2j4_BDT.weights.xml";
   b2j4HutSTreader->BookMVA("BDTG method",b2j4HutSTweightsFile.c_str());

   b3j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   b3j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("LepCharge", &LepCharge);
   b3j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b3j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b3j4_BDT.weights.xml";
   b3j4HutSTreader->BookMVA("BDTG method",b3j4HutSTweightsFile.c_str());

   b4j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   b4j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD", &TopHadMass_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("LepCharge", &LepCharge);
   b4j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB", &HiggsEta_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB", &TopLepPt_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB", &TopLepEta_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   b4j4HutSTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b4j4HutSTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b4j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_all_b4j4_BDT.weights.xml";
   b4j4HutSTreader->BookMVA("BDTG method",b4j4HutSTweightsFile.c_str());

   b0j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   b0j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b0j3HutTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b0j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b0j3_BDT.weights.xml";
   b0j3HutTTreader->BookMVA("BDTG method",b0j3HutTTweightsFile.c_str());

   b1j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   b1j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b1j3HutTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b1j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b1j3_BDT.weights.xml";
   b1j3HutTTreader->BookMVA("BDTG method",b1j3HutTTweightsFile.c_str());

   b2j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   b2j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b2j3HutTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b2j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b2j3_BDT.weights.xml";
   b2j3HutTTreader->BookMVA("BDTG method",b2j3HutTTweightsFile.c_str());

   b3j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   b3j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB", &HiggsMass_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("MVA_TOPHLEPBB", &MVA_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepMass_TOPHLEPBB", &TopLepMass_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB", &HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB", &TopLepHiggsDr_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB", &HiggsBJet1CSVv2_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB", &HiggsBJet2CSVv2_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepBJetCSVv2_TOPHLEPBB", &TopLepBJetCSVv2_TOPHLEPBB);
   std::string b3j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b3j3_BDT.weights.xml";
   b3j3HutTTreader->BookMVA("BDTG method",b3j3HutTTweightsFile.c_str());

   b0j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   b0j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b0j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b0j4HutTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b0j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b0j4_BDT.weights.xml";
   b0j4HutTTreader->BookMVA("BDTG method",b0j4HutTTweightsFile.c_str());

   b1j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   b1j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b1j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b1j4HutTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b1j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b1j4_BDT.weights.xml";
   b1j4HutTTreader->BookMVA("BDTG method",b1j4HutTTweightsFile.c_str());

   b2j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   b2j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b2j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b2j4_BDT.weights.xml";
   b2j4HutTTreader->BookMVA("BDTG method",b2j4HutTTweightsFile.c_str());

   b3j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   b3j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b3j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b3j4_BDT.weights.xml";
   b3j4HutTTreader->BookMVA("BDTG method",b3j4HutTTweightsFile.c_str());

   b4j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   b4j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB", &HiggsMass_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB", &MVA_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD", &MVA_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB", &TopLepMass_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD", &TopLepMass_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB", &HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB", &TopLepHiggsDr_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD", &TopLepTopHadDr_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB", &HiggsBJet1CSVv2_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB", &HiggsBJet2CSVv2_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB", &TopLepBJetCSVv2_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB", &TopHadNonBJetCSVv2_TOPTOPLEPHBB);
   b4j4HutTTreader->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD", &TopLepBJetCSVv2_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD", &TopHadBJetCSVv2_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
   b4j4HutTTreader->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD", &TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
   std::string b4j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_all_b4j4_BDT.weights.xml";
   b4j4HutTTreader->BookMVA("BDTG method",b4j4HutTTweightsFile.c_str());

   //std::cout << "CIRKOVIC EXIT" << std::endl;    exit(0);

   std::cout << "ApplyMVA initialisation done" << std::endl;
}

double ApplyMVA::run(std::string chan)
{
   _disc = 10E+10;
   
   /*
   if( chan == "b3j4HutST" ) _disc = b3j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HctST" ) _disc = b3j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HutTT" ) _disc = b3j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HctTT" ) _disc = b3j4HctTTreader->EvaluateMVA("BDTG method");

   else if( chan == "b3j3HutST" ) _disc = b3j3HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HctST" ) _disc = b3j3HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HutTT" ) _disc = b3j3HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HctTT" ) _disc = b3j3HctTTreader->EvaluateMVA("BDTG method");

   else if( chan == "b2j4HutST" ) _disc = b2j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HctST" ) _disc = b2j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HutTT" ) _disc = b2j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HctTT" ) _disc = b2j4HctTTreader->EvaluateMVA("BDTG method");
   */

   if( chan == "b0j3HctST" ) _disc = b0j3HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j3HctST" ) _disc = b1j3HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j3HctST" ) _disc = b2j3HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HctST" ) _disc = b3j3HctSTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j4HctST" ) _disc = b0j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j4HctST" ) _disc = b1j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HctST" ) _disc = b2j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HctST" ) _disc = b3j4HctSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b4j4HctST" ) _disc = b4j4HctSTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j3HctTT" ) _disc = b0j3HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j3HctTT" ) _disc = b1j3HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j3HctTT" ) _disc = b2j3HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HctTT" ) _disc = b3j3HctTTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j4HctTT" ) _disc = b0j4HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j4HctTT" ) _disc = b1j4HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HctTT" ) _disc = b2j4HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HctTT" ) _disc = b3j4HctTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b4j4HctTT" ) _disc = b4j4HctTTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j3HutST" ) _disc = b0j3HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j3HutST" ) _disc = b1j3HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j3HutST" ) _disc = b2j3HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HutST" ) _disc = b3j3HutSTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j4HutST" ) _disc = b0j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j4HutST" ) _disc = b1j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HutST" ) _disc = b2j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HutST" ) _disc = b3j4HutSTreader->EvaluateMVA("BDTG method");
   else if( chan == "b4j4HutST" ) _disc = b4j4HutSTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j3HutTT" ) _disc = b0j3HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j3HutTT" ) _disc = b1j3HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j3HutTT" ) _disc = b2j3HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j3HutTT" ) _disc = b3j3HutTTreader->EvaluateMVA("BDTG method");

   else if( chan == "b0j4HutTT" ) _disc = b0j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b1j4HutTT" ) _disc = b1j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b2j4HutTT" ) _disc = b2j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b3j4HutTT" ) _disc = b3j4HutTTreader->EvaluateMVA("BDTG method");
   else if( chan == "b4j4HutTT" ) _disc = b4j4HutTTreader->EvaluateMVA("BDTG method");
   
   return _disc;
}

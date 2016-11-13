#include "../include/ApplyMVA.h"

ApplyMVA::ApplyMVA(std::string home)
{
   _home = home;
}

ApplyMVA::~ApplyMVA()
{
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
}

void ApplyMVA::init()
{
   b3j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("LepCharge",&LepCharge);
   b3j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b3j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
   b3j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b3j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b3j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
   b3j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b3j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b3j4HutTTreader->AddVariable("TopLepPt_TOPTOPLEPHBB",&TopLepPt_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
   b3j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b3j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b3j4HctTTreader->AddVariable("TopLepPt_TOPTOPLEPHBB",&TopLepPt_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
   b3j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
   
   std::string b3j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b3j4_BDT.weights.xml";
   b3j4HutSTreader->BookMVA("BDTG method",b3j4HutSTweightsFile.c_str());
   std::string b3j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b3j4_BDT.weights.xml";
   b3j4HctSTreader->BookMVA("BDTG method",b3j4HctSTweightsFile.c_str());
   std::string b3j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b3j4_BDT.weights.xml";
   b3j4HutTTreader->BookMVA("BDTG method",b3j4HutTTweightsFile.c_str());
   std::string b3j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b3j4_BDT.weights.xml";
   b3j4HctTTreader->BookMVA("BDTG method",b3j4HctTTweightsFile.c_str());

   b3j3HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("LepCharge",&LepCharge);
   b3j3HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);

   b3j3HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);

   b3j3HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HutTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HutTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);

   b3j3HctTTreader = new TMVA::Reader("!Color:!Silent");
   
   b3j3HctTTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b3j3HctTTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
   
   std::string b3j3HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b3j3_BDT.weights.xml";
   b3j3HutSTreader->BookMVA("BDTG method",b3j3HutSTweightsFile.c_str());
   std::string b3j3HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b3j3_BDT.weights.xml";
   b3j3HctSTreader->BookMVA("BDTG method",b3j3HctSTweightsFile.c_str());
   std::string b3j3HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b3j3_BDT.weights.xml";
   b3j3HutTTreader->BookMVA("BDTG method",b3j3HutTTweightsFile.c_str());
   std::string b3j3HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b3j3_BDT.weights.xml";
   b3j3HctTTreader->BookMVA("BDTG method",b3j3HctTTweightsFile.c_str());

   b2j4HutSTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HutSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("LepCharge",&LepCharge);
   b2j4HutSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b2j4HutSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
   b2j4HutSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b2j4HctSTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HctSTreader->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b2j4HctSTreader->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB);
   b2j4HctSTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b2j4HutTTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HutTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b2j4HutTTreader->AddVariable("TopLepPt_TOPTOPLEPHBB",&TopLepPt_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HutTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
   b2j4HutTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);

   b2j4HctTTreader = new TMVA::Reader("!Color:!Silent");
   
   b2j4HctTTreader->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("HiggsEta_TOPTOPLEPHBB",&HiggsEta_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
   b2j4HctTTreader->AddVariable("TopLepPt_TOPTOPLEPHBB",&TopLepPt_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepEta_TOPTOPLEPHBB",&TopLepEta_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB);
   b2j4HctTTreader->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
   b2j4HctTTreader->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
   
   std::string b2j4HutSTweightsFile = _home+"/MVA/weights/TMVA_HutST_b2j4_BDT.weights.xml";
   b2j4HutSTreader->BookMVA("BDTG method",b2j4HutSTweightsFile.c_str());
   std::string b2j4HctSTweightsFile = _home+"/MVA/weights/TMVA_HctST_b2j4_BDT.weights.xml";
   b2j4HctSTreader->BookMVA("BDTG method",b2j4HctSTweightsFile.c_str());
   std::string b2j4HutTTweightsFile = _home+"/MVA/weights/TMVA_HutTT_b2j4_BDT.weights.xml";
   b2j4HutTTreader->BookMVA("BDTG method",b2j4HutTTweightsFile.c_str());
   std::string b2j4HctTTweightsFile = _home+"/MVA/weights/TMVA_HctTT_b2j4_BDT.weights.xml";
   b2j4HctTTreader->BookMVA("BDTG method",b2j4HctTTweightsFile.c_str());

   std::cout << "ApplyMVA initialisation done" << std::endl;
}

double ApplyMVA::run(std::string chan)
{
   _disc = 10E+10;
   
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
   
   return _disc;
}

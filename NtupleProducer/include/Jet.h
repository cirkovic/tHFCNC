#ifndef JET_H
#define JET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Jet : public Base
{
 public:
   Jet();
   virtual ~Jet();

   static bool sortPtPredicate(Jet lhs, Jet rhs)
     {return (lhs.pt() > rhs.pt());};

   static bool sortCSVv2Predicate(Jet lhs, Jet rhs)
     {return (lhs.CSVv2() > rhs.CSVv2());};

   static bool sortcMVAv2Predicate(Jet lhs, Jet rhs)
     {return (lhs.cMVAv2() > rhs.cMVAv2());};
   
   int ID()    {return _ID;};
   
   void sel();
   
   // kinematics
   float E()         {return _E;};
   float pt()        {return _pt;};
   float eta()       {return _eta;};
   float phi()       {return _phi;};
   float m()         {return _m;};

   float gen_E()         {return _gen_E;};
   float gen_pt()        {return _gen_pt;};
   float gen_eta()       {return _gen_eta;};
   float gen_phi()       {return _gen_phi;};
   float gen_m()         {return _gen_m;};
   int gen_status()      {return _gen_status;};
   int gen_id()          {return _gen_id;};
   
   TLorentzVector p4()  {return _p4;};
   
   float CSVv2()         {return _CSVv2;};
   float cMVAv2()         {return _cMVAv2;};
   
   float CharmCvsL()         {return _CharmCvsL;};
   float CharmCvsB()         {return _CharmCvsB;};
   
   bool isLoose()         {return _isLoose;};
   bool isTight()         {return _isTight;};
   
   bool isBTag()         {return _isBTag;};
   
   bool isLooseJetID()         {return _isLooseJetID;};
   bool isTightJetID()         {return _isTightJetID;};
   bool passElecOverlap()      {return _passElecOverlap;};
   bool passMuonOverlap()      {return _passMuonOverlap;};
   
   int hadronFlavour()         {return _hadronFlavour;};
   
   bool hasGenJet()         {return _hasGenJet;};
   float genJetPt()         {return _genJetPt;};
   float genJetEta()        {return _genJetEta;};
   float genJetPhi()        {return _genJetPhi;};
   float genJetE()          {return _genJetE;};
   
   float SfIterativeFitCentral()         {return _SfIterativeFitCentral;};
   float SfIterativeFitJesUp()         {return _SfIterativeFitJesUp;};
   float SfIterativeFitJesDown()         {return _SfIterativeFitJesDown;};
   float SfIterativeFitLfUp()         {return _SfIterativeFitLfUp;};
   float SfIterativeFitLfDown()         {return _SfIterativeFitLfDown;};
   float SfIterativeFitHfstats1Up()         {return _SfIterativeFitHfstats1Up;};
   float SfIterativeFitHfstats1Down()         {return _SfIterativeFitHfstats1Down;};
   float SfIterativeFitHfstats2Up()         {return _SfIterativeFitHfstats2Up;};
   float SfIterativeFitHfstats2Down()         {return _SfIterativeFitHfstats2Down;};
   float SfIterativeFitCferr1Up()         {return _SfIterativeFitCferr1Up;};
   float SfIterativeFitCferr1Down()         {return _SfIterativeFitCferr1Down;};
   float SfIterativeFitCferr2Up()         {return _SfIterativeFitCferr2Up;};
   float SfIterativeFitCferr2Down()         {return _SfIterativeFitCferr2Down;};

   float r_AK4PF_pt()             {return _r_AK4PF_pt;};
   float r_AK4PchsF_pt()          {return _r_AK4PchsF_pt;};
   float r_AK8PF_pt()             {return _r_AK8PF_pt;};
   float r_AK8PFchs_pt()          {return _r_AK8PFchs_pt;};
   float r_AK4PF_phi()            {return _r_AK4PF_phi;};
   float r_AK4PchsF_phi()         {return _r_AK4PchsF_phi;};
   float r_AK8PF_phi()            {return _r_AK8PF_phi;};
   float r_AK8PFchs_phi()         {return _r_AK8PFchs_phi;};
   
   void read();
   void init();
	
 protected:

   int _ID;
	
   float _E;
   float _pt;
   float _eta;
   float _phi;
   float _m;

   float _gen_E;
   float _gen_pt;
   float _gen_eta;
   float _gen_phi;
   float _gen_m;
   int _gen_status;
   int _gen_id;
   
   TLorentzVector _p4;

   float _CSVv2;
   float _cMVAv2;
   
   float _CharmCvsL;
   float _CharmCvsB;

   bool _isLoose;
   bool _isTight;
   
   bool _isBTag;

   bool _isLooseJetID;
   bool _isTightJetID;
   
   int _hadronFlavour;
   
   bool _hasGenJet;
   float _genJetPt;
   float _genJetEta;
   float _genJetPhi;
   float _genJetE;
   
   bool _passElecOverlap;
   bool _passMuonOverlap;
   
   float _SfIterativeFitCentral;
   float _SfIterativeFitJesUp;
   float _SfIterativeFitJesDown;
   float _SfIterativeFitLfUp;
   float _SfIterativeFitLfDown;
   float _SfIterativeFitHfstats1Up;
   float _SfIterativeFitHfstats1Down;
   float _SfIterativeFitHfstats2Up;
   float _SfIterativeFitHfstats2Down;
   float _SfIterativeFitCferr1Up;
   float _SfIterativeFitCferr1Down;
   float _SfIterativeFitCferr2Up;
   float _SfIterativeFitCferr2Down;

   float _r_AK4PF_pt;
   float _r_AK4PchsF_pt;
   float _r_AK8PF_pt;
   float _r_AK8PFchs_pt;
   float _r_AK4PF_phi;
   float _r_AK4PchsF_phi;
   float _r_AK8PF_phi;
   float _r_AK8PFchs_phi;

   ClassDef(Jet,1)
};

#endif

#include "include/NtupleProducer.h"

ClassImp(Jet)
    
Jet::Jet()
{
}

Jet::~Jet()
{
}

void Jet::read()
{
   _ID = idx;
   
   _E   = ntP->jet_E->at(idx);
   _pt  = ntP->jet_pt->at(idx);
   _eta = ntP->jet_eta->at(idx);
   _phi = ntP->jet_phi->at(idx);
   _m   = ntP->jet_m->at(idx);

   if( CHECK(ntP->jet_genJet_E) )                _gen_E   = ntP->jet_genJet_E->at(idx);
   if( CHECK(ntP->jet_genJet_pt) )               _gen_pt  = ntP->jet_genJet_pt->at(idx);
   if( CHECK(ntP->jet_genJet_eta) )              _gen_eta = ntP->jet_genJet_eta->at(idx);   
   if( CHECK(ntP->jet_genJet_phi) )              _gen_phi = ntP->jet_genJet_phi->at(idx);
   if( CHECK(ntP->jet_genJet_m) )                _gen_m   = ntP->jet_genJet_m->at(idx);
   if( CHECK(ntP->jet_genJet_status) )           _gen_status   = ntP->jet_genJet_status->at(idx);
   if( CHECK(ntP->jet_genJet_id) )               _gen_id   = ntP->jet_genJet_id->at(idx);
   
   if( CHECK(ntP->jet_CSVv2) )                _CSVv2   = ntP->jet_CSVv2->at(idx);
   if( CHECK(ntP->jet_cMVAv2) )                _cMVAv2   = ntP->jet_cMVAv2->at(idx);
   
   if( CHECK(ntP->jet_CharmCvsL) )                _CharmCvsL   = ntP->jet_CharmCvsL->at(idx);
   if( CHECK(ntP->jet_CharmCvsB) )                _CharmCvsB   = ntP->jet_CharmCvsB->at(idx);
   
   if( CHECK(ntP->jet_looseJetID) )               _isLooseJetID   = ntP->jet_looseJetID->at(idx);
   if( CHECK(ntP->jet_tightJetID) )               _isTightJetID   = ntP->jet_tightJetID->at(idx);
   
   if( CHECK(ntP->jet_hadronFlavour) )               _hadronFlavour   = ntP->jet_hadronFlavour->at(idx);
   
   if( CHECK(ntP->jet_hasGenJet) )               _hasGenJet   = ntP->jet_hasGenJet->at(idx);
   if( CHECK(ntP->jet_genJet_pt) )               _genJetPt   = ntP->jet_genJet_pt->at(idx);
   if( CHECK(ntP->jet_genJet_eta) )              _genJetEta   = ntP->jet_genJet_eta->at(idx);
   if( CHECK(ntP->jet_genJet_phi) )              _genJetPhi   = ntP->jet_genJet_phi->at(idx);
   if( CHECK(ntP->jet_genJet_E) )                _genJetE   = ntP->jet_genJet_E->at(idx);
   
   if(_pt > 0.) _p4.SetPtEtaPhiE(_pt,_eta,_phi,_E);

   if( CHECK(ntP->jet_r_AK4PF_pt) )                    _r_AK4PF_pt     = ntP->jet_r_AK4PF_pt->at(idx);
   if( CHECK(ntP->jet_r_AK4PchsF_pt) )                 _r_AK4PchsF_pt  = ntP->jet_r_AK4PchsF_pt->at(idx);
   if( CHECK(ntP->jet_r_AK8PF_pt) )                    _r_AK8PF_pt     = ntP->jet_r_AK8PF_pt->at(idx);
   if( CHECK(ntP->jet_r_AK8PFchs_pt) )                 _r_AK8PFchs_pt  = ntP->jet_r_AK8PFchs_pt->at(idx);
   if( CHECK(ntP->jet_r_AK4PF_phi) )                   _r_AK4PF_phi    = ntP->jet_r_AK4PF_phi->at(idx);
   if( CHECK(ntP->jet_r_AK4PchsF_phi) )                _r_AK4PchsF_phi = ntP->jet_r_AK4PchsF_phi->at(idx);
   if( CHECK(ntP->jet_r_AK8PF_phi) )                   _r_AK8PF_phi    = ntP->jet_r_AK8PF_phi->at(idx);
   if( CHECK(ntP->jet_r_AK8PFchs_phi) )                _r_AK8PFchs_phi = ntP->jet_r_AK8PFchs_phi->at(idx);
}

void Jet::init()
{
   _E        = VDEF;
   _pt       = VDEF;
   _eta      = VDEF;
   _phi      = VDEF;
   _m        = VDEF;

   _gen_E        = VDEF;
   _gen_pt       = VDEF;
   _gen_eta      = VDEF;
   _gen_phi      = VDEF;
   _gen_m        = VDEF;
   _gen_status        = VDEF;
   _gen_id        = VDEF;
   
   _CSVv2        = VDEF;
   _cMVAv2        = VDEF;
   
   _CharmCvsL        = VDEF;
   _CharmCvsB        = VDEF;

   _isLoose   = 0;
   _isTight   = 0;
   
   _isBTag   = 0;
   
   _isLooseJetID   = 0;
   _isTightJetID   = 0;
   
   _hadronFlavour   = VDEF;
   
   _hasGenJet = 0;
   _genJetPt = VDEF;
   _genJetEta = VDEF;
   _genJetPhi = VDEF;
   _genJetE = VDEF;
   
   _passElecOverlap   = 0;
   _passMuonOverlap   = 0;
   
   _SfIterativeFitCentral = 1.;
   _SfIterativeFitJesUp = 1.;
   _SfIterativeFitJesDown = 1.;
   _SfIterativeFitLfUp = 1.;
   _SfIterativeFitLfDown = 1.;
   _SfIterativeFitHfstats1Up = 1.;
   _SfIterativeFitHfstats1Down = 1.;
   _SfIterativeFitHfstats2Up = 1.;
   _SfIterativeFitHfstats2Down = 1.;
   _SfIterativeFitCferr1Up = 1.;
   _SfIterativeFitCferr1Down = 1.;
   _SfIterativeFitCferr2Up = 1.;
   _SfIterativeFitCferr2Down = 1.;

   _r_AK4PF_pt     = VDEF;
   _r_AK4PchsF_pt  = VDEF;
   _r_AK8PF_pt     = VDEF;
   _r_AK8PFchs_pt  = VDEF;
   _r_AK4PF_phi    = VDEF;
   _r_AK4PchsF_phi = VDEF;
   _r_AK8PF_phi    = VDEF;
   _r_AK8PFchs_phi = VDEF;
}

void Jet::sel()
{   
   bool passPt = (_pt > 30.);
   bool passEta = (fabs(_eta) < 5.2);

   TLorentzVector vjet;
   vjet.SetPtEtaPhiE(_pt,_eta,_phi,_E);
   
   bool passMuOverlap = 1;
   int nMuonLoose = nt->NtMuonLoose->size();
   for(int im=0;im<nMuonLoose;im++)
     {
	float dr = GetDeltaR(_eta,_phi,nt->NtMuonLoose->at(im).eta(),nt->NtMuonLoose->at(im).phi());
	if( dr < 0.4 && nt->NtMuonLoose->at(im).pt() > 10. ) passMuOverlap = 0;
     }  

   bool passElOverlap = 1;
   int nElectronLoose = nt->NtElectronLoose->size();
   for(int ie=0;ie<nElectronLoose;ie++)
     {
	float dr = GetDeltaR(_eta,_phi,nt->NtElectronLoose->at(ie).eta(),nt->NtElectronLoose->at(ie).phi());
	if( dr < 0.4 && nt->NtElectronLoose->at(ie).pt() > 10. ) passElOverlap = 0;
     }  

   _passElecOverlap = passElOverlap;
   _passMuonOverlap = passMuOverlap;

   if( !_isdata )
     {	
	if( abs(_hadronFlavour) == 5 )
	  {	
	     _SfIterativeFitCentral = reader_iterativefit->eval_auto_bounds("central",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitJesUp = reader_iterativefit->eval_auto_bounds("up_jes",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitJesDown = reader_iterativefit->eval_auto_bounds("down_jes",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitLfUp = reader_iterativefit->eval_auto_bounds("up_lf",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitLfDown = reader_iterativefit->eval_auto_bounds("down_lf",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitHfstats1Up = reader_iterativefit->eval_auto_bounds("up_hfstats1",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitHfstats1Down = reader_iterativefit->eval_auto_bounds("down_hfstats1",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitHfstats2Up = reader_iterativefit->eval_auto_bounds("up_hfstats2",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitHfstats2Down = reader_iterativefit->eval_auto_bounds("down_hfstats2",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitCferr1Up = reader_iterativefit->eval_auto_bounds("up_cferr1",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitCferr1Down = reader_iterativefit->eval_auto_bounds("down_cferr1",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitCferr2Up = reader_iterativefit->eval_auto_bounds("up_cferr2",BTagEntry::FLAV_B,_eta,_pt);
	     _SfIterativeFitCferr2Down = reader_iterativefit->eval_auto_bounds("down_cferr2",BTagEntry::FLAV_B,_eta,_pt);
	  }
	else if( abs(_hadronFlavour) == 4 )
	  {
	     _SfIterativeFitCentral = reader_iterativefit->eval_auto_bounds("central",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitJesUp = reader_iterativefit->eval_auto_bounds("up_jes",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitJesDown = reader_iterativefit->eval_auto_bounds("down_jes",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitLfUp = reader_iterativefit->eval_auto_bounds("up_lf",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitLfDown = reader_iterativefit->eval_auto_bounds("down_lf",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitHfstats1Up = reader_iterativefit->eval_auto_bounds("up_hfstats1",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitHfstats1Down = reader_iterativefit->eval_auto_bounds("down_hfstats1",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitHfstats2Up = reader_iterativefit->eval_auto_bounds("up_hfstats2",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitHfstats2Down = reader_iterativefit->eval_auto_bounds("down_hfstats2",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitCferr1Up = reader_iterativefit->eval_auto_bounds("up_cferr1",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitCferr1Down = reader_iterativefit->eval_auto_bounds("down_cferr1",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitCferr2Up = reader_iterativefit->eval_auto_bounds("up_cferr2",BTagEntry::FLAV_C,_eta,_pt);
	     _SfIterativeFitCferr2Down = reader_iterativefit->eval_auto_bounds("down_cferr2",BTagEntry::FLAV_C,_eta,_pt);
	  }
	else
	  {
	     _SfIterativeFitCentral = reader_iterativefit->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitJesUp = reader_iterativefit->eval_auto_bounds("up_jes",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitJesDown = reader_iterativefit->eval_auto_bounds("down_jes",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitLfUp = reader_iterativefit->eval_auto_bounds("up_lf",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitLfDown = reader_iterativefit->eval_auto_bounds("down_lf",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitHfstats1Up = reader_iterativefit->eval_auto_bounds("up_hfstats1",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitHfstats1Down = reader_iterativefit->eval_auto_bounds("down_hfstats1",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitHfstats2Up = reader_iterativefit->eval_auto_bounds("up_hfstats2",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitHfstats2Down = reader_iterativefit->eval_auto_bounds("down_hfstats2",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitCferr1Up = reader_iterativefit->eval_auto_bounds("up_cferr1",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitCferr1Down = reader_iterativefit->eval_auto_bounds("down_cferr1",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitCferr2Up = reader_iterativefit->eval_auto_bounds("up_cferr2",BTagEntry::FLAV_UDSG,_eta,_pt);
	     _SfIterativeFitCferr2Down = reader_iterativefit->eval_auto_bounds("down_cferr2",BTagEntry::FLAV_UDSG,_eta,_pt);
	  }   
     }
   
   _isLoose = (
	       passPt &&
	       passEta &&
	       passMuOverlap &&
	       passElOverlap
	      );

   _isTight = (
	       _isLoose &&
	       _isLooseJetID
	      );
   
   _isBTag = (_CSVv2 > 0.800);

   for(int id=0;id<evdebug->size();id++)
     {	
	if( nt->NtEvent->at(0).id() == evdebug->at(id) )
	  {
	     std::cout << "Jet #" << _ID << std::endl;
	     std::cout << "   pt=" << _pt << " eta=" << _eta << " phi=" << _phi << std::endl;
	     std::cout << "   isTight=" << _isTight << std::endl;
	     std::cout << "   passElOverlap=" << passElOverlap << std::endl;
	     std::cout << "   passMuOverlap=" << passMuOverlap << std::endl;
	  }	
     }
}

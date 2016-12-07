#ifndef NTUPLE_H
#define NTUPLE_H

#include "Electron.h"
#include "Muon.h"
#include "Event.h"
#include "Jet.h"
#include "Truth.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class Ntuple
{
 public:
   Ntuple();
   virtual ~Ntuple();
   
   std::vector<Electron>* NtElectronLoose;
   std::vector<Electron>* NtElectronTight;
   
   std::vector<Muon>* NtMuonLoose;
   std::vector<Muon>* NtMuonTight;

   std::vector<Jet>* NtJetLoose;
   std::vector<Jet>* NtJetTight;
   std::vector<Jet>* NtBJetTight;
   
   std::vector<Event>* NtEvent;
   std::vector<Truth>* NtTruth;
   
   void Init(std::string fname);
   
   void setBranchAddress();
   void createVar();
   void clearVar();
   void fill();
   void write();
   void close();
   
   TFile*  m_file;
   TTree*  m_tree;
   
 private:

   TChain* m_chain;
};

#endif

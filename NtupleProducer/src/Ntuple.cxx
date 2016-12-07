#include "include/Ntuple.h"

Ntuple::Ntuple()
{
}

Ntuple::~Ntuple()
{   
   m_file->Write();
   m_file->Close();
}

void Ntuple::Init(std::string fname)
{
   std::string foutName = fname+".root";
   m_file = new TFile(foutName.c_str(),"RECREATE");
   
   m_tree = new TTree("Nt"," Ntuple");
}

void Ntuple::setBranchAddress()
{
   m_tree->Branch("Electron","std::vector<Electron>",(NtElectronTight),32000,1);
   m_tree->Branch("Muon","std::vector<Muon>",(NtMuonTight),32000,1);
   m_tree->Branch("Event","std::vector<Event>",(NtEvent),32000,1);
   m_tree->Branch("Jet","std::vector<Jet>",(NtJetTight),32000,1);
   m_tree->Branch("Truth","std::vector<Truth>",(NtTruth),32000,1);
}

void Ntuple::createVar()
{
   NtElectronLoose           = new std::vector<Electron>;
   NtElectronTight           = new std::vector<Electron>;
   
   NtMuonLoose               = new std::vector<Muon>;
   NtMuonTight               = new std::vector<Muon>;

   NtJetLoose               = new std::vector<Jet>;
   NtJetTight               = new std::vector<Jet>;
   NtBJetTight               = new std::vector<Jet>;
   
   NtEvent               = new std::vector<Event>;
   NtTruth               = new std::vector<Truth>;
}

void Ntuple::clearVar()
{
   NtElectronLoose->clear();
   NtElectronTight->clear();
   
   NtMuonLoose->clear();
   NtMuonTight->clear();

   NtJetLoose->clear();
   NtJetTight->clear();
   NtBJetTight->clear();
   
   NtEvent->clear();
   NtTruth->clear();
}

void Ntuple::fill()
{
   m_tree->Fill();
}

void Ntuple::write()
{
   m_file->Write();
}

void Ntuple::close()
{
   m_file->Close();
}


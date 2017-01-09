#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Tree {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ev_run;
   Int_t           ev_id;
   Int_t           ev_lumi;
   Float_t         ev_rho;
   Float_t         met_px;
   Float_t         met_py;
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         met_sumet;
   Double_t        met_sig;
   Double_t        met_cov00;
   Double_t        met_cov10;
   Double_t        met_cov01;
   Double_t        met_cov11;
   Float_t         metGen_px;
   Float_t         metGen_py;
   Float_t         metGen_pt;
   Float_t         metGen_phi;
   Float_t         metGen_sumet;
   Float_t         metGen_MuonEt;
   Float_t         met_shiftedPx_JetEnUp;
   Float_t         met_shiftedPx_JetEnDown;
   Float_t         met_shiftedPx_JetResUp;
   Float_t         met_shiftedPx_JetResDown;
   Float_t         met_shiftedPx_MuonEnUp;
   Float_t         met_shiftedPx_MuonEnDown;
   Float_t         met_shiftedPx_ElectronEnUp;
   Float_t         met_shiftedPx_ElectronEnDown;
   Float_t         met_shiftedPx_TauEnUp;
   Float_t         met_shiftedPx_TauEnDown;
   Float_t         met_shiftedPx_UnclusteredEnUp;
   Float_t         met_shiftedPx_UnclusteredEnDown;
   Float_t         met_shiftedPx_NoShift;
   Float_t         met_shiftedPx_PhotonEnUp;
   Float_t         met_shiftedPx_PhotonEnDown;
   Float_t         met_shiftedPy_JetEnUp;
   Float_t         met_shiftedPy_JetEnDown;
   Float_t         met_shiftedPy_JetResUp;
   Float_t         met_shiftedPy_JetResDown;
   Float_t         met_shiftedPy_MuonEnUp;
   Float_t         met_shiftedPy_MuonEnDown;
   Float_t         met_shiftedPy_ElectronEnUp;
   Float_t         met_shiftedPy_ElectronEnDown;
   Float_t         met_shiftedPy_TauEnUp;
   Float_t         met_shiftedPy_TauEnDown;
   Float_t         met_shiftedPy_UnclusteredEnUp;
   Float_t         met_shiftedPy_UnclusteredEnDown;
   Float_t         met_shiftedPy_NoShift;
   Float_t         met_shiftedPy_PhotonEnUp;
   Float_t         met_shiftedPy_PhotonEnDown;
   Float_t         met_shiftedPhi_JetEnUp;
   Float_t         met_shiftedPhi_JetEnDown;
   Float_t         met_shiftedPhi_JetResUp;
   Float_t         met_shiftedPhi_JetResDown;
   Float_t         met_shiftedPhi_MuonEnUp;
   Float_t         met_shiftedPhi_MuonEnDown;
   Float_t         met_shiftedPhi_ElectronEnUp;
   Float_t         met_shiftedPhi_ElectronEnDown;
   Float_t         met_shiftedPhi_TauEnUp;
   Float_t         met_shiftedPhi_TauEnDown;
   Float_t         met_shiftedPhi_UnclusteredEnUp;
   Float_t         met_shiftedPhi_UnclusteredEnDown;
   Float_t         met_shiftedPhi_NoShift;
   Float_t         met_shiftedPhi_PhotonEnUp;
   Float_t         met_shiftedPhi_PhotonEnDown;
   Float_t         met_shiftedSumEt_JetEnUp;
   Float_t         met_shiftedSumEt_JetEnDown;
   Float_t         met_shiftedSumEt_JetResUp;
   Float_t         met_shiftedSumEt_JetResDown;
   Float_t         met_shiftedSumEt_MuonEnUp;
   Float_t         met_shiftedSumEt_MuonEnDown;
   Float_t         met_shiftedSumEt_ElectronEnUp;
   Float_t         met_shiftedSumEt_ElectronEnDown;
   Float_t         met_shiftedSumEt_TauEnUp;
   Float_t         met_shiftedSumEt_TauEnDown;
   Float_t         met_shiftedSumEt_UnclusteredEnUp;
   Float_t         met_shiftedSumEt_UnclusteredEnDown;
   Float_t         met_shiftedSumEt_NoShift;
   Float_t         met_shiftedSumEt_PhotonEnUp;
   Float_t         met_shiftedSumEt_PhotonEnDown;
   Float_t         metNoHF_pt;
   Float_t         metNoHF_phi;
   Float_t         metNoHF_sumet;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Float_t         pv_xError;
   Float_t         pv_yError;
   Float_t         pv_zError;
   Float_t         pv_chi2;
   Int_t           pv_ndof;
   Float_t         pv_rho;
   Int_t           pv_isFake;
   Float_t         mc_weight;
   Int_t           mc_id;
   Int_t           mc_f1;
   Int_t           mc_f2;
   Float_t         mc_x1;
   Float_t         mc_x2;
   Float_t         mc_scale;
   Float_t         mc_ptHat;
   Float_t         weight_originalXWGTUP;
   Float_t         weight_scale_muF0p5;
   Float_t         weight_scale_muF2;
   Float_t         weight_scale_muR0p5;
   Float_t         weight_scale_muR2;
   vector<float>   *mc_pdfweights;
   vector<string>  *mc_pdfweightIds;
   Int_t           mc_pu_intime_NumInt;
   Int_t           mc_pu_trueNumInt;
   Int_t           mc_pu_before_npu;
   Int_t           mc_pu_after_npu;
   Int_t           mc_pu_Npvi;
   vector<int>     *mc_pu_Nzpositions;
   vector<int>     *mc_pu_BunchCrossing;
   vector<vector<float> > *mc_pu_zpositions;
   vector<vector<float> > *mc_pu_sumpT_lowpT;
   vector<vector<float> > *mc_pu_sumpT_highpT;
   vector<vector<int> > *mc_pu_ntrks_lowpT;
   vector<vector<int> > *mc_pu_ntrks_highpT;
   Int_t           trigger_n;
   vector<int>     *trigger;
   vector<string>  *trigger_name;
   vector<bool>    *trigger_pass;
   vector<int>     *trigger_prescale;
   vector<int>     *trigger_HLTprescale;
   vector<int>     *trigger_L1prescale;
   Int_t           triggerobject_n;
   vector<float>   *triggerobject_pt;
   vector<float>   *triggerobject_eta;
   vector<float>   *triggerobject_phi;
   vector<string>  *triggerobject_collection;
   vector<int>     *triggerobject_filterIds_n;
   vector<int>     *triggerobject_filterIds;
   vector<bool>    *triggerobject_isTriggerL1Mu;
   vector<bool>    *triggerobject_isTriggerL1NoIsoEG;
   vector<bool>    *triggerobject_isTriggerL1IsoEG;
   vector<bool>    *triggerobject_isTriggerL1CenJet;
   vector<bool>    *triggerobject_isTriggerL1ForJet;
   vector<bool>    *triggerobject_isTriggerL1TauJet;
   vector<bool>    *triggerobject_isTriggerL1ETM;
   vector<bool>    *triggerobject_isTriggerL1ETT;
   vector<bool>    *triggerobject_isTriggerL1HTT;
   vector<bool>    *triggerobject_isTriggerL1HTM;
   vector<bool>    *triggerobject_isTriggerL1JetCounts;
   vector<bool>    *triggerobject_isTriggerL1HfBitCounts;
   vector<bool>    *triggerobject_isTriggerL1HfRingEtSums;
   vector<bool>    *triggerobject_isTriggerL1TechTrig;
   vector<bool>    *triggerobject_isTriggerL1Castor;
   vector<bool>    *triggerobject_isTriggerL1BPTX;
   vector<bool>    *triggerobject_isTriggerL1GtExternal;
   vector<bool>    *triggerobject_isHLT_TriggerPhoton;
   vector<bool>    *triggerobject_isHLT_TriggerElectron;
   vector<bool>    *triggerobject_isHLT_TriggerMuon;
   vector<bool>    *triggerobject_isHLT_TriggerTau;
   vector<bool>    *triggerobject_isHLT_TriggerJet;
   vector<bool>    *triggerobject_isHLT_TriggerBJet;
   vector<bool>    *triggerobject_isHLT_TriggerMET;
   vector<bool>    *triggerobject_isHLT_TriggerTET;
   vector<bool>    *triggerobject_isHLT_TriggerTHT;
   vector<bool>    *triggerobject_isHLT_TriggerMHT;
   vector<bool>    *triggerobject_isHLT_TriggerTrack;
   vector<bool>    *triggerobject_isHLT_TriggerCluster;
   vector<bool>    *triggerobject_isHLT_TriggerMETSig;
   vector<bool>    *triggerobject_isHLT_TriggerELongit;
   vector<bool>    *triggerobject_isHLT_TriggerMHTSig;
   vector<bool>    *triggerobject_isHLT_TriggerHLongit;
   vector<int>     *triggerobject_filterLabels_n;
   vector<string>  *triggerobject_filterLabels;
   vector<int>     *triggerobject_pathNamesAll_n;
   vector<string>  *triggerobject_pathNamesAll;
   vector<bool>    *triggerobject_pathNamesAll_isBoth;
   vector<bool>    *triggerobject_pathNamesAll_isL3;
   vector<bool>    *triggerobject_pathNamesAll_isLF;
   vector<bool>    *triggerobject_pathNamesAll_isNone;
   Int_t           nvertex;
   Int_t           el_n;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_m;
   vector<float>   *el_E;
   vector<int>     *el_id;
   vector<int>     *el_charge;
   vector<int>     *el_passConversionVeto;
   vector<int>     *el_isGsfCtfScPixChargeConsistent;
   vector<int>     *el_isGsfScPixChargeConsistent;
   vector<float>   *el_ecalEnergy;
   vector<float>   *el_correctedEcalEnergy;
   vector<float>   *el_trackMomentumError;
   vector<float>   *el_ip3d;
   vector<float>   *el_ip3dErr;
   vector<float>   *el_ip2d;
   vector<float>   *el_ip2dErr;
   vector<float>   *el_neutralHadronIso;
   vector<float>   *el_chargedHadronIso;
   vector<float>   *el_puChargedHadronIso;
   vector<float>   *el_ecalIso;
   vector<float>   *el_hcalIso;
   vector<float>   *el_particleIso;
   vector<float>   *el_photonIso;
   vector<float>   *el_trackIso;
   vector<float>   *el_ecalPFClusterIso;
   vector<float>   *el_hcalPFClusterIso;
   vector<float>   *el_pfIso_sumChargedHadronPt;
   vector<float>   *el_pfIso_sumNeutralHadronEt;
   vector<float>   *el_pfIso_sumPhotonEt;
   vector<float>   *el_pfIso_sumPUPt;
   vector<float>   *el_vx;
   vector<float>   *el_vy;
   vector<float>   *el_vz;
   vector<bool>    *el_hasGsfTrack;
   vector<float>   *el_gsfTrack_d0;
   vector<float>   *el_gsfTrack_z0;
   vector<float>   *el_gsfTrack_d0Error;
   vector<float>   *el_gsfTrack_z0Error;
   vector<float>   *el_gsfTrack_PV_dxy;
   vector<float>   *el_gsfTrack_PV_dz;
   vector<float>   *el_gsfTrack_dxyError;
   vector<float>   *el_gsfTrack_dzError;
   vector<float>   *el_gsfTrack_normalizedChi2;
   vector<float>   *el_superCluster_eta;
   vector<float>   *el_superCluster_phi;
   vector<float>   *el_superCluster_energy;
   vector<float>   *el_superCluster_etaWidth;
   vector<float>   *el_superCluster_phiWidth;
   vector<float>   *el_e5x5;
   vector<float>   *el_sigmaEtaEta;
   vector<float>   *el_sigmaIetaIeta;
   vector<float>   *el_sigmaIphiIphi;
   vector<float>   *el_sigmaIetaIphi;
   vector<float>   *el_full5x5_sigmaIphiIphi;
   vector<float>   *el_full5x5_sigmaEtaEta;
   vector<float>   *el_full5x5_sigmaIetaIeta;
   vector<float>   *el_full5x5_sigmaIetaIphi;
   vector<float>   *el_full5x5_r9;
   vector<float>   *el_full5x5_e1x5;
   vector<float>   *el_full5x5_e5x5;
   vector<float>   *el_full5x5_e2x5Max;
   vector<int>     *el_expectedMissingOuterHits;
   vector<int>     *el_numberOfValidPixelHits;
   vector<int>     *el_numberOfLostPixelHits;
   vector<int>     *el_trackerLayersWithMeasurement;
   vector<int>     *el_pixelLayersWithMeasurement;
   vector<int>     *el_trackerLayersWithoutMeasurement;
   vector<int>     *el_numberOfHits;
   vector<int>     *el_numberOfValidHits;
   vector<float>   *el_hadronicOverEm;
   vector<int>     *el_numberOfLostHits;
   vector<int>     *el_numberOfLostHitsDefault;
   vector<float>   *el_fbrem;
   vector<float>   *el_kf_normalizedChi2;
   vector<float>   *el_gsf_normalizedChi2;
   vector<float>   *el_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *el_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *el_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *el_deltaPhiSeedClusterTrackAtCalo;
   vector<float>   *el_eSuperClusterOverP;
   vector<float>   *el_IoEmIoP;
   vector<float>   *el_ooEmooP;
   vector<float>   *el_eleEoPout;
   vector<float>   *el_mvaNonTrigV0;
   vector<float>   *el_mvaNonTrigCat;
   vector<bool>    *el_vetoCBId;
   vector<bool>    *el_looseCBId;
   vector<bool>    *el_mediumCBId;
   vector<bool>    *el_tightCBId;
   vector<bool>    *el_heepCBId;
   vector<bool>    *el_mediumMVAId;
   vector<bool>    *el_tightMVAId;
   vector<int>     *el_hasMCMatch;
   vector<float>   *el_gen_pt;
   vector<float>   *el_gen_eta;
   vector<float>   *el_gen_phi;
   vector<float>   *el_gen_m;
   vector<int>     *el_gen_status;
   vector<int>     *el_gen_id;
   vector<int>     *el_gen_charge;
   vector<float>   *el_gen_dr;
   vector<int>     *el_hasMCMatchPAT;
   vector<float>   *el_genPAT_pt;
   vector<float>   *el_genPAT_eta;
   vector<float>   *el_genPAT_phi;
   vector<float>   *el_genPAT_m;
   vector<int>     *el_genPAT_status;
   vector<int>     *el_genPAT_id;
   vector<int>     *el_genPAT_charge;
   vector<bool>    *el_hasMatchedConversion;
   vector<int>     *el_expectedMissingInnerHits;
   Int_t           mu_n;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_m;
   vector<float>   *mu_E;
   vector<int>     *mu_id;
   vector<int>     *mu_charge;
   vector<float>   *mu_ip3d;
   vector<float>   *mu_ip3dErr;
   vector<float>   *mu_ip2d;
   vector<float>   *mu_ip2dErr;
   vector<float>   *mu_ip3dBS;
   vector<float>   *mu_ip3dBSErr;
   vector<float>   *mu_neutralHadronIso;
   vector<float>   *mu_chargedHadronIso;
   vector<float>   *mu_puChargedHadronIso;
   vector<float>   *mu_ecalIso;
   vector<float>   *mu_hcalIso;
   vector<float>   *mu_photonIso;
   vector<float>   *mu_trackIso;
   vector<float>   *mu_pfIso03_sumChargedHadronPt;
   vector<float>   *mu_pfIso03_sumChargedParticlePt;
   vector<float>   *mu_pfIso03_sumNeutralHadronEt;
   vector<float>   *mu_pfIso03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIso03_sumPhotonEt;
   vector<float>   *mu_pfIso03_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIso03_sumPUPt;
   vector<float>   *mu_pfIso04_sumChargedHadronPt;
   vector<float>   *mu_pfIso04_sumChargedParticlePt;
   vector<float>   *mu_pfIso04_sumNeutralHadronEt;
   vector<float>   *mu_pfIso04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIso04_sumPhotonEt;
   vector<float>   *mu_pfIso04_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIso04_sumPUPt;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isCaloMuon;
   vector<int>     *mu_isPFMuon;
   vector<int>     *mu_isRPCMuon;
   vector<float>   *mu_vx;
   vector<float>   *mu_vy;
   vector<float>   *mu_vz;
   vector<int>     *mu_numberOfMatches;
   vector<int>     *mu_numberOfMatchedStations;
   vector<float>   *mu_segmentCompatibility;
   vector<float>   *mu_caloCompatibility;
   vector<float>   *mu_combinedQuality_chi2LocalPosition;
   vector<float>   *mu_combinedQuality_trkKink;
   vector<bool>    *mu_isLooseMuon;
   vector<bool>    *mu_isMediumMuon;
   vector<bool>    *mu_isTightMuon;
   vector<bool>    *mu_isSoftMuon;
   vector<bool>    *mu_isHighPtMuon;
   vector<int>     *mu_hasGlobalTrack;
   vector<float>   *mu_globalTrack_d0;
   vector<float>   *mu_globalTrack_z0;
   vector<float>   *mu_globalTrack_d0Error;
   vector<float>   *mu_globalTrack_z0Error;
   vector<float>   *mu_globalTrack_PV_dxy;
   vector<float>   *mu_globalTrack_PV_dz;
   vector<float>   *mu_globalTrack_RP_dxy;
   vector<float>   *mu_globalTrack_RP_dz;
   vector<float>   *mu_globalTrack_dxyError;
   vector<float>   *mu_globalTrack_dzError;
   vector<float>   *mu_globalTrack_normalizedChi2;
   vector<int>     *mu_globalTrack_numberOfValidHits;
   vector<int>     *mu_globalTrack_numberOfValidMuonHits;
   vector<int>     *mu_globalTrack_numberOfLostHits;
   vector<float>   *mu_globalTrack_pt;
   vector<float>   *mu_globalTrack_eta;
   vector<float>   *mu_globalTrack_phi;
   vector<float>   *mu_globalTrack_ptError;
   vector<float>   *mu_globalTrack_etaError;
   vector<float>   *mu_globalTrack_phiError;
   vector<float>   *mu_globalTrack_vx;
   vector<float>   *mu_globalTrack_vy;
   vector<float>   *mu_globalTrack_vz;
   vector<float>   *mu_globalTrack_qoverp;
   vector<float>   *mu_globalTrack_qoverpError;
   vector<int>     *mu_globalTrack_charge;
   vector<int>     *mu_bestTrackType;
   vector<int>     *mu_hasBestTrack;
   vector<float>   *mu_bestTrack_d0;
   vector<float>   *mu_bestTrack_z0;
   vector<float>   *mu_bestTrack_d0Error;
   vector<float>   *mu_bestTrack_z0Error;
   vector<float>   *mu_bestTrack_PV_dxy;
   vector<float>   *mu_bestTrack_PV_dz;
   vector<float>   *mu_bestTrack_dxyError;
   vector<float>   *mu_bestTrack_dzError;
   vector<float>   *mu_bestTrack_normalizedChi2;
   vector<int>     *mu_bestTrack_numberOfValidHits;
   vector<int>     *mu_bestTrack_numberOfLostHits;
   vector<float>   *mu_bestTrack_pt;
   vector<float>   *mu_bestTrack_eta;
   vector<float>   *mu_bestTrack_phi;
   vector<float>   *mu_bestTrack_ptError;
   vector<float>   *mu_bestTrack_etaError;
   vector<float>   *mu_bestTrack_phiError;
   vector<float>   *mu_bestTrack_vx;
   vector<float>   *mu_bestTrack_vy;
   vector<float>   *mu_bestTrack_vz;
   vector<float>   *mu_bestTrack_qoverp;
   vector<float>   *mu_bestTrack_qoverpError;
   vector<int>     *mu_bestTrack_charge;
   vector<int>     *mu_hasInnerTrack;
   vector<float>   *mu_innerTrack_d0;
   vector<float>   *mu_innerTrack_z0;
   vector<float>   *mu_innerTrack_d0Error;
   vector<float>   *mu_innerTrack_z0Error;
   vector<float>   *mu_innerTrack_PV_dxy;
   vector<float>   *mu_innerTrack_PV_dz;
   vector<float>   *mu_innerTrack_dxyError;
   vector<float>   *mu_innerTrack_dzError;
   vector<float>   *mu_innerTrack_normalizedChi2;
   vector<int>     *mu_innerTrack_numberOfValidHits;
   vector<int>     *mu_innerTrack_numberOfLostHits;
   vector<float>   *mu_innerTrack_pt;
   vector<float>   *mu_innerTrack_eta;
   vector<float>   *mu_innerTrack_phi;
   vector<float>   *mu_innerTrack_ptError;
   vector<float>   *mu_innerTrack_etaError;
   vector<float>   *mu_innerTrack_phiError;
   vector<float>   *mu_innerTrack_vx;
   vector<float>   *mu_innerTrack_vy;
   vector<float>   *mu_innerTrack_vz;
   vector<float>   *mu_innerTrack_qoverp;
   vector<float>   *mu_innerTrack_qoverpError;
   vector<int>     *mu_innerTrack_charge;
   vector<int>     *mu_innerTrack_trackerLayersWithMeasurement;
   vector<int>     *mu_innerTrack_numberOfValidPixelHits;
   vector<int>     *mu_innerTrack_numberOfLostPixelHits;
   vector<float>   *mu_innerTrack_validFraction;
   vector<int>     *mu_type;
   vector<int>     *mu_hasMCMatch;
   vector<float>   *mu_gen_pt;
   vector<float>   *mu_gen_eta;
   vector<float>   *mu_gen_phi;
   vector<float>   *mu_gen_m;
   vector<int>     *mu_gen_status;
   vector<int>     *mu_gen_id;
   vector<int>     *mu_gen_charge;
   vector<float>   *mu_gen_dr;
   vector<int>     *mu_hasMCMatchPAT;
   vector<float>   *mu_genPAT_pt;
   vector<float>   *mu_genPAT_eta;
   vector<float>   *mu_genPAT_phi;
   vector<float>   *mu_genPAT_m;
   vector<int>     *mu_genPAT_status;
   vector<int>     *mu_genPAT_id;
   vector<int>     *mu_genPAT_charge;
   Int_t           tau_n;
   vector<float>   *tau_pt;
   vector<float>   *tau_eta;
   vector<float>   *tau_phi;
   vector<float>   *tau_m;
   vector<float>   *tau_E;
   vector<int>     *tau_id;
   vector<int>     *tau_charge;
   vector<bool>    *tau_hasLeadChargedHadrCand;
   vector<float>   *tau_leadingTrackPt;
   vector<int>     *tau_leadingTrackCharge;
   vector<float>   *tau_leadingTrackDz;
   vector<float>   *tau_leadingTrackDxy;
   vector<int>     *tau_decayMode;
   vector<float>   *tau_decayModeFindingOldDMs;
   vector<float>   *tau_decayModeFindingNewDMs;
   vector<float>   *tau_puCorrPtSum;
   vector<float>   *tau_neutralIsoPtSum;
   vector<float>   *tau_chargedIsoPtSum;
   vector<float>   *tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<float>   *tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_againstMuonLoose3;
   vector<float>   *tau_againstMuonTight3;
   vector<float>   *tau_pfEssential_jet_pt;
   vector<float>   *tau_pfEssential_jet_eta;
   vector<float>   *tau_pfEssential_jet_phi;
   vector<float>   *tau_pfEssential_jet_m;
   vector<float>   *tau_pfEssential_jetCorr_pt;
   vector<float>   *tau_pfEssential_jetCorr_eta;
   vector<float>   *tau_pfEssential_jetCorr_phi;
   vector<float>   *tau_pfEssential_jetCorr_m;
   vector<bool>    *tau_pfEssential_hasSV;
   vector<float>   *tau_pfEssential_sv_x;
   vector<float>   *tau_pfEssential_sv_y;
   vector<float>   *tau_pfEssential_sv_z;
   vector<float>   *tau_pfEssential_flightLengthSig;
   vector<float>   *tau_pfEssential_dxy;
   vector<float>   *tau_pfEssential_dxy_error;
   vector<float>   *tau_pfEssential_dxy_Sig;
   Int_t           jet_n;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_m;
   vector<float>   *jet_E;
   vector<int>     *jet_ntrk;
   vector<float>   *jet_JBP;
   vector<float>   *jet_JP;
   vector<float>   *jet_TCHP;
   vector<float>   *jet_TCHE;
   vector<float>   *jet_SSVHE;
   vector<float>   *jet_SSVHP;
   vector<float>   *jet_CMVA;
   vector<float>   *jet_CSVv2;
   vector<float>   *jet_cMVAv2;
   vector<float>   *jet_CharmCvsL;
   vector<float>   *jet_CharmCvsB;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_neutralHadronEnergy;
   vector<float>   *jet_neutralEmEnergy;
   vector<float>   *jet_chargedHadronEnergy;
   vector<float>   *jet_chargedEmEnergy;
   vector<float>   *jet_electronEnergy;
   vector<float>   *jet_muonEnergy;
   vector<float>   *jet_photonEnergy;
   vector<int>     *jet_chargedMultiplicity;
   vector<int>     *jet_neutralMultiplicity;
   vector<float>   *jet_jetArea;
   vector<float>   *jet_jecFactorUncorrected;
   vector<float>   *jet_jecFactorL1FastJet;
   vector<float>   *jet_jecFactorL2Relative;
   vector<float>   *jet_jecFactorL3Absolute;
   vector<float>   *jet_neutralHadronEnergyFraction;
   vector<float>   *jet_neutralEmEnergyFraction;
   vector<float>   *jet_chargedHadronEnergyFraction;
   vector<float>   *jet_muonEnergyFraction;
   vector<float>   *jet_chargedEmEnergyFraction;
   vector<float>   *jet_Unc;
   vector<float>   *jet_pileupJetId;
   vector<bool>    *jet_looseJetID;
   vector<bool>    *jet_tightJetID;
   vector<bool>    *jet_hasGenJet;
   vector<float>   *jet_genJet_pt;
   vector<float>   *jet_genJet_eta;
   vector<float>   *jet_genJet_phi;
   vector<float>   *jet_genJet_m;
   vector<float>   *jet_genJet_E;
   vector<int>     *jet_genJet_status;
   vector<int>     *jet_genJet_id;
   vector<bool>    *jet_hasGenParton;
   vector<float>   *jet_genParton_pt;
   vector<float>   *jet_genParton_eta;
   vector<float>   *jet_genParton_phi;
   vector<float>   *jet_genParton_m;
   vector<float>   *jet_genParton_E;
   vector<int>     *jet_genParton_status;
   vector<int>     *jet_genParton_id;
   vector<float>   *jet_r_AK4PF_pt;
   vector<float>   *jet_r_AK4PchsF_pt;
   vector<float>   *jet_r_AK8PF_pt;
   vector<float>   *jet_r_AK8PFchs_pt;
   vector<float>   *jet_r_AK4PF_phi;
   vector<float>   *jet_r_AK4PchsF_phi;
   vector<float>   *jet_r_AK8PF_phi;
   vector<float>   *jet_r_AK8PFchs_phi;
   Int_t           genJet_n;
   vector<float>   *genJet_pt;
   vector<float>   *genJet_eta;
   vector<float>   *genJet_phi;
   vector<float>   *genJet_m;
   vector<float>   *genJet_E;
   vector<float>   *genJet_emEnergy;
   vector<float>   *genJet_hadEnergy;
   vector<float>   *genJet_invisibleEnergy;
   vector<float>   *genJet_auxiliaryEnergy;
   vector<int>     *genJet_flavour;
   Int_t           gen_n;
   vector<float>   *gen_pt;
   vector<float>   *gen_eta;
   vector<float>   *gen_phi;
   vector<float>   *gen_m;
   vector<float>   *gen_E;
   vector<int>     *gen_status;
   vector<int>     *gen_id;
   vector<int>     *gen_charge;
   vector<int>     *gen_index;
   vector<int>     *gen_mother_index;
   vector<int>     *gen_daughter_n;
   vector<vector<int> > *gen_daughter_index;
   
   // List of branches
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_id;   //!
   TBranch        *b_ev_lumi;   //!
   TBranch        *b_ev_rho;   //!
   TBranch        *b_met_px;   //!
   TBranch        *b_met_py;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_sig;   //!
   TBranch        *b_met_cov00;   //!
   TBranch        *b_met_cov10;   //!
   TBranch        *b_met_cov01;   //!
   TBranch        *b_met_cov11;   //!
   TBranch        *b_metGen_px;   //!
   TBranch        *b_metGen_py;   //!
   TBranch        *b_metGen_pt;   //!
   TBranch        *b_metGen_phi;   //!
   TBranch        *b_metGen_sumet;   //!
   TBranch        *b_metGen_MuonEt;   //!
   TBranch        *b_met_shiftedPx_JetEnUp;   //!
   TBranch        *b_met_shiftedPx_JetEnDown;   //!
   TBranch        *b_met_shiftedPx_JetResUp;   //!
   TBranch        *b_met_shiftedPx_JetResDown;   //!
   TBranch        *b_met_shiftedPx_MuonEnUp;   //!
   TBranch        *b_met_shiftedPx_MuonEnDown;   //!
   TBranch        *b_met_shiftedPx_ElectronEnUp;   //!
   TBranch        *b_met_shiftedPx_ElectronEnDown;   //!
   TBranch        *b_met_shiftedPx_TauEnUp;   //!
   TBranch        *b_met_shiftedPx_TauEnDown;   //!
   TBranch        *b_met_shiftedPx_UnclusteredEnUp;   //!
   TBranch        *b_met_shiftedPx_UnclusteredEnDown;   //!
   TBranch        *b_met_shiftedPx_NoShift;   //!
   TBranch        *b_met_shiftedPx_PhotonEnUp;   //!
   TBranch        *b_met_shiftedPx_PhotonEnDown;   //!
   TBranch        *b_met_shiftedPy_JetEnUp;   //!
   TBranch        *b_met_shiftedPy_JetEnDown;   //!
   TBranch        *b_met_shiftedPy_JetResUp;   //!
   TBranch        *b_met_shiftedPy_JetResDown;   //!
   TBranch        *b_met_shiftedPy_MuonEnUp;   //!
   TBranch        *b_met_shiftedPy_MuonEnDown;   //!
   TBranch        *b_met_shiftedPy_ElectronEnUp;   //!
   TBranch        *b_met_shiftedPy_ElectronEnDown;   //!
   TBranch        *b_met_shiftedPy_TauEnUp;   //!
   TBranch        *b_met_shiftedPy_TauEnDown;   //!
   TBranch        *b_met_shiftedPy_UnclusteredEnUp;   //!
   TBranch        *b_met_shiftedPy_UnclusteredEnDown;   //!
   TBranch        *b_met_shiftedPy_NoShift;   //!
   TBranch        *b_met_shiftedPy_PhotonEnUp;   //!
   TBranch        *b_met_shiftedPy_PhotonEnDown;   //!
   TBranch        *b_met_shiftedPhi_JetEnUp;   //!
   TBranch        *b_met_shiftedPhi_JetEnDown;   //!
   TBranch        *b_met_shiftedPhi_JetResUp;   //!
   TBranch        *b_met_shiftedPhi_JetResDown;   //!
   TBranch        *b_met_shiftedPhi_MuonEnUp;   //!
   TBranch        *b_met_shiftedPhi_MuonEnDown;   //!
   TBranch        *b_met_shiftedPhi_ElectronEnUp;   //!
   TBranch        *b_met_shiftedPhi_ElectronEnDown;   //!
   TBranch        *b_met_shiftedPhi_TauEnUp;   //!
   TBranch        *b_met_shiftedPhi_TauEnDown;   //!
   TBranch        *b_met_shiftedPhi_UnclusteredEnUp;   //!
   TBranch        *b_met_shiftedPhi_UnclusteredEnDown;   //!
   TBranch        *b_met_shiftedPhi_NoShift;   //!
   TBranch        *b_met_shiftedPhi_PhotonEnUp;   //!
   TBranch        *b_met_shiftedPhi_PhotonEnDown;   //!
   TBranch        *b_met_shiftedSumEt_JetEnUp;   //!
   TBranch        *b_met_shiftedSumEt_JetEnDown;   //!
   TBranch        *b_met_shiftedSumEt_JetResUp;   //!
   TBranch        *b_met_shiftedSumEt_JetResDown;   //!
   TBranch        *b_met_shiftedSumEt_MuonEnUp;   //!
   TBranch        *b_met_shiftedSumEt_MuonEnDown;   //!
   TBranch        *b_met_shiftedSumEt_ElectronEnUp;   //!
   TBranch        *b_met_shiftedSumEt_ElectronEnDown;   //!
   TBranch        *b_met_shiftedSumEt_TauEnUp;   //!
   TBranch        *b_met_shiftedSumEt_TauEnDown;   //!
   TBranch        *b_met_shiftedSumEt_UnclusteredEnUp;   //!
   TBranch        *b_met_shiftedSumEt_UnclusteredEnDown;   //!
   TBranch        *b_met_shiftedSumEt_NoShift;   //!
   TBranch        *b_met_shiftedSumEt_PhotonEnUp;   //!
   TBranch        *b_met_shiftedSumEt_PhotonEnDown;   //!
   TBranch        *b_metNoHF_pt;   //!
   TBranch        *b_metNoHF_phi;   //!
   TBranch        *b_metNoHF_sumet;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xError;   //!
   TBranch        *b_pv_yError;   //!
   TBranch        *b_pv_zError;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_rho;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_f1;   //!
   TBranch        *b_mc_f2;   //!
   TBranch        *b_mc_x1;   //!
   TBranch        *b_mc_x2;   //!
   TBranch        *b_mc_scale;   //!
   TBranch        *b_mc_ptHat;   //!
   TBranch        *b_weight_originalXWGTUP;   //!
   TBranch        *b_weight_scale_muF0p5;   //!
   TBranch        *b_weight_scale_muF2;   //!
   TBranch        *b_weight_scale_muR0p5;   //!
   TBranch        *b_weight_scale_muR2;   //!
   TBranch        *b_mc_pdfweights;   //!
   TBranch        *b_mc_pdfweightIds;   //!
   TBranch        *b_mc_pu_intime_NumInt;   //!
   TBranch        *b_mc_pu_trueNumInt;   //!
   TBranch        *b_mc_pu_before_npu;   //!
   TBranch        *b_mc_pu_after_npu;   //!
   TBranch        *b_mc_pu_Npvi;   //!
   TBranch        *b_mc_pu_Nzpositions;   //!
   TBranch        *b_mc_pu_BunchCrossing;   //!
   TBranch        *b_mc_pu_zpositions;   //!
   TBranch        *b_mc_pu_sumpT_lowpT;   //!
   TBranch        *b_mc_pu_sumpT_highpT;   //!
   TBranch        *b_mc_pu_ntrks_lowpT;   //!
   TBranch        *b_mc_pu_ntrks_highpT;   //!
   TBranch        *b_trigger_n;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_pass;   //!
   TBranch        *b_trigger_prescale;   //!
   TBranch        *b_trigger_HLTprescale;   //!
   TBranch        *b_trigger_L1prescale;   //!
   TBranch        *b_triggerobject_n;   //!
   TBranch        *b_triggerobject_pt;   //!
   TBranch        *b_triggerobject_eta;   //!
   TBranch        *b_triggerobject_phi;   //!
   TBranch        *b_triggerobject_collection;   //!
   TBranch        *b_triggerobject_filterIds_n;   //!
   TBranch        *b_triggerobject_filterIds;   //!
   TBranch        *b_triggerobject_isTriggerL1Mu;   //!
   TBranch        *b_triggerobject_isTriggerL1NoIsoEG;   //!
   TBranch        *b_triggerobject_isTriggerL1IsoEG;   //!
   TBranch        *b_triggerobject_isTriggerL1CenJet;   //!
   TBranch        *b_triggerobject_isTriggerL1ForJet;   //!
   TBranch        *b_triggerobject_isTriggerL1TauJet;   //!
   TBranch        *b_triggerobject_isTriggerL1ETM;   //!
   TBranch        *b_triggerobject_isTriggerL1ETT;   //!
   TBranch        *b_triggerobject_isTriggerL1HTT;   //!
   TBranch        *b_triggerobject_isTriggerL1HTM;   //!
   TBranch        *b_triggerobject_isTriggerL1JetCounts;   //!
   TBranch        *b_triggerobject_isTriggerL1HfBitCounts;   //!
   TBranch        *b_triggerobject_isTriggerL1HfRingEtSums;   //!
   TBranch        *b_triggerobject_isTriggerL1TechTrig;   //!
   TBranch        *b_triggerobject_isTriggerL1Castor;   //!
   TBranch        *b_triggerobject_isTriggerL1BPTX;   //!
   TBranch        *b_triggerobject_isTriggerL1GtExternal;   //!
   TBranch        *b_triggerobject_isHLT_TriggerPhoton;   //!
   TBranch        *b_triggerobject_isHLT_TriggerElectron;   //!
   TBranch        *b_triggerobject_isHLT_TriggerMuon;   //!
   TBranch        *b_triggerobject_isHLT_TriggerTau;   //!
   TBranch        *b_triggerobject_isHLT_TriggerJet;   //!
   TBranch        *b_triggerobject_isHLT_TriggerBJet;   //!
   TBranch        *b_triggerobject_isHLT_TriggerMET;   //!
   TBranch        *b_triggerobject_isHLT_TriggerTET;   //!
   TBranch        *b_triggerobject_isHLT_TriggerTHT;   //!
   TBranch        *b_triggerobject_isHLT_TriggerMHT;   //!
   TBranch        *b_triggerobject_isHLT_TriggerTrack;   //!
   TBranch        *b_triggerobject_isHLT_TriggerCluster;   //!
   TBranch        *b_triggerobject_isHLT_TriggerMETSig;   //!
   TBranch        *b_triggerobject_isHLT_TriggerELongit;   //!
   TBranch        *b_triggerobject_isHLT_TriggerMHTSig;   //!
   TBranch        *b_triggerobject_isHLT_TriggerHLongit;   //!
   TBranch        *b_triggerobject_filterLabels_n;   //!
   TBranch        *b_triggerobject_filterLabels;   //!
   TBranch        *b_triggerobject_pathNamesAll_n;   //!
   TBranch        *b_triggerobject_pathNamesAll;   //!
   TBranch        *b_triggerobject_pathNamesAll_isBoth;   //!
   TBranch        *b_triggerobject_pathNamesAll_isL3;   //!
   TBranch        *b_triggerobject_pathNamesAll_isLF;   //!
   TBranch        *b_triggerobject_pathNamesAll_isNone;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_m;   //!
   TBranch        *b_el_E;   //!
   TBranch        *b_el_id;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_passConversionVeto;   //!
   TBranch        *b_el_isGsfCtfScPixChargeConsistent;   //!
   TBranch        *b_el_isGsfScPixChargeConsistent;   //!
   TBranch        *b_el_ecalEnergy;   //!
   TBranch        *b_el_correctedEcalEnergy;   //!
   TBranch        *b_el_trackMomentumError;   //!
   TBranch        *b_el_ip3d;   //!
   TBranch        *b_el_ip3dErr;   //!
   TBranch        *b_el_ip2d;   //!
   TBranch        *b_el_ip2dErr;   //!
   TBranch        *b_el_neutralHadronIso;   //!
   TBranch        *b_el_chargedHadronIso;   //!
   TBranch        *b_el_puChargedHadronIso;   //!
   TBranch        *b_el_ecalIso;   //!
   TBranch        *b_el_hcalIso;   //!
   TBranch        *b_el_particleIso;   //!
   TBranch        *b_el_photonIso;   //!
   TBranch        *b_el_trackIso;   //!
   TBranch        *b_el_ecalPFClusterIso;   //!
   TBranch        *b_el_hcalPFClusterIso;   //!
   TBranch        *b_el_pfIso_sumChargedHadronPt;   //!
   TBranch        *b_el_pfIso_sumNeutralHadronEt;   //!
   TBranch        *b_el_pfIso_sumPhotonEt;   //!
   TBranch        *b_el_pfIso_sumPUPt;   //!
   TBranch        *b_el_miniIso;   //!
   TBranch        *b_el_miniIsoTTH;   //!
   TBranch        *b_el_vx;   //!
   TBranch        *b_el_vy;   //!
   TBranch        *b_el_vz;   //!
   TBranch        *b_el_hasGsfTrack;   //!
   TBranch        *b_el_gsfTrack_d0;   //!
   TBranch        *b_el_gsfTrack_z0;   //!
   TBranch        *b_el_gsfTrack_d0Error;   //!
   TBranch        *b_el_gsfTrack_z0Error;   //!
   TBranch        *b_el_gsfTrack_PV_dxy;   //!
   TBranch        *b_el_gsfTrack_PV_dz;   //!
   TBranch        *b_el_gsfTrack_dxyError;   //!
   TBranch        *b_el_gsfTrack_dzError;   //!
   TBranch        *b_el_gsfTrack_normalizedChi2;   //!
   TBranch        *b_el_superCluster_eta;   //!
   TBranch        *b_el_superCluster_phi;   //!
   TBranch        *b_el_superCluster_energy;   //!
   TBranch        *b_el_superCluster_etaWidth;   //!
   TBranch        *b_el_superCluster_phiWidth;   //!
   TBranch        *b_el_e5x5;   //!
   TBranch        *b_el_sigmaEtaEta;   //!
   TBranch        *b_el_sigmaIetaIeta;   //!
   TBranch        *b_el_sigmaIphiIphi;   //!
   TBranch        *b_el_sigmaIetaIphi;   //!
   TBranch        *b_el_full5x5_sigmaIphiIphi;   //!
   TBranch        *b_el_full5x5_sigmaEtaEta;   //!
   TBranch        *b_el_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_el_full5x5_sigmaIetaIphi;   //!
   TBranch        *b_el_full5x5_r9;   //!
   TBranch        *b_el_full5x5_e1x5;   //!
   TBranch        *b_el_full5x5_e5x5;   //!
   TBranch        *b_el_full5x5_e2x5Max;   //!
   TBranch        *b_el_expectedMissingOuterHits;   //!
   TBranch        *b_el_numberOfValidPixelHits;   //!
   TBranch        *b_el_numberOfLostPixelHits;   //!
   TBranch        *b_el_trackerLayersWithMeasurement;   //!
   TBranch        *b_el_pixelLayersWithMeasurement;   //!
   TBranch        *b_el_trackerLayersWithoutMeasurement;   //!
   TBranch        *b_el_numberOfHits;   //!
   TBranch        *b_el_numberOfValidHits;   //!
   TBranch        *b_el_hadronicOverEm;   //!
   TBranch        *b_el_numberOfLostHits;   //!
   TBranch        *b_el_numberOfLostHitsDefault;   //!
   TBranch        *b_el_fbrem;   //!
   TBranch        *b_el_kf_normalizedChi2;   //!
   TBranch        *b_el_gsf_normalizedChi2;   //!
   TBranch        *b_el_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_el_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_el_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_el_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_el_eSuperClusterOverP;   //!
   TBranch        *b_el_IoEmIoP;   //!
   TBranch        *b_el_ooEmooP;   //!
   TBranch        *b_el_eleEoPout;   //!
   TBranch        *b_el_mvaNonTrigV0;   //!
   TBranch        *b_el_mvaNonTrigCat;   //!
   TBranch        *b_el_vetoCBId;   //!
   TBranch        *b_el_looseCBId;   //!
   TBranch        *b_el_mediumCBId;   //!
   TBranch        *b_el_tightCBId;   //!
   TBranch        *b_el_heepCBId;   //!
   TBranch        *b_el_mediumMVAId;   //!
   TBranch        *b_el_tightMVAId;   //!
   TBranch        *b_el_hasMCMatch;   //!
   TBranch        *b_el_gen_pt;   //!
   TBranch        *b_el_gen_eta;   //!
   TBranch        *b_el_gen_phi;   //!
   TBranch        *b_el_gen_m;   //!
   TBranch        *b_el_gen_status;   //!
   TBranch        *b_el_gen_id;   //!
   TBranch        *b_el_gen_charge;   //!
   TBranch        *b_el_gen_dr;   //!
   TBranch        *b_el_hasMCMatchPAT;   //!
   TBranch        *b_el_genPAT_pt;   //!
   TBranch        *b_el_genPAT_eta;   //!
   TBranch        *b_el_genPAT_phi;   //!
   TBranch        *b_el_genPAT_m;   //!
   TBranch        *b_el_genPAT_status;   //!
   TBranch        *b_el_genPAT_id;   //!
   TBranch        *b_el_genPAT_charge;   //!
   TBranch        *b_el_hasMatchedConversion;   //!
   TBranch        *b_el_expectedMissingInnerHits;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_m;   //!
   TBranch        *b_mu_E;   //!
   TBranch        *b_mu_id;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_ip3d;   //!
   TBranch        *b_mu_ip3dErr;   //!
   TBranch        *b_mu_ip2d;   //!
   TBranch        *b_mu_ip2dErr;   //!
   TBranch        *b_mu_ip3dBS;   //!
   TBranch        *b_mu_ip3dBSErr;   //!
   TBranch        *b_mu_neutralHadronIso;   //!
   TBranch        *b_mu_chargedHadronIso;   //!
   TBranch        *b_mu_puChargedHadronIso;   //!
   TBranch        *b_mu_ecalIso;   //!
   TBranch        *b_mu_hcalIso;   //!
   TBranch        *b_mu_photonIso;   //!
   TBranch        *b_mu_trackIso;   //!
   TBranch        *b_mu_pfIso03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIso03_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIso03_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIso03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIso03_sumPhotonEt;   //!
   TBranch        *b_mu_pfIso03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIso03_sumPUPt;   //!
   TBranch        *b_mu_pfIso04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIso04_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIso04_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIso04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIso04_sumPhotonEt;   //!
   TBranch        *b_mu_pfIso04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIso04_sumPUPt;   //!
   TBranch        *b_mu_miniIso;   //!
   TBranch        *b_mu_miniIsoTTH;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isCaloMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isRPCMuon;   //!
   TBranch        *b_mu_vx;   //!
   TBranch        *b_mu_vy;   //!
   TBranch        *b_mu_vz;   //!
   TBranch        *b_mu_numberOfMatches;   //!
   TBranch        *b_mu_numberOfMatchedStations;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_caloCompatibility;   //!
   TBranch        *b_mu_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_mu_combinedQuality_trkKink;   //!
   TBranch        *b_mu_isLooseMuon;   //!
   TBranch        *b_mu_isMediumMuon;   //!
   TBranch        *b_mu_isTightMuon;   //!
   TBranch        *b_mu_isSoftMuon;   //!
   TBranch        *b_mu_isHighPtMuon;   //!
   TBranch        *b_mu_hasGlobalTrack;   //!
   TBranch        *b_mu_globalTrack_d0;   //!
   TBranch        *b_mu_globalTrack_z0;   //!
   TBranch        *b_mu_globalTrack_d0Error;   //!
   TBranch        *b_mu_globalTrack_z0Error;   //!
   TBranch        *b_mu_globalTrack_PV_dxy;   //!
   TBranch        *b_mu_globalTrack_PV_dz;   //!
   TBranch        *b_mu_globalTrack_RP_dxy;   //!
   TBranch        *b_mu_globalTrack_RP_dz;   //!
   TBranch        *b_mu_globalTrack_dxyError;   //!
   TBranch        *b_mu_globalTrack_dzError;   //!
   TBranch        *b_mu_globalTrack_normalizedChi2;   //!
   TBranch        *b_mu_globalTrack_numberOfValidHits;   //!
   TBranch        *b_mu_globalTrack_numberOfValidMuonHits;   //!
   TBranch        *b_mu_globalTrack_numberOfLostHits;   //!
   TBranch        *b_mu_globalTrack_pt;   //!
   TBranch        *b_mu_globalTrack_eta;   //!
   TBranch        *b_mu_globalTrack_phi;   //!
   TBranch        *b_mu_globalTrack_ptError;   //!
   TBranch        *b_mu_globalTrack_etaError;   //!
   TBranch        *b_mu_globalTrack_phiError;   //!
   TBranch        *b_mu_globalTrack_vx;   //!
   TBranch        *b_mu_globalTrack_vy;   //!
   TBranch        *b_mu_globalTrack_vz;   //!
   TBranch        *b_mu_globalTrack_qoverp;   //!
   TBranch        *b_mu_globalTrack_qoverpError;   //!
   TBranch        *b_mu_globalTrack_charge;   //!
   TBranch        *b_mu_bestTrackType;   //!
   TBranch        *b_mu_hasBestTrack;   //!
   TBranch        *b_mu_bestTrack_d0;   //!
   TBranch        *b_mu_bestTrack_z0;   //!
   TBranch        *b_mu_bestTrack_d0Error;   //!
   TBranch        *b_mu_bestTrack_z0Error;   //!
   TBranch        *b_mu_bestTrack_PV_dxy;   //!
   TBranch        *b_mu_bestTrack_PV_dz;   //!
   TBranch        *b_mu_bestTrack_dxyError;   //!
   TBranch        *b_mu_bestTrack_dzError;   //!
   TBranch        *b_mu_bestTrack_normalizedChi2;   //!
   TBranch        *b_mu_bestTrack_numberOfValidHits;   //!
   TBranch        *b_mu_bestTrack_numberOfLostHits;   //!
   TBranch        *b_mu_bestTrack_pt;   //!
   TBranch        *b_mu_bestTrack_eta;   //!
   TBranch        *b_mu_bestTrack_phi;   //!
   TBranch        *b_mu_bestTrack_ptError;   //!
   TBranch        *b_mu_bestTrack_etaError;   //!
   TBranch        *b_mu_bestTrack_phiError;   //!
   TBranch        *b_mu_bestTrack_vx;   //!
   TBranch        *b_mu_bestTrack_vy;   //!
   TBranch        *b_mu_bestTrack_vz;   //!
   TBranch        *b_mu_bestTrack_qoverp;   //!
   TBranch        *b_mu_bestTrack_qoverpError;   //!
   TBranch        *b_mu_bestTrack_charge;   //!
   TBranch        *b_mu_hasInnerTrack;   //!
   TBranch        *b_mu_innerTrack_d0;   //!
   TBranch        *b_mu_innerTrack_z0;   //!
   TBranch        *b_mu_innerTrack_d0Error;   //!
   TBranch        *b_mu_innerTrack_z0Error;   //!
   TBranch        *b_mu_innerTrack_PV_dxy;   //!
   TBranch        *b_mu_innerTrack_PV_dz;   //!
   TBranch        *b_mu_innerTrack_dxyError;   //!
   TBranch        *b_mu_innerTrack_dzError;   //!
   TBranch        *b_mu_innerTrack_normalizedChi2;   //!
   TBranch        *b_mu_innerTrack_numberOfValidHits;   //!
   TBranch        *b_mu_innerTrack_numberOfLostHits;   //!
   TBranch        *b_mu_innerTrack_pt;   //!
   TBranch        *b_mu_innerTrack_eta;   //!
   TBranch        *b_mu_innerTrack_phi;   //!
   TBranch        *b_mu_innerTrack_ptError;   //!
   TBranch        *b_mu_innerTrack_etaError;   //!
   TBranch        *b_mu_innerTrack_phiError;   //!
   TBranch        *b_mu_innerTrack_vx;   //!
   TBranch        *b_mu_innerTrack_vy;   //!
   TBranch        *b_mu_innerTrack_vz;   //!
   TBranch        *b_mu_innerTrack_qoverp;   //!
   TBranch        *b_mu_innerTrack_qoverpError;   //!
   TBranch        *b_mu_innerTrack_charge;   //!
   TBranch        *b_mu_innerTrack_trackerLayersWithMeasurement;   //!
   TBranch        *b_mu_innerTrack_numberOfValidPixelHits;   //!
   TBranch        *b_mu_innerTrack_numberOfLostPixelHits;   //!
   TBranch        *b_mu_innerTrack_validFraction;   //!
   TBranch        *b_mu_type;   //!
   TBranch        *b_mu_hasMCMatch;   //!
   TBranch        *b_mu_gen_pt;   //!
   TBranch        *b_mu_gen_eta;   //!
   TBranch        *b_mu_gen_phi;   //!
   TBranch        *b_mu_gen_m;   //!
   TBranch        *b_mu_gen_status;   //!
   TBranch        *b_mu_gen_id;   //!
   TBranch        *b_mu_gen_charge;   //!
   TBranch        *b_mu_gen_dr;   //!
   TBranch        *b_mu_hasMCMatchPAT;   //!
   TBranch        *b_mu_genPAT_pt;   //!
   TBranch        *b_mu_genPAT_eta;   //!
   TBranch        *b_mu_genPAT_phi;   //!
   TBranch        *b_mu_genPAT_m;   //!
   TBranch        *b_mu_genPAT_status;   //!
   TBranch        *b_mu_genPAT_id;   //!
   TBranch        *b_mu_genPAT_charge;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_m;   //!
   TBranch        *b_tau_E;   //!
   TBranch        *b_tau_id;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_hasLeadChargedHadrCand;   //!
   TBranch        *b_tau_leadingTrackPt;   //!
   TBranch        *b_tau_leadingTrackCharge;   //!
   TBranch        *b_tau_leadingTrackDz;   //!
   TBranch        *b_tau_leadingTrackDxy;   //!
   TBranch        *b_tau_decayMode;   //!
   TBranch        *b_tau_decayModeFindingOldDMs;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_puCorrPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSum;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_againstMuonLoose3;   //!
   TBranch        *b_tau_againstMuonTight3;   //!
   TBranch        *b_tau_pfEssential_jet_pt;   //!
   TBranch        *b_tau_pfEssential_jet_eta;   //!
   TBranch        *b_tau_pfEssential_jet_phi;   //!
   TBranch        *b_tau_pfEssential_jet_m;   //!
   TBranch        *b_tau_pfEssential_jetCorr_pt;   //!
   TBranch        *b_tau_pfEssential_jetCorr_eta;   //!
   TBranch        *b_tau_pfEssential_jetCorr_phi;   //!
   TBranch        *b_tau_pfEssential_jetCorr_m;   //!
   TBranch        *b_tau_pfEssential_hasSV;   //!
   TBranch        *b_tau_pfEssential_sv_x;   //!
   TBranch        *b_tau_pfEssential_sv_y;   //!
   TBranch        *b_tau_pfEssential_sv_z;   //!
   TBranch        *b_tau_pfEssential_flightLengthSig;   //!
   TBranch        *b_tau_pfEssential_dxy;   //!
   TBranch        *b_tau_pfEssential_dxy_error;   //!
   TBranch        *b_tau_pfEssential_dxy_Sig;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_ntrk;   //!
   TBranch        *b_jet_JBP;   //!
   TBranch        *b_jet_JP;   //!
   TBranch        *b_jet_TCHP;   //!
   TBranch        *b_jet_TCHE;   //!
   TBranch        *b_jet_SSVHE;   //!
   TBranch        *b_jet_SSVHP;   //!
   TBranch        *b_jet_CMVA;   //!
   TBranch        *b_jet_CSVv2;   //!
   TBranch        *b_jet_cMVAv2;   //!
   TBranch        *b_jet_CharmCvsL;   //!
   TBranch        *b_jet_CharmCvsB;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_neutralHadronEnergy;   //!
   TBranch        *b_jet_neutralEmEnergy;   //!
   TBranch        *b_jet_chargedHadronEnergy;   //!
   TBranch        *b_jet_chargedEmEnergy;   //!
   TBranch        *b_jet_electronEnergy;   //!
   TBranch        *b_jet_muonEnergy;   //!
   TBranch        *b_jet_photonEnergy;   //!
   TBranch        *b_jet_chargedMultiplicity;   //!
   TBranch        *b_jet_neutralMultiplicity;   //!
   TBranch        *b_jet_jetArea;   //!
   TBranch        *b_jet_jecFactorUncorrected;   //!
   TBranch        *b_jet_jecFactorL1FastJet;   //!
   TBranch        *b_jet_jecFactorL2Relative;   //!
   TBranch        *b_jet_jecFactorL3Absolute;   //!
   TBranch        *b_jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_jet_muonEnergyFraction;   //!
   TBranch        *b_jet_chargedEmEnergyFraction;   //!
   TBranch        *b_jet_Unc;   //!
   TBranch        *b_jet_pileupJetId;   //!
   TBranch        *b_jet_looseJetID;   //!
   TBranch        *b_jet_tightJetID;   //!
   TBranch        *b_jet_hasGenJet;   //!
   TBranch        *b_jet_genJet_pt;   //!
   TBranch        *b_jet_genJet_eta;   //!
   TBranch        *b_jet_genJet_phi;   //!
   TBranch        *b_jet_genJet_m;   //!
   TBranch        *b_jet_genJet_E;   //!
   TBranch        *b_jet_genJet_status;   //!
   TBranch        *b_jet_genJet_id;   //!
   TBranch        *b_jet_hasGenParton;   //!
   TBranch        *b_jet_genParton_pt;   //!
   TBranch        *b_jet_genParton_eta;   //!
   TBranch        *b_jet_genParton_phi;   //!
   TBranch        *b_jet_genParton_m;   //!
   TBranch        *b_jet_genParton_E;   //!
   TBranch        *b_jet_genParton_status;   //!
   TBranch        *b_jet_genParton_id;   //!
   TBranch        *b_jet_r_AK4PF_pt;   //!
   TBranch        *b_jet_r_AK4PchsF_pt;   //!
   TBranch        *b_jet_r_AK8PF_pt;   //!
   TBranch        *b_jet_r_AK8PFchs_pt;   //!
   TBranch        *b_jet_r_AK4PF_phi;   //!
   TBranch        *b_jet_r_AK4PchsF_phi;   //!
   TBranch        *b_jet_r_AK8PF_phi;   //!
   TBranch        *b_jet_r_AK8PFchs_phi;   //!
   TBranch        *b_genJet_n;   //!
   TBranch        *b_genJet_pt;   //!
   TBranch        *b_genJet_eta;   //!
   TBranch        *b_genJet_phi;   //!
   TBranch        *b_genJet_m;   //!
   TBranch        *b_genJet_E;   //!
   TBranch        *b_genJet_emEnergy;   //!
   TBranch        *b_genJet_hadEnergy;   //!
   TBranch        *b_genJet_invisibleEnergy;   //!
   TBranch        *b_genJet_auxiliaryEnergy;   //!
   TBranch        *b_genJet_flavour;   //!
   TBranch        *b_gen_n;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_m;   //!
   TBranch        *b_gen_E;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_id;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_index;   //!
   TBranch        *b_gen_mother_index;   //!
   TBranch        *b_gen_daughter_n;   //!
   TBranch        *b_gen_daughter_index;   //!
   
   Tree(TChain *tree=0,std::string fname="output.root",std::string treename="FlatTree/tree");
   virtual ~Tree();
   virtual Int_t    GetEntry(Long64_t entry);
//   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *ch);
   virtual void     registerInputBranches(TChain *ch);
};

#endif

#include "PlotStyle.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"

#include <vector>
#include <iostream>

void addbin(TH1F*);

int main()
{
   gROOT->SetBatch(1);
   
   SetPlotStyle();

   gStyle->SetOptFit();

   TFile *fout = new TFile("hist.root","RECREATE");
   
   std::vector<std::string> ttbar;
   std::vector<std::string> stop;
   std::vector<std::string> higgs;
   std::vector<std::string> ttv;
   std::vector<std::string> thq;
   std::vector<std::string> tzq;
   std::vector<std::string> wjets;
   std::vector<std::string> zjets;
   
   std::vector<std::string> sig_hut_stop;
   std::vector<std::string> sig_hct_stop;
   std::vector<std::string> sig_gut;
   std::vector<std::string> sig_gct;

   std::vector<std::string> sig_hut_ttbar;
   std::vector<std::string> sig_hct_ttbar;
   
   std::vector<std::string> sig_hut;
   std::vector<std::string> sig_hct;
   
   ttbar.push_back("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola");
   
   stop.push_back("TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");
   stop.push_back("TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola");
   stop.push_back("T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
   stop.push_back("Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
   stop.push_back("TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");
   stop.push_back("TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola");
   
   higgs.push_back("VBF_HToBB_M-125_13TeV-powheg-pythia6");
   higgs.push_back("WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp");
////   higgs.push_back("ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp");
   
   ttv.push_back("TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola");
   ttv.push_back("TTWJets_Tune4C_13TeV-madgraph-tauola");
   ttv.push_back("TTZJets_Tune4C_13TeV-madgraph-tauola");
   
   thq.push_back("THQtoTBBQ_1L");
   
   tzq.push_back("TZQtoTQQQ_1L");
   
   wjets.push_back("WJetsToLNu_13TeV-madgraph-pythia8-tauola");
   
   zjets.push_back("DYJetsToLL_M-50_13TeV-madgraph-pythia8");
   
   sig_hut_stop.push_back("tHToBB_1L_Kappa_hut");
   sig_hct_stop.push_back("tHToBB_1L_Kappa_hct");
   
   sig_gut.push_back("tHToBB_1L_Kappa_gut");
   sig_gct.push_back("tHToBB_1L_Kappa_gct");
   
   sig_hut_ttbar.push_back("TTtoTHToBB_1L_Kappa_hut");
   sig_hct_ttbar.push_back("TTtoTHToBB_1L_Kappa_hct");

   sig_hut.push_back("tHToBB_1L_Kappa_hut");
   sig_hut.push_back("TTtoTHToBB_1L_Kappa_hut");
   
   sig_hct.push_back("tHToBB_1L_Kappa_hct");
   sig_hct.push_back("TTtoTHToBB_1L_Kappa_hct");   
   
   const int ttbar_n = ttbar.size();
   const int stop_n = stop.size();
   const int ttv_n = ttv.size();
   const int zjets_n = zjets.size();
   const int wjets_n = wjets.size();
   const int higgs_n = higgs.size();
   const int thq_n = thq.size();
   const int tzq_n = tzq.size();
   const int sig_hut_n = sig_hut.size();
   const int sig_hut_stop_n = sig_hut_stop.size();
   const int sig_hut_ttbar_n = sig_hut_ttbar.size();

   float ksig = 1;
   
   const int nvar = 28;
   std::string var[nvar] =
     {
	"H_m",
	"H_pt",
	"H_eta",
	"top_m",
	"top_pt",
	"top_eta",
	"W_m",
	"W_pt",
	"W_eta",
	"l_pt",
	"l_eta",
	"nu_pt",
	"nu_eta",
	"Hb_pt",
	"Hb_eta",
	"topb_pt",
	"topb_eta",
	"Hb1_Hb2_dr",
	"H_top_dr",
	"l_nu_dr",
	"W_topb_dr",
	"Hb1_Hb2_costheta",
	"H_top_costheta",
	"l_nu_costheta",
	"W_topb_costheta",
	"HT",
	"MET",
	"njet"
     };   

   std::string xtit[nvar] =
     {
	"m(H) [GeV]",
	"p_{T} (H) [GeV]",
	"#eta (H)",
	"m(top) [GeV]",
	"p_{T} (top) [GeV]",
	"#eta (top)",
	"m(W) [GeV]",
	"p_{T} (W) [GeV]",
	"#eta (W)",
	"p_{T} (l) [GeV]",
	"#eta (l)",
	"p_{T} (#nu) [GeV]",
	"#eta (#nu)",
	"p_{T} (b_{H}) [GeV]",
	"#eta (b_{H})",
	"p_{T} (b_{top}) [GeV]",
	"#eta (b_{top})",
	"#Delta R (b_{H1},b_{H2})",
	"#Delta R (H,top)",
	"#Delta R (l,#nu)",
	"#Delta R (W,b_{top})",
	"cos #theta (b_{H1},b_{H2})",
	"cos #theta (H,top)",
	"cos #theta (l,#nu)",
	"cos #theta (W,b_{top})",
	"H_{T} [GeV]",
	"E_{T}^{miss} [GeV]",
	"Number of jets"
     };   
   
   const int nchan = 2;
   std::string chan[nchan] = 
     {
	"e",
	"m"
     };   

   const int njet = 1;
   std::string jet[njet] =
     {
	"bjge3"
     };   

   const int nsel = 1;
   std::string sel[nsel] =
     {
	"nosel"
     };   

   std::string hname[1000];
   int varid[1000];
   
   int idx = 0;
   for(int iv=0;iv<nvar;iv++)
     {
	for(int ic=0;ic<nchan;ic++)
	  {
	     for(int ij=0;ij<njet;ij++)
	       {
		  for(int is=0;is<nsel;is++)
		    {
		       hname[idx] = "h_"+var[iv]+"_"+chan[ic]+"_"+jet[ij]+"_"+sel[is]+"_ALL";
		       varid[idx] = iv;
		       idx++;
		    }
	       }
	  }	
     }
   const int hn = idx;
   
   TH1F *h_ttbar[hn];   
   TFile *f_ttbar[ttbar_n];   

   TH1F *h_stop[hn];
   TFile *f_stop[stop_n];   

   TH1F *h_ttv[hn];
   TFile *f_ttv[ttv_n];   
   
   TH1F *h_zjets[hn];
   TFile *f_zjets[zjets_n];   

   TH1F *h_wjets[hn];
   TFile *f_wjets[wjets_n];   

   TH1F *h_higgs[hn];
   TFile *f_higgs[higgs_n];   
   
   TH1F *h_thq[hn];
   TFile *f_thq[thq_n];   

   TH1F *h_tzq[hn];
   TFile *f_tzq[tzq_n];   

   TH1F *h_sig_hut[hn];
   TFile *f_sig_hut[sig_hut_n];

   TH1F *h_sig_hut_stop[hn];
   TFile *f_sig_hut_stop[sig_hut_stop_n];

   TH1F *h_sig_hut_ttbar[hn];
   TFile *f_sig_hut_ttbar[sig_hut_ttbar_n];
   
   for(int i=0;i<ttbar_n;i++)
     {	
	std::string fname = "../histMERGE/"+ttbar[i]+".root";
	f_ttbar[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {	     
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_ttbar[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_ttbar[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_ttbar[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<zjets_n;i++)
     {	
	std::string fname = "../histMERGE/"+zjets[i]+".root";
	f_zjets[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_zjets[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_zjets[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_zjets[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<wjets_n;i++)
     {	
	std::string fname = "../histMERGE/"+wjets[i]+".root";
	f_wjets[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_wjets[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_wjets[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_wjets[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<stop_n;i++)
     {	
	std::string fname = "../histMERGE/"+stop[i]+".root";
	f_stop[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_stop[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_stop[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_stop[ih]->Add(h);
	       }	
	  }	
     }
   
   for(int i=0;i<ttv_n;i++)
     {	
	std::string fname = "../histMERGE/"+ttv[i]+".root";
	f_ttv[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_ttv[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_ttv[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_ttv[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<higgs_n;i++)
     {	
	std::string fname = "../histMERGE/"+higgs[i]+".root";
	f_higgs[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_higgs[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_higgs[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_higgs[ih]->Add(h);
	       }	
	  }	
     }
   
   for(int i=0;i<thq_n;i++)
     {	
	std::string fname = "../histMERGE/"+thq[i]+".root";
	f_thq[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_thq[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_thq[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_thq[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<tzq_n;i++)
     {	
	std::string fname = "../histMERGE/"+tzq[i]+".root";
	f_tzq[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_tzq[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_tzq[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_tzq[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<sig_hut_n;i++)
     {	
	std::string fname = "../histMERGE/"+sig_hut[i]+".root";
	f_sig_hut[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_sig_hut[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_sig_hut[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_sig_hut[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<sig_hut_stop_n;i++)
     {	
	std::string fname = "../histMERGE/"+sig_hut_stop[i]+".root";
	f_sig_hut_stop[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_sig_hut_stop[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_sig_hut_stop[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_sig_hut_stop[ih]->Add(h);
	       }	
	  }	
     }

   for(int i=0;i<sig_hut_ttbar_n;i++)
     {	
	std::string fname = "../histMERGE/"+sig_hut_ttbar[i]+".root";
	f_sig_hut_ttbar[i] = TFile::Open(fname.c_str());
	for(int ih=0;ih<hn;ih++)
	  {
	     std::string hname_clone = hname[ih]+"_clone";
	     TH1F *h = (TH1F*)f_sig_hut_ttbar[i]->Get(hname[ih].c_str());
	     if( i == 0 )
	       {
		  h_sig_hut_ttbar[ih] = (TH1F*)h->Clone(hname_clone.c_str());
	       }
	     else
	       {
		  h_sig_hut_ttbar[ih]->Add(h);
	       }	
	  }	
     }

   fout->cd();
   
   for(int ih=0;ih<hn;ih++)
     {
	TCanvas *c1 = new TCanvas();
	
	addbin(h_ttbar[ih]);
	addbin(h_stop[ih]);
	addbin(h_ttv[ih]);
	addbin(h_higgs[ih]);
	addbin(h_zjets[ih]);
	addbin(h_wjets[ih]);
	addbin(h_thq[ih]);
	addbin(h_tzq[ih]);
	addbin(h_sig_hut[ih]);
	addbin(h_sig_hut_stop[ih]);
	addbin(h_sig_hut_ttbar[ih]);

        h_zjets[ih]->SetMarkerSize(0);
	h_zjets[ih]->SetMarkerColor(94);
	h_zjets[ih]->SetLineColor(94);
	h_zjets[ih]->SetFillColor(94);
	h_zjets[ih]->SetLineStyle(1);
	
	h_ttbar[ih]->SetMarkerSize(0);
	h_ttbar[ih]->SetMarkerColor(kAzure+1);
	h_ttbar[ih]->SetLineColor(kAzure+1);
	h_ttbar[ih]->SetFillColor(kAzure+1);
	h_ttbar[ih]->SetLineStyle(1);

	h_stop[ih]->SetMarkerSize(0);
	h_stop[ih]->SetMarkerColor(kGreen);
	h_stop[ih]->SetLineColor(kGreen);
	h_stop[ih]->SetFillColor(kGreen);
	h_stop[ih]->SetLineStyle(1);
	
	h_ttv[ih]->SetMarkerSize(0);
	h_ttv[ih]->SetMarkerColor(kOrange);
	h_ttv[ih]->SetLineColor(kOrange);
	h_ttv[ih]->SetFillColor(kOrange);
	h_ttv[ih]->SetLineStyle(1);

	h_higgs[ih]->SetMarkerSize(0);
	h_higgs[ih]->SetMarkerColor(kYellow);
	h_higgs[ih]->SetLineColor(kYellow);
	h_higgs[ih]->SetFillColor(kYellow);
	h_higgs[ih]->SetLineStyle(1);
	
	h_wjets[ih]->SetMarkerSize(0);
	h_wjets[ih]->SetMarkerColor(kMagenta-7);
	h_wjets[ih]->SetLineColor(kMagenta-7);
	h_wjets[ih]->SetFillColor(kMagenta-7);
	h_wjets[ih]->SetLineStyle(1);

	h_thq[ih]->SetMarkerSize(0);
	h_thq[ih]->SetMarkerColor(kGreen+1);
	h_thq[ih]->SetLineColor(kGreen+1);
	h_thq[ih]->SetFillColor(kGreen+1);
	h_thq[ih]->SetLineStyle(1);

	h_tzq[ih]->SetMarkerSize(0);
	h_tzq[ih]->SetMarkerColor(kGray);
	h_tzq[ih]->SetLineColor(kGray);
	h_tzq[ih]->SetFillColor(kGray);
	h_tzq[ih]->SetLineStyle(1);
	
	h_sig_hut[ih]->SetMarkerSize(0);
	h_sig_hut[ih]->SetMarkerColor(kRed);
	h_sig_hut[ih]->SetLineColor(kRed);
	h_sig_hut[ih]->SetFillColor(0);
	h_sig_hut[ih]->SetLineWidth(2);

	h_sig_hut_stop[ih]->SetMarkerSize(0);
	h_sig_hut_stop[ih]->SetMarkerColor(kRed+1);
	h_sig_hut_stop[ih]->SetLineColor(kRed+1);
	h_sig_hut_stop[ih]->SetFillColor(0);
	h_sig_hut_stop[ih]->SetLineWidth(2);

	h_sig_hut_ttbar[ih]->SetMarkerSize(0);
	h_sig_hut_ttbar[ih]->SetMarkerColor(kRed-7);
	h_sig_hut_ttbar[ih]->SetLineColor(kRed-7);
	h_sig_hut_ttbar[ih]->SetFillColor(0);
	h_sig_hut_ttbar[ih]->SetLineWidth(2);
	
	float thq_int = h_thq[ih]->Integral(0,h_thq[ih]->GetXaxis()->GetNbins()+1);
	float tzq_int = h_tzq[ih]->Integral(0,h_tzq[ih]->GetXaxis()->GetNbins()+1);
	float ttbar_int = h_ttbar[ih]->Integral(0,h_ttbar[ih]->GetXaxis()->GetNbins()+1);
	float stop_int = h_stop[ih]->Integral(0,h_stop[ih]->GetXaxis()->GetNbins()+1);
	float ttv_int = h_ttv[ih]->Integral(0,h_ttv[ih]->GetXaxis()->GetNbins()+1);
	float zjets_int = h_zjets[ih]->Integral(0,h_zjets[ih]->GetXaxis()->GetNbins()+1);
	float wjets_int = h_wjets[ih]->Integral(0,h_wjets[ih]->GetXaxis()->GetNbins()+1);
	float higgs_int = h_higgs[ih]->Integral(0,h_higgs[ih]->GetXaxis()->GetNbins()+1);
	
	float bg_int = thq_int+tzq_int+ttbar_int+stop_int+ttv_int+zjets_int+wjets_int+higgs_int;
	
	std::string h_thq_name = "thq_"+hname[ih];
	TH1F * h_thq_orig = (TH1F*)h_thq[ih]->Clone(h_thq_name.c_str());
	std::string h_thq_name_norm = "thq_"+hname[ih]+"_norm";
	TH1F * h_thq_norm = (TH1F*)h_thq[ih]->Clone(h_thq_name_norm.c_str());
	if(bg_int > 0) h_thq_norm->Scale(1./bg_int);
	
	std::string h_tzq_name = "tzq_"+hname[ih];
	TH1F * h_tzq_orig = (TH1F*)h_tzq[ih]->Clone(h_tzq_name.c_str());
	std::string h_tzq_name_norm = "tzq_"+hname[ih]+"_norm";
	TH1F * h_tzq_norm = (TH1F*)h_tzq[ih]->Clone(h_tzq_name_norm.c_str());
	if(bg_int > 0) h_tzq_norm->Scale(1./bg_int);
	
	std::string h_ttbar_name = "ttbar_"+hname[ih];
	TH1F * h_ttbar_orig = (TH1F*)h_ttbar[ih]->Clone(h_ttbar_name.c_str());
	std::string h_ttbar_name_norm = "ttbar_"+hname[ih]+"_norm";
	TH1F * h_ttbar_norm = (TH1F*)h_ttbar[ih]->Clone(h_ttbar_name_norm.c_str());
	if(bg_int > 0) h_ttbar_norm->Scale(1./bg_int);
	
	std::string h_stop_name = "stop_"+hname[ih];
	TH1F * h_stop_orig = (TH1F*)h_stop[ih]->Clone(h_stop_name.c_str());
	std::string h_stop_name_norm = "stop_"+hname[ih]+"_norm";
	TH1F * h_stop_norm = (TH1F*)h_stop[ih]->Clone(h_stop_name_norm.c_str());
	if(bg_int > 0) h_stop_norm->Scale(1./bg_int);
	
	std::string h_ttv_name = "ttv_"+hname[ih];
	TH1F * h_ttv_orig = (TH1F*)h_ttv[ih]->Clone(h_ttv_name.c_str());
	std::string h_ttv_name_norm = "ttv_"+hname[ih]+"_norm";
	TH1F * h_ttv_norm = (TH1F*)h_ttv[ih]->Clone(h_ttv_name_norm.c_str());
	if(bg_int > 0) h_ttv_norm->Scale(1./bg_int);
	
	std::string h_zjets_name = "zjets_"+hname[ih];
	TH1F * h_zjets_orig = (TH1F*)h_zjets[ih]->Clone(h_zjets_name.c_str());
	std::string h_zjets_name_norm = "zjets_"+hname[ih]+"_norm";
	TH1F * h_zjets_norm = (TH1F*)h_zjets[ih]->Clone(h_zjets_name_norm.c_str());
	if(bg_int > 0) h_zjets_norm->Scale(1./bg_int);
	
	std::string h_wjets_name = "wjets_"+hname[ih];
	TH1F * h_wjets_orig = (TH1F*)h_wjets[ih]->Clone(h_wjets_name.c_str());
	std::string h_wjets_name_norm = "wjets_"+hname[ih]+"_norm";
	TH1F * h_wjets_norm = (TH1F*)h_wjets[ih]->Clone(h_wjets_name_norm.c_str());
	if(bg_int > 0) h_wjets_norm->Scale(1./bg_int);
	
	std::string h_higgs_name = "higgs_"+hname[ih];
	TH1F * h_higgs_orig = (TH1F*)h_higgs[ih]->Clone(h_higgs_name.c_str());
	std::string h_higgs_name_norm = "higgs_"+hname[ih]+"_norm";
	TH1F * h_higgs_norm = (TH1F*)h_higgs[ih]->Clone(h_higgs_name_norm.c_str());
	if(bg_int > 0) h_higgs_norm->Scale(1./bg_int);

	float sig_hut_stop_int = h_sig_hut_stop[ih]->Integral(0,h_sig_hut_stop[ih]->GetXaxis()->GetNbins()+1);
	float sig_hut_ttbar_int = h_sig_hut_ttbar[ih]->Integral(0,h_sig_hut_ttbar[ih]->GetXaxis()->GetNbins()+1);
	
	std::string h_sig_hut_stop_name = "sig_hut_stop_"+hname[ih];
	TH1F * h_sig_hut_stop_orig = (TH1F*)h_sig_hut_stop[ih]->Clone(h_sig_hut_stop_name.c_str());
	std::string h_sig_hut_stop_name_norm = "sig_hut_stop_"+hname[ih]+"_norm";
	TH1F * h_sig_hut_stop_norm = (TH1F*)h_sig_hut_stop[ih]->Clone(h_sig_hut_stop_name_norm.c_str());
	if(sig_hut_stop_int > 0) h_sig_hut_stop_norm->Scale(1./sig_hut_stop_int);

	std::string h_sig_hut_ttbar_name = "sig_hut_ttbar_"+hname[ih];
	TH1F * h_sig_hut_ttbar_orig = (TH1F*)h_sig_hut_ttbar[ih]->Clone(h_sig_hut_ttbar_name.c_str());
	std::string h_sig_hut_ttbar_name_norm = "sig_hut_ttbar_"+hname[ih]+"_norm";
	TH1F * h_sig_hut_ttbar_norm = (TH1F*)h_sig_hut_ttbar[ih]->Clone(h_sig_hut_ttbar_name_norm.c_str());
	if(sig_hut_ttbar_int > 0) h_sig_hut_ttbar_norm->Scale(1./sig_hut_ttbar_int);

	std::string h_bg_name_merged = "bg_"+hname[ih]+"_merged";
	TH1F *h_bg_merged = (TH1F*)h_thq[ih]->Clone(h_bg_name_merged.c_str());
	h_bg_merged->Add(h_tzq[ih]);
	h_bg_merged->Add(h_ttbar[ih]);
	h_bg_merged->Add(h_stop[ih]);
	h_bg_merged->Add(h_ttv[ih]);
	h_bg_merged->Add(h_zjets[ih]);
	h_bg_merged->Add(h_wjets[ih]);
	h_bg_merged->Add(h_higgs[ih]);
	
	THStack *h_bg_norm = new THStack();
	h_bg_norm->Add(h_thq_norm);
	h_bg_norm->Add(h_tzq_norm);
	h_bg_norm->Add(h_ttbar_norm);
	h_bg_norm->Add(h_stop_norm);
	h_bg_norm->Add(h_ttv_norm);
	h_bg_norm->Add(h_zjets_norm);
	h_bg_norm->Add(h_wjets_norm);
	h_bg_norm->Add(h_higgs_norm);
	
	THStack *h_bg = new THStack();
	h_bg->Add(h_thq[ih]);
	h_bg->Add(h_tzq[ih]);
	h_bg->Add(h_ttbar[ih]);
	h_bg->Add(h_stop[ih]);
	h_bg->Add(h_ttv[ih]);
	h_bg->Add(h_zjets[ih]);
	h_bg->Add(h_wjets[ih]);
	h_bg->Add(h_higgs[ih]);

	int itit = varid[ih];
	std::string xtitle = xtit[itit];
	
	// IL normalised
	h_bg->Draw("hist e1");
	h_bg->GetXaxis()->SetTitle(xtitle.c_str());
	h_bg->GetYaxis()->SetTitle("Number of events");

	h_sig_hut_stop[ih]->Scale(ksig);
	h_sig_hut_stop[ih]->Draw("hist e1 same");

	h_sig_hut_ttbar[ih]->Scale(ksig);
	h_sig_hut_ttbar[ih]->Draw("hist e1 same");

	std::string h_sig_name_merged = "sig_"+hname[ih]+"_merged";
	TH1F *h_sig_merged = (TH1F*)h_sig_hut_stop[ih]->Clone(h_sig_name_merged.c_str());
	h_sig_merged->Add(h_sig_hut_ttbar[ih]);
	
	TLegend *leg = new TLegend(0.70,0.92,0.92,0.60);
	leg->SetFillColor(253);
	leg->SetBorderSize(0);
	leg->AddEntry(h_ttbar[ih],"ttbar","f");
	leg->AddEntry(h_stop[ih],"Single top","f");
	leg->AddEntry(h_zjets[ih],"Z+jets","f");
	leg->AddEntry(h_ttv[ih],"ttbarV","f");
	leg->AddEntry(h_wjets[ih],"W+jets","f");
	leg->AddEntry(h_higgs[ih],"Higgs","f");
	leg->AddEntry(h_thq[ih],"tHq","f");
	leg->AddEntry(h_tzq[ih],"tZq","f");
	leg->AddEntry(h_sig_hut_stop[ih],"FCNC #kappa_{Hut}","lp");
	leg->AddEntry(h_sig_hut_ttbar[ih],"FCNC #kappa_{Hut} (ttbar)","lp");
	leg->Draw();

	c1->SetLogy(1);	
	h_bg->SetMinimum(0.001);
	
	std::string pname = "pics/"+hname[ih]+".eps";
	c1->Print(pname.c_str());
	c1->Clear();

	// unity normalised
	h_bg_norm->Draw("hist e1");
	h_bg_norm->GetXaxis()->SetTitle(xtitle.c_str());
	h_bg_norm->GetYaxis()->SetTitle("");

	h_sig_hut_stop_norm->Draw("hist e1 same");
	h_sig_hut_ttbar_norm->Draw("hist e1 same");

	float bg_norm_max = h_bg_norm->GetMaximum();
	float sig_hut_stop_max = h_sig_hut_stop_norm->GetMaximum();

	float max = std::max(bg_norm_max,sig_hut_stop_max);
	h_bg_norm->SetMaximum(1.2*max);
	
	leg->Draw();

	c1->SetLogy(0);
	
	std::string pname_norm = "pics/"+hname[ih]+"_norm.eps";
	c1->Print(pname_norm.c_str());
	c1->Clear();
	
	delete leg;
	delete c1;
     }

   fout->Write();   
   fout->Close();
   
   for(int i=0;i<ttbar_n;i++) f_ttbar[i]->Close();
   for(int i=0;i<zjets_n;i++) f_zjets[i]->Close();
}

void addbin(TH1F *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   // Set overflow and underflow bins to 0
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

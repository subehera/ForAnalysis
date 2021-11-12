#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TCanvas.h>
#include "TClonesArray.h"
#include <TDirectory.h>
#include <TFile.h>
#include "TH1F.h"
#include <TLatex.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "TRandom.h"
#include <TStyle.h>
#include <TSystem.h>
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include "TROOT.h"
#include "TLine.h"
#include "TF1.h"
#include <iomanip>

const double muonPtCut = 0.0; // For Upsilon candidates


void Accptance()
{
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  TFile *f = new TFile("../coherent_accEff_pPb8TeV.root");  //---coherent process
  //TFile *f = new TFile("../incoherent_accEff_pPb8TeV.root"); //--- incoherent process
  TDirectory *dir = (TDirectory*)f->Get("hionia");
  dir->ls();
  fChain = (TTree*)dir->Get("myTree");
  fChain->ls();
  cout<<fChain->GetEntries()<<"   "<<endl;


  //------------------------                                                                                                                
  TString outfileName = "Acceptance_Coherent_JPsi";  //---coherent process
  //TString outfileName = "Acceptance_Incoherent_JPsi"; //-- Incoherent process
  outfileName        += ".root";
  TFile *fout = new TFile(outfileName,"recreate");
  TTree *tr = new TTree("tr","recreate");

  
  // Declaring variables to be read from the data	
  Float_t        muMiDxy;
  Float_t        muMiDz;
  Int_t          muMiNPxlLayers;
  Int_t          muMiNTrkLayers;
  Float_t        muPlDxy;
  Float_t        muPlDz;
  Int_t          muPlNPxlLayers;
  Int_t          muPlNTrkLayers;
  Float_t        vProb;
  
  Bool_t         muPlGoodMu;
  Bool_t         muMiGoodMu;
  
  Bool_t	 muPlHighPurity;
  Bool_t	 muPlTrkMuArb;
  Bool_t	 muPlTMOneStaTight;
  Bool_t         muMiHighPurity;
  Bool_t         muMiTrkMuArb;
  Bool_t         muMiTMOneStaTight;
  
  Int_t          Centrality;
  Int_t          HLTriggers;

  ULong64_t	  HLTriggers_pp;
  Int_t           Reco_QQ_size;
  Int_t           Reco_QQ_sign[45];    
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  TClonesArray    *Gen_mu_4mom;
  Int_t           Reco_QQ_trig[45];
  ULong64_t       Reco_QQ_trig_pp[45]; 
  Float_t         Reco_QQ_VtxProb[45];    
  Int_t           Reco_QQ_mupl_nPixWMea[45];    
  Int_t           Reco_QQ_mumi_nPixWMea[45];    
  Int_t           Reco_QQ_mupl_nTrkWMea[45];    
  Int_t           Reco_QQ_mumi_nTrkWMea[45];    
  Float_t         Reco_QQ_mupl_dxy[45];    
  Float_t         Reco_QQ_mumi_dxy[45];    
  Float_t         Reco_QQ_mupl_dz[45];    
  Float_t         Reco_QQ_mumi_dz[45];    
  
  Bool_t          Reco_QQ_mupl_isGoodMuon[45];
  Bool_t          Reco_QQ_mumi_isGoodMuon[45];
  
  Bool_t          Reco_QQ_mupl_isHighPurity[45];
  Bool_t          Reco_QQ_mupl_TrkMuArb[45];
  Bool_t          Reco_QQ_mupl_TMOneStaTight[45];
  Bool_t          Reco_QQ_mumi_isHighPurity[45];
  Bool_t          Reco_QQ_mumi_TrkMuArb[45];
  Bool_t          Reco_QQ_mumi_TMOneStaTight[45];
  
  Int_t           Gen_QQ_size;
  Int_t           Gen_mu_size;
  Int_t           Gen_QQ_sign[45];    
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  

  // Setting pointers for branches
  TBranch        	*b_Centrality;    
  TBranch        	*b_HLTriggers;
  TBranch               *b_HLTriggers_pp;    
  TBranch        	*b_Reco_QQ_size;    
  TBranch        	*b_Reco_QQ_sign;    
  TBranch        	*b_Reco_QQ_4mom;    
  TBranch        	*b_Reco_QQ_mupl_4mom;    
  TBranch        	*b_Reco_QQ_mumi_4mom;
  TBranch               *b_Gen_mu_4mom;   //!
  TBranch        	*b_Reco_QQ_trig; 
  TBranch               *b_Reco_QQ_trig_pp;
  TBranch        	*b_Reco_QQ_VtxProb;    
  TBranch        	*b_Reco_QQ_mupl_nPixWMea;    
  TBranch        	*b_Reco_QQ_mumi_nPixWMea;    
  TBranch        	*b_Reco_QQ_mupl_nTrkWMea;    
  TBranch        	*b_Reco_QQ_mumi_nTrkWMea;    
  TBranch        	*b_Reco_QQ_mupl_dxy;    
  TBranch        	*b_Reco_QQ_mumi_dxy;    
  TBranch        	*b_Reco_QQ_mupl_dz;    
  TBranch        	*b_Reco_QQ_mumi_dz;    
  
  TBranch        	*b_Reco_QQ_mupl_isGoodMuon;
  TBranch        	*b_Reco_QQ_mumi_isGoodMuon;
  
  TBranch		*b_Reco_QQ_mupl_isHighPurity;
  TBranch		*b_Reco_QQ_mupl_TrkMuArb;
  TBranch		*b_Reco_QQ_mupl_TMOneStaTight;
  TBranch        	*b_Reco_QQ_mumi_isHighPurity;
  TBranch        	*b_Reco_QQ_mumi_TrkMuArb;
  TBranch        	*b_Reco_QQ_mumi_TMOneStaTight;
  
  TBranch        	*b_Gen_QQ_size;
  TBranch        	*b_Gen_mu_size;  
  TBranch        	*b_Gen_QQ_4mom;    
  TBranch        	*b_Gen_QQ_mupl_4mom;    
  TBranch        	*b_Gen_QQ_mumi_4mom;


  // Function Initializations
  bool PtCut(TLorentzVector* Muon);
  bool IsAccept(TLorentzVector* Muon);
  bool MassCut(TLorentzVector* DiMuon, double LowM, double HighM);
  double PtReweight(TLorentzVector* DiMuon, TF1 *Pt_ReWeights);

  float rapLowRpANeg;
  float rapLowRpAPos;
  float rapHighRpANeg;
  float rapHighRpAPos;
  float lowpTLow;
  float highpTLow;
  float lowpTHigh;
  float highpTHigh;
  
  // Setting mass ranges for Upsilon states
  float       massLow = 2.6;
  float       massHigh = 3.5;
  float       ptReweight = 0.0;
  float       ptReweight_XS = 0.0;
  
  rapLowRpANeg = 0.0;
  rapLowRpAPos = 0.0;
  //rapHighRpANeg = 1.93;
  rapHighRpAPos = 2.4;
  rapLowRpANeg = -2.4;
  rapHighRpANeg = 2.4;
  rapLowRpAPos = -2.4;
  rapHighRpAPos = 2.4;
  lowpTLow = 0.0;
  lowpTHigh = 1.0;
  highpTLow = -.0;
  highpTHigh = 1.0;
  
  // Grabbing data from trees. 
  // Set object pointers and Initialize.
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  Gen_mu_4mom = 0;
  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;


  fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  
  fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  
  fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  
  fChain->SetBranchAddress("Reco_QQ_mupl_isHighPurity", Reco_QQ_mupl_isHighPurity, &b_Reco_QQ_mupl_isHighPurity);
  fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
  fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  fChain->SetBranchAddress("Reco_QQ_mumi_isHighPurity", Reco_QQ_mumi_isHighPurity, &b_Reco_QQ_mumi_isHighPurity);
  fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
  fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);


  fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
  fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

    
  // Declaring histograms
  TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);
  TH1D  *RecoEventsPtRpArapNeg = new TH1D("RecoEventsPtRpArapNeg", "Reconstructed", 10, 0., 1.0);
  TH1D  *RecoEventsRapRpAlowpT = new TH1D("RecoEventsRapRpAlowpT", "Reconstructed", 35, -4., 4.);
  TH1D  *GenEventsPtRpArapPos = new TH1D("GenEventsPtRpArapPos", "Generated", 5, 0., 1.0);
  TH1D  *GenEventsRapRpAlowpT= new TH1D("GenEventsRapRpAlowpT", "Generated", 35, -4., 4.);

  TH1D  *GenEventsPTBeforeCut = new TH1D("GenEventsPTBeforeCut", "Generated-nocut", 5, 0., 1.0);
  TH1D  *GenEventsRapBeforeCut= new TH1D("GenEventsRapBeforeCut", "Generated-nocut", 35, -4., 4.);

  TH2D  *GenEventsAcc = new TH2D("GenEventsAcc", "", 5, 0., 1.0, 35, -4., 4.);
  
  TLorentzVector muon1;
  TLorentzVector muon2;
  TLorentzVector dimuon_gen_4mom;
  TLorentzVector muonvec;
  
  // Get total number of event records in tree
  Long64_t nentries = fChain->GetEntries();
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
      fChain->GetEntry(jentry);
      if(jentry%100000 == 0){
	cout<<"--Processing Event: "<<jentry<<endl;
      }
      
      int nsize = Gen_mu_size;
      for (int iQQ = 0; iQQ < nsize; iQQ++)
	{
	  hCrossCheck->Fill(0);
	  TLorentzVector *g_qq4mom = (TLorentzVector*)Gen_mu_4mom->At(iQQ);
	  TLorentzVector *mu1_gen_4mom = (TLorentzVector*)Gen_mu_4mom->At(0);
	  TLorentzVector *mu2_gen_4mom = (TLorentzVector*)Gen_mu_4mom->At(1);
	  Double_t Pxp=mu1_gen_4mom->Px();
	  Double_t Pyp=mu1_gen_4mom->Py();
	  Double_t Pzp=mu1_gen_4mom->Pz();
	  
	  Double_t Pxm=mu2_gen_4mom->Px();
	  Double_t Pym=mu2_gen_4mom->Py();
	  Double_t Pzm=mu2_gen_4mom->Pz();
	  
	  muon1.SetXYZM(Pxp, Pyp, Pzp,.105);
	  muon2.SetXYZM(Pxm, Pym, Pzm,.105);

	  float etaplus = 0;
	  float	etaminus = 0;
	  etaplus = muon1.Eta();
	  etaminus = muon2.Eta();
	  //cout<<"etaplus: " << etaplus << " etaplus:"<< etaplus << endl;
	  muonvec = muon1 + muon2;

	  float ptGen = 0;
	  float rapGen = 0;
	  float rapGenCM = 0;
	  float Phi = 0;
	  float Eta;
	  float Mass;
	  
	  Mass =  muonvec.M();
	  ptGen =  muonvec.Pt();
	  rapGen = muonvec.Rapidity();
	  Phi = muonvec.Phi();
	  Eta = muonvec.Eta();
	  rapGenCM = (-1.*rapGen)-0.47;
	  
	  bool acceptMu = 0;
	  bool PtCutPass = 0;
	  bool MassCutPass = 0;
	  GenEventsRapBeforeCut->Fill(rapGen);
	  GenEventsPTBeforeCut->Fill(ptGen);
	  
	  if (PtCut(mu1_gen_4mom) && PtCut(mu2_gen_4mom)){ PtCutPass = 1; }
	  if (IsAccept(mu1_gen_4mom) && IsAccept(mu2_gen_4mom)){ acceptMu = 1; }
	  MassCutPass = MassCut(g_qq4mom, massLow, massHigh);
	  if ((fabs(rapGen) < 2.4)&&fabs(etaplus)< 2.4 && fabs(etaminus)< 2.4 && ptGen < 1. ){
	    
	    GenEventsPtRpArapPos->Fill(ptGen);
	    GenEventsRapRpAlowpT->Fill(rapGen);
	    
	    }
	  
	} //---QQ GEN loop
    }//---event loop

  TGraphAsymmErrors *EffPtRpArapNeg = new TGraphAsymmErrors;

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->SetRightMargin(1);
  c2->cd();

  EffPtRpArapNeg->BayesDivide(GenEventsRapRpAlowpT,GenEventsRapBeforeCut);
  EffPtRpArapNeg->SetName("EffPtRpArapNeg");

  EffPtRpArapNeg->SetMarkerSize(1.0);
  EffPtRpArapNeg->SetMarkerColor(kRed);
  EffPtRpArapNeg->SetMarkerStyle(20);

  EffPtRpArapNeg->SetTitle("");
  EffPtRpArapNeg->GetXaxis()->SetTitle("Rapidity");
  EffPtRpArapNeg->GetYaxis()->SetTitle("Acceptance");
  EffPtRpArapNeg->GetYaxis()->SetRangeUser(0,1);
  EffPtRpArapNeg->GetXaxis()->SetRangeUser(-30.0, 30.0);
  EffPtRpArapNeg->GetXaxis()->CenterTitle();
  EffPtRpArapNeg->GetYaxis()->CenterTitle();
  EffPtRpArapNeg->GetXaxis()->SetTitleOffset(0.92);
  EffPtRpArapNeg->GetYaxis()->SetTitleOffset(0.92);
  EffPtRpArapNeg->GetXaxis()->SetLabelSize(0.04);
  EffPtRpArapNeg->GetYaxis()->SetLabelSize(0.04);
  EffPtRpArapNeg->Draw("AP");

  //---------------------PT---------
  TGraphAsymmErrors *EffRapRpAlowpT = new TGraphAsymmErrors;
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  c3->SetRightMargin(1);
  c3->cd();

  EffRapRpAlowpT->BayesDivide(GenEventsPtRpArapPos,GenEventsPTBeforeCut);
  EffRapRpAlowpT->SetName("EffRapRpAlowpT");

  EffRapRpAlowpT->SetMarkerSize(1.0);
  EffRapRpAlowpT->SetMarkerColor(kRed);
  EffRapRpAlowpT->SetMarkerStyle(20);
  EffRapRpAlowpT->GetXaxis()->SetRangeUser(0.0, 1.0);
  EffRapRpAlowpT->GetYaxis()->SetRangeUser(0,1.);
  EffRapRpAlowpT->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  EffRapRpAlowpT->GetYaxis()->SetTitle("Acceptance");
  EffRapRpAlowpT->GetXaxis()->CenterTitle();
  EffRapRpAlowpT->GetYaxis()->CenterTitle();
  EffRapRpAlowpT->GetXaxis()->SetTitleOffset(0.92);
  EffRapRpAlowpT->GetYaxis()->SetTitleOffset(0.92);
  EffRapRpAlowpT->GetXaxis()->SetLabelSize(0.04);
  EffRapRpAlowpT->GetYaxis()->SetLabelSize(0.04);
  EffRapRpAlowpT->Draw("AP");

  //GenEventsAcc->Fill(EffPtRpArapNeg, EffRapRpAlowpT);
  //------------Write the tree
  
  hCrossCheck->Write();
  GenEventsRapBeforeCut->Write();
  GenEventsPTBeforeCut->Write();
  GenEventsPtRpArapPos->Write();
  GenEventsRapRpAlowpT->Write();
  EffPtRpArapNeg->Write();
  EffRapRpAlowpT->Write();
  GenEventsAcc->Write();
}//--E O main fun


// Functions
bool PtCut(TLorentzVector* Muon){
        if (Muon->Pt() < muonPtCut){ return false; }
        else return true;
}

bool IsAccept(TLorentzVector* Muon){
    return ((Muon->Pt() > muonPtCut) && (fabs(Muon->Eta())<2.4));
}


bool MassCut(TLorentzVector* DiMuon, double LowM, double HighM){
        if (DiMuon->M() < LowM){ return false; }
        if (DiMuon->M() > HighM){ return false; }
        return true;
}

double PtReweight(TLorentzVector* DiMuon, TF1 *Pt_ReWeights){
        double pT = (DiMuon->Pt());
        return Pt_ReWeights->Eval(pT);
}


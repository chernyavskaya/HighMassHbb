#include <ctype.h>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm> 
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>
#include <TProfile.h>
#include "/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/preselection_highmass.C"
#include "math.h"


Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t erf2( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]))+ (1.-par[0]);
}
#define SWAP2(A, B) { TLorentzVector t = A; A = B; B = t; }
void SortByEta(std::vector<TLorentzVector> &jets){
  int i, j;
	int n=jets.size();
  for (i = n - 1; i >= 0; i--){
    for (j = 0; j < i; j++){
      if (jets[j].Eta() < jets[j + 1].Eta() ){
        SWAP2( jets[j], jets[j + 1] );
		}
    }
	}
}

typedef std::map<double, int> JetList;
const int njets = 300;

typedef struct {
   Float_t eta[njets];
   Float_t pt[njets];
   Float_t JEC_corr[njets];
   Float_t JEC_corr_up[njets];
   Float_t JEC_corr_down[njets];
   Float_t JER_corr[njets];
   Float_t JER_corr_up[njets];
   Float_t JER_corr_down[njets];
   Float_t phi[njets];
	Float_t mass[njets];
	Float_t btag[njets];
	Int_t nsoft;
	Float_t soft_pt[njets];
	Float_t soft_eta[njets];
	Float_t soft_mass[njets];
	Float_t qgl[njets];
	Int_t nsoft2;
	Int_t nsoft5;
	Int_t nsoft10;
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
	Float_t qgl1_VBF[njets];
	Float_t qgl2_VBF[njets];
} Jets;

using namespace std;


int main(int argc, char* argv[]){


TString file_name = std::string(argv[1]);
TString file_tag = std::string(argv[2]);
TString region = std::string(argv[3]); 
int data = atoi(argv[4]);
int applyTrigWeight = atoi(argv[5]);
TString trigWeight_str = std::string(argv[6]);
TString heppyVersion = std::string(argv[7]);
TString postfix = std::string(argv[8]);
TString output = std::string(argv[9]);


std::map <TString, float> xsec;
xsec["BTagCSV"] = 1.;
xsec["VBFHToBB_M-125"] =  2.20;
xsec["GluGluHToBB_M-125"] = 25.69;
xsec["QCD_HT100to200"] = 2.75E07; 
xsec["QCD_HT200to300"] = 1.74E06 ; 
xsec["QCD_HT300to500"] = 3.67E05; 
xsec["QCD_HT500to700"] = 2.94E04; 
xsec["QCD_HT700to1000"] = 6.52E03;
xsec["QCD_HT1000to1500"] = 1.064E03;
xsec["QCD_HT1500to2000"] =  121.5;
xsec["QCD_HT2000toInf"] = 2.54E01;
xsec["TT"] = 831.76;
xsec["TT_madgraph"] = 831.76;
xsec["TT_powheg"] = 831.76;
xsec["ST_tW"] = 71.7 ;			//inclusive decays
xsec["ST_s-channel"] = 3.36; //leptonic decays
xsec["ST_t-channel_top_4f_inclusiveDecays"] = 136.;
xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 81.;
xsec["DYJetsToQQ"] = 1461.02;
xsec["WJetsToQQ"] = 3539.25 ;
xsec["DYJetsToLL"] = 6025.2;
xsec["WJetsToLNu"]  = 61526.7;
xsec["ttHTobb_M125"]  = .295;
xsec["ttHToNonbb_M125"]  = .2117 ;
xsec["ZH_HToBB"] = 3.093083e-01; //to qq
xsec["WplusH_HToBB"] = 3.057600e-01; 
xsec["WminusH_HToBB"] = 1.940120e-01;
xsec["WW"] = 118.7 ;
xsec["ZZ"] = 16.523 ;
xsec["WZ"] = 47.13;
xsec["bbHSusy120"] = 0.284;
xsec["bbHToBB_yb2"] = 0.3104;
xsec["bbHToBB_ybyt"] = 0.0262;
xsec["VBFHToBB_M-125_amc"] = 2.20;
xsec["GluGluHToBB_M-125_amc"] = 25.69;
xsec["VBFHToBB_M-125_ueps_up"] = 2.20;
xsec["VBFHToBB_M-125_ueps_down"] = 2.20;
xsec["VBFHToBB_M-125_herwig"] = 2.20;
xsec["GluGluHToBB_M-125_ueps_up"] = 25.69;
xsec["GluGluHToBB_M-125_ueps_down"] = 25.69;
xsec["GluGluHToBB_M-125_herwig"] = 25.69;
xsec["VBFSpin0ToBBbar_W_1p0_M_300"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_375"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_450"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_525"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_600"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_675"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_750"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_825"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_900"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_975"] = 1.;
xsec["VBFSpin0ToBBbar_W_1p0_M_1050"] = 1.;

 int counter=0;


int whichTrigWeight;
if ((trigWeight_str.CompareTo("none")==0)||(trigWeight_str.CompareTo("nom")==0)) whichTrigWeight=0;
if (trigWeight_str.CompareTo("up")==0) whichTrigWeight=1;
if (trigWeight_str.CompareTo("down")==0) whichTrigWeight=2;

    
float gen_pos=0; 
float gen_neg=0; 
float gen_pos_weight=0; 
float gen_neg_weight=0; 


	
	Float_t presel=0;
	Float_t presel_vtype[10] = {0,0,0,0,0,0,0,0,0};
	Float_t pos_puweight=0;
	Float_t all_puweight=0.;
	Float_t puweight;
	Float_t PU=1.;
	Float_t genweight;
	Float_t bTagWeight;
	Float_t genweight0;
	float  trigWeight_tree;
	Int_t global_counter = 0;
	Int_t HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200;
	Int_t HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460;
	TFile *file_initial;
	TChain *tree_initial;

//	TString path = "/shome/nchernya/Hbb/skim_trees/v14/";
//	file_initial = TFile::Open(path+file_names[files]+"_v14"+dataset_type[set_type]+"/"+file_names[files]+file_postfix[set_type]+dataset_type[set_type]+".root");

/////////////////qgd//////////////
	TString path ;
//	path = "~/Hbb/slim_trees/v14_slimmed_bdt/";
//	file_initial =TFile::Open(path+file_names[files]+dataset_type[set_type]+".root");
	
//	path = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/Hbb/v14/main_tmva/v14_last/main_mva_v14_";
	//path = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/Hbb/v14/main_tmva/new5jet_float_final/main_mva_v14_";
//	path= "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/Hbb/v14/qgd/";
//	if (set_type==0)	file_initial = TFile::Open(path+file_names[files]+dataset_type[set_type]+".root");
//	if (set_type==1) file_initial = TFile::Open(path+"qgd_"+file_names[files]+dataset_type[set_type]+"_2jets"+".root");
//	file_initial=TFile::Open(path+file_names[files]+dataset_type[set_type]+".root");
//	if ((files==9)||(files==8)) file_initial=TFile::Open("/shome/nchernya/Hbb/slim_trees/SysUnc/JECR/main_mva_v14_"+file_names[files]+dataset_type[set_type]+".root");



	//if (files==18) file_initial=TFile::Open("/shome/nchernya/Hbb/skim_trees/v14/vbf_76"+dataset_type[set_type]+"/"+file_names[files]+"_v14"+dataset_type[set_type]+".root");


//	path= "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/Hbb/v14/qgd/";
//	if (set_type==0)	file_initial = TFile::Open(path+file_names[files]+dataset_type[set_type]+".root");
//	if (set_type==1) file_initial = TFile::Open(path+"qgd_"+file_names[files]+dataset_type[set_type]+"_2jets"+".root");
////////////////////////////

	file_initial = TFile::Open(file_name);


	
	tree_initial = (TChain*)file_initial->Get("tree");
	Int_t events_generated;
	TH1F *countPos;
	TH1F *countNeg;
	TH1F *countLHEScale;
	TH1F *countLHEPdf;
	TH1F *countWeighted;
	if ((data!=1)){
		countPos = (TH1F*)file_initial->Get("CountPosWeight");
		countPos = (TH1F*)file_initial->Get("CountPosWeight");
 		countWeighted = (TH1F*)file_initial->Get("CountWeighted");
 		countLHEScale = (TH1F*)file_initial->Get("CountWeightedLHEWeightScale");
		countLHEPdf=	(TH1F*)file_initial->Get("CountWeightedLHEWeightPdf");
 		events_generated = countWeighted->GetEntries();
	} else events_generated = 1;
    Jets Jet;
    Float_t v_type;
    Float_t wrong_type=0.;
    Int_t nJets;
	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	

	Float_t bdt;
	Float_t met_pt;
	Float_t met_phi;

	Jets GenHiggsSisters;

	int pos_weight_presel=0;
 	Int_t selLeptons_tightId[20];
	Float_t selLeptons_relIso03[20] , selLeptons_chargedHadRelIso03[20], selLeptons_pfRelIso03[20];
	Float_t vLeptons_dz[20], vLeptons_edz[20];


Float_t LHE_weights_pdf_wgt[103];
Float_t LHE_weights_scale_wgt[10];

    tree_initial->SetBranchAddress("Vtype",&v_type);
    tree_initial->SetBranchAddress("rho",&rho);
    tree_initial->SetBranchAddress("nJet",&nJets);
    tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
    tree_initial->SetBranchAddress("Jet_corr_JECUp",Jet.JEC_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JECDown",Jet.JEC_corr_down);
    tree_initial->SetBranchAddress("Jet_corr",Jet.JEC_corr);
    tree_initial->SetBranchAddress("Jet_corr_JERUp",Jet.JER_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JERDown",Jet.JER_corr_down);
    tree_initial->SetBranchAddress("Jet_corr_JER",Jet.JER_corr);
    tree_initial->SetBranchAddress("Jet_eta",Jet.eta);
    tree_initial->SetBranchAddress("Jet_phi",Jet.phi);
	tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
	tree_initial->SetBranchAddress("Jet_btagCSV",Jet.btag);
	tree_initial->SetBranchAddress("Jet_blike_VBF",Jet.blike_VBF);
	tree_initial->SetBranchAddress("Jet_id",Jet.id);	
	tree_initial->SetBranchAddress("Jet_puId",Jet.puId);
 	tree_initial->SetBranchAddress("Jet_leadTrackPt",Jet.leadTrackPt);
	tree_initial->SetBranchAddress("met_pt",&met_pt);
	tree_initial->SetBranchAddress("met_phi",&met_phi);
	
	tree_initial->SetBranchAddress("softActivityJets_pt",Jet.soft_pt);
	tree_initial->SetBranchAddress("softActivityJets_eta",Jet.soft_eta);
	tree_initial->SetBranchAddress("softActivityJets_mass",Jet.soft_mass);
	tree_initial->SetBranchAddress("softActivity_HT",&Jet.HTsoft);
	tree_initial->SetBranchAddress("softActivity_njets2",&Jet.nsoft2);
	tree_initial->SetBranchAddress("softActivity_njets5",&Jet.nsoft5);
	tree_initial->SetBranchAddress("softActivity_njets10",&Jet.nsoft10);
	tree_initial->SetBranchAddress("Jet_qgl",Jet.qgl);
	tree_initial->SetBranchAddress("genWeight",&genweight);
	tree_initial->SetBranchAddress("bTagWeight",&bTagWeight);
	tree_initial->SetBranchAddress("puWeight",&puweight);
	tree_initial->SetBranchAddress("nPVs",&nPVs);
	tree_initial->SetBranchAddress("Jet_ptd",Jet.ptd);
	tree_initial->SetBranchAddress("Jet_axis2",Jet.axis2);
	tree_initial->SetBranchAddress("Jet_mult",Jet.mult);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",&HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",&HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460);
	tree_initial->SetBranchAddress("Jet_pt_regVBF",Jet.pt_regVBF);
	tree_initial->SetBranchAddress("json",&JSON);
    
	tree_initial->SetBranchAddress("GenHiggsSisters_pt",GenHiggsSisters.pt);
    tree_initial->SetBranchAddress("GenHiggsSisters_eta",GenHiggsSisters.eta);
    tree_initial->SetBranchAddress("GenHiggsSisters_phi",GenHiggsSisters.phi);
	tree_initial->SetBranchAddress("GenHiggsSisters_mass",GenHiggsSisters.mass);

	tree_initial->SetBranchAddress("selLeptons_tightId",selLeptons_tightId);
	tree_initial->SetBranchAddress("selLeptons_relIso03",selLeptons_relIso03);
	tree_initial->SetBranchAddress("selLeptons_chargedHadRelIso03",selLeptons_chargedHadRelIso03);
	tree_initial->SetBranchAddress("selLeptons_pfRelIso03",selLeptons_pfRelIso03);
	tree_initial->SetBranchAddress("vLeptons_dz",vLeptons_dz);
	tree_initial->SetBranchAddress("vLeptons_edz",vLeptons_edz);

	tree_initial->SetBranchAddress("Jet_qgl1_VBF",Jet.qgl1_VBF);
	tree_initial->SetBranchAddress("Jet_qgl2_VBF",Jet.qgl2_VBF);

	tree_initial->SetBranchAddress("BDT_VBF",&bdt);
	tree_initial->SetBranchAddress("trigWeight",&trigWeight_tree);
	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

 	
    TH1F *hJet1_pt_bin = new TH1F("hJet1_pt_bin", "", 50, 90., 140.);
    hJet1_pt_bin->GetXaxis()->SetTitle("1^{st} Jet p_{T} (GeV)");
    TH1F *hJet2_pt_bin = new TH1F("hJet2_pt_bin", "", 40, 70., 110.);
    hJet2_pt_bin->GetXaxis()->SetTitle("2^{nd} Jet p_{T} (GeV)");
    TH1F *hJet3_pt_bin = new TH1F("hJet3_pt_bin", "", 30, 60., 90.);
    hJet3_pt_bin->GetXaxis()->SetTitle("3^{rd} Jet p_{T} (GeV)");
    TH1F *hJet4_pt_bin = new TH1F("hJet4_pt_bin", "", 40, 30., 70.);
    hJet4_pt_bin->GetXaxis()->SetTitle("4^{th} Jet p_{T} (GeV)");
	TH1F *hMqq_bin = new TH1F("hMqq_bin","",50.,200.,250.);
	hMqq_bin->GetXaxis()->SetTitle("m_{qq} (GeV)");
   
    TH1F *hJet1_pt = new TH1F("hJet1_pt", "", 30, 0., 600.);
    hJet1_pt->GetXaxis()->SetTitle("1^{st} Jet p_{T} (GeV)");
    TH1F *hJet2_pt = new TH1F("hJet2_pt", "", 30, 0., 600.);
    hJet2_pt->GetXaxis()->SetTitle("2^{nd} Jet p_{T} (GeV)");
    TH1F *hJet3_pt = new TH1F("hJet3_pt", "", 20, 0., 400.);
    hJet3_pt->GetXaxis()->SetTitle("3^{rd} Jet p_{T} (GeV)");
    TH1F *hJet4_pt = new TH1F("hJet4_pt", "", 20, 0., 400.);
    hJet4_pt->GetXaxis()->SetTitle("4^{th} Jet p_{T} (GeV)");
    
    TH1F *hJet5_pt = new TH1F("hJet5_pt", "", 25, 0., 500.);
    hJet5_pt->GetXaxis()->SetTitle("5^{th} Jet p_{T} (GeV)");
    
	TH1F *hJet1_eta = new TH1F("hJet1_eta", "", 80, -5., 5.);
    hJet1_eta->GetXaxis()->SetTitle("1^{st} Jet #eta");
    TH1F *hJet2_eta = new TH1F("hJet2_eta", "", 80, -5., 5.);
    hJet2_eta->GetXaxis()->SetTitle("2^{nd} Jet #eta");
    TH1F *hJet3_eta = new TH1F("hJet3_eta", "", 80, -5., 5.);
    hJet3_eta->GetXaxis()->SetTitle("3^{rd} Jet #eta");
    TH1F *hJet4_eta = new TH1F("hJet4_eta", "", 80, -5., 5.);
    hJet4_eta->GetXaxis()->SetTitle("4^{th} Jet #eta");
    
    TH1F *hJet1_phi = new TH1F("hJet1_phi", "", 32,-3.2,3.2);
    hJet1_phi->GetXaxis()->SetTitle("1^{st} Jet #phi");
    TH1F *hJet2_phi = new TH1F("hJet2_phi", "", 32,-3.2,3.2);
    hJet2_phi->GetXaxis()->SetTitle("2^{nd} Jet #phi");
    TH1F *hJet3_phi = new TH1F("hJet3_phi", "", 32,-3.2,3.2);
    hJet3_phi->GetXaxis()->SetTitle("3^{rd} Jet #phi");
    TH1F *hJet4_phi = new TH1F("hJet4_phi", "", 32,-3.2,3.2);
    hJet4_phi->GetXaxis()->SetTitle("4^{th} Jet #phi");
    TH1F *hVtype = new TH1F("hVtype","", 6,-1.,6.);
    hVtype->GetXaxis()->SetTitle("vtype");

	TH1F *hMqq = new TH1F("hMqq","",100.,0.,3000.);
	hMqq->GetXaxis()->SetTitle("m_{qq} (GeV)");
	TH1F *hMqq_pt = new TH1F("hMqq_pt","",100.,0.,3000.);
	hMqq_pt->GetXaxis()->SetTitle("m_{qq}, chosen by p_{T} (GeV)");
   
	TH1F *hEtaQQ = new TH1F("hEtaQQ","",90,0.,9.);
	hEtaQQ->GetXaxis()->SetTitle("|#Delta#eta_{qq}|");

	TH1F *hPhiQQ = new TH1F("hPhiQQ","",32,0.,3.2);
	hPhiQQ->GetXaxis()->SetTitle("|#Delta#phi_{qq}|");
    
	TH1F *hPhiBB = new TH1F("hPhiBB","",32,0.,3.2);
	hPhiBB->GetXaxis()->SetTitle("|#Delta#phi_{bb}|");
	
	TH1F *hEtaSoftJets = new TH1F("hEtaSoftJets","",12,-3.,3.);
	hEtaSoftJets->GetXaxis()->SetTitle("|#eta^{soft}|");
    
	
	TH1F *hMassSoftJets = new TH1F("hMassSoftJets","",10,0.,100.);
	hMassSoftJets->GetXaxis()->SetTitle("m^{soft}");

	TH1F *hHTsoft = new TH1F("hHTsoft","",60,0.,600.);
	hHTsoft->GetXaxis()->SetTitle("H_{T}^{soft} (GeV)" );
	TH1F *hSoft_n2 = new TH1F("hSoft_n2","",25,0.,25.);
	hSoft_n2->GetXaxis()->SetTitle("N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5 = new TH1F("hSoft_n5","",10,0.,10.);
	hSoft_n5->GetXaxis()->SetTitle("N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10 = new TH1F("hSoft_n10","",6,0.,6.);
	hSoft_n10->GetXaxis()->SetTitle("N soft jets, p_{T} > 10 GeV");
	
	TH1F* hMbb = new TH1F("hMbb","",150,0.,1500);
	hMbb->GetXaxis()->SetTitle("m_{bb} (GeV)");
	TH1F* hqgl = new TH1F("hqgl","",20.,0.,1.);
	hqgl->GetXaxis()->SetTitle("QGL 1^{st} q-jet");
	
	TH1F* hqgl2 = new TH1F("hqgl2","",20.,0.,1.);
	hqgl2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");

	TH1F *hbtag = new TH1F("hbtag","",110,0.,1.1);
	hbtag->GetXaxis()->SetTitle("CSV 1^{st} b-jet");
 
	TH1F *hbtag2 = new TH1F("hbtag2","",110,0.,1.1);
	hbtag2->GetXaxis()->SetTitle("CSV 2^{nd} b-jet");
	
	TH1F *hPtSoftJets = new TH1F("hPtSoftJets","",30,0.,300);
	hPtSoftJets->GetXaxis()->SetTitle("p_{T}^{soft} (GeV)");
    TH1F *hPtSoftJets2 = new TH1F("hPtSoftJets2", "", 20, 0., 200.);
    hPtSoftJets2->GetXaxis()->SetTitle("2nd Soft Jet p_{T} (GeV)");
    TH1F *hPtSoftJets3 = new TH1F("hPtSoftJets3", "", 20, 0., 200.);
    hPtSoftJets3->GetXaxis()->SetTitle("3rd Soft Jet p_{T} (GeV)");
	
	TH1F *hcosOqqbb = new TH1F("hcosOqqbb","",100,-1.,1.);
	hcosOqqbb->GetXaxis()->SetTitle("cos(#theta_{bb_qq})");
	TH1F *hEtaQB1 = new TH1F("hEtaQB1","",160.,-8.,8.);
	hEtaQB1->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}");
	TH1F *hEtaQB2 = new TH1F("hEtaQB2","",160.,-8.,8.);
	hEtaQB2->GetXaxis()->SetTitle("#Delta#eta_{qb}^{backward}");
	TH1F *hPhiQB1 = new TH1F("hPhiQB1","",32,0.,3.2);
	hPhiQB1->GetXaxis()->SetTitle("#Delta#phi_{qb}^{forward}");
	TH1F *hPhiQB2 = new TH1F("hPhiQB2","",32,0.,3.2);
	hPhiQB2->GetXaxis()->SetTitle("#Delta#phi_{qb}^{backward}");
	TH1F *hx1 = new TH1F("hx1","",100.,0.,1.);
	hx1->GetXaxis()->SetTitle("x_{1}");
	TH1F *hx2 = new TH1F("hx2","",100.,0.,1.);
	hx2->GetXaxis()->SetTitle("x_{2}");
	TH1F *hVB1_mass = new TH1F("hVB1_mass","",100,0.,1000.);
	hVB1_mass->GetXaxis()->SetTitle("M_{W'_{1}} (GeV)");
	TH1F *hVB2_mass = new TH1F("hVB2_mass","",100.,0.,1000.);
	hVB2_mass->GetXaxis()->SetTitle("M_{W'_{2}} (GeV)");

	TH1F* hEtot = new TH1F("hEtot","",150.,0.,6000.);
	hEtot->GetXaxis()->SetTitle("E^{tot} (GeV)");
	TH1F* hPxtot= new TH1F("hPxtot","",100,-500.,500.);
	hPxtot->GetXaxis()->SetTitle("P_{x}^{tot} (GeV)");
	TH1F* hPytot= new TH1F("hPytot","",100,-500.,500.);
	hPytot->GetXaxis()->SetTitle("P_{y}^{tot} (GeV)");
	TH1F* hPztot= new TH1F("hPztot","",100,-5000.,5000);
	hPztot->GetXaxis()->SetTitle("P_{z}^{tot} (GeV)");

	
	TH1F *hPtqqbb = new TH1F("hPtqqbb","",50.,0.,500.);
	hPtqqbb->GetXaxis()->SetTitle("p_{T} of qqbb system (GeV)");
	TH1F *hPhiqqbb = new TH1F("hPhiqqbb","",32,-3.2,3.2);
	hPhiqqbb->GetXaxis()->SetTitle("-#phi of qqbb system");
	TH1F *hEtaqqbb = new TH1F("hEtaqqbb","",160,0,8);
	hEtaqqbb->GetXaxis()->SetTitle("#eta of qqbb system");

	TH1F *hnPVs = new TH1F("hPVs","",50,0,50);
	hnPVs->GetXaxis()->SetTitle("nPVs");
	
	TH1F* hMbb_fsr = new TH1F("hMbb_fsr","",150,0.,1500.);
	hMbb_fsr->GetXaxis()->SetTitle("m_{bb} [FSR] (GeV)");

	TH1F* hMbb_met_fsr = new TH1F("hMbb_met_fsr","",220,0.,2200.);
	hMbb_met_fsr->GetXaxis()->SetTitle("m_{bb} [FSR + MET] (GeV)");
	TH1F* hMbb_met_fsr_bg = new TH1F("hMbb_met_fsr_bg","",220,0.,2200.);
	hMbb_met_fsr_bg->GetXaxis()->SetTitle("m_{bb} [FSR + MET] (GeV), -0.6 < BDT < 0.6");
	TH1F* hMbb_met_fsr_sg = new TH1F("hMbb_met_fsr_sg","",220,0.,2200.);
	hMbb_met_fsr_sg->GetXaxis()->SetTitle("m_{bb} [FSR + MET], (GeV), BDT > 0.8");
	
	TH1F* hMbb_met_fsr_dphibb = new TH1F("hMbb_met_fsr_dphibb","",150,0.,1500.);
	hMbb_met_fsr_dphibb->GetXaxis()->SetTitle("m_{bb} [FSR + MET], #Delta#phi_{bb} < 2.4 (GeV)");
	
	TH1F* hMbb_met_fsr_long = new TH1F("hMbb_met_fsr_long","",260,0.,13000.);
	hMbb_met_fsr_long->GetXaxis()->SetTitle("m_{bb} [FSR + MET] (GeV)");
	

	TH1F *hJet1q_pt = new TH1F("hJet1q_pt","",17,30,200.);
	hJet1q_pt->GetXaxis()->SetTitle("p_{T} 1^{st} q-jet");
	TH1F *hJet1q_eta = new TH1F("hJet1q_eta","",20,-5,5);
	hJet1q_eta->GetXaxis()->SetTitle("#eta 1^{st} q-jet");
	TH1F *hJet1q_ptd = new TH1F("hJet1q_ptd","",100,0,1);
	hJet1q_ptd->GetXaxis()->SetTitle("ptd 1^{st} q-jet");
	TH1F *hJet1q_axis2= new TH1F("hJet1q_axis2","",80,0.,0.16);
	hJet1q_axis2->GetXaxis()->SetTitle("#sigma_{2} 1^{st} q-jet");
	TH1F *hJet1q_mult= new TH1F("hJet1q_mult","",30,0,30.);
	hJet1q_mult->GetXaxis()->SetTitle("N 1^{st} q-jet");
	TH1F *hJet1q_leadTrackPt= new TH1F("hJet1q_leadTrackPt","",20,0,100.);
	hJet1q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 1^{st} q-jet");


	TH1F *hJet2q_pt = new TH1F("hJet2q_pt","",17,30,200);
	hJet2q_pt->GetXaxis()->SetTitle("p_{T} 2^{nd} q-jet");
	TH1F *hJet2q_eta = new TH1F("hJet2q_eta","",20,-5,5);
	hJet2q_eta->GetXaxis()->SetTitle("#eta 2^{nd} q-jet");
	TH1F *hJet2q_ptd = new TH1F("hJet2q_ptd","",100,0,1);
	hJet2q_ptd->GetXaxis()->SetTitle("ptd 2^{nd} q-jet");
	TH1F *hJet2q_axis2= new TH1F("hJet2q_axis2","",80,0.,0.16);
	hJet2q_axis2->GetXaxis()->SetTitle("#sigma_{2} 2^{nd} q-jet");

	TH1F *hist_bins = new TH1F("bins","",80,0.,0.16);
	


	TH1F *hJet2q_mult= new TH1F("hJet2q_mult","",30,0,30.);
	hJet2q_mult->GetXaxis()->SetTitle("N 2^{nd} q-jet");
	TH1F *hJet2q_leadTrackPt= new TH1F("hJet2q_leadTrackPt","",20,0,100.);
	hJet2q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 2^{nd} q-jet");

	TH1F *hblike1 = new TH1F("hblike1","",100,-1,1);
	hblike1->GetXaxis()->SetTitle("b-likelihood of 1^{st} b-jet");
	TH1F *hblike2 = new TH1F("hblike2","",100,-1,1);
	hblike2->GetXaxis()->SetTitle("b-likelihood of 2^{nd} b-jet");

	TH1F *hmet = new TH1F("hmet","",100,0.,1000.);
	hmet->GetXaxis()->SetTitle("MET p_{T} (GeV)");
	TH1F *hrho = new TH1F("hrho","",60,0.,30.);
	hrho->GetXaxis()->SetTitle("rho");


	TH1F *hselLeptons_tightId = new TH1F("hselLeptons_tightId","",8,-2.5,5.5 );
	hselLeptons_tightId->GetXaxis()->SetTitle("selLeptons_tightId[0]");
	TH1F *hselLeptons_relIso03= new TH1F("hselLeptons_relIso03","",50,-5.5,-0.5 );
	hselLeptons_relIso03->GetXaxis()->SetTitle("selLeptons_relIso03[0]");
	TH1F *hselLeptons_chargedHadRelIso03 = new TH1F("hselLeptons_chargedHadRelIso03","",50,-5.5,-0.5 );
	hselLeptons_chargedHadRelIso03->GetXaxis()->SetTitle("selLeptons_chargedHadRelIso03[0]");
	TH1F *hselLeptons_pfRelIso03= new TH1F("hselLeptons_pfRelIso03","",50,-5.5,-0.5 );
	hselLeptons_pfRelIso03->GetXaxis()->SetTitle("selLeptons_pfRelIso03[0]");

	TH1F *hqgl1_VBF = new TH1F("hqgl1_VBF","",50,0,1);
	hqgl1_VBF->GetXaxis()->SetTitle("VBF QGD 1^{st} q-jet");
	TH1F *hqgl2_VBF = new TH1F("hqgl2_VBF","",50,0,1);
	hqgl2_VBF->GetXaxis()->SetTitle("VBF QGD 2^{nd} q-jet");

	TH1F *hbdt = new TH1F("hbdt","",100,-1,1);
	hbdt->GetXaxis()->SetTitle("BDT output");


	TH1F *hbdt_first = new TH1F("hbdt_first","",1000,-1,1);
	hbdt_first->GetXaxis()->SetTitle("BDT output");

	TH1F *hbdt_second = new TH1F("hbdt_second","",1000,-1,1);
	hbdt_second->GetXaxis()->SetTitle("BDT output");

	
	TH1F *hbtag_log = new TH1F("hbtag_log","",100,0.,10);
	hbtag_log->GetXaxis()->SetTitle("-log(1-CSV1)");
	TH1F *hbtag2_log = new TH1F("hbtag2_log","",100,0.,10);
	hbtag2_log->GetXaxis()->SetTitle("-log(1-CSV2)");

 	TProfile *hprof  = new TProfile("hprof","",1000,-1.,1.,0.,1500.);
	hprof->GetXaxis()->SetTitle("BDT output");
	hprof->GetYaxis()->SetTitle("<m_{bb}> (GeV)");
	
	TH1F *hMbbcat3 = new TH1F("hMbbcat3","",300.,200.,3200.);
	hMbbcat3->GetXaxis()->SetTitle("m_{qq} (GeV)");
	TH1F* hMbbcat2 = new TH1F("hMbbcat2","",300,200.,3200.);
	hMbbcat2->GetXaxis()->SetTitle("m_{bb} (GeV)");
	TH1F* hMbbcat1 = new TH1F("hMbbcat1","",300,200.,3200.);
	hMbbcat1->GetXaxis()->SetTitle("m_{bb} (GeV)");



	TH1F *hvLeptons_dz = new TH1F("hvLeptons_dz","",100,0.,0.2);
	hvLeptons_dz->GetXaxis()->SetTitle("vLeptons_dz");
	TH1F *hvLeptons_edz = new TH1F("hvLeptons_edz","",100,0.,0.2);
	hvLeptons_edz->GetXaxis()->SetTitle("vLeptons_edz");
	TH1F *hvLeptons_sig = new TH1F("hvLeptons_sig","",200,0.,10.);
	hvLeptons_sig->GetXaxis()->SetTitle("vLeptons_significance");




		
    TH1F *hJet1b_pt = new TH1F("hJet1b_pt", "", 30, 0., 600.);
    hJet1b_pt->GetXaxis()->SetTitle("1^{st}b-jet p_{T} (GeV)");
    TH1F *hJet2b_pt = new TH1F("hJet2b_pt", "", 30, 0., 600.);
    hJet2b_pt->GetXaxis()->SetTitle("2^{nd} b-jet p_{T} (GeV)");
    TH1F *hJetbb_pt = new TH1F("hJetbb_pt", "", 30, 0., 600.);
    hJetbb_pt->GetXaxis()->SetTitle("bb p_{T} (GeV)");
    TH1F *hJetqqbb_pz = new TH1F("hJetqqbb_pz", "", 200, -4000, 4000.);
    hJetqqbb_pz->GetXaxis()->SetTitle("qqbb p_{z} (GeV)");
    
    
	TH1F *hJet1b_eta = new TH1F("hJet1b_eta", "", 80, -5., 5.);
    hJet1b_eta->GetXaxis()->SetTitle("1^{st} b-jet #eta");
    TH1F *hJet2b_eta = new TH1F("hJet2b_eta", "", 80, -5., 5.);
    hJet2b_eta->GetXaxis()->SetTitle("2^{nd} b-jet #eta");
	TH1F *hEtaQBplus = new TH1F("hEtaQBplus","",160.,-8.,8.);
	hEtaQBplus->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}+#Delta#eta_{qb}^{backward}");
	TH1F *hEtaQBminus = new TH1F("hEtaQBminus","",160.,-8.,8.);
	hEtaQBminus->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}-#Delta#eta_{qb}^{backward}");

 	TProfile *hprof_etot  = new TProfile("hprof_etot","E_{tot} (GeV)",150,0.,6000.,0.,1500.);
 	TProfile *hprof_mqq  = new TProfile("hprof_mqq","M_{qq} (GeV)",150,0.,3000.,0.,1500.);
 	TProfile *hprof_etaqq  = new TProfile("hprof_etaqq","#Delta#eta_{qq}",90,0.,9.,0.,1500.);
 	TProfile *hprof_qqbb_pz  = new TProfile("hprof_qqbb_pz","p_{z} of qqbb system",200,-4000.,4000.,0.,1500.);
 	TProfile *hprof_qqbb_pt  = new TProfile("hprof_qqbb_pt","p_{T} of qqbb system",30,0.,150.,0.,1500.);
 	TProfile *hprof_x1  = new TProfile("hprof_x1","x_{1}",100,0,1.,0.,1500.);
 	TProfile *hprof_x2  = new TProfile("hprof_x2","x_{2}",100,0,1.,0,1500.);
 	TProfile *hprof_btag_log1  = new TProfile("hprof_btag_log1","-log(1-CSV_{1})",100,0,10.,0,1500.);
 	TProfile *hprof_btag_log2  = new TProfile("hprof_btag_log2","-log(1-CSV_{2})",100,0,10.,0,1500.);
 	TProfile *hprof_etaqb  = new TProfile("hprof_etaqb","#Delta#eta_{qb}^{forward}+#Delta#eta_{qb}^{backward}",160,-8,8.,0,1500.);
 	TProfile *hprof_softn2  = new TProfile("hprof_softn2","N soft jets, p_{T} > 2 GeV",25,0,25,0,1500.);
 	TProfile *hprof_phiqq  = new TProfile("hprof_phiqq","#Delta#phi_{qq}",32,0,3.2,0,1500.);
 	TProfile *hprof_axis2_1  = new TProfile("hprof_axis2_1","#sigma_{2}, 1^{st} q-jet",80,0,0.16,0,1500.);
 	TProfile *hprof_axis2_2  = new TProfile("hprof_axis2_2","#sigma_{2}, 2^{nd} q-jet",80,0,0.16,0,1500.);
 	TProfile *hprof_mult_1  = new TProfile("hprof_mult_1","N, 1^{st} q-jet",80,0,30,0,1500.);
 	TProfile *hprof_mult_2  = new TProfile("hprof_mult_2","N, 2^{nd} q-jet",80,0,30,0,1500.);

			const int numProfiles = 17;
			TProfile *histProfile[numProfiles] = {hprof, hprof_etot, hprof_mqq, hprof_etaqq, hprof_qqbb_pz, hprof_qqbb_pt, hprof_x1, hprof_x2, hprof_btag_log1, hprof_btag_log2, hprof_etaqb, hprof_softn2, hprof_phiqq, hprof_axis2_1, hprof_axis2_2, hprof_mult_1, hprof_mult_2};

   		const int numArray= 98; 
   		TH1F* histArray[numArray] = {hJet1_pt,hJet2_pt,hJet3_pt,hJet4_pt,  hJet1_eta,hJet2_eta,hJet3_eta,hJet4_eta,  hJet1_phi,hJet2_phi,hJet3_phi,hJet4_phi, hMqq, hEtaQQ, hPhiBB, hEtaSoftJets, hPtSoftJets,hMassSoftJets,hHTsoft,hSoft_n2,hSoft_n5,hSoft_n10,hMbb,hqgl,hbtag,hqgl2,hbtag2,hPtSoftJets2,hPtSoftJets3,hcosOqqbb,hEtaQB1, hEtaQB2, hPhiQB1, hPhiQB2,hx1,hx2,hVB1_mass, hVB2_mass, hEtot, hPxtot, hPytot, hPztot, hJet5_pt,hPtqqbb, hPhiqqbb, hEtaqqbb, hJet1_pt_bin,hJet2_pt_bin,hJet3_pt_bin,hJet4_pt_bin, hMqq_bin, hnPVs, hMbb_fsr, hMbb_met_fsr,hMbb_met_fsr_long, hMbb_met_fsr_dphibb, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult,hblike1,hblike2, hmet,  hselLeptons_tightId , hselLeptons_relIso03 , hselLeptons_chargedHadRelIso03, hselLeptons_pfRelIso03, hqgl1_VBF,hqgl2_VBF, hPhiQQ,hbdt, hbtag_log, hbtag2_log,   hvLeptons_dz, hvLeptons_edz, hvLeptons_sig, hrho, hJet1q_leadTrackPt, hJet2q_leadTrackPt, hJet1b_pt, hJet2b_pt, hJet1b_eta, hJet2b_eta, hJetbb_pt, hJetqqbb_pz, hEtaQBplus, hEtaQBminus, hbdt_first, hbdt_second , hMqq_pt, hMbb_met_fsr_bg, hMbb_met_fsr_sg};
			for (int i=0;i<numArray;i++){
				histArray[i]->Sumw2();
			}
	
	float qq_matching = 0;
	float qq_matching_all = 0;
	

	int nentries = tree_initial->GetEntries() ;
	
//int scale_counter=5;

	TF1 *func_r = new TF1("erffunc",erf,0.,1000.,3);

	for (int entry=0; entry<nentries;++entry){
//	for (int entry=0; entry<1000000;++entry){
        tree_initial->GetEntry(entry);

		if (JSON!=1) {
			continue;
		}


		if (region.CompareTo("analysis")==0) {
			if (!((v_type==-1)||(v_type>3))) continue;
	
		} 
		else if (region.CompareTo("controlTop")==0) {
			if (!((v_type==2)||(v_type==3))) continue;
			if (met_pt<50) continue;
			if (!(selLeptons_relIso03[0]<0.08)) continue;
	//		if (!(selLeptons_chargedHadRelIso03[0]<0.05)) continue;
		}
		else if (region.CompareTo("controlDY")==0) {
			if (!((v_type==0)||(v_type==1))) continue;
			if (!((selLeptons_tightId[0]==1)||(selLeptons_tightId[0]==3))) continue;
			if (!((selLeptons_tightId[1]==1)||(selLeptons_tightId[1]==3))) continue;
		}
		

		if (data==1) PU=1.;
		else PU=puweight;
		genweight0 = genweight/TMath::Abs(genweight);
		genweight=genweight/TMath::Abs(genweight)*PU*bTagWeight;
		genweight/=events_generated/xsec[file_tag]; 
	

		int btag_max1_number = -1;
		int btag_max2_number = -1;
		int eta_max1_number = -1;
		int eta_max2_number = -1;
		int pt_num1 = -1;
		int pt_num2 = -1;
		TLorentzVector Bjet1;
		TLorentzVector Bjet2;
		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector Qjet1_pt;
		TLorentzVector Qjet2_pt;
		TLorentzVector qq;
		int good_jets = 0;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		vector<float> jets_btag;
		vector<float> jets_pt;
		TLorentzVector met;
		vector<TLorentzVector> FSRjet;
		///////////////////////
		//preselection/////
		//////////////////////

		if  (preselection(nJets, Jet.pt,Jet.eta, Jet.phi, Jet.mass, Jet.btag, Jet.id,Jet.puId, btag_max1_number, btag_max2_number, eta_max1_number, eta_max2_number,HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200, Bjet1, Bjet2, Qjet1, Qjet2, qq,jets_pv, jets_btag, jets_pt, good_jets, jets_indices, met_pt, met_phi, met, FSRjet) != 0)  continue; 


		Float_t Mqq = qq.M();
		Float_t bbDeltaPhi = TMath::Abs(Bjet1.DeltaPhi(Bjet2));
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));

///////////////////////////////////////////////////////////////////////////////
//////////////////////////Q jets chosen by pT/////////////////////////////////
    	jets_pt[btag_max1_number] = 0; 
    	jets_pt[btag_max2_number] = 0; 
      pt_num1 =  max_element(jets_pt.begin(), jets_pt.end()) - jets_pt.begin() ;
    	jets_pt[pt_num1] = 0; 
      pt_num2 =  max_element(jets_pt.begin(), jets_pt.end()) - jets_pt.begin();
	
		Qjet1_pt = jets_pv[pt_num1];
		Qjet2_pt = jets_pv[pt_num2];
		TLorentzVector qq_pt = Qjet1_pt + Qjet2_pt;	
		Float_t Mqq_pt = qq_pt.M();
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



		presel+=genweight;

		presel_vtype[(int)(v_type+1)]+=genweight;



		///////TRIGGER WEIGHTS////////////////////////// 
		float trigCor0_nom_double[10] = {1.,1.,1., 1. , 1., 1.31914e+00  };
		float trigCor1_nom_double[10] = {7.60461e+01,  6.92756e+01 , 4.90344e+01,-9.29334e+00,-4.39848e+01, -6.02564e+01  };
		float trigCor2_nom_double[10] = { 1.00000e+02,  2.05669e+01 ,2.29954e+01,3.86211e+01,1.01285e+01,1.41746e+01  };
		float trigCor0_up_double[10] = {1.,1.,1.,1.,1. ,   1.45106e+00 };
		float trigCor1_up_double[10] = { 5.70346e+01,6.23481e+01 , 4.32436e+01,-7.43467e+00,  -1.33079e+01,-1.93149e+01   };
		float trigCor2_up_double[10] = { 1.05000e+02,2.15952e+01, 2.41452e+01 ,  4.05522e+01,     1.21542e+01,  1.70095e+01  };
		float trigCor0_down_double[10] = {1.,1.,1., 1.,1. ,  1.18723e+00};
		float trigCor1_down_double[10] = {9.50576e+01, 7.62032e+01 ,  5.39379e+01,-1.11520e+01, -6.31126e+01,  -7.37354e+01 };
		float trigCor2_down_double[10] = {9.50000e+01,1.64535e+01  ,2.18456e+01 ,  3.66901e+01, 8.41821e+00, 1.27518e+01  };

		float trigWeight=1;
		func_r->FixParameter(0,1.);
		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[0] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[0] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[0] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[0] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[0] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[0] );
			trigWeight*=func_r->Eval(jets_pv[0].Pt());
		}
		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[1] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[1] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[1] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[1] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[1] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[1] );
			trigWeight*=func_r->Eval(jets_pv[1].Pt());
		}
		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[2] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[2] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[2] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[2] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[2] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[2] );
			trigWeight*=func_r->Eval(jets_pv[2].Pt());
		}

		//////////////////////forth jet, not used////////////////////
/*		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[3] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[3] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[3] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[3] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[3] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[3] );
			trigWeight*=func_r->Eval(jets_pv[3].Pt());
		}
*/

	if (data!=1){
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[4] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[4] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[4] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[4] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[4] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[4] );
			trigWeight*=func_r->Eval(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max1_number]]));
		}
		if (data!=1){
			if (whichTrigWeight==0) func_r->FixParameter(0,trigCor0_nom_double[5] );
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[5] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[5] );
			if (whichTrigWeight==1) func_r->FixParameter(0,trigCor0_up_double[5] );
			if (whichTrigWeight==1) func_r->FixParameter(1,trigCor1_up_double[5] );
			if (whichTrigWeight==1) func_r->FixParameter(2,trigCor2_up_double[5] );
			if (whichTrigWeight==2) func_r->FixParameter(0,trigCor0_down_double[5] );
			if (whichTrigWeight==2) func_r->FixParameter(1,trigCor1_down_double[5] );
			if (whichTrigWeight==2) func_r->FixParameter(2,trigCor2_down_double[5] );
			trigWeight*=func_r->Eval(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max2_number]]));
		}


		if (applyTrigWeight==1) genweight*=trigWeight;
		
/////////////////////////////////////////////////////////


		TLorentzVector bb;
		bb = Bjet1+Bjet2;
		Float_t Mbb = bb.M();


		TLorentzVector bb_fsr = Bjet1+Bjet2;
		for (int i=0;i<FSRjet.size();i++)
			bb_fsr+= FSRjet[i];
		
		float Mbb_fsr=bb_fsr.M();
		Float_t Mbb_met_fsr = (bb_fsr+met).M(); 

		TLorentzVector bbqq;
		bbqq = bb_fsr + qq;

		Float_t cosObbqq =TMath::Cos( ( ( Bjet1.Vect() ).Cross(Bjet2.Vect()) ).Angle( ( Qjet1.Vect() ).Cross(Qjet2.Vect()) ) );	

		Float_t EtaBQ1;
	 	Float_t EtaBQ2;
		Float_t PhiBQ1; 	
		Float_t PhiBQ2;
		 if (Qjet1.Eta() >= Qjet2.Eta()) {
			if (Bjet1.Eta() >= Bjet2.Eta())  {
				EtaBQ1 = Qjet1.Eta()-Bjet1.Eta();
				PhiBQ1 = TMath::Abs(Bjet1.DeltaPhi(Qjet1));		
			}
			else {
				EtaBQ1 = Qjet1.Eta()-Bjet2.Eta();
				PhiBQ1 = TMath::Abs(Bjet2.DeltaPhi(Qjet1));	
			}	
		} else if (Bjet1.Eta() >= Bjet2.Eta()) {
				EtaBQ1 = Qjet2.Eta()-Bjet1.Eta();
				PhiBQ1 = TMath::Abs(Bjet1.DeltaPhi(Qjet2));	
				
				}
			else {
				EtaBQ1 = Qjet2.Eta()-Bjet2.Eta();
				PhiBQ1 = TMath::Abs(Bjet2.DeltaPhi(Qjet2));	
			}


		 if (Qjet1.Eta() <= Qjet2.Eta()) {
			if (Bjet1.Eta() <= Bjet2.Eta())  {
				EtaBQ2 = Qjet1.Eta()-Bjet1.Eta();
				PhiBQ2 = TMath::Abs(Bjet1.DeltaPhi(Qjet1));		
			}
			else {
				EtaBQ2 = Qjet1.Eta()-Bjet2.Eta();
				PhiBQ2 = TMath::Abs(Bjet2.DeltaPhi(Qjet1));	
			}	
		} else if (Bjet1.Eta() <= Bjet2.Eta()) {
				EtaBQ2 = Qjet2.Eta()-Bjet1.Eta();
				PhiBQ2 = TMath::Abs(Bjet1.DeltaPhi(Qjet2));	
				
				}
			else {
				EtaBQ2 = Qjet2.Eta()-Bjet2.Eta();
				PhiBQ2 = TMath::Abs(Bjet2.DeltaPhi(Qjet2));	
			}
		


		Float_t Etot = Bjet1.E()+Bjet2.E()+Qjet1.E()+Qjet2.E();
		Float_t PzTot = Bjet1.Pz()+Bjet2.Pz()+Qjet1.Pz()+Qjet2.Pz();
		Float_t PxTot = Bjet1.Px()+Bjet2.Px()+Qjet1.Px()+Qjet2.Px();
		Float_t PyTot = Bjet1.Py()+Bjet2.Py()+Qjet1.Py()+Qjet2.Py();
	
		Float_t x1 = 0.;
		Float_t x2 = 0.;
		x1 = (Etot + PzTot)/2./13000.;
		x2 = (Etot - PzTot)/2./13000.;

		TLorentzVector q1,q2,q1_after,q2_after, VB1, VB2;
		Float_t VB1_mass = -1;
		Float_t VB2_mass = -1.;
		
		q1.SetPxPyPzE(0.,0.,13000./2.*x1,13000./2.*x1);
		q2.SetPxPyPzE(0.,0.,-13000./2.*x2,13000./2.*x2);
		q1_after.SetPxPyPzE(Qjet1.Px(),Qjet1.Py(),Qjet1.Pz(),Qjet1.E());
		q2_after.SetPxPyPzE(Qjet2.Px(),Qjet2.Py(),Qjet2.Pz(),Qjet2.E());
		if (q1_after.Eta()>=0.) {
			VB1 = -q1_after+q1;
			VB2 = -q2_after+q2;
		} else {
			VB1 = -q2_after+q1;
			VB2 = -q1_after+q2;
		} 
		VB1_mass = TMath::Abs(VB1.M());
		VB2_mass = TMath::Abs(VB2.M()); 

		Float_t Jet5pt = -1;
		for (int i=0;i<good_jets;i++){
			if ((i==btag_max1_number)||(i==btag_max2_number)||(i==eta_max1_number)||(i==eta_max2_number)||(Jet.id[i]<=2)||(Jet.puId[i]==0)) continue;
			Jet5pt=jets_pv[i].Pt();
			break;
		}



			counter++;

            hJet1_pt->Fill(jets_pv[0].Pt(),genweight);
            hJet1_eta->Fill(jets_pv[0].Eta(),genweight);
            hJet1_phi->Fill(jets_pv[0].Phi(),genweight);
            hJet2_pt->Fill(jets_pv[1].Pt(),genweight);
            hJet2_eta->Fill(jets_pv[1].Eta(),genweight);
            hJet2_phi->Fill(jets_pv[1].Phi(),genweight);
            hJet3_pt->Fill(jets_pv[2].Pt(),genweight);
            hJet3_eta->Fill(jets_pv[2].Eta(),genweight);
            hJet3_phi->Fill(jets_pv[2].Phi(),genweight);
            hJet4_pt->Fill(jets_pv[3].Pt(),genweight);
            hJet4_eta->Fill(jets_pv[3].Eta(),genweight);
            hJet4_phi->Fill(jets_pv[3].Phi(),genweight);
            hVtype->Fill(v_type,genweight);
            hJet5_pt->Fill(Jet5pt,genweight);
			hMqq->Fill(Mqq,genweight);
			hMqq_pt->Fill(Mqq_pt,genweight);
			hEtaQQ->Fill(qqDeltaEta,genweight);
			hPhiQQ->Fill(qqDeltaPhi,genweight);
			hPhiBB->Fill(bbDeltaPhi,genweight);
			hEtaSoftJets->Fill(Jet.soft_eta[0],genweight);
			hPtSoftJets->Fill(Jet.soft_pt[0],genweight);
			hPtSoftJets2->Fill(Jet.soft_pt[1],genweight);
			hPtSoftJets3->Fill(Jet.soft_pt[2],genweight);
			hMassSoftJets->Fill(Jet.soft_mass[0],genweight);
			hHTsoft->Fill(Jet.HTsoft,genweight);
			hSoft_n2->Fill(Jet.nsoft2, genweight);
			hSoft_n5->Fill(Jet.nsoft5, genweight);
			hSoft_n10->Fill(Jet.nsoft10, genweight);
			hMbb->Fill(Mbb,genweight);
			hqgl->Fill(Jet.qgl[jets_indices[eta_max1_number]],genweight);
			hbtag->Fill(Jet.btag[jets_indices[btag_max1_number]],genweight);
			hqgl2->Fill(Jet.qgl[jets_indices[eta_max2_number]],genweight);
			hbtag2->Fill(Jet.btag[jets_indices[btag_max2_number]],genweight);
			hcosOqqbb->Fill(cosObbqq,genweight);
			hEtaQB1->Fill(EtaBQ1,genweight);
			hEtaQB2->Fill(EtaBQ2,genweight);
			hPhiQB1->Fill(PhiBQ1,genweight);
			hPhiQB2->Fill(PhiBQ2,genweight);
			hx1->Fill(x1,genweight);
			hx2->Fill(x2,genweight);
			hVB1_mass->Fill(VB1_mass,genweight);
			hVB2_mass->Fill(VB2_mass,genweight);
			hEtot->Fill(Etot,genweight);
			hPxtot->Fill(PxTot,genweight);
			hPytot->Fill(PyTot,genweight);
			hPztot->Fill(PzTot,genweight);
			hPtqqbb->Fill(bbqq.Pt(),genweight);
			hPhiqqbb->Fill((-1)*bbqq.Phi(),genweight);
			hEtaqqbb->Fill(TMath::Abs(bbqq.Eta()), genweight);
            hJet1_pt_bin->Fill(jets_pv[0].Pt(),genweight);
            hJet2_pt_bin->Fill(jets_pv[1].Pt(),genweight);
            hJet3_pt_bin->Fill(jets_pv[2].Pt(),genweight);
            hJet4_pt_bin->Fill(jets_pv[3].Pt(),genweight);
				hMqq_bin->Fill(Mqq,genweight);
				hnPVs->Fill(nPVs,genweight);
				hMbb_fsr->Fill(Mbb_fsr,genweight);
				hMbb_met_fsr->Fill(Mbb_met_fsr,genweight);
				if ((bdt>-0.6)&&(bdt<0.6)) hMbb_met_fsr_bg->Fill(Mbb_met_fsr,genweight);
		//		if (data!=1) if ((bdt>0.8)) hMbb_met_fsr_sg->Fill(Mbb_met_fsr,genweight);
			 if ((bdt>0.8)) hMbb_met_fsr_sg->Fill(Mbb_met_fsr,genweight);
				hMbb_met_fsr_long->Fill(Mbb_met_fsr,genweight);
				if (bbDeltaPhi < 2.4) hMbb_met_fsr_dphibb->Fill(Mbb_met_fsr,genweight);


				hJet1q_pt->Fill(jets_pv[eta_max1_number].Pt(),genweight);
				hJet1q_eta->Fill(jets_pv[eta_max1_number].Eta(),genweight);
				hJet1q_ptd->Fill(Jet.ptd[jets_indices[eta_max1_number]],genweight);
				hJet1q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[eta_max1_number]]),genweight);
				hJet1q_mult->Fill(Jet.mult[jets_indices[eta_max1_number]],genweight);
				hJet1q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[eta_max1_number]],genweight);
				hJet2q_pt->Fill(jets_pv[eta_max2_number].Pt(),genweight);
				hJet2q_eta->Fill(jets_pv[eta_max2_number].Eta(),genweight);
				hJet2q_ptd->Fill(Jet.ptd[jets_indices[eta_max2_number]],genweight);
				hJet2q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[eta_max2_number]]),genweight);
				hJet2q_mult->Fill(Jet.mult[jets_indices[eta_max2_number]],genweight);
				hJet2q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[eta_max2_number]],genweight);
				hblike1->Fill(Jet.blike_VBF[jets_indices[btag_max1_number]],genweight);
				hblike2->Fill(Jet.blike_VBF[jets_indices[btag_max2_number]],genweight);
				hmet->Fill(met_pt,genweight);
				hselLeptons_tightId->Fill(selLeptons_tightId[0],genweight);

				hJet1b_pt->Fill(jets_pv[btag_max1_number].Pt(),genweight);
				hJet1b_eta->Fill(jets_pv[btag_max1_number].Eta(),genweight);
				hJet2b_pt->Fill(jets_pv[btag_max2_number].Pt(),genweight);
				hJet2b_eta->Fill(jets_pv[btag_max2_number].Eta(),genweight);
				hJetqqbb_pz->Fill(bbqq.Pz(),genweight);
				hEtaQBplus->Fill(EtaBQ1+EtaBQ2 , genweight);
				hEtaQBminus->Fill(EtaBQ1-EtaBQ2, genweight);
				hJetbb_pt->Fill(bb_fsr.Pt() , genweight);
			
				Float_t logIso  = -5;
				if (selLeptons_relIso03[0] > TMath::Exp(-5)) logIso = TMath::Log(selLeptons_relIso03[0]);
				hselLeptons_relIso03->Fill(logIso,genweight);
				logIso  = -5;
			//	if (selLeptons_chargedHadRelIso03[0] > TMath::Exp(-5)) logIso = TMath::Log(selLeptons_chargedHadRelIso03[0]);
			//	hselLeptons_chargedHadRelIso03->Fill(logIso,genweight);
				logIso  = -5;
				if (selLeptons_pfRelIso03[0] > TMath::Exp(-5)) logIso = TMath::Log(selLeptons_pfRelIso03[0]);
				hselLeptons_pfRelIso03->Fill(logIso,genweight);
 
			hbdt->Fill(bdt,genweight);



			hbtag_log->Fill(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max1_number]]),genweight);
			hbtag2_log->Fill(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max2_number]]),genweight);
			hrho->Fill(rho,genweight);
			hprof->Fill(bdt,Mbb_met_fsr,genweight);
			hprof_etot->Fill(Etot, Mbb_met_fsr, genweight);
			hprof_mqq->Fill(Mqq, Mbb_met_fsr, genweight);
			hprof_etaqq->Fill(qqDeltaEta, Mbb_met_fsr, genweight);
			hprof_qqbb_pz->Fill(bbqq.Pz(), Mbb_met_fsr, genweight);
			hprof_qqbb_pt->Fill(bbqq.Pt(), Mbb_met_fsr, genweight);
			hprof_x1->Fill(x1, Mbb_met_fsr, genweight);
			hprof_x2->Fill(x2, Mbb_met_fsr, genweight);
			hprof_btag_log1->Fill(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max1_number]]), Mbb_met_fsr, genweight);
			hprof_btag_log2->Fill(-1.*TMath::Log(1.-Jet.btag[jets_indices[btag_max2_number]]), Mbb_met_fsr, genweight);
			hprof_etaqb->Fill(EtaBQ1+EtaBQ2, Mbb_met_fsr, genweight);
			hprof_softn2->Fill(Jet.nsoft2, Mbb_met_fsr, genweight);
			hprof_phiqq->Fill(qqDeltaPhi, Mbb_met_fsr, genweight);
			hprof_axis2_1->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[eta_max1_number]]), Mbb_met_fsr, genweight);
			hprof_axis2_2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[eta_max2_number]]), Mbb_met_fsr, genweight);
			hprof_mult_1->Fill(Jet.mult[jets_indices[eta_max1_number]], Mbb_met_fsr, genweight);
			hprof_mult_2->Fill(Jet.mult[jets_indices[eta_max2_number]], Mbb_met_fsr, genweight);



		
			if (genweight>0) pos_weight_presel++;



		float mcweight=genweight*events_generated/xsec[file_tag]; 
		
		if (genweight>0) gen_pos_weight+=mcweight;
		if (genweight<0) gen_neg_weight+=mcweight;
		if (genweight>0) gen_pos+=genweight0;
		if (genweight<0) gen_neg+=genweight0;

				
			global_counter++;
        }

		cout<<counter<<endl;
		TFile file(output+"/"+file_tag+"_"+region+"_"+heppyVersion+"_"+postfix+".root","recreate");
    
		for (int i=0;i<numArray;++i){
    	    	histArray[i]->SetLineWidth(2);
    	   	histArray[i]->GetYaxis()->SetTitle("N_{events}");
       		histArray[i]->GetYaxis()->SetTitleFont(42);
       		histArray[i]->GetYaxis()->SetTitleSize(0.060);
        		histArray[i]->GetYaxis()->SetTitleOffset(0.8);
        		histArray[i]->SetLineColor(kBlue);
        		histArray[i]->Draw();
        		histArray[i]->Write();
   		} 
			for (int i=0;i<numProfiles;i++){
				histProfile[i]->SetLineWidth(2);
      	  	histProfile[i]->SetLineColor(kBlue);
       	 	histProfile[i]->Draw();
       	 	histProfile[i]->Write();
			}
    		file.Write();
    		file.Close();
	 ofstream out(output+"/"+file_tag+"_"+region+"_"+heppyVersion+"_"+postfix+".txt");
	out<< "positive pure selected = "<<gen_pos<<"  , positive weighted selected =  "<<gen_pos_weight<<" , negative pure selected = "<<gen_neg<< ", negative weighted selected = "<<gen_neg_weight<< ", all evetns in the begining = "<<events_generated<<" , xsec = "<<xsec[file_tag]<<endl;
	out<<"positive weight in so many events : "<<  pos_weight_presel<<endl;

return 0;
    
}

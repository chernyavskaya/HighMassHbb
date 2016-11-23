#include <stdio.h>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>
#include "TROOT.h"
#include "TPad.h"
#include "TGaxis.h"
#include <boost/format.hpp>
#include <array>
#include <algorithm>
#include <iterator>

#define SWAP(A, B) { Float_t t = A; A = B; B = t; }
#define SWAP2(A, B) { B.swap(A); }

void bubblesort(Float_t *a, std::string *str, int n)
{
  int i, j;
  for (i = n - 1; i >= 0; i--)
  {
    for (j = 0; j < i; j++)
    {
      if (a[j] < a[j + 1]){
        SWAP( a[j], a[j + 1] );
	//	cout<<str[j]<<"   "<<str[j+1]<<endl;
		SWAP2( str[j], str[j+1]);
	//	cout<<str[j]<<"   "<<str[j+1]<<endl;
		}
    }
  }
}


Float_t Discr(TH1F *h1, TH1F *h2){  
  Int_t h1s = h1->GetNbinsX();
  Int_t h2s = h2->GetNbinsX();
  if (h1s != h2s) {
  	printf("h1 size %d   h2 size %d \n",h1s,h2s);
  	return -1;
  }
  if (h1->Integral()!=0 ) h1->Scale(1./h1->Integral());
  if (h2->Integral()!=0 ) h2->Scale(1./h2->Integral());
  Float_t adiff = 0;
  for (Int_t i=0;i<h1s;i++){
	adiff +=  TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i));
	//cout << "h1 bin " << i << " " << h1->GetBinContent(i) << endl;
	}
  return adiff/2;
}

using namespace std;

int main(int argc, char* argv[]){
//int analyzer_stackRatio(){
gROOT->ProcessLine(".x /afs/cern.ch/work/n/nchernya/setTDRStyle.C");
int region_type=atoi(argv[1]); // 0 - analysis, 1 - control region , top
 
const int nfiles  = 22;
TString leg_names[nfiles] = {"Data","VBF, m(H) = 750 GeV","","","","Z(ll) + jets","W(l#nu) + jets","Single top","Single t,t-top","Single t,t-antitop","Single t","Z + jets","W + jets","t#bar{t}","QCD","QCD, H_{T}>200 GeV","QCD, H_{T}>300 GeV","QCD, H_{T}>500 GeV","QCD, H_{T}>700 GeV","QCD, H_{T}>1000 GeV","QCD, H_{T}>1500 GeV","QCD, H_{T}>2000 GeV"};
TString set_names[2] = {"DoubleB","SingleB"}; 


TString file_names[nfiles] = { "BTagCSV","VBFSpin0ToBBbar_W_1p0_M_750","VBFSpin0ToBBbar_W_1p0_M_750","VBFSpin0ToBBbar_W_1p0_M_750","VBFSpin0ToBBbar_W_1p0_M_750","DYJetsToLL","WJetsToLNu","ST_s-channel","ST_t-channel_top_4f_inclusiveDecays","ST_t-channel_antitop_4f_inclusiveDecays","ST_tW","DYJetsToQQ","WJetsToQQ",/*"TT","TT_madgraph"*/"TT_powheg","QCD_HT2000toInf","QCD_HT1500to2000","QCD_HT1000to1500","QCD_HT100to200","QCD_HT700to1000","QCD_HT200to300","QCD_HT500to700","QCD_HT300to500"};
int bg_begin;
int qcd_begin=14;
if (region_type==0) bg_begin=8;
else bg_begin=5;

int FILLCOLOR[nfiles] = {1,kRed,kBlue+0,kViolet+2,kGreen+2,kYellow-4,kYellow-8,kYellow-9,kYellow-9,kYellow-9,kYellow-9,kOrange-9,kYellow-10,kOrange-2,kRed-9,kRed-9,kRed-9, kRed-9,kRed-9,kRed-9,kRed-9,kRed-9};
//int LINECOLOR[nfiles] = {1,kRed,kBlue+0,kViolet+2,kGreen+2,1,1,kYellow-9,1,1,1,1,1,1,1,1,1,1,1,1};
int LINECOLOR[nfiles] = {1,kRed,kBlue+0,kViolet+2,kGreen+2,1,1,kYellow-9,kYellow-9,kYellow-9,1,1,1,1,kRed-9,kRed-9,kRed-9, kRed-9,kRed-9,kRed-9,kRed-9,1};
int LINESTYLE[nfiles] = {1,1,2,3,8,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
int LINEWIDTH[nfiles] = {1,3,3,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001};
	
int order[nfiles] = {0,1,2,3,4, 5, 6, 7, 8, 9,10,11,12,13,  21,20,19,14,18,15,17,16};// number in file_names array as they should appear in the legend
//int order[nfiles] = {0,1,2,3,4, 5, 6, 7, 8,9,10,  18,17,16,11,15,12,14,13};// number in file_names array as they should appear in the legend
int order_legend[nfiles]; 
for (int i=0;i<nfiles;i++){
	order_legend[order[i]]=i;
}
	

TString btag[2] = {"v14, DoubleBtag","v14, SingleBtagBlike"}; 
TString region[3]={"_analysis","_controlTop","_controlDY"}; // 0 - analysis, 1 - control region , top; 2 - control region DY

for (int i=0;i<nfiles;i++){
	file_names[i].Prepend("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/HighMassHbb/plotterOutput/v21/");
	file_names[i].Append(region[region_type]);
//	file_names[i].Append("_v21_oldTrigWtMbbcut200");
//	file_names[i].Append("_v21_NoTrigWt_bbMbbcut250");
//	file_names[i].Append("_v21_WoTriggerFinal");
	file_names[i].Append("_v21_TriggerFinal");
//	file_names[i].Append("_v21_TrigWtFinal");
	file_names[i].Append(".root");
}
//TString dir_name= "plots_mbbcut300";
TString dir_name= "plots_triggerFinal";
//TString dir_name= "plots_wo_trigger";
dir_name.Append(region[region_type]);
dir_name.Append("/");
//TString dir_name = "plots_amc/";
Float_t lumi = 2320;

	


TLegend *leg = new TLegend(0.77,0.45,0.92,0.9); //without writing about SF
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetTextSize(0.025);
TLegend *leg2 = new TLegend(0.77,0.35,0.92,0.4); //with  writing about SF
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->SetTextFont(42);
leg2->SetTextSize(0.025);

const int nhistos = 86; //40//52
//TString hist_names[nhistos]={"hJet1_pt","hJet2_pt","hJet3_pt","hJet4_pt","hJet1_eta","hJet2_eta","hJet3_eta","hJet4_eta","hJet1_phi","hJet2_phi","hJet3_phi","hJet4_phi","hMqq", "hEtaQQ", "hPhiBB","hMbb","hbtag","hbtag2","hcosOqqbb","hEtaQB1", "hEtaQB2", "hPhiQB1", "hPhiQB2","hx1","hx2","hVB1_mass","hVB2_mass","hEtot","hPxtot","hPytot","hPztot","hJet5_pt","hPtqqbb","hEtaqqbb","hPhiqqbb","hJet1_pt_bin","hJet2_pt_bin","hJet3_pt_bin","hJet4_pt_bin", "hMqq_bin","hEtaSoftJets", "hPtSoftJets","hMassSoftJets","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hqgl","hqgl2", "hPtSoftJets2","hPtSoftJets3","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult","hMbb_regVBF","hMbb_regVBF_fsr","hblike1","hblike2", "hmet", "hqgl1_VBF", "hqgl2_VBF", "hbdt", "hPhiQQ", "hMbb_regVBF_fsr_high", "hMbb_regVBF_fsr_high_cat","hselLeptons_relIso03"};
TString hist_names[nhistos]={"hJet1_pt","hJet2_pt","hJet3_pt","hJet4_pt","hJet1_eta","hJet2_eta","hJet3_eta","hJet4_eta","hJet1_phi","hJet2_phi","hJet3_phi","hJet4_phi","hMqq", "hEtaQQ", "hPhiBB","hMbb","hbtag","hbtag2","hcosOqqbb","hEtaQB1", "hEtaQB2", "hPhiQB1", "hPhiQB2","hx1","hx2","hVB1_mass","hVB2_mass","hEtot","hPxtot","hPytot","hPztot","hJet5_pt","hPtqqbb","hEtaqqbb","hPhiqqbb","hJet1_pt_bin","hJet2_pt_bin","hJet3_pt_bin","hJet4_pt_bin", "hMqq_bin","hEtaSoftJets", "hPtSoftJets","hMassSoftJets","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hqgl","hqgl2", "hPtSoftJets2","hPtSoftJets3","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult","hMbb_fsr","hMbb_met_fsr","hMbb_met_fsr_long","hblike1","hblike2", "hmet", "hmet", "hmet", "hmet", "hPhiQQ", "hselLeptons_relIso03", "hbtag_log","hbtag2_log","hrho","hJet1q_leadTrackPt","hJet2q_leadTrackPt" ,"hJet1b_pt", "hJet2b_pt", "hJet1b_eta", "hJet2b_eta", "hJetbb_pt", "hJetqqbb_pz", "hEtaQBplus", "hEtaQBminus"};
int UNITS[100]={1,1,1,1 ,0,0,0,0 ,0,0,0,0, 1 ,0,0 ,1, 0,0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1 ,0,0, 1,1,1,1,1, 0, 1,1,1 ,0,0,0,0,0, 1,1, 0, 1,0,0,0,0, 1,0,0,0,0, 1,1 ,0,0, 1, 0,0,0,0 ,1,1, 0,1,1,1,1,1,1,1,1,1,1};
std::array<int,100> LOGY_array = {71,72,51,52,53,54,55,56,57,58,59,60,61};

TString hist_names_sum[nhistos]={};
TString sum_histos_names[nhistos]={};
std::string hist_names_sort[nhistos];
for (int i=0;i<nhistos;i++){
	hist_names_sort[i] = hist_names[i];
	hist_names_sum[i] = hist_names[i];
	hist_names_sum[i].Append("_sum");
	sum_histos_names[i] = hist_names[i];
	sum_histos_names[i].Append("_sum0");
}
Float_t discriminators[nhistos];
TString stacks_names[nhistos];
for (int i=0;i<nhistos;i++){
	stacks_names[i] = hist_names[i];
	stacks_names[i].Prepend("s");
}
TString output_names[nhistos];
TString output_names2[nhistos];
for (int i=0;i<nhistos;i++){
	output_names[i] = hist_names[i];
	output_names[i].Prepend(dir_name);
	output_names[i].Append(region[region_type]);
	output_names2[i]  = output_names[i];
	output_names[i].Append("_highmass_v21.pdf");
	output_names2[i].Append("_highmass_v21.png");
}

TH1F *data_histos[nhistos];
TH1F *data_histos2[nhistos];
TH1F *data_histosTrig[nhistos];
TH1F *signal_histos[nhistos];
TH1F *signal_histos2[nhistos];
TH1F *tthbb_sum_histos[nhistos];
TH1F *tthbb_histos[nhistos];
TH1F *tthnbb_histos[nhistos];
TH1F *gf_histos[nhistos];
TH1F *sum_histos[nhistos];
TH1F *sum_histosUp[nhistos];
TH1F *sum_histosDown[nhistos];
TH1F *histos_forUnc[nhistos];
TH1F *histos_for_legened[nhistos];
TH1F *discr_histos[nhistos];//Mqq,delta eta, delta phi, qgl, btag //12,13,14,21,22
TH1F *hBkgVis[nhistos];
TH1F *hBkgUncUp[nhistos];
TH1F *hBkgUncLo[nhistos];
TH1F *hBkgUncUpTrig[nhistos];
TH1F *hBkgUncLoTrig[nhistos];


int files=0; 
THStack *stacks[nhistos];
for (int i=0;i<nhistos;++i){
	stacks[i] = new THStack(stacks_names[i],"");
}
Double_t totalBG=0.;
Double_t totalQCD=0.;
Double_t totalMC=0.;
Double_t totalData=0.;
Double_t totalDataQCD=0.;
ofstream out_efficiency;
ofstream out_discrimination;
out_efficiency.open(dir_name+"efficiency.txt"); 

//Float_t MC_data = 1/1.17;////////// without dynamical scale factors , cut on Mbb > 250 GeV. 
//Float_t MC_data = 1./1.12;////////////for efficiency for mbb 200 GeV cut 
//Float_t MC_data = 1./1.39;////////////first iteration of dynamical scale factors, cut on Mbb > 250 GeV. TrigWt1 
Float_t MC_data = 1./1.;////////////second iteration of dynamical scale factors, cut on Mbb > 250 GeV. TrigWt1 
//Float_t MC_data = 1./1.03;////////////second iteration of dynamical scale factors, cut on Mbb > 250 GeV. TrigWt1 
Float_t MC_data_top_1l = 1./0.745;
Float_t MC_data_top_2l = 1;
Float_t lumi_qcd=lumi/MC_data;
Float_t lumi_top;
lumi_top=lumi/MC_data_top_1l;
files=0;

do{
	TFile *file_initial;
  	file_initial = TFile::Open(file_names[files]);
	TH1F *histos[100];
	for (int hist=0;hist<nhistos;++hist){
		int divisor=1;
		for (int i=3;i<8;i++){
			if ((((TH1F*)file_initial->Get(hist_names[hist]))->GetNbinsX())%i==0) {
				divisor=i;
				break;
			}
		}
		if ((region_type==1)&&(hist_names[hist].CompareTo("hMbb_met_fsr")==0)) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(4);
//		if (region_type==1) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2);
//		if (region_type==1) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(4);
//		if (region_type==2) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(divisor);
//		if ((region_type!=0)&&(hist==87)) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2);
		if (region_type==2) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2);
//		if ((hist==71)||(hist==72)) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2);
//		if (hist==12) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2.);
		histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("h");
		if (files==0) data_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("data");
		if ((files>0)&&(files<qcd_begin)&&(files!=13))histos[hist]->Scale(lumi);  //top
		if (files>=qcd_begin)histos[hist]->Scale(lumi_qcd);
		if (files==13) histos[hist]->Scale(lumi_top);
		

		if (files==bg_begin) 	sum_histos[hist] = (TH1F*)histos[hist]->Clone(sum_histos_names[hist]);
		if (files>bg_begin)	sum_histos[hist]->Add(histos[hist]); 



		if (files>0) histos[hist]->Sumw2(kFALSE);
		if (hist==4) cout<<files<<"   "<<histos[4]->Integral() <<endl;

		if (files==1) {
			signal_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"newhist");
			signal_histos[hist]->Scale(lumi);
			signal_histos[hist]->Sumw2(kFALSE);
			signal_histos[hist]->SetLineColor(LINECOLOR[files]);
			signal_histos[hist]->SetLineStyle(LINESTYLE[files]);
			signal_histos2[hist]=(TH1F*)signal_histos[hist]->Clone("signalHist2");
		}
		if (files==2) {
			gf_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"tthnb");
			gf_histos[hist]->Scale(lumi);
			gf_histos[hist]->Sumw2(kFALSE);
			gf_histos[hist]->SetLineColor(LINECOLOR[files]);
			gf_histos[hist]->SetLineStyle(LINESTYLE[files]);
		}
		if (files==3) {
			tthbb_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"tthb");
			tthbb_histos[hist]->Scale(lumi);
			tthbb_histos[hist]->Sumw2(kFALSE);
			tthbb_histos[hist]->SetLineColor(LINECOLOR[files]);
			tthbb_histos[hist]->SetLineStyle(LINESTYLE[files]);
		}
		if (files==4) {
			tthnbb_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"tthnb");
			tthnbb_histos[hist]->Scale(lumi);
			tthnbb_histos[hist]->Sumw2(kFALSE);
			tthnbb_histos[hist]->SetLineColor(LINECOLOR[files]);
			tthnbb_histos[hist]->SetLineStyle(LINESTYLE[files]);
			tthbb_sum_histos[hist] = (TH1F*)tthbb_histos[hist]->Clone("sum_tth");
			tthbb_sum_histos[hist]->Add(tthnbb_histos[hist]);
			tthbb_sum_histos[hist]->SetLineColor(LINECOLOR[files]);
			tthbb_sum_histos[hist]->SetLineStyle(LINESTYLE[files-1]);
		}
		histos[hist]->SetLineColor(LINECOLOR[files]);
		histos[hist]->SetLineStyle(LINESTYLE[files]);
		histos[hist]->SetLineWidth(LINEWIDTH[files]);
		histos[hist]->SetFillStyle(FILLSTYLE[files]);
		if ((files!=0)) histos[hist]->SetFillColor(FILLCOLOR[files]);
		if ((region_type!=0)&&(files>=qcd_begin)) histos[hist]->SetFillColor(kRed+2);
		
		if (files==0) {
			histos[hist]->SetMarkerStyle(20);
			data_histos[hist]->SetLineColor(1);
			data_histos[hist]->SetMarkerStyle(20);
			data_histos[hist]->SetMarkerSize(.8);
		}
	 	if (files>=bg_begin) stacks[hist]->Add(histos[hist]);
	//	if (files==bg_begin) 	sum_histos[hist] = (TH1F*)histos[hist]->Clone(sum_histos_names[hist]);
	//	if (files>bg_begin)	sum_histos[hist]->Add(histos[hist]); 
	//
	//	if (files==nfiles-1)	sum_histos[hist] = (TH1F*)stacks[hist]->GetStack()->Last()->Clone("sum_histos"); 
	//
		if (hist==0) histos_for_legened[files] = (TH1F*)histos[0]->Clone("newd");
		if (files==bg_begin)	discr_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("discr");
		if (files>bg_begin)	discr_histos[hist]->Add(histos[hist]); 
   	}
		if (files>=bg_begin) totalBG+=histos[4]->Integral();
		if (files>=qcd_begin) totalQCD+=histos[4]->Integral();
		if (files>0) totalMC+=histos[4]->Integral();
		if (files==0) {totalData+=histos[4]->Integral(); totalDataQCD+=histos[4]->Integral();}
		if ((files>0)&&(files<qcd_begin)) {totalDataQCD-=histos[4]->Integral();}
		if (files==0) out_efficiency<<"Sample  \t\t\t yield(per "<< lumi<<" pb^-1)"<<endl;
		if (files==0) out_efficiency<<leg_names[order[files]]<<"\t \t \t"<< std::setprecision(5)<<histos[4]->Integral() <<endl;
		else out_efficiency<<leg_names[order[files]]<<"\t\t\t  "<<std::setprecision(5)<<histos[4]->Integral()<<endl;
		if (files==nfiles-1) out_efficiency<<"Total BG"<<"\t \t \t  "<<std::setprecision(5)<<totalBG<<endl;
		if (files==nfiles-1) out_efficiency<<"Total MC"<<"\t \t \t  "<<std::setprecision(5)<<totalMC<<endl;
		if (files==nfiles-1) out_efficiency<<"Data/MC"<<"\t \t \t  "<<std::setprecision(3)<<totalData/totalMC<<endl;
		if (files==nfiles-1) out_efficiency<<"DataQCD/QCD"<<"\t \t \t  "<<std::setprecision(3)<<totalDataQCD/totalQCD<<endl;
		if (files==nfiles-1) cout<<"DataQCD/QCD"<<"\t \t \t  "<<std::setprecision(3)<<totalDataQCD/totalQCD<<endl;
files++;
}while (files<nfiles);
out_efficiency.close();

for (int hist=0;hist<nhistos;hist++){
	hBkgUncUp[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncUp");
	hBkgUncLo[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncLo");
  	hBkgVis[hist]   = (TH1F*)sum_histos[hist]->Clone("hbkgVis");
  	for(int i=0;i<hBkgUncUp[hist]->GetNbinsX();i++) {
  		float e = 0.0;
    	if (sum_histos[hist]->GetBinContent(i+1) != 0) {
      	e = sum_histos[hist]->GetBinError(i+1)/sum_histos[hist]->GetBinContent(i+1);
   	}
    	hBkgUncUp[hist]->SetBinContent(i+1,e);
    	hBkgUncLo[hist]->SetBinContent(i+1,-e);
	}
  	hBkgVis[hist]->SetMarkerSize(1.);
 	hBkgVis[hist]->SetFillColor(kBlack);
 	hBkgVis[hist]->SetFillStyle(3004);
 	hBkgUncUp[hist]->SetLineColor(kBlack);
 	hBkgUncUp[hist]->SetLineWidth(1);
 	hBkgUncUp[hist]->SetFillColor(kBlack);
 	hBkgUncUp[hist]->SetFillStyle(3004);
	hBkgUncLo[hist]->SetLineColor(kBlack);
 	hBkgUncLo[hist]->SetLineWidth(1);
  	hBkgUncLo[hist]->SetFillColor(kBlack);
	hBkgUncLo[hist]->SetFillStyle(3004);
}

//Float_t TSF = 1.07;   
Float_t TSF = 1.;   
Float_t kfactors = 0;
kfactors = 1./(MC_data * TSF);
if (region_type==1) kfactors = 1./(MC_data_top_1l * TSF );
if (region_type==2) kfactors = 1./(MC_data_top_2l  * TSF );
cout<<"kfactor = "<<kfactors<<endl;
cout<<"Data/MC = "<< MC_data<<endl;


for (int i=0;i<bg_begin;i++){
	if ((i==7)||(i==8)||(i==9)) continue;
	if (i==0) leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"P");
	if (region_type==0) if (i==1)  leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"L");
}
TString k_factor_str;
TString temp_str;
temp_str.Form("%2.2f",kfactors);
if (region_type!=0) k_factor_str= temp_str.Prepend("t#bar{t} MC #times  ");

if (region_type==0) leg->AddEntry(histos_for_legened[nfiles-1],k_factor_str.Format("QCD (#times %2.2f)",kfactors),"F");
else leg->AddEntry(histos_for_legened[nfiles-1],"QCD","F");
//leg->AddEntry(histos_for_legened[nfiles-1],k_factor_str.Format("QCD",kfactors),"F");

for (int i=nfiles-1;i>=bg_begin;i--){
	if ((i==7)||(i==8)||(i==9)) continue;
	if ((i>=bg_begin)&&(i<qcd_begin)) leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"F");
}
leg->AddEntry(hBkgUncUp[0],"MC stat. unc.","F");
leg2->AddEntry(hBkgUncUp[0],"without dynamical","");
leg2->AddEntry(hBkgUncUp[0],"trigger SF","");



for (int d=0;d<nhistos;d++){
	discriminators[d] = Discr(discr_histos[d],signal_histos2[d]);
}

bubblesort(discriminators, hist_names_sort,nhistos);

//out_discrimination.open("Aftertriggercorr2/"+trigger[set_type]+dir_name+"discrimination.txt");
out_discrimination.open(dir_name+"discrimination.txt");
for (int d=0;d<nhistos;d++){
	if (d==0) out_discrimination<<"Variable &\t d"<<endl;
	out_discrimination<<"$"<<hist_names_sort[d]<<"$"<<" & \t "<< std::setprecision(2)<< discriminators[d]<<endl;
}
out_discrimination.close();

TLatex* tex = new TLatex(0.75,0.95,"2.32 fb^{-1} (13 TeV)");
tex->SetNDC();
tex->SetTextAlign(35);
tex->SetTextFont(42);
tex->SetTextSize(0.035);
tex->SetLineWidth(2);
TLatex *tex1 = new TLatex(0.17,0.95,"CMS");
tex1->SetNDC();
tex1->SetTextAlign(20);
tex1->SetTextFont(61);
tex1->SetTextSize(0.04);
tex1->SetLineWidth(2);
TLatex* tex2 = new TLatex(0.32,0.941,"Preliminary");
tex2->SetNDC();
tex2->SetTextAlign(20);
tex2->SetTextFont(52);
tex2->SetTextSize(0.04);
tex2->SetLineWidth(2);
// CMS info
	float left2 = gStyle->GetPadLeftMargin();
	float right2 = gStyle->GetPadRightMargin();
	float top2 = gStyle->GetPadTopMargin();
	float bottom2 = gStyle->GetPadBottomMargin();
	TPaveText pCMS1(left2,1.-top2,0.4,1.,"NDC");
	pCMS1.SetTextFont(62);
	pCMS1.SetTextSize(top2*0.75);
	pCMS1.SetTextAlign(12);
	pCMS1.SetFillStyle(-1);
	pCMS1.SetBorderSize(0);
	pCMS1.AddText("CMS");
	TPaveText pCMS12(left2+0.08,1.-top2*1.12,0.6,1.,"NDC");
	pCMS12.SetTextFont(52);
	pCMS12.SetTextSize(top2*0.75);
	pCMS12.SetTextAlign(12);
	pCMS12.SetFillStyle(-1);
	pCMS12.SetBorderSize(0);
	pCMS12.AddText("Preliminary");
	TPaveText pCMSset(0.57,1.-top2*2.,0.67,0.92,"NDC");
	pCMSset.SetTextFont(42);
	pCMSset.SetTextSize(top2*0.75);
	pCMSset.SetTextAlign(12);
	pCMSset.SetFillStyle(-1);
	pCMSset.SetBorderSize(0);
//	pCMSset.AddText(set_names);
TLatex* tex_k = new TLatex(0.63,0.89,k_factor_str);
tex_k->SetNDC();
tex_k->SetTextAlign(20);
tex_k->SetTextFont(42);
tex_k->SetTextSize(0.03);
tex_k->SetLineWidth(2);

float axis2_array1[500];
float axis2_array2[500];
int ar1N;
int ar2N;

for (int i=0;i<nhistos;i++){
//	for (int i=0;i<4;i++){
		//temp_str.Form("%2.2f",Discr(discr_histos[i],signal_histos2[i]));
		Float_t d_value = Discr((TH1F*)data_histos[i]->Clone("data2"),signal_histos2[i]);
		temp_str.Form("%2.2f",d_value);
		TString disc_value = temp_str.Prepend(" d = ");
		TLatex *disc_value_text = new TLatex(0.668,0.86,disc_value);
      disc_value_text->SetNDC();
     	disc_value_text->SetTextAlign(20);
      disc_value_text->SetTextFont(42);
      disc_value_text->SetTextSize(0.03);
      disc_value_text->SetLineWidth(2);
		
		temp_str.Form("%2d",i);
		TString can_name="c1";
		can_name.Append(temp_str);
		TCanvas *c1 = new TCanvas(can_name,"",900,750);
		c1->cd();
		gPad->SetLogy(0);
		c1->SetBottomMargin(.3);
		c1->SetRightMargin(.25);
	
		bool LOGY=true;
		if ((region_type==0)&&(i==51)) LOGY=false;	
	//	bool LOGY=std::find(LOGY_array.begin(),LOGY_array.end(),i) != LOGY_array.end();
  	//	if ((region_type==0) && (LOGY==true)) gPad->SetLogy();
  		if ((region_type==0) && (LOGY==true)) gPad->SetLogy();
		
		
		Double_t xmin = signal_histos[i]->GetBinCenter(1)-signal_histos[i]->GetBinWidth(1)/2.;
		Double_t xmax = signal_histos[i]->GetBinCenter(signal_histos[i]->GetNbinsX())+signal_histos[i]->GetBinWidth(signal_histos[i]->GetNbinsX())/2.;
		if ((region_type!=0 ) && (i==66)) {xmin=0.;xmax=400;}
		if ((region_type!=0 ) && (i==87)) {xmin=-1.;xmax=1.;}
		if (hist_names[i].CompareTo("hPVs")==0) xmax=30;
		TH1F *frame = new TH1F("frame","",1,xmin,xmax);
		TGaxis::SetExponentOffset(-0.07,0,"xy");
		frame->Reset();
		frame->SetMinimum(1e-3);
      frame->SetMaximum(1e8);
		if ((region_type!=0) || (LOGY==false))	frame->SetMaximum(hBkgVis[i]->GetMaximum()*1.4);
      frame->GetXaxis()->SetTitleOffset(0.91);
      frame->SetStats(0);
		frame->GetYaxis()->SetNdivisions(505);
		if (hist_names[i].CompareTo("hSoft_n5")==0)  frame->GetXaxis()->SetNdivisions(18,-2);
	 	frame->GetXaxis()->SetLabelSize(0.0);
		char name[1000];
		if (UNITS[i]==0) {
			if (data_histos[i]->GetBinWidth(1)>1) sprintf(name,"Events / %1.0f",data_histos[i]->GetBinWidth(1));
			else 
				if (data_histos[i]->GetBinWidth(1)<0.01)  sprintf(name,"Events / %1.3f",data_histos[i]->GetBinWidth(1));
					else sprintf(name,"Events / %1.2f",data_histos[i]->GetBinWidth(1));
		} else 	sprintf(name,"Events / %1.0f",data_histos[i]->GetBinWidth(1),"GeV");
		frame->GetYaxis()->SetTitle(name);

      frame->Draw();
		tex->Draw();
	//	tex1->Draw();
//		tex2->Draw();
	pCMS1.Draw("same");
//	pCMSset.Draw("same");
	pCMS12.Draw("same");
		if (region_type==1) tex_k->Draw();
//		tex_k->Draw();
		if ((region_type==0)&&(d_value>0.1)&&(d_value<1.)) disc_value_text->Draw();
    	stacks[i]->Draw("same");	
//		stacks[i]->GetStack()->Last()->Draw("E2same"); 
		if  (region_type==0) signal_histos[i]->Draw("same");
	//	if  (region_type==0)tthbb_histos[i]->Draw("same");
	//	if  (region_type==0)tthnbb_histos[i]->Draw("same");
	//	if  (region_type==0)tthbb_sum_histos[i]->Draw("same");
	//	if  (region_type==0)gf_histos[i]->Draw("same");
/////////////////cross check of stupid THStack//////
//
//		sum_histos[i]->SetLineColor(kCyan);
//		sum_histos[i]->Scale(lumi);
//		sum_histos[i]->Draw("Lsame");
//
///////////////////////////////////////////////////
		data_histos[i]->Draw("Psame");
		hBkgVis[i]->Draw("same E2");
	//	if (hist_names[i].CompareTo("hMbb_regVBF_fsr")==0) 
		//	cout<<" Mbb  "<<data_histos[i]->GetBinContent(data_histos[i]->GetXaxis()->GetNbins()+1)<<"   , mc  "<<hBkgVis[i]->GetBinContent(hBkgVis[i]->GetXaxis()->GetNbins()+1) <<endl;

		

		leg->Draw("same");
	//	leg2->Draw("same");
  		gPad->RedrawAxis();
	
  		TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);
 		pad2->SetTopMargin(0.73);
 		pad2->SetRightMargin(0.25);
 		pad2->SetFillColor(0);
  		pad2->SetFillStyle(0);
  		pad2->Draw();
  		pad2->cd(0);
  		gPad->SetGridy();

		TH1F *frame2 = new TH1F("frame2","",1,xmin,xmax);
		frame2->SetMinimum(-1.);
      frame2->SetMaximum(1.);
      frame2->SetStats(0);
      frame2->SetTitleFont(42,"x");
		frame2->SetTitleFont(42,"y");
      frame2->SetTitleSize(0.13, "XYZ");
		frame2->GetYaxis()->SetNdivisions(505);
		if (hist_names[i].CompareTo("hSoft_n5")==0)  frame2->GetXaxis()->SetNdivisions(18,-2);
		if (hist_names[i].CompareTo("hSoft_n5")==0)  frame2->GetXaxis()->CenterLabels(kTRUE);
 		frame2->GetYaxis()->SetTickLength(0.06);
  		frame2->GetYaxis()->SetTitleSize(0.04);
  		frame2->GetYaxis()->SetTitleOffset(1.5);
 		frame2->GetYaxis()->SetLabelSize(0.03);
  		frame2->GetYaxis()->CenterTitle(kTRUE);
  		frame2->GetXaxis()->SetTitleSize(0.05);
  		frame2->GetXaxis()->SetLabelSize(0.04);
		frame2->SetXTitle(signal_histos[i]->GetXaxis()->GetTitle());
		frame2->SetYTitle("Data / MC - 1");
		frame2->Draw();	

 		Double_t aa[2] = {xmin,xmax};
   	Double_t bb[2] = {0,0};
   	TGraph *cons = new TGraph(2,aa,bb);
    	cons->SetLineStyle(2);
		cons->Draw("Lsame");

		data_histos2[i] = (TH1F*)data_histos[i]->Clone("new");
		data_histos2[i]->Add(sum_histos[i],-1);
		data_histos2[i]->Divide(sum_histos[i]);
		data_histos2[i]->Draw("PEsame");
		if (hist_names[i].CompareTo("hJet1q_axis2")==0) {
			ar1N=data_histos2[i]->GetXaxis()->GetNbins();
			for (int k=1;k<data_histos2[i]->GetXaxis()->GetNbins();k++)
					axis2_array1[k-1] = 1+data_histos2[i]->GetBinContent(k);
		}
		if (hist_names[i].CompareTo("hJet2q_axis2")==0) {
			ar2N=data_histos2[i]->GetXaxis()->GetNbins();
			for (int k=1;k<data_histos2[i]->GetXaxis()->GetNbins();k++)
					axis2_array2[k-1] = 1+ data_histos2[i]->GetBinContent(k);
		}
		hBkgUncUp[i]->Draw("HIST same");
		hBkgUncLo[i]->Draw("HIST same");
		pad2->RedrawAxis();
		c1->Print(output_names[i]);
		c1->Print(output_names2[i]);
		c1->Delete();
	}
		cout<<endl;
		cout<<endl;
		cout<<endl;
		cout<<endl;
		cout<<endl;

//		for (int k=1;k<ar1N;k++)
	//			cout<<axis2_array1[k] <<" , ";
//		cout<<endl;
//		for (int k=1;k<ar2N;k++)
//				cout<<axis2_array2[k]<<" , ";
		cout<<endl;
		cout<<endl;
		cout<<endl;


return 0;
}

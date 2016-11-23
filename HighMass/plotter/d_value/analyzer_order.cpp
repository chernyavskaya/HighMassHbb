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
int set_type = atoi(argv[1]); //0 - double, 1 - single
int region_type=atoi(argv[2]); // 0 - analysis, 1 - control region , top
 
const int nfiles  = 3;
TString leg_names[nfiles] = {"VBF, v14","VBF, v20", "VBF, v21"};


TString file_names[nfiles] = { "VBFHToBB_M-125_13TeV_powheg","VBFHToBB_M-125_13TeV_powheg",  "VBFHToBB_M-125_13TeV_powheg"}; 
int bg_begin;
int qcd_begin=10;
if (region_type==0) bg_begin=7;
else bg_begin=5;

	

TString btag[2] = {"DoubleB","SingleB"}; 
TString region[3]={"_analysis","_controlTop","_controlDY"}; // 0 - analysis, 1 - control region , top; 2 - control region DY

TString set[2] = {"_double","_single"};
for (int i=0;i<nfiles;i++){
	file_names[i].Prepend(set[set_type]);
	file_names[i].Prepend(region[region_type]);
	file_names[i].Prepend("skimmed_tree");
	file_names[i].Prepend("../../output_hist/v14/golden/");
	if (i==0) file_names[i].Append("_v14_final_74cmssw");
	if (i==1) file_names[i].Append("_v20");
	if (i==2) file_names[i].Append("_v21first");
//	file_names[i].Append("_v14_bdt");
//	if (i!=0) file_names[i].Append("_v14_qgl");
//	else file_names[i].Append("_v14_qgl2");
	file_names[i].Append(".root");
}
TString trigger[2] = {"RatioDoubleBtag_blike_", "RatioSingleBtag_blike_"};
//TString trigger[2] = {"trigWeightRatioDoubleBtag_", "trigWeightSilvioRatioSingleBtag_"};
//TString dir_name= "plots_powheg_130/";
TString dir_name= "plots_powheg_125_exclusive";
dir_name.Append(region[region_type]+"/");
//TString dir_name = "plots_amc/";
Float_t lumi = 2190;

	


TLegend *leg = new TLegend(0.77,0.80,0.92,0.9);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetTextSize(0.025);

const int nhistos = 70; //40//52
TString hist_names[nhistos]={"hJet1_pt","hJet2_pt","hJet3_pt","hJet4_pt","hJet1_eta","hJet2_eta","hJet3_eta","hJet4_eta","hJet1_phi","hJet2_phi","hJet3_phi","hJet4_phi","hMqq", "hEtaQQ", "hPhiBB","hMbb","hbtag","hbtag2","hcosOqqbb","hEtaQB1", "hEtaQB2", "hPhiQB1", "hPhiQB2","hx1","hx2","hVB1_mass","hVB2_mass","hEtot","hPxtot","hPytot","hPztot","hJet5_pt","hPtqqbb","hEtaqqbb","hPhiqqbb","hJet1_pt_bin","hJet2_pt_bin","hJet3_pt_bin","hJet4_pt_bin", "hMqq_bin","hEtaSoftJets", "hPtSoftJets","hMassSoftJets","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hqgl","hqgl2", "hPtSoftJets2","hPtSoftJets3","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult","hMbb_regVBF","hMbb_regVBF_fsr","hblike1","hblike2", "hmet",  "hPhiQQ", "hMbb_regVBF_fsr_high", "hMbb_regVBF_fsr_high_cat"};
int UNITS[nhistos]={1,1,1,1 ,0,0,0,0 ,0,0,0,0, 1 ,0,0 ,1, 0,0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1 ,0,0, 1,1,1,1,1, 0, 1,1,1 ,0,0,0,0,0, 1,1, 0, 1,0,0,0,0, 1,0,0,0,0, 1,1 ,0,0, 1, 0,1,1};
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
for (int i=0;i<nhistos;i++){
	output_names[i] = hist_names[i];
	output_names[i].Prepend("plots"+set[set_type]+"/VBF");
	output_names[i].Append(set[set_type]+"_CMSSW76_woweight.png");
}

TH1F *data_histos[nhistos];
TH1F *data_histos2[nhistos];
TH1F *signal_histos[nhistos];
TH1F *signal_histos2[nhistos];
TH1F *tthbb_histos[nhistos];
TH1F *tthnbb_histos[nhistos];
TH1F *gf_histos[nhistos];
TH1F *sum_histos[nhistos];
TH1F *histos_forUnc[nhistos];
TH1F *histos_for_legened[nhistos];
TH1F *discr_histos[nhistos];//Mqq,delta eta, delta phi, qgl, btag //12,13,14,21,22
TH1F *hBkgVis[nhistos];
TH1F *hBkgUncUp[nhistos];
TH1F *hBkgUncLo[nhistos];


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
out_efficiency.open(trigger[set_type]+dir_name+"efficiency.txt"); 

do{
   TFile *file_initial = new TFile(file_names[files]);
	TH1F *histos[100];
	for (int hist=0;hist<nhistos;++hist){
		if ((hist==71)||(hist==72)) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(2);

		if (files==0) {
			signal_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"newhist");
			signal_histos[hist]->Scale(1./signal_histos[hist]->Integral());
			signal_histos[hist]->Sumw2(kFALSE);
			signal_histos[hist]->SetLineColor(kBlue);
			signal_histos[hist]->SetLineStyle(2);
			signal_histos2[hist]=(TH1F*)signal_histos[hist]->Clone("signalHist2");
		}
		if (files==1) {
			gf_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"tthnb");
			gf_histos[hist]->Scale(1./gf_histos[hist]->Integral());
			gf_histos[hist]->Sumw2(kFALSE);
			gf_histos[hist]->SetLineColor(kRed);
			gf_histos[hist]->SetLineStyle(1);
		}
		if (files==2) {
			tthbb_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"tthnvv");
			tthbb_histos[hist]->Scale(1./tthbb_histos[hist]->Integral());
			tthbb_histos[hist]->Sumw2(kFALSE);
			tthbb_histos[hist]->SetLineColor(kBlack);
			tthbb_histos[hist]->SetLineStyle(1);
		}
	}
files++;
}while (files<nfiles);
out_efficiency.close();

leg->AddEntry(signal_histos[0],btag[set_type],"");
leg->AddEntry(signal_histos[0],leg_names[0],"L");
leg->AddEntry(gf_histos[0],leg_names[1],"L");
leg->AddEntry(tthbb_histos[0],leg_names[2],"L");




TLatex* tex = new TLatex(0.75,0.95,"13 TeV");
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
TLatex* tex2 = new TLatex(0.27,0.89,"Work in progress");
tex2->SetNDC();
tex2->SetTextAlign(20);
tex2->SetTextFont(52);
tex2->SetTextSize(0.035);
tex2->SetLineWidth(2);
TString temp_str;
for (int i=0;i<nhistos;i++){
//	for (int i=0;i<1;i++){
		
		temp_str.Form("%2d",i);
		TString can_name="c1";
		can_name.Append(temp_str);
		TCanvas *c1 = new TCanvas(can_name,"",900,750);
		c1->cd();
		gPad->SetLogy(0);
		c1->SetBottomMargin(.13);
		c1->SetRightMargin(.25);
	
		Double_t xmin = signal_histos[i]->GetBinCenter(0);
		Double_t xmax = signal_histos[i]->GetBinCenter(signal_histos[i]->GetNbinsX())*1.02;
		if ((region_type!=0 ) && (i==66)) {xmin=0.;xmax=400;}
		TH1F *frame = new TH1F("frame","",1,xmin,xmax);
		frame->Reset();
		frame->SetMinimum(1e-3);
      frame->SetMaximum(gf_histos[i]->GetMaximum()*1.2);
      frame->GetXaxis()->SetTitleOffset(0.91);
      frame->SetStats(0);
		frame->GetYaxis()->SetNdivisions(505);
	 	frame->GetXaxis()->SetLabelSize(0.0);
  		frame->GetXaxis()->SetTitleSize(0.05);
  		frame->GetXaxis()->SetLabelSize(0.04);
		frame->SetXTitle(signal_histos[i]->GetXaxis()->GetTitle());
		char name[1000];
		if (UNITS[i]==0) {
			if (signal_histos[i]->GetBinWidth(1)>1) sprintf(name,"Events / %1.0f",signal_histos[i]->GetBinWidth(1));
			else sprintf(name,"Events / %1.2f",signal_histos[i]->GetBinWidth(1));
		} else {
      	sprintf(name,"Events / %1.0f %s",signal_histos[i]->GetBinWidth(1),"GeV");
		}
		frame->GetYaxis()->SetTitle(name);

      frame->Draw();
		tex->Draw();
		tex1->Draw();
		tex2->Draw();
		signal_histos[i]->Draw("same");
		gf_histos[i]->Draw("same");
		tthbb_histos[i]->Draw("same");
		leg->Draw("same");
  		gPad->RedrawAxis();
	
		c1->Print(output_names[i]);
		c1->Delete();
	}



return 0;
}

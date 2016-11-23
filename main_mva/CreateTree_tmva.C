#define CreateTree_tmva_cxx
#include "CreateTree_tmva.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/preselection_highmass.C"
#include "math.h"


Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
const int njets = 300;

typedef struct {
	Float_t CSV1;
	Float_t CSV2;
	Float_t logCSV1;
	Float_t logCSV2;
	Float_t Mqq;
	Float_t Mbb;
	Float_t bb_pt;
	Float_t DeltaEtaQQ;
	Float_t DeltaPhiQQ;
	Int_t SoftN5;
	Int_t SoftN2;
	Float_t HTsoft;
	Float_t DeltaEtaQB1;
	Float_t DeltaEtaQB2;
	Float_t DeltaEtaQB;
	Float_t DeltaEtaQB_minus;
	Float_t cosOqqbb;
	Float_t qgl1_VBF;
	Float_t qgl2_VBF;
	Float_t x1;
	Float_t x2;
	Float_t VB1;
	Float_t VB2;
	Float_t Jet5_pt;
	Float_t Etot;
	Float_t axis2_jet1;
	Float_t axis2_jet2;
	Float_t qqbb_pt;
	Float_t qqbb_eta;
	Float_t qqbb_pz;
	Int_t jet1q_mult; 
	Int_t jet2q_mult; 
}TMVAstruct;


using namespace std;

void CreateTree_tmva::Loop(TString input_filename,TString output_dir, int data)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TMVAstruct TMVA;


	int events_saved=0;

	float weight;	
 

	TFile file("main_tmva_tree_"+output_dir+"_v21.root","recreate");
	TTree *tree0 = new TTree("TMVA","TMVA");
	tree0->Branch("CSV1",&TMVA.CSV1,"CSV1/F");
	tree0->Branch("CSV2",&TMVA.CSV2,"CSV2/F");
	tree0->Branch("logCSV1",&TMVA.logCSV1,"logCSV1/F");
	tree0->Branch("logCSV2",&TMVA.logCSV2,"logCSV2/F");
	tree0->Branch("Mqq",&TMVA.Mqq,"Mqq/F");
	tree0->Branch("Mbb",&TMVA.Mbb,"Mbb/F");
	tree0->Branch("bb_pt",&TMVA.bb_pt,"bb_pt/F");
	tree0->Branch("DeltaEtaQQ",&TMVA.DeltaEtaQQ,"DeltaEtaQQ/F");
	tree0->Branch("DeltaPhiQQ",&TMVA.DeltaPhiQQ,"DeltaPhiQQ/F");
	tree0->Branch("SoftN5",&TMVA.SoftN5,"SoftN5/I");
	tree0->Branch("SoftN2",&TMVA.SoftN2,"SoftN2/I");
	tree0->Branch("HTsoft",&TMVA.HTsoft,"HTsoft/F");
	tree0->Branch("axis2_jet1",&TMVA.axis2_jet1,"axis2_jet1/F");
	tree0->Branch("axis2_jet2",&TMVA.axis2_jet2,"axis2_jet2/F");
	tree0->Branch("DeltaEtaQB1",&TMVA.DeltaEtaQB1,"DeltaEtaQB1/F");
	tree0->Branch("DeltaEtaQB2",&TMVA.DeltaEtaQB2,"DeltaEtaQB2/F");
	tree0->Branch("DeltaEtaQB",&TMVA.DeltaEtaQB,"DeltaEtaQB/F");
	tree0->Branch("DeltaEtaQB_minus",&TMVA.DeltaEtaQB_minus,"DeltaEtaQB_minus/F");
	tree0->Branch("cosOqqbb",&TMVA.cosOqqbb,"cosOqqbb/F");
	tree0->Branch("qqbb_pt",&TMVA.qqbb_pt,"qqbb_pt/F");
	tree0->Branch("qqbb_eta",&TMVA.qqbb_eta,"qqbb_eta/F");
	tree0->Branch("qqbb_pz",&TMVA.qqbb_pz,"qqbb_pz/F");
//	tree0->Branch("qgl1_VBF",&TMVA.qgl1_VBF,"qgl1_VBF/F");
//	tree0->Branch("qgl2_VBF",&TMVA.qgl2_VBF,"qgl2_VBF/F");
	tree0->Branch("jet1q_mult",&TMVA.jet1q_mult,"jet1q_mult/I");
	tree0->Branch("jet2q_mult",&TMVA.jet2q_mult,"jet2q_mult/I");
	tree0->Branch("Jet5_pt",&TMVA.Jet5_pt,"Jet5_pt/F");
	tree0->Branch("x1",&TMVA.x1,"x1/F");
	tree0->Branch("x2",&TMVA.x2,"x2/F");
//	tree0->Branch("VB1",&TMVA.VB1,"VB1/F");
//	tree0->Branch("VB2",&TMVA.VB2,"VB2/F");
	tree0->Branch("weight",&weight);


	TF1 *func_r = new TF1("erffunc",erf,0.,1000.,3);
	
	cout<<nentries<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;

		if (data==1) if ((jentry%50!=0) && (jentry%49!=0)&& (jentry%48!=0)) continue;

		if (data==1) genWeight = 1.;
		if (genWeight <0) continue;
		if (json!=1) continue;

		weight = bTagWeight*puWeight*TMath::Abs(genWeight)/genWeight;

		if (!((Vtype==-1)||(Vtype>3))) continue;

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
		TLorentzVector bb;
		vector<TLorentzVector> FSRjet;
		///////////////////////
		//preselection/////
		//////////////////////
		
		
		if  (preselection(nJet, Jet_pt,Jet_eta, Jet_phi, Jet_mass, Jet_btagCSV, Jet_id,Jet_puId, btag_max1_number, btag_max2_number, eta_max1_number, eta_max2_number,HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v, Bjet1, Bjet2, Qjet1, Qjet2, qq,jets_pv, jets_btag, jets_pt, good_jets, jets_indices, met_pt, met_phi, met, FSRjet) != 0)  continue; 

		bb = Bjet1+Bjet2;

		Float_t Jet5_pt = -1;
		for (int i=0;i<good_jets;i++){
			if ((i==btag_max1_number)||(i==btag_max2_number)||(i==eta_max1_number)||(i==eta_max2_number)||(Jet_id[i]<=2)||(Jet_puId[i]==0)) continue;
			Jet5_pt=jets_pv[i].Pt();
			break;
		}

		for (int i=0;i<FSRjet.size();i++)
			bb+= FSRjet[i];

		bb+=met;

	
		///////TRIGGER WEIGHTS////////////////////////// 
		float trigCor0_nom_double[10] = {1.,1.,1., 1. , 1., 1.31914e+00  };
		float trigCor1_nom_double[10] = {7.60461e+01,  6.92756e+01 , 4.90344e+01,-9.29334e+00,-4.39848e+01, -6.02564e+01  };
		float trigCor2_nom_double[10] = { 1.00000e+02,  2.05669e+01 ,2.29954e+01,3.86211e+01,1.01285e+01,1.41746e+01  };

		float trigWeight=1;
		func_r->FixParameter(0,1.);
		int whichTrigWeight =0;		

		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[0] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[0] );
			trigWeight*=func_r->Eval(jets_pv[0].Pt());
		}
		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[1] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[1] );
			trigWeight*=func_r->Eval(jets_pv[1].Pt());
		}
		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[2] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[2] );
			trigWeight*=func_r->Eval(jets_pv[2].Pt());
		}

		//////////////////////forth jet, not used////////////////////
/*		if (data!=1) {
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[3] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[3] );
			trigWeight*=func_r->Eval(jets_pv[3].Pt());
		}
*/

	if (data!=1){
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[4] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[4] );
			trigWeight*=func_r->Eval(-1.*TMath::Log(1.-Jet_btagCSV[jets_indices[btag_max1_number]]));
		}
		if (data!=1){
			if (whichTrigWeight==0) func_r->FixParameter(0,trigCor0_nom_double[5] );
			if (whichTrigWeight==0) func_r->FixParameter(1,trigCor1_nom_double[5] );
			if (whichTrigWeight==0) func_r->FixParameter(2,trigCor2_nom_double[5] );
			trigWeight*=func_r->Eval(-1.*TMath::Log(1.-Jet_btagCSV[jets_indices[btag_max2_number]]));
		}


	weight*=trigWeight;
		
/////////////////////////////////////////////////////////


		Float_t Mbb = bb.M();

		TLorentzVector bbqq;
		bbqq = bb  + qq;

		Float_t cosObbqq =TMath::Cos( ( ( Bjet1.Vect() ).Cross(Bjet2.Vect()) ).Angle( ( Qjet1.Vect() ).Cross(Qjet2.Vect()) ) );	

		Float_t EtaBQ1;
	 	Float_t EtaBQ2;
		Float_t PhiBQ1; 	
		Float_t PhiBQ2;

		
		Float_t Mqq = qq.M();
		Float_t bbDeltaPhi = TMath::Abs(Bjet1.DeltaPhi(Bjet2));
		Double_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t cosOqqbb =TMath::Cos( ( ( Bjet1.Vect() ).Cross(Bjet2.Vect()) ).Angle( ( Qjet1.Vect() ).Cross(Qjet2.Vect()) ) );	

		
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


		for (int i=0;i<nJet;i++){
			if (Jet_btagCSV[i]>1) Jet_btagCSV[i]=1.;
			if (Jet_btagCSV[i]<0) Jet_btagCSV[i]=0.;
		}

		TMVA.Mqq = Mqq;
		TMVA.CSV1 = Jet_btagCSV[jets_indices[btag_max1_number]];	
		TMVA.CSV2 = Jet_btagCSV[jets_indices[btag_max2_number]];
		TMVA.logCSV1 = -1.*TMath::Log(1.-Jet_btagCSV[jets_indices[btag_max1_number]]);	
		TMVA.logCSV2 = -1.*TMath::Log(1.-Jet_btagCSV[jets_indices[btag_max2_number]]);
		TMVA.DeltaEtaQQ = qqDeltaEta;
		TMVA.DeltaPhiQQ = qqDeltaPhi;
		TMVA.SoftN5 = softActivity_njets5;
		TMVA.SoftN2 = softActivity_njets2;
		TMVA.HTsoft = softActivity_HT;
		TMVA.DeltaEtaQB1 = EtaBQ1;
		TMVA.DeltaEtaQB2 = EtaBQ2;
		TMVA.DeltaEtaQB = (EtaBQ1 + EtaBQ2);
		TMVA.cosOqqbb = cosOqqbb;
		TMVA.Jet5_pt = Jet5_pt;
		TMVA.x1 = x1;
		TMVA.x2 = x2;
		TMVA.Mbb = Mbb;
		TMVA.axis2_jet1 = Jet_axis2[jets_indices[eta_max1_number]];
		TMVA.axis2_jet2 = Jet_axis2[jets_indices[eta_max2_number]];
		TMVA.bb_pt = bb.Pt();
		TMVA.DeltaEtaQB_minus = (EtaBQ1 - EtaBQ2);
		TMVA.qqbb_pt = (bb +qq).Pt();
		TMVA.qqbb_eta = (bb +qq).Eta();
		TMVA.qqbb_pz = (bb+qq).Pz();
		TMVA.jet1q_mult = Jet_mult[jets_indices[eta_max1_number]];
		TMVA.jet2q_mult = Jet_mult[jets_indices[eta_max2_number]];




		tree0->Fill();	
		events_saved++;		
	 	if ((data==1)&& (events_saved>=100000)) break;
		if ((data==0)&&(events_saved>=10000)) break; 

	}  

	file.Write();
	file.Close();

}


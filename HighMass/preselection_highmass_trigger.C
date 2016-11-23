#include <vector>
#include <TLorentzVector.h>

int preselection_double(Int_t nJets, Float_t Jet_pt[300], Float_t Jet_eta[300], Float_t Jet_phi[300], Float_t Jet_mass[300], Float_t Jet_btagCSV[300], Int_t Jet_id[300],Int_t Jet_puId[300], Int_t& btag_max1_number, Int_t& btag_max2_number, Int_t& eta_max1_number, Int_t& eta_max2_number, Float_t trigger, TLorentzVector& Bjet1,TLorentzVector& Bjet2, TLorentzVector& Qjet1, TLorentzVector& Qjet2,TLorentzVector& qq, Int_t& cut_count, std::vector<TLorentzVector> &jets_pv, std::vector<float> &jets_btag){
	
	btag_max1_number = -1;
	btag_max2_number = -1;
	eta_max1_number = -1;
	eta_max2_number = -1;
	cut_count = 0;
	int not_pass=0;


		int good_jets = 0;
		std::vector<float> jets_pt;
		std::vector<float> jets_btag_times_pt;
		for (int i=0;i<nJets;i++){
			TLorentzVector jet0;
			if (!((Jet_id[i]>2)&&(Jet_puId[i]>0))) continue;
			jet0.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
			jets_pv.push_back(jet0);
			good_jets++;
			jets_btag.push_back(Jet_btagCSV[i]);
			jets_pt.push_back(Jet_pt[i]);
			jets_btag_times_pt.push_back(Jet_btagCSV[i]*Jet_pt[i]);
		}

		if (good_jets<4) return -1;


		if (!(jets_pv[0].Pt()>92.)) {not_pass= -2; cut_count++;}
		if (!(jets_pv[1].Pt()>76.)) {not_pass= -3; cut_count++;}
		if (!(jets_pv[2].Pt()>64.)) {not_pass= -4; cut_count++;}
		if (!(jets_pv[3].Pt()>30.))  {not_pass= -5; cut_count++;}

      btag_max1_number =  max_element(jets_btag_times_pt.begin(), jets_btag_times_pt.end() )-  jets_btag_times_pt.begin();
    	jets_btag_times_pt[btag_max1_number] = 0; 
      btag_max2_number =   max_element( jets_btag_times_pt.begin(), jets_btag_times_pt.end()) -  jets_btag_times_pt.begin();

		Bjet1 = jets_pv[btag_max1_number];
		Bjet2 = jets_pv[btag_max2_number];
		double maxdEta = -1;
		for (int i=0;i<std::min(7,good_jets);i++){
			if ((i==btag_max1_number) || (i==btag_max2_number)) continue;
			for (int j=i+1;j<std::min(7,good_jets);j++){
				if ((j==btag_max1_number) || (j==btag_max2_number)) continue;
				double d = TMath::Abs(jets_pv[i].Eta() - jets_pv[j].Eta());
				if (d> maxdEta) {
					maxdEta = d;
					eta_max1_number = i;
					eta_max2_number = j;
				}
			}
		}
		Qjet1 = jets_pv[eta_max1_number];
		Qjet2 = jets_pv[eta_max2_number];
			

		qq = Qjet1 + Qjet2;
		TLorentzVector bb = Bjet1 + Bjet2;	
		Float_t Mbb = bb.M();
		int fsr_bjet=0;
		std::vector<TLorentzVector> FSRjet;
		for (int i=0;i<good_jets;i++){	
			if ((i!=btag_max1_number)&&(i!=btag_max2_number)&&(i!=eta_max1_number)&&(i!=eta_max2_number)){
				TLorentzVector jet = jets_pv[i]; 
				if ((jet.DeltaR(Bjet1)<0.8) || (jet.DeltaR(Bjet2)<0.8)) {
					fsr_bjet++;
					FSRjet.push_back(jet);
				} 
			}
		}
		TLorentzVector bb_fsr = Bjet1+Bjet2;
		for (int i=0;i<FSRjet.size();i++)
			bb_fsr+= FSRjet[i];
		 float Mbb_fsr=bb_fsr.M();


		Float_t bbDeltaPhi = TMath::Abs(Bjet1.DeltaPhi(Bjet2));
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));
		Float_t Mqq = qq.M();	

		if (!(Mqq>200)) {not_pass=-6; cut_count++;}
		if (!(qqDeltaEta>1.2)) {not_pass = -7; cut_count++;}
		if (!(Mbb_fsr>200)) {not_pass=-8; cut_count++;}
		if (trigger!=1) not_pass =  -9;

	return not_pass;
}

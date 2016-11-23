import ROOT
from ROOT import gROOT
from ROOT import gStyle

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]



gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/workflow/VBFSpin0ToBBbar_W_1p0_M_1050_TuneCUEP8M1_13TeV_pythia8_v21.root")
#f = ROOT.TFile.Open('/afs/cern.ch/work/n/nchernya/tmva/TMVA-v4.2.0/test/HighMassQlikelihood/treeWithQlike750GeV.root')
tree = f.Get("tree")

hist_count = f.Get('CountWeighted')
events_generated = hist_count.GetEntries()

hist_pt_ord_pt = []
hist_pt_ord_eta = []
hist_eta_ord_pt = []
hist_eta_ord_eta = []
hist_flavour_ord_pt = []
hist_flavour_ord_eta = []

hist_num = 12 
hist_mqq = ROOT.TH1F("mqq","m_{qq} (GeV)",150,0,3000)
hist_etaqq = ROOT.TH1F("etaqq","#eta_{qq}",16,0,10)
hist_etabb = ROOT.TH1F("etabb","#eta_{bb}",16,0,10)
hist_etaqqbb = ROOT.TH1F("etaqqbb","#Delta#eta_{qb}^{forward}+#Delta#eta_{qb}^{backward}",160,-8,8)
hist_phiqq = ROOT.TH1F("phiqq","|#phi_{qq}|",32,0,3.2)
hist_phibb = ROOT.TH1F("phibb","|#phi_{bb}|",32,0,3.2)
hist_ptqqbb = ROOT.TH1F("ptqqbb","p_{T} qqbb (GeV)",30,0,150)
hist_ptqq = ROOT.TH1F("ptqq","p_{T} qq (GeV)",50,0,500)
hist_ptbb = ROOT.TH1F("ptbb","p_{T} bb (GeV)",50,0,500)
hist_pzqqbb = ROOT.TH1F("pzqqbb","p_{z} qqbb (GeV)",200,-4000,4000)
hist_pzqq = ROOT.TH1F("pzqq","p_{z} qq (GeV)",100,-2000,2000)
hist_pzbb = ROOT.TH1F("pzbb","p_{z} bb (GeV)",100,-2000,2000)

for i in range(0,4):
	hist_pt_ord_pt.append(ROOT.TH1F("pt_ord_pt_%i" % i,"pt_ord_pt_%i" % i,120,0,600))
	hist_pt_ord_eta.append(ROOT.TH1F("pt_ord_eta_%i" % i,"pt_ord_eta_%i" % i,32,-8,8))
	hist_eta_ord_pt.append(ROOT.TH1F("eta_ord_pt_%i" % i,"eta_ord_pt_%i" % i,120,0,600))
	hist_eta_ord_eta.append(ROOT.TH1F("eta_ord_eta_%i" % i,"eta_ord_eta_%i" % i,32,-8,8))
	hist_flavour_ord_pt.append(ROOT.TH1F("flavour_ord_pt_%i" % i,"flavour_ord_pt_%i" % i,120,0,600))
	hist_flavour_ord_eta.append(ROOT.TH1F("flavour_ord_eta_%i" % i,"flavour_ord_eta_%i" % i,32,-8,8))

triggers = {'DoubleB200':0, 'DoubleB240':0,'SingleB460':0, 'SingleB500':0, 'SingleBorDoubleB200_460':0,'SingleBorDoubleB240_500':0}

presel = 0
presel_trigger = 0
eff_choosing_jets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_bjets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_qjets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_jets_csv_pt = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_jets_csv_times_pt_pt = {'b1':0,'b2':0,'q1':0,'q2':0,'bb':0,'qq':0}
eff_choosing_jets_csv_times_pt_deta = {'b1':0,'b2':0,'q1':0,'q2':0,'bb':0,'qq':0}
eff_choosing_jets_csv_times_pt_qlike = {'b1':0,'b2':0,'q1':0,'q2':0,'bb':0,'qq':0}
eff_choosing_jets_csv_deta = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_ref = 0

hist_mqq_pt = ROOT.TH1F("mqq_pt","m_{qq} (GeV), chosen by p_{T}",150,0,3000)
hist_mqq_deta = ROOT.TH1F("mqq_deta","m_{qq} (GeV), chosen by #Delta#eta",150,0,3000)


num_jets_considered = 5
for entry in range(1,tree.GetEntries()):
#for entry in range(1,20):
	if entry%10000==0 : print entry
	tree.GetEntry(entry)
		
#	if tree.json != 1 : continue
#	weight = 1
	
	weight = abs(tree.genWeight)/tree.genWeight * tree.bTagWeight * tree.puWeight

#	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==1) and (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v==0) : triggers['DoubleB200']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==1) : triggers['DoubleB200']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v==1) and (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v==0) : triggers['DoubleB240']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==0) and (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v==1) : triggers['SingleB460']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v==0) and (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v==1) : triggers['SingleB500']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v==1) or (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v==1) : triggers['SingleBorDoubleB240_500']+=weight
	if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==1) or (tree.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v==1) : triggers['SingleBorDoubleB200_460']+=weight

	jets = []
	jets.append(ROOT.TLorentzVector())
	jets[0].SetPtEtaPhiM(tree.GenBQuarkFromH_pt[0],tree.GenBQuarkFromH_eta[0],tree.GenBQuarkFromH_phi[0],tree.GenBQuarkFromH_mass[0])		
	jets.append(ROOT.TLorentzVector())
	jets[1].SetPtEtaPhiM(tree.GenBQuarkFromH_pt[1],tree.GenBQuarkFromH_eta[1],tree.GenBQuarkFromH_phi[1],tree.GenBQuarkFromH_mass[1])	
	jets.append(ROOT.TLorentzVector())	
	jets[2].SetPtEtaPhiM(tree.GenHiggsSisters_pt[0],tree.GenHiggsSisters_eta[0],tree.GenHiggsSisters_phi[0],tree.GenHiggsSisters_mass[0])	
	jets.append(ROOT.TLorentzVector())
	jets[3].SetPtEtaPhiM(tree.GenHiggsSisters_pt[1],tree.GenHiggsSisters_eta[1],tree.GenHiggsSisters_phi[1],tree.GenHiggsSisters_mass[1])	

	good_jets = []
	good_jets_array_idx = []
	num_good_jets = 0 
	for i in range (0,tree.nJet):
		if (tree.Jet_puId[i] > 0) and (tree.Jet_id[i] > 2) :
			good_jets.append(ROOT.TLorentzVector())
			good_jets[num_good_jets].SetPtEtaPhiM(tree.Jet_pt[i],tree.Jet_eta[i],tree.Jet_phi[i],tree.Jet_mass[i])
			num_good_jets+=1
			good_jets_array_idx.append(i)
	
	if len(good_jets) >= 4 :
		if (good_jets[0].Pt() > 92) and (good_jets[1].Pt() >76 ) and (good_jets[2].Pt() > 64) and (good_jets[3].Pt() >30 ) and ((tree.Vtype==-1) or (tree.Vtype>3)) :
			btag = []
			btag_pt = []
			pt = []
			qlike = []
			for i in range(0,min(len(good_jets),15)) :
				btag.append(tree.Jet_btagCSV[good_jets_array_idx[i]])
				btag_pt.append(tree.Jet_btagCSV[good_jets_array_idx[i]]*good_jets[i].Pt())
				pt.append(tree.Jet_pt[good_jets_array_idx[i]])
		#		qlike.append(tree.Jet_qlikelihood[i])
			chosen_bjet1 = ROOT.TLorentzVector()	
			chosen_bjet2 = ROOT.TLorentzVector()
			chosen_bjet1_csvpt = ROOT.TLorentzVector()	
			chosen_bjet2_csvpt = ROOT.TLorentzVector()
			chosen_bjet1 = good_jets[btag.index(max(btag))]
			bjet1_idx = btag.index(max(btag)) 
			pt[btag.index(max(btag))] = -99
			btag[btag.index(max(btag))] = -99
			chosen_bjet2 = good_jets[btag.index(max(btag))]
			bjet2_idx = btag.index(max(btag))
			bjets_idc = bjet1_idx, bjet2_idx 
			jets_for_deta = [i for j, i in enumerate(good_jets) if j not in bjets_idc]
			pt[btag.index(max(btag))] = -99
			btag[btag.index(max(btag))] = -99
			chosen_qjet1 = ROOT.TLorentzVector()	
			chosen_qjet2 = ROOT.TLorentzVector()
			chosen_qjet1 = good_jets[pt.index(max(pt))]
			pt[pt.index(max(pt))] = -99
			chosen_qjet2 = good_jets[pt.index(max(pt))]
			if ((chosen_qjet1+chosen_qjet2).M() > 200) and (abs(chosen_qjet1.Eta()-chosen_qjet2.Eta())>1.2) : presel+=weight	
			if ((chosen_qjet1+chosen_qjet2).M() > 200) and (abs(chosen_qjet1.Eta()-chosen_qjet2.Eta())>1.2) and (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==1): presel_trigger+=weight	

			
			chosen_bjet1_csvpt = good_jets[btag_pt.index(max(btag_pt))]
			bjet1_idx_csvpt = btag_pt.index(max(btag_pt)) 
			btag_pt[btag_pt.index(max(btag_pt))] = 0
			chosen_bjet2_csvpt = good_jets[btag_pt.index(max(btag_pt))]
			bjet2_idx_csvpt = btag_pt.index(max(btag_pt)) 
			bjets_idc_csvpt = bjet1_idx_csvpt, bjet2_idx_csvpt
			jets_for_csvpt_pt = [i for j, i in enumerate(good_jets) if j not in bjets_idc_csvpt]		
	#		jets_for_csvpt_qlike = [i for j, i in enumerate(good_jets) if j not in bjets_idc_csvpt]		
			chosen_qjet1_csvpt = jets_for_csvpt_pt[0]
			chosen_qjet2_csvpt = jets_for_csvpt_pt[1]

			chosen_qq_csvpt = chosen_qjet1_csvpt + chosen_qjet2_csvpt

	# event interpretation part
	
			if (tree.HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v==1):
				eff_choosing_ref+=weight
				if (chosen_bjet1.DeltaR(jets[0])<0.8) or (chosen_bjet1.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_pt['b1']+=weight
				if (chosen_bjet2.DeltaR(jets[0])<0.8) or (chosen_bjet2.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_pt['b2']+=weight
				if (chosen_qjet1.DeltaR(jets[2])<0.8) or (chosen_qjet1.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_pt['q1']+=weight
				if (chosen_qjet2.DeltaR(jets[2])<0.8) or (chosen_qjet2.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_pt['q2']+=weight

				if (chosen_bjet1.DeltaR(jets[0])<0.8) or (chosen_bjet1.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_deta['b1']+=weight
				if (chosen_bjet2.DeltaR(jets[0])<0.8) or (chosen_bjet2.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_deta['b2']+=weight

				maxEta1, maxEta2 = None, None
				maxDeltaEta = -1
				chosen_qjet1_deta = ROOT.TLorentzVector()
				chosen_qjet2_deta = ROOT.TLorentzVector()
				for i in range(min(num_jets_considered,len(jets_for_deta))):
					for j in range(i + 1, min(num_jets_considered,len(jets_for_deta))):
						d = abs(jets_for_deta[i].Eta() - jets_for_deta[j].Eta())
						if d > maxDeltaEta:	
							maxEta1, maxEta2, maxDeltaEta = i, j, d
				if (maxEta1>=0) and (maxEta2>=0) :
					chosen_qjet1_deta = jets_for_deta[maxEta1]		
					chosen_qjet2_deta = jets_for_deta[maxEta2]

				maxEta1, maxEta2 = None, None
				maxDeltaEta = -1
				chosen_qjet1_csvpt_deta = ROOT.TLorentzVector()
				chosen_qjet2_csvpt_deta = ROOT.TLorentzVector()
				for i in range(min(num_jets_considered,len(jets_for_csvpt_pt))):
				#	print jets_for_csvpt_pt[i].Eta()
					for j in range(i + 1, min(num_jets_considered,len(jets_for_csvpt_pt))):
						d = abs(jets_for_csvpt_pt[i].Eta() - jets_for_csvpt_pt[j].Eta())
						if d > maxDeltaEta:	
							maxEta1, maxEta2, maxDeltaEta = i, j, d
				if (maxEta1>=0) and (maxEta2>=0) :
					chosen_qjet1_csvpt_deta = jets_for_csvpt_pt[maxEta1]		
					chosen_qjet2_csvpt_deta = jets_for_csvpt_pt[maxEta2]
				#	print chosen_bjet1_csvpt.Pt(), chosen_bjet2_csvpt.Pt()
				#	print chosen_qjet1_csvpt_deta.Eta(), chosen_qjet2_csvpt_deta.Eta()

				chosen_qq_csvpt_deta = chosen_qjet1_csvpt_deta + chosen_qjet2_csvpt_deta
				hist_mqq_deta.Fill(chosen_qq_csvpt_deta.M(),weight)
				hist_mqq_pt.Fill(chosen_qq_csvpt.M(),weight)


	#			chosen_qjet1_csvpt_qlike = ROOT.TLorentzVector()
#				chosen_qjet2_csvpt_qlike = ROOT.TLorentzVector()
#				good_jets_for_qlike = [i for j, i in enumerate(good_jets) if j not in bjets_idc_csvpt]
#				qlike_clean = [i for j, i in enumerate(qlike) if j not in bjets_idc_csvpt]
#				qlike_jets = good_jets_for_qlike[:min(len(good_jets_for_qlike),4)]
#				chosen_qjet1_csvpt_qlike = qlike_jets[qlike_clean.index(max(qlike_clean))]
#				qlike_clean[qlike_clean.index(max(qlike_clean))] = -3
#				chosen_qjet2_csvpt_qlike = qlike_jets[qlike_clean.index(max(qlike_clean))]
					
					

				if (chosen_qjet1_deta.DeltaR(jets[2])<0.8) or (chosen_qjet1_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_deta['q1']+=weight
				if (chosen_qjet2_deta.DeltaR(jets[2])<0.8) or (chosen_qjet2_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_deta['q2']+=weight

				if (chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_pt['b1']+=weight
				if (chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_pt['b2']+=weight
				if (chosen_qjet1_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_pt['q1']+=weight
				if (chosen_qjet2_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_pt['q2']+=weight
				if ((chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8)) and ((chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) ) : eff_choosing_jets_csv_times_pt_pt['bb']+=weight
				if ((chosen_qjet1_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt.DeltaR(jets[3])<0.8)) and ((chosen_qjet2_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt.DeltaR(jets[3])<0.8) ) : eff_choosing_jets_csv_times_pt_pt['qq']+=weight
				
				if (chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_deta['b1']+=weight
				if (chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_deta['b2']+=weight
				if (chosen_qjet1_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_deta['q1']+=weight
				if (chosen_qjet2_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_deta['q2']+=weight
				if ((chosen_qjet1_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt_deta.DeltaR(jets[3])<0.8)) and ((chosen_qjet2_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt_deta.DeltaR(jets[3])<0.8) ) : eff_choosing_jets_csv_times_pt_deta['qq']+=weight

		#		if (chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_qlike['b1']+=weight
		#		if (chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_qlike['b2']+=weight
		#		if (chosen_qjet1_csvpt_qlike.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt_qlike.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_qlike['q1']+=weight
	#			if (chosen_qjet2_csvpt_qlike.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt_qlike.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_qlike['q2']+=weight
#				if ((chosen_qjet1_csvpt_qlike.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt_qlike.DeltaR(jets[3])<0.8)) and ((chosen_qjet2_csvpt_qlike.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt_qlike.DeltaR(jets[3])<0.8) ) : eff_choosing_jets_csv_times_pt_qlike['qq']+=weight


				for i,jet in enumerate(good_jets):
					if (jet.DeltaR(jets[0])<0.8) or (jet.DeltaR(jets[1])<0.8) or (jet.DeltaR(jets[2])<0.8) or (jet.DeltaR(jets[3])<0.8) : eff_choosing_jets_pt[str(i+1)]+=weight
					if (jet.DeltaR(jets[0])<0.8) or (jet.DeltaR(jets[1])<0.8)  : eff_choosing_bjets_pt[str(i+1)]+=weight
					if (jet.DeltaR(jets[2])<0.8) or (jet.DeltaR(jets[3])<0.8)  : eff_choosing_qjets_pt[str(i+1)]+=weight	
					if i==6 : break
			
				
			
			
						





	pteta_jets = [[jets[0].Pt(),jets[0].Eta()], [jets[1].Pt(),jets[1].Eta()], [jets[2].Pt(),jets[2].Eta()],[jets[3].Pt(),jets[3].Eta()]]
	jets_ptordered = sorted(pteta_jets,key=getPt)
	for i in  range(len(jets_ptordered)):
		hist_pt_ord_pt[i].Fill(jets_ptordered[i][0],weight)
		hist_pt_ord_eta[i].Fill(jets_ptordered[i][1],weight)
	jets_etaordered = sorted(pteta_jets,key=getEta)
	for i in  range(len(jets_etaordered)):
		hist_eta_ord_pt[i].Fill(jets_etaordered[i][0],weight)
		hist_eta_ord_eta[i].Fill(jets_etaordered[i][1],weight)
	for i in  range(len(jets)):
		hist_flavour_ord_pt[i].Fill(jets[i].Pt(),weight)
		hist_flavour_ord_eta[i].Fill(jets[i].Eta(),weight)
	qq_jet = jets[2] + jets[3]
	bb_jet = jets[0] + jets[1]
	qqbb_jet = jets[0] + jets[1] + jets[2] + jets[3]
	hist_mqq.Fill(qq_jet.M(),weight)
	hist_etaqq.Fill(abs(jets[2].Eta() - jets[3].Eta()),weight)
	hist_etabb.Fill(abs(jets[0].Eta() - jets[1].Eta()),weight)
	hist_phibb.Fill(abs(jets[0].DeltaPhi(jets[1])),weight)
	hist_phiqq.Fill(abs(jets[2].DeltaPhi(jets[3])),weight)
	hist_ptqq.Fill(qq_jet.Pt(),weight)
	hist_pzqq.Fill(qq_jet.Pz(),weight)
	hist_ptbb.Fill(bb_jet.Pt(),weight)
	hist_pzbb.Fill(bb_jet.Pz(),weight)
	hist_ptqqbb.Fill(qqbb_jet.Pt(),weight)
	hist_pzqqbb.Fill(qqbb_jet.Pz(),weight)
	hist_etaqqbb.Fill( (jets[2].Eta() + jets[3].Eta() - jets[0].Eta() - jets[1].Eta())  ,weight)

	




print 'only DoubleB200 = %.3f  '%(triggers['DoubleB200']/events_generated), 'presel + doubleB200 = ', presel_trigger/events_generated,'presel = ', presel/events_generated , 'DoubleB240 = ', triggers['DoubleB240']/events_generated, 'SingleB460 = ', triggers['SingleB460']/events_generated, 'SingleB500 = ', triggers['SingleB500']/events_generated, 'SingleBorDoubleB200_460 = ', triggers['SingleBorDoubleB200_460']/events_generated, 'SingleBorDoubleB240_500 = ', triggers['SingleBorDoubleB240_500']/events_generated,'\n\n'

#for i in range(0,7):
#	print 'matching eff all jets for %i jet = '%(i+1), eff_choosing_jets_pt[str(i+1)]/eff_choosing_ref,'\n'
#for i in range(0,7):
#	print 'matching eff b-jets for %i jet = '%(i+1), eff_choosing_bjets_pt[str(i+1)]/eff_choosing_ref,'\n'
#for i in range(0,7):
#	print 'matching eff q-jets for %i jet = '%(i+1), eff_choosing_qjets_pt[str(i+1)]/eff_choosing_ref,'\n'

print 'matching efficiency for csv + pt for number of jets considered %i : \n'%num_jets_considered
print 'b1 =  %.3f '%(eff_choosing_jets_csv_pt['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_pt['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_pt['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_pt['q2']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv + deta for number of jets considered %i : \n'%num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_deta['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_deta['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_deta['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_deta['q2']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv*pt and pt for number of jets considered %i : \n'%num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['q2']/eff_choosing_ref), '\n'
print 'bb = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['bb']/eff_choosing_ref), '\n'
print 'qq = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['qq']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv*pt and deta for number of jets considered %i : \n'%num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['q2']/eff_choosing_ref), '\n'
print 'qq = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['qq']/eff_choosing_ref), '\n\n'
#print 'matching efficiency for csv*pt and qlikelihood: \n',num_jets_considered
#print 'b1 = %.3f  '%(eff_choosing_jets_csv_times_pt_qlike['b1']/eff_choosing_ref), '\n'
#print 'b2 = %.3f  '%(eff_choosing_jets_csv_times_pt_qlike['b2']/eff_choosing_ref), '\n'
#print 'q1 = %.3f  '%(eff_choosing_jets_csv_times_pt_qlike['q1']/eff_choosing_ref), '\n'
#print 'q2 = %.3f  '%(eff_choosing_jets_csv_times_pt_qlike['q2']/eff_choosing_ref), '\n'
#print 'qq = %.3f  '%(eff_choosing_jets_csv_times_pt_qlike['qq']/eff_choosing_ref), '\n\n'
	
all_hist = []
all_hist.append(hist_pt_ord_pt)
all_hist.append(hist_pt_ord_eta)
all_hist.append(hist_eta_ord_pt)
all_hist.append(hist_eta_ord_eta)
all_hist.append(hist_flavour_ord_pt)
all_hist.append(hist_flavour_ord_eta)

for i in range(0,6):
	all_hist[i][0].SetLineColor(ROOT.kMagenta)
	all_hist[i][1].SetLineColor(ROOT.kBlue)
	all_hist[i][2].SetLineColor(ROOT.kGreen)
	all_hist[i][3].SetLineColor(ROOT.kOrange-2)
	for j in range(0,4):
		all_hist[i][j].SetLineWidth(2)
	

right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
left,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadBottomMargin()
pCMS1 = ROOT.TPaveText(left,1.-top,0.4,1.,"NDC")
pCMS1.SetTextFont(62)
pCMS1.SetTextSize(top*0.75)
pCMS1.SetTextAlign(12)
pCMS1.SetFillStyle(-1)
pCMS1.SetBorderSize(0)
pCMS1.AddText("CMS")

pCMS12 = ROOT.TPaveText(left+0.1,1.-top*1.13,0.57,1.,"NDC")
pCMS12.SetTextFont(52)
pCMS12.SetTextSize(top*0.73)
pCMS12.SetTextAlign(12)
pCMS12.SetFillStyle(-1)
pCMS12.SetBorderSize(0)
pCMS12.AddText("Simulation")



pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("(13 TeV)")

list_qq = []
list_qq.append(hist_mqq)
list_qq.append(hist_etaqq)
list_qq.append(hist_etabb)
list_qq.append(hist_etaqqbb)
list_qq.append(hist_phiqq)
list_qq.append(hist_phibb)
list_qq.append(hist_ptqqbb)
list_qq.append(hist_ptqq)
list_qq.append(hist_ptbb)
list_qq.append(hist_pzqqbb)
list_qq.append(hist_pzqq)
list_qq.append(hist_pzbb)



for num in range(0,hist_num):
	c = ROOT.TCanvas("c","c",900,900)
	list_qq[num].Scale(1./list_qq[num].Integral())
	list_qq[num].GetYaxis().SetTitle("Events")
	list_qq[num].GetXaxis().SetTitle(list_qq[num].GetTitle())
	list_qq[num].SetStats(0)
	list_qq[num].SetLineWidth(2)
	list_qq[num].GetYaxis().SetRangeUser(0.,list_qq[num].GetMaximum()*1.2)
	list_qq[num].Draw("HIST")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	pMeanRMS = ROOT.TPaveText(0.7,0.79,0.85,0.9,"NDC")
	pMeanRMS.SetTextFont(42)
	pMeanRMS.SetTextSize(top*0.5)
	pMeanRMS.SetTextAlign(11)
	pMeanRMS.SetFillStyle(-1)
	pMeanRMS.SetBorderSize(0)
	pMeanRMS.AddText("Mean = %2.2f"%(list_qq[num].GetMean()))
	pMeanRMS.AddText("RMS = %2.2f"%(list_qq[num].GetRMS() ))
	pMeanRMS.Draw()
#	c.SaveAs("plots/plot_%s.pdf" %list_qq[num].GetName() )



for num in range(0,6):
	c = ROOT.TCanvas("c","c",900,900)
	for num_in in range(0,4):
		all_hist[num][num_in].Scale(1./all_hist[num][num_in].Integral())
	all_hist[num][0].GetYaxis().SetTitle("Events")
	if (num % 2 == 0) : all_hist[num][0].GetXaxis().SetTitle("p_{T} (GeV)")
	if (num % 2 == 1) : all_hist[num][0].GetXaxis().SetTitle("#eta")
	
	all_hist[num][0].SetStats(0)
	if (num % 2 == 0) : ymax = 0.12
	if (num % 2 == 1) : ymax = 0.25
	if (num % 4 == 0) : ymax = 0.08
	if (num == 0) : ymax = 0.12
	all_hist[num][0].GetYaxis().SetRangeUser(0.,ymax)
	all_hist[num][0].Draw("HIST")
	for num_in in range(1,4):
		all_hist[num][num_in].Draw("HISTsame")
		all_hist[num][num_in].SetStats(0)
	

	leg = ROOT.TLegend(0.35,0.73,0.65,0.9)
	for num_in in range(0,4):
		leg.AddEntry(all_hist[num][num_in],"Parton[%i], Mean = %2.2f, RMS = %2.2f"%(num_in, all_hist[num][num_in].GetMean(),all_hist[num][num_in].GetRMS() ) ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()


	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()


hist_fout = ROOT.TFile('hist_sig.root','recreate')
hist_mqq_pt.Write()
hist_mqq_deta.Write()
#	c.SaveAs("plots/plot_%s.pdf" %all_hist[num][0].GetName() )



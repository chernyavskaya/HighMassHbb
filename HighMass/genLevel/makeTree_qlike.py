import ROOT
from ROOT import gROOT
from ROOT import gStyle
from array import array

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

f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/workflow/VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8_v21.root")
tree = f.Get("tree")

hist_count = f.Get('CountWeighted')
events_generated = hist_count.GetEntries()


qlike_fout = ROOT.TFile('qlikelihood_intput_tree.root','recreate')

tree_qlike = ROOT.TTree('Jet_tree_q','Jet_tree_q')
qlike_jet_pt = array('f',[0])
qlike_jet_pt_idx = array('i',[0])
qlike_jet_eta = array('f',[0])
qlike_jet_eta_idx = array('i',[0])
q_matched = array('i',[0])
tree_qlike.Branch('Jet_pt',qlike_jet_pt,'qlike_jet_pt[1]/F')
tree_qlike.Branch('Jet_pt_idx',qlike_jet_pt_idx,'qlike_jet_pt_idx[1]/I')
tree_qlike.Branch('Jet_eta',qlike_jet_eta,'qlike_jet_eta[1]/F')
tree_qlike.Branch('Jet_eta_idx',qlike_jet_eta_idx,'qlike_jet_eta_idx[1]/I')
tree_qlike.Branch('q_matched',q_matched,'q_matched[1]/I')


presel=0
presel_trigger=0

eff_choosing_jets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_bjets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_qjets_pt = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
eff_choosing_jets_csv_pt = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_jets_csv_times_pt_pt = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_jets_csv_times_pt_deta = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_jets_csv_deta = {'b1':0,'b2':0,'q1':0,'q2':0}
eff_choosing_ref = 0

num_jets_considered = 5
for entry in range(1,tree.GetEntries()):
#for entry in range(1,100):
	if entry%10000==0 : print entry
	tree.GetEntry(entry)

	weight = abs(tree.genWeight)/tree.genWeight * tree.bTagWeight * tree.puWeight


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
		if (tree.Jet_puId[i] > 0) and (tree.Jet_id[i] > 2) and (tree.Jet_pt[i]>20.):
			good_jets.append(ROOT.TLorentzVector())
			good_jets[num_good_jets].SetPtEtaPhiM(tree.Jet_pt[i],tree.Jet_eta[i],tree.Jet_phi[i],tree.Jet_mass[i])
			num_good_jets+=1
			good_jets_array_idx.append(i)
	
	if len(good_jets) >= 4 :
		if (good_jets[0].Pt() > 92) and (good_jets[1].Pt() >76 ) and (good_jets[2].Pt() > 64) and (good_jets[3].Pt() >30 ) :
			btag = []
			btag_pt = []
			pt = []
			for i in range(0,len(good_jets)) :
				btag.append(tree.Jet_btagCSV[good_jets_array_idx[i]])
				btag_pt.append(tree.Jet_btagCSV[good_jets_array_idx[i]]*good_jets[i].Pt())
				pt.append(tree.Jet_pt[good_jets_array_idx[i]])
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
			chosen_qjet1_csvpt = jets_for_csvpt_pt[0]
			chosen_qjet2_csvpt = jets_for_csvpt_pt[1]

###########################################################
###########################################################
###########################################################
			
			eta_qlike_jets = []
			for jet in jets_for_csvpt_pt:
				eta_qlike_jets.append(abs(jet.Eta()))
	
			eta_idx_qlike_jets = sorted(range(min(len(eta_qlike_jets),4)), key=lambda k: eta_qlike_jets[k])	
		

			for i in range(0,min(len(jets_for_csvpt_pt),4)) :
				q_matched[0] = 0	
				qlike_jet_pt[0] = jets_for_csvpt_pt[i].Pt()
				qlike_jet_pt_idx[0] = i
				qlike_jet_eta[0] = abs(jets_for_csvpt_pt[i].Eta())
				qlike_jet_eta_idx[0] = eta_idx_qlike_jets[i]
				if (jets_for_csvpt_pt[i].DeltaR(jets[2])<0.8) or (jets_for_csvpt_pt[i].DeltaR(jets[3])<0.8) : q_matched[0] = 1
				tree_qlike.Fill()
			
		
###########################################################
###########################################################
###########################################################
					

	# matching efficiency part
	
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
					for j in range(i + 1, min(num_jets_considered,len(jets_for_csvpt_pt))):
						d = abs(jets_for_csvpt_pt[i].Eta() - jets_for_csvpt_pt[j].Eta())
						if d > maxDeltaEta:	
							maxEta1, maxEta2, maxDeltaEta = i, j, d
				if (maxEta1>=0) and (maxEta2>=0) :
					chosen_qjet1_csvpt_deta = jets_for_csvpt_pt[maxEta1]		
					chosen_qjet2_csvpt_deta = jets_for_csvpt_pt[maxEta2]

				if (chosen_qjet1_deta.DeltaR(jets[2])<0.8) or (chosen_qjet1_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_deta['q1']+=weight
				if (chosen_qjet2_deta.DeltaR(jets[2])<0.8) or (chosen_qjet2_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_deta['q2']+=weight

				if (chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_pt['b1']+=weight
				if (chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_pt['b2']+=weight
				if (chosen_qjet1_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_pt['q1']+=weight
				if (chosen_qjet2_csvpt.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_pt['q2']+=weight
				
				if (chosen_bjet1_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet1_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_deta['b1']+=weight
				if (chosen_bjet2_csvpt.DeltaR(jets[0])<0.8) or (chosen_bjet2_csvpt.DeltaR(jets[1])<0.8) : eff_choosing_jets_csv_times_pt_deta['b2']+=weight
				if (chosen_qjet1_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet1_csvpt_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_deta['q1']+=weight
				if (chosen_qjet2_csvpt_deta.DeltaR(jets[2])<0.8) or (chosen_qjet2_csvpt_deta.DeltaR(jets[3])<0.8) : eff_choosing_jets_csv_times_pt_deta['q2']+=weight


				for i,jet in enumerate(good_jets):
					if (jet.DeltaR(jets[0])<0.8) or (jet.DeltaR(jets[1])<0.8) or (jet.DeltaR(jets[2])<0.8) or (jet.DeltaR(jets[3])<0.8) : eff_choosing_jets_pt[str(i+1)]+=weight
					if (jet.DeltaR(jets[0])<0.8) or (jet.DeltaR(jets[1])<0.8)  : eff_choosing_bjets_pt[str(i+1)]+=weight
					if (jet.DeltaR(jets[2])<0.8) or (jet.DeltaR(jets[3])<0.8)  : eff_choosing_qjets_pt[str(i+1)]+=weight	
					if i==6 : break
			
				
			
			
					
qlike_fout.Write()
qlike_fout.Close()



print 'matching efficiency for csv + pt for number of jets considered %i : \n',num_jets_considered
print 'b1 =  %.3f '%(eff_choosing_jets_csv_pt['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_pt['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_pt['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_pt['q2']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv + deta for number of jets considered %i : \n',num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_deta['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_deta['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_deta['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_deta['q2']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv*pt and pt for number of jets considered %i : \n',num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_times_pt_pt['q2']/eff_choosing_ref), '\n\n'
print 'matching efficiency for csv*pt and deta for number of jets considered %i : \n',num_jets_considered
print 'b1 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['b1']/eff_choosing_ref), '\n'
print 'b2 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['b2']/eff_choosing_ref), '\n'
print 'q1 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['q1']/eff_choosing_ref), '\n'
print 'q2 = %.3f  '%(eff_choosing_jets_csv_times_pt_deta['q2']/eff_choosing_ref), '\n\n'



	

import math
import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import RooRealVar
from ROOT import RooBukinPdf
from ROOT import RooDataHist
from ROOT import RooArgList
from ROOT import RooFit
from ROOT import RooFormulaVar
from ROOT import RooCBShape
from ROOT import RooBernstein
from ROOT import RooAddPdf



def getPt(item):
	return item[0]
def getEta(item):
	return item[1]



gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadLeftMargin(0.18)

#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/HighMassHbb/workflow/VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8_v21.root")
f = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/genLevel/VBFSpin0ToBBbar_W_1p0_M_750_TuneCUEP8M1_13TeV_pythia8_v21.root")
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


hist_mbb = ROOT.TH1F("mbb","m_{bb} (GeV)",150,0,1500)
hist_mbb_fsr = ROOT.TH1F("mbb_fsr","m_{bb} + FSR (GeV)",150,0,1500)
hist_mbb_fsr_smreg = ROOT.TH1F("mbb_fsr_reg_smreg","m_{bb} + FSR + SMReg (GeV)",150,0,1500)
hist_mbb_fsr_met = ROOT.TH1F("mbb_fsr_met","m_{bb} + FSR + MET (GeV)",150,0,1500)
hist_mbb.Sumw2(ROOT.kTRUE)
hist_mbb_fsr.Sumw2(ROOT.kTRUE)
hist_mbb_fsr_met.Sumw2(ROOT.kTRUE)
hist_mbb_fsr_smreg.Sumw2(ROOT.kTRUE)

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

num_jets_considered = 5
for entry in range(1,tree.GetEntries()):
#for entry in range(1,10000):
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
		if (tree.Jet_puId[i] > 0) and (tree.Jet_id[i] > 2) :
			good_jets.append(ROOT.TLorentzVector())
			good_jets[num_good_jets].SetPtEtaPhiM(tree.Jet_pt[i],tree.Jet_eta[i],tree.Jet_phi[i],tree.Jet_mass[i])
			num_good_jets+=1
			good_jets_array_idx.append(i)
	
	if len(good_jets) >= 4 :
		if (good_jets[0].Pt() > 92) and (good_jets[1].Pt() >76 ) and (good_jets[2].Pt() > 64) and (good_jets[3].Pt() >30 ) :
			btag = []
			btag_pt = []
			pt = []
			for i in range(0,min(len(good_jets),15)) :
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
			qjets_idc_csvpt_pt = 0, 1
			all4jets_idc_csvpt_pt = 0, 1, bjet1_idx_csvpt, bjet2_idx_csvpt


		
			######################################################
			################## Mbb part ##########################
			######################################################

			bb_jet_csvpt = chosen_bjet1_csvpt + chosen_bjet2_csvpt
			hist_mbb.Fill(bb_jet_csvpt.M(),weight)
		

			chosen_bjet1_csvpt_smreg = ROOT.TLorentzVector()
			chosen_bjet2_csvpt_smreg = ROOT.TLorentzVector()
			chosen_bjet1_csvpt_smreg.SetPtEtaPhiM(tree.Jet_pt_reg[good_jets_array_idx[bjet1_idx_csvpt]], good_jets[bjet1_idx_csvpt].Eta(),good_jets[bjet1_idx_csvpt].Phi(),good_jets[bjet1_idx_csvpt].M())
			chosen_bjet2_csvpt_smreg.SetPtEtaPhiM(tree.Jet_pt_reg[good_jets_array_idx[bjet2_idx_csvpt]], good_jets[bjet2_idx_csvpt].Eta(),good_jets[bjet2_idx_csvpt].Phi(),good_jets[bjet2_idx_csvpt].M())
			bb_jet_csvpt_smreg = chosen_bjet1_csvpt_smreg + chosen_bjet2_csvpt_smreg

			good_jets_for_fsr = [i for j, i in enumerate(good_jets) if j not in all4jets_idc_csvpt_pt]
			num_fsr = 0		
			fsr_jets = []
			for jet in good_jets_for_fsr:
				if ( jet.DeltaR(chosen_bjet1_csvpt) < 0.8 ) or (jet.DeltaR(chosen_bjet2_csvpt) < 0.8) : 
					fsr_jets.append(jet)
					num_fsr+=1
			
				
			bb_jets_csvpt_fsr = bb_jet_csvpt 
			if (num_fsr > 0) :
				for jet in fsr_jets:
					bb_jets_csvpt_fsr += jet
					bb_jet_csvpt_smreg +=jet 

			hist_mbb_fsr.Fill(bb_jets_csvpt_fsr.M(), weight)
			hist_mbb_fsr_smreg.Fill(bb_jet_csvpt_smreg.M(), weight)


			bb_jets_csvpt_fsr_met = bb_jets_csvpt_fsr
			met = ROOT.TLorentzVector()
		#	met.SetPtEtaPhiM(tree.met_pt, tree.met_eta, tree.met_phi, tree.met_mass)
		#	if (met.DeltaPhi(chosen_bjet1_csvpt) < 0.1) or (met.DeltaPhi(chosen_bjet2_csvpt) < 0.1) :  #the best
		# 	if (met.DeltaPhi(chosen_bjet1_csvpt) < 0.2) or (met.DeltaPhi(chosen_bjet2_csvpt) < 0.2) :
		#	if (met.DeltaR(chosen_bjet1_csvpt) < 0.8) or (met.DeltaR(chosen_bjet2_csvpt) < 0.8) :  #also very good 
		#		bb_jets_csvpt_fsr_met += met

			met_dphi1 = abs(tree.met_phi - chosen_bjet1_csvpt.Phi())
			met_dphi2 = abs(tree.met_phi - chosen_bjet2_csvpt.Phi()) 
			if (met_dphi1 < met_dphi2):
				met.SetPtEtaPhiM(tree.met_pt * math.cos(min(met_dphi1,met_dphi2)),chosen_bjet1_csvpt.Eta(),chosen_bjet1_csvpt.Phi(),0)
			else : met.SetPtEtaPhiM(tree.met_pt * math.cos(min(met_dphi1,met_dphi2)),chosen_bjet2_csvpt.Eta(),chosen_bjet2_csvpt.Phi(),0)
			bb_jets_csvpt_fsr_met += met
			hist_mbb_fsr_met.Fill(bb_jets_csvpt_fsr_met.M(), weight)
			
			######################################################
			######################################################
			######################################################



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

	




all_hist = []
all_hist.append(hist_pt_ord_pt)
all_hist.append(hist_pt_ord_eta)
all_hist.append(hist_eta_ord_pt)
all_hist.append(hist_eta_ord_eta)
all_hist.append(hist_flavour_ord_pt)
all_hist.append(hist_flavour_ord_eta)

list_mbb = []
list_mbb.append(hist_mbb)
list_mbb.append(hist_mbb_fsr)
list_mbb.append(hist_mbb_fsr_met)
#list_mbb.append(hist_mbb_fsr_smreg)
mbb_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
for j,i in enumerate(mbb_colors):
	list_mbb[j].SetLineColor(i)

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
	list_qq[num].GetYaxis().SetRangeUser(0.,list_qq[num].GetMaximum()*1.1)
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



c = ROOT.TCanvas("c","c",900,900)
pMeanRMS = []
pMSig = []
for num in range(0,len(list_mbb)):
	list_mbb[num].Scale(1./list_mbb[num].Integral())
	list_mbb[num].SetStats(0)
	list_mbb[num].SetLineWidth(2)


	x   = RooRealVar("mbb%i"%num,"mbb%i"%num,200,1200)
	roo_mbb_hist = RooDataHist("roohist_fit_mbb_%i"%num,"roohist_fit_mbb_%i"%num,RooArgList(x),list_mbb[num])
 	Xp = RooRealVar("Xp_%i"%num,"Xp_%i"%num,750,650, 770)
 	sigp = RooRealVar("sigp_%i"%num,"sigp_%i"%num,80,50, 120)
 	xi = RooRealVar("xi_%i"%num,"xi_%i"%num,-0.25,-1, 1)
 	rho1 = RooRealVar("rho1_%i"%num,"rho1_%i"%num,0.15,-1, 1)
 	rho2 = RooRealVar("rho2_%i"%num,"rho2_%i"%num,0.22,-1, 1)
	sig = RooBukinPdf("signal_bukin_%i"%num,"signal_bukin_%i"%num,x,Xp,sigp,xi,rho1,rho2)


#	m             = RooRealVar("mean_%i"%num,"mean_%i"%num,750,620,850)
#	s             = RooRealVar("sigma_%i"%num,"sigma_%i"%num,70,50,150)
#	width         = RooRealVar("fwhm_%i"%num,"fwhm_%i"%num,100,50,200)
#	a             = RooRealVar("alpha_%i"%num,"alpha_%i"%num,1,0,10)
#	n             = RooRealVar("exp_%i"%num,"exp_%i"%num,1,0,10)
#	fsig          = RooRealVar("fsig_%i"%num,"fsig_%i"%num,0.03,0.,1.)
#	sigSig        = RooCBShape("signal_gauss_%i"%num,"signal_gauss_%i"%num,x,m,s,a,n)

#	b0,b1,b2,b3   = [RooRealVar("b%d_%i"%(num,num),"b%d_%d"%(num,num),0.5,0.,3.) for i in range(4)]
  ### Bkg part: Bernstein     
#	bkg           = RooBernstein("signal_bkg_%d"%num,"signal_bkg_%d"%num,x,RooArgList(b0,b1,b2))
#	sig         = RooAddPdf("signal_model_%d"%num,"signal_model_%d"%num,RooArgList(sigSig,bkg),RooArgList(fsig))

	res = sig.fitTo(roo_mbb_hist,RooFit.SumW2Error(ROOT.kFALSE))
	for o in [Xp,sigp,xi,rho1,rho2]: o.setConstant(ROOT.kTRUE)
#	for o in [m,s,width,a,n,fsig]: o.setConstant(ROOT.kTRUE)


	frame = x.frame()
	roo_mbb_hist.plotOn(frame,RooFit.DrawOption("Psame"), RooFit.LineWidth(2), RooFit.LineColor(mbb_colors[num]), RooFit.MarkerColor(mbb_colors[num]))
	sig.plotOn(frame,RooFit.LineColor(mbb_colors[num]))
#	bkg.plotOn(frame,RooFit.LineColor(mbb_colors[num]),RooFit.LineStyle(2))
	if (num==0): 
		frame.GetYaxis().SetRangeUser(0.,list_mbb[0].GetMaximum()*1.4)
		frame.GetXaxis().SetNdivisions(5,False)
		frame.GetYaxis().SetTitleOffset(1.5)
		frame.GetYaxis().SetTitle("Events")
		frame.GetXaxis().SetTitle(list_mbb[0].GetTitle())
		frame.Draw()
	else :  frame.Draw("same")


	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	pMeanRMS.append(ROOT.TPaveText(0.2,0.59+0.1*num,0.35,0.69+0.1*num,"NDC"))
	pMeanRMS[num].SetTextFont(42)
	pMeanRMS[num].SetTextSize(top*0.5)
	pMeanRMS[num].SetTextAlign(11)
	pMeanRMS[num].SetFillStyle(-1)
	pMeanRMS[num].SetBorderSize(0)
	pMeanRMS[num].SetTextColor(mbb_colors[num])
	pMeanRMS[num].AddText("Mean = %2.2f"%(list_mbb[num].GetMean()))
	pMeanRMS[num].AddText("RMS = %2.2f"%(list_mbb[num].GetRMS() ))
	pMeanRMS[num].Draw("same")

	pMSig.append(ROOT.TPaveText(0.65,0.59+0.1*num,0.8,0.69+0.1*num,"NDC"))
	pMSig[num].SetTextFont(42)
	pMSig[num].SetTextSize(top*0.5)
	pMSig[num].SetTextAlign(11)
	pMSig[num].SetFillStyle(-1)
	pMSig[num].SetBorderSize(0)
	pMSig[num].SetTextColor(mbb_colors[num])
	pMSig[num].AddText("M = %2.2f #pm %.2f"%(Xp.getVal(),Xp.getError()))
	pMSig[num].AddText("#sigma = %2.2f #pm %.2f"%(sigp.getVal(),sigp.getError()))
#	pMSig[num].AddText("M = %2.2f #pm %.2f"%(m.getVal(),m.getError()))
#	pMSig[num].AddText("#sigma = %2.2f #pm %.2f"%(s.getVal(),s.getError()))
	pMSig[num].Draw("same")
c.SaveAs("plots/plot_%s_fitBukin_project.pdf" %list_mbb[len(list_mbb)-1].GetName() )


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


#	c.SaveAs("plots/plot_%s.pdf" %all_hist[num][0].GetName() )



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

f = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFHToBB_M-125_13TeV_powheg_Py8_weightfix__fall15MAv2-pu25ns15v1_76r2as.root")
tree = f.Get("tree")
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

for entry in range(1,tree.GetEntries()):
#for entry in range(1,10000):
	if entry%10000==0 : print entry
	tree.GetEntry(entry)

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

	if (qq_jet.M()>450) : presel+=weight


print 'DoubleB200 = ', triggers['DoubleB200']/100000., 'presel + doubleB200 = ', presel_trigger/100000.,'presel = ', presel/100000. , 'DoubleB240 = ', triggers['DoubleB240']/100000., 'SingleB460 = ', triggers['SingleB460']/100000., 'SingleB500 = ', triggers['SingleB500']/100000., 'SingleBorDoubleB200_460 = ', triggers['SingleBorDoubleB200_460']/100000., 'SingleBorDoubleB240_500 = ', triggers['SingleBorDoubleB240_500']/100000
	
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
	c.SaveAs("plots_125/plot_%s.pdf" %list_qq[num].GetName() )



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
		all_hist[num][num_in].SetStats(0);
	

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


	c.SaveAs("plots_125/plot_%s.pdf" %all_hist[num][0].GetName() )



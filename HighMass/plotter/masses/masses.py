import ROOT
from ROOT import gROOT
from ROOT import gStyle

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]
def getD(item):
	return item[1]

def getDvalue(h1, h2):
	h1s = h1.GetNbinsX()
	h2s = h2.GetNbinsX()
	if (h1s!=h2s) :
		print "hists have different N bins"
		return -1
	if (h1.Integral()!=0) : h1.Scale(1./h1.Integral()) 
	if (h2.Integral()!=0) : h2.Scale(1./h2.Integral()) 
	adiff = 0
	for i in range(0,h1s):
		adiff+=abs(h1.GetBinContent(i) - h2.GetBinContent(i))
	return adiff/2

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.30)
gStyle.SetPadLeftMargin(0.15)


path="dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/HighMassHbb/plotterOutput/v21/"
#path="/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/plotter/"
files = [] 
files.append(ROOT.TFile.Open(path+"BTagCSV_analysis_v21_TriggerFinal.root"))
num_masses = 11
masses = ['375','450','525','600','675','750','825','900','975','1050']
for i in range(0,len(masses)):
	name = path+'VBFSpin0ToBBbar_W_1p0_M_' + masses[i] + '_analysis_v21_TriggerFinal.root'
	files.append(ROOT.TFile.Open(name))


colors = [ROOT.kBlack,ROOT.kRed, ROOT.kRed-7, ROOT.kOrange, ROOT.kOrange-1, ROOT.kYellow-6, ROOT.kSpring+10, ROOT.kGreen+1, ROOT.kCyan,ROOT.kBlue, ROOT.kViolet+1, ROOT.kMagenta  ]
names = ["BTagCSV","VBF Hbb, 375 GeV","VBF Hbb, 450 GeV","VBF Hbb, 525 GeV","VBF Hbb, 600 GeV","VBF Hbb, 675 GeV","VBF Hbb, 750 GeV","VBF Hbb, 825 GeV","VBF Hbb, 900 GeV","VBF Hbb, 975 GeV","VBF Hbb, 1050 GeV"]
#colors = [ROOT.kRed, ROOT.kRed-7, ROOT.kOrange, ROOT.kOrange-1, ROOT.kYellow-6, ROOT.kSpring+10, ROOT.kGreen+1, ROOT.kCyan,ROOT.kBlue, ROOT.kViolet+1, ROOT.kMagenta  ]
#names = ["VBF Hbb, 300 GeV","VBF Hbb, 375 GeV","VBF Hbb, 450 GeV","VBF Hbb, 525 GeV","VBF Hbb, 600 GeV","VBF Hbb, 675 GeV","VBF Hbb, 750 GeV","VBF Hbb, 825 GeV","VBF Hbb, 900 GeV","VBF Hbb, 975 GeV","VBF Hbb, 1050 GeV"]


right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
left,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadBottomMargin()
pCMS1 = ROOT.TPaveText(left,1.-top,0.4,1.,"NDC")
pCMS1.SetTextFont(62)
pCMS1.SetTextSize(top*0.75)
pCMS1.SetTextAlign(12)
pCMS1.SetFillStyle(-1)
pCMS1.SetBorderSize(0)
pCMS1.AddText("CMS")

pCMS12 = ROOT.TPaveText(left+0.1,1.-top*1.13,0.48,1.,"NDC")
pCMS12.SetTextFont(52)
pCMS12.SetTextSize(top*0.73)
pCMS12.SetTextAlign(12)
pCMS12.SetFillStyle(-1)
pCMS12.SetBorderSize(0)
pCMS12.AddText("Work in progress")

pCMS2 = ROOT.TPaveText(0.6,1.-top,0.7,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("(13 TeV)")


hist_names = ["hMqq", "hMqq_pt","hJet1_pt","hJet2_pt","hJet3_pt","hJet4_pt","hJet1_eta","hJet2_eta","hJet3_eta","hJet4_eta","hJet1_phi","hJet2_phi","hJet3_phi","hJet4_phi", "hEtaQQ", "hPhiBB","hMbb","hbtag","hbtag2","hcosOqqbb","hEtaQB1", "hEtaQB2", "hPhiQB1", "hPhiQB2","hx1","hx2","hVB1_mass","hVB2_mass","hEtot","hPxtot","hPytot","hPztot","hJet5_pt","hPtqqbb","hEtaqqbb","hPhiqqbb","hJet1_pt_bin","hJet2_pt_bin","hJet3_pt_bin","hJet4_pt_bin", "hMqq_bin","hEtaSoftJets", "hPtSoftJets","hMassSoftJets","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hqgl","hqgl2", "hPtSoftJets2","hPtSoftJets3","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult","hMbb","hMbb_met_fsr","hblike1","hblike2", "hmet", "hmet", "hmet", "hmet", "hPhiQQ","hselLeptons_relIso03", "hbtag_log","hbtag2_log","hrho","hJet1q_leadTrackPt","hJet2q_leadTrackPt" ,"hJet1b_pt", "hJet2b_pt", "hJet1b_eta", "hJet2b_eta", "hJetbb_pt", "hJetqqbb_pz", "hEtaQBplus", "hEtaQBminus", "hbdt"]


c = ROOT.TCanvas("c","c",1200,900)
for num in range(0,len(hist_names)):
	leg = ROOT.TLegend(0.72,0.5,0.92,0.9)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.027)
	hists = []
	for n_files,o in enumerate(files):
		hists.append(o.Get(hist_names[num]))
		hists[n_files].SetLineColor(colors[n_files])
		hists[n_files].SetLineWidth(2)
		hists[0].GetYaxis().SetTitle("Events")
		hists[0].GetXaxis().SetTitle(hists[0].GetXaxis().GetTitle())
		hists[0].SetStats(0)
		hists[n_files].Scale(1./hists[n_files].Integral())
		ymax =  hists[0].GetBinContent(hists[0].GetMaximumBin())*1.5
		hists[0].GetYaxis().SetRangeUser(0,ymax)
		if n_files==0 : hists[n_files].Draw("HIST")
		else : hists[n_files].Draw("HISTsame")
		pCMS1.Draw()
		pCMS12.Draw()
		pCMS2.Draw()
		leg.AddEntry(hists[n_files],names[n_files],"L")
	leg.Draw()
	c.SaveAs("plots_new/plot_%s.pdf" %hists[n_files].GetName() )
	c.SaveAs("plots_new/plot_%s.png" %hists[n_files].GetName() )


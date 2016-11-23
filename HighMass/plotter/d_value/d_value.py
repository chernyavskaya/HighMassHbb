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
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

f_sig = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/plotter/VBFSpin0ToBBbar_W_1p0_M_750_analysis_v21_presel.root")
f_data = ROOT.TFile.Open('/afs/cern.ch/work/n/nchernya/HighMassHbb/HighMass/plotter/BTagCSV_analysis_v21_presel.root')

hists_data  = []
hists_sig  = []
hist_names = ["hMqq", "hMqq_pt","hJet1_pt","hJet2_pt","hJet3_pt","hJet4_pt","hJet1_eta","hJet2_eta","hJet3_eta","hJet4_eta","hJet1_phi","hJet2_phi","hJet3_phi","hJet4_phi", "hEtaQQ", "hPhiBB","hMbb","hbtag","hbtag2","hcosOqqbb","hEtaQB1", "hEtaQB2", "hPhiQB1", "hPhiQB2","hx1","hx2","hVB1_mass","hVB2_mass","hEtot","hPxtot","hPytot","hPztot","hJet5_pt","hPtqqbb","hEtaqqbb","hPhiqqbb","hJet1_pt_bin","hJet2_pt_bin","hJet3_pt_bin","hJet4_pt_bin", "hMqq_bin","hEtaSoftJets", "hPtSoftJets","hMassSoftJets","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hqgl","hqgl2", "hPtSoftJets2","hPtSoftJets3","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult","hMbb_regVBF","hMbb_regVBF_fsr","hblike1","hblike2", "hmet", "hmet", "hmet", "hmet", "hPhiQQ","hselLeptons_relIso03", "hbtag_log","hbtag2_log","hrho","hJet1q_leadTrackPt","hJet2q_leadTrackPt" ,"hJet1b_pt", "hJet2b_pt", "hJet1b_eta", "hJet2b_eta", "hJetbb_pt", "hJetqqbb_pz", "hEtaQBplus", "hEtaQBminus"]

for o in hist_names:
	hists_data.append(f_data.Get(o))
	hists_sig.append(f_sig.Get(o))




for i in range(0,len(hist_names)):
	hists_data[i].SetLineColor(ROOT.kBlack)
	hists_data[i].SetLineWidth(2)
	hists_sig[i].SetLineWidth(2)
	hists_sig[i].SetLineColor(ROOT.kRed)
	

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
pCMS12.AddText("Work in progress")



pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("(13 TeV)")


d_values = []

for num in range(0,len(hist_names)):
	c = ROOT.TCanvas("c","c",900,900)
	hists_sig[num].GetYaxis().SetTitle("Events")
	hists_sig[num].GetXaxis().SetTitle(hists_data[num].GetXaxis().GetTitle())
	hists_sig[num].SetStats(0)
	hists_sig[num].SetLineWidth(2)
	hists_sig[num].Scale(1./hists_sig[num].Integral())
	hists_data[num].Scale(1./hists_data[num].Integral())
	ymax =  max(hists_sig[num].GetBinContent(hists_sig[num].GetMaximumBin()), hists_data[num].GetBinContent(hists_data[num].GetMaximumBin()))*1.2
	hists_sig[num].GetYaxis().SetRangeUser(0,ymax)
	hists_sig[num].Draw("HIST")
	hists_data[num].Draw("HISTsame")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	dValue = ROOT.TPaveText(0.6,0.79,0.85,0.9,"NDC")
	dValue.SetTextFont(42)
	dValue.SetTextSize(top*0.5)
	dValue.SetTextAlign(11)
	dValue.SetFillStyle(-1)
	dValue.SetBorderSize(0)
	h1 = hists_data[num].Clone("newdata")
	h2 = hists_sig[num].Clone("newsig")
	dval = getDvalue(h1,h2)
	dValue.AddText("d = %2.2f"%(dval))
	dValue.Draw()
	leg = ROOT.TLegend(0.6,0.7,0.85,0.83)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.AddEntry(hists_sig[num],"VBF Hbb, 750 GeV","L")
	leg.AddEntry(hists_data[num],"BTagCSV","L")
	leg.Draw()
	d_values.append( (hists_sig[num].GetName(),dval)  )
	c.SaveAs("plots/plot_%s.pdf" %hists_sig[num].GetName() )
	c.SaveAs("plots/plot_%s.png" %hists_sig[num].GetName() )

d_values_sorted = sorted(d_values,key=getD, reverse=True)
for o in d_values_sorted:
	print o, "\n" 

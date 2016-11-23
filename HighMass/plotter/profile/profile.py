import ROOT
from ROOT import gROOT
from ROOT import gStyle

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]
def getD(item):
	return item[1]


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
hist_names = ["hprof_etot", "hprof_mqq", "hprof_etaqq", "hprof_qqbb_pz", "hprof_qqbb_pt", "hprof_x1", "hprof_x2", "hprof_btag_log1", "hprof_btag_log2", "hprof_etaqb", "hprof_softn2", "hprof_phiqq", "hprof_axis2_1", "hprof_axis2_2"]

for o in hist_names:
	hists_data.append(f_data.Get(o))
#	hists_sig.append(f_sig.Get(o))




for i in range(0,len(hist_names)):
	hists_data[i].SetLineColor(ROOT.kBlack)
	hists_data[i].SetLineWidth(2)
#	hists_sig[i].SetLineWidth(2)
#	hists_sig[i].SetLineColor(ROOT.kRed)
	

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
	hists_data[num].GetYaxis().SetTitle("M_{bb}")
	hists_data[num].GetXaxis().SetTitle(hists_data[num].GetTitle())
	hists_data[num].SetStats(0)
	hists_data[num].SetLineWidth(2)
	ymax = hists_data[num].GetBinContent(hists_data[num].GetMaximumBin())*1.2
	hists_data[num].GetYaxis().SetRangeUser(0,ymax)
#	hists_data[num].Draw("HIST")
	hists_data[num].Draw("PE")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	leg = ROOT.TLegend(0.6,0.85,0.85,0.9)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
#	leg.AddEntry(hists_data[num],"VBF Hbb, 750 GeV","L")
	leg.AddEntry(hists_data[num],"BTagCSV","L")
	leg.Draw()
	c.SaveAs("plots/plot_%s.pdf" %hists_data[num].GetName() )
	c.SaveAs("plots/plot_%s.png" %hists_data[num].GetName() )


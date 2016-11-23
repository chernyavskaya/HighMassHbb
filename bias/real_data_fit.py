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

#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/HighMassHbb/plotterOutput/v21/BTagCSV_analysis_v21_TrigWtFinal.root")
f = ROOT.TFile.Open("BTagCSV_analysis_v21_TrigWtFinal.root")
data_hist = f.Get("hMbb_met_fsr")
data_hist.SetLineWidth(2)
	

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



sqrts=1.3e+04

canvas = ROOT.TCanvas("canvas","canvas",900,900)
x   = RooRealVar("mbb","mbb",310,1200)
roo_mbb_hist = RooDataHist("roohist_fit_mbb","roohist_fit_mbb",RooArgList(x),data_hist)

#func = ROOT.TF1("func","[0]*exp(-[1]*log(x)-[1]*[2]*log(x)*log(x))",310.,1400.)
#func.SetParameters(10,0.1,0.5)
#data_hist.Draw()
#data_hist.Fit(func,"R")


pdf = None
coeff = ROOT.RooArgList()
pdf_name = "dijet"
n_param = 2
#formula = ("TMath::Exp(-@0*TMath::Log(x/%E))" % (sqrts))
#formula = ("TMath::Exp(-@0*TMath::Log(x))")
#formula = ("TMath::Exp(-@0*TMath::Log(x/%E)-@0*@1*TMath::Log(x/%E)*TMath::Log(x/%E))" % (sqrts, sqrts, sqrts))

formula = "TMath::Exp(-@0*TMath::Log(@2)-@0*@1*TMath::Log(@2)*TMath::Log(@2))*(1+@3*(@2))"

#c  = RooRealVar("c", "", -5,0);
#c2 = RooRealVar("c2", "", -5,0);
#c3 = RooRealVar("c3", "", -1,1);
c  = RooRealVar("c", "", 0,10);
c2 = RooRealVar("c2", "", 0,10);
c3 = RooRealVar("c3", "", 0,10);

#for p in range(0,n_param):
#	p_name = ("a%d" % (p))
#	p_min = -5
#	p_max = 0	
#	param = ROOT.RooRealVar( p_name, "", p_min, p_max)
#	gcs.append(param)
#	coeff.add(param)
#	print 'her'
print formula
#coeff.add(x)
#coeff.Print()
pdf = ROOT.RooGenericPdf( ("%s" % (pdf_name)), "", "TMath::Exp(-@0*TMath::Log(@2/13000)-@0*@1*TMath::Log(@2/13000)*TMath::Log(@2/13000))*(1+@3*(@2/13000))", RooArgList(c,c2,x,c3) )
#pdf = ROOT.RooGenericPdf( ("%s" % (pdf_name)), "", "TMath::Exp(-@0*TMath::Log(@2)-@0*@1*TMath::Log(@2)*TMath::Log(@2))", RooArgList(c,c2,x) )
res = pdf.fitTo(roo_mbb_hist, RooFit.SumW2Error(ROOT.kFALSE), 
RooFit.Warnings(ROOT.kFALSE))



frame = x.frame()
roo_mbb_hist.plotOn(frame,RooFit.DrawOption("Psame"), RooFit.LineWidth(2))
pdf.plotOn(frame,RooFit.LineColor(2))
frame.GetYaxis().SetRangeUser(0.,data_hist.GetMaximum()*1.1)
frame.GetXaxis().SetNdivisions(5,False)
frame.GetYaxis().SetTitleOffset(1.5)
frame.GetYaxis().SetTitle("Events")
frame.GetXaxis().SetTitle(data_hist.GetTitle())
chi2 = frame.chiSquare( )
print chi2
frame.Draw()
print 'res minNll  =  ', res.minNll()



pCMS1.Draw()
pCMS12.Draw()
pCMS2.Draw()
canvas.SaveAs("plots/plot_test.pdf" )


